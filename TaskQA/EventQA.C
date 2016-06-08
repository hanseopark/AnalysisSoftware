#include "QA.h"

void EventQA(
                Int_t nSetsIn,
                TString fEnergyFlag,
                TString* DataSets,
                TString* plotDataSets,
                TString* pathDataSets,
                Int_t mode              = 2,            // standard mode for analysis
                Int_t cutNr             = -1,           // if -1: you have to choose number at runtime
                Int_t doExtQA           = 2,            // 0: switched off, 1: normal extQA, 2: with Cell level plots
                TString suffix          = "eps",        // output format of plots
                TString labelData       = "Data",       // Label for data    
                Bool_t addSubfolder     = kFALSE        // flag to enable subdirectory creation for primary cut
            )
{
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "EventQA" << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

    gROOT->Reset();
    TH1::AddDirectory(kFALSE);

    const Int_t nSets = nSetsIn;

    StyleSettingsThesis();
    SetPlotStyle();
    
    //**************************************************************************************************************
    TString fDate = ReturnDateString();
    TString fTextMeasurement = Form("#pi^{0} #rightarrow #gamma#gamma");
    TString fTextMeasurementEta = Form("#eta #rightarrow #gamma#gamma");
    TString fTextMeasurementMeson[2]={fTextMeasurement,fTextMeasurementEta};

    const Int_t maxSets = 12;
    //nSets == 0 is always data!

    if(nSets>maxSets){
        cout << "Maximum hardcoded number of Data Sets in EventQA: " << maxSets << endl;
        cout << "You have chosen: " << nSets << ", returning!" << endl;
        return;
    }

    Color_t colorCompare[maxSets] = {kBlack, kRed+1, kMagenta+2, 807, 800, kGreen+2, kCyan+2, kBlue+1, kOrange+2, kAzure, kViolet, kGray+1};
    TString nameMainDir[maxSets];

    Int_t fMode             = mode;
    Bool_t isCalo           = kFALSE;
    Bool_t isMergedCalo     = kFALSE;
    Bool_t isConv           = kFALSE;
    // mode:	0 // new output PCM-PCM
    //			1 // new output PCM dalitz
    //			2 // new output PCM-EMCal
    //			3 // new output PCM-PHOS
    //          4 // new output EMCal-EMCal
    //          5 // new output PHOS-PHOS
    //			9 // old output PCM-PCM
    //			10 // merged EMCal
    //			11 // merged PHOS
    if (fMode == 0 || fMode == 1 || fMode == 2 || fMode == 3 || fMode == 9) 
        isConv         = kTRUE;
    if (fMode == 2 || fMode == 3 || fMode == 4 || fMode == 5 || fMode == 10 || fMode == 11) 
        isCalo          = kTRUE;
    if (fMode == 10 || fMode == 11)
        isMergedCalo    = kTRUE;

    TString* fCutSelection          = new TString[nSets];
    TString* fEventCutSelection     = new TString[nSets];
    TString* fGammaCutSelection     = new TString[nSets];
    TString* fClusterCutSelection   = new TString[nSets];
    TString* fElectronCutSelection  = new TString[nSets];
    TString* fMesonCutSelection     = new TString[nSets];

    //*****************************************************************************************************
    //*****************************************************************************************************
    //****************************** Determine which cut to process ***************************************
    //*****************************************************************************************************
    //*****************************************************************************************************

    vector <TString> cutsTemp;
    map<TString,Int_t> mapCuts;

    for(Int_t i=0; i<nSets; i++){
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
        fCutSelection[i] = cuts.at(cutNr);
        fEventCutSelection[i] = ""; fGammaCutSelection[i] = ""; fClusterCutSelection[i] = ""; fElectronCutSelection[i] = ""; fMesonCutSelection[i] = "";
        ReturnSeparatedCutNumberAdvanced(fCutSelection[i], fEventCutSelection[i], fGammaCutSelection[i], fClusterCutSelection[i], fElectronCutSelection[i], fMesonCutSelection[i], fMode);
    }
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

    cout << "Obtaining trigger - ";
    TString* fTrigger = new TString[nSets];
    for(Int_t iT=0;iT<nSets;iT++){fTrigger[iT]="";}
    TString fTriggerCut = fEventCutSelection[0](3,2);
    fTrigger[0] = AnalyseSpecialTriggerCut(fTriggerCut.Atoi(), DataSets[0].Data());
    cout  << "'" << fTrigger[0].Data() << "' - was found!" << endl;
    if(fTrigger[0].Contains("not defined")){
        fTrigger[0] = "";
        cout << "INFO: Trigger cut not defined!" << endl;
    }

    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

    TString calo = "";
    TString fClusters = "";
    Int_t iCalo = 0;
    Int_t nCaloModules = 0;
    Int_t nCaloCells = 0;

    if(isCalo){
        if(fClusterCutSelection[0].BeginsWith('1')){
            calo="EMCal"; iCalo=1;
            fClusters = Form("%s clusters", calo.Data());
            nCaloModules = 10;
            nCaloCells = 11520;
            if(DataSets[0].Contains("LHC10")){
                nCaloModules = 4;
                nCaloCells = 4608;
            }
        } else if(fClusterCutSelection[0].BeginsWith('2')){
            calo="PHOS"; iCalo=2;
            fClusters = Form("%s clusters", calo.Data());
            nCaloModules = 5;
            nCaloCells = 6000;
            if(DataSets[0].Contains("LHC10")){
                nCaloModules = 5;
                nCaloCells = 6000;
            }
        } else {cout << "No correct calorimeter type found: " << calo.Data() << ", returning..." << endl; return;}
    }

    TString fCollisionSystem = ReturnFullCollisionsSystem(fEnergyFlag);
    if (fCollisionSystem.CompareTo("") == 0){
        cout << "No correct collision system specification, has been given" << endl;
        return;
    }

    TString fDetectionProcess = ReturnFullTextReconstructionProcess(fMode);

    TString outputDir = Form("%s/%s/EventQA/%s",cuts.at(cutNr).Data(),fEnergyFlag.Data(),suffix.Data());
    if(addSubfolder) outputDir+=Form("/%s",DataSets[0].Data());

    gSystem->Exec("mkdir -p "+outputDir);
    gSystem->Exec("mkdir -p "+outputDir+"/Comparison");
    gSystem->Exec("mkdir -p "+outputDir+"/Comparison/Ratios");

    //*****************************************************************************************************
    //*****************************************************************************************************
    //*****************************************************************************************************

    fstream fLog;
    fLog.open(Form("%s/A-%s.log",outputDir.Data(),DataSets[0].Data()), ios::out);
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

    std::vector<TH1D*> vecNEvents;
    std::vector<TH1D*> vecVertexZ;
    std::vector<TH1D*> vecNGoodTracks;
    std::vector<TH1D*> vecGammaCandidates;
    std::vector<TH1D*> vecV0Mult;

    std::vector<TH1D*> vecInvMassBeforeAfter[2];

    std::vector<TH1D*> signalPi0;
    std::vector<TH1D*> signalEta;

    std::vector<TH1D*> vecPi0Pt;
    std::vector<TH1D*> vecPi0Alpha;
    std::vector<TH1D*> vecPi0Y;
    std::vector<TH1D*> vecPi0OpenAngle;

    std::vector<TH1D*> vecEtaPt;
    std::vector<TH1D*> vecEtaAlpha;
    std::vector<TH1D*> vecEtaY;
    std::vector<TH1D*> vecEtaOpenAngle;

    std::vector<TH2D*> vecESDMother;
    std::vector<TH2D*> vecESDBackground;

    Double_t* nEventsAll = new Double_t[nSets];
    Double_t* nEvents = new Double_t[nSets];

    Int_t minB = 0; Int_t maxB = 0;
    Int_t minYB = 0; Int_t maxYB = 0;

    Int_t minB_SPD = 0; Int_t maxB_SPD = 0;
    Int_t minYB_SPD = 0; Int_t maxYB_SPD = 0;
    Bool_t isMinMaxSPD = kTRUE;

    MesonFit fitter;

    //*****************************************************************************************************
    //*****************************************************************************************************
    //****************************** Looping over DataSets ************************************************
    //*****************************************************************************************************
    //*****************************************************************************************************

    TCanvas* canvas = new TCanvas("canvas","",10,10,750,500);  // gives the page size
    TCanvas* cvsQuadratic = new TCanvas("cvsQuadratic","",10,10,500,500);  // gives the page size
    Double_t leftMargin = 0.09; Double_t rightMargin = 0.02; Double_t topMargin = 0.04; Double_t bottomMargin = 0.09;
    DrawGammaCanvasSettings(canvas,leftMargin,rightMargin,topMargin,bottomMargin);

    for(Int_t i=0; i<nSets; i++) {
        TFile *fFile = new TFile(pathDataSets[i].Data(),"READ");
        if(fFile->IsZombie()) {
            cout << "ERROR: ROOT file '" << pathDataSets[i].Data() << "' could not be openend, returning!" << endl;
            fLog << "ERROR: ROOT file '" << pathDataSets[i].Data() << "' could not be openend, returning!" << endl;
            return;
        }
        //-------------------------------------------------------------------------------------------------------------------------------
        TList* TopDir = (TList*) fFile->Get(nameMainDir[i].Data());
            if(TopDir == NULL) {cout << "ERROR: TopDir not Found"<<endl; return;}
            else TopDir->SetOwner(kTRUE);
        TList* TopContainer = (TList*) TopDir->FindObject(Form("Cut Number %s",fCutSelection[i].Data()));
            if(TopContainer == NULL) {cout << "ERROR: " << Form("Cut Number %s",fCutSelection[i].Data()) << " not found in File" << endl; return;}
            else TopContainer->SetOwner(kTRUE);
        TList* ESDContainer = (TList*) TopContainer->FindObject(Form("%s ESD histograms",fCutSelection[i].Data()));
            if(ESDContainer == NULL) {cout << "ERROR: " << Form("%s ESD histograms",fCutSelection[i].Data()) << " not found in File" << endl; return;}
            else ESDContainer->SetOwner(kTRUE);
        TList* ConvEventCutsContainer = (TList*) TopContainer->FindObject(Form("ConvEventCuts_%s",fEventCutSelection[i].Data()));
            if(ConvEventCutsContainer == NULL) {cout << "ERROR: " << Form("ConvEventCuts_%s",fEventCutSelection[i].Data()) << " not found in File" << endl; return;}
            else if(ConvEventCutsContainer) ConvEventCutsContainer->SetOwner(kTRUE);
        TList* ConvCutsContainer = (TList*) TopContainer->FindObject(Form("ConvCuts_%s",fGammaCutSelection[i].Data()));
            if(isConv && ConvCutsContainer == NULL) {cout << "ERROR: " << Form("ConvCuts_%s",fGammaCutSelection[i].Data()) << " not found in File" << endl; return;}
            else if(ConvCutsContainer) ConvCutsContainer->SetOwner(kTRUE);
        TList* CaloCutsContainer = (TList*) TopContainer->FindObject(Form("CaloCuts_%s",fClusterCutSelection[i].Data()));
            if(isCalo && CaloCutsContainer == NULL) {cout << "ERROR: " << Form("CaloCuts_%s",fClusterCutSelection[i].Data()) << " not found in File" << endl; return;}
            else if(CaloCutsContainer) CaloCutsContainer->SetOwner(kTRUE);
        TList* MesonCutsContainer = (TList*) TopContainer->FindObject(Form("ConvMesonCuts_%s",fMesonCutSelection[i].Data()));
            if(MesonCutsContainer == NULL) {cout << "ERROR: " << Form("ConvMesonCuts_%s",fMesonCutSelection[i].Data()) << " not found in File" << endl; return;}
            else if(MesonCutsContainer) MesonCutsContainer->SetOwner(kTRUE);
        TList* TrueContainer = (TList*) TopContainer->FindObject(Form("%s True histograms",fCutSelection[i].Data()));
            if(TrueContainer == NULL) {cout << "INFO: " << Form("%s True histograms",fCutSelection[i].Data()) << " not found in File, processing data?" << endl;}
            else TrueContainer->SetOwner(kTRUE);

        TList* TopContainerGamma = NULL;
        TString ContainerGammaCut = "";
            if(isConv){
                TList *listCuts=NULL; TString name=""; Int_t j=0;
                do{ listCuts = (TList*)TopDir->At(j); name = listCuts->GetName(); } while(!name.BeginsWith("ConvCuts_") && ++j<TopDir->GetSize());
                if(j>=TopDir->GetSize()) TopContainerGamma = NULL;
                else TopContainerGamma = listCuts;

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
        nEvents[i] = 0;
        nEventsAll[i] = 0;
        TH1D* fHistNEvents = (TH1D*)ESDContainer->FindObject("NEvents");
        //if(plotDataSets[i].Contains("JetJet") || plotDataSets[i].Contains("jetjet")) fHistNEvents = (TH1D*)ESDContainer->FindObject("NEventsWOWeight");
        if(fHistNEvents){
            nEvents[i] = (Double_t) GetNEvents(fHistNEvents,kFALSE);
            nEventsAll[i] = fHistNEvents->GetEntries() - fHistNEvents->GetBinContent(4);
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
        const char* nameOutput = Form("%s/%s/EventQA/EventQA_%s.root",fCutSelection[i].Data(),fEnergyFlag.Data(),DataSets[i].Data());
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
        //----------------------------------------------- Event properties --------------------------------------------------------------
        //-------------------------------------------------------------------------------------------------------------------------------
        // Event counting histos
        if(fHistNEvents){
            AdjustHistRange(fHistNEvents,1,10,kTRUE,1,0.1);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                fHistNEvents,"","","# of Entries",1,1,
                                0.82,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHistNEvents, Form("%s/NEvents_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
            vecNEvents.push_back(new TH1D(*fHistNEvents));
        }
        //-------------------------------------------------------------------------------------------------------------------------------
        // vertex Z distribution
        TH1D* fHistVertexZ = (TH1D*)ESDContainer->FindObject("VertexZ");
        if(fHistVertexZ){
            GetMinMaxBin(fHistVertexZ,minB,maxB);
            SetXRange(fHistVertexZ,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kFALSE,kFALSE,
                                fHistVertexZ,"","z-Vertex (cm)","#frac{dN}{dz}",1,1,
                                0.82,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            WriteHistogram(fHistVertexZ);
            vecVertexZ.push_back(new TH1D(*fHistVertexZ));
        } else cout << "INFO: Object |fHistVertexZ| could not be found! Skipping Draw..." << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        // number of good tracks in central acceptance (|\eta| < 0.8)
        TH1D* fHistNGoodTracks = (TH1D*)ESDContainer->FindObject("GoodESDTracks");
        if(fHistNGoodTracks){
            GetMinMaxBin(fHistNGoodTracks,minB,maxB);
            SetXRange(fHistNGoodTracks,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                fHistNGoodTracks,"","N_{Good Tracks}","#frac{dN}{dN_{Good Tracks}}",1,1,
                                0.82,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            WriteHistogram(fHistNGoodTracks);
            vecNGoodTracks.push_back(new TH1D(*fHistNGoodTracks));
        } else cout << "INFO: Object |fHistNGoodTracks| could not be found! Skipping Draw..." << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        // VZERO multiplicity
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        TH1D* fHistV0Mult = (TH1D*)ESDContainer->FindObject("V0 Multiplicity");
        if(fHistV0Mult){
            GetMinMaxBin(fHistV0Mult,minB,maxB);
            SetXRange(fHistV0Mult,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                fHistV0Mult,"","V0 Multiplicity","#frac{dN}{dV0Mult}",1,1,
                                0.82,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            WriteHistogram(fHistV0Mult);
            vecV0Mult.push_back(new TH1D(*fHistV0Mult));
        } else cout << "INFO: Object |V0 Multiplicity| could not be found! Skipping Draw..." << endl;
        TGaxis::SetExponentOffset(0, 0, "x");
        //-------------------------------------------------------------------------------------------------------------------------------
        // number of gamma candidates 
        TH1D* fHistGammaCandidates = (TH1D*)ESDContainer->FindObject("GammaCandidates");
        if(fHistGammaCandidates){
            GetMinMaxBin(fHistGammaCandidates,minB,maxB);
            SetXRange(fHistGammaCandidates,1,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                fHistGammaCandidates,"","N_{GammaCandidates}","#frac{dN}{dN_{Cand}}",1,1,
                                0.82,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            WriteHistogram(fHistGammaCandidates);
            vecGammaCandidates.push_back(new TH1D(*fHistGammaCandidates));
        } else cout << "INFO: Object |GammaCandidates| could not be found! Skipping Draw..." << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        // SPD tracklets vs SPD cluster before cuts
        TH2D* fHistSPDtracklets_clusters_before = (TH2D*)TopContainerEvent->FindObject(Form("SPD tracklets vs SPD clusters %s",ContainerEventCut.Data()));
        if(fHistSPDtracklets_clusters_before){
            if(isMinMaxSPD){
                GetMinMaxBin(fHistSPDtracklets_clusters_before,minB_SPD,maxB_SPD); maxB_SPD+=5;
                GetMinMaxBinY(fHistSPDtracklets_clusters_before,minYB_SPD,maxYB_SPD); maxYB_SPD-=1;
                isMinMaxSPD=kFALSE;
            }
            SetXRange(fHistSPDtracklets_clusters_before,1,maxB_SPD);
            SetYRange(fHistSPDtracklets_clusters_before,1,maxYB_SPD);
            SetZMinMaxTH2(fHistSPDtracklets_clusters_before,1,maxB_SPD,1,maxYB_SPD);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                                fHistSPDtracklets_clusters_before,"",
                                "N_{SPD tracklets}","N_{SPD clusters}",1,1.4,
                                0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            TF1 *cut = new TF1("cut","65. + 4 * x",fHistSPDtracklets_clusters_before->GetXaxis()->GetBinLowEdge(1),fHistSPDtracklets_clusters_before->GetXaxis()->GetBinUpEdge(maxB_SPD+5));
            cut->SetLineColor(kRed);
            cut->SetLineStyle(2);
            cut->SetLineWidth(4);
            cut->Draw("SAME");
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistSPDtracklets_clusters_before, Form("%s/SPD_TrackletsVsClusters_Before_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
            delete cut;
        } else cout << Form("INFO: Object |SPD tracklets vs SPD clusters %s| could not be found! Skipping Draw...",ContainerEventCut.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        // SPD tracklets vs SPD cluster before pileup cuts
        TH2D* fHistSPDtracklets_clusters_beforePileUp = (TH2D*)ConvEventCutsContainer->FindObject(Form("SPD tracklets vs SPD clusters %s before Pileup Cut",fEventCutSelection[i].Data()));
        if(fHistSPDtracklets_clusters_beforePileUp){
            SetXRange(fHistSPDtracklets_clusters_beforePileUp,1,maxB_SPD);
            SetYRange(fHistSPDtracklets_clusters_beforePileUp,1,maxYB_SPD);
            SetZMinMaxTH2(fHistSPDtracklets_clusters_beforePileUp,1,maxB_SPD,1,maxYB_SPD);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                                fHistSPDtracklets_clusters_beforePileUp,"",
                                "N_{SPD tracklets}","N_{SPD clusters}",1,1.4,
                                0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            TF1 *cut = new TF1("cut","65. + 4 * x",fHistSPDtracklets_clusters_beforePileUp->GetXaxis()->GetBinLowEdge(1),fHistSPDtracklets_clusters_beforePileUp->GetXaxis()->GetBinUpEdge(maxB_SPD+5));
            cut->SetLineColor(kRed);
            cut->SetLineStyle(2);
            cut->SetLineWidth(4);
            cut->Draw("SAME");
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistSPDtracklets_clusters_beforePileUp, Form("%s/SPD_TrackletsVsClusters_BeforePileupCut_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
            delete cut;
        } else cout << Form("INFO: Object |SPD tracklets vs SPD clusters %s before Pileup Cut| could not be found! Skipping Draw...",ContainerEventCut.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        // SPD tracklets vs SPD cluster after cuts
        TH2D* fHistSPDtracklets_clusters = (TH2D*)ESDContainer->FindObject("SPD tracklets vs SPD clusters");
        if(fHistSPDtracklets_clusters){
            SetXRange(fHistSPDtracklets_clusters,1,maxB_SPD);
            SetYRange(fHistSPDtracklets_clusters,1,maxYB_SPD);
            SetZMinMaxTH2(fHistSPDtracklets_clusters,1,maxB_SPD,1,maxYB_SPD);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                                fHistSPDtracklets_clusters,"",
                                "N_{SPD tracklets}","N_{SPD clusters}",1,1.4,
                                0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            TF1 *cut = new TF1("cut","65. + 4 * x",fHistSPDtracklets_clusters->GetXaxis()->GetBinLowEdge(1),fHistSPDtracklets_clusters->GetXaxis()->GetBinUpEdge(maxB_SPD+5));
            cut->SetLineColor(kRed);
            cut->SetLineStyle(2);
            cut->SetLineWidth(4);
            cut->Draw("SAME");
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistSPDtracklets_clusters, Form("%s/SPD_TrackletsVsClusters_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
            delete cut;
        } else cout << "INFO: Object |SPD tracklets vs SPD clusters| could not be found! Skipping Draw..." << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        // number of gamma candidates vs track candidates in central acceptance
        TH2D* fHistTracksVsCandidates = (TH2D*)ESDContainer->FindObject("GoodESDTracksVsGammaCandidates");
        if(fHistTracksVsCandidates){
            GetMinMaxBin(fHistTracksVsCandidates,minB,maxB);
            SetXRange(fHistTracksVsCandidates,1,maxB+1);
            GetMinMaxBinY(fHistTracksVsCandidates,minYB,maxYB);
            SetYRange(fHistTracksVsCandidates,1,maxYB+1);
            SetZMinMaxTH2(fHistTracksVsCandidates,1,maxB+1,1,maxYB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                                fHistTracksVsCandidates,"",
                                "N_{Good Tracks}","N_{GammaCandidates}",1,1.4,
                                0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistTracksVsCandidates, Form("%s/GoodESDTracksVsGammaCandidates_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
        } else cout << "INFO: Object |GoodESDTracksVsGammaCandidates| could not be found! Skipping Draw..." << endl;
        
        //-------------------------------------------------------------------------------------------------------------------------------

        Float_t nEventsBin1 = fHistNEvents->GetBinContent(1);
        Float_t nEventsBin2 = fHistNEvents->GetBinContent(2);
        Float_t nEventsBin3 = fHistNEvents->GetBinContent(3);
        Float_t nEventsBin4 = fHistNEvents->GetBinContent(4);
        Float_t nEventsBin5 = fHistNEvents->GetBinContent(5);
        Float_t nEventsBin6 = fHistNEvents->GetBinContent(6);
        Float_t nEventsBin7 = fHistNEvents->GetBinContent(7);
        Float_t nEventsBin12 = 0;
        if(fHistNEvents->GetNbinsX()==12) nEventsBin12 = fHistNEvents->GetBinContent(12);

        //-------------------------------------------------------------------------------------------------------------------------------

        Float_t nEventsAllEvt = nEventsBin1+nEventsBin2+nEventsBin3+nEventsBin4+nEventsBin5+nEventsBin6+nEventsBin7;
        Float_t nEventsAllEvtErr = sqrt(nEventsAllEvt);
        Float_t nEventsMinBiasEvt = nEventsBin1+nEventsBin2+nEventsBin5+nEventsBin6+nEventsBin7;
        Float_t nEventsMinBiasEvtErr = sqrt(nEventsMinBiasEvt);
        Float_t nEventsNormEvt = nEventsBin1+(nEventsBin1/(nEventsBin1+nEventsBin5))*nEventsBin6;
        Float_t nEventsNormEvtErr = sqrt(nEventsNormEvt);

        cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        cout << "nEventsAllEvt: \t\t\t(" << nEventsAllEvt << " +- " << nEventsAllEvtErr << ")" << endl;
        cout << "nEventsMinBiasEvt: \t\t(" << nEventsMinBiasEvt << " +- " << nEventsMinBiasEvtErr << ")" << endl;
        cout << "nEventsNormEvt: \t\t(" << nEventsNormEvt << " +- " << nEventsNormEvtErr << ")" << endl;
        fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        fLog << "nEventsAllEvt: \t\t\t(" << nEventsAllEvt << " +- " << nEventsAllEvtErr << ")" << endl;
        fLog << "nEventsMinBiasEvt: \t\t(" << nEventsMinBiasEvt << " +- " << nEventsMinBiasEvtErr << ")" << endl;
        fLog << "nEventsNormEvt: \t\t(" << nEventsNormEvt << " +- " << nEventsNormEvtErr << ")" << endl;

        //-------------------------------------------------------------------------------------------------------------------------------
        //-------------------------------------------------------------------------------------------------------------------------------

        Float_t ratioWOVtxEvt = nEventsBin6/nEventsMinBiasEvt;
        Float_t ratioWOVtxEvtErr = sqrt( pow(sqrt(nEventsBin6)/nEventsMinBiasEvt,2)  + pow( nEventsMinBiasEvtErr*nEventsBin6/pow(nEventsMinBiasEvt,2),2) );
        Float_t ratioWVtxOutside10cmEvt = nEventsBin5/nEventsMinBiasEvt;
        Float_t ratioWVtxOutside10cmEvtErr = sqrt( pow(sqrt(nEventsBin5)/nEventsMinBiasEvt,2)  + pow( nEventsMinBiasEvtErr*nEventsBin5/pow(nEventsMinBiasEvt,2),2) );
        Float_t ratioPileUpEvt = nEventsBin7/nEventsMinBiasEvt;
        Float_t ratioPileUpEvtErr = sqrt( pow(sqrt(nEventsBin7)/nEventsMinBiasEvt,2)  + pow( nEventsMinBiasEvtErr*nEventsBin7/pow(nEventsMinBiasEvt,2),2) );
        Float_t ratioSPDClusTrackEvt = nEventsBin12/nEventsMinBiasEvt;
        Float_t ratioSPDClusTrackEvtErr = sqrt( pow(sqrt(nEventsBin12)/nEventsMinBiasEvt,2)  + pow( nEventsMinBiasEvtErr*nEventsBin12/pow(nEventsMinBiasEvt,2),2) );

        Float_t ratioGoodEventsEvt = nEventsBin1/nEventsMinBiasEvt;
        Float_t ratioGoodEventsEvtErr = sqrt( pow(sqrt(nEventsBin1)/nEventsMinBiasEvt,2)  + pow( nEventsMinBiasEvtErr*nEventsBin1/pow(nEventsMinBiasEvt,2),2) );
        Float_t ratioNormAllEvt = nEventsNormEvt/nEventsAllEvt;
        Float_t ratioNormAllEvtErr = sqrt( pow(sqrt(nEventsNormEvt)/nEventsAllEvt,2)  + pow( nEventsAllEvtErr*nEventsNormEvt/pow(nEventsAllEvt,2),2) );
        Float_t ratioNormMinBiasEvt = nEventsNormEvt/nEventsMinBiasEvt;
        Float_t ratioNormMinBiasEvtErr = sqrt( pow(sqrt(nEventsNormEvt)/nEventsMinBiasEvt,2)  + pow( nEventsMinBiasEvtErr*nEventsNormEvt/pow(nEventsMinBiasEvt,2),2) );

        cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        cout << "ratioWOVtxEvt: \t\t\t(" << ratioWOVtxEvt << " +- " << ratioWOVtxEvtErr << ")" << endl;
        cout << "ratioWVtxOutside10cmEvt: \t(" << ratioWVtxOutside10cmEvt << " +- " << ratioWVtxOutside10cmEvtErr << ")" << endl;
        cout << "ratioPileUpEvt: \t\t(" << ratioPileUpEvt << " +- " << ratioPileUpEvtErr << ")" << endl;
        cout << "ratioSPDClusTrackEvt: \t\t(" << ratioSPDClusTrackEvt << " +- " << ratioSPDClusTrackEvtErr << ")" << endl;
        cout << "ratioGoodEventsEvt: \t\t(" << ratioGoodEventsEvt << " +- " << ratioGoodEventsEvtErr << ")" << endl;
        cout << "ratioNormAllEvt: \t\t(" << ratioNormAllEvt << " +- " << ratioNormAllEvtErr << ")" << endl;
        cout << "ratioNormMinBiasEvt: \t\t(" << ratioNormMinBiasEvt << " +- " << ratioNormMinBiasEvtErr << ")" << endl;
        fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        fLog << "ratioWOVtxEvt: \t\t\t(" << ratioWOVtxEvt << " +- " << ratioWOVtxEvtErr << ")" << endl;
        fLog << "ratioWVtxOutside10cmEvt: \t(" << ratioWVtxOutside10cmEvt << " +- " << ratioWVtxOutside10cmEvtErr << ")" << endl;
        fLog << "ratioPileUpEvt: \t\t(" << ratioPileUpEvt << " +- " << ratioPileUpEvtErr << ")" << endl;
        fLog << "ratioSPDClusTrackEvt: \t\t(" << ratioSPDClusTrackEvt << " +- " << ratioSPDClusTrackEvtErr << ")" << endl;
        fLog << "ratioGoodEventsEvt: \t\t(" << ratioGoodEventsEvt << " +- " << ratioGoodEventsEvtErr << ")" << endl;
        fLog << "ratioNormAllEvt: \t\t(" << ratioNormAllEvt << " +- " << ratioNormAllEvtErr << ")" << endl;
        fLog << "ratioNormMinBiasEvt: \t\t(" << ratioNormMinBiasEvt << " +- " << ratioNormMinBiasEvtErr << ")" << endl;
        
        
        //-------------------------------------------------------------------------------------------------------------------------------
        //---------------------------------------------- Cluster properties -------------------------------------------------------------
        //-------------------------------------------------------------------------------------------------------------------------------
        TH1D* fHistAcceptance = 0x0;
        TH1D* fHistCutIndex = 0x0;
        if(isCalo){
            fHistAcceptance = (TH1D*)CaloCutsContainer->FindObject(Form("AcceptanceCuts %s", fClusterCutSelection[i].Data()));
            if(fHistAcceptance){
                DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kFALSE,kFALSE,
                                    fHistAcceptance,"","","# of Entries",1,1,
                                    0.82,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
                SaveCanvasAndWriteHistogram(canvas, fHistAcceptance, Form("%s/AcceptanceCuts_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
            } else cout << "INFO: Object |AcceptanceCuts| could not be found! Skipping Draw..." << endl;
            //-------------------------------------------------------------------------------------------------------------------------------
            fHistCutIndex = (TH1D*)CaloCutsContainer->FindObject(Form("IsPhotonSelected %s", fClusterCutSelection[i].Data()));
            if(fHistCutIndex){
                DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kFALSE,kFALSE,
                                    fHistCutIndex,"","","# of Entries",1,1,
                                    0.82,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
                SaveCanvasAndWriteHistogram(canvas, fHistCutIndex, Form("%s/IsPhotonSelected_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
            } else cout << "INFO: Object |fHistCutIndex| could not be found! Skipping Draw..." << endl;
            //-------------------------------------------------------------------------------------------------------------------------------
            TH1D* fHistClusterIdentificationCuts = (TH1D*)CaloCutsContainer->FindObject(Form("ClusterQualityCuts %s", fClusterCutSelection[i].Data()));
            TH2D* fHistClusterIdentificationCuts2D = (TH2D*)CaloCutsContainer->FindObject(Form("ClusterQualityCuts vs E %s", fClusterCutSelection[i].Data()));
            if(fHistClusterIdentificationCuts){
                DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kFALSE,kFALSE,
                                    fHistClusterIdentificationCuts,"","","# of Entries",1,1,
                                    0.82,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
                SaveCanvasAndWriteHistogram(canvas, fHistClusterIdentificationCuts, Form("%s/ClusterQualityCuts_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
            } else if(fHistClusterIdentificationCuts2D && fHistClusterIdentificationCuts2D->IsA()==TH2F::Class()){
                GetMinMaxBinY(fHistClusterIdentificationCuts2D,minB,maxB);
                SetYRange(fHistClusterIdentificationCuts2D,1,maxB+1);
                SetZMinMaxTH2(fHistClusterIdentificationCuts2D,1,fHistClusterIdentificationCuts2D->GetNbinsX(),1,maxB+1);
                DrawPeriodQAHistoTH2(canvas,leftMargin,0.1,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                                    fHistClusterIdentificationCuts2D,Form("%s - %s - %s",fCollisionSystem.Data(), plotDataSets[i].Data(), fClusters.Data()),
                                    "","Cluster Energy (GeV)",0.9,0.8);
                SaveCanvasAndWriteHistogram(canvas, fHistClusterIdentificationCuts2D, Form("%s/ClusterQualityCuts_vs_E_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));

                PlotCutHistoReasons(canvas,leftMargin,rightMargin,topMargin,bottomMargin, fHistClusterIdentificationCuts2D, "#it{E}_{Cluster}", "#frac{d#it{E}_{Cluster}}{dN}",
                                    5,10,0,0);
                SaveCanvas(canvas, Form("%s/ClusterQualityCuts_Projected_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()), kTRUE, kTRUE);
            } else cout << "INFO: Object |fHistClusterIdentificationCuts (TH2 vs pT)| could not be found! Skipping Draw..." << endl;

            //-------------------------------------------------------------------------------------------------------------------------------
            //---------------------------------------------- statistics for Calo photons -----------------------------------------------------
            //-------------------------------------------------------------------------------------------------------------------------------        
            if(fHistAcceptance && fHistCutIndex){
                Double_t CaloNClusters = fHistAcceptance->GetBinContent(1);
                Double_t CaloNClustersQA = fHistCutIndex->GetBinContent(5);
                CaloNClusters/=nEvents[i];
                CaloNClustersQA/=nEvents[i];
                cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
                cout << "CaloNClusters: \t\t\t(" << CaloNClusters << " +- " << sqrt(CaloNClusters) << ")" << endl;
                cout << "CaloNClustersQA: \t\t(" << CaloNClustersQA << " +- " << sqrt(CaloNClustersQA) << ")" << endl;
                fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
                fLog << "CaloNClusters: \t\t\t(" << CaloNClusters << " +- " << sqrt(CaloNClusters) << ")" << endl;
                fLog << "CaloNClustersQA: \t\t(" << CaloNClustersQA << " +- " << sqrt(CaloNClustersQA) << ")" << endl;
            }
            
        }

        
        //-------------------------------------------------------------------------------------------------------------------------------
        //------------------------------------------- PCM photon properties -------------------------------------------------------------
        //-------------------------------------------------------------------------------------------------------------------------------
        TH1D* fHistConvGamma = 0x0;
        if(isConv){
            fHistConvGamma = (TH1D*)ESDContainer->FindObject("ESD_ConvGamma_Pt");
            if(!fHistConvGamma) {
                cout << Form("ERROR: Object |ESD_ConvGamma_Pt| could not be found, although running mode: '%i'! Returning...",fMode) << endl;
                return;
            }
            //-------------------------------------------------------------------------------------------------------------------------------
            TH2D* fHistPhotonCuts2DPreSel = (TH2D*)TopContainerGamma->FindObject(Form("PhotonCuts %s", ContainerGammaCut.Data()));
            if(fHistPhotonCuts2DPreSel && fHistPhotonCuts2DPreSel->IsA()==TH2F::Class()){
                GetMinMaxBinY(fHistPhotonCuts2DPreSel,minB,maxB);
                SetYRange(fHistPhotonCuts2DPreSel,1,maxB+1);
                SetZMinMaxTH2(fHistPhotonCuts2DPreSel,1,fHistPhotonCuts2DPreSel->GetNbinsX(),1,maxB+1);
                DrawPeriodQAHistoTH2(canvas,leftMargin,0.1,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                                    fHistPhotonCuts2DPreSel,Form("%s - %s - %s",fCollisionSystem.Data(), plotDataSets[i].Data(), fClusters.Data()),
                                    "","#it{p}_{T, #gamma}",0.9,0.8);
                SaveCanvasAndWriteHistogram(canvas, fHistPhotonCuts2DPreSel, Form("%s/Preselect_PhotonCuts_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));

                PlotCutHistoReasons(canvas,leftMargin,rightMargin,topMargin,bottomMargin, fHistPhotonCuts2DPreSel, "#it{p}_{T, #gamma}", "#frac{d#it{p}_{T, #gamma}}{dN}",
                                    5,10,0,0);
                SaveCanvas(canvas, Form("%s/Preselect_PhotonCuts_Projected_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()), kTRUE, kTRUE);
            } else cout << Form("INFO: Object |Preselect_PhotonCuts %s (TH2 vs pT)| could not be found! Skipping Draw...", ContainerGammaCut.Data()) << endl;
            //-------------------------------------------------------------------------------------------------------------------------------
            TH2D* fHistPhotonAccCuts2DPreSel = (TH2D*)TopContainerGamma->FindObject(Form("PhotonAcceptanceCuts %s", ContainerGammaCut.Data()));
            if(fHistPhotonAccCuts2DPreSel && fHistPhotonAccCuts2DPreSel->IsA()==TH2F::Class()){
                GetMinMaxBinY(fHistPhotonAccCuts2DPreSel,minB,maxB);
                SetYRange(fHistPhotonAccCuts2DPreSel,1,maxB+1);
                SetZMinMaxTH2(fHistPhotonAccCuts2DPreSel,1,fHistPhotonAccCuts2DPreSel->GetNbinsX(),1,maxB+1);
                DrawPeriodQAHistoTH2(canvas,leftMargin,0.1,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                                    fHistPhotonAccCuts2DPreSel,Form("%s - %s - %s",fCollisionSystem.Data(), plotDataSets[i].Data(), fClusters.Data()),
                                    "","#it{p}_{T, #gamma}",0.9,0.8);
                SaveCanvasAndWriteHistogram(canvas, fHistPhotonAccCuts2DPreSel, Form("%s/Preselect_PhotonAcceptanceCuts_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));

                PlotCutHistoReasons(canvas,leftMargin,rightMargin,topMargin,bottomMargin, fHistPhotonAccCuts2DPreSel, "#it{p}_{T, #gamma}", "#frac{d#it{p}_{T, #gamma}}{dN}",
                                    5,10,0,0);
                SaveCanvas(canvas, Form("%s/Preselect_PhotonAcceptanceCuts_Projected_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()), kTRUE, kTRUE);
            } else cout << Form("INFO: Object |Preselect_PhotonAcceptanceCuts %s (TH2 vs pT)| could not be found! Skipping Draw...", ContainerGammaCut.Data()) << endl;
            //-------------------------------------------------------------------------------------------------------------------------------
            TH2D* fHistPhotonCuts2D = (TH2D*)ConvCutsContainer->FindObject(Form("PhotonCuts %s", fGammaCutSelection[i].Data()));
            if(fHistPhotonCuts2D && fHistPhotonCuts2D->IsA()==TH2F::Class()){
                GetMinMaxBinY(fHistPhotonCuts2D,minB,maxB);
                SetYRange(fHistPhotonCuts2D,1,maxB+1);
                SetZMinMaxTH2(fHistPhotonCuts2D,1,fHistPhotonCuts2D->GetNbinsX(),1,maxB+1);
                DrawPeriodQAHistoTH2(canvas,leftMargin,0.1,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                                    fHistPhotonCuts2D,Form("%s - %s - %s",fCollisionSystem.Data(), plotDataSets[i].Data(), fClusters.Data()),
                                    "","#it{p}_{T, #gamma}",0.9,0.8);
                SaveCanvasAndWriteHistogram(canvas, fHistPhotonCuts2D, Form("%s/PhotonCuts_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));

                PlotCutHistoReasons(canvas,leftMargin,rightMargin,topMargin,bottomMargin, fHistPhotonCuts2D, "#it{p}_{T, #gamma}", "#frac{d#it{p}_{T, #gamma}}{dN}",
                                    5,10,0,0);
                SaveCanvas(canvas, Form("%s/PhotonCuts_Projected_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()), kTRUE, kTRUE);
            } else cout << Form("INFO: Object |PhotonCuts %s (TH2 vs pT)| could not be found! Skipping Draw...", fGammaCutSelection[i].Data()) << endl;
            //-------------------------------------------------------------------------------------------------------------------------------
            TH2D* fHistPhotonAccCuts2D = (TH2D*)ConvCutsContainer->FindObject(Form("PhotonAcceptanceCuts %s", fGammaCutSelection[i].Data()));
            if(fHistPhotonAccCuts2D && fHistPhotonAccCuts2D->IsA()==TH2F::Class()){
                GetMinMaxBinY(fHistPhotonAccCuts2D,minB,maxB);
                SetYRange(fHistPhotonAccCuts2D,1,maxB+1);
                SetZMinMaxTH2(fHistPhotonAccCuts2D,1,fHistPhotonAccCuts2D->GetNbinsX(),1,maxB+1);
                DrawPeriodQAHistoTH2(canvas,leftMargin,0.1,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                                    fHistPhotonAccCuts2D,Form("%s - %s - %s",fCollisionSystem.Data(), plotDataSets[i].Data(), fClusters.Data()),
                                    "","#it{p}_{T, #gamma}",0.9,0.8);
                SaveCanvasAndWriteHistogram(canvas, fHistPhotonAccCuts2D, Form("%s/PhotonAcceptanceCuts_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));

                PlotCutHistoReasons(canvas,leftMargin,rightMargin,topMargin,bottomMargin, fHistPhotonAccCuts2D, "#it{p}_{T, #gamma}", "#frac{d#it{p}_{T, #gamma}}{dN}",
                                    5,10,0,0);
                SaveCanvas(canvas, Form("%s/PhotonAcceptanceCuts_Projected_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()), kTRUE, kTRUE);
            } else cout << Form("INFO: Object |PhotonAcceptanceCuts %s (TH2 vs pT)| could not be found! Skipping Draw...", fGammaCutSelection[i].Data()) << endl;
            //-------------------------------------------------------------------------------------------------------------------------------
            TH2D* fHistGammaAlphaQt = (TH2D*)TopContainerGamma->FindObject(Form("Armenteros_before %s",ContainerGammaCut.Data()));
            if(fHistGammaAlphaQt){
                fHistGammaAlphaQt->Sumw2();
                fHistGammaAlphaQt->Scale(1./fHistGammaAlphaQt->GetEntries());
                fHistGammaAlphaQt->GetYaxis()->SetRangeUser(0,0.25);
                SetZMinMaxTH2(fHistGammaAlphaQt,1,fHistGammaAlphaQt->GetNbinsX(),1,fHistGammaAlphaQt->GetXaxis()->FindBin(0.25),kTRUE);
                DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                                    fHistGammaAlphaQt,"",
                                    "#alpha = (#it{p}^{+}_{L}-#it{p}^{-}_{L})/(#it{p}^{+}_{L}+#it{p}^{-}_{L})","#it{q}_{T} (GeV/#it{c})",1,1.4,
                                    0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
                DrawArmenterosLines(kFALSE);
                SaveCanvasAndWriteHistogram(cvsQuadratic, fHistGammaAlphaQt, Form("%s/Armenteros_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
            } else cout << Form("INFO: Object |Armenteros_before %s| could not be found! Skipping Draw...",ContainerGammaCut.Data()) << endl;
            //-------------------------------------------------------------------------------------------------------------------------------
            TString str[2]={"before","after"};
            Int_t sigmadEdxMinMax[2]={0,0};
            Int_t sigmaTOFMinMax[2]={0,0};
            for(Int_t iBefore=0; iBefore<2; iBefore++){
                TList* list = TopContainerGamma;
                TString iCutNumber = ContainerGammaCut;
                if(iBefore==1){
                    list = ConvCutsContainer;
                    iCutNumber = fGammaCutSelection[i];
                }

                TH2D* histodEdx = (TH2D*)list->FindObject(Form("Gamma_dEdx_%s %s",str[iBefore].Data(),iCutNumber.Data()));
                if(histodEdx){
                    histodEdx->Sumw2();
                    histodEdx->Scale(1./histodEdx->GetEntries());
                    GetMinMaxBin(histodEdx,minB,maxB);
                    SetXRange(histodEdx,minB,maxB);
                    SetZMinMaxTH2(histodEdx,1,histodEdx->GetNbinsX(),1,histodEdx->GetNbinsY(),kTRUE);
                    DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                        histodEdx,"",
                                        "#it{p} (GeV/#it{c})","d#it{E}/d#it{x}",1,1.4,
                                        0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
                    SaveCanvasAndWriteHistogram(cvsQuadratic, histodEdx, Form("%s/dEdx_Candidates_%s_%s.%s", outputDir.Data(), str[iBefore].Data(),DataSets[i].Data(), suffix.Data()));
                } else cout << Form("INFO: Object |Gamma_dEdx_%s %s| could not be found! Skipping Draw...",str[iBefore].Data(),iCutNumber.Data()) << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                TH2D* histonSigmadEdx = (TH2D*)list->FindObject(Form("Gamma_dEdxSig_%s %s",str[iBefore].Data(),iCutNumber.Data()));
                if(histonSigmadEdx){
                    histonSigmadEdx->Sumw2();
                    histonSigmadEdx->Scale(1./histonSigmadEdx->GetEntries());
                    GetMinMaxBin(histonSigmadEdx,minB,maxB);
                    SetXRange(histonSigmadEdx,minB,maxB);
                    if(iBefore==0) GetMinMaxBinY(histonSigmadEdx,sigmadEdxMinMax[0],sigmadEdxMinMax[1]);
                    SetYRange(histonSigmadEdx,sigmadEdxMinMax[0]-1,sigmadEdxMinMax[1]+1);
                    SetZMinMaxTH2(histonSigmadEdx,1,histonSigmadEdx->GetNbinsX(),1,histonSigmadEdx->GetNbinsY(),kTRUE);
                    DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                        histonSigmadEdx,"",
                                        "#it{p} (GeV/#it{c})","#it{n} #sigma_{e^{#pm}} d#it{E}/d#it{x}",1,1.4,
                                        0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
                        SaveCanvasAndWriteHistogram(cvsQuadratic, histonSigmadEdx, Form("%s/nSigma_dEdx_Candidates_%s_%s.%s", outputDir.Data(), str[iBefore].Data(),DataSets[i].Data(), suffix.Data()));
                } else cout << Form("INFO: Object |Gamma_dEdxSig_%s %s| could not be found! Skipping Draw...",str[iBefore].Data(),iCutNumber.Data()) << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                TH2D* histoTOF = (TH2D*)list->FindObject(Form("Gamma_TOF_%s %s",str[iBefore].Data(),iCutNumber.Data()));
                if(iBefore==0){
                    if(histoTOF){
                    histoTOF->Sumw2();
                    histoTOF->Scale(1./histoTOF->GetEntries());
                    GetMinMaxBin(histoTOF,minB,maxB);
                    SetXRange(histoTOF,minB,maxB);
                    SetZMinMaxTH2(histoTOF,1,histoTOF->GetNbinsX(),1,histoTOF->GetNbinsY(),kTRUE);
                    DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                        histoTOF,"",
                                        "#it{p} (GeV/#it{c})","#it{t}_{measured}-#it{t}_{expected}",1,1.4,
                                        0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
                    SaveCanvasAndWriteHistogram(cvsQuadratic, histoTOF, Form("%s/TOF_Candidates_%s_%s.%s", outputDir.Data(),str[iBefore].Data(), DataSets[i].Data(), suffix.Data()));
                    } else cout << Form("INFO: Object |Gamma_TOF_%s %s| could not be found! Skipping Draw...",str[iBefore].Data(),iCutNumber.Data()) << endl;
                }
                //-------------------------------------------------------------------------------------------------------------------------------
                TH2D* histonSigmaTOF = (TH2D*)list->FindObject(Form("Gamma_TOFSig_%s %s",str[iBefore].Data(),iCutNumber.Data()));
                if(histonSigmaTOF){
                    histonSigmaTOF->Sumw2();
                    histonSigmaTOF->Scale(1./histonSigmaTOF->GetEntries());
                    GetMinMaxBin(histonSigmaTOF,minB,maxB);
                    SetXRange(histonSigmaTOF,minB,maxB);
                    if(iBefore==0) GetMinMaxBinY(histonSigmaTOF,sigmaTOFMinMax[0],sigmaTOFMinMax[1]);
                    SetYRange(histonSigmaTOF,sigmaTOFMinMax[0]-1,sigmaTOFMinMax[1]+1);
                    SetZMinMaxTH2(histonSigmaTOF,1,histonSigmadEdx->GetNbinsX(),1,histonSigmadEdx->GetNbinsY(),kTRUE);
                    DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                        histonSigmaTOF,"",
                                        "#it{p} (GeV/#it{c})","#it{n} #sigma_{e^{#pm}} TOF",1,1.4,
                                        0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
                    SaveCanvasAndWriteHistogram(cvsQuadratic, histonSigmaTOF, Form("%s/nSigma_TOF_Candidates_%s_%s.%s", outputDir.Data(), str[iBefore].Data(),DataSets[i].Data(), suffix.Data()));
                } else cout << Form("INFO: Object |Gamma_TOFSig_%s %s| could not be found! Skipping Draw...",str[iBefore].Data(),iCutNumber.Data()) << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                TH1D* fHistInvMass = (TH1D*)list->FindObject(Form("InvMass_%s %s",str[iBefore].Data(),iCutNumber.Data()));
                if(fHistInvMass){
                    fHistInvMass->Sumw2();
                    fHistInvMass->Scale(1./fHistConvGamma->Integral());
                    DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                        fHistInvMass,"","M_{e^{+}e^{-}} (GeV/#it{c}^{2})","#frac{d#it{N}_{#gamma}}{#it{N} d#it{M}_{e^{+}e^{-}}} (GeV/#it{c}^{2})^{-1}",1,1,
                                        0.82,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
                    WriteHistogram(fHistInvMass);
                    vecInvMassBeforeAfter[iBefore].push_back(new TH1D(*fHistInvMass));
                } else cout << Form("INFO: Object |InvMass_%s %s| could not be found! Skipping Draw...",str[iBefore].Data(),iCutNumber.Data()) << endl;
            }
            //-------------------------------------------------------------------------------------------------------------------------------
            //---------------------------------------------- statistics for PCM photons -----------------------------------------------------
            //-------------------------------------------------------------------------------------------------------------------------------

            if( fHistConvGamma){
                Float_t nGamma = fHistConvGamma->Integral();
                Float_t nGammaErr = sqrt(nGamma);
                Float_t ConvNCandidatesQA = (Float_t) nGamma/nEvents[i];
                Float_t ConvNCandidatesQAErr = sqrt( pow(nGammaErr/nEvents[i],2)  + pow( sqrt(nEvents[i])*nGamma/pow(nEvents[i],2),2) );
                cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
                cout << "ConvNCandidatesQA: \t\t(" << ConvNCandidatesQA << " +- " << ConvNCandidatesQAErr << ")" << endl;
                fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
                fLog << "ConvNCandidatesQA: \t\t(" << ConvNCandidatesQA << " +- " << ConvNCandidatesQAErr << ")" << endl;
            }

            
        }
        
        //-------------------------------------------------------------------------------------------------------------------------------
        //----------------------------------------------- Meson properties --------------------------------------------------------------
        //-------------------------------------------------------------------------------------------------------------------------------
        // meson cut tracking histo
        TH2D* fHistMesonCuts = (TH2D*) MesonCutsContainer->FindObject(Form("MesonCuts %s", fMesonCutSelection[i].Data()));
        if(fHistMesonCuts && fHistMesonCuts->IsA()==TH2F::Class()){
            GetMinMaxBinY(fHistMesonCuts,minB,maxB);
            SetYRange(fHistMesonCuts,1,maxB+1);
            SetZMinMaxTH2(fHistMesonCuts,1,fHistMesonCuts->GetNbinsX(),1,maxB+1);
            DrawPeriodQAHistoTH2(canvas,leftMargin,0.1,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                                fHistMesonCuts,Form("%s - %s - %s",fCollisionSystem.Data(), plotDataSets[i].Data(), fClusters.Data()),
                                "","#it{p}_{T}",0.9,0.8);
            SaveCanvasAndWriteHistogram(canvas, fHistMesonCuts, Form("%s/Meson_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));

            PlotCutHistoReasons(canvas,leftMargin,rightMargin,topMargin,bottomMargin, fHistMesonCuts, "#it{p}_{T}", "#frac{d#it{p}_{T}}{dN}",
                                5,10,0,0);
            SaveCanvas(canvas, Form("%s/Meson_Projected_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()), kTRUE, kTRUE);
        } else cout << Form("INFO: Object |MesonCuts %s (TH2 vs pT)| could not be found! Skipping Draw...", fMesonCutSelection[i].Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        // meson BG cut tracking histo
        TH2D* fHistMesonBGCuts = (TH2D*)MesonCutsContainer->FindObject(Form("MesonBGCuts %s", fMesonCutSelection[i].Data()));
        if(fHistMesonBGCuts && fHistMesonBGCuts->IsA()==TH2F::Class()){
            GetMinMaxBinY(fHistMesonBGCuts,minB,maxB);
            SetYRange(fHistMesonBGCuts,1,maxB+1);
            SetZMinMaxTH2(fHistMesonBGCuts,1,fHistMesonBGCuts->GetNbinsX(),1,maxB+1);
            DrawPeriodQAHistoTH2(canvas,leftMargin,0.1,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                                fHistMesonBGCuts,Form("%s - %s - %s",fCollisionSystem.Data(), plotDataSets[i].Data(), fClusters.Data()),
                                "","#it{p}_{T}",0.9,0.8);
            SaveCanvasAndWriteHistogram(canvas, fHistMesonBGCuts, Form("%s/MesonBG_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));

            PlotCutHistoReasons(canvas,leftMargin,rightMargin,topMargin,bottomMargin, fHistMesonBGCuts, "#it{p}_{T}", "#frac{d#it{p}_{T}}{dN}",
                                5,10,0,0);
            SaveCanvas(canvas, Form("%s/MesonBGCuts_Projected_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()), kTRUE, kTRUE);
        } else cout << Form("INFO: Object |MesonBGCuts %s (TH2 vs pT)| could not be found! Skipping Draw...", fMesonCutSelection[i].Data()) << endl;
        
        //-------------------------------------------------------------------------------------------------------------------------------
        // define name for meson property histograms
        TString namePi0MesonY       = "ESD_MotherPi0_Pt_Y";
        TString namePi0MesonAlpha   = "ESD_MotherPi0_Pt_Alpha";
        TString namePi0MesonOpen    = "ESD_MotherPi0_Pt_OpenAngle";
        if ( isMergedCalo ){
            namePi0MesonY           = "ESD_Mother_Pt_Y";
            namePi0MesonAlpha       = "ESD_Mother_Pt_Alpha";
            namePi0MesonOpen        = "ESD_Mother_Pt_OpenAngle";
        } 
        //-------------------------------------------------------------------------------------------------------------------------------
        // neutral pion pt vs rapidity
        TH2D* Pi0PtY = (TH2D*) ESDContainer->FindObject(namePi0MesonY.Data());
        if(Pi0PtY){
            TH1D* Pi0Pt = (TH1D*) Pi0PtY->ProjectionX("Pi0_Pt");
            GetMinMaxBin(Pi0Pt,minB,maxB);
            SetXRange(Pi0Pt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                Pi0Pt,"","#it{p}_{T}^{#pi^{0}}","#frac{dN^{#pi^{0}}}{d#it{p}_{T}}",1,1,
                                0.82,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            WriteHistogram(Pi0Pt);
            vecPi0Pt.push_back(new TH1D(*Pi0Pt));

            TH1D* Pi0Y = (TH1D*) Pi0PtY->ProjectionY("Pi0_Y");
            GetMinMaxBin(Pi0Y,minB,maxB);
            SetXRange(Pi0Y,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                Pi0Y,"","Y^{#pi^{0}}","#frac{dN^{#pi^{0}}}{dY}",1,1,
                                0.82,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            WriteHistogram(Pi0Y);
            vecPi0Y.push_back(new TH1D(*Pi0Y));
        } else cout << "INFO: Object |"<< namePi0MesonY.Data()<<"| could not be found! Skipping Draw..." << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        // neutral pion pt vs alpha
        TH2D* Pi0PtAlpha = (TH2D*) ESDContainer->FindObject(namePi0MesonAlpha.Data());
        if(Pi0PtAlpha){
            TH1D* Pi0Alpha = (TH1D*) Pi0PtAlpha->ProjectionY("Pi0_Alpha");
            GetMinMaxBin(Pi0Alpha,minB,maxB);
            SetXRange(Pi0Alpha,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                Pi0Alpha,"","#alpha^{#pi^{0}}","#frac{dN^{#pi^{0}}}{d#alpha}",1,1,
                                0.82,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            WriteHistogram(Pi0Alpha);
            vecPi0Alpha.push_back(new TH1D(*Pi0Alpha));
        }else cout << "INFO: Object |"<< namePi0MesonAlpha.Data()<<"| could not be found! Skipping Draw..." << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        // neutral pion pt vs opening angle
        TH2D* Pi0PtOpenAngle = (TH2D*) ESDContainer->FindObject(namePi0MesonOpen.Data());
        if(Pi0PtOpenAngle){
            TH1D* Pi0OpenAngle = (TH1D*) Pi0PtOpenAngle->ProjectionY("Pi0_Alpha");
            GetMinMaxBin(Pi0OpenAngle,minB,maxB);
            SetXRange(Pi0OpenAngle,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                Pi0OpenAngle,"","#theta_{#pi^{0}}","#frac{dN^{#pi^{0}}}{d#theta}",1,1,
                                0.82,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            WriteHistogram(Pi0OpenAngle);
            vecPi0OpenAngle.push_back(new TH1D(*Pi0OpenAngle));
        }else cout << "INFO: Object |"<< namePi0MesonOpen.Data()<<"| could not be found! Skipping Draw..." << endl;
        
        //-------------------------------------------------------------------------------------------------------------------------------
        // Eta properties 
        if ( !isMergedCalo ){
            //-------------------------------------------------------------------------------------------------------------------------------
            // eta meson pt vs rapidity
            TH2D* EtaPtY = (TH2D*) ESDContainer->FindObject("ESD_MotherEta_Pt_Y");
            if(EtaPtY){
                TH1D* EtaPt = (TH1D*) EtaPtY->ProjectionX("Eta_Pt");
                GetMinMaxBin(EtaPt,minB,maxB);
                SetXRange(EtaPt,minB,maxB);
                DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                    EtaPt,"","#it{p}_{T}^{#eta}","#frac{dN^{#eta}}{d#it{p}_{T}}",1,1,
                                    0.82,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
                WriteHistogram(EtaPt);
                vecEtaPt.push_back(new TH1D(*EtaPt));

                TH1D* EtaY = (TH1D*) EtaPtY->ProjectionY("Eta_Y");
                GetMinMaxBin(EtaY,minB,maxB);
                SetXRange(EtaY,minB,maxB);
                DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                    EtaY,"","Y^{#eta}","#frac{dN^{#eta}}{dY}",1,1,
                                    0.82,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
                WriteHistogram(EtaY);
                vecEtaY.push_back(new TH1D(*EtaY));
            }else cout << "INFO: Object |ESD_MotherEta_Pt_Y| could not be found! Skipping Draw..." << endl;
            //-------------------------------------------------------------------------------------------------------------------------------
            // eta meson pt vs alpha
            TH2D* EtaPtAlpha = (TH2D*) ESDContainer->FindObject("ESD_MotherEta_Pt_Alpha");
            if(EtaPtAlpha){
                TH1D* EtaAlpha = (TH1D*) EtaPtAlpha->ProjectionY("Eta_Alpha");
                GetMinMaxBin(EtaAlpha,minB,maxB);
                SetXRange(EtaAlpha,minB,maxB);
                DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                    EtaAlpha,"","#alpha^{#eta}","#frac{dN^{#eta}}{d#alpha}",1,1,
                                    0.82,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
                WriteHistogram(EtaAlpha);
                vecEtaAlpha.push_back(new TH1D(*EtaAlpha));
            }else cout << "INFO: Object |ESD_MotherEta_Pt_Alpha| could not be found! Skipping Draw..." << endl;
            //-------------------------------------------------------------------------------------------------------------------------------
            // eta meson pt vs opening angle
            TH2D* EtaPtOpenAngle = (TH2D*) ESDContainer->FindObject("ESD_MotherEta_Pt_OpenAngle");
            if(EtaPtOpenAngle){
                TH1D* EtaOpenAngle = (TH1D*) EtaPtOpenAngle->ProjectionY("Eta_Alpha");
                GetMinMaxBin(EtaOpenAngle,minB,maxB);
                SetXRange(EtaOpenAngle,minB,maxB);
                DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                    EtaOpenAngle,"","#theta_{#eta}","#frac{dN^{#eta}}{d#theta}",1,1,
                                    0.82,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
                WriteHistogram(EtaOpenAngle);
                vecEtaOpenAngle.push_back(new TH1D(*EtaOpenAngle));
            }else cout << "INFO: Object |ESD_MotherEta_Pt_OpenAngle| could not be found! Skipping Draw..." << endl;
        }    
        //-------------------------------------------------------------------------------------------------------------------------------
        // Invariant mass of meson candidates
        TH2D* ESDMother = (TH2D*) ESDContainer->FindObject("ESD_Mother_InvMass_Pt");
        if(ESDMother){
            WriteHistogram(ESDMother);
            vecESDMother.push_back(new TH2D(*ESDMother));
        }else {cout << "ERROR: Object |ESD_Mother_InvMass_Pt| could not be found! Skipping Draw & return..." << endl; return;}
        //-------------------------------------------------------------------------------------------------------------------------------
        // Invariant mass of BG candidates
        TH2D* ESDBackground = NULL;
        if (mode != 10 && mode != 11){
            TH2D* ESDBackground = (TH2D*) ESDContainer->FindObject("ESD_Background_InvMass_Pt");
            if(ESDBackground){
                WriteHistogram(ESDBackground);
                vecESDBackground.push_back(new TH2D(*ESDBackground));
            }else {cout << "ERROR: Object |ESD_Background_InvMass_Pt| could not be found! Skipping Draw & return..." << endl; return;}
        }   
        //-------------------------------------------------------------------------------------------------------------------------------
        //--------------------------------- statistics for meson properties -------------------------------------------------------------
        //-------------------------------------------------------------------------------------------------------------------------------
        if ( !isMergedCalo ){
            Bool_t kScs = fitter.DoFitting((TH2D*)vecESDMother.at(i), (TH2D*)vecESDBackground.at(i), nEventsBin1, fMode, outputDir.Data(), DataSets[i],kFALSE,kTRUE,fLog);

            if(kScs){
            TH1D* tempPi = new TH1D(*fitter.GetSignalPi0());
            tempPi->GetXaxis()->SetRangeUser(0.,0.25);
            signalPi0.push_back(tempPi);

            TH1D* tempEta = new TH1D(*fitter.GetSignalEta());
            tempEta->GetXaxis()->SetRangeUser(0.4,0.7);
            signalEta.push_back(tempEta);
            }
            Double_t ratioPi0 = 0; Double_t ratioPi0Err = 0;
            Double_t ratioEta = 0; Double_t ratioEtaErr = 0;

            fitter.GetMesonRatios(kTRUE,ratioPi0,ratioPi0Err);
            fitter.GetMesonRatios(kFALSE,ratioEta,ratioEtaErr);

            cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
            cout << "ratioPi0: \t\t\t(" << ratioPi0 << " +- " << ratioPi0Err << ")" << endl;
            cout << "ratioEta: \t\t\t(" << ratioEta << " +- " << ratioEtaErr << ")" << endl;
            cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
            cout << endl;
            fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
            fLog << "ratioPi0: \t\t\t(" << ratioPi0 << " +- " << ratioPi0Err << ")" << endl;
            fLog << "ratioEta: \t\t\t(" << ratioEta << " +- " << ratioEtaErr << ")" << endl;
            fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
            fLog << endl;
        }

        //-------------------------------------------------------------------------------------------------------------------------------
        //----------------------------------- clean up ----------------------------------------------------------------------------------
        //-------------------------------------------------------------------------------------------------------------------------------
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

    std::vector<TH1D*> temp;
    
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "Drawing Comparison Histograms" << endl;
    //-------------------------------------------------------------------------------------------------------------------------------
    //---------------------------------------- Event histogram comparisons ----------------------------------------------------------
    //-------------------------------------------------------------------------------------------------------------------------------
    // event statistics overview
    for(Int_t iVec=0; iVec<(Int_t)vecNEvents.size(); iVec++){
        TH1D* temp = vecNEvents.at(iVec);
        temp->Sumw2();
    }
    AdjustHistRange(vecNEvents,1,10,kTRUE,1,0.1);
    DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                vecNEvents,"","","N_{Events}",1,1.1,
                                labelData, colorCompare, kFALSE, 1.1, 1.1, kTRUE,
                                0.82,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
    SaveCanvas(canvas, Form("%s/Comparison/NEvents.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);

    DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                    vecNEvents,"","","N_{Events}",1,1.1,
                                    labelData, colorCompare, kTRUE, 1.1, 1.1, kTRUE,
                                    0.82,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
    SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_NEvents.%s", outputDir.Data(), suffix.Data()));
    //-------------------------------------------------------------------------------------------------------------------------------
    // z-vertex distribution
    GetMinMaxBin(vecVertexZ,minB,maxB);
    for(Int_t iVec=0; iVec<(Int_t)vecVertexZ.size(); iVec++){
        TH1D* temp = vecVertexZ.at(iVec);
        temp->Sumw2();
        temp->Scale(1./temp->Integral());
        //temp->Scale(1./nEvents[iVec]);
        SetXRange(temp,minB,maxB);
    }
    DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                        vecVertexZ,"","Vertex z (cm)","#frac{1}{N} #frac{dN}{dz}",1,1.1,
                        labelData, colorCompare, kTRUE, 1.1, 1.1, kTRUE,
                        0.82,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
    SaveCanvas(canvas, Form("%s/Comparison/Vertex_Z.%s", outputDir.Data(), suffix.Data()));

    DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                        vecVertexZ,"","Vertex z (cm)","#frac{1}{N} #frac{dN}{dz}",1,1.1,
                        labelData, colorCompare, kTRUE, 1.1, 1.1, kTRUE,
                        0.82,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
    SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Vertex_Z.%s", outputDir.Data(), suffix.Data()));
    //-------------------------------------------------------------------------------------------------------------------------------
    // number of good tracks
    GetMinMaxBin(vecNGoodTracks,minB,maxB);
    for(Int_t iVec=0; iVec<(Int_t)vecNGoodTracks.size(); iVec++){
        TH1D* temp = vecNGoodTracks.at(iVec);
        temp->Sumw2();
        temp->Scale(1./temp->Integral());
        //temp->Scale(1./nEvents[iVec]);
        SetXRange(temp,1,maxB);
    }
    DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                vecNGoodTracks,"","Number of Good Tracks","#frac{1}{N} #frac{dN}{dNTracks}",1,1.1,
                                labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                0.82,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
    SaveCanvas(canvas, Form("%s/Comparison/NGoodTracks.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);

    DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                    vecNGoodTracks,"","Number of Good Tracks","#frac{1}{N} #frac{dN}{dNTracks}",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                    0.82,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
    SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_NGoodTracks.%s", outputDir.Data(), suffix.Data()));
    //-------------------------------------------------------------------------------------------------------------------------------
    // VZERO multiplicity
    GetMinMaxBin(vecV0Mult,minB,maxB);
    for(Int_t iVec=0; iVec<(Int_t)vecV0Mult.size(); iVec++){
        TH1D* temp = vecV0Mult.at(iVec);
        temp->Sumw2();
        temp->Scale(1./temp->Integral());
        //temp->Scale(1./nEvents[iVec]);
        SetXRange(temp,1,maxB);
    }
    TGaxis::SetExponentOffset(-0.03, -0.04, "x");
    DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                vecV0Mult,"","V0 Multiplicity","#frac{1}{N} #frac{dN}{dV0Mult}",1,1.1,
                                labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                0.82,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
    SaveCanvas(canvas, Form("%s/Comparison/V0Mult.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
    TGaxis::SetExponentOffset(0, 0, "x");

    DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                vecV0Mult,"","V0 Multiplicity","#frac{1}{N} #frac{dN}{dV0Mult}",1,1.1,
                                labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                0.82,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
    SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_V0Mult.%s", outputDir.Data(), suffix.Data()));
    //-------------------------------------------------------------------------------------------------------------------------------
    // Gamma candidates
    GetMinMaxBin(vecGammaCandidates,minB,maxB);
    for(Int_t iVec=0; iVec<(Int_t)vecGammaCandidates.size(); iVec++){
        TH1D* temp = vecGammaCandidates.at(iVec);
        temp->Sumw2();
        temp->Scale(1./nEvents[iVec]);
        SetXRange(temp,minB,maxB);
    }
    DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                vecGammaCandidates,"","Number of GammaCandidates","#frac{1}{N_{Events}} #frac{dN}{dNCand}",1,1.1,
                                labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                0.82,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
    SaveCanvas(canvas, Form("%s/Comparison/GammaCandidates.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);

    DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                    vecGammaCandidates,"","Number of GammaCandidates","#frac{1}{N_{Events}} #frac{dN}{dNCand}",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                    0.82,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
    SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_GammaCandidates.%s", outputDir.Data(), suffix.Data()));

    //-------------------------------------------------------------------------------------------------------------------------------
    //---------------------------------------- Meson histogram comparisons ----------------------------------------------------------
    //-------------------------------------------------------------------------------------------------------------------------------
    if ( !isMergedCalo ){
        //-------------------------------------------------------------------------------------------------------------------------------
        // subtracted signal pi0
        for(Int_t iVec=0; iVec<(Int_t)signalPi0.size(); iVec++){
            TH1D* temp = signalPi0.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
        }
        AdjustHistRange(signalPi0,1.1,1.1,kTRUE,1,0);
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                    signalPi0,"","M_{#gamma#gamma}","#frac{1}{N} #frac{1}{N} #frac{dN}{dM_{#gamma#gamma}}",1,1.1,
                                    labelData, colorCompare, kFALSE, 1, 1, kFALSE,
                                    0.82,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
        SaveCanvas(canvas, Form("%s/Comparison/Pi0.%s", outputDir.Data(), suffix.Data()));

        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                        signalPi0,"","M_{#gamma#gamma}","#frac{1}{N} #frac{1}{N} #frac{dN}{dM_{#gamma#gamma}}",1,1.1,
                                        labelData, colorCompare, kFALSE, 1, 1, kTRUE,
                                        0.82,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Pi0.%s", outputDir.Data(), suffix.Data()));
        //-------------------------------------------------------------------------------------------------------------------------------
        // subtracted signal eta
        for(Int_t iVec=0; iVec<(Int_t)signalEta.size(); iVec++){
            TH1D* temp = signalEta.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
        }
        AdjustHistRange(signalEta,1.1,1.1,kTRUE,1,0);
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                    signalEta,"","M_{#gamma#gamma}","#frac{1}{N} #frac{1}{N} #frac{dN}{dM_{#gamma#gamma}}",1,1.1,
                                    labelData, colorCompare, kFALSE, 1, 1, kFALSE,
                                    0.82,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
        SaveCanvas(canvas, Form("%s/Comparison/Eta.%s", outputDir.Data(), suffix.Data()));

        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                        signalEta,"","M_{#gamma#gamma}","#frac{1}{N} #frac{1}{N} #frac{dN}{dM_{#gamma#gamma}}",1,1.1,
                                        labelData, colorCompare, kFALSE, 1, 1, kTRUE,
                                        0.82,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Eta.%s", outputDir.Data(), suffix.Data()));
    }    
  
    //-------------------------------------------------------------------------------------------------------------------------------
    // pt distribution pi0
    GetMinMaxBin(vecPi0Pt,minB,maxB);
    for(Int_t iVec=0; iVec<(Int_t)vecPi0Pt.size(); iVec++){
        TH1D* temp = vecPi0Pt.at(iVec);
        temp->Sumw2();
        temp->Scale(1./temp->Integral());
        //temp->Scale(1./nEvents[iVec]);
        SetXRange(temp,minB,maxB);
    }
    TGaxis::SetExponentOffset(-0.03, -0.04, "x");
    DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kTRUE,kTRUE,kFALSE,
                                vecPi0Pt,"","#it{p}_{T}^{#pi^{0}}","#frac{1}{N} #frac{dN^{#pi^{0}}}{d#it{p}_{T}}",1,1.1,
                                labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                0.82,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
    SaveCanvas(canvas, Form("%s/Comparison/Pi0_Pt.%s", outputDir.Data(), suffix.Data()), kTRUE, kTRUE);
    TGaxis::SetExponentOffset(0, 0, "x");

    DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                    vecPi0Pt,"","#it{p}_{T}^{#pi^{0}}","#frac{1}{N} #frac{dN^{#pi^{0}}}{d#it{p}_{T}}",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                    0.82,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
    SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Pi0_Pt.%s", outputDir.Data(), suffix.Data()));
    //-------------------------------------------------------------------------------------------------------------------------------
    // alpha distribution Pi0
    GetMinMaxBin(vecPi0Alpha,minB,maxB);
    for(Int_t iVec=0; iVec<(Int_t)vecPi0Alpha.size(); iVec++){
        TH1D* temp = vecPi0Alpha.at(iVec);
        temp->Sumw2();
        temp->Scale(1./temp->Integral());
        //temp->Scale(1./nEvents[iVec]);
        SetXRange(temp,minB,maxB);
    }
    TGaxis::SetExponentOffset(-0.03, -0.04, "x");
    DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                vecPi0Alpha,"","#alpha^{#pi^{0}}","#frac{1}{N} #frac{dN^{#pi^{0}}}{d#alpha}",1,1.1,
                                labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                0.82,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
    SaveCanvas(canvas, Form("%s/Comparison/Pi0_Alpha.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
    TGaxis::SetExponentOffset(0, 0, "x");

    DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecPi0Alpha,"","#alpha^{#pi^{0}}","#frac{1}{N} #frac{dN^{#pi^{0}}}{d#alpha}",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                    0.82,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
    SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Pi0_Alpha.%s", outputDir.Data(), suffix.Data()));
    //-------------------------------------------------------------------------------------------------------------------------------
    // y distribution Pi0
    GetMinMaxBin(vecPi0Y,minB,maxB);
    for(Int_t iVec=0; iVec<(Int_t)vecPi0Y.size(); iVec++){
        TH1D* temp = vecPi0Y.at(iVec);
        temp->Sumw2();
        temp->Scale(1./temp->Integral());
        //temp->Scale(1./nEvents[iVec]);
        SetXRange(temp,minB,maxB);
    }
    TGaxis::SetExponentOffset(-0.03, -0.04, "x");
    DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                vecPi0Y,"","Y^{#pi^{0}}","#frac{1}{N} #frac{dN^{#pi^{0}}}{dY}",1,1.1,
                                labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                0.82,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
    SaveCanvas(canvas, Form("%s/Comparison/Pi0_Y.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
    TGaxis::SetExponentOffset(0, 0, "x");

    DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                    vecPi0Y,"","Y^{#pi^{0}}","#frac{1}{N} #frac{dN^{#pi^{0}}}{dY}",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                    0.82,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
    SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Pi0_Y.%s", outputDir.Data(), suffix.Data()));
    //-------------------------------------------------------------------------------------------------------------------------------
    // opening angle Pi0
    GetMinMaxBin(vecPi0OpenAngle,minB,maxB);
    for(Int_t iVec=0; iVec<(Int_t)vecPi0OpenAngle.size(); iVec++){
        TH1D* temp = vecPi0OpenAngle.at(iVec);
        temp->Sumw2();
        temp->Scale(1./temp->Integral());
        //temp->Scale(1./nEvents[iVec]);
        SetXRange(temp,minB,maxB);
    }
    TGaxis::SetExponentOffset(-0.03, -0.04, "x");
    DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                vecPi0OpenAngle,"","#theta_{#pi^{0}}","#frac{1}{N} #frac{dN^{#pi^{0}}}{d#theta}",1,1.1,
                                labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                0.82,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
    SaveCanvas(canvas, Form("%s/Comparison/Pi0_OpenAngle.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
    TGaxis::SetExponentOffset(0, 0, "x");

    DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                    vecPi0OpenAngle,"","#theta_{#pi^{0}}","#frac{1}{N} #frac{dN^{#pi^{0}}}{d#theta}",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                    0.82,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
    SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Pi0_OpenAngle.%s", outputDir.Data(), suffix.Data()));
    
    if (!isMergedCalo){
        //-------------------------------------------------------------------------------------------------------------------------------
        // pt distribution eta
        GetMinMaxBin(vecEtaPt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecEtaPt.size(); iVec++){
            TH1D* temp = vecEtaPt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kTRUE,kTRUE,kFALSE,
                                    vecEtaPt,"","#it{p}_{T}^{#eta}","#frac{1}{N} #frac{dN^{#eta}}{d#it{p}_{T}}",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.82,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
        SaveCanvas(canvas, Form("%s/Comparison/Eta_Pt.%s", outputDir.Data(), suffix.Data()), kTRUE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");

        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                        vecEtaPt,"","#it{p}_{T}^{#eta}","#frac{1}{N} #frac{dN^{#eta}}{d#it{p}_{T}}",1,1.1,
                                        labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                        0.82,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Eta_Pt.%s", outputDir.Data(), suffix.Data()));
        //-------------------------------------------------------------------------------------------------------------------------------
        // alpha distribution eta
        GetMinMaxBin(vecEtaAlpha,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecEtaAlpha.size(); iVec++){
            TH1D* temp = vecEtaAlpha.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecEtaAlpha,"","#alpha^{#eta}","#frac{1}{N} #frac{dN^{#eta}}{d#alpha}",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.82,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
        SaveCanvas(canvas, Form("%s/Comparison/Eta_Alpha.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");

        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                        vecEtaAlpha,"","#alpha^{#eta}","#frac{1}{N} #frac{dN^{#eta}}{d#alpha}",1,1.1,
                                        labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                        0.82,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Eta_Alpha.%s", outputDir.Data(), suffix.Data()));
        //-------------------------------------------------------------------------------------------------------------------------------
        // y distribution eta
        GetMinMaxBin(vecEtaY,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecEtaY.size(); iVec++){
            TH1D* temp = vecEtaY.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecEtaY,"","Y^{#eta}","#frac{1}{N} #frac{dN^{#eta}}{dY}",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.82,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
        SaveCanvas(canvas, Form("%s/Comparison/Eta_Y.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");

        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                        vecEtaY,"","Y^{#eta}","#frac{1}{N} #frac{dN^{#eta}}{dY}",1,1.1,
                                        labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                        0.82,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Eta_Y.%s", outputDir.Data(), suffix.Data()));
        //-------------------------------------------------------------------------------------------------------------------------------
        // opening angle distribution eta
        GetMinMaxBin(vecEtaOpenAngle,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecEtaOpenAngle.size(); iVec++){
            TH1D* temp = vecEtaOpenAngle.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecEtaOpenAngle,"","#theta_{#eta}","#frac{1}{N} #frac{dN^{#eta}}{d#theta}",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.82,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
        SaveCanvas(canvas, Form("%s/Comparison/Eta_OpenAngle.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");

        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                        vecEtaOpenAngle,"","#theta_{#eta}","#frac{1}{N} #frac{dN^{#eta}}{d#theta}",1,1.1,
                                        labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                        0.82,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Eta_OpenAngle.%s", outputDir.Data(), suffix.Data()));
    }
    
    
    //*****************************************************************************************************
    //*****************************************************************************************************
    //****************************** Special Histograms ***************************************************
    //*****************************************************************************************************
    //*****************************************************************************************************

    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "Drawing Special Histograms" << endl;

    if (isConv){
        for(Int_t i=0; i<nSets; i++){
            temp.clear();
            if((Int_t)vecInvMassBeforeAfter[0].size()>0 && (Int_t)vecInvMassBeforeAfter[1].size()>0){
                temp.push_back(vecInvMassBeforeAfter[0].at(i));
                temp.push_back(vecInvMassBeforeAfter[1].at(i));

                AdjustHistRange(temp,5,5,kFALSE);

                DrawGammaSetMarker(vecInvMassBeforeAfter[0].at(i),20,0.8, kBlack , kBlack);
                DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                    vecInvMassBeforeAfter[0].at(i),"","M_{e^{+}e^{-}} (GeV/#it{c}^{2})","#frac{d#it{N}_{#gamma}}{#it{N}_{#gamma} d#it{M}_{e^{+}e^{-}}} (GeV/#it{c}^{2})^{-1}",1,1,
                                    0.82,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
                DrawGammaSetMarker(vecInvMassBeforeAfter[1].at(i),20,0.8, kRed , kRed);
                vecInvMassBeforeAfter[1].at(i)->Draw("e,histsame");

                TLegend *legAdditional = new TLegend(0.12,0.9,0.8,0.94);
                legAdditional->SetNColumns(2);
                legAdditional->SetLineColor(0);
                legAdditional->SetTextSize(0.03);
                legAdditional->SetFillColor(0);
                legAdditional->AddEntry(vecInvMassBeforeAfter[0].at(i),("before cuts"));
                legAdditional->AddEntry(vecInvMassBeforeAfter[1].at(i),("after cuts"));
                legAdditional->Draw();

                SaveCanvas(canvas, Form("%s/Comparison/Photon_InvMass_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()),kFALSE,kTRUE);
                delete legAdditional;
            }
        }

        DeleteVecTH1D(vecInvMassBeforeAfter[0]);
        DeleteVecTH1D(vecInvMassBeforeAfter[1]);
    }
    


    //*****************************************************************************************************
    //***************************** Cleanup vectors *******************************************************
    //*****************************************************************************************************
    //*****************************************************************************************************

    fLog.close();

    delete cvsQuadratic;
    delete canvas;

    DeleteVecTH1D(vecNEvents);
    DeleteVecTH1D(vecVertexZ);
    DeleteVecTH1D(vecNGoodTracks);
    DeleteVecTH1D(vecGammaCandidates);
    DeleteVecTH1D(vecV0Mult);

    DeleteVecTH1D(signalPi0);
    DeleteVecTH1D(signalEta);

    DeleteVecTH1D(vecPi0Pt);
    DeleteVecTH1D(vecPi0Alpha);
    DeleteVecTH1D(vecPi0Y);
    DeleteVecTH1D(vecPi0OpenAngle);

    DeleteVecTH1D(vecEtaPt);
    DeleteVecTH1D(vecEtaAlpha);
    DeleteVecTH1D(vecEtaY);
    DeleteVecTH1D(vecEtaOpenAngle);

    DeleteVecTH2D(vecESDMother);
    DeleteVecTH2D(vecESDBackground);

    delete[] nEvents;
    delete[] nEventsAll;

    delete[] fTrigger;
    delete[] fCutSelection;
    delete[] fEventCutSelection;
    delete[] fGammaCutSelection;
    delete[] fClusterCutSelection;
    delete[] fElectronCutSelection;
    delete[] fMesonCutSelection;

    TH1::AddDirectory(kTRUE);

    cout << "Done with EventQA" << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

    return;

}//end

