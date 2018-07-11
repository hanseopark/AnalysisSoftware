/*******************************************************************************
 ******  provided by Gamma Conversion Group, PWGGA,                        *****
 ******     Daniel Muehlheim, d.muehlheim@cern.ch                          *****
 ******     Nicolas Schmidt, nicolas.schmidt@cern.ch                       *****
 ******     Friederike Bock, fbock@cern.ch                                 *****
 *******************************************************************************/

#include "QA.h"
#include "../TaskQA/BuildHistogramsTriggerQA.C"

void TriggerQA_Runwise(
                        Int_t nSetsIn,
                        Int_t nDataIn,
                        TString fEnergyFlag,
                        TString filePath,
                        TString fileName,
                        TString* DataSets,
                        TString* plotDataSets,
                        Int_t mode                      = 2,
                        Int_t cutNr                     = 0,               // if -1: you have to choose number at runtime, if 0 the first or only number is chosen
                        Int_t doExtQA                   = 2,                // 0: switched off, 1: normal extQA, 2: with Cell level plots
                        Bool_t doEquidistantXaxis       = kFALSE,
                        Bool_t doTrigger                = kTRUE,
                        Bool_t doHistsForEverySet       = kTRUE,
                        Bool_t addSubFolder             = kFALSE,
                        Bool_t useDataRunListForMC      = kFALSE,
                        Size_t markerSize               = 1,
                        TString suffix                  = "eps",
                        TString folderRunlists          = "",
                        TString addLabelRunList         = "",
                        TString cutTreeProjection       = ""
                    ){
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "TriggerQA_Runwise" << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

    gROOT->Reset();
    TH1::AddDirectory(kFALSE);

    const Int_t nSets = nSetsIn;
    const Int_t nData = nDataIn;

    StyleSettingsThesis();
    SetPlotStyle();
    //**************************************************************************************************************
    TString fDate = ReturnDateString();
    TString fTextMeasurement = Form("#pi^{0} #rightarrow #gamma#gamma");

    const Int_t maxSets = 20;
    if(nSets>maxSets){
        cout << "Maximum hardcoded number of Data Sets: " << maxSets << endl;
        cout << "You have chosen: " << nSets << ", returning!" << endl;
        return;
    }

    std::vector<TString> vecDataSet;
    Style_t hMarkerStyle[maxSets];
    Size_t hMarkerSize[maxSets];
    Color_t hMarkerColor[maxSets];
    Color_t hLineColor[maxSets];

//**************************************************************************************************************
    for(Int_t i=0; i<nSets; i++){
        vecDataSet.push_back(DataSets[i].Data());
        hMarkerSize[i]=markerSize;
    }

    Int_t fMode = mode;
    // mode:    0 // new output PCM-PCM
    //          1 // new output PCM dalitz
    //          2 // new output PCM-EMCal
    //          3 // new output PCM-PHOS
    //          4 // new output EMCal-EMCal
    //          5 // new output PHOS-PHOS
    //          9 // old output PCM-PCM
    //          10 // merged EMCal
    //          11 // merged PHOS
    if(fMode == 4 || fMode == 5 || fMode == 10 || fMode == 11){ cout << "Returning, given mode contains no PCM information: " << fMode << endl; return;}

    std::vector<TString> vecRuns;
    TString fileRuns[maxSets];
    for(Int_t i=0; i<nSets; i++){
        fileRuns[i]             = Form("%s/runNumbers%s%s.txt", folderRunlists.Data(), (vecDataSet.at(i)).Data(),addLabelRunList.Data());
        if(useDataRunListForMC && i>=nData) {
            fileRuns[i]         = Form("%s/runNumbers%s%s-%s.txt", folderRunlists.Data(), vecDataSet.at(i).Data(), addLabelRunList.Data(),vecDataSet.at(0).Data());
            cout << "Switch useDataRunListForMC is true, reading runs from: " << fileRuns[i].Data() << endl;
        }
        cout << "trying to read: " << fileRuns[i].Data() << endl;
        if(!readin(fileRuns[i], vecRuns, kFALSE)) {cout << "ERROR, no Run Numbers could be found! Returning..." << endl; return;}
    }

    //******************************************************************************
    cout << "creating: " << Form("%s/%s/%s/%s", filePath.Data(), ((TString)vecDataSet.at(0)).Data(), ((TString)vecRuns.at(0)).Data(),fileName.Data()) << endl;
    TFile* fTriggerQAFile = new TFile(Form("%s/%s/%s/%s", filePath.Data(), ((TString)vecDataSet.at(0)).Data(), ((TString)vecRuns.at(0)).Data(),fileName.Data()));
    if(fTriggerQAFile->IsZombie()) {
        cout << "ERROR: ROOT file '" << Form("%s/%s/%s/%s", filePath.Data(), ((TString)vecDataSet.at(0)).Data(), ((TString)vecRuns.at(0)).Data(), fileName.Data()) << "' could not be openend, testing stages output!" << endl;
        delete fTriggerQAFile;

        // check wether there are subfiles
        TString listForTesting  = Form("%s/%s/%s/Stage*/*/%s", filePath.Data(), ((TString)vecDataSet.at(0)).Data(), ((TString)vecRuns.at(0)).Data(), fileName.Data());
        TString listName        = Form("fileList%s.txt", ((TString)vecRuns.at(0)).Data());
        gSystem->Exec("ls "+listForTesting+ " > " + listName);

        ifstream fileListName;
        fileListName.open(listName,ios_base::in);
        if (!fileListName) {
            cout << "ERROR: settings " << listName.Data() << " not found!" << endl;
            return;
        }

        Bool_t hasStagedFile    = kFALSE ;
        for( TString tempLine; tempLine.ReadLine(fileListName, kTRUE) && !hasStagedFile ; ) {
            cout << tempLine.Data() << endl;
            fTriggerQAFile = new TFile(tempLine.Data(),"READ");
            if(!fTriggerQAFile->IsZombie()) {
                hasStagedFile = kTRUE;
                cout << "found corresponding file in Stages output" << endl;
            }
        }
        if (!hasStagedFile){
            cout << "could not find any merges of Stages either, return!" << endl;
            return;
        }
    }

    TString nameCutsPQA;        // full cut number
    TString nameCutsPQAshort;   // event cut number
    vector <TString> cutsPQA;   // in case there are several different event cut selections (eg centralities) saved in one AnalysisResults.root file
    vector <TString> cutsPQAshort;
    TString nameMainDir;
    TKey *keyPQA;
    TIter nextPQA(fTriggerQAFile->GetListOfKeys());
    while ((keyPQA=(TKey*)nextPQA())){
      cout << Form("Found TopDir: '%s' ",keyPQA->GetName()) << endl;
        nameMainDir = keyPQA->GetName();
        if(nameMainDir.Contains("GammaConvV1_QA_")){
            nameCutsPQA = Form("%s",nameMainDir.Data());
            nameCutsPQAshort = Form("%s",nameMainDir.Data());
            nameCutsPQA.Replace(0,15,"");
            nameCutsPQAshort.Replace(0,15,"");
            nameCutsPQAshort.Replace(8,35,"");
            cutsPQA.push_back(nameCutsPQA);
            cutsPQAshort.push_back(nameCutsPQAshort);
        }
    }
    cout << "The following cuts are available:" << endl;
    for(Int_t i = 0; i < (Int_t) cutsPQA.size(); i++) {
      cout << Form("(%i) -- %s", i, cutsPQA[i].Data()) << endl;
    }
    if(cutNr == -1){
        do{ cin >> cutNr;}
        while( (cutNr < 0) || (cutNr > (Int_t) cutsPQA.size()) );
    }
    cout << "Processing Cut Number: " << cutNr << endl;
    nameCutsPQA      = cutsPQA.at(cutNr);
    nameCutsPQAshort = cutsPQAshort.at(cutNr);

    TString fCentralityFromCut = GetCentralityString(nameCutsPQAshort);
    cout << "Corresponding to centrality: " << fCentralityFromCut << endl;
    // plotting settings according to default for energy, dataset and centrality
    for(Int_t i=0; i<nSets; i++){
        hMarkerStyle[i]=GetDefaultMarkerStyle(fEnergyFlag.Data(),DataSets[i].Data(),fCentralityFromCut.Data());
        hMarkerColor[i]=GetColorDefaultColor(fEnergyFlag.Data(),DataSets[i].Data(),fCentralityFromCut.Data());
        hLineColor[i]=GetColorDefaultColor(fEnergyFlag.Data(),DataSets[i].Data(),fCentralityFromCut.Data());
    }

    Bool_t cutTreeStdCut = kFALSE;
    if (cutTreeProjection.CompareTo("") == 0 ){
        cutTreeProjection   = "0004314141";
        cutTreeStdCut       = kTRUE;
    }
    //choosing cut to project from tree
    if(fEnergyFlag.CompareTo("8TeV")==0){
        cutTreeProjection   = "0005314140";
        cutTreeStdCut       = kTRUE;
    } else if(fEnergyFlag.CompareTo("PbPb_2.76TeV")==0){
        cutTreeProjection   = "0005310040";
        cutTreeStdCut       = kTRUE;
    }
    // if(fEnergyFlag.CompareTo("pPb_5.023TeV")==0){cutTreeProjection = "0000004140"; cutTreeStdCut = kFALSE;} //without chi2 and psipair

    cout << endl;
    cout << "long cutnumber for TriggerQA: '" << nameCutsPQA << "'\nevent cutnumber for TriggerQA: '"<< nameCutsPQAshort << "'" << endl;
    cout << "****************************************************************************" << endl;
    cout << "****************************************************************************" << endl;
    cout << "****************************************************************************\n" << endl;
    if (cutTreeStdCut) cout << "WARNING: TriggerQA hist will be build with a preset cut, make sure its what you wanted!!" << endl;
    cout << "using cutnumber: '" << cutTreeProjection << "' to project out of tree, double check with BuildHistogramsTriggerQA.C!!!" << endl;
    cout << "\n****************************************************************************" << endl;
    cout << "****************************************************************************" << endl;
    cout << "****************************************************************************" << endl;

    fTriggerQAFile->Close();
    delete fTriggerQAFile;
    if(nameMainDir.IsNull() || !nameMainDir.BeginsWith("Gamma")){cout << "ERROR, Unable to obtain valid name of MainDir:|" << nameMainDir.Data() << "|, running in mode: " << fMode << endl; return;}

    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "Processing Cut Number: " << nameCutsPQA << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "Pictures are saved as " << suffix.Data() << "!" << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

    TString fCollisionSystem = ReturnFullCollisionsSystem(fEnergyFlag);
    if (fCollisionSystem.CompareTo("") == 0){
        cout << "No correct collision system specification, has been given" << endl;
        return;
    }
    if (fCentralityFromCut.CompareTo("pp")!=0 && !fCentralityFromCut.Contains("0-100%") ){
      fCollisionSystem    = Form("%s %s", fCentralityFromCut.Data(), fCollisionSystem.Data());
    }

    TString fDetectionProcess = ReturnFullTextReconstructionProcess(fMode);

    TString outputDir = Form("%s_%s/%s/TriggerQA/%s/Runwise",nameCutsPQA.Data(),cutTreeProjection.Data(), fEnergyFlag.Data(),suffix.Data());
    if(addSubFolder){
        outputDir += "/";
        outputDir += DataSets[0];
    }
    gSystem->Exec("mkdir -p "+outputDir);

//****************************** Histograms ************************************************

    std::vector<TString> globalRuns;

    for(Int_t i=0; i<nSets; i++)
    {
        vecRuns.clear();
        if(!readin(fileRuns[i], vecRuns, kFALSE)) {cout << "ERROR, no Run Numbers could be found! Returning..." << endl; return;}

        for(Int_t j=0; j<(Int_t) vecRuns.size();j++)
        {
            if( i==0 ) globalRuns.push_back(vecRuns.at(j));
            else
            {
                Bool_t bFound = kFALSE;
                for(Int_t k=0; k<(Int_t) globalRuns.size();k++){ if(globalRuns.at(k)==vecRuns.at(j)) bFound=kTRUE;}
                if(!bFound) globalRuns.push_back(vecRuns.at(j));
            }
        }
    }

    selection_sort(globalRuns.begin(), globalRuns.end());

    map<TString,Int_t> mapBin;

    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "Processing following list of " << globalRuns.size() << " Runs:";
    for(Int_t i=0; i<(Int_t) globalRuns.size(); i++) {
        mapBin[globalRuns.at(i)]=i+1;
        if(i%10==0) cout << endl;
        cout << globalRuns.at(i) << ", ";
    }
    cout << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;


    TH1D* hGammaN[nSets];
    TH1D* hGammaPt[nSets];

    TH1D* hGammaPtRMS[nSets];
    std::vector<TH1D*>* vecGammaPt = new std::vector<TH1D*>[nSets];

    Int_t hFBin;
    Int_t hLBin;
    Int_t hNBin;

    if(doEquidistantXaxis) {
        hFBin = 0;
        hLBin = globalRuns.size();
        hNBin = globalRuns.size();
    } else {
        if(nSets>2){
            hFBin = ((TString) globalRuns.at(0)).Atoi() - 500;
            hLBin = ((TString) globalRuns.at(globalRuns.size()-1)).Atoi()  + 500;
            hNBin = hLBin - hFBin;
        }else{
            hFBin = ((TString) globalRuns.at(0)).Atoi() - 25;
            hLBin = ((TString) globalRuns.at(globalRuns.size()-1)).Atoi()  + 25;
            hNBin = hLBin - hFBin;
        }
    }

    std::vector<TH1D*>* vecHistos = new std::vector<TH1D*>[nSets];
    std::vector<TString> vecHistosName;
    TString histoName;

    for(Int_t i=0; i<nSets; i++){
        histoName = "nGamma";
        if(i==0) vecHistosName.push_back(histoName);
        hGammaN[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"nGamma; Run Number ; Number of Photons",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hGammaN[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hGammaN[i]);

        histoName = "hGammaPt";
        if(i==0) vecHistosName.push_back(histoName);
        hGammaPt[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hGammaPt; Run Number ; Mean Photon #it{p}_{T}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hGammaPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hGammaPt[i]);

        histoName = "hGammaPtRMS";
        if(i==0) vecHistosName.push_back(histoName);
        hGammaPtRMS[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hGammaPtRMS; Run Number ; #sigma_{Photon #it{p}_{T}}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hGammaPtRMS[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hGammaPtRMS[i]);
    }

//****************************** Looping over DataSets ************************************************

    fstream* fLog = new fstream[nSets];
    for(Int_t iStr=0; iStr<nSets; iStr++){
        if(useDataRunListForMC && iStr>=nData) fLog[iStr].open(Form("%s/A-%s-%s.log",outputDir.Data(),DataSets[iStr].Data(),DataSets[0].Data()), ios::out);
        else fLog[iStr].open(Form("%s/A-%s.log",outputDir.Data(),DataSets[iStr].Data()), ios::out);
        fLog[iStr] << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        fLog[iStr] << DataSets[iStr].Data() << endl;
        fLog[iStr] << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        fLog[iStr] << fCollisionSystem.Data() << endl;
        fLog[iStr] << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        fLog[iStr] << "processed cut: " << nameCutsPQA.Data() << endl;
        fLog[iStr] << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    }

    TString fRootFile;
    TString fDataSet;
    TString fRunNumber;
    Int_t bin = -1;

    std::vector<TString>* vecMissingRuns = new std::vector<TString>[nSets];

    for(Int_t i=0; i<nSets; i++) {
        vecRuns.clear();
        fDataSet = vecDataSet.at(i);
        if(!readin(fileRuns[i], vecRuns)) {cout << "ERROR, no Run Numbers could be found! Returning..." << endl; return;}

        Bool_t doMergeOutput = kFALSE;
        TString mergeCommand = Form("hadd -f -k %s/%s/TriggerQA_%s.root",filePath.Data(),fDataSet.Data(),fDataSet.Data());
        //check if file exists, otherwise merge
        TFile* mergedFile = new TFile(Form("%s/%s/TriggerQA_%s.root",filePath.Data(),fDataSet.Data(),fDataSet.Data()),"READ");
        if(mergedFile->IsZombie()){doMergeOutput = kTRUE; cout << "--- Switched on merging of Runs! ---" << endl;}
        mergedFile->Close();
        delete mergedFile;

        //****************************** Looping over Runs ************************************************
        cout << endl;
        cout << "\t----------------------------------------------------------------------------" << endl;
        cout << "\tLooping over Runs of DataSet |" << (vecDataSet.at(i)).Data() << "|" << endl;
        cout << "\t----------------------------------------------------------------------------" << endl;
        cout << endl;
        fLog[i] << "Looping over Runs:" << endl;

        vecMissingRuns[i].clear();
        for(Int_t j=0; j<(Int_t) vecRuns.size(); j++) {
            fRunNumber = vecRuns.at(j);
            fRootFile = Form("%s/%s/%s/TriggerQA_Data.root", filePath.Data(), fDataSet.Data(), fRunNumber.Data());

            cout << endl;
            cout << "\t\t----------------------------------------------------------------------------" << endl;
            cout << Form("\t\tRun %s", fRunNumber.Data()) << endl;
            cout << "\t\tProcessing file: " << fRootFile.Data() << endl;
            cout << "\t\t----------------------------------------------------------------------------" << endl;
            cout << endl;

            Bool_t isTriggerQAFileOpen   = kTRUE;
            Bool_t needToBuildQAFile    = kFALSE;
            TFile* fTriggerQAFile = new TFile(fRootFile.Data(),"READ");
            if(fTriggerQAFile->IsZombie()){
                cout << endl;
                cout << "\t\t----------------------------------------------------------------------------" << endl;
                cout << "\t\t" << fRootFile << ", could not be opened." << endl;
                cout << endl;
                isTriggerQAFileOpen      = kFALSE;
                fTriggerQAFile->Close();
                delete fTriggerQAFile;

                needToBuildQAFile       = kTRUE;

            } else {
                TDirectory* TopDir = (TDirectory*) fTriggerQAFile->Get(Form("GammaConvV1_QA_%s_%s",nameCutsPQA.Data(),cutTreeProjection.Data()));
                if(TopDir == NULL) {
                    cout << "ERROR: TopDir not Found"<<endl;
                    cout << endl;
                    cout << "\t\t----------------------------------------------------------------------------" << endl;
                    cout << endl;
                    isTriggerQAFileOpen      = kFALSE;
                    fTriggerQAFile->Close();
                    delete fTriggerQAFile;
                    needToBuildQAFile       = kTRUE;
                }
            }

            if (needToBuildQAFile){
                TString fAnalysisResultsFile = Form("%s/%s/%s/%s", filePath.Data(), fDataSet.Data(), fRunNumber.Data(), fileName.Data());
                cout << "\t\tTrying to open " << fAnalysisResultsFile << "..." << endl;
                TFile* RootFileTriggerQA = new TFile(fAnalysisResultsFile.Data(),"READ");
                if(RootFileTriggerQA->IsZombie()) {
                    vecMissingRuns[i].push_back(fRunNumber);
                    cout << "INFO: ROOT file '" << fAnalysisResultsFile.Data() << "' could not be openend, continue!" << endl;
                    cout << "\t\t----------------------------------------------------------------------------" << endl;
                    RootFileTriggerQA->Close();
                    delete RootFileTriggerQA;

                    // check wether there are subfiles
                    TString listForTesting  = Form("%s/%s/%s/Stage*/*/%s", filePath.Data(), fDataSet.Data(), fRunNumber.Data(), fileName.Data());
                    TString listName        = Form("fileList%s.txt", fRunNumber.Data());
                    gSystem->Exec("ls "+listForTesting+ " > " + listName);

                    TString fileNamesInList[100];
                    Int_t nFilesInList                  = 0;

                    ifstream fileListName;
                    fileListName.open(listName,ios_base::in);
                    if (!fileListName) {
                        cout << "ERROR: settings " << listName.Data() << " not found!" << endl;
                        return;
                    }

                    for( TString tempLine; tempLine.ReadLine(fileListName, kTRUE); ) {
                        cout << tempLine.Data() << endl;
                        TFile* RootFileTriggerQA = new TFile(tempLine.Data(),"READ");
                        if(!RootFileTriggerQA->IsZombie()) {
                            fileNamesInList[nFilesInList]   = tempLine;
                            nFilesInList++;
                        }
                        RootFileTriggerQA->Close();
                        delete RootFileTriggerQA;
                    }

                    for (Int_t readFile = 0; readFile < nFilesInList; readFile++ ){
                        doMergeOutput = kTRUE;
                        if (readFile == 0){
                            cout << "\t\tCalling BuildHistogramsTriggerQA for " << fileNamesInList[readFile].Data() << endl;
                            BuildHistogramsTriggerQA(fileNamesInList[readFile].Data(),nameCutsPQA,cutTreeProjection,(i>=nData),0,Form("%s/%s/%s/TriggerQA",filePath.Data(), fDataSet.Data(), fRunNumber.Data()),kTRUE);
                            cout << "\t\tdone!" << endl;
                            cout << "\t\t----------------------------------------------------------------------------" << endl;
                        } else {
                            cout << "\t\tCalling BuildHistogramsTriggerQA... to add " << fileNamesInList[readFile].Data() << endl;
                            BuildHistogramsTriggerQA(fileNamesInList[readFile].Data(),nameCutsPQA,cutTreeProjection,(i>=nData),1,Form("%s/%s/%s/TriggerQA",filePath.Data(), fDataSet.Data(), fRunNumber.Data()),kTRUE);
                            cout << "\t\tdone!" << endl;
                            cout << "\t\t----------------------------------------------------------------------------" << endl;
                        }
                    }
                } else {
                    RootFileTriggerQA->Close();
                    delete RootFileTriggerQA;
                    cout << "\t\tCalling BuildHistogramsTriggerQA..." << endl;
                    BuildHistogramsTriggerQA(fAnalysisResultsFile.Data(),nameCutsPQA,cutTreeProjection,(i>=nData),0,Form("%s/%s/%s/TriggerQA",filePath.Data(), fDataSet.Data(), fRunNumber.Data()),kTRUE);
                    cout << "\t\tdone!" << endl;
                    cout << "\t\t----------------------------------------------------------------------------" << endl;
                    doMergeOutput = kTRUE;
                    cout << "--- Switched on merging of Runs! ---" << endl;
                }
            }

            mergeCommand += Form(" %s",fRootFile.Data());
            if(!isTriggerQAFileOpen){
                fTriggerQAFile = new TFile(fRootFile.Data(),"READ");
                if(fTriggerQAFile->IsZombie()) {vecMissingRuns[i].push_back(fRunNumber); cout << "INFO: ROOT file '" << fRootFile.Data() << "' could not be openend, continue!" << endl; fTriggerQAFile->Close(); delete fTriggerQAFile; continue;}
            }

            TDirectory* TopDir = (TDirectory*) fTriggerQAFile->Get(Form("GammaConvV1_QA_%s_%s",nameCutsPQA.Data(), cutTreeProjection.Data()));
                if(TopDir == NULL) {cout << "ERROR: TopDir not Found"<<endl; return;}

            //---------
            if(doEquidistantXaxis) bin = mapBin[fRunNumber];
            else bin = fRunNumber.Atoi() - hFBin;

            Double_t nEvents = 0;
            TH1D* GOODESD = (TH1D*) TopDir->Get("histoGoodESDTracks");
            if(GOODESD){
            nEvents = (Double_t) GOODESD->Integral();
            delete GOODESD;
            }

            TH2D* GammaEtaPt = (TH2D*) TopDir->Get("histoGammaEtaPt");
            Double_t nGammas = 1.;
            if(GammaEtaPt){
            nGammas = GammaEtaPt->Integral();
            hGammaN[i]->SetBinContent(bin, nGammas);
            }
            //---------
            fLog[i] << "----------------------------------------------------------------------------" << endl;
            fLog[i] << "Processing file: " << fRootFile.Data() << endl;
            fLog[i] << "Run: " << fRunNumber.Data() << ", with NEvents: " << nEvents << endl;
            fLog[i] << "----------------------------------------------------------------------------" << endl;
            //---------
            if( nEvents <= 1. ){cout << "Warning: number of accepted events in run: " << nEvents << "! Setting nEvents to 1..." << endl; nEvents = 1.;}
            if( nGammas <= 1. ){cout << "Warning: number of gammas: " << nGammas << "! Setting nGammas to 1..." << endl; nGammas = 1.;}
            //---------
    //---------
            TH1D* GammaPhi = (TH1D*) TopDir->Get("histoGammaPhi");
            if(GammaPhi){
            hGammaPhi[i]->SetBinContent(bin, GammaPhi->GetMean());
            hGammaPhi[i]->SetBinError(bin, GammaPhi->GetMeanError());
            hGammaPhiRMS[i]->SetBinContent(bin, GammaPhi->GetRMS());
            hGammaPhiRMS[i]->SetBinError(bin, GammaPhi->GetRMSError());
            hGammaPhiLow[i]->SetBinContent(bin, GetHistogramIntegral(GammaPhi,0,TMath::Pi())/GammaPhi->Integral());
            hGammaPhiLow[i]->SetBinError(bin, GetHistogramIntegralError(GammaPhi,0,TMath::Pi())/GammaPhi->Integral());
            hGammaPhiHigh[i]->SetBinContent(bin, GetHistogramIntegral(GammaPhi,TMath::Pi(),2*TMath::Pi())/GammaPhi->Integral());
            hGammaPhiHigh[i]->SetBinError(bin, GetHistogramIntegralError(GammaPhi,TMath::Pi(),2*TMath::Pi())/GammaPhi->Integral());
            GammaPhi->GetXaxis()->SetTitle("Photon #phi");
            GammaPhi->GetYaxis()->SetTitle("#frac{1}{N_{#gamma}} #frac{dN}{d#phi}");
            GammaPhi->Sumw2();
            GammaPhi->Scale(1 / nGammas);
            vecGammaPhi[i].push_back(GammaPhi);
            }else cout << "INFO: Object |histoGammaPhi| could not be found! Skipping Fill..." << endl;

            TopDir->Clear();
            delete TopDir;

            fTriggerQAFile->Close();
            delete fTriggerQAFile;
        }

        if(doMergeOutput) {
            cout << "Flag doMergeOutput = kTRUE, merging TriggerQA output..." << endl;
            gSystem->Exec(mergeCommand.Data());
            cout << "done!" << endl;
        }
    }

//****************************** Drawing Histograms ************************************************
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "Drawing Histograms" << endl;

    TCanvas* canvas2D = new TCanvas("canvas2D","",10,10,750*6,500*6);  // gives the page size
    Double_t leftMar = 0.09; Double_t rightMar = 0.025; Double_t topMargin = 0.04; Double_t bottomMargin = 0.09;
    DrawGammaCanvasSettings(canvas2D, leftMar, rightMar, topMargin, bottomMargin);
    TCanvas* canvas = new TCanvas("canvas","",10,10,750,500);  // gives the page size
    DrawGammaCanvasSettings(canvas, leftMar, rightMar, topMargin, bottomMargin);

    Float_t xPosLabel = 0.8;
    if(fEnergyFlag.Contains("PbPb")) xPosLabel = 0.75;

    if(doHistsForEverySet) {
        TBox *boxLabel = new TBox(1.37,0.7,1.78,0.83);
        boxLabel->SetFillStyle(0);boxLabel->SetFillColor(0);boxLabel->SetLineColor(1);boxLabel->SetLineWidth(0.6);
        TBox *boxLabel2 = new TBox(-0.4,51,6.5,56.5);
        boxLabel2->SetFillStyle(0);boxLabel2->SetFillColor(0);boxLabel2->SetLineColor(1);boxLabel2->SetLineWidth(0.6);

        TString outputDirDataSet;

        for(Int_t i=0; i<nSets; i++) {
            cout << "DataSet: " << DataSets[i].Data() << endl;
            outputDirDataSet = Form("%s/%s",outputDir.Data(), DataSets[i].Data());
            if(useDataRunListForMC && i>=nData) {
                outputDirDataSet = Form("%s/%s-%s", outputDir.Data(), DataSets[i].Data(),DataSets[0].Data());
                cout << "Switch useDataRunListForMC is true, output to: " << outputDirDataSet.Data() << endl;
            }
            gSystem->Exec("mkdir -p "+outputDirDataSet);

            //--------
            for(Int_t h=0; h<(Int_t) vecHistos[i].size(); h++){
                (((TH1D*) vecHistos[i].at(h)))->SetTitle("");
                AdjustHistRange(((TH1D*) vecHistos[i].at(h)),1.1,1.1,kTRUE);
                if(((TString)vecHistosName.at(h)).CompareTo("nGamma")==0) AdjustHistRange(((TH1D*) vecHistos[i].at(h)),10,10,kTRUE);
                ((TH1D*) vecHistos[i].at(h))->Draw("px0e1");

                if(doTrigger) PutProcessLabelAndEnergyOnPlot(xPosLabel, 0.92, 0.03, fCollisionSystem.Data(), plotDataSets[i].Data(), fTrigger.Data());
                else PutProcessLabelAndEnergyOnPlot(xPosLabel, 0.92, 0.03, fCollisionSystem.Data(), plotDataSets[i].Data(), "");

                if(((TString)vecHistosName.at(h)).CompareTo("nGamma")==0) SaveCanvas(canvas, Form("%s/%s.%s", outputDirDataSet.Data(), vecHistosName.at(h).Data(), suffix.Data()), kFALSE, kTRUE);
                else SaveCanvas(canvas, Form("%s/%s.%s", outputDirDataSet.Data(), vecHistosName.at(h).Data(), suffix.Data()));
            }
        }
        delete boxLabel;
        delete boxLabel2;
//****************************** Drawing Runwise Histograms ************************************************
        cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        cout << "Drawing Runwise Histograms" << endl;

        TCanvas* canvasRunwise;
        TLegend *legendRuns;

        for(Int_t i=0; i<nSets; i++) {
            vecRuns.clear();
            fDataSet = vecDataSet.at(i);
            outputDirDataSet = Form("%s/%s",outputDir.Data(), DataSets[i].Data());
            if(useDataRunListForMC && i>=nData) {
                outputDirDataSet = Form("%s/%s-%s", outputDir.Data(), DataSets[i].Data(),DataSets[0].Data());
            }
            if(!readin(fileRuns[i], vecRuns, kFALSE)) {cout << "ERROR, no Run Numbers could be found! Returning..." << endl; return;}

            for(Int_t iRun=0; iRun<(Int_t)vecMissingRuns[i].size(); iRun++){
                vecRuns.erase(std::remove(vecRuns.begin(), vecRuns.end(), vecMissingRuns[i].at(iRun)), vecRuns.end());
            }

            Int_t NColumns = ((Int_t) vecRuns.size() / 31 ) + 1;

            canvasRunwise = new TCanvas("canvasRunwise","",10,10,1350+(NColumns*108),900);  // gives the page size
            DrawGammaCanvasSettings(canvasRunwise, 130.5/(1350.+(NColumns*108.)), (40.5+(NColumns*108.))/(1350.+(NColumns*108.)), topMargin, bottomMargin);
            canvasRunwise->cd();

            Double_t addRight = ((Double_t)NColumns*108)/(1350+((Double_t)NColumns*108));
            legendRuns = new TLegend(0.98-addRight,0.09,0.98,0.94);
            legendRuns->SetNColumns(NColumns);
            legendRuns->SetFillColor(0);
            legendRuns->SetLineColor(0);
            legendRuns->SetTextSize(0.03);
            legendRuns->SetTextFont(42);

            cout << "DataSet: " << plotDataSets[i].Data() << endl;
            TString fTrigger = "";
            if(i<nData){
                TString fTriggerCut = nameCutsPQAshort(3,2);
                fTrigger = AnalyseSpecialTriggerCut(fTriggerCut.Atoi(), plotDataSets[i]);
                if(fTrigger.Contains("not defined")) fTrigger = "";
            }

            //---------
            DrawVectorRunwiseTH1D(	canvasRunwise, legendRuns, vecGammaPt[i], vecRuns, 5, 5, kTRUE, addRight, xPosLabel, 0.94, 0.03, xPosLabel, 0.8, 0.03,
                                    doTrigger, fTrigger, (Bool_t)(i<nData), outputDirDataSet, "GammaPt", plotDataSets[i], kFALSE,
                                    fCollisionSystem, "", suffix, kTRUE, kTRUE, kFALSE);
            //--------
            delete legendRuns;
            delete canvasRunwise;
        }
        canvas->cd();
    }

//****************************** Combined Trending Histograms ************************************************
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "Drawing Trending Histograms" << endl;
//---------
    if(useDataRunListForMC) cout << "WARNING: useDataRunListForMC is true, overwriting histograms for DataSet!" << endl;

    TLegend *legend = new TLegend(0.15,0.95,0.95,0.98);
    legend->SetNColumns(nSets);
    legend->SetFillColor(0);
    legend->SetLineColor(0);
    legend->SetTextSize(0.03);
    legend->SetTextFont(42);
//---------
    if(nSets > 1)
    {
        for(Int_t h=0; h<(Int_t) vecHistos[0].size(); h++)
        {
            TString fTrigger = "";
            TString fTriggerCut = nameCutsPQAshort(3,2);
            for(Int_t i=0; i<nSets; i++)
            {
                if(doTrigger && i<nData && i==0){
                    fTrigger = AnalyseSpecialTriggerCut(fTriggerCut.Atoi(), plotDataSets[i]);
                    if(fTrigger.Contains("not defined")) fTrigger = "";
                }
                legend->AddEntry(((TH1D*) vecHistos[i].at(h)),plotDataSets[i].Data(),"p");
            }
            AdjustHistRange(vecHistos,1.1,1.1,h,nSets,kTRUE);
            if(((TString)vecHistosName.at(h)).CompareTo("hGammaN")==0) AdjustHistRange(vecHistos,10,10,h,nSets,kTRUE);
            for(Int_t i=nSets-1; i>=0; i--)
            {
                TString draw = (i==nSets-1)?"px0e1":"px0e1, same";
                ((TH1D*) vecHistos[i].at(h))->SetTitle("");
                ((TH1D*) vecHistos[i].at(h))->Draw(draw.Data());
            }
            legend->Draw();

            PutProcessLabelAndEnergyOnPlot(xPosLabel, 0.92, 0.03, fCollisionSystem.Data(), fTrigger.Data(), "");

            if(canvas->GetTopMargin()!=0.06) canvas->SetTopMargin(0.06);
            if(useDataRunListForMC && !addSubFolder){
                if( ((TString)vecHistosName.at(h)).CompareTo("nGamma")==0 ) SaveCanvas(canvas, Form("%s/%s/%s.%s", outputDir.Data(), DataSets[0].Data(),vecHistosName.at(h).Data(),suffix.Data()), kFALSE, kTRUE);
                else SaveCanvas(canvas, Form("%s/%s/%s.%s", outputDir.Data(), DataSets[0].Data(),vecHistosName.at(h).Data(),suffix.Data()));
            }else{
                if(((TString)vecHistosName.at(h)).CompareTo("nGamma")==0) SaveCanvas(canvas, Form("%s/%s.%s", outputDir.Data(), vecHistosName.at(h).Data(), suffix.Data()), kFALSE, kTRUE);
                else SaveCanvas(canvas, Form("%s/%s.%s", outputDir.Data(), vecHistosName.at(h).Data(), suffix.Data()));
            }
            legend->Clear();
        }
    }

//****************************** Combined Ratio Trending Histograms ************************************************
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "Drawing Ratio Trending Histograms" << endl;

    if(doHistsForEverySet)
    {
        TString* ratioSets[nSets-nData];
        TH1D* ratio[nSets-nData];
        TString outputDirDataSet;

        if(nSets>1 && nSets>nData)
        {
            legend->SetNColumns(nData*(nSets-nData));
            Int_t markerStyles[14]={2,4,5,20,21,22,23,24,25,26,27,28,29,30};

            for(Int_t h=0; h<(Int_t) vecHistos[0].size(); h++)
            {
                for(Int_t i=0; i<nData; i++)
                {
                    for(Int_t j=0; j<nSets-nData; j++)
                    {
                        ratioSets[j] = new TString(Form("%s / %s", plotDataSets[i].Data(), plotDataSets[j+nData].Data()));
                        ratio[j] = new TH1D(Form("%s%i%i",((TH1D*) vecHistos[i].at(h))->GetName(),i,j),
                                            Form("%s%i%i;%s;%s - Ratio: Data / MC",((TH1D*) vecHistos[i].at(h))->GetTitle(),i,j,((TH1D*) vecHistos[i].at(h))->GetXaxis()->GetTitle(),((TH1D*) vecHistos[i].at(h))->GetYaxis()->GetTitle()),
                                            hNBin,hFBin,hLBin);
                        EditTH1(globalRuns, doEquidistantXaxis, ratio[j], markerStyles[j % 14], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
                    }

                    for(Int_t b=1; b<hNBin+1; b++)
                    {
                        Double_t hData = ((TH1D*) vecHistos[i].at(h))->GetBinContent(b);
                        for(Int_t j=0; j<nSets-nData; j++)
                        {
                            Double_t hMC = ((TH1D*) vecHistos[j+nData].at(h))->GetBinContent(b);
                            if(hMC!=0) {
                                if(hData/hMC>1.98) ratio[j]->SetBinContent(b,1.98);
                                else if(hData/hMC<0.02) ratio[j]->SetBinContent(b,0.02);
                                else ratio[j]->SetBinContent(b,hData/hMC);
                            }else ratio[j]->SetBinContent(b,1.98);
                        }
                    }

                    for(Int_t j=0; j<nSets-nData; j++)
                    {
                        TString draw = (i==0&&j==0)?"p":"p, same";
                        ratio[j]->SetTitle("");
                        ratio[j]->GetYaxis()->SetRangeUser(0,2);
                        ratio[j]->Draw(draw.Data());
                        legend->AddEntry(ratio[j],ratioSets[j]->Data(),"p");
                    }

                    legend->Draw();
                    outputDirDataSet = Form("%s/TrendingRatios",outputDir.Data());
                    gSystem->Exec("mkdir -p "+outputDirDataSet);

                    if(doTrigger){
                    TString fTrigger = "";
                    TString fTriggerCut = nameCutsPQAshort(3,2);
                    fTrigger = AnalyseSpecialTriggerCut(fTriggerCut.Atoi(), plotDataSets[i]);
                    if(fTrigger.Contains("not defined")) fTrigger = "";

                    PutProcessLabelAndEnergyOnPlot(xPosLabel, 0.92, 0.03, fCollisionSystem.Data(), fTrigger.Data(), "");
                    }else PutProcessLabelAndEnergyOnPlot(xPosLabel, 0.92, 0.03, fCollisionSystem.Data(), "", "");

                    SaveCanvas(canvas, Form("%s/%s.%s", outputDirDataSet.Data(),Form("%s",((TH1D*) vecHistos[i].at(h))->GetName()),suffix.Data()));
                    legend->Clear();
                    for(Int_t j=0; j<nSets-nData; j++)
                    {
                    delete ratio[j];
                    delete ratioSets[j];
                    }
                }
            }
        }else cout << "...skipped due to nSets<=1 or nSets==nData!" << endl;

    }
    delete legend;
    delete canvas;

//****************************** Create Output ROOT-File ************************************************

    const char* nameOutput;

    for(Int_t i=0; i<nSets; i++){
        fDataSet = vecDataSet.at(i);

        if(useDataRunListForMC && i>=nData) nameOutput = Form("%s_%s/%s/TriggerQA/%s-%s_TriggerQARunwise.root",nameCutsPQA.Data(), cutTreeProjection.Data() ,fEnergyFlag.Data(),fDataSet.Data(),vecDataSet.at(0).Data());
        else nameOutput = Form("%s_%s/%s/TriggerQA/%s_TriggerQARunwise.root",nameCutsPQA.Data(),cutTreeProjection.Data(),fEnergyFlag.Data(),fDataSet.Data());

        TFile* fOutput = new TFile(nameOutput,"RECREATE");
        cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        cout << "Output file: " << nameOutput << endl;
        cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        fLog[i] << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        fLog[i] << "Output file: " << nameOutput << endl;
        fLog[i] << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

        for(Int_t h=0; h<(Int_t) vecHistos[i].size(); h++) WriteHistogram(((TH1D*) vecHistos[i].at(h)));

        DeleteVecTH1D(vecGammaPt[i]);

        fOutput->Write();
        fOutput->Close();
        delete fOutput;
        fLog[i].close();
    }

    delete[] vecHistos;
    delete[] vecGammaPt;
    delete[] vecMissingRuns;
    delete[] fLog;

    TH1::AddDirectory(kTRUE);

    cout << "Done with TriggerQA_Runwise" << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

    return;

}//end
