/*******************************************************************************
******  provided by Gamma Conversion Group, PWGGA,                        *****
******     Daniel Muehlheim, d.muehlheim@cern.ch                          *****
******     Adrian Mechler, mechler@ikf.uni-frnakfurt.de                   *****
*******************************************************************************/

#include "QA.h"

void ClusterQA_Runwise_V2(
    Int_t nSetsIn,                                  // number of sets to be analysed
    Int_t nDataIn,                                  // number of real data sets to be analysed
    TString fEnergyFlag,                            // energy flag
    TString filePath,                               // path to the data
    TString fileName,                               // file name of the data
    TString fileNameMC,
    TString* DataSets,                              // technical names of data sets for output
    TString* plotDataSets,                          // labels of data sets in plots
    Int_t mode                      = 2,            // standard mode for analysis
    Int_t cutNr                     = -1,           // if -1: you have to choose number at runtime
    Int_t doExtQA                   = 2,            // 0: switched off, 1: normal extQA, 2: with Cell level plots
    Bool_t doEquidistantXaxis       = kFALSE,       // kTRUE: each run in runlist corresponds to 1 bin in X in histogram,
    // kFALSE: histograms contain the complete specified run number range, where each run represents a bin - even if it is not specified
    Bool_t doTrigger                = kTRUE,        // enables trigger analysis
    Bool_t doHistsForEverySet       = kTRUE,        // kTRUE: output done for each set separately as well
    // kFALSE: only full run range output is produced
    Bool_t addSubFolder             = kFALSE,       // kTRUE: adds another subfolder for QA output fo reach DataSet[i]
    // kFALSE: stores the runwise output all together
    Bool_t useDataRunListForMC      = kFALSE,       // kTRUE: use the same run list for data and MC
    // kFALSE: use specified
    Size_t markerSize               = 1,            // how large should the markers be?
    TString suffix                  = "eps",        // output format of plots
    TString folderRunlists          = "",           // path to the runlists
    TString addLabelRunList         = "",           // additional name for runlist
    Bool_t runMergedClust           = kFALSE,        // flag to switch to merged cluster cuts instead of standard cluster cut
    Int_t NMaxRuns                  = 50,           // Maximum of Runs to be analysed at once (depending on availible RAM)
    Bool_t onlytrending                  = kFALSE           // Maximum of Runs to be analysed at once (depending on availible RAM)
)
{
    cout << endl << endl << endl << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<< endl;
    cout << "ClusterQA_Runwise_V2("<< endl;
    cout << "nSetsIn:               \t" << nSetsIn << endl;
    cout << "nDataIn:               \t" << nDataIn << endl;
    cout << "fEnergyFlag:           \t" << fEnergyFlag << endl;
    cout << "filePath:              \t" << filePath << endl;
    cout << "fileName:              \t" << fileName << endl;
    cout << "fileNameMC:            \t" << fileNameMC << endl;
    cout << "DataSets:              \t" << DataSets << endl;
    cout << "plotDataSets:          \t" << plotDataSets << endl;
    cout << "mode:                  \t" << mode << endl;
    cout << "cutNr:                 \t" << cutNr << endl;
    cout << "doExtQA:               \t" << doExtQA << endl;
    cout << "doEquidistantXaxis:    \t" << doEquidistantXaxis << endl;
    cout << "kFALSE:                \t" << kFALSE << endl;
    cout << "doTrigger:             \t" << doTrigger << endl;
    cout << "doHistsForEverySet:    \t" << doHistsForEverySet << endl;
    cout << "kFALSE:                \t" << kFALSE << endl;
    cout << "addSubFolder:          \t" << addSubFolder << endl;
    cout << "kFALSE:                \t" << kFALSE << endl;
    cout << "useDataRunListForMC:   \t" << useDataRunListForMC << endl;
    cout << "kFALSE:                \t" << kFALSE << endl;
    cout << "markerSize:            \t" << markerSize << endl;
    cout << "suffix:                \t" << suffix << endl;
    cout << "folderRunlists:        \t" << folderRunlists << endl;
    cout << "addLabelRunList:       \t" << addLabelRunList << endl;
    cout << "runMergedClust:        \t" << runMergedClust << endl;
    cout << "NMaxRuns:              \t" << NMaxRuns << endl;
    cout << "onlytrending:              \t" << NMaxRuns << endl;
    cout << ")"<< endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<< endl;

    // **************************************************************************************************************
    // **************************** Setting general plotting style **************************************************
    // **************************************************************************************************************

    gROOT->Reset();
    TH1::AddDirectory(kFALSE);
    StyleSettingsThesis();
    SetPlotStyle();

    // **************************************************************************************************************
    // ****************************** Setting common variables ******************************************************
    // **************************************************************************************************************


    TString fileNameData        = fileName;

    TString fDate               = ReturnDateString();
    TString fTextMeasurement    = Form("#pi^{0} #rightarrow #gamma#gamma");
    TString nameMainDir, nameMainDirData, nameMainDirMC;

    Float_t xPosLabel           = 0.8;
    if (fEnergyFlag.CompareTo("pPb_5.023TeV") == 0)
    xPosLabel               = 0.75;


    Int_t fMode                 = mode;
    Bool_t isPCMCalo            = kTRUE;
    Bool_t isMerged             = kFALSE;
    // mode:    0 // new output PCM-PCM
    //          1 // new output PCM dalitz
    //          2 // new output PCM-EMCal
    //          3 // new output PCM-PHOS
    //          4 // new output EMCal-EMCal
    //          5 // new output PHOS-PHOS
    //          9 // old output PCM-PCM
    //          10 // merged EMCal
    //          11 // merged EMCal
    //          14 // new output PCM-EDC (EMCal + DCal)
    //          15 // new output EDC (EMCal + DCal)
    if(fMode == 0 || fMode == 1 || fMode == 9){ printf("Returning, given mode contains no calo information: %i\n", fMode); return;}
    if(fMode != 2 && fMode != 3 && fMode != 14)
    isPCMCalo               = kFALSE;
    if(fMode == 10 || fMode == 11)
    isMerged                = kTRUE;

    TString fCollisionSystem = ReturnFullCollisionsSystem(fEnergyFlag);
    if (fCollisionSystem.CompareTo("") == 0){
        printf("No correct collision system specification, has been given\n");
        return;
    }
    TString fDetectionProcess = ReturnFullTextReconstructionProcess(fMode);


    // *************************************************************************************************************
    // runNumbers
    std::vector<TString> vecRuns;
    std::vector<TString> fileRuns;
    std::vector<TString> vecDataSet;
    std::vector<TString> vecDataSetIn;
    std::vector<TString> vecplotDataSets;
    std::vector<Int_t> vecSubRunlists;
    std::vector<Int_t> vecNSplits;

    Int_t nSetsInTmp = 0;
    Int_t nDataInTmp = 0;
    //
    printf("\n\n-----------------------------------------------\n");
    printf("Setup: %i Sets, %i Data, %i MC\n#\ttype\tName\tSplits\tNRuns\n", nSetsIn, nDataIn, nSetsIn-nDataIn);
    for(Int_t i=0; i<nSetsIn; i++){
        TString fileRunstmp             = Form("%s/runNumbers%s%s.txt", folderRunlists.Data(), (DataSets[i]).Data(),addLabelRunList.Data());
        if(useDataRunListForMC && i>=nDataIn) {
            fileRunstmp         = Form("%s/runNumbers%s%s-%s.txt", folderRunlists.Data(), DataSets[i].Data(), addLabelRunList.Data(),DataSets[0].Data());
            printf("Switch useDataRunListForMC is true, reading runs from: %s\n", fileRunstmp.Data());
        }
        // printf("trying to read: " << fileRunstmp.Data() << endl;
        vecRuns.clear();
        if(!readin(fileRunstmp, vecRuns, kFALSE)) {printf("\033[0;31mERROR\033[0m, no Run Numbers could be found! Returning...\n"); return;}
        Int_t NRuns=vecRuns.size();
        Int_t NSplits = (Int_t)ceil(((Double_t)(NRuns))/((Double_t)NMaxRuns));
        nSetsInTmp += (NSplits-1);
        TString type = "MC";
        if(i < nDataIn) {
            nDataInTmp += (NSplits-1);
            type = "Data";
        }
        printf("%i\t%s\t%s\t%i\t%i\n",i,type.Data(), DataSets[i].Data(), NSplits, NRuns);

        for (Int_t iSubRunlists = 0; iSubRunlists < NSplits; iSubRunlists++) {
            vecNSplits.push_back(NSplits);
            vecSubRunlists.push_back(iSubRunlists);
            if( NSplits > 1 ) {
                vecDataSet.push_back(Form("%s_%i",DataSets[i].Data(),iSubRunlists));
            } else {
                vecDataSet.push_back(Form("%s",DataSets[i].Data()));
            }
            vecplotDataSets.push_back(Form("%s",plotDataSets[i].Data()));
            vecDataSetIn.push_back(Form("%s",DataSets[i].Data()));
            fileRuns.push_back(fileRunstmp);
        }
    }
    // const Int_t nSets           = nSetsIn + nSetsInTmp;
    if( nSetsIn + nSetsInTmp != (Int_t)vecDataSet.size()) printf("#######################################################\n");
    const Int_t nSets           = (Int_t)vecDataSet.size();
    const Int_t nData           = nDataIn + nDataInTmp;
    printf("\n->\t%i Sets, %i Data, %i MC\n", nSets, nData, nSets-nData);
    printf("-----------------------------------------------\n\n");





    // *************************************************************************************************************
    // ****************************** Determine which cut to process ***********************************************
    // *************************************************************************************************************
    UInt_t ActualRunIndexInVector=0;
    TFile* fCutFile             = NULL;
    vecRuns.clear();
    if(!readin(fileRuns.at(0), vecRuns, kFALSE)) {printf("\033[0;31mERROR\033[0m, no Run Numbers could be found! Returning...\n"); return;}
    while ( (fCutFile==NULL)&&(ActualRunIndexInVector<vecRuns.size()) ){
        fCutFile             = new TFile(Form("%s/%s/%s/%s", filePath.Data(), ((TString)vecDataSetIn.at(0)).Data(), ((TString)vecRuns.at(ActualRunIndexInVector)).Data(), fileName.Data()));
        if(fCutFile->IsZombie()) {
            printf("\033[0;31mERROR\033[0m: ROOT file '%s/%s/%s/%s' could not be openend, return!\n", filePath.Data(), ((TString)vecDataSetIn.at(0)).Data(), ((TString)vecRuns.at(ActualRunIndexInVector)).Data(), fileName.Data());
            fCutFile->Close();
            delete fCutFile;
            fCutFile=NULL;
        }
        ActualRunIndexInVector++;
    }
    if(fCutFile==NULL) {
        printf("\033[0;31mERROR\033[0m: no ROOT file; last tried file: '%s/%s/%s/%s'; return!\n", filePath.Data(), ((TString)vecDataSet.at(0)).Data(), ((TString)vecRuns.at(ActualRunIndexInVector-1)).Data(), fileName.Data());
        return;
    }
    TKey *key;
    TIter next(fCutFile->GetListOfKeys());
    while ((key=(TKey*)next())){
        cout << Form("Found TopDir: '%s' ",key->GetName())<<endl;
        nameMainDir             = key->GetName();
        if (nameMainDir.Contains("Gamma")) break;

    }
    nameMainDirData=nameMainDir;
    if (fileNameData.CompareTo(fileNameMC.Data())==0){
        nameMainDirMC=nameMainDirData;
    }
    else {
        cout<<"fileNameData and fileNameMC differ => Changing corresponding name variable nameMainDirMC"<<endl;
        ActualRunIndexInVector=0;
        TFile* fCutFileMC=NULL;
        while ( (fCutFileMC==NULL)&&(ActualRunIndexInVector<vecRuns.size()) ){
            fCutFileMC             = new TFile(Form("%s/%s/%s/%s", filePath.Data(), ((TString)vecDataSetIn.at(nData)).Data(), ((TString)vecRuns.at(ActualRunIndexInVector)).Data(), fileNameMC.Data()));
            if(fCutFileMC->IsZombie()) {
                printf("\033[0;31mERROR\033[0m: MC ROOT file %s/%s/%s/%s' could not be openend, return!\n", filePath.Data(), ((TString)vecDataSetIn.at(nData)).Data(), ((TString)vecRuns.at(ActualRunIndexInVector)).Data(), fileNameMC.Data());
                fCutFileMC->Close();
                delete fCutFileMC;
                fCutFileMC=NULL;
            }
            ActualRunIndexInVector++;
        }
        if(fCutFileMC==NULL) {
            printf("\033[0;31mERROR\033[0m: no MC ROOT file; last tried File: %s/%s/%s/%s'; return!\n", filePath.Data(), ((TString)vecDataSet.at(nData)).Data(), ((TString)vecRuns.at(ActualRunIndexInVector-1)).Data(), fileNameMC.Data());
            return;
        }
        TKey *keyMC;
        TIter nextMC(fCutFileMC->GetListOfKeys());
        while ((keyMC=(TKey*)nextMC())){
            cout << Form("Found TopDir for MC: '%s' ",keyMC->GetName())<<endl;
            nameMainDirMC             = keyMC->GetName();
            if (nameMainDirMC.Contains("Gamma")) break;
        }
        cout<<"nameMainDirMC changed to: "<<nameMainDirMC.Data()<<endl;
        delete keyMC;
        fCutFileMC->Close();
        delete fCutFileMC;
    }
    cout << endl;
    printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
    if(nameMainDir.IsNull() || !nameMainDir.BeginsWith("Gamma")){printf("\033[0;31mERROR\033[0m, Unable to obtain valid name of MainDir:|%s|, running in mode: %i\n", nameMainDir.Data(), fMode); return;}

    TList *listInput            = (TList*)fCutFile->Get(nameMainDir.Data());
    listInput->SetOwner(kTRUE);
    if(!listInput) {printf("\033[0;31mERROR\033[0m: Could not find main dir: %s in file! Returning...\n", nameMainDir.Data()); return;}
    vector <TString> cuts;
    for(Int_t i = 0; i<listInput->GetSize(); i++){
        TList *listCuts         = (TList*)listInput->At(i);
        TString nameCuts        = listCuts->GetName();
        if(nameCuts.BeginsWith("Cut Number")){
            nameCuts.Replace(0,11,"");
            cuts.push_back(nameCuts);
        }
    }
    delete listInput;

    printf("The following cuts are available:\n");
    for(Int_t i = 0; i < (Int_t) cuts.size(); i++) {cout << Form("(%i) -- %s", i, cuts[i].Data()) << endl;}

    if(cutNr == -1){
        do{ cin >> cutNr;}
        while( (cutNr < 0) || (cutNr > (Int_t) cuts.size()) );
    }

    fCutFile->Close();
    delete fCutFile;

    printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
    printf("Processing Cut Number: %i\n", cutNr);
    cout << cuts.at(cutNr) << endl;
    printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
    printf("Pictures are saved as %s!\n", suffix.Data());
    printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");

    TString fCutSelection            = cuts.at(cutNr);
    TString fEventCutSelection       = "";
    TString fGammaCutSelection       = "";
    TString fClusterCutSelection     = "";
    TString fMClusterCutSelection    = "";
    TString fElectronCutSelection    = "";
    TString fMesonCutSelection       = "";
    if (!isMerged){
        ReturnSeparatedCutNumberAdvanced(fCutSelection, fEventCutSelection, fGammaCutSelection, fClusterCutSelection, fElectronCutSelection, fMesonCutSelection, fMode);
    } else {
        ReturnSeparatedCutNumberAdvanced(fCutSelection, fEventCutSelection, fClusterCutSelection, fMClusterCutSelection, fElectronCutSelection, fMesonCutSelection, fMode);
    }
    if (runMergedClust){
        fClusterCutSelection        = fMClusterCutSelection;
    }


    // *****************************************************************************************************
    // **************************** Determine whether to run track matching or not *************************
    // *****************************************************************************************************
    Bool_t isTrackMatching          = kTRUE;
    TString trackMatchingCut(fClusterCutSelection(GetClusterTrackMatchingCutPosition(fClusterCutSelection),1));
    Int_t trackMatch                = trackMatchingCut.Atoi();
    if(trackMatch == 0 || runMergedClust ){
        printf("\033[0;33mINFO\033[0m: TrackMatching cut found to be '0' for '%s', deactivating generation of track matching histograms!\n",(vecDataSet.at(0)).Data());
        isTrackMatching             = kFALSE;
    }

    // *****************************************************************************************************
    // ************************** Set proper cluster nomenclature ******************************************
    // *****************************************************************************************************
    TString calo                    = "";
    Int_t iCalo                     = 0;
    Int_t nCaloModules              = 0;
    Int_t nCaloCells                = 0;
    if(fClusterCutSelection.BeginsWith('1')){
        calo                        = "EMCal";
        iCalo                       = 1;
        nCaloModules                = 10;
        nCaloCells                  = 11520;
        if((vecDataSet.at(0)).Contains("LHC10")){
            nCaloModules            = 4;
            nCaloCells              = 4608;
        }
    } else if(fClusterCutSelection.BeginsWith('2')){
        calo                        = "PHOS";
        iCalo                       = 2;
        nCaloModules                = 4;
        nCaloCells                  = 10752;
        if((vecDataSet.at(0)).Contains("LHC15") || (vecDataSet.at(0)).Contains("LHC16") || (vecDataSet.at(0)).Contains("LHC17") || (vecDataSet.at(0)).Contains("LHC18")){
            nCaloModules            = 5;
            nCaloCells              = 14336;
        }
    } else if(fClusterCutSelection.BeginsWith('3')){
        calo                        = "DCal";
        iCalo                       = 3;
        nCaloModules                = 18;
        nCaloCells                  = 18000;
    } else if(fClusterCutSelection.BeginsWith('4')){
        calo                        = "EMCAL+DCAL";
        iCalo                       = 1;
        nCaloModules                = 20;
        nCaloCells                  = 19000;
    } else {printf("No correct calorimeter type found: %s, returning...",calo.Data()); return;}


    // *****************************************************************************************************
    // ************************** Define output directories*************************************************
    // *****************************************************************************************************
    TString outputDir   = "";
    if (!runMergedClust){
        outputDir       = Form("%s/%s/ClusterQA/%s/Runwise",cuts.at(cutNr).Data(),fEnergyFlag.Data(),suffix.Data());
    } else {
        outputDir       = Form("%s/%s/MergedClusterQA/%s/Runwise",cuts.at(cutNr).Data(),fEnergyFlag.Data(),suffix.Data());
    }
    if(addSubFolder)
    outputDir       +=Form("/%s",vecDataSetIn.at(0).Data());
    gSystem->Exec("mkdir -p "+outputDir);

    // *****************************************************************************************************
    // **************************** Determine global run list **********************************************
    // *****************************************************************************************************
    std::vector<TString> globalRuns;

    for(Int_t i=0; i<nSets; i++) {
        vecRuns.clear();
        if(!readin(fileRuns.at(i), vecRuns, kFALSE)) {printf("\033[0;31mERROR\033[0m, no Run Numbers could be found! Returning...\n"); return;}
        for(Int_t j=0; j<(Int_t) vecRuns.size();j++) {

            if( i==0 )
            globalRuns.push_back(vecRuns.at(j));
            else {
                Bool_t bFound = kFALSE;
                for(Int_t k=0; k<(Int_t) globalRuns.size();k++){ if(globalRuns.at(k)==vecRuns.at(j)) bFound=kTRUE;}
                if(!bFound) globalRuns.push_back(vecRuns.at(j));
            }
        }
    }
    selection_sort(globalRuns.begin(), globalRuns.end());

    map<TString,Int_t> mapBin;

    printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
    printf("Processing following list of %i Runs:\n",(Int_t)globalRuns.size());
    for(Int_t i=0; i<(Int_t) globalRuns.size(); i++) {
        mapBin[globalRuns.at(i)]=i+1;
        if(i%10==0) cout << endl;
        cout << globalRuns.at(i) << ", ";
    }
    cout << endl;
    printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");






    std::vector<TH1D*> vecHistos;
    std::vector<TString>  vecHistosName;
    TString histoName;

    Int_t hFBin;
    Int_t hLBin;
    Int_t hNBin;

    if(doEquidistantXaxis) {
        hFBin       = 0;
        hLBin       = globalRuns.size();
        hNBin       = globalRuns.size();
    } else {
        if(nSets>2){
            hFBin   = ((TString) globalRuns.at(0)).Atoi() - 500;
            hLBin   = ((TString) globalRuns.at(globalRuns.size()-1)).Atoi()  + 500;
            hNBin   = hLBin - hFBin;
        } else{
            hFBin   = ((TString) globalRuns.at(0)).Atoi() - 25;
            hLBin   = ((TString) globalRuns.at(globalRuns.size()-1)).Atoi()  + 25;
            hNBin   = hLBin - hFBin;
        }
    }

    // *****************************************************************************************************
    // ******************************* define transverse momentum binning **********************************
    // *****************************************************************************************************

    const Int_t fNBinsClusterPt                 = 60;
    Double_t fBinsClusterPt[fNBinsClusterPt+1]  = { 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.4, 7.8, 8.2, 8.6, 9.0, 9.5, 10, 11, 12, 14, 16, 18, 20, 25, 30};

    TH1D* fDeltaPt              = new TH1D("deltaPt","",fNBinsClusterPt,fBinsClusterPt);
    for(Int_t iPt=1;iPt<fNBinsClusterPt+1;iPt++){
        fDeltaPt->SetBinContent(iPt,fBinsClusterPt[iPt]-fBinsClusterPt[iPt-1]);
        fDeltaPt->SetBinError(iPt,0);
    }

    // *****************************************************************************************************
    // ********************************** Supporting variables *********************************************
    // *****************************************************************************************************
    Double_t minT_Energy[5]     = {0, 0, 1.2, 4, 16};
    Double_t maxT_Energy[5]     = {40, 1.2, 4, 16, 40};

    TString sT[2]               = {"Pi0", "Eta"};
    TString sL[2]               = {"#pi^{0}", "#eta"};

    TString fRootFile;
    TString fRunNumber;
    TString fileCells;
    Int_t bin                   = -1;

    std::vector<TString>* vecMissingRuns = new std::vector<TString>[nSets];

    Int_t width = (Int_t)((Double_t)vecRuns.size()*2.);
    TCanvas* canvas                 = new TCanvas("canvas","",10,10,750,500);  // gives the page size
    TCanvas* canvaswide                 = new TCanvas("canvaswide","",10,10,width,500);  // gives the page size
    TCanvas* canvaslarge                 = new TCanvas("canvaslarge","",10,10,4000,2666);  // gives the page size
    Double_t leftMar                = 0.09;
    Double_t rightMar               = 0.025;
    Double_t topMargin              = 0.04;
    Double_t bottomMargin           = 0.09;
    DrawGammaCanvasSettings(canvas, leftMar, rightMar, topMargin, bottomMargin);
    DrawGammaCanvasSettings(canvaswide, 0.04, 0.01, topMargin, bottomMargin);

    TBox *boxLabel              = new TBox(1.37,0.7,1.78,0.83);
    boxLabel->SetFillStyle(0);boxLabel->SetFillColor(0);boxLabel->SetLineColor(1);boxLabel->SetLineWidth(1);
    TBox *boxLabel2             = new TBox(-0.4,51,6.5,56.5);
    boxLabel2->SetFillStyle(0);boxLabel2->SetFillColor(0);boxLabel2->SetLineColor(1);boxLabel2->SetLineWidth(1);

    TString outputDirDataSet;

    TCanvas* canvasRunwise;
    TLegend *legendRuns;

    const char* nameOutput;
    // /*
    // *****************************************************************************************************
    // ****************************** Looping over DataSets ************************************************
    // *****************************************************************************************************
    for(Int_t i=0; i<nSets; i++) {
        vecRuns.clear();
        if(!readin(fileRuns.at(i), vecRuns, kTRUE, kFALSE, NMaxRuns*vecSubRunlists.at(i), NMaxRuns*(vecSubRunlists.at(i)+1))){printf("\033[0;31mERROR\033[0m, no Run Numbers could be found! Returning...\n"); return;}
        // if(!readin(fileRuns.at(i), vecRuns)) {printf("\033[0;31mERROR\033[0m, no Run Numbers could be found! Returning...\n"); return;}

        // *****************************************************************************************************
        // ******************************* create log files *****************************************************
        // *****************************************************************************************************
        fstream fLog;
        fstream fLogRunwiseHotCells;
        fstream fLogRunwiseColdCells;
        if(useDataRunListForMC && i>=nData) {
            fLog.open(Form("%s/A-%s-%s.log",outputDir.Data(),vecDataSet.at(i).Data(),vecDataSet.at(0).Data()), ios::out);
            if(doExtQA==2) fLogRunwiseHotCells.open(Form("%s/HotCellsRunwise-%s-%s.log",outputDir.Data(),vecDataSet.at(i).Data(),vecDataSet.at(0).Data()), ios::out);
            if(doExtQA==2) fLogRunwiseColdCells.open(Form("%s/ColdCellsRunwise-%s-%s.log",outputDir.Data(),vecDataSet.at(i).Data(),vecDataSet.at(0).Data()), ios::out);
        } else{
            fLog.open(Form("%s/A-%s.log",outputDir.Data(),vecDataSet.at(i).Data()), ios::out);
            if(doExtQA==2) fLogRunwiseHotCells.open(Form("%s/HotCellsRunwise-%s.log",outputDir.Data(),vecDataSet.at(i).Data()), ios::out);
            if(doExtQA==2) fLogRunwiseColdCells.open(Form("%s/ColdCellsRunwise-%s.log",outputDir.Data(),vecDataSet.at(i).Data()), ios::out);
        }
        if( vecNSplits.at(i) > 1 && vecSubRunlists.at(i) == 0) {
            fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
            fLog << vecDataSetIn.at(i).Data() << endl;
            fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
            fLog << fCollisionSystem.Data() << endl;
            fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
            fLog << "processed cut: " << fCutSelection.Data() << endl;
            fLog << calo.Data() << ", Modules: " << nCaloModules << ", Cells: " << nCaloCells << endl;
            fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
            cout << endl << endl << endl;
            printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
            cout << vecDataSetIn.at(i).Data() << endl;
            printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
            cout << fCollisionSystem.Data() << endl;
            printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
            cout << "processed cut: " << fCutSelection.Data() << endl;
            cout << calo.Data() << ", Modules: " << nCaloModules << ", Cells: " << nCaloCells << endl;
            printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
        }

        // **************************************************************************************************************
        // ****************************** Define plotting settings ******************************************************
        // **************************************************************************************************************

        Style_t hMarkerStyle;
        Size_t hMarkerSize;
        Color_t hMarkerColor;
        Color_t hLineColor;

        if (i < nData){
            hMarkerStyle  = GetDefaultMarkerStyle(fEnergyFlag.Data(),vecDataSetIn.at(i).Data(),"");
            hMarkerColor  = GetColorDefaultColor(fEnergyFlag.Data(),vecDataSetIn.at(i).Data(),"",kFALSE,kFALSE);
            hLineColor    = GetColorDefaultColor(fEnergyFlag.Data(),vecDataSetIn.at(i).Data(),"",kFALSE,kFALSE);
            hMarkerSize   = markerSize;
        } else {
            hMarkerStyle  = GetDefaultMarkerStyle(fEnergyFlag.Data(),vecplotDataSets.at(i).Data(),"");
            hMarkerColor  = GetColorDefaultColor(fEnergyFlag.Data(),vecplotDataSets.at(i).Data(),"",kFALSE,kFALSE);
            hLineColor    = GetColorDefaultColor(fEnergyFlag.Data(),vecplotDataSets.at(i).Data(),"",kFALSE,kFALSE);
            hMarkerSize   = markerSize;
        }

        // *****************************************************************************************************
        // ****************************** Vectors for Histograms ***********************************************
        // *****************************************************************************************************

        std::vector<TString> vecBadCells;
        std::vector<TH2D*> vecClusterEtaPhi;
        std::vector<TH2D*> vecCellTimeID;
        std::vector<TH2D*> vecClusterEnergyTime;
        std::vector<TH2D*> vecClusterEVsNCells;
        std::vector<TH1D*> vecClusterIncludedCells;
        std::vector<TH1D*> vecClusterIncludedCellsBefore;
        std::vector<TH1D*> vecClusterEFracCells;
        std::vector<TH1D*> vecClusterEFracCellsBefore;
        std::vector<TH2D*> vecClusterEnergyMeanSigma;
        std::vector<TH2D*> vecClusterTimeMeanSigma;
        std::vector<TH2D*> vecClusterFiredCellIDs;
        // std::vector<TH2D*> vecClusterMissingCellIDs;
        std::vector<TProfile*> vecClusterBadCells;
        std::vector<TH1D*> vecClusterEnergy;
        std::vector<TH1D*> vecClusterM02;
        std::vector<TH1D*> vecClusterM20;
        std::vector<TH1D*> vecClusterNCells;
        std::vector<TH1D*> vecClusterDispersion;
        std::vector<TH1D*> vecClusterTime[5];
        std::vector<TH1D*> vecClusterDeltaEta;
        std::vector<TH1D*> vecClusterDeltaPhi;
        std::vector<TH1D*> vecClusterPi0ConvPhotonEta;
        std::vector<TH1D*> vecClusterPi0ConvPhotonPhi;
        std::vector<TH1D*> vecClusterEtaConvPhotonEta;
        std::vector<TH1D*> vecClusterEtaConvPhotonPhi;
        std::vector<TH1D*> vecBadCellsEnergy;
        std::vector<TH1D*> vecBadCellsTime;
        std::vector<TH1D*> vecClusterEVsModule[nCaloModules];
        std::vector<TH1D*> vecModuleEVsModule[nCaloModules];
        std::vector<TH1D*> vecClusterNCells100VsModule[nCaloModules];
        std::vector<TH1D*> vecClusterNCells1500VsModule[nCaloModules];

        // *****************************************************************************************************
        // ********************* Create histograms for plotting ************************************************
        // *****************************************************************************************************


        histoName                   = "hNEvents";
        if(i==0) vecHistosName.push_back(histoName);
        TH1D* hNEvents                 = new TH1D(Form("%s_%s", histoName.Data(), vecDataSet.at(i).Data()),"hNEvents; Run Number ; #Ge of Events",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hNEvents, hMarkerStyle, hMarkerSize, hMarkerColor, hLineColor);
        vecHistos.push_back(hNEvents);


        histoName                   = "hCaloNClusters";
        if(i==0) vecHistosName.push_back(histoName);
        TH1D* hCaloNClusters           = new TH1D(Form("%s_%s", histoName.Data(), vecDataSet.at(i).Data()),"hCaloNClusters; Run Number ; #frac{1}{N_{Events}} N_{Clusters before Cuts}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hCaloNClusters, hMarkerStyle, hMarkerSize, hMarkerColor, hLineColor);
        vecHistos.push_back(hCaloNClusters);


        histoName                   = "hCaloNClustersQA";
        if(i==0) vecHistosName.push_back(histoName);
        TH1D* hCaloNClustersQA         = new TH1D(Form("%s_%s", histoName.Data(), vecDataSet.at(i).Data()),"hCaloNClustersQA; Run Number ; #frac{1}{N_{Events}} N_{Clusters}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hCaloNClustersQA, hMarkerStyle, hMarkerSize, hMarkerColor, hLineColor);
        vecHistos.push_back(hCaloNClustersQA);


        histoName                   = "hClusterEnergy-Mean";
        if(i==0) vecHistosName.push_back(histoName);
        TH1D* hClusterMeanEnergy       = new TH1D(Form("%s_%s", histoName.Data(), vecDataSet.at(i).Data()),"hClusterMeanEnergy; Run Number ; #bar{#lower[0.1]{E}}_{Cluster} (GeV)",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hClusterMeanEnergy, hMarkerStyle, hMarkerSize, hMarkerColor, hLineColor);
        vecHistos.push_back(hClusterMeanEnergy);


        histoName                   = "hClusterEnergy-RMS";
        if(i==0) vecHistosName.push_back(histoName);
        TH1D* hClusterRMSEnergy        = new TH1D(Form("%s_%s", histoName.Data(), vecDataSet.at(i).Data()),"hClusterRMSEnergy; Run Number ; #sigma_{E_{Cluster}} (GeV)",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hClusterRMSEnergy, hMarkerStyle, hMarkerSize, hMarkerColor, hLineColor);
        vecHistos.push_back(hClusterRMSEnergy);


        histoName                   = "hClusterEnergy-01";
        if(i==0) vecHistosName.push_back(histoName);
        TH1D* hClusterEnergy01         = new TH1D(Form("%s_%s", histoName.Data(), vecDataSet.at(i).Data()),"hClusterEnergy01; Run Number ; #frac{1}{N_{Events}} N_{Clusters, E<1 GeV}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hClusterEnergy01, hMarkerStyle, hMarkerSize, hMarkerColor, hLineColor);
        vecHistos.push_back(hClusterEnergy01);


        histoName                   = "hClusterEnergy-14";
        if(i==0) vecHistosName.push_back(histoName);
        TH1D* hClusterEnergy14         = new TH1D(Form("%s_%s", histoName.Data(), vecDataSet.at(i).Data()),"hClusterEnergy14; Run Number ; #frac{1}{N_{Events}} N_{Clusters, 1<E<4 GeV}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hClusterEnergy14, hMarkerStyle, hMarkerSize, hMarkerColor, hLineColor);
        vecHistos.push_back(hClusterEnergy14);


        histoName                   = "hClusterEnergy-4";
        if(i==0) vecHistosName.push_back(histoName);
        TH1D* hClusterEnergy4          = new TH1D(Form("%s_%s", histoName.Data(), vecDataSet.at(i).Data()),"hClusterEnergy4; Run Number ; #frac{1}{N_{Events}} N_{Clusters, E>4 GeV}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hClusterEnergy4, hMarkerStyle, hMarkerSize, hMarkerColor, hLineColor);
        vecHistos.push_back(hClusterEnergy4);


        histoName                   = "hClusterNCells-Mean";
        if(i==0) vecHistosName.push_back(histoName);
        TH1D* hClusterMeanNCells       = new TH1D(Form("%s_%s", histoName.Data(), vecDataSet.at(i).Data()),"hClusterMeanNCells; Run Number ; #bar{#lower[0.1]{N}}_{Cells, Cluster}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hClusterMeanNCells, hMarkerStyle, hMarkerSize, hMarkerColor, hLineColor);
        vecHistos.push_back(hClusterMeanNCells);


        histoName                   = "hClusterNCells-RMS";
        if(i==0) vecHistosName.push_back(histoName);
        TH1D* hClusterRMSNCells        = new TH1D(Form("%s_%s", histoName.Data(), vecDataSet.at(i).Data()),"hClusterRMSNCells; Run Number ; #sigma_{N_{Cells, Cluster}}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hClusterRMSNCells, hMarkerStyle, hMarkerSize, hMarkerColor, hLineColor);
        vecHistos.push_back(hClusterRMSNCells);


        histoName                   = "hClusterDispersion-Mean";
        if(i==0) vecHistosName.push_back(histoName);
        TH1D* hClusterMeanDispersion   = new TH1D(Form("%s_%s", histoName.Data(), vecDataSet.at(i).Data()),"hClusterMeanDispersion; Run Number ; #bar{#lower[0.1]{Dispersion}}_{Cluster}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hClusterMeanDispersion, hMarkerStyle, hMarkerSize, hMarkerColor, hLineColor);
        vecHistos.push_back(hClusterMeanDispersion);


        histoName                   = "hClusterDispersion-RMS";
        if(i==0) vecHistosName.push_back(histoName);
        TH1D* hClusterRMSDispersion    = new TH1D(Form("%s_%s", histoName.Data(), vecDataSet.at(i).Data()),"hClusterRMSDispersion; Run Number ; #sigma_{Dispersion_{Cluster}}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hClusterRMSDispersion, hMarkerStyle, hMarkerSize, hMarkerColor, hLineColor);
        vecHistos.push_back(hClusterRMSDispersion);


        histoName                   = "hClusterM02-Mean";
        if(i==0) vecHistosName.push_back(histoName);
        TH1D* hClusterMeanM02          = new TH1D(Form("%s_%s", histoName.Data(), vecDataSet.at(i).Data()),"hClusterMeanM02; Run Number ; #bar{#lower[0.1]{#lambda_{0}^{2}}}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hClusterMeanM02, hMarkerStyle, hMarkerSize, hMarkerColor, hLineColor);
        vecHistos.push_back(hClusterMeanM02);


        histoName                   = "hClusterM02-RMS";
        if(i==0) vecHistosName.push_back(histoName);
        TH1D* hClusterRMSM02           = new TH1D(Form("%s_%s", histoName.Data(), vecDataSet.at(i).Data()),"hClusterRMSM02; Run Number ; #sigma_{#lambda_{0}^{2}}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hClusterRMSM02, hMarkerStyle, hMarkerSize, hMarkerColor, hLineColor);
        vecHistos.push_back(hClusterRMSM02);


        histoName                   = "hClusterR-Mean";
        if(i==0) vecHistosName.push_back(histoName);
        TH1D* hClusterMeanR            = new TH1D(Form("%s_%s", histoName.Data(), vecDataSet.at(i).Data()),"hClusterMeanR; Run Number ; #bar{#lower[0.1]{R}}_{Cluster}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hClusterMeanR, hMarkerStyle, hMarkerSize, hMarkerColor, hLineColor);
        vecHistos.push_back(hClusterMeanR);


        histoName                   = "hClusterR-RMS";
        if(i==0) vecHistosName.push_back(histoName);
        TH1D* hClusterRMSR             = new TH1D(Form("%s_%s", histoName.Data(), vecDataSet.at(i).Data()),"hClusterRMSR; Run Number ; #sigma_{R_{Cluster}}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hClusterRMSR, hMarkerStyle, hMarkerSize, hMarkerColor, hLineColor);
        vecHistos.push_back(hClusterRMSR);


        histoName                   = "hClusterTime-Mean";
        if(i==0) vecHistosName.push_back(histoName);
        TH1D* hClusterMeanTime         = new TH1D(Form("%s_%s", histoName.Data(), vecDataSet.at(i).Data()),"hClusterMeanTime; Run Number ; #bar{#lower[0.1]{t}}_{Cluster} (s)",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hClusterMeanTime, hMarkerStyle, hMarkerSize, hMarkerColor, hLineColor);
        vecHistos.push_back(hClusterMeanTime);


        histoName                   = "hClusterTime-RMS";
        if(i==0) vecHistosName.push_back(histoName);
        TH1D* hClusterRMSTime          = new TH1D(Form("%s_%s", histoName.Data(), vecDataSet.at(i).Data()),"hClusterRMSTime; Run Number ; #sigma_{t_{Cluster}} (s)",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hClusterRMSTime, hMarkerStyle, hMarkerSize, hMarkerColor, hLineColor);
        vecHistos.push_back(hClusterRMSTime);


        histoName                   = "hCluster-FractionMatches";
        if(i==0) vecHistosName.push_back(histoName);
        TH1D* hClusterFractionMatches  = new TH1D(Form("%s_%s", histoName.Data(), vecDataSet.at(i).Data()),"hClusterFractionMatches; Run Number ; #frac{N_{Removed Mother Candidates by Cluster-Track Matching}}{N_{Mother Candidates}}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hClusterFractionMatches, hMarkerStyle, hMarkerSize, hMarkerColor, hLineColor);
        vecHistos.push_back(hClusterFractionMatches);


        histoName                   = "hCluster-FractionMatchesS";
        if(i==0) vecHistosName.push_back(histoName);
        TH1D* hClusterFractionMatchesS = new TH1D(Form("%s_%s", histoName.Data(), vecDataSet.at(i).Data()),"hClusterFractionMatchesS; Run Number ; #frac{N_{Removed Mother Candidates by Cluster-Track Matching}}{N_{Mother Candidates}} #scale[0.8]{for #it{p}_{T,Pair}<1.5 GeV/#it{c}}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hClusterFractionMatchesS, hMarkerStyle, hMarkerSize, hMarkerColor, hLineColor);
        vecHistos.push_back(hClusterFractionMatchesS);


        histoName                   = "hCluster-FractionMatchesM";
        if(i==0) vecHistosName.push_back(histoName);
        TH1D* hClusterFractionMatchesM = new TH1D(Form("%s_%s", histoName.Data(), vecDataSet.at(i).Data()),"hClusterFractionMatchesM; Run Number ; #frac{N_{Removed Mother Candidates by Cluster-Track Matching}}{N_{Mother Candidates}} #scale[0.8]{for 1.5<#it{p}_{T,Pair}<2.5 GeV/#it{c}}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hClusterFractionMatchesM, hMarkerStyle, hMarkerSize, hMarkerColor, hLineColor);
        vecHistos.push_back(hClusterFractionMatchesM);


        histoName                   = "hCluster-FractionMatchesH";
        if(i==0) vecHistosName.push_back(histoName);
        TH1D* hClusterFractionMatchesH = new TH1D(Form("%s_%s", histoName.Data(), vecDataSet.at(i).Data()),"hClusterFractionMatchesH; Run Number ; #frac{N_{Removed Mother Candidates by Cluster-Track Matching}}{N_{Mother Candidates}} #scale[0.8]{for #it{p}_{T,Pair}>2.5 GeV/#it{c}}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hClusterFractionMatchesH, hMarkerStyle, hMarkerSize, hMarkerColor, hLineColor);
        vecHistos.push_back(hClusterFractionMatchesH);


        histoName                   = "hClusterDeltaEta-Mean";
        if(i==0) vecHistosName.push_back(histoName);
        TH1D* hClusterMeanDeltaEta     = new TH1D(Form("%s_%s", histoName.Data(), vecDataSet.at(i).Data()),"hClusterMeanDeltaEta; Run Number ; #bar{#lower[0.1]{#Delta#eta}}_{Cluster, Tracks}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hClusterMeanDeltaEta, hMarkerStyle, hMarkerSize, hMarkerColor, hLineColor);
        vecHistos.push_back(hClusterMeanDeltaEta);


        histoName                   = "hClusterDeltaPhi-Mean";
        if(i==0) vecHistosName.push_back(histoName);
        TH1D* hClusterMeanDeltaPhi     = new TH1D(Form("%s_%s", histoName.Data(), vecDataSet.at(i).Data()),"hClusterMeanDeltaPhi; Run Number ; #bar{#lower[0.1]{#Delta#phi}}_{Cluster, Tracks}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hClusterMeanDeltaPhi, hMarkerStyle, hMarkerSize, hMarkerColor, hLineColor);
        vecHistos.push_back(hClusterMeanDeltaPhi);


        histoName                   = "hClusterDeltaEta-RMS";
        if(i==0) vecHistosName.push_back(histoName);
        TH1D* hClusterRMSDeltaEta      = new TH1D(Form("%s_%s", histoName.Data(), vecDataSet.at(i).Data()),"hClusterRMSDeltaEta; Run Number ; #sigma_{#Delta#eta_{Cluster, Tracks}}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hClusterRMSDeltaEta, hMarkerStyle, hMarkerSize, hMarkerColor, hLineColor);
        vecHistos.push_back(hClusterRMSDeltaEta);


        histoName                   = "hClusterDeltaPhi-RMS";
        if(i==0) vecHistosName.push_back(histoName);
        TH1D* hClusterRMSDeltaPhi      = new TH1D(Form("%s_%s", histoName.Data(), vecDataSet.at(i).Data()),"hClusterRMSDeltaPhi; Run Number ; #sigma_{#Delta#phi_{Cluster, Tracks}}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hClusterRMSDeltaPhi, hMarkerStyle, hMarkerSize, hMarkerColor, hLineColor);
        vecHistos.push_back(hClusterRMSDeltaPhi);


        histoName                   = "hClusterM20-Mean";
        if(i==0) vecHistosName.push_back(histoName);
        TH1D* hClusterMeanM20          = new TH1D(Form("%s_%s", histoName.Data(), vecDataSet.at(i).Data()),"hClusterMeanM20; Run Number ; #bar{#lower[0.1]{#lambda_{1}^{2}}}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hClusterMeanM20, hMarkerStyle, hMarkerSize, hMarkerColor, hLineColor);
        vecHistos.push_back(hClusterMeanM20);


        histoName                   = "hClusterM20-RMS";
        if(i==0) vecHistosName.push_back(histoName);
        TH1D* hClusterRMSM20           = new TH1D(Form("%s_%s", histoName.Data(), vecDataSet.at(i).Data()),"hClusterRMSM20; Run Number ; #sigma_{#lower[0.1]{#lambda_{1}^{2}}}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hClusterRMSM20, hMarkerStyle, hMarkerSize, hMarkerColor, hLineColor);
        vecHistos.push_back(hClusterRMSM20);


        histoName                           = "hClusterConvPhotonPi0_Eta-Mean";
        if(i==0) vecHistosName.push_back(histoName);
        TH1D* hClusterMeanConvPhotonPi0_Eta    = new TH1D(Form("%s_%s", histoName.Data(), vecDataSet.at(i).Data()),"hClusterConvPhotonPi0_Eta; Run Number ; #bar{#lower[0.1]{#eta}}_{#gamma_{conv} under #pi^{0}-peak}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hClusterMeanConvPhotonPi0_Eta, hMarkerStyle, hMarkerSize, hMarkerColor, hLineColor);
        vecHistos.push_back(hClusterMeanConvPhotonPi0_Eta);


        histoName                           = "hClusterConvPhotonPi0_Phi-Mean";
        if(i==0) vecHistosName.push_back(histoName);
        TH1D* hClusterMeanConvPhotonPi0_Phi    = new TH1D(Form("%s_%s", histoName.Data(), vecDataSet.at(i).Data()),"hClusterConvPhotonPi0_Phi; Run Number ; #bar{#lower[0.1]{#phi}}_{#gamma_{conv} under #pi^{0}-peak}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hClusterMeanConvPhotonPi0_Phi, hMarkerStyle, hMarkerSize, hMarkerColor, hLineColor);
        vecHistos.push_back(hClusterMeanConvPhotonPi0_Phi);


        histoName                           = "hClusterConvPhotonEta_Eta-Mean";
        if(i==0) vecHistosName.push_back(histoName);
        TH1D* hClusterMeanConvPhotonEta_Eta    = new TH1D(Form("%s_%s", histoName.Data(), vecDataSet.at(i).Data()),"hClusterConvPhotonEta_Eta; Run Number ; #bar{#lower[0.1]{#eta}}_{#gamma_{conv} under #eta-peak}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hClusterMeanConvPhotonEta_Eta, hMarkerStyle, hMarkerSize, hMarkerColor, hLineColor);
        vecHistos.push_back(hClusterMeanConvPhotonEta_Eta);


        histoName                           = "hClusterConvPhotonEta_Phi-Mean";
        if(i==0) vecHistosName.push_back(histoName);
        TH1D* hClusterMeanConvPhotonEta_Phi    = new TH1D(Form("%s_%s", histoName.Data(), vecDataSet.at(i).Data()),"hClusterConvPhotonEta_Phi; Run Number ; #bar{#lower[0.1]{#phi}}_{#gamma_{conv} under #eta-peak}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hClusterMeanConvPhotonEta_Phi, hMarkerStyle, hMarkerSize, hMarkerColor, hLineColor);
        vecHistos.push_back(hClusterMeanConvPhotonEta_Phi);


        histoName                           = "hClusterConvPhotonPi0_Eta-RMS";
        if(i==0) vecHistosName.push_back(histoName);
        TH1D* hClusterRMSConvPhotonPi0_Eta     = new TH1D(Form("%s_%s", histoName.Data(), vecDataSet.at(i).Data()),"hClusterConvPhotonPi0_Eta; Run Number ; #sigma_{#eta_{#gamma_{conv} under #pi^{0}-peak}}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hClusterRMSConvPhotonPi0_Eta, hMarkerStyle, hMarkerSize, hMarkerColor, hLineColor);
        vecHistos.push_back(hClusterRMSConvPhotonPi0_Eta);


        histoName                           = "hClusterConvPhotonPi0_Phi-RMS";
        if(i==0) vecHistosName.push_back(histoName);
        TH1D* hClusterRMSConvPhotonPi0_Phi     = new TH1D(Form("%s_%s", histoName.Data(), vecDataSet.at(i).Data()),"hClusterConvPhotonPi0_Phi; Run Number ; #sigma_{#phi_{#gamma_{conv} under #pi^{0}-peak}}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hClusterRMSConvPhotonPi0_Phi, hMarkerStyle, hMarkerSize, hMarkerColor, hLineColor);
        vecHistos.push_back(hClusterRMSConvPhotonPi0_Phi);


        histoName                           = "hClusterConvPhotonEta_Eta-RMS";
        if(i==0) vecHistosName.push_back(histoName);
        TH1D* hClusterRMSConvPhotonEta_Eta     = new TH1D(Form("%s_%s", histoName.Data(), vecDataSet.at(i).Data()),"hClusterConvPhotonEta_Eta; Run Number ; #sigma_{#eta_{#gamma_{conv} under #eta-peak}}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hClusterRMSConvPhotonEta_Eta, hMarkerStyle, hMarkerSize, hMarkerColor, hLineColor);
        vecHistos.push_back(hClusterRMSConvPhotonEta_Eta);


        histoName                           = "hClusterConvPhotonEta_Phi-RMS";
        if(i==0) vecHistosName.push_back(histoName);
        TH1D* hClusterRMSConvPhotonEta_Phi     = new TH1D(Form("%s_%s", histoName.Data(), vecDataSet.at(i).Data()),"hClusterConvPhotonEta_Phi; Run Number ; #sigma_{#phi_{#gamma_{conv} under #eta-peak}}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hClusterRMSConvPhotonEta_Phi, hMarkerStyle, hMarkerSize, hMarkerColor, hLineColor);
        vecHistos.push_back(hClusterRMSConvPhotonEta_Phi);


        histoName                   = "hClusterNLM-Mean";
        if(i==0) vecHistosName.push_back(histoName);
        TH1D* hClusterMeanNLM          = new TH1D(Form("%s_%s", histoName.Data(), vecDataSet.at(i).Data()),"hClusterMeanNLM; Run Number ; #bar{#lower[0.1]{NLM}}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hClusterMeanNLM, hMarkerStyle, hMarkerSize, hMarkerColor, hLineColor);
        vecHistos.push_back(hClusterMeanNLM);


        histoName                   = "hClusterNLM-RMS";
        if(i==0) vecHistosName.push_back(histoName);
        TH1D* hClusterRMSNLM           = new TH1D(Form("%s_%s", histoName.Data(), vecDataSet.at(i).Data()),"hClusterRMSNLM; Run Number ; #sigma_{NLM}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hClusterRMSNLM, hMarkerStyle, hMarkerSize, hMarkerColor, hLineColor);
        vecHistos.push_back(hClusterRMSNLM);

        fileCells               = Form("%s/%s/ClusterQA/%s/%s/Cells/%s.log",fCutSelection.Data(),fEnergyFlag.Data(),suffix.Data(),vecDataSet.at(i).Data(),vecDataSet.at(i).Data());
        if (runMergedClust)
        fileCells               = Form("%s/%s/MergedClusterQA/%s/%s/Cells/%s.log",fCutSelection.Data(),fEnergyFlag.Data(),suffix.Data(),vecDataSet.at(i).Data(),vecDataSet.at(i).Data());
        readin(fileCells, vecBadCells, kTRUE, kTRUE);
        fLog << "Read in " << vecBadCells.size() << " bad cells: ";
        for(Int_t iBad=0; iBad<(Int_t) vecBadCells.size(); iBad++)
        fLog << vecBadCells.at(iBad) << ", ";
        fLog << endl;

        if(!onlytrending){
            // *****************************************************************************************************
            // ****************************** Looping over Runs ****************************************************
            // *****************************************************************************************************
            printf("\n\t----------------------------------------------------------------------------\n\tLooping over Runs of DataSet |%s|\n\t----------------------------------------------------------------------------\n", (vecDataSetIn.at(i)).Data());
            fLog << "Looping over Runs:" << endl;




            if( vecNSplits.at(i) > 1) {
                printf("\t---------------------------------\n\tSplit: %i/%i | DataSet: %s | Runs: %i-%i \n\t---------------------------------\n",vecSubRunlists.at(i)+1, vecNSplits.at(i), (vecDataSet.at(i)).Data(), NMaxRuns*vecSubRunlists.at(i), (Int_t) (NMaxRuns*vecSubRunlists.at(i)+vecRuns.size()) );
            }
            Int_t NRunsProcessedInTotal=0;
            for(Int_t j=0; j<(Int_t) vecRuns.size(); j++) {
                NRunsProcessedInTotal++;
                fRunNumber                  = vecRuns.at(j);
                if ((fileNameData!=fileNameMC)&&(j==0)){
                    if (i>=nData){
                        if (i==nData){
                            cout<<"Switching File Name and nameMainDir to MC"<<endl;
                            cout<<"fileNameData: "<<fileNameData.Data()<<"; nameMainDirData: "<<nameMainDirData<<endl;
                            cout<<"fileNameMC: "<<fileNameMC.Data()<<"; nameMainDirMC: "<<nameMainDirMC<<endl;

                        }
                        fileName=fileNameMC;
                        nameMainDir=nameMainDirMC;
                    }
                    else{
                        fileName=fileNameData;
                        nameMainDir=nameMainDirData;
                    }
                }
                fRootFile                   = Form("%s/%s/%s/%s", filePath.Data(), vecDataSetIn.at(i).Data(), fRunNumber.Data(), fileName.Data());
                TFile* RootFile             = new TFile(fRootFile.Data(),"READ");
                if(RootFile->IsZombie()) {vecMissingRuns[i].push_back(fRunNumber); printf("\033[0;33mINFO\033[0m: ROOT file '%s' could not be openend, continue!\n",fRootFile.Data()); continue;}

                cout << endl;
                printf("\t\t----------------------------------------------------------------------------\n");
                cout << Form("%i/%i \t\tRun %s", NRunsProcessedInTotal, (Int_t) (vecRuns.size()), fRunNumber.Data()) << endl;
                printf("\t\tProcessing file: %s\n", fRootFile.Data() );
                printf("\t\t----------------------------------------------------------------------------\n");
                cout << endl;

                TList* TopDir               = (TList*) RootFile->Get(nameMainDir.Data());
                if(TopDir == NULL) {printf("\033[0;31mERROR\033[0m: TopDir not Found\n"); return;}
                else TopDir->SetOwner(kTRUE);
                TList* TopContainer         = (TList*) TopDir->FindObject(Form("Cut Number %s",fCutSelection.Data()));
                if(TopContainer == NULL) {printf("\033[0;31mERROR\033[0m: Cut Number %s not found in File",fCutSelection.Data()); return;}
                else TopContainer->SetOwner(kTRUE);
                TList* ESDContainer         = (TList*) TopContainer->FindObject(Form("%s ESD histograms",fCutSelection.Data()));
                if(ESDContainer == NULL) {printf("\033[0;31mERROR\033[0m: %s ESD histograms not found in File",fCutSelection.Data()); return;}
                else ESDContainer->SetOwner(kTRUE);
                TList* CaloCutsContainer    = (TList*) TopContainer->FindObject(Form("CaloCuts_%s",fClusterCutSelection.Data()));
                if(CaloCutsContainer == NULL) {printf("\033[0;31mERROR\033[0m: CaloCuts_%s not found in File",fClusterCutSelection.Data()); return;}
                else CaloCutsContainer->SetOwner(kTRUE);
                TList* ConvCutsContainer    = (TList*) TopContainer->FindObject(Form("ConvCuts_%s",fGammaCutSelection.Data()));
                if(isPCMCalo && ConvCutsContainer == NULL) {printf("\033[0;31mERROR\033[0m: ConvCuts_%s not found in File",fGammaCutSelection.Data()); return;}
                else if(ConvCutsContainer) ConvCutsContainer->SetOwner(kTRUE);
                TList* CaloExtQAContainer   = (TList*) TopContainer->FindObject(Form("CaloExtQA_%s",fClusterCutSelection.Data()));
                if(CaloExtQAContainer == NULL) {
                    printf("\033[0;35mWARNING\033[0m: CaloExtQA_%s not found in File, using CaloCuts-Container",fClusterCutSelection.Data());
                    CaloExtQAContainer  = CaloCutsContainer;
                } else CaloExtQAContainer->SetOwner(kTRUE);
                //--------------------------------------------------------------------------------------------------------
                //if(extendedRunwiseQA) ClusterQA(fRootFile,"","",nameMainDir,fEnergyFlag,vecplotDataSets.at(i),"","",fCutSelection,fRunNumber,suffix,mode);
                //--------------------------------------------------------------------------------------------------------
                if(doEquidistantXaxis)
                bin                     = mapBin[fRunNumber];
                else
                bin                     = fRunNumber.Atoi() - hFBin;
                Double_t nEvents            = 0;
                if ((TH1*) ESDContainer->FindObject("NEventsWOWeights")){
                    if (fEnergyFlag.Contains("PbPb") || fEnergyFlag.Contains("pPb"))
                    nEvents         = ((TH1*) ESDContainer->FindObject("NEventsWOWeights"))->GetBinContent(1);
                    else
                    nEvents         = GetNEvents((TH1*) ESDContainer->FindObject("NEventsWOWeights"),kFALSE);
                    printf("\033[0;33mINFO\033[0m: Output contains event weights\n");
                } else {
                    if (fEnergyFlag.Contains("PbPb") || fEnergyFlag.Contains("pPb"))
                    nEvents        = ((TH1*) ESDContainer->FindObject("NEvents"))->GetBinContent(1);
                    else
                    nEvents        = GetNEvents((TH1*) ESDContainer->FindObject("NEvents"),kFALSE);
                }
                //--------------------------------------------------------------------------------------------------------
                fLog << "----------------------------------------------------------------------------" << endl;
                fLog << "Processing file: " << fRootFile.Data() << endl;
                fLog << "Run: " << fRunNumber.Data() << ", with NEvents: " << nEvents << endl;
                fLog << "----------------------------------------------------------------------------" << endl;
                if(doExtQA==2) fLogRunwiseHotCells << "Run-" << fRunNumber.Data() << endl;
                if(doExtQA==2) fLogRunwiseColdCells << "Run-" << fRunNumber.Data() << endl;
                //--------------------------------------------------------------------------------------------------------
                if( nEvents < 1 ){printf("\033[0;35mWARNING\033[0m: number of accepted events in run: %f! Setting nEvents to 1...\n",nEvents); nEvents = 1;}
                //--------------------------------------------------------------------------------------------------------
                if(hNEvents){
                    hNEvents->SetBinContent(bin, nEvents);
                } else printf("\033[0;33mINFO\033[0m: Object |NEvents| could not be found! Skipping Fill...\n");

                //--------------------------------------------------------------------------------------------------------
                //------------------------------ Cluster properties: cluster energy --------------------------------------
                //--------------------------------------------------------------------------------------------------------
                TH1F* EOC                   = (TH1F*) CaloCutsContainer->FindObject(Form("EnergyOfCluster_afterClusterQA %s", fClusterCutSelection.Data()));
                if(EOC){
                    hClusterMeanEnergy->SetBinContent(bin, EOC->GetMean());
                    hClusterMeanEnergy->SetBinError(bin, EOC->GetMeanError());
                    hClusterRMSEnergy->SetBinContent(bin, EOC->GetRMS());
                    hClusterRMSEnergy->SetBinError(bin, EOC->GetRMSError());
                } else printf("\033[0;33mINFO\033[0m: Object |EnergyOfCluster_afterClusterQA| could not be found! Skipping Fill...\n");
                if(EOC){
                    Double_t eps = 0.0001;
                    Double_t Energy01 = EOC->Integral(1,EOC->GetXaxis()->FindBin(1-eps));
                    Double_t Energy14 = EOC->Integral(EOC->GetXaxis()->FindBin(1+eps),EOC->GetXaxis()->FindBin(4-eps));
                    Double_t Energy4 = EOC->Integral(EOC->GetXaxis()->FindBin(4+eps), EOC->GetXaxis()->GetLast());
                    hClusterEnergy01->SetBinContent(bin, Energy01 / nEvents);
                    hClusterEnergy01->SetBinError(bin, sqrt(Energy01) / nEvents);
                    hClusterEnergy14->SetBinContent(bin, Energy14 / nEvents);
                    hClusterEnergy14->SetBinError(bin, sqrt(Energy14) / nEvents);
                    hClusterEnergy4->SetBinContent(bin, Energy4 / nEvents);
                    hClusterEnergy4->SetBinError(bin, sqrt(Energy4) / nEvents);
                } else printf("\033[0;33mINFO\033[0m: Object |EnergyOfCluster_afterClusterQA| could not be found! Skipping Fill...\n");

                TH1D* ClusEnergy = (TH1D*) CaloCutsContainer->FindObject(Form("EnergyOfCluster_afterClusterQA %s", fClusterCutSelection.Data()));
                if(ClusEnergy){
                    TH1D* ClusEnergyBin = (TH1D*) ClusEnergy->Rebin(fNBinsClusterPt, "name", fBinsClusterPt);
                    TH1D* tempClusEnergy = new TH1D(*ClusEnergyBin);
                    delete ClusEnergyBin;
                    tempClusEnergy->Divide(fDeltaPt);
                    tempClusEnergy->Sumw2();
                    tempClusEnergy->Scale(1 / nEvents);
                    tempClusEnergy->GetXaxis()->SetTitle("Cluster Energy (GeV)");
                    tempClusEnergy->GetYaxis()->SetTitle("#frac{1}{N_{Events}} #frac{dN}{dE}");
                    vecClusterEnergy.push_back(tempClusEnergy);
                } else printf("\033[0;33mINFO\033[0m: Object |EnergyOfCluster_afterClusterQA| could not be found! Skipping...\n");

                //--------------------------------------------------------------------------------------------------------
                //----------------------------- Cluster properties: distribution on detector -----------------------------
                //--------------------------------------------------------------------------------------------------------
                TH2D* ClusEtaPhi = (TH2D*) CaloCutsContainer->FindObject(Form("EtaPhi_afterClusterQA %s", fClusterCutSelection.Data()));
                if(ClusEtaPhi){
                    TH2D* tempClusEtaPhi = new TH2D(*ClusEtaPhi);
                    tempClusEtaPhi->Scale(1 / nEvents);
                    tempClusEtaPhi->Scale(1 / GetMeanTH2(tempClusEtaPhi));
                    tempClusEtaPhi->SetName(Form("%s_EtaPhi_%s",fRunNumber.Data(),vecDataSet.at(i).Data()));
                    tempClusEtaPhi->SetTitle(fRunNumber);
                    tempClusEtaPhi->GetXaxis()->SetTitle("#phi");
                    tempClusEtaPhi->GetYaxis()->SetTitle("#eta");
                    tempClusEtaPhi->GetZaxis()->SetRangeUser(0, 2);
                    vecClusterEtaPhi.push_back(tempClusEtaPhi);
                } else printf("\033[0;33mINFO\033[0m: Object |EtaPhi_afterClusterQA| could not be found! Skipping...\n");


                if (i<nData){
                    TH2D* CellTimeVsCellID     = (TH2D*)CaloExtQAContainer->FindObject(Form("CellTimeVsCellID %s", fClusterCutSelection.Data()));
                    if (CellTimeVsCellID){
                        TH2D* tempCellTimeVsCellId = new TH2D(*CellTimeVsCellID);
                        tempCellTimeVsCellId->SetName(Form("%s_CellTimeVsId_%s",fRunNumber.Data(),vecDataSet.at(i).Data()));
                        tempCellTimeVsCellId->SetTitle(fRunNumber);
                        tempCellTimeVsCellId->GetXaxis()->SetTitle("Cell time (#mu s)");
                        tempCellTimeVsCellId->GetYaxis()->SetTitle("Cell ID");
                        vecCellTimeID.push_back(tempCellTimeVsCellId);

                    }
                }
                //--------------------------------------------------------------------------------------------------------
                //------------------------------ Cluster properties: cluster NCells --------------------------------------
                //--------------------------------------------------------------------------------------------------------
                TH1F* NCells = (TH1F*) CaloCutsContainer->FindObject(Form("NCellPerCluster_afterClusterQA %s", fClusterCutSelection.Data()));
                if(NCells){
                    hClusterMeanNCells->SetBinContent(bin, NCells->GetMean());
                    hClusterMeanNCells->SetBinError(bin, NCells->GetMeanError());
                    hClusterRMSNCells->SetBinContent(bin, NCells->GetRMS());
                    hClusterRMSNCells->SetBinError(bin, NCells->GetRMSError());
                } else printf("\033[0;33mINFO\033[0m: Object |NCellPerCluster_afterClusterQA| could not be found! Skipping Fill...\n");
                TH1D* ClusNCells = (TH1D*) CaloCutsContainer->FindObject(Form("NCellPerCluster_afterClusterQA %s", fClusterCutSelection.Data()));
                if(ClusNCells){
                    TH1D* tempClusNCells = new TH1D(*ClusNCells);
                    tempClusNCells->Sumw2();
                    tempClusNCells->Scale(1 / nEvents);
                    tempClusNCells->GetXaxis()->SetTitle("N_{Cells} per Cluster");
                    tempClusNCells->GetYaxis()->SetTitle("#frac{1}{N_{Events}} #frac{dN}{dN_{Cells}}");
                    vecClusterNCells.push_back(tempClusNCells);
                } else printf("\033[0;33mINFO\033[0m: Object |NCellPerCluster_afterClusterQA| could not be found! Skipping...\n");

                if(doExtQA>0){
                    TH2D* ClusEVsNCells = (TH2D*) CaloExtQAContainer->FindObject(Form("ClusterEnergyVsNCells_afterQA %s", fClusterCutSelection.Data()));
                    if(ClusEVsNCells){
                        TH2D* tempClusEVsNCells = new TH2D(*ClusEVsNCells);
                        tempClusEVsNCells->SetName(Form("%s_ClusterEnergyVsNCells_%s",fRunNumber.Data(),vecDataSet.at(i).Data()));
                        tempClusEVsNCells->SetTitle(fRunNumber);
                        tempClusEVsNCells->GetXaxis()->SetTitle("Cluster Energy (GeV)");
                        tempClusEVsNCells->GetYaxis()->SetTitle("#it{N}_{Cells} per Cluster");
                        vecClusterEVsNCells.push_back(tempClusEVsNCells);
                    } else printf("\033[0;33mINFO\033[0m: Object |ClusterEnergyVsNCells_afterQA| could not be found! Skipping...\n");
                }

                //--------------------------------------------------------------------------------------------------------
                //------------------------------ Cluster properties: cluster shape ---------------------------------------
                //--------------------------------------------------------------------------------------------------------
                TH1F* Dispersion = (TH1F*) CaloCutsContainer->FindObject(Form("Dispersion_afterClusterQA %s", fClusterCutSelection.Data()));
                if(Dispersion){
                    hClusterMeanDispersion->SetBinContent(bin, Dispersion->GetMean());
                    hClusterMeanDispersion->SetBinError(bin, Dispersion->GetMeanError());
                    hClusterRMSDispersion->SetBinContent(bin, Dispersion->GetRMS());
                    hClusterRMSDispersion->SetBinError(bin, Dispersion->GetRMSError());
                } else printf("\033[0;33mINFO\033[0m: Object |Dispersion_afterClusterQA| could not be found! Skipping Fill...\n");
                TH1F* M02 = (TH1F*) CaloCutsContainer->FindObject(Form("M02_afterClusterQA %s", fClusterCutSelection.Data()));
                if(M02){
                    hClusterMeanM02->SetBinContent(bin, M02->GetMean());
                    hClusterMeanM02->SetBinError(bin, M02->GetMeanError());
                    hClusterRMSM02->SetBinContent(bin, M02->GetRMS());
                    hClusterRMSM02->SetBinError(bin, M02->GetRMSError());
                } else printf("\033[0;33mINFO\033[0m: Object |M02_afterClusterQA| could not be found! Skipping Fill...\n");
                TH1F* M20 = (TH1F*) CaloCutsContainer->FindObject(Form("M20_afterClusterQA %s", fClusterCutSelection.Data()));
                if(M20){
                    hClusterMeanM20->SetBinContent(bin, M20->GetMean());
                    hClusterMeanM20->SetBinError(bin, M20->GetMeanError());
                    hClusterRMSM20->SetBinContent(bin, M20->GetRMS());
                    hClusterRMSM20->SetBinError(bin, M20->GetRMSError());
                } else printf("\033[0;33mINFO\033[0m: Object |M20_afterClusterQA| could not be found! Skipping Fill...\n");

                TH1D* ClusM02 = (TH1D*) CaloCutsContainer->FindObject(Form("M02_afterClusterQA %s", fClusterCutSelection.Data()));
                if(ClusM02){
                    TH1D* tempClusM02 = new TH1D(*ClusM02);
                    tempClusM02->Sumw2();
                    tempClusM02->Scale(1 / nEvents);
                    tempClusM02->GetXaxis()->SetTitle("#lambda_{0}^{2}");
                    tempClusM02->GetYaxis()->SetTitle("#frac{1}{N_{Events}} #frac{dN}{d#lambda_{0}^{2}}");
                    vecClusterM02.push_back(tempClusM02);
                } else printf("\033[0;33mINFO\033[0m: Object |M02_afterClusterQA| could not be found! Skipping...\n");
                TH1D* ClusM20 = (TH1D*) CaloCutsContainer->FindObject(Form("M20_afterClusterQA %s", fClusterCutSelection.Data()));
                if(ClusM20){
                    TH1D* tempClusM20 = new TH1D(*ClusM20);
                    tempClusM20->Sumw2();
                    tempClusM20->Scale(1 / nEvents);
                    tempClusM20->GetXaxis()->SetTitle("#lambda_{1}^{2}");
                    tempClusM20->GetYaxis()->SetTitle("#frac{1}{N_{Events}} #frac{dN}{d#lambda_{1}^{2}}");
                    vecClusterM20.push_back(tempClusM20);
                } else printf("\033[0;33mINFO\033[0m: Object |M20_afterClusterQA| could not be found! Skipping...\n");
                TH1D* ClusDispersion = (TH1D*) CaloCutsContainer->FindObject(Form("Dispersion_afterClusterQA %s", fClusterCutSelection.Data()));
                if(ClusDispersion){
                    TH1D* tempClusDispersion = new TH1D(*ClusDispersion);
                    tempClusDispersion->Sumw2();
                    tempClusDispersion->Scale(1 / nEvents);
                    tempClusDispersion->GetXaxis()->SetTitle("Dispersion of Cluster");
                    tempClusDispersion->GetYaxis()->SetTitle("#frac{1}{N_{Events}} #frac{dN}{dDisp}");
                    vecClusterDispersion.push_back(tempClusDispersion);
                } else printf("\033[0;33mINFO\033[0m: Object |Dispersion_afterClusterQA| could not be found! Skipping...\n");

                //--------------------------------------------------------------------------------------------------------
                //----------------------------- conversion photon bookkeeping --------------------------------------------
                //--------------------------------------------------------------------------------------------------------
                TH1F* CaloACC = (TH1F*) CaloCutsContainer->FindObject(Form("AcceptanceCuts %s", fClusterCutSelection.Data()));
                TH1F* CaloIPS = (TH1F*) CaloCutsContainer->FindObject(Form("IsPhotonSelected %s", fClusterCutSelection.Data()));
                if(CaloACC && CaloIPS){
                    Double_t CaloNClusters = CaloACC->GetBinContent(1);
                    Double_t CaloNClustersQA = CaloIPS->GetBinContent(5);
                    hCaloNClusters->SetBinContent(bin, CaloNClusters / nEvents);
                    hCaloNClusters->SetBinError(bin, sqrt(CaloNClusters) / nEvents);
                    hCaloNClustersQA->SetBinContent(bin, CaloNClustersQA / nEvents);
                    hCaloNClustersQA->SetBinError(bin, sqrt(CaloNClustersQA) / nEvents);
                } else printf("\033[0;33mINFO\033[0m: Object |AcceptanceCuts| or |IsPhotonSelected| could not be found! Skipping Fill...\n");

                //--------------------------------------------------------------------------------------------------------
                //------------------------------ Cluster properties: cluster timing --------------------------------------
                //--------------------------------------------------------------------------------------------------------
                TH2F* Time = (TH2F*) CaloCutsContainer->FindObject(Form("ClusterTimeVsE_afterClusterQA %s", fClusterCutSelection.Data()));
                if(Time){
                    hClusterMeanTime->SetBinContent(bin, Time->GetMean());
                    hClusterMeanTime->SetBinError(bin, Time->GetMeanError());
                    hClusterRMSTime->SetBinContent(bin, Time->GetRMS());
                    hClusterRMSTime->SetBinError(bin, Time->GetRMSError());
                } else printf("\033[0;33mINFO\033[0m: Object |ClusterTimeVsE_afterClusterQA| could not be found! Skipping Fill...\n");

                if(i<nData){
                    TH2D* ClusEnergyTime = (TH2D*) CaloCutsContainer->FindObject(Form("ClusterTimeVsE_afterClusterQA %s", fClusterCutSelection.Data()));
                    if(ClusEnergyTime){
                        TH2D* tempClusEnergyTime = new TH2D(*ClusEnergyTime);
                        //tempClusEnergyTime->Scale(1 / nEvents);
                        tempClusEnergyTime->SetName(Form("%s_EnergyTime_%s",fRunNumber.Data(),vecDataSet.at(i).Data()));
                        tempClusEnergyTime->SetTitle(fRunNumber);
                        Int_t min = 0, max = 0;
                        GetMinMaxBin(tempClusEnergyTime,min,max);
                        SetXRange(tempClusEnergyTime,min,max);
                        tempClusEnergyTime->GetXaxis()->SetTitle("#it{t}_{Cluster} (s)");
                        tempClusEnergyTime->GetYaxis()->SetTitle("#it{E}_{Cluster} (GeV)");
                        vecClusterEnergyTime.push_back(tempClusEnergyTime);
                    } else printf("\033[0;33mINFO\033[0m: Object |ClusterTimeVsE_afterClusterQA| could not be found! Skipping...\n");
                }
                TH2D* ClusTimeE = (TH2D*) CaloCutsContainer->FindObject(Form("ClusterTimeVsE_afterClusterQA %s", fClusterCutSelection.Data()));
                if(ClusTimeE){
                    for(Int_t iT=0; iT<5; iT++){
                        TH1D* ClusTime = (TH1D*) ClusTimeE->ProjectionX(Form("ProjectionClusTime_%f-%f",(minT_Energy[iT]*2.5)+1,(maxT_Energy[iT]*2.5)),(minT_Energy[iT]*2.5)+1,(maxT_Energy[iT]*2.5));
                        TH1D* tempClusTime = new TH1D(*ClusTime);
                        delete ClusTime;
                        tempClusTime->Sumw2();
                        tempClusTime->Scale(1 / nEvents);
                        tempClusTime->GetXaxis()->SetTitle("t_{Cluster} (ns)");
                        tempClusTime->GetYaxis()->SetTitle("#frac{1}{N_{Events}} #frac{dN}{dt}");
                        vecClusterTime[iT].push_back(tempClusTime);
                    }
                } else printf("\033[0;33mINFO\033[0m: Object |ClusterTimeVsE_afterClusterQA| could not be found! Skipping...\n");

                //--------------------------------------------------------------------------------------------------------
                //------------------------- Cluster properties: number of local maxima in cluster ------------------------
                //--------------------------------------------------------------------------------------------------------
                TH1F* NLM = (TH1F*) CaloCutsContainer->FindObject(Form("NLM_afterClusterQA %s", fClusterCutSelection.Data()));
                if(NLM){
                    hClusterMeanNLM->SetBinContent(bin, NLM->GetMean());
                    hClusterMeanNLM->SetBinError(bin, NLM->GetMeanError());
                    if(NLM->GetRMS()>0){
                        hClusterRMSNLM->SetBinContent(bin, NLM->GetRMS());
                        hClusterRMSNLM->SetBinError(bin, NLM->GetRMSError());
                    } else{
                        hClusterRMSNLM->SetBinContent(bin, 1);
                        hClusterRMSNLM->SetBinError(bin, 0);
                    }
                } else printf("\033[0;33mINFO\033[0m: Object |NLM_afterClusterQA| could not be found! Skipping Fill...\n");

                //--------------------------------------------------------------------------------------------------------
                //-------------------- Neutral mesons vs neutral mesons rejected due to V0 matching ----------------------
                //--------------------------------------------------------------------------------------------------------
                if(isPCMCalo){
                    TH2F* ESD_Mother = (TH2F*) ESDContainer->FindObject("ESD_Mother_InvMass_Pt");
                    TH2F* ESD_Mother_Matched = (TH2F*) ESDContainer->FindObject("ESD_MotherMatched_InvMass_Pt");
                    if(ESD_Mother && ESD_Mother_Matched){
                        CalculateFractionMatches(hClusterFractionMatches, ESD_Mother, ESD_Mother_Matched, bin, 0, 30);
                        CalculateFractionMatches(hClusterFractionMatchesS, ESD_Mother, ESD_Mother_Matched, bin, 0, 2.25);
                        CalculateFractionMatches(hClusterFractionMatchesM, ESD_Mother, ESD_Mother_Matched, bin, 2.25, 4);
                        CalculateFractionMatches(hClusterFractionMatchesH, ESD_Mother, ESD_Mother_Matched, bin, 4, 30);
                    } else printf("\033[0;33mINFO\033[0m: Object |ESD_Mother_InvMass_Pt| or |ESD_MotherMatched_InvMass_Pt| could not be found! Skipping Fill...\n");

                    for(Int_t iT=0; iT<2; iT++) {
                        TH2D* ConvPhotonEtaPhi = (TH2D*) ESDContainer->FindObject(Form("ESD_Mother%sConvPhoton_Eta_Phi", sT[iT].Data()));
                        if(ConvPhotonEtaPhi){
                            if(iT==0){
                                hClusterMeanConvPhotonPi0_Eta->SetBinContent(bin, ConvPhotonEtaPhi->GetMean(2));
                                hClusterMeanConvPhotonPi0_Eta->SetBinError(bin, ConvPhotonEtaPhi->GetMeanError(2));
                                hClusterMeanConvPhotonPi0_Phi->SetBinContent(bin, ConvPhotonEtaPhi->GetMean(1));
                                hClusterMeanConvPhotonPi0_Phi->SetBinError(bin, ConvPhotonEtaPhi->GetMeanError(1));
                                hClusterRMSConvPhotonPi0_Eta->SetBinContent(bin, ConvPhotonEtaPhi->GetRMS(2));
                                hClusterRMSConvPhotonPi0_Eta->SetBinError(bin, ConvPhotonEtaPhi->GetRMSError(2));
                                hClusterRMSConvPhotonPi0_Phi->SetBinContent(bin, ConvPhotonEtaPhi->GetRMS(1));
                                hClusterRMSConvPhotonPi0_Phi->SetBinError(bin, ConvPhotonEtaPhi->GetRMSError(1));
                            } else{
                                hClusterMeanConvPhotonEta_Eta->SetBinContent(bin, ConvPhotonEtaPhi->GetMean(2));
                                hClusterMeanConvPhotonEta_Eta->SetBinError(bin, ConvPhotonEtaPhi->GetMeanError(2));
                                hClusterMeanConvPhotonEta_Phi->SetBinContent(bin, ConvPhotonEtaPhi->GetMean(1));
                                hClusterMeanConvPhotonEta_Phi->SetBinError(bin, ConvPhotonEtaPhi->GetMeanError(1));
                                hClusterRMSConvPhotonEta_Eta->SetBinContent(bin, ConvPhotonEtaPhi->GetRMS(2));
                                hClusterRMSConvPhotonEta_Eta->SetBinError(bin, ConvPhotonEtaPhi->GetRMSError(2));
                                hClusterRMSConvPhotonEta_Phi->SetBinContent(bin, ConvPhotonEtaPhi->GetRMS(1));
                                hClusterRMSConvPhotonEta_Phi->SetBinError(bin, ConvPhotonEtaPhi->GetRMSError(1));
                            }
                            TH1D* ConvPhotonEta = (TH1D*) ConvPhotonEtaPhi->ProjectionY("ProjectionConvPhotonEta",1,200);
                            TH1D* ConvPhotonPhi = (TH1D*) ConvPhotonEtaPhi->ProjectionX("ProjectionConvPhotonPhi",1,600);
                            TH1D* tempConvPhotonEta = new TH1D(*ConvPhotonEta);
                            TH1D* tempConvPhotonPhi = new TH1D(*ConvPhotonPhi);
                            delete ConvPhotonEta;
                            delete ConvPhotonPhi;
                            //tempConvPhotonEta->GetXaxis()->SetRangeUser(-1,1);
                            tempConvPhotonEta->GetXaxis()->SetTitle(Form("#eta_{#gamma_{conv} under %s-peak}",sL[iT].Data()));
                            tempConvPhotonEta->GetYaxis()->SetTitle("#frac{1}{N_{Events}} #frac{dN}{d#eta}");
                            tempConvPhotonEta->SetTitle(Form("%s, nEvents: %.2e",fRunNumber.Data(),nEvents));
                            tempConvPhotonPhi->GetXaxis()->SetTitle(Form("#phi_{#gamma_{conv} under %s-peak}",sL[iT].Data()));
                            tempConvPhotonPhi->GetYaxis()->SetTitle("#frac{1}{N_{Events}} #frac{dN}{d#phi}");
                            tempConvPhotonPhi->SetTitle(Form("%s, nEvents: %.2e",fRunNumber.Data(),nEvents));
                            tempConvPhotonEta->Sumw2();
                            tempConvPhotonPhi->Sumw2();
                            tempConvPhotonEta->Scale(1 / nEvents);
                            tempConvPhotonPhi->Scale(1 / nEvents);
                            if(iT==0){
                                vecClusterPi0ConvPhotonEta.push_back(tempConvPhotonEta);
                                vecClusterPi0ConvPhotonPhi.push_back(tempConvPhotonPhi);
                            } else{
                                vecClusterEtaConvPhotonEta.push_back(tempConvPhotonEta);
                                vecClusterEtaConvPhotonPhi.push_back(tempConvPhotonPhi);
                            }
                        } else cout << Form("Info: Object |ESD_Mother%sConvPhoton_Eta_Phi| could not be found! Skipping...", sT[iT].Data()) << endl;
                    }
                }

                //--------------------------------------------------------------------------------------------------------
                //----------------------------- Cluster properties: track matching to cluster ----------------------------
                //--------------------------------------------------------------------------------------------------------
                if(isTrackMatching){
                    TH1F* R = (TH1F*) CaloCutsContainer->FindObject(Form("R_Cluster_afterClusterQA %s", fClusterCutSelection.Data()));
                    if(R){
                        hClusterMeanR->SetBinContent(bin, R->GetMean());
                        hClusterMeanR->SetBinError(bin, R->GetMeanError());
                        hClusterRMSR->SetBinContent(bin, R->GetRMS());
                        hClusterRMSR->SetBinError(bin, R->GetRMSError());
                    } else printf("\033[0;33mINFO\033[0m: Object |R_Cluster_afterClusterQA| could not be found! Skipping Fill...\n");

                    TH2F* DeltaEtaPhi = (TH2F*) CaloCutsContainer->FindObject(Form("dEtaVsdPhi_beforeClusterQA %s", fClusterCutSelection.Data()));
                    if(DeltaEtaPhi){
                        hClusterMeanDeltaEta->SetBinContent(bin, DeltaEtaPhi->GetMean(1));
                        hClusterMeanDeltaEta->SetBinError(bin, DeltaEtaPhi->GetMeanError(1));
                        hClusterMeanDeltaPhi->SetBinContent(bin, DeltaEtaPhi->GetMean(2));
                        hClusterMeanDeltaPhi->SetBinError(bin, DeltaEtaPhi->GetMeanError(2));
                        hClusterRMSDeltaEta->SetBinContent(bin, DeltaEtaPhi->GetRMS(1));
                        hClusterRMSDeltaEta->SetBinError(bin, DeltaEtaPhi->GetRMSError(1));
                        hClusterRMSDeltaPhi->SetBinContent(bin, DeltaEtaPhi->GetRMS(2));
                        hClusterRMSDeltaPhi->SetBinError(bin, DeltaEtaPhi->GetRMSError(2));
                    } else printf("\033[0;33mINFO\033[0m: Object |dEtaVsdPhi_beforeClusterQA| could not be found! Skipping Fill...\n");

                    TH2D* ClusdEta_dPhi_before = (TH2D*)CaloCutsContainer->FindObject(Form("dEtaVsdPhi_beforeClusterQA %s", fClusterCutSelection.Data()));
                    TH2D* ClusdEta_dPhi_after = (TH2D*)CaloCutsContainer->FindObject(Form("dEtaVsdPhi_afterClusterQA %s", fClusterCutSelection.Data()));
                    if(ClusdEta_dPhi_before && ClusdEta_dPhi_after){
                        TH2D* ClusdEta_dPhi_matched = (TH2D*) ClusdEta_dPhi_before->Clone();
                        ClusdEta_dPhi_matched->Add(ClusdEta_dPhi_after,-1);
                        TH1D* ClusdEta_matched = (TH1D*) ClusdEta_dPhi_matched->ProjectionX("ProjectionClusdEta",1,240);
                        TH1D* ClusdPhi_matched = (TH1D*) ClusdEta_dPhi_matched->ProjectionY("ProjectionClusdPhi",1,240);
                        TH1D* tempClusdEta_matched = new TH1D(*ClusdEta_matched);
                        TH1D* tempClusdPhi_matched = new TH1D(*ClusdPhi_matched);
                        delete ClusdEta_matched;
                        delete ClusdPhi_matched;
                        tempClusdEta_matched->GetXaxis()->SetTitle("#Delta#eta");
                        tempClusdEta_matched->GetYaxis()->SetTitle("#frac{1}{N} #frac{dN}{d#Delta#eta}");
                        tempClusdPhi_matched->GetXaxis()->SetTitle("#Delta#phi");
                        tempClusdPhi_matched->GetYaxis()->SetTitle("#frac{1}{N} #frac{dN}{d#Delta#phi}");
                        tempClusdEta_matched->Sumw2();
                        tempClusdPhi_matched->Sumw2();
                        if(tempClusdEta_matched->Integral() > 0) tempClusdEta_matched->Scale(1 / ((Double_t)tempClusdEta_matched->Integral()));
                        if(tempClusdPhi_matched->Integral() > 0) tempClusdPhi_matched->Scale(1 / ((Double_t)tempClusdPhi_matched->Integral()));

                        vecClusterDeltaEta.push_back(tempClusdEta_matched);
                        vecClusterDeltaPhi.push_back(tempClusdPhi_matched);
                    } else printf("\033[0;33mINFO\033[0m: Object |dEtaVsdPhi_beforeClusterQA| or |dEtaVsdPhi_afterClusterQA| could not be found! Skipping...\n");
                }

                //--------------------------------------------------------------------------------------------------------
                //------------------------------ Super module properties -------------------------------------------------
                //--------------------------------------------------------------------------------------------------------
                if(doExtQA>0){
                    TH2D* ClusterEVsModule = (TH2D*) CaloExtQAContainer->FindObject(Form("ClusterEnergyVsModule_afterClusterQA %s", fClusterCutSelection.Data()));
                    TH2D* ModuleEnergyEVsModule = (TH2D*) CaloExtQAContainer->FindObject(Form("ModuleEnergyVsModule %s", fClusterCutSelection.Data()));
                    TH2D* NCells100VsModule = (TH2D*) CaloExtQAContainer->FindObject(Form("NCellsAbove100VsModule %s", fClusterCutSelection.Data()));
                    TH2D* NCells1500VsModule = (TH2D*) CaloExtQAContainer->FindObject(Form("NCellsAbove1500VsModule %s", fClusterCutSelection.Data()));

                    for(Int_t iModule=0; iModule<nCaloModules; iModule++) {
                        if(ClusterEVsModule){
                            TH1D* ClusterEForModule = (TH1D*) ClusterEVsModule->ProjectionX(Form("projectClusterE_%i",iModule),iModule+1,iModule+1);
                            TH1D* tempClusterEForModule = new TH1D(*ClusterEForModule);
                            delete ClusterEForModule;
                            tempClusterEForModule->GetXaxis()->SetTitle(Form("Cluster Energy of SM%i (GeV)",iModule));
                            tempClusterEForModule->GetYaxis()->SetTitle("#frac{dE}{dN}");
                            tempClusterEForModule->SetTitle(Form("%s, nEvents: %.2e",fRunNumber.Data(),nEvents));
                            tempClusterEForModule->Sumw2();
                            tempClusterEForModule->Scale(1 / nEvents);
                            vecClusterEVsModule[iModule].push_back(tempClusterEForModule);
                        } else printf("\033[0;33mINFO\033[0m: Object |ClusterEnergyVsModule_afterClusterQA| could not be found! Skipping...\n");
                        if(ModuleEnergyEVsModule){
                            TH1D* ModuleEForModule = (TH1D*) ModuleEnergyEVsModule->ProjectionX(Form("projectModuleE_%i",iModule),iModule+1,iModule+1);
                            TH1D* tempModuleEForModule = new TH1D(*ModuleEForModule);
                            delete ModuleEForModule;
                            tempModuleEForModule->GetXaxis()->SetTitle(Form("Total SuperModule Energy per Event for SM%i(GeV)",iModule));
                            tempModuleEForModule->GetYaxis()->SetTitle("#frac{dE}{dN}");
                            tempModuleEForModule->SetTitle(Form("%s, nEvents: %.2e",fRunNumber.Data(),nEvents));
                            tempModuleEForModule->Sumw2();
                            tempModuleEForModule->Scale(1 / nEvents);
                            vecModuleEVsModule[iModule].push_back(tempModuleEForModule);
                        } else printf("\033[0;33mINFO\033[0m: Object |ModuleEnergyVsModule_afterClusterQA| could not be found! Skipping...\n");
                        if(NCells100VsModule){
                            TH1D* NCells100ForModule = (TH1D*) NCells100VsModule->ProjectionX(Form("projectNCells100_%i",iModule),iModule+1,iModule+1);
                            TH1D* tempNCells100ForModule = new TH1D(*NCells100ForModule);
                            delete NCells100ForModule;
                            tempNCells100ForModule->GetXaxis()->SetTitle(Form("#it{N}_{Cells}>100 MeV in SM%i per Event",iModule));
                            tempNCells100ForModule->GetYaxis()->SetTitle("#frac{d#it{N}_{Cells}>100 MeV}{dN}");
                            tempNCells100ForModule->SetTitle(Form("%s, nEvents: %.2e",fRunNumber.Data(),nEvents));
                            tempNCells100ForModule->Sumw2();
                            tempNCells100ForModule->Scale(1 / nEvents);
                            vecClusterNCells100VsModule[iModule].push_back(tempNCells100ForModule);
                        } else printf("\033[0;33mINFO\033[0m: Object |NCellsAbove100VsModule| could not be found! Skipping...\n");
                        if(NCells1500VsModule){
                            TH1D* NCells1500ForModule = (TH1D*) NCells1500VsModule->ProjectionX(Form("projectNCells1500_%i",iModule),iModule+1,iModule+1);
                            TH1D* tempNCells1500ForModule = new TH1D(*NCells1500ForModule);
                            delete NCells1500ForModule;
                            tempNCells1500ForModule->GetXaxis()->SetTitle(Form("#it{N}_{Cells}>1500 MeV in SM%i per Event",iModule));
                            tempNCells1500ForModule->GetYaxis()->SetTitle("#frac{d#it{N}_{Cells}>1500 MeV}{dN}");
                            tempNCells1500ForModule->SetTitle(Form("%s, nEvents: %.2e",fRunNumber.Data(),nEvents));
                            tempNCells1500ForModule->Sumw2();
                            tempNCells1500ForModule->Scale(1 / nEvents);
                            vecClusterNCells1500VsModule[iModule].push_back(tempNCells1500ForModule);
                        } else printf("\033[0;33mINFO\033[0m: Object |NCellsAbove1500VsModule| could not be found! Skipping...\n");
                    }

                    //--------------------------------------------------------------------------------------------------------
                    //------------------------- Cell QA: cell frequency, timing and energy distribution ----------------------
                    //--------------------------------------------------------------------------------------------------------
                    TH1D* ClusIncludedCells = (TH1D*)CaloExtQAContainer->FindObject(Form("ClusterIncludedCells_afterClusterQA %s", fClusterCutSelection.Data()));
                    if(ClusIncludedCells){
                        TH1D* tempClusIncludedCells = new TH1D(*ClusIncludedCells);
                        tempClusIncludedCells->SetName(Form("%s_ClusterIncludedCells_%s",fRunNumber.Data(),vecDataSet.at(i).Data()));
                        tempClusIncludedCells->SetTitle(fRunNumber);
                        tempClusIncludedCells->GetXaxis()->SetTitle("Cell ID in accepted Clusters");
                        tempClusIncludedCells->GetYaxis()->SetTitle("d#it{N}/d#it{CellID}");
                        tempClusIncludedCells->Sumw2();
                        vecClusterIncludedCells.push_back(tempClusIncludedCells);
                    } else printf("\033[0;33mINFO\033[0m: Object |ClusterIncludedCells_afterClusterQA| could not be found! Skipping...\n");
                    //--------------------------------------------------------------------------------------------------------
                    if(doExtQA==2){
                        TH1D* DeadCellsClusIncludedCellsBefore = (TH1D*)CaloExtQAContainer->FindObject(Form("ClusterIncludedCells_beforeClusterQA %s", fClusterCutSelection.Data()));
                        if(DeadCellsClusIncludedCellsBefore){
                            TH2D* tempDeadCellsRunwise = CompareDeadCellsRunwise(DeadCellsClusIncludedCellsBefore, nCaloCells, Form("%s, nEvents: %.2e",fRunNumber.Data(),nEvents));
                            tempDeadCellsRunwise->SetName(Form("%s_ClusterIncludedCells_%s",fRunNumber.Data(),vecDataSet.at(i).Data()));
                            tempDeadCellsRunwise->GetYaxis()->SetTitle("Cell ID");
                            vecClusterFiredCellIDs.push_back(tempDeadCellsRunwise);
                        } else printf("\033[0;33mINFO\033[0m: Object |ClusterIncludedCells_beforeClusterQA| could not be found! Skipping...\n");
                    }
                    //--------------------------------------------------------------------------------------------------------
                    TH1D* ClusIncludedCellsBefore = (TH1D*)CaloExtQAContainer->FindObject(Form("ClusterIncludedCells_beforeClusterQA %s", fClusterCutSelection.Data()));
                    if(ClusIncludedCellsBefore){
                        TH1D* tempClusIncludedCellBefore = new TH1D(*ClusIncludedCellsBefore);
                        //					if(doExtQA==2) CheckHotCellsRunwise(fLogRunwiseBadCells[i],tempClusIncludedCellBefore,nCaloCells,kFALSE);
                        //					if(doExtQA==2) CheckDeadCellsRunwise(fLogRunwiseDeadCells[i],tempClusIncludedCellBefore,nCaloCells);
                        tempClusIncludedCellBefore->SetName(Form("%s_ClusterIncludedCells_beforeClusterQA_%s",fRunNumber.Data(),vecDataSet.at(i).Data()));
                        tempClusIncludedCellBefore->SetTitle(fRunNumber);
                        tempClusIncludedCellBefore->GetXaxis()->SetTitle("Cell ID in all Clusters");
                        tempClusIncludedCellBefore->GetYaxis()->SetTitle("d#it{N}/d#it{CellID}");
                        tempClusIncludedCellBefore->Sumw2();
                        vecClusterIncludedCellsBefore.push_back(tempClusIncludedCellBefore);
                    } else printf("\033[0;33mINFO\033[0m: Object |ClusterIncludedCells_beforeClusterQA| could not be found! Skipping...\n");
                    //--------------------------------------------------------------------------------------------------------
                    TProfile* ClusBadCells = (TProfile*)CaloExtQAContainer->FindObject(Form("%s - Bad Channels",calo.Data()));
                    if(ClusBadCells){
                        TProfile* tempClusBadCells = new TProfile(*ClusBadCells);
                        tempClusBadCells->SetName(Form("%s_BadCells_%s",fRunNumber.Data(),vecDataSet.at(i).Data()));
                        tempClusBadCells->SetTitle(fRunNumber);
                        tempClusBadCells->GetXaxis()->SetTitle("Cell ID");
                        tempClusBadCells->GetYaxis()->SetTitle("Cell ID Bad in Fraction of Events");
                        vecClusterBadCells.push_back(tempClusBadCells);
                    } else cout << Form("Info: Object |%s - Bad Channels| could not be found! Skipping...",calo.Data()) << endl;
                    //--------------------------------------------------------------------------------------------------------
                    TH1D* ClusEFracCells = (TH1D*)CaloExtQAContainer->FindObject(Form("ClusterEnergyFracCells_afterClusterQA %s", fClusterCutSelection.Data()));
                    if(ClusEFracCells){
                        TH1D* tempClusEFracCells = new TH1D(*ClusEFracCells);
                        tempClusEFracCells->SetName(Form("%s_EFracCells_%s",fRunNumber.Data(),vecDataSet.at(i).Data()));
                        tempClusEFracCells->SetTitle(fRunNumber);
                        tempClusEFracCells->GetXaxis()->SetTitle("Cell ID in accepted Clusters");
                        tempClusEFracCells->GetYaxis()->SetTitle("#sum^{events} E-Frac of Cell");
                        tempClusEFracCells->Sumw2();
                        vecClusterEFracCells.push_back(tempClusEFracCells);
                    } else printf("\033[0;33mINFO\033[0m: Object |ClusterEnergyFracCells_afterClusterQA| could not be found! Skipping...\n");
                    //--------------------------------------------------------------------------------------------------------
                    TH1D* ClusEFracCellsBefore = (TH1D*)CaloExtQAContainer->FindObject(Form("ClusterEnergyFracCells_beforeClusterQA %s", fClusterCutSelection.Data()));
                    if(ClusEFracCellsBefore){
                        TH1D* tempClusEFracCellBefore = new TH1D(*ClusEFracCellsBefore);
                        if(doExtQA==2) CheckHotAndColdCellsEFracRunwise(fLogRunwiseHotCells,fLogRunwiseColdCells,tempClusEFracCellBefore,ClusBadCells,nCaloCells,kFALSE);
                        tempClusEFracCellBefore->SetName(Form("%s_EFracCells_beforeQA_%s",fRunNumber.Data(),vecDataSet.at(i).Data()));
                        tempClusEFracCellBefore->SetTitle(fRunNumber);
                        tempClusEFracCellBefore->GetXaxis()->SetTitle("Cell ID in all Clusters");
                        tempClusEFracCellBefore->GetYaxis()->SetTitle("#sum^{events} E-Frac of Cell");
                        tempClusEFracCellBefore->Sumw2();
                        vecClusterEFracCellsBefore.push_back(tempClusEFracCellBefore);
                    } else printf("\033[0;33mINFO\033[0m: Object |ClusterEnergyFracCells_beforeClusterQA| could not be found! Skipping...\n");
                    //--------------------------------------------------------------------------------------------------------
                    if(doExtQA==2){
                        TH2D* fHistCellEnergyVsCellID = (TH2D*)CaloExtQAContainer->FindObject(Form("CellEnergyVsCellID %s", fClusterCutSelection.Data()));
                        TH2D* fHistCellTimeVsCellID = (TH2D*)CaloExtQAContainer->FindObject(Form("CellTimeVsCellID %s", fClusterCutSelection.Data()));
                        if(fHistCellEnergyVsCellID && fHistCellTimeVsCellID){
                            TH2D** tempClusCell = PlotCellMeanVsSigmaForRunwise(nCaloCells,fHistCellEnergyVsCellID,fHistCellTimeVsCellID, "Mean Cell Energy (GeV)","#sigma_{Cell Energy} (GeV)","Mean Cell Time (s)","#sigma_{Cell Time} (s)",(i>=nData));
                            tempClusCell[0]->SetName(Form("%s_ClusterEMeanVsSigma_%s",fRunNumber.Data(),vecDataSet.at(i).Data()));
                            tempClusCell[0]->SetTitle(fRunNumber);
                            vecClusterEnergyMeanSigma.push_back(tempClusCell[0]);
                            tempClusCell[1]->SetName(Form("%s_ClusterTimeMeanVsSigma_%s",fRunNumber.Data(),vecDataSet.at(i).Data()));
                            tempClusCell[1]->SetTitle(fRunNumber);
                            vecClusterTimeMeanSigma.push_back(tempClusCell[1]);
                        } else printf("\033[0;33mINFO\033[0m: Object |CellEnergyVsCellID or CellTimeVsCellID| could not be found! Skipping...\n");

                        for(Int_t j=0; j<(Int_t) vecBadCells.size(); j++) {
                            TH2D* fHistCellEnergyVsCellID = (TH2D*)CaloExtQAContainer->FindObject(Form("CellEnergyVsCellID %s", fClusterCutSelection.Data()));
                            TH2D* fHistCellTimeVsCellID = (TH2D*)CaloExtQAContainer->FindObject(Form("CellTimeVsCellID %s", fClusterCutSelection.Data()));
                            if(fHistCellEnergyVsCellID && fHistCellTimeVsCellID){
                                TH1D* tempEnergy;
                                TH1D* tempTime;
                                tempEnergy = (TH1D*) fHistCellEnergyVsCellID->ProjectionX("energy",((TString)vecBadCells.at(j)).Atoi(),((TString)vecBadCells.at(j)).Atoi());
                                tempTime = (TH1D*) fHistCellTimeVsCellID->ProjectionX("time",((TString)vecBadCells.at(j)).Atoi(),((TString)vecBadCells.at(j)).Atoi());

                                TH1D* tempEnergyCell = new TH1D(*tempEnergy);
                                delete tempEnergy;
                                tempEnergyCell->GetXaxis()->SetTitle(Form("Cell Energy of ID %i (GeV)",((TString)vecBadCells.at(j)).Atoi()));
                                tempEnergyCell->GetYaxis()->SetTitle("#frac{dE}{dN}");
                                tempEnergyCell->Sumw2();
                                tempEnergyCell->Scale(1 / nEvents);
                                vecBadCellsEnergy.push_back(tempEnergyCell);

                                TH1D* tempTimeCell = new TH1D(*tempTime);
                                delete tempTime;
                                tempTimeCell->GetXaxis()->SetTitle(Form("Cell Time of ID %i (s)",((TString)vecBadCells.at(j)).Atoi()));
                                tempTimeCell->GetYaxis()->SetTitle("#frac{dTime}{dN}");
                                tempTimeCell->Sumw2();
                                tempTimeCell->Scale(1 / nEvents);
                                vecBadCellsTime.push_back(tempTimeCell);
                            }
                        }
                    }
                }

                //--------------------------------------------------------------------------------------------------------
                //-------------------------------- clean up --------------------------------------------------------------
                //--------------------------------------------------------------------------------------------------------
                delete TopDir;

                RootFile->Close();
                delete RootFile;
            }

            // ****************************** Drawing Histograms ************************************************
            printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
            printf("Drawing Histograms\n");



            if(doHistsForEverySet) {

                printf("DataSet: %s\n", vecDataSet.at(i).Data());
                outputDirDataSet        = Form("%s/%s",outputDir.Data(), vecDataSet.at(i).Data());

                if(useDataRunListForMC && i>=nData) {
                    outputDirDataSet    = Form("%s/%s-%s", outputDir.Data(), vecDataSet.at(i).Data(),vecDataSet.at(0).Data());
                    printf("Switch useDataRunListForMC is true, output to: %s\n", outputDirDataSet.Data());
                }
                gSystem->Exec("mkdir -p "+outputDirDataSet);
                gSystem->Exec("mkdir -p "+outputDirDataSet+"/ExtQA");
                gSystem->Exec("mkdir -p "+outputDirDataSet+"/ExtQA/EFrac");
                if(doExtQA>0) gSystem->Exec("mkdir -p "+outputDirDataSet+"/ExtQA/IncludedCells");
                if(doExtQA==2) gSystem->Exec("mkdir -p "+outputDirDataSet+"/ExtQA/MissingCells");
                gSystem->Exec("mkdir -p "+outputDirDataSet+"/ExtQA/ConvPhotonEtaPhi");
                if(doExtQA==2) gSystem->Exec("mkdir -p "+outputDirDataSet+"/BadCells");

                TString fTrigger        = "";
                if(i<nData){
                    TString fTriggerCut = fEventCutSelection(3,2);
                    fTrigger            = AnalyseSpecialTriggerCut(fTriggerCut.Atoi(), vecplotDataSets.at(i));
                    printf("Trigger: '%s'\n",fTrigger.Data());
                    if(fTrigger.Contains("not defined")){
                        fTrigger        = "";
                        printf("\033[0;33mINFO\033[0m: Trigger cut not defined!\n");
                    }
                }

                DrawVectorOverviewTH2D( canvaslarge, vecClusterEtaPhi, "hClusterEtaPhi_scaledNEventsAndMean", outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, boxLabel, kFALSE, kFALSE);
                if(i<nData)
                DrawVectorOverviewTH2D( canvaslarge, vecCellTimeID, "hCellTimeVsId", outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, boxLabel, kFALSE, kTRUE);

                TGaxis::SetExponentOffset(0, -0.1, "x");
                if(i<nData)
                DrawVectorOverviewTH2D( canvaslarge, vecClusterEnergyTime, "ExtQA/hClusterEnergyTime", outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.14, 0.8, 0.8, 0.12, 0.93, 0x0, kTRUE, kTRUE);
                TGaxis::SetExponentOffset(0, 0, "x");

                if(doExtQA>0){
                    DrawVectorOverviewTH2D( canvaslarge, vecClusterEVsNCells, "ExtQA/hClusterEVsNCells", outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.14, 0.8, 0.8, 0.12, 0.93, boxLabel2, kFALSE, kTRUE);
                    DrawVectorOverviewTH1D( canvaslarge, vecClusterIncludedCells, "ExtQA/IncludedCells/hClusterIncludedCells", outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.14, 0.8, 0.9, 0.12, 0.93, 0x0, kFALSE, kTRUE);
                    DrawVectorOverviewTH1D(canvaslarge, vecClusterIncludedCellsBefore, "ExtQA/IncludedCells/hClusterIncludedCellsBefore", outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.14, 0.8, 0.9, 0.12, 0.93, 0x0, kFALSE, kTRUE);
                }


                DrawVectorOverviewTH1D( canvaslarge, vecClusterEFracCellsBefore, "ExtQA/EFrac/hClusterEFracCellsBefore", outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.14, 0.8, 0.9, 0.12, 0.93, 0x0, kFALSE, kTRUE);
                DrawVectorOverviewTH1D( canvaslarge, vecClusterEFracCells, "ExtQA/EFrac/hClusterEFracCells", outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.14, 0.8, 0.9, 0.12, 0.93, 0x0, kFALSE, kTRUE);

                if(isPCMCalo){
                    DrawVectorOverviewTH1D( canvas, vecClusterPi0ConvPhotonEta, "ExtQA/ConvPhotonEtaPhi/hClusterPi0ConvPhotonEta_Runwise", outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.14, 0.8, 0.9, 0.2, 0.93, 0x0, kFALSE, kTRUE);
                    DrawVectorOverviewTH1D( canvas, vecClusterPi0ConvPhotonPhi, "ExtQA/ConvPhotonEtaPhi/hClusterPi0ConvPhotonPhi_Runwise", outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.14, 0.8, 0.9, 0.2, 0.93, 0x0, kFALSE, kTRUE);
                    DrawVectorOverviewTH1D( canvas, vecClusterEtaConvPhotonEta, "ExtQA/ConvPhotonEtaPhi/hClusterEtaConvPhotonEta_Runwise", outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.14, 0.8, 0.9, 0.2, 0.93, 0x0, kFALSE, kTRUE);
                    DrawVectorOverviewTH1D( canvas, vecClusterEtaConvPhotonPhi, "ExtQA/ConvPhotonEtaPhi/hClusterEtaConvPhotonPhi_Runwise", outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.14, 0.8, 0.9, 0.2, 0.93, 0x0, kFALSE, kTRUE);
                }

                //--------------------------------------------------------------------------------------------------------
                for(Int_t h=0; h<(Int_t) vecHistos.size(); h++) {
                    if(!vecHistos.at(h)) {
                        printf("\033[0;31mERROR\033[0m*********   Vector to Histo %s is empty\n",vecHistosName.at(h).Data());
                        continue;
                    }
                    vecHistos.at(h)->SetTitle("");
                    if( ((TString)vecHistosName.at(h)).CompareTo("hNEvents")==0 ) {
                        AdjustHistRange(((TH1D*) vecHistos.at(h)),10.,10.,kTRUE, 0, 0.);
                        ((TH1D*) vecHistos.at(h))->Draw("p");
                    } else{
                        AdjustHistRange(((TH1D*) vecHistos.at(h)),1.1,1.1,kTRUE, 0, 0.);
                        ((TH1D*) vecHistos.at(h))->Draw("px0e1");
                    }

                    if(doTrigger && i<nData){
                        PutProcessLabelAndEnergyOnPlot(xPosLabel, 0.92, 0.03, fCollisionSystem.Data(), vecplotDataSets.at(i).Data(), fTrigger.Data());
                        PutProcessLabelAndEnergyOnPlot(xPosLabel, 0.81, 0.03, Form("%s clusters", calo.Data()), "", "");
                    } else{
                        PutProcessLabelAndEnergyOnPlot(xPosLabel, 0.92, 0.03, fCollisionSystem.Data(), vecplotDataSets.at(i).Data(), Form("%s clusters", calo.Data()));
                    }

                    if( ((TString)vecHistosName.at(h)).CompareTo("hNEvents")==0 )
                    SaveCanvas(canvas, Form("%s/%s.%s", outputDirDataSet.Data(), vecHistosName.at(h).Data(), suffix.Data()), kFALSE, kTRUE);
                    else
                    SaveCanvas(canvas, Form("%s/%s.%s", outputDirDataSet.Data(), vecHistosName.at(h).Data(), suffix.Data()));
                }




                // **********************************************************************************************************
                // ****************************** Drawing Runwise Histograms ************************************************
                // **********************************************************************************************************
                printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
                printf("Drawing Runwise Histograms\n");

                vecRuns.clear();
                // if(!readin(fileRuns.at(i), vecRuns, kFALSE)) {printf("\033[0;31mERROR\033[0m, no Run Numbers could be found! Returning...\n"); return;}
                if(!readin(fileRuns.at(i), vecRuns, kTRUE, kFALSE, NMaxRuns*vecSubRunlists.at(i), NMaxRuns*(vecSubRunlists.at(i)+1))){printf("\033[0;31mERROR\033[0m, no Run Numbers could be found! Returning...\n"); return;}

                for(Int_t iRun=0; iRun<(Int_t)vecMissingRuns[i].size(); iRun++){
                    vecRuns.erase(std::remove(vecRuns.begin(), vecRuns.end(), vecMissingRuns[i].at(iRun)), vecRuns.end());
                }

                Int_t NColumns              = ((Int_t) vecRuns.size() / 31 ) + 1;

                //extending nHist for testing
                //for(Int_t i=123567; i<123667; i++) {globalRuns.push_back(Form("%i",i)); vecClusterEnergy[0].push_back(vecClusterEnergy[0].at(0));}

                canvasRunwise               = new TCanvas("canvasRunwise","",10,10,1350+(NColumns*108),900);  // gives the page size
                DrawGammaCanvasSettings(canvasRunwise, 130.5/(1350.+(NColumns*108.)), (40.5+(NColumns*108.))/(1350.+(NColumns*108.)), topMargin, bottomMargin);
                canvasRunwise->cd();

                Double_t addRight           = ((Double_t)NColumns*108)/(1350+((Double_t)NColumns*108));
                legendRuns                  = new TLegend(0.98-addRight,0.09,0.98,0.94);
                legendRuns->SetNColumns(NColumns);
                legendRuns->SetFillColor(0);
                legendRuns->SetLineColor(0);
                legendRuns->SetTextSize(0.03);
                legendRuns->SetTextFont(42);

                printf("DataSet: %s\n", vecplotDataSets.at(i).Data());
                // TString fTrigger            = "";
                if(doTrigger && i<nData){
                    TString fTriggerCut     = fEventCutSelection(3,2);
                    fTrigger                = AnalyseSpecialTriggerCut(fTriggerCut.Atoi(), vecplotDataSets.at(i));
                    if(fTrigger.Contains("not defined"))
                    fTrigger            = "";
                }
                //--------------------------------------------------------------------------------------------------------
                DrawVectorRunwiseTH1D(	canvasRunwise, legendRuns, vecClusterEnergy, vecRuns, 5, 5, kFALSE, addRight, 0.8, 0.92, 0.03, 0.8, 0.81, 0.03, doTrigger, fTrigger, (Bool_t)(i<nData), outputDirDataSet, "ClusterEnergy_Runwise", vecplotDataSets.at(i),kFALSE, fCollisionSystem, Form("%s clusters", calo.Data()), suffix, kTRUE, kTRUE, kFALSE);
                //--------------------------------------------------------------------------------------------------------
                DrawVectorRunwiseTH1D(	canvasRunwise, legendRuns, vecClusterM02, vecRuns, 5, 5, kTRUE, addRight, 0.8, 0.92, 0.03, 0.8, 0.81, 0.03, doTrigger, fTrigger, (Bool_t)(i<nData), outputDirDataSet, "ClusterM02_Runwise", vecplotDataSets.at(i), kFALSE, fCollisionSystem, Form("%s clusters", calo.Data()), suffix, kFALSE, kTRUE, kFALSE);
                //--------------------------------------------------------------------------------------------------------
                DrawVectorRunwiseTH1D(	canvasRunwise, legendRuns, vecClusterM20, vecRuns, 5, 5, kTRUE, addRight, 0.8, 0.92, 0.03, 0.8, 0.81, 0.03, doTrigger, fTrigger, (Bool_t)(i<nData), outputDirDataSet, "ClusterM20_Runwise", vecplotDataSets.at(i), kFALSE, fCollisionSystem, Form("%s clusters", calo.Data()), suffix, kFALSE, kTRUE, kFALSE);
                //--------------------------------------------------------------------------------------------------------
                DrawVectorRunwiseTH1D(	canvasRunwise, legendRuns, vecClusterNCells, vecRuns, 5, 5, kTRUE, addRight, 0.8, 0.92, 0.03, 0.8, 0.81, 0.03, doTrigger, fTrigger, (Bool_t)(i<nData), outputDirDataSet, "ClusterNCells_Runwise", vecplotDataSets.at(i), kFALSE, fCollisionSystem, Form("%s clusters", calo.Data()), suffix, kTRUE, kTRUE, kFALSE);
                //--------------------------------------------------------------------------------------------------------
                DrawVectorRunwiseTH1D(	canvasRunwise, legendRuns, vecClusterDispersion, vecRuns, 5, 5, kTRUE, addRight, 0.8, 0.92, 0.03, 0.8, 0.81, 0.03, doTrigger, fTrigger, (Bool_t)(i<nData), outputDirDataSet, "ClusterDispersion_Runwise", vecplotDataSets.at(i), kFALSE, fCollisionSystem, Form("%s clusters", calo.Data()), suffix, kFALSE, kTRUE, kFALSE);
                //--------------------------------------------------------------------------------------------------------
                TGaxis::SetExponentOffset(0.5, 0, "x");
                DrawVectorRunwiseTH1D(	canvasRunwise, legendRuns, vecClusterTime[0], vecRuns, 5, 5, kTRUE, addRight, 0.8, 0.92, 0.03, 0.8, 0.81, 0.03, doTrigger, fTrigger, (Bool_t)(i<nData), outputDirDataSet, "ClusterTime_Runwise", vecplotDataSets.at(i), kTRUE, fCollisionSystem, Form("%s clusters", calo.Data()), suffix, kFALSE, kTRUE, kFALSE);
                for(Int_t iT=1; iT<5; iT++){
                    DrawVectorRunwiseTH1D(	canvasRunwise, legendRuns, vecClusterTime[iT], vecRuns, 5, 5, kTRUE, addRight, 0.8, 0.92, 0.03, 0.8, 0.81, 0.03, doTrigger, fTrigger, (Bool_t)(i<nData), outputDirDataSet, Form("ClusterTime_Runwise_%.01f-%.01f",minT_Energy[iT], maxT_Energy[iT]), vecplotDataSets.at(i), kTRUE, fCollisionSystem, Form("%s clusters", calo.Data()), suffix, kFALSE, kTRUE, kFALSE);
                }
                TGaxis::SetExponentOffset(0, 0, "x");
                //--------------------------------------------------------------------------------------------------------
                if(isTrackMatching){
                    TGaxis::SetExponentOffset(0.013, -0.0285, "x");
                    DrawVectorRunwiseTH1D(	canvasRunwise, legendRuns, vecClusterDeltaEta, vecRuns, 5, 5, kTRUE, addRight, 0.8, 0.92, 0.03, 0.8, 0.81, 0.03, doTrigger, fTrigger, (Bool_t)(i<nData), outputDirDataSet, "ClusterDeltaEta_Runwise", vecplotDataSets.at(i), kTRUE, fCollisionSystem, Form("%s clusters", calo.Data()), suffix, kFALSE, kTRUE, kFALSE);
                    DrawVectorRunwiseTH1D(	canvasRunwise, legendRuns, vecClusterDeltaPhi, vecRuns, 5, 5, kTRUE, addRight, 0.8, 0.92, 0.03, 0.8, 0.81, 0.03, doTrigger, fTrigger, (Bool_t)(i<nData), outputDirDataSet, "ClusterDeltaPhi_Runwise", vecplotDataSets.at(i), kTRUE, fCollisionSystem, Form("%s clusters", calo.Data()), suffix, kFALSE, kTRUE, kFALSE);
                    TGaxis::SetExponentOffset(0, 0, "x");
                }
                //--------------------------------------------------------------------------------------------------------
                if(doExtQA>0){
                    for(Int_t iModule=0; iModule<nCaloModules; iModule++){
                        DrawVectorRunwiseTH1D(	canvasRunwise, legendRuns, vecClusterEVsModule[iModule], vecRuns, 5, 5, kTRUE, addRight, 0.8, 0.92, 0.03, 0.8, 0.81, 0.03, doTrigger, fTrigger, (Bool_t)(i<nData), outputDirDataSet, Form("ExtQA/ClusterEVsModule%i_Runwise",iModule), vecplotDataSets.at(i), kFALSE, fCollisionSystem, Form("%s clusters", calo.Data()), suffix, kTRUE, kTRUE, kFALSE);
                        DrawVectorRunwiseTH1D(	canvasRunwise, legendRuns, vecModuleEVsModule[iModule], vecRuns, 5, 5, kTRUE, addRight, 0.8, 0.92, 0.03, 0.8, 0.81, 0.03, doTrigger, fTrigger, (Bool_t)(i<nData), outputDirDataSet, Form("ExtQA/ModuleEVsModule%i_Runwise",iModule), vecplotDataSets.at(i), kFALSE, fCollisionSystem, Form("%s clusters", calo.Data()), suffix, kTRUE, kTRUE, kFALSE);
                        DrawVectorRunwiseTH1D(	canvasRunwise, legendRuns, vecClusterNCells100VsModule[iModule], vecRuns, 5, 5, kTRUE, addRight, 0.8, 0.92, 0.03, 0.8, 0.81, 0.03, doTrigger, fTrigger, (Bool_t)(i<nData), outputDirDataSet, Form("ExtQA/NCells100VsModule%i_Runwise",iModule), vecplotDataSets.at(i), kFALSE, fCollisionSystem, Form("%s clusters", calo.Data()), suffix, kFALSE, kTRUE, kFALSE);
                        DrawVectorRunwiseTH1D(	canvasRunwise, legendRuns, vecClusterNCells1500VsModule[iModule], vecRuns, 5, 5, kTRUE, addRight, 0.8, 0.92, 0.03, 0.8, 0.81, 0.03, doTrigger, fTrigger, (Bool_t)(i<nData), outputDirDataSet, Form("ExtQA/NCells1500VsModule%i_Runwise",iModule), vecplotDataSets.at(i), kFALSE, fCollisionSystem, Form("%s clusters", calo.Data()), suffix, kFALSE, kTRUE, kFALSE);
                    }
                }

                if(doExtQA>1){
                    for(Int_t iBad=0; iBad<(Int_t)vecBadCells.size(); iBad++){
                        DrawVectorRunwiseBadCells(canvasRunwise, legendRuns, vecBadCellsEnergy, (Int_t)vecBadCells.size(), iBad, vecRuns, 5, 5, kTRUE, addRight, 0.8, 0.92, 0.03, 0.8, 0.81, 0.03, doTrigger, fTrigger, (Bool_t)(i<nData), outputDirDataSet, Form("BadCells/Cell%i_Energy",((TString)vecBadCells.at(iBad)).Atoi()), vecplotDataSets.at(i), kTRUE, fCollisionSystem, Form("%s clusters", calo.Data()), suffix, kFALSE, kTRUE, kFALSE);
                        DrawVectorRunwiseBadCells(canvasRunwise, legendRuns, vecBadCellsTime, (Int_t)vecBadCells.size(), iBad, vecRuns, 5, 5, kTRUE, addRight, 0.8, 0.92, 0.03, 0.8, 0.81, 0.03, doTrigger, fTrigger, (Bool_t)(i<nData), outputDirDataSet, Form("BadCells/Cell%i_Time",((TString)vecBadCells.at(iBad)).Atoi()), vecplotDataSets.at(i), kTRUE, fCollisionSystem, Form("%s clusters", calo.Data()), suffix, kFALSE, kTRUE, kFALSE);
                    }
                }
                //--------------------------------------------------------------------------------------------------------
                delete legendRuns;
                delete canvasRunwise;
            }

            canvas->cd();


            // ****************************** Create Output ROOT-File ************************************************

            if (!runMergedClust){
                if(useDataRunListForMC && i>=nData) nameOutput = Form("%s/%s/ClusterQA/%s-%s_ClusterQARunwise.root",fCutSelection.Data(),fEnergyFlag.Data(),vecDataSet.at(i).Data(),vecDataSet.at(0).Data());
                else nameOutput = Form("%s/%s/ClusterQA/%s_ClusterQARunwise.root",fCutSelection.Data(),fEnergyFlag.Data(),vecDataSet.at(i).Data());
            } else {
                if(useDataRunListForMC && i>=nData) nameOutput = Form("%s/%s/MergedClusterQA/%s-%s_ClusterQARunwise.root",fCutSelection.Data(),fEnergyFlag.Data(),vecDataSet.at(i).Data(),vecDataSet.at(0).Data());
                else nameOutput = Form("%s/%s/MergedClusterQA/%s_ClusterQARunwise.root",fCutSelection.Data(),fEnergyFlag.Data(),vecDataSet.at(i).Data());
            }

            TFile* fOutput = new TFile(nameOutput,"RECREATE");
            printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
            printf("Output file: %s\n", nameOutput);
            printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
            fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
            fLog << "Output file: " << nameOutput << endl;
            fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

            for(Int_t h=0; h<(Int_t) vecHistos.size(); h++) WriteHistogram(((TH1D*) vecHistos.at(h)));

            WriteHistogramTH2DVec(fOutput,vecClusterEtaPhi,"EtaVsPhi");
            if(i<nData) WriteHistogramTH2DVec(fOutput,vecClusterEnergyTime,"EnergyVsTime");

            if(doExtQA>0){
                WriteHistogramTH2DVec(fOutput,vecCellTimeID,"CellTimeVsCellId");
                WriteHistogramTH2DVec(fOutput,vecClusterEVsNCells,"ClusterEVsNCells");
                WriteHistogramTH1DVec(fOutput,vecClusterIncludedCells,"IncludedCells");
                WriteHistogramTH1DVec(fOutput,vecClusterIncludedCellsBefore,"IncludedCells_beforeQA");
                WriteHistogramTProfileVec(fOutput,vecClusterBadCells,"BadCells");
                // if(i<nData) WriteHistogramTH2DVec(fOutput,vecClusterMissingCellIDs,"MissingCells");
                WriteHistogramTH1DVec(fOutput,vecClusterEFracCells,"EFracCells");
                WriteHistogramTH1DVec(fOutput,vecClusterEFracCellsBefore,"EFracCells_beforeQA");
                if(doExtQA==2){
                    WriteHistogramTH2DVec(fOutput,vecClusterEnergyMeanSigma,"ClusterEMeanVsSigma");
                    WriteHistogramTH2DVec(fOutput,vecClusterTimeMeanSigma,"ClusterTimeMeanVsSigma");
                    WriteHistogramTH2DVec(fOutput,vecClusterFiredCellIDs,"ClusterFiredCellIDs");
                }
            }

            for(Int_t iCM=0; iCM<nCaloModules; iCM++){
                WriteHistogramTH1DVec(fOutput,vecClusterEVsModule[iCM],Form("Mod%02i/ClusterE",iCM));
                WriteHistogramTH1DVec(fOutput,vecModuleEVsModule[iCM],Form("Mod%02i/ModuleE",iCM));
                WriteHistogramTH1DVec(fOutput,vecClusterNCells100VsModule[iCM],Form("Mod%02i/NCells100",iCM));
                WriteHistogramTH1DVec(fOutput,vecClusterNCells1500VsModule[iCM],Form("Mod%02i/NCells1500",iCM));
            }

            DeleteVecTH1D(vecClusterEnergy);
            DeleteVecTH1D(vecClusterM02);
            DeleteVecTH1D(vecClusterM20);
            DeleteVecTH1D(vecClusterNCells);
            DeleteVecTH1D(vecClusterDispersion);
            for(Int_t iT=0; iT<5; iT++) DeleteVecTH1D(vecClusterTime[iT]);
            DeleteVecTH1D(vecClusterDeltaEta);
            DeleteVecTH1D(vecClusterDeltaPhi);
            DeleteVecTH1D(vecClusterPi0ConvPhotonEta);
            DeleteVecTH1D(vecClusterPi0ConvPhotonPhi);
            DeleteVecTH1D(vecClusterEtaConvPhotonEta);
            DeleteVecTH1D(vecClusterEtaConvPhotonPhi);

            DeleteVecTH1D(vecBadCellsEnergy);
            DeleteVecTH1D(vecBadCellsTime);


            DeleteVecTH1D(vecHistos);
            vecHistos.clear();

            fOutput->Write();
            fOutput->Close();
            delete fOutput;

            fLog.close();
            if(doExtQA==2) fLogRunwiseHotCells.close();
            if(doExtQA==2) fLogRunwiseColdCells.close();
        } else printf("Skipping.... | Only Trending\n");

        if(doExtQA==2){
            std::vector<TString>  vec;
            std::map <Int_t,Int_t> ma;
            TH1D* nCellsRun;
            Int_t nRuns;
            TCanvas* cvsRun             = new TCanvas("canvas","",10,10,750,500);  // gives the page size
            DrawGammaCanvasSettings(cvsRun, leftMar, rightMar, 0.06, bottomMargin);

            nRuns             = 0;
            fLogRunwiseHotCells.open(Form("%s/HotCellsRunwise-%s.log",outputDir.Data(),vecDataSet.at(i).Data()), ios::in);
            if(fLogRunwiseHotCells.good()) {
                fLogRunwiseHotCells.seekg(0L, ios::beg);
                TString fVar;
                while(!fLogRunwiseHotCells.eof()) {
                    fLogRunwiseHotCells >> fVar;
                    if(fVar.BeginsWith("Run-")) {
                        nRuns++;
                        vec.push_back(fVar);
                    } else if(fVar.BeginsWith("NoNoisy")||fVar.BeginsWith("NotEnough")){
                        continue;
                    } else if(fVar.Sizeof()>1) {
                        TObjArray *rNumber = fVar.Tokenize("-");
                        TObjString* rString = (TObjString*)rNumber->At(0);
                        TString vecString = rString->GetString();
                        vec.push_back(vecString);
                        if( ma.find(vecString.Atoi()) != ma.end() ) ma[vecString.Atoi()] += 1;
                        else ma[vecString.Atoi()] = 1;
                    }
                }
            }

            if((Int_t)ma.size()==0){
                printf("No Bad Cells found for: %s\n", vecDataSet.at(i).Data());
                continue;
            }

            Int_t plotting = 0;
            Int_t nPlot = 100;
            Int_t nPart = 0;
            map<Int_t, Int_t>::iterator it = ma.begin();
            do {
                plotting+=nPlot;
                if(plotting>(Int_t)ma.size()){
                    plotting -= nPlot;
                    nPlot = ma.size()-plotting;
                    plotting = ma.size();
                }
                nCellsRun = new TH1D(Form("%s, %i runs, %s clusters",vecDataSet.at(i).Data(), nRuns, calo.Data()),Form("%s, %i runs, Runwise Hot Cells for %s; CellID; Number of Runs",vecDataSet.at(i).Data(), nRuns, calo.Data()),nPlot,0,nPlot);
                OnlyEditTH1(nCellsRun, 20, 1, kBlack, kBlack);

                for (Int_t h=1; h<=nPlot; it++,h++){
                    if(it == ma.end()) {printf("\033[0;31mERROR\033[0m while plotting HotCellsRunwiseOverview\n");break;}
                    nCellsRun->GetXaxis()->SetBinLabel(h,Form("%i",it->first));
                    nCellsRun->SetBinContent(h,it->second);
                }
                nCellsRun->SetFillStyle(3004);
                nCellsRun->SetFillColor(1);
                nCellsRun->Draw();
                if(useDataRunListForMC && i>=nData) SaveCanvas(cvsRun, Form("%s/%s-%s/%s_%i.%s", outputDir.Data(), vecDataSet.at(i).Data(),vecDataSet.at(0).Data(), "HotCells_FiredInNRuns",nPart++,suffix.Data()));
                else SaveCanvas(cvsRun, Form("%s/%s/%s_%i.%s", outputDir.Data(), vecDataSet.at(i).Data(), "HotCells_FiredInNRuns",nPart++,suffix.Data()));
                delete nCellsRun;
            } while(plotting < (Int_t)ma.size());

            plotting = 0;
            nPlot = 100;
            nPart = 0;
            Int_t posInVec = 0;
            std::vector<TString>::iterator itVec = vec.begin();
            do {
                plotting+=nPlot;
                if(plotting>nRuns){
                    plotting -= nPlot;
                    nPlot = nRuns-plotting;
                    plotting = nRuns;
                }
                nCellsRun = new TH1D(Form("%s, %i runs, %s clusters",vecDataSet.at(i).Data(), nRuns, calo.Data()),Form("%s, %i runs, Runwise Hot Cells for %s; Run Number; Number of Hot Cells",vecDataSet.at(i).Data(), nRuns, calo.Data()),nPlot,0,nPlot);
                OnlyEditTH1(nCellsRun, 20, 1, kBlack, kBlack);

                for (Int_t h=1; h<=nPlot; h++){
                    Int_t nBadCells = 0;
                    do{
                        if((vec.at(posInVec)).BeginsWith("Run-")){
                            TString tempStr = (vec.at(posInVec)).Remove(0,4);
                            nCellsRun->GetXaxis()->SetBinLabel(h,Form("%s",tempStr.Data()));
                        } else nBadCells++;
                    }while(++itVec!=vec.end() && !((vec.at(++posInVec)).BeginsWith("Run-")));
                    nCellsRun->SetBinContent(h,nBadCells);
                }
                nCellsRun->SetFillStyle(3004);
                nCellsRun->SetFillColor(1);
                nCellsRun->Draw();
                if(useDataRunListForMC && i>=nData) SaveCanvas(cvsRun, Form("%s/%s-%s/%s_%i.%s", outputDir.Data(), vecDataSet.at(i).Data(),vecDataSet.at(0).Data(), "HotCells_Runwise",nPart++,suffix.Data()));
                else SaveCanvas(cvsRun, Form("%s/%s/%s_%i.%s", outputDir.Data(), vecDataSet.at(i).Data(), "HotCells_Runwise",nPart++,suffix.Data()));
                delete nCellsRun;
            } while (plotting < nRuns);
            fLogRunwiseHotCells.close();

            delete cvsRun;
            // delete nRuns;
            vec.clear();
            // ma.clear();
        }


        vecMissingRuns[i].clear();


    }
    // */
    delete boxLabel;
    delete boxLabel2;
    delete[] vecMissingRuns;

    delete fDeltaPt;
    TH1::AddDirectory(kTRUE);


    std::vector<TH2D*> vecTMParr[nSets];
    std::vector<TH2D*> vecClusterMissingCellIDstmp;
    if(doExtQA==2){
        for(Int_t i=0; i<nData; i++) {
            if (!runMergedClust){
                if(useDataRunListForMC && i>=nData) nameOutput = Form("%s/%s/ClusterQA/%s-%s_ClusterQARunwise.root",fCutSelection.Data(),fEnergyFlag.Data(),vecDataSet.at(i).Data(),vecDataSet.at(0).Data());
                else nameOutput = Form("%s/%s/ClusterQA/%s_ClusterQARunwise.root",fCutSelection.Data(),fEnergyFlag.Data(),vecDataSet.at(i).Data());
            } else {
                if(useDataRunListForMC && i>=nData) nameOutput = Form("%s/%s/MergedClusterQA/%s-%s_ClusterQARunwise.root",fCutSelection.Data(),fEnergyFlag.Data(),vecDataSet.at(i).Data(),vecDataSet.at(0).Data());
                else nameOutput = Form("%s/%s/MergedClusterQA/%s_ClusterQARunwise.root",fCutSelection.Data(),fEnergyFlag.Data(),vecDataSet.at(i).Data());
            }
            TFile* fOutput = new TFile(nameOutput,"READ");
            fOutput->cd();
            TString dir = "ClusterFiredCellIDs";
            TDirectory *tmpDir = (TDirectory *)fOutput->Get(dir.Data());
            if(!readin(fileRuns.at(i), vecRuns)) {printf("\033[0;31mERROR\033[0m, no Run Numbers could be found! Returning...\n"); return;}
            for(Int_t j=0; j<(Int_t) vecRuns.size(); j++) {
                TH2D* tmpDeadCellsRunwise = (TH2D*)tmpDir->Get(Form("%s_ClusterIncludedCells_%s",vecRuns.at(j).Data(),vecDataSet.at(i).Data()));
                vecTMParr[i].push_back(tmpDeadCellsRunwise);
            }
            fOutput->Close();
        }
        for(Int_t i=0; i<nData; i++) {
            if(i < nData){
                for(Int_t iMC=nData; iMC<nData; iMC++){
                    if(vecTMParr[i].size() == vecTMParr[iMC].size()){
                        DrawVectorOverviewMissingCells(canvaslarge, vecTMParr[i], vecTMParr[iMC], vecClusterMissingCellIDstmp, Form("ExtQA/MissingCells/hCellsMissingData_%s",vecDataSet.at(iMC).Data()), outputDirDataSet, suffix, vecDataSet.at(iMC));
                    } else {
                        printf("\033[0;31mERROR\033[0m! MissingCell overviews will not be drawn, because vecTMParr[i].size() != vecTMParr[iMC].size() \n");
                    }
                }
            }
        }
    }




    // ************************************************************************************************************
    // ****************************** Combined Trending Histograms ************************************************
    // ************************************************************************************************************
    printf("\n\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
    printf("Drawing Trending Histograms\n");


    //--------------------------------------------------------------------------------------------------------
    if(useDataRunListForMC) printf("\033[0;35mWARNING\033[0m: useDataRunListForMC is true, overwriting histograms for DataSet!\n");

    TLegend *legend = new TLegend(0.15,0.95,0.95,0.98);
    legend->SetNColumns(nSets);
    legend->SetFillColor(0);
    legend->SetLineColor(0);
    legend->SetTextSize(0.03);
    legend->SetTextFont(42);
    //--------------------------------------------------------------------------------------------------------
    std::vector< TH1D* > vecHTMP;
    if(nSets > 1){
        TH1D* hTMP;
        for(Int_t h=0; h<(Int_t) vecHistosName.size(); h++){
            TString fTrigger        = "";
            TString fTriggerCut     = fEventCutSelection(3,2);
            printf("%s\n", vecHistosName.at(h).Data());
            for(Int_t i=0; i<nSets; i++){
            // for(Int_t i=(nSets)-1; i>=0; i--){
                if (!runMergedClust){
                    if(useDataRunListForMC && i>=nData) nameOutput = Form("%s/%s/ClusterQA/%s-%s_ClusterQARunwise.root",fCutSelection.Data(),fEnergyFlag.Data(),vecDataSet.at(i).Data(),vecDataSet.at(0).Data());
                    else nameOutput = Form("%s/%s/ClusterQA/%s_ClusterQARunwise.root",fCutSelection.Data(),fEnergyFlag.Data(),vecDataSet.at(i).Data());
                } else {
                    if(useDataRunListForMC && i>=nData) nameOutput = Form("%s/%s/MergedClusterQA/%s-%s_ClusterQARunwise.root",fCutSelection.Data(),fEnergyFlag.Data(),vecDataSet.at(i).Data(),vecDataSet.at(0).Data());
                    else nameOutput = Form("%s/%s/MergedClusterQA/%s_ClusterQARunwise.root",fCutSelection.Data(),fEnergyFlag.Data(),vecDataSet.at(i).Data());
                }
                // printf("%s\n", nameOutput);
                TFile* fOutput = new TFile(nameOutput,"READ");
                fOutput->cd();
                hTMP = (TH1D*)(fOutput->Get(Form("%s_%s",vecHistosName.at(h).Data(),vecDataSet.at(i).Data())))->Clone(Form("%s_%s_%i",vecHistosName.at(h).Data(),vecDataSet.at(i).Data(),i));

                if(!hTMP) {
                    printf("\n%s-%s\t\033[0;31mERROR\033[0m not found \n", vecDataSet.at(i).Data(), vecHistosName.at(h).Data());
                    continue;
                }

                hTMP->SetDirectory(0);
                vecHTMP.push_back(hTMP);
                printf("%s, ", vecDataSet.at(i).Data());

                fOutput->Close();
            }

            if( ((TString)vecHistosName.at(h)).CompareTo("hNEvents")==0 || ((TString)vecHistosName.at(h)).CompareTo("hClusterTime-Mean")==0 )
            AdjustHistRange(vecHTMP,10.,10.,nSets,kTRUE);
            else AdjustHistRange(vecHTMP,1.1,1.1,nSets,kTRUE);

            for(Int_t i=(nSets)-1; i>=0; i--){
                if(vecSubRunlists.at(i) == 0) legend->AddEntry(((TH1D*) vecHTMP.at(i)),vecplotDataSets.at(i).Data(),"p");
                TString draw;
                if(h==0) draw = (i==(nSets)-1)?"p":"p, same";
                else draw = (i==(nSets)-1)?"px0e1":"px0e1, same";
                ((TH1D*) vecHTMP.at(i))->SetTitle("");
                canvaswide->cd();
                ((TH1D*) vecHTMP.at(i))->Draw(draw.Data());
                canvas->cd();
                ((TH1D*) vecHTMP.at(i))->Draw(draw.Data());
            }

            canvaswide->cd();
            legend->Draw();
            canvas->cd();
            legend->Draw();

            if(doTrigger && 0<nData){
                fTrigger        = AnalyseSpecialTriggerCut(fTriggerCut.Atoi(), vecplotDataSets.at(0));
                if(fTrigger.Contains("not defined")) fTrigger = "";
            }

            if(doTrigger) PutProcessLabelAndEnergyOnPlot(xPosLabel, 0.92, 0.03, fCollisionSystem.Data(), fTrigger.Data(), Form("%s clusters", calo.Data()));
            else PutProcessLabelAndEnergyOnPlot(xPosLabel, 0.92, 0.03, fCollisionSystem.Data(), Form("%s clusters", calo.Data()), "");

            if(canvas->GetTopMargin()!=0.06) canvas->SetTopMargin(0.06);
            if(canvaswide->GetTopMargin()!=0.06) canvaswide->SetTopMargin(0.06);

            if(useDataRunListForMC && !addSubFolder){
                if( ((TString)vecHistosName.at(h)).CompareTo("hNEvents")==0 ||
                ((TString)vecHistosName.at(h)).CompareTo("hClusterTime-Mean")==0 ) {
                    SaveCanvas(canvas, Form("%s/%s/%s.%s", outputDir.Data(), vecDataSet.at(0).Data(),vecHistosName.at(h).Data(),suffix.Data()), kFALSE, kTRUE);
                    SaveCanvas(canvaswide, Form("%s/%s/%s_wide.%s", outputDir.Data(), vecDataSet.at(0).Data(),vecHistosName.at(h).Data(),suffix.Data()), kFALSE, kTRUE);
                }
                else {
                    SaveCanvas(canvas, Form("%s/%s/%s.%s", outputDir.Data(), vecDataSet.at(0).Data(),vecHistosName.at(h).Data(),suffix.Data()));
                    SaveCanvas(canvaswide, Form("%s/%s/%s_wide.%s", outputDir.Data(), vecDataSet.at(0).Data(),vecHistosName.at(h).Data(),suffix.Data()));
                }
            } else{
                if( ((TString)vecHistosName.at(h)).CompareTo("hNEvents")==0 ||
                ((TString)vecHistosName.at(h)).CompareTo("hClusterTime-Mean")==0 ){
                    SaveCanvas(canvas, Form("%s/%s.%s", outputDir.Data(), vecHistosName.at(h).Data(),suffix.Data()), kFALSE, kTRUE);
                    SaveCanvas(canvaswide, Form("%s/%s_wide.%s", outputDir.Data(), vecHistosName.at(h).Data(),suffix.Data()), kFALSE, kTRUE);
                }
                else {
                    SaveCanvas(canvas, Form("%s/%s.%s", outputDir.Data(), vecHistosName.at(h).Data(),suffix.Data()));
                    SaveCanvas(canvaswide, Form("%s/%s_wide.%s", outputDir.Data(), vecHistosName.at(h).Data(),suffix.Data()));
                }
            }
            legend->Clear();
            printf("\n----------------------------\n");

            DeleteVecTH1D(vecHTMP);
            vecHTMP.clear();
        }
    }
    // /*
    // // ******************************************************************************************************************
    // // ****************************** Combined Ratio Trending Histograms ************************************************
    // // ******************************************************************************************************************
    printf("\n\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
    printf("Drawing Ratio Trending Histograms\n");

    if(doHistsForEverySet) {
        TString* ratioSets[nSets-nData];
        TH1D* ratio[nSets-nData];
        TString outputDirDataSet;
        TH1D* hTMP;
        if(nSets>1 && nSets>nData) {

            legend->SetNColumns(nData*(nSets-nData));
            Int_t markerStyles[14]={2,4,5,20,21,22,23,24,25,26,27,28,29,30};

            for(Int_t h=0; h<(Int_t) vecHistosName.size(); h++) {
                printf("%s\n", vecHistosName.at(h).Data());
                for(Int_t i=0; i<nSets; i++){
                    if (!runMergedClust){
                        if(useDataRunListForMC && i>=nData) nameOutput = Form("%s/%s/ClusterQA/%s-%s_ClusterQARunwise.root",fCutSelection.Data(),fEnergyFlag.Data(),vecDataSet.at(i).Data(),vecDataSet.at(0).Data());
                        else nameOutput = Form("%s/%s/ClusterQA/%s_ClusterQARunwise.root",fCutSelection.Data(),fEnergyFlag.Data(),vecDataSet.at(i).Data());
                    } else {
                        if(useDataRunListForMC && i>=nData) nameOutput = Form("%s/%s/MergedClusterQA/%s-%s_ClusterQARunwise.root",fCutSelection.Data(),fEnergyFlag.Data(),vecDataSet.at(i).Data(),vecDataSet.at(0).Data());
                        else nameOutput = Form("%s/%s/MergedClusterQA/%s_ClusterQARunwise.root",fCutSelection.Data(),fEnergyFlag.Data(),vecDataSet.at(i).Data());
                    }
                    // printf("%s\n", nameOutput);

                    TFile* fOutput = new TFile(nameOutput,"READ");
                    fOutput->cd();
                    hTMP = (TH1D*)(fOutput->Get(Form("%s_%s",vecHistosName.at(h).Data(),vecDataSet.at(i).Data())))->Clone(Form("%s_%s_%i",vecHistosName.at(h).Data(),vecDataSet.at(i).Data(),i));
                    if(!hTMP) {
                        printf("\n%s-%s\t\033[0;31mERROR\033[0m not found \n", vecDataSet.at(i).Data(), vecHistosName.at(h).Data());
                        continue;
                    }
                    hTMP->SetDirectory(0);
                    vecHTMP.push_back(hTMP);
                    printf("%s, ", vecDataSet.at(i).Data());


                    fOutput->Close();

                }

                if( ((TString)vecHistosName.at(h)).CompareTo("hNEvents")==0 || ((TString)vecHistosName.at(h)).CompareTo("hClusterTime-Mean")==0 )
                AdjustHistRange(vecHTMP,10.,10.,nSets,kTRUE);
                else AdjustHistRange(vecHTMP,1.1,1.1,nSets,kTRUE);

                for(Int_t i=0; i<nData; i++) {
                    for(Int_t j=0; j<nSets-nData; j++) {
                        ratioSets[j] = new TString(Form("%s / %s", vecplotDataSets.at(i).Data(), vecplotDataSets.at(j+nData).Data()));
                        ratio[j] = new TH1D(Form("%s%i%i",((TH1D*) vecHTMP.at(i))->GetName(),i,j),
                        Form("%s%i%i;%s;%s - Ratio: Data / MC",((TH1D*) vecHTMP.at(i))->GetTitle(),i,j,((TH1D*) vecHTMP.at(i))->GetXaxis()->GetTitle(),((TH1D*) vecHTMP.at(i))->GetYaxis()->GetTitle()),
                        hNBin,hFBin,hLBin);
                        EditTH1(globalRuns, doEquidistantXaxis, ratio[j], GetDefaultMarkerStyle(fEnergyFlag.Data(),vecDataSet.at(j % 14).Data(),""), markerSize, GetColorDefaultColor(fEnergyFlag.Data(),vecDataSetIn.at(i).Data(),"",kFALSE,kFALSE), GetColorDefaultColor(fEnergyFlag.Data(),vecDataSetIn.at(i).Data(),"",kFALSE,kFALSE));
                    }

                    for(Int_t b=1; b<hNBin+1; b++) {
                        Double_t hData = ((TH1D*) vecHTMP.at(i))->GetBinContent(b);
                        for(Int_t j=0; j<nSets-nData; j++) {
                            Double_t hMC = ((TH1D*) vecHTMP.at(j+nData))->GetBinContent(b);
                            if(hMC!=0) {
                                if(hData/hMC>1.98) ratio[j]->SetBinContent(b,1.98);
                                else if(hData/hMC<0.02) ratio[j]->SetBinContent(b,0.02);
                                else ratio[j]->SetBinContent(b,hData/hMC);
                            } else ratio[j]->SetBinContent(b,1.98);
                        }
                    }


                    for(Int_t j=0; j<nSets-nData; j++) {
                        TString draw = (i==0&&j==0)?"p":"p, same";
                        ratio[j]->SetTitle("");
                        ratio[j]->GetYaxis()->SetRangeUser(0,2);
                        canvaswide->cd();
                        ratio[j]->Draw(draw.Data());
                        canvas->cd();
                        ratio[j]->Draw(draw.Data());
                        legend->AddEntry(ratio[j],ratioSets[j]->Data(),"p");
                    }

                    canvaswide->cd();
                    legend->Draw();
                    canvas->cd();
                    legend->Draw();
                    outputDirDataSet = Form("%s/%s",outputDir.Data(), vecDataSet.at(i).Data());
                    gSystem->Exec("mkdir -p "+outputDirDataSet+"/TrendingRatios");

                    if(doTrigger){
                        TString fTrigger            = "";
                        TString fTriggerCut         = fEventCutSelection(3,2);
                        fTrigger                    = AnalyseSpecialTriggerCut(fTriggerCut.Atoi(), vecplotDataSets.at(i));
                        if(fTrigger.Contains("not defined"))
                        fTrigger                = "";
                        canvaswide->cd();
                        PutProcessLabelAndEnergyOnPlot(xPosLabel, 0.92, 0.03, fCollisionSystem.Data(), fTrigger.Data(), Form("%s clusters", calo.Data()));
                        canvas->cd();
                        PutProcessLabelAndEnergyOnPlot(xPosLabel, 0.92, 0.03, fCollisionSystem.Data(), fTrigger.Data(), Form("%s clusters", calo.Data()));
                    } else
                    canvaswide->cd();
                    PutProcessLabelAndEnergyOnPlot(xPosLabel, 0.92, 0.03, fCollisionSystem.Data(), Form("%s clusters", calo.Data()), "");
                    canvas->cd();
                    PutProcessLabelAndEnergyOnPlot(xPosLabel, 0.92, 0.03, fCollisionSystem.Data(), Form("%s clusters", calo.Data()), "");


                    SaveCanvas(canvas, Form("%s/TrendingRatios/%s.%s", outputDirDataSet.Data(),Form("%s",((TH1D*) vecHTMP.at(i))->GetName()),suffix.Data()));
                    SaveCanvas(canvaswide, Form("%s/TrendingRatios/%s_wide.%s", outputDirDataSet.Data(),Form("%s",((TH1D*) vecHTMP.at(i))->GetName()),suffix.Data()));
                    legend->Clear();
                    for(Int_t j=0; j<nSets-nData; j++){
                        delete ratio[j];
                        delete ratioSets[j];
                    }

                }
                printf("\n----------------------------\n");
                DeleteVecTH1D(vecHTMP);
                vecHTMP.clear();
            }
        } else printf("...skipped due to nSets<=1 or nSets==nData!\n");
    }
    // */
    printf("\n\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
    printf("Done with ClusterQA_Runwise_V2\n");
    printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
    delete legend;
    delete canvas;
    delete canvaslarge;
    return;

}//end
