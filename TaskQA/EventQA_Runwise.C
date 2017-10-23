/*******************************************************************************
 ******  provided by Gamma Conversion Group, PWGGA,                        *****
 ******     Daniel Muehlheim, d.muehlheim@cern.ch                          *****
 ******     Friederike Bock, fbock@cern.ch                                 *****
 *******************************************************************************/

#include "QA.h"

//**************************************************************************************************************
//***************************** Main routine *******************************************************************
//**************************************************************************************************************
void EventQA_Runwise(
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
			Int_t *nSigmasBadRun            = NULL              // array of 8 integers 
                    )
{
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "EventQA_Runwise" << endl;
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
    if(nSets>maxSets){
        cout << "Maximum hardcoded number of Data Sets: " << maxSets << endl;
        cout << "You have chosen: " << nSets << ", returning!" << endl;
        return;
    }

    Int_t fMode         = mode;
    Bool_t isCalo       = kFALSE;
    Bool_t isMerged     = kFALSE;
    Bool_t isConv       = kFALSE;
    // mode:    0 // new output PCM-PCM
    //          1 // new output PCM dalitz
    //          2 // new output PCM-EMCal
    //          3 // new output PCM-PHOS
    //          4 // new output EMCal-EMCal
    //          5 // new output PHOS-PHOS
    //          9 // old output PCM-PCM
    //          10 // merged EMCal
    //          11 // merged PHOS
    if (fMode == 0 || fMode == 1 || fMode == 2 || fMode == 3 || fMode == 9)
        isConv          = kTRUE;
    if (fMode == 2 || fMode == 3 || fMode == 4 || fMode == 5 || fMode == 10 || fMode == 11)
        isCalo          = kTRUE;
    if (fMode == 10 || fMode == 11 )
        isMerged        = kTRUE;


    TString fDate               = ReturnDateString();
    TString fTextMeasurement    = Form("#pi^{0} #rightarrow #gamma#gamma");
    TString fCentrality[30];
    for(Int_t i=0; i<nSets; i++) {
      if(fEnergyFlag.Contains("PbPb")){
        if(plotDataSets[i].Contains("0-10%")) fCentrality[i] = "0-10%";
        else if(plotDataSets[i].Contains("10-20%")) fCentrality[i] = "10-20%";
        else if(plotDataSets[i].Contains("20-50%")) fCentrality[i] = "20-50%";
        else if(plotDataSets[i].Contains("50-90%")) fCentrality[i] = "50-90%";
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
    if(fEnergyFlag.Contains("PbPb")){
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

    //*************************************************************************************************************
    //****************************** Determine which cut to process ***********************************************
    //*************************************************************************************************************
    TString fileRuns            = Form("%s/runNumbers%s.txt", folderRunlists.Data(), (vecDataSet.at(0)).Data());
    if(!readin(fileRuns, vecRuns, kFALSE)) {cout << "ERROR, no Run Numbers could be found! Returning..., tried to find: " << fileRuns.Data() << endl; return;}
    TFile* fCutFile             = new TFile(Form("%s/%s/%s/%s", filePath.Data(), ((TString)vecDataSet.at(0)).Data(), ((TString)vecRuns.at(0)).Data(), fileName.Data()));
    if(fCutFile->IsZombie()) {cout << "ERROR: ROOT file '" << Form("%s/%s/%s/%s", filePath.Data(), ((TString)vecDataSet.at(0)).Data(), ((TString)vecRuns.at(0)).Data(), fileName.Data()) << "' could not be openend, return!" << endl; return;}

    TKey *key;
    TIter next(fCutFile->GetListOfKeys());
    while ((key=(TKey*)next())){
        cout << Form("Found TopDir: '%s' ",key->GetName());
        nameMainDir             = key->GetName();
    }
    cout << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    if(nameMainDir.IsNull() || !nameMainDir.BeginsWith("Gamma")){cout << "ERROR, Unable to obtain valid name of MainDir:|" << nameMainDir.Data() << "|, running in mode: " << fMode << endl; return;}

    TList *listInput            = (TList*)fCutFile->Get(nameMainDir.Data());
    if(!listInput) {cout << "ERROR: Could not find main dir: " << nameMainDir.Data() << " in file! Returning..." << endl; return;}
        listInput->SetOwner(kTRUE);
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

    cout << "The following cuts are available:" << endl;
    for(Int_t i = 0; i < (Int_t) cuts.size(); i++) {cout << Form("(%i) -- %s", i, cuts[i].Data()) << endl;}

    if(cutNr == -1){
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


    //*****************************************************************************************************
    //************************** Set proper cluster nomenclature ******************************************
    //*****************************************************************************************************

    TString calo                = "";
    TString fClusters           = "";

    if(isCalo){
        if(fClusterCutSelection.BeginsWith('1')){
            calo                = "EMCal";
            fClusters           = Form("%s clusters", calo.Data());
        }else if(fClusterCutSelection.BeginsWith('2')){
            calo                = "PHOS";
            fClusters           = Form("%s clusters", calo.Data());
        }else {cout << "No correct calorimeter type found: " << calo.Data() << ", returning..." << endl; return;}
    }

    //*****************************************************************************************************
    //************************** Define output directories*************************************************
    //*****************************************************************************************************
    TString outputDir           = Form("%s/%s/EventQA/%s/Runwise",fCutSelection.Data(),fEnergyFlag.Data(),suffix.Data());
    if(addSubFolder){
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
        fileRuns = Form("%s/runNumbers%s.txt", folderRunlists.Data(), vecDataSet.at(i).Data());
        if(useDataRunListForMC && i>=nData) {
            fileRuns = Form("%s/runNumbers%s-%s.txt", folderRunlists.Data(), vecDataSet.at(i).Data(),vecDataSet.at(0).Data());
            cout << "Switch useDataRunListForMC is true, reading runs from: " << fileRuns.Data() << endl;
        }
        if(!readin(fileRuns, vecRuns, kFALSE)) {cout << "ERROR, no Run Numbers could be found! Returning..." << endl; return;}

        for(Int_t j=0; j<(Int_t) vecRuns.size();j++){
            if( i==0 ) globalRuns.push_back(vecRuns.at(j));
            else {
                Bool_t bFound = kFALSE;
                for(Int_t k=0; k<(Int_t) globalRuns.size();k++){ if(globalRuns.at(k)==vecRuns.at(j)) bFound=kTRUE;}
                if(!bFound) globalRuns.push_back(vecRuns.at(j));
            }
        }

        if( !doEquidistantXaxis && ((Int_t) vecRuns.size())>0 ){
            if(nSets>2){
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
        if(i%10==0) cout << endl;
        cout << globalRuns.at(i) << ", ";
    }
    cout << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

    TH1D* hNEvents[nSets];
    TH1D* hNEventsMinBias[nSets];
    TH1D* hNEventsAll[nSets];
    TH1D* hNEventsFracGoodEvents[nSets];
    TH1D* hNEventsFracNormAll[nSets];
    TH1D* hNEventsFracMinBias[nSets];

    TH1D* hTracksMeanGood[nSets];
    TH1D* hTracksRMSGood[nSets];

    TH1D* hVertexZMean[nSets];
    TH1D* hVertexZRMS[nSets];
    TH1D* hCentralityMean[nSets];
    TH1D* hEventPlaneAngleMean[nSets];
    TH1D* hFracWVtxOutside10cm[nSets];
    TH1D* hFracWOVtx[nSets];
    TH1D* hFracPileUp[nSets];
    TH1D* hFracSPDClusTrack[nSets];

    TH1D* hConvNCandidates[nSets];
    TH1D* hConvNCandidatesQA[nSets];

    TH1D* hCaloNClusters[nSets];
    TH1D* hCaloNClustersQA[nSets];

    TH1D* hCaloMergedNClusters[nSets];
    TH1D* hCaloMergedNClustersQA[nSets];

    TH1D* hPi0Frac[nSets];
    TH1D* hPi0Mass[nSets];
    TH1D* hPi0Width[nSets];
    TH1D* hPi0Pt[nSets];
    TH1D* hPi0PtRMS[nSets];
    TH1D* hPi0Alpha[nSets];
    TH1D* hPi0AlphaRMS[nSets];
    TH1D* hPi0Y[nSets];
    TH1D* hPi0YRMS[nSets];
    TH1D* hPi0OpenAngle[nSets];
    TH1D* hPi0OpenAngleRMS[nSets];

    TH1D* hEtaFrac[nSets];
    TH1D* hEtaMass[nSets];
    TH1D* hEtaWidth[nSets];
    TH1D* hEtaPt[nSets];
    TH1D* hEtaPtRMS[nSets];
    TH1D* hEtaAlpha[nSets];
    TH1D* hEtaAlphaRMS[nSets];
    TH1D* hEtaY[nSets];
    TH1D* hEtaYRMS[nSets];
    TH1D* hEtaOpenAngle[nSets];
    TH1D* hEtaOpenAngleRMS[nSets];

    //*****************************************************************************************************
    //*****************************************************************************************************
    //****************************** Vectors for Histograms ***********************************************
    //*****************************************************************************************************
    //*****************************************************************************************************

    std::vector<TH1D*>* vecVertexZ      = new std::vector<TH1D*>[nSets];
    std::vector<TH1D*>* vecPi0Signal    = new std::vector<TH1D*>[nSets];
    std::vector<TH1D*>* vecEtaSignal    = new std::vector<TH1D*>[nSets];
    std::vector<TF1*>* vecPi0Fit        = new std::vector<TF1*>[nSets];
    std::vector<TF1*>* vecEtaFit        = new std::vector<TF1*>[nSets];

    std::vector<TH1D*>* vecVertexZRatio = new std::vector<TH1D*>[nSets-1];
    std::map<Int_t,Int_t> mapVertexRatio;
    std::vector<TH1D*>* vecHistos       = new std::vector<TH1D*>[nSets];
    std::vector<TString> vecHistosName;
    TString histoName;


    Int_t hFBin;
    Int_t hLBin;
    Int_t hNBin;

    if(doEquidistantXaxis)    {
        hFBin       = 0;
        hLBin       = globalRuns.size();
        hNBin       = globalRuns.size();
    } else {
        if(nSets>2){
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
        histoName               = "hNEvents";
        if(i==0) vecHistosName.push_back(histoName);
        hNEvents[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hNEvents; Run Number ; N^{Evt}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hNEvents[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hNEvents[i]);

        histoName                   = "hNEventsMinBias";
        if(i==0) vecHistosName.push_back(histoName);
        hNEventsMinBias[i]          = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hNEventsMinBias; Run Number ; N^{MinBias}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hNEventsMinBias[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hNEventsMinBias[i]);

        histoName                   = "hNEventsAll";
        if(i==0) vecHistosName.push_back(histoName);
        hNEventsAll[i]              = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hNEventsAll; Run Number ; N^{All}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hNEventsAll[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hNEventsAll[i]);

        histoName                   = "hNEventsFracGoodEvents";
        if(i==0) vecHistosName.push_back(histoName);
        hNEventsFracGoodEvents[i]   = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hNEventsFracGoodEvents; Run Number ; N^{good Evt}/N^{MinBias}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hNEventsFracGoodEvents[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hNEventsFracGoodEvents[i]);

        histoName                   = "hNEventsFracNormAll";
        if(i==0) vecHistosName.push_back(histoName);
        hNEventsFracNormAll[i]      = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hNEventsFracNormAll; Run Number ; N_{Norm}^{Evt}/N^{All}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hNEventsFracNormAll[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hNEventsFracNormAll[i]);

        histoName                   = "hNEventsFracMinBias";
        if(i==0) vecHistosName.push_back(histoName);
        hNEventsFracMinBias[i]      = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hNEventsFracMinBias; Run Number ; N_{Norm}^{Evt}/N^{MinBias}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hNEventsFracMinBias[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hNEventsFracMinBias[i]);

        histoName                   = "hTracksGood-Mean";
        if(i==0) vecHistosName.push_back(histoName);
        hTracksMeanGood[i]          = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hTracksMeanGood; Run Number ; #bar{#lower[0.1]{N}}_{Good Tracks}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hTracksMeanGood[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hTracksMeanGood[i]);

        histoName                   = "hTracksGood-RMS";
        if(i==0) vecHistosName.push_back(histoName);
        hTracksRMSGood[i]           = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hTracksRMSGood; Run Number ; #sigma_{N_{Good Tracks}}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hTracksRMSGood[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hTracksRMSGood[i]);

        histoName                   = "hVertexZ-Mean";
        if(i==0) vecHistosName.push_back(histoName);
        hVertexZMean[i]             = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hVertexZ-Mean; Run Number ; #bar{#lower[0.1]{z}}_{Vertex} (cm)",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hVertexZMean[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hVertexZMean[i]);

        histoName                   = "hVertexZ-RMS";
        if(i==0) vecHistosName.push_back(histoName);
        hVertexZRMS[i]              = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hVertexZ-RMS; Run Number ; #sigma_{z-vertex} (cm)",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hVertexZRMS[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hVertexZRMS[i]);

	if(fEnergyFlag.Contains("PbPb")){
	  histoName                   = "hCentrality-Mean";
	  if(i==0) vecHistosName.push_back(histoName);
	  hCentralityMean[i]             = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hCentrality-Mean; Run Number ; #bar{cent} (%)",hNBin,hFBin,hLBin);
	  EditTH1(globalRuns, doEquidistantXaxis, hCentralityMean[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
	  vecHistos[i].push_back(hCentralityMean[i]);

	  histoName                   = "hEventPlaneAngle-Mean";
	  if(i==0) vecHistosName.push_back(histoName);
	  hEventPlaneAngleMean[i]     = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hEventPlaneAngle-Mean; Run Number ; #bar{e.p.a.} (rad)",hNBin,hFBin,hLBin);
	  EditTH1(globalRuns, doEquidistantXaxis, hEventPlaneAngleMean[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
	  vecHistos[i].push_back(hEventPlaneAngleMean[i]);
	}

        histoName                   = "hFracWVtxOutside10cm";
        if(i==0) vecHistosName.push_back(histoName);
        hFracWVtxOutside10cm[i]     = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hFracWVtxOutside10cm; Run Number ; N^{zVTX>10cm}/N^{Evt}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hFracWVtxOutside10cm[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hFracWVtxOutside10cm[i]);

        histoName                   = "hFracWOVtx";
        if(i==0) vecHistosName.push_back(histoName);
        hFracWOVtx[i]               = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hFracWOVtx; Run Number ; N^{without VTX}/N^{Evt}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hFracWOVtx[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hFracWOVtx[i]);

        histoName                   = "hFracPileUp";
        if(i==0) vecHistosName.push_back(histoName);
        hFracPileUp[i]              = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hFracPileUp; Run Number ; N^{pileup}/N^{Evt}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hFracPileUp[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hFracPileUp[i]);

        histoName                   = "hFracSPDClusTrack";
        if(i==0) vecHistosName.push_back(histoName);
        hFracSPDClusTrack[i]        = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hFracSPDClusTrack; Run Number ; N^{background}_{SPD Cluster vs Tracklet}/N^{Evt}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hFracSPDClusTrack[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hFracSPDClusTrack[i]);

        if (isConv){
            histoName               = "hConvNCandidates";
            if(i==0) vecHistosName.push_back(histoName);
            hConvNCandidates[i]     = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hConvNCandidates; Run Number ; #frac{1}{N_{Events}} N_{#gamma_{conv}, Candidates before Cuts}",hNBin,hFBin,hLBin);
            EditTH1(globalRuns, doEquidistantXaxis, hConvNCandidates[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
            vecHistos[i].push_back(hConvNCandidates[i]);

            histoName               = "hConvNCandidatesQA";
            if(i==0) vecHistosName.push_back(histoName);
            hConvNCandidatesQA[i]   = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hConvNCandidatesQA; Run Number ; #frac{1}{N_{Events}} N_{#gamma_{conv}, Candidates}",hNBin,hFBin,hLBin);
            EditTH1(globalRuns, doEquidistantXaxis, hConvNCandidatesQA[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
            vecHistos[i].push_back(hConvNCandidatesQA[i]);
        }

        histoName                   = "hCaloNClusters";
        if(i==0) vecHistosName.push_back(histoName);
        hCaloNClusters[i]           = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hCaloNClusters; Run Number ; #frac{1}{N_{Events}} N_{Clusters before Cuts}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hCaloNClusters[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hCaloNClusters[i]);

        histoName                   = "hCaloNClustersQA";
        if(i==0) vecHistosName.push_back(histoName);
        hCaloNClustersQA[i]         = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hCaloNClustersQA; Run Number ; #frac{1}{N_{Events}} N_{Clusters}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hCaloNClustersQA[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hCaloNClustersQA[i]);

        if (isMerged){
            histoName               = "hCaloMergedNClusters";
            if(i==0) vecHistosName.push_back(histoName);
            hCaloMergedNClusters[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hCaloMergedNClusters; Run Number ; #frac{1}{N_{Events}} N_{merged clusters before Cuts}",hNBin,hFBin,hLBin);
            EditTH1(globalRuns, doEquidistantXaxis, hCaloMergedNClusters[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
            vecHistos[i].push_back(hCaloMergedNClusters[i]);

            histoName               = "hCaloMergedNClustersQA";
            if(i==0) vecHistosName.push_back(histoName);
            hCaloMergedNClustersQA[i]= new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hCaloMergedNClustersQA; Run Number ; #frac{1}{N_{Events}} N_{merged clusters}",hNBin,hFBin,hLBin);
            EditTH1(globalRuns, doEquidistantXaxis, hCaloMergedNClustersQA[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
            vecHistos[i].push_back(hCaloMergedNClustersQA[i]);
        }

        if (!isMerged){
            histoName               = "hPi0Frac";
            if(i==0) vecHistosName.push_back(histoName);
            hPi0Frac[i]             = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hPi0Frac; Run Number ; N^{#pi^{0}}/N^{Evt}",hNBin,hFBin,hLBin);
            EditTH1(globalRuns, doEquidistantXaxis, hPi0Frac[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
            vecHistos[i].push_back(hPi0Frac[i]);

            histoName               = "hPi0Mass";
            if(i==0) vecHistosName.push_back(histoName);
            hPi0Mass[i]             = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hPi0Mass; Run Number ; m_{#pi^{0}}",hNBin,hFBin,hLBin);
            EditTH1(globalRuns, doEquidistantXaxis, hPi0Mass[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
            vecHistos[i].push_back(hPi0Mass[i]);

            histoName               = "hPi0Width";
            if(i==0) vecHistosName.push_back(histoName);
            hPi0Width[i]            = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hPi0Width; Run Number ; #sigma_{#pi^{0}}",hNBin,hFBin,hLBin);
            EditTH1(globalRuns, doEquidistantXaxis, hPi0Width[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
            vecHistos[i].push_back(hPi0Width[i]);

        }
        histoName               = "hPi0Pt";
        if(i==0) vecHistosName.push_back(histoName);
        hPi0Pt[i]               = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hPi0Pt; Run Number ; Mean #it{p}_{T}^{#pi^{0}}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hPi0Pt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hPi0Pt[i]);

        histoName               = "hPi0PtRMS";
        if(i==0) vecHistosName.push_back(histoName);
        hPi0PtRMS[i]            = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hPi0PtRMS; Run Number ; #sigma_{#it{p}_{T}^{#pi^{0}}}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hPi0PtRMS[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hPi0PtRMS[i]);

        histoName                   = "hPi0Alpha";
        if(i==0) vecHistosName.push_back(histoName);
        hPi0Alpha[i]                = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hPi0Alpha; Run Number ; Mean #alpha_{#pi^{0}}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hPi0Alpha[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hPi0Alpha[i]);

        histoName                   = "hPi0AlphaRMS";
        if(i==0) vecHistosName.push_back(histoName);
        hPi0AlphaRMS[i]             = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hPi0AlphaRMS; Run Number ; #sigma_{#alpha_{#pi^{0}}}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hPi0AlphaRMS[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hPi0AlphaRMS[i]);

        histoName                   = "hPi0Y";
        if(i==0) vecHistosName.push_back(histoName);
        hPi0Y[i]                    = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hPi0Y; Run Number ; Mean Y_{#pi^{0}}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hPi0Y[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hPi0Y[i]);

        histoName                   = "hPi0YRMS";
        if(i==0) vecHistosName.push_back(histoName);
        hPi0YRMS[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hPi0YRMS; Run Number ; #sigma_{Y_{#pi^{0}}}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hPi0YRMS[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hPi0YRMS[i]);

        histoName                   = "hPi0OpenAngle";
        if(i==0) vecHistosName.push_back(histoName);
        hPi0OpenAngle[i]            = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hPi0OpenAngle; Run Number ; Mean #theta_{#pi^{0}}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hPi0OpenAngle[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hPi0OpenAngle[i]);

        histoName                   = "hPi0OpenAngleRMS";
        if(i==0) vecHistosName.push_back(histoName);
        hPi0OpenAngleRMS[i]         = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hPi0OpenAngleRMS; Run Number ; #sigma_{#theta_{#pi^{0}}}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hPi0OpenAngleRMS[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hPi0OpenAngleRMS[i]);

        if (!isMerged){
            histoName               = "hEtaFrac";
            if(i==0) vecHistosName.push_back(histoName);
            hEtaFrac[i]             = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hEtaFrac; Run Number ; N^{#eta}/N^{Evt}",hNBin,hFBin,hLBin);
            EditTH1(globalRuns, doEquidistantXaxis, hEtaFrac[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
            vecHistos[i].push_back(hEtaFrac[i]);

            histoName               = "hEtaMass";
            if(i==0) vecHistosName.push_back(histoName);
            hEtaMass[i]             = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hEtaMass; Run Number ; m_{#eta}",hNBin,hFBin,hLBin);
            EditTH1(globalRuns, doEquidistantXaxis, hEtaMass[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
            vecHistos[i].push_back(hEtaMass[i]);

            histoName               = "hEtaWidth";
            if(i==0) vecHistosName.push_back(histoName);
            hEtaWidth[i]            = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hEtaWidth; Run Number ; #sigma_{#eta}",hNBin,hFBin,hLBin);
            EditTH1(globalRuns, doEquidistantXaxis, hEtaWidth[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
            vecHistos[i].push_back(hEtaWidth[i]);

            histoName               = "hEtaPt";
            if(i==0) vecHistosName.push_back(histoName);
            hEtaPt[i]               = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hEtaPt; Run Number ; Mean #it{p}_{T}^{#eta}",hNBin,hFBin,hLBin);
            EditTH1(globalRuns, doEquidistantXaxis, hEtaPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
            vecHistos[i].push_back(hEtaPt[i]);

            histoName               = "hEtaPtRMS";
            if(i==0) vecHistosName.push_back(histoName);
            hEtaPtRMS[i]            = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hEtaPtRMS; Run Number ; #sigma_{#it{p}_{T}^{#eta}}",hNBin,hFBin,hLBin);
            EditTH1(globalRuns, doEquidistantXaxis, hEtaPtRMS[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
            vecHistos[i].push_back(hEtaPtRMS[i]);

            histoName               = "hEtaAlpha";
            if(i==0) vecHistosName.push_back(histoName);
            hEtaAlpha[i]            = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hEtaAlpha; Run Number ; Mean #alpha_{#eta}",hNBin,hFBin,hLBin);
            EditTH1(globalRuns, doEquidistantXaxis, hEtaAlpha[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
            vecHistos[i].push_back(hEtaAlpha[i]);

            histoName               = "hEtaAlphaRMS";
            if(i==0) vecHistosName.push_back(histoName);
            hEtaAlphaRMS[i]         = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hEtaAlphaRMS; Run Number ; #sigma_{#alpha_{#eta}}",hNBin,hFBin,hLBin);
            EditTH1(globalRuns, doEquidistantXaxis, hEtaAlphaRMS[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
            vecHistos[i].push_back(hEtaAlphaRMS[i]);

            histoName               = "hEtaY";
            if(i==0) vecHistosName.push_back(histoName);
            hEtaY[i]                = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hEtaY; Run Number ; Mean Y_{#eta}",hNBin,hFBin,hLBin);
            EditTH1(globalRuns, doEquidistantXaxis, hEtaY[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
            vecHistos[i].push_back(hEtaY[i]);

            histoName               = "hEtaYRMS";
            if(i==0) vecHistosName.push_back(histoName);
            hEtaYRMS[i]             = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hEtaYRMS; Run Number ; #sigma_{Y_{#eta}}",hNBin,hFBin,hLBin);
            EditTH1(globalRuns, doEquidistantXaxis, hEtaYRMS[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
            vecHistos[i].push_back(hEtaYRMS[i]);

            histoName               = "hEtaOpenAngle";
            if(i==0) vecHistosName.push_back(histoName);
            hEtaOpenAngle[i]        = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hEtaOpenAngle; Run Number ; Mean #theta_{#eta}",hNBin,hFBin,hLBin);
            EditTH1(globalRuns, doEquidistantXaxis, hEtaOpenAngle[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
            vecHistos[i].push_back(hEtaOpenAngle[i]);

            histoName               = "hEtaOpenAngleRMS";
            if(i==0) vecHistosName.push_back(histoName);
            hEtaOpenAngleRMS[i]     = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hEtaOpenAngleRMS; Run Number ; #sigma_{#theta_{#eta}}",hNBin,hFBin,hLBin);
            EditTH1(globalRuns, doEquidistantXaxis, hEtaOpenAngleRMS[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
            vecHistos[i].push_back(hEtaOpenAngleRMS[i]);
        }
    } // end of loop over datasets

    //*****************************************************************************************************
    //******************************* create log files *****************************************************
    //*****************************************************************************************************
    fstream* fLog           = new fstream[nSets];
    fstream* fEventLog      = new fstream[nSets];
    for(Int_t iStr=0; iStr<nSets; iStr++){
        if(useDataRunListForMC && iStr>=nData) fLog[iStr].open(Form("%s/A-%s-%s.log",outputDir.Data(),DataSets[iStr].Data(),DataSets[0].Data()), ios::out);
        else fLog[iStr].open(Form("%s/A-%s.log",outputDir.Data(),DataSets[iStr].Data()), ios::out);
        fLog[iStr] << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        fLog[iStr] << DataSets[iStr].Data() << endl;
        fLog[iStr] << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        fLog[iStr] << fCollisionSystem.Data() << endl;
        fLog[iStr] << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        fLog[iStr] << "processed cut: " << fCutSelection.Data() << endl;
        fLog[iStr] << calo.Data() << endl;
        fLog[iStr] << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

        if(useDataRunListForMC && iStr>=nData) fEventLog[iStr].open(Form("%s/A-NEvents-%s-%s.log",outputDir.Data(),DataSets[iStr].Data(),DataSets[0].Data()), ios::out);
        else fEventLog[iStr].open(Form("%s/A-NEvents-%s.log",outputDir.Data(),DataSets[iStr].Data()), ios::out);
        fEventLog[iStr] << Form("Run\t%s",DataSets[iStr].Data()) << endl;
    }

    //*****************************************************************************************************
    //********************************** Supporting variables *********************************************
    //*****************************************************************************************************
    MesonFit fitter;

    TString fRootFile;
    TString fDataSet;
    TString fRunNumber;
    Int_t bin                               = -1;

    Double_t** fitValues                    = new Double_t*[nSets];
    Double_t* sumEvents                     = new Double_t[nSets];
    Double_t* sumEventsAll                  = new Double_t[nSets];

    // bad QA runs
    std::vector<TString> vecRunsBad;
    TString fileRunsBad                     = "";

    std::vector<TString>* vecMissingRuns    = new std::vector<TString>[nSets];

    Bool_t isNullNSigmas = kFALSE;
    if(nSigmasBadRun == NULL){
      isNullNSigmas = kTRUE;
      const Int_t nQuantities = 8;
      nSigmasBadRun = (Int_t*)calloc(nQuantities,sizeof(Int_t));
      for(Int_t i=0; i<nQuantities; i++) nSigmasBadRun[i] = 2;
    }


    //*****************************************************************************************************
    //****************************** Looping over DataSets ************************************************
    //*****************************************************************************************************
    for(Int_t i=0; i<nSets; i++) {
        fitValues[i]        = new Double_t[(Int_t)vecHistos[i].size()];
        vecRuns.clear();
	vecRunsBad.clear();
	badRunCalc sVertexZMean       = {nSigmasBadRun[0],0,0.,0.,0.,hVertexZMean[i]};        // calculate mean for every dataset separately
	badRunCalc sCentralityMean    = {nSigmasBadRun[1],0,0.,0.,0.,hCentralityMean[i]};
	badRunCalc sFracWOVtx         = {nSigmasBadRun[2],0,0.,0.,0.,hFracWOVtx[i]};
	badRunCalc sTracksMeanGood    = {nSigmasBadRun[3],0,0.,0.,0.,hTracksMeanGood[i]};
	badRunCalc sConvNCandidatesQA = {nSigmasBadRun[4],0,0.,0.,0.,hConvNCandidatesQA[i]};
	badRunCalc sPi0Frac           = {nSigmasBadRun[5],0,0.,0.,0.,hPi0Frac[i]};
	badRunCalc sPi0Mass           = {nSigmasBadRun[6],0,0.,0.,0.,hPi0Mass[i]};
	badRunCalc sPi0Width          = {nSigmasBadRun[7],0,0.,0.,0.,hPi0Width[i]};
        fDataSet            = vecDataSet.at(i);
        fileRuns            = Form("%s/runNumbers%s.txt", folderRunlists.Data(), fDataSet.Data());
        fileRunsBad         = Form("%s/runNumbers%sBadQA.txt", folderRunlists.Data(), fDataSet.Data());
        if(useDataRunListForMC && i>=nData) {
            fileRuns        = Form("%s/runNumbers%s-%s.txt", folderRunlists.Data(), vecDataSet.at(i).Data(),vecDataSet.at(0).Data());
            cout << "Switch useDataRunListForMC is true, reading runs from: " << fileRuns.Data() << endl;
        }
        if(!readin(fileRuns, vecRuns)) {cout << "ERROR, no Run Numbers could be found! Returning..." << endl; return;}

        //*****************************************************************************************************
        //****************************** Looping over Runs ****************************************************
        //*****************************************************************************************************
        cout << endl;
        cout << "\t----------------------------------------------------------------------------" << endl;
        cout << "\tLooping over Runs of DataSet |" << (vecDataSet.at(i)).Data() << "|" << endl;
        cout << "\t----------------------------------------------------------------------------" << endl;
        cout << endl;
        fLog[i] << "Looping over Runs:" << endl;
        sumEvents[i] = 0;
        sumEventsAll[i] = 0;

        vecMissingRuns[i].clear();
        for(Int_t j=0; j<(Int_t) vecRuns.size(); j++){
            //--------------------------------------------------------------------------------------------------------
            //------------------------- Read in individual files -----------------------------------------------------
            //--------------------------------------------------------------------------------------------------------
            fRunNumber                  = vecRuns.at(j);
            fRootFile                   = Form("%s/%s/%s/%s", filePath.Data(), fDataSet.Data(), fRunNumber.Data(), fileName.Data());
            TFile* RootFile             = new TFile(fRootFile.Data(),"READ");
            if(RootFile->IsZombie()) {vecMissingRuns[i].push_back(fRunNumber); cout << "INFO: ROOT file '" << fRootFile.Data() << "' could not be openend, continue!" << endl; continue;}

            cout << endl;
            cout << "\t\t----------------------------------------------------------------------------" << endl;
            cout << Form("\t\tRun %s", fRunNumber.Data()) << endl;
            cout << "\t\tProcessing file: " << fRootFile.Data() << endl;
            cout << "\t\t----------------------------------------------------------------------------" << endl;
            cout << endl;

            TList* TopDir               = (TList*) RootFile->Get(nameMainDir.Data());
                if(TopDir == NULL) {cout << "ERROR: TopDir not Found"<<endl; return;}
                else TopDir->SetOwner(kTRUE);
            TList* TopContainer         = (TList*) TopDir->FindObject(Form("Cut Number %s",fCutSelection.Data()));
                if(TopContainer == NULL) {cout << "ERROR: " << Form("Cut Number %s",fCutSelection.Data()) << " not found in File" << endl; return;}
                else TopContainer->SetOwner(kTRUE);
            TList* ESDContainer         = (TList*) TopContainer->FindObject(Form("%s ESD histograms",fCutSelection.Data()));
                if(ESDContainer == NULL) {cout << "ERROR: " << Form("%s ESD histograms",fCutSelection.Data()) << " not found in File" << endl; return;}
                else ESDContainer->SetOwner(kTRUE);
            TList* CaloCutsContainer    = (TList*) TopContainer->FindObject(Form("CaloCuts_%s",fClusterCutSelection.Data()));
                if(isCalo && CaloCutsContainer == NULL) {cout << "ERROR: " << Form("CaloCuts_%s",fClusterCutSelection.Data()) << " not found in File" << endl; return;}
                else if(CaloCutsContainer) CaloCutsContainer->SetOwner(kTRUE);
            TList* CaloMergedCutsContainer = (TList*) TopContainer->FindObject(Form("CaloCuts_%s",fMClusterCutSelection.Data()));
                if(isMerged && CaloMergedCutsContainer == NULL) {cout << "ERROR: " << Form("CaloCuts_%s",fMClusterCutSelection.Data()) << " not found in File" << endl; return;}
                else if(CaloMergedCutsContainer) CaloMergedCutsContainer->SetOwner(kTRUE);
            TList* ConvCutsContainer    = (TList*) TopContainer->FindObject(Form("ConvCuts_%s",fGammaCutSelection.Data()));
                if(isConv && ConvCutsContainer == NULL) {cout << "ERROR: " << Form("ConvCuts_%s",fGammaCutSelection.Data()) << " not found in File" << endl; return;}
                else if(ConvCutsContainer) ConvCutsContainer->SetOwner(kTRUE);
            TList* ConvEventCutsContainer    = (TList*) TopContainer->FindObject(Form("ConvEventCuts_%s",fEventCutSelection.Data()));
                if(ConvEventCutsContainer == NULL) {cout << "ERROR: " << Form("ConvEventCuts_%s",fEventCutSelection.Data()) << " not found in File" << endl; return;}
                else ConvEventCutsContainer->SetOwner(kTRUE);
            //--------------------------------------------------------------------------------------------------------
            if(doEquidistantXaxis) bin  = mapBin[fRunNumber];
            else bin = fRunNumber.Atoi() - hFBin;

            //--------------------------------------------------------------------------------------------------------
            //----------------------------- Event quality histograms processing --------------------------------------
            //--------------------------------------------------------------------------------------------------------
            TH1D* EVENTS    = NULL;
            EVENTS          = (TH1D*) ESDContainer->FindObject("NEventsWOWeight");
            if (EVENTS){
                cout << "INFO: Output contains event weights" << endl;
            } else {
                EVENTS      = (TH1D*) ESDContainer->FindObject("NEvents");
            }

            if(!EVENTS){
                cout << "ERROR: Object |NEvents| could not be found! Returning!" << endl;
                return;
            }
            Double_t nEvents;
            if(fEnergyFlag.Contains("PbPb") || fEnergyFlag.Contains("pPb"))
                nEvents = EVENTS->GetBinContent(1);
            else
                nEvents = GetNEvents((TH1*)EVENTS,kFALSE);
            Double_t nEventsAll         = EVENTS->GetEntries() - EVENTS->GetBinContent(4);
            sumEventsAll[i]             += nEventsAll;
            sumEvents[i]                += nEvents;
            //--------------------------------------------------------------------------------------------------------
            Float_t nEventsBin1         = EVENTS->GetBinContent(1);
            if(nEventsBin1==0)
                nEventsBin1             = 1;
            Float_t nEventsBin2         = EVENTS->GetBinContent(2);
            Float_t nEventsBin3         = EVENTS->GetBinContent(3);
            Float_t nEventsBin4         = EVENTS->GetBinContent(4);
            Float_t nEventsBin5         = EVENTS->GetBinContent(5);
            Float_t nEventsBin6         = EVENTS->GetBinContent(6);
            Float_t nEventsBin7         = EVENTS->GetBinContent(7);
            Float_t nEventsBin12 = 0;
            if(EVENTS->GetNbinsX()==12)
                nEventsBin12            = EVENTS->GetBinContent(12);

            Float_t nEventsAllEvt       = nEventsBin1+nEventsBin2+nEventsBin3+nEventsBin4+nEventsBin5+nEventsBin6+nEventsBin7;
            Float_t nEventsAllEvtErr    = sqrt(nEventsAllEvt);
            Float_t nEventsMinBiasEvt   = nEventsBin1+nEventsBin2+nEventsBin5+nEventsBin6+nEventsBin7;
            Float_t nEventsMinBiasEvtErr        = sqrt(nEventsMinBiasEvt);
            Float_t nEventsNormEvt      = nEventsBin1+(nEventsBin1/(nEventsBin1+nEventsBin5))*nEventsBin6;
            Float_t nEventsNormEvtErr   = sqrt(nEventsNormEvt);

            Float_t ratioWOVtxEvt       = nEventsBin6/nEventsMinBiasEvt;
            Float_t ratioWOVtxEvtErr    = sqrt( pow(sqrt(nEventsBin6)/nEventsMinBiasEvt,2)  + pow( nEventsMinBiasEvtErr*nEventsBin6/pow(nEventsMinBiasEvt,2),2) );
            Float_t ratioWVtxOutside10cmEvt     = nEventsBin5/nEventsMinBiasEvt;
            Float_t ratioWVtxOutside10cmEvtErr  = sqrt( pow(sqrt(nEventsBin5)/nEventsMinBiasEvt,2)  + pow( nEventsMinBiasEvtErr*nEventsBin5/pow(nEventsMinBiasEvt,2),2) );
            Float_t ratioPileUpEvt      = nEventsBin7/nEventsMinBiasEvt;
            Float_t ratioPileUpEvtErr   = sqrt( pow(sqrt(nEventsBin7)/nEventsMinBiasEvt,2)  + pow( nEventsMinBiasEvtErr*nEventsBin7/pow(nEventsMinBiasEvt,2),2) );
            Float_t ratioSPDClusTrackEvt        = nEventsBin12/nEventsMinBiasEvt;
            Float_t ratioSPDClusTrackEvtErr     = sqrt( pow(sqrt(nEventsBin12)/nEventsMinBiasEvt,2)  + pow( nEventsMinBiasEvtErr*nEventsBin12/pow(nEventsMinBiasEvt,2),2) );

            Float_t ratioGoodEventsEvt  = nEventsBin1/nEventsMinBiasEvt;
            Float_t ratioGoodEventsEvtErr       = sqrt( pow(sqrt(nEventsBin1)/nEventsMinBiasEvt,2)  + pow( nEventsMinBiasEvtErr*nEventsBin1/pow(nEventsMinBiasEvt,2),2) );
            Float_t ratioNormAllEvt     = nEventsNormEvt/nEventsAllEvt;
            Float_t ratioNormAllEvtErr  = sqrt( pow(sqrt(nEventsNormEvt)/nEventsAllEvt,2)  + pow( nEventsAllEvtErr*nEventsNormEvt/pow(nEventsAllEvt,2),2) );
            Float_t ratioNormMinBiasEvt = nEventsNormEvt/nEventsMinBiasEvt;
            Float_t ratioNormMinBiasEvtErr      = sqrt( pow(sqrt(nEventsNormEvt)/nEventsMinBiasEvt,2)  + pow( nEventsMinBiasEvtErr*nEventsNormEvt/pow(nEventsMinBiasEvt,2),2) );
            //--------------------------------------------------------------------------------------------------------
            // write to log file
            fLog[i] << "----------------------------------------------------------------------------" << endl;
            fLog[i] << "Processing file: " << fRootFile.Data() << endl;
            fLog[i] << "Run: " << fRunNumber.Data() << endl;
            fLog[i] << "NEventsAll: '" << nEventsAll << "', NEvents (for normalization): '" << nEvents << "'" << endl;
            fLog[i] << "----------------------------------------------------------------------------" << endl;
            fEventLog[i] << Form("%s\t%0.0f",fRunNumber.Data(),nEventsAll) << endl;
            //--------------------------------------------------------------------------------------------------------
            if( nEvents < 1 ){cout << "Warning: number of accepted events in run: " << nEvents << "! Setting nEvents to 1..." << endl; nEvents = 1;}

            //--------------------------------------------------------------------------------------------------------
            // fill trending histos for events
            hNEvents[i]->SetBinContent(bin, nEvents);
            hNEventsMinBias[i]->SetBinContent(bin, nEventsMinBiasEvt);
            hNEventsAll[i]->SetBinContent(bin, nEventsAllEvt);
            if(nEventsBin1==1 && ratioGoodEventsEvt==1){
                hNEventsFracGoodEvents[i]->SetBinContent(bin,0);
                hNEventsFracGoodEvents[i]->SetBinError(bin,0);
            }else{
                hNEventsFracGoodEvents[i]->SetBinContent(bin,ratioGoodEventsEvt);
                hNEventsFracGoodEvents[i]->SetBinError(bin,ratioGoodEventsEvtErr);
            }
            hNEventsFracNormAll[i]->SetBinContent(bin,ratioNormAllEvt);
            hNEventsFracNormAll[i]->SetBinError(bin,ratioNormAllEvtErr);
            hNEventsFracMinBias[i]->SetBinContent(bin,ratioNormMinBiasEvt);
            hNEventsFracMinBias[i]->SetBinError(bin,ratioNormMinBiasEvtErr);

            hFracWVtxOutside10cm[i]->SetBinContent(bin,ratioWVtxOutside10cmEvt);
            hFracWVtxOutside10cm[i]->SetBinError(bin,ratioWVtxOutside10cmEvtErr);
            hFracWOVtx[i]->SetBinContent(bin,ratioWOVtxEvt);
            hFracWOVtx[i]->SetBinError(bin,ratioWOVtxEvtErr);
	    sFracWOVtx.mean += ratioWOVtxEvt;
            hFracPileUp[i]->SetBinContent(bin,ratioPileUpEvt);
            hFracPileUp[i]->SetBinError(bin,ratioPileUpEvtErr);
            hFracSPDClusTrack[i]->SetBinContent(bin,ratioSPDClusTrackEvt);
            hFracSPDClusTrack[i]->SetBinError(bin,ratioSPDClusTrackEvtErr);

            //--------------------------------------------------------------------------------------------------------
            //----------------------- number of reference tracks in the TPC (|eta| < 0.8) ----------------------------
            //--------------------------------------------------------------------------------------------------------
            TH1F* GOODESD               = (TH1F*) ESDContainer->FindObject("GoodESDTracks");
            if(GOODESD){
                hTracksMeanGood[i]->SetBinContent(bin, GOODESD->GetMean());
                hTracksMeanGood[i]->SetBinError(bin, GOODESD->GetMeanError());
		sTracksMeanGood.mean += GOODESD->GetMean();
                hTracksRMSGood[i]->SetBinContent(bin, GOODESD->GetRMS());
                hTracksRMSGood[i]->SetBinError(bin, GOODESD->GetRMSError());
            }else cout << "INFO: Object |GoodESDTracks| could not be found! Skipping Fill..." << endl;

            //--------------------------------------------------------------------------------------------------------
            //------------------------- primary vertex distribution in Z ---------------------------------------------
            //--------------------------------------------------------------------------------------------------------
            TH1D* ZVertex               = (TH1D*) ESDContainer->FindObject("VertexZ");
            if(ZVertex){
                hVertexZMean[i]->SetBinContent(bin, ZVertex->GetMean());
                hVertexZMean[i]->SetBinError(bin, ZVertex->GetMeanError());
                hVertexZRMS[i]->SetBinContent(bin, ZVertex->GetRMS());
                hVertexZRMS[i]->SetBinError(bin, ZVertex->GetRMSError());
		sVertexZMean.mean +=  ZVertex->GetMean();
                TH1D* tempVertexZ       = new TH1D(*ZVertex);
                tempVertexZ->GetXaxis()->SetTitle("z-Vertex (cm)");
                tempVertexZ->GetYaxis()->SetTitle("#frac{1}{N_{Events}} #frac{dN}{dz}");
                tempVertexZ->SetTitle(Form("%s",fRunNumber.Data()));
                tempVertexZ->Sumw2();
                tempVertexZ->Scale(1 / tempVertexZ->GetEntries());
                tempVertexZ->SetDirectory(0);
                vecVertexZ[i].push_back(tempVertexZ);
                if(i==0) mapVertexRatio[j]=bin;
                else if( (mapVertexRatio.find(j) != mapVertexRatio.end()) && mapVertexRatio[j]>=0 && mapVertexRatio[j]<(Int_t)vecVertexZ[0].size() ){
                    TH1D* tempVertexRatio   = new TH1D(*tempVertexZ);
                    tempVertexRatio->Divide(tempVertexRatio,vecVertexZ[0].at(mapVertexRatio[j]),1,1,"B");
                    tempVertexRatio->SetName(Form("%s_z-VertexRatio_%s",fRunNumber.Data(),fDataSet.Data()));
                    tempVertexRatio->GetXaxis()->SetRangeUser(-10,10);
                    tempVertexRatio->GetYaxis()->SetTitle("ratio");
                    tempVertexRatio->GetYaxis()->SetRangeUser(0,2);
                    vecVertexZRatio[i-1].push_back(tempVertexRatio);
                }
            }else cout << "INFO: Object |VertexZ| could not be found! Skipping Fill..." << endl;

            if(fEnergyFlag.Contains("PbPb")){
            //--------------------------------------------------------------------------------------------------------
            //--------------------------------------- centrality -----------------------------------------------------
            //--------------------------------------------------------------------------------------------------------
            TH1D* Centrality               = (TH1D*) ESDContainer->FindObject("Centrality");
            if(Centrality){
                hCentralityMean[i]->SetBinContent(bin, Centrality->GetMean());
                hCentralityMean[i]->SetBinError(bin, Centrality->GetMeanError());
		sCentralityMean.mean += Centrality->GetMean();
            }else cout << "INFO: Object |Centrality| could not be found! Skipping Fill..." << endl;

            //--------------------------------------------------------------------------------------------------------
            //--------------------------------------- event plane angle ----------------------------------------------
            //--------------------------------------------------------------------------------------------------------
            TH1D* EventPlaneAngle          = (TH1D*) ConvEventCutsContainer->FindObject(Form("EventPlaneAngle %s",fEventCutSelection.Data()));
            if(EventPlaneAngle){
                hEventPlaneAngleMean[i]->SetBinContent(bin, EventPlaneAngle->GetMean());
                hEventPlaneAngleMean[i]->SetBinError(bin, EventPlaneAngle->GetMeanError());
            } else cout << "INFO: Object |EventPlaneAngle| could not be found! Skipping Fill..." << endl;
            }

            //--------------------------------------------------------------------------------------------------------
            //------------------------- Calorimeter selection histograms ---------------------------------------------
            //--------------------------------------------------------------------------------------------------------
            TH1F* CaloACC;
            TH1F* CaloIPS;
            if(isCalo){
                CaloACC                 = (TH1F*) CaloCutsContainer->FindObject(Form("AcceptanceCuts %s", fClusterCutSelection.Data()));
                CaloIPS                 = (TH1F*) CaloCutsContainer->FindObject(Form("IsPhotonSelected %s", fClusterCutSelection.Data()));
                if(CaloACC && CaloIPS){
                    Double_t CaloNClusters      = CaloACC->GetBinContent(1);
                    Double_t CaloNClustersQA    = CaloIPS->GetBinContent(5);
                    hCaloNClusters[i]->SetBinContent(bin, CaloNClusters / nEvents);
                    hCaloNClusters[i]->SetBinError(bin, sqrt(CaloNClusters) / nEvents);
                    hCaloNClustersQA[i]->SetBinContent(bin, CaloNClustersQA / nEvents);
                    hCaloNClustersQA[i]->SetBinError(bin, sqrt(CaloNClustersQA) / nEvents);
                }else cout << "INFO: Object |AcceptanceCuts| or |IsPhotonSelected| could not be found! Skipping Fill..." << endl;
            }

            //--------------------------------------------------------------------------------------------------------
            //------------------------- Calorimeter selection histograms ---------------------------------------------
            //--------------------------------------------------------------------------------------------------------
            TH1F* CaloMergedACC;
            TH1F* CaloMergedIPS;
            if(isMerged){
                CaloMergedACC           = (TH1F*) CaloMergedCutsContainer->FindObject(Form("AcceptanceCuts %s", fMClusterCutSelection.Data()));
                CaloMergedIPS           = (TH1F*) CaloMergedCutsContainer->FindObject(Form("IsPhotonSelected %s", fMClusterCutSelection.Data()));
                if(CaloMergedACC && CaloMergedIPS){
                    Double_t CaloNClusters      = CaloMergedACC->GetBinContent(1);
                    Double_t CaloNClustersQA    = CaloMergedIPS->GetBinContent(5);
                    hCaloMergedNClusters[i]->SetBinContent(bin, CaloNClusters / nEvents);
                    hCaloMergedNClusters[i]->SetBinError(bin, sqrt(CaloNClusters) / nEvents);
                    hCaloMergedNClustersQA[i]->SetBinContent(bin, CaloNClustersQA / nEvents);
                    hCaloMergedNClustersQA[i]->SetBinError(bin, sqrt(CaloNClustersQA) / nEvents);
                }else cout << "INFO: Object |AcceptanceCuts| or |IsPhotonSelected| could not be found! Skipping Fill..." << endl;
            }


            //--------------------------------------------------------------------------------------------------------
            //-------------------------- Conversion selection histograms ---------------------------------------------
            //--------------------------------------------------------------------------------------------------------
            if(isConv){
                TH1F* ConvIPS           = (TH1F*) ConvCutsContainer->FindObject(Form("IsPhotonSelected %s", fGammaCutSelection.Data()));
                if(ConvIPS){
                    Double_t ConvNCandidates    = ConvIPS->GetBinContent(1);
                    Double_t ConvNCandidatesQA  = 0;
                    TString ConvIPSbin9         = ConvIPS->GetXaxis()->GetBinLabel(9);
                    TString ConvIPSbin10        = ConvIPS->GetXaxis()->GetBinLabel(10);
                    if(ConvIPSbin9.CompareTo("out")==0) ConvNCandidatesQA = ConvIPS->GetBinContent(9);
                    else if(ConvIPSbin10.CompareTo("out")==0) ConvNCandidatesQA = ConvIPS->GetBinContent(10);
                    else {cout << "------EventQA_Runwise: Could not determine ConvNCandidatesQA, setting to zero------" << endl;}
                    hConvNCandidates[i]->SetBinContent(bin, ConvNCandidates / nEvents);
                    hConvNCandidates[i]->SetBinError(bin, sqrt(ConvNCandidates) / nEvents);
                    hConvNCandidatesQA[i]->SetBinContent(bin, ConvNCandidatesQA / nEvents);
                    hConvNCandidatesQA[i]->SetBinError(bin, sqrt(ConvNCandidatesQA) / nEvents);
		    sConvNCandidatesQA.mean += ConvNCandidatesQA / nEvents;
                }else cout << "INFO: Object |IsPhotonSelected| could not be found! Skipping Fill..." << endl;
            }
            //--------------------------------------------------------------------------------------------------------
            //---------------------------- Neutral meson properties --------------------------------------------------
            //--------------------------------------------------------------------------------------------------------
            if (!isMerged){
                TH2D* ESD_Mother        = (TH2D*) ESDContainer->FindObject("ESD_Mother_InvMass_Pt");
                TH2D* ESD_Background    = (TH2D*) ESDContainer->FindObject("ESD_Background_InvMass_Pt");
                if( ESD_Mother && ESD_Background ){
                    if( nEvents > 10 ){
                        TString outputDirDataSet    = Form("%s/%s",outputDir.Data(), DataSets[i].Data());
                        if(useDataRunListForMC && i>=nData) {
                            outputDirDataSet        = Form("%s/%s-%s", outputDir.Data(), DataSets[i].Data(),DataSets[0].Data());
                        }
                        gSystem->Exec("mkdir -p "+outputDirDataSet);
                        gSystem->Exec("mkdir -p "+outputDirDataSet+"/ExtQA");

                        Bool_t kScs = fitter.DoFitting(ESD_Mother, ESD_Background, nEventsBin1, fMode, Form("%s/ExtQA",outputDirDataSet.Data()), fRunNumber,kTRUE,kTRUE,fLog[i]);

                        Double_t widthPi    = 0; Double_t widthPiErr    = 0; Double_t massPi    = 0; Double_t massPiErr     = 0;
                        Double_t ratioPi0   = 0; Double_t ratioPi0Err   = 0;
                        Double_t widthEta   = 0; Double_t widthEtaErr   = 0; Double_t massEta   = 0; Double_t massEtaErr    = 0;
                        Double_t ratioEta   = 0; Double_t ratioEtaErr   = 0;

                        fitter.GetMeson(kTRUE,widthPi,widthPiErr,massPi,massPiErr);
                        fitter.GetMesonRatios(kTRUE,ratioPi0,ratioPi0Err);
                        fitter.GetMeson(kFALSE,widthEta,widthEtaErr,massEta,massEtaErr);
                        fitter.GetMesonRatios(kFALSE,ratioEta,ratioEtaErr);

                        hPi0Frac[i]->SetBinContent(bin,ratioPi0);
                        hPi0Frac[i]->SetBinError(bin,ratioPi0Err);
			sPi0Frac.mean += ratioPi0;
                        hPi0Mass[i]->SetBinContent(bin,massPi);
                        hPi0Mass[i]->SetBinError(bin,massPiErr);
			sPi0Mass.mean += massPi;
                        hPi0Width[i]->SetBinContent(bin,widthPi/2.35);
                        hPi0Width[i]->SetBinError(bin,widthPiErr/2.35);
			sPi0Width.mean += widthPi/2.35;
                        hEtaFrac[i]->SetBinContent(bin,ratioEta);
                        hEtaFrac[i]->SetBinError(bin,ratioEtaErr);
                        hEtaMass[i]->SetBinContent(bin,massEta);
                        hEtaMass[i]->SetBinError(bin,massEtaErr);
                        hEtaWidth[i]->SetBinContent(bin,widthEta/2.35);
                        hEtaWidth[i]->SetBinError(bin,widthEtaErr/2.35);

                        if(kScs){
                            TH1D* tempPi0   = new TH1D(*fitter.GetSignalPi0());
                            tempPi0->GetXaxis()->SetRangeUser(0.02,0.22);
                            tempPi0->GetXaxis()->SetTitle("M_{#gamma#gamma} (GeV/#it{c}^2)");
                            tempPi0->GetYaxis()->SetTitle("#frac{1}{N_{Events}} #frac{dN}{dM_{#gamma#gamma}}");
                            tempPi0->SetTitle(Form("%s, nEvents: %.2e",fRunNumber.Data(),nEvents));
                            tempPi0->Sumw2();
                            tempPi0->Scale(1 / nEvents);
                            tempPi0->SetDirectory(0);
                            vecPi0Signal[i].push_back(tempPi0);

                            TH1D* tempEta   = new TH1D(*fitter.GetSignalEta());
                            tempEta->GetXaxis()->SetRangeUser(0.45,0.65);
                            tempEta->GetXaxis()->SetTitle("M_{#gamma#gamma} (GeV/#it{c}^2)");
                            tempEta->GetYaxis()->SetTitle("#frac{1}{N_{Events}} #frac{dN}{dM_{#gamma#gamma}}");
                            tempEta->SetTitle(Form("%s, nEvents: %.2e",fRunNumber.Data(),nEvents));
                            tempEta->Sumw2();
                            tempEta->Scale(1 / nEvents);
                            tempEta->SetDirectory(0);
                            vecEtaSignal[i].push_back(tempEta);

                            vecPi0Fit[i].push_back(new TF1(*fitter.GetFitPi0()));
                            vecEtaFit[i].push_back(new TF1(*fitter.GetFitEta()));
                        }
                    }else{
                        cout << "INFO: nEvents<10, skipping mass fits...!" << endl;
                        hPi0Frac[i]->SetBinContent(bin,0.);
                        hPi0Frac[i]->SetBinError(bin,0.);
                        hPi0Mass[i]->SetBinContent(bin,0.);
                        hPi0Mass[i]->SetBinError(bin,0.);
                        hPi0Width[i]->SetBinContent(bin,0.);
                        hPi0Width[i]->SetBinError(bin,0.);

                        hEtaFrac[i]->SetBinContent(bin,0.);
                        hEtaFrac[i]->SetBinError(bin,0.);
                        hEtaMass[i]->SetBinContent(bin,0.);
                        hEtaMass[i]->SetBinError(bin,0.);
                        hEtaWidth[i]->SetBinContent(bin,0.);
                        hEtaWidth[i]->SetBinError(bin,0.);
                    }
                }else{ cout << "ERROR: Object |ESD_Mother_InvMass_Pt| or |ESD_BG_InvMass| could not be found! Returning..." << endl; return;}
            }
            //--------------------------------------------------------------------------------------------------------
            //--------------------------- Pion properties ------------------------------------------------------------
            //--------------------------------------------------------------------------------------------------------
            // define name for meson property histograms
            TString namePi0MesonY       = "ESD_MotherPi0_Pt_Y";
            TString namePi0MesonAlpha   = "ESD_MotherPi0_Pt_Alpha";
            TString namePi0MesonOpen    = "ESD_MotherPi0_Pt_OpenAngle";
            if ( isMerged ){
                namePi0MesonY           = "ESD_Mother_Pt_Y";
                namePi0MesonAlpha       = "ESD_Mother_Pt_Alpha";
                namePi0MesonOpen        = "ESD_Mother_Pt_OpenAngle";
            }

            TH2D* Pi0PtY                = (TH2D*) ESDContainer->FindObject(namePi0MesonY.Data());
            if(Pi0PtY){
                TH1D* Pi0Pt             = (TH1D*) Pi0PtY->ProjectionX("Pi0_Pt");
                TH1D* Pi0Y              = (TH1D*) Pi0PtY->ProjectionY("Pi0_Y");
                hPi0Pt[i]->SetBinContent(bin, Pi0Pt->GetMean());
                hPi0Pt[i]->SetBinError(bin, Pi0Pt->GetMeanError());
                hPi0PtRMS[i]->SetBinContent(bin, Pi0Pt->GetRMS());
                hPi0PtRMS[i]->SetBinError(bin, Pi0Pt->GetRMSError());
                hPi0Y[i]->SetBinContent(bin, Pi0Y->GetMean());
                hPi0Y[i]->SetBinError(bin, Pi0Y->GetMeanError());
                hPi0YRMS[i]->SetBinContent(bin, Pi0Y->GetRMS());
                hPi0YRMS[i]->SetBinError(bin, Pi0Y->GetRMSError());
                delete Pi0Y;
                delete Pi0Pt;
            }else cout << "INFO: Object |ESD_MotherPi0_Pt_Y| could not be found! Skipping Fill..." << endl;

            TH2D* Pi0PtAlpha            = (TH2D*) ESDContainer->FindObject(namePi0MesonAlpha.Data());
            if(Pi0PtAlpha){
                TH1D* Pi0Alpha          = (TH1D*) Pi0PtAlpha->ProjectionY("Pi0_Alpha");
                hPi0Alpha[i]->SetBinContent(bin, Pi0Alpha->GetMean());
                hPi0Alpha[i]->SetBinError(bin, Pi0Alpha->GetMeanError());
                hPi0AlphaRMS[i]->SetBinContent(bin, Pi0Alpha->GetRMS());
                hPi0AlphaRMS[i]->SetBinError(bin, Pi0Alpha->GetRMSError());
                delete Pi0Alpha;
            }else cout << "INFO: Object |ESD_MotherPi0_Pt_Alpha| could not be found! Skipping Fill..." << endl;

            TH2D* Pi0PtOpenAngle        = (TH2D*) ESDContainer->FindObject(namePi0MesonOpen.Data());
            if(Pi0PtOpenAngle){
                TH1D* Pi0OpenAngle      = (TH1D*) Pi0PtOpenAngle->ProjectionY("Pi0_OpenAngle");
                hPi0OpenAngle[i]->SetBinContent(bin, Pi0OpenAngle->GetMean());
                hPi0OpenAngle[i]->SetBinError(bin, Pi0OpenAngle->GetMeanError());
                hPi0OpenAngleRMS[i]->SetBinContent(bin, Pi0OpenAngle->GetRMS());
                hPi0OpenAngleRMS[i]->SetBinError(bin, Pi0OpenAngle->GetRMSError());
                delete Pi0OpenAngle;
            }else cout << "INFO: Object |ESD_MotherPi0_Pt_OpenAngle| could not be found! Skipping Fill..." << endl;


            //--------------------------------------------------------------------------------------------------------
            //---------------------------- Eta properties ------------------------------------------------------------
            //--------------------------------------------------------------------------------------------------------
            if (!isMerged){
                TH2D* EtaPtY            = (TH2D*) ESDContainer->FindObject("ESD_MotherEta_Pt_Y");
                if(EtaPtY){
                    TH1D* EtaPt         = (TH1D*) EtaPtY->ProjectionX("Eta_Pt");
                    TH1D* EtaY          = (TH1D*) EtaPtY->ProjectionY("Eta_Y");
                    hEtaPt[i]->SetBinContent(bin, EtaPt->GetMean());
                    hEtaPt[i]->SetBinError(bin, EtaPt->GetMeanError());
                    hEtaPtRMS[i]->SetBinContent(bin, EtaPt->GetRMS());
                    hEtaPtRMS[i]->SetBinError(bin, EtaPt->GetRMSError());
                    hEtaY[i]->SetBinContent(bin, EtaY->GetMean());
                    hEtaY[i]->SetBinError(bin, EtaY->GetMeanError());
                    hEtaYRMS[i]->SetBinContent(bin, EtaY->GetRMS());
                    hEtaYRMS[i]->SetBinError(bin, EtaY->GetRMSError());
                    delete EtaY;
                    delete EtaPt;
                }else cout << "INFO: Object |ESD_MotherEta_Pt_Y| could not be found! Skipping Fill..." << endl;

                TH2D* EtaPtAlpha        = (TH2D*) ESDContainer->FindObject("ESD_MotherEta_Pt_Alpha");
                if(EtaPtAlpha){
                    TH1D* EtaAlpha      = (TH1D*) EtaPtAlpha->ProjectionY("Eta_Alpha");
                    hEtaAlpha[i]->SetBinContent(bin, EtaAlpha->GetMean());
                    hEtaAlpha[i]->SetBinError(bin, EtaAlpha->GetMeanError());
                    hEtaAlphaRMS[i]->SetBinContent(bin, EtaAlpha->GetRMS());
                    hEtaAlphaRMS[i]->SetBinError(bin, EtaAlpha->GetRMSError());
                    delete EtaAlpha;
                }else cout << "INFO: Object |ESD_MotherEta_Pt_Alpha| could not be found! Skipping Fill..." << endl;

                TH2D* EtaPtOpenAngle    = (TH2D*) ESDContainer->FindObject("ESD_MotherEta_Pt_OpenAngle");
                if(EtaPtOpenAngle){
                    TH1D* EtaOpenAngle  = (TH1D*) EtaPtOpenAngle->ProjectionY("Eta_OpenAngle");
                    hEtaOpenAngle[i]->SetBinContent(bin, EtaOpenAngle->GetMean());
                    hEtaOpenAngle[i]->SetBinError(bin, EtaOpenAngle->GetMeanError());
                    hEtaOpenAngleRMS[i]->SetBinContent(bin, EtaOpenAngle->GetRMS());
                    hEtaOpenAngleRMS[i]->SetBinError(bin, EtaOpenAngle->GetRMSError());
                    delete EtaOpenAngle;
                }else cout << "INFO: Object |ESD_MotherEta_Pt_OpenAngle| could not be found! Skipping Fill..." << endl;
            }

            //--------------------------------------------------------------------------------------------------------
            //---------------------------- Cleanup -------------------------------------------------------------------
            //--------------------------------------------------------------------------------------------------------
            delete TopDir;

            RootFile->Close();
            delete RootFile;
        } // end of loop over runs

        //--------------------------------------------------------------------------------------------------------
        //--------------------------------------- Mean values ----------------------------------------------------
        //--------------------------------------------------------------------------------------------------------
	if(doExtQA==3){
	  Bool_t isRunBadFlag = kFALSE;
	  sVertexZMean.mean /= vecRuns.size();  // divide the sum of bin contents by number of runs in the current dataset
	  sCentralityMean.mean /= vecRuns.size();
	  sFracWOVtx.mean /= vecRuns.size();
	  sTracksMeanGood.mean /= vecRuns.size();
	  sConvNCandidatesQA.mean /= vecRuns.size();
	  sPi0Frac.mean /= vecRuns.size();
	  sPi0Mass.mean /= vecRuns.size();
	  sPi0Width.mean /= vecRuns.size();
	  for(Int_t j=0; j<(Int_t) vecRuns.size(); j++){  // loop over runs j of this dataset
	    fRunNumber = vecRuns.at(j);
	    if(doEquidistantXaxis) bin  = mapBin[fRunNumber]; else bin = fRunNumber.Atoi() - hFBin;     // bin: run number in global run list, starting from 1
	    if(isBadRun(&sVertexZMean,bin))       isRunBadFlag = kTRUE;
	    if(isBadRun(&sCentralityMean,bin))    isRunBadFlag = kTRUE;
	    if(isBadRun(&sFracWOVtx,bin))         isRunBadFlag = kTRUE;
	    if(isBadRun(&sTracksMeanGood,bin))    isRunBadFlag = kTRUE;
	    if(isBadRun(&sConvNCandidatesQA,bin)) isRunBadFlag = kTRUE;
	    if(isBadRun(&sPi0Frac,bin))           isRunBadFlag = kTRUE;
	    if(isBadRun(&sPi0Mass,bin))           isRunBadFlag = kTRUE;
	    if(isBadRun(&sPi0Width,bin))          isRunBadFlag = kTRUE;
	    if(isRunBadFlag) vecRunsBad.push_back(fRunNumber);
	  } // end of loop over runs
	  cout << "INFO: hVertexZMean: "       << sVertexZMean.nRunsBad       << " runs of " <<  vecRuns.size() << " runs deviate by more than " << sVertexZMean.nSigmaBad       << " error bars from the mean over all runs" << endl;
	  cout << "INFO: hCentralityMean: "    << sCentralityMean.nRunsBad    << " runs of " <<  vecRuns.size() << " runs deviate by more than " << sCentralityMean.nSigmaBad    << " error bars from the mean over all runs" << endl;
	  cout << "INFO: hFracWOVtx: "         << sFracWOVtx.nRunsBad         << " runs of " <<  vecRuns.size() << " runs deviate by more than " << sFracWOVtx.nSigmaBad         << " error bars from the mean over all runs" << endl;
	  cout << "INFO: hTracksMeanGood: "    << sTracksMeanGood.nRunsBad    << " runs of " <<  vecRuns.size() << " runs deviate by more than " << sTracksMeanGood.nSigmaBad    << " error bars from the mean over all runs" << endl;
	  cout << "INFO: hConvNCandidatesQA: " << sConvNCandidatesQA.nRunsBad << " runs of " <<  vecRuns.size() << " runs deviate by more than " << sConvNCandidatesQA.nSigmaBad << " error bars from the mean over all runs" << endl;
	  cout << "INFO: hPi0Frac: "           << sPi0Frac.nRunsBad           << " runs of " <<  vecRuns.size() << " runs deviate by more than " << sPi0Frac.nSigmaBad           << " error bars from the mean over all runs" << endl;
	  cout << "INFO: hPi0Mass: "           << sPi0Mass.nRunsBad           << " runs of " <<  vecRuns.size() << " runs deviate by more than " << sPi0Mass.nSigmaBad << " error bars from the mean over all runs" << endl;
	  cout << "INFO: hPi0Width: "          << sPi0Width.nRunsBad          << " runs of " <<  vecRuns.size() << " runs deviate by more than " << sPi0Width.nSigmaBad << " error bars from the mean over all runs" << endl;
	  if(!vecRunsBad.empty()) writeout(fileRunsBad, vecRunsBad, kTRUE);
	}

        //--------------------------------------------------------------------------------------------------------
        //---------------------------------------- Fitters -------------------------------------------------------
        //--------------------------------------------------------------------------------------------------------
        TF1 *tfFit;
        for(Int_t iHist=0; iHist<(Int_t)vecHistos[i].size(); iHist++){
            TH1D* temp          = vecHistos[i].at(iHist);
            if(doEquidistantXaxis)
                tfFit           = new TF1("tfFit","[0]",0,temp->GetNbinsX());
            else
                tfFit           = new TF1("tfFit","[0]",rangesRuns[i][0],rangesRuns[i][1]);
            temp->Fit(tfFit,"QRME0");
            fitValues[i][iHist] = tfFit->GetParameter(0);
            delete tfFit;
            tfFit               = 0x0;
        }
    } // end of loop over datasets
    
    if(isNullNSigmas) free(nSigmasBadRun);

    //**************************************************************************************************
    //****************************** Drawing Histograms ************************************************
    //**************************************************************************************************
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "Drawing Histograms" << endl;
    TCanvas* canvas         = new TCanvas("canvas","",10,10,750,500);  // gives the page size
    Double_t leftMar        = 0.09;
    Double_t rightMar       = 0.025;
    Double_t topMargin      = 0.04;
    Double_t bottomMargin   = 0.09;
    DrawGammaCanvasSettings(canvas, leftMar, rightMar, topMargin, bottomMargin);
    Double_t yMax;  // min of y range for vertical lines
    Double_t yMin;  // max of y range for vertical lines
    Bool_t adjustedRange = kFALSE;

    // Plot single periods as well
    if(doHistsForEverySet) {
        TString outputDirDataSet;

        for(Int_t i=0; i<nSets; i++) {
            cout << "DataSet: " << DataSets[i].Data() << endl;
            outputDirDataSet = Form("%s/%s",outputDir.Data(), DataSets[i].Data());
            if(useDataRunListForMC && i>=nData) {
                outputDirDataSet = Form("%s/%s-%s", outputDir.Data(), DataSets[i].Data(),DataSets[0].Data());
                cout << "Switch useDataRunListForMC is true, output to: " << outputDirDataSet.Data() << endl;
            }
            gSystem->Exec("mkdir -p "+outputDirDataSet);
            gSystem->Exec("mkdir -p "+outputDirDataSet+"/ExtQA");

            TString fTrigger = "";
            if(doTrigger && i<nData){
                TString fTriggerCut = fEventCutSelection(3,2);
                fTrigger = AnalyseSpecialTriggerCut(fTriggerCut.Atoi(), plotDataSets[i]);
                cout << "Trigger: '" << fTrigger.Data() << "'" << endl;
                if(fTrigger.Contains("not defined")){
                    fTrigger = "";
                    cout << "INFO: Trigger cut not defined!" << endl;
                }
            }

            if(i>0){
                DrawVectorOverviewTH1D( canvas, vecVertexZRatio[i-1], "hRatioVertexZ", outputDirDataSet, suffix,
                                        0.13, 0.15, 0.1, 0.14, 0.8, 0.9, 0.12, 0.93, 0x0, kFALSE, kTRUE);
            }

            //--------------------------------------------------------------------------------------------------------
            for(Int_t h=0; h<(Int_t) vecHistos[i].size(); h++) {
                ((TH1D*) vecHistos[i].at(h))->SetTitle("");
                if( ((TString)vecHistosName.at(h)).CompareTo("hNEvents")==0 ||
                    ((TString)vecHistosName.at(h)).CompareTo("hNEventsMinBias")==0 ||
                    ((TString)vecHistosName.at(h)).CompareTo("hNEventsAll")==0 ) {
                    AdjustHistRange(((TH1D*) vecHistos[i].at(h)),10,10,kTRUE);
                    ((TH1D*) vecHistos[i].at(h))->Draw("p");
                } else {
                    AdjustHistRange(((TH1D*) vecHistos[i].at(h)),1.1,1.1,kTRUE);
                    ((TH1D*) vecHistos[i].at(h))->Draw("px0e1");
                }
                if(i<nData) DrawFit(((TH1D*) vecHistos[i].at(h)),i,fitValues[i][h],rangesRuns[i],DataSets[i],plotDataSets[i],0.15,0.9,0.03,1);

                if(doTrigger && i<nData){
                    PutProcessLabelAndEnergyOnPlot(xPosLabel, 0.92, 0.03, fCollisionSystem.Data(), plotDataSets[i].Data(), fTrigger.Data());
                    if(((TString)vecHistosName.at(h)).CompareTo("hCaloNClusters")==0 ||
                        ((TString)vecHistosName.at(h)).CompareTo("hCaloNClustersQA")==0 ||
                        ((TString)vecHistosName.at(h)).CompareTo("hCaloMergedNClusters")==0 ||
                        ((TString)vecHistosName.at(h)).CompareTo("hCaloMergedNClustersQA")==0 )
                        PutProcessLabelAndEnergyOnPlot(xPosLabel, 0.81, 0.03, fClusters.Data(), "", "");
                } else{
                    TString temp="";
                    if(((TString)vecHistosName.at(h)).CompareTo("hCaloNClusters")==0 ||
                        ((TString)vecHistosName.at(h)).CompareTo("hCaloNClustersQA")==0 ||
                        ((TString)vecHistosName.at(h)).CompareTo("hCaloMergedNClusters")==0 ||
                        ((TString)vecHistosName.at(h)).CompareTo("hCaloMergedNClustersQA")==0 )
                        temp = fClusters;
                    PutProcessLabelAndEnergyOnPlot(xPosLabel, 0.92, 0.03, fCollisionSystem.Data(), plotDataSets[i].Data(), temp.Data());

                }

                if( ((TString)vecHistosName.at(h)).CompareTo("hNEvents")==0 ||
                    ((TString)vecHistosName.at(h)).CompareTo("hNEventsMinBias")==0 ||
                    ((TString)vecHistosName.at(h)).CompareTo("hNEventsAll")==0 )
                    SaveCanvas(canvas, Form("%s/%s.%s", outputDirDataSet.Data(), vecHistosName.at(h).Data(), suffix.Data()), kFALSE, kTRUE);
                else SaveCanvas(canvas, Form("%s/%s.%s", outputDirDataSet.Data(), vecHistosName.at(h).Data(), suffix.Data()));
            }
        }

        //**********************************************************************************************************
        //****************************** Drawing Runwise Histograms ************************************************
        //**********************************************************************************************************
        cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        cout << "Drawing Runwise Histograms" << endl;

        TCanvas* canvasRunwise;
        TLegend *legendRuns;

        for(Int_t i=0; i<nSets; i++) {
            vecRuns.clear();
            fDataSet                = vecDataSet.at(i);
            fileRuns                = Form("%s/runNumbers%s.txt", folderRunlists.Data(), fDataSet.Data());
            outputDirDataSet        = Form("%s/%s",outputDir.Data(), DataSets[i].Data());

            if(useDataRunListForMC && i>=nData) {
                fileRuns            = Form("%s/runNumbers%s-%s.txt", folderRunlists.Data(), vecDataSet.at(i).Data(),vecDataSet.at(0).Data());
                outputDirDataSet    = Form("%s/%s-%s", outputDir.Data(), DataSets[i].Data(),DataSets[0].Data());
                cout << "Switch useDataRunListForMC is true, reading runs from: " << fileRuns.Data() << endl;
            }
            if(!readin(fileRuns, vecRuns, kFALSE)) {cout << "ERROR, no Run Numbers could be found! Returning..." << endl; return;}

            for(Int_t iRun=0; iRun<(Int_t)vecMissingRuns[i].size(); iRun++){
                vecRuns.erase(std::remove(vecRuns.begin(), vecRuns.end(), vecMissingRuns[i].at(iRun)), vecRuns.end());
            }

            Int_t NColumns          = ((Int_t) vecRuns.size() / 31 ) + 1;

            canvasRunwise           = new TCanvas("canvasRunwise","",10,10,1350+(NColumns*108),900);  // gives the page size
            DrawGammaCanvasSettings(canvasRunwise, 130.5/(1350.+(NColumns*108.)), (40.5+(NColumns*108.))/(1350.+(NColumns*108.)), topMargin, bottomMargin);
            canvasRunwise->cd();

            Double_t addRight       = ((Double_t)NColumns*108)/(1350+((Double_t)NColumns*108));
            legendRuns              = new TLegend(0.98-addRight,0.09,0.98,0.94);
            legendRuns->SetNColumns(NColumns);
            legendRuns->SetFillColor(0);
            legendRuns->SetLineColor(0);
            legendRuns->SetTextSize(0.03);
            legendRuns->SetTextFont(42);

            cout << "DataSet: " << DataSets[i].Data() << endl;
            TString fTrigger        = "";
            if(doTrigger && i<nData){
                TString fTriggerCut = fEventCutSelection(3,2);
                fTrigger            = AnalyseSpecialTriggerCut(fTriggerCut.Atoi(), plotDataSets[i]);
                if(fTrigger.Contains("not defined"))
                    fTrigger        = "";
            }
            //--------------------------------------------------------------------------------------------------------
            DrawVectorRunwiseTH1D(	canvasRunwise, legendRuns, vecVertexZ[i], vecRuns, 2, 1.1, kTRUE, addRight, 0.8, 0.94, 0.03, 0.8, 0.83, 0.03,
                                    doTrigger, fTrigger, (Bool_t)(i<nData), outputDirDataSet, "VertexZ_Runwise", plotDataSets[i], kFALSE,
                                    fCollisionSystem, "", suffix, kFALSE, kFALSE, kFALSE);
            if (!isMerged){
                //--------------------------------------------------------------------------------------------------------
                DrawVectorRunwiseTH1D(	canvasRunwise, legendRuns, vecPi0Signal[i], vecRuns, 1, 1.1, kTRUE, addRight, 0.8, 0.94, 0.03, 0.8, 0.83, 0.03,
                                        doTrigger, fTrigger, (Bool_t)(i<nData), outputDirDataSet, "Pi0Signal_Runwise", plotDataSets[i], kFALSE,
                                        fCollisionSystem, "", suffix, kFALSE, kFALSE, kFALSE, kFALSE);
                //--------------------------------------------------------------------------------------------------------
                DrawVectorRunwiseTH1D(	canvasRunwise, legendRuns, vecEtaSignal[i], vecRuns, 1, 1.1, kTRUE, addRight, 0.8, 0.94, 0.03, 0.8, 0.83, 0.03,
                                        doTrigger, fTrigger, (Bool_t)(i<nData), outputDirDataSet, "EtaSignal_Runwise", plotDataSets[i], kFALSE,
                                        fCollisionSystem, "", suffix, kFALSE, kFALSE, kFALSE, kFALSE);
            }

            delete legendRuns;
            delete canvasRunwise;
        }
        canvas->cd();
    }

    //**********************************************************************************************************
    //****************************** Combined Trending Histograms **********************************************
    //**********************************************************************************************************
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "Drawing Trending Histograms" << endl;
    //--------------------------------------------------------------------------------------------------------
    if(useDataRunListForMC) cout << "WARNING: useDataRunListForMC is true, overwriting histograms for DataSet!" << endl;

    TLegend *legend = new TLegend(0.15,0.95,0.95,0.98);
    legend->SetNColumns(nSets);
    legend->SetFillColor(0);
    legend->SetLineColor(0);
    legend->SetTextSize(0.03);
    legend->SetTextFont(42);
    //--------------------------------------------------------------------------------------------------------
    if(nSets > 1) {
        for(Int_t h=0; h<(Int_t) vecHistos[0].size(); h++) {
            TString fTrigger = "";
            TString fTriggerCut = fEventCutSelection(3,2);
            for(Int_t i=0; i<nSets; i++) {
                if(doTrigger && i<nData && i==0){
                    fTrigger = AnalyseSpecialTriggerCut(fTriggerCut.Atoi(), plotDataSets[i]);
                    if(fTrigger.Contains("not defined")) fTrigger = "";
                }
                legend->AddEntry(((TH1D*) vecHistos[i].at(h)),plotDataSets[i].Data(),"p");
            }
            if( ((TString)vecHistosName.at(h)).CompareTo("hNEvents")==0 ||
                ((TString)vecHistosName.at(h)).CompareTo("hNEventsMinBias")==0 ||
                ((TString)vecHistosName.at(h)).CompareTo("hNEventsAll")==0 )
	      adjustedRange = AdjustHistRange(vecHistos,10,10,h,nSets,kTRUE, &yMin, &yMax);
            else adjustedRange = AdjustHistRange(vecHistos,1.1,1.1,h,nSets,kTRUE, &yMin, &yMax);
            for(Int_t i=nSets-1; i>=0; i--)            {
                TString draw;
                if(h==0) draw = (i==nSets-1)?"p":"p, same";
                else draw = (i==nSets-1)?"px0e1":"px0e1, same";
                ((TH1D*) vecHistos[i].at(h))->SetTitle("");
                ((TH1D*) vecHistos[i].at(h))->Draw(draw.Data());
                if(i<nData) DrawFit(((TH1D*) vecHistos[i].at(h)),i,fitValues[i][h],rangesRuns[i],DataSets[i],plotDataSets[i],0.15,0.9,0.03,1);
            }
            if (drawVerticalLines){
                if(!adjustedRange){
                    canvas->Update();
                    yMax = canvas->GetUymax();
                    yMin = canvas->GetUymin();
                }
                for(Int_t lineBin=0; lineBin<nLines; lineBin++){
                    verticalLines[lineBin] = new TLine(runRanges[lineBin],yMin,runRanges[lineBin],yMax);
                    verticalLines[lineBin]->SetLineWidth(1);
                    verticalLines[lineBin]->SetLineStyle(2);
                    verticalLines[lineBin]->Draw("same");
                }
            }
            legend->Draw();

            if(doTrigger){
                if( CheckForCaloHist((TString)vecHistosName.at(h)) ) PutProcessLabelAndEnergyOnPlot(xPosLabel, 0.92, 0.03, fCollisionSystem.Data(), fTrigger.Data(), fClusters.Data());
                else PutProcessLabelAndEnergyOnPlot(xPosLabel, 0.92, 0.03, fCollisionSystem.Data(), fTrigger.Data(), "");
            }else{
                if( CheckForCaloHist((TString)vecHistosName.at(h)) ) PutProcessLabelAndEnergyOnPlot(xPosLabel, 0.92, 0.03, fCollisionSystem.Data(), fClusters.Data(), "");
                else PutProcessLabelAndEnergyOnPlot(xPosLabel, 0.92, 0.03, fCollisionSystem.Data(), "", "");
            }

            if(canvas->GetTopMargin()!=0.06) canvas->SetTopMargin(0.06);
            if(useDataRunListForMC && !addSubFolder){
                if( ((TString)vecHistosName.at(h)).CompareTo("hNEvents")==0 ||
                    ((TString)vecHistosName.at(h)).CompareTo("hNEventsMinBias")==0 ||
                    ((TString)vecHistosName.at(h)).CompareTo("hNEventsAll")==0 )
                        SaveCanvas(canvas, Form("%s/%s/%s.%s", outputDir.Data(), DataSets[0].Data(),vecHistosName.at(h).Data(),suffix.Data()), kFALSE, kTRUE);
                else SaveCanvas(canvas, Form("%s/%s/%s.%s", outputDir.Data(), DataSets[0].Data(),vecHistosName.at(h).Data(),suffix.Data()));
            } else {
                if( ((TString)vecHistosName.at(h)).CompareTo("hNEvents")==0 ||
                    ((TString)vecHistosName.at(h)).CompareTo("hNEventsMinBias")==0 ||
                    ((TString)vecHistosName.at(h)).CompareTo("hNEventsAll")==0 )
                     SaveCanvas(canvas, Form("%s/%s.%s", outputDir.Data(), vecHistosName.at(h).Data(),suffix.Data()), kFALSE, kTRUE);
                else SaveCanvas(canvas, Form("%s/%s.%s", outputDir.Data(), vecHistosName.at(h).Data(),suffix.Data()));
            }
            legend->Clear();
	    if(drawVerticalLines){
	      for(Int_t lineBin=0; lineBin<nLines; lineBin++){
		delete verticalLines[lineBin];
	      }
	    }
        }
    }

    //**********************************************************************************************************
    //****************************** Combined Ratio Trending Histograms ****************************************
    //**********************************************************************************************************
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "Drawing Ratio Trending Histograms" << endl;

    if(doHistsForEverySet) {
        TString* ratioSets[nSets-nData];
        TH1D* ratio[nSets-nData];
        TString outputDirDataSet;

        if(nSets>1 && nSets>nData) {
            legend->SetNColumns(nData*(nSets-nData));
            Int_t markerStyles[14]={2,4,5,20,21,22,23,24,25,26,27,28,29,30};

            for(Int_t h=0; h<(Int_t) vecHistos[0].size(); h++) {
                for(Int_t i=0; i<nData; i++) {
                    for(Int_t j=0; j<nSets-nData; j++){
                        ratioSets[j] = new TString(Form("%s / %s", plotDataSets[i].Data(), plotDataSets[j+nData].Data()));
                        ratio[j] = new TH1D(Form("%s%i%i",((TH1D*) vecHistos[i].at(h))->GetName(),i,j),
                                            Form("%s%i%i;%s;%s - Ratio: Data / MC",((TH1D*) vecHistos[i].at(h))->GetTitle(),i,j,((TH1D*) vecHistos[i].at(h))->GetXaxis()->GetTitle(),((TH1D*) vecHistos[i].at(h))->GetYaxis()->GetTitle()),
                                            hNBin,hFBin,hLBin);
                        EditTH1(globalRuns, doEquidistantXaxis, ratio[j], markerStyles[j % 14], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
                    }

                    for(Int_t b=1; b<hNBin+1; b++){
                        Double_t hData = ((TH1D*) vecHistos[i].at(h))->GetBinContent(b);
			Double_t hErrData = ((TH1D*) vecHistos[i].at(h))->GetBinError(b);
                        for(Int_t j=0; j<nSets-nData; j++){
                            Double_t hMC = ((TH1D*) vecHistos[j+nData].at(h))->GetBinContent(b);
                            Double_t hErrMC = ((TH1D*) vecHistos[j+nData].at(h))->GetBinError(b);
                            if(hMC!=0) {
                                if(hData/hMC>1.98) ratio[j]->SetBinContent(b,1.98);
                                else if(hData/hMC<0.02) ratio[j]->SetBinContent(b,0.02);
                                else {
				  ratio[j]->SetBinContent(b,hData/hMC);
				  Double_t hErrRatio = hData*sqrt((hErrData/hData)*(hErrData/hData)+(hErrMC/hMC)*(hErrMC/hMC))/hMC;
				  ratio[j]->SetBinError(b,hErrRatio);
				}
                            }else ratio[j]->SetBinContent(b,1.98);
                        }
                    }

                    for(Int_t j=0; j<nSets-nData; j++) {
                        TString draw = (i==0&&j==0)?"p":"p, same";
                        ratio[j]->SetTitle("");
                        ratio[j]->GetYaxis()->SetRangeUser(0,2);
                        ratio[j]->Draw(draw.Data());
                        legend->AddEntry(ratio[j],ratioSets[j]->Data(),"p");
                    }
                    if (drawVerticalLines){
                        canvas->Update();
                        yMax = canvas->GetUymax();
                        yMin = canvas->GetUymin();
                        for(Int_t lineBin=0; lineBin<nLines; lineBin++){
                            verticalLines[lineBin] = new TLine(runRanges[lineBin],yMin,runRanges[lineBin],yMax);
                            verticalLines[lineBin]->SetLineWidth(1);
                            verticalLines[lineBin]->SetLineStyle(2);
                            verticalLines[lineBin]->Draw("same");
                        }
                    }
                    legend->Draw();
                    outputDirDataSet = Form("%s/%s",outputDir.Data(), DataSets[i].Data());
                    gSystem->Exec("mkdir -p "+outputDirDataSet+"/TrendingRatios");

                    if(doTrigger){
                        TString fTrigger    = "";
                        TString fTriggerCut = fEventCutSelection(3,2);
                        fTrigger            = AnalyseSpecialTriggerCut(fTriggerCut.Atoi(), plotDataSets[i]);
                        if(fTrigger.Contains("not defined"))
                            fTrigger        = "";

                        if( CheckForCaloHist((TString)vecHistosName.at(h)) ) PutProcessLabelAndEnergyOnPlot(xPosLabel, 0.92, 0.03, fCollisionSystem.Data(), fTrigger.Data(), fClusters.Data());
                        else PutProcessLabelAndEnergyOnPlot(xPosLabel, 0.92, 0.03, fCollisionSystem.Data(), fTrigger.Data(), "");
                    }else{
                        if( CheckForCaloHist((TString)vecHistosName.at(h)) ) PutProcessLabelAndEnergyOnPlot(xPosLabel, 0.92, 0.03, fCollisionSystem.Data(), fClusters.Data(), "");
                        else PutProcessLabelAndEnergyOnPlot(xPosLabel, 0.92, 0.03, fCollisionSystem.Data(), "", "");
                    }

                    if(canvas->GetTopMargin()!=0.06) canvas->SetTopMargin(0.06);
                    SaveCanvas(canvas, Form("%s/TrendingRatios/%s.%s", outputDirDataSet.Data(),Form("%s",((TH1D*) vecHistos[i].at(h))->GetName()),suffix.Data()));
                    legend->Clear();
                    if(drawVerticalLines){
                        for(Int_t lineBin=0; lineBin<nLines; lineBin++){
                            delete verticalLines[lineBin];
                        }
                    }
                    for(Int_t j=0; j<nSets-nData; j++) {
                        delete ratio[j];
                        delete ratioSets[j];
                    }
                }
            }
        }else cout << "...skipped due to nSets<=1 or nSets==nData!" << endl;

    }
    delete legend;
    delete canvas;

    //**********************************************************************************************************
    //****************************** Create Output ROOT-File ***************************************************
    //**********************************************************************************************************
    const char* nameOutput;

    for(Int_t i=0; i<nSets; i++){
        fDataSet = vecDataSet.at(i);

        if(useDataRunListForMC && i>=nData) nameOutput = Form("%s/%s/EventQA/%s-%s_EventQARunwise.root",fCutSelection.Data(),fEnergyFlag.Data(),fDataSet.Data(),vecDataSet.at(0).Data());
        else nameOutput = Form("%s/%s/EventQA/%s_EventQARunwise.root",fCutSelection.Data(),fEnergyFlag.Data(),fDataSet.Data());

        TFile* fOutput = new TFile(nameOutput,"RECREATE");
        cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        cout << "Output file: " << nameOutput << endl;
        cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        cout << "Summed NEvents for period: '" << sumEventsAll[i] << "'" << endl;
        cout << "Summed NEvents (for normalization) for period: '" << sumEvents[i] << "'" << endl;
        cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        fLog[i] << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        fLog[i] << "Output file: " << nameOutput << endl;
        fLog[i] << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        fLog[i] << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        fLog[i] << "Summed NEvents for period: '" << sumEventsAll[i] << "'" << endl;
        fLog[i] << "Summed NEvents (for normalization) for period: '" << sumEvents[i] << "'" << endl;
        fLog[i] << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

        for(Int_t h=0; h<(Int_t) vecHistos[i].size(); h++) WriteHistogram(((TH1D*) vecHistos[i].at(h)));

        WriteHistogramTH1DVec(fOutput,vecVertexZ[i],"VertexZ");

        if(i>0) WriteHistogramTH1DVec(fOutput,vecVertexZRatio[i-1],"VertexZ-Ratios");

        WriteHistogramTH1DVec(fOutput,vecPi0Signal[i],"Pi0Signal");
        WriteHistogramTH1DVec(fOutput,vecEtaSignal[i],"EtaSignal");

        WriteHistogramTF1Vec(fOutput,vecPi0Fit[i],"Pi0Fit");
        WriteHistogramTF1Vec(fOutput,vecEtaFit[i],"EtaFit");

        delete[] fitValues[i];

        fOutput->Write();
        fOutput->Close();
        delete fOutput;
        fLog[i].close();
        fEventLog[i].close();
    }

    delete[] vecMissingRuns;

    delete[] sumEvents;
    delete[] sumEventsAll;
    delete[] fLog;
    delete[] fEventLog;

    delete[] vecHistos;
    delete[] vecVertexZ;
    delete[] vecVertexZRatio;
    delete[] vecPi0Signal;
    delete[] vecEtaSignal;
    delete[] vecPi0Fit;
    delete[] vecEtaFit;

    delete[] fitValues;

    TH1::AddDirectory(kTRUE);

    cout << "Done with EventQA_Runwise" << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

    return;

}//end
