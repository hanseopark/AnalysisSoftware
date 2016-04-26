#include "QA.h"

void ClusterQA_Runwise(
                Int_t nSetsIn,
                Int_t nDataIn,
				TString fEnergyFlag,
				TString filePath,
				TString fileName,
				TString* DataSets,
				TString* plotDataSets,
				Int_t mode = 2,
				Int_t cutNr = -1,				// if -1: you have to choose number at runtime
				Int_t doExtQA = 2,				// 0: switched off, 1: normal extQA, 2: with Cell level plots
				Bool_t doEquidistantXaxis= kFALSE,
				Bool_t doTrigger = kTRUE,
				Bool_t doHistsForEverySet = kTRUE,
				Bool_t addSubFolder = kFALSE,
				Bool_t useDataRunListForMC = kFALSE,
				Size_t markerSize = 1,
				TString suffix = "eps"
				)
{
	cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
	cout << "ClusterQA_Runwise" << endl;
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
	TString nameMainDir;

    Float_t xPosLabel = 0.8;
    if (fEnergyFlag.CompareTo("pPb_5.023TeV") == 0)
        xPosLabel = 0.75;    
    
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
		hMarkerStyle[i]=GetDefaultMarkerStyle(fEnergyFlag.Data(),DataSets[i].Data(),"");
		hMarkerColor[i]=GetColorDefaultColor(fEnergyFlag.Data(),DataSets[i].Data(),"");
		hLineColor[i]=GetColorDefaultColor(fEnergyFlag.Data(),DataSets[i].Data(),"");
		hMarkerSize[i]=markerSize;
	}

	Int_t fMode = mode;
	Bool_t isPCMCalo = kTRUE;
	Bool_t isMerged = kFALSE;
	// mode:	0 // new output PCM-PCM
	//			1 // new output PCM dalitz
	//			2 // new output PCM-EMCal
	//			3 // new output PCM-PHOS
	//          4 // new output EMCal-EMCal
	//          5 // new output PHOS-PHOS
	//			9 // old output PCM-PCM
	//			10 // merged EMCal
	//			11 // merged EMCal
	if(fMode == 0 || fMode == 1 || fMode == 9){ cout << "Returning, given mode contains no calo information: " << fMode << endl; return;}
    if(fMode != 2 && fMode != 3) isPCMCalo = kFALSE;
	if(fMode == 10 || fMode == 11) isMerged = kTRUE;

	std::vector<TString> vecRuns;

//****************************** Determine which cut to process ************************************************
	TString fileRuns = Form("runNumbers%s.txt", (vecDataSet.at(0)).Data());
	if(!readin(fileRuns, vecRuns, kFALSE)) {cout << "ERROR, no Run Numbers could be found! Returning..." << endl; return;}
	TFile* fCutFile = new TFile(Form("%s/%s/%s/%s", filePath.Data(), ((TString)vecDataSet.at(0)).Data(), ((TString)vecRuns.at(0)).Data(), fileName.Data()));
	if(fCutFile->IsZombie()) {cout << "ERROR: ROOT file '" << Form("%s/%s/%s/%s", filePath.Data(), ((TString)vecDataSet.at(0)).Data(), ((TString)vecRuns.at(0)).Data(), fileName.Data()) << "' could not be openend, return!" << endl; return;}

	TKey *key;
	TIter next(fCutFile->GetListOfKeys());
	while ((key=(TKey*)next())){
		cout << Form("Found TopDir: '%s' ",key->GetName());
		nameMainDir = key->GetName();
	}
	cout << endl;
	cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
	if(nameMainDir.IsNull() || !nameMainDir.BeginsWith("Gamma")){cout << "ERROR, Unable to obtain valid name of MainDir:|" << nameMainDir.Data() << "|, running in mode: " << fMode << endl; return;}

	TList *listInput =(TList*)fCutFile->Get(nameMainDir.Data());
		listInput->SetOwner(kTRUE);
	if(!listInput) {cout << "ERROR: Could not find main dir: " << nameMainDir.Data() << " in file! Returning..." << endl; return;}
	vector <TString> cuts;
	for(Int_t i = 0; i<listInput->GetSize(); i++){
		TList *listCuts = (TList*)listInput->At(i);
		TString nameCuts = listCuts->GetName();
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

	TString fCutSelection = cuts.at(cutNr);
	TString fEventCutSelection, fGammaCutSelection, fClusterCutSelection, fElectronCutSelection, fMesonCutSelection;
	ReturnSeparatedCutNumberAdvanced(fCutSelection, fEventCutSelection, fGammaCutSelection, fClusterCutSelection, fElectronCutSelection, fMesonCutSelection, fMode);

	Bool_t isTrackMatching = kTRUE;
	TString trackMatchingCut(fClusterCutSelection(GetClusterTrackMatchingCutPosition(fClusterCutSelection),1));
	Int_t trackMatch = trackMatchingCut.Atoi();
	if(trackMatch == 0){
		cout << "INFO: TrackMatching cut found to be '0' for '" << (vecDataSet.at(0)).Data() << "', deactivating generation of track matching histograms!" << endl;
		isTrackMatching = kFALSE;
	}

	TString calo;
	Int_t iCalo = 0;
	Int_t nCaloModules = 0;
	Int_t nCaloCells = 0;
	if(fClusterCutSelection.BeginsWith('1')){
		calo="EMCal"; iCalo=1;
		nCaloModules = 10;
		nCaloCells = 11520;
		if((vecDataSet.at(0)).Contains("LHC10")){
			nCaloModules = 4;
			nCaloCells = 4608;
		}
	}else if(fClusterCutSelection.BeginsWith('2')){
		calo="PHOS"; iCalo=2;
		nCaloModules = 5;
		nCaloCells = 6000;
		if((vecDataSet.at(0)).Contains("LHC10")){
			nCaloModules = 5;
			nCaloCells = 6000;
		}
	}else {cout << "No correct calorimeter type found: " << calo.Data() << ", returning..." << endl; return;}

	TString fCollisionSystem = ReturnFullCollisionsSystem(fEnergyFlag);
	if (fCollisionSystem.CompareTo("") == 0){
		cout << "No correct collision system specification, has been given" << endl;
		return;
	}

	TString fDetectionProcess = ReturnFullTextReconstructionProcess(fMode);

	TString outputDir = Form("%s/%s/ClusterQA/%s/Runwise",fCutSelection.Data(),fEnergyFlag.Data(),suffix.Data());
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
		fileRuns = Form("runNumbers%s.txt", vecDataSet.at(i).Data());
		if(useDataRunListForMC && i>=nData) {
			fileRuns = Form("runNumbers%s-%s.txt", vecDataSet.at(i).Data(),vecDataSet.at(0).Data());
			cout << "Switch useDataRunListForMC is true, reading runs from: " << fileRuns.Data() << endl;
		}
		if(!readin(fileRuns, vecRuns, kFALSE)) {cout << "ERROR, no Run Numbers could be found! Returning..." << endl; return;}

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

	TH1D* hNEvents[nSets];

	TH1D* hCaloNClusters[nSets];
	TH1D* hCaloNClustersQA[nSets];

	TH1D* hClusterMeanEnergy[nSets];
	TH1D* hClusterRMSEnergy[nSets];

	TH1D* hClusterEnergy01[nSets];
	TH1D* hClusterEnergy14[nSets];
	TH1D* hClusterEnergy4[nSets];

	TH1D* hClusterMeanNCells[nSets];
	TH1D* hClusterRMSNCells[nSets];

	TH1D* hClusterMeanDispersion[nSets];
	TH1D* hClusterRMSDispersion[nSets];

	TH1D* hClusterMeanM02[nSets];
	TH1D* hClusterRMSM02[nSets];

	TH1D* hClusterMeanR[nSets];
	TH1D* hClusterRMSR[nSets];

	TH1D* hClusterMeanTime[nSets];
	TH1D* hClusterRMSTime[nSets];

	TH1D* hClusterFractionMatches[nSets];
	TH1D* hClusterFractionMatchesS[nSets];
	TH1D* hClusterFractionMatchesM[nSets];
	TH1D* hClusterFractionMatchesH[nSets];

	TH1D* hClusterMeanDeltaEta[nSets];
	TH1D* hClusterMeanDeltaPhi[nSets];
	TH1D* hClusterRMSDeltaEta[nSets];
	TH1D* hClusterRMSDeltaPhi[nSets];

	TH1D* hClusterMeanM20[nSets];
	TH1D* hClusterRMSM20[nSets];

	TH1D* hClusterMeanConvPhotonPi0_Eta[nSets];
	TH1D* hClusterMeanConvPhotonPi0_Phi[nSets];
	TH1D* hClusterMeanConvPhotonEta_Eta[nSets];
	TH1D* hClusterMeanConvPhotonEta_Phi[nSets];
	TH1D* hClusterRMSConvPhotonPi0_Eta[nSets];
	TH1D* hClusterRMSConvPhotonPi0_Phi[nSets];
	TH1D* hClusterRMSConvPhotonEta_Eta[nSets];
	TH1D* hClusterRMSConvPhotonEta_Phi[nSets];

    TH1D* hClusterMeanNLM[nSets];
    TH1D* hClusterRMSNLM[nSets];

    std::vector<TString>* vecBadCells = new std::vector<TString>[nSets];

    std::vector<TH2D*>* vecClusterEtaPhi = new std::vector<TH2D*>[nSets];
    std::vector<TH2D*>* vecClusterEnergyTime = new std::vector<TH2D*>[nSets];
    std::vector<TH2D*>* vecClusterEVsNCells = new std::vector<TH2D*>[nSets];
    std::vector<TH1D*>* vecClusterIncludedCells = new std::vector<TH1D*>[nSets];
    std::vector<TH1D*>* vecClusterIncludedCellsBefore = new std::vector<TH1D*>[nSets];
    std::vector<TH1D*>* vecClusterEFracCells = new std::vector<TH1D*>[nSets];
    std::vector<TH1D*>* vecClusterEFracCellsBefore = new std::vector<TH1D*>[nSets];
    std::vector<TH2D*>* vecClusterEnergyMeanSigma = new std::vector<TH2D*>[nSets];
    std::vector<TH2D*>* vecClusterTimeMeanSigma = new std::vector<TH2D*>[nSets];
    std::vector<TH2D*>* vecClusterFiredCellIDs = new std::vector<TH2D*>[nSets];
    std::vector<TH2D*>* vecClusterMissingCellIDs = new std::vector<TH2D*>[nSets];
    std::vector<TProfile*>* vecClusterBadCells = new std::vector<TProfile*>[nSets];

    std::vector<TH1D*>* vecClusterEnergy = new std::vector<TH1D*>[nSets];
    std::vector<TH1D*>* vecClusterM02 = new std::vector<TH1D*>[nSets];
    std::vector<TH1D*>* vecClusterM20 = new std::vector<TH1D*>[nSets];
    std::vector<TH1D*>* vecClusterNCells = new std::vector<TH1D*>[nSets];
    std::vector<TH1D*>* vecClusterDispersion = new std::vector<TH1D*>[nSets];

    std::vector<TH1D*>** vecClusterTime = new std::vector<TH1D*>*[nSets];
    for (Int_t i=0; i<nSets; i++) {vecClusterTime[i] = new std::vector<TH1D*>[5];}

    std::vector<TH1D*>* vecClusterDeltaEta = new std::vector<TH1D*>[nSets];
    std::vector<TH1D*>* vecClusterDeltaPhi = new std::vector<TH1D*>[nSets];
    std::vector<TH1D*>* vecClusterPi0ConvPhotonEta = new std::vector<TH1D*>[nSets];
    std::vector<TH1D*>* vecClusterPi0ConvPhotonPhi = new std::vector<TH1D*>[nSets];
    std::vector<TH1D*>* vecClusterEtaConvPhotonEta = new std::vector<TH1D*>[nSets];
    std::vector<TH1D*>* vecClusterEtaConvPhotonPhi = new std::vector<TH1D*>[nSets];

    std::vector<TH1D*>* vecBadCellsEnergy = new std::vector<TH1D*>[nSets];
    std::vector<TH1D*>* vecBadCellsTime = new std::vector<TH1D*>[nSets];

    std::vector<TH1D*>** vecClusterEVsModule = new std::vector<TH1D*>*[nSets];
    std::vector<TH1D*>** vecModuleEVsModule = new std::vector<TH1D*>*[nSets];
    std::vector<TH1D*>** vecClusterNCells100VsModule = new std::vector<TH1D*>*[nSets];
    std::vector<TH1D*>** vecClusterNCells1500VsModule = new std::vector<TH1D*>*[nSets];
    for (Int_t i=0; i<nSets; i++) {
      vecClusterEVsModule[i] = new std::vector<TH1D*>[nCaloModules];
      vecModuleEVsModule[i] = new std::vector<TH1D*>[nCaloModules];
      vecClusterNCells100VsModule[i] = new std::vector<TH1D*>[nCaloModules];
      vecClusterNCells1500VsModule[i] = new std::vector<TH1D*>[nCaloModules];
    }

	Int_t hFBin;
	Int_t hLBin;
	Int_t hNBin;

	if(doEquidistantXaxis)
	{
		hFBin = 0;
		hLBin = globalRuns.size();
		hNBin = globalRuns.size();
	}
	else
	{
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
    std::vector<TString>  vecHistosName;
	TString histoName;

	for(Int_t i=0; i<nSets; i++)
	{
		histoName = "hNEvents";
        if(i==0) vecHistosName.push_back(histoName);
		hNEvents[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hNEvents; Run Number ; # of Events",hNBin,hFBin,hLBin);
		EditTH1(globalRuns, doEquidistantXaxis, hNEvents[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
		vecHistos[i].push_back(hNEvents[i]);

		histoName = "hCaloNClusters";
        if(i==0) vecHistosName.push_back(histoName);
		hCaloNClusters[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hCaloNClusters; Run Number ; #frac{1}{N_{Events}} N_{Clusters before Cuts}",hNBin,hFBin,hLBin);
		EditTH1(globalRuns, doEquidistantXaxis, hCaloNClusters[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
		vecHistos[i].push_back(hCaloNClusters[i]);

		histoName = "hCaloNClustersQA";
        if(i==0) vecHistosName.push_back(histoName);
		hCaloNClustersQA[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hCaloNClustersQA; Run Number ; #frac{1}{N_{Events}} N_{Clusters}",hNBin,hFBin,hLBin);
		EditTH1(globalRuns, doEquidistantXaxis, hCaloNClustersQA[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
		vecHistos[i].push_back(hCaloNClustersQA[i]);

		histoName = "hClusterEnergy-Mean";
        if(i==0) vecHistosName.push_back(histoName);
		hClusterMeanEnergy[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hClusterMeanEnergy; Run Number ; #bar{#lower[0.1]{E}}_{Cluster} (GeV)",hNBin,hFBin,hLBin);
		EditTH1(globalRuns, doEquidistantXaxis, hClusterMeanEnergy[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
		vecHistos[i].push_back(hClusterMeanEnergy[i]);

		histoName = "hClusterEnergy-RMS";
        if(i==0) vecHistosName.push_back(histoName);
		hClusterRMSEnergy[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hClusterRMSEnergy; Run Number ; #sigma_{E_{Cluster}} (GeV)",hNBin,hFBin,hLBin);
		EditTH1(globalRuns, doEquidistantXaxis, hClusterRMSEnergy[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
		vecHistos[i].push_back(hClusterRMSEnergy[i]);

		histoName = "hClusterEnergy-01";
        if(i==0) vecHistosName.push_back(histoName);
		hClusterEnergy01[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hClusterEnergy01; Run Number ; #frac{1}{N_{Events}} N_{Clusters, E<1 GeV}",hNBin,hFBin,hLBin);
		EditTH1(globalRuns, doEquidistantXaxis, hClusterEnergy01[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
		vecHistos[i].push_back(hClusterEnergy01[i]);

		histoName = "hClusterEnergy-14";
        if(i==0) vecHistosName.push_back(histoName);
		hClusterEnergy14[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hClusterEnergy14; Run Number ; #frac{1}{N_{Events}} N_{Clusters, 1<E<4 GeV}",hNBin,hFBin,hLBin);
		EditTH1(globalRuns, doEquidistantXaxis, hClusterEnergy14[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
		vecHistos[i].push_back(hClusterEnergy14[i]);

		histoName = "hClusterEnergy-4";
        if(i==0) vecHistosName.push_back(histoName);
		hClusterEnergy4[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hClusterEnergy4; Run Number ; #frac{1}{N_{Events}} N_{Clusters, E>4 GeV}",hNBin,hFBin,hLBin);
		EditTH1(globalRuns, doEquidistantXaxis, hClusterEnergy4[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
		vecHistos[i].push_back(hClusterEnergy4[i]);

		histoName = "hClusterNCells-Mean";
        if(i==0) vecHistosName.push_back(histoName);
		hClusterMeanNCells[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hClusterMeanNCells; Run Number ; #bar{#lower[0.1]{N}}_{Cells, Cluster}",hNBin,hFBin,hLBin);
		EditTH1(globalRuns, doEquidistantXaxis, hClusterMeanNCells[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
		vecHistos[i].push_back(hClusterMeanNCells[i]);

		histoName = "hClusterNCells-RMS";
        if(i==0) vecHistosName.push_back(histoName);
		hClusterRMSNCells[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hClusterRMSNCells; Run Number ; #sigma_{N_{Cells, Cluster}}",hNBin,hFBin,hLBin);
		EditTH1(globalRuns, doEquidistantXaxis, hClusterRMSNCells[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
		vecHistos[i].push_back(hClusterRMSNCells[i]);

		histoName = "hClusterDispersion-Mean";
        if(i==0) vecHistosName.push_back(histoName);
		hClusterMeanDispersion[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hClusterMeanDispersion; Run Number ; #bar{#lower[0.1]{Dispersion}}_{Cluster}",hNBin,hFBin,hLBin);
		EditTH1(globalRuns, doEquidistantXaxis, hClusterMeanDispersion[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
		vecHistos[i].push_back(hClusterMeanDispersion[i]);

		histoName = "hClusterDispersion-RMS";
        if(i==0) vecHistosName.push_back(histoName);
		hClusterRMSDispersion[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hClusterRMSDispersion; Run Number ; #sigma_{Dispersion_{Cluster}}",hNBin,hFBin,hLBin);
		EditTH1(globalRuns, doEquidistantXaxis, hClusterRMSDispersion[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
		vecHistos[i].push_back(hClusterRMSDispersion[i]);

		histoName = "hClusterM02-Mean";
        if(i==0) vecHistosName.push_back(histoName);
		hClusterMeanM02[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hClusterMeanM02; Run Number ; #bar{#lower[0.1]{#lambda_{0}^{2}}}",hNBin,hFBin,hLBin);
		EditTH1(globalRuns, doEquidistantXaxis, hClusterMeanM02[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
		vecHistos[i].push_back(hClusterMeanM02[i]);

		histoName = "hClusterM02-RMS";
        if(i==0) vecHistosName.push_back(histoName);
		hClusterRMSM02[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hClusterRMSM02; Run Number ; #sigma_{#lambda_{0}^{2}}",hNBin,hFBin,hLBin);
		EditTH1(globalRuns, doEquidistantXaxis, hClusterRMSM02[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
		vecHistos[i].push_back(hClusterRMSM02[i]);

		histoName = "hClusterR-Mean";
        if(i==0) vecHistosName.push_back(histoName);
		hClusterMeanR[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hClusterMeanR; Run Number ; #bar{#lower[0.1]{R}}_{Cluster}",hNBin,hFBin,hLBin);
		EditTH1(globalRuns, doEquidistantXaxis, hClusterMeanR[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
		vecHistos[i].push_back(hClusterMeanR[i]);

		histoName = "hClusterR-RMS";
        if(i==0) vecHistosName.push_back(histoName);
		hClusterRMSR[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hClusterRMSR; Run Number ; #sigma_{R_{Cluster}}",hNBin,hFBin,hLBin);
		EditTH1(globalRuns, doEquidistantXaxis, hClusterRMSR[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
		vecHistos[i].push_back(hClusterRMSR[i]);

		histoName = "hClusterTime-Mean";
        if(i==0) vecHistosName.push_back(histoName);
		hClusterMeanTime[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hClusterMeanTime; Run Number ; #bar{#lower[0.1]{t}}_{Cluster} (s)",hNBin,hFBin,hLBin);
		EditTH1(globalRuns, doEquidistantXaxis, hClusterMeanTime[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
		vecHistos[i].push_back(hClusterMeanTime[i]);

		histoName = "hClusterTime-RMS";
        if(i==0) vecHistosName.push_back(histoName);
		hClusterRMSTime[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hClusterRMSTime; Run Number ; #sigma_{t_{Cluster}} (s)",hNBin,hFBin,hLBin);
		EditTH1(globalRuns, doEquidistantXaxis, hClusterRMSTime[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
		vecHistos[i].push_back(hClusterRMSTime[i]);

		histoName = "hCluster-FractionMatches";
        if(i==0) vecHistosName.push_back(histoName);
		hClusterFractionMatches[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hClusterFractionMatches; Run Number ; #frac{N_{Removed Mother Candidates by Cluster-Track Matching}}{N_{Mother Candidates}}",hNBin,hFBin,hLBin);
		EditTH1(globalRuns, doEquidistantXaxis, hClusterFractionMatches[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
		vecHistos[i].push_back(hClusterFractionMatches[i]);

		histoName = "hCluster-FractionMatchesS";
        if(i==0) vecHistosName.push_back(histoName);
		hClusterFractionMatchesS[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hClusterFractionMatchesS; Run Number ; #frac{N_{Removed Mother Candidates by Cluster-Track Matching}}{N_{Mother Candidates}} #scale[0.8]{for #it{p}_{T,Pair}<1.5 GeV/#it{c}}",hNBin,hFBin,hLBin);
		EditTH1(globalRuns, doEquidistantXaxis, hClusterFractionMatchesS[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
		vecHistos[i].push_back(hClusterFractionMatchesS[i]);

		histoName = "hCluster-FractionMatchesM";
        if(i==0) vecHistosName.push_back(histoName);
		hClusterFractionMatchesM[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hClusterFractionMatchesM; Run Number ; #frac{N_{Removed Mother Candidates by Cluster-Track Matching}}{N_{Mother Candidates}} #scale[0.8]{for 1.5<#it{p}_{T,Pair}<2.5 GeV/#it{c}}",hNBin,hFBin,hLBin);
		EditTH1(globalRuns, doEquidistantXaxis, hClusterFractionMatchesM[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
		vecHistos[i].push_back(hClusterFractionMatchesM[i]);

		histoName = "hCluster-FractionMatchesH";
        if(i==0) vecHistosName.push_back(histoName);
		hClusterFractionMatchesH[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hClusterFractionMatchesH; Run Number ; #frac{N_{Removed Mother Candidates by Cluster-Track Matching}}{N_{Mother Candidates}} #scale[0.8]{for #it{p}_{T,Pair}>2.5 GeV/#it{c}}",hNBin,hFBin,hLBin);
		EditTH1(globalRuns, doEquidistantXaxis, hClusterFractionMatchesH[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
		vecHistos[i].push_back(hClusterFractionMatchesH[i]);

		histoName = "hClusterDeltaEta-Mean";
        if(i==0) vecHistosName.push_back(histoName);
		hClusterMeanDeltaEta[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hClusterMeanDeltaEta; Run Number ; #bar{#lower[0.1]{#Delta#eta}}_{Cluster, Tracks}",hNBin,hFBin,hLBin);
		EditTH1(globalRuns, doEquidistantXaxis, hClusterMeanDeltaEta[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
		vecHistos[i].push_back(hClusterMeanDeltaEta[i]);

		histoName = "hClusterDeltaPhi-Mean";
        if(i==0) vecHistosName.push_back(histoName);
		hClusterMeanDeltaPhi[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hClusterMeanDeltaPhi; Run Number ; #bar{#lower[0.1]{#Delta#phi}}_{Cluster, Tracks}",hNBin,hFBin,hLBin);
		EditTH1(globalRuns, doEquidistantXaxis, hClusterMeanDeltaPhi[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
		vecHistos[i].push_back(hClusterMeanDeltaPhi[i]);

		histoName = "hClusterDeltaEta-RMS";
        if(i==0) vecHistosName.push_back(histoName);
		hClusterRMSDeltaEta[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hClusterRMSDeltaEta; Run Number ; #sigma_{#Delta#eta_{Cluster, Tracks}}",hNBin,hFBin,hLBin);
		EditTH1(globalRuns, doEquidistantXaxis, hClusterRMSDeltaEta[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
		vecHistos[i].push_back(hClusterRMSDeltaEta[i]);

		histoName = "hClusterDeltaPhi-RMS";
        if(i==0) vecHistosName.push_back(histoName);
		hClusterRMSDeltaPhi[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hClusterRMSDeltaPhi; Run Number ; #sigma_{#Delta#phi_{Cluster, Tracks}}",hNBin,hFBin,hLBin);
		EditTH1(globalRuns, doEquidistantXaxis, hClusterRMSDeltaPhi[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
		vecHistos[i].push_back(hClusterRMSDeltaPhi[i]);

		histoName = "hClusterM20-Mean";
        if(i==0) vecHistosName.push_back(histoName);
		hClusterMeanM20[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hClusterMeanM20; Run Number ; #bar{#lower[0.1]{M20}}_{Cluster}",hNBin,hFBin,hLBin);
		EditTH1(globalRuns, doEquidistantXaxis, hClusterMeanM20[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
		vecHistos[i].push_back(hClusterMeanM20[i]);

		histoName = "hClusterM20-RMS";
        if(i==0) vecHistosName.push_back(histoName);
		hClusterRMSM20[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hClusterRMSM20; Run Number ; #sigma_{M20_{Cluster}}",hNBin,hFBin,hLBin);
		EditTH1(globalRuns, doEquidistantXaxis, hClusterRMSM20[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
		vecHistos[i].push_back(hClusterRMSM20[i]);

		histoName = "hClusterConvPhotonPi0_Eta-Mean";
        if(i==0) vecHistosName.push_back(histoName);
		hClusterMeanConvPhotonPi0_Eta[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hClusterConvPhotonPi0_Eta; Run Number ; #bar{#lower[0.1]{#eta}}_{#gamma_{conv} under #pi^{0}-peak}",hNBin,hFBin,hLBin);
		EditTH1(globalRuns, doEquidistantXaxis, hClusterMeanConvPhotonPi0_Eta[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
		vecHistos[i].push_back(hClusterMeanConvPhotonPi0_Eta[i]);

		histoName = "hClusterConvPhotonPi0_Phi-Mean";
        if(i==0) vecHistosName.push_back(histoName);
		hClusterMeanConvPhotonPi0_Phi[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hClusterConvPhotonPi0_Phi; Run Number ; #bar{#lower[0.1]{#phi}}_{#gamma_{conv} under #pi^{0}-peak}",hNBin,hFBin,hLBin);
		EditTH1(globalRuns, doEquidistantXaxis, hClusterMeanConvPhotonPi0_Phi[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
		vecHistos[i].push_back(hClusterMeanConvPhotonPi0_Phi[i]);

		histoName = "hClusterConvPhotonEta_Eta-Mean";
        if(i==0) vecHistosName.push_back(histoName);
		hClusterMeanConvPhotonEta_Eta[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hClusterConvPhotonEta_Eta; Run Number ; #bar{#lower[0.1]{#eta}}_{#gamma_{conv} under #eta-peak}",hNBin,hFBin,hLBin);
		EditTH1(globalRuns, doEquidistantXaxis, hClusterMeanConvPhotonEta_Eta[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
		vecHistos[i].push_back(hClusterMeanConvPhotonEta_Eta[i]);

		histoName = "hClusterConvPhotonEta_Phi-Mean";
        if(i==0) vecHistosName.push_back(histoName);
		hClusterMeanConvPhotonEta_Phi[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hClusterConvPhotonEta_Phi; Run Number ; #bar{#lower[0.1]{#phi}}_{#gamma_{conv} under #eta-peak}",hNBin,hFBin,hLBin);
		EditTH1(globalRuns, doEquidistantXaxis, hClusterMeanConvPhotonEta_Phi[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
		vecHistos[i].push_back(hClusterMeanConvPhotonEta_Phi[i]);

		histoName = "hClusterConvPhotonPi0_Eta-RMS";
        if(i==0) vecHistosName.push_back(histoName);
		hClusterRMSConvPhotonPi0_Eta[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hClusterConvPhotonPi0_Eta; Run Number ; #sigma_{#eta_{#gamma_{conv} under #pi^{0}-peak}}",hNBin,hFBin,hLBin);
		EditTH1(globalRuns, doEquidistantXaxis, hClusterRMSConvPhotonPi0_Eta[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
		vecHistos[i].push_back(hClusterRMSConvPhotonPi0_Eta[i]);

		histoName = "hClusterConvPhotonPi0_Phi-RMS";
        if(i==0) vecHistosName.push_back(histoName);
		hClusterRMSConvPhotonPi0_Phi[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hClusterConvPhotonPi0_Phi; Run Number ; #sigma_{#phi_{#gamma_{conv} under #pi^{0}-peak}}",hNBin,hFBin,hLBin);
		EditTH1(globalRuns, doEquidistantXaxis, hClusterRMSConvPhotonPi0_Phi[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
		vecHistos[i].push_back(hClusterRMSConvPhotonPi0_Phi[i]);

		histoName = "hClusterConvPhotonEta_Eta-RMS";
        if(i==0) vecHistosName.push_back(histoName);
		hClusterRMSConvPhotonEta_Eta[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hClusterConvPhotonEta_Eta; Run Number ; #sigma_{#eta_{#gamma_{conv} under #eta-peak}}",hNBin,hFBin,hLBin);
		EditTH1(globalRuns, doEquidistantXaxis, hClusterRMSConvPhotonEta_Eta[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
		vecHistos[i].push_back(hClusterRMSConvPhotonEta_Eta[i]);

		histoName = "hClusterConvPhotonEta_Phi-RMS";
        if(i==0) vecHistosName.push_back(histoName);
		hClusterRMSConvPhotonEta_Phi[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hClusterConvPhotonEta_Phi; Run Number ; #sigma_{#phi_{#gamma_{conv} under #eta-peak}}",hNBin,hFBin,hLBin);
		EditTH1(globalRuns, doEquidistantXaxis, hClusterRMSConvPhotonEta_Phi[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
		vecHistos[i].push_back(hClusterRMSConvPhotonEta_Phi[i]);

        histoName = "hClusterNLM-Mean";
        if(i==0) vecHistosName.push_back(histoName);
        hClusterMeanNLM[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hClusterMeanNLM; Run Number ; #bar{#lower[0.1]{NLM}}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hClusterMeanNLM[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hClusterMeanNLM[i]);

        histoName = "hClusterNLM-RMS";
        if(i==0) vecHistosName.push_back(histoName);
        hClusterRMSNLM[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hClusterRMSNLM; Run Number ; #sigma_{NLM}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hClusterRMSNLM[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hClusterRMSNLM[i]);
	}

//****************************** Looping over DataSets ************************************************

    fstream* fLog = new fstream[nSets];
    fstream* fLogRunwiseHotCells = new fstream[nSets];
    fstream* fLogRunwiseColdCells = new fstream[nSets];
	for(Int_t iStr=0; iStr<nSets; iStr++){
		if(useDataRunListForMC && iStr>=nData) {
			fLog[iStr].open(Form("%s/A-%s-%s.log",outputDir.Data(),DataSets[iStr].Data(),DataSets[0].Data()), ios::out);
			if(doExtQA==2) fLogRunwiseHotCells[iStr].open(Form("%s/HotCellsRunwise-%s-%s.log",outputDir.Data(),DataSets[iStr].Data(),DataSets[0].Data()), ios::out);
			if(doExtQA==2) fLogRunwiseColdCells[iStr].open(Form("%s/ColdCellsRunwise-%s-%s.log",outputDir.Data(),DataSets[iStr].Data(),DataSets[0].Data()), ios::out);
		}else{
			fLog[iStr].open(Form("%s/A-%s.log",outputDir.Data(),DataSets[iStr].Data()), ios::out);
			if(doExtQA==2) fLogRunwiseHotCells[iStr].open(Form("%s/HotCellsRunwise-%s.log",outputDir.Data(),DataSets[iStr].Data()), ios::out);
			if(doExtQA==2) fLogRunwiseColdCells[iStr].open(Form("%s/ColdCellsRunwise-%s.log",outputDir.Data(),DataSets[iStr].Data()), ios::out);
		}
		fLog[iStr] << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
		fLog[iStr] << DataSets[iStr].Data() << endl;
		fLog[iStr] << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
		fLog[iStr] << fCollisionSystem.Data() << endl;
		fLog[iStr] << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
		fLog[iStr] << "processed cut: " << fCutSelection.Data() << endl;
		fLog[iStr] << calo.Data() << ", Modules: " << nCaloModules << ", Cells: " << nCaloCells << endl;
		fLog[iStr] << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
	}

	const Int_t fNBinsClusterPt 			= 60;
	Double_t fBinsClusterPt[fNBinsClusterPt+1] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
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

	Double_t minT_Energy[5]={0,0,1.2,4,16};
	Double_t maxT_Energy[5]={40,1.2,4,16,40};

	TString sT[2]={"Pi0","Eta"};
	TString sL[2]={"#pi^{0}","#eta"};

	TString fRootFile;
	TString fDataSet;
	TString fRunNumber;
	TString fileCells;
	Int_t bin = -1;

    std::vector<TString>* vecMissingRuns = new std::vector<TString>[nSets];

	for(Int_t i=0; i<nSets; i++)
	{
		vecRuns.clear();
		fDataSet = vecDataSet.at(i);
		fileRuns = Form("runNumbers%s.txt", fDataSet.Data());
		if(useDataRunListForMC && i>=nData) {
			fileRuns = Form("runNumbers%s-%s.txt", vecDataSet.at(i).Data(),vecDataSet.at(0).Data());
			cout << "Switch useDataRunListForMC is true, reading runs from: " << fileRuns.Data() << endl;
		}
		if(!readin(fileRuns, vecRuns)) {cout << "ERROR, no Run Numbers could be found! Returning..." << endl; return;}

		fileCells = Form("%s/%s/ClusterQA/%s/%s/Cells/%s.log",fCutSelection.Data(),fEnergyFlag.Data(),suffix.Data(),DataSets[i].Data(),DataSets[i].Data());
		readin(fileCells, vecBadCells[i], kTRUE, kTRUE);
		fLog[i] << "Read in " << vecBadCells[i].size() << " bad cells: ";
		for(Int_t iBad=0; iBad<(Int_t) vecBadCells[i].size(); iBad++) fLog[i] << vecBadCells[i].at(iBad) << ", ";
		fLog[i] << endl;

//****************************** Looping over Runs ************************************************
		cout << endl;
		cout << "\t----------------------------------------------------------------------------" << endl;
		cout << "\tLooping over Runs of DataSet |" << (vecDataSet.at(i)).Data() << "|" << endl;
		cout << "\t----------------------------------------------------------------------------" << endl;
		cout << endl;
		fLog[i] << "Looping over Runs:" << endl;

		vecMissingRuns[i].clear();
		for(Int_t j=0; j<(Int_t) vecRuns.size(); j++)
		{
			fRunNumber = vecRuns.at(j);
			fRootFile = Form("%s/%s/%s/%s", filePath.Data(), fDataSet.Data(), fRunNumber.Data(), fileName.Data());
			TFile* RootFile = new TFile(fRootFile.Data(),"READ");
			if(RootFile->IsZombie()) {vecMissingRuns[i].push_back(fRunNumber); cout << "INFO: ROOT file '" << fRootFile.Data() << "' could not be openend, continue!" << endl; continue;}

			cout << endl;
			cout << "\t\t----------------------------------------------------------------------------" << endl;
			cout << Form("\t\tRun %s", fRunNumber.Data()) << endl;
			cout << "\t\tProcessing file: " << fRootFile.Data() << endl;
			cout << "\t\t----------------------------------------------------------------------------" << endl;
			cout << endl;

			TList* TopDir = (TList*) RootFile->Get(nameMainDir.Data());
				if(TopDir == NULL) {cout << "ERROR: TopDir not Found"<<endl; return;}
				else TopDir->SetOwner(kTRUE);
			TList* TopContainer = (TList*) TopDir->FindObject(Form("Cut Number %s",fCutSelection.Data()));
				if(TopContainer == NULL) {cout << "ERROR: " << Form("Cut Number %s",fCutSelection.Data()) << " not found in File" << endl; return;}
				else TopContainer->SetOwner(kTRUE);
			TList* ESDContainer = (TList*) TopContainer->FindObject(Form("%s ESD histograms",fCutSelection.Data()));
				if(ESDContainer == NULL) {cout << "ERROR: " << Form("%s ESD histograms",fCutSelection.Data()) << " not found in File" << endl; return;}
				else ESDContainer->SetOwner(kTRUE);
			TList* CaloCutsContainer = (TList*) TopContainer->FindObject(Form("CaloCuts_%s",fClusterCutSelection.Data()));
				if(CaloCutsContainer == NULL) {cout << "ERROR: " << Form("CaloCuts_%s",fClusterCutSelection.Data()) << " not found in File" << endl; return;}
				else CaloCutsContainer->SetOwner(kTRUE);
			TList* ConvCutsContainer = (TList*) TopContainer->FindObject(Form("ConvCuts_%s",fGammaCutSelection.Data()));
				if(isPCMCalo && ConvCutsContainer == NULL) {cout << "ERROR: " << Form("ConvCuts_%s",fGammaCutSelection.Data()) << " not found in File" << endl; return;}
                else if(ConvCutsContainer) ConvCutsContainer->SetOwner(kTRUE);
			TList* CaloExtQAContainer = (TList*) TopContainer->FindObject(Form("CaloExtQA_%s",fClusterCutSelection.Data()));
				if(CaloExtQAContainer == NULL) {
					cout << "WARNING: " << Form("CaloExtQA_%s",fClusterCutSelection.Data()) << " not found in File, using CaloCuts-Container" << endl;
					CaloExtQAContainer = CaloCutsContainer;
				}else CaloExtQAContainer->SetOwner(kTRUE);
	//---------
			//if(extendedRunwiseQA) ClusterQA(fRootFile,"","",nameMainDir,fEnergyFlag,plotDataSets[i],"","",fCutSelection,fRunNumber,suffix,mode);
	//---------
			if(doEquidistantXaxis) bin = mapBin[fRunNumber];
			else bin = fRunNumber.Atoi() - hFBin;
			Double_t nEvents = GetNEvents((TH1*) ESDContainer->FindObject("NEvents"),kFALSE);
	//---------
			fLog[i] << "----------------------------------------------------------------------------" << endl;
			fLog[i] << "Processing file: " << fRootFile.Data() << endl;
			fLog[i] << "Run: " << fRunNumber.Data() << ", with NEvents: " << nEvents << endl;
			fLog[i] << "----------------------------------------------------------------------------" << endl;
			if(doExtQA==2) fLogRunwiseHotCells[i] << "Run-" << fRunNumber.Data() << endl;
			if(doExtQA==2) fLogRunwiseColdCells[i] << "Run-" << fRunNumber.Data() << endl;
	//---------
			if( nEvents < 1 ){cout << "Warning: number of accepted events in run: " << nEvents << "! Setting nEvents to 1..." << endl; nEvents = 1;}
	//---------
			if(hNEvents){
				hNEvents[i]->SetBinContent(bin, nEvents);
			}else cout << "Info: Object |NEvents| could not be found! Skipping Fill..." << endl;
	//---------
			TH1F* EOC = (TH1F*) CaloCutsContainer->FindObject(Form("EnergyOfCluster_afterClusterQA %s", fClusterCutSelection.Data()));
			if(EOC){
				hClusterMeanEnergy[i]->SetBinContent(bin, EOC->GetMean());
				hClusterMeanEnergy[i]->SetBinError(bin, EOC->GetMeanError());
				hClusterRMSEnergy[i]->SetBinContent(bin, EOC->GetRMS());
				hClusterRMSEnergy[i]->SetBinError(bin, EOC->GetRMSError());
			}else cout << "Info: Object |EnergyOfCluster_afterClusterQA| could not be found! Skipping Fill..." << endl;
	//---------
			if(EOC){
				Double_t eps = 0.0001;
				Double_t Energy01 = EOC->Integral(1,EOC->GetXaxis()->FindBin(1-eps));
				Double_t Energy14 = EOC->Integral(EOC->GetXaxis()->FindBin(1+eps),EOC->GetXaxis()->FindBin(4-eps));
				Double_t Energy4 = EOC->Integral(EOC->GetXaxis()->FindBin(4+eps), EOC->GetXaxis()->GetLast());
				hClusterEnergy01[i]->SetBinContent(bin, Energy01 / nEvents);
				hClusterEnergy01[i]->SetBinError(bin, sqrt(Energy01) / nEvents);
				hClusterEnergy14[i]->SetBinContent(bin, Energy14 / nEvents);
				hClusterEnergy14[i]->SetBinError(bin, sqrt(Energy14) / nEvents);
				hClusterEnergy4[i]->SetBinContent(bin, Energy4 / nEvents);
				hClusterEnergy4[i]->SetBinError(bin, sqrt(Energy4) / nEvents);
			}else cout << "Info: Object |EnergyOfCluster_afterClusterQA| could not be found! Skipping Fill..." << endl;
	//---------
			TH1F* NCells = (TH1F*) CaloCutsContainer->FindObject(Form("NCellPerCluster_afterClusterQA %s", fClusterCutSelection.Data()));
			if(NCells){
				hClusterMeanNCells[i]->SetBinContent(bin, NCells->GetMean());
				hClusterMeanNCells[i]->SetBinError(bin, NCells->GetMeanError());
				hClusterRMSNCells[i]->SetBinContent(bin, NCells->GetRMS());
				hClusterRMSNCells[i]->SetBinError(bin, NCells->GetRMSError());
			}else cout << "Info: Object |NCellPerCluster_afterClusterQA| could not be found! Skipping Fill..." << endl;
	//---------
			TH1F* Dispersion = (TH1F*) CaloCutsContainer->FindObject(Form("Dispersion_afterClusterQA %s", fClusterCutSelection.Data()));
			if(Dispersion){
				hClusterMeanDispersion[i]->SetBinContent(bin, Dispersion->GetMean());
				hClusterMeanDispersion[i]->SetBinError(bin, Dispersion->GetMeanError());
				hClusterRMSDispersion[i]->SetBinContent(bin, Dispersion->GetRMS());
				hClusterRMSDispersion[i]->SetBinError(bin, Dispersion->GetRMSError());
			}else cout << "Info: Object |Dispersion_afterClusterQA| could not be found! Skipping Fill..." << endl;
	//---------
			TH1F* M02 = (TH1F*) CaloCutsContainer->FindObject(Form("M02_afterClusterQA %s", fClusterCutSelection.Data()));
			if(M02){
				hClusterMeanM02[i]->SetBinContent(bin, M02->GetMean());
				hClusterMeanM02[i]->SetBinError(bin, M02->GetMeanError());
				hClusterRMSM02[i]->SetBinContent(bin, M02->GetRMS());
				hClusterRMSM02[i]->SetBinError(bin, M02->GetRMSError());
			}else cout << "Info: Object |M02_afterClusterQA| could not be found! Skipping Fill..." << endl;
	//---------
			TH1F* R = (TH1F*) CaloCutsContainer->FindObject(Form("R_Cluster_afterClusterQA %s", fClusterCutSelection.Data()));
			if(R){
				hClusterMeanR[i]->SetBinContent(bin, R->GetMean());
				hClusterMeanR[i]->SetBinError(bin, R->GetMeanError());
				hClusterRMSR[i]->SetBinContent(bin, R->GetRMS());
				hClusterRMSR[i]->SetBinError(bin, R->GetRMSError());
			}else cout << "Info: Object |R_Cluster_afterClusterQA| could not be found! Skipping Fill..." << endl;
	//---------
			TH1F* CaloACC = (TH1F*) CaloCutsContainer->FindObject(Form("AcceptanceCuts %s", fClusterCutSelection.Data()));
			TH1F* CaloIPS = (TH1F*) CaloCutsContainer->FindObject(Form("IsPhotonSelected %s", fClusterCutSelection.Data()));
			if(CaloACC && CaloIPS){
				Double_t CaloNClusters = CaloACC->GetBinContent(1);
				Double_t CaloNClustersQA = CaloIPS->GetBinContent(5);
				hCaloNClusters[i]->SetBinContent(bin, CaloNClusters / nEvents);
				hCaloNClusters[i]->SetBinError(bin, sqrt(CaloNClusters) / nEvents);
				hCaloNClustersQA[i]->SetBinContent(bin, CaloNClustersQA / nEvents);
				hCaloNClustersQA[i]->SetBinError(bin, sqrt(CaloNClustersQA) / nEvents);
			}else cout << "Info: Object |AcceptanceCuts| or |IsPhotonSelected| could not be found! Skipping Fill..." << endl;
	//---------
			TH2F* Time = (TH2F*) CaloCutsContainer->FindObject(Form("ClusterTimeVsE_afterClusterQA %s", fClusterCutSelection.Data()));
			if(Time){
				hClusterMeanTime[i]->SetBinContent(bin, Time->GetMean());
				hClusterMeanTime[i]->SetBinError(bin, Time->GetMeanError());
				hClusterRMSTime[i]->SetBinContent(bin, Time->GetRMS());
				hClusterRMSTime[i]->SetBinError(bin, Time->GetRMSError());
			}else cout << "Info: Object |ClusterTimeVsE_afterClusterQA| could not be found! Skipping Fill..." << endl;
	//---------
			if(isPCMCalo){
				TH2F* ESD_Mother = (TH2F*) ESDContainer->FindObject("ESD_Mother_InvMass_Pt");
				TH2F* ESD_Mother_Matched = (TH2F*) ESDContainer->FindObject("ESD_MotherMatched_InvMass_Pt");
				if(ESD_Mother && ESD_Mother_Matched){
					CalculateFractionMatches(hClusterFractionMatches[i], ESD_Mother, ESD_Mother_Matched, bin, 0, 30);
                    CalculateFractionMatches(hClusterFractionMatchesS[i], ESD_Mother, ESD_Mother_Matched, bin, 0, 2.25);
                    CalculateFractionMatches(hClusterFractionMatchesM[i], ESD_Mother, ESD_Mother_Matched, bin, 2.25, 4);
                    CalculateFractionMatches(hClusterFractionMatchesH[i], ESD_Mother, ESD_Mother_Matched, bin, 4, 30);
				}else cout << "Info: Object |ESD_Mother_InvMass_Pt| or |ESD_MotherMatched_InvMass_Pt| could not be found! Skipping Fill..." << endl;
			}
	//---------
			if(isTrackMatching){
				TH2F* DeltaEtaPhi = (TH2F*) CaloCutsContainer->FindObject(Form("dEtaVsdPhi_beforeClusterQA %s", fClusterCutSelection.Data()));
				if(DeltaEtaPhi){
					hClusterMeanDeltaEta[i]->SetBinContent(bin, DeltaEtaPhi->GetMean(1));
					hClusterMeanDeltaEta[i]->SetBinError(bin, DeltaEtaPhi->GetMeanError(1));
					hClusterMeanDeltaPhi[i]->SetBinContent(bin, DeltaEtaPhi->GetMean(2));
					hClusterMeanDeltaPhi[i]->SetBinError(bin, DeltaEtaPhi->GetMeanError(2));
					hClusterRMSDeltaEta[i]->SetBinContent(bin, DeltaEtaPhi->GetRMS(1));
					hClusterRMSDeltaEta[i]->SetBinError(bin, DeltaEtaPhi->GetRMSError(1));
					hClusterRMSDeltaPhi[i]->SetBinContent(bin, DeltaEtaPhi->GetRMS(2));
					hClusterRMSDeltaPhi[i]->SetBinError(bin, DeltaEtaPhi->GetRMSError(2));
				}else cout << "Info: Object |dEtaVsdPhi_beforeClusterQA| could not be found! Skipping Fill..." << endl;
			}
	//---------
			TH1F* M20 = (TH1F*) CaloCutsContainer->FindObject(Form("M20_afterClusterQA %s", fClusterCutSelection.Data()));
			if(M20){
				hClusterMeanM20[i]->SetBinContent(bin, M20->GetMean());
				hClusterMeanM20[i]->SetBinError(bin, M20->GetMeanError());
				hClusterRMSM20[i]->SetBinContent(bin, M20->GetRMS());
				hClusterRMSM20[i]->SetBinError(bin, M20->GetRMSError());
			}else cout << "Info: Object |M20_afterClusterQA| could not be found! Skipping Fill..." << endl;
	//---------
			TH2D* ClusEtaPhi = (TH2D*) CaloCutsContainer->FindObject(Form("EtaPhi_afterClusterQA %s", fClusterCutSelection.Data()));
			if(ClusEtaPhi){
				TH2D* tempClusEtaPhi = new TH2D(*ClusEtaPhi);
				tempClusEtaPhi->Scale(1 / nEvents);
				tempClusEtaPhi->Scale(1 / GetMeanTH2(tempClusEtaPhi));
				tempClusEtaPhi->SetName(Form("%s_EtaPhi_%s",fRunNumber.Data(),fDataSet.Data()));
				tempClusEtaPhi->SetTitle(fRunNumber);
				tempClusEtaPhi->GetXaxis()->SetTitle("#phi");
				tempClusEtaPhi->GetYaxis()->SetTitle("#eta");
				tempClusEtaPhi->GetZaxis()->SetRangeUser(0, 2);
				vecClusterEtaPhi[i].push_back(tempClusEtaPhi);
			}else cout << "Info: Object |EtaPhi_afterClusterQA| could not be found! Skipping..." << endl;
	//---------
			if(i<nData)
			{
				TH2D* ClusEnergyTime = (TH2D*) CaloCutsContainer->FindObject(Form("ClusterTimeVsE_afterClusterQA %s", fClusterCutSelection.Data()));
				if(ClusEnergyTime){
					TH2D* tempClusEnergyTime = new TH2D(*ClusEnergyTime);
					//tempClusEnergyTime->Scale(1 / nEvents);
					tempClusEnergyTime->SetName(Form("%s_EnergyTime_%s",fRunNumber.Data(),fDataSet.Data()));
					tempClusEnergyTime->SetTitle(fRunNumber);
					Int_t min = 0, max = 0;
					GetMinMaxBin(tempClusEnergyTime,min,max);
                    SetXRange(tempClusEnergyTime,min,max);
					tempClusEnergyTime->GetXaxis()->SetTitle("#it{t}_{Cluster} (s)");
					tempClusEnergyTime->GetYaxis()->SetTitle("#it{E}_{Cluster} (GeV)");
					vecClusterEnergyTime[i].push_back(tempClusEnergyTime);
				}else cout << "Info: Object |ClusterTimeVsE_afterClusterQA| could not be found! Skipping..." << endl;
			}
	//---------
			if(doExtQA>0){
				TH2D* ClusEVsNCells = (TH2D*) CaloExtQAContainer->FindObject(Form("ClusterEnergyVsNCells_afterQA %s", fClusterCutSelection.Data()));
				if(ClusEVsNCells){
					TH2D* tempClusEVsNCells = new TH2D(*ClusEVsNCells);
					tempClusEVsNCells->SetName(Form("%s_ClusterEnergyVsNCells_%s",fRunNumber.Data(),fDataSet.Data()));
					tempClusEVsNCells->SetTitle(fRunNumber);
					tempClusEVsNCells->GetXaxis()->SetTitle("Cluster Energy (GeV)");
					tempClusEVsNCells->GetYaxis()->SetTitle("#it{N}_{Cells} per Cluster");
					vecClusterEVsNCells[i].push_back(tempClusEVsNCells);
				}else cout << "Info: Object |ClusterEnergyVsNCells_afterQA| could not be found! Skipping..." << endl;
			}
//---------
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
				vecClusterEnergy[i].push_back(tempClusEnergy);
			}else cout << "Info: Object |EnergyOfCluster_afterClusterQA| could not be found! Skipping..." << endl;
	//---------
			TH1D* ClusM02 = (TH1D*) CaloCutsContainer->FindObject(Form("M02_afterClusterQA %s", fClusterCutSelection.Data()));
			if(ClusM02){
				TH1D* tempClusM02 = new TH1D(*ClusM02);
				tempClusM02->Sumw2();
				tempClusM02->Scale(1 / nEvents);
				tempClusM02->GetXaxis()->SetTitle("#lambda_{0}^{2}");
				tempClusM02->GetYaxis()->SetTitle("#frac{1}{N_{Events}} #frac{dN}{d#lambda_{0}^{2}}");
				vecClusterM02[i].push_back(tempClusM02);
			}else cout << "Info: Object |M02_afterClusterQA| could not be found! Skipping..." << endl;
	//---------
			TH1D* ClusM20 = (TH1D*) CaloCutsContainer->FindObject(Form("M20_afterClusterQA %s", fClusterCutSelection.Data()));
			if(ClusM20){
				TH1D* tempClusM20 = new TH1D(*ClusM20);
				tempClusM20->Sumw2();
				tempClusM20->Scale(1 / nEvents);
				tempClusM20->GetXaxis()->SetTitle("#lambda_{1}^{2}");
				tempClusM20->GetYaxis()->SetTitle("#frac{1}{N_{Events}} #frac{dN}{d#lambda_{1}^{2}}");
				vecClusterM20[i].push_back(tempClusM20);
			}else cout << "Info: Object |M20_afterClusterQA| could not be found! Skipping..." << endl;
	//---------
			TH1D* ClusNCells = (TH1D*) CaloCutsContainer->FindObject(Form("NCellPerCluster_afterClusterQA %s", fClusterCutSelection.Data()));
			if(ClusNCells){
				TH1D* tempClusNCells = new TH1D(*ClusNCells);
				tempClusNCells->Sumw2();
				tempClusNCells->Scale(1 / nEvents);
				tempClusNCells->GetXaxis()->SetTitle("N_{Cells} per Cluster");
				tempClusNCells->GetYaxis()->SetTitle("#frac{1}{N_{Events}} #frac{dN}{dN_{Cells}}");
				vecClusterNCells[i].push_back(tempClusNCells);
			}else cout << "Info: Object |NCellPerCluster_afterClusterQA| could not be found! Skipping..." << endl;
	//---------
			TH1D* ClusDispersion = (TH1D*) CaloCutsContainer->FindObject(Form("Dispersion_afterClusterQA %s", fClusterCutSelection.Data()));
			if(ClusDispersion){
				TH1D* tempClusDispersion = new TH1D(*ClusDispersion);
				tempClusDispersion->Sumw2();
				tempClusDispersion->Scale(1 / nEvents);
				tempClusDispersion->GetXaxis()->SetTitle("Dispersion of Cluster");
				tempClusDispersion->GetYaxis()->SetTitle("#frac{1}{N_{Events}} #frac{dN}{dDisp}");
				vecClusterDispersion[i].push_back(tempClusDispersion);
			}else cout << "Info: Object |Dispersion_afterClusterQA| could not be found! Skipping..." << endl;
	//---------
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
                    vecClusterTime[i][iT].push_back(tempClusTime);
				}
			}else cout << "Info: Object |ClusterTimeVsE_afterClusterQA| could not be found! Skipping..." << endl;
    //---------
            TH1F* NLM = (TH1F*) CaloCutsContainer->FindObject(Form("NLM_afterClusterQA %s", fClusterCutSelection.Data()));
            if(NLM){
              hClusterMeanNLM[i]->SetBinContent(bin, NLM->GetMean());
              hClusterMeanNLM[i]->SetBinError(bin, NLM->GetMeanError());
              if(NLM->GetRMS()>0){
                hClusterRMSNLM[i]->SetBinContent(bin, NLM->GetRMS());
                hClusterRMSNLM[i]->SetBinError(bin, NLM->GetRMSError());
              }else{
                hClusterRMSNLM[i]->SetBinContent(bin, 1);
                hClusterRMSNLM[i]->SetBinError(bin, 0);
              }
            }else cout << "Info: Object |NLM_afterClusterQA| could not be found! Skipping Fill..." << endl;
	//---------
			if(isTrackMatching){
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

					vecClusterDeltaEta[i].push_back(tempClusdEta_matched);
					vecClusterDeltaPhi[i].push_back(tempClusdPhi_matched);
				}else cout << "Info: Object |dEtaVsdPhi_beforeClusterQA| or |dEtaVsdPhi_afterClusterQA| could not be found! Skipping..." << endl;
			}
	//---------
			if(isPCMCalo){
				for(Int_t iT=0; iT<2; iT++)
				{
					TH2D* ConvPhotonEtaPhi = (TH2D*) ESDContainer->FindObject(Form("ESD_Mother%sConvPhoton_Eta_Phi", sT[iT].Data()));
					if(ConvPhotonEtaPhi){
						if(iT==0){
							hClusterMeanConvPhotonPi0_Eta[i]->SetBinContent(bin, ConvPhotonEtaPhi->GetMean(2));
							hClusterMeanConvPhotonPi0_Eta[i]->SetBinError(bin, ConvPhotonEtaPhi->GetMeanError(2));
							hClusterMeanConvPhotonPi0_Phi[i]->SetBinContent(bin, ConvPhotonEtaPhi->GetMean(1));
							hClusterMeanConvPhotonPi0_Phi[i]->SetBinError(bin, ConvPhotonEtaPhi->GetMeanError(1));
							hClusterRMSConvPhotonPi0_Eta[i]->SetBinContent(bin, ConvPhotonEtaPhi->GetRMS(2));
							hClusterRMSConvPhotonPi0_Eta[i]->SetBinError(bin, ConvPhotonEtaPhi->GetRMSError(2));
							hClusterRMSConvPhotonPi0_Phi[i]->SetBinContent(bin, ConvPhotonEtaPhi->GetRMS(1));
							hClusterRMSConvPhotonPi0_Phi[i]->SetBinError(bin, ConvPhotonEtaPhi->GetRMSError(1));
						}else{
							hClusterMeanConvPhotonEta_Eta[i]->SetBinContent(bin, ConvPhotonEtaPhi->GetMean(2));
							hClusterMeanConvPhotonEta_Eta[i]->SetBinError(bin, ConvPhotonEtaPhi->GetMeanError(2));
							hClusterMeanConvPhotonEta_Phi[i]->SetBinContent(bin, ConvPhotonEtaPhi->GetMean(1));
							hClusterMeanConvPhotonEta_Phi[i]->SetBinError(bin, ConvPhotonEtaPhi->GetMeanError(1));
							hClusterRMSConvPhotonEta_Eta[i]->SetBinContent(bin, ConvPhotonEtaPhi->GetRMS(2));
							hClusterRMSConvPhotonEta_Eta[i]->SetBinError(bin, ConvPhotonEtaPhi->GetRMSError(2));
							hClusterRMSConvPhotonEta_Phi[i]->SetBinContent(bin, ConvPhotonEtaPhi->GetRMS(1));
							hClusterRMSConvPhotonEta_Phi[i]->SetBinError(bin, ConvPhotonEtaPhi->GetRMSError(1));
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
							vecClusterPi0ConvPhotonEta[i].push_back(tempConvPhotonEta);
							vecClusterPi0ConvPhotonPhi[i].push_back(tempConvPhotonPhi);
						}else{
							vecClusterEtaConvPhotonEta[i].push_back(tempConvPhotonEta);
							vecClusterEtaConvPhotonPhi[i].push_back(tempConvPhotonPhi);
						}
					}else cout << Form("Info: Object |ESD_Mother%sConvPhoton_Eta_Phi| could not be found! Skipping...", sT[iT].Data()) << endl;
				}
			}
	//---------
			if(doExtQA>0){
				TH2D* ClusterEVsModule = (TH2D*) CaloExtQAContainer->FindObject(Form("ClusterEnergyVsModule_afterClusterQA %s", fClusterCutSelection.Data()));
				TH2D* ModuleEnergyEVsModule = (TH2D*) CaloExtQAContainer->FindObject(Form("ModuleEnergyVsModule %s", fClusterCutSelection.Data()));
				TH2D* NCells100VsModule = (TH2D*) CaloExtQAContainer->FindObject(Form("NCellsAbove100VsModule %s", fClusterCutSelection.Data()));
				TH2D* NCells1500VsModule = (TH2D*) CaloExtQAContainer->FindObject(Form("NCellsAbove1500VsModule %s", fClusterCutSelection.Data()));

				for(Int_t iModule=0; iModule<nCaloModules; iModule++)
				{
					if(ClusterEVsModule){
						TH1D* ClusterEForModule = (TH1D*) ClusterEVsModule->ProjectionX(Form("projectClusterE_%i",iModule),iModule+1,iModule+1);
						TH1D* tempClusterEForModule = new TH1D(*ClusterEForModule);
						delete ClusterEForModule;
						tempClusterEForModule->GetXaxis()->SetTitle(Form("Cluster Energy of SM%i (GeV)",iModule));
						tempClusterEForModule->GetYaxis()->SetTitle("#frac{dE}{dN}");
                        tempClusterEForModule->SetTitle(Form("%s, nEvents: %.2e",fRunNumber.Data(),nEvents));
						tempClusterEForModule->Sumw2();
						tempClusterEForModule->Scale(1 / nEvents);
						vecClusterEVsModule[i][iModule].push_back(tempClusterEForModule);
					}else cout << "Info: Object |ClusterEnergyVsModule_afterClusterQA| could not be found! Skipping..." << endl;
					if(ModuleEnergyEVsModule){
						TH1D* ModuleEForModule = (TH1D*) ModuleEnergyEVsModule->ProjectionX(Form("projectModuleE_%i",iModule),iModule+1,iModule+1);
						TH1D* tempModuleEForModule = new TH1D(*ModuleEForModule);
						delete ModuleEForModule;
						tempModuleEForModule->GetXaxis()->SetTitle(Form("Total SuperModule Energy per Event for SM%i(GeV)",iModule));
						tempModuleEForModule->GetYaxis()->SetTitle("#frac{dE}{dN}");
                        tempModuleEForModule->SetTitle(Form("%s, nEvents: %.2e",fRunNumber.Data(),nEvents));
						tempModuleEForModule->Sumw2();
						tempModuleEForModule->Scale(1 / nEvents);
						vecModuleEVsModule[i][iModule].push_back(tempModuleEForModule);
					}else cout << "Info: Object |ModuleEnergyVsModule_afterClusterQA| could not be found! Skipping..." << endl;
					if(NCells100VsModule){
						TH1D* NCells100ForModule = (TH1D*) NCells100VsModule->ProjectionX(Form("projectNCells100_%i",iModule),iModule+1,iModule+1);
						TH1D* tempNCells100ForModule = new TH1D(*NCells100ForModule);
						delete NCells100ForModule;
						tempNCells100ForModule->GetXaxis()->SetTitle(Form("#it{N}_{Cells}>100 MeV in SM%i per Event",iModule));
						tempNCells100ForModule->GetYaxis()->SetTitle("#frac{d#it{N}_{Cells}>100 MeV}{dN}");
                        tempNCells100ForModule->SetTitle(Form("%s, nEvents: %.2e",fRunNumber.Data(),nEvents));
						tempNCells100ForModule->Sumw2();
						tempNCells100ForModule->Scale(1 / nEvents);
						vecClusterNCells100VsModule[i][iModule].push_back(tempNCells100ForModule);
					}else cout << "Info: Object |NCellsAbove100VsModule| could not be found! Skipping..." << endl;
					if(NCells1500VsModule){
						TH1D* NCells1500ForModule = (TH1D*) NCells1500VsModule->ProjectionX(Form("projectNCells1500_%i",iModule),iModule+1,iModule+1);
						TH1D* tempNCells1500ForModule = new TH1D(*NCells1500ForModule);
						delete NCells1500ForModule;
						tempNCells1500ForModule->GetXaxis()->SetTitle(Form("#it{N}_{Cells}>1500 MeV in SM%i per Event",iModule));
						tempNCells1500ForModule->GetYaxis()->SetTitle("#frac{d#it{N}_{Cells}>1500 MeV}{dN}");
                        tempNCells1500ForModule->SetTitle(Form("%s, nEvents: %.2e",fRunNumber.Data(),nEvents));
						tempNCells1500ForModule->Sumw2();
						tempNCells1500ForModule->Scale(1 / nEvents);
						vecClusterNCells1500VsModule[i][iModule].push_back(tempNCells1500ForModule);
					}else cout << "Info: Object |NCellsAbove1500VsModule| could not be found! Skipping..." << endl;
				}
		//---------
				TH1D* ClusIncludedCells = (TH1D*)CaloExtQAContainer->FindObject(Form("ClusterIncludedCells_afterClusterQA %s", fClusterCutSelection.Data()));
				if(ClusIncludedCells){
					TH1D* tempClusIncludedCells = new TH1D(*ClusIncludedCells);
					tempClusIncludedCells->SetName(Form("%s_ClusterIncludedCells_%s",fRunNumber.Data(),fDataSet.Data()));
					tempClusIncludedCells->SetTitle(fRunNumber);
					tempClusIncludedCells->GetXaxis()->SetTitle("Cell ID in accepted Clusters");
					tempClusIncludedCells->GetYaxis()->SetTitle("d#it{N}/d#it{CellID}");
					tempClusIncludedCells->Sumw2();
					vecClusterIncludedCells[i].push_back(tempClusIncludedCells);
				}else cout << "Info: Object |ClusterIncludedCells_afterClusterQA| could not be found! Skipping..." << endl;
		//---------
				if(doExtQA==2){
						TH1D* DeadCellsClusIncludedCellsBefore = (TH1D*)CaloExtQAContainer->FindObject(Form("ClusterIncludedCells_beforeClusterQA %s", fClusterCutSelection.Data()));
						if(DeadCellsClusIncludedCellsBefore){
							TH2D* tempDeadCellsRunwise = CompareDeadCellsRunwise(DeadCellsClusIncludedCellsBefore, nCaloCells, Form("%s, nEvents: %.2e",fRunNumber.Data(),nEvents));
							tempDeadCellsRunwise->SetName(Form("%s_ClusterIncludedCells_%s",fRunNumber.Data(),fDataSet.Data()));
							tempDeadCellsRunwise->GetYaxis()->SetTitle("Cell ID");
							vecClusterFiredCellIDs[i].push_back(tempDeadCellsRunwise);
						}else cout << "Info: Object |ClusterIncludedCells_beforeClusterQA| could not be found! Skipping..." << endl;
				}
		//---------
				TH1D* ClusIncludedCellsBefore = (TH1D*)CaloExtQAContainer->FindObject(Form("ClusterIncludedCells_beforeClusterQA %s", fClusterCutSelection.Data()));
				if(ClusIncludedCellsBefore){
					TH1D* tempClusIncludedCellBefore = new TH1D(*ClusIncludedCellsBefore);
//					if(doExtQA==2) CheckHotCellsRunwise(fLogRunwiseBadCells[i],tempClusIncludedCellBefore,nCaloCells,kFALSE);
//					if(doExtQA==2) CheckDeadCellsRunwise(fLogRunwiseDeadCells[i],tempClusIncludedCellBefore,nCaloCells);
					tempClusIncludedCellBefore->SetName(Form("%s_ClusterIncludedCells_beforeClusterQA_%s",fRunNumber.Data(),fDataSet.Data()));
					tempClusIncludedCellBefore->SetTitle(fRunNumber);
					tempClusIncludedCellBefore->GetXaxis()->SetTitle("Cell ID in all Clusters");
					tempClusIncludedCellBefore->GetYaxis()->SetTitle("d#it{N}/d#it{CellID}");
					tempClusIncludedCellBefore->Sumw2();
					vecClusterIncludedCellsBefore[i].push_back(tempClusIncludedCellBefore);
				}else cout << "Info: Object |ClusterIncludedCells_beforeClusterQA| could not be found! Skipping..." << endl;
				//---------
				TProfile* ClusBadCells = (TProfile*)CaloExtQAContainer->FindObject(Form("%s - Bad Channels",calo.Data()));
				if(ClusBadCells){
					TProfile* tempClusBadCells = new TProfile(*ClusBadCells);
					tempClusBadCells->SetName(Form("%s_BadCells_%s",fRunNumber.Data(),fDataSet.Data()));
					tempClusBadCells->SetTitle(fRunNumber);
					tempClusBadCells->GetXaxis()->SetTitle("Cell ID");
					tempClusBadCells->GetYaxis()->SetTitle("Cell ID Bad in Fraction of Events");
					vecClusterBadCells[i].push_back(tempClusBadCells);
				}else cout << Form("Info: Object |%s - Bad Channels| could not be found! Skipping...",calo.Data()) << endl;
		//---------
				TH1D* ClusEFracCells = (TH1D*)CaloExtQAContainer->FindObject(Form("ClusterEnergyFracCells_afterClusterQA %s", fClusterCutSelection.Data()));
				if(ClusEFracCells){
					TH1D* tempClusEFracCells = new TH1D(*ClusEFracCells);
					tempClusEFracCells->SetName(Form("%s_EFracCells_%s",fRunNumber.Data(),fDataSet.Data()));
					tempClusEFracCells->SetTitle(fRunNumber);
					tempClusEFracCells->GetXaxis()->SetTitle("Cell ID in accepted Clusters");
					tempClusEFracCells->GetYaxis()->SetTitle("#sum^{events} E-Frac of Cell");
					tempClusEFracCells->Sumw2();
					vecClusterEFracCells[i].push_back(tempClusEFracCells);
				}else cout << "Info: Object |ClusterEnergyFracCells_afterClusterQA| could not be found! Skipping..." << endl;
		//---------
				TH1D* ClusEFracCellsBefore = (TH1D*)CaloExtQAContainer->FindObject(Form("ClusterEnergyFracCells_beforeClusterQA %s", fClusterCutSelection.Data()));
				if(ClusEFracCellsBefore){
					TH1D* tempClusEFracCellBefore = new TH1D(*ClusEFracCellsBefore);
					if(doExtQA==2) CheckHotAndColdCellsEFracRunwise(fLogRunwiseHotCells[i],fLogRunwiseColdCells[i],tempClusEFracCellBefore,ClusBadCells,nCaloCells,kFALSE);
					tempClusEFracCellBefore->SetName(Form("%s_EFracCells_beforeQA_%s",fRunNumber.Data(),fDataSet.Data()));
					tempClusEFracCellBefore->SetTitle(fRunNumber);
					tempClusEFracCellBefore->GetXaxis()->SetTitle("Cell ID in all Clusters");
					tempClusEFracCellBefore->GetYaxis()->SetTitle("#sum^{events} E-Frac of Cell");
					tempClusEFracCellBefore->Sumw2();
					vecClusterEFracCellsBefore[i].push_back(tempClusEFracCellBefore);
				}else cout << "Info: Object |ClusterEnergyFracCells_beforeClusterQA| could not be found! Skipping..." << endl;
		//---------
				if(doExtQA==2){
					TH2D* fHistCellEnergyVsCellID = (TH2D*)CaloExtQAContainer->FindObject(Form("CellEnergyVsCellID %s", fClusterCutSelection.Data()));
					TH2D* fHistCellTimeVsCellID = (TH2D*)CaloExtQAContainer->FindObject(Form("CellTimeVsCellID %s", fClusterCutSelection.Data()));
					if(fHistCellEnergyVsCellID && fHistCellTimeVsCellID){
						TH2D** tempClusCell = PlotCellMeanVsSigmaForRunwise(nCaloCells,fHistCellEnergyVsCellID,fHistCellTimeVsCellID,
													"Mean Cell Energy (GeV)","#sigma_{Cell Energy} (GeV)","Mean Cell Time (s)","#sigma_{Cell Time} (s)",(i>=nData));
						tempClusCell[0]->SetName(Form("%s_ClusterEMeanVsSigma_%s",fRunNumber.Data(),fDataSet.Data()));
						tempClusCell[0]->SetTitle(fRunNumber);
						vecClusterEnergyMeanSigma[i].push_back(tempClusCell[0]);
						tempClusCell[1]->SetName(Form("%s_ClusterTimeMeanVsSigma_%s",fRunNumber.Data(),fDataSet.Data()));
						tempClusCell[1]->SetTitle(fRunNumber);
						vecClusterTimeMeanSigma[i].push_back(tempClusCell[1]);
					}else cout << "Info: Object |CellEnergyVsCellID or CellTimeVsCellID| could not be found! Skipping..." << endl;

					for(Int_t j=0; j<(Int_t) vecBadCells[i].size(); j++)
					{
						TH2D* fHistCellEnergyVsCellID = (TH2D*)CaloExtQAContainer->FindObject(Form("CellEnergyVsCellID %s", fClusterCutSelection.Data()));
						TH2D* fHistCellTimeVsCellID = (TH2D*)CaloExtQAContainer->FindObject(Form("CellTimeVsCellID %s", fClusterCutSelection.Data()));
						if(fHistCellEnergyVsCellID && fHistCellTimeVsCellID){
							TH1D* tempEnergy;
							TH1D* tempTime;
							tempEnergy = (TH1D*) fHistCellEnergyVsCellID->ProjectionX("energy",((TString)vecBadCells[i].at(j)).Atoi(),((TString)vecBadCells[i].at(j)).Atoi());
							tempTime = (TH1D*) fHistCellTimeVsCellID->ProjectionX("time",((TString)vecBadCells[i].at(j)).Atoi(),((TString)vecBadCells[i].at(j)).Atoi());

							TH1D* tempEnergyCell = new TH1D(*tempEnergy);
							delete tempEnergy;
							tempEnergyCell->GetXaxis()->SetTitle(Form("Cell Energy of ID %i (GeV)",((TString)vecBadCells[i].at(j)).Atoi()));
							tempEnergyCell->GetYaxis()->SetTitle("#frac{dE}{dN}");
							tempEnergyCell->Sumw2();
							tempEnergyCell->Scale(1 / nEvents);
							vecBadCellsEnergy[i].push_back(tempEnergyCell);

							TH1D* tempTimeCell = new TH1D(*tempTime);
							delete tempTime;
							tempTimeCell->GetXaxis()->SetTitle(Form("Cell Time of ID %i (s)",((TString)vecBadCells[i].at(j)).Atoi()));
							tempTimeCell->GetYaxis()->SetTitle("#frac{dTime}{dN}");
							tempTimeCell->Sumw2();
							tempTimeCell->Scale(1 / nEvents);
							vecBadCellsTime[i].push_back(tempTimeCell);
						}
					}
				}
			}
	//---------
			delete TopDir;

			RootFile->Close();
			delete RootFile;
		}

	}

//****************************** Drawing Histograms ************************************************
	cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
	cout << "Drawing Histograms" << endl;

    TCanvas* canvas = new TCanvas("canvas","",10,10,750,500);  // gives the page size
    Double_t leftMar = 0.09; Double_t rightMar = 0.025; Double_t topMargin = 0.04; Double_t bottomMargin = 0.09;
	DrawGammaCanvasSettings(canvas, leftMar, rightMar, topMargin, bottomMargin);

	if(doHistsForEverySet)
	{
		TBox *boxLabel = new TBox(1.37,0.7,1.78,0.83);
		boxLabel->SetFillStyle(0);boxLabel->SetFillColor(0);boxLabel->SetLineColor(1);boxLabel->SetLineWidth(0.6);
		TBox *boxLabel2 = new TBox(-0.4,51,6.5,56.5);
		boxLabel2->SetFillStyle(0);boxLabel2->SetFillColor(0);boxLabel2->SetLineColor(1);boxLabel2->SetLineWidth(0.6);

		TString outputDirDataSet;

		for(Int_t i=0; i<nSets; i++)
		{
			cout << "DataSet: " << DataSets[i].Data() << endl;
			outputDirDataSet = Form("%s/%s",outputDir.Data(), DataSets[i].Data());

			if(useDataRunListForMC && i>=nData) {
				outputDirDataSet = Form("%s/%s-%s", outputDir.Data(), DataSets[i].Data(),DataSets[0].Data());
				cout << "Switch useDataRunListForMC is true, output to: " << outputDirDataSet.Data() << endl;
			}
			gSystem->Exec("mkdir -p "+outputDirDataSet);
			gSystem->Exec("mkdir -p "+outputDirDataSet+"/ExtQA");
			gSystem->Exec("mkdir -p "+outputDirDataSet+"/ExtQA/EFrac");
			if(doExtQA>0) gSystem->Exec("mkdir -p "+outputDirDataSet+"/ExtQA/IncludedCells");
			if(doExtQA==2) gSystem->Exec("mkdir -p "+outputDirDataSet+"/ExtQA/MissingCells");
			gSystem->Exec("mkdir -p "+outputDirDataSet+"/ExtQA/ConvPhotonEtaPhi");
			if(doExtQA==2) gSystem->Exec("mkdir -p "+outputDirDataSet+"/BadCells");

			TString fTrigger = "";
			if(i<nData){
				TString fTriggerCut = fEventCutSelection(3,2);
				fTrigger = AnalyseSpecialTriggerCut(fTriggerCut.Atoi(), plotDataSets[i]);
				cout << "Trigger: '" << fTrigger.Data() << "'" << endl;
				if(fTrigger.Contains("not defined")){
					fTrigger = "";
					cout << "INFO: Trigger cut not defined!" << endl;
				}
			}

			DrawVectorOverviewTH2D(canvas, vecClusterEtaPhi[i], "hClusterEtaPhi_scaledNEventsAndMean", outputDirDataSet, suffix,
										0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, boxLabel, kFALSE, kFALSE);

			TGaxis::SetExponentOffset(0, -0.1, "x");
			if(i<nData) DrawVectorOverviewTH2D(canvas, vecClusterEnergyTime[i], "ExtQA/hClusterEnergyTime", outputDirDataSet, suffix,
										0.13, 0.15, 0.1, 0.14, 0.8, 0.8, 0.12, 0.93, 0x0, kTRUE, kTRUE);
			TGaxis::SetExponentOffset(0, 0, "x");

			if(doExtQA>0){
				DrawVectorOverviewTH2D(canvas, vecClusterEVsNCells[i], "ExtQA/hClusterEVsNCells", outputDirDataSet, suffix,
										0.13, 0.15, 0.1, 0.14, 0.8, 0.8, 0.12, 0.93, boxLabel2, kFALSE, kTRUE);
				DrawVectorOverviewTH1D(canvas, vecClusterIncludedCells[i], "ExtQA/IncludedCells/hClusterIncludedCells", outputDirDataSet, suffix,
										0.13, 0.15, 0.1, 0.14, 0.8, 0.9, 0.12, 0.93, 0x0, kFALSE, kTRUE);
				DrawVectorOverviewTH1D(canvas, vecClusterIncludedCellsBefore[i], "ExtQA/IncludedCells/hClusterIncludedCellsBefore", outputDirDataSet, suffix,
										0.13, 0.15, 0.1, 0.14, 0.8, 0.9, 0.12, 0.93, 0x0, kFALSE, kTRUE);
			}

			if(i < nData && doExtQA==2){
				for(Int_t iMC=nData; iMC<nSets; iMC++){
					if(vecClusterFiredCellIDs[i].size() == vecClusterFiredCellIDs[iMC].size() && vecClusterMissingCellIDs[i].size() == vecClusterFiredCellIDs[iMC].size()){
						DrawVectorOverviewMissingCells(canvas, vecClusterFiredCellIDs[i], vecClusterFiredCellIDs[iMC], vecClusterMissingCellIDs[i],
															Form("ExtQA/MissingCells/hCellsMissingData_%s",DataSets[iMC].Data()), outputDirDataSet, suffix, DataSets[iMC]);
					}
				}
			}

			DrawVectorOverviewTH1D(canvas, vecClusterEFracCellsBefore[i], "ExtQA/EFrac/hClusterEFracCellsBefore", outputDirDataSet, suffix,
									0.13, 0.15, 0.1, 0.14, 0.8, 0.9, 0.12, 0.93, 0x0, kFALSE, kTRUE);
			DrawVectorOverviewTH1D(canvas, vecClusterEFracCells[i], "ExtQA/EFrac/hClusterEFracCells", outputDirDataSet, suffix,
									0.13, 0.15, 0.1, 0.14, 0.8, 0.9, 0.12, 0.93, 0x0, kFALSE, kTRUE);

			if(isPCMCalo){
				DrawVectorOverviewTH1D(canvas, vecClusterPi0ConvPhotonEta[i], "ExtQA/ConvPhotonEtaPhi/hClusterPi0ConvPhotonEta_Runwise", outputDirDataSet, suffix,
										0.13, 0.15, 0.1, 0.14, 0.8, 0.9, 0.2, 0.93, 0x0, kFALSE, kTRUE);
				DrawVectorOverviewTH1D(canvas, vecClusterPi0ConvPhotonPhi[i], "ExtQA/ConvPhotonEtaPhi/hClusterPi0ConvPhotonPhi_Runwise", outputDirDataSet, suffix,
										0.13, 0.15, 0.1, 0.14, 0.8, 0.9, 0.2, 0.93, 0x0, kFALSE, kTRUE);
				DrawVectorOverviewTH1D(canvas, vecClusterEtaConvPhotonEta[i], "ExtQA/ConvPhotonEtaPhi/hClusterEtaConvPhotonEta_Runwise", outputDirDataSet, suffix,
										0.13, 0.15, 0.1, 0.14, 0.8, 0.9, 0.2, 0.93, 0x0, kFALSE, kTRUE);
				DrawVectorOverviewTH1D(canvas, vecClusterEtaConvPhotonPhi[i], "ExtQA/ConvPhotonEtaPhi/hClusterEtaConvPhotonPhi_Runwise", outputDirDataSet, suffix,
										0.13, 0.15, 0.1, 0.14, 0.8, 0.9, 0.2, 0.93, 0x0, kFALSE, kTRUE);
			}

		//--------
			for(Int_t h=0; h<(Int_t) vecHistos[i].size(); h++)
			{
				(((TH1D*) vecHistos[i].at(h)))->SetTitle("");
                if( ((TString)vecHistosName.at(h)).CompareTo("hNEvents")==0 ) {
					AdjustHistRange(((TH1D*) vecHistos[i].at(h)),10,10,kTRUE);
					((TH1D*) vecHistos[i].at(h))->Draw("p");
				}else{
					AdjustHistRange(((TH1D*) vecHistos[i].at(h)),1.1,1.1,kTRUE);
					((TH1D*) vecHistos[i].at(h))->Draw("px0e1");
				}

                if(doTrigger && i<nData){
                  PutProcessLabelAndEnergyOnPlot(xPosLabel, 0.92, 0.03, fCollisionSystem.Data(), plotDataSets[i].Data(), fTrigger.Data());
                  PutProcessLabelAndEnergyOnPlot(xPosLabel, 0.81, 0.03, Form("%s clusters", calo.Data()), "", "");
                }else{
                  PutProcessLabelAndEnergyOnPlot(xPosLabel, 0.92, 0.03, fCollisionSystem.Data(), plotDataSets[i].Data(), Form("%s clusters", calo.Data()));
                }

                if( ((TString)vecHistosName.at(h)).CompareTo("hNEvents")==0 ) SaveCanvas(canvas, Form("%s/%s.%s", outputDirDataSet.Data(), vecHistosName.at(h).Data(), suffix.Data()), kFALSE, kTRUE);
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

		for(Int_t i=0; i<nSets; i++)
		{
			vecRuns.clear();
			fDataSet = vecDataSet.at(i);
			fileRuns = Form("runNumbers%s.txt", fDataSet.Data());
			outputDirDataSet = Form("%s/%s",outputDir.Data(), DataSets[i].Data());

			if(useDataRunListForMC && i>=nData) {
				fileRuns = Form("runNumbers%s-%s.txt", vecDataSet.at(i).Data(),vecDataSet.at(0).Data());
				outputDirDataSet = Form("%s/%s-%s", outputDir.Data(), DataSets[i].Data(),DataSets[0].Data());
				cout << "Switch useDataRunListForMC is true, reading runs from: " << fileRuns.Data() << endl;
			}
			if(!readin(fileRuns, vecRuns, kFALSE)) {cout << "ERROR, no Run Numbers could be found! Returning..." << endl; return;}

			for(Int_t iRun=0; iRun<(Int_t)vecMissingRuns[i].size(); iRun++){
				vecRuns.erase(std::remove(vecRuns.begin(), vecRuns.end(), vecMissingRuns[i].at(iRun)), vecRuns.end());
			}

			Int_t NColumns = ((Int_t) vecRuns.size() / 31 ) + 1;

			//extending nHist for testing
			//for(Int_t i=123567; i<123667; i++) {globalRuns.push_back(Form("%i",i)); vecClusterEnergy[0].push_back(vecClusterEnergy[0].at(0));}

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
            if(doTrigger && i<nData){
				TString fTriggerCut = fEventCutSelection(3,2);
				fTrigger = AnalyseSpecialTriggerCut(fTriggerCut.Atoi(), plotDataSets[i]);
				if(fTrigger.Contains("not defined")) fTrigger = "";
			}
		//---------
            DrawVectorRunwiseTH1D(	canvasRunwise, legendRuns, vecClusterEnergy[i], vecRuns, 5, 5, kFALSE, addRight, 0.8, 0.92, 0.03, 0.8, 0.81, 0.03,
									doTrigger, fTrigger, (Bool_t)(i<nData), outputDirDataSet, "ClusterEnergy_Runwise", plotDataSets[i],kFALSE,
                                    fCollisionSystem, Form("%s clusters", calo.Data()), suffix, kTRUE, kTRUE, kFALSE);
		//---------
            DrawVectorRunwiseTH1D(	canvasRunwise, legendRuns, vecClusterM02[i], vecRuns, 5, 5, kTRUE, addRight, 0.8, 0.92, 0.03, 0.8, 0.81, 0.03,
									doTrigger, fTrigger, (Bool_t)(i<nData), outputDirDataSet, "ClusterM02_Runwise", plotDataSets[i], kFALSE,
                                    fCollisionSystem, Form("%s clusters", calo.Data()), suffix, kFALSE, kTRUE, kFALSE);
		//---------
            DrawVectorRunwiseTH1D(	canvasRunwise, legendRuns, vecClusterM20[i], vecRuns, 5, 5, kTRUE, addRight, 0.8, 0.92, 0.03, 0.8, 0.81, 0.03,
									doTrigger, fTrigger, (Bool_t)(i<nData), outputDirDataSet, "ClusterM20_Runwise", plotDataSets[i], kFALSE,
                                    fCollisionSystem, Form("%s clusters", calo.Data()), suffix, kFALSE, kTRUE, kFALSE);
		//---------
            DrawVectorRunwiseTH1D(	canvasRunwise, legendRuns, vecClusterNCells[i], vecRuns, 5, 5, kTRUE, addRight, 0.8, 0.92, 0.03, 0.8, 0.81, 0.03,
									doTrigger, fTrigger, (Bool_t)(i<nData), outputDirDataSet, "ClusterNCells_Runwise", plotDataSets[i], kFALSE,
                                    fCollisionSystem, Form("%s clusters", calo.Data()), suffix, kTRUE, kTRUE, kFALSE);
		//---------
            DrawVectorRunwiseTH1D(	canvasRunwise, legendRuns, vecClusterDispersion[i], vecRuns, 5, 5, kTRUE, addRight, 0.8, 0.92, 0.03, 0.8, 0.81, 0.03,
									doTrigger, fTrigger, (Bool_t)(i<nData), outputDirDataSet, "ClusterDispersion_Runwise", plotDataSets[i], kFALSE,
                                    fCollisionSystem, Form("%s clusters", calo.Data()), suffix, kFALSE, kTRUE, kFALSE);
		//---------
			TGaxis::SetExponentOffset(0.5, 0, "x");
            DrawVectorRunwiseTH1D(	canvasRunwise, legendRuns, vecClusterTime[i][0], vecRuns, 5, 5, kTRUE, addRight, 0.8, 0.92, 0.03, 0.8, 0.81, 0.03,
									doTrigger, fTrigger, (Bool_t)(i<nData), outputDirDataSet, "ClusterTime_Runwise", plotDataSets[i], kTRUE,
                                    fCollisionSystem, Form("%s clusters", calo.Data()), suffix, kFALSE, kTRUE, kFALSE);
            for(Int_t iT=1; iT<5; iT++){
                DrawVectorRunwiseTH1D(	canvasRunwise, legendRuns, vecClusterTime[i][iT], vecRuns, 5, 5, kTRUE, addRight, 0.8, 0.92, 0.03, 0.8, 0.81, 0.03,
                                        doTrigger, fTrigger, (Bool_t)(i<nData), outputDirDataSet, Form("ClusterTime_Runwise_%.01f-%.01f",minT_Energy[iT], maxT_Energy[iT]), plotDataSets[i], kTRUE,
                                        fCollisionSystem, Form("%s clusters", calo.Data()), suffix, kFALSE, kTRUE, kFALSE);
            }
            TGaxis::SetExponentOffset(0, 0, "x");
		//---------
			if(isTrackMatching){
				TGaxis::SetExponentOffset(0.013, -0.0285, "x");
                DrawVectorRunwiseTH1D(	canvasRunwise, legendRuns, vecClusterDeltaEta[i], vecRuns, 5, 5, kTRUE, addRight, 0.8, 0.92, 0.03, 0.8, 0.81, 0.03,
										doTrigger, fTrigger, (Bool_t)(i<nData), outputDirDataSet, "ClusterDeltaEta_Runwise", plotDataSets[i], kTRUE,
										fCollisionSystem, Form("%s clusters", calo.Data()), suffix, kFALSE, kTRUE, kFALSE);
                DrawVectorRunwiseTH1D(	canvasRunwise, legendRuns, vecClusterDeltaPhi[i], vecRuns, 5, 5, kTRUE, addRight, 0.8, 0.92, 0.03, 0.8, 0.81, 0.03,
										doTrigger, fTrigger, (Bool_t)(i<nData), outputDirDataSet, "ClusterDeltaPhi_Runwise", plotDataSets[i], kTRUE,
										fCollisionSystem, Form("%s clusters", calo.Data()), suffix, kFALSE, kTRUE, kFALSE);
				TGaxis::SetExponentOffset(0, 0, "x");
            }
		//--------
			if(doExtQA>0){
				for(Int_t iModule=0; iModule<nCaloModules; iModule++){
                    DrawVectorRunwiseTH1D(	canvasRunwise, legendRuns, vecClusterEVsModule[i][iModule], vecRuns, 5, 5, kTRUE, addRight, 0.8, 0.92, 0.03, 0.8, 0.81, 0.03,
											doTrigger, fTrigger, (Bool_t)(i<nData), outputDirDataSet, Form("ExtQA/ClusterEVsModule%i_Runwise",iModule), plotDataSets[i], kFALSE,
											fCollisionSystem, Form("%s clusters", calo.Data()), suffix, kTRUE, kTRUE, kFALSE);
                    DrawVectorRunwiseTH1D(	canvasRunwise, legendRuns, vecModuleEVsModule[i][iModule], vecRuns, 5, 5, kTRUE, addRight, 0.8, 0.92, 0.03, 0.8, 0.81, 0.03,
											doTrigger, fTrigger, (Bool_t)(i<nData), outputDirDataSet, Form("ExtQA/ModuleEVsModule%i_Runwise",iModule), plotDataSets[i], kFALSE,
											fCollisionSystem, Form("%s clusters", calo.Data()), suffix, kTRUE, kTRUE, kFALSE);
                    DrawVectorRunwiseTH1D(	canvasRunwise, legendRuns, vecClusterNCells100VsModule[i][iModule], vecRuns, 5, 5, kTRUE, addRight, 0.8, 0.92, 0.03, 0.8, 0.81, 0.03,
											doTrigger, fTrigger, (Bool_t)(i<nData), outputDirDataSet, Form("ExtQA/NCells100VsModule%i_Runwise",iModule), plotDataSets[i], kFALSE,
											fCollisionSystem, Form("%s clusters", calo.Data()), suffix, kFALSE, kTRUE, kFALSE);
                    DrawVectorRunwiseTH1D(	canvasRunwise, legendRuns, vecClusterNCells1500VsModule[i][iModule], vecRuns, 5, 5, kTRUE, addRight, 0.8, 0.92, 0.03, 0.8, 0.81, 0.03,
											doTrigger, fTrigger, (Bool_t)(i<nData), outputDirDataSet, Form("ExtQA/NCells1500VsModule%i_Runwise",iModule), plotDataSets[i], kFALSE,
											fCollisionSystem, Form("%s clusters", calo.Data()), suffix, kFALSE, kTRUE, kFALSE);
				}
            }

            if(doExtQA>1){
				for(Int_t iBad=0; iBad<(Int_t)vecBadCells[i].size(); iBad++){
                    DrawVectorRunwiseBadCells(canvasRunwise, legendRuns, vecBadCellsEnergy[i], (Int_t)vecBadCells[i].size(), iBad, vecRuns, 5, 5, kTRUE, addRight, 0.8, 0.92, 0.03, 0.8, 0.81, 0.03,
											   doTrigger, fTrigger, (Bool_t)(i<nData), outputDirDataSet, Form("BadCells/Cell%i_Energy",((TString)vecBadCells[i].at(iBad)).Atoi()), plotDataSets[i], kTRUE,
											   fCollisionSystem, Form("%s clusters", calo.Data()), suffix, kFALSE, kTRUE, kFALSE);
                    DrawVectorRunwiseBadCells(canvasRunwise, legendRuns, vecBadCellsTime[i], (Int_t)vecBadCells[i].size(), iBad, vecRuns, 5, 5, kTRUE, addRight, 0.8, 0.92, 0.03, 0.8, 0.81, 0.03,
											   doTrigger, fTrigger, (Bool_t)(i<nData), outputDirDataSet, Form("BadCells/Cell%i_Time",((TString)vecBadCells[i].at(iBad)).Atoi()), plotDataSets[i], kTRUE,
											   fCollisionSystem, Form("%s clusters", calo.Data()), suffix, kFALSE, kTRUE, kFALSE);
				}
            }
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
			TString fTriggerCut = fEventCutSelection(3,2);
			for(Int_t i=0; i<nSets; i++)
			{
				if(doTrigger && i<nData && i==0){
					fTrigger = AnalyseSpecialTriggerCut(fTriggerCut.Atoi(), plotDataSets[i]);
					if(fTrigger.Contains("not defined")) fTrigger = "";
				}
				legend->AddEntry(((TH1D*) vecHistos[i].at(h)),plotDataSets[i].Data(),"p");
			}
            if( ((TString)vecHistosName.at(h)).CompareTo("hNEvents")==0 ||
                ((TString)vecHistosName.at(h)).CompareTo("hClusterTime-Mean")==0 )
                 AdjustHistRange(vecHistos,10,10,h,nSets,kTRUE);
			else AdjustHistRange(vecHistos,1.1,1.1,h,nSets,kTRUE);
			for(Int_t i=nSets-1; i>=0; i--)
			{
				TString draw;
				if(h==0) draw = (i==nSets-1)?"p":"p, same";
				else draw = (i==nSets-1)?"px0e1":"px0e1, same";
				((TH1D*) vecHistos[i].at(h))->SetTitle("");
				((TH1D*) vecHistos[i].at(h))->Draw(draw.Data());
			}
			legend->Draw();

            if(doTrigger) PutProcessLabelAndEnergyOnPlot(xPosLabel, 0.92, 0.03, fCollisionSystem.Data(), fTrigger.Data(), Form("%s clusters", calo.Data()));
            else PutProcessLabelAndEnergyOnPlot(xPosLabel, 0.92, 0.03, fCollisionSystem.Data(), Form("%s clusters", calo.Data()), "");

            if(canvas->GetTopMargin()!=0.06) canvas->SetTopMargin(0.06);
			if(useDataRunListForMC && !addSubFolder){
               if( ((TString)vecHistosName.at(h)).CompareTo("hNEvents")==0 ||
                   ((TString)vecHistosName.at(h)).CompareTo("hClusterTime-Mean")==0 )
                     SaveCanvas(canvas, Form("%s/%s/%s.%s", outputDir.Data(), DataSets[0].Data(),vecHistosName.at(h).Data(),suffix.Data()), kFALSE, kTRUE);
				else SaveCanvas(canvas, Form("%s/%s/%s.%s", outputDir.Data(), DataSets[0].Data(),vecHistosName.at(h).Data(),suffix.Data()));
			}else{
              if( ((TString)vecHistosName.at(h)).CompareTo("hNEvents")==0 ||
                  ((TString)vecHistosName.at(h)).CompareTo("hClusterTime-Mean")==0 )
                     SaveCanvas(canvas, Form("%s/%s.%s", outputDir.Data(), vecHistosName.at(h).Data(),suffix.Data()), kFALSE, kTRUE);
				else SaveCanvas(canvas, Form("%s/%s.%s", outputDir.Data(), vecHistosName.at(h).Data(),suffix.Data()));
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
					outputDirDataSet = Form("%s/%s",outputDir.Data(), DataSets[i].Data());
					gSystem->Exec("mkdir -p "+outputDirDataSet+"/TrendingRatios");

                    if(doTrigger){
                      TString fTrigger = "";
                      TString fTriggerCut = fEventCutSelection(3,2);
                      fTrigger = AnalyseSpecialTriggerCut(fTriggerCut.Atoi(), plotDataSets[i]);
                      if(fTrigger.Contains("not defined")) fTrigger = "";

                      PutProcessLabelAndEnergyOnPlot(xPosLabel, 0.92, 0.03, fCollisionSystem.Data(), fTrigger.Data(), Form("%s clusters", calo.Data()));
                    }else PutProcessLabelAndEnergyOnPlot(xPosLabel, 0.92, 0.03, fCollisionSystem.Data(), Form("%s clusters", calo.Data()), "");

					SaveCanvas(canvas, Form("%s/TrendingRatios/%s.%s", outputDirDataSet.Data(),Form("%s",((TH1D*) vecHistos[i].at(h))->GetName()),suffix.Data()));
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

	for(Int_t i=0; i<nSets; i++)
	{
		fDataSet = vecDataSet.at(i);

		if(useDataRunListForMC && i>=nData) nameOutput = Form("%s/%s/ClusterQA/%s-%s_ClusterQARunwise.root",fCutSelection.Data(),fEnergyFlag.Data(),fDataSet.Data(),vecDataSet.at(0).Data());
		else nameOutput = Form("%s/%s/ClusterQA/%s_ClusterQARunwise.root",fCutSelection.Data(),fEnergyFlag.Data(),fDataSet.Data());

		TFile* fOutput = new TFile(nameOutput,"RECREATE");
		cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
		cout << "Output file: " << nameOutput << endl;
		cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
		fLog[i] << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
		fLog[i] << "Output file: " << nameOutput << endl;
		fLog[i] << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

		for(Int_t h=0; h<(Int_t) vecHistos[i].size(); h++) WriteHistogram(((TH1D*) vecHistos[i].at(h)));

		WriteHistogramTH2DVec(fOutput,vecClusterEtaPhi[i],"EtaVsPhi");
		if(i<nData) WriteHistogramTH2DVec(fOutput,vecClusterEnergyTime[i],"EnergyVsTime");

		if(doExtQA>0){
			WriteHistogramTH2DVec(fOutput,vecClusterEVsNCells[i],"ClusterEVsNCells");
			WriteHistogramTH1DVec(fOutput,vecClusterIncludedCells[i],"IncludedCells");
			WriteHistogramTH1DVec(fOutput,vecClusterIncludedCellsBefore[i],"IncludedCells_beforeQA");
			WriteHistogramTProfileVec(fOutput,vecClusterBadCells[i],"BadCells");
			if(i<nData) WriteHistogramTH2DVec(fOutput,vecClusterMissingCellIDs[i],"MissingCells");
			WriteHistogramTH1DVec(fOutput,vecClusterEFracCells[i],"EFracCells");
			WriteHistogramTH1DVec(fOutput,vecClusterEFracCellsBefore[i],"EFracCells_beforeQA");
			if(doExtQA==2){
				WriteHistogramTH2DVec(fOutput,vecClusterEnergyMeanSigma[i],"ClusterEMeanVsSigma");
				WriteHistogramTH2DVec(fOutput,vecClusterTimeMeanSigma[i],"ClusterTimeMeanVsSigma");
			}
		}

		DeleteVecTH2D(vecClusterFiredCellIDs[i]);

		DeleteVecTH1D(vecClusterEnergy[i]);
		DeleteVecTH1D(vecClusterM02[i]);
		DeleteVecTH1D(vecClusterM20[i]);
		DeleteVecTH1D(vecClusterNCells[i]);
		DeleteVecTH1D(vecClusterDispersion[i]);
		for(Int_t iT=0; iT<5; iT++) DeleteVecTH1D(vecClusterTime[i][iT]);
		DeleteVecTH1D(vecClusterDeltaEta[i]);
		DeleteVecTH1D(vecClusterDeltaPhi[i]);
		DeleteVecTH1D(vecClusterPi0ConvPhotonEta[i]);
		DeleteVecTH1D(vecClusterPi0ConvPhotonPhi[i]);
		DeleteVecTH1D(vecClusterEtaConvPhotonEta[i]);
		DeleteVecTH1D(vecClusterEtaConvPhotonPhi[i]);

		DeleteVecTH1D(vecBadCellsEnergy[i]);
		DeleteVecTH1D(vecBadCellsTime[i]);

		for(Int_t iCM=0; iCM<nCaloModules; iCM++){
            WriteHistogramTH1DVec(fOutput,vecClusterEVsModule[i][iCM],Form("Mod%02i/ClusterE",iCM));
            WriteHistogramTH1DVec(fOutput,vecModuleEVsModule[i][iCM],Form("Mod%02i/ModuleE",iCM));
            WriteHistogramTH1DVec(fOutput,vecClusterNCells100VsModule[i][iCM],Form("Mod%02i/NCells100",iCM));
            WriteHistogramTH1DVec(fOutput,vecClusterNCells1500VsModule[i][iCM],Form("Mod%02i/NCells1500",iCM));
		}

		fOutput->Write();
		fOutput->Close();
		delete fOutput;
		fLog[i].close();
		if(doExtQA==2) fLogRunwiseHotCells[i].close();
		if(doExtQA==2) fLogRunwiseColdCells[i].close();
	}

	if(doExtQA==2){
        std::vector< std::vector<TString> > vec(nSets);
        std::vector< std::map <Int_t,Int_t> > ma(nSets);
		TH1D* nCellsRun[nSets];
        Int_t* nRuns = new Int_t[nSets];
        TCanvas* cvsRun = new TCanvas("canvas","",10,10,750,500);  // gives the page size
        DrawGammaCanvasSettings(cvsRun, leftMar, rightMar, 0.06, bottomMargin);

		for(Int_t iStr=0; iStr<nSets; iStr++){
			nRuns[iStr]=0;
			fLogRunwiseHotCells[iStr].open(Form("%s/HotCellsRunwise-%s.log",outputDir.Data(),DataSets[iStr].Data()), ios::in);
			if(fLogRunwiseHotCells[iStr].good())
			{
				fLogRunwiseHotCells[iStr].seekg(0L, ios::beg);
				TString fVar;
				while(!fLogRunwiseHotCells[iStr].eof())
				{
					fLogRunwiseHotCells[iStr] >> fVar;
					if(fVar.BeginsWith("Run-")) {
						nRuns[iStr]++;
						vec[iStr].push_back(fVar);
					}
					else if(fVar.BeginsWith("NoNoisy")||fVar.BeginsWith("NotEnough")){
						continue;
					}
					else if(fVar.Sizeof()>1) {
						TObjArray *rNumber = fVar.Tokenize("-");
						TObjString* rString = (TObjString*)rNumber->At(0);
						TString vecString = rString->GetString();
						vec[iStr].push_back(vecString);
						if( ma[iStr].find(vecString.Atoi()) != ma[iStr].end() ) ma[iStr][vecString.Atoi()] += 1;
						else ma[iStr][vecString.Atoi()] = 1;
					}
				}
			}

			if((Int_t)ma[iStr].size()==0){
				cout << "No Bad Cells found for: " << DataSets[iStr].Data() << endl;
				continue;
			}

			Int_t plotting = 0;
			Int_t nPlot = 100;
			Int_t nPart = 0;
			map<Int_t, Int_t>::iterator it = ma[iStr].begin();
			do
			{
				plotting+=nPlot;
				if(plotting>(Int_t)ma[iStr].size()){
					plotting -= nPlot;
					nPlot = ma[iStr].size()-plotting;
					plotting = ma[iStr].size();
				}
				nCellsRun[iStr] = new TH1D(Form("%s, %i runs, %s clusters",DataSets[iStr].Data(), nRuns[iStr], calo.Data()),Form("%s, %i runs, Runwise Hot Cells for %s; CellID; Number of Runs",DataSets[iStr].Data(), nRuns[iStr], calo.Data()),nPlot,0,nPlot);
				OnlyEditTH1(nCellsRun[iStr], 20, 1, kBlack, kBlack);

				for (Int_t h=1; h<=nPlot; it++,h++){
					if(it == ma[iStr].end()) {cout << "ERROR while plotting HotCellsRunwiseOverview" << endl;break;}
					nCellsRun[iStr]->GetXaxis()->SetBinLabel(h,Form("%i",it->first));
					nCellsRun[iStr]->SetBinContent(h,it->second);
				}
				nCellsRun[iStr]->SetFillStyle(3004);
				nCellsRun[iStr]->SetFillColor(1);
				nCellsRun[iStr]->Draw();
				if(useDataRunListForMC && iStr>=nData) SaveCanvas(cvsRun, Form("%s/%s-%s/%s_%i.%s", outputDir.Data(), DataSets[iStr].Data(),DataSets[0].Data(), "HotCells_FiredInNRuns",nPart++,suffix.Data()));
				else SaveCanvas(cvsRun, Form("%s/%s/%s_%i.%s", outputDir.Data(), DataSets[iStr].Data(), "HotCells_FiredInNRuns",nPart++,suffix.Data()));
				delete nCellsRun[iStr];
			}while(plotting < (Int_t)ma[iStr].size());

			plotting = 0;
			nPlot = 100;
			nPart = 0;
			Int_t posInVec = 0;
			std::vector<TString>::iterator itVec = vec[iStr].begin();
			do
			{
				plotting+=nPlot;
				if(plotting>nRuns[iStr]){
					plotting -= nPlot;
					nPlot = nRuns[iStr]-plotting;
					plotting = nRuns[iStr];
				}
				nCellsRun[iStr] = new TH1D(Form("%s, %i runs, %s clusters",DataSets[iStr].Data(), nRuns[iStr], calo.Data()),Form("%s, %i runs, Runwise Hot Cells for %s; Run Number; Number of Hot Cells",DataSets[iStr].Data(), nRuns[iStr], calo.Data()),nPlot,0,nPlot);
				OnlyEditTH1(nCellsRun[iStr], 20, 1, kBlack, kBlack);

				for (Int_t h=1; h<=nPlot; h++){
					Int_t nBadCells = 0;
					do{
						if((vec[iStr].at(posInVec)).BeginsWith("Run-")){
							TString tempStr = (vec[iStr].at(posInVec)).Remove(0,4);
							nCellsRun[iStr]->GetXaxis()->SetBinLabel(h,Form("%s",tempStr.Data()));
						}else nBadCells++;
					}while(++itVec!=vec[iStr].end() && !((vec[iStr].at(++posInVec)).BeginsWith("Run-")));
					nCellsRun[iStr]->SetBinContent(h,nBadCells);
				}
				nCellsRun[iStr]->SetFillStyle(3004);
				nCellsRun[iStr]->SetFillColor(1);
				nCellsRun[iStr]->Draw();
				if(useDataRunListForMC && iStr>=nData) SaveCanvas(cvsRun, Form("%s/%s-%s/%s_%i.%s", outputDir.Data(), DataSets[iStr].Data(),DataSets[0].Data(), "HotCells_Runwise",nPart++,suffix.Data()));
				else SaveCanvas(cvsRun, Form("%s/%s/%s_%i.%s", outputDir.Data(), DataSets[iStr].Data(), "HotCells_Runwise",nPart++,suffix.Data()));
				delete nCellsRun[iStr];
			}while(plotting < nRuns[iStr]);
			fLogRunwiseHotCells[iStr].close();
		}
		delete cvsRun;
        delete[] nRuns;
        vec.clear();
        ma.clear();
	}

    delete[] vecBadCells;

    delete[] vecClusterEtaPhi;
    delete[] vecClusterEnergyTime;
    delete[] vecClusterEVsNCells;
    delete[] vecClusterIncludedCells;
    delete[] vecClusterIncludedCellsBefore;
    delete[] vecClusterEFracCells;
    delete[] vecClusterEFracCellsBefore;
    delete[] vecClusterEnergyMeanSigma;
    delete[] vecClusterTimeMeanSigma;
    delete[] vecClusterFiredCellIDs;
    delete[] vecClusterMissingCellIDs;
    delete[] vecClusterBadCells;

    delete[] vecClusterEnergy;
    delete[] vecClusterM02;
    delete[] vecClusterM20;
    delete[] vecClusterNCells;
    delete[] vecClusterDispersion;

    delete[] vecClusterDeltaEta;
    delete[] vecClusterDeltaPhi;
    delete[] vecClusterPi0ConvPhotonEta;
    delete[] vecClusterPi0ConvPhotonPhi;
    delete[] vecClusterEtaConvPhotonEta;
    delete[] vecClusterEtaConvPhotonPhi;

    delete[] vecBadCellsEnergy;
    delete[] vecBadCellsTime;

    for (Int_t i=0; i<nSets; i++) {
      delete[] vecClusterTime[i];
      delete[] vecClusterEVsModule[i];
      delete[] vecModuleEVsModule[i];
      delete[] vecClusterNCells100VsModule[i];
      delete[] vecClusterNCells1500VsModule[i];
    }
    delete[] vecClusterTime;
    delete[] vecClusterEVsModule;
    delete[] vecModuleEVsModule;
    delete[] vecClusterNCells100VsModule;
    delete[] vecClusterNCells1500VsModule;

    delete[] vecHistos;
    delete[] vecMissingRuns;

    delete[] fLog;
    delete[] fLogRunwiseHotCells;
    delete[] fLogRunwiseColdCells;

	delete fDeltaPt;
	TH1::AddDirectory(kTRUE);

	cout << "Done with ClusterQA_Runwise" << endl;
	cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

	return;

}//end
