#include "QA.h"
#include "../TaskV1/BuildHistogramsForGammaQAAdvV3.C"

void PhotonQA_Runwise(
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
	cout << "PhotonQA_Runwise" << endl;
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
		hMarkerStyle[i]=GetDefaultMarkerStyle(fEnergyFlag.Data(),DataSets[i].Data(),"");
		hMarkerColor[i]=GetColorDefaultColor(fEnergyFlag.Data(),DataSets[i].Data(),"");
		hLineColor[i]=GetColorDefaultColor(fEnergyFlag.Data(),DataSets[i].Data(),"");
		hMarkerSize[i]=markerSize;
	}

	Int_t fMode = mode;
	// mode:	0 // new output PCM-PCM
	//			1 // new output PCM dalitz
	//			2 // new output PCM-EMCal
	//			3 // new output PCM-PHOS
	//          4 // new output EMCal-EMCal
	//          5 // new output PHOS-PHOS
	//			9 // old output PCM-PCM
    //			10 // merged EMCal
    //			11 // merged PHOS
    if(fMode == 4 || fMode == 5 || fMode == 10 || fMode == 11){ cout << "Returning, given mode contains no PCM information: " << fMode << endl; return;}

	std::vector<TString> vecRuns;

//******************************************************************************
	TString fileRuns = Form("runNumbers%s.txt", (vecDataSet.at(0)).Data());
	if(!readin(fileRuns, vecRuns, kFALSE)) {cout << "ERROR, no Run Numbers could be found! Returning..." << endl; return;}

	TFile* fPhotonQAFile = new TFile(Form("%s/%s/%s/AnalysisResults.root", filePath.Data(), ((TString)vecDataSet.at(0)).Data(), ((TString)vecRuns.at(0)).Data()));
	if(fPhotonQAFile->IsZombie()) {cout << "ERROR: ROOT file '" << Form("%s/%s/%s/AnalysisResults.root", filePath.Data(), ((TString)vecDataSet.at(0)).Data(), ((TString)vecRuns.at(0)).Data()) << "' could not be openend, return!" << endl; return;}

	TString nameCutsPQA;
	TString nameCutsPQAshort;
	TString nameMainDir;

	TKey *keyPQA;
	TIter nextPQA(fPhotonQAFile->GetListOfKeys());
	while ((keyPQA=(TKey*)nextPQA())){
		cout << Form("Found TopDir: '%s' ",keyPQA->GetName());
		nameMainDir = keyPQA->GetName();
	}
	nameCutsPQA = Form("%s",nameMainDir.Data());
	nameCutsPQAshort = Form("%s",nameMainDir.Data());
	nameCutsPQA.Replace(0,15,"");
	nameCutsPQAshort.Replace(0,15,"");
	nameCutsPQAshort.Replace(8,35,"");

	cout << endl;
	cout << "long cutnumber for PhotonQA: '" << nameCutsPQA << "'\nshort cutnumber for PhotonQA: '"<< nameCutsPQAshort << "'" << endl;

	fPhotonQAFile->Close();
	delete fPhotonQAFile;
	cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
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

	TString fDetectionProcess = ReturnFullTextReconstructionProcess(fMode);

	TString outputDir = Form("%s/%s/PhotonQA/%s/Runwise",nameCutsPQA.Data(),fEnergyFlag.Data(),suffix.Data());
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

    TString lepton[2] = {"Electron","Positron"};
    TString charge[2] = {"-","+"};

    TH1D* hGammaN[nSets];
	TH1D* hGammaPt[nSets];
	TH1D* hGammaAlpha[nSets];
	TH1D* hGammaEta[nSets];
	TH1D* hGammaEtaNeg[nSets];
	TH1D* hGammaEtaPos[nSets];
    TH1D* hGammaPhi[nSets];
    TH1D* hGammaPhiLow[nSets];
    TH1D* hGammaPhiHigh[nSets];
    TH1D* hGammaPhiEtaNeg[nSets];
    TH1D* hGammaPhiEtaPos[nSets];
    TH1D* hGammaChi2[nSets];
    TH1D* hGammaPsiPair[nSets];
    TH1D* hGammaCosPoint[nSets];
    TH1D* hGammaInvMass[nSets];
    TH1D* hGammaQt[nSets];
    TH1D* hGammaAsym[nSets];
    TH1D* hGammaR[nSets];
    TH1D* hNSdEdx[nSets][2];
    TH1D* hEta[nSets][2];
    TH1D* hPt[nSets][2];
    TH1D* hTPCclusters[nSets][2];
    TH1D* hfTPCclusters[nSets][2];

    TH1D* hGammaPtRMS[nSets];
    TH1D* hGammaAlphaRMS[nSets];
    TH1D* hGammaEtaRMS[nSets];
    TH1D* hGammaPhiRMS[nSets];
    TH1D* hGammaChi2RMS[nSets];
    TH1D* hGammaPsiPairRMS[nSets];
    TH1D* hGammaCosPointRMS[nSets];
    TH1D* hGammaInvMassRMS[nSets];
    TH1D* hGammaQtRMS[nSets];
    TH1D* hGammaAsymRMS[nSets];
    TH1D* hGammaRRMS[nSets];
    TH1D* hNSdEdxRMS[nSets][2];
    TH1D* hEtaRMS[nSets][2];
    TH1D* hPtRMS[nSets][2];
    TH1D* hTPCclustersRMS[nSets][2];
    TH1D* hfTPCclustersRMS[nSets][2];

    std::vector<TH1D*>* vecGammaPt = new std::vector<TH1D*>[nSets];
    std::vector<TH1D*>* vecGammaAlpha = new std::vector<TH1D*>[nSets];
    std::vector<TH1D*>* vecGammaEta = new std::vector<TH1D*>[nSets];
    std::vector<TH1D*>* vecGammaPhi = new std::vector<TH1D*>[nSets];
    std::vector<TH1D*>* vecGammaChi2 = new std::vector<TH1D*>[nSets];
    std::vector<TH1D*>* vecGammaPsiPair = new std::vector<TH1D*>[nSets];
    std::vector<TH1D*>* vecGammaCosPoint = new std::vector<TH1D*>[nSets];
    std::vector<TH1D*>* vecGammaInvMass = new std::vector<TH1D*>[nSets];

    std::vector<TH2D*>* vecInvMassR = new std::vector<TH2D*>[nSets];
    std::vector<TH2D*>* vecEtaR = new std::vector<TH2D*>[nSets];
    std::vector<TH2D*>* vecPhiR = new std::vector<TH2D*>[nSets];
    std::vector<TH2D*>* vecAlphaR = new std::vector<TH2D*>[nSets];
    std::vector<TH2D*>* vecPsiPairR = new std::vector<TH2D*>[nSets];
    std::vector<TH2D*>* vecAsymR = new std::vector<TH2D*>[nSets];

    std::vector<TH2D*>** vecNSdEdx = new std::vector<TH2D*>*[nSets];
    std::vector<TH2D*>** vecTPCClusR = new std::vector<TH2D*>*[nSets];
    for (Int_t i=0; i<nSets; i++) {
      vecNSdEdx[i] = new std::vector<TH2D*>[2];
      vecTPCClusR[i] = new std::vector<TH2D*>[2];
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
	std::vector<TString> vecHistosName;
	TString histoName;

	for(Int_t i=0; i<nSets; i++)
	{
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

		histoName = "hGammaAlpha";
        if(i==0) vecHistosName.push_back(histoName);
        hGammaAlpha[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hGammaAlpha; Run Number ; Mean Photon #alpha",hNBin,hFBin,hLBin);
		EditTH1(globalRuns, doEquidistantXaxis, hGammaAlpha[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
		vecHistos[i].push_back(hGammaAlpha[i]);

        histoName = "hGammaAlphaRMS";
        if(i==0) vecHistosName.push_back(histoName);
        hGammaAlphaRMS[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hGammaAlphaRMS; Run Number ; #sigma_{Photon #alpha}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hGammaAlphaRMS[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hGammaAlphaRMS[i]);

		histoName = "hGammaEta";
        if(i==0) vecHistosName.push_back(histoName);
        hGammaEta[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hGammaEta; Run Number ; Mean Photon #eta",hNBin,hFBin,hLBin);
		EditTH1(globalRuns, doEquidistantXaxis, hGammaEta[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
		vecHistos[i].push_back(hGammaEta[i]);

        histoName = "hGammaEtaRMS";
        if(i==0) vecHistosName.push_back(histoName);
        hGammaEtaRMS[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hGammaEtaRMS; Run Number ; #sigma_{Photon #eta}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hGammaEtaRMS[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hGammaEtaRMS[i]);

		histoName = "hGammaEtaNeg";
        if(i==0) vecHistosName.push_back(histoName);
        hGammaEtaNeg[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hGammaEtaNeg; Run Number ; #gamma_{#eta < 0}/#gamma_{total}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hGammaEtaNeg[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hGammaEtaNeg[i]);

		histoName = "hGammaEtaPos";
        if(i==0) vecHistosName.push_back(histoName);
		hGammaEtaPos[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hGammaEtaPos; Run Number ; #gamma_{#eta > 0}/#gamma_{total}",hNBin,hFBin,hLBin);
		EditTH1(globalRuns, doEquidistantXaxis, hGammaEtaPos[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
		vecHistos[i].push_back(hGammaEtaPos[i]);

        histoName = "hGammaPhi";
        if(i==0) vecHistosName.push_back(histoName);
        hGammaPhi[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hGammaPhi; Run Number ; Mean Photon #phi",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hGammaPhi[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hGammaPhi[i]);

        histoName = "hGammaPhiRMS";
        if(i==0) vecHistosName.push_back(histoName);
        hGammaPhiRMS[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hGammaPhiRMS; Run Number ; #sigma_{Photon #phi}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hGammaPhiRMS[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hGammaPhiRMS[i]);

        histoName = "hGammaPhiLow";
        if(i==0) vecHistosName.push_back(histoName);
        hGammaPhiLow[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hGammaPhiLow; Run Number ; #gamma_{#phi < #pi}/#gamma_{total}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hGammaPhiLow[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hGammaPhiLow[i]);

        histoName = "hGammaPhiHigh";
        if(i==0) vecHistosName.push_back(histoName);
        hGammaPhiHigh[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hGammaPhiHigh; Run Number ; #gamma_{#phi > #pi}/#gamma_{total}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hGammaPhiHigh[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hGammaPhiHigh[i]);

        histoName = "hGammaPhiEtaNeg";
        if(i==0) vecHistosName.push_back(histoName);
        hGammaPhiEtaNeg[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hGammaPhiEtaNeg; Run Number ; Mean Photon #phi C side",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hGammaPhiEtaNeg[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hGammaPhiEtaNeg[i]);

        histoName = "hGammaPhiEtaPos";
        if(i==0) vecHistosName.push_back(histoName);
        hGammaPhiEtaPos[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hGammaPhiEtaPos; Run Number ; Mean Photon #phi A side",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hGammaPhiEtaPos[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hGammaPhiEtaPos[i]);

        histoName = "hGammaChi2";
        if(i==0) vecHistosName.push_back(histoName);
        hGammaChi2[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hGammaChi2; Run Number ; Mean #chi^{2}_{#gamma}/ndf",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hGammaChi2[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hGammaChi2[i]);

        histoName = "hGammaChi2RMS";
        if(i==0) vecHistosName.push_back(histoName);
        hGammaChi2RMS[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hGammaChi2RMS; Run Number ; #sigma_{#chi^{2}_{#gamma}/ndf}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hGammaChi2RMS[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hGammaChi2RMS[i]);

        histoName = "hGammaPsiPair";
        if(i==0) vecHistosName.push_back(histoName);
        hGammaPsiPair[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hGammaPsiPair; Run Number ; Mean #psi_{pair}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hGammaPsiPair[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hGammaPsiPair[i]);

        histoName = "hGammaPsiPairRMS";
        if(i==0) vecHistosName.push_back(histoName);
        hGammaPsiPairRMS[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hGammaPsiPairRMS; Run Number ; #sigma_{#psi_{pair}}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hGammaPsiPairRMS[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hGammaPsiPairRMS[i]);

        histoName = "hGammaCosPoint";
        if(i==0) vecHistosName.push_back(histoName);
        hGammaCosPoint[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hGammaCosPoint; Run Number ; Mean cos(#theta_{point})",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hGammaCosPoint[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hGammaCosPoint[i]);

        histoName = "hGammaCosPointRMS";
        if(i==0) vecHistosName.push_back(histoName);
        hGammaCosPointRMS[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hGammaCosPointRMS; Run Number ; #sigma_{cos(#theta_{point})}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hGammaCosPointRMS[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hGammaCosPointRMS[i]);

        histoName = "hGammaInvMass";
        if(i==0) vecHistosName.push_back(histoName);
        hGammaInvMass[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hGammaInvMass; Run Number ; Mean M_{e^{+}e^{-}} (GeV/#it{c}^{2})",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hGammaInvMass[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hGammaInvMass[i]);

        histoName = "hGammaInvMassRMS";
        if(i==0) vecHistosName.push_back(histoName);
        hGammaInvMassRMS[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hGammaInvMassRMS; Run Number ; #sigma_{M_{e^{+}e^{-}}} (GeV/#it{c}^{2})",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hGammaInvMassRMS[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hGammaInvMassRMS[i]);

        histoName = "hGammaQt";
        if(i==0) vecHistosName.push_back(histoName);
        hGammaQt[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hGammaQt; Run Number ; Mean Photon #it{q}_{T}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hGammaQt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hGammaQt[i]);

        histoName = "hGammaQtRMS";
        if(i==0) vecHistosName.push_back(histoName);
        hGammaQtRMS[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hGammaQtRMS; Run Number ; #sigma_{Photon #it{q}_{T}}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hGammaQtRMS[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hGammaQtRMS[i]);

        histoName = "hGammaAsym";
        if(i==0) vecHistosName.push_back(histoName);
        hGammaAsym[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hGammaAsym; Run Number ; Mean Asym",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hGammaAsym[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hGammaAsym[i]);

        histoName = "hGammaAsymRMS";
        if(i==0) vecHistosName.push_back(histoName);
        hGammaAsymRMS[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hGammaAsymRMS; Run Number ; #sigma_{Asym}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hGammaAsymRMS[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hGammaAsymRMS[i]);

        histoName = "hGammaR";
        if(i==0) vecHistosName.push_back(histoName);
        hGammaR[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hGammaR; Run Number ; Mean Photon R",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hGammaR[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hGammaR[i]);

        histoName = "hGammaRRMS";
        if(i==0) vecHistosName.push_back(histoName);
        hGammaRRMS[i] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),"hGammaRRMS; Run Number ; #sigma_{Photon R}",hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hGammaRRMS[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hGammaRRMS[i]);

        //-----------------------------------
        //-----------------------------------
        //-- Leptons
        //-----------------------------------
        //-----------------------------------

        for(Int_t iL=0; iL<2; iL++){
          //iL==0, Electron - iL==1, Positron
          histoName = Form("h%sEta",lepton[iL].Data());
          if(i==0) vecHistosName.push_back(histoName);
          hEta[i][iL] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%sEta; Run Number ; Mean %s #eta",lepton[iL].Data(),lepton[iL].Data()),hNBin,hFBin,hLBin);
          EditTH1(globalRuns, doEquidistantXaxis, hEta[i][iL], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
          vecHistos[i].push_back(hEta[i][iL]);

          histoName = Form("h%sEtaRMS",lepton[iL].Data());
          if(i==0) vecHistosName.push_back(histoName);
          hEtaRMS[i][iL] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%sEtaRMS; Run Number ; #sigma_{%s #eta}",lepton[iL].Data(),lepton[iL].Data()),hNBin,hFBin,hLBin);
          EditTH1(globalRuns, doEquidistantXaxis, hEtaRMS[i][iL], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
          vecHistos[i].push_back(hEtaRMS[i][iL]);

          histoName = Form("h%sPt",lepton[iL].Data());
          if(i==0) vecHistosName.push_back(histoName);
          hPt[i][iL] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%sPt; Run Number ; Mean %s #it{p}_{T}",lepton[iL].Data(),lepton[iL].Data()),hNBin,hFBin,hLBin);
          EditTH1(globalRuns, doEquidistantXaxis, hPt[i][iL], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
          vecHistos[i].push_back(hPt[i][iL]);

          histoName = Form("h%sPtRMS",lepton[iL].Data());
          if(i==0) vecHistosName.push_back(histoName);
          hPtRMS[i][iL] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%sPtRMS; Run Number ; #sigma_{%s #it{p}_{T}}",lepton[iL].Data(),lepton[iL].Data()),hNBin,hFBin,hLBin);
          EditTH1(globalRuns, doEquidistantXaxis, hPtRMS[i][iL], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
          vecHistos[i].push_back(hPtRMS[i][iL]);

          histoName = Form("hNSdEdx%s",lepton[iL].Data());
          if(i==0) vecHistosName.push_back(histoName);
          hNSdEdx[i][iL] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("hNSdEdx%s; Run Number ; Mean #it{n}#sigma_{e^{%s}} d#it{E}/d#it{x} TPC",lepton[iL].Data(),charge[iL].Data()),hNBin,hFBin,hLBin);
          EditTH1(globalRuns, doEquidistantXaxis, hNSdEdx[i][iL], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
          vecHistos[i].push_back(hNSdEdx[i][iL]);

          histoName = Form("hNSdEdxRMS%s",lepton[iL].Data());
          if(i==0) vecHistosName.push_back(histoName);
          hNSdEdxRMS[i][iL] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("hNSdEdxRMS%s; Run Number ; #sigma_{#it{n}#sigma_{e^{%s}} d#it{E}/d#it{x} TPC}",lepton[iL].Data(),charge[iL].Data()),hNBin,hFBin,hLBin);
          EditTH1(globalRuns, doEquidistantXaxis, hNSdEdxRMS[i][iL], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
          vecHistos[i].push_back(hNSdEdxRMS[i][iL]);

          histoName = Form("hTPCclusters%s",lepton[iL].Data());
          if(i==0) vecHistosName.push_back(histoName);
          hTPCclusters[i][iL] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("hTPCclusters%s; Run Number ; Mean number of TPC clusters in e^{%s} track",lepton[iL].Data(),charge[iL].Data()),hNBin,hFBin,hLBin);
          EditTH1(globalRuns, doEquidistantXaxis, hTPCclusters[i][iL], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
          vecHistos[i].push_back(hTPCclusters[i][iL]);

          histoName = Form("hTPCclustersRMS%s",lepton[iL].Data());
          if(i==0) vecHistosName.push_back(histoName);
          hTPCclustersRMS[i][iL] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("hTPCclustersRMS%s; Run Number ; #sigma_{number of TPC clusters in e^{%s} track}",lepton[iL].Data(),charge[iL].Data()),hNBin,hFBin,hLBin);
          EditTH1(globalRuns, doEquidistantXaxis, hTPCclustersRMS[i][iL], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
          vecHistos[i].push_back(hTPCclustersRMS[i][iL]);

          histoName = Form("hfTPCclusters%s",lepton[iL].Data());
          if(i==0) vecHistosName.push_back(histoName);
          hfTPCclusters[i][iL] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("hfTPCclusters%s; Run Number ; Mean number of TPC clusters in e^{%s} track / findable",lepton[iL].Data(),charge[iL].Data()),hNBin,hFBin,hLBin);
          EditTH1(globalRuns, doEquidistantXaxis, hfTPCclusters[i][iL], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
          vecHistos[i].push_back(hfTPCclusters[i][iL]);

          histoName = Form("hfTPCclustersRMS%s",lepton[iL].Data());
          if(i==0) vecHistosName.push_back(histoName);
          hfTPCclustersRMS[i][iL] = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("hfTPCclustersRMS%s; Run Number ; #sigma_{number of TPC clusters in e^{%s} track / findable}",lepton[iL].Data(),charge[iL].Data()),hNBin,hFBin,hLBin);
          EditTH1(globalRuns, doEquidistantXaxis, hfTPCclustersRMS[i][iL], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
          vecHistos[i].push_back(hfTPCclustersRMS[i][iL]);
        }
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

		Bool_t doMergeOutput = kFALSE;
        TString mergeCommand = Form("hadd -f -k %s/%s/PhotonQA_%s.root",filePath.Data(),fDataSet.Data(),fDataSet.Data());
        //check if file exists, otherwise merge
        TFile* mergedFile = new TFile(Form("%s/%s/PhotonQA_%s.root",filePath.Data(),fDataSet.Data(),fDataSet.Data()),"READ");
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
		for(Int_t j=0; j<(Int_t) vecRuns.size(); j++)
        {
			fRunNumber = vecRuns.at(j);
			if(i>=nData) fRootFile = Form("%s/%s/%s/PhotonQA_MC.root", filePath.Data(), fDataSet.Data(), fRunNumber.Data());
			else fRootFile = Form("%s/%s/%s/PhotonQA_Data.root", filePath.Data(), fDataSet.Data(), fRunNumber.Data());

			cout << endl;
			cout << "\t\t----------------------------------------------------------------------------" << endl;
			cout << Form("\t\tRun %s", fRunNumber.Data()) << endl;
			cout << "\t\tProcessing file: " << fRootFile.Data() << endl;
			cout << "\t\t----------------------------------------------------------------------------" << endl;
			cout << endl;

            Bool_t isPhotonQAFileOpen = kTRUE;
			TFile* fPhotonQAFile = new TFile(fRootFile.Data(),"READ");
            if(fPhotonQAFile->IsZombie()){
				cout << endl;
				cout << "\t\t----------------------------------------------------------------------------" << endl;
				cout << "\t\t" << fRootFile << ", could not be opnened." << endl;
				cout << endl;
                isPhotonQAFileOpen = kFALSE;
				fPhotonQAFile->Close();
				delete fPhotonQAFile;
				TString fAnalysisResultsFile = Form("%s/%s/%s/AnalysisResults.root", filePath.Data(), fDataSet.Data(), fRunNumber.Data());
				cout << "\t\tTrying to open " << fAnalysisResultsFile << "..." << endl;
				TFile* RootFilePQA = new TFile(fAnalysisResultsFile.Data(),"READ");
				if(RootFilePQA->IsZombie()) {
                  vecMissingRuns[i].push_back(fRunNumber);
                  cout << "INFO: ROOT file '" << fAnalysisResultsFile.Data() << "' could not be openend, continue!" << endl;
                  cout << "\t\t----------------------------------------------------------------------------" << endl;
                  RootFilePQA->Close();
                  delete RootFilePQA;
                  continue;
				}
				RootFilePQA->Close();
				delete RootFilePQA;

				cout << "\t\tCalling BuildHistogramsForGammaQAAdvV3..." << endl;
				BuildHistogramsForGammaQAAdvV3(fAnalysisResultsFile.Data(),nameCutsPQA,nameCutsPQAshort,(i>=nData),0,Form("%s/%s/%s/PhotonQA",filePath.Data(), fDataSet.Data(), fRunNumber.Data()),kTRUE);
				cout << "\t\tdone!" << endl;
				cout << "\t\t----------------------------------------------------------------------------" << endl;
				doMergeOutput = kTRUE;
                cout << "--- Switched on merging of Runs! ---" << endl;
            }

			mergeCommand += Form(" %s",fRootFile.Data());
            if(!isPhotonQAFileOpen){
              fPhotonQAFile = new TFile(fRootFile.Data(),"READ");
                if(fPhotonQAFile->IsZombie()) {vecMissingRuns[i].push_back(fRunNumber); cout << "INFO: ROOT file '" << fRootFile.Data() << "' could not be openend, continue!" << endl; fPhotonQAFile->Close(); delete fPhotonQAFile; continue;}
            }

			TDirectory* TopDir = (TDirectory*) fPhotonQAFile->Get(Form("GammaConvV1_QA_%s",nameCutsPQA.Data()));
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
			if(GammaEtaPt){
				TH1D * GammaEta = (TH1D*) GammaEtaPt->ProjectionX("GammaEta");
				hGammaEta[i]->SetBinContent(bin, GammaEta->GetMean());
				hGammaEta[i]->SetBinError(bin, GammaEta->GetMeanError());
                hGammaEtaRMS[i]->SetBinContent(bin, GammaEta->GetRMS());
                hGammaEtaRMS[i]->SetBinError(bin, GammaEta->GetRMSError());
				hGammaEtaNeg[i]->SetBinContent(bin, GetHistogramIntegral(GammaEta,-1,0)/GammaEta->Integral());
				hGammaEtaNeg[i]->SetBinError(bin, GetHistogramIntegralError(GammaEta,-1,0)/GammaEta->Integral());
				hGammaEtaPos[i]->SetBinContent(bin, GetHistogramIntegral(GammaEta,0,1)/GammaEta->Integral());
				hGammaEtaPos[i]->SetBinError(bin, GetHistogramIntegralError(GammaEta,0,1)/GammaEta->Integral());

				TH1D * GammaPt = (TH1D*) GammaEtaPt->ProjectionY("GammaPt");
				hGammaPt[i]->SetBinContent(bin, GammaPt->GetMean());
				hGammaPt[i]->SetBinError(bin, GammaPt->GetMeanError());
                hGammaPtRMS[i]->SetBinContent(bin, GammaPt->GetRMS());
                hGammaPtRMS[i]->SetBinError(bin, GammaPt->GetRMSError());
                GammaPt->GetXaxis()->SetTitle("Photon #it{p}_{T}");
                GammaPt->GetYaxis()->SetTitle("#frac{1}{N_{#gamma}} #frac{dN}{d#it{p}_{T}}");
                GammaPt->Sumw2();
                GammaPt->Scale(1 / nGammas);
				vecGammaPt[i].push_back(GammaPt);

				GammaEta->GetXaxis()->SetTitle("Photon #eta");
				GammaEta->GetXaxis()->SetRangeUser(-1,1);
                GammaEta->GetYaxis()->SetTitle("#frac{1}{N_{#gamma}} #frac{dN}{d#eta}");
				GammaEta->Sumw2();
                GammaEta->Scale(1 / nGammas);
				vecGammaEta[i].push_back(GammaEta);
                delete GammaEtaPt;
            }else cout << "INFO: Object |histoGammaEtaPt| could not be found! Skipping Fill..." << endl;
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

            TH1D* GammaPhiEtaNeg = (TH1D*) TopDir->Get("histoGammaPhiEtaNeg");
            if(GammaPhiEtaNeg){
                  hGammaPhiEtaNeg[i]->SetBinContent(bin, GammaPhiEtaNeg->GetMean());
                  hGammaPhiEtaNeg[i]->SetBinError(bin, GammaPhiEtaNeg->GetMeanError());
                  delete GammaPhiEtaNeg;
            }else cout << "INFO: Object |histoGammaPhiEtaNeg| could not be found! Skipping Fill..." << endl;

            TH1D* GammaPhiEtaPos = (TH1D*) TopDir->Get("histoGammaPhiEtaPos");
            if(GammaPhiEtaPos){
                  hGammaPhiEtaPos[i]->SetBinContent(bin, GammaPhiEtaPos->GetMean());
                  hGammaPhiEtaPos[i]->SetBinError(bin, GammaPhiEtaPos->GetMeanError());
                  delete GammaPhiEtaPos;
            }else cout << "INFO: Object |histoGammaPhiEtaPos| could not be found! Skipping Fill..." << endl;

    //---------
            TH2D* GammaAlphaQt = (TH2D*) TopDir->Get("histoGammaAlphaQt");
			if(GammaAlphaQt){
				TH1D * GammaAlpha = (TH1D*) GammaAlphaQt->ProjectionX("GammaAlpha");
				hGammaAlpha[i]->SetBinContent(bin, GammaAlpha->GetMean());
				hGammaAlpha[i]->SetBinError(bin, GammaAlpha->GetMeanError());
                hGammaAlphaRMS[i]->SetBinContent(bin, GammaAlpha->GetRMS());
                hGammaAlphaRMS[i]->SetBinError(bin, GammaAlpha->GetRMSError());
                GammaAlpha->GetXaxis()->SetTitle("Photon #alpha");
                GammaAlpha->GetYaxis()->SetTitle("#frac{1}{N_{#gamma}} #frac{dN}{d#alpha}");
                GammaAlpha->Sumw2();
                GammaAlpha->Scale(1 / nGammas);
				vecGammaAlpha[i].push_back(GammaAlpha);
                delete GammaAlphaQt;
            }else cout << "INFO: Object |histoGammaAlphaQt| could not be found! Skipping Fill..." << endl;

            TH2D* GammaChi2Pt = (TH2D*) TopDir->Get("histoGammaChi2NDFPt");
			if(GammaChi2Pt){
				TH1D * GammaChi2 = (TH1D*) GammaChi2Pt->ProjectionX("GammaChi2");
				hGammaChi2[i]->SetBinContent(bin, GammaChi2->GetMean());
				hGammaChi2[i]->SetBinError(bin, GammaChi2->GetMeanError());
                hGammaChi2RMS[i]->SetBinContent(bin, GammaChi2->GetRMS());
                hGammaChi2RMS[i]->SetBinError(bin, GammaChi2->GetRMSError());
                GammaChi2->GetXaxis()->SetTitle("Photon #chi^{2}_{#gamma}/ndf");
                GammaChi2->GetYaxis()->SetTitle("#frac{1}{N_{#gamma}} #frac{dN}{d#chi^{2}}");
                GammaChi2->Sumw2();
                GammaChi2->Scale(1 / nGammas);
                vecGammaChi2[i].push_back(GammaChi2);
                delete GammaChi2Pt;
            }else cout << "INFO: Object |histoGammaChi2NDFPt| could not be found! Skipping Fill..." << endl;

            TH2D* GammaPsiPairPt = (TH2D*) TopDir->Get("histoGammaPsiPairPt");
			if(GammaPsiPairPt){
				TH1D * GammaPsiPair = (TH1D*) GammaPsiPairPt->ProjectionX("GammaPsiPair");
				hGammaPsiPair[i]->SetBinContent(bin, GammaPsiPair->GetMean());
				hGammaPsiPair[i]->SetBinError(bin, GammaPsiPair->GetMeanError());
                hGammaPsiPairRMS[i]->SetBinContent(bin, GammaPsiPair->GetRMS());
                hGammaPsiPairRMS[i]->SetBinError(bin, GammaPsiPair->GetRMSError());
                GammaPsiPair->GetXaxis()->SetTitle("Photon #psi_{pair}");
                GammaPsiPair->GetYaxis()->SetTitle("#frac{1}{N_{#gamma}} #frac{dN}{d#psi_{pair}}");
                GammaPsiPair->Sumw2();
                GammaPsiPair->Scale(1 / nGammas);
                vecGammaPsiPair[i].push_back(GammaPsiPair);
                delete GammaPsiPairPt;
            }else cout << "INFO: Object |histoGammaPsiPairPt| could not be found! Skipping Fill..." << endl;

           TH2D* GammaCosPointPt = (TH2D*) TopDir->Get("histoGammaCosPointPt");
			if(GammaCosPointPt){
				TH1D * GammaCosPoint = (TH1D*) GammaCosPointPt->ProjectionX("GammaCosPoint");
				hGammaCosPoint[i]->SetBinContent(bin, GammaCosPoint->GetMean());
				hGammaCosPoint[i]->SetBinError(bin, GammaCosPoint->GetMeanError());
                hGammaCosPointRMS[i]->SetBinContent(bin, GammaCosPoint->GetRMS());
                hGammaCosPointRMS[i]->SetBinError(bin, GammaCosPoint->GetRMSError());
                GammaCosPoint->GetXaxis()->SetTitle("Photon cos(#theta_{point})");
                GammaCosPoint->GetYaxis()->SetTitle("#frac{1}{N_{#gamma}} #frac{dN}{d cos(#theta_{point}})");
                GammaCosPoint->Sumw2();
                GammaCosPoint->Scale(1 / nGammas);
                vecGammaCosPoint[i].push_back(GammaCosPoint);
                delete GammaCosPointPt;
            }else cout << "INFO: Object |histoGammaCosPointPt| could not be found! Skipping Fill..." << endl;

            TH2D* GammaInvMassPt = (TH2D*) TopDir->Get("histoGammaInvMassPt");
			if(GammaInvMassPt){
				TH1D * GammaInvMass = (TH1D*) GammaInvMassPt->ProjectionX("GammaInvMass");
				hGammaInvMass[i]->SetBinContent(bin, GammaInvMass->GetMean());
				hGammaInvMass[i]->SetBinError(bin, GammaInvMass->GetMeanError());
                hGammaInvMassRMS[i]->SetBinContent(bin, GammaInvMass->GetRMS());
                hGammaInvMassRMS[i]->SetBinError(bin, GammaInvMass->GetRMSError());
                GammaInvMass->GetXaxis()->SetTitle("M_{e^{+}e^{-}}");
                GammaInvMass->GetYaxis()->SetTitle("#frac{1}{N_{#gamma}} #frac{dN}{dM_{e^{+}e^{-}}}");
                GammaInvMass->Sumw2();
                GammaInvMass->Scale(1 / nGammas);
                vecGammaInvMass[i].push_back(GammaInvMass);
                delete GammaInvMassPt;
            }else cout << "INFO: Object |histoGammaInvMassPt| could not be found! Skipping Fill..." << endl;

            TH2D* GammaQtPt = (TH2D*) TopDir->Get("histoGammaQtPt");
			if(GammaQtPt){
				TH1D * GammaQt = (TH1D*) GammaQtPt->ProjectionX("GammaQt");
				hGammaQt[i]->SetBinContent(bin, GammaQt->GetMean());
				hGammaQt[i]->SetBinError(bin, GammaQt->GetMeanError());
                hGammaQtRMS[i]->SetBinContent(bin, GammaQt->GetRMS());
                hGammaQtRMS[i]->SetBinError(bin, GammaQt->GetRMSError());
				delete GammaQt;
                delete GammaQtPt;
            }else cout << "INFO: Object |histoGammaQtPt| could not be found! Skipping Fill..." << endl;

            TH2D* GammaInvMassR = (TH2D*) TopDir->Get("histoGammaInvMassR");
            if(GammaInvMassR){
              GammaInvMassR->SetTitle(fRunNumber);
              vecInvMassR[i].push_back(GammaInvMassR);
            }else cout << "INFO: Object |histoGammaInvMassR| could not be found! Skipping Fill..." << endl;

            TH2D* GammaEtaR = (TH2D*) TopDir->Get("histoGammaEtaR");
            if(GammaEtaR){
              GammaEtaR->SetTitle(fRunNumber);
              vecEtaR[i].push_back(GammaEtaR);
            }else cout << "INFO: Object |histoGammaEtaR| could not be found! Skipping Fill..." << endl;

            TH2D* GammaPhiR = (TH2D*) TopDir->Get("histoGammaPhiR");
            if(GammaPhiR){
              GammaPhiR->SetTitle(fRunNumber);
              vecPhiR[i].push_back(GammaPhiR);
            }else cout << "INFO: Object |histoGammaPhiR| could not be found! Skipping Fill..." << endl;

            TH2D* GammaAlphaR = (TH2D*) TopDir->Get("histoGammaAlphaR");
            if(GammaAlphaR){
              GammaAlphaR->SetTitle(fRunNumber);
              vecAlphaR[i].push_back(GammaAlphaR);
            }else cout << "INFO: Object |histoGammaAlphaR| could not be found! Skipping Fill..." << endl;

            TH2D* GammaPsiPairR = (TH2D*) TopDir->Get("histoGammaPsiPairR");
            if(GammaPsiPairR){
              GammaPsiPairR->SetTitle(fRunNumber);
              vecPsiPairR[i].push_back(GammaPsiPairR);
            }else cout << "INFO: Object |histoGammaPsiPairR| could not be found! Skipping Fill..." << endl;

            TH2D* GammaAsymR = (TH2D*) TopDir->Get("histoGammaAsymR");
            if(GammaAsymR){
              TH1D * GammaAsym = (TH1D*) GammaAsymR->ProjectionX("GammaAsym");
              hGammaAsym[i]->SetBinContent(bin, GammaAsym->GetMean());
              hGammaAsym[i]->SetBinError(bin, GammaAsym->GetMeanError());
              hGammaAsymRMS[i]->SetBinContent(bin, GammaAsym->GetRMS());
              hGammaAsymRMS[i]->SetBinError(bin, GammaAsym->GetRMSError());
              delete GammaAsym;

              TH1D * GammaR = (TH1D*) GammaAsymR->ProjectionY("GammaR");
              hGammaR[i]->SetBinContent(bin, GammaR->GetMean());
              hGammaR[i]->SetBinError(bin, GammaR->GetMeanError());
              hGammaRRMS[i]->SetBinContent(bin, GammaR->GetRMS());
              hGammaRRMS[i]->SetBinError(bin, GammaR->GetRMSError());
              delete GammaR;

              GammaAsymR->SetTitle(fRunNumber);
              vecAsymR[i].push_back(GammaAsymR);
            }else cout << "INFO: Object |histoGammaAsymR| could not be found! Skipping Fill..." << endl;

            //-----------------------------------
            //-----------------------------------
            //-- Leptons
            //-----------------------------------
            //-----------------------------------

            for(Int_t iL=0; iL<2; iL++){
              //iL==0, Electron - iL==1, Positron

              TH3D* NSigmadEdxEtaP = (TH3D*) TopDir->Get(Form("histo%sNSigmadEdxEtaP",lepton[iL].Data()));
              if(NSigmadEdxEtaP){
                  TH1D * NSigmadEdx = (TH1D*) NSigmadEdxEtaP->Project3D("x");
                  hNSdEdx[i][iL]->SetBinContent(bin, NSigmadEdx->GetMean(1));
                  hNSdEdx[i][iL]->SetBinError(bin, NSigmadEdx->GetMeanError(1));
                  hNSdEdxRMS[i][iL]->SetBinContent(bin, NSigmadEdx->GetRMS(1));
                  hNSdEdxRMS[i][iL]->SetBinError(bin, NSigmadEdx->GetRMSError(1));
                  delete NSigmadEdx;
                  TH1D * Eta = (TH1D*) NSigmadEdxEtaP->Project3D("y");
                  hEta[i][iL]->SetBinContent(bin, Eta->GetMean(1));
                  hEta[i][iL]->SetBinError(bin, Eta->GetMeanError(1));
                  hEtaRMS[i][iL]->SetBinContent(bin, Eta->GetRMS(1));
                  hEtaRMS[i][iL]->SetBinError(bin, Eta->GetRMSError(1));
                  delete Eta;
                  TH1D * Pt = (TH1D*) NSigmadEdxEtaP->Project3D("z");
                  hPt[i][iL]->SetBinContent(bin, Pt->GetMean(1));
                  hPt[i][iL]->SetBinError(bin, Pt->GetMeanError(1));
                  hPtRMS[i][iL]->SetBinContent(bin, Pt->GetRMS(1));
                  hPtRMS[i][iL]->SetBinError(bin, Pt->GetRMSError(1));
                  delete Pt;
                  delete NSigmadEdxEtaP;
              }else cout << Form("INFO: Object |histo%sNSigmadEdxEtaP| could not be found! Skipping Fill...",lepton[iL].Data()) << endl;

              TH2D* TPCclustersPt= (TH2D*) TopDir->Get(Form("histo%sClPt",lepton[iL].Data()));
              if(TPCclustersPt){
                TH1D * TPCclusters = (TH1D*) TPCclustersPt->ProjectionX("TPCclusters");
                hTPCclusters[i][iL]->SetBinContent(bin, TPCclusters->GetMean());
                hTPCclusters[i][iL]->SetBinError(bin, TPCclusters->GetMeanError());
                hTPCclustersRMS[i][iL]->SetBinContent(bin, TPCclusters->GetRMS());
                hTPCclustersRMS[i][iL]->SetBinError(bin, TPCclusters->GetRMSError());
                delete TPCclusters;
                delete TPCclustersPt;
              }else cout << Form("INFO: Object |histo%sClPt| could not be found! Skipping Fill...", lepton[iL].Data()) << endl;

              TH2D* fTPCclustersPt= (TH2D*) TopDir->Get(Form("histo%sFClPt",lepton[iL].Data()));
              if(fTPCclustersPt){
                TH1D * fTPCclusters = (TH1D*) fTPCclustersPt->ProjectionX("fTPCclusters");
                hfTPCclusters[i][iL]->SetBinContent(bin, fTPCclusters->GetMean());
                hfTPCclusters[i][iL]->SetBinError(bin, fTPCclusters->GetMeanError());
                hfTPCclustersRMS[i][iL]->SetBinContent(bin, fTPCclusters->GetRMS());
                hfTPCclustersRMS[i][iL]->SetBinError(bin, fTPCclusters->GetRMSError());
                delete fTPCclusters;
                delete fTPCclustersPt;
              }else cout << Form("INFO: Object |histo%sFClPt| could not be found! Skipping Fill...", lepton[iL].Data()) << endl;
            //---------
              TH2D* NSdEdx = (TH2D*) TopDir->Get(Form("histo%sNSigmadEdxPhi",lepton[iL].Data()));
              if(NSdEdx){
                  vecNSdEdx[i][iL].push_back(NSdEdx);
              }else cout << Form("INFO: Object |histo%sNSigmadEdxPhi| could not be found! Skipping Fill...",lepton[iL].Data()) << endl;

              TH2D* TPCclustersR = (TH2D*) TopDir->Get(Form("histo%sClR",lepton[iL].Data()));
              if(TPCclustersR){
                  vecTPCClusR[i][iL].push_back(TPCclustersR);
              }else cout << Form("INFO: Object |histo%sClR| could not be found! Skipping Fill...",lepton[iL].Data()) << endl;

            }

			TopDir->Clear();
			delete TopDir;

			fPhotonQAFile->Close();
			delete fPhotonQAFile;
		}

		if(doMergeOutput) {
			cout << "Flag doMergeOutput = kTRUE, merging PhotonQA output..." << endl;
			gSystem->Exec(mergeCommand.Data());
			cout << "done!" << endl;
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

			TString fTrigger = "";
            if(doTrigger && i<nData){
				cout << "Obtaining trigger - ";
                TString fTriggerCut = nameCutsPQAshort(3,2);
				fTrigger = AnalyseSpecialTriggerCut(fTriggerCut.Atoi(), plotDataSets[i]);
				cout << "'" << fTrigger.Data() << "' - was found!" << endl;
				if(fTrigger.Contains("not defined")){
					fTrigger = "";
					cout << "INFO: Trigger cut not defined!" << endl;
				}
			}

            DrawVectorOverviewTH2D(canvas, vecInvMassR[i], "oInvMassR", outputDirDataSet, suffix,
                                        0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            DrawVectorOverviewTH2D(canvas, vecEtaR[i], "oEtaR", outputDirDataSet, suffix,
                                        0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            DrawVectorOverviewTH2D(canvas, vecPhiR[i], "oPhiR", outputDirDataSet, suffix,
                                        0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            DrawVectorOverviewTH2D(canvas, vecAlphaR[i], "oAlphaR", outputDirDataSet, suffix,
                                        0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            DrawVectorOverviewTH2D(canvas, vecPsiPairR[i], "oPsiPairR", outputDirDataSet, suffix,
                                        0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            DrawVectorOverviewTH2D(canvas, vecAsymR[i], "oAsymR", outputDirDataSet, suffix,
                                        0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);

            DrawVectorOverviewTH2D(canvas, vecNSdEdx[i][0], "oNSdEdx_Electron", outputDirDataSet, suffix,
                                        0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            DrawVectorOverviewTH2D(canvas, vecNSdEdx[i][1], "oNSdEdx_Positron", outputDirDataSet, suffix,
                                        0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            DrawVectorOverviewTH2D(canvas, vecTPCClusR[i][0], "oTPCClusR_Electron", outputDirDataSet, suffix,
                                        0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            DrawVectorOverviewTH2D(canvas, vecTPCClusR[i][1], "oTPCClusR_Positron", outputDirDataSet, suffix,
                                        0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);

		//--------
			for(Int_t h=0; h<(Int_t) vecHistos[i].size(); h++)
			{
				(((TH1D*) vecHistos[i].at(h)))->SetTitle("");
                AdjustHistRange(((TH1D*) vecHistos[i].at(h)),1.1,1.1,kTRUE);
                if(((TString)vecHistosName.at(h)).CompareTo("hGammaCosPoint")==0) AdjustHistRange(((TH1D*) vecHistos[i].at(h)),1.001,1.001,kTRUE);
                if(((TString)vecHistosName.at(h)).CompareTo("nGamma")==0) AdjustHistRange(((TH1D*) vecHistos[i].at(h)),10,10,kTRUE);
                ((TH1D*) vecHistos[i].at(h))->Draw("px0e1");

                if(doTrigger) PutProcessLabelAndEnergyOnPlot(0.8, 0.92, 0.03, fCollisionSystem.Data(), plotDataSets[i].Data(), fTrigger.Data());
                else PutProcessLabelAndEnergyOnPlot(0.8, 0.92, 0.03, fCollisionSystem.Data(), plotDataSets[i].Data(), "");

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

            DrawVectorRunwiseTH1D(	canvasRunwise, legendRuns, vecGammaPt[i], vecRuns, 5, 5, kTRUE, addRight, 0.8, 0.94, 0.03, 0.8, 0.8, 0.03,
                                    doTrigger, fTrigger, (Bool_t)(i<nData), outputDirDataSet, "GammaPt", plotDataSets[i], kFALSE,
                                    fCollisionSystem, "", suffix, kTRUE, kTRUE, kFALSE);

            DrawVectorRunwiseTH1D(	canvasRunwise, legendRuns, vecGammaAlpha[i], vecRuns, 2, 1.1, kTRUE, addRight, 0.8, 0.94, 0.03, 0.8, 0.8, 0.03,
                                    doTrigger, fTrigger, (Bool_t)(i<nData), outputDirDataSet, "GammaAlpha", plotDataSets[i], kFALSE,
                                    fCollisionSystem, "", suffix, kFALSE, kFALSE, kFALSE);

            DrawVectorRunwiseTH1D(	canvasRunwise, legendRuns, vecGammaEta[i], vecRuns, 2, 1.1, kTRUE, addRight, 0.8, 0.94, 0.03, 0.8, 0.8, 0.03,
                                    doTrigger, fTrigger, (Bool_t)(i<nData), outputDirDataSet, "GammaEta", plotDataSets[i], kFALSE,
                                    fCollisionSystem, "", suffix, kFALSE, kFALSE, kFALSE);

            DrawVectorRunwiseTH1D(	canvasRunwise, legendRuns, vecGammaPhi[i], vecRuns, 2, 1.1, kTRUE, addRight, 0.8, 0.94, 0.03, 0.8, 0.8, 0.03,
                                    doTrigger, fTrigger, (Bool_t)(i<nData), outputDirDataSet, "GammaPhi", plotDataSets[i], kFALSE,
                                    fCollisionSystem, "", suffix, kFALSE, kFALSE, kFALSE);

            DrawVectorRunwiseTH1D(	canvasRunwise, legendRuns, vecGammaChi2[i], vecRuns, 5, 5, kTRUE, addRight, 0.8, 0.94, 0.03, 0.8, 0.8, 0.03,
                                    doTrigger, fTrigger, (Bool_t)(i<nData), outputDirDataSet, "GammaChi2", plotDataSets[i], kFALSE,
                                    fCollisionSystem, "", suffix, kFALSE, kTRUE, kFALSE);

            DrawVectorRunwiseTH1D(	canvasRunwise, legendRuns, vecGammaPsiPair[i], vecRuns, 2, 1.1, kTRUE, addRight, 0.8, 0.94, 0.03, 0.8, 0.8, 0.03,
                                    doTrigger, fTrigger, (Bool_t)(i<nData), outputDirDataSet, "GammaPsiPair", plotDataSets[i], kFALSE,
                                    fCollisionSystem, "", suffix, kFALSE, kFALSE, kFALSE);

            DrawVectorRunwiseTH1D(	canvasRunwise, legendRuns, vecGammaCosPoint[i], vecRuns, 5, 5, kTRUE, addRight, 0.8, 0.94, 0.03, 0.8, 0.8, 0.03,
                                    doTrigger, fTrigger, (Bool_t)(i<nData), outputDirDataSet, "GammaCosPoint", plotDataSets[i], kFALSE,
                                    fCollisionSystem, "", suffix, kTRUE, kTRUE, kFALSE);

            DrawVectorRunwiseTH1D(	canvasRunwise, legendRuns, vecGammaInvMass[i], vecRuns, 5, 5, kTRUE, addRight, 0.8, 0.94, 0.03, 0.8, 0.8, 0.03,
                                    doTrigger, fTrigger, (Bool_t)(i<nData), outputDirDataSet, "GammaInvMass", plotDataSets[i], kFALSE,
                                    fCollisionSystem, "", suffix, kFALSE, kTRUE, kFALSE);

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
            if(((TString)vecHistosName.at(h)).CompareTo("hGammaCosPoint")==0) AdjustHistRange(vecHistos,1.001,1.001,h,nSets,kTRUE);
            if(((TString)vecHistosName.at(h)).CompareTo("hGammaN")==0) AdjustHistRange(vecHistos,10,10,h,nSets,kTRUE);
			for(Int_t i=nSets-1; i>=0; i--)
			{
                TString draw = (i==nSets-1)?"px0e1":"px0e1, same";
				((TH1D*) vecHistos[i].at(h))->SetTitle("");
				((TH1D*) vecHistos[i].at(h))->Draw(draw.Data());
			}
			legend->Draw();

            if(doTrigger) PutProcessLabelAndEnergyOnPlot(0.8, 0.92, 0.03, fCollisionSystem.Data(), fTrigger.Data(), "");
            else PutProcessLabelAndEnergyOnPlot(0.8, 0.92, 0.03, fCollisionSystem.Data(), "", "");

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

                      PutProcessLabelAndEnergyOnPlot(0.8, 0.92, 0.03, fCollisionSystem.Data(), fTrigger.Data(), "");
                    }else PutProcessLabelAndEnergyOnPlot(0.8, 0.92, 0.03, fCollisionSystem.Data(), "", "");

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

	for(Int_t i=0; i<nSets; i++)
	{
		fDataSet = vecDataSet.at(i);

		if(useDataRunListForMC && i>=nData) nameOutput = Form("%s/%s/PhotonQA/%s-%s_PhotonQARunwise.root",nameCutsPQA.Data(),fEnergyFlag.Data(),fDataSet.Data(),vecDataSet.at(0).Data());
		else nameOutput = Form("%s/%s/PhotonQA/%s_PhotonQARunwise.root",nameCutsPQA.Data(),fEnergyFlag.Data(),fDataSet.Data());

		TFile* fOutput = new TFile(nameOutput,"RECREATE");
		cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
		cout << "Output file: " << nameOutput << endl;
		cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
		fLog[i] << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
		fLog[i] << "Output file: " << nameOutput << endl;
		fLog[i] << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

        for(Int_t h=0; h<(Int_t) vecHistos[i].size(); h++) WriteHistogram(((TH1D*) vecHistos[i].at(h)));

		DeleteVecTH1D(vecGammaPt[i]);
        DeleteVecTH1D(vecGammaAlpha[i]);
		DeleteVecTH1D(vecGammaEta[i]);
        DeleteVecTH1D(vecGammaPhi[i]);
        DeleteVecTH1D(vecGammaChi2[i]);
        DeleteVecTH1D(vecGammaPsiPair[i]);
        DeleteVecTH1D(vecGammaCosPoint[i]);
        DeleteVecTH1D(vecGammaInvMass[i]);

        WriteHistogramTH2DVec(fOutput,vecAlphaR[i],"GammaAlphaR");
        WriteHistogramTH2DVec(fOutput,vecEtaR[i],"GammaEtaR");
        WriteHistogramTH2DVec(fOutput,vecPhiR[i],"GammaPhiR");
        WriteHistogramTH2DVec(fOutput,vecInvMassR[i],"GammaInvMassR");
        WriteHistogramTH2DVec(fOutput,vecPsiPairR[i],"GammaPsiPairR");
        WriteHistogramTH2DVec(fOutput,vecAsymR[i],"GammaAsymR");

        WriteHistogramTH2DVec(fOutput,vecNSdEdx[i][0],"ElectronNSdEdx");
        WriteHistogramTH2DVec(fOutput,vecNSdEdx[i][1],"PositronNSdEdx");
        WriteHistogramTH2DVec(fOutput,vecTPCClusR[i][0],"ElectronTPCClusR");
        WriteHistogramTH2DVec(fOutput,vecTPCClusR[i][1],"PositronTPCClusR");

		fOutput->Write();
		fOutput->Close();
		delete fOutput;
		fLog[i].close();
	}

    delete[] vecHistos;

    delete[] vecGammaPt;
    delete[] vecGammaAlpha;
    delete[] vecGammaEta;
    delete[] vecGammaPhi;
    delete[] vecGammaChi2;
    delete[] vecGammaPsiPair;
    delete[] vecGammaCosPoint;
    delete[] vecGammaInvMass;

    delete[] vecInvMassR;
    delete[] vecEtaR;
    delete[] vecPhiR;
    delete[] vecAlphaR;
    delete[] vecPsiPairR;
    delete[] vecAsymR;
    for (Int_t i=0; i<nSets; i++) {
      delete[] vecNSdEdx[i];
      delete[] vecTPCClusR[i];
    }
    delete[] vecNSdEdx;
    delete[] vecTPCClusR;

    delete[] vecMissingRuns;

    delete[] fLog;

	TH1::AddDirectory(kTRUE);

	cout << "Done with PhotonQA_Runwise" << endl;
	cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

	return;

}//end

