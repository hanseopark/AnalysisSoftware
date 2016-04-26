#include "QA.h"

void ClusterQA_Compare(TString suffix = "eps", Int_t mode = 2)
{
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "ClusterQA_Compare" << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

	const Int_t maxHistosRange = 25;
	const Int_t maxHistosRunwiseRange = 50;
	Bool_t doPlotsByTrigger = kTRUE;

    gROOT->Reset();

    StyleSettingsThesis();
    SetPlotStyle();

//choose which data sets to process
//**************************************************************************************************************

    const Int_t nSets = 4;
    Size_t constMarkerSize = 1;
    TString outputDir ="LHC13g";
    TString fEnergyFlag = "2.76TeV";
    TString DataSets[nSets]={"LHC13g","LHC13g-kEMC7","LHC13g-kEMCEG1","LHC13g-kEMCEG2"};
    TString plotDataSet = "LHC13g";
    TString cutsDataSets[nSets]={
        "00000113_00200009327000008250400000_1111121063032230000_0163103100000010",
        "00052013_00200009327000008250400000_1111121063032230000_0163103100000010",
        "00083013_00200009327000008250400000_1111121063032230000_0163103100000010",
        "00085013_00200009327000008250400000_1111121063032230000_0163103100000010"
        };
    mode = 2;

	//**************************************************************************************************************

//	const Int_t nSets = 4;
//	Size_t constMarkerSize = 1.5;
//	TString outputDir ="LHC12";
//	TString fEnergyFlag = "8TeV";
//	TString DataSets[nSets]={"LHC12","LHC12-kEMC7","LHC12-kEMCEGA","LHC12-kEMCEJE"};
//	TString plotDataSet = "LHC12";
//	TString cutsDataSets[nSets]={
//		"0000011_00200009327000008250400000_10000063032230000_0163103100000000",
//		"0005211_00200009327000008250400000_10000063032230000_0163103100000000",
//		"0008111_00200009327000008250400000_10000063032230000_0163103100000000",
//		"0009111_00200009327000008250400000_10000063032230000_0163103100000000"
//		};

//**************************************************************************************************************

//	const Int_t nSets = 8;
//	doPlotsByTrigger = kFALSE;
//	Size_t constMarkerSize = 1.5;
//	TString outputDir ="LHC12x";
//	TString fEnergyFlag = "8TeV";
//	TString DataSets[nSets]={"LHC12a","LHC12b","LHC12c","LHC12d","LHC12f","LHC12g","LHC12h","LHC12i"};
//	TString plotDataSet = "LHC12";
//	TString cutsDataSets[nSets]={
//		"0000011_00200009327000008250400000_11111063032230000_0163103100000000",
//		"0000011_00200009327000008250400000_11111063032230000_0163103100000000",
//		"0000011_00200009327000008250400000_11111063032230000_0163103100000000",
//		"0000011_00200009327000008250400000_11111063032230000_0163103100000000",
//		"0000011_00200009327000008250400000_11111063032230000_0163103100000000",
//		"0000011_00200009327000008250400000_11111063032230000_0163103100000000",
//		"0000011_00200009327000008250400000_11111063032230000_0163103100000000",
//		"0000011_00200009327000008250400000_11111063032230000_0163103100000000"
//		};


//**************************************************************************************************************

//    const Int_t nSets = 2;
//    Size_t constMarkerSize = 1;
//    TString outputDir ="LHC11a";
//    TString fEnergyFlag = "2.76TeV";
//    TString DataSets[nSets]={"LHC11a_p4_wSDD","LHC11a_p4_wSDD-kEMC1"};
//    TString plotDataSet = "LHC11a";
//    TString cutsDataSets[nSets]={
//        "00003113_1111121053032220000_0163103100000050",
//        "00051013_1111121053032220000_0163103100000050"
//    };
//    mode = 4;

//**************************************************************************************************************
	Style_t hMarkerStyle[nSets];
	Size_t hMarkerSize[nSets];
	Color_t hMarkerColor[nSets];
	Color_t hLineColor[nSets];

	TString plotDataSets[nSets];
    TString pathDataSets[nSets];
	Int_t markerStyles[14]={2,4,5,20,21,22,23,24,25,26,27,28,29,30};

    for(Int_t i=0; i<nSets; i++)
    {
        pathDataSets[i]=Form("%s/%s/ClusterQA", cutsDataSets[i].Data(),fEnergyFlag.Data());

		if(doPlotsByTrigger){
			TString tempCutSelection = cutsDataSets[i];
			TString fTriggerCut = tempCutSelection(3,2);
			TString fTrigger = AnalyseSpecialTriggerCut(fTriggerCut.Atoi(), DataSets[i]);
			plotDataSets[i] = fTrigger;
			hMarkerStyle[i] = GetDefaultTriggerMarker(DataSets[i],fTriggerCut.Atoi());
			hMarkerSize[i] = constMarkerSize;
			hMarkerColor[i] = GetDefaultTriggerColor(DataSets[i],fTriggerCut.Atoi());
			hLineColor[i] = GetDefaultTriggerColor(DataSets[i],fTriggerCut.Atoi());
			cout << "Found Trigger: " << fTrigger.Data() << " determined from cut: " << fTriggerCut.Atoi() << endl;
		}else{
			plotDataSets[i] = DataSets[i];
			hMarkerStyle[i] = markerStyles[i%14];
			hMarkerSize[i] = constMarkerSize;
			hMarkerColor[i] = 1+i;
			hLineColor[i] = 1+i;
		}
    }

	std::vector<TString> vecDataSet;
	std::vector<TString> vecPlotDataSet;
	std::vector<TString> vecPathDataSet;
    for(Int_t i=0; i<nSets; i++){vecDataSet.push_back(DataSets[i].Data());}
    for(Int_t i=0; i<nSets; i++){vecPlotDataSet.push_back(plotDataSets[i].Data());}
    for(Int_t i=0; i<nSets; i++){vecPathDataSet.push_back(pathDataSets[i].Data());}

	TString fCutSelection[nSets];
	TString fEventCutSelection[nSets], fGammaCutSelection[nSets], fClusterCutSelection[nSets], fElectronCutSelection[nSets], fMesonCutSelection[nSets];
	TString calo;
	Int_t iCalo = 0;

	for(Int_t i=0; i<nSets; i++){
		fCutSelection[i] = cutsDataSets[i];
		ReturnSeparatedCutNumberAdvanced(fCutSelection[i], fEventCutSelection[i], fGammaCutSelection[i], fClusterCutSelection[i], fElectronCutSelection[i], fMesonCutSelection[i], mode);
		if(fClusterCutSelection[i].BeginsWith('1')){
			if(iCalo == 0){ calo="EMCal"; iCalo=1;}
			else if(iCalo != 1) {cout << "ERROR: Different calorimeters compared, returning..." << endl; return;}
		}else if(fClusterCutSelection[i].BeginsWith('2')){
			if(iCalo == 0){ calo="PHOS"; iCalo=2;}
			else if(iCalo != 2) {cout << "ERROR: Different calorimeters compared, returning..." << endl; return;}
		}else {cout << "No correct calorimeter type found: " << calo.Data() << ", returning..." << endl; return;}
	}

	cout << "Analyzing " << calo.Data() << " clusters..." << endl;

    TString fCollisionSystem = ReturnFullCollisionsSystem(fEnergyFlag);
    if (fCollisionSystem.CompareTo("") == 0){
        cout << "No correct collision system specification, has been given" << endl;
        return;
    }

    TString fDetectionProcess = ReturnFullTextReconstructionProcess(mode);

    std::vector<TH1D*> vecHistos[nSets];
	std::vector<TString> vecHistosName[nSets];
	std::vector<TString> vecHistosNameForSaving[nSets];

    std::vector<TH1D*> vecHistosRunwise[nSets];
	std::vector<TString> vecHistosRunwiseName[nSets];

	TH1D* hcDeltaEtaMatch[nSets];
	TH1D* hcDeltaPhiMatch[nSets];
	TH1D* hcDeltaEta[nSets];
	TH1D* hcDeltaPhi[nSets];
    TH1D* hcDispersion[nSets];
    TH1D* hcClusterEnergy[nSets];
    TH1D* hcM02[nSets];
    TH1D* hcM20[nSets];
    TH1D* hcNCells[nSets];
	TH1D* hcTimeBefore[nSets];
	TH1D* hcTime[nSets];

	TH1D* hcMotherPi0Conv_Phi[nSets];
	TH1D* hcMotherPi0Conv_Eta[nSets];
	TH1D* hcMotherEtaConv_Phi[nSets];
	TH1D* hcMotherEtaConv_Eta[nSets];

//---------
    TH1D* hNEvents[nSets];

//    TH1D* hTracksMeanGood[nSets];
//    TH1D* hTracksRMSGood[nSets];

//    TH1D* hConvNCandidates[nSets];
//    TH1D* hConvNCandidatesQA[nSets];

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

	Int_t nRange = 0;
	Int_t nRangeRunwise = 0;

    for(Int_t i=0; i<nSets; i++)
    {
		vecHistos[i].push_back(hcDeltaEtaMatch[i]);
		vecHistosName[i].push_back("fHistClusterdEtaMatchedAfterQA");
		vecHistosNameForSaving[i].push_back("dEtaMatched");
		if(i==0)nRange++;

		vecHistos[i].push_back(hcDeltaPhiMatch[i]);
		vecHistosName[i].push_back("fHistClusterdPhiMatchedAfterQA");
		vecHistosNameForSaving[i].push_back("dPhiMatched");
		if(i==0)nRange++;

        vecHistos[i].push_back(hcDeltaEta[i]);
		vecHistosName[i].push_back("fHistClusterdEtaBeforeQA");
		vecHistosNameForSaving[i].push_back("dEta");
		if(i==0)nRange++;

        vecHistos[i].push_back(hcDeltaPhi[i]);
		vecHistosName[i].push_back("fHistClusterdPhiBeforeQA");
		vecHistosNameForSaving[i].push_back("dPhi");
		if(i==0)nRange++;

        vecHistos[i].push_back(hcDispersion[i]);
		vecHistosName[i].push_back(Form("Dispersion_afterClusterQA %s",fClusterCutSelection[i].Data()));
		vecHistosNameForSaving[i].push_back("Dispersion");
		if(i==0)nRange++;

        vecHistos[i].push_back(hcClusterEnergy[i]);
		vecHistosName[i].push_back("energyAfter");
		vecHistosNameForSaving[i].push_back("ClusterEnergy");
		if(i==0)nRange++;

        vecHistos[i].push_back(hcM02[i]);
		vecHistosName[i].push_back(Form("M02_afterClusterQA %s",fClusterCutSelection[i].Data()));
		vecHistosNameForSaving[i].push_back("M02");
		if(i==0)nRange++;

        vecHistos[i].push_back(hcM20[i]);
		vecHistosName[i].push_back(Form("M20_afterClusterQA %s",fClusterCutSelection[i].Data()));
		vecHistosNameForSaving[i].push_back("M20");
		if(i==0)nRange++;

        vecHistos[i].push_back(hcNCells[i]);
		vecHistosName[i].push_back(Form("NCellPerCluster_afterClusterQA %s",fClusterCutSelection[i].Data()));
		vecHistosNameForSaving[i].push_back("NCells");
		if(i==0)nRange++;

		vecHistos[i].push_back(hcTimeBefore[i]);
		vecHistosName[i].push_back("fHistClusterTimeBeforeQA");
		vecHistosNameForSaving[i].push_back("Time_Before");
		if(i==0)nRange++;

		vecHistos[i].push_back(hcTime[i]);
		vecHistosName[i].push_back("fHistClusterTimeAfterQA");
		vecHistosNameForSaving[i].push_back("Time");
		if(i==0)nRange++;

        if(mode != 4){
          vecHistos[i].push_back(hcMotherPi0Conv_Phi[i]);
          vecHistosName[i].push_back("MotherPi0Conv_Phi");
          vecHistosNameForSaving[i].push_back("MotherPi0Conv_Phi");
          if(i==0)nRange++;

          vecHistos[i].push_back(hcMotherPi0Conv_Eta[i]);
          vecHistosName[i].push_back("MotherPi0Conv_Eta");
          vecHistosNameForSaving[i].push_back("MotherPi0Conv_Eta");
          if(i==0)nRange++;

          vecHistos[i].push_back(hcMotherEtaConv_Phi[i]);
          vecHistosName[i].push_back("MotherEtaConv_Phi");
          vecHistosNameForSaving[i].push_back("MotherEtaConv_Phi");
          if(i==0)nRange++;

          vecHistos[i].push_back(hcMotherEtaConv_Eta[i]);
          vecHistosName[i].push_back("MotherEtaConv_Eta");
          vecHistosNameForSaving[i].push_back("MotherEtaConv_Eta");
          if(i==0)nRange++;
        }

//		vecHistos[i].push_back(hClusterTime[i]);
//		vecHistosName[i].push_back("MotherEtaConv_Eta");
//---------
        vecHistosRunwise[i].push_back(hNEvents[i]);
		vecHistosRunwiseName[i].push_back("hNEvents");
		if(i==0)nRangeRunwise++;

//        vecHistosRunwise[i].push_back(hTracksMeanGood[i]);
//		vecHistosRunwiseName[i].push_back("hTracksGood-Mean");
//		if(i==0)nRangeRunwise++;

//        vecHistosRunwise[i].push_back(hTracksRMSGood[i]);
//		vecHistosRunwiseName[i].push_back("hTracksGood-RMS");
//		if(i==0)nRangeRunwise++;

//        vecHistosRunwise[i].push_back(hConvNCandidates[i]);
//		vecHistosRunwiseName[i].push_back("hConvNCandidates");
//		if(i==0)nRangeRunwise++;

//        vecHistosRunwise[i].push_back(hConvNCandidatesQA[i]);
//		vecHistosRunwiseName[i].push_back("hConvNCandidatesQA");
//		if(i==0)nRangeRunwise++;

		vecHistosRunwise[i].push_back(hCaloNClusters[i]);
		vecHistosRunwiseName[i].push_back("hCaloNClusters");
		if(i==0)nRangeRunwise++;

		vecHistosRunwise[i].push_back(hCaloNClustersQA[i]);
		vecHistosRunwiseName[i].push_back("hCaloNClustersQA");
		if(i==0)nRangeRunwise++;

        vecHistosRunwise[i].push_back(hClusterMeanEnergy[i]);
		vecHistosRunwiseName[i].push_back("hClusterEnergy-Mean");
		if(i==0)nRangeRunwise++;

        vecHistosRunwise[i].push_back(hClusterRMSEnergy[i]);
		vecHistosRunwiseName[i].push_back("hClusterEnergy-RMS");
		if(i==0)nRangeRunwise++;

        vecHistosRunwise[i].push_back(hClusterEnergy01[i]);
		vecHistosRunwiseName[i].push_back("hClusterEnergy-01");
		if(i==0)nRangeRunwise++;

        vecHistosRunwise[i].push_back(hClusterEnergy14[i]);
		vecHistosRunwiseName[i].push_back("hClusterEnergy-14");
		if(i==0)nRangeRunwise++;

        vecHistosRunwise[i].push_back(hClusterEnergy4[i]);
		vecHistosRunwiseName[i].push_back("hClusterEnergy-4");
		if(i==0)nRangeRunwise++;

        vecHistosRunwise[i].push_back(hClusterMeanNCells[i]);
		vecHistosRunwiseName[i].push_back("hClusterNCells-Mean");
		if(i==0)nRangeRunwise++;

        vecHistosRunwise[i].push_back(hClusterRMSNCells[i]);
		vecHistosRunwiseName[i].push_back("hClusterNCells-RMS");
		if(i==0)nRangeRunwise++;

        vecHistosRunwise[i].push_back(hClusterMeanDispersion[i]);
		vecHistosRunwiseName[i].push_back("hClusterDispersion-Mean");
		if(i==0)nRangeRunwise++;

        vecHistosRunwise[i].push_back(hClusterRMSDispersion[i]);
		vecHistosRunwiseName[i].push_back("hClusterDispersion-RMS");
		if(i==0)nRangeRunwise++;

        vecHistosRunwise[i].push_back(hClusterMeanM02[i]);
		vecHistosRunwiseName[i].push_back("hClusterM02-Mean");
		if(i==0)nRangeRunwise++;

        vecHistosRunwise[i].push_back(hClusterRMSM02[i]);
		vecHistosRunwiseName[i].push_back("hClusterM02-RMS");
		if(i==0)nRangeRunwise++;

        vecHistosRunwise[i].push_back(hClusterMeanR[i]);
		vecHistosRunwiseName[i].push_back("hClusterR-Mean");
		if(i==0)nRangeRunwise++;

        vecHistosRunwise[i].push_back(hClusterRMSR[i]);
		vecHistosRunwiseName[i].push_back("hClusterR-RMS");
		if(i==0)nRangeRunwise++;

        vecHistosRunwise[i].push_back(hClusterMeanTime[i]);
		vecHistosRunwiseName[i].push_back("hClusterTime-Mean");
		if(i==0)nRangeRunwise++;

        vecHistosRunwise[i].push_back(hClusterRMSTime[i]);
		vecHistosRunwiseName[i].push_back("hClusterTime-RMS");
		if(i==0)nRangeRunwise++;

        vecHistosRunwise[i].push_back(hClusterFractionMatches[i]);
		vecHistosRunwiseName[i].push_back("hCluster-FractionMatches");
		if(i==0)nRangeRunwise++;

		vecHistosRunwise[i].push_back(hClusterFractionMatchesS[i]);
		vecHistosRunwiseName[i].push_back("hCluster-FractionMatchesS");
		if(i==0)nRangeRunwise++;

		vecHistosRunwise[i].push_back(hClusterFractionMatchesM[i]);
		vecHistosRunwiseName[i].push_back("hCluster-FractionMatchesM");
		if(i==0)nRangeRunwise++;

		vecHistosRunwise[i].push_back(hClusterFractionMatchesH[i]);
		vecHistosRunwiseName[i].push_back("hCluster-FractionMatchesH");
		if(i==0)nRangeRunwise++;

        vecHistosRunwise[i].push_back(hClusterMeanDeltaEta[i]);
		vecHistosRunwiseName[i].push_back("hClusterDeltaEta-Mean");
		if(i==0)nRangeRunwise++;

        vecHistosRunwise[i].push_back(hClusterMeanDeltaPhi[i]);
		vecHistosRunwiseName[i].push_back("hClusterDeltaPhi-Mean");
		if(i==0)nRangeRunwise++;

        vecHistosRunwise[i].push_back(hClusterRMSDeltaEta[i]);
		vecHistosRunwiseName[i].push_back("hClusterDeltaEta-RMS");
		if(i==0)nRangeRunwise++;

        vecHistosRunwise[i].push_back(hClusterRMSDeltaPhi[i]);
		vecHistosRunwiseName[i].push_back("hClusterDeltaPhi-RMS");
		if(i==0)nRangeRunwise++;
    }

	Double_t nEvents[nSets];

    for(Int_t i=0; i<nSets; i++)
    {
        TString fFile = Form("%s/ClusterQA_%s.root", pathDataSets[i].Data(), DataSets[i].Data());
		TString fFileRunwise = Form("%s/%s_ClusterQARunwise.root", pathDataSets[i].Data(), DataSets[i].Data());
        TFile* File = new TFile(fFile.Data(),"READ");
        TFile* FileRunwise = new TFile(fFileRunwise.Data(),"READ");

		if(File->IsZombie()) {cout << "Warning: ROOT file '" << fFile.Data() << "' could not be openend!" << endl;}
		if(FileRunwise->IsZombie()) {cout << "Warning: ROOT file '" << fFileRunwise.Data() << "' could not be openend!" << endl;}
		if(File->IsZombie()&&FileRunwise->IsZombie()) {cout << "ERROR: Only Zombies around, return..." << endl; return;}

		nEvents[i] = GetNEvents((TH1*) File->Get("NEvents"),kFALSE);

        cout << endl;
        cout << "\t\t----------------------------------------------------------------------------" << endl;
        cout << "\t\tDataSet: " << DataSets[i].Data() << ", " << plotDataSets[i].Data() << ":" << endl;
        cout << "\t\tTotal NEvents: " << nEvents[i] << endl;
        cout << "\t\tProcessing file: " << fFile.Data() << endl;
        cout << "\t\tProcessing runwise file: " << fFileRunwise.Data() << endl;
        cout << "\t\t----------------------------------------------------------------------------" << endl;
        cout << endl;

		cout << "\t\t----------------------------------------------------------------------------" << endl;
		cout << "\t\t\tNumber of Runwise Histograms being processed: " << nRangeRunwise
			 << " of max. possible: " << maxHistosRunwiseRange << endl;
		cout << "\t\t----------------------------------------------------------------------------" << endl;

		for(Int_t j=0; j<nRangeRunwise; j++)
        {
			cout << Form("%s_%s", vecHistosRunwiseName[i].at(j).Data(), DataSets[i].Data());
			TH1D* temp = (TH1D*) FileRunwise->Get(Form("%s_%s", vecHistosRunwiseName[i].at(j).Data(), DataSets[i].Data()));
            cout << " - "<< temp << ", ";
			OnlyEditTH1(temp, hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
			temp->SetTitle("");
            vecHistosRunwise[i].at(j) = new TH1D(*temp);
        }

		cout << "\n\t\t----------------------------------------------------------------------------" << endl;
		cout << "\t\t\tNumber of Histograms being processed: " << nRange
			 << " of max. possible: " << maxHistosRange << endl;
		cout << "\t\t----------------------------------------------------------------------------" << endl;

		for(Int_t j=0; j<nRange; j++)
        {
			cout << ((TString)vecHistosName[i].at(j)).Data();
			TH1D* temp = (TH1D*) File->Get(((TString)vecHistosName[i].at(j)).Data());
            cout << " - "<< temp << ", ";
			temp->Sumw2();
			if( j >= 4 ) {
				temp->Scale(1./nEvents[i]);
				temp->GetYaxis()->SetTitle(Form("#frac{1}{N_{Events}} %s",temp->GetYaxis()->GetTitle()));
			}
			else if( j == 0 || j == 1 || j == 2 || j == 3){
				Double_t tempEntries = temp->GetEntries();
				temp->Scale(1./tempEntries);
				temp->GetYaxis()->SetTitle(Form("#frac{1}{N} %s",temp->GetYaxis()->GetTitle()));
			}

			OnlyEditTH1(temp, hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
			if(j==9) temp->GetXaxis()->SetNoExponent(kFALSE);
			temp->GetXaxis()->SetTitleOffset(1);
			temp->GetYaxis()->SetTitleOffset(1.1);
			temp->SetTitle("");
            vecHistos[i].at(j) = new TH1D(*temp);
        }

        cout << "done!" << endl;
    }

	gSystem->Exec("mkdir -p ClusterQA_Compare/"+outputDir);
	gSystem->Exec("mkdir -p ClusterQA_Compare/"+outputDir+"/QA");

	TString nameOutput = Form("ClusterQA_Compare/%s/ClusterQA_Compare.root",outputDir.Data());
    TFile* fOutput = new TFile(nameOutput,"RECREATE");
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "Output file: " << nameOutput << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

	Int_t nTopRows = (nSets-1)/(3) + 1;
	TCanvas* canvas = new TCanvas("canvas","",200,10,1350,(846.+nTopRows*54.));  // gives the page size
	Double_t leftMar = 0.09; Double_t rightMar = 0.02; Double_t topMargin = (nTopRows*54.)/(846.+nTopRows*54.); Double_t bottomMargin = 0.09;
    DrawGammaCanvasSettings(canvas, leftMar, rightMar, topMargin, bottomMargin);

//****************************** Combined Trending Histograms ************************************************

	TLegend *legend = new TLegend(0.12,0.99 - (nTopRows*0.04),0.95,0.98);
	legend->SetNColumns(3);
    legend->SetFillColor(0);
    legend->SetLineColor(0);
    legend->SetTextSize(0.03);
    legend->SetTextFont(42);
//---------
	if(doPlotsByTrigger){
		for(Int_t h=0; h<(Int_t)vecHistosRunwise[0].size(); h++)
		{
			cout << h << ", ";
			AdjustHistRange(vecHistosRunwise,1.2,1.2,h,nSets,kTRUE);
			for(Int_t i=0; i<nSets; i++)
			{
				TString draw;
				if(h==0) draw = (i==0)?"p":"p, same";
				else draw = (i==0)?"px0e1":"px0e1, same";
				((TH1D*) vecHistosRunwise[i].at(h))->Draw(draw.Data());
				legend->AddEntry(((TH1D*) vecHistosRunwise[i].at(h)),plotDataSets[i].Data(),"p");
			}
			legend->Draw();
			PutProcessLabelAndEnergyOnPlot(0.8, 0.97-nTopRows*0.06, 0.03, fCollisionSystem.Data(), Form("%s clusters", calo.Data()), plotDataSet.Data());
			if(h==0) SaveWriteCanvas(canvas, Form("ClusterQA_Compare/%s/%s.%s", outputDir.Data(), vecHistosRunwiseName[0].at(h).Data(),suffix.Data()), kFALSE, kTRUE);
			else SaveWriteCanvas(canvas, Form("ClusterQA_Compare/%s/%s.%s", outputDir.Data(), vecHistosRunwiseName[0].at(h).Data(),suffix.Data()));
			legend->Clear();
		}
		cout << endl;
	}

//****************************** Combined Comparison Histograms ************************************************
	for(Int_t h=0; h<(Int_t)vecHistos[0].size(); h++)
    {
		cout << h << ", ";
		if(h==9) TGaxis::SetExponentOffset(0.5, 0, "x");
		AdjustHistRange(vecHistos,5,5,h,nSets,kFALSE);
        for(Int_t i=0; i<nSets; i++)
        {
            TString draw;
            if(h==0) draw = (i==0)?"p":"p, same";
            else draw = (i==0)?"px0e1":"px0e1, same";
			((TH1D*) vecHistos[i].at(h))->Draw(draw.Data());
			legend->AddEntry(((TH1D*) vecHistos[i].at(h)),plotDataSets[i].Data(),"p");
        }
        legend->Draw();
		PutProcessLabelAndEnergyOnPlot(0.8, 0.97-nTopRows*0.06, 0.03, fCollisionSystem.Data(), Form("%s clusters", calo.Data()), plotDataSet.Data());
		SaveWriteCanvas(canvas, Form("ClusterQA_Compare/%s/QA/%s.%s", outputDir.Data(), vecHistosNameForSaving[0].at(h).Data(),suffix.Data()), kFALSE, kTRUE);
        legend->Clear();
		if(h==9) TGaxis::SetExponentOffset(0, 0, "x");
    }
	cout << endl;

    fOutput->Write();
    fOutput->Close();

    cout << "Done with ClusterQA_Compare" << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;


    return;
}
