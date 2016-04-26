#include "QA.h"
#include "../TaskV1/BuildHistogramsForGammaQAAdvV3.C"

void PhotonQA(
                Int_t nSetsIn,
				TString fEnergyFlag,
				TString* DataSets,
				TString* plotDataSets,
				TString* pathDataSets,
				Int_t mode = 2,
				Int_t cutNr = -1,			// if -1: you have to choose number at runtime
				Int_t doExtQA = 2,			// 0: switched off, 1: normal extQA, 2: with Cell level plots
				TString suffix = "eps",
				TString labelData = "Data",
				Bool_t addSubfolder = kFALSE
				)
{
	cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
	cout << "PhotonQA" << endl;
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
		cout << "Maximum hardcoded number of Data Sets in PhotonQA: " << maxSets << endl;
		cout << "You have chosen: " << nSets << ", returning!" << endl;
		return;
	}

    Color_t colorCompare[maxSets] = {kBlack, kRed+1, kMagenta+2, 807, 800, kGreen+2, kCyan+2, kBlue+1, kOrange+2, kAzure, kViolet, kGray+1};
	TString nameMainDir[maxSets];
	TString nameCutsPQA[maxSets];
	TString nameCutsPQAshort[maxSets];

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

	for(Int_t i=0; i<nSets; i++){
		TFile* fFile = new TFile(pathDataSets[i].Data(),"READ");
		if(fFile->IsZombie()){cout << "ERROR: File " << pathDataSets[i].Data() << " could not be openend! Returning..." << endl; return;}
		else{
			cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
			cout << "Processing file: " << pathDataSets[i].Data();
			TKey *key;
			TIter next(fFile->GetListOfKeys());
			while ((key=(TKey*)next())){
				cout << Form(" - found TopDir: %s",key->GetName());
				nameMainDir[i] = key->GetName();
			}
			nameCutsPQA[i] = Form("%s",nameMainDir[i].Data());
			nameCutsPQAshort[i] = Form("%s",nameMainDir[i].Data());
			nameCutsPQA[i].Replace(0,15,"");
			nameCutsPQAshort[i].Replace(0,15,"");
			nameCutsPQAshort[i].Replace(8,35,"");

			cout << endl;
			cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
			cout << "nameMainDir:" << nameMainDir[i].Data() << endl;
			cout << "long cutnumber for PhotonQA: " << nameCutsPQA[i] << "\nshort cutnumber for PhotonQA: "<< nameCutsPQAshort[i] << endl;
			cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
			cout << endl;
			if(nameMainDir[i].IsNull() || !nameMainDir[i].BeginsWith("Gamma")){cout << "ERROR, Unable to obtain valid name of MainDir:|" << nameMainDir[i].Data() << "|, running in mode: " << fMode << endl; return;}
		}
		fFile->Close();
		delete fFile;
	}

//*****************************************************************************************************
//*****************************************************************************************************

	cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
	cout << "Pictures are saved as " << suffix.Data() << "!" << endl;
	cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

	cout << "Obtaining trigger - ";
    TString* fTrigger = new TString[nSets];
    for(Int_t iT=0;iT<nSets;iT++){fTrigger[iT]="";}
    TString fTriggerCut = nameCutsPQAshort[0](3,2);
    fTrigger[0] = AnalyseSpecialTriggerCut(fTriggerCut.Atoi(), DataSets[0].Data());
    cout  << "'" << fTrigger[0].Data() << "' - was found!" << endl;
    if(fTrigger[0].Contains("not defined")){
        fTrigger[0] = "";
        cout << "INFO: Trigger cut not defined!" << endl;
    }

	cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
	TString fCollisionSystem = ReturnFullCollisionsSystem(fEnergyFlag);
	if (fCollisionSystem.CompareTo("") == 0){
		cout << "No correct collision system specification, has been given" << endl;
		return;
	}

	TString fDetectionProcess = ReturnFullTextReconstructionProcess(fMode);

	TString outputDir = Form("%s/%s/PhotonQA/%s",nameCutsPQA[0].Data(),fEnergyFlag.Data(),suffix.Data());
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
	for(Int_t i=0; i<nSets; i++){
		fLog << "Using file: " << pathDataSets[i].Data() << endl;
		fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
		fLog << DataSets[i].Data() << endl;
		fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
		fLog << "processed: " << nameCutsPQA[i].Data() << endl;
	}
	fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
	fLog << "Energy: " << fEnergyFlag.Data() << endl;
	fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
	fLog << fCollisionSystem.Data() << endl;
	fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

//*****************************************************************************************************
//*****************************************************************************************************
//****************************** Vectors for Histograms ***********************************************
//*****************************************************************************************************
//*****************************************************************************************************

	std::vector<TH1D*> vecGammaPt;
	std::vector<TH1D*> vecGammaEta;
	std::vector<TH1D*> vecGammaAlpha;
	std::vector<TH1D*> vecGammaPhi;
	std::vector<TH1D*> vecGammaPhiEtaNeg;
	std::vector<TH1D*> vecGammaPhiEtaPos;

	std::vector<TH1D*> vecEta[2];
	std::vector<TH1D*> vecPt[2];
	std::vector<TH1D*> vecNSigmadEdx[2];
    std::vector<TH1D*> vecFCL[2];
    std::vector<TH1D*> vecTPCCL[2];
    std::vector<TH1D*> vecChi2[5];
    std::vector<TH1D*> vecPsi[5];
    std::vector<TH1D*> vecCos[5];
    std::vector<TH1D*> vecGamma[5];
    std::vector<TH1D*> vecInvMass[4];

    std::vector<TH1D*> vecITSAllclsR[2];
    std::vector<TH1D*> vecITSclsR[2][7];
    std::vector<TH1D*> vecTrueITSclsR[2][7];

    Double_t* nEvents = new Double_t[nSets];
    Int_t minB = 0; Int_t maxB = 0;
    Int_t minYB = 0; Int_t maxYB = 0;

	TString lepton[2] = {"Electron","Positron"};
	TString charge[2] = {"-","+"};

//*****************************************************************************************************
//*****************************************************************************************************
//****************************** Looping over DataSets ************************************************
//*****************************************************************************************************
//*****************************************************************************************************

    TCanvas* canvas = new TCanvas("canvas","",10,10,750,500);  // gives the page size
    TCanvas* cvsQuadratic = new TCanvas("cvsQuadratic","",10,10,500,500);  // gives the page size
    Double_t leftMargin = 0.09; Double_t rightMargin = 0.02; Double_t topMargin = 0.04; Double_t bottomMargin = 0.09;
	DrawGammaCanvasSettings(canvas,leftMargin,rightMargin,topMargin,bottomMargin);

	for(Int_t i=0; i<nSets; i++)
	{
	//-----------------------------------
		TString path = pathDataSets[i];
		if(!path.EndsWith(Form("/PhotonQA_%s.root",DataSets[i].Data()))){
			cout << "ERROR: Path to PhotonQA output does not end with '" << Form("/PhotonQA_%s.root",DataSets[i].Data()) << "' for DataSet '" << DataSets[i].Data() << "'" << endl;
			return;
		}

		cout << endl;
		cout << "\t\t----------------------------------------------------------------------------" << endl;
		cout << Form("\t\tSet %s", DataSets[i].Data()) << endl;
		cout << "\t\tProcessing file: " << path.Data() << endl;
		cout << "\t\t----------------------------------------------------------------------------" << endl;
		cout << endl;

		TFile* fPhotonQAFile = new TFile(path.Data(),"READ");
		if(fPhotonQAFile->IsZombie()) {
			cout << "ERROR: ROOT file '" << pathDataSets[i].Data() << "' could not be openend, returning!" << endl;
			fLog << "ERROR: ROOT file '" << pathDataSets[i].Data() << "' could not be openend, returning!" << endl;
			return;
		}

		TDirectory* TopDir = (TDirectory*) fPhotonQAFile->Get(Form("GammaConvV1_QA_%s",nameCutsPQA[i].Data()));
			if(TopDir == NULL) {cout << "ERROR: TopDir not Found"<<endl; return;}
	//-----------------------------------
		nEvents[i] = 0;
		TH1D* GOODESD = (TH1D*) TopDir->Get("histoGoodESDTracks");
		if(GOODESD) nEvents[i] = (Double_t) GOODESD->Integral();
	//-----------------------------------
		fLog << "----------------------------------------------------------------------------" << endl;
		fLog << "Processing file: " << pathDataSets[i].Data() << endl;
		fLog << "Set: " << plotDataSets[i].Data() << " - NEvents: " << nEvents[i] << endl;
		fLog << "----------------------------------------------------------------------------" << endl;
	//-----------------------------------
		if( nEvents[i] < 1. ){cout << "ERROR: number of accepted events in data set is <1: " << nEvents[i] << "! Returning..." << endl; return;}
	//-----------------------------------
		const char* nameOutput = Form("%s/%s/PhotonQA/PhotonQA_%s.root",nameCutsPQA[i].Data(),fEnergyFlag.Data(),DataSets[i].Data());
		TFile* fOutput = new TFile(nameOutput,"RECREATE");
		cout << endl;
		cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
		cout << "Output file: " << nameOutput << endl;
		cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
		cout << endl;
		fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
		fLog << "Output file: " << nameOutput << endl;
		fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;


	//-----------------------------------
	//-----------------------------------
	//-- Gammas
	//-----------------------------------
	//-----------------------------------

	//-----------------------------------
		Double_t nGammas = 1.;
		TH2D* fHistGammaAlphaQt = (TH2D*)TopDir->Get("histoGammaAlphaQt");
		if(fHistGammaAlphaQt){
            nGammas = fHistGammaAlphaQt->Integral();
			fHistGammaAlphaQt->Sumw2();
			fHistGammaAlphaQt->Scale(1./nGammas);
            SetZMinMaxTH2(fHistGammaAlphaQt,1,fHistGammaAlphaQt->GetNbinsX(),1,fHistGammaAlphaQt->GetNbinsY(),kTRUE);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                                fHistGammaAlphaQt,"",
                                 "#alpha = (#it{p}^{+}_{L}-#it{p}^{-}_{L})/(#it{p}^{+}_{L}+#it{p}^{-}_{L})","#it{q}_{T} (GeV/#it{c})",1,1.4,
                                 0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
			DrawArmenterosLines();
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistGammaAlphaQt, Form("%s/Armenteros_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));

            TH1D* fHistGammaAlpha = (TH1D*)fHistGammaAlphaQt->ProjectionX("fHistGammaAlpha",1,fHistGammaAlphaQt->GetNbinsY());
			GetMinMaxBin(fHistGammaAlpha,minB,maxB);
            SetXRange(fHistGammaAlpha,minB,maxB);
			DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kFALSE,kFALSE,
                                 fHistGammaAlpha,"","Photon #alpha","# of Entries",1,1,
                                 0.82,0.94,0.03,fCollisionSystem,plotDataSets[i],"");
            SaveCanvasAndWriteHistogram(canvas, fHistGammaAlpha, Form("%s/Photon_Alpha_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
			vecGammaAlpha.push_back(fHistGammaAlpha);
        }else cout << "INFO: Object |histoGammaAlphaQt| could not be found! Skipping Draw..." << endl;
    //-----------------------------------
        if(nGammas <= 1.){
          cout << "********************************************************" << endl;
          cout << "WARNING: nGammas <= 1! Scaling will not work properly..." << endl;
          cout << "********************************************************" << endl;
        }
	//-----------------------------------
		TH2D* fHistTruePrimGammaAlphaQt = (TH2D*)TopDir->Get("histoTruePrimGammaAlphaQt");
		if(fHistTruePrimGammaAlphaQt){
            CheckNEntries(fHistTruePrimGammaAlphaQt);
			fHistTruePrimGammaAlphaQt->Sumw2();
			fHistTruePrimGammaAlphaQt->Scale(1./nGammas);
            SetZMinMaxTH2(fHistTruePrimGammaAlphaQt,1,fHistTruePrimGammaAlphaQt->GetNbinsX(),1,fHistTruePrimGammaAlphaQt->GetNbinsY(),kTRUE);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                                fHistGammaAlphaQt,"",
                                 "#alpha = (#it{p}^{+}_{L}-#it{p}^{-}_{L})/(#it{p}^{+}_{L}+#it{p}^{-}_{L})", "#it{q}_{T} (GeV/#it{c})",1,1.4,
                                 0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
			DrawArmenterosLines();
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistTruePrimGammaAlphaQt, Form("%s/Armenteros_TruePrimGamma_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
        }else cout << "INFO: Object |histoTruePrimGammaAlphaQt| could not be found! Skipping Draw..." << endl;
	//-----------------------------------
		TH2D* fHistTrueSecGammaAlphaQt = (TH2D*)TopDir->Get("histoTrueSecGammaAlphaQt");
		if(fHistTrueSecGammaAlphaQt){
            CheckNEntries(fHistTrueSecGammaAlphaQt);
			fHistTrueSecGammaAlphaQt->Sumw2();
			fHistTrueSecGammaAlphaQt->Scale(1./nGammas);
            SetZMinMaxTH2(fHistTrueSecGammaAlphaQt,1,fHistTrueSecGammaAlphaQt->GetNbinsX(),1,fHistTrueSecGammaAlphaQt->GetNbinsY(),kTRUE);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                                fHistTrueSecGammaAlphaQt,"",
                                 "#alpha = (#it{p}^{+}_{L}-#it{p}^{-}_{L})/(#it{p}^{+}_{L}+#it{p}^{-}_{L})", "#it{q}_{T} (GeV/#it{c})",1,1.4,
                                 0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
			DrawArmenterosLines();
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistTrueSecGammaAlphaQt, Form("%s/Armenteros_TrueSecGamma_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
        }else cout << "INFO: Object |histoTrueSecGammaAlphaQt| could not be found! Skipping Draw..." << endl;
	//-----------------------------------
		TH2D* fHistTrueDalitzGammaAlphaQt = (TH2D*)TopDir->Get("histoTrueDalitzGammaAlphaQt");
		if(fHistTrueDalitzGammaAlphaQt){
            CheckNEntries(fHistTrueDalitzGammaAlphaQt);
			fHistTrueDalitzGammaAlphaQt->Sumw2();
			fHistTrueDalitzGammaAlphaQt->Scale(1./nGammas);
            SetZMinMaxTH2(fHistTrueDalitzGammaAlphaQt,1,fHistTrueDalitzGammaAlphaQt->GetNbinsX(),1,fHistTrueDalitzGammaAlphaQt->GetNbinsY(),kTRUE);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                                fHistTrueDalitzGammaAlphaQt,"",
                                 "#alpha = (#it{p}^{+}_{L}-#it{p}^{-}_{L})/(#it{p}^{+}_{L}+#it{p}^{-}_{L})", "#it{q}_{T} (GeV/#it{c})",1,1.4,
                                 0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
			DrawArmenterosLines();
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistTrueDalitzGammaAlphaQt, Form("%s/Armenteros_TrueDalitzGamma_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
        }else cout << "INFO: Object |histoTrueDalitzGammaAlphaQt| could not be found! Skipping Draw..." << endl;
	//-----------------------------------
		if(fHistGammaAlphaQt && fHistTruePrimGammaAlphaQt && fHistTrueSecGammaAlphaQt && fHistTrueDalitzGammaAlphaQt){
			TH2D* fHistTrueBGAlphaQt = (TH2D*) fHistGammaAlphaQt->Clone("histoMonteCarloBGAlphaQtMonteCarlo");
            fHistTrueBGAlphaQt->Sumw2();
			fHistTrueBGAlphaQt->Add(fHistTruePrimGammaAlphaQt,-1);
			fHistTrueBGAlphaQt->Add(fHistTrueSecGammaAlphaQt,-1);
			fHistTrueBGAlphaQt->Add(fHistTrueDalitzGammaAlphaQt,-1);
			fHistTrueBGAlphaQt->Scale(1./nGammas);
            SetZMinMaxTH2(fHistTrueBGAlphaQt,1,fHistTrueBGAlphaQt->GetNbinsX(),1,fHistTrueBGAlphaQt->GetNbinsY(),kTRUE);
            CheckForNegativeEntries(fHistTrueBGAlphaQt);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                                fHistTrueBGAlphaQt,"",
                                 "#alpha = (#it{p}^{+}_{L}-#it{p}^{-}_{L})/(#it{p}^{+}_{L}+#it{p}^{-}_{L})", "#it{q}_{T} (GeV/#it{c})",1,1.4,
                                 0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
			DrawArmenterosLines();
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistTrueBGAlphaQt, Form("%s/Armenteros_BG_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
			delete fHistTrueBGAlphaQt;
		}
	//-----------------------------------
		TH2D* fHistGammaEtaPt = (TH2D*)TopDir->Get("histoGammaEtaPt");
		if(fHistGammaEtaPt){
            TH1D* fHistGammaPt = (TH1D*)fHistGammaEtaPt->ProjectionY("fHistGammaPt",1,fHistGammaEtaPt->GetNbinsX());
            fHistGammaPt->Sumw2();
            fHistGammaPt->Scale(1./nGammas);
            GetMinMaxBin(fHistGammaPt,minB,maxB);
            SetXRange(fHistGammaPt,fHistGammaPt->GetXaxis()->FindBin(0.1),maxB+5);
			DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kTRUE,kTRUE,kFALSE,
                                 fHistGammaPt,"","#it{p}_{T} (GeV/#it{c})","#frac{1}{N_{#gamma}} #frac{dN}{d#it{p}_{T}}",1,1,
                                 0.82,0.94,0.03,fCollisionSystem,plotDataSets[i],"");
            SaveCanvasAndWriteHistogram(canvas, fHistGammaPt, Form("%s/Photon_Pt_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
			vecGammaPt.push_back(fHistGammaPt);
			//-----------------------------------
            TH1D* fHistGammaEta = (TH1D*)fHistGammaEtaPt->ProjectionX("fHistGammaEta",1,fHistGammaEtaPt->GetNbinsY());
            fHistGammaEta->Sumw2();
            fHistGammaEta->Scale(1./nGammas);
			GetMinMaxBin(fHistGammaEta,minB,maxB);
            SetXRange(fHistGammaEta,minB-5,maxB+5);
			DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kFALSE,kFALSE,
                                 fHistGammaEta,"","Photon #eta","#frac{1}{N_{#gamma}} #frac{dN}{d#eta}",1,1,
                                 0.82,0.94,0.03,fCollisionSystem,plotDataSets[i],"");
            SaveCanvasAndWriteHistogram(canvas, fHistGammaEta, Form("%s/Photon_Eta_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
			vecGammaEta.push_back(fHistGammaEta);
        }else cout << "INFO: Object |histoGammaEtaPt| could not be found! Skipping Draw..." << endl;
	//-----------------------------------
		TH1D* fHistGammaPhi = (TH1D*)TopDir->Get("histoGammaPhi");
		if(fHistGammaPhi){
            fHistGammaPhi->Sumw2();
            fHistGammaPhi->Scale(1./nGammas);
			GetMinMaxBin(fHistGammaPhi,minB,maxB);
            SetXRange(fHistGammaPhi,minB,maxB);
			DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kFALSE,kFALSE,
                                 fHistGammaPhi,"","Photon #phi","#frac{1}{N_{#gamma}} #frac{dN}{d#phi}",1,1,
                                 0.82,0.94,0.03,fCollisionSystem,plotDataSets[i],"");
            SaveCanvasAndWriteHistogram(canvas, fHistGammaPhi, Form("%s/Photon_Phi_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
			vecGammaPhi.push_back(new TH1D(*fHistGammaPhi));
        }else cout << "INFO: Object |histoGammaPhi| could not be found! Skipping Draw..." << endl;
	//-----------------------------------
	TH1D* fHistGammaPhiEtaNeg = (TH1D*)TopDir->Get("histoGammaPhiEtaNeg");
	if(fHistGammaPhiEtaNeg){
            fHistGammaPhiEtaNeg->Sumw2();
            fHistGammaPhiEtaNeg->Scale(1./nGammas);
			GetMinMaxBin(fHistGammaPhiEtaNeg,minB,maxB);
            SetXRange(fHistGammaPhiEtaNeg,minB,maxB);
			DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kFALSE,kFALSE,
                                 fHistGammaPhiEtaNeg,"","Photon #phi C side","#frac{1}{N_{#gamma}} #frac{dN}{d#phi}",1,1,
                                 0.82,0.94,0.03,fCollisionSystem,plotDataSets[i],"");
            SaveCanvasAndWriteHistogram(canvas, fHistGammaPhiEtaNeg, Form("%s/Photon_Phi_EtaNeg_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
			vecGammaPhiEtaNeg.push_back(new TH1D(*fHistGammaPhiEtaNeg));
        }else cout << "INFO: Object |histoGammaPhiEtaNeg| could not be found! Skipping Draw..." << endl;
	//-----------------------------------
	TH1D* fHistGammaPhiEtaPos = (TH1D*)TopDir->Get("histoGammaPhiEtaPos");
	if(fHistGammaPhiEtaPos){
            fHistGammaPhiEtaPos->Sumw2();
            fHistGammaPhiEtaPos->Scale(1./nGammas);
			GetMinMaxBin(fHistGammaPhiEtaPos,minB,maxB);
            SetXRange(fHistGammaPhiEtaPos,minB,maxB);
			DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kFALSE,kFALSE,
                                 fHistGammaPhiEtaPos,"","Photon #phi A side","#frac{1}{N_{#gamma}} #frac{dN}{d#phi}",1,1,
                                 0.82,0.94,0.03,fCollisionSystem,plotDataSets[i],"");
            SaveCanvasAndWriteHistogram(canvas, fHistGammaPhiEtaPos, Form("%s/Photon_Phi_EtaPos_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
			vecGammaPhiEtaPos.push_back(new TH1D(*fHistGammaPhiEtaPos));
        }else cout << "INFO: Object |histoGammaPhiEtaPos| could not be found! Skipping Draw..." << endl;
    //-----------------------------------
        TH2D* fHistXY = (TH2D*)TopDir->Get("histoGammaXY");
        if(fHistXY){
          fHistXY->Sumw2();
          fHistXY->Scale(1./nGammas);
          GetMinMaxBin(fHistXY,minB,maxB); minB-=5; maxB+=5;
          GetMinMaxBinY(fHistXY,minYB,maxYB); minYB-=5; maxYB+=5;
          SetXRange(fHistXY,minB,maxB);
          SetYRange(fHistXY,minYB,maxYB);
          SetZMinMaxTH2(fHistXY,minB,maxB,minYB,maxYB,kTRUE);
          DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,0.12,0.12,kFALSE,kFALSE,kTRUE,
                               fHistXY,"",
                               "X (cm)", "Y (cm)",0.9,1.3,
                               0.65,0.85,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
          SaveCanvasAndWriteHistogram(cvsQuadratic, fHistXY, Form("%s/Photon_XY_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
        }else cout << "INFO: Object |histoGammaXY| could not be found! Skipping Draw..." << endl;
    //-----------------------------------
        TH2D* fHistZR = (TH2D*)TopDir->Get("histoGammaZR");
        if(fHistZR){
          fHistZR->Sumw2();
          fHistZR->Scale(1./nGammas);
          GetMinMaxBin(fHistZR,minB,maxB); minB-=1; maxB+=1;
          GetMinMaxBinY(fHistZR,minYB,maxYB);
          SetXRange(fHistZR,minB,maxB);
          SetYRange(fHistZR,minYB,maxYB);
          SetZMinMaxTH2(fHistZR,minB,maxB,minYB,maxYB,kTRUE);
          DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                               fHistZR,"",
                               "Z (cm)", "R (cm)",0.9,1.3,
                               0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
          SaveCanvasAndWriteHistogram(cvsQuadratic, fHistZR, Form("%s/Photon_ZR_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
        }else cout << "INFO: Object |histoGammaZR| could not be found! Skipping Draw..." << endl;

  //-----------------------------------
  //-------------------------------------- Photon Chi2 ----------------------------------------------
  //-----------------------------------

        TH2D* fHistGammaChi2NDFPt = (TH2D*)TopDir->Get("histoGammaChi2NDFPt");
        if(fHistGammaChi2NDFPt){
          fHistGammaChi2NDFPt->Sumw2();
          fHistGammaChi2NDFPt->Scale(1./nGammas);
          GetMinMaxBin(fHistGammaChi2NDFPt,minB,maxB); maxB+=5;
          GetMinMaxBinY(fHistGammaChi2NDFPt,minYB,maxYB); maxYB+=5;
          SetXRange(fHistGammaChi2NDFPt,minB,maxB);
          SetYRange(fHistGammaChi2NDFPt,fHistGammaChi2NDFPt->GetYaxis()->FindBin(0.1),maxYB);
          SetZMinMaxTH2(fHistGammaChi2NDFPt,minB,maxB,fHistGammaChi2NDFPt->GetYaxis()->FindBin(0.1),maxYB,kTRUE);
          //-----------------------------------
          TH1D* histoGammaChi2 = (TH1D*)fHistGammaChi2NDFPt->ProjectionX("histoGammaChi2");
          ConvGammaRebinWithBinCorrection(histoGammaChi2,1);
          CheckForNegativeEntries(histoGammaChi2);
          vecChi2[0].push_back(histoGammaChi2);
          //-----------------------------------
          DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                               fHistGammaChi2NDFPt,"",
                               "#chi^{2}_{#gamma}/ndf", "#it{p}_{T,#gamma} (GeV/#it{c})",1,1.4,
                               0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
          SaveCanvasAndWriteHistogram(cvsQuadratic, fHistGammaChi2NDFPt, Form("%s/Chi2NDF_Pt_Photon_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
        }else cout << "INFO: Object |histoGammaChi2NDFPt| could not be found! Skipping Draw..." << endl;
    //-----------------------------------
        TH2D* fHistTruePrimGammaChi2NDFPt = (TH2D*)TopDir->Get("histoTruePrimGammaChi2NDFPt");
        if(fHistTruePrimGammaChi2NDFPt){
          CheckNEntries(fHistTruePrimGammaChi2NDFPt);
          fHistTruePrimGammaChi2NDFPt->Sumw2();
          fHistTruePrimGammaChi2NDFPt->Scale(1./nGammas);
          SetXRange(fHistTruePrimGammaChi2NDFPt,minB,maxB);
          SetYRange(fHistTruePrimGammaChi2NDFPt,fHistTruePrimGammaChi2NDFPt->GetYaxis()->FindBin(0.1),maxYB);
          SetZMinMaxTH2(fHistTruePrimGammaChi2NDFPt,minB,maxB,fHistTruePrimGammaChi2NDFPt->GetYaxis()->FindBin(0.1),maxYB,kTRUE);
          //-----------------------------------
          TH1D* histoMonteCarloTruePrimGammaChi2 = (TH1D*)fHistTruePrimGammaChi2NDFPt->ProjectionX("histoMonteCarloTruePrimGammaChi2");
          ConvGammaRebinWithBinCorrection(histoMonteCarloTruePrimGammaChi2,1);
          vecChi2[1].push_back(histoMonteCarloTruePrimGammaChi2);
          //-----------------------------------
          DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                               fHistTruePrimGammaChi2NDFPt,"",
                               "#chi^{2}_{#gamma}/ndf", "#it{p}_{T,#gamma} (GeV/#it{c})",1,1.4,
                               0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
          SaveCanvasAndWriteHistogram(cvsQuadratic, fHistTruePrimGammaChi2NDFPt, Form("%s/Chi2NDF_Pt_TruePrimGamma_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
        }else cout << "INFO: Object |histoTruePrimGammaChi2NDFPt| could not be found! Skipping Draw..." << endl;
    //-----------------------------------
        TH2D* fHistTrueSecGammaChi2NDFPt = (TH2D*)TopDir->Get("histoTrueSecGammaChi2NDFPt");
        if(fHistTrueSecGammaChi2NDFPt){
          CheckNEntries(fHistTrueSecGammaChi2NDFPt);
          fHistTrueSecGammaChi2NDFPt->Sumw2();
          fHistTrueSecGammaChi2NDFPt->Scale(1./nGammas);
          SetXRange(fHistTrueSecGammaChi2NDFPt,minB,maxB);
          SetYRange(fHistTrueSecGammaChi2NDFPt,fHistTrueSecGammaChi2NDFPt->GetYaxis()->FindBin(0.1),maxYB);
          SetZMinMaxTH2(fHistTrueSecGammaChi2NDFPt,minB,maxB,fHistTrueSecGammaChi2NDFPt->GetYaxis()->FindBin(0.1),maxYB,kTRUE);
          //-----------------------------------
          TH1D* histoMonteCarloTrueSecGammaChi2 = (TH1D*)fHistTrueSecGammaChi2NDFPt->ProjectionX("histoMonteCarloTrueSecGammaChi2");
          ConvGammaRebinWithBinCorrection(histoMonteCarloTrueSecGammaChi2,1);
          vecChi2[2].push_back(histoMonteCarloTrueSecGammaChi2);
          //-----------------------------------
          DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                               fHistTrueSecGammaChi2NDFPt,"",
                               "#chi^{2}_{#gamma}/ndf", "#it{p}_{T,#gamma} (GeV/#it{c})",1,1.4,
                               0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
          SaveCanvasAndWriteHistogram(cvsQuadratic, fHistTrueSecGammaChi2NDFPt, Form("%s/Chi2NDF_Pt_TrueSecGamma_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
        }else cout << "INFO: Object |histoTrueSecGammaChi2NDFPt| could not be found! Skipping Draw..." << endl;
    //-----------------------------------
        TH2D* fHistTrueDalitzGammaChi2NDFPt = (TH2D*)TopDir->Get("histoTrueDalitzGammaChi2NDFPt");
        if(fHistTrueDalitzGammaChi2NDFPt){
          CheckNEntries(fHistTrueDalitzGammaChi2NDFPt);
          fHistTrueDalitzGammaChi2NDFPt->Sumw2();
          fHistTrueDalitzGammaChi2NDFPt->Scale(1./nGammas);
          SetXRange(fHistTrueDalitzGammaChi2NDFPt,minB,maxB);
          SetYRange(fHistTrueDalitzGammaChi2NDFPt,fHistTrueDalitzGammaChi2NDFPt->GetYaxis()->FindBin(0.1),maxYB);
          SetZMinMaxTH2(fHistTrueDalitzGammaChi2NDFPt,minB,maxB,fHistTrueDalitzGammaChi2NDFPt->GetYaxis()->FindBin(0.1),maxYB,kTRUE);
          //-----------------------------------
          TH1D* histoMonteCarloTrueDalitzGammaChi2 = (TH1D*)fHistTrueDalitzGammaChi2NDFPt->ProjectionX("histoMonteCarloTrueDalitzGammaChi2");
          ConvGammaRebinWithBinCorrection(histoMonteCarloTrueDalitzGammaChi2,1);
          vecChi2[3].push_back(histoMonteCarloTrueDalitzGammaChi2);
          //-----------------------------------
          DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                               fHistTrueDalitzGammaChi2NDFPt,"",
                               "#chi^{2}_{#gamma}/ndf", "#it{p}_{T,#gamma} (GeV/#it{c})",1,1.4,
                               0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
          SaveCanvasAndWriteHistogram(cvsQuadratic, fHistTrueDalitzGammaChi2NDFPt, Form("%s/Chi2NDF_Pt_TrueDalitz_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
        }else cout << "INFO: Object |histoTrueDalitzGammaChi2NDFPt| could not be found! Skipping Draw..." << endl;
    //-----------------------------------
        if(fHistGammaChi2NDFPt && fHistTruePrimGammaChi2NDFPt && fHistTrueSecGammaChi2NDFPt && fHistTrueDalitzGammaChi2NDFPt){
          TH2D* fHistBGChi2NDFPt = (TH2D*) fHistGammaChi2NDFPt->Clone("histoMonteCarloBGChi2NDFPt");
          fHistBGChi2NDFPt->Sumw2();
          fHistBGChi2NDFPt->Add(fHistTruePrimGammaChi2NDFPt,-1.);
          fHistBGChi2NDFPt->Add(fHistTrueSecGammaChi2NDFPt,-1.);
          fHistBGChi2NDFPt->Add(fHistTrueDalitzGammaChi2NDFPt,-1.);
          SetXRange(fHistBGChi2NDFPt,minB,maxB);
          SetYRange(fHistBGChi2NDFPt,fHistBGChi2NDFPt->GetYaxis()->FindBin(0.1),maxYB);
          SetZMinMaxTH2(fHistBGChi2NDFPt,minB,maxB,fHistBGChi2NDFPt->GetYaxis()->FindBin(0.1),maxYB,kTRUE);
          CheckForNegativeEntries(fHistBGChi2NDFPt);
          //-----------------------------------
          TH1D* histoMonteCarloBGChi2 = (TH1D*)fHistBGChi2NDFPt->ProjectionX("histoMonteCarloBGChi2");
          ConvGammaRebinWithBinCorrection(histoMonteCarloBGChi2,1);
          CheckForNegativeEntries(histoMonteCarloBGChi2);
          vecChi2[4].push_back(histoMonteCarloBGChi2);
          //-----------------------------------
          DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                               fHistBGChi2NDFPt,"",
                               "#chi^{2}_{#gamma}/ndf", "#it{p}_{T,#gamma} (GeV/#it{c})",1,1.4,
                               0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
          SaveCanvasAndWriteHistogram(cvsQuadratic, fHistBGChi2NDFPt, Form("%s/Chi2NDF_Pt_BG_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));

          delete fHistBGChi2NDFPt;
        }

    //-----------------------------------
    //-------------------------------------- Photon PsiPair ----------------------------------------------
    //-----------------------------------

        TH2D* fHistGammaPsiPairPt = (TH2D*)TopDir->Get("histoGammaPsiPairPt");
        if(fHistGammaPsiPairPt){
          fHistGammaPsiPairPt->Sumw2();
          fHistGammaPsiPairPt->Scale(1./nGammas);
          GetMinMaxBin(fHistGammaPsiPairPt,minB,maxB); minB-=5; maxB+=5;
          GetMinMaxBinY(fHistGammaPsiPairPt,minYB,maxYB); maxYB+=5;
          SetXRange(fHistGammaPsiPairPt,minB,maxB);
          SetYRange(fHistGammaPsiPairPt,fHistGammaPsiPairPt->GetYaxis()->FindBin(0.1),maxYB);
          SetZMinMaxTH2(fHistGammaPsiPairPt,minB,maxB,fHistGammaPsiPairPt->GetYaxis()->FindBin(0.1),maxYB,kTRUE);
          //-----------------------------------
          TH1D* histoGammaPsiPair = (TH1D*)fHistGammaPsiPairPt->ProjectionX("histoGammaPsiPair");
          ConvGammaRebinWithBinCorrection(histoGammaPsiPair,1);
          vecPsi[0].push_back(histoGammaPsiPair);
          //-----------------------------------
          DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                               fHistGammaPsiPairPt,"",
                               "#psi_{pair}", "#it{p}_{T,#gamma} (GeV/#it{c})",1,1.4,
                               0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
          SaveCanvasAndWriteHistogram(cvsQuadratic, fHistGammaPsiPairPt, Form("%s/PsiPair_Pt_Photon_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
        }else cout << "INFO: Object |histoGammaPsiPairPt| could not be found! Skipping Draw..." << endl;
    //-----------------------------------
        TH2D* fHistTruePrimGammaPsiPairPt = (TH2D*)TopDir->Get("histoTruePrimGammaPsiPairPt");
        if(fHistTruePrimGammaPsiPairPt){
          CheckNEntries(fHistTruePrimGammaPsiPairPt);
          fHistTruePrimGammaPsiPairPt->Sumw2();
          fHistTruePrimGammaPsiPairPt->Scale(1./nGammas);
          SetXRange(fHistTruePrimGammaPsiPairPt,minB,maxB);
          SetYRange(fHistTruePrimGammaPsiPairPt,fHistTruePrimGammaPsiPairPt->GetYaxis()->FindBin(0.1),maxYB);
          SetZMinMaxTH2(fHistTruePrimGammaPsiPairPt,minB,maxB,fHistTruePrimGammaPsiPairPt->GetYaxis()->FindBin(0.1),maxYB,kTRUE);
          //-----------------------------------
          TH1D* histoMonteCarloTruePrimGammaPsiPair = (TH1D*)fHistTruePrimGammaPsiPairPt->ProjectionX("histoMonteCarloTruePrimGammaPsiPair");
          ConvGammaRebinWithBinCorrection(histoMonteCarloTruePrimGammaPsiPair,1);
          vecPsi[1].push_back(histoMonteCarloTruePrimGammaPsiPair);
          //-----------------------------------
          DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                               fHistTruePrimGammaPsiPairPt,"",
                               "#psi_{pair}", "#it{p}_{T,#gamma} (GeV/#it{c})",1,1.4,
                               0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
          SaveCanvasAndWriteHistogram(cvsQuadratic, fHistTruePrimGammaPsiPairPt, Form("%s/PsiPair_Pt_TruePrimGamma_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
        }else cout << "INFO: Object |histoTruePrimGammaPsiPairPt| could not be found! Skipping Draw..." << endl;
    //-----------------------------------
        TH2D* fHistTrueSecGammaPsiPairPt = (TH2D*)TopDir->Get("histoTrueSecGammaPsiPairPt");
        if(fHistTrueSecGammaPsiPairPt){
          CheckNEntries(fHistTrueSecGammaPsiPairPt);
          fHistTrueSecGammaPsiPairPt->Sumw2();
          fHistTrueSecGammaPsiPairPt->Scale(1./nGammas);
          SetXRange(fHistTrueSecGammaPsiPairPt,minB,maxB);
          SetYRange(fHistTrueSecGammaPsiPairPt,fHistTrueSecGammaPsiPairPt->GetYaxis()->FindBin(0.1),maxYB);
          SetZMinMaxTH2(fHistTrueSecGammaPsiPairPt,minB,maxB,fHistTrueSecGammaPsiPairPt->GetYaxis()->FindBin(0.1),maxYB,kTRUE);
          //-----------------------------------
          TH1D* histoMonteCarloTrueSecGammaPsiPair = (TH1D*)fHistTrueSecGammaPsiPairPt->ProjectionX("histoMonteCarloTrueSecGammaPsiPair");
          ConvGammaRebinWithBinCorrection(histoMonteCarloTrueSecGammaPsiPair,1);
          vecPsi[2].push_back(histoMonteCarloTrueSecGammaPsiPair);
          //-----------------------------------
          DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                               fHistTrueSecGammaPsiPairPt,"",
                               "#psi_{pair}", "#it{p}_{T,#gamma} (GeV/#it{c})",1,1.4,
                               0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
          SaveCanvasAndWriteHistogram(cvsQuadratic, fHistTrueSecGammaPsiPairPt, Form("%s/PsiPair_Pt_TrueSecGamma_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
        }else cout << "INFO: Object |histoTrueSecGammaPsiPairPt| could not be found! Skipping Draw..." << endl;
    //-----------------------------------
        TH2D* fHistTrueDalitzGammaPsiPairPt = (TH2D*)TopDir->Get("histoTrueDalitzGammaPsiPairPt");
        if(fHistTrueDalitzGammaPsiPairPt){
          CheckNEntries(fHistTrueDalitzGammaPsiPairPt);
          fHistTrueDalitzGammaPsiPairPt->Sumw2();
          fHistTrueDalitzGammaPsiPairPt->Scale(1./nGammas);
          SetXRange(fHistTrueDalitzGammaPsiPairPt,minB,maxB);
          SetYRange(fHistTrueDalitzGammaPsiPairPt,fHistTrueDalitzGammaPsiPairPt->GetYaxis()->FindBin(0.1),maxYB);
          SetZMinMaxTH2(fHistTrueDalitzGammaPsiPairPt,minB,maxB,fHistTrueDalitzGammaPsiPairPt->GetYaxis()->FindBin(0.1),maxYB,kTRUE);
          //-----------------------------------
          TH1D* histoMonteCarloTrueDalitzGammaPsiPair = (TH1D*)fHistTrueDalitzGammaPsiPairPt->ProjectionX("histoMonteCarloTrueDalitzGammaPsiPair");
          ConvGammaRebinWithBinCorrection(histoMonteCarloTrueDalitzGammaPsiPair,1);
          vecPsi[3].push_back(histoMonteCarloTrueDalitzGammaPsiPair);
          //-----------------------------------
          DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                               fHistTrueDalitzGammaPsiPairPt,"",
                               "#psi_{pair}", "#it{p}_{T,#gamma} (GeV/#it{c})",1,1.4,
                               0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
          SaveCanvasAndWriteHistogram(cvsQuadratic, fHistTrueDalitzGammaPsiPairPt, Form("%s/PsiPair_Pt_TrueDalitz_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
        }else cout << "INFO: Object |histoTrueDalitzGammaPsiPairPt| could not be found! Skipping Draw..." << endl;
    //-----------------------------------
        if(fHistGammaPsiPairPt && fHistTruePrimGammaPsiPairPt && fHistTrueSecGammaPsiPairPt && fHistTrueDalitzGammaPsiPairPt){
          TH2D* fHistBGPsiPairPt = (TH2D*) fHistGammaPsiPairPt->Clone("histoBGPsiPairMonteCarloPt");
          fHistBGPsiPairPt->Sumw2();
          fHistBGPsiPairPt->Add(fHistTruePrimGammaPsiPairPt,-1.);
          fHistBGPsiPairPt->Add(fHistTrueSecGammaPsiPairPt,-1.);
          fHistBGPsiPairPt->Add(fHistTrueDalitzGammaPsiPairPt,-1.);
          SetXRange(fHistBGPsiPairPt,minB,maxB);
          SetYRange(fHistBGPsiPairPt,fHistBGPsiPairPt->GetYaxis()->FindBin(0.1),maxYB);
          SetZMinMaxTH2(fHistBGPsiPairPt,minB,maxB,fHistBGPsiPairPt->GetYaxis()->FindBin(0.1),maxYB,kTRUE);
          CheckForNegativeEntries(fHistBGPsiPairPt);
          //-----------------------------------
          TH1D* histoMonteCarloBGPsiPair = (TH1D*)fHistBGPsiPairPt->ProjectionX("histoMonteCarloBGPsiPair");
          ConvGammaRebinWithBinCorrection(histoMonteCarloBGPsiPair,1);
          CheckForNegativeEntries(histoMonteCarloBGPsiPair);
          vecPsi[4].push_back(histoMonteCarloBGPsiPair);
          //-----------------------------------
          DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                               fHistBGPsiPairPt,"",
                               "#psi_{pair}", "#it{p}_{T,#gamma} (GeV/#it{c})",1,1.4,
                               0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
          SaveCanvasAndWriteHistogram(cvsQuadratic, fHistBGPsiPairPt, Form("%s/PsiPair_Pt_BG_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));

          delete fHistBGPsiPairPt;
        }
  //-----------------------------------
  //-------------------------------------- Photon CosPoint ----------------------------------------------
  //-----------------------------------

        TH2D* fHistGammaCosPointPt = (TH2D*)TopDir->Get("histoGammaCosPointPt");
        if(fHistGammaCosPointPt){
          fHistGammaCosPointPt->Sumw2();
          fHistGammaCosPointPt->Scale(1./nGammas);
          GetMinMaxBin(fHistGammaCosPointPt,minB,maxB); minB-=5;
          GetMinMaxBinY(fHistGammaCosPointPt,minYB,maxYB); maxYB+=5;
          SetXRange(fHistGammaCosPointPt,minB,maxB);
          SetYRange(fHistGammaCosPointPt,fHistGammaCosPointPt->GetYaxis()->FindBin(0.1),maxYB);
          SetZMinMaxTH2(fHistGammaCosPointPt,minB,maxB,fHistGammaCosPointPt->GetYaxis()->FindBin(0.1),maxYB,kTRUE);
          //-----------------------------------
          TH1D* histoGammaCosPoint = (TH1D*)fHistGammaCosPointPt->ProjectionX("histoGammaCosPoint");
          ConvGammaRebinWithBinCorrection(histoGammaCosPoint,1);
          vecCos[0].push_back(histoGammaCosPoint);
          //-----------------------------------
          TH1D* histoGammaPt = (TH1D*)fHistGammaCosPointPt->ProjectionY("histoGammaPt");
          ConvGammaRebinWithBinCorrection(histoGammaPt,1);
          vecGamma[0].push_back(histoGammaPt);
          //-----------------------------------
          DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kTRUE,kTRUE,
                               fHistGammaCosPointPt,"",
                               "cos(#theta_{point})", "#it{p}_{T,#gamma} (GeV/#it{c})",1,1.4,
                               0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
          SaveCanvasAndWriteHistogram(cvsQuadratic, fHistGammaCosPointPt, Form("%s/CosPoint_Pt_Photon_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
        }else cout << "INFO: Object |histoGammaCosPointPt| could not be found! Skipping Draw..." << endl;
    //-----------------------------------
        TH2D* fHistTruePrimGammaCosPointPt = (TH2D*)TopDir->Get("histoTruePrimGammaCosPointPt");
        if(fHistTruePrimGammaCosPointPt){
          CheckNEntries(fHistTruePrimGammaCosPointPt);
          fHistTruePrimGammaCosPointPt->Sumw2();
          fHistTruePrimGammaCosPointPt->Scale(1./nGammas);
          SetXRange(fHistTruePrimGammaCosPointPt,minB,maxB);
          SetYRange(fHistTruePrimGammaCosPointPt,fHistTruePrimGammaCosPointPt->GetYaxis()->FindBin(0.1),maxYB);
          SetZMinMaxTH2(fHistTruePrimGammaCosPointPt,minB,maxB,fHistTruePrimGammaCosPointPt->GetYaxis()->FindBin(0.1),maxYB,kTRUE);
          //-----------------------------------
          TH1D* histoMonteCarloTruePrimGammaCosPoint = (TH1D*)fHistTruePrimGammaCosPointPt->ProjectionX("histoMonteCarloTruePrimGammaCosPoint");
          ConvGammaRebinWithBinCorrection(histoMonteCarloTruePrimGammaCosPoint,1);
          vecCos[1].push_back(histoMonteCarloTruePrimGammaCosPoint);
          //-----------------------------------
          TH1D* histoTruePrimGammaPt = (TH1D*)fHistTruePrimGammaCosPointPt->ProjectionY("histoTruePrimGammaPt");
          ConvGammaRebinWithBinCorrection(histoTruePrimGammaPt,1);
          vecGamma[1].push_back(histoTruePrimGammaPt);
          //-----------------------------------
          DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kTRUE,kTRUE,
                               fHistTruePrimGammaCosPointPt,"",
                               "cos(#theta_{point})", "#it{p}_{T,#gamma} (GeV/#it{c})",1,1.4,
                               0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
          SaveCanvasAndWriteHistogram(cvsQuadratic, fHistTruePrimGammaCosPointPt, Form("%s/CosPoint_Pt_TruePrimGamma_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
        }else cout << "INFO: Object |histoTruePrimGammaCosPointPt| could not be found! Skipping Draw..." << endl;
    //-----------------------------------
        TH2D* fHistTrueSecGammaCosPointPt = (TH2D*)TopDir->Get("histoTrueSecGammaCosPointPt");
        if(fHistTrueSecGammaCosPointPt){
          CheckNEntries(fHistTrueSecGammaCosPointPt);
          fHistTrueSecGammaCosPointPt->Sumw2();
          fHistTrueSecGammaCosPointPt->Scale(1./nGammas);
          SetXRange(fHistTrueSecGammaCosPointPt,minB,maxB);
          SetYRange(fHistTrueSecGammaCosPointPt,fHistTrueSecGammaCosPointPt->GetYaxis()->FindBin(0.1),maxYB);
          SetZMinMaxTH2(fHistTrueSecGammaCosPointPt,minB,maxB,fHistTrueSecGammaCosPointPt->GetYaxis()->FindBin(0.1),maxYB,kTRUE);
          //-----------------------------------
          TH1D* histoMonteCarloTrueSecGammaCosPoint = (TH1D*)fHistTrueSecGammaCosPointPt->ProjectionX("histoMonteCarloTrueSecGammaCosPoint");
          ConvGammaRebinWithBinCorrection(histoMonteCarloTrueSecGammaCosPoint,1);
          vecCos[2].push_back(histoMonteCarloTrueSecGammaCosPoint);
          //-----------------------------------
          TH1D* histoTrueSecGammaPt = (TH1D*)fHistTrueSecGammaCosPointPt->ProjectionY("histoTrueSecGammaPt");
          ConvGammaRebinWithBinCorrection(histoTrueSecGammaPt,1);
          vecGamma[2].push_back(histoTrueSecGammaPt);
          //-----------------------------------
          DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kTRUE,kTRUE,
                               fHistTrueSecGammaCosPointPt,"",
                               "cos(#theta_{point})", "#it{p}_{T,#gamma} (GeV/#it{c})",1,1.4,
                               0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
          SaveCanvasAndWriteHistogram(cvsQuadratic, fHistTrueSecGammaCosPointPt, Form("%s/CosPoint_Pt_TrueSecGamma_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
        }else cout << "INFO: Object |histoTrueSecGammaCosPointPt| could not be found! Skipping Draw..." << endl;
    //-----------------------------------
        TH2D* fHistTrueDalitzGammaCosPointPt = (TH2D*)TopDir->Get("histoTrueDalitzGammaCosPointPt");
        if(fHistTrueDalitzGammaCosPointPt){
          CheckNEntries(fHistTrueDalitzGammaCosPointPt);
          fHistTrueDalitzGammaCosPointPt->Sumw2();
          fHistTrueDalitzGammaCosPointPt->Scale(1./nGammas);
          SetXRange(fHistTrueDalitzGammaCosPointPt,minB,maxB);
          SetYRange(fHistTrueDalitzGammaCosPointPt,fHistTrueDalitzGammaCosPointPt->GetYaxis()->FindBin(0.1),maxYB);
          SetZMinMaxTH2(fHistTrueDalitzGammaCosPointPt,minB,maxB,fHistTrueDalitzGammaCosPointPt->GetYaxis()->FindBin(0.1),maxYB,kTRUE);
          //-----------------------------------
          TH1D* histoMonteCarloTrueDalitzGammaCosPoint = (TH1D*)fHistTrueDalitzGammaCosPointPt->ProjectionX("histoMonteCarloTrueDalitzGammaCosPoint");
          ConvGammaRebinWithBinCorrection(histoMonteCarloTrueDalitzGammaCosPoint,1);
          vecCos[3].push_back(histoMonteCarloTrueDalitzGammaCosPoint);
          //-----------------------------------
          TH1D* histoTrueDalitzGammaPt = (TH1D*)fHistTrueDalitzGammaCosPointPt->ProjectionY("histoTrueDalitzGammaPt");
          ConvGammaRebinWithBinCorrection(histoTrueDalitzGammaPt,1);
          vecGamma[3].push_back(histoTrueDalitzGammaPt);
          //-----------------------------------
          DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kTRUE,kTRUE,
                               fHistTrueDalitzGammaCosPointPt,"",
                               "cos(#theta_{point})", "#it{p}_{T,#gamma} (GeV/#it{c})",1,1.4,
                               0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
          SaveCanvasAndWriteHistogram(cvsQuadratic, fHistTrueDalitzGammaCosPointPt, Form("%s/CosPoint_Pt_TrueDalitz_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
        }else cout << "INFO: Object |histoTrueDalitzGammaCosPointPt| could not be found! Skipping Draw..." << endl;
    //-----------------------------------
        if(fHistGammaCosPointPt && fHistTruePrimGammaCosPointPt && fHistTrueSecGammaCosPointPt && fHistTrueDalitzGammaCosPointPt){
          TH2D* fHistBGCosPointPt = (TH2D*) fHistGammaCosPointPt->Clone("histoBGCosPointMonteCarloPt");
          fHistBGCosPointPt->Sumw2();
          fHistBGCosPointPt->Add(fHistTruePrimGammaCosPointPt,-1.);
          fHistBGCosPointPt->Add(fHistTrueSecGammaCosPointPt,-1.);
          fHistBGCosPointPt->Add(fHistTrueDalitzGammaCosPointPt,-1.);
          SetXRange(fHistBGCosPointPt,minB,maxB);
          SetYRange(fHistBGCosPointPt,fHistBGCosPointPt->GetYaxis()->FindBin(0.1),maxYB);
          SetZMinMaxTH2(fHistBGCosPointPt,minB,maxB,minYB,maxYB,kTRUE);
          CheckForNegativeEntries(fHistBGCosPointPt);
          //-----------------------------------
          TH1D* histoMonteCarloBGCosPoint = (TH1D*)fHistBGCosPointPt->ProjectionX("histoMonteCarloBGCosPoint");
          ConvGammaRebinWithBinCorrection(histoMonteCarloBGCosPoint,1);
          vecCos[4].push_back(histoMonteCarloBGCosPoint);
          //-----------------------------------
          TH1D* histoBGGammaPt = (TH1D*)fHistBGCosPointPt->ProjectionY("histoBGGammaPt");
          ConvGammaRebinWithBinCorrection(histoBGGammaPt,1);
          CheckForNegativeEntries(histoBGGammaPt);
          vecGamma[4].push_back(histoBGGammaPt);
          //-----------------------------------
          DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kTRUE,kTRUE,
                               fHistBGCosPointPt,"",
                               "cos(#theta_{point})", "#it{p}_{T,#gamma} (GeV/#it{c})",1,1.4,
                               0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
          SaveCanvasAndWriteHistogram(cvsQuadratic, fHistBGCosPointPt, Form("%s/CosPoint_Pt_BG_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));

          delete fHistBGCosPointPt;
        }

//-----------------------------------
//-------------------------------------- Photon CosPoint vs Chi2 ----------------------------------------------
//-----------------------------------

        TH2D* fHistGammaCosPointChi2 = (TH2D*)TopDir->Get("histoGammaCosPointChi2");
        if(fHistGammaCosPointChi2){
          fHistGammaCosPointChi2->Sumw2();
          fHistGammaCosPointChi2->Scale(1./nGammas);
          GetMinMaxBin(fHistGammaCosPointChi2,minB,maxB); minB-=5;
          GetMinMaxBinY(fHistGammaCosPointChi2,minYB,maxYB); maxYB+=5;
          SetXRange(fHistGammaCosPointChi2,minB,maxB);
          SetYRange(fHistGammaCosPointChi2,minYB,maxYB);
          SetZMinMaxTH2(fHistGammaCosPointChi2,minB,maxB,minYB,maxYB,kTRUE);
          //-----------------------------------
          DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                               fHistGammaCosPointChi2,"",
                               "cos(#theta_{point})", "#chi^{2}_{#gamma}/ndf",1,1.4,
                               0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
          SaveCanvasAndWriteHistogram(cvsQuadratic, fHistGammaCosPointChi2, Form("%s/CosPoint_Chi2_Photon_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
        }else cout << "INFO: Object |histoGammaCosPointChi2| could not be found! Skipping Draw..." << endl;
    //-----------------------------------
        TH2D* fHistTrueGammaCosPointChi2 = (TH2D*)TopDir->Get("histoTrueGammaCosPointChi2");
        if(fHistTrueGammaCosPointChi2){
          CheckNEntries(fHistTrueGammaCosPointChi2);
          fHistTrueGammaCosPointChi2->Sumw2();
          fHistTrueGammaCosPointChi2->Scale(1./nGammas);
          SetXRange(fHistTrueGammaCosPointChi2,minB,maxB);
          SetYRange(fHistTrueGammaCosPointChi2,minYB,maxYB);
          SetZMinMaxTH2(fHistTrueGammaCosPointChi2,minB,maxB,minYB,maxYB,kTRUE);
          //-----------------------------------
          DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                               fHistTrueGammaCosPointChi2,"",
                               "cos(#theta_{point})", "#chi^{2}_{#gamma}/ndf",1,1.4,
                               0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
          SaveCanvasAndWriteHistogram(cvsQuadratic, fHistTrueGammaCosPointChi2, Form("%s/CosPoint_Chi2_TrueGamma_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
        }else cout << "INFO: Object |histoTrueGammaCosPointChi2| could not be found! Skipping Draw..." << endl;
    //-----------------------------------
        TH2D* fHistTrueDalitzGammaCosPointChi2 = (TH2D*)TopDir->Get("histoTrueDalitzGammaCosPointChi2");
        if(fHistTrueDalitzGammaCosPointChi2){
          CheckNEntries(fHistTrueDalitzGammaCosPointChi2);
          fHistTrueDalitzGammaCosPointChi2->Sumw2();
          fHistTrueDalitzGammaCosPointChi2->Scale(1./nGammas);
          SetXRange(fHistTrueDalitzGammaCosPointChi2,minB,maxB);
          SetYRange(fHistTrueDalitzGammaCosPointChi2,minYB,maxYB);
          SetZMinMaxTH2(fHistTrueDalitzGammaCosPointChi2,minB,maxB,minYB,maxYB,kTRUE);
          //-----------------------------------
          DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                               fHistTrueDalitzGammaCosPointChi2,"",
                               "cos(#theta_{point})", "#chi^{2}_{#gamma}/ndf",1,1.4,
                               0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
          SaveCanvasAndWriteHistogram(cvsQuadratic, fHistTrueDalitzGammaCosPointChi2, Form("%s/CosPoint_Chi2_TrueDalitz_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
        }else cout << "INFO: Object |histoTrueDalitzGammaCosPointChi2| could not be found! Skipping Draw..." << endl;
    //-----------------------------------
        if(fHistGammaCosPointChi2 && fHistTrueGammaCosPointChi2 && fHistTrueDalitzGammaCosPointChi2){
          TH2D* fHistBGCosPointChi2 = (TH2D*) fHistGammaCosPointChi2->Clone("histoBGCosPointChi2MonteCarlo");
          fHistBGCosPointChi2->Sumw2();
          fHistBGCosPointChi2->Add(fHistTrueGammaCosPointChi2,-1.);
          fHistBGCosPointChi2->Add(fHistTrueDalitzGammaCosPointChi2,-1.);
          SetXRange(fHistBGCosPointChi2,minB,maxB);
          SetYRange(fHistBGCosPointChi2,minYB,maxYB);
          SetZMinMaxTH2(fHistBGCosPointChi2,minB,maxB,minYB,maxYB,kTRUE);
          CheckForNegativeEntries(fHistBGCosPointChi2);
          //-----------------------------------
          DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                               fHistBGCosPointChi2,"",
                               "cos(#theta_{point})", "#chi^{2}_{#gamma}/ndf",1,1.4,
                               0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
          SaveCanvasAndWriteHistogram(cvsQuadratic, fHistBGCosPointChi2, Form("%s/CosPoint_Chi2_BG_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));

          delete fHistBGCosPointChi2;
        }

//-----------------------------------
//-------------------------------------- Photon Chi2 vs PsiPair ----------------------------------------------
//-----------------------------------

        TH2D* fHistGammaChi2PsiPair = (TH2D*)TopDir->Get("histoGammaChi2PsiPair");
        if(fHistGammaChi2PsiPair){
          fHistGammaChi2PsiPair->Sumw2();
          fHistGammaChi2PsiPair->Scale(1./nGammas);
          GetMinMaxBin(fHistGammaChi2PsiPair,minB,maxB); maxB+=5;
          GetMinMaxBinY(fHistGammaChi2PsiPair,minYB,maxYB); minYB-=5; maxYB+=5;
          SetXRange(fHistGammaChi2PsiPair,minB,maxB);
          SetYRange(fHistGammaChi2PsiPair,minYB,maxYB);
          SetZMinMaxTH2(fHistGammaChi2PsiPair,minB,maxB,minYB,maxYB,kTRUE);
          //-----------------------------------
          DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                               fHistGammaChi2PsiPair,"",
                               "#chi^{2}_{#gamma}/ndf", "#psi_{pair}",1,1.4,
                               0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
          SaveCanvasAndWriteHistogram(cvsQuadratic, fHistGammaChi2PsiPair, Form("%s/Chi2_PsiPair_Photon_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
        }else cout << "INFO: Object |histoGammaChi2PsiPair| could not be found! Skipping Draw..." << endl;
    //-----------------------------------
        TH2D* fHistTrueGammaChi2PsiPair = (TH2D*)TopDir->Get("histoTrueGammaChi2PsiPair");
        if(fHistTrueGammaChi2PsiPair){
          fHistTrueGammaChi2PsiPair->Sumw2();
          fHistTrueGammaChi2PsiPair->Scale(1./nGammas);
          SetXRange(fHistTrueGammaChi2PsiPair,minB,maxB);
          SetYRange(fHistTrueGammaChi2PsiPair,minYB,maxYB);
          SetZMinMaxTH2(fHistTrueGammaChi2PsiPair,minB,maxB,minYB,maxYB,kTRUE);
          //-----------------------------------
          DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                               fHistTrueGammaChi2PsiPair,"",
                               "#chi^{2}_{#gamma}/ndf", "#psi_{pair}",1,1.4,
                               0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
          SaveCanvasAndWriteHistogram(cvsQuadratic, fHistTrueGammaChi2PsiPair, Form("%s/Chi2_PsiPair_TrueGamma_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
        }else cout << "INFO: Object |histoTrueGammaChi2PsiPair| could not be found! Skipping Draw..." << endl;
    //-----------------------------------
        TH2D* fHistTrueDalitzGammaChi2PsiPair = (TH2D*)TopDir->Get("histoTrueDalitzGammaChi2PsiPair");
        if(fHistTrueDalitzGammaChi2PsiPair){
          fHistTrueDalitzGammaChi2PsiPair->Sumw2();
          fHistTrueDalitzGammaChi2PsiPair->Scale(1./nGammas);
          SetXRange(fHistTrueDalitzGammaChi2PsiPair,minB,maxB);
          SetYRange(fHistTrueDalitzGammaChi2PsiPair,minYB,maxYB);
          SetZMinMaxTH2(fHistTrueDalitzGammaChi2PsiPair,minB,maxB,minYB,maxYB,kTRUE);
          //-----------------------------------
          DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                               fHistTrueDalitzGammaChi2PsiPair,"",
                               "#chi^{2}_{#gamma}/ndf", "#psi_{pair}",1,1.4,
                               0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
          SaveCanvasAndWriteHistogram(cvsQuadratic, fHistTrueDalitzGammaChi2PsiPair, Form("%s/Chi2_PsiPair_TrueDalitz_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
        }else cout << "INFO: Object |histoTrueDalitzGammaChi2PsiPair| could not be found! Skipping Draw..." << endl;
    //-----------------------------------
        if(fHistGammaChi2PsiPair && fHistTrueGammaChi2PsiPair && fHistTrueDalitzGammaChi2PsiPair){
          TH2D* fHistBGChi2PsiPair = (TH2D*) fHistGammaChi2PsiPair->Clone("histoBGChi2PsiPairMonteCarlo");
          fHistBGChi2PsiPair->Sumw2();
          fHistBGChi2PsiPair->Add(fHistTrueGammaChi2PsiPair,-1.);
          fHistBGChi2PsiPair->Add(fHistTrueDalitzGammaChi2PsiPair,-1.);
          SetXRange(fHistBGChi2PsiPair,minB,maxB);
          SetYRange(fHistBGChi2PsiPair,minYB,maxYB);
          SetZMinMaxTH2(fHistBGChi2PsiPair,minB,maxB,minYB,maxYB,kTRUE);
          CheckForNegativeEntries(fHistBGChi2PsiPair);
          //-----------------------------------
          DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                               fHistBGChi2PsiPair,"",
                               "#chi^{2}_{#gamma}/ndf", "#psi_{pair}",1,1.4,
                               0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
          SaveCanvasAndWriteHistogram(cvsQuadratic, fHistBGChi2PsiPair, Form("%s/Chi2_PsiPair_BG_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));

          delete fHistBGChi2PsiPair;
        }

//-----------------------------------
//-------------------------------------- Photon CosPoint vs PsiPair ----------------------------------------------
//-----------------------------------

        TH2D* fHistGammaCosPointPsiPair = (TH2D*)TopDir->Get("histoGammaCosPointPsiPair");
        if(fHistGammaCosPointPsiPair){
          fHistGammaCosPointPsiPair->Sumw2();
          fHistGammaCosPointPsiPair->Scale(1./nGammas);
          GetMinMaxBin(fHistGammaCosPointPsiPair,minB,maxB); minB-=5;
          GetMinMaxBinY(fHistGammaCosPointPsiPair,minYB,maxYB); minYB-=5; maxYB+=5;
          SetXRange(fHistGammaCosPointPsiPair,minB,maxB);
          SetYRange(fHistGammaCosPointPsiPair,minYB,maxYB);
          SetZMinMaxTH2(fHistGammaCosPointPsiPair,minB,maxB,minYB,maxYB,kTRUE);
          //-----------------------------------
          DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                               fHistGammaCosPointPsiPair,"",
                               "cos(#theta_{point})", "#psi_{pair}",1,1.4,
                               0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
          SaveCanvasAndWriteHistogram(cvsQuadratic, fHistGammaCosPointPsiPair, Form("%s/CosPoint_PsiPair_Photon_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
        }else cout << "INFO: Object |histoGammaCosPointPsiPair| could not be found! Skipping Draw..." << endl;
    //-----------------------------------
        TH2D* fHistTrueGammaCosPointPsiPair = (TH2D*)TopDir->Get("histoTrueGammaCosPointPsiPair");
        if(fHistTrueGammaCosPointPsiPair){
          CheckNEntries(fHistTrueGammaCosPointPsiPair);
          fHistTrueGammaCosPointPsiPair->Sumw2();
          fHistTrueGammaCosPointPsiPair->Scale(1./nGammas);
          SetXRange(fHistTrueGammaCosPointPsiPair,minB,maxB);
          SetYRange(fHistTrueGammaCosPointPsiPair,minYB,maxYB);
          SetZMinMaxTH2(fHistTrueGammaCosPointPsiPair,minB,maxB,minYB,maxYB,kTRUE);
          //-----------------------------------
          DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                               fHistTrueGammaCosPointPsiPair,"",
                               "cos(#theta_{point})", "#psi_{pair}",1,1.4,
                               0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
          SaveCanvasAndWriteHistogram(cvsQuadratic, fHistTrueGammaCosPointPsiPair, Form("%s/CosPoint_PsiPair_TrueGamma_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
        }else cout << "INFO: Object |histoTrueGammaCosPointPsiPair| could not be found! Skipping Draw..." << endl;
    //-----------------------------------
        TH2D* fHistTrueDalitzGammaCosPointPsiPair = (TH2D*)TopDir->Get("histoTrueDalitzGammaCosPointPsiPair");
        if(fHistTrueDalitzGammaCosPointPsiPair){
          CheckNEntries(fHistTrueDalitzGammaCosPointPsiPair);
          fHistTrueDalitzGammaCosPointPsiPair->Sumw2();
          fHistTrueDalitzGammaCosPointPsiPair->Scale(1./nGammas);
          SetXRange(fHistTrueDalitzGammaCosPointPsiPair,minB,maxB);
          SetYRange(fHistTrueDalitzGammaCosPointPsiPair,minYB,maxYB);
          SetZMinMaxTH2(fHistTrueDalitzGammaCosPointPsiPair,minB,maxB,minYB,maxYB,kTRUE);
          //-----------------------------------
          DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                               fHistTrueDalitzGammaCosPointPsiPair,"",
                               "cos(#theta_{point})", "#psi_{pair}",1,1.4,
                               0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
          SaveCanvasAndWriteHistogram(cvsQuadratic, fHistTrueDalitzGammaCosPointPsiPair, Form("%s/CosPoint_PsiPair_TrueDalitz_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
        }else cout << "INFO: Object |histoTrueDalitzGammaCosPointPsiPair| could not be found! Skipping Draw..." << endl;
    //-----------------------------------
        if(fHistGammaCosPointPsiPair && fHistTrueGammaCosPointPsiPair && fHistTrueDalitzGammaCosPointPsiPair){
          TH2D* fHistBGCosPointPsiPair = (TH2D*) fHistGammaCosPointPsiPair->Clone("histoBGCosPointPsiPairMonteCarlo");
          fHistBGCosPointPsiPair->Sumw2();
          fHistBGCosPointPsiPair->Add(fHistTrueGammaCosPointPsiPair,-1.);
          fHistBGCosPointPsiPair->Add(fHistTrueDalitzGammaCosPointPsiPair,-1.);
          SetXRange(fHistBGCosPointPsiPair,minB,maxB);
          SetYRange(fHistBGCosPointPsiPair,minYB,maxYB);
          SetZMinMaxTH2(fHistBGCosPointPsiPair,minB,maxB,minYB,maxYB,kTRUE);
          CheckForNegativeEntries(fHistBGCosPointPsiPair);
          //-----------------------------------
          DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                               fHistBGCosPointPsiPair,"",
                               "cos(#theta_{point})", "#psi_{pair}",1,1.4,
                               0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
          SaveCanvasAndWriteHistogram(cvsQuadratic, fHistBGCosPointPsiPair, Form("%s/CosPoint_PsiPair_BG_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));

          delete fHistBGCosPointPsiPair;
        }

//-----------------------------------
//-------------------------------------- Photon Asymmetrie vs. P ----------------------------------------------
//-----------------------------------

        TH2D* fHistGammaAsymP = (TH2D*)TopDir->Get("histoGammaAsymP");
        if(fHistGammaAsymP){
          fHistGammaAsymP->Sumw2();
          fHistGammaAsymP->Scale(1./nGammas);
          GetMinMaxBin(fHistGammaAsymP,minB,maxB);
          GetMinMaxBinY(fHistGammaAsymP,minYB,maxYB);
          SetXRange(fHistGammaAsymP,minB,maxB);
          SetYRange(fHistGammaAsymP,fHistGammaAsymP->GetYaxis()->FindBin(0.1),maxYB);
          SetZMinMaxTH2(fHistGammaAsymP,minB,maxB,fHistGammaAsymP->GetYaxis()->FindBin(0.1),maxYB,kTRUE);
          //-----------------------------------
          DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kTRUE,kTRUE,
                               fHistGammaAsymP,"",
                               "#alpha", "#it{p}_{#gamma} (GeV/#it{c})",1,1.4,
                               0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
          SaveCanvasAndWriteHistogram(cvsQuadratic, fHistGammaAsymP, Form("%s/Asymmetry_P_Photon_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
        }else cout << "INFO: Object |histoGammaAsymP| could not be found! Skipping Draw..." << endl;
        //-----------------------------------
        TH2D* fHistTruePrimGammaAsymP = (TH2D*)TopDir->Get("histoTruePrimGammaAsymP");
        if(fHistTruePrimGammaAsymP){
          CheckNEntries(fHistTruePrimGammaAsymP);
          fHistTruePrimGammaAsymP->Sumw2();
          fHistTruePrimGammaAsymP->Scale(1./nGammas);
          SetXRange(fHistTruePrimGammaAsymP,minB,maxB);
          SetYRange(fHistTruePrimGammaAsymP,fHistTruePrimGammaAsymP->GetYaxis()->FindBin(0.1),maxYB);
          SetZMinMaxTH2(fHistTruePrimGammaAsymP,minB,maxB,fHistTruePrimGammaAsymP->GetYaxis()->FindBin(0.1),maxYB,kTRUE);
          //-----------------------------------
          DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kTRUE,kTRUE,
                               fHistTruePrimGammaAsymP,"",
                               "#alpha", "#it{p}_{#gamma} (GeV/#it{c})",1,1.4,
                               0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
          SaveCanvasAndWriteHistogram(cvsQuadratic, fHistTruePrimGammaAsymP, Form("%s/Asymmetry_P_TruePrimGamma_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
        }else cout << "INFO: Object |histoTruePrimGammaAsymP| could not be found! Skipping Draw..." << endl;
        //-----------------------------------
        TH2D* fHistTrueSecGammaAsymP = (TH2D*)TopDir->Get("histoTrueSecGammaAsymP");
        if(fHistTrueSecGammaAsymP){
          CheckNEntries(fHistTrueSecGammaAsymP);
          fHistTrueSecGammaAsymP->Sumw2();
          fHistTrueSecGammaAsymP->Scale(1./nGammas);
          SetXRange(fHistTrueSecGammaAsymP,minB,maxB);
          SetYRange(fHistTrueSecGammaAsymP,fHistTrueSecGammaAsymP->GetYaxis()->FindBin(0.1),maxYB);
          SetZMinMaxTH2(fHistTrueSecGammaAsymP,minB,maxB,fHistTrueSecGammaAsymP->GetYaxis()->FindBin(0.1),maxYB,kTRUE);
          //-----------------------------------
          DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kTRUE,kTRUE,
                               fHistTrueSecGammaAsymP,"",
                               "#alpha", "#it{p}_{#gamma} (GeV/#it{c})",1,1.4,
                               0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
          SaveCanvasAndWriteHistogram(cvsQuadratic, fHistTrueSecGammaAsymP, Form("%s/Asymmetry_P_TrueSecGamma_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
        }else cout << "INFO: Object |histoTrueSecGammaAsymP| could not be found! Skipping Draw..." << endl;
        //-----------------------------------
        TH2D* fHistTrueDalitzGammaAsymP = (TH2D*)TopDir->Get("histoTrueDalitzGammaAsymP");
        if(fHistTrueDalitzGammaAsymP){
          CheckNEntries(fHistTrueDalitzGammaAsymP);
          fHistTrueDalitzGammaAsymP->Sumw2();
          fHistTrueDalitzGammaAsymP->Scale(1./nGammas);
          SetXRange(fHistTrueDalitzGammaAsymP,minB,maxB);
          SetYRange(fHistTrueDalitzGammaAsymP,fHistTrueDalitzGammaAsymP->GetYaxis()->FindBin(0.1),maxYB);
          SetZMinMaxTH2(fHistTrueDalitzGammaAsymP,minB,maxB,fHistTrueDalitzGammaAsymP->GetYaxis()->FindBin(0.1),maxYB,kTRUE);
          //-----------------------------------
          DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kTRUE,kTRUE,
                               fHistTrueDalitzGammaAsymP,"",
                               "#alpha", "#it{p}_{#gamma} (GeV/#it{c})",1,1.4,
                               0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
          SaveCanvasAndWriteHistogram(cvsQuadratic, fHistTrueDalitzGammaAsymP, Form("%s/Asymmetry_P_TrueDalitz_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
        }else cout << "INFO: Object |histoTrueDalitzGammaAsymP| could not be found! Skipping Draw..." << endl;
        //-----------------------------------
        if(fHistGammaAsymP && fHistTruePrimGammaAsymP && fHistTrueSecGammaAsymP && fHistTrueDalitzGammaAsymP){
          TH2D* fHistBGAsymP = (TH2D*) fHistGammaAsymP->Clone("histoBGAsymPMonteCarlo");
          fHistBGAsymP->Sumw2();
          fHistBGAsymP->Add(fHistTruePrimGammaAsymP,-1.);
          fHistBGAsymP->Add(fHistTrueSecGammaAsymP,-1.);
          fHistBGAsymP->Add(fHistTrueDalitzGammaAsymP,-1.);
          SetXRange(fHistBGAsymP,minB,maxB);
          SetYRange(fHistBGAsymP,fHistBGAsymP->GetYaxis()->FindBin(0.1),maxYB);
          SetZMinMaxTH2(fHistBGAsymP,minB,maxB,fHistBGAsymP->GetYaxis()->FindBin(0.1),maxYB,kTRUE);
          CheckForNegativeEntries(fHistBGAsymP);
          //-----------------------------------
          DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kTRUE,kTRUE,
                               fHistBGAsymP,"",
                               "#alpha", "#it{p}_{#gamma} (GeV/#it{c})",1,1.4,
                               0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
          SaveCanvasAndWriteHistogram(cvsQuadratic, fHistBGAsymP, Form("%s/Asymmetry_P_BG_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));

          delete fHistBGAsymP;
        }

//-----------------------------------
//-------------------------------------- Photon InvMass vs Pt ----------------------------------------------
//-----------------------------------

        TH2D* fHistGammaInvMassPt = (TH2D*)TopDir->Get("histoGammaInvMassPt");
        if(fHistGammaInvMassPt){
          fHistGammaInvMassPt->Sumw2();
          fHistGammaInvMassPt->Scale(1./nGammas);
          GetMinMaxBin(fHistGammaInvMassPt,minB,maxB);
          GetMinMaxBinY(fHistGammaInvMassPt,minYB,maxYB); maxYB+=5;
          SetXRange(fHistGammaInvMassPt,minB,maxB);
          SetYRange(fHistGammaInvMassPt,fHistGammaInvMassPt->GetYaxis()->FindBin(0.1),maxYB);
          SetZMinMaxTH2(fHistGammaInvMassPt,minB,maxB,fHistGammaInvMassPt->GetYaxis()->FindBin(0.1),maxYB,kTRUE);
          //-----------------------------------
          TH1D* histoGammaInvMass = (TH1D*)fHistGammaInvMassPt->ProjectionX("histoGammaInvMass");
          ConvGammaRebinWithBinCorrection(histoGammaInvMass,1);
          vecInvMass[0].push_back(histoGammaInvMass);
          //-----------------------------------
          DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kTRUE,kTRUE,
                               fHistGammaInvMassPt,"",
                               "M_{e^{+}e^{-}} (GeV/#it{c}^{2})", "#it{p}_{T,#gamma} (GeV/#it{c})",1,1.4,
                               0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
          SaveCanvasAndWriteHistogram(cvsQuadratic, fHistGammaInvMassPt, Form("%s/InvMass_Pt_Photon_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
        }else cout << "INFO: Object |histoGammaInvMassPt| could not be found! Skipping Draw..." << endl;
    //-----------------------------------
        TH2D* fHistTrueGammaInvMassPt = (TH2D*)TopDir->Get("histoTrueGammaInvMassPt");
        if(fHistTrueGammaInvMassPt){
          CheckNEntries(fHistTrueGammaInvMassPt);
          fHistTrueGammaInvMassPt->Sumw2();
          fHistTrueGammaInvMassPt->Scale(1./nGammas);
          SetXRange(fHistTrueGammaInvMassPt,minB,maxB);
          SetYRange(fHistTrueGammaInvMassPt,fHistTrueGammaInvMassPt->GetYaxis()->FindBin(0.1),maxYB);
          SetZMinMaxTH2(fHistTrueGammaInvMassPt,minB,maxB,fHistTrueGammaInvMassPt->GetYaxis()->FindBin(0.1),maxYB,kTRUE);
          //-----------------------------------
          TH1D* histoMonteCarloTrueGammaInvMass = (TH1D*)fHistTrueGammaInvMassPt->ProjectionX("histoMonteCarloTrueGammaInvMass");
          ConvGammaRebinWithBinCorrection(histoMonteCarloTrueGammaInvMass,1);
          vecInvMass[1].push_back(histoMonteCarloTrueGammaInvMass);
          //-----------------------------------
          DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kTRUE,kTRUE,
                               fHistTrueGammaInvMassPt,"",
                               "M_{e^{+}e^{-}} (GeV/#it{c}^{2})", "#it{p}_{T,#gamma} (GeV/#it{c})",1,1.4,
                               0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
          SaveCanvasAndWriteHistogram(cvsQuadratic, fHistTrueGammaInvMassPt, Form("%s/InvMass_Pt_TrueGamma_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
        }else cout << "INFO: Object |histoTrueGammaInvMassPt| could not be found! Skipping Draw..." << endl;
    //-----------------------------------
        TH2D* fHistTrueDalitzGammaInvMassPt = (TH2D*)TopDir->Get("histoTrueDalitzGammaInvMassPt");
        if(fHistTrueDalitzGammaInvMassPt){
          CheckNEntries(fHistTrueDalitzGammaInvMassPt);
          fHistTrueDalitzGammaInvMassPt->Sumw2();
          fHistTrueDalitzGammaInvMassPt->Scale(1./nGammas);
          SetXRange(fHistTrueDalitzGammaInvMassPt,minB,maxB);
          SetYRange(fHistTrueDalitzGammaInvMassPt,fHistTrueDalitzGammaInvMassPt->GetYaxis()->FindBin(0.1),maxYB);
          SetZMinMaxTH2(fHistTrueDalitzGammaInvMassPt,minB,maxB,fHistTrueDalitzGammaInvMassPt->GetYaxis()->FindBin(0.1),maxYB,kTRUE);
          //-----------------------------------
          TH1D* histoMonteCarloTrueDalitzGammaInvMass = (TH1D*)fHistTrueDalitzGammaInvMassPt->ProjectionX("histoMonteCarloTrueDalitzGammaInvMass");
          ConvGammaRebinWithBinCorrection(histoMonteCarloTrueDalitzGammaInvMass,1);
          vecInvMass[2].push_back(histoMonteCarloTrueDalitzGammaInvMass);
          //-----------------------------------
          DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kTRUE,kTRUE,
                               fHistTrueDalitzGammaInvMassPt,"",
                               "M_{e^{+}e^{-}} (GeV/#it{c}^{2})", "#it{p}_{T,#gamma} (GeV/#it{c})",1,1.4,
                               0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
          SaveCanvasAndWriteHistogram(cvsQuadratic, fHistTrueDalitzGammaInvMassPt, Form("%s/InvMass_Pt_TrueDalitz_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
        }else cout << "INFO: Object |histoTrueDalitzGammaInvMassPt| could not be found! Skipping Draw..." << endl;
    //-----------------------------------
        if(fHistGammaInvMassPt && fHistTrueGammaInvMassPt && fHistTrueDalitzGammaInvMassPt){
          TH2D* fHistBGInvMassPt = (TH2D*) fHistGammaInvMassPt->Clone("histoBGCosPointMonteCarloPt");
          fHistBGInvMassPt->Sumw2();
          fHistBGInvMassPt->Add(fHistTrueGammaInvMassPt,-1.);
          fHistBGInvMassPt->Add(fHistTrueDalitzGammaInvMassPt,-1.);
          SetXRange(fHistBGInvMassPt,minB,maxB);
          SetYRange(fHistBGInvMassPt,fHistBGInvMassPt->GetYaxis()->FindBin(0.1),maxYB);
          SetZMinMaxTH2(fHistBGInvMassPt,minB,maxB,fHistBGInvMassPt->GetYaxis()->FindBin(0.1),maxYB,kTRUE);
          CheckForNegativeEntries(fHistBGInvMassPt);
          //-----------------------------------
          TH1D* histoMonteCarloBGInvMass = (TH1D*)fHistBGInvMassPt->ProjectionX("histoMonteCarloBGInvMass");
          ConvGammaRebinWithBinCorrection(histoMonteCarloBGInvMass,1);
          CheckForNegativeEntries(histoMonteCarloBGInvMass);
          vecInvMass[3].push_back(histoMonteCarloBGInvMass);
          //-----------------------------------
          DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kTRUE,kTRUE,
                               fHistBGInvMassPt,"",
                               "M_{e^{+}e^{-}} (GeV/#it{c}^{2})", "#it{p}_{T,#gamma} (GeV/#it{c})",1,1.4,
                               0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
          SaveCanvasAndWriteHistogram(cvsQuadratic, fHistBGInvMassPt, Form("%s/InvMass_Pt_BG_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));

          delete fHistBGInvMassPt;
        }

//-----------------------------------
//-------------------------------------- Photon Qt ----------------------------------------------
//-----------------------------------

        TH2D* fHistGammaQtPt = (TH2D*)TopDir->Get("histoGammaQtPt");
        if(fHistGammaQtPt){
          fHistGammaQtPt->Sumw2();
          fHistGammaQtPt->Scale(1./nGammas);
          GetMinMaxBin(fHistGammaQtPt,minB,maxB);
          GetMinMaxBinY(fHistGammaQtPt,minYB,maxYB); maxYB+=5;
          SetXRange(fHistGammaQtPt,minB,maxB);
          SetYRange(fHistGammaQtPt,fHistGammaQtPt->GetYaxis()->FindBin(0.1),maxYB);
          SetZMinMaxTH2(fHistGammaQtPt,minB,maxB,fHistGammaQtPt->GetYaxis()->FindBin(0.1),maxYB,kTRUE);
          //-----------------------------------
          DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kTRUE,kTRUE,
                               fHistGammaQtPt,"",
                               "#it{q}_{T,#gamma}", "#it{p}_{T,#gamma} (GeV/#it{c})",1,1.4,
                               0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
          SaveCanvasAndWriteHistogram(cvsQuadratic, fHistGammaQtPt, Form("%s/Qt_Pt_Photon_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
        }else cout << "INFO: Object |histoGammaQtPt| could not be found! Skipping Draw..." << endl;
        //-----------------------------------
        TH2D* fHistTruePrimGammaQtPt = (TH2D*)TopDir->Get("histoTruePrimGammaQtPt");
        if(fHistTruePrimGammaQtPt){
          CheckNEntries(fHistTruePrimGammaQtPt);
          fHistTruePrimGammaQtPt->Sumw2();
          fHistTruePrimGammaQtPt->Scale(1./nGammas);
          SetXRange(fHistTruePrimGammaQtPt,minB,maxB);
          SetYRange(fHistTruePrimGammaQtPt,fHistTruePrimGammaQtPt->GetYaxis()->FindBin(0.1),maxYB);
          SetZMinMaxTH2(fHistTruePrimGammaQtPt,minB,maxB,fHistTruePrimGammaQtPt->GetYaxis()->FindBin(0.1),maxYB,kTRUE);
          //-----------------------------------
          DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kTRUE,kTRUE,
                               fHistTruePrimGammaQtPt,"",
                               "#it{q}_{T,#gamma}", "#it{p}_{T,#gamma} (GeV/#it{c})",1,1.4,
                               0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
          SaveCanvasAndWriteHistogram(cvsQuadratic, fHistTruePrimGammaQtPt, Form("%s/Qt_Pt_TruePrimGamma_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
        }else cout << "INFO: Object |histoTruePrimGammaQtPt| could not be found! Skipping Draw..." << endl;
        //-----------------------------------
        TH2D* fHistTrueSecGammaQtPt = (TH2D*)TopDir->Get("histoTrueSecGammaQtPt");
        if(fHistTrueSecGammaQtPt){
          CheckNEntries(fHistTrueSecGammaQtPt);
          fHistTrueSecGammaQtPt->Sumw2();
          fHistTrueSecGammaQtPt->Scale(1./nGammas);
          SetXRange(fHistTrueSecGammaQtPt,minB,maxB);
          SetYRange(fHistTrueSecGammaQtPt,fHistTrueSecGammaQtPt->GetYaxis()->FindBin(0.1),maxYB);
          SetZMinMaxTH2(fHistTrueSecGammaQtPt,minB,maxB,fHistTrueSecGammaQtPt->GetYaxis()->FindBin(0.1),maxYB,kTRUE);
          //-----------------------------------
          DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kTRUE,kTRUE,
                               fHistTrueSecGammaQtPt,"",
                               "#it{q}_{T,#gamma}", "#it{p}_{T,#gamma} (GeV/#it{c})",1,1.4,
                               0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
          SaveCanvasAndWriteHistogram(cvsQuadratic, fHistTrueSecGammaQtPt, Form("%s/Qt_Pt_TrueSecGamma_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
        }else cout << "INFO: Object |histoTrueSecGammaQtPt| could not be found! Skipping Draw..." << endl;
        //-----------------------------------
        TH2D* fHistTrueDalitzGammaQtPt = (TH2D*)TopDir->Get("histoTrueDalitzGammaQtPt");
        if(fHistTrueDalitzGammaQtPt){
          CheckNEntries(fHistTrueDalitzGammaQtPt);
          fHistTrueDalitzGammaQtPt->Sumw2();
          fHistTrueDalitzGammaQtPt->Scale(1./nGammas);
          SetXRange(fHistTrueDalitzGammaQtPt,minB,maxB);
          SetYRange(fHistTrueDalitzGammaQtPt,fHistTrueDalitzGammaQtPt->GetYaxis()->FindBin(0.1),maxYB);
          SetZMinMaxTH2(fHistTrueDalitzGammaQtPt,minB,maxB,fHistTrueDalitzGammaQtPt->GetYaxis()->FindBin(0.1),maxYB,kTRUE);
          //-----------------------------------
          DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kTRUE,kTRUE,
                               fHistTrueDalitzGammaQtPt,"",
                               "#it{q}_{T,#gamma}", "#it{p}_{T,#gamma} (GeV/#it{c})",1,1.4,
                               0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
          SaveCanvasAndWriteHistogram(cvsQuadratic, fHistTrueDalitzGammaQtPt, Form("%s/Qt_Pt_TrueDalitz_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
        }else cout << "INFO: Object |histoTrueDalitzGammaQtPt| could not be found! Skipping Draw..." << endl;
        //-----------------------------------
        if(fHistGammaQtPt && fHistTruePrimGammaQtPt && fHistTrueSecGammaQtPt && fHistTrueDalitzGammaQtPt){
          TH2D* fHistBGQtPt = (TH2D*) fHistGammaQtPt->Clone("histoBGQtMonteCarloPt");
          fHistBGQtPt->Sumw2();
          fHistBGQtPt->Add(fHistTruePrimGammaQtPt,-1.);
          fHistBGQtPt->Add(fHistTrueSecGammaQtPt,-1.);
          fHistBGQtPt->Add(fHistTrueDalitzGammaQtPt,-1.);
          SetXRange(fHistBGQtPt,minB,maxB);
          SetYRange(fHistBGQtPt,fHistBGQtPt->GetYaxis()->FindBin(0.1),maxYB);
          SetZMinMaxTH2(fHistBGQtPt,minB,maxB,fHistBGQtPt->GetYaxis()->FindBin(0.1),maxYB,kTRUE);
          CheckForNegativeEntries(fHistBGQtPt);
          //-----------------------------------
          DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kTRUE,kTRUE,
                               fHistBGQtPt,"",
                               "#it{q}_{T,#gamma}", "#it{p}_{T,#gamma} (GeV/#it{c})",1,1.4,
                               0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
          SaveCanvasAndWriteHistogram(cvsQuadratic, fHistBGQtPt, Form("%s/Qt_Pt_BG_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));

          delete fHistBGQtPt;
        }

	//-----------------------------------
	//-----------------------------------
	//-- Leptons
	//-----------------------------------
	//-----------------------------------

	for(Int_t iL=0; iL<2; iL++){
		TH3D* dEdxEtaP = (TH3D*)TopDir->Get(Form("histo%sdEdxEtaP",lepton[iL].Data()));
        TString labelNLeptons = Form("#frac{1}{N_{e^{%s}}}",charge[iL].Data());
        Double_t nLeptons = 1.;
		if(dEdxEtaP){
            nLeptons = dEdxEtaP->Integral();
			TH2D* histodEdxP = (TH2D*)dEdxEtaP->Project3D("xz");
			histodEdxP->Sumw2();
            histodEdxP->Scale(1./nLeptons);
			GetMinMaxBin(histodEdxP,minB,maxB);
            SetXRange(histodEdxP,minB,maxB);
            SetZMinMaxTH2(histodEdxP,minB,maxB,1,histodEdxP->GetNbinsY(),kTRUE);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                histodEdxP,"",
                                 Form("#it{p}_{e^{%s}} (GeV/#it{c})",charge[iL].Data()),Form("d#it{E}_{e^{%s}-cand} /d#it{x} TPC",charge[iL].Data()),1,1.4,
                                 0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, histodEdxP, Form("%s/dEdxTPC_%s_%s.%s", outputDir.Data(), lepton[iL].Data(), DataSets[i].Data(), suffix.Data()));
			delete histodEdxP;
			//-----------------------------------
			TH1D* NSigmadEdx = (TH1D*)dEdxEtaP->Project3D("x");
			NSigmadEdx->Sumw2();
            NSigmadEdx->Scale(1./nLeptons);
			GetMinMaxBin(NSigmadEdx,minB,maxB);
            SetXRange(NSigmadEdx,minB,maxB);
			DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                 NSigmadEdx,"",Form("n#sigma dEdx %s",lepton[iL].Data()),Form("%s #frac{dN}{dn#sigma dEdx}",labelNLeptons.Data()),1,1,
                                 0.82,0.94,0.03,fCollisionSystem,plotDataSets[i],"");
            SaveCanvasAndWriteHistogram(canvas, NSigmadEdx, Form("%s/nSigmadEdxTPC_%s_%s.%s", outputDir.Data(), lepton[iL].Data(), DataSets[i].Data(), suffix.Data()));
			vecNSigmadEdx[iL].push_back(NSigmadEdx);
        }else cout << Form("INFO: Object |histo%sdEdxEtaP| could not be found! Skipping Draw...",lepton[iL].Data()) << endl;
     //-----------------------------------
        if(nLeptons <= 1.){
          cout << "********************************************************" << endl;
          cout << "WARNING: nLeptons <= 1! Scaling will not work properly..." << endl;
          cout << "********************************************************" << endl;
        }
	//-----------------------------------
		TH3D* NSigmadEdxEtaP = (TH3D*)TopDir->Get(Form("histo%sNSigmadEdxEtaP",lepton[iL].Data()));
		if(NSigmadEdxEtaP){
			TH2D* histoSigmadEdxP = (TH2D*)NSigmadEdxEtaP->Project3D("xz");
			histoSigmadEdxP->Sumw2();
            histoSigmadEdxP->Scale(1./nLeptons);
            GetMinMaxBin(histoSigmadEdxP,minB,maxB);
            SetXRange(histoSigmadEdxP,minB,maxB);
            SetZMinMaxTH2(histoSigmadEdxP,minB,maxB,1,histoSigmadEdxP->GetNbinsY(),kTRUE);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                histoSigmadEdxP,"",
                                 Form("#it{p}_{e^{%s}} (GeV/#it{c})",charge[iL].Data()),Form("#it{n} #sigma_{e^{%s}} d#it{E}/d#it{x} TPC",charge[iL].Data()),1,1.4,
                                 0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, histoSigmadEdxP, Form("%s/nSigma_dEdxTPC_%s_%s.%s", outputDir.Data(), lepton[iL].Data(), DataSets[i].Data(), suffix.Data()));
			delete histoSigmadEdxP;
        }else cout << Form("INFO: Object |histo%sNSigmadEdxEtaP| could not be found! Skipping Draw...",lepton[iL].Data()) << endl;
    //-----------------------------------
        TH2D* TrueNSigmadEdxTPCP = (TH2D*)TopDir->Get(Form("histoTrue%sNSigmadEdxTPCP",lepton[iL].Data()));
        if(TrueNSigmadEdxTPCP){
          TrueNSigmadEdxTPCP->Sumw2();
          TrueNSigmadEdxTPCP->Scale(1./nLeptons);
          GetMinMaxBin(TrueNSigmadEdxTPCP,minB,maxB);
          SetXRange(TrueNSigmadEdxTPCP,minB,maxB);
          SetZMinMaxTH2(TrueNSigmadEdxTPCP,minB,maxB,1,TrueNSigmadEdxTPCP->GetNbinsY(),kTRUE);
          DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                               TrueNSigmadEdxTPCP,"",
                               Form("#it{p}_{e^{%s}} (GeV/#it{c})",charge[iL].Data()),Form("#it{n} #sigma_{e^{%s}} d#it{E}/d#it{x} TPC",charge[iL].Data()),1,1.4,
                               0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
          SaveCanvasAndWriteHistogram(cvsQuadratic, TrueNSigmadEdxTPCP, Form("%s/nSigma_dEdxTPC_True%s_%s.%s", outputDir.Data(), lepton[iL].Data(), DataSets[i].Data(), suffix.Data()));
        }else cout << Form("INFO: Object |histoTrue%sNSigmadEdxTPCP| could not be found! Skipping Draw...",lepton[iL].Data()) << endl;
	//-----------------------------------
		TH3D* ITSdEdxEtaP = (TH3D*)TopDir->Get(Form("histo%sITSdEdxEtaP",lepton[iL].Data()));
		if(ITSdEdxEtaP){
			TH2D* histoITSdEdxEtaP = (TH2D*)ITSdEdxEtaP->Project3D("xz");
			histoITSdEdxEtaP->Sumw2();
            histoITSdEdxEtaP->Scale(1./nLeptons);
			GetMinMaxBin(histoITSdEdxEtaP,minB,maxB);
            SetXRange(histoITSdEdxEtaP,minB,maxB);
            SetZMinMaxTH2(histoITSdEdxEtaP,minB,maxB,1,histoITSdEdxEtaP->GetNbinsY(),kTRUE);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                 histoITSdEdxEtaP,"",
                                 Form("#it{p}_{e^{%s}} (GeV/#it{c})",charge[iL].Data()),Form("d#it{E}_{e^{%s}-cand} /d#it{x} ITS",charge[iL].Data()),1,1.4,
                                 0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, histoITSdEdxEtaP, Form("%s/dEdxITS_%s_%s.%s", outputDir.Data(), lepton[iL].Data(), DataSets[i].Data(), suffix.Data()));
			delete histoITSdEdxEtaP;
        }else cout << Form("INFO: Object |histo%sITSdEdxEtaP| could not be found! Skipping Draw...",lepton[iL].Data()) << endl;
	//-----------------------------------
		TH3D* NSigmaITSEtaP = (TH3D*)TopDir->Get(Form("histo%sNSigmaITSEtaP",lepton[iL].Data()));
		if(NSigmaITSEtaP){
			TH2D* histoNSigmaITSEtaP = (TH2D*)NSigmaITSEtaP->Project3D("xz");
			histoNSigmaITSEtaP->Sumw2();
            histoNSigmaITSEtaP->Scale(1./nLeptons);
			GetMinMaxBin(histoNSigmaITSEtaP,minB,maxB);
            SetXRange(histoNSigmaITSEtaP,minB,maxB);
            SetZMinMaxTH2(histoNSigmaITSEtaP,minB,maxB,1,histoNSigmaITSEtaP->GetNbinsY(),kTRUE);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                 histoNSigmaITSEtaP,"",
                                 Form("#it{p}_{e^{%s}} (GeV/#it{c})",charge[iL].Data()),Form("#it{n} #sigma_{e^{%s}} d#it{E}/d#it{x} ITS",charge[iL].Data()),1,1.4,
                                 0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, histoNSigmaITSEtaP, Form("%s/nSigma_dEdxITS_%s_%s.%s", outputDir.Data(), lepton[iL].Data(), DataSets[i].Data(), suffix.Data()));
			delete histoNSigmaITSEtaP;
        }else cout << Form("INFO: Object |histo%sNSigmaITSEtaP| could not be found! Skipping Draw...",lepton[iL].Data()) << endl;
    //-----------------------------------
        TH2D* TrueNSigmadEdxITSP = (TH2D*)TopDir->Get(Form("histoTrue%sNSigmadEdxITSP",lepton[iL].Data()));
        if(TrueNSigmadEdxITSP){
          TrueNSigmadEdxITSP->Sumw2();
          TrueNSigmadEdxITSP->Scale(1./nLeptons);
          GetMinMaxBin(TrueNSigmadEdxITSP,minB,maxB);
          SetXRange(TrueNSigmadEdxITSP,minB,maxB);
          SetZMinMaxTH2(TrueNSigmadEdxITSP,minB,maxB,1,TrueNSigmadEdxITSP->GetNbinsY(),kTRUE);
          DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                               TrueNSigmadEdxITSP,"",
                               Form("#it{p}_{e^{%s}} (GeV/#it{c})",charge[iL].Data()),Form("#it{n} #sigma_{e^{%s}} d#it{E}/d#it{x} ITS",charge[iL].Data()),1,1.4,
                               0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
          SaveCanvasAndWriteHistogram(cvsQuadratic, TrueNSigmadEdxITSP, Form("%s/nSigma_dEdxITS_True%s_%s.%s", outputDir.Data(), lepton[iL].Data(), DataSets[i].Data(), suffix.Data()));
        }else cout << Form("INFO: Object |histoTrue%sNSigmadEdxITSP| could not be found! Skipping Draw...",lepton[iL].Data()) << endl;
	//-----------------------------------
		TH3D* TOFEtaP = (TH3D*)TopDir->Get(Form("histo%sTOFEtaP",lepton[iL].Data()));
		if(TOFEtaP){
			TH2D* histoSigmadEdxP = (TH2D*)TOFEtaP->Project3D("xz");
			histoSigmadEdxP->Sumw2();
            histoSigmadEdxP->Scale(1./nLeptons);
			GetMinMaxBin(histoSigmadEdxP,minB,maxB);
            SetXRange(histoSigmadEdxP,minB,maxB);
            SetZMinMaxTH2(histoSigmadEdxP,minB,maxB,1,histoSigmadEdxP->GetNbinsY());
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kFALSE,kFALSE,
                                 histoSigmadEdxP,"",
                                 Form("#it{p}_{e^{%s}} (GeV/#it{c})",charge[iL].Data()),Form("#it{t}_{measured}-#it{t}_{expected} e^{%s}",charge[iL].Data()),1,1.4,
                                 0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, histoSigmadEdxP, Form("%s/TOF_%s_%s.%s", outputDir.Data(), lepton[iL].Data(), DataSets[i].Data(), suffix.Data()));
			delete histoSigmadEdxP;
        }else cout << Form("INFO: Object |histo%sTOFEtaP| could not be found! Skipping Draw...",lepton[iL].Data()) << endl;
	//-----------------------------------
		TH3D* NSigmaTOFEtaP = (TH3D*)TopDir->Get(Form("histo%sNSigmaTOFEtaP",lepton[iL].Data()));
		if(NSigmaTOFEtaP){
			TH2D* histoSigmadEdxP = (TH2D*)NSigmaTOFEtaP->Project3D("xz");
			histoSigmadEdxP->Sumw2();
            histoSigmadEdxP->Scale(1./nLeptons);
			GetMinMaxBin(histoSigmadEdxP,minB,maxB);
            SetXRange(histoSigmadEdxP,minB,maxB);
            SetZMinMaxTH2(histoSigmadEdxP,minB,maxB,1,histoSigmadEdxP->GetNbinsY());
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kFALSE,kFALSE,
                                 histoSigmadEdxP,"",
                                 Form("#it{p}_{e^{%s}} (GeV/#it{c})",charge[iL].Data()),Form("#it{n} #sigma_{e^{%s}} TOF",charge[iL].Data()),1,1.4,
                                 0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, histoSigmadEdxP, Form("%s/nSigma_TOF_%s_%s.%s", outputDir.Data(), lepton[iL].Data(), DataSets[i].Data(), suffix.Data()));
			delete histoSigmadEdxP;
        }else cout << Form("INFO: Object |histo%sNSigmaTOFEtaP| could not be found! Skipping Draw...",lepton[iL].Data()) << endl;
    //-----------------------------------
        TH2D* TrueNSigmadEdxTOFP = (TH2D*)TopDir->Get(Form("histoTrue%sNSigmaTOFP",lepton[iL].Data()));
        if(TrueNSigmadEdxTOFP){
          TrueNSigmadEdxTOFP->Sumw2();
          TrueNSigmadEdxTOFP->Scale(1./nLeptons);
          GetMinMaxBin(TrueNSigmadEdxTOFP,minB,maxB); maxB+=5;
          SetXRange(TrueNSigmadEdxTOFP,TrueNSigmadEdxTOFP->GetXaxis()->FindBin(0.05),maxB);
          SetZMinMaxTH2(TrueNSigmadEdxTOFP,TrueNSigmadEdxTOFP->GetYaxis()->FindBin(0.05),maxB,1,TrueNSigmadEdxTOFP->GetNbinsY(),kTRUE);
          DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                               TrueNSigmadEdxTOFP,"",
                               Form("#it{p}_{e^{%s}} (GeV/#it{c})",charge[iL].Data()),Form("#it{n} #sigma_{e^{%s}} TOF",charge[iL].Data()),1,1.4,
                               0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
          SaveCanvasAndWriteHistogram(cvsQuadratic, TrueNSigmadEdxTOFP, Form("%s/nSigma_dEdxTOF_True%s_%s.%s", outputDir.Data(), lepton[iL].Data(), DataSets[i].Data(), suffix.Data()));
        }else cout << Form("INFO: Object |histoTrue%sNSigmadEdxTOFP| could not be found! Skipping Draw...",lepton[iL].Data()) << endl;
    //-----------------------------------
        TH2D* ITSclsPt = (TH2D*)TopDir->Get(Form("histo%sITSClPt",lepton[iL].Data()));
        if(ITSclsPt){
          ITSclsPt->Sumw2();
          ITSclsPt->Scale(1./nLeptons);
          GetMinMaxBin(ITSclsPt,minB,maxB);
          GetMinMaxBinY(ITSclsPt,minYB,maxYB); maxYB+=5;
          SetXRange(ITSclsPt,minB,maxB);
          SetYRange(ITSclsPt,ITSclsPt->GetYaxis()->FindBin(0.05),maxYB);
          SetZMinMaxTH2(ITSclsPt,minB,maxB,ITSclsPt->GetYaxis()->FindBin(0.05),maxYB,kTRUE);
          DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kTRUE,kTRUE,
                               ITSclsPt,"","# ITS cl","#it{p}_{T,e^{-}} (GeV/#it{c})",1,1.4,
                               0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],"");
          SaveCanvasAndWriteHistogram(cvsQuadratic, ITSclsPt, Form("%s/nITSCl_%s_%s.%s", outputDir.Data(), lepton[iL].Data(), DataSets[i].Data(), suffix.Data()));
        }else cout << Form("INFO: Object |histo%sITSClPt| could not be found! Skipping Draw...",lepton[iL].Data()) << endl;
	//-----------------------------------
		TH2D* ITSclsR = (TH2D*)TopDir->Get(Form("histo%sITSClR",lepton[iL].Data()));
		if(ITSclsR){
            ITSclsR->Sumw2();
            TH1D* AllRcls = (TH1D*) ITSclsR->ProjectionX("Rcls");
            AllRcls->Sumw2();
            ConvGammaRebinWithBinCorrection(AllRcls,1);
            AllRcls->Scale(1./nLeptons);
            GetMinMaxBin(AllRcls,minB,maxB);
            SetXRange(AllRcls,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 AllRcls,"","# ITS cl ",Form("%s #frac{dN}{dITS Cl}",labelNLeptons.Data()),1,1,
                                 0.82,0.94,0.03,fCollisionSystem,plotDataSets[i],"");
            SaveCanvasAndWriteHistogram(canvas, AllRcls, Form("%s/%s_ITScl_%s.%s", outputDir.Data(), lepton[iL].Data(), DataSets[i].Data(), suffix.Data()));
            vecITSAllclsR[iL].push_back(AllRcls);
            //------------------------------------------
            for(Int_t iB=0; iB<7; iB++){
                TH1D* Rcls = (TH1D*) ITSclsR->ProjectionY(Form("Rcls%i",iB),iB+1,iB+1);
                Rcls->Sumw2();
                ConvGammaRebinWithBinCorrection(Rcls,1);
                Rcls->Scale(1./nLeptons);
				GetMinMaxBin(Rcls,minB,maxB);
                SetXRange(Rcls,minB,maxB);
				DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                     Rcls,"","R (cm)",Form("#frac{d#it{N}_{e^{%s}}}{#it{N}_{e^{%s}} d#it{R}}",charge[iL].Data(),charge[iL].Data()),1,1,
                                     0.82,0.94,0.03,fCollisionSystem,plotDataSets[i],"");
                SaveCanvasAndWriteHistogram(canvas, Rcls, Form("%s/%s_%iITSClvsR_%s.%s", outputDir.Data(), lepton[iL].Data(), iB, DataSets[i].Data(), suffix.Data()));
                vecITSclsR[iL][iB].push_back(Rcls);
			}
        }else cout << Form("INFO: Object |histo%sITSClR| could not be found! Skipping Draw...",lepton[iL].Data()) << endl;
    //-----------------------------------
        TH2D* TrueITSclsR = (TH2D*)TopDir->Get(Form("histoTrue%sITSClR",lepton[iL].Data()));
        if(TrueITSclsR){
          TrueITSclsR->Sumw2();
          for(Int_t iB=0; iB<7; iB++){
            TH1D* TrueRcls = (TH1D*) TrueITSclsR->ProjectionY(Form("TrueRcls%i",iB),iB+1,iB+1);
            TrueRcls->Sumw2();
            ConvGammaRebinWithBinCorrection(TrueRcls,1);
            TrueRcls->Scale(1./nLeptons);
            GetMinMaxBin(TrueRcls,minB,maxB);
            SetXRange(TrueRcls,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 TrueRcls,"","R (cm)",Form("#frac{d#it{N}_{e^{%s}}}{#it{N}_{e^{%s}} d#it{R}}",charge[iL].Data(),charge[iL].Data()),1,1,
                                 0.82,0.94,0.03,fCollisionSystem,plotDataSets[i],"");
            SaveCanvasAndWriteHistogram(canvas, TrueRcls, Form("%s/True%s_%iITSClvsR_%s.%s", outputDir.Data(), lepton[iL].Data(), iB, DataSets[i].Data(), suffix.Data()));
            vecTrueITSclsR[iL][iB].push_back(TrueRcls);
          }
        }else cout << Form("INFO: Object |histoTrue%sITSClR| could not be found! Skipping Draw...",lepton[iL].Data()) << endl;
	//-----------------------------------
		TH2D* EtaPt = (TH2D*)TopDir->Get(Form("histo%sEtaPt",lepton[iL].Data()));
		if(EtaPt){
            TH1D* Eta = (TH1D*)EtaPt->ProjectionX(Form("%sEta",lepton[iL].Data()),1,EtaPt->GetNbinsY());
            Eta->Sumw2();
            Eta->Scale(1./nLeptons);
			GetMinMaxBin(Eta,minB,maxB);
            SetXRange(Eta,minB-1,maxB+1);
			DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kFALSE,kFALSE,
                                 Eta,"",Form("%s #eta",lepton[iL].Data()),Form("%s #frac{dN}{d#eta}",labelNLeptons.Data()),1,1,
                                 0.82,0.94,0.03,fCollisionSystem,plotDataSets[i],"");
			SaveCanvasAndWriteHistogram(canvas, Eta, Form("%s/%sEta_%s.%s", outputDir.Data(), lepton[iL].Data(), DataSets[i].Data(), suffix.Data()));
			vecEta[iL].push_back(Eta);
			//-----------------------------------
            TH1D* Pt = (TH1D*)EtaPt->ProjectionY(Form("%sPt",lepton[iL].Data()),1,EtaPt->GetNbinsX());
            Pt->Sumw2();
            Pt->Scale(1./nLeptons);
			GetMinMaxBin(Pt,minB,maxB);
            SetXRange(Pt,Pt->GetXaxis()->FindBin(0.05),maxB);
			DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kTRUE,kTRUE,kFALSE,
                                 Pt,"",Form("%s #it{p}_{T} (GeV/#it{c})" ,lepton[iL].Data()),Form("%s #frac{dN}{d#it{p}_{T}}",labelNLeptons.Data()),1,1,
                                 0.82,0.94,0.03,fCollisionSystem,plotDataSets[i],"");
			SaveCanvasAndWriteHistogram(canvas, Pt, Form("%s/%sPt_%s.%s", outputDir.Data(), lepton[iL].Data(), DataSets[i].Data(), suffix.Data()));
			vecPt[iL].push_back(Pt);
        }else cout << Form("INFO: Object |histo%sEtaPt| could not be found! Skipping Draw...",lepton[iL].Data()) << endl;
    //-----------------------------------
        TH2D* FCLPt = (TH2D*)TopDir->Get(Form("histo%sFClPt",lepton[iL].Data()));
        if(FCLPt){
          FCLPt->Sumw2();
          FCLPt->Scale(1./nLeptons);
          GetMinMaxBin(FCLPt,minB,maxB); minB-=2;
          SetXRange(FCLPt,minB,maxB);
          GetMinMaxBinY(FCLPt,minYB,maxYB); minYB-=2;
          SetYRange(FCLPt,minYB,maxYB);
          SetZMinMaxTH2(FCLPt,minB,maxB,minYB,maxYB,kTRUE);
          DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kTRUE,kTRUE,
                               FCLPt,"",
                               Form("TPC Clusters/Findable Clusters e^{%s}",charge[iL].Data()),Form("#it{p}_{T,e^{%s}} (GeV/#it{c})",charge[iL].Data()),1,1.4,
                               0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
          SaveCanvasAndWriteHistogram(cvsQuadratic, FCLPt, Form("%s/%sFindClusterTPC_Pt_%s.%s", outputDir.Data(), lepton[iL].Data(), DataSets[i].Data(), suffix.Data()));
          //-----------------------------------
          SetXRange(FCLPt,1,FCLPt->GetNbinsX());
          TH1D* FCL = (TH1D*) FCLPt->ProjectionX("FCL",1,FCLPt->GetNbinsY());
          GetMinMaxBin(FCL,minB,maxB); minB-=2;
          SetXRange(FCL,minB,maxB);
          DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kFALSE,kFALSE,
                               FCL,"",Form("TPC Clusters/Findable Clusters e^{%s}",charge[iL].Data()),Form("%s #frac{dN}{dRatio Cl}",labelNLeptons.Data()),1,1,
                               0.82,0.94,0.03,fCollisionSystem,plotDataSets[i],"");
          SaveCanvasAndWriteHistogram(canvas, FCL, Form("%s/%sFindTPCClusters_%s.%s", outputDir.Data(), lepton[iL].Data(), DataSets[i].Data(), suffix.Data()));
          vecFCL[iL].push_back(FCL);
        }else cout << Form("INFO: Object |histo%sFClPt| could not be found! Skipping Draw...",lepton[iL].Data()) << endl;
	//-----------------------------------
        TH2D* TPCCLPt = (TH2D*)TopDir->Get(Form("histo%sClPt",lepton[iL].Data()));
        if(TPCCLPt){
          TPCCLPt->Sumw2();
          TPCCLPt->Scale(1./nLeptons);
          GetMinMaxBin(TPCCLPt,minB,maxB); minB-=5;
          SetXRange(TPCCLPt,minB,maxB);
          GetMinMaxBinY(TPCCLPt,minYB,maxYB); minYB-=2;
          SetYRange(TPCCLPt,minYB,maxYB);
          SetZMinMaxTH2(TPCCLPt,minB,maxB,minYB,maxYB,kTRUE);
          DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kTRUE,kTRUE,
                               TPCCLPt,"",
                               Form("TPC Clusters e^{%s}",charge[iL].Data()),Form("#it{p}_{T,e^{%s}} (GeV/#it{c})",charge[iL].Data()),1,1.4,
                               0.65,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
          SaveCanvasAndWriteHistogram(cvsQuadratic, TPCCLPt, Form("%s/%sClusterTPC_Pt_%s.%s", outputDir.Data(), lepton[iL].Data(), DataSets[i].Data(), suffix.Data()));
          //-----------------------------------
          SetXRange(TPCCLPt,1,TPCCLPt->GetNbinsX());
          TH1D* TPCCL = (TH1D*) TPCCLPt->ProjectionX("FCL",1,FCLPt->GetNbinsY());
          GetMinMaxBin(TPCCL,minB,maxB);
          SetXRange(TPCCL,minB-5,maxB+5);
          DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kFALSE,kFALSE,
                               TPCCL,"",Form("TPC Clusters e^{%s}",charge[iL].Data()),Form("%s #frac{dN}{dTPC Cl}",labelNLeptons.Data()),1,1,
                               0.82,0.94,0.03,fCollisionSystem,plotDataSets[i],"");
          SaveCanvasAndWriteHistogram(canvas, TPCCL, Form("%s/%sTPCClusters_%s.%s", outputDir.Data(), lepton[iL].Data(), DataSets[i].Data(), suffix.Data()));
          vecTPCCL[iL].push_back(TPCCL);
        }else cout << Form("INFO: Object |histo%sClPt| could not be found! Skipping Draw...",lepton[iL].Data()) << endl;
  //-----------------------------------

	}
		delete TopDir;
		fOutput->Delete();
		delete fOutput;
	}
//*****************************************************************************************************
//*****************************************************************************************************
//****************************** Drawing Special Histograms *******************************************
//*****************************************************************************************************
//*****************************************************************************************************

    std::vector<TH1D*> temp;

    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "Drawing Special Histograms" << endl;
    
    // compare photon phi on A and C side
    if(vecGammaPhiEtaNeg.size()>0 && vecGammaPhiEtaPos.size()>0){
        Int_t i = 0; //only for data
	temp.push_back(vecGammaPhiEtaNeg.at(i));
	temp.push_back(vecGammaPhiEtaPos.at(i));
	AdjustHistRange(temp,2.0,1.2,kFALSE);
	DrawGammaSetMarker(vecGammaPhiEtaNeg.at(i),20,0.8, kBlack , kBlack);
    DrawPeriodQAHistoTH1(canvas,0.11, 0.02, 0.06, 0.11,kFALSE,kFALSE,kFALSE,
                           vecGammaPhiEtaNeg.at(i),"","Photon #phi","#frac{1}{N_{#gamma}} #frac{dN}{d#phi}",1,1,
               0.82,0.92,0.03,fCollisionSystem,plotDataSets[i],"");
	DrawGammaSetMarker(vecGammaPhiEtaPos.at(i),20,0.8, kRed , kRed);
	vecGammaPhiEtaPos.at(i)->Draw("e,hist,same");
	TLegend *legPhotonPhiSpecial = new TLegend(0.3,0.95,0.95,0.99);
	legPhotonPhiSpecial->SetNColumns(2);
	legPhotonPhiSpecial->SetLineColor(0);
	legPhotonPhiSpecial->SetTextSize(0.03);
	legPhotonPhiSpecial->SetFillColor(0);
	legPhotonPhiSpecial->AddEntry(vecGammaPhiEtaPos.at(i),"A side");
	legPhotonPhiSpecial->AddEntry(vecGammaPhiEtaNeg.at(i),"C side");
	legPhotonPhiSpecial->Draw();
	SaveCanvasOnly(canvas, Form("%s/Comparison/Photon_Phi_EtaNegPos_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
	temp.clear(); 
    }

    if(nSets>1){

      if((Int_t)vecChi2[0].size()>0 && (Int_t)vecChi2[1].size()>0 && (Int_t)vecChi2[2].size()>0 && (Int_t)vecChi2[3].size()>0 && (Int_t)vecChi2[4].size()>0){
        for(Int_t i=0; i<nSets; i++){
          temp.push_back(vecChi2[0].at(i));
          if(i!=0){
            temp.push_back(vecChi2[1].at(i-1));
            temp.push_back(vecChi2[2].at(i-1));
            temp.push_back(vecChi2[3].at(i-1));
            temp.push_back(vecChi2[4].at(i-1));
          }
        }
        AdjustHistRange(temp,5,5,kFALSE);
        for(Int_t i=1; i<nSets; i++){
          DrawGammaSetMarker(vecChi2[0].at(0),20,0.8, kBlack , kBlack);
          DrawPeriodQAHistoTH1(canvas,0.11, 0.02, 0.06, 0.11,kFALSE,kTRUE,kFALSE,
                               vecChi2[0].at(0),"","#chi^{2}_{#gamma}/ndf ","#frac{d#it{N}_{#gamma}}{#it{N}_{#gamma} d#chi^{2}}",1,1,
                               0.82,0.92,0.03,fCollisionSystem,plotDataSets[i],"");
          DrawGammaSetMarker(vecChi2[0].at(i),20,0.8, kRed , kRed);
          vecChi2[0].at(i)->Draw("histsame");
          DrawGammaSetMarker(vecChi2[1].at(i-1),20,0.8, kBlue+2 , kBlue+2);
          vecChi2[1].at(i-1)->Draw("histsame");
          DrawGammaSetMarker(vecChi2[2].at(i-1),20,0.8, kGreen+2 , kGreen+2);
          vecChi2[2].at(i-1)->Draw("histsame");
          DrawGammaSetMarker(vecChi2[3].at(i-1),20,0.8, kOrange+2 , kOrange+2);
          vecChi2[3].at(i-1)->Draw("histsame");
          DrawGammaSetMarker(vecChi2[4].at(i-1),20,0.8, kViolet , kViolet);
          vecChi2[4].at(i-1)->Draw("histsame");

          TLegend *legAdditionalInfoChi2 = new TLegend(0.12,0.95,0.95,0.99);
          legAdditionalInfoChi2->SetNColumns(6);
          legAdditionalInfoChi2->SetLineColor(0);
          //TLegend* legAdditionalInfoChi2 = new TLegend( 0.65,0.65,0.97,0.85);
          legAdditionalInfoChi2->SetTextSize(0.03);
          legAdditionalInfoChi2->SetFillColor(0);
          legAdditionalInfoChi2->AddEntry(vecChi2[0].at(0),("Data"));
          legAdditionalInfoChi2->AddEntry(vecChi2[0].at(i),("MC"));
          legAdditionalInfoChi2->AddEntry(vecChi2[1].at(i-1),("True Prim #gamma"),"l");
          legAdditionalInfoChi2->AddEntry(vecChi2[2].at(i-1),("True Sec #gamma"),"l");
          legAdditionalInfoChi2->AddEntry(vecChi2[3].at(i-1),("True Dalitz"),"l");
          legAdditionalInfoChi2->AddEntry(vecChi2[4].at(i-1),("BG"),"l");
          legAdditionalInfoChi2->Draw();

          SaveCanvas(canvas, Form("%s/Comparison/Photon_Chi2_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()),kFALSE,kTRUE);
          delete legAdditionalInfoChi2;
        }
        DeleteVecTH1D(vecChi2[1]);
        DeleteVecTH1D(vecChi2[2]);
        DeleteVecTH1D(vecChi2[3]);
        DeleteVecTH1D(vecChi2[4]);
      }else cout << "INFO: Skipping drawing of: Photon_Chi2, vector 'vecChi2' is missing some histograms!" << endl;
  //-----------------------------------
  //-----------------------------------
      temp.clear();
      if((Int_t)vecPsi[0].size()>0 && (Int_t)vecPsi[1].size()>0 && (Int_t)vecPsi[2].size()>0 && (Int_t)vecPsi[3].size()>0 && (Int_t)vecPsi[4].size()>0){
        for(Int_t i=0; i<nSets; i++){
          temp.push_back(vecPsi[0].at(i));
          if(i!=0){
            temp.push_back(vecPsi[1].at(i-1));
            temp.push_back(vecPsi[2].at(i-1));
            temp.push_back(vecPsi[3].at(i-1));
            temp.push_back(vecPsi[4].at(i-1));
          }
        }
        AdjustHistRange(temp,5,5,kFALSE);
        for(Int_t i=1; i<nSets; i++){
          DrawGammaSetMarker(vecPsi[0].at(0),20,0.8, kBlack , kBlack);
          DrawPeriodQAHistoTH1(canvas,0.11, 0.02, 0.06, 0.11,kFALSE,kTRUE,kFALSE,
                               vecPsi[0].at(0),"","#psi_{pair} ","#frac{d#it{N}_{#gamma}}{#it{N}_{#gamma} d#psi_{pair}}",1,1,
                               0.82,0.92,0.03,fCollisionSystem,plotDataSets[i],"");
          DrawGammaSetMarker(vecPsi[0].at(i),20,0.8, kRed , kRed);
          vecPsi[0].at(i)->Draw("histsame");
          DrawGammaSetMarker(vecPsi[1].at(i-1),20,0.8, kBlue+2 , kBlue+2);
          vecPsi[1].at(i-1)->Draw("histsame");
          DrawGammaSetMarker(vecPsi[2].at(i-1),20,0.8, kGreen+2 , kGreen+2);
          vecPsi[2].at(i-1)->Draw("histsame");
          DrawGammaSetMarker(vecPsi[3].at(i-1),20,0.8, kOrange+2 , kOrange+2);
          vecPsi[3].at(i-1)->Draw("histsame");
          DrawGammaSetMarker(vecPsi[4].at(i-1),20,0.8, kViolet , kViolet);
          vecPsi[4].at(i-1)->Draw("histsame");

          TLegend *legAdditionalInfoPsi = new TLegend(0.12,0.95,0.95,0.99);
          legAdditionalInfoPsi->SetNColumns(6);
          legAdditionalInfoPsi->SetLineColor(0);
          //TLegend* legAdditionalInfoPsi = new TLegend( 0.65,0.65,0.97,0.85);
          legAdditionalInfoPsi->SetTextSize(0.03);
          legAdditionalInfoPsi->SetFillColor(0);
          legAdditionalInfoPsi->AddEntry(vecPsi[0].at(0),("Data"));
          legAdditionalInfoPsi->AddEntry(vecPsi[0].at(i),("MC"));
          legAdditionalInfoPsi->AddEntry(vecPsi[1].at(i-1),("True Prim #gamma"),"l");
          legAdditionalInfoPsi->AddEntry(vecPsi[2].at(i-1),("True Sec #gamma"),"l");
          legAdditionalInfoPsi->AddEntry(vecPsi[3].at(i-1),("True Dalitz"),"l");
          legAdditionalInfoPsi->AddEntry(vecPsi[4].at(i-1),("BG"),"l");
          legAdditionalInfoPsi->Draw();

          SaveCanvas(canvas, Form("%s/Comparison/Photon_PsiPair_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()),kFALSE,kTRUE);
          delete legAdditionalInfoPsi;
        }
        DeleteVecTH1D(vecPsi[1]);
        DeleteVecTH1D(vecPsi[2]);
        DeleteVecTH1D(vecPsi[3]);
        DeleteVecTH1D(vecPsi[4]);
      }else cout << "INFO: Skipping drawing of: Photon_PsiPair, vector 'vecPsi' is missing some histograms!" << endl;
  //-----------------------------------
  //-----------------------------------
      temp.clear();
      if((Int_t)vecCos[0].size()>0 && (Int_t)vecCos[1].size()>0 && (Int_t)vecCos[2].size()>0 && (Int_t)vecCos[3].size()>0 && (Int_t)vecCos[4].size()>0){
        for(Int_t i=0; i<nSets; i++){
          temp.push_back(vecCos[0].at(i));
          if(i!=0){
            temp.push_back(vecCos[1].at(i-1));
            temp.push_back(vecCos[2].at(i-1));
            temp.push_back(vecCos[3].at(i-1));
            temp.push_back(vecCos[4].at(i-1));
          }
        }
        AdjustHistRange(temp,5,5,kFALSE);
        for(Int_t i=1; i<nSets; i++){
          DrawGammaSetMarker(vecCos[0].at(0),20,0.8, kBlack , kBlack);
          DrawPeriodQAHistoTH1(canvas,0.11, 0.02, 0.06, 0.11,kFALSE,kTRUE,kFALSE,
                               vecCos[0].at(0),"","cos(#theta_{point})","#frac{d#it{N}_{#gamma}}{#it{N}_{#gamma} dcos(#theta_{point})}",1,1,
              0.82,0.92,0.03,fCollisionSystem,plotDataSets[i],"");
          DrawGammaSetMarker(vecCos[0].at(i),20,0.8, kRed , kRed);
          vecCos[0].at(i)->Draw("histsame");
          DrawGammaSetMarker(vecCos[1].at(i-1),20,0.8, kBlue+2 , kBlue+2);
          vecCos[1].at(i-1)->Draw("histsame");
          DrawGammaSetMarker(vecCos[2].at(i-1),20,0.8, kGreen+2 , kGreen+2);
          vecCos[2].at(i-1)->Draw("histsame");
          DrawGammaSetMarker(vecCos[3].at(i-1),20,0.8, kOrange+2 , kOrange+2);
          vecCos[3].at(i-1)->Draw("histsame");
          DrawGammaSetMarker(vecCos[4].at(i-1),20,0.8, kViolet , kViolet);
          vecCos[4].at(i-1)->Draw("histsame");

          TLegend *legAdditionalInfoCos = new TLegend(0.12,0.95,0.95,0.99);
          legAdditionalInfoCos->SetNColumns(6);
          legAdditionalInfoCos->SetLineColor(0);
          //TLegend* legAdditionalInfoCos = new TLegend( 0.65,0.65,0.97,0.85);
          legAdditionalInfoCos->SetTextSize(0.03);
          legAdditionalInfoCos->SetFillColor(0);
          legAdditionalInfoCos->AddEntry(vecCos[0].at(0),("Data"));
          legAdditionalInfoCos->AddEntry(vecCos[0].at(i),("MC"));
          legAdditionalInfoCos->AddEntry(vecCos[1].at(i-1),("True Prim #gamma"),"l");
          legAdditionalInfoCos->AddEntry(vecCos[2].at(i-1),("True Sec #gamma"),"l");
          legAdditionalInfoCos->AddEntry(vecCos[3].at(i-1),("True Dalitz"),"l");
          legAdditionalInfoCos->AddEntry(vecCos[4].at(i-1),("BG"),"l");
          legAdditionalInfoCos->Draw();

          SaveCanvas(canvas, Form("%s/Comparison/Photon_CosPoint_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()),kFALSE,kTRUE);
          delete legAdditionalInfoCos;
        }
        DeleteVecTH1D(vecCos[1]);
        DeleteVecTH1D(vecCos[2]);
        DeleteVecTH1D(vecCos[3]);
        DeleteVecTH1D(vecCos[4]);
      }else cout << "INFO: Skipping drawing of: Photon_CosPoint, vector 'vecCos is missing some histograms!" << endl;
  //-----------------------------------
  //-----------------------------------
      temp.clear();
      if((Int_t)vecGamma[0].size()>0 && (Int_t)vecGamma[1].size()>0 && (Int_t)vecGamma[2].size()>0 && (Int_t)vecGamma[3].size()>0 && (Int_t)vecGamma[4].size()>0){
        for(Int_t i=0; i<nSets; i++){
          temp.push_back(vecGamma[0].at(i));
        }
        AdjustHistRange(temp,5,5,kFALSE);
        for(Int_t i=1; i<nSets; i++){
          DrawGammaSetMarker(vecGamma[0].at(0),20,0.8, kBlack , kBlack);
          DrawPeriodQAHistoTH1(canvas,0.11, 0.02, 0.06, 0.11,kTRUE,kTRUE,kFALSE,
                               vecGamma[0].at(0),"","#it{p}_{T} (GeV/#it{c})","#frac{d#it{N}_{#gamma}}{#it{N}_{#gamma} d#it{p}_{T}} (GeV/#it{c})^{-1}",1,1,
              0.82,0.92,0.03,fCollisionSystem,plotDataSets[i],"");
          DrawGammaSetMarker(vecGamma[0].at(i),20,0.8, kRed , kRed);
          vecGamma[0].at(i)->Draw("histsame");
          DrawGammaSetMarker(vecGamma[1].at(i-1),20,0.8, kBlue+2 , kBlue+2);
          vecGamma[1].at(i-1)->Draw("histsame");
          DrawGammaSetMarker(vecGamma[2].at(i-1),20,0.8, kGreen+2 , kGreen+2);
          vecGamma[2].at(i-1)->Draw("histsame");
          DrawGammaSetMarker(vecGamma[3].at(i-1),20,0.8, kOrange+2 , kOrange+2);
          vecGamma[3].at(i-1)->Draw("histsame");
          DrawGammaSetMarker(vecGamma[4].at(i-1),20,0.8, kViolet , kViolet);
          vecGamma[4].at(i-1)->Draw("histsame");

          TLegend *legAdditionalInfoGamma = new TLegend(0.12,0.95,0.95,0.99);
          legAdditionalInfoGamma->SetNColumns(6);
          legAdditionalInfoGamma->SetLineColor(0);
          //TLegend* legAdditionalInfoGamma = new TLegend( 0.65,0.65,0.97,0.85);
          legAdditionalInfoGamma->SetTextSize(0.03);
          legAdditionalInfoGamma->SetFillColor(0);
          legAdditionalInfoGamma->AddEntry(vecGamma[0].at(0),("Data"));
          legAdditionalInfoGamma->AddEntry(vecGamma[0].at(i),("MC"));
          legAdditionalInfoGamma->AddEntry(vecGamma[1].at(i-1),("True Prim #gamma"),"l");
          legAdditionalInfoGamma->AddEntry(vecGamma[2].at(i-1),("True Sec #gamma"),"l");
          legAdditionalInfoGamma->AddEntry(vecGamma[3].at(i-1),("True Dalitz"),"l");
          legAdditionalInfoGamma->AddEntry(vecGamma[4].at(i-1),("BG"),"l");
          legAdditionalInfoGamma->Draw();

          SaveCanvas(canvas, Form("%s/Comparison/Photon_Pt_AddMCInfo_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()),kTRUE, kTRUE);
          delete legAdditionalInfoGamma;
        }

        for(Int_t i=1; i<nSets; i++){
          TH1D* histoTruePhotonPt = (TH1D*) vecGamma[1].at(i-1)->Clone("histoTruePhotonPt");
          histoTruePhotonPt->Add(vecGamma[2].at(i-1));
          TH1D* histoPurity = (TH1D*) histoTruePhotonPt->Clone("histoPurity");
          histoPurity->Divide(histoPurity,vecGamma[0].at(i),1,1,"B");

          DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kFALSE,kFALSE,
                               histoPurity,"","#it{p}_{T} (GeV/#it{c})","#epsilon_{pur}",1,1,
                               0.82,0.94,0.03,fCollisionSystem,plotDataSets[i],"");
          SaveCanvas(canvas, Form("%s/Comparison/Photon_Purity_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()), kTRUE);

          histoPurity->GetXaxis()->SetRangeUser(0.05,60.);
          histoPurity->GetYaxis()->SetRangeUser(0.9,1.05);
          DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kFALSE,kFALSE,
                               histoPurity,"","#it{p}_{T} (GeV/#it{c})","#epsilon_{pur}",1,1,
                               0.82,0.94,0.03,fCollisionSystem,plotDataSets[i],"");
          SaveCanvas(canvas, Form("%s/Comparison/Photon_Purity_FixRange_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()), kTRUE);

          delete histoTruePhotonPt;
          delete histoPurity;
        }
        DeleteVecTH1D(vecGamma[0]);
        DeleteVecTH1D(vecGamma[1]);
        DeleteVecTH1D(vecGamma[2]);
        DeleteVecTH1D(vecGamma[3]);
        DeleteVecTH1D(vecGamma[4]);
      }else cout << "INFO: Skipping drawing of: Photon_Pt, vector 'vecGamma' is missing some histograms!" << endl;
//-----------------------------------
//-----------------------------------
      temp.clear();
      if((Int_t)vecInvMass[0].size()>0 && (Int_t)vecInvMass[1].size()>0 && (Int_t)vecInvMass[2].size()>0 && (Int_t)vecInvMass[3].size()>0){
        for(Int_t i=0; i<nSets; i++){
          temp.push_back(vecInvMass[0].at(i));
          if(i!=0){
            temp.push_back(vecInvMass[1].at(i-1));
            temp.push_back(vecInvMass[2].at(i-1));
            temp.push_back(vecInvMass[3].at(i-1));
          }
        }
        AdjustHistRange(temp,5,5,kFALSE);
        for(Int_t i=1; i<nSets; i++){
          DrawGammaSetMarker(vecInvMass[0].at(0),20,0.8, kBlack , kBlack);
          DrawPeriodQAHistoTH1(canvas,0.11, 0.02, 0.06, 0.11,kFALSE,kTRUE,kFALSE,
                               vecInvMass[0].at(0),"","M_{e^{+}e^{-}} (GeV/#it{c}^{2})","#frac{d#it{N}_{#gamma}}{#it{N}_{#gamma} d#it{M}_{e^{+}e^{-}}} (GeV/#it{c}^{2})^{-1}",1,1,
              0.82,0.94,0.03,fCollisionSystem,plotDataSets[i],"");
          DrawGammaSetMarker(vecInvMass[0].at(i),20,0.8, kRed , kRed);
          vecInvMass[0].at(i)->Draw("histsame");
          DrawGammaSetMarker(vecInvMass[1].at(i-1),20,0.8, kBlue+2 , kBlue+2);
          vecInvMass[1].at(i-1)->Draw("histsame");
          DrawGammaSetMarker(vecInvMass[2].at(i-1),20,0.8, kOrange+2 , kOrange+2);
          vecInvMass[2].at(i-1)->Draw("histsame");
          DrawGammaSetMarker(vecInvMass[3].at(i-1),20,0.8, kViolet , kViolet);
          vecInvMass[3].at(i-1)->Draw("histsame");

          TLegend *legAdditionalInfoInvMass = new TLegend(0.12,0.95,0.95,0.99);
          legAdditionalInfoInvMass->SetNColumns(6);
          legAdditionalInfoInvMass->SetLineColor(0);
          //TLegend* legAdditionalInfoInvMass = new TLegend( 0.65,0.65,0.97,0.85);
          legAdditionalInfoInvMass->SetTextSize(0.03);
          legAdditionalInfoInvMass->SetFillColor(0);
          legAdditionalInfoInvMass->AddEntry(vecInvMass[0].at(0),("Data"));
          legAdditionalInfoInvMass->AddEntry(vecInvMass[0].at(i),("MC"));
          legAdditionalInfoInvMass->AddEntry(vecInvMass[1].at(i-1),("True #gamma"),"l");
          legAdditionalInfoInvMass->AddEntry(vecInvMass[2].at(i-1),("True Dalitz"),"l");
          legAdditionalInfoInvMass->AddEntry(vecInvMass[3].at(i-1),("BG"),"l");
          legAdditionalInfoInvMass->Draw();

          SaveCanvas(canvas, Form("%s/Comparison/Photon_InvMass_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()),kFALSE,kTRUE);
          delete legAdditionalInfoInvMass;
        }
        DeleteVecTH1D(vecInvMass[1]);
        DeleteVecTH1D(vecInvMass[2]);
        DeleteVecTH1D(vecInvMass[3]);
      }else cout << "INFO: Skipping drawing of: Photon_InvMass, vector 'vecInvMass' is missing some histograms!" << endl;
//-----------------------------------
//-----------------------------------
      for(Int_t iC=0; iC<2; iC++){
        if((Int_t)vecITSclsR[iC][0].size()>0 && (Int_t)vecITSclsR[iC][1].size()>0 && (Int_t)vecITSclsR[iC][2].size()>0 &&
           (Int_t)vecITSclsR[iC][3].size()>0 && (Int_t)vecITSclsR[iC][4].size()>0 && (Int_t)vecITSclsR[iC][5].size()>0 &&
           (Int_t)vecITSclsR[iC][6].size()>0){

          for(Int_t i=0; i<nSets; i++){
            temp.clear();
            for(Int_t iB=0; iB<7; iB++) temp.push_back(vecITSclsR[iC][iB].at(i));
            AdjustHistRange(temp,5,5,kFALSE);
            DrawGammaSetMarker(vecITSclsR[iC][0].at(i),20,0.8, kBlack , kBlack);
            DrawPeriodQAHistoTH1(canvas,0.11, 0.02, 0.06, 0.11,kFALSE,kTRUE,kFALSE,
                                 vecITSclsR[iC][0].at(i),"","R (cm)",Form("#frac{d#it{N}_{e^{%s}}}{#it{N}_{e^{%s}} d#it{R}}",charge[iC].Data(),charge[iC].Data()),1,1,
                0.82,0.94,0.03,fCollisionSystem,plotDataSets[i],"");
            DrawGammaSetMarker(vecITSclsR[iC][1].at(i),20,0.8, kRed , kRed);
            vecITSclsR[iC][1].at(i)->Draw("histsame");
            DrawGammaSetMarker(vecITSclsR[iC][2].at(i),20,0.8, kBlue+2 , kBlue+2);
            vecITSclsR[iC][2].at(i)->Draw("histsame");
            DrawGammaSetMarker(vecITSclsR[iC][3].at(i),20,0.8, kOrange+2 , kOrange+2);
            vecITSclsR[iC][3].at(i)->Draw("histsame");
            DrawGammaSetMarker(vecITSclsR[iC][4].at(i),20,0.8, kViolet+2 , kViolet+2);
            vecITSclsR[iC][4].at(i)->Draw("histsame");
            DrawGammaSetMarker(vecITSclsR[iC][5].at(i),20,0.8, kGreen+2 , kGreen+2);
            vecITSclsR[iC][5].at(i)->Draw("histsame");
            DrawGammaSetMarker(vecITSclsR[iC][6].at(i),20,0.8, kPink+2 , kPink+2);
            vecITSclsR[iC][6].at(i)->Draw("histsame");

            TLegend *legAdditionalInfoITSClvsR = new TLegend(0.12,0.95,0.95,0.99);
            legAdditionalInfoITSClvsR->SetNColumns(7);
            legAdditionalInfoITSClvsR->SetLineColor(0);
            legAdditionalInfoITSClvsR->SetTextSize(0.03);
            legAdditionalInfoITSClvsR->SetFillColor(0);
            legAdditionalInfoITSClvsR->AddEntry(vecITSclsR[iC][0].at(i),("0 ITS cl"),"l");
            legAdditionalInfoITSClvsR->AddEntry(vecITSclsR[iC][1].at(i),("1 ITS cl"),"l");
            legAdditionalInfoITSClvsR->AddEntry(vecITSclsR[iC][2].at(i),("2 ITS cl"),"l");
            legAdditionalInfoITSClvsR->AddEntry(vecITSclsR[iC][3].at(i),("3 ITS cl"),"l");
            legAdditionalInfoITSClvsR->AddEntry(vecITSclsR[iC][4].at(i),("4 ITS cl"),"l");
            legAdditionalInfoITSClvsR->AddEntry(vecITSclsR[iC][5].at(i),("5 ITS cl"),"l");
            legAdditionalInfoITSClvsR->AddEntry(vecITSclsR[iC][6].at(i),("6 ITS cl"),"l");
            legAdditionalInfoITSClvsR->Draw();

            SaveCanvas(canvas, Form("%s/Comparison/%s_ITSClvsR_%s.%s", outputDir.Data(), lepton[iC].Data(), DataSets[i].Data(), suffix.Data()),kFALSE,kTRUE);
            delete legAdditionalInfoITSClvsR;
          }
        }else cout << Form("INFO: Skipping drawing of: %s_ITSClvsR, vector 'vecITSclsR' is missing some histograms!",lepton[iC].Data()) << endl;

        if((Int_t)vecTrueITSclsR[iC][0].size()>0 && (Int_t)vecTrueITSclsR[iC][1].size()>0 && (Int_t)vecTrueITSclsR[iC][2].size()>0 &&
           (Int_t)vecTrueITSclsR[iC][3].size()>0 && (Int_t)vecTrueITSclsR[iC][4].size()>0 && (Int_t)vecTrueITSclsR[iC][5].size()>0 &&
           (Int_t)vecTrueITSclsR[iC][6].size()>0){
          for(Int_t i=1; i<nSets; i++){
            temp.clear();
            Int_t iMC = i-1;
            for(Int_t iB=0; iB<7; iB++) temp.push_back(vecTrueITSclsR[iC][iB].at(iMC));
            AdjustHistRange(temp,5,5,kFALSE);
            DrawGammaSetMarker(vecTrueITSclsR[iC][0].at(iMC),20,0.8, kBlack , kBlack);
            DrawPeriodQAHistoTH1(canvas,0.11, 0.02, 0.06, 0.11,kFALSE,kTRUE,kFALSE,
                                 vecTrueITSclsR[iC][0].at(iMC),"","R (cm)",Form("#frac{d#it{N}_{e^{%s}}}{#it{N}_{e^{%s}} d#it{R}}",charge[iC].Data(),charge[iC].Data()),1,1,
                0.82,0.94,0.03,fCollisionSystem,plotDataSets[i],"");
            DrawGammaSetMarker(vecTrueITSclsR[iC][1].at(iMC),20,0.8, kRed , kRed);
            vecTrueITSclsR[iC][1].at(iMC)->Draw("histsame");
            DrawGammaSetMarker(vecTrueITSclsR[iC][2].at(iMC),20,0.8, kBlue+2 , kBlue+2);
            vecTrueITSclsR[iC][2].at(iMC)->Draw("histsame");
            DrawGammaSetMarker(vecTrueITSclsR[iC][3].at(iMC),20,0.8, kOrange+2 , kOrange+2);
            vecTrueITSclsR[iC][3].at(iMC)->Draw("histsame");
            DrawGammaSetMarker(vecTrueITSclsR[iC][4].at(iMC),20,0.8, kViolet+2 , kViolet+2);
            vecTrueITSclsR[iC][4].at(iMC)->Draw("histsame");
            DrawGammaSetMarker(vecTrueITSclsR[iC][5].at(iMC),20,0.8, kGreen+2 , kGreen+2);
            vecTrueITSclsR[iC][5].at(iMC)->Draw("histsame");
            DrawGammaSetMarker(vecTrueITSclsR[iC][6].at(iMC),20,0.8, kPink+2 , kPink+2);
            vecTrueITSclsR[iC][6].at(iMC)->Draw("histsame");

            TLegend *legAdditionalInfoITSClvsR = new TLegend(0.12,0.95,0.95,0.99);
            legAdditionalInfoITSClvsR->SetNColumns(7);
            legAdditionalInfoITSClvsR->SetLineColor(0);
            legAdditionalInfoITSClvsR->SetTextSize(0.03);
            legAdditionalInfoITSClvsR->SetFillColor(0);
            legAdditionalInfoITSClvsR->AddEntry(vecTrueITSclsR[iC][0].at(iMC),("0 ITS cl"),"l");
            legAdditionalInfoITSClvsR->AddEntry(vecTrueITSclsR[iC][1].at(iMC),("1 ITS cl"),"l");
            legAdditionalInfoITSClvsR->AddEntry(vecTrueITSclsR[iC][2].at(iMC),("2 ITS cl"),"l");
            legAdditionalInfoITSClvsR->AddEntry(vecTrueITSclsR[iC][3].at(iMC),("3 ITS cl"),"l");
            legAdditionalInfoITSClvsR->AddEntry(vecTrueITSclsR[iC][4].at(iMC),("4 ITS cl"),"l");
            legAdditionalInfoITSClvsR->AddEntry(vecTrueITSclsR[iC][5].at(iMC),("5 ITS cl"),"l");
            legAdditionalInfoITSClvsR->AddEntry(vecTrueITSclsR[iC][6].at(iMC),("6 ITS cl"),"l");
            legAdditionalInfoITSClvsR->Draw();

            SaveCanvas(canvas, Form("%s/Comparison/True%s_ITSClvsR_%s.%s", outputDir.Data(), lepton[iC].Data(), DataSets[i].Data(), suffix.Data()),kFALSE,kTRUE);
            delete legAdditionalInfoITSClvsR;
          }
        }else cout << Form("INFO: Skipping drawing of: True%s_ITSClvsR, vector 'vecTrueITSclsR' is missing some histograms!", lepton[iC].Data()) << endl;
      }
    }

//*****************************************************************************************************
//*****************************************************************************************************
//****************************** Comparison Histograms ************************************************
//*****************************************************************************************************
//*****************************************************************************************************

    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "Drawing Comparison Histograms" << endl;

	//-----------------------------------
	//-----------------------------------
	//-- Gammas
	//-----------------------------------
	//-----------------------------------

	//-----------------------------------
		DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kTRUE,kTRUE,kFALSE,
                             vecGammaPt,"","Photon #it{p}_{T} (GeV/#it{c})","#frac{1}{N_{#gamma}} #frac{dN}{d#it{p}_{T}}",1,1.1,
							 labelData, colorCompare, kTRUE, 2, 5, kTRUE,
                             0.82,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
        SaveCanvas(canvas, Form("%s/Comparison/Photon_Pt.%s", outputDir.Data(), suffix.Data()),kTRUE,kTRUE);

        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kTRUE,kFALSE,kFALSE,
                             vecGammaPt,"","Photon #it{p}_{T} (GeV/#it{c})","#frac{1}{N_{#gamma}} #frac{dN}{d#it{p}_{T}}",1,1.1,
							 labelData, colorCompare, kTRUE, 1.1, 1.1, kTRUE,
                             0.82,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Photon_Pt.%s", outputDir.Data(), suffix.Data()),kTRUE);
		DeleteVecTH1D(vecGammaPt);
	//-----------------------------------
		DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                             vecGammaEta,"","Photon #eta","#frac{1}{N_{#gamma}} #frac{dN}{d#eta}",1,1.1,
							 labelData, colorCompare, kTRUE, 1.1, 1.1, kTRUE,
                             0.82,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
        SaveCanvas(canvas, Form("%s/Comparison/Photon_Eta.%s", outputDir.Data(), suffix.Data()));

		DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                             vecGammaEta,"","Photon #eta","#frac{1}{N_{#gamma}} #frac{dN}{d#eta}",1,1.1,
							 labelData, colorCompare, kTRUE, 1.1, 1.1, kTRUE,
                             0.82,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Photon_Eta.%s", outputDir.Data(), suffix.Data()));
		DeleteVecTH1D(vecGammaEta);
	//-----------------------------------
		DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                             vecGammaAlpha,"","Photon #alpha","#frac{1}{N_{#gamma}} #frac{dN}{d#alpha}",1,1.1,
							 labelData, colorCompare, kTRUE, 1.1, 1.1, kTRUE,
                             0.82,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
        SaveCanvas(canvas, Form("%s/Comparison/Photon_Alpha.%s", outputDir.Data(), suffix.Data()));

		DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                             vecGammaAlpha,"","Photon #alpha","#frac{1}{N_{#gamma}} #frac{dN}{d#alpha}",1,1.1,
							 labelData, colorCompare, kTRUE, 1.1, 1.1, kTRUE,
                             0.82,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Photon_Alpha.%s", outputDir.Data(), suffix.Data()));
		DeleteVecTH1D(vecGammaAlpha);
	//-----------------------------------
		DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                             vecGammaPhi,"","Photon #phi","#frac{1}{N_{#gamma}} #frac{dN}{d#phi}",1,1.1,
							 labelData, colorCompare, kTRUE, 1.1, 1.1, kTRUE,
                             0.82,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
        SaveCanvas(canvas, Form("%s/Comparison/Photon_Phi.%s", outputDir.Data(), suffix.Data()));

		DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                             vecGammaPhi,"","Photon #phi","#frac{1}{N_{#gamma}} #frac{dN}{d#phi}",1,1.1,
							 labelData, colorCompare, kTRUE, 1.1, 1.1, kTRUE,
                             0.82,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Photon_Phi.%s", outputDir.Data(), suffix.Data()));
		DeleteVecTH1D(vecGammaPhi);
	//-----------------------------------
		DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                             vecGammaPhiEtaNeg,"","Photon #phi C side","#frac{1}{N_{#gamma}} #frac{dN}{d#phi}",1,1.1,
							 labelData, colorCompare, kTRUE, 1.1, 1.1, kTRUE,
                             0.82,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
        SaveCanvas(canvas, Form("%s/Comparison/Photon_Phi_EtaNeg.%s", outputDir.Data(), suffix.Data()));

		DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                             vecGammaPhiEtaNeg,"","Photon #phi C side","#frac{1}{N_{#gamma}} #frac{dN}{d#phi}",1,1.1,
							 labelData, colorCompare, kTRUE, 1.1, 1.1, kTRUE,
                             0.82,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Photon_Phi_EtaNeg.%s", outputDir.Data(), suffix.Data()));
		DeleteVecTH1D(vecGammaPhiEtaNeg);
	//-----------------------------------
		DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                             vecGammaPhiEtaPos,"","Photon #phi A side","#frac{1}{N_{#gamma}} #frac{dN}{d#phi}",1,1.1,
							 labelData, colorCompare, kTRUE, 1.1, 1.1, kTRUE,
                             0.82,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
        SaveCanvas(canvas, Form("%s/Comparison/Photon_Phi_EtaPos.%s", outputDir.Data(), suffix.Data()));

		DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                             vecGammaPhiEtaPos,"","Photon #phi A side","#frac{1}{N_{#gamma}} #frac{dN}{d#phi}",1,1.1,
							 labelData, colorCompare, kTRUE, 1.1, 1.1, kTRUE,
                             0.82,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Photon_Phi_EtaPos.%s", outputDir.Data(), suffix.Data()));
		DeleteVecTH1D(vecGammaPhiEtaPos);
    //-----------------------------------
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecChi2[0],"","#chi^{2}_{#gamma}/ndf","#frac{1}{N_{#gamma}} #frac{dN}{d#chi^{2}}",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                    0.82,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
        SaveCanvas(canvas, Form("%s/Comparison/Chi2NDF.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);

        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                         vecChi2[0],"","#chi^{2}_{#gamma}/ndf","#frac{1}{N_{#gamma}} #frac{dN}{d#chi^{2}}",1,1.1,
                                         labelData, colorCompare, kTRUE, 1.1, 1.1, kTRUE,
                                         0.82,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Chi2NDF.%s", outputDir.Data(), suffix.Data()));
        DeleteVecTH1D(vecChi2[0]);
    //-----------------------------------
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecPsi[0],"","#psi_{pair}","#frac{1}{N_{#gamma}} #frac{dN}{d#psi_{pair}}",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                    0.82,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
        SaveCanvas(canvas, Form("%s/Comparison/PsiPair.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);

        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                         vecPsi[0],"","#psi_{pair}","#frac{1}{N_{#gamma}} #frac{dN}{d#psi_{pair}}",1,1.1,
                                         labelData, colorCompare, kTRUE, 1.1, 1.1, kTRUE,
                                         0.82,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_PsiPair.%s", outputDir.Data(), suffix.Data()));
        DeleteVecTH1D(vecPsi[0]);
    //-----------------------------------
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecCos[0],"","cos(#theta_{point})","#frac{1}{N_{#gamma}} #frac{d#it{N}_{#gamma}}{dcos(#theta_{point})}",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.82,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
        SaveCanvas(canvas, Form("%s/Comparison/CosPoint.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);

        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                         vecCos[0],"","cos(#theta_{point})","#frac{1}{N_{#gamma}} #frac{d#it{N}_{#gamma}}{dcos(#theta_{point})}",1,1.1,
                                         labelData, colorCompare, kTRUE, 1.1, 1.1, kFALSE,
                                         0.82,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_CosPoint.%s", outputDir.Data(), suffix.Data()));
        DeleteVecTH1D(vecCos[0]);
    //-----------------------------------
        GetMinMaxBin(vecInvMass[0].at(0),minB,maxB);
        for(Int_t i=0; i<nSets; i++) SetXRange(vecInvMass[0].at(i),minB,maxB+5);
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecInvMass[0],"","M_{e^{+}e^{-}} (GeV/#it{c}^{2})", "#frac{d#it{N}_{#gamma}}{#it{N}_{#gamma} d#it{M}_{e^{+}e^{-}}} (GeV/#it{c}^{2})^{-1}",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.82,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
        SaveCanvas(canvas, Form("%s/Comparison/InvMass.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);

        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                         vecInvMass[0],"","M_{e^{+}e^{-}} (GeV/#it{c}^{2})", "#frac{d#it{N}_{#gamma}}{#it{N}_{#gamma} d#it{M}_{e^{+}e^{-}}} (GeV/#it{c}^{2})^{-1}",1,1.1,
                                         labelData, colorCompare, kTRUE, 1.1, 1.1, kFALSE,
                                         0.82,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_InvMass.%s", outputDir.Data(), suffix.Data()));
        DeleteVecTH1D(vecInvMass[0]);
	//-----------------------------------
	//-----------------------------------
	//-- Leptons
	//-----------------------------------
	//-----------------------------------

	for(Int_t iL=0; iL<2; iL++){
        TString labelNLeptons = Form("#frac{1}{N_{e^{%s}}}",charge[iL].Data());
		DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                             vecEta[iL],"",Form("%s #eta",lepton[iL].Data()),Form("%s #frac{dN}{d#eta}",labelNLeptons.Data()),1,1.1,
							 labelData, colorCompare, kTRUE, 1.1, 1.1, kTRUE,
                             0.82,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
        SaveCanvas(canvas, Form("%s/Comparison/%s_Eta.%s", outputDir.Data(),lepton[iL].Data(), suffix.Data()));

		DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                             vecEta[iL],"",Form("%s #eta",lepton[iL].Data()),Form("%s #frac{dN}{d#eta}",labelNLeptons.Data()),1,1.1,
							 labelData, colorCompare, kTRUE, 1.1, 1.1, kTRUE,
                             0.82,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_%s_Eta.%s", outputDir.Data(),lepton[iL].Data(), suffix.Data()));
		DeleteVecTH1D(vecEta[iL]);
	//-----------------------------------
		DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kTRUE,kTRUE,kFALSE,
                             vecPt[iL],"",Form("%s #it{p}_{T} (GeV/#it{c})",lepton[iL].Data()),Form("%s #frac{dN}{d#it{p}_{T}}",labelNLeptons.Data()),1,1.1,
							 labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                             0.82,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
        SaveCanvas(canvas, Form("%s/Comparison/%s_Pt.%s", outputDir.Data(),lepton[iL].Data(), suffix.Data()),kTRUE,kTRUE);

        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kTRUE,kFALSE,kFALSE,
                             vecPt[iL],"",Form("%s #it{p}_{T} (GeV/#it{c})",lepton[iL].Data()),Form("%s #frac{dN}{d#it{p}_{T}}",labelNLeptons.Data()),1,1.1,
							 labelData, colorCompare, kTRUE, 1.1, 1.1, kTRUE,
                             0.82,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_%s_Pt.%s", outputDir.Data(),lepton[iL].Data(), suffix.Data()),kTRUE);
		DeleteVecTH1D(vecPt[iL]);
	//-----------------------------------
		DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                             vecNSigmadEdx[iL],"",Form("%s n#sigma dEdx",lepton[iL].Data()),Form("%s #frac{dN}{dn#sigma dEdx}",labelNLeptons.Data()),1,1.1,
							 labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                             0.82,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
        SaveCanvas(canvas, Form("%s/Comparison/%s_NSigdEdx.%s", outputDir.Data(),lepton[iL].Data(), suffix.Data()), kFALSE, kTRUE);

		DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                             vecNSigmadEdx[iL],"",Form("%s n#sigma dEdx",lepton[iL].Data()),Form("%s #frac{dN}{dn#sigma dEdx}",labelNLeptons.Data()),1,1.1,
							 labelData, colorCompare, kTRUE, 1.1, 1.1, kTRUE,
                             0.82,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_%s_NSigdEdx.%s", outputDir.Data(),lepton[iL].Data(), suffix.Data()));
		DeleteVecTH1D(vecNSigmadEdx[iL]);
    //-----------------------------------
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                             vecITSAllclsR[iL],"","# ITS cl ",Form("%s #frac{dN}{dITS Cl}",labelNLeptons.Data()),1,1.1,
                             labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                             0.82,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
        SaveCanvas(canvas, Form("%s/Comparison/%s_ITSCl.%s", outputDir.Data(),lepton[iL].Data(), suffix.Data()), kFALSE, kTRUE);

        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                             vecITSAllclsR[iL],"","# ITS cl ",Form("%s #frac{dN}{dITS Cl}",labelNLeptons.Data()),1,1.1,
                             labelData, colorCompare, kTRUE, 1.1, 1.1, kTRUE,
                             0.82,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_%s_ITSCl.%s", outputDir.Data(),lepton[iL].Data(), suffix.Data()));
        DeleteVecTH1D(vecITSAllclsR[iL]);
	//-----------------------------------
        for(Int_t iB=0; iB<7; iB++){
			DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                 vecITSclsR[iL][iB],"","R (cm)",Form("%s #frac{d#it{N}_{e^{%s}}}{d#it{R}}",labelNLeptons.Data(),charge[iL].Data()),1,1.1,
								 labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                 0.82,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
            SaveCanvas(canvas, Form("%s/Comparison/%s_%iITSClvsR.%s", outputDir.Data(),lepton[iL].Data(), iB, suffix.Data()), kFALSE, kTRUE);

			DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                 vecITSclsR[iL][iB],"","R (cm)",Form("%s #frac{d#it{N}_{e^{%s}}}{d#it{R}}",labelNLeptons.Data(),charge[iL].Data()),1,1.1,
								 labelData, colorCompare, kTRUE, 1.1, 1.1, kTRUE,
                                 0.82,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
            SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_%s_%iITSClvsR.%s", outputDir.Data(),lepton[iL].Data(), iB, suffix.Data()));
            DeleteVecTH1D(vecITSclsR[iL][iB]);
		}
	//-----------------------------------
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                             vecFCL[iL],"",Form("TPC Clusters/Findable Clusters e^{%s}",charge[iL].Data()),Form("%s #frac{dN}{dRatio Cl}",labelNLeptons.Data()),1,1.1,
                             labelData, colorCompare, kTRUE, 1.1, 1.1, kTRUE,
                             0.15,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
        SaveCanvas(canvas, Form("%s/Comparison/%s_FindTPCClusters.%s", outputDir.Data(),lepton[iL].Data(), suffix.Data()));

        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                             vecFCL[iL],"",Form("TPC Clusters/Findable Clusters e^{%s}",charge[iL].Data()),Form("%s #frac{dN}{dRatio Cl}",labelNLeptons.Data()),1,1.1,
                             labelData, colorCompare, kTRUE, 1.1, 1.1, kTRUE,
                             0.82,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_%s_FindTPCClusters.%s", outputDir.Data(),lepton[iL].Data(), suffix.Data()));
        DeleteVecTH1D(vecFCL[iL]);
   //-----------------------------------
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                    vecTPCCL[iL],"",Form("TPC Clusters e^{%s}",charge[iL].Data()),Form("%s #frac{dN}{dTPC Cl}",labelNLeptons.Data()),1,1.1,
                                    labelData, colorCompare, kTRUE, 1.1, 1.1, kTRUE,
                                    0.15,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
        SaveCanvas(canvas, Form("%s/Comparison/%s_TPCClusters.%s", outputDir.Data(),lepton[iL].Data(), suffix.Data()));

        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecTPCCL[iL],"",Form("TPC Clusters e^{%s}",charge[iL].Data()),Form("%s #frac{dN}{dTPC Cl}",labelNLeptons.Data()),1,1.1,
                                         labelData, colorCompare, kTRUE, 1.1, 1.1, kTRUE,
                                         0.82,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_%s_TPCClusters.%s", outputDir.Data(),lepton[iL].Data(), suffix.Data()));
        DeleteVecTH1D(vecTPCCL[iL]);
	}

	fLog.close();

    delete[] fTrigger;

    delete[] nEvents;

	delete cvsQuadratic;
	delete canvas;

	TH1::AddDirectory(kTRUE);

	cout << "Done with PhotonQA" << endl;
	cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

	return;

}//end
