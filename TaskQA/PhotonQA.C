/*******************************************************************************
 ******  provided by Gamma Conversion Group, PWGGA,                        *****
 ******     Daniel Muehlheim, d.muehlheim@cern.ch                          *****
 ******     Nicolas Schmidt, nicolas.schmidt@cern.ch                       *****
 ******     Friederike Bock, fbock@cern.ch                                 *****
 *******************************************************************************/

#include "QA.h"
#include "../TaskQA/BuildHistogramsForGammaQAAdvV3.C"


void PalBW(){
    static Int_t  colors[50];
    static Bool_t initialized = kFALSE;

    Double_t stopsBW[5] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t redBW[5]   = { 1.00, 0.84, 0.61, 0.34, 0.00 };
    Double_t greenBW[5] = { 1.00, 0.84, 0.61, 0.34, 0.00 };
    Double_t blueBW[5]  = { 1.00, 0.84, 0.61, 0.34, 0.00 };

    if(!initialized){
        Int_t FI = TColor::CreateGradientColorTable(5, stopsBW, redBW, greenBW, blueBW, 50);
        for (int i=0; i<50; i++) colors[i] = FI+i;
        initialized = kTRUE;
        return;
    }
    gStyle->SetPalette(50,colors);
}

void PhotonQA(
                Int_t nSetsIn,
                TString fEnergyFlag,
                TString* DataSets,
                TString* plotDataSets,
                TString* pathDataSets,
                Int_t mode                      = 2,
                Int_t cutNr                     = 0,         // if -1: you have to choose number at runtime, if 0: the first or only number is chosen
                TString suffix                  = "eps",
                TString labelData               = "Data",
                Bool_t addSubfolder             = kFALSE,
                TString cutTreeProjection       = "",
                Bool_t useConsistentCut         = kFALSE
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
    TString fDate                   = ReturnDateString();
    TString fTextMeasurement        = Form("#pi^{0} #rightarrow #gamma#gamma");
    TString fTextMeasurementEta     = Form("#eta #rightarrow #gamma#gamma");
    TString fTextMeasurementMeson[2]={fTextMeasurement,fTextMeasurementEta};

    const Int_t maxSets = 12;
    //nSets == 0 is always data!

    if(nSets>maxSets){
        cout << "Maximum hardcoded number of Data Sets in PhotonQA: " << maxSets << endl;
        cout << "You have chosen: " << nSets << ", returning!" << endl;
        return;
    }
    Int_t textSizeLabelsPixel                   = 1200*0.04;
    Color_t colorCompare[maxSets] = {kBlack, kRed+1, kMagenta+2, 807, 800, kGreen+2, kCyan+2, kBlue+1, kOrange+2, kAzure, kViolet, kGray+1};
    TString nameMainDir[maxSets];
    TString nameCutsPQA[maxSets];
    TString nameCutsPQAEvent[maxSets];
    TString nameCutsPQAmain[maxSets];
    vector <TString> mainDir[maxSets];
    vector <TString> cutsPQA[maxSets];
    Bool_t isMC [maxSets];

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
    if(fMode == 4 || fMode == 5 || fMode == 10 || fMode == 11 || fMode == 12){ cout << "Returning, given mode contains no PCM information: " << fMode << endl; return;}

    for(Int_t i=0; i<nSets; i++){
        TFile* fFile = new TFile(pathDataSets[i].Data(),"READ");
        if(fFile->IsZombie()){cout << "ERROR: File " << pathDataSets[i].Data() << " could not be openend! Returning..." << endl; return;}
        else{
            cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
            cout << "Processing file: " << pathDataSets[i].Data() << endl;
            TKey *key;
            TIter next(fFile->GetListOfKeys());
            while ((key=(TKey*)next())){
                cout << Form(" - found TopDir: %s",key->GetName()) << endl;
                nameMainDir[i] = key->GetName();
                if(nameMainDir[i].Contains("GammaConvV1_QA_")){
                    nameCutsPQA[i] = Form("%s",nameMainDir[i].Data());
                    nameCutsPQA[i].Replace(0,15,"");
                    mainDir[i].push_back(nameMainDir[i]);
                    cutsPQA[i].push_back(nameCutsPQA[i]);
                }
            }
            cout << "Found " << cutsPQA[i].size() << " different cuts. The following cuts are available:" << endl;
            for(Int_t j = 0; j < (Int_t) cutsPQA[i].size(); j++) {
                cout << Form("(%i) -- %s", j, cutsPQA[i].at(j).Data()) << endl;
            }
            if(cutNr == -1){
                do{ cin >> cutNr;}
                while( (cutNr < 0) || (cutNr > (Int_t) cutsPQA[i].size()) );
            }
            if (i == 0 || !useConsistentCut){
              cout << "Processing Cut Number: " << cutNr << endl;
              nameMainDir[i]      = mainDir[i].at(cutNr);
              nameCutsPQA[i]      = cutsPQA[i].at(cutNr);
              nameCutsPQAEvent[i] = Form("%s",nameMainDir[i].Data());
              nameCutsPQAEvent[i].Replace(0,15,"");
              nameCutsPQAEvent[i].Replace(8,39,"");
              nameCutsPQAmain[i] = Form("%s",nameMainDir[i].Data());
              nameCutsPQAmain[i].Replace(0,15,"");
              nameCutsPQAmain[i].Replace(35,12,"");
            } else {
              nameMainDir[i]      = nameMainDir[0];
              nameCutsPQA[i]      = nameCutsPQA[0];
              nameCutsPQAEvent[i] = nameCutsPQAEvent[0];
              nameCutsPQAmain[i]  = nameCutsPQAmain[0];
            }
            cout << endl;
            cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
            cout << "nameMainDir:" << nameMainDir[i].Data() << endl;
            cout << "long cutnumber for PhotonQA: " << nameCutsPQA[i] << "\nevent cutnumber for PhotonQA: "<< nameCutsPQAEvent[i] <<  "\nbase cutnumber for PhotonQA: "<< nameCutsPQAmain[i]<< endl;
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
    TString fTriggerCut = nameCutsPQAEvent[0](3,2);
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
    TString centralityString    = GetCentralityString(nameCutsPQAEvent[0]);
    if (centralityString.CompareTo("pp")!=0 && !centralityString.Contains("0-100%") ){
      fCollisionSystem    = Form("%s %s", centralityString.Data(), fCollisionSystem.Data());
    }

    TString fDetectionProcess   = ReturnFullTextReconstructionProcess(fMode);
    TString labelALICEforPlots  = "ALICE";

    TString outputDir = Form("%s/%s/PhotonQA/%s",nameCutsPQA[0].Data(),fEnergyFlag.Data(),suffix.Data());
    if(addSubfolder) outputDir+=Form("/%s",DataSets[0].Data());
    TString outputDirSp = Form("%s/%s/PhotonQA/%s/Special",nameCutsPQA[0].Data(),fEnergyFlag.Data(),suffix.Data());
    if(addSubfolder) outputDirSp+=Form("/%s",DataSets[0].Data());

    gSystem->Exec("mkdir -p "+outputDir);
    gSystem->Exec("mkdir -p "+outputDir+"/Comparison");
    gSystem->Exec("mkdir -p "+outputDir+"/Comparison/Ratios");
    gSystem->Exec("mkdir -p "+outputDir + "/Special/AlphaQt");
    gSystem->Exec("mkdir -p "+outputDir + "/Special/Chi2PsiPair");


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

    Double_t projBinnin[16] = { 0.0, 0.1, 0.2, 0.4, 0.7,
                                1.0, 1.5, 2.0, 3.0, 4.0,
                                6.0, 10., 15., 20., 30.,
                                50.};

    TCanvas* canvas = new TCanvas("canvas","",10,10,750,500);  // gives the page size
    TCanvas* cvsQuadratic = new TCanvas("cvsQuadratic","",10,10,500,500);  // gives the page size
    Double_t leftMargin = 0.09; Double_t rightMargin = 0.02; Double_t topMargin = 0.04; Double_t bottomMargin = 0.09;
    DrawGammaCanvasSettings(canvas,leftMargin,rightMargin,topMargin,bottomMargin);
    Double_t xPosLabel1D = 0.95;
    Double_t xPosLabel2D = 0.85;
    TLatex *labelColor                  = NULL;
    TLatex *labelBW                     = NULL;
    TLatex *labelEnergy                 = NULL;
    TLatex *labelProcess                = NULL;
    TExec *ex1                          = new TExec("ex1","PalColor();");
    TExec *ex2                          = new TExec("ex2","gStyle->SetPalette(9);");
    TLatex* labelpTrange                = NULL;
    for(Int_t i=0; i<nSets; i++){
        //-----------------------------------
        TString path = pathDataSets[i];
        if(!path.Contains("/PhotonQA")){
            cout << "ERROR: Path to PhotonQA output does contin with '/PhotonQA' for DataSet '" << DataSets[i].Data() << "'" << endl;
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

        TString topdirName = Form("GammaConvV1_QA_%s_%s",nameCutsPQAmain[i].Data(), cutTreeProjection.Data());
        TDirectory* TopDir = (TDirectory*) fPhotonQAFile->Get(topdirName.Data());
        if(TopDir == NULL) {cout << "ERROR: TopDir not Found:"<< topdirName.Data() <<endl; return;}
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
        TH3F* fHistGammaAlphaQtPt                   = (TH3F*)TopDir->Get("histoGammaAlphaQtPt");
        TH3F* fHistGammaAlphaQtPtCopy               = NULL;
        TH2D* fHistGammaAlphaQt                     = NULL;
        TH2D* fHistGammaAlphaPt                     = NULL;
        TH2D* fHistGammaQtPt                        = NULL;
        TH2D* fHistGammaAlphaQtPtSliced[15]         = {NULL};
        TH3F* fHistTrueMCKindAlphaQtPt[20]          = {NULL};
        TH3F* fHistTrueMCKindAlphaQtPtCopy[20]      = {NULL};
        TH2D* fHistTrueMCKindQtPt[20]               = {NULL};
        TH2D* fHistTrueMCKindAlphaPt[20]            = {NULL};
        TH2D* fHistTrueMCKindAlphaQt[20]            = {NULL};
        TH2D* fHistTrueMCKindAlphaQtPtSliced[20][15]= {NULL};
        TH2D* fHistTrueDalitzGammaAlphaPt           = NULL;
        TH2D* fHistTrueDalitzGammaAlphaQt           = NULL;
        TH2D* fHistTrueDalitzGammaQtPt              = NULL;

        textSizeLabelsPixel                   = 1200*0.04;
        TCanvas* canvasAlphaQtPlots = new TCanvas("canvasAlphaQtPlots","",200,10,1350,1200);  // gives the page size
        DrawGammaCanvasSettings( canvasAlphaQtPlots, 0.1, 0.02, 0.02, 0.09);
        // canvasAlphaQtPlots->SetLogy();
        canvasAlphaQtPlots->SetLogz();
        TH2F * histo2DAlphaQtDummy = new TH2F("histo2DAlphaQtDummy","histo2DAlphaQtDummy",200,-1.05,1.05,200,0.,0.15);
        SetStyleHistoTH2ForGraphs(histo2DAlphaQtDummy, "#alpha^{#gamma} = (#it{p}^{+}_{L}-#it{p}^{-}_{L})/(#it{p}^{+}_{L}+#it{p}^{-}_{L})", "#it{q}_{T}^{#gamma}",
                                  0.035,0.04, 0.035,0.04, 0.98,1.2);

        isMC[i] = kFALSE;

        if (fHistGammaAlphaQtPt){
            fHistGammaAlphaQt   = (TH2D*)fHistGammaAlphaQtPt->Project3D("yx");
            nGammas             = fHistGammaAlphaQt->Integral();
            fHistGammaAlphaQt->Sumw2();
            fHistGammaAlphaQt->Scale(1./nGammas);
            fHistGammaQtPt      = (TH2D*)fHistGammaAlphaQtPt->Project3D("zy");
            fHistGammaQtPt->Sumw2();
            fHistGammaQtPt->Scale(1./nGammas);
            fHistGammaAlphaPt      = (TH2D*)fHistGammaAlphaQtPt->Project3D("zx");
            fHistGammaAlphaPt->Sumw2();
            fHistGammaAlphaPt->Scale(1./nGammas);


            fHistGammaAlphaQtPtCopy                 = (TH3F*)fHistGammaAlphaQtPt->Clone("fHistGammaAlphaQtPtCopy");
            for(Int_t j=0; j<15; j++){
                fHistGammaAlphaQtPtCopy->GetZaxis()->SetRangeUser(projBinnin[j],projBinnin[j+1]);
                fHistGammaAlphaQtPtSliced[j]       = (TH2D*)fHistGammaAlphaQtPtCopy->Project3D(Form("yx%d",j));
                fHistGammaAlphaQtPtSliced[j]->Sumw2();
                fHistGammaAlphaQtPtSliced[j]->Scale(1./nGammas);
            }
        }

        for(Int_t k=0;k<20;k++){
            if(k != 0 && k != 3 && k !=5 && k != 4 && k!=10 && k!=11 && k!=13) continue;
            fHistTrueMCKindAlphaQtPt[k]             = (TH3F*)TopDir->Get(Form("histoTrueMCKindAlphaQtPt_kind%d",k));
            if (fHistTrueMCKindAlphaQtPt[k]){
                isMC[i]     = kTRUE;
            } else {
                cout << "Couldn't find: " << Form("histoTrueMCKindAlphaQtPt_kind%d",k) << endl;
                continue;
            }
            fHistTrueMCKindQtPt[k]                  = (TH2D*)fHistTrueMCKindAlphaQtPt[k]->Project3D(Form("zy%d",k));
            CheckNEntries(fHistTrueMCKindQtPt[k]);
            fHistTrueMCKindQtPt[k]->Sumw2();
            fHistTrueMCKindQtPt[k]->Scale(1./nGammas);
            fHistTrueMCKindAlphaPt[k]               = (TH2D*)fHistTrueMCKindAlphaQtPt[k]->Project3D(Form("zx%d",k));
            CheckNEntries(fHistTrueMCKindAlphaPt[k]);
            fHistTrueMCKindAlphaPt[k]->Sumw2();
            fHistTrueMCKindAlphaPt[k]->Scale(1./nGammas);
            fHistTrueMCKindAlphaQt[k]               = (TH2D*)fHistTrueMCKindAlphaQtPt[k]->Project3D(Form("yx%d",k));
            CheckNEntries(fHistTrueMCKindAlphaQt[k]);
            fHistTrueMCKindAlphaQt[k]->Sumw2();
            fHistTrueMCKindAlphaQt[k]->Scale(1./nGammas);

            fHistTrueMCKindAlphaQtPtCopy[k]         = (TH3F*)fHistTrueMCKindAlphaQtPt[k]->Clone(Form("histoTrueMCKindAlphaQtPtCopy%d",k));
            for(Int_t j=0; j<15; j++){
                fHistTrueMCKindAlphaQtPtCopy[k]->GetZaxis()->SetRangeUser(projBinnin[j],projBinnin[j+1]);
                fHistTrueMCKindAlphaQtPtSliced[k][j]       = (TH2D*)fHistTrueMCKindAlphaQtPtCopy[k]->Project3D(Form("yx%d",j));
                fHistTrueMCKindAlphaQtPtSliced[k][j]->Sumw2();
                fHistTrueMCKindAlphaQtPtSliced[k][j]->Scale(1./nGammas);
            }
        }
        if (isMC[i]){
            cout << "adding hists for Dalitz and true gamma" << endl;
            fHistTrueDalitzGammaAlphaQt         = (TH2D*)fHistTrueMCKindAlphaQt[3]->Clone("fHistTrueDalitzGammaAlphaQt");
            fHistTrueDalitzGammaAlphaQt->Sumw2();
            fHistTrueDalitzGammaAlphaQt->Add(fHistTrueMCKindAlphaQt[4]);
            fHistTrueDalitzGammaAlphaPt         = (TH2D*)fHistTrueMCKindAlphaPt[3]->Clone("fHistTrueDalitzGammaAlphaPt");
            fHistTrueDalitzGammaAlphaPt->Sumw2();
            fHistTrueDalitzGammaAlphaPt->Add(fHistTrueMCKindAlphaPt[4]);
            fHistTrueDalitzGammaQtPt         = (TH2D*)fHistTrueMCKindQtPt[3]->Clone("fHistTrueDalitzGammaQtPt");
            fHistTrueDalitzGammaQtPt->Sumw2();
            fHistTrueDalitzGammaQtPt->Add(fHistTrueMCKindQtPt[4]);
        }

        if (isMC[i])  cout << "I found MC" << endl;

        //    ____ _______                   _     ____  _    _
        //   / __ \__   __|           /\    | |   |  _ \| |  | |    /\
        //  | |  | | | |             /  \   | |   | |_| | |__| |   /  \
        //  | |  | | | |            / /\ \  | |   |  __/|  __  |  / /\ \
        //  | |__| | | |           / ____ \ | |___| |   | |  | | / ____ \
        //   \___\_\ |_|          /_/    \_\|_____|_|   |_|  |_|/_/    \_\
        //
        SetPlotStyle();
        if(fHistGammaAlphaQtPt){
            SetZMinMaxTH2(fHistGammaAlphaQt,1,fHistGammaAlphaQt->GetNbinsX(),1,fHistGammaAlphaQt->GetNbinsY(),kTRUE);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.117,0.117,0.014,0.092,kFALSE,kFALSE,kTRUE,
                                fHistGammaAlphaQt,"",
                                "#alpha = (#it{p}^{+}_{L}-#it{p}^{-}_{L})/(#it{p}^{+}_{L}+#it{p}^{-}_{L})","#it{q}_{T} (GeV/#it{c})",1,1.4,
                                xPosLabel2D,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            DrawArmenterosLines();
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistGammaAlphaQt, Form("%s/Armenteros_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));

            TH1D* fHistGammaAlpha = (TH1D*)fHistGammaAlphaQt->ProjectionX("fHistGammaAlpha",1,fHistGammaAlphaQt->GetNbinsY());
            GetMinMaxBin(fHistGammaAlpha,minB,maxB);
            SetXRange(fHistGammaAlpha,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kFALSE,kFALSE,
                                fHistGammaAlpha,"","Photon #alpha","# of Entries",1,1,
                                xPosLabel1D,0.94,0.03,fCollisionSystem,plotDataSets[i],"");
            SaveCanvasAndWriteHistogram(canvas, fHistGammaAlpha, Form("%s/Photon_Alpha_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
            vecGammaAlpha.push_back(fHistGammaAlpha);


            // special plotting
            histo2DAlphaQtDummy->GetZaxis()->SetRangeUser(2./nGammas,fHistGammaAlphaQt->GetMaximum());
            canvasAlphaQtPlots->cd();

            // RECONSTRUCTED PLOT
            histo2DAlphaQtDummy->Draw("copy");
            fHistGammaAlphaQt->Draw("col,same");
            ex1 = new TExec("ex1","PalColor();");
            ex1->Draw();
            fHistGammaAlphaQt->Draw("col,same");
            labelEnergy                         = new TLatex(0.95,0.90,Form("%s, %s",labelALICEforPlots.Data(),fCollisionSystem.Data()));
            SetStyleTLatex( labelEnergy, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,31);
            labelEnergy->Draw();
            labelProcess                        = new TLatex(0.95,0.86,isMC[i] ? "#gamma candidates (MC rec.)" : "#gamma candidates (data)");
            SetStyleTLatex( labelProcess, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,31);
            labelProcess->Draw();
            histo2DAlphaQtDummy->Draw("axis,same");
            canvasAlphaQtPlots->SaveAs(Form("%s/AlphaQt_AllRec_%s.%s",outputDirSp.Data(),DataSets[i].Data(),suffix.Data()));

            for(Int_t bin=0; bin<15; bin++){
                // RECONSTRUCTED PLOT in pT slices
                histo2DAlphaQtDummy->Draw("copy");
                fHistGammaAlphaQtPtSliced[bin]->Draw("col,same");
//                 ex1 = new TExec("ex1","PalColor();");
                ex1->Draw();
                fHistGammaAlphaQtPtSliced[bin]->Draw("col,same");
                labelEnergy->Draw();
                labelProcess->Draw();
                labelpTrange                     = new TLatex(0.95,0.82,Form("%2.1f < #it{p}_{T} < %2.1f GeV/#it{c}",projBinnin[bin],projBinnin[bin+1]));
                SetStyleTLatex( labelpTrange, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,31);
                labelpTrange->Draw();
                histo2DAlphaQtDummy->Draw("axis,same");
                canvasAlphaQtPlots->SaveAs(Form("%s/AlphaQt/%dbin_AlphaQt_vs_Pt_AllRec_%s.%s",outputDirSp.Data(), bin, DataSets[i].Data(),suffix.Data()));
            }
        }else cout << "INFO: Object |fHistGammaAlphaQtPt| could not be found! Skipping Draw..." << endl;
        //-----------------------------------
        if(nGammas <= 1.){
            cout << "********************************************************" << endl;
            cout << "WARNING: nGammas <= 1! Scaling will not work properly..." << endl;
            cout << "********************************************************" << endl;
        }
        if(isMC[i]){
            SetZMinMaxTH2(fHistTrueMCKindAlphaQt[0],1,fHistTrueMCKindAlphaQt[0]->GetNbinsX(),1,fHistTrueMCKindAlphaQt[0]->GetNbinsY(),kTRUE);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                                fHistGammaAlphaQt,"",
                                "#alpha = (#it{p}^{+}_{L}-#it{p}^{-}_{L})/(#it{p}^{+}_{L}+#it{p}^{-}_{L})", "#it{q}_{T} (GeV/#it{c})",1,1.4,
                                xPosLabel2D,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            DrawArmenterosLines();
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistTrueMCKindAlphaQt[0], Form("%s/Armenteros_TruePrimGamma_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));

            //-----------------------------------
            SetZMinMaxTH2(fHistTrueMCKindAlphaQt[5],1,fHistTrueMCKindAlphaQt[5]->GetNbinsX(),1,fHistTrueMCKindAlphaQt[5]->GetNbinsY(),kTRUE);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                                 fHistTrueMCKindAlphaQt[5],"",
                                 "#alpha = (#it{p}^{+}_{L}-#it{p}^{-}_{L})/(#it{p}^{+}_{L}+#it{p}^{-}_{L})", "#it{q}_{T} (GeV/#it{c})",1,1.4,
                                 xPosLabel2D,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            DrawArmenterosLines();
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistTrueMCKindAlphaQt[5], Form("%s/Armenteros_TrueSecGamma_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
            //-----------------------------------

            fHistTrueDalitzGammaAlphaQt         = (TH2D*)fHistTrueMCKindAlphaQt[3]->Clone("fHistTrueDalitzGammaAlphaQt");
            fHistTrueDalitzGammaAlphaQt->Sumw2();
            fHistTrueDalitzGammaAlphaQt->Add(fHistTrueMCKindAlphaQt[4]);

            SetZMinMaxTH2(fHistTrueDalitzGammaAlphaQt,1,fHistTrueDalitzGammaAlphaQt->GetNbinsX(),1,fHistTrueDalitzGammaAlphaQt->GetNbinsY(),kTRUE);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                                 fHistTrueDalitzGammaAlphaQt,"",
                                 "#alpha = (#it{p}^{+}_{L}-#it{p}^{-}_{L})/(#it{p}^{+}_{L}+#it{p}^{-}_{L})", "#it{q}_{T} (GeV/#it{c})",1,1.4,
                                 xPosLabel2D,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            DrawArmenterosLines();
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistTrueDalitzGammaAlphaQt, Form("%s/Armenteros_TrueDalitzGamma_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));

            if(fHistGammaAlphaQt && fHistTrueMCKindAlphaQt[0] && fHistTrueMCKindAlphaQt[5] && fHistTrueDalitzGammaAlphaQt){

                TH2D* fHistTrueBGAlphaQt = (TH2D*) fHistGammaAlphaQt->Clone("histoMonteCarloBGAlphaQtMonteCarlo");
                fHistTrueBGAlphaQt->Sumw2();
                fHistTrueBGAlphaQt->Add(fHistTrueMCKindAlphaQt[5],-1);
                fHistTrueBGAlphaQt->Add(fHistTrueMCKindAlphaQt[0],-1);
                fHistTrueBGAlphaQt->Add(fHistTrueDalitzGammaAlphaQt,-1);
                SetZMinMaxTH2(fHistTrueBGAlphaQt,1,fHistTrueBGAlphaQt->GetNbinsX(),1,fHistTrueBGAlphaQt->GetNbinsY(),kTRUE);
                CheckForNegativeEntries(fHistTrueBGAlphaQt);
                DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                                    fHistTrueBGAlphaQt,"",
                                    "#alpha = (#it{p}^{+}_{L}-#it{p}^{-}_{L})/(#it{p}^{+}_{L}+#it{p}^{-}_{L})", "#it{q}_{T} (GeV/#it{c})",1,1.4,
                                    xPosLabel2D,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
                DrawArmenterosLines();
                SaveCanvasAndWriteHistogram(cvsQuadratic, fHistTrueBGAlphaQt, Form("%s/Armenteros_BG_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
                delete fHistTrueBGAlphaQt;
            }

            // TRUE GAMMAS PLOT
            canvasAlphaQtPlots->cd();
            histo2DAlphaQtDummy->GetZaxis()->SetRangeUser(2./nGammas,fHistGammaAlphaQt->GetMaximum());
            histo2DAlphaQtDummy->Draw("copy");
            // draw true gamma->ee
            fHistTrueMCKindAlphaQt[0]->Draw("col,same");
//             ex1 = new TExec("ex1","PalColor();");
            ex1->Draw();
            fHistTrueMCKindAlphaQt[0]->Draw("col,same");
            labelColor                  = new TLatex(0.12,0.90, "#gamma rec. from e^{+}e^{-} (color)" );
            SetStyleTLatex( labelColor, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,11);
            labelColor->Draw();
            labelProcess                     = new TLatex(0.95,0.86, "#gamma candidates (MC true)");
            SetStyleTLatex( labelProcess, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,31);
            labelEnergy->Draw();
            labelProcess->Draw();
            histo2DAlphaQtDummy->Draw("axis,same");
            canvasAlphaQtPlots->SaveAs(Form("%s/AlphaQt_trueGamma_woCuts_%s.%s",outputDirSp.Data(),DataSets[i].Data(),suffix.Data()));

            labelBW                     = new TLatex(0.12,0.86,"#gamma rec. from e^{#pm}#pi^{#pm} (BAW)");
            SetStyleTLatex( labelBW, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,11);

            for(Int_t bin=0; bin<15; bin++){
                // RECONSTRUCTED PLOT
                histo2DAlphaQtDummy->GetZaxis()->SetRangeUser(2./nGammas,fHistGammaAlphaQt->GetMaximum());
                histo2DAlphaQtDummy->Draw("copy");
                    fHistTrueMCKindAlphaQtPtSliced[0][bin]->Draw("col,same");
                    ex1 = new TExec("ex1","PalColor();");
                    ex1->Draw();
                    fHistTrueMCKindAlphaQtPtSliced[0][bin]->Draw("col,same");
                    labelColor->Draw();
                    labelEnergy->Draw();
                    labelProcess->Draw();
                    labelpTrange    = new TLatex(0.95,0.82,Form("%2.1f < #it{p}_{T} < %2.1f GeV/#it{c}",projBinnin[bin],projBinnin[bin+1]));
                    SetStyleTLatex( labelpTrange, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,31);
                    labelpTrange->Draw();
                histo2DAlphaQtDummy->Draw("axis,same");
                canvasAlphaQtPlots->SaveAs(Form("%s/AlphaQt/%dbin_AlphaQt_vs_Pt_trueGamma_%s.%s",outputDirSp.Data(), bin, DataSets[i].Data() ,suffix.Data()));

                // RECONSTRUCTED PLOT
                histo2DAlphaQtDummy->GetZaxis()->SetRangeUser(1./nGammas,fHistTrueMCKindAlphaQtPtSliced[11][bin]->GetMaximum()*10);
                histo2DAlphaQtDummy->Draw("copy");
                    fHistTrueMCKindAlphaQtPtSliced[0][bin]->Draw("col,same");
                    ex1 = new TExec("ex1","PalColor();");
                    ex1->Draw();
                    fHistTrueMCKindAlphaQtPtSliced[0][bin]->Draw("col,same");
                    fHistTrueMCKindAlphaQtPtSliced[11][bin]->Draw("col,same");
                    ex2 = new TExec("ex2","gStyle->SetPalette(9);");
                    ex2->Draw();
                    fHistTrueMCKindAlphaQtPtSliced[11][bin]->Draw("col,same");
                    labelEnergy->Draw();
                    labelProcess->Draw();
                    labelColor->Draw();
                    labelBW->Draw();
                    labelpTrange->Draw();
                histo2DAlphaQtDummy->Draw("axis,same");
                canvasAlphaQtPlots->SaveAs(Form("%s/AlphaQt/%dbin_AlphaQt_vs_Pt_trueGammaAndComb_%s.%s",outputDirSp.Data(), bin,DataSets[i].Data() ,suffix.Data()));
                ex1->Draw();
            }
        }
        //-----------------------------------


        //    ____ _______      ___  _______
        //   / __ \__   __|    |  _ \__   __|
        //  | |  | | | |       | |_| | | |
        //  | |  | | | |       |  __/  | |
        //  | |__| | | |       | |     | |
        //   \___\_\ |_|       |_|     |_|

        TCanvas* canvasQtPlots = new TCanvas("canvasQtPlots","",200,10,1350,1200);  // gives the page size
        DrawGammaCanvasSettings( canvasQtPlots, 0.08, 0.02, 0.02, 0.09);
        canvasQtPlots->SetLogy();
        canvasQtPlots->SetLogz();

        TH2F * histo2DQtDummy = new TH2F("histo2DQtDummy","histo2DQtDummy",100,0,0.1,1000,0.04,40);
        SetStyleHistoTH2ForGraphs(histo2DQtDummy, "#it{q}_{T}^{#gamma}","#it{p}_{T} (GeV/#it{c})",0.035,0.04, 0.035,0.04, 0.98,0.9);
        // initialize cuts for QT vs pT
        Double_t funcParamQtPt[4][2]    = {{0.11,0.04},{0.125,0.05},{0.14,0.06},{0.16,0.07}};
        if (fEnergyFlag.CompareTo("XeXe_5.44TeV") == 0 || fEnergyFlag.CompareTo("13TeVLowB") == 0){
            funcParamQtPt[0][0]     = 0.18;
            funcParamQtPt[0][1]     = 0.032;
            funcParamQtPt[1][0]     = 0.2;
            funcParamQtPt[1][1]     = 0.035;
            funcParamQtPt[2][0]     = 0.25;
            funcParamQtPt[2][1]     = 0.04;
            funcParamQtPt[3][0]     = 0.3;
            funcParamQtPt[3][1]     = 0.045;
        }
        Color_t colorQtFunc[5]      = {kMagenta+2, kMagenta+4, kMagenta-4, kMagenta-8, kMagenta+1};
        Style_t styleQtFunc[5]      = {1,7,8,9,6};
        TF1* funcPtDepQtCut[4]      = {NULL};
        TF1 *funcOldCutDummy = new TF1("funcOldCutDummy","x/[0]",0.005,1.);
        DrawGammaSetMarkerTF1( funcOldCutDummy, 9, 3, kGreen-8);
        TLegend* legendQtPlotFits  = GetAndSetLegend2(0.58, 0.13, 0.95, 0.13+(0.035*5*1.15), 0.75*textSizeLabelsPixel,1,"",43,0.1);
        legendQtPlotFits->AddEntry(funcOldCutDummy, "#it{q}_{T}^{max} = 0.05","l");
        for (Int_t k = 0; k< 4; k++){
            legendQtPlotFits->AddEntry(funcPtDepQtCut[k], Form("#it{q}_{T}^{max} = %0.3f#it{p}_{T}, #it{q}_{T}^{max} = %0.3f",funcParamQtPt[k][0],funcParamQtPt[k][1] ),"l");
        }
        SetPlotStyle();

        if(fHistGammaAlphaQtPt){
            GetMinMaxBin(fHistGammaQtPt,minB,maxB);
            GetMinMaxBinY(fHistGammaQtPt,minYB,maxYB); maxYB+=5;
            SetXRange(fHistGammaQtPt,minB,maxB);
            SetYRange(fHistGammaQtPt,fHistGammaQtPt->GetYaxis()->FindBin(0.1),maxYB);
            SetZMinMaxTH2(fHistGammaQtPt,minB,maxB,fHistGammaQtPt->GetYaxis()->FindBin(0.1),maxYB,kTRUE);
            //-----------------------------------
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kTRUE,kTRUE,
                                 fHistGammaQtPt,"",
                                 "#it{q}_{T,#gamma}", "#it{p}_{T,#gamma} (GeV/#it{c})",1,1.4,
                                 xPosLabel2D,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistGammaQtPt, Form("%s/Qt_Pt_Photon_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));

            histo2DQtDummy->GetZaxis()->SetRangeUser(6/nGammas,fHistGammaQtPt->GetMaximum());

            //-----------------------------------
            canvasQtPlots->cd();
            histo2DQtDummy->Draw("copy");

            fHistGammaQtPt->Draw("col,same");
            fHistGammaQtPt->Draw("col,same");
            DrawGammaLines(0.05, 0.05, 0.04, 8, 3, kGreen-8, 9);
            for (Int_t k = 0; k< 4; k++){
                funcPtDepQtCut[k] = new TF1(Form("funcPtDepQtCut_%d",k),"x/[0]",0.005,funcParamQtPt[k][1]);
                DrawGammaSetMarkerTF1( funcPtDepQtCut[k], styleQtFunc[k+1], 2, colorQtFunc[k+1]);
                funcPtDepQtCut[k]->SetParameter(0,funcParamQtPt[k][0]);
                funcPtDepQtCut[k]->Draw("same");
                DrawGammaLines(funcParamQtPt[k][1], funcParamQtPt[k][1], funcParamQtPt[k][1]/funcParamQtPt[k][0], 8, 2, colorQtFunc[k+1], styleQtFunc[k+1]);
            }

            labelEnergy                  = new TLatex(0.95,0.90,Form("%s, %s",labelALICEforPlots.Data(),fCollisionSystem.Data()));
            SetStyleTLatex( labelEnergy, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,31);
            labelEnergy->Draw();
            labelProcess                     = new TLatex(0.95,0.86,isMC[i] ? "#gamma candidates (MC rec.)" : "#gamma candidates (data)");
            SetStyleTLatex( labelProcess, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,31);
            labelProcess->Draw();
            histo2DQtDummy->Draw("axis,same");
            canvasQtPlots->SaveAs(Form("%s/QtPt_withCuts_%s.%s",outputDirSp.Data(),DataSets[i].Data(),suffix.Data()));

        }else cout << "INFO: Object |fHistGammaQtPt| could not be found! Skipping Draw..." << endl;
        //-----------------------------------
        if(isMC[i]){
            SetXRange(fHistTrueMCKindQtPt[0],minB,maxB);
            SetYRange(fHistTrueMCKindQtPt[0],fHistTrueMCKindQtPt[0]->GetYaxis()->FindBin(0.1),maxYB);
            SetZMinMaxTH2(fHistTrueMCKindQtPt[0],minB,maxB,fHistTrueMCKindQtPt[0]->GetYaxis()->FindBin(0.1),maxYB,kTRUE);
            //-----------------------------------
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kTRUE,kTRUE,
                                 fHistTrueMCKindQtPt[0],"",
                                 "#it{q}_{T,#gamma}", "#it{p}_{T,#gamma} (GeV/#it{c})",1,1.4,
                                 xPosLabel2D,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistTrueMCKindQtPt[0], Form("%s/Qt_Pt_TruePrimGamma_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));

            //-----------------------------------
            SetXRange(fHistTrueMCKindQtPt[5],minB,maxB);
            SetYRange(fHistTrueMCKindQtPt[5],fHistTrueMCKindQtPt[5]->GetYaxis()->FindBin(0.1),maxYB);
            SetZMinMaxTH2(fHistTrueMCKindQtPt[5],minB,maxB,fHistTrueMCKindQtPt[5]->GetYaxis()->FindBin(0.1),maxYB,kTRUE);
            //-----------------------------------
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kTRUE,kTRUE,
                                 fHistTrueMCKindQtPt[5],"",
                                 "#it{q}_{T,#gamma}", "#it{p}_{T,#gamma} (GeV/#it{c})",1,1.4,
                                 xPosLabel2D,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistTrueMCKindQtPt[5], Form("%s/Qt_Pt_TrueSecGamma_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
            //-----------------------------------
            SetXRange(fHistTrueDalitzGammaQtPt,minB,maxB);
            SetYRange(fHistTrueDalitzGammaQtPt,fHistTrueDalitzGammaQtPt->GetYaxis()->FindBin(0.1),maxYB);
            SetZMinMaxTH2(fHistTrueDalitzGammaQtPt,minB,maxB,fHistTrueDalitzGammaQtPt->GetYaxis()->FindBin(0.1),maxYB,kTRUE);
            //-----------------------------------
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kTRUE,kTRUE,
                                 fHistTrueDalitzGammaQtPt,"",
                                 "#it{q}_{T,#gamma}", "#it{p}_{T,#gamma} (GeV/#it{c})",1,1.4,
                                 xPosLabel2D,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistTrueDalitzGammaQtPt, Form("%s/Qt_Pt_TrueDalitz_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));

            //-----------------------------------
            if(fHistGammaQtPt && fHistTrueMCKindQtPt[0] && fHistTrueMCKindQtPt[5] && fHistTrueDalitzGammaQtPt){
                ex1->Draw();
                TH2D* fHistBGQtPt = (TH2D*) fHistGammaQtPt->Clone("histoBGQtMonteCarloPt");
                fHistBGQtPt->Sumw2();
                fHistBGQtPt->Add(fHistTrueMCKindQtPt[0],-1.);
                fHistBGQtPt->Add(fHistTrueMCKindQtPt[5],-1.);
                fHistBGQtPt->Add(fHistTrueDalitzGammaQtPt,-1.);
                SetXRange(fHistBGQtPt,minB,maxB);
                SetYRange(fHistBGQtPt,fHistBGQtPt->GetYaxis()->FindBin(0.1),maxYB);
                SetZMinMaxTH2(fHistBGQtPt,minB,maxB,fHistBGQtPt->GetYaxis()->FindBin(0.1),maxYB,kTRUE);
                CheckForNegativeEntries(fHistBGQtPt);
                //-----------------------------------
                DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kTRUE,kTRUE,
                                     fHistBGQtPt,"",
                                     "#it{q}_{T,#gamma}", "#it{p}_{T,#gamma} (GeV/#it{c})",1,1.4,
                                     xPosLabel2D,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
                SaveCanvasAndWriteHistogram(cvsQuadratic, fHistBGQtPt, Form("%s/Qt_Pt_BG_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
                delete fHistBGQtPt;
            }

            histo2DQtDummy->GetZaxis()->SetRangeUser(6./nGammas,fHistTrueMCKindQtPt[11]->GetMaximum()*1000);
            canvasQtPlots->cd();
            histo2DQtDummy->Draw("copy");
                // draw true gamma->ee
                fHistTrueMCKindQtPt[0]->Draw("col,same");
                ex1 = new TExec("ex1","PalColor();");
                ex1->Draw();
                fHistTrueMCKindQtPt[0]->Draw("col,same");
                // draw ee combinatorics
                fHistTrueMCKindQtPt[11]->Draw("col,same");
                ex2 = new TExec("ex2","gStyle->SetPalette(9);");
                ex2->Draw();
                fHistTrueMCKindQtPt[11]->Draw("col,same");
                // draw pipi combinatorics
                fHistTrueMCKindQtPt[13]->Draw("col,same");
                ex2 = new TExec("ex2","gStyle->SetPalette(9);");
                ex2->Draw();
                fHistTrueMCKindQtPt[13]->Draw("col,same");
                DrawGammaLines(0.05, 0.05, 0.04, 8, 3, kGreen-8, 9);
                for (Int_t k = 0; k< 4; k++){
                    funcPtDepQtCut[k]->Draw("same");
                    DrawGammaLines(funcParamQtPt[k][1], funcParamQtPt[k][1], funcParamQtPt[k][1]/funcParamQtPt[k][0], 8, 2, colorQtFunc[k+1], styleQtFunc[k+1]);
                }

                labelColor                  = new TLatex(0.12,0.90, "#gamma rec. from e^{+}e^{-} (color)" );
                SetStyleTLatex( labelColor, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,11);
                labelColor->Draw();
                labelBW                     = new TLatex(0.12,0.86, "#gamma rec. from e^{#pm}#pi^{#pm}/#pi^{+}#pi^{-} (BAW)" );
                SetStyleTLatex( labelBW, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,11);
                labelBW->Draw();
                labelEnergy->Draw();
                labelProcess->Draw();
            histo2DQtDummy->Draw("axis,same");
            canvasQtPlots->SaveAs(Form("%s/QtPt_MCSep_withCuts_%s.%s",outputDirSp.Data(),DataSets[i].Data(),suffix.Data()));

            histo2DQtDummy->GetZaxis()->SetRangeUser(6/nGammas,fHistGammaQtPt->GetMaximum());
            histo2DQtDummy->Draw("copy");
                fHistTrueMCKindQtPt[0]->Draw("col,same");
                ex1 = new TExec("ex1","PalColor();");
                ex1->Draw();
                fHistTrueMCKindQtPt[0]->Draw("col,same");
                DrawGammaLines(0.05, 0.05, 0.04, 8, 3, kGreen-8, 9);
                for (Int_t k = 0; k< 4; k++){
                    funcPtDepQtCut[k]->Draw("same");
                    DrawGammaLines(funcParamQtPt[k][1], funcParamQtPt[k][1], funcParamQtPt[k][1]/funcParamQtPt[k][0], 8, 2, colorQtFunc[k+1], styleQtFunc[k+1]);
                }
                legendQtPlotFits->Draw();
                labelColor->Draw();
                labelEnergy->Draw();
            labelProcess->Draw();
            histo2DQtDummy->Draw("axis,same");
            canvasQtPlots->SaveAs(Form("%s/QtPt_trueGamma_wCuts_%s.%s",outputDirSp.Data(),DataSets[i].Data(),suffix.Data()));

            histo2DQtDummy->GetZaxis()->SetRangeUser(6./nGammas,fHistTrueMCKindQtPt[11]->GetMaximum()*1000);
            histo2DQtDummy->Draw("copy");
                fHistTrueMCKindQtPt[0]->Draw("col,same");
                ex1 = new TExec("ex1","PalColor();");
                ex1->Draw();
                fHistTrueMCKindQtPt[0]->Draw("col,same");
                fHistTrueMCKindQtPt[11]->Draw("col,same");
                ex2 = new TExec("ex2","gStyle->SetPalette(9);");
                ex2->Draw();
                fHistTrueMCKindQtPt[11]->Draw("col,same");

                labelBW                     = new TLatex(0.12,0.86,"#gamma rec. from #pi^{+}#pi^{-} (BAW)");
                SetStyleTLatex( labelBW, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,11);
                labelBW->Draw();
                labelColor->Draw();
                labelEnergy->Draw();
                labelProcess->Draw();
            histo2DQtDummy->Draw("axis,same");
            canvasQtPlots->SaveAs(Form("%s/QtPt_trueGammaAndComb_woCuts_%s.%s",outputDirSp.Data(),DataSets[i].Data(),suffix.Data()));

            histo2DQtDummy->Draw("copy");
                fHistTrueMCKindQtPt[0]->Draw("col,same");
                ex1 = new TExec("ex1","PalColor();");
                ex1->Draw();
                fHistTrueMCKindQtPt[0]->Draw("col,same");
                fHistTrueMCKindQtPt[11]->Draw("col,same");
                ex2 = new TExec("ex2","gStyle->SetPalette(9);");
                ex2->Draw();
                fHistTrueMCKindQtPt[11]->Draw("col,same");
                fHistTrueMCKindQtPt[13]->Draw("col,same");
                ex2 = new TExec("ex2","gStyle->SetPalette(9);");
                ex2->Draw();
                fHistTrueMCKindQtPt[13]->Draw("col,same");
                labelBW                     = new TLatex(0.12,0.86,"#gamma rec. from e^{#pm}#pi^{#pm}/#pi^{+}#pi^{-} (BAW)" );
                SetStyleTLatex( labelBW, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,11);
                labelBW->Draw();
                labelColor->Draw();
                labelEnergy->Draw();
                labelProcess->Draw();
            histo2DQtDummy->Draw("axis,same");
            canvasQtPlots->SaveAs(Form("%s/QtPt_trueGammaAndCombPlus_woCuts_%s.%s",outputDirSp.Data(),DataSets[i].Data(),suffix.Data()));
        }

        //             _     ____  _    _                    ___  _______
        //      /\    | |   |  _ \| |  | |    /\            |  _ \__   __|
        //     /  \   | |   | |_| | |__| |   /  \           | |_| | | |
        //    / /\ \  | |   |  __/|  __  |  / /\ \          |  __/  | |
        //   / ____ \ | |___| |   | |  | | / ____ \         | |     | |
        //  /_/    \_\|_____|_|   |_|  |_|/_/    \_\        |_|     |_|
        //

        TCanvas* canvasAsymPlots = new TCanvas("canvasAsymPlots","",200,10,1350,1200);  // gives the page size
        DrawGammaCanvasSettings( canvasAsymPlots, 0.08, 0.02, 0.02, 0.09);
        canvasAsymPlots->SetLogy();
        canvasAsymPlots->SetLogz();
        TH2F * histo2DAlphaDummy = new TH2F("histo2DAlphaDummy","histo2DAlphaDummy",100,-1.07,1.07,1000,0.04,200);
        SetStyleHistoTH2ForGraphs(histo2DAlphaDummy, "#alpha^{#gamma} = (#it{p}^{+}_{L}-#it{p}^{-}_{L})/(#it{p}^{+}_{L}+#it{p}^{-}_{L})","#it{p}_{T} (GeV/#it{c})",0.035,0.04, 0.035,0.04, 0.98,0.9);
        SetPlotStyle();
        if(fHistGammaAlphaQtPt){
            histo2DAlphaDummy->GetZaxis()->SetRangeUser(1/nGammas,fHistGammaAlphaPt->GetMaximum());
            canvasAsymPlots->cd();

                // RECONSTRUCTED PLOT
                histo2DAlphaDummy->Draw("copy");
                fHistGammaAlphaPt->Draw("col,same");
                labelEnergy                  = new TLatex(0.95,0.90,Form("%s, %s",labelALICEforPlots.Data(),fCollisionSystem.Data()));
                SetStyleTLatex( labelEnergy, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,31);
                labelEnergy->Draw();
                labelProcess                     = new TLatex(0.95,0.86,isMC[i] ? "#gamma candidates (MC rec.)" : "#gamma candidates (data)");
                SetStyleTLatex( labelProcess, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,31);
                labelProcess->Draw();
            histo2DAlphaDummy->Draw("axis,same");
            canvasAsymPlots->SaveAs(Form("%s/AlphaPt_AllRec_%s.%s",outputDirSp.Data(),DataSets[i].Data(),suffix.Data()));

                DrawGammaLines(-0.95, -0.95, 0.04, 50, 3, kGreen+3, 9);
                DrawGammaLines( 0.95,  0.95, 0.04, 50, 3, kGreen+3, 9);
                funcOldCutDummy = new TF1("funcOldCutDummy","x/[0]",0.005,1.);
                DrawGammaSetMarkerTF1( funcOldCutDummy, 9, 3, kGreen+3);
                TLegend* legendAlphaPlotFits  = GetAndSetLegend2(0.67, 0.13, 0.95, 0.13+(0.035*1*1.35), 0.85*textSizeLabelsPixel);
                legendAlphaPlotFits->AddEntry(funcOldCutDummy, "|#alpha^{max}| = 0.95","l");
                legendAlphaPlotFits->Draw();
                TLatex* labelHighpTSignal    = new TLatex(0.90,0.78,"high #it{p}_{T} signal #rightarrow");
                SetStyleTLatex( labelHighpTSignal, 0.95*textSizeLabelsPixel,4,kBlue+2,43,kTRUE,31);
                labelHighpTSignal->Draw();
                labelEnergy->Draw();
                labelProcess->Draw();
            histo2DAlphaDummy->Draw("axis,same");
            canvasAsymPlots->SaveAs(Form("%s/AlphaPt_Final_withCuts_%s.%s",outputDirSp.Data() ,DataSets[i].Data(),suffix.Data()));
        }

        if (isMC[i]){
            // TRUE GAMMAS PLOT
            canvasAsymPlots->cd();
            histo2DAlphaDummy->Draw("copy");
                fHistTrueMCKindAlphaPt[0]->Draw("col,same");
                ex1 = new TExec("ex1","PalColor();");
                ex1->Draw();
                fHistTrueMCKindAlphaPt[0]->Draw("col,same");
                labelColor                  = new TLatex(0.12,0.90, "#gamma rec. from e^{+}e^{-} (color)" );
                SetStyleTLatex( labelColor, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,11);
                labelColor->Draw();
                labelProcess                     = new TLatex(0.95,0.86,"#gamma candidates (MC true)");
                SetStyleTLatex( labelProcess, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,31);
                labelEnergy->Draw();
                labelProcess->Draw();
            histo2DAlphaDummy->Draw("axis,same");
            canvasAsymPlots->SaveAs(Form("%s/AlphaPt_trueGamma_woCuts_%s.%s",outputDirSp.Data(),DataSets[i].Data(),suffix.Data()));


            // TRUE GAMMAS AND COMBINATORIAL GAMMA CANDIDATES PLOT
            histo2DAlphaDummy->GetZaxis()->SetRangeUser(1/nGammas,fHistTrueMCKindAlphaPt[13]->GetMaximum()*1000);
            histo2DAlphaDummy->Draw("copy");
                // draw true gamma->ee
                fHistTrueMCKindAlphaPt[0]->Draw("col,same");
                ex1 = new TExec("ex1","PalColor();");
                ex1->Draw();
                fHistTrueMCKindAlphaPt[0]->Draw("col,same");
                // draw ee combinatorics
                fHistTrueMCKindAlphaPt[11]->Draw("col,same");
                ex2 = new TExec("ex2","gStyle->SetPalette(9);");
                ex2->Draw();
                fHistTrueMCKindAlphaPt[11]->Draw("col,same");
                // draw pipi combinatorics
                fHistTrueMCKindAlphaPt[13]->Draw("col,same");
                ex2 = new TExec("ex2","gStyle->SetPalette(9);");
                ex2->Draw();
                fHistTrueMCKindAlphaPt[13]->Draw("col,same");
                labelColor->Draw();
                labelBW                     = new TLatex(0.12,0.86,"#gamma rec. from e^{#pm}#pi^{#pm}/#pi^{+}#pi^{-} (BAW)" );
                SetStyleTLatex( labelBW, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,11);
                labelBW->Draw();
                labelEnergy->Draw();
                labelProcess->Draw();
            histo2DAlphaDummy->Draw("axis,same");
            canvasAsymPlots->SaveAs(Form("%s/AlphaPt_trueGammaAndCombPlus_woCuts_%s.%s",outputDirSp.Data(),DataSets[i].Data(),suffix.Data()));

            // TRUE GAMMAS AND COMBINATORIAL GAMMA CANDIDATES PLOT
            histo2DAlphaDummy->Draw("copy");
                // draw true gamma->ee
                fHistTrueMCKindAlphaPt[0]->Draw("col,same");
                ex1 = new TExec("ex1","PalColor();");
                ex1->Draw();
                fHistTrueMCKindAlphaPt[0]->Draw("col,same");
                // draw ee combinatorics
                fHistTrueMCKindAlphaPt[11]->Draw("col,same");
                ex2 = new TExec("ex2","gStyle->SetPalette(9);");
                ex2->Draw();
                fHistTrueMCKindAlphaPt[11]->Draw("col,same");
                // draw pipi combinatorics
                fHistTrueMCKindAlphaPt[13]->Draw("col,same");
                ex2 = new TExec("ex2","gStyle->SetPalette(9);");
                ex2->Draw();
                fHistTrueMCKindAlphaPt[13]->Draw("col,same");
                DrawGammaLines(-0.95, -0.95, 0.04, 50, 3, kGreen+3, 9);
                DrawGammaLines( 0.95,  0.95, 0.04, 50, 3, kGreen+3, 9);
                funcOldCutDummy = new TF1("funcOldCutDummy","x/[0]",0.005,1.);
                DrawGammaSetMarkerTF1( funcOldCutDummy, 9, 3, kGreen+3);
                TLegend* legendAlphaPlotFits  = GetAndSetLegend2(0.67, 0.13, 0.95, 0.13+(0.035*1*1.35), 0.85*textSizeLabelsPixel);
                legendAlphaPlotFits->AddEntry(funcOldCutDummy, "|#alpha^{max}| = 0.95","l");
                legendAlphaPlotFits->Draw();
                TLatex* labelHighpTSignal    = new TLatex(0.90,0.78,"high #it{p}_{T} signal #rightarrow");
                SetStyleTLatex( labelHighpTSignal, 0.95*textSizeLabelsPixel,4,kBlue+2,43,kTRUE,31);
                labelHighpTSignal->Draw();
                labelColor->Draw();
                labelBW->Draw();
                labelEnergy->Draw();
                labelProcess->Draw();
            histo2DAlphaDummy->Draw("axis,same");
            canvasAsymPlots->SaveAs(Form("%s/AlphaPt_MCSep_withCuts_%s.%s",outputDirSp.Data(),DataSets[i].Data(),suffix.Data()));
        }


        //   _____   _____ _____   _____        _____ _____         _____ _    _ _____ ___        ____ _______
        //  |  __ \ / ____|_   _| |  __ \ /\   |_   _|  __ \       / ____| |  | |_   _|__ \      |  _ \ _   __|
        //  | |__) | (___   | |   | |__) /  \    | | | |__) | ___ | |    | |__| | | |    ) | ___ | |_| | | |
        //  |  ___/ \___ \  | |   |  ___/ /\ \   | | |  _  /  ___ | |    |  __  | | |   / /  ___ |  __/  | |
        //  | |     ____) |_| |_  | |  / ____ \ _| |_| | \ \      | |____| |  | |_| |_ / /_      | |     | |
        //  |_|    |_____/|_____| |_| /_/    \_\_____|_|  \_\      \_____|_|  |_|_____|____|     |_|     |_|

        TH3F* fHistGammaChi2PsiPairPt               = (TH3F*)TopDir->Get("histoGammaChi2PsiPairPt");
        TH3F* fHistGammaChi2PsiPairPtCopy           = NULL;
        TH2D* fHistGammaPsiPairPt                   = NULL;
        TH2D* fHistGammaChi2PsiPair                 = NULL;
        TH2D* fHistGammaChi2PsiPairPtSliced[15]     = {NULL};
        TH2D* fHistGammaChi2Pt                      = NULL;

        TH3F* fHistTrueMCKindChi2PsiPairPt[20]      = {NULL};
        TH3F* fHistTrueMCKindChi2PsiPairPtCopy[20]  = {NULL};
        TH2D* fHistTrueMCKindPsiPairPt[20]          = {NULL};
        TH2D* fHistTrueMCKindChi2Pt[20]             = {NULL};
        TH2D* fHistTrueMCKindChi2PsiPair[20]        = {NULL};
        TH2D* fHistTrueMCKindChi2PsiPairPtSliced[20][15]    = {NULL};
        TH2D* histoTrueGammaChi2PsiPair             = NULL;
        TH2D* fHistTrueDalitzGammaPsiPairPt         = NULL;
        TH2D* fHistTrueDalitzGammaChi2PsiPair       = NULL;
        TH2D* fHistTrueDalitzGammaChi2Pt            = NULL;

        if (fHistGammaChi2PsiPairPt){
            fHistGammaPsiPairPt         = (TH2D*)fHistGammaChi2PsiPairPt->Project3D("zy");
            fHistGammaPsiPairPt->Sumw2();
            fHistGammaPsiPairPt->Scale(1./nGammas);
            fHistGammaChi2Pt            = (TH2D*)fHistGammaChi2PsiPairPt->Project3D("zx");
            fHistGammaChi2Pt->Sumw2();
            fHistGammaChi2Pt->Scale(1./nGammas);
            fHistGammaChi2PsiPair       = (TH2D*)fHistGammaChi2PsiPairPt->Project3D("yx");
            fHistGammaChi2PsiPair->Sumw2();
            fHistGammaChi2PsiPair->Scale(1./nGammas);

            fHistGammaChi2PsiPairPtCopy             = (TH3F*)fHistGammaChi2PsiPairPt->Clone("histoGammaChi2PsiPairPtCopy");
            for(Int_t j=0; j<15; j++){
                fHistGammaChi2PsiPairPtCopy->GetZaxis()->SetRangeUser(projBinnin[j],projBinnin[j+1]);
                fHistGammaChi2PsiPairPtSliced[j]   = (TH2D*)fHistGammaChi2PsiPairPtCopy->Project3D(Form("yx%d",j));
                fHistGammaChi2PsiPairPtSliced[j]->Sumw2();
                fHistGammaChi2PsiPairPtSliced[j]->Scale(1./nGammas);
            }

        }

        for(Int_t k=0;k<20;k++){
            if(k != 0 && k != 3 && k !=5 && k != 4 && k!=10 && k!=11 && k!=13) continue;
            fHistTrueMCKindChi2PsiPairPt[k] = (TH3F*)TopDir->Get(Form("histoTrueMCKindChi2PsiPairPt_kind%d",k));
            if (!fHistTrueMCKindChi2PsiPairPt[k]){
                cout << "Couldn't find: " << Form("histoTrueMCKindChi2PsiPairPt_kind%d",k) << endl;
                continue;
            }

            fHistTrueMCKindPsiPairPt[k]     = (TH2D*)fHistTrueMCKindChi2PsiPairPt[k]->Project3D(Form("zy%d",k));
            CheckNEntries(fHistTrueMCKindPsiPairPt[k]);
            fHistTrueMCKindPsiPairPt[k]->Sumw2();
            fHistTrueMCKindPsiPairPt[k]->Scale(1./nGammas);
            fHistTrueMCKindChi2Pt[k]        = (TH2D*)fHistTrueMCKindChi2PsiPairPt[k]->Project3D(Form("zx%d",k));
            CheckNEntries(fHistTrueMCKindChi2Pt[k]);
            fHistTrueMCKindChi2Pt[k]->Sumw2();
            fHistTrueMCKindChi2Pt[k]->Scale(1./nGammas);
            fHistTrueMCKindChi2PsiPair[k]   = (TH2D*)fHistTrueMCKindChi2PsiPairPt[k]->Project3D(Form("yx%d",k));
            CheckNEntries(fHistTrueMCKindChi2PsiPair[k]);
            fHistTrueMCKindChi2PsiPair[k]->Sumw2();
            fHistTrueMCKindChi2PsiPair[k]->Scale(1./nGammas);

            fHistTrueMCKindChi2PsiPairPtCopy[k] = (TH3F*)fHistTrueMCKindChi2PsiPairPt[k]->Clone(Form("histoTrueMCKindChi2PsiPairPtCopy%d",k));
            for(Int_t j=0; j<15; j++){
                fHistTrueMCKindChi2PsiPairPtCopy[k]->GetZaxis()->SetRangeUser(projBinnin[j],projBinnin[j+1]);
                fHistTrueMCKindChi2PsiPairPtSliced[k][j]   = (TH2D*)fHistTrueMCKindChi2PsiPairPtCopy[k]->Project3D(Form("yx%d",j));
                fHistTrueMCKindChi2PsiPairPtSliced[k][j]->Sumw2();
                fHistTrueMCKindChi2PsiPairPtSliced[k][j]->Scale(1./nGammas);
            }
        }

        if (isMC[i]){
            fHistTrueDalitzGammaPsiPairPt      = (TH2D*)fHistTrueMCKindPsiPairPt[3]->Clone("fHistTrueDalitzGammaPsiPairPt");
            fHistTrueDalitzGammaPsiPairPt->Sumw2();
            fHistTrueDalitzGammaPsiPairPt->Add(fHistTrueMCKindPsiPairPt[4]);
            fHistTrueDalitzGammaChi2PsiPair    = (TH2D*)fHistTrueMCKindChi2PsiPair[3]->Clone("fHistTrueDalitzGammaChi2PsiPair");
            fHistTrueDalitzGammaChi2PsiPair->Sumw2();
            fHistTrueDalitzGammaChi2PsiPair->Add(fHistTrueMCKindChi2PsiPair[4]);
            fHistTrueDalitzGammaChi2Pt         = (TH2D*)fHistTrueMCKindChi2Pt[3]->Clone("fHistTrueDalitzGammaChi2Pt");
            fHistTrueDalitzGammaChi2Pt->Sumw2();
            fHistTrueDalitzGammaChi2Pt->Add(fHistTrueMCKindChi2Pt[4]);

            histoTrueGammaChi2PsiPair           = (TH2D*)fHistTrueMCKindChi2PsiPair[0]->Clone("fHistTrueGammaChi2PsiPair");
            histoTrueGammaChi2PsiPair->Sumw2();
            histoTrueGammaChi2PsiPair->Add(fHistTrueMCKindChi2PsiPair[5]);
        }

        //   _____   _____ _____   _____        _____ _____       ____ _______
        //  |  __ \ / ____|_   _| |  __ \ /\   |_   _|  __ \     |  _ \ _   __|
        //  | |__) | (___   | |   | |__) /  \    | | | |__) |    | |_| | | |
        //  |  ___/ \___ \  | |   |  ___/ /\ \   | | |  _  /     |  __/  | |
        //  | |     ____) |_| |_  | |  / ____ \ _| |_| | \ \     | |     | |
        //  |_|    |_____/|_____| |_| /_/    \_\_____|_|  \_\    |_|     |_|

        TCanvas* canvasPsiPairPlots = new TCanvas("canvasPsiPairPlots","",200,10,1350,1200);  // gives the page size
        DrawGammaCanvasSettings( canvasPsiPairPlots, 0.08, 0.02, 0.02, 0.09);
        canvasPsiPairPlots->SetLogy();
        canvasPsiPairPlots->SetLogz();
        TH2F * histo2DPsiPairDummy = new TH2F("histo2DPsiPairDummy","histo2DPsiPairDummy",200,-0.16,0.16,1000,0.04,200);
        SetStyleHistoTH2ForGraphs(histo2DPsiPairDummy, "#Psi_{pair}","#it{p}_{T} (GeV/#it{c})",0.035,0.04, 0.035,0.04, 0.98,0.9);
        TF1* funcOldCutPsiPairDummy = new TF1("funcOldCutDummy","x/[0]",0.005,1.);
        DrawGammaSetMarkerTF1( funcOldCutPsiPairDummy, 9, 3, kMagenta+2);
        TF1* funcOldCutPsiPairDummy1 = new TF1("funcOldCutDummy","x/[0]",0.005,1.);
        DrawGammaSetMarkerTF1( funcOldCutPsiPairDummy1, 7, 3, kMagenta+1);
        TLegend* legendPsiPairPtPlotFits  = GetAndSetLegend2(0.74, 0.75, 0.95, 0.75+(0.035*2*1.35), 0.85*textSizeLabelsPixel);
        legendPsiPairPtPlotFits->AddEntry(funcOldCutPsiPairDummy, "|#Psi_{pair}| = 0.1","l");
        legendPsiPairPtPlotFits->AddEntry(funcOldCutPsiPairDummy1, "|#Psi_{pair}| = 0.15","l");
        SetPlotStyle();

        if(fHistGammaChi2PsiPairPt){
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
                                 xPosLabel2D,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistGammaPsiPairPt, Form("%s/PsiPair_Pt_Photon_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));

            histo2DPsiPairDummy->GetZaxis()->SetRangeUser(1/nGammas,fHistGammaPsiPairPt->GetMaximum());
            canvasPsiPairPlots->cd();

            // RECONSTRUCTED PLOT
            histo2DPsiPairDummy->Draw("copy");
            fHistGammaPsiPairPt->Draw("col,same");
                fHistGammaPsiPairPt->Draw("col,same");
                labelEnergy                  = new TLatex(0.95,0.90,Form("%s, %s",labelALICEforPlots.Data(),fCollisionSystem.Data()));
                SetStyleTLatex( labelEnergy, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,31);
                labelEnergy->Draw();
                labelProcess                     = new TLatex(0.95,0.86,isMC[i] ? "#gamma candidates (MC rec.)" : "#gamma candidates (data)");
                SetStyleTLatex( labelProcess, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,31);
                labelProcess->Draw();
            histo2DPsiPairDummy->Draw("axis,same");
            canvasPsiPairPlots->SaveAs(Form("%s/PsiPairPt_AllRec_%s.%s",outputDirSp.Data(), DataSets[i].Data() ,suffix.Data()));

            histo2DPsiPairDummy->Draw("copy");
            fHistGammaPsiPairPt->Draw("col,same");

                DrawGammaLines( 0.1,  0.1, 0.04, 15, 3, kMagenta+2, 9);
                DrawGammaLines( -0.1,  -0.1, 0.04, 15, 3, kMagenta+2, 9);
                DrawGammaLines(0.15,  0.15, 0.04, 15, 3, kMagenta+1, 7);
                DrawGammaLines(-0.15,  -0.15, 0.04, 15, 3, kMagenta+1, 7);
                labelEnergy->Draw();
                labelProcess->Draw();
                legendPsiPairPtPlotFits->Draw();
            histo2DPsiPairDummy->Draw("axis,same");
            canvasPsiPairPlots->SaveAs(Form("%s/PsiPairPt_Final_withCuts_%s.%s",outputDirSp.Data(), DataSets[i].Data(),suffix.Data()));
            ex1->Draw();

        }else cout << "INFO: Object |histoGammaPsiPairPt| could not be found! Skipping Draw..." << endl;
        //-----------------------------------
        if(isMC[i]){
            SetXRange(fHistTrueMCKindPsiPairPt[0],minB,maxB);
            SetYRange(fHistTrueMCKindPsiPairPt[0],fHistTrueMCKindPsiPairPt[0]->GetYaxis()->FindBin(0.1),maxYB);
            SetZMinMaxTH2(fHistTrueMCKindPsiPairPt[0],minB,maxB,fHistTrueMCKindPsiPairPt[0]->GetYaxis()->FindBin(0.1),maxYB,kTRUE);
            //-----------------------------------
            TH1D* histoMonteCarloTruePrimGammaPsiPair = (TH1D*)fHistTrueMCKindPsiPairPt[0]->ProjectionX("histoMonteCarloTruePrimGammaPsiPair");
            ConvGammaRebinWithBinCorrection(histoMonteCarloTruePrimGammaPsiPair,1);
            vecPsi[1].push_back(histoMonteCarloTruePrimGammaPsiPair);
            //-----------------------------------
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                                 fHistTrueMCKindPsiPairPt[0],"",
                                 "#psi_{pair}", "#it{p}_{T,#gamma} (GeV/#it{c})",1,1.4,
                                 xPosLabel2D,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistTrueMCKindPsiPairPt[0], Form("%s/PsiPair_Pt_TruePrimGamma_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));

            SetXRange(fHistTrueMCKindPsiPairPt[5],minB,maxB);
            SetYRange(fHistTrueMCKindPsiPairPt[5],fHistTrueMCKindPsiPairPt[5]->GetYaxis()->FindBin(0.1),maxYB);
            SetZMinMaxTH2(fHistTrueMCKindPsiPairPt[5],minB,maxB,fHistTrueMCKindPsiPairPt[5]->GetYaxis()->FindBin(0.1),maxYB,kTRUE);
            //-----------------------------------
            TH1D* histoMonteCarloTrueSecGammaPsiPair = (TH1D*)fHistTrueMCKindPsiPairPt[5]->ProjectionX("histoMonteCarloTrueSecGammaPsiPair");
            ConvGammaRebinWithBinCorrection(histoMonteCarloTrueSecGammaPsiPair,1);
            vecPsi[2].push_back(histoMonteCarloTrueSecGammaPsiPair);
            //-----------------------------------
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                                 fHistTrueMCKindPsiPairPt[5],"",
                                 "#psi_{pair}", "#it{p}_{T,#gamma} (GeV/#it{c})",1,1.4,
                                 xPosLabel2D,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistTrueMCKindPsiPairPt[5], Form("%s/PsiPair_Pt_TrueSecGamma_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));

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
                                 xPosLabel2D,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistTrueDalitzGammaPsiPairPt, Form("%s/PsiPair_Pt_TrueDalitz_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));

            //-----------------------------------
            if(fHistGammaPsiPairPt && fHistTrueMCKindPsiPairPt[0] && fHistTrueMCKindPsiPairPt[5] && fHistTrueDalitzGammaPsiPairPt){
                TH2D* fHistBGPsiPairPt = (TH2D*) fHistGammaPsiPairPt->Clone("histoBGPsiPairMonteCarloPt");
                fHistBGPsiPairPt->Sumw2();
                fHistBGPsiPairPt->Add(fHistTrueMCKindPsiPairPt[0],-1.);
                fHistBGPsiPairPt->Add(fHistTrueMCKindPsiPairPt[5],-1.);
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
                                     xPosLabel2D,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
                SaveCanvasAndWriteHistogram(cvsQuadratic, fHistBGPsiPairPt, Form("%s/PsiPair_Pt_BG_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));

                delete fHistBGPsiPairPt;
            }

            // TRUE GAMMAS PLOT
            histo2DPsiPairDummy->Draw("copy");
                // draw true gamma->ee
                fHistTrueMCKindPsiPairPt[0]->Draw("col,same");
                ex1 = new TExec("ex1","PalColor();");
                ex1->Draw();
                fHistTrueMCKindPsiPairPt[0]->Draw("col,same");
                labelColor                  = new TLatex(0.12,0.90, isMC[i] ? "#gamma rec. from e^{+}e^{-} (color)" : "");
                SetStyleTLatex( labelColor, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,11);
                labelColor->Draw();
                labelProcess                     = new TLatex(0.95,0.86,"#gamma candidates (MC true)" );
                SetStyleTLatex( labelProcess, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,31);
                labelEnergy->Draw();
                labelProcess->Draw();
            histo2DPsiPairDummy->Draw("axis,same");
            canvasPsiPairPlots->SaveAs(Form("%s/PsiPairPt_trueGamma_woCuts_%s.%s",outputDirSp.Data(), DataSets[i].Data(), suffix.Data()));

            // TRUE GAMMAS AND COMBINATORIAL GAMMA CANDIDATES PLOT
            histo2DPsiPairDummy->Draw("copy");
                // draw true gamma->ee
                fHistTrueMCKindPsiPairPt[0]->Draw("col,same");
                ex1 = new TExec("ex1","PalColor();");
                ex1->Draw();
                fHistTrueMCKindPsiPairPt[0]->Draw("col,same");
                fHistTrueMCKindPsiPairPt[13]->Draw("col,same");
                ex2 = new TExec("ex2","gStyle->SetPalette(9);");
                ex2->Draw();
                fHistTrueMCKindPsiPairPt[13]->Draw("col,same");

                labelColor->Draw();
                labelBW                     = new TLatex(0.12,0.86,"#gamma rec. from #pi^{+}#pi^{-} (BAW)");
                SetStyleTLatex( labelBW, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,11);
                labelBW->Draw();
                labelEnergy->Draw();
                labelProcess->Draw();
            histo2DPsiPairDummy->Draw("axis,same");
            canvasPsiPairPlots->SaveAs(Form("%s/PsiPairPt_trueGammaAndCombPlus_woCuts_%s.%s",outputDirSp.Data(), DataSets[i].Data(), suffix.Data()));

            // TRUE GAMMAS AND COMBINATORIAL GAMMA CANDIDATES PLOT
            histo2DPsiPairDummy->Draw("copy");
                // draw true gamma->ee
                fHistTrueMCKindPsiPairPt[0]->Draw("col,same");
                ex1 = new TExec("ex1","PalColor();");
                ex1->Draw();
                fHistTrueMCKindPsiPairPt[0]->Draw("col,same");

                DrawGammaLines( 0.1,  0.1, 0.04, 15, 3, kMagenta+2, 9);
                DrawGammaLines( -0.1,  -0.1, 0.04, 15, 3, kMagenta+2, 9);
                DrawGammaLines(0.15,  0.15, 0.04, 15, 3, kMagenta+1, 7);
                DrawGammaLines(-0.15,  -0.15, 0.04, 15, 3, kMagenta+1, 7);
                labelColor->Draw();
                labelEnergy->Draw();
                labelProcess->Draw();
                legendPsiPairPtPlotFits->Draw();
            histo2DPsiPairDummy->Draw("axis,same");
            canvasPsiPairPlots->SaveAs(Form("%s/PsiPairPt_TrueGamma_withCuts_%s.%s",outputDirSp.Data(), DataSets[i].Data(), suffix.Data()));
        }


        //    _____ _    _ _____ ___       ____ _______
        //   / ____| |  | |_   _|__ \     |  _ \ _   __|
        //  | |    | |__| | | |    ) |    | |_| | | |
        //  | |    |  __  | | |   / /     |  __/  | |
        //  | |____| |  | |_| |_ / /_     | |     | |
        //   \_____|_|  |_|_____|____|    |_|     |_|
        TCanvas* canvasChi2Plots = new TCanvas("canvasChi2Plots","",200,10,1350,1200);  // gives the page size
        DrawGammaCanvasSettings( canvasChi2Plots, 0.08, 0.02, 0.02, 0.09);
        canvasChi2Plots->SetLogy();
        canvasChi2Plots->SetLogz();
        TH2F * histo2DChi2Dummy = new TH2F("histo2DChi2Dummy","histo2DChi2Dummy",200,0.00,50,1000,0.04,200);
        SetStyleHistoTH2ForGraphs(histo2DChi2Dummy, "#chi^{2}/NDF","#it{p}_{T} (GeV/#it{c})",0.035,0.04, 0.035,0.04, 0.98,0.9);
        canvasChi2Plots->cd();
        TF1* funcOldCutChi2PtDummy = new TF1("funcOldCutChi2PtDummy","x/[0]",0.005,1.);
        DrawGammaSetMarkerTF1( funcOldCutChi2PtDummy, 9, 3, kMagenta+2);
        TF1* funcOldCutChi2PtDummy1 = new TF1("funcOldCutChi2PtDummy1","x/[0]",0.005,1.);
        DrawGammaSetMarkerTF1( funcOldCutChi2PtDummy1, 7, 3, kMagenta+1);
        TLegend* legendChi2PlotFits  = GetAndSetLegend2(0.74, 0.75, 0.95, 0.75+(0.035*2*1.35), 0.85*textSizeLabelsPixel);
        legendChi2PlotFits->AddEntry(funcOldCutChi2PtDummy, "#chi^{2}_{max} = 30","l");
        legendChi2PlotFits->AddEntry(funcOldCutChi2PtDummy1, "#chi^{2}_{max} = 40","l");
        SetPlotStyle();

        if(fHistGammaChi2PsiPairPt){
            GetMinMaxBin(fHistGammaChi2Pt,minB,maxB); maxB+=5;
            GetMinMaxBinY(fHistGammaChi2Pt,minYB,maxYB); maxYB+=5;
            SetXRange(fHistGammaChi2Pt,minB,maxB);
            SetYRange(fHistGammaChi2Pt,fHistGammaChi2Pt->GetYaxis()->FindBin(0.1),maxYB);
            SetZMinMaxTH2(fHistGammaChi2Pt,minB,maxB,fHistGammaChi2Pt->GetYaxis()->FindBin(0.1),maxYB,kTRUE);
            //-----------------------------------
            TH1D* histoGammaChi2 = (TH1D*)fHistGammaChi2Pt->ProjectionX("histoGammaChi2");
            ConvGammaRebinWithBinCorrection(histoGammaChi2,1);
            CheckForNegativeEntries(histoGammaChi2);
            vecChi2[0].push_back(histoGammaChi2);
            //-----------------------------------
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                                 fHistGammaChi2Pt,"",
                                 "#chi^{2}_{#gamma}/ndf", "#it{p}_{T,#gamma} (GeV/#it{c})",1,1.4,
                                 xPosLabel2D,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistGammaChi2Pt, Form("%s/Chi2_Pt_Photon_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));

            canvasChi2Plots->cd();
            histo2DChi2Dummy->GetZaxis()->SetRangeUser(1/nGammas,fHistGammaChi2Pt->GetMaximum());

            // RECONSTRUCTED PLOT
            histo2DChi2Dummy->Draw("copy");
                fHistGammaChi2Pt->Draw("col,same");
                labelEnergy                  = new TLatex(0.95,0.90,Form("%s, %s",labelALICEforPlots.Data(),fCollisionSystem.Data()));
                SetStyleTLatex( labelEnergy, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,31);
                labelEnergy->Draw();
                labelProcess                     = new TLatex(0.95,0.86,isMC[i] ? "#gamma candidates (MC rec.)" : "#gamma candidates (data)");
                SetStyleTLatex( labelProcess, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,31);
                labelProcess->Draw();
            histo2DChi2Dummy->Draw("axis,same");
            canvasChi2Plots->SaveAs(Form("%s/Chi2Pt_AllRec_%s.%s",outputDir.Data(), DataSets[i].Data(), suffix.Data()));

            histo2DChi2Dummy->Draw("copy");
                fHistGammaChi2Pt->Draw("col,same");
                DrawGammaLines( 30.,  30., 0.04, 15, 3, kMagenta+2, 9);
                DrawGammaLines(40.,  40., 0.04, 15, 3, kMagenta+1, 7);
                legendChi2PlotFits->Draw();
                labelEnergy->Draw();
                labelProcess->Draw();
            histo2DChi2Dummy->Draw("axis,same");
            canvasChi2Plots->SaveAs(Form("%s/Chi2Pt_Final_withCuts_%s.%s",outputDirSp.Data(), DataSets[i].Data(), suffix.Data()));
        }else cout << "INFO: Object |histoGammaChi2Pt| could not be found! Skipping Draw..." << endl;
        //-----------------------------------
        if(isMC[i]){
            SetXRange(fHistTrueMCKindChi2Pt[0],minB,maxB);
            SetYRange(fHistTrueMCKindChi2Pt[0],fHistTrueMCKindChi2Pt[0]->GetYaxis()->FindBin(0.1),maxYB);
            SetZMinMaxTH2(fHistTrueMCKindChi2Pt[0],minB,maxB,fHistTrueMCKindChi2Pt[0]->GetYaxis()->FindBin(0.1),maxYB,kTRUE);
            //-----------------------------------
            TH1D* histoMonteCarloTruePrimGammaChi2 = (TH1D*)fHistTrueMCKindChi2Pt[0]->ProjectionX("histoMonteCarloTruePrimGammaChi2");
            ConvGammaRebinWithBinCorrection(histoMonteCarloTruePrimGammaChi2,1);
            vecChi2[1].push_back(histoMonteCarloTruePrimGammaChi2);
            //-----------------------------------
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                                 fHistTrueMCKindChi2Pt[0],"",
                                 "#chi^{2}_{#gamma}/ndf", "#it{p}_{T,#gamma} (GeV/#it{c})",1,1.4,
                                 xPosLabel2D,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistTrueMCKindChi2Pt[0], Form("%s/Chi2_Pt_TruePrimGamma_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));

            SetXRange(fHistTrueMCKindChi2Pt[5],minB,maxB);
            SetYRange(fHistTrueMCKindChi2Pt[5],fHistTrueMCKindChi2Pt[5]->GetYaxis()->FindBin(0.1),maxYB);
            SetZMinMaxTH2(fHistTrueMCKindChi2Pt[5],minB,maxB,fHistTrueMCKindChi2Pt[5]->GetYaxis()->FindBin(0.1),maxYB,kTRUE);
            //-----------------------------------
            TH1D* histoMonteCarloTrueSecGammaChi2 = (TH1D*)fHistTrueMCKindChi2Pt[5]->ProjectionX("histoMonteCarloTrueSecGammaChi2");
            ConvGammaRebinWithBinCorrection(histoMonteCarloTrueSecGammaChi2,1);
            vecChi2[2].push_back(histoMonteCarloTrueSecGammaChi2);
            //-----------------------------------
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                                 fHistTrueMCKindChi2Pt[5],"",
                                 "#chi^{2}_{#gamma}/ndf", "#it{p}_{T,#gamma} (GeV/#it{c})",1,1.4,
                                 xPosLabel2D,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistTrueMCKindChi2Pt[5], Form("%s/Chi2_Pt_TrueSecGamma_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));

            SetXRange(fHistTrueDalitzGammaChi2Pt,minB,maxB);
            SetYRange(fHistTrueDalitzGammaChi2Pt,fHistTrueDalitzGammaChi2Pt->GetYaxis()->FindBin(0.1),maxYB);
            SetZMinMaxTH2(fHistTrueDalitzGammaChi2Pt,minB,maxB,fHistTrueDalitzGammaChi2Pt->GetYaxis()->FindBin(0.1),maxYB,kTRUE);
            //-----------------------------------
            TH1D* histoMonteCarloTrueDalitzGammaChi2 = (TH1D*)fHistTrueDalitzGammaChi2Pt->ProjectionX("histoMonteCarloTrueDalitzGammaChi2");
            ConvGammaRebinWithBinCorrection(histoMonteCarloTrueDalitzGammaChi2,1);
            vecChi2[3].push_back(histoMonteCarloTrueDalitzGammaChi2);
            //-----------------------------------
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                                 fHistTrueDalitzGammaChi2Pt,"",
                                 "#chi^{2}_{#gamma}/ndf", "#it{p}_{T,#gamma} (GeV/#it{c})",1,1.4,
                                 xPosLabel2D,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistTrueDalitzGammaChi2Pt, Form("%s/Chi2_Pt_TrueDalitz_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));


            //-----------------------------------
            if(fHistGammaChi2Pt && fHistTrueMCKindChi2Pt[0] && fHistTrueMCKindChi2Pt[5] && fHistTrueDalitzGammaChi2Pt){
                TH2D* fHistBGChi2Pt = (TH2D*) fHistGammaChi2Pt->Clone("histoMonteCarloBGChi2Pt");
                fHistBGChi2Pt->Sumw2();
                fHistBGChi2Pt->Add(fHistTrueMCKindChi2Pt[0],-1.);
                fHistBGChi2Pt->Add(fHistTrueMCKindChi2Pt[5],-1.);
                fHistBGChi2Pt->Add(fHistTrueDalitzGammaChi2Pt,-1.);
                SetXRange(fHistBGChi2Pt,minB,maxB);
                SetYRange(fHistBGChi2Pt,fHistBGChi2Pt->GetYaxis()->FindBin(0.1),maxYB);
                SetZMinMaxTH2(fHistBGChi2Pt,minB,maxB,fHistBGChi2Pt->GetYaxis()->FindBin(0.1),maxYB,kTRUE);
                CheckForNegativeEntries(fHistBGChi2Pt);
                //-----------------------------------
                TH1D* histoMonteCarloBGChi2 = (TH1D*)fHistBGChi2Pt->ProjectionX("histoMonteCarloBGChi2");
                ConvGammaRebinWithBinCorrection(histoMonteCarloBGChi2,1);
                CheckForNegativeEntries(histoMonteCarloBGChi2);
                vecChi2[4].push_back(histoMonteCarloBGChi2);
                //-----------------------------------
                DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                                     fHistBGChi2Pt,"",
                                     "#chi^{2}_{#gamma}/ndf", "#it{p}_{T,#gamma} (GeV/#it{c})",1,1.4,
                                     xPosLabel2D,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
                SaveCanvasAndWriteHistogram(cvsQuadratic, fHistBGChi2Pt, Form("%s/Chi2_Pt_BG_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));

                delete fHistBGChi2Pt;
            }

            // TRUE GAMMAS PLOT
            histo2DChi2Dummy->GetZaxis()->SetRangeUser(1/nGammas,fHistGammaChi2Pt->GetMaximum());
            histo2DChi2Dummy->Draw("copy");
                // draw true gamma->ee
                fHistTrueMCKindChi2Pt[0]->Draw("col,same");
                ex1 = new TExec("ex1","PalColor();");
                ex1->Draw();
                fHistTrueMCKindChi2Pt[0]->Draw("col,same");
                labelColor                  = new TLatex(0.12,0.90, "#gamma rec. from e^{+}e^{-} (color)");
                SetStyleTLatex( labelColor, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,11);
                labelColor->Draw();
                labelProcess                     = new TLatex(0.95,0.86,isMC[i] ? "#gamma candidates (MC true)" : "#gamma candidates (data)");
                SetStyleTLatex( labelProcess, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,31);
                labelEnergy->Draw();
                labelProcess->Draw();
            histo2DChi2Dummy->Draw("axis,same");
            canvasChi2Plots->SaveAs(Form("%s/Chi2Pt_trueGamma_woCut_%ss.%s",outputDirSp.Data(), DataSets[i].Data(), suffix.Data()));


            // draw true gamma->ee
            fHistTrueMCKindChi2Pt[0]->Draw("col,same");
                ex1 = new TExec("ex1","PalColor();");
                ex1->Draw();
                fHistTrueMCKindChi2Pt[0]->Draw("col,same");
                DrawGammaLines( 30.,  30., 0.04, 15, 3, kMagenta+2, 9);
                DrawGammaLines(40.,  40., 0.04, 15, 3, kMagenta+1, 7);
                legendChi2PlotFits->Draw();
                labelColor->Draw();
                labelEnergy->Draw();
                labelProcess->Draw();
                histo2DChi2Dummy->Draw("axis,same");
            canvasChi2Plots->SaveAs(Form("%s/Chi2Pt_trueGamma_withCuts_%s.%s",outputDirSp.Data(), DataSets[i].Data(), suffix.Data()));


            // TRUE GAMMAS AND COMBINATORIAL GAMMA CANDIDATES PLOT
            histo2DChi2Dummy->GetZaxis()->SetRangeUser(1/nGammas,fHistTrueMCKindChi2Pt[13]->GetMaximum()*1000);
            histo2DChi2Dummy->Draw("copy");

                // draw true gamma->ee
                fHistTrueMCKindChi2Pt[0]->Draw("col,same");
                ex1 = new TExec("ex1","PalColor();");
                ex1->Draw();
                fHistTrueMCKindChi2Pt[0]->Draw("col,same");
                fHistTrueMCKindChi2Pt[13]->Draw("col,same");
                ex2 = new TExec("ex2","gStyle->SetPalette(9);");
                ex2->Draw();
                fHistTrueMCKindChi2Pt[13]->Draw("col,same");
                labelColor->Draw();
                labelBW                     = new TLatex(0.12,0.86,"#gamma rec. from #pi^{+}#pi^{-} (BAW)");
                SetStyleTLatex( labelBW, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,11);
                labelBW->Draw();
                labelEnergy->Draw();
                labelProcess->Draw();
            histo2DChi2Dummy->Draw("axis,same");
            canvasChi2Plots->SaveAs(Form("%s/Chi2Pt_trueGammaAndCombPlus_woCuts_%s.%s",outputDirSp.Data(), DataSets[i].Data(), suffix.Data()));
        }


        //   _____   _____ _____   _____        _____ _____         _____ _    _ _____ ___
        //  |  __ \ / ____|_   _| |  __ \ /\   |_   _|  __ \       / ____| |  | |_   _|__ \
        //  | |__) | (___   | |   | |__) /  \    | | | |__) | ___ | |    | |__| | | |    ) |
        //  |  ___/ \___ \  | |   |  ___/ /\ \   | | |  _  /  ___ | |    |  __  | | |   / /
        //  | |     ____) |_| |_  | |  / ____ \ _| |_| | \ \      | |____| |  | |_| |_ / /_
        //  |_|    |_____/|_____| |_| /_/    \_\_____|_|  \_\      \_____|_|  |_|_____|____|

        TCanvas* canvasChi2PsiPairPlots = new TCanvas("canvasChi2PsiPairPlots","",200,10,1350,1200);  // gives the page size
        DrawGammaCanvasSettings( canvasChi2PsiPairPlots, 0.1, 0.02, 0.02, 0.09);
        // canvasChi2PsiPairPlots->SetLogy();
        canvasChi2PsiPairPlots->SetLogz();
        TH2F * histo2DChi2PsiPairDummy = new TH2F("histo2DChi2PsiPairDummy","histo2DChi2PsiPairDummy",200,0.0,50,200,-0.215,0.215);
        SetStyleHistoTH2ForGraphs(histo2DChi2PsiPairDummy, "#chi^{2}/NDF", "#Psi_{pair}",0.035,0.04, 0.035,0.04, 0.98,1.2);
        canvasChi2PsiPairPlots->cd();

        Double_t funcParamChi2PsiPairLin[2][2]  = {{0.1,30}, {0.15, 50}};
        Double_t funcParamChi2PsiPair[3][2]     = {{0.15,-0.065}, {0.18, -0.055}, { 0.20, -0.050}};
        if (fEnergyFlag.CompareTo("XeXe_5.44TeV") == 0 || fEnergyFlag.CompareTo("13TeVLowB") == 0){
            funcParamChi2PsiPair[0][0]      = 0.3;
            funcParamChi2PsiPair[0][1]      = -0.085;
            funcParamChi2PsiPair[1][0]      = 0.35;
            funcParamChi2PsiPair[1][1]      = -0.075;
            funcParamChi2PsiPair[2][0]      = 0.4;
            funcParamChi2PsiPair[2][1]      = -0.065;
        }
        Color_t colorTriFunc[2]                 = {kMagenta+2, kMagenta+4};
        Color_t colorModFunc[3]                 = {kMagenta-4, kMagenta-8, kMagenta+1};

        TF1 *funcChi2PsiPairTri[4]          = {NULL};
        for (Int_t k = 0; k< 2; k++){
            funcChi2PsiPairTri[k*2] =  new TF1(Form("funcChi2PsiPairTri_%d",k*2),Form("-%0.2f+(%0.2f/%f*x)",funcParamChi2PsiPairLin[k][0],funcParamChi2PsiPairLin[k][0],funcParamChi2PsiPairLin[k][1]),0,funcParamChi2PsiPairLin[k][1]);
            DrawGammaSetMarkerTF1( funcChi2PsiPairTri[k*2], 2, 3, colorTriFunc[k]);
            funcChi2PsiPairTri[k*2+1] =   new TF1(Form("funcChi2PsiPairTri_%d",k*2+1),Form("%0.2f-(%0.2f/%f*x)",funcParamChi2PsiPairLin[k][0],funcParamChi2PsiPairLin[k][0],funcParamChi2PsiPairLin[k][1]),0,funcParamChi2PsiPairLin[k][1]);
            DrawGammaSetMarkerTF1( funcChi2PsiPairTri[k*2+1], 2, 3, colorTriFunc[k]);
        }
        TF1 *funcChi2PsiPairExp[6]      = {NULL};
        for (Int_t k = 0; k< 3; k++){
            funcChi2PsiPairExp[k*2] =  new TF1(Form("funcChi2PsiPairExp_%d",k*2),Form("%0.2f*TMath::Exp(%0.3f*x)",funcParamChi2PsiPair[k][0],funcParamChi2PsiPair[k][1]),0,50);
            DrawGammaSetMarkerTF1( funcChi2PsiPairExp[k*2], 2, 3, colorModFunc[k]);
            funcChi2PsiPairExp[k*2+1] =  new TF1(Form("funcChi2PsiPairExp_%d",k*2+1),Form("-%0.2f*TMath::Exp(%0.3f*x)",funcParamChi2PsiPair[k][0],funcParamChi2PsiPair[k][1]),0,50);
            DrawGammaSetMarkerTF1( funcChi2PsiPairExp[k*2+1], 2, 3, colorModFunc[k]);
        }
        TLegend* legendChi2PsiPairPlotFits  = GetAndSetLegend2(0.65, 0.13, 0.95, 0.13+(0.035*5*1.15), 0.75*textSizeLabelsPixel,1, "", 43, 0.1);
        for (Int_t k = 0; k< 2; k++){
            legendChi2PsiPairPlotFits->AddEntry(funcChi2PsiPairTri[k*2], Form("|#Psi_{pair}| < %0.2f/%.0f#chi^{2}+%0.2f",funcParamChi2PsiPairLin[k][0], funcParamChi2PsiPairLin[k][1],funcParamChi2PsiPairLin[k][0]),"l");
        }
        for (Int_t k = 0; k< 3; k++){
            legendChi2PsiPairPlotFits->AddEntry(funcChi2PsiPairExp[k*2], Form("|#Psi_{pair}| < %0.2fexp(%0.3f#chi^{2})",funcParamChi2PsiPair[k][0],funcParamChi2PsiPair[k][1]),"l");
        }

        SetPlotStyle();
        if(fHistGammaChi2PsiPairPt){
            GetMinMaxBin(fHistGammaChi2PsiPair,minB,maxB);
            GetMinMaxBinY(fHistGammaChi2PsiPair,minYB,maxYB); minYB-=5; maxYB+=5;
            SetXRange(fHistGammaChi2PsiPair,minB,maxB);
            SetYRange(fHistGammaChi2PsiPair,minYB,maxYB);
            SetZMinMaxTH2(fHistGammaChi2PsiPair,minB,maxB,minYB,maxYB,kTRUE);
            //-----------------------------------
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                                 fHistGammaChi2PsiPair,"",
                                 "#chi^{2}_{#gamma}/ndf", "#psi_{pair}",1,1.4,
                                 xPosLabel2D,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistGammaChi2PsiPair, Form("%s/Chi2_PsiPair_Photon_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));

            histo2DChi2PsiPairDummy->GetZaxis()->SetRangeUser(8./nGammas,fHistGammaChi2PsiPair->GetMaximum());


            // RECONSTRUCTED PLOT
            canvasChi2PsiPairPlots->cd();
            histo2DChi2PsiPairDummy->Draw("copy");
            fHistGammaChi2PsiPair->Draw("col,same");
            labelEnergy                  = new TLatex(0.95,0.90,Form("%s, %s",labelALICEforPlots.Data(),fCollisionSystem.Data()));
            SetStyleTLatex( labelEnergy, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,31);
            labelEnergy->Draw();
            labelProcess                     = new TLatex(0.95,0.86,isMC[i] ? "#gamma candidates (MC rec.)" : "#gamma candidates (data)");
            SetStyleTLatex( labelProcess, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,31);
            labelProcess->Draw();
            histo2DChi2PsiPairDummy->Draw("axis,same");
            canvasChi2PsiPairPlots->SaveAs(Form("%s/Chi2PsiPair_AllRec_%s.%s",outputDirSp.Data(),DataSets[i].Data(),suffix.Data()));

                for (Int_t k = 0; k< 2; k++){
                    funcChi2PsiPairTri[k*2]->Draw("same");
                    funcChi2PsiPairTri[k*2+1]->Draw("same");
                }
                for (Int_t k = 0; k< 3; k++){
                    funcChi2PsiPairExp[k*2]->Draw("same");
                    funcChi2PsiPairExp[k*2+1]->Draw("same");
                }
                legendChi2PsiPairPlotFits->Draw();
                labelEnergy->Draw();
                labelProcess->Draw();
                histo2DChi2PsiPairDummy->Draw("axis,same");
            canvasChi2PsiPairPlots->SaveAs(Form("%s/Chi2PsiPair_AllRec_withCuts_%s.%s",outputDirSp.Data(),  DataSets[i].Data(), suffix.Data()));

            for(Int_t bin=0; bin<15; bin++){
                histo2DChi2PsiPairDummy->GetZaxis()->SetRangeUser(1/nGammas,fHistGammaChi2PsiPair->GetMaximum());
                histo2DChi2PsiPairDummy->Draw("copy");
                    fHistGammaChi2PsiPairPtSliced[bin]->Draw("col,same");
                    labelpTrange    = new TLatex(0.95,0.82,Form("%2.1f < #it{p}_{T} < %2.1f GeV/#it{c}",projBinnin[bin],projBinnin[bin+1]));
                    SetStyleTLatex( labelpTrange, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,31);
                    labelpTrange->Draw();
                    labelEnergy->Draw();
                    labelProcess->Draw();
                histo2DChi2PsiPairDummy->Draw("axis,same");
                canvasChi2PsiPairPlots->SaveAs(Form("%s/Chi2PsiPair/%dbin_Chi2PsiPair_vs_Pt_AllRec_%s.%s",outputDirSp.Data(),bin,DataSets[i].Data(),suffix.Data()));

            }
            ex1->Draw();
            histo2DChi2PsiPairDummy->GetZaxis()->SetRangeUser(8/nGammas,fHistGammaChi2PsiPair->GetMaximum());
        }else cout << "INFO: Object |histoGammaChi2PsiPair| could not be found! Skipping Draw..." << endl;
        //-----------------------------------
        if(isMC[i]){
            SetXRange(fHistTrueMCKindChi2PsiPair[0],minB,maxB);
            SetYRange(fHistTrueMCKindChi2PsiPair[0],minYB,maxYB);
            SetZMinMaxTH2(fHistTrueMCKindChi2PsiPair[0],minB,maxB,minYB,maxYB,kTRUE);
            //-----------------------------------
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                                 fHistTrueMCKindChi2PsiPair[0],"",
                                 "#chi^{2}_{#gamma}/ndf", "#psi_{pair}",1,1.4,
                                 xPosLabel2D,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistTrueMCKindChi2PsiPair[0], Form("%s/Chi2_PsiPair_TruePrimGamma_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));

            SetXRange(fHistTrueMCKindChi2PsiPair[5],minB,maxB);
            SetYRange(fHistTrueMCKindChi2PsiPair[5],minYB,maxYB);
            SetZMinMaxTH2(fHistTrueMCKindChi2PsiPair[5],minB,maxB,minYB,maxYB,kTRUE);
            //-----------------------------------
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                                 fHistTrueMCKindChi2PsiPair[5],"",
                                 "#chi^{2}_{#gamma}/ndf", "#psi_{pair}",1,1.4,
                                 xPosLabel2D,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistTrueMCKindChi2PsiPair[5], Form("%s/Chi2_PsiPair_TrueSecGamma_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));

            SetXRange(fHistTrueDalitzGammaChi2PsiPair,minB,maxB);
            SetYRange(fHistTrueDalitzGammaChi2PsiPair,minYB,maxYB);
            SetZMinMaxTH2(fHistTrueDalitzGammaChi2PsiPair,minB,maxB,minYB,maxYB,kTRUE);
            //-----------------------------------
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                                 fHistTrueDalitzGammaChi2PsiPair,"",
                                 "#chi^{2}_{#gamma}/ndf", "#psi_{pair}",1,1.4,
                                 xPosLabel2D,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistTrueDalitzGammaChi2PsiPair, Form("%s/Chi2_PsiPair_TrueDalitz_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));

            //-----------------------------------
            if(fHistGammaChi2PsiPair && fHistTrueMCKindChi2PsiPair[0] && fHistTrueMCKindChi2PsiPair[5] &&fHistTrueDalitzGammaChi2PsiPair){
                TH2D* fHistBGChi2PsiPair = (TH2D*) fHistGammaChi2PsiPair->Clone("histoBGChi2PsiPairMonteCarlo");
                fHistBGChi2PsiPair->Sumw2();
                fHistBGChi2PsiPair->Add(fHistTrueMCKindChi2PsiPair[0],-1.);
                fHistBGChi2PsiPair->Add(fHistTrueMCKindChi2PsiPair[5],-1.);
                fHistBGChi2PsiPair->Add(fHistTrueDalitzGammaChi2PsiPair,-1.);
                SetXRange(fHistBGChi2PsiPair,minB,maxB);
                SetYRange(fHistBGChi2PsiPair,minYB,maxYB);
                SetZMinMaxTH2(fHistBGChi2PsiPair,minB,maxB,minYB,maxYB,kTRUE);
                CheckForNegativeEntries(fHistBGChi2PsiPair);
                //-----------------------------------
                DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                                     fHistBGChi2PsiPair,"",
                                     "#chi^{2}_{#gamma}/ndf", "#psi_{pair}",1,1.4,
                                     xPosLabel2D,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
                SaveCanvasAndWriteHistogram(cvsQuadratic, fHistBGChi2PsiPair, Form("%s/Chi2_PsiPair_BG_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));

                delete fHistBGChi2PsiPair;
            }

            canvasChi2PsiPairPlots->cd();
            // TRUE GAMMAS PLOT
                histo2DChi2PsiPairDummy->Draw("copy");
                // draw true gamma->ee
                histoTrueGammaChi2PsiPair->Draw("col,same");
                ex1 = new TExec("ex1","PalColor();");
                ex1->Draw();
                histoTrueGammaChi2PsiPair->Draw("col,same");
                labelColor                  = new TLatex(0.12,0.90, "#gamma rec. from e^{+}e^{-} (color)");
                SetStyleTLatex( labelColor, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,11);
                labelColor->Draw();
                labelProcess                     = new TLatex(0.95,0.86, "#gamma candidates (MC true)");
                SetStyleTLatex( labelProcess, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,31);
                labelEnergy->Draw();
                labelProcess->Draw();
            histo2DChi2PsiPairDummy->Draw("axis,same");
            canvasChi2PsiPairPlots->SaveAs(Form("%s/Chi2PsiPair_trueGamma_woCuts.%s",outputDirSp.Data(),suffix.Data()));

                // draw true gamma->ee
                for (Int_t k = 0; k< 2; k++){
                    funcChi2PsiPairTri[k*2]->Draw("same");
                    funcChi2PsiPairTri[k*2+1]->Draw("same");
                }
                for (Int_t k = 0; k< 3; k++){
                    funcChi2PsiPairExp[k*2]->Draw("same");
                    funcChi2PsiPairExp[k*2+1]->Draw("same");
                }
                legendChi2PsiPairPlotFits->Draw();
                labelEnergy->Draw();
                labelProcess->Draw();
            histo2DChi2PsiPairDummy->Draw("axis,same");
            canvasChi2PsiPairPlots->SaveAs(Form("%s/Chi2PsiPair_trueGamma_%s.%s",outputDirSp.Data(), DataSets[i].Data(), suffix.Data()));

            for(Int_t bin=0; bin<15; bin++){
                histo2DChi2PsiPairDummy->Draw("copy");
                // draw true gamma->ee
                fHistTrueMCKindChi2PsiPairPtSliced[0][bin]->Draw("col,same");
                ex1 = new TExec("ex1","PalColor();");
                ex1->Draw();
                fHistTrueMCKindChi2PsiPairPtSliced[0][bin]->Draw("col,same");
                labelColor->Draw();
                labelEnergy->Draw();
                labelProcess->Draw();
                labelpTrange    = new TLatex(0.95,0.82,Form("%2.1f < #it{p}_{T} < %2.1f GeV/#it{c}",projBinnin[bin],projBinnin[bin+1]));
                SetStyleTLatex( labelpTrange, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,31);
                labelpTrange->Draw();
                histo2DChi2PsiPairDummy->Draw("axis,same");
                canvasChi2PsiPairPlots->SaveAs(Form("%s/Chi2PsiPair/%dbin_Chi2PsiPair_vs_Pt_trueGamma_%s.%s",outputDirSp.Data(),bin, DataSets[i].Data(), suffix.Data()));
            }

            labelBW                     = new TLatex(0.12,0.86,"#gamma rec. from e^{#pm}#pi^{#pm} (BAW)");
            SetStyleTLatex( labelBW, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,11);

            for(Int_t bin=0; bin<15; bin++){
                histo2DChi2PsiPairDummy->GetZaxis()->SetRangeUser(1/nGammas,fHistTrueMCKindChi2PsiPairPtSliced[11][bin]->GetMaximum()*10);
                histo2DChi2PsiPairDummy->Draw("copy");
                // draw true gamma->ee
                fHistTrueMCKindChi2PsiPairPtSliced[0][bin]->Draw("col,same");
                ex1 = new TExec("ex1","PalColor();");
                ex1->Draw();
                fHistTrueMCKindChi2PsiPairPtSliced[0][bin]->Draw("col,same");
                // draw ee combinatorics
                fHistTrueMCKindChi2PsiPairPtSliced[11][bin]->Draw("col,same");
                ex2 = new TExec("ex2","gStyle->SetPalette(9);");
                ex2->Draw();
                fHistTrueMCKindChi2PsiPairPtSliced[11][bin]->Draw("col,same");
                labelColor->Draw();
                labelBW->Draw();
                labelEnergy->Draw();
                labelProcess->Draw();
                labelpTrange    = new TLatex(0.95,0.82,Form("%2.1f < #it{p}_{T} < %2.1f GeV/#it{c}",projBinnin[bin],projBinnin[bin+1]));
                SetStyleTLatex( labelpTrange, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,31);
                labelpTrange->Draw();
                histo2DChi2PsiPairDummy->Draw("axis,same");
                canvasChi2PsiPairPlots->SaveAs(Form("%s/Chi2PsiPair/%dbin_Chi2PsiPair_vs_Pt_trueGammaAndComb_%s.%s",outputDirSp.Data(),bin, DataSets[i].Data(), suffix.Data()));
            }

            histo2DChi2PsiPairDummy->GetZaxis()->SetRangeUser(1/nGammas,fHistTrueMCKindChi2PsiPair[11]->GetMaximum()*1000);
            // TRUE GAMMAS AND COMBINATORICS PLOT
            histo2DChi2PsiPairDummy->Draw("copy");
                // draw true gamma->ee
                histoTrueGammaChi2PsiPair->Draw("col,same");
                ex1 = new TExec("ex1","PalColor();");
                ex1->Draw();
                histoTrueGammaChi2PsiPair->Draw("col,same");
                // draw ee combinatorics
                fHistTrueMCKindChi2PsiPair[11]->Draw("col,same");
                ex2 = new TExec("ex2","gStyle->SetPalette(9);");
                ex2->Draw();
                fHistTrueMCKindChi2PsiPair[11]->Draw("col,same");
                // draw pipi combinatorics
                fHistTrueMCKindChi2PsiPair[13]->Draw("col,same");
                ex2 = new TExec("ex2","gStyle->SetPalette(9);");
                ex2->Draw();
                fHistTrueMCKindChi2PsiPair[13]->Draw("col,same");
                labelColor                  = new TLatex(0.12,0.90, "#gamma rec. from e^{+}e^{-} (color)");
                SetStyleTLatex( labelColor, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,11);
                labelColor->Draw();
                labelBW                     = new TLatex(0.12,0.86,"#gamma rec. from e^{#pm}#pi^{#pm} (BAW)");
                SetStyleTLatex( labelBW, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,11);
                labelBW->Draw();
                labelProcess                     = new TLatex(0.95,0.86,"#gamma candidates (MC true)");
                SetStyleTLatex( labelProcess, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,31);
                labelEnergy->Draw();
                labelProcess->Draw();
            histo2DChi2PsiPairDummy->Draw("axis,same");
            canvasChi2PsiPairPlots->SaveAs(Form("%s/Chi2PsiPair_trueGammaAndComb_woCuts_%s.%s",outputDirSp.Data(), DataSets[i].Data(), suffix.Data()));

                for (Int_t k = 0; k< 2; k++){
                    funcChi2PsiPairTri[k*2]->Draw("same");
                    funcChi2PsiPairTri[k*2+1]->Draw("same");
                }
                for (Int_t k = 0; k< 3; k++){
                    funcChi2PsiPairExp[k*2]->Draw("same");
                    funcChi2PsiPairExp[k*2+1]->Draw("same");
                }
                legendChi2PsiPairPlotFits->Draw();
                labelEnergy->Draw();
                labelProcess->Draw();
            histo2DChi2PsiPairDummy->Draw("axis,same");
            canvasChi2PsiPairPlots->SaveAs(Form("%s/Chi2PsiPair_MCSep_withCuts_%s.%s",outputDirSp.Data(),  DataSets[i].Data(), suffix.Data()));
        }else if (isMC[i]) cout << "INFO: Object |histoTrueDalitzGammaChi2PsiPair| could not be found! Skipping Draw..." << endl;


        //    ___ _______ _   _ _____ ____
        //   / _ \__   __| | | |  ___|  _ \
        //  | | | | | |  | |_| | |_  | |_| |
        //  | | | | | |  |  _  |  _| |    /
        //  | |_| | | |  | | | | |___| |\ \
        //   \___/  |_|  |_| |_|_____|_| \_\

        SetPlotStyle();

        TH2D* fHistGammaEtaPt = (TH2D*)TopDir->Get("histoGammaEtaPt");
        if(fHistGammaEtaPt){
            TH1D* fHistGammaPt = (TH1D*)fHistGammaEtaPt->ProjectionY("fHistGammaPt",1,fHistGammaEtaPt->GetNbinsX());
            fHistGammaPt->Sumw2();
            fHistGammaPt->Scale(1./nGammas);
            GetMinMaxBin(fHistGammaPt,minB,maxB);
            SetXRange(fHistGammaPt,fHistGammaPt->GetXaxis()->FindBin(0.1),maxB+5);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kTRUE,kTRUE,kFALSE,
                                fHistGammaPt,"","#it{p}_{T} (GeV/#it{c})","#frac{1}{N_{#gamma}} #frac{dN}{d#it{p}_{T}}",1,1,
                                xPosLabel1D,0.94,0.03,fCollisionSystem,plotDataSets[i],"");
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
                                xPosLabel1D,0.94,0.03,fCollisionSystem,plotDataSets[i],"");
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
                                xPosLabel1D,0.94,0.03,fCollisionSystem,plotDataSets[i],"");
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
                                    xPosLabel1D,0.94,0.03,fCollisionSystem,plotDataSets[i],"");
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
                                xPosLabel1D,0.94,0.03,fCollisionSystem,plotDataSets[i],"");
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
                                xPosLabel2D,0.85,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
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
                                xPosLabel2D,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistZR, Form("%s/Photon_ZR_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
        }else cout << "INFO: Object |histoGammaZR| could not be found! Skipping Draw..." << endl;


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
                                xPosLabel2D,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
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
                                xPosLabel2D,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistTruePrimGammaCosPointPt, Form("%s/CosPoint_Pt_TruePrimGamma_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
        }else if (isMC[i]) cout << "INFO: Object |histoTruePrimGammaCosPointPt| could not be found! Skipping Draw..." << endl;
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
                                xPosLabel2D,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistTrueSecGammaCosPointPt, Form("%s/CosPoint_Pt_TrueSecGamma_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
        }else if (isMC[i]) cout << "INFO: Object |histoTrueSecGammaCosPointPt| could not be found! Skipping Draw..." << endl;
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
                                xPosLabel2D,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistTrueDalitzGammaCosPointPt, Form("%s/CosPoint_Pt_TrueDalitz_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
        }else if (isMC[i]) cout << "INFO: Object |histoTrueDalitzGammaCosPointPt| could not be found! Skipping Draw..." << endl;
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
                                xPosLabel2D,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
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
                                xPosLabel2D,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
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
                                xPosLabel2D,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistTrueGammaCosPointChi2, Form("%s/CosPoint_Chi2_TrueGamma_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
        }else if (isMC[i]) cout << "INFO: Object |histoTrueGammaCosPointChi2| could not be found! Skipping Draw..." << endl;
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
                                xPosLabel2D,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistTrueDalitzGammaCosPointChi2, Form("%s/CosPoint_Chi2_TrueDalitz_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
        }else if (isMC[i]) cout << "INFO: Object |histoTrueDalitzGammaCosPointChi2| could not be found! Skipping Draw..." << endl;
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
                                xPosLabel2D,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistBGCosPointChi2, Form("%s/CosPoint_Chi2_BG_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));

            delete fHistBGCosPointChi2;
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
                                xPosLabel2D,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
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
                                xPosLabel2D,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistTrueGammaCosPointPsiPair, Form("%s/CosPoint_PsiPair_TrueGamma_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
        }else if (isMC[i]) cout << "INFO: Object |histoTrueGammaCosPointPsiPair| could not be found! Skipping Draw..." << endl;
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
                                xPosLabel2D,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistTrueDalitzGammaCosPointPsiPair, Form("%s/CosPoint_PsiPair_TrueDalitz_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
        }else if (isMC[i]) cout << "INFO: Object |histoTrueDalitzGammaCosPointPsiPair| could not be found! Skipping Draw..." << endl;
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
                                xPosLabel2D,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
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
                                xPosLabel2D,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
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
                                xPosLabel2D,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistTruePrimGammaAsymP, Form("%s/Asymmetry_P_TruePrimGamma_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
        }else if (isMC[i]) cout << "INFO: Object |histoTruePrimGammaAsymP| could not be found! Skipping Draw..." << endl;
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
                                xPosLabel2D,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistTrueSecGammaAsymP, Form("%s/Asymmetry_P_TrueSecGamma_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
        }else if (isMC[i]) cout << "INFO: Object |histoTrueSecGammaAsymP| could not be found! Skipping Draw..." << endl;
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
                                xPosLabel2D,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistTrueDalitzGammaAsymP, Form("%s/Asymmetry_P_TrueDalitz_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
        }else if (isMC[i]) cout << "INFO: Object |histoTrueDalitzGammaAsymP| could not be found! Skipping Draw..." << endl;
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
                                xPosLabel2D,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
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
                                xPosLabel2D,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
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
                                xPosLabel2D,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistTrueGammaInvMassPt, Form("%s/InvMass_Pt_TrueGamma_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
        }else if (isMC[i]) cout << "INFO: Object |histoTrueGammaInvMassPt| could not be found! Skipping Draw..." << endl;
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
                                xPosLabel2D,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistTrueDalitzGammaInvMassPt, Form("%s/InvMass_Pt_TrueDalitz_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
        }else if (isMC[i]) cout << "INFO: Object |histoTrueDalitzGammaInvMassPt| could not be found! Skipping Draw..." << endl;
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
                                xPosLabel2D,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistBGInvMassPt, Form("%s/InvMass_Pt_BG_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));

            delete fHistBGInvMassPt;
        }

        //-----------------------------------
        //-------------------------------------- Photon Qt ----------------------------------------------
        //-----------------------------------

        SetPlotStyle();

        //-----------------------------------
        //-----------------------------------
        //-- Leptons
        //-----------------------------------
        //-----------------------------------
        TH2D* histodEdx_bothLeptons = 0x0;
        TH2D* NSigmadEdxEta_bothLeptons = 0x0;

        for(Int_t iL=0; iL<2; iL++){
            TH3D* dEdxEtaP = (TH3D*)TopDir->Get(Form("histo%sdEdxEtaP",lepton[iL].Data()));
            TString labelNLeptons = Form("#frac{1}{N_{e^{%s}}}",charge[iL].Data());
            Double_t nLeptons = 1.;
            if(dEdxEtaP){
                nLeptons = dEdxEtaP->Integral();
                TH2D* histodEdxP = (TH2D*)dEdxEtaP->Project3D("xz");
                histodEdxP->Sumw2();
                if(iL==0) histodEdx_bothLeptons = (TH2D*) histodEdxP->Clone("histodEdx_bothLeptons");
                else if(histodEdx_bothLeptons) histodEdx_bothLeptons->Add(histodEdxP,1.);
                histodEdxP->Scale(1./nLeptons);
                GetMinMaxBin(histodEdxP,minB,maxB);
                //SetXRange(histodEdxP,minB,maxB);
                SetZMinMaxTH2(histodEdxP,minB,maxB,1,histodEdxP->GetNbinsY(),kTRUE);
                DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                    histodEdxP,"",
                                    Form("#it{p}_{e^{%s}} (GeV/#it{c})",charge[iL].Data()),Form("d#it{E}_{e^{%s}-cand} /d#it{x}",charge[iL].Data()),1,1.4,
                                    xPosLabel2D,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
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
                                    xPosLabel1D,0.94,0.03,fCollisionSystem,plotDataSets[i],"");
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
                if(iL==0) NSigmadEdxEta_bothLeptons = (TH2D*) histoSigmadEdxP->Clone("NSigmadEdxEta_bothLeptons");
                else if(NSigmadEdxEta_bothLeptons) NSigmadEdxEta_bothLeptons->Add(histoSigmadEdxP,1.);
                histoSigmadEdxP->Scale(1./nLeptons);
                GetMinMaxBin(histoSigmadEdxP,minB,maxB);
                //SetXRange(histoSigmadEdxP,minB,maxB);
                SetZMinMaxTH2(histoSigmadEdxP,minB,maxB,1,histoSigmadEdxP->GetNbinsY(),kTRUE);
                DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                    histoSigmadEdxP,"",
                                    Form("#it{p}_{e^{%s}} (GeV/#it{c})",charge[iL].Data()),Form("#it{n} #sigma_{e^{%s}} d#it{E}/d#it{x}",charge[iL].Data()),1,1.4,
                                    xPosLabel2D,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
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
                                xPosLabel2D,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, TrueNSigmadEdxTPCP, Form("%s/nSigma_dEdxTPC_True%s_%s.%s", outputDir.Data(), lepton[iL].Data(), DataSets[i].Data(), suffix.Data()));
            }else if (isMC[i]) cout << Form("INFO: Object |histoTrue%sNSigmadEdxTPCP| could not be found! Skipping Draw...",lepton[iL].Data()) << endl;
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
                                    xPosLabel2D,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
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
                                    xPosLabel2D,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
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
                                xPosLabel2D,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, TrueNSigmadEdxITSP, Form("%s/nSigma_dEdxITS_True%s_%s.%s", outputDir.Data(), lepton[iL].Data(), DataSets[i].Data(), suffix.Data()));
            }else if (isMC[i]) cout << Form("INFO: Object |histoTrue%sNSigmadEdxITSP| could not be found! Skipping Draw...",lepton[iL].Data()) << endl;
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
                                    xPosLabel2D,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
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
                                    xPosLabel2D,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
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
                                    xPosLabel2D,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
                SaveCanvasAndWriteHistogram(cvsQuadratic, TrueNSigmadEdxTOFP, Form("%s/nSigma_dEdxTOF_True%s_%s.%s", outputDir.Data(), lepton[iL].Data(), DataSets[i].Data(), suffix.Data()));
            }else if (isMC[i]) cout << Form("INFO: Object |histoTrue%sNSigmadEdxTOFP| could not be found! Skipping Draw...",lepton[iL].Data()) << endl;
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
                                    xPosLabel2D,0.95,0.03,fCollisionSystem,plotDataSets[i],"");
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
                                    xPosLabel1D,0.94,0.03,fCollisionSystem,plotDataSets[i],"");
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
                                        xPosLabel1D,0.94,0.03,fCollisionSystem,plotDataSets[i],"");
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
                                    xPosLabel1D,0.94,0.03,fCollisionSystem,plotDataSets[i],"");
                SaveCanvasAndWriteHistogram(canvas, TrueRcls, Form("%s/True%s_%iITSClvsR_%s.%s", outputDir.Data(), lepton[iL].Data(), iB, DataSets[i].Data(), suffix.Data()));
                vecTrueITSclsR[iL][iB].push_back(TrueRcls);
            }
            }else if (isMC[i]) cout << Form("INFO: Object |histoTrue%sITSClR| could not be found! Skipping Draw...",lepton[iL].Data()) << endl;
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
                                    xPosLabel1D,0.94,0.03,fCollisionSystem,plotDataSets[i],"");
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
                                    xPosLabel1D,0.94,0.03,fCollisionSystem,plotDataSets[i],"");
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
                                    xPosLabel2D,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
                SaveCanvasAndWriteHistogram(cvsQuadratic, FCLPt, Form("%s/%sFindClusterTPC_Pt_%s.%s", outputDir.Data(), lepton[iL].Data(), DataSets[i].Data(), suffix.Data()));
                //-----------------------------------
                SetXRange(FCLPt,1,FCLPt->GetNbinsX());
                TH1D* FCL = (TH1D*) FCLPt->ProjectionX("FCL",1,FCLPt->GetNbinsY());
                GetMinMaxBin(FCL,minB,maxB); minB-=2;
                SetXRange(FCL,minB,maxB);
                DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kFALSE,kFALSE,
                                    FCL,"",Form("TPC Clusters/Findable Clusters e^{%s}",charge[iL].Data()),Form("%s #frac{dN}{dRatio Cl}",labelNLeptons.Data()),1,1,
                                    0.15,0.94,0.03,fCollisionSystem,plotDataSets[i],"",11);
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
                                    xPosLabel2D,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
                SaveCanvasAndWriteHistogram(cvsQuadratic, TPCCLPt, Form("%s/%sClusterTPC_Pt_%s.%s", outputDir.Data(), lepton[iL].Data(), DataSets[i].Data(), suffix.Data()));
                //-----------------------------------
                SetXRange(TPCCLPt,1,TPCCLPt->GetNbinsX());
                TH1D* TPCCL = (TH1D*) TPCCLPt->ProjectionX("FCL",1,FCLPt->GetNbinsY());
                GetMinMaxBin(TPCCL,minB,maxB);
                SetXRange(TPCCL,minB-5,maxB+5);
                DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kFALSE,kFALSE,
                                    TPCCL,"",Form("TPC Clusters e^{%s}",charge[iL].Data()),Form("%s #frac{dN}{dTPC Cl}",labelNLeptons.Data()),1,1,
                                    0.15,0.94,0.03,fCollisionSystem,plotDataSets[i],"",11);
                SaveCanvasAndWriteHistogram(canvas, TPCCL, Form("%s/%sTPCClusters_%s.%s", outputDir.Data(), lepton[iL].Data(), DataSets[i].Data(), suffix.Data()));
                vecTPCCL[iL].push_back(TPCCL);
            }else cout << Form("INFO: Object |histo%sClPt| could not be found! Skipping Draw...",lepton[iL].Data()) << endl;
    //-----------------------------------

        }

        //-----------------------------------
        if(histodEdx_bothLeptons){
            histodEdx_bothLeptons->Sumw2();
            histodEdx_bothLeptons->Scale(1./histodEdx_bothLeptons->GetEntries());
            SetZMinMaxTH2(histodEdx_bothLeptons,1,histodEdx_bothLeptons->GetNbinsX(),1,histodEdx_bothLeptons->GetNbinsY(),kTRUE);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.117,0.117,0.014,0.092,kTRUE,kFALSE,kTRUE,
                                histodEdx_bothLeptons,"",
                                "#it{p} (GeV/#it{c})","d#it{E}_{e^{#pm}-cand} /d#it{x}",1,1.4,
                                xPosLabel2D,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
            SaveCanvasAndWriteHistogram(cvsQuadratic, histodEdx_bothLeptons, Form("%s/dEdxTPC_bothLeptons_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
            delete histodEdx_bothLeptons;
        }else cout << "INFO: Object |histodEdx_bothLeptons| could not be found! Skipping Draw..." << endl;

        //-----------------------------------
        if(NSigmadEdxEta_bothLeptons){
            NSigmadEdxEta_bothLeptons->Sumw2();
            NSigmadEdxEta_bothLeptons->Scale(1./NSigmadEdxEta_bothLeptons->GetEntries());
            SetZMinMaxTH2(NSigmadEdxEta_bothLeptons,1,NSigmadEdxEta_bothLeptons->GetNbinsX(),1,NSigmadEdxEta_bothLeptons->GetNbinsY(),kTRUE);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.117,0.117,0.014,0.092,kTRUE,kFALSE,kTRUE,
                                NSigmadEdxEta_bothLeptons,"",
                                "#it{p} (GeV/#it{c})","#it{n} #sigma_{e^{#pm}} d#it{E}/d#it{x}",1,1.4,
                                xPosLabel2D,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
            SaveCanvasAndWriteHistogram(cvsQuadratic, NSigmadEdxEta_bothLeptons, Form("%s/nSigma_dEdxTPC_bothLeptons_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
            delete NSigmadEdxEta_bothLeptons;
        }else cout << "INFO: Object |NSigmadEdxEta_bothLeptons| could not be found! Skipping Draw..." << endl;


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
                            xPosLabel1D,0.92,0.03,fCollisionSystem,plotDataSets[i],"");
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

    if(nSets>1 ){

        if((Int_t)vecChi2[0].size()>0 && (Int_t)vecChi2[1].size()>0 && (Int_t)vecChi2[2].size()>0 && (Int_t)vecChi2[3].size()>0 && (Int_t)vecChi2[4].size()>0){
            for(Int_t i=0; i<nSets; i++){
                temp.push_back(vecChi2[0].at(i));
                if(i!=0 && isMC[i]){
                    temp.push_back(vecChi2[1].at(i-1));
                    temp.push_back(vecChi2[2].at(i-1));
                    temp.push_back(vecChi2[3].at(i-1));
                    temp.push_back(vecChi2[4].at(i-1));
                }
            }
            AdjustHistRange(temp,5,5,kFALSE);
            for(Int_t i=1; i<nSets; i++){
                if (!isMC[i]) continue;
                DrawGammaSetMarker(vecChi2[0].at(0),20,0.8, kBlack , kBlack);
                DrawPeriodQAHistoTH1(canvas,0.11, 0.02, 0.06, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecChi2[0].at(0),"","#chi^{2}_{#gamma}/ndf ","#frac{d#it{N}_{#gamma}}{#it{N}_{#gamma} d#chi^{2}}",1,1,
                                    xPosLabel1D,0.92,0.03,fCollisionSystem,plotDataSets[i],"");
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
                if(i!=0 && isMC[i]){
                    temp.push_back(vecPsi[1].at(i-1));
                    temp.push_back(vecPsi[2].at(i-1));
                    temp.push_back(vecPsi[3].at(i-1));
                    temp.push_back(vecPsi[4].at(i-1));
                }
            }
            AdjustHistRange(temp,5,5,kFALSE);
            for(Int_t i=1; i<nSets; i++){
                if (!isMC[i]) continue;
                DrawGammaSetMarker(vecPsi[0].at(0),20,0.8, kBlack , kBlack);
                DrawPeriodQAHistoTH1(canvas,0.11, 0.02, 0.06, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecPsi[0].at(0),"","#psi_{pair} ","#frac{d#it{N}_{#gamma}}{#it{N}_{#gamma} d#psi_{pair}}",1,1,
                                    xPosLabel1D,0.92,0.03,fCollisionSystem,plotDataSets[i],"");
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
                if(i!=0 && isMC[i]){
                    temp.push_back(vecCos[1].at(i-1));
                    temp.push_back(vecCos[2].at(i-1));
                    temp.push_back(vecCos[3].at(i-1));
                    temp.push_back(vecCos[4].at(i-1));
                }
            }
            AdjustHistRange(temp,5,5,kFALSE);
            for(Int_t i=1; i<nSets; i++){
                if (!isMC[i]) continue;
                DrawGammaSetMarker(vecCos[0].at(0),20,0.8, kBlack , kBlack);
                DrawPeriodQAHistoTH1(canvas,0.11, 0.02, 0.06, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecCos[0].at(0),"","cos(#theta_{point})","#frac{d#it{N}_{#gamma}}{#it{N}_{#gamma} dcos(#theta_{point})}",1,1,
                    xPosLabel1D,0.92,0.03,fCollisionSystem,plotDataSets[i],"");
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
                if (!isMC[i]) continue;
                temp.push_back(vecGamma[0].at(i));
            }
            AdjustHistRange(temp,5,5,kFALSE);
            for(Int_t i=1; i<nSets; i++){
                if (!isMC[i]) continue;
                DrawGammaSetMarker(vecGamma[0].at(0),20,0.8, kBlack , kBlack);
                DrawPeriodQAHistoTH1(canvas,0.11, 0.02, 0.06, 0.11,kTRUE,kTRUE,kFALSE,
                                    vecGamma[0].at(0),"","#it{p}_{T} (GeV/#it{c})","#frac{d#it{N}_{#gamma}}{#it{N}_{#gamma} d#it{p}_{T}} (GeV/#it{c})^{-1}",1,1,
                    xPosLabel1D,0.92,0.03,fCollisionSystem,plotDataSets[i],"");
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
                if (!isMC[i]) continue;
                TH1D* histoTruePhotonPt = (TH1D*) vecGamma[1].at(i-1)->Clone("histoTruePhotonPt");
                histoTruePhotonPt->Add(vecGamma[2].at(i-1));
                TH1D* histoPurity = (TH1D*) histoTruePhotonPt->Clone("histoPurity");
                histoPurity->Divide(histoPurity,vecGamma[0].at(i),1,1,"B");

                DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kFALSE,kFALSE,
                                    histoPurity,"","#it{p}_{T} (GeV/#it{c})","#epsilon_{pur}",1,1,
                                    xPosLabel1D,0.94,0.03,fCollisionSystem,plotDataSets[i],"");
                SaveCanvas(canvas, Form("%s/Comparison/Photon_Purity_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()), kTRUE);

                histoPurity->GetXaxis()->SetRangeUser(0.05,60.);
                histoPurity->GetYaxis()->SetRangeUser(0.9,1.05);
                DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kFALSE,kFALSE,
                                    histoPurity,"","#it{p}_{T} (GeV/#it{c})","#epsilon_{pur}",1,1,
                                    xPosLabel1D,0.94,0.03,fCollisionSystem,plotDataSets[i],"");
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
                if(i!=0 && isMC[i]){
                    temp.push_back(vecInvMass[1].at(i-1));
                    temp.push_back(vecInvMass[2].at(i-1));
                    temp.push_back(vecInvMass[3].at(i-1));
                }
            }
            AdjustHistRange(temp,5,5,kFALSE);
            for(Int_t i=1; i<nSets; i++){
                if (!isMC[i]) continue;
                DrawGammaSetMarker(vecInvMass[0].at(0),20,0.8, kBlack , kBlack);
                DrawPeriodQAHistoTH1(canvas,0.11, 0.02, 0.06, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecInvMass[0].at(0),"","M_{e^{+}e^{-}} (GeV/#it{c}^{2})","#frac{d#it{N}_{#gamma}}{#it{N}_{#gamma} d#it{M}_{e^{+}e^{-}}} (GeV/#it{c}^{2})^{-1}",1,1,
                    xPosLabel1D,0.94,0.03,fCollisionSystem,plotDataSets[i],"");
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
                        xPosLabel1D,0.94,0.03,fCollisionSystem,plotDataSets[i],"");
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
                    if (!isMC[i]) continue;
                    Int_t iMC = i-1;
                    for(Int_t iB=0; iB<7; iB++) temp.push_back(vecTrueITSclsR[iC][iB].at(iMC));
                    AdjustHistRange(temp,5,5,kFALSE);
                    DrawGammaSetMarker(vecTrueITSclsR[iC][0].at(iMC),20,0.8, kBlack , kBlack);
                    DrawPeriodQAHistoTH1(canvas,0.11, 0.02, 0.06, 0.11,kFALSE,kTRUE,kFALSE,
                                        vecTrueITSclsR[iC][0].at(iMC),"","R (cm)",Form("#frac{d#it{N}_{e^{%s}}}{#it{N}_{e^{%s}} d#it{R}}",charge[iC].Data(),charge[iC].Data()),1,1,
                        xPosLabel1D,0.94,0.03,fCollisionSystem,plotDataSets[i],"");
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
                        xPosLabel1D,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
    SaveCanvas(canvas, Form("%s/Comparison/Photon_Pt.%s", outputDir.Data(), suffix.Data()),kTRUE,kTRUE);

    DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kTRUE,kFALSE,kFALSE,
                        vecGammaPt,"","Photon #it{p}_{T} (GeV/#it{c})","#frac{1}{N_{#gamma}} #frac{dN}{d#it{p}_{T}}",1,1.1,
                        labelData, colorCompare, kTRUE, 1.1, 1.1, kTRUE,
                        xPosLabel1D,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
    SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Photon_Pt.%s", outputDir.Data(), suffix.Data()),kTRUE);
    DeleteVecTH1D(vecGammaPt);
    //-----------------------------------
    DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                        vecGammaEta,"","Photon #eta","#frac{1}{N_{#gamma}} #frac{dN}{d#eta}",1,1.1,
                        labelData, colorCompare, kTRUE, 1.1, 1.1, kTRUE,
                        xPosLabel1D,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
    SaveCanvas(canvas, Form("%s/Comparison/Photon_Eta.%s", outputDir.Data(), suffix.Data()));

    DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                        vecGammaEta,"","Photon #eta","#frac{1}{N_{#gamma}} #frac{dN}{d#eta}",1,1.1,
                        labelData, colorCompare, kTRUE, 1.1, 1.1, kTRUE,
                        xPosLabel1D,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
    SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Photon_Eta.%s", outputDir.Data(), suffix.Data()));
    DeleteVecTH1D(vecGammaEta);
    //-----------------------------------
    DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                        vecGammaAlpha,"","Photon #alpha","#frac{1}{N_{#gamma}} #frac{dN}{d#alpha}",1,1.1,
                        labelData, colorCompare, kTRUE, 1.1, 1.1, kTRUE,
                        xPosLabel1D,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
    SaveCanvas(canvas, Form("%s/Comparison/Photon_Alpha.%s", outputDir.Data(), suffix.Data()));

    DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                        vecGammaAlpha,"","Photon #alpha","#frac{1}{N_{#gamma}} #frac{dN}{d#alpha}",1,1.1,
                        labelData, colorCompare, kTRUE, 1.1, 1.1, kTRUE,
                        xPosLabel1D,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
    SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Photon_Alpha.%s", outputDir.Data(), suffix.Data()));
    DeleteVecTH1D(vecGammaAlpha);
    //-----------------------------------
    DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                        vecGammaPhi,"","Photon #phi","#frac{1}{N_{#gamma}} #frac{dN}{d#phi}",1,1.1,
                        labelData, colorCompare, kTRUE, 1.1, 1.1, kTRUE,
                        xPosLabel1D,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
    SaveCanvas(canvas, Form("%s/Comparison/Photon_Phi.%s", outputDir.Data(), suffix.Data()));

    DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                        vecGammaPhi,"","Photon #phi","#frac{1}{N_{#gamma}} #frac{dN}{d#phi}",1,1.1,
                        labelData, colorCompare, kTRUE, 1.1, 1.1, kTRUE,
                        xPosLabel1D,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
    SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Photon_Phi.%s", outputDir.Data(), suffix.Data()));
    DeleteVecTH1D(vecGammaPhi);
    //-----------------------------------
    DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                        vecGammaPhiEtaNeg,"","Photon #phi C side","#frac{1}{N_{#gamma}} #frac{dN}{d#phi}",1,1.1,
                        labelData, colorCompare, kTRUE, 1.1, 1.1, kTRUE,
                        xPosLabel1D,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
    SaveCanvas(canvas, Form("%s/Comparison/Photon_Phi_EtaNeg.%s", outputDir.Data(), suffix.Data()));

    DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                        vecGammaPhiEtaNeg,"","Photon #phi C side","#frac{1}{N_{#gamma}} #frac{dN}{d#phi}",1,1.1,
                        labelData, colorCompare, kTRUE, 1.1, 1.1, kTRUE,
                        xPosLabel1D,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
    SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Photon_Phi_EtaNeg.%s", outputDir.Data(), suffix.Data()));
    DeleteVecTH1D(vecGammaPhiEtaNeg);
    //-----------------------------------
    DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                        vecGammaPhiEtaPos,"","Photon #phi A side","#frac{1}{N_{#gamma}} #frac{dN}{d#phi}",1,1.1,
                        labelData, colorCompare, kTRUE, 1.1, 1.1, kTRUE,
                        xPosLabel1D,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
    SaveCanvas(canvas, Form("%s/Comparison/Photon_Phi_EtaPos.%s", outputDir.Data(), suffix.Data()));

    DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                        vecGammaPhiEtaPos,"","Photon #phi A side","#frac{1}{N_{#gamma}} #frac{dN}{d#phi}",1,1.1,
                        labelData, colorCompare, kTRUE, 1.1, 1.1, kTRUE,
                        xPosLabel1D,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
    SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Photon_Phi_EtaPos.%s", outputDir.Data(), suffix.Data()));
    DeleteVecTH1D(vecGammaPhiEtaPos);
    //-----------------------------------
    DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                vecChi2[0],"","#chi^{2}_{#gamma}/ndf","#frac{1}{N_{#gamma}} #frac{dN}{d#chi^{2}}",1,1.1,
                                labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                xPosLabel1D,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
    SaveCanvas(canvas, Form("%s/Comparison/Chi2.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);

    DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecChi2[0],"","#chi^{2}_{#gamma}/ndf","#frac{1}{N_{#gamma}} #frac{dN}{d#chi^{2}}",1,1.1,
                                    labelData, colorCompare, kTRUE, 1.1, 1.1, kTRUE,
                                    xPosLabel1D,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
    SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Chi2.%s", outputDir.Data(), suffix.Data()));
    DeleteVecTH1D(vecChi2[0]);
    //-----------------------------------
    DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                vecPsi[0],"","#psi_{pair}","#frac{1}{N_{#gamma}} #frac{dN}{d#psi_{pair}}",1,1.1,
                                labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                xPosLabel1D,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
    SaveCanvas(canvas, Form("%s/Comparison/PsiPair.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);

    DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecPsi[0],"","#psi_{pair}","#frac{1}{N_{#gamma}} #frac{dN}{d#psi_{pair}}",1,1.1,
                                    labelData, colorCompare, kTRUE, 1.1, 1.1, kTRUE,
                                    xPosLabel1D,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
    SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_PsiPair.%s", outputDir.Data(), suffix.Data()));
    DeleteVecTH1D(vecPsi[0]);
    //-----------------------------------
    DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                vecCos[0],"","cos(#theta_{point})","#frac{1}{N_{#gamma}} #frac{d#it{N}_{#gamma}}{dcos(#theta_{point})}",1,1.1,
                                labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                xPosLabel1D,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
    SaveCanvas(canvas, Form("%s/Comparison/CosPoint.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);

    DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecCos[0],"","cos(#theta_{point})","#frac{1}{N_{#gamma}} #frac{d#it{N}_{#gamma}}{dcos(#theta_{point})}",1,1.1,
                                    labelData, colorCompare, kTRUE, 1.1, 1.1, kFALSE,
                                    xPosLabel1D,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
    SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_CosPoint.%s", outputDir.Data(), suffix.Data()));
    DeleteVecTH1D(vecCos[0]);
    //-----------------------------------
    GetMinMaxBin(vecInvMass[0].at(0),minB,maxB);
    for(Int_t i=0; i<nSets; i++) SetXRange(vecInvMass[0].at(i),minB,maxB+5);
    DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                vecInvMass[0],"","M_{e^{+}e^{-}} (GeV/#it{c}^{2})", "#frac{d#it{N}_{#gamma}}{#it{N}_{#gamma} d#it{M}_{e^{+}e^{-}}} (GeV/#it{c}^{2})^{-1}",1,1.1,
                                labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                xPosLabel1D,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
    SaveCanvas(canvas, Form("%s/Comparison/InvMass.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);

    DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecInvMass[0],"","M_{e^{+}e^{-}} (GeV/#it{c}^{2})", "#frac{d#it{N}_{#gamma}}{#it{N}_{#gamma} d#it{M}_{e^{+}e^{-}}} (GeV/#it{c}^{2})^{-1}",1,1.1,
                                    labelData, colorCompare, kTRUE, 1.1, 1.1, kFALSE,
                                    xPosLabel1D,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
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
                            xPosLabel1D,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
        SaveCanvas(canvas, Form("%s/Comparison/%s_Eta.%s", outputDir.Data(),lepton[iL].Data(), suffix.Data()));

        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                            vecEta[iL],"",Form("%s #eta",lepton[iL].Data()),Form("%s #frac{dN}{d#eta}",labelNLeptons.Data()),1,1.1,
                            labelData, colorCompare, kTRUE, 1.1, 1.1, kTRUE,
                            xPosLabel1D,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_%s_Eta.%s", outputDir.Data(),lepton[iL].Data(), suffix.Data()));
        DeleteVecTH1D(vecEta[iL]);
        //-----------------------------------
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kTRUE,kTRUE,kFALSE,
                            vecPt[iL],"",Form("%s #it{p}_{T} (GeV/#it{c})",lepton[iL].Data()),Form("%s #frac{dN}{d#it{p}_{T}}",labelNLeptons.Data()),1,1.1,
                            labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                            xPosLabel1D,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
        SaveCanvas(canvas, Form("%s/Comparison/%s_Pt.%s", outputDir.Data(),lepton[iL].Data(), suffix.Data()),kTRUE,kTRUE);

        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kTRUE,kFALSE,kFALSE,
                            vecPt[iL],"",Form("%s #it{p}_{T} (GeV/#it{c})",lepton[iL].Data()),Form("%s #frac{dN}{d#it{p}_{T}}",labelNLeptons.Data()),1,1.1,
                            labelData, colorCompare, kTRUE, 1.1, 1.1, kTRUE,
                            xPosLabel1D,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_%s_Pt.%s", outputDir.Data(),lepton[iL].Data(), suffix.Data()),kTRUE);
        DeleteVecTH1D(vecPt[iL]);
        //-----------------------------------
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                            vecNSigmadEdx[iL],"",Form("%s n#sigma dEdx",lepton[iL].Data()),Form("%s #frac{dN}{dn#sigma dEdx}",labelNLeptons.Data()),1,1.1,
                            labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                            xPosLabel1D,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
        SaveCanvas(canvas, Form("%s/Comparison/%s_NSigdEdx.%s", outputDir.Data(),lepton[iL].Data(), suffix.Data()), kFALSE, kTRUE);

        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                            vecNSigmadEdx[iL],"",Form("%s n#sigma dEdx",lepton[iL].Data()),Form("%s #frac{dN}{dn#sigma dEdx}",labelNLeptons.Data()),1,1.1,
                            labelData, colorCompare, kTRUE, 1.1, 1.1, kTRUE,
                            xPosLabel1D,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_%s_NSigdEdx.%s", outputDir.Data(),lepton[iL].Data(), suffix.Data()));
        DeleteVecTH1D(vecNSigmadEdx[iL]);
        //-----------------------------------
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                            vecITSAllclsR[iL],"","# ITS cl ",Form("%s #frac{dN}{dITS Cl}",labelNLeptons.Data()),1,1.1,
                            labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                            xPosLabel1D,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
        SaveCanvas(canvas, Form("%s/Comparison/%s_ITSCl.%s", outputDir.Data(),lepton[iL].Data(), suffix.Data()), kFALSE, kTRUE);

        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                            vecITSAllclsR[iL],"","# ITS cl ",Form("%s #frac{dN}{dITS Cl}",labelNLeptons.Data()),1,1.1,
                            labelData, colorCompare, kTRUE, 1.1, 1.1, kTRUE,
                            xPosLabel1D,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_%s_ITSCl.%s", outputDir.Data(),lepton[iL].Data(), suffix.Data()));
        DeleteVecTH1D(vecITSAllclsR[iL]);
        //-----------------------------------
        for(Int_t iB=0; iB<7; iB++){
            DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                vecITSclsR[iL][iB],"","R (cm)",Form("%s #frac{d#it{N}_{e^{%s}}}{d#it{R}}",labelNLeptons.Data(),charge[iL].Data()),1,1.1,
                                labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                xPosLabel1D,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
            SaveCanvas(canvas, Form("%s/Comparison/%s_%iITSClvsR.%s", outputDir.Data(),lepton[iL].Data(), iB, suffix.Data()), kFALSE, kTRUE);

            DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                vecITSclsR[iL][iB],"","R (cm)",Form("%s #frac{d#it{N}_{e^{%s}}}{d#it{R}}",labelNLeptons.Data(),charge[iL].Data()),1,1.1,
                                labelData, colorCompare, kTRUE, 1.1, 1.1, kTRUE,
                                xPosLabel1D,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0]);
            SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_%s_%iITSClvsR.%s", outputDir.Data(),lepton[iL].Data(), iB, suffix.Data()));
            DeleteVecTH1D(vecITSclsR[iL][iB]);
        }
        //-----------------------------------
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                            vecFCL[iL],"",Form("TPC Clusters/Findable Clusters e^{%s}",charge[iL].Data()),Form("%s #frac{dN}{dRatio Cl}",labelNLeptons.Data()),1,1.1,
                            labelData, colorCompare, kTRUE, 1.1, 1.1, kTRUE,
                            0.15,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0], 11);
        SaveCanvas(canvas, Form("%s/Comparison/%s_FindTPCClusters.%s", outputDir.Data(),lepton[iL].Data(), suffix.Data()));

        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                            vecFCL[iL],"",Form("TPC Clusters/Findable Clusters e^{%s}",charge[iL].Data()),Form("%s #frac{dN}{dRatio Cl}",labelNLeptons.Data()),1,1.1,
                            labelData, colorCompare, kTRUE, 1.1, 1.1, kTRUE,
                            xPosLabel1D,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0], 31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_%s_FindTPCClusters.%s", outputDir.Data(),lepton[iL].Data(), suffix.Data()));
        DeleteVecTH1D(vecFCL[iL]);
        //-----------------------------------
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                    vecTPCCL[iL],"",Form("TPC Clusters e^{%s}",charge[iL].Data()),Form("%s #frac{dN}{dTPC Cl}",labelNLeptons.Data()),1,1.1,
                                    labelData, colorCompare, kTRUE, 1.1, 1.1, kTRUE,
                                    0.15,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],11);
        SaveCanvas(canvas, Form("%s/Comparison/%s_TPCClusters.%s", outputDir.Data(),lepton[iL].Data(), suffix.Data()));

        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                        vecTPCCL[iL],"",Form("TPC Clusters e^{%s}",charge[iL].Data()),Form("%s #frac{dN}{dTPC Cl}",labelNLeptons.Data()),1,1.1,
                                        labelData, colorCompare, kTRUE, 1.1, 1.1, kTRUE,
                                        xPosLabel1D,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
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
