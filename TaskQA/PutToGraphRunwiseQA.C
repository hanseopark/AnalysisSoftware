/********************************************************************************************************
****** 		provided by Gamma Conversion Group,  													*****
******		Lucia Leardini , lucia.leardini@cern.ch													*****
********************************************************************************************************/
#include <Riostream.h>
#include <fstream>
#include "TMath.h"
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TPaveLabel.h>
#include <TSystem.h>
#include <TFrame.h>
#include <TStyle.h>
#include <TString.h>
#include "TGaxis.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TF1.h"
#include "TVirtualFitter.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TDatabasePDG.h"
#include "TMinuit.h"
#include "TLatex.h"
#include "TTree.h"
#include "TASImage.h"
#include "TPostScript.h"
#include "TGraphErrors.h"
#include "TArrow.h"
#include "TMarker.h"
#include "TGraphAsymmErrors.h"
#include "PutToGraphRunwiseQA.h"

void PutToGraphRunwiseQA(TString mode = "data o MC", TString dEdxFile = ".txt", TString GammaTPCSectorsFile = ".txt", TString GammaTPCSectorsErrorFile = ".txt", TString centrality = "", TString outputDir = "", TString suffix = "eps")
{
	gROOT->Reset();	
	gROOT->SetStyle("Plain");
	
	ifstream filedEdxFile;
 	filedEdxFile.open(dEdxFile.Data(),ios_base::in);
	cout << dEdxFile.Data() << endl;

	ifstream fileTPCSectors;
 	fileTPCSectors.open(GammaTPCSectorsFile.Data(),ios_base::in);
	cout << GammaTPCSectorsFile.Data() << endl;

	ifstream fileTPCSectorsError;
 	fileTPCSectorsError.open(GammaTPCSectorsErrorFile.Data(),ios_base::in);
	cout << GammaTPCSectorsErrorFile.Data() << endl;


   	Int_t run = 0;
	while(!filedEdxFile.eof()){
	   	ex[run] = 0;
		ey[run] = 0;

		filedEdxFile >> nRun[run] >> Events[run] >> meanNSigmadEdxElectron[run] >> meanErrNSigmadEdxElectron[run] >> meanNSigmadEdxPositron[run] >> meanErrNSigmadEdxPositron[run] >> 
		meanelNSigmaProtonCut[run] >> meanErrelNSigmaProtonCut[run] >> meanposNSigmaProtonCut[run] >> meanErrposNSigmaProtonCut[run] >> sigmaelNSigmaElectronCut[run] >> sigmaErrelNSigmaElectronCut[run] >> sigmaposNSigmaElectronCut[run] >> sigmaErrposNSigmaElectronCut[run] >> widthelNSigmaElectronCut[run] >> widthErrelNSigmaElectronCut[run] >> widthposNSigmaElectronCut[run] >> widthErrposNSigmaElectronCut[run] >> sigmaelNSigmaProtonCut[run] >> sigmaErrelNSigmaProtonCut[run] >> sigmaposNSigmaProtonCut[run] >> sigmaErrposNSigmaProtonCut[run] >> widthelNSigmaProtonCut[run] >> widthErrelNSigmaProtonCut[run] >> widthposNSigmaProtonCut[run] >> widthErrposNSigmaProtonCut[run] >> meanEtaElectron[run] >> meanErrEtaElectron[run] >> meanEtaPositron[run] >> meanErrEtaPositron[run] >> meanEtaNegElectron[run] >> meanErrEtaNegElectron[run] >> meanEtaNegPositron[run] >> meanErrEtaNegPositron[run] >> meanEtaPosElectron[run] >> meanErrEtaPosElectron[run] >> meanEtaPosPositron[run] >> meanErrEtaPosPositron[run] >> meanGammaEta[run] >> meanErrGammaEta[run] >> meanGammaEtaNeg[run] >> meanErrGammaEtaNeg[run] >> meanGammaEtaPos[run] >> meanErrGammaEtaPos[run] >> meanGammaPt[run] >> meanErrGammaPt[run];

		//TPC sector Nev only
		fileTPCSectors >> a[run] >> Events[run] >> NumGammaAllSectorEtaNeg[run] >> NumGammaAllSectorEtaPos[run] >> NumGammaSector0EtaNeg[run] >> NumGammaSector0EtaPos[run] >> NumGammaSector1EtaNeg[run] >> NumGammaSector1EtaPos[run] >> NumGammaSector2EtaNeg[run] >> NumGammaSector2EtaPos[run] >> NumGammaSector3EtaNeg[run] >> NumGammaSector3EtaPos[run] >> NumGammaSector4EtaNeg[run] >> NumGammaSector4EtaPos[run] >> NumGammaSector5EtaNeg[run] >> NumGammaSector5EtaPos[run] >> NumGammaSector6EtaNeg[run] >> NumGammaSector6EtaPos[run] >> NumGammaSector7EtaNeg[run] >> NumGammaSector7EtaPos[run] >> NumGammaSector8EtaNeg[run] >> NumGammaSector8EtaPos[run] >> NumGammaSector9EtaNeg[run] >> NumGammaSector9EtaPos[run] >> NumGammaSector10EtaNeg[run] >> NumGammaSector10EtaPos[run] >> NumGammaSector11EtaNeg[run] >> NumGammaSector11EtaPos[run] >> NumGammaSector12EtaNeg[run] >> NumGammaSector12EtaPos[run] >> NumGammaSector13EtaNeg[run] >> NumGammaSector13EtaPos[run] >> NumGammaSector14EtaNeg[run] >> NumGammaSector14EtaPos[run] >> NumGammaSector15EtaNeg[run] >> NumGammaSector15EtaPos[run] >> NumGammaSector16EtaNeg[run] >> NumGammaSector16EtaPos[run] >> NumGammaSector17EtaNeg[run] >> NumGammaSector17EtaPos[run];

		fileTPCSectorsError >> b[run] >> Events[run] >> ErrNumGammaAllSectorEtaNeg[run] >> ErrNumGammaAllSectorEtaPos[run] >> ErrNumGammaSector0EtaNeg[run] >> ErrNumGammaSector0EtaPos[run] >> ErrNumGammaSector1EtaNeg[run] >> ErrNumGammaSector1EtaPos[run] >> ErrNumGammaSector2EtaNeg[run] >> ErrNumGammaSector2EtaPos[run] >> ErrNumGammaSector3EtaNeg[run] >> ErrNumGammaSector3EtaPos[run] >> ErrNumGammaSector4EtaNeg[run] >> ErrNumGammaSector4EtaPos[run] >> ErrNumGammaSector5EtaNeg[run] >> ErrNumGammaSector5EtaPos[run] >> ErrNumGammaSector6EtaNeg[run] >> ErrNumGammaSector6EtaPos[run] >> ErrNumGammaSector7EtaNeg[run] >> ErrNumGammaSector7EtaPos[run] >> ErrNumGammaSector8EtaNeg[run] >> ErrNumGammaSector8EtaPos[run] >> ErrNumGammaSector9EtaNeg[run] >> ErrNumGammaSector9EtaPos[run] >> ErrNumGammaSector10EtaNeg[run] >> ErrNumGammaSector10EtaPos[run] >> ErrNumGammaSector11EtaNeg[run] >> ErrNumGammaSector11EtaPos[run] >> ErrNumGammaSector12EtaNeg[run] >> ErrNumGammaSector12EtaPos[run] >> ErrNumGammaSector13EtaNeg[run] >> ErrNumGammaSector13EtaPos[run] >> ErrNumGammaSector14EtaNeg[run] >> ErrNumGammaSector14EtaPos[run] >> ErrNumGammaSector15EtaNeg[run] >> ErrNumGammaSector15EtaPos[run] >> ErrNumGammaSector16EtaNeg[run] >> ErrNumGammaSector16EtaPos[run] >> ErrNumGammaSector17EtaNeg[run] >> ErrNumGammaSector17EtaPos[run];	
		run++;
	}
	cout << "here" << endl;
	cout << "=====================================================================================" << endl;
	filedEdxFile.close();
	fileTPCSectors.close();
	fileTPCSectorsError.close();
	run = run-1;

	
	//Draw variables run wise
	TGraphErrors* EventsPerRun = new TGraphErrors(run,nRun,Events,ex,ey);	

	TGraphErrors* NSigmadEdxElectron = new TGraphErrors(run,nRun,meanNSigmadEdxElectron,ex,meanErrNSigmadEdxElectron);
	TGraphErrors* NSigmadEdxPositron = new TGraphErrors(run,nRun,meanNSigmadEdxPositron,ex,meanErrNSigmadEdxPositron);
	TGraphErrors* sigmaNSigmadEdxElectron = new TGraphErrors(run,nRun,sigmaelNSigmaElectronCut,ex,sigmaErrelNSigmaElectronCut);
	TGraphErrors* sigmaNSigmadEdxPositron = new TGraphErrors(run,nRun,sigmaposNSigmaElectronCut,ex,sigmaErrposNSigmaElectronCut);
	TGraphErrors* widthNSigmadEdxElectron = new TGraphErrors(run,nRun,widthelNSigmaElectronCut,ex,widthErrelNSigmaElectronCut);
	TGraphErrors* widthNSigmadEdxPositron = new TGraphErrors(run,nRun,widthposNSigmaElectronCut,ex,widthErrposNSigmaElectronCut);

	TGraphErrors* NSigmadEdxElectronProtonCut = new TGraphErrors(run,nRun,meanelNSigmaProtonCut,ex,meanErrelNSigmaProtonCut);
	TGraphErrors* NSigmadEdxPositronProtonCut = new TGraphErrors(run,nRun,meanposNSigmaProtonCut,ex,meanErrposNSigmaProtonCut);
	TGraphErrors* sigmaNSigmadEdxElectronProtonCut = new TGraphErrors(run,nRun,sigmaelNSigmaProtonCut,ex,sigmaErrelNSigmaProtonCut);
	TGraphErrors* sigmaNSigmadEdxPositronProtonCut = new TGraphErrors(run,nRun,sigmaposNSigmaProtonCut,ex,sigmaErrposNSigmaProtonCut);
	TGraphErrors* widthNSigmadEdxElectronProtonCut = new TGraphErrors(run,nRun,widthelNSigmaProtonCut,ex,widthErrelNSigmaProtonCut);
	TGraphErrors* widthNSigmadEdxPositronProtonCut = new TGraphErrors(run,nRun,widthposNSigmaProtonCut,ex,widthErrposNSigmaProtonCut);

	//eta
	TGraphErrors* ElectronEta = new TGraphErrors(run,nRun,meanEtaElectron,ex,meanErrEtaElectron);
	TGraphErrors* PositronEta = new TGraphErrors(run,nRun,meanEtaPositron,ex,meanErrEtaPositron);
	TGraphErrors* ElectronEtaNeg = new TGraphErrors(run,nRun,meanEtaNegElectron,ex,meanErrEtaNegElectron);
	TGraphErrors* PositronEtaNeg = new TGraphErrors(run,nRun,meanEtaNegPositron,ex,meanErrEtaNegPositron);
	TGraphErrors* ElectronEtaPos = new TGraphErrors(run,nRun,meanEtaPosElectron,ex,meanErrEtaPosElectron);
	TGraphErrors* PositronEtaPos = new TGraphErrors(run,nRun,meanEtaPosPositron,ex,meanErrEtaPosPositron);

	TGraphErrors* PhotonPt = new TGraphErrors(run,nRun,meanGammaPt,ex,meanErrGammaPt);
	
	TGraphErrors* PhotonEta = new TGraphErrors(run,nRun,meanGammaEta,ex,meanErrGammaEta);
	TGraphErrors* PhotonEtaNeg = new TGraphErrors(run,nRun,meanGammaEtaNeg,ex,meanErrGammaEtaNeg);
	TGraphErrors* PhotonEtaPos = new TGraphErrors(run,nRun,meanGammaEtaPos,ex,meanErrGammaEtaPos);

	TGraphErrors* PhotonAllSectorEtaNeg = new TGraphErrors(run,nRun,NumGammaAllSectorEtaNeg,ex,ErrNumGammaAllSectorEtaNeg);
	TGraphErrors* PhotonAllSectorEtaPos = new TGraphErrors(run,nRun,NumGammaAllSectorEtaPos,ex,ErrNumGammaAllSectorEtaPos);
	TGraphErrors* PhotonSector0EtaNeg = new TGraphErrors(run,nRun,NumGammaSector0EtaNeg,ex,ErrNumGammaSector0EtaNeg);
	TGraphErrors* PhotonSector0EtaPos = new TGraphErrors(run,nRun,NumGammaSector0EtaPos,ex,ErrNumGammaSector0EtaPos);
	TGraphErrors* PhotonSector1EtaNeg = new TGraphErrors(run,nRun,NumGammaSector1EtaNeg,ex,ErrNumGammaSector1EtaNeg);
	TGraphErrors* PhotonSector1EtaPos = new TGraphErrors(run,nRun,NumGammaSector1EtaPos,ex,ErrNumGammaSector1EtaPos);
	TGraphErrors* PhotonSector2EtaNeg = new TGraphErrors(run,nRun,NumGammaSector2EtaNeg,ex,ErrNumGammaSector2EtaNeg);
	TGraphErrors* PhotonSector2EtaPos = new TGraphErrors(run,nRun,NumGammaSector2EtaPos,ex,ErrNumGammaSector2EtaPos);
	TGraphErrors* PhotonSector3EtaNeg = new TGraphErrors(run,nRun,NumGammaSector3EtaNeg,ex,ErrNumGammaSector3EtaNeg);
	TGraphErrors* PhotonSector3EtaPos = new TGraphErrors(run,nRun,NumGammaSector3EtaPos,ex,ErrNumGammaSector3EtaPos);
	TGraphErrors* PhotonSector4EtaNeg = new TGraphErrors(run,nRun,NumGammaSector4EtaNeg,ex,ErrNumGammaSector4EtaNeg);
	TGraphErrors* PhotonSector4EtaPos = new TGraphErrors(run,nRun,NumGammaSector4EtaPos,ex,ErrNumGammaSector4EtaPos);
	TGraphErrors* PhotonSector5EtaNeg = new TGraphErrors(run,nRun,NumGammaSector5EtaNeg,ex,ErrNumGammaSector5EtaNeg);
	TGraphErrors* PhotonSector5EtaPos = new TGraphErrors(run,nRun,NumGammaSector5EtaPos,ex,ErrNumGammaSector5EtaPos);
	TGraphErrors* PhotonSector6EtaNeg = new TGraphErrors(run,nRun,NumGammaSector6EtaNeg,ex,ErrNumGammaSector6EtaNeg);
	TGraphErrors* PhotonSector6EtaPos = new TGraphErrors(run,nRun,NumGammaSector6EtaPos,ex,ErrNumGammaSector6EtaPos);
	TGraphErrors* PhotonSector7EtaNeg = new TGraphErrors(run,nRun,NumGammaSector7EtaNeg,ex,ErrNumGammaSector7EtaNeg);
	TGraphErrors* PhotonSector7EtaPos = new TGraphErrors(run,nRun,NumGammaSector7EtaPos,ex,ErrNumGammaSector7EtaPos);
	TGraphErrors* PhotonSector8EtaNeg = new TGraphErrors(run,nRun,NumGammaSector8EtaNeg,ex,ErrNumGammaSector8EtaNeg);
	TGraphErrors* PhotonSector8EtaPos = new TGraphErrors(run,nRun,NumGammaSector8EtaPos,ex,ErrNumGammaSector8EtaPos);
	TGraphErrors* PhotonSector9EtaNeg = new TGraphErrors(run,nRun,NumGammaSector9EtaNeg,ex,ErrNumGammaSector9EtaNeg);
	TGraphErrors* PhotonSector9EtaPos = new TGraphErrors(run,nRun,NumGammaSector9EtaPos,ex,ErrNumGammaSector9EtaPos);
	TGraphErrors* PhotonSector10EtaNeg = new TGraphErrors(run,nRun,NumGammaSector10EtaNeg,ex,ErrNumGammaSector10EtaNeg);
	TGraphErrors* PhotonSector10EtaPos = new TGraphErrors(run,nRun,NumGammaSector10EtaPos,ex,ErrNumGammaSector10EtaPos);
	TGraphErrors* PhotonSector11EtaNeg = new TGraphErrors(run,nRun,NumGammaSector11EtaNeg,ex,ErrNumGammaSector11EtaNeg);
	TGraphErrors* PhotonSector11EtaPos = new TGraphErrors(run,nRun,NumGammaSector11EtaPos,ex,ErrNumGammaSector11EtaPos);
	TGraphErrors* PhotonSector12EtaNeg = new TGraphErrors(run,nRun,NumGammaSector12EtaNeg,ex,ErrNumGammaSector12EtaNeg);
	TGraphErrors* PhotonSector12EtaPos = new TGraphErrors(run,nRun,NumGammaSector12EtaPos,ex,ErrNumGammaSector12EtaPos);
	TGraphErrors* PhotonSector13EtaNeg = new TGraphErrors(run,nRun,NumGammaSector13EtaNeg,ex,ErrNumGammaSector13EtaNeg);
	TGraphErrors* PhotonSector13EtaPos = new TGraphErrors(run,nRun,NumGammaSector13EtaPos,ex,ErrNumGammaSector13EtaPos);
	TGraphErrors* PhotonSector14EtaNeg = new TGraphErrors(run,nRun,NumGammaSector14EtaNeg,ex,ErrNumGammaSector14EtaNeg);
	TGraphErrors* PhotonSector14EtaPos = new TGraphErrors(run,nRun,NumGammaSector14EtaPos,ex,ErrNumGammaSector14EtaPos);
	TGraphErrors* PhotonSector15EtaNeg = new TGraphErrors(run,nRun,NumGammaSector15EtaNeg,ex,ErrNumGammaSector15EtaNeg);
	TGraphErrors* PhotonSector15EtaPos = new TGraphErrors(run,nRun,NumGammaSector15EtaPos,ex,ErrNumGammaSector15EtaPos);
	TGraphErrors* PhotonSector16EtaNeg = new TGraphErrors(run,nRun,NumGammaSector16EtaNeg,ex,ErrNumGammaSector16EtaNeg);
	TGraphErrors* PhotonSector16EtaPos = new TGraphErrors(run,nRun,NumGammaSector16EtaPos,ex,ErrNumGammaSector16EtaPos);
	TGraphErrors* PhotonSector17EtaNeg = new TGraphErrors(run,nRun,NumGammaSector17EtaNeg,ex,ErrNumGammaSector17EtaNeg);
	TGraphErrors* PhotonSector17EtaPos = new TGraphErrors(run,nRun,NumGammaSector17EtaPos,ex,ErrNumGammaSector17EtaPos);

	
	cout << "Saving to File" << endl;
	TFile *fileoutput = new TFile(Form("%s/OutputAnalysisQA_%s.root",outputDir.Data(), mode.Data()),"UPDATE");
	fileoutput->cd();

	EventsPerRun->Write(Form("EventsPerRun_%s",centrality.Data()));

	NSigmadEdxElectron->Write(Form("NSigmadEdxElectron_%s",centrality.Data()));
	NSigmadEdxPositron->Write(Form("NSigmadEdxPositron_%s",centrality.Data()));

	NSigmadEdxElectron->Write(Form("meanNSigmadEdxElectron_%s",centrality.Data()));
	NSigmadEdxPositron->Write(Form("meanNSigmadEdxPositron_%s",centrality.Data()));
	sigmaNSigmadEdxElectron->Write(Form("sigmaelNSigmaElectronCut_%s",centrality.Data()));
	sigmaNSigmadEdxPositron->Write(Form("sigmaposNSigmaElectronCut_%s",centrality.Data()));
	widthNSigmadEdxElectron->Write(Form("widthelNSigmaElectronCut_%s",centrality.Data()));
	widthNSigmadEdxPositron->Write(Form("widthposNSigmaElectronCut_%s",centrality.Data()));

	NSigmadEdxElectronProtonCut->Write(Form("meanelNSigmaProtonCut_%s",centrality.Data()));
	NSigmadEdxPositronProtonCut->Write(Form("meanposNSigmaProtonCut_%s",centrality.Data()));
	sigmaNSigmadEdxElectronProtonCut->Write(Form("sigmaelNSigmaProtonCut_%s",centrality.Data()));
	sigmaNSigmadEdxPositronProtonCut->Write(Form("sigmaposNSigmaProtonCut_%s",centrality.Data()));
	widthNSigmadEdxElectronProtonCut->Write(Form("widthelNSigmaProtonCut_%s",centrality.Data()));
	widthNSigmadEdxPositronProtonCut->Write(Form("widthposNSigmaProtonCut_%s",centrality.Data()));


	ElectronEta->Write(Form("ElectronEta_%s",centrality.Data()));
	PositronEta->Write(Form("PositronEta_%s",centrality.Data()));
	PhotonEta->Write(Form("PhotonEta_%s",centrality.Data()));
	ElectronEtaNeg->Write(Form("ElectronEtaNeg_%s",centrality.Data()));
	PositronEtaNeg->Write(Form("PositronEtaNeg_%s",centrality.Data()));
	PhotonEtaNeg->Write(Form("PhotonEtaNeg_%s",centrality.Data()));
	ElectronEtaPos->Write(Form("ElectronEtaPos_%s",centrality.Data()));
	PositronEtaPos->Write(Form("PositronEtaPos_%s",centrality.Data()));
	PhotonEtaPos->Write(Form("PhotonEtaPos_%s",centrality.Data()));
	PhotonPt->Write(Form("PhotonPt_%s",centrality.Data()));

	PhotonAllSectorEtaNeg->Write(Form("PhotonAllSectorsn_%s",centrality.Data()));
	PhotonAllSectorEtaPos->Write(Form("PhotonAllSectorsp_%s",centrality.Data()));
	PhotonSector0EtaNeg->Write(Form("PhotonSector0n_%s",centrality.Data()));
	PhotonSector0EtaPos->Write(Form("PhotonSector0p_%s",centrality.Data()));
	PhotonSector1EtaNeg->Write(Form("PhotonSector1n_%s",centrality.Data()));
	PhotonSector1EtaPos->Write(Form("PhotonSector1p_%s",centrality.Data()));
	PhotonSector2EtaNeg->Write(Form("PhotonSector2n_%s",centrality.Data()));
	PhotonSector2EtaPos->Write(Form("PhotonSector2p_%s",centrality.Data()));
	PhotonSector3EtaNeg->Write(Form("PhotonSector3n_%s",centrality.Data()));
	PhotonSector3EtaPos->Write(Form("PhotonSector3p_%s",centrality.Data()));
	PhotonSector4EtaNeg->Write(Form("PhotonSector4n_%s",centrality.Data()));
	PhotonSector4EtaPos->Write(Form("PhotonSector4p_%s",centrality.Data()));
	PhotonSector5EtaNeg->Write(Form("PhotonSector5n_%s",centrality.Data()));
	PhotonSector5EtaPos->Write(Form("PhotonSector5p_%s",centrality.Data()));
	PhotonSector6EtaNeg->Write(Form("PhotonSector6n_%s",centrality.Data()));
	PhotonSector6EtaPos->Write(Form("PhotonSector6p_%s",centrality.Data()));
	PhotonSector7EtaNeg->Write(Form("PhotonSector7n_%s",centrality.Data()));
	PhotonSector7EtaPos->Write(Form("PhotonSector7p_%s",centrality.Data()));
	PhotonSector8EtaNeg->Write(Form("PhotonSector8n_%s",centrality.Data()));
	PhotonSector8EtaPos->Write(Form("PhotonSector8p_%s",centrality.Data()));
	PhotonSector9EtaNeg->Write(Form("PhotonSector9n_%s",centrality.Data()));
	PhotonSector9EtaPos->Write(Form("PhotonSector9p_%s",centrality.Data()));
	PhotonSector10EtaNeg->Write(Form("PhotonSector10n_%s",centrality.Data()));
	PhotonSector10EtaPos->Write(Form("PhotonSector10p_%s",centrality.Data()));
	PhotonSector11EtaNeg->Write(Form("PhotonSector11n_%s",centrality.Data()));
	PhotonSector11EtaPos->Write(Form("PhotonSector11p_%s",centrality.Data()));
	PhotonSector12EtaNeg->Write(Form("PhotonSector12n_%s",centrality.Data()));
	PhotonSector12EtaPos->Write(Form("PhotonSector12p_%s",centrality.Data()));
	PhotonSector13EtaNeg->Write(Form("PhotonSector13n_%s",centrality.Data()));
	PhotonSector13EtaPos->Write(Form("PhotonSector13p_%s",centrality.Data()));
	PhotonSector14EtaNeg->Write(Form("PhotonSector14n_%s",centrality.Data()));
	PhotonSector14EtaPos->Write(Form("PhotonSector14p_%s",centrality.Data()));
	PhotonSector15EtaNeg->Write(Form("PhotonSector15n_%s",centrality.Data()));
	PhotonSector15EtaPos->Write(Form("PhotonSector15p_%s",centrality.Data()));
	PhotonSector16EtaNeg->Write(Form("PhotonSector16n_%s",centrality.Data()));
	PhotonSector16EtaPos->Write(Form("PhotonSector16p_%s",centrality.Data()));
	PhotonSector17EtaNeg->Write(Form("PhotonSector17n_%s",centrality.Data()));
	PhotonSector17EtaPos->Write(Form("PhotonSector17p_%s",centrality.Data()));

	fileoutput->Write();
	fileoutput->Close();	


}



