/********************************************************************************************************
****** 		provided by Gamma Conversion Group,  													*****
******		Lucia Leardini , lucia.leardini@cern.ch													*****
********************************************************************************************************/
#include <Riostream.h>
#include <iostream>
#include "TMath.h"
#include <stdio.h>
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
#include "TASImage.h"
#include "TPostScript.h"
#include "TGraphErrors.h"
#include "TArrow.h"
#include "TMarker.h"
#include "TGraphAsymmErrors.h" 
#include "/home/admin1/leardini/photonconv/AnalysisSoftware/CommonHeaders/PlottingGammaConversionHistos.h"
#include "/home/admin1/leardini/photonconv/AnalysisSoftware/CommonHeaders/PlottingGammaConversionAdditional.h"
#include "/home/admin1/leardini/photonconv/AnalysisSoftware/CommonHeaders/FittingGammaConversion.h"
#include "/home/admin1/leardini/photonconv/AnalysisSoftware/CommonHeaders/ConversionFunctions.h"
#include "TVectorD.h"
#include "TTreeStream.h"
#include "ExtractRunwiseQA.h"


using namespace std;

void ExtractRunwiseQA(TString mode="Data",TString fileName = ".root",TString cutNumber="",TString cutNumberQA="",TString fileDirectory="",TString outputDirectory="",Int_t run=167902,TString suffix= "pdf", TString cent = "010")
{
	fstream EventFileQA(Form("%s/EventsRunwise_%s_%s.txt",outputDirectory.Data(),cent.Data(),mode.Data()),ios::out|ios::app);
	if(!EventFileQA.is_open()){
	   cout<<"Problem opening file events"<<endl;
	   return;
	}
 
	fstream outputFile(Form("%s/dEdxRunwise_%s_%s.txt",outputDirectory.Data(),cent.Data(),mode.Data()),ios::out|ios::app);
	if(!outputFile.is_open()){
	   cout<<"Problem opening file 1"<<endl;
	   return;
	}
	fstream outputFile2(Form("%s/GammaTPCSectors_%s_%s.txt",outputDirectory.Data(),cent.Data(), mode.Data()),ios::out|ios::app);
	if(!outputFile2.is_open()){
	   cout<<"Problem opening file 2"<<endl;
	   return;
	}
	fstream outputFileErr2(Form("%s/GammaTPCSectorsErrors_%s_%s.txt",outputDirectory.Data(),cent.Data(), mode.Data()),ios::out|ios::app);
	if(!outputFileErr2.is_open()){
	   cout<<"Problem opening file 2"<<endl;
	   return;
	}

	fstream MissingFiles(Form("%s/MissingFiles_%s_%s.txt",outputDirectory.Data(),cent.Data(),mode.Data()),ios::out|ios::app);
	if(!MissingFiles.is_open()){
	   cout<<"Problem opening file events"<<endl;
	   return;
	}
	
	cout << run << endl;
	SetPlotStyle();
	
	
 	File = new TFile(Form("%s/%d/%s", fileDirectory.Data(), run, fileName.Data()));
	if (File->IsZombie()) { 
		MissingFiles << run << "\t" << Form("%s/%d/%s", fileDirectory.Data(), run, fileName.Data()) << endl; 
		return;
	} else cout << "File is there" << endl;

	folderInFile = (TDirectoryFile*) File->Get(Form("GammaConv_%s",cutNumberQA.Data()));
	//get the number of events	
	VertexZ = (TH1F*) folderInFile->Get("histoVertexZ");
	Events = VertexZ->GetEntries();
	cout << "Number of events for run number " << run << ": " << Events << endl;

	//get the histograms		
	NSigmadEdxElectron = (TH3F*) folderInFile->Get("histoElectronNSigmadEdxEtaP")->Clone("ElectronNSigmadEdxEtaP");
	NSigmadEdxPositron = (TH3F*) folderInFile->Get("histoPositronNSigmadEdxEtaP")->Clone("PositronNSigmadEdxEtaP");          
	GammaEtaPt = (TH2F*) folderInFile->Get("histoGammaEtaPt")->Clone("GammaEtaPt");

  	//projection of NSigmadEdx
	projxNSigmadEdxElectron = (TH2D*) NSigmadEdxElectron->Project3D("xz");
	projxNSigmadEdxPositron = (TH2D*) NSigmadEdxPositron->Project3D("xz");

	//cut below kaons line ("only" electrons) //11 invece che 12 per la proj.
	elNSigmaElectronCut = (TH1D*) projxNSigmadEdxElectron->ProjectionY("elNSigmaElectronCut", 1,13,"e"); 
	meanelNSigmaElectronCut = elNSigmaElectronCut->GetMean(1);
	meanErrelNSigmaElectronCut = elNSigmaElectronCut->GetRMSError(1);
	posNSigmaElectronCut = (TH1D*) projxNSigmadEdxPositron->ProjectionY("posNSigmaElectronCut", 1,13, "e");
	meanposNSigmaElectronCut = posNSigmaElectronCut->GetMean(1);
	meanErrposNSigmaElectronCut = posNSigmaElectronCut->GetRMSError(1);
	
	//get sigma
	sigmaelNSigmaElectronCut = elNSigmaElectronCut->GetStdDev();
	sigmaErrelNSigmaElectronCut = elNSigmaElectronCut->GetRMSError(1);
	sigmaposNSigmaElectronCut = posNSigmaElectronCut->GetStdDev();
	sigmaErrposNSigmaElectronCut = posNSigmaElectronCut->GetRMSError(1);
	
	//get width
	bin1el = elNSigmaElectronCut->FindFirstBinAbove(elNSigmaElectronCut->GetMaximum()/2);
	bin1elErr = elNSigmaElectronCut->GetBinError(bin1el);
   	bin2el = elNSigmaElectronCut->FindLastBinAbove(elNSigmaElectronCut->GetMaximum()/2);
	bin2elErr = elNSigmaElectronCut->GetBinError(bin2el);
	widthelNSigmaElectronCut = elNSigmaElectronCut->GetBinCenter(bin2el) - elNSigmaElectronCut->GetBinCenter(bin1el);
	widthErrelNSigmaElectronCut = elNSigmaElectronCut->GetRMSError(1); //bin1elErr + bin2elErr;
	bin1pos = posNSigmaElectronCut->FindFirstBinAbove(posNSigmaElectronCut->GetMaximum()/2);
	bin1posErr = posNSigmaElectronCut->GetBinError(bin1pos);
   	bin2pos = posNSigmaElectronCut->FindLastBinAbove(posNSigmaElectronCut->GetMaximum()/2);
	bin2posErr = posNSigmaElectronCut->GetBinError(bin2pos);
	widthposNSigmaElectronCut = posNSigmaElectronCut->GetBinCenter(bin2pos) - posNSigmaElectronCut->GetBinCenter(bin1pos);
	widthErrposNSigmaElectronCut =  posNSigmaElectronCut->GetRMSError(1); //bin1posErr + bin2posErr;
	
	//which is also used for the final plotting -> avoid pions contamination
	meanNSigmadEdxElectron = meanelNSigmaElectronCut;
	meanErrNSigmadEdxElectron = meanErrelNSigmaElectronCut;
	meanNSigmadEdxPositron = meanposNSigmaElectronCut;
	meanErrNSigmadEdxPositron = meanErrposNSigmaElectronCut;

	//cut ABOVE proton line
	elNSigmaProtonCut = (TH1D*) projxNSigmadEdxElectron->ProjectionY("elNSigmaProtonCut", 50,99,"e"); 
	meanelNSigmaProtonCut = elNSigmaProtonCut->GetMean(1);
	meanErrelNSigmaProtonCut = elNSigmaProtonCut->GetRMSError(1);
	posNSigmaProtonCut = (TH1D*) projxNSigmadEdxPositron->ProjectionY("posNSigmaProtonCut", 50,99, "e");
	meanposNSigmaProtonCut = posNSigmaProtonCut->GetMean(1);
	meanErrposNSigmaProtonCut = posNSigmaProtonCut->GetRMSError(1);
	
	//get sigma
	sigmaelNSigmaProtonCut = elNSigmaProtonCut->GetStdDev();
	sigmaErrelNSigmaProtonCut = elNSigmaProtonCut->GetRMSError(1);
	sigmaposNSigmaProtonCut = posNSigmaProtonCut->GetStdDev();
	sigmaErrposNSigmaProtonCut = posNSigmaProtonCut->GetRMSError(1);

	//get width
	bin3el = elNSigmaProtonCut->FindFirstBinAbove(elNSigmaProtonCut->GetMaximum()/2);
	bin3elErr = elNSigmaProtonCut->GetBinError(bin3el);
   	bin4el = elNSigmaProtonCut->FindLastBinAbove(elNSigmaProtonCut->GetMaximum()/2);
	bin4elErr = elNSigmaProtonCut->GetBinError(bin4el);
	widthelNSigmaProtonCut = elNSigmaProtonCut->GetBinCenter(bin4el) - elNSigmaProtonCut->GetBinCenter(bin3el);
	widthErrelNSigmaProtonCut = elNSigmaProtonCut->GetRMSError(1); //bin3elErr + bin4elErr;
	bin3pos = posNSigmaProtonCut->FindFirstBinAbove(posNSigmaProtonCut->GetMaximum()/2);
	bin3posErr = posNSigmaProtonCut->GetBinError(bin3pos);
   	bin4pos = posNSigmaProtonCut->FindLastBinAbove(posNSigmaProtonCut->GetMaximum()/2);
	bin4posErr = posNSigmaProtonCut->GetBinError(bin4pos);
	widthposNSigmaProtonCut = posNSigmaProtonCut->GetBinCenter(bin4pos) - posNSigmaProtonCut->GetBinCenter(bin3pos);
	widthErrposNSigmaProtonCut = posNSigmaProtonCut->GetRMSError(1); // = bin3posErr + bin4posErr;

	//Eta total 
	projyNSigmadEdxElectron = (TH1D*) NSigmadEdxElectron->Project3D("y");
	meanEtaElectron = projyNSigmadEdxElectron->GetMean(1);
	meanErrEtaElectron = projyNSigmadEdxElectron->GetRMSError(1);
	projyNSigmadEdxPositron = (TH1D*) NSigmadEdxPositron->Project3D("y");
	meanEtaPositron = projyNSigmadEdxPositron->GetMean(1);
	meanErrEtaPositron = projyNSigmadEdxPositron->GetRMSError(1);
	GammaEta = (TH1D*) GammaEtaPt->ProjectionX("GammaEta");
	meanGammaEta = GammaEta->GetMean(1);
	meanErrGammaEta = GammaEta->GetRMSError(1);

	//Eta negative 
	projyNSigmadEdxElectron->GetXaxis()->SetRangeUser(-0.9,0.);
	projyNSigmadEdxPositron->GetXaxis()->SetRangeUser(-0.9,0.);
	GammaEta->GetXaxis()->SetRangeUser(-0.9,0.);
	meanEtaNegElectron = projyNSigmadEdxElectron->GetMean(1);
	meanErrEtaNegElectron = projyNSigmadEdxElectron->GetRMSError(1);
	meanEtaNegPositron = projyNSigmadEdxPositron->GetMean(1);
	meanErrEtaNegPositron = projyNSigmadEdxPositron->GetRMSError(1);

	meanGammaEtaNeg = GammaEta->GetMean(1);
	meanErrGammaEtaNeg = GammaEta->GetRMSError(1);
	//Eta positive 
	projyNSigmadEdxElectron->GetXaxis()->SetRangeUser(0.,0.9);
	projyNSigmadEdxPositron->GetXaxis()->SetRangeUser(0.,0.9);
	GammaEta->GetXaxis()->SetRangeUser(0.,0.9);
	meanEtaPosElectron = projyNSigmadEdxElectron->GetMean(1);
	meanErrEtaPosElectron = projyNSigmadEdxElectron->GetRMSError(1);
	meanEtaPosPositron = projyNSigmadEdxPositron->GetMean(1);
	meanErrEtaPosPositron = projyNSigmadEdxPositron->GetRMSError(1);
	meanGammaEtaPos = GammaEta->GetMean(1);
	meanErrGammaEtaPos = GammaEta->GetRMSError(1);

	//Gamma Pt
	GammaPt = (TH1D*) GammaEtaPt->ProjectionY("GammaPt");
	meanGammaPt = GammaPt->GetMean(1);
	meanErrGammaPt = GammaPt->GetRMSError(1);

	//Gamma Phi
	TH1F *GammaPhiEtaNeg = (TH1F*) folderInFile->Get("histoGammaPhiEtaNeg")->Clone("GammaPhiEtaNeg");
	TH1F *GammaPhiEtaPos = (TH1F*) folderInFile->Get("histoGammaPhiEtaPos")->Clone("GammaPhiEtaPos");
	Double_t pi = TMath::Pi();
	Double_t rangePhi = pi/9.0;
	Int_t nBinsEtaNeg = GammaPhiEtaNeg->GetXaxis()->GetNbins();
	Int_t nBinsEtaPos = GammaPhiEtaPos->GetXaxis()->GetNbins();
	GammaPhiEtaNeg->Sumw2();	
	GammaPhiEtaPos->Sumw2();
	GammaPhiEtaNeg->Scale(1./Events);	
	GammaPhiEtaPos->Scale(1./Events);

	//Divided only
// 	Double_t ErrNumGammaAllSectorEtaPos;
// 	Double_t ErrNumGammaAllSectorEtaNeg;
	Double_t NumGammaAllSectorEtaNeg = GammaPhiEtaNeg->IntegralAndError(1,nBinsEtaNeg,ErrNumGammaAllSectorEtaNeg);
	Double_t NumGammaAllSectorEtaPos = GammaPhiEtaPos->IntegralAndError(1,nBinsEtaPos,ErrNumGammaAllSectorEtaPos);

	//Events File
	EventFileQA << "&" << run << "\t &" << Events << endl;

	//analysis QA
	outputFile << run << "\t" << Events << "\t" << meanNSigmadEdxElectron << "\t" << meanErrNSigmadEdxElectron << "\t" << meanNSigmadEdxPositron << "\t" << meanErrNSigmadEdxPositron << "\t" << 
meanelNSigmaProtonCut << "\t" << meanErrelNSigmaProtonCut << "\t" << meanposNSigmaProtonCut << "\t" << meanErrposNSigmaProtonCut << "\t" << sigmaelNSigmaElectronCut << "\t" << sigmaErrelNSigmaElectronCut  << "\t" << sigmaposNSigmaElectronCut  << "\t" << sigmaErrposNSigmaElectronCut  << "\t" << widthelNSigmaElectronCut << "\t"  << widthErrelNSigmaElectronCut << "\t" << widthposNSigmaElectronCut << "\t" << widthErrposNSigmaElectronCut << "\t" << sigmaelNSigmaProtonCut << "\t" << sigmaErrelNSigmaProtonCut  << "\t" << sigmaposNSigmaProtonCut  << "\t" << sigmaErrposNSigmaProtonCut  << "\t" << widthelNSigmaProtonCut << "\t" << widthErrelNSigmaProtonCut << "\t" << widthposNSigmaProtonCut << "\t" << widthErrposNSigmaProtonCut << "\t" << meanEtaElectron << "\t" << meanErrEtaElectron << "\t" << meanEtaPositron << "\t" << meanErrEtaPositron << "\t" << meanEtaNegElectron << "\t" << meanErrEtaNegElectron << "\t" << meanEtaNegPositron << "\t" << meanErrEtaNegPositron << "\t" << meanEtaPosElectron << "\t" << meanErrEtaPosElectron << "\t" << meanEtaPosPositron << "\t" << meanErrEtaPosPositron << "\t" << meanGammaEta << "\t" << meanErrGammaEta << "\t" << meanGammaEtaNeg << "\t" << meanErrGammaEtaNeg << "\t" << meanGammaEtaPos << "\t" << meanErrGammaEtaPos << "\t" << meanGammaPt << "\t" << meanErrGammaPt << endl;
 	outputFile.close();
	
	Double_t NumGammaSectorEtaPos[18];
	Double_t NumGammaSectorEtaNeg[18];
	Double_t NumGammaSectorEtaPosError[18];
	Double_t NumGammaSectorEtaNegError[18];
	
// 	TVectorD NumGammaSectorEtaPosV(18);
// 	TVectorD NumGammaSectorEtaNegV(18);
// 	TVectorD NumGammaSectorEtaPosErrorV(18);
// 	TVectorD NumGammaSectorEtaNegErrorV(18);

	outputFile2 << run << "\t" << Events << "\t" << NumGammaAllSectorEtaNeg << "\t" << NumGammaAllSectorEtaPos << "\t";
	outputFileErr2 << run << "\t" << Events << "\t" << ErrNumGammaAllSectorEtaNeg << "\t" << ErrNumGammaAllSectorEtaPos << "\t";
	
	//loop on phi sectors
	for(Int_t i = 0; i <18; i++){
	
		//filling normal double array
		NumGammaSectorEtaPos[i] = GammaPhiEtaNeg->IntegralAndError(GammaPhiEtaNeg->GetXaxis()->FindBin((rangePhi*i)+0.001), GammaPhiEtaNeg->GetXaxis()->FindBin((rangePhi*(i+1))-0.0001),NumGammaSectorEtaPosError[i]);
		NumGammaSectorEtaNeg[i] = GammaPhiEtaPos->IntegralAndError(GammaPhiEtaPos->GetXaxis()->FindBin((rangePhi*i)+0.001), GammaPhiEtaPos->GetXaxis()->FindBin((rangePhi*(i+1))-0.0001),NumGammaSectorEtaNegError[i]);

		//TPC sector Nev only
		outputFile2 << NumGammaSectorEtaNeg[i] << "\t" << NumGammaSectorEtaPos[i] << endl;
		outputFileErr2 << NumGammaSectorEtaNegError[i] << "\t" << NumGammaSectorEtaPosError[i] << endl;
		
		//filling TVector
// 		NumGammaSectorEtaPosV[i] = 0;
// 		NumGammaSectorEtaNegV[i] = 0;
// 		NumGammaSectorEtaPosErrorV[i] = 0;
// 		NumGammaSectorEtaPosErrorV[i] = 0;
// 		
// 		NumGammaSectorEtaPosV[i] = GammaPhiEtaNeg->IntegralAndError(GammaPhiEtaNeg->GetXaxis()->FindBin((rangePhi*i)+0.001), GammaPhiEtaNeg->GetXaxis()->FindBin((rangePhi*(i+1))-0.0001),NumGammaSectorEtaPosErrorV[i]);
// 		NumGammaSectorEtaNegV[i] = GammaPhiEtaPos->IntegralAndError(GammaPhiEtaPos->GetXaxis()->FindBin((rangePhi*i)+0.001), GammaPhiEtaPos->GetXaxis()->FindBin((rangePhi*(i+1))-0.0001),NumGammaSectorEtaPosErrorV[i]);

	}
	
	outputFile2.close();
	outputFileErr2.close();
	
	//saving variables in trendign file for general QA (for Marian)
// 	TString outfile = Form("%s/trending_%s.root",outputDirectory.Data(),cent.Data());
// 	TTreeSRedirector* pcstream = 0;
// 	pcstream = new TTreeSRedirector(outfile, "UPDATE");
// 	if (!pcstream) return;
// 	
// 	(*pcstream)<<Form("PWGGApcmQA_%d",run)<<  
//     "run="<<run<<
// 	"meanNSigmadEdxElectron="<< meanNSigmadEdxElectron<<
// 	"meanNSigmadEdxPositron="<< meanNSigmadEdxPositron<<
// 	"NumGammaSectorEtaPosV.="<< &NumGammaSectorEtaPosV<<
// 	"NumGammaSectorEtaNegV.="<< &NumGammaSectorEtaNegV;
// 	"NumGammaSectorEtaPosErrorV.="<< &NumGammaSectorEtaPosErrorV<<
// 	"NumGammaSectorEtaNegErrorV.="<< &NumGammaSectorEtaNegErrorV;
// 	
// 	 (*pcstream)<<Form("PWGGApcmQA_%d",run)<<"\n";
// 	
//     if (!pcstream) return;
//     if (pcstream) { delete pcstream; pcstream = 0; }    
// 	 
	
	
}//end macro

