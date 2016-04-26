/****************************************************************************************************************************
****** 		provided by Gamma Conversion Group, PWGGA, 													*****
******		Friederike Bock, friederike.bock@cern.ch													*****
*****************************************************************************************************************************/

#include <Riostream.h>
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
#include "TBenchmark.h"
#include "TRandom.h"
#include "TLatex.h"
#include "TASImage.h"
#include "TPostScript.h"
#include "TGraphErrors.h"
#include "TArrow.h"
#include "TGraphAsymmErrors.h" 
#include "TGaxis.h"
#include "TMarker.h"
#include "Math/WrappedTF1.h"
#include "Math/BrentRootFinder.h"
#include "CommonHeaders/PlottingGammaConversionHistos.h"
#include "CommonHeaders/PlottingGammaConversionAdditional.h"
#include "CommonHeaders/FittingGammaConversion.h"
#include "CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "CommonHeaders/ConversionFunctions.h"

extern TRandom*	gRandom;
extern TBenchmark*	gBenchmark;
extern TSystem*	gSystem;
extern TMinuit*  	gMinuit;

void ProduceHEPDataFileForRAAPaper2014( TString nameFilePP = "", TString nameFilePbPb = "CombinedResultsPbPb_15_May_2013.root"){

	gROOT->Reset();	
	gROOT->SetStyle("Plain");
	
	StyleSettingsThesis();	
	SetPlotStyle();

   TString dateForOutput = ReturnDateStringForOutput();
   
	Double_t xSection2760GeVpp = 		55.416*1e-3;
	Double_t xSection2760GeVErrpp = 	3.9;
	Double_t xSection2760GeVppINEL = 	62.8*1e9;
	Double_t xSection900GeVppINEL = 	52.5*1e9;
	Double_t xSection7TeVppINEL = 		73.2*1e9;	
	Double_t recalcBarn = 				1e12; //NLO in pbarn!!!!

	TString collisionSystemPP = "pp #sqrt{#it{s}} = 2.76 TeV";		
	TString collisionSystemCent0 = "0-5% Pb-Pb #sqrt{s_{_{NN}}} = 2.76 TeV";		
	TString collisionSystemCent1 = "5-10% Pb-Pb #sqrt{s_{_{NN}}} = 2.76 TeV";		
	TString collisionSystemCent2 = "10-20% Pb-Pb #sqrt{s_{_{NN}}} = 2.76 TeV";		
	TString collisionSystemCent = "0-20% Pb-Pb #sqrt{s_{_{NN}}} = 2.76 TeV";		
	TString collisionSystemSemiCent = "20-40% Pb-Pb #sqrt{s_{_{NN}}} = 2.76 TeV";		
	TString collisionSystemSemiPer = "40-60% Pb-Pb #sqrt{s_{_{NN}}} = 2.76 TeV";		
	TString collisionSystemPer = "60-80% Pb-Pb #sqrt{s_{_{NN}}} = 2.76 TeV";		

	Size_t markerSizeComparison = 1.5;
	
	TFile* fileNeutralPionCombDataPP = new TFile(nameFilePP.Data());
	TGraphAsymmErrors* graphInvXSectionPi0Comb7TeVStatErr= (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0Comb7TeVStatErr");
	TGraphAsymmErrors* graphInvXSectionPi0Comb7TeVSysErr= (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0Comb7TeVSysErr");	
	TGraphAsymmErrors* graphInvXSectionPi0Comb2760GeVStatErr= (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0Comb2760GeVStatErr");
	TGraphAsymmErrors* graphInvXSectionPi0Comb2760GeVSysErr= (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0Comb2760GeVSysErr");
	TGraphAsymmErrors* graphInvXSectionPi0Comb900GeVStatErr= (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0Comb900GeVStatErr");
	TGraphAsymmErrors* graphInvXSectionPi0Comb900GeVSysErr= (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0Comb900GeVSysErr");

	TFile* fCombResults= new TFile(nameFilePbPb.Data());
		TGraphAsymmErrors*	graphYieldCombStatPi02760GeV = (TGraphAsymmErrors*)fCombResults->Get("InvYieldPPComb2760GeV_StatErr");
		TGraphAsymmErrors*	graphYieldCombSysPi02760GeV = (TGraphAsymmErrors*)fCombResults->Get("InvYieldPPComb2760GeV_SysErr");
		TGraphAsymmErrors*	graphYieldPi0CombPbPb0005StatErr = (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbStatErr_0005");
		TGraphAsymmErrors*	graphYieldPi0CombPbPb0005SysErr = (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbSysErr_0005");
		TGraphAsymmErrors*	graphYieldPi0CombPbPb0010StatErr = (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbStatErr_0010");
		TGraphAsymmErrors*	graphYieldPi0CombPbPb0010SysErr = (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbSysErr_0010");
		TGraphAsymmErrors*	graphYieldPi0CombPbPb0510StatErr = (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbStatErr_0510");
		TGraphAsymmErrors*	graphYieldPi0CombPbPb0510SysErr = (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbSysErr_0510");
		TGraphAsymmErrors*	graphYieldPi0CombPbPb1020StatErr = (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbStatErr_1020");
		TGraphAsymmErrors*	graphYieldPi0CombPbPb1020SysErr = (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbSysErr_1020");
		TGraphAsymmErrors*	graphYieldPi0CombPbPb2040StatErr = (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbStatErr_2040");
		TGraphAsymmErrors*	graphYieldPi0CombPbPb2040SysErr = (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbSysErr_2040");
		TGraphAsymmErrors*	graphYieldPi0CombPbPb4060StatErr = (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbStatErr_4060");
		TGraphAsymmErrors*	graphYieldPi0CombPbPb4060SysErr = (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbSysErr_4060");
		TGraphAsymmErrors*	graphYieldPi0CombPbPb6080StatErr = (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbStatErr_6080");
		TGraphAsymmErrors*	graphYieldPi0CombPbPb6080SysErr = (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbSysErr_6080");

		TGraphAsymmErrors*	graphRAAPi0CombPbPb0005StatErr = (TGraphAsymmErrors*)fCombResults->Get("graphRAAStatErr_0005");
		TGraphAsymmErrors*	graphRAAPi0CombPbPb0005SysErr = (TGraphAsymmErrors*)fCombResults->Get("graphRAASysErr_0005");
		TGraphAsymmErrors*	graphRAAPi0CombPbPb0010StatErr = (TGraphAsymmErrors*)fCombResults->Get("graphRAAStatErr_0010");
		TGraphAsymmErrors*	graphRAAPi0CombPbPb0010SysErr = (TGraphAsymmErrors*)fCombResults->Get("graphRAASysErr_0010");
		TGraphAsymmErrors*	graphRAAPi0CombPbPb0510StatErr = (TGraphAsymmErrors*)fCombResults->Get("graphRAAStatErr_0510");
		TGraphAsymmErrors*	graphRAAPi0CombPbPb0510SysErr = (TGraphAsymmErrors*)fCombResults->Get("graphRAASysErr_0510");
		TGraphAsymmErrors*	graphRAAPi0CombPbPb1020StatErr = (TGraphAsymmErrors*)fCombResults->Get("graphRAAStatErr_1020");
		TGraphAsymmErrors*	graphRAAPi0CombPbPb1020SysErr = (TGraphAsymmErrors*)fCombResults->Get("graphRAASysErr_1020");
		TGraphAsymmErrors*	graphRAAPi0CombPbPb2040StatErr = (TGraphAsymmErrors*)fCombResults->Get("graphRAAStatErr_2040");
		TGraphAsymmErrors*	graphRAAPi0CombPbPb2040SysErr = (TGraphAsymmErrors*)fCombResults->Get("graphRAASysErr_2040");
		TGraphAsymmErrors*	graphRAAPi0CombPbPb4060StatErr = (TGraphAsymmErrors*)fCombResults->Get("graphRAAStatErr_4060");
		TGraphAsymmErrors*	graphRAAPi0CombPbPb4060SysErr = (TGraphAsymmErrors*)fCombResults->Get("graphRAASysErr_4060");
		TGraphAsymmErrors*	graphRAAPi0CombPbPb6080StatErr = (TGraphAsymmErrors*)fCombResults->Get("graphRAAStatErr_6080");
		TGraphAsymmErrors*	graphRAAPi0CombPbPb6080SysErr = (TGraphAsymmErrors*)fCombResults->Get("graphRAASysErr_6080");

		TGraphErrors*	graphRAAPi0CombPbPbvsNPartStatErr = (TGraphErrors*)fCombResults->Get("graphRAAvsNPart_StatErr");
		TGraphErrors*	graphRAAPi0CombPbPbvsNPartSysErr = (TGraphErrors*)fCombResults->Get("graphRAAvsNPart_SysErr");
	
	TString headerPaper = "*author: ABELEV \n *reference: ARXIV:1405.3794 : 2014\n *reference:  CERN-PH-EP-2014-091 : 2014\n *reference: Eur. Phys. J. C (2014) 74-3108 : 2014\n *status: Encoded 19 Nov 2014 by F. Bock\n *experiment: CERN-LHC-ALICE\n *detector: ALICE\n *inspireId: 1296306\n *cdsId: 1702183\n *title: Neutral pion production at midrapidity in pp and Pb–Pb collisions at \sqrt{s_{NN}}=2.76TeV\n *comment: CERN-LHC.  Invariant yields of neutral pions at midrapidity in the transverse momentum range 0.6<pT<12 GeV/c measured in Pb-Pb collisions at \sqrt{s_{NN}} = 2.76 TeV are presented for six centrality classes. The pp reference spectrum was measured in the range 0.4<pT<10 GeV/c at the same center-of-mass energy. The nuclear modification factor, RAA, shows a suppression of neutral pions in central Pb-Pb collisions by a factor of up to about 8−10 for 5\leq pT \leq 7 GeV/c. The presented measurements are compared with results at lower center-of-mass energies and with theoretical calculations.";	
		
	cout << headerPaper.Data() << endl;
	
	TString header = "*location: Fig. 5\n *dscomment: Invariant differential yields of PI0 produced in inelastic pp collisions at center-of-mass energy 2.76 TeV\n *reackey: P P --> PI0 X\n *obskey: DN/DPT**2/DYRAP\n *qual: RE: P P --> PI0 X\n *qual: YRAP : 0.0\n *qual: SQRT(S) IN GEV : 2760.0\n *yheader: 1/(2*PI*Nevt)*(d**2N)/(PT*dPT*dY) IN (GEV/c)**-2\n *xheader: PT IN GEV/c\n *data: x : y";
	ProduceHEPDataFile(graphYieldCombStatPi02760GeV, graphYieldCombSysPi02760GeV, header);

	TString headerXSec = "*location: Fig. 5\n *dscomment: Invariant differential cross section of PI0 produced in inelastic pp collisions at center-of-mass energy 2.76 TeV, the uncertainty of \sigma_{inel} of 3.9% is not included in the systematic error \n *reackey: P P --> PI0 X\n *obskey: E*D3SIG/D3P\n *qual: RE: P P --> PI0 X\n *qual: YRAP : 0.0\n *qual: SQRT(S) IN GEV : 2760.0\n *yheader: E d**3sigma/dp**3 IN (pb GEV**-2 c**3)\n *xheader: PT IN GEV/c\n *data: x : y";
	ProduceHEPDataFile(graphInvXSectionPi0Comb2760GeVStatErr, graphInvXSectionPi0Comb2760GeVSysErr, headerXSec);

	TString header0005 = "*location: Fig. 5\n *dscomment: Invariant differential yields of PI0 produced in 0-5% central inelastic PbPb collisions at center-of-mass energy per nucleon 2.76 TeV\n *reackey: PB PB --> PI0 X\n *obskey:  DN/DPT**2/DYRAP\n *qual: RE: PB PB --> PI0 X\n *qual: YRAP : 0.0\n *qual: CENTRALITY : 0.0 TO 5.0\n *qual: SQRT(S)/NUCLEON IN GEV : 2760.0\n *yheader: 1/(2*PI*Nevt)*(d**2N)/(PT*dPT*dY) IN (GEV/c)**-2\n *xheader: PT IN GEV/c\n *data: x : y";
	ProduceHEPDataFile(graphYieldPi0CombPbPb0005StatErr, graphYieldPi0CombPbPb0005SysErr, header0005);

	TString header0010 = "*location: Fig. 5\n *dscomment: Invariant differential yields of PI0 produced in 0-10% central inelastic PbPb collisions at center-of-mass energy per nucleon 2.76 TeV\n *reackey: PB PB --> PI0 X\n *obskey: DN/DPT**2/DYRAP\n *qual: RE: PB PB --> PI0 X\n *qual: YRAP : 0.0\n *qual: CENTRALITY : 0.0 TO 10.0\n *qual: SQRT(S)/NUCLEON IN GEV : 2760.0\n *yheader: 1/(2*PI*Nevt)*(d**2N)/(PT*dPT*dY) IN (GEV/c)**-2\n *xheader: PT IN GEV/c\n *data: x : y";
	ProduceHEPDataFile(graphYieldPi0CombPbPb0010StatErr, graphYieldPi0CombPbPb0010SysErr, header0010);
	
	TString header0510 = "*location: Fig. 5\n *dscomment: Invariant differential yields of PI0 produced in 5-10% central inelastic PbPb collisions at center-of-mass energy per nucleon 2.76 TeV\n *reackey: PB PB --> PI0 X\n *obskey: DN/DPT**2/DYRAP\n *qual: RE: PB PB --> PI0 X\n *qual: YRAP : 0.0\n *qual: CENTRALITY : 5.0 TO 10.0\n *qual: SQRT(S)/NUCLEON IN GEV : 2760.0\n *yheader: 1/(2*PI*Nevt)*(d**2N)/(PT*dPT*dY) IN (GEV/c)**-2\n *xheader: PT IN GEV/c\n *data: x : y";
	ProduceHEPDataFile(graphYieldPi0CombPbPb0510StatErr, graphYieldPi0CombPbPb0510SysErr, header0510);

	TString header1020 = "*location: Fig. 5\n *dscomment: Invariant differential yields of PI0 produced in 10-20% central inelastic PbPb collisions at center-of-mass energy per nucleon 2.76 TeV\n *reackey: PB PB --> PI0 X\n *obskey: DN/DPT**2/DYRAP\n *qual: RE: PB PB --> PI0 X\n *qual: YRAP : 0.0\n *qual: CENTRALITY : 10.0 TO 20.0\n *qual: SQRT(S)/NUCLEON IN GEV : 2760.0\n *yheader: 1/(2*PI*Nevt)*(d**2N)/(PT*dPT*dY) IN (GEV/c)**-2\n *xheader: PT IN GEV/c\n *data: x : y";
	ProduceHEPDataFile(graphYieldPi0CombPbPb1020StatErr, graphYieldPi0CombPbPb1020SysErr, header1020);

	TString header2040 = "*location: Fig. 5\n *dscomment: Invariant differential yields of PI0 produced in 20-40% central inelastic PbPb collisions at center-of-mass energy per nucleon 2.76 TeV\n *reackey: PB PB --> PI0 X\n *obskey: DN/DPT**2/DYRAP\n *qual: RE: PB PB --> PI0 X\n *qual: YRAP : 0.0\n *qual: CENTRALITY : 20.0 TO 40.0\n *qual: SQRT(S)/NUCLEON IN GEV : 2760.0\n *yheader: 1/(2*PI*Nevt)*(d**2N)/(PT*dPT*dY) IN (GEV/c)**-2\n *xheader: PT IN GEV/c\n *data: x : y";
	ProduceHEPDataFile(graphYieldPi0CombPbPb2040StatErr, graphYieldPi0CombPbPb2040SysErr, header2040);

	TString header4060 = "*location: Fig. 5\n *dscomment: Invariant differential yields of PI0 produced in 40-60% central inelastic PbPb collisions at center-of-mass energy per nucleon 2.76 TeV\n *reackey: PB PB --> PI0 X\n *obskey: DN/DPT**2/DYRAP\n *qual: RE: PB PB --> PI0 X\n *qual: YRAP : 0.0\n *qual: CENTRALITY : 40.0 TO 60.0\n *qual: SQRT(S)/NUCLEON IN GEV : 2760.0\n *yheader: 1/(2*PI*Nevt)*(d**2N)/(PT*dPT*dY) IN (GEV/c)**-2\n *xheader: PT IN GEV/c\n *data: x : y";
	ProduceHEPDataFile(graphYieldPi0CombPbPb4060StatErr, graphYieldPi0CombPbPb4060SysErr, header4060);

	TString header6080 = "*location: Fig. 5\n *dscomment: Invariant differential yields of PI0 produced in 60-80% central inelastic PbPb collisions at center-of-mass energy per nucleon 2.76 TeV\n *reackey: PB PB --> PI0 X\n *obskey: DN/DPT**2/DYRAP\n *qual: RE: PB PB --> PI0 X\n *qual: YRAP : 0.0\n *qual: CENTRALITY : 60.0 TO 80.0\n *qual: SQRT(S)/NUCLEON IN GEV : 2760.0\n *yheader: 1/(2*PI*Nevt)*(d**2N)/(PT*dPT*dY) IN (GEV/c)**-2\n *xheader: PT IN GEV/c\n *data: x : y";
	ProduceHEPDataFile(graphYieldPi0CombPbPb6080StatErr, graphYieldPi0CombPbPb6080SysErr, header6080);

	TString headerRAA0005 = "*location: Fig. 6,7,10 \n *dscomment: Nuclear Suppression factor RAA of PI0 produced in 0-5% central inelastic PbPb collisions at center-of-mass energy per nucleon 2.76 TeV, TAA uncertainty of 3.2% not and uncertainty of \sigma_{inel} of 3.9% included in systematics\n *reackey: PB PB --> PI0 X\n *obskey: RAA\n *qual: RE: PB PB --> PI0 X\n *qual: YRAP : 0.0\n *qual: CENTRALITY : 0.0 TO 5.0\n *qual: SQRT(S)/NUCLEON IN GEV : 2760.0\n *yheader: RAA\n *xheader: PT IN GEV/c\n *data: x : y";
	ProduceHEPDataFile(graphRAAPi0CombPbPb0005StatErr, graphRAAPi0CombPbPb0005SysErr, headerRAA0005);
	
	TString headerRAA0010 = "*location: Fig. 6,7,10 \n *dscomment: Nuclear Suppression factor RAA of PI0 produced in 0-10% central inelastic PbPb collisions at center-of-mass energy per nucleon 2.76 TeV, TAA uncertainty of 3.25% not and uncertainty of \sigma_{inel} of 3.9% included in systematics\n *reackey: PB PB --> PI0 X\n *obskey: RAA\n *qual: RE: PB PB --> PI0 X\n *qual: YRAP : 0.0\n *qual: CENTRALITY : 0.0 TO 10.0\n *qual: SQRT(S)/NUCLEON IN GEV : 2760.0\n *yheader:  RAA\n *xheader: PT IN GEV/c\n *data: x : y";
	ProduceHEPDataFile(graphRAAPi0CombPbPb0010StatErr, graphRAAPi0CombPbPb0010SysErr, headerRAA0010);
	
	TString headerRAA0510 = "*location: Fig. 6,7,10 \n *dscomment: Nuclear Suppression factor RAA of PI0 produced in 5-10% central inelastic PbPb collisions at center-of-mass energy per nucleon 2.76 TeV, TAA uncertainty of 3.3% not and uncertainty of \sigma_{inel} of 3.9% included in systematics\n *reackey: PB PB --> PI0 X\n *obskey: RAA\n *qual: RE: PB PB --> PI0 X\n *qual: YRAP : 0.0\n *qual: CENTRALITY : 5.0 TO 10.0\n *qual: SQRT(S)/NUCLEON IN GEV : 2760.0\n *yheader:  RAA\n *xheader: PT IN GEV/c\n *data: x : y";
	ProduceHEPDataFile(graphRAAPi0CombPbPb0510StatErr, graphRAAPi0CombPbPb0510SysErr, headerRAA0510);

	TString headerRAA1020 = "*location: Fig. 6,7,10 \n *dscomment: Nuclear Suppression factor RAA of PI0 produced in 10-20% central inelastic PbPb collisions at center-of-mass energy per nucleon 2.76 TeV, TAA uncertainty of 3.1% not and uncertainty of \sigma_{inel} of 3.9% included in systematics\n *reackey: PB PB --> PI0 X\n *obskey: RAA\n *qual: RE: PB PB --> PI0 X\n *qual: YRAP : 0.0\n *qual: CENTRALITY : 10.0 TO 20.0\n *qual: SQRT(S)/NUCLEON IN GEV : 2760.0\n *yheader:  RAA\n *xheader: PT IN GEV/c\n *data: x : y";
	ProduceHEPDataFile(graphRAAPi0CombPbPb1020StatErr, graphRAAPi0CombPbPb1020SysErr, headerRAA1020);

	TString headerRAA2040 = "*location: Fig. 6,7,10 \n *dscomment: Nuclear Suppression factor RAA of PI0 produced in 20-40% central inelastic PbPb collisions at center-of-mass energy per nucleon 2.76 TeV, TAA uncertainty of 3.3% not and uncertainty of \sigma_{inel} of 3.9% included in systematics\n *reackey: PB PB --> PI0 X\n *obskey: RAA\n *qual: RE: PB PB --> PI0 X\n *qual: YRAP : 0.0\n *qual: CENTRALITY : 20.0 TO 40.0\n *qual: SQRT(S)/NUCLEON IN GEV : 2760.0\n *yheader:  RAA\n *xheader: PT IN GEV/c\n *data: x : y";
	ProduceHEPDataFile(graphRAAPi0CombPbPb2040StatErr, graphRAAPi0CombPbPb2040SysErr, headerRAA2040);

	TString headerRAA4060 = "*location: Fig. 6,7,10 \n *dscomment: Nuclear Suppression factor RAA of PI0 produced in 40-60% central inelastic PbPb collisions at center-of-mass energy per nucleon 2.76 TeV, TAA uncertainty of 4.9% not and uncertainty of \sigma_{inel} of 3.9% included in systematics\n *reackey: PB PB --> PI0 X\n *obskey: RAA\n *qual: RE: PB PB --> PI0 X\n *qual: YRAP : 0.0\n *qual: CENTRALITY : 40.0 TO 60.0\n *qual: SQRT(S)/NUCLEON IN GEV : 2760.0\n *yheader:  RAA\n *xheader: PT IN GEV/c\n *data: x : y";
	ProduceHEPDataFile(graphRAAPi0CombPbPb4060StatErr, graphRAAPi0CombPbPb4060SysErr, headerRAA4060);

	TString headerRAA6080 = "*location: Fig. 6,7,10 \n *dscomment: Nuclear Suppression factor RAA of PI0 produced in 60-80% central inelastic PbPb collisions at center-of-mass energy per nucleon 2.76 TeV, TAA uncertainty of 6.2% not and uncertainty of \sigma_{inel} of 3.9% included in systematics\n *reackey: PB PB --> PI0 X\n *obskey: RAA\n *qual: RE: PB PB --> PI0 X\n *qual: YRAP : 0.0\n *qual: CENTRALITY : 60.0 TO 80.0\n *qual: SQRT(S)/NUCLEON IN GEV : 2760.0\n *yheader:  RAA\n *xheader: PT IN GEV/c\n *data: x : y";
	ProduceHEPDataFile(graphRAAPi0CombPbPb6080StatErr, graphRAAPi0CombPbPb6080SysErr, headerRAA6080);

	TString headerRAANPart = "*location: Fig. 9 \n *dscomment: Nuclear Suppression factor RAA at PT = 7 GeV/c of PI0 produced in inelastic PbPb collisions at center-of-mass energy per nucleon 2.76 TeV vs NPart, NPart uncertainty indicated by bin width in x \n *reackey: PB PB --> PI0 X\n *obskey: RAA\n *qual: RE: PB PB --> PI0 X\n *qual: YRAP : 0.0\n *qual: PT : 7 GEV/c\n *qual: CENTRALITY : 0.0 TO 80.0\n *qual: SQRT(S)/NUCLEON IN GEV : 2760.0\n *yheader:  RAA\n *xheader: NPART\n *data: x : y";
	ProduceHEPDataFileNPart(graphRAAPi0CombPbPbvsNPartStatErr, graphRAAPi0CombPbPbvsNPartSysErr, headerRAANPart);
	
	
	
}