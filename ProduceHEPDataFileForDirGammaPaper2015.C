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

void ProduceHEPDataFileForDirGammaPaper2015( TString nameFilePbPb = "Gamma_CombResults_PbPb_2.76TeV.root", Bool_t simpleErrors = kTRUE){

	gROOT->Reset();	
	gROOT->SetStyle("Plain");
	
	StyleSettingsThesis();	
	SetPlotStyle();

   TString dateForOutput = ReturnDateStringForOutput();
   
	Double_t xSection2760GeVpp 			= 55.416*1e-3;
	Double_t xSection2760GeVErrpp 		= 3.9;
	Double_t xSection2760GeVppINEL 		= 62.8*1e9;
	Double_t xSection900GeVppINEL	 	= 52.5*1e9;
	Double_t xSection7TeVppINEL 		= 73.2*1e9;	
	Double_t recalcBarn 				= 1e12; //NLO in pbarn!!!!

	TString collisionSystemCent 		= "0-20% Pb-Pb #sqrt{s_{_{NN}}} = 2.76 TeV";		
	TString collisionSystemSemiCent 	= "20-40% Pb-Pb #sqrt{s_{_{NN}}} = 2.76 TeV";		
	TString collisionSystemPer 			= "40-80% Pb-Pb #sqrt{s_{_{NN}}} = 2.76 TeV";		

	Size_t markerSizeComparison 		= 1.5;
	
	Double_t newBinsComb[21] 			= {	0.9, 1.1, 1.3, 1.5, 1.7, 
											1.9, 2.1, 2.3, 2.5, 2.7, 
											3.0, 3.3, 3.7, 4.1, 4.6, 
											5.4, 6.2, 7.0, 8.0, 11.0, 
											14.0};
	Int_t nPointsTot					= 20;
	
	TFile* fCombResults= new TFile(nameFilePbPb.Data());
		TDirectory* 		directoryCombGamma0020 					= (TDirectory*)fCombResults->Get("Gamma_PbPb_2.76TeV_0-20%"); 
		TGraphAsymmErrors*	graphDirGammaYieldStat0020 				= (TGraphAsymmErrors*)directoryCombGamma0020->Get("DirGammaSpec_comb_StatErr");
		TGraphAsymmErrors*	graphDirGammaYieldSys0020 				= (TGraphAsymmErrors*)directoryCombGamma0020->Get("DirGammaSpec_comb_SysErr");
		TGraphAsymmErrors*	graphDirGammaYieldSysA0020 				= (TGraphAsymmErrors*)directoryCombGamma0020->Get("DirGammaSpec_comb_SysAErr");
		TGraphAsymmErrors*	graphDirGammaYieldSysB0020 				= (TGraphAsymmErrors*)directoryCombGamma0020->Get("DirGammaSpec_comb_SysBErr");
		TGraphAsymmErrors*	graphDirGammaYieldSysC0020 				= (TGraphAsymmErrors*)directoryCombGamma0020->Get("DirGammaSpec_comb_SysCErr");
		TGraphAsymmErrors*	graphDRStat0020 						= (TGraphAsymmErrors*)directoryCombGamma0020->Get("DR_comb_StatErr");
		TGraphAsymmErrors*	graphDRSys0020 							= (TGraphAsymmErrors*)directoryCombGamma0020->Get("DR_comb_SysErr");
		TGraphAsymmErrors*	graphDRSysA0020 						= (TGraphAsymmErrors*)directoryCombGamma0020->Get("DR_comb_SysAErr");
		TGraphAsymmErrors*	graphDRSysB0020 						= (TGraphAsymmErrors*)directoryCombGamma0020->Get("DR_comb_SysBErr")	;
		TGraphAsymmErrors*	graphDRSysC0020 						= (TGraphAsymmErrors*)directoryCombGamma0020->Get("DR_comb_SysCErr");
		TGraphAsymmErrors*	graphIncGammaYieldStat0020 				= (TGraphAsymmErrors*)directoryCombGamma0020->Get("IncGammaSpec_comb_StatErr");
		TGraphAsymmErrors*	graphIncGammaYieldSys0020 				= (TGraphAsymmErrors*)directoryCombGamma0020->Get("IncGammaSpec_comb_SysErr");
		TGraphAsymmErrors*	graphIncGammaYieldSysA0020 				= (TGraphAsymmErrors*)directoryCombGamma0020->Get("IncGammaSpec_comb_SysAErr");
		TGraphAsymmErrors*	graphIncGammaYieldSysB0020 				= (TGraphAsymmErrors*)directoryCombGamma0020->Get("IncGammaSpec_comb_SysBErr");
		TGraphAsymmErrors*	graphIncGammaYieldSysC0020 				= (TGraphAsymmErrors*)directoryCombGamma0020->Get("IncGammaSpec_comb_SysCErr");

		TDirectory* 		directoryCombGamma2040 					= (TDirectory*)fCombResults->Get("Gamma_PbPb_2.76TeV_20-40%"); 
		TGraphAsymmErrors*	graphDirGammaYieldStat2040 				= (TGraphAsymmErrors*)directoryCombGamma2040->Get("DirGammaSpec_comb_StatErr");
		TGraphAsymmErrors*	graphDirGammaYieldSys2040 				= (TGraphAsymmErrors*)directoryCombGamma2040->Get("DirGammaSpec_comb_SysErr");
		TGraphAsymmErrors*	graphDirGammaYieldSysA2040 				= (TGraphAsymmErrors*)directoryCombGamma2040->Get("DirGammaSpec_comb_SysAErr");
		TGraphAsymmErrors*	graphDirGammaYieldSysB2040 				= (TGraphAsymmErrors*)directoryCombGamma2040->Get("DirGammaSpec_comb_SysBErr");
		TGraphAsymmErrors*	graphDirGammaYieldSysC2040 				= (TGraphAsymmErrors*)directoryCombGamma2040->Get("DirGammaSpec_comb_SysCErr");
		TGraphAsymmErrors*	graphDirGammaYieldUpper2040 			= (TGraphAsymmErrors*)directoryCombGamma2040->Get("DirGammaSpec_comb_upperLimits");
		TGraphAsymmErrors*	graphDRStat2040 						= (TGraphAsymmErrors*)directoryCombGamma2040->Get("DR_comb_StatErr");
		TGraphAsymmErrors*	graphDRSys2040 							= (TGraphAsymmErrors*)directoryCombGamma2040->Get("DR_comb_SysErr");
		TGraphAsymmErrors*	graphDRSysA2040 						= (TGraphAsymmErrors*)directoryCombGamma2040->Get("DR_comb_SysAErr");
		TGraphAsymmErrors*	graphDRSysB2040 						= (TGraphAsymmErrors*)directoryCombGamma2040->Get("DR_comb_SysBErr");
		TGraphAsymmErrors*	graphDRSysC2040 						= (TGraphAsymmErrors*)directoryCombGamma2040->Get("DR_comb_SysCErr");
		TGraphAsymmErrors*	graphIncGammaYieldStat2040 				= (TGraphAsymmErrors*)directoryCombGamma2040->Get("IncGammaSpec_comb_StatErr");
		TGraphAsymmErrors*	graphIncGammaYieldSys2040 				= (TGraphAsymmErrors*)directoryCombGamma2040->Get("IncGammaSpec_comb_SysErr");
		TGraphAsymmErrors*	graphIncGammaYieldSysA2040 				= (TGraphAsymmErrors*)directoryCombGamma2040->Get("IncGammaSpec_comb_SysAErr");
		TGraphAsymmErrors*	graphIncGammaYieldSysB2040 				= (TGraphAsymmErrors*)directoryCombGamma2040->Get("IncGammaSpec_comb_SysBErr");
		TGraphAsymmErrors*	graphIncGammaYieldSysC2040 				= (TGraphAsymmErrors*)directoryCombGamma2040->Get("IncGammaSpec_comb_SysCErr");
		
		TDirectory* 		directoryCombGamma4080 					= (TDirectory*)fCombResults->Get("Gamma_PbPb_2.76TeV_40-80%"); 
		TGraphAsymmErrors*	graphDirGammaYieldStat4080 				= (TGraphAsymmErrors*)directoryCombGamma4080->Get("DirGammaSpec_comb_StatErr");
		TGraphAsymmErrors*	graphDirGammaYieldSys4080 				= (TGraphAsymmErrors*)directoryCombGamma4080->Get("DirGammaSpec_comb_SysErr");
		TGraphAsymmErrors*	graphDirGammaYieldSysA4080 				= (TGraphAsymmErrors*)directoryCombGamma4080->Get("DirGammaSpec_comb_SysAErr");
		TGraphAsymmErrors*	graphDirGammaYieldSysB4080 				= (TGraphAsymmErrors*)directoryCombGamma4080->Get("DirGammaSpec_comb_SysBErr");
		TGraphAsymmErrors*	graphDirGammaYieldSysC4080 				= (TGraphAsymmErrors*)directoryCombGamma4080->Get("DirGammaSpec_comb_SysCErr");
		TGraphAsymmErrors*	graphDirGammaYieldUpper4080 			= (TGraphAsymmErrors*)directoryCombGamma4080->Get("DirGammaSpec_comb_upperLimits");
		TGraphAsymmErrors*	graphDRStat4080 						= (TGraphAsymmErrors*)directoryCombGamma4080->Get("DR_comb_StatErr");
		TGraphAsymmErrors*	graphDRSys4080 							= (TGraphAsymmErrors*)directoryCombGamma4080->Get("DR_comb_SysErr");
		TGraphAsymmErrors*	graphDRSysA4080 						= (TGraphAsymmErrors*)directoryCombGamma4080->Get("DR_comb_SysAErr");
		TGraphAsymmErrors*	graphDRSysB4080 						= (TGraphAsymmErrors*)directoryCombGamma4080->Get("DR_comb_SysBErr");
		TGraphAsymmErrors*	graphDRSysC4080 						= (TGraphAsymmErrors*)directoryCombGamma4080->Get("DR_comb_SysCErr");
		TGraphAsymmErrors*	graphIncGammaYieldStat4080 				= (TGraphAsymmErrors*)directoryCombGamma4080->Get("IncGammaSpec_comb_StatErr");
		TGraphAsymmErrors*	graphIncGammaYieldSys4080 				= (TGraphAsymmErrors*)directoryCombGamma4080->Get("IncGammaSpec_comb_SysErr");
		TGraphAsymmErrors*	graphIncGammaYieldSysA4080 				= (TGraphAsymmErrors*)directoryCombGamma4080->Get("IncGammaSpec_comb_SysAErr");
		TGraphAsymmErrors*	graphIncGammaYieldSysB4080 				= (TGraphAsymmErrors*)directoryCombGamma4080->Get("IncGammaSpec_comb_SysBErr");
		TGraphAsymmErrors*	graphIncGammaYieldSysC4080 				= (TGraphAsymmErrors*)directoryCombGamma4080->Get("IncGammaSpec_comb_SysCErr");

		
	TString headerPaper = " *author: ADAM \n *reference: CERN-PH-EP-2015-254 : 2015\n *status: Encoded 15 Sep 2015 by F. Bock\n *experiment: CERN-LHC-ALICE\n *detector: ALICE\n \n *title: Direct photon production in Pb-Pb collisions at \sqrt{s_{NN}}=2.76TeV\n *comment: CERN-LHC.  Direct photon production at mid-rapidity in Pb-Pb collisions at \sqrt{s_{NN}}=2.76TeV TeV was studied in the transverse momentum range 0.9<pT<14 GeV/c. Photons were detected with the highly segmented electromagnetic calorimeter PHOS and via conversions in the ALICE detector material with the $e^+e^-$ pair reconstructed in the central tracking system. The results of the two methods were combined and direct photon spectra were measured for the 0--20%, 20--40%, and 40--80% centrality classes. For all three classes, agreement was found with perturbative QCD calculations for pT \gtrsim 5 GeV/c. Direct photon spectra down to pT \approx 1 GeV/c could be extracted for the 20--40% and 0--20% centrality classes. The significance of the direct photon signal for 0.9<pT<2.1 GeV/c is 2.6\sigma for the 0--20% class. The spectrum in this pT range and centrality class can be described by an exponential with an inverse slope parameter of (297 \pm 12 (stat) \pm 41 (syst) ) MeV. State-of-the-art models for photon production in heavy-ion collisions agree with the data within uncertainties.\n";	
	
	if (simpleErrors){
		cout << headerPaper.Data() << endl;
		TString header0020 = " *location: Fig. 5\n *dscomment: Invariant differential yields of direct GAMMA produced in 0-20% central inelastic PbPb collisions at center-of-mass energy per nucleon 2.76 TeV\n *reackey: PB PB --> GAMMA X\n *obskey: DN/DPT**2/DYRAP\n *qual: RE: PB PB --> GAMMA X\n *qual: YRAP : 0.0\n *qual: CENTRALITY : 0.0 TO 20.0\n *qual: SQRT(S)/NUCLEON IN GEV : 2760.0\n *yheader: 1/(2*PI*Nevt)*(d**2N)/(PT*dPT*dY) IN (GEV/c)**-2 \n *xheader: PT IN GEV/c\n *data: x : y";
		ProduceHEPDataFileWithUpperLimits(graphDirGammaYieldStat0020, graphDirGammaYieldSys0020, NULL, newBinsComb, nPointsTot, header0020);
			
		cout << endl;	
		TString header2040 = " *location: Fig. 5\n *dscomment: Invariant differential yields of direct GAMMA produced in 20-40% central inelastic PbPb collisions at center-of-mass energy per nucleon 2.76 TeV\n *reackey: PB PB --> GAMMA X\n *obskey: DN/DPT**2/DYRAP\n *qual: RE: PB PB --> GAMMA X\n *qual: YRAP : 0.0\n *qual: CENTRALITY : 20.0 TO 40.0\n *qual: SQRT(S)/NUCLEON IN GEV : 2760.0\n *yheader: 1/(2*PI*Nevt)*(d**2N)/(PT*dPT*dY) IN (GEV/c)**-2\n *xheader: PT IN GEV/c\n *data: x : y";
		ProduceHEPDataFileWithUpperLimits(graphDirGammaYieldStat2040, graphDirGammaYieldSys2040, graphDirGammaYieldUpper2040, newBinsComb, nPointsTot, header2040);
		
		cout << endl;
		TString header4080 = " *location: Fig. 5\n *dscomment: Invariant differential yields of direct GAMMA produced in 40-80% central inelastic PbPb collisions at center-of-mass energy per nucleon 2.76 TeV\n *reackey: PB PB --> GAMMA X\n *obskey: DN/DPT**2/DYRAP\n *qual: RE: PB PB --> GAMMA X\n *qual: YRAP : 0.0\n *qual: CENTRALITY : 40.0 TO 80.0\n *qual: SQRT(S)/NUCLEON IN GEV : 2760.0\n *yheader: 1/(2*PI*Nevt)*(d**2N)/(PT*dPT*dY) IN (GEV/c)**-2\n *xheader: PT IN GEV/c\n *data: x : y";
		ProduceHEPDataFileWithUpperLimits(graphDirGammaYieldStat4080, graphDirGammaYieldSys4080, graphDirGammaYieldUpper4080, newBinsComb, nPointsTot, header4080);

		cout << endl;		
		TString headerRAA0020 = " *location: Fig. 4\n *dscomment: Double Ratio RGAMMA  in 0-20% central inelastic PbPb collisions at center-of-mass energy per nucleon 2.76 TeV\n *reackey: PB PB --> GAMMA X\n *obskey: RGAMMA\n *qual: RE: PB PB --> GAMMA X\n *qual: YRAP : 0.0\n *qual: CENTRALITY : 0.0 TO 20.0\n *qual: SQRT(S)/NUCLEON IN GEV : 2760.0\n *yheader: RGAMMA \n *xheader: PT IN GEV/c\n *data: x : y";
		ProduceHEPDataFileWithUpperLimits(graphDRStat0020, graphDRSys0020, NULL, newBinsComb, nPointsTot, headerRAA0020);
			
		cout << endl;	
		TString headerRAA2040 = " *location: Fig. 4\n *dscomment: Double Ratio RGAMMA in 20-40% central inelastic PbPb collisions at center-of-mass energy per nucleon 2.76 TeV\n *reackey: PB PB --> GAMMA X\n *obskey: RGAMMA\n *qual: RE: PB PB --> GAMMA X\n *qual: YRAP : 0.0\n *qual: CENTRALITY : 20.0 TO 40.0\n *qual: SQRT(S)/NUCLEON IN GEV : 2760.0\n *yheader: RGAMMA\n *xheader: PT IN GEV/c\n *data: x : y";
		ProduceHEPDataFileWithUpperLimits(graphDRStat2040, graphDRSys2040, NULL, newBinsComb, nPointsTot, headerRAA2040);
		
		cout << endl;
		TString headerRAA4080 = " *location: Fig. 4\n *dscomment: Double Ratio RGAMMA in 40-80% central inelastic PbPb collisions at center-of-mass energy per nucleon 2.76 TeV\n *reackey: PB PB --> GAMMA X\n *obskey: RGAMMA\n *qual: RE: PB PB --> GAMMA X\n *qual: YRAP : 0.0\n *qual: CENTRALITY : 40.0 TO 80.0\n *qual: SQRT(S)/NUCLEON IN GEV : 2760.0\n *yheader: RGAMMA\n *xheader: PT IN GEV/c\n *data: x : y";
		ProduceHEPDataFileWithUpperLimits(graphDRStat4080, graphDRSys4080, NULL, newBinsComb, nPointsTot, headerRAA4080);

		TString headerInc0020 = " *location: not in paper\n *dscomment: Invariant differential yields of inclusive GAMMA produced in 0-20% central inelastic PbPb collisions at center-of-mass energy per nucleon 2.76 TeV\n *reackey: PB PB --> GAMMA X\n *obskey: DN/DPT**2/DYRAP\n *qual: RE: PB PB --> GAMMA X\n *qual: YRAP : 0.0\n *qual: CENTRALITY : 0.0 TO 20.0\n *qual: SQRT(S)/NUCLEON IN GEV : 2760.0\n *yheader: 1/(2*PI*Nevt)*(d**2N)/(PT*dPT*dY) IN (GEV/c)**-2 \n *xheader: PT IN GEV/c\n *data: x : y";
		ProduceHEPDataFileWithUpperLimits(graphIncGammaYieldStat0020, graphIncGammaYieldSys0020, NULL, newBinsComb, nPointsTot, headerInc0020);
			
		cout << endl;	
		TString headerInc2040 = " *location: not in paper\n *dscomment: Invariant differential yields of inclusive GAMMA produced in 20-40% central inelastic PbPb collisions at center-of-mass energy per nucleon 2.76 TeV\n *reackey: PB PB --> GAMMA X\n *obskey: DN/DPT**2/DYRAP\n *qual: RE: PB PB --> GAMMA X\n *qual: YRAP : 0.0\n *qual: CENTRALITY : 20.0 TO 40.0\n *qual: SQRT(S)/NUCLEON IN GEV : 2760.0\n *yheader: 1/(2*PI*Nevt)*(d**2N)/(PT*dPT*dY) IN (GEV/c)**-2\n *xheader: PT IN GEV/c\n *data: x : y";
		ProduceHEPDataFileWithUpperLimits(graphIncGammaYieldStat2040, graphIncGammaYieldSys2040, NULL, newBinsComb, nPointsTot, headerInc2040);
		
		cout << endl;
		TString headerInc4080 = " *location: not in paper\n *dscomment: Invariant differential yields of inclusive GAMMA produced in 40-80% central inelastic PbPb collisions at center-of-mass energy per nucleon 2.76 TeV\n *reackey: PB PB --> GAMMA X\n *obskey: DN/DPT**2/DYRAP\n *qual: RE: PB PB --> GAMMA X\n *qual: YRAP : 0.0\n *qual: CENTRALITY : 40.0 TO 80.0\n *qual: SQRT(S)/NUCLEON IN GEV : 2760.0\n *yheader: 1/(2*PI*Nevt)*(d**2N)/(PT*dPT*dY) IN (GEV/c)**-2\n *xheader: PT IN GEV/c\n *data: x : y";
		ProduceHEPDataFileWithUpperLimits(graphIncGammaYieldStat4080, graphIncGammaYieldSys4080, NULL, newBinsComb, nPointsTot, headerInc4080);
		
	} else {
		cout << headerPaper.Data() << endl;
		TString header0020 = " *location: Fig. 5\n *dscomment: Invariant differential yields of direct GAMMA produced in 0-20% central inelastic PbPb collisions at center-of-mass energy per nucleon 2.76 TeV\n *reackey: PB PB --> GAMMA X\n *obskey: DN/DPT**2/DYRAP\n *qual: RE: PB PB --> GAMMA X\n *qual: YRAP : 0.0\n *qual: CENTRALITY : 0.0 TO 20.0\n *qual: SQRT(S)/NUCLEON IN GEV : 2760.0\n *yheader: 1/(2*PI*Nevt)*(d**2N)/(PT*dPT*dY) IN (GEV/c)**-2 \n *xheader: PT IN GEV/c\n *data: x : y";
		ProduceHEPDataFileWithUpperLimitsErrorsSplit(graphDirGammaYieldStat0020, graphDirGammaYieldSysA0020, graphDirGammaYieldSysB0020, graphDirGammaYieldSysC0020, NULL, newBinsComb, nPointsTot, header0020);
			
		cout << endl;	
		TString header2040 = " *location: Fig. 5\n *dscomment: Invariant differential yields of direct GAMMA produced in 20-40% central inelastic PbPb collisions at center-of-mass energy per nucleon 2.76 TeV\n *reackey: PB PB --> GAMMA X\n *obskey: DN/DPT**2/DYRAP\n *qual: RE: PB PB --> GAMMA X\n *qual: YRAP : 0.0\n *qual: CENTRALITY : 20.0 TO 40.0\n *qual: SQRT(S)/NUCLEON IN GEV : 2760.0\n *yheader: 1/(2*PI*Nevt)*(d**2N)/(PT*dPT*dY) IN (GEV/c)**-2\n *xheader: PT IN GEV/c\n *data: x : y";
		ProduceHEPDataFileWithUpperLimitsErrorsSplit(graphDirGammaYieldStat2040, graphDirGammaYieldSysA2040, graphDirGammaYieldSysC2040, graphDirGammaYieldSysC2040, graphDirGammaYieldUpper2040, newBinsComb, nPointsTot, header2040);
		
		cout << endl;
		TString header4080 = " *location: Fig. 5\n *dscomment: Invariant differential yields of direct GAMMA produced in 40-80% central inelastic PbPb collisions at center-of-mass energy per nucleon 2.76 TeV\n *reackey: PB PB --> GAMMA X\n *obskey: DN/DPT**2/DYRAP\n *qual: RE: PB PB --> GAMMA X\n *qual: YRAP : 0.0\n *qual: CENTRALITY : 40.0 TO 80.0\n *qual: SQRT(S)/NUCLEON IN GEV : 2760.0\n *yheader: 1/(2*PI*Nevt)*(d**2N)/(PT*dPT*dY) IN (GEV/c)**-2\n *xheader: PT IN GEV/c\n *data: x : y";
		ProduceHEPDataFileWithUpperLimitsErrorsSplit(graphDirGammaYieldStat4080, graphDirGammaYieldSysA4080, graphDirGammaYieldSysB4080, graphDirGammaYieldSysC4080, graphDirGammaYieldUpper4080, newBinsComb, nPointsTot, header4080);

		cout << endl;		
		TString headerRAA0020 = " *location: Fig. 4\n *dscomment: Double Ratio RGAMMA  in 0-20% central inelastic PbPb collisions at center-of-mass energy per nucleon 2.76 TeV\n *reackey: PB PB --> GAMMA X\n *obskey: RGAMMA\n *qual: RE: PB PB --> GAMMA X\n *qual: YRAP : 0.0\n *qual: CENTRALITY : 0.0 TO 20.0\n *qual: SQRT(S)/NUCLEON IN GEV : 2760.0\n *yheader: RGAMMA \n *xheader: PT IN GEV/c\n *data: x : y";
		ProduceHEPDataFileWithUpperLimitsErrorsSplit(graphDRStat0020, graphDRSysA0020, graphDRSysB0020, graphDRSysC0020, NULL, newBinsComb, nPointsTot, headerRAA0020);
			
		cout << endl;	
		TString headerRAA2040 = " *location: Fig. 4\n *dscomment: Double Ratio RGAMMA in 20-40% central inelastic PbPb collisions at center-of-mass energy per nucleon 2.76 TeV\n *reackey: PB PB --> GAMMA X\n *obskey: RGAMMA\n *qual: RE: PB PB --> GAMMA X\n *qual: YRAP : 0.0\n *qual: CENTRALITY : 20.0 TO 40.0\n *qual: SQRT(S)/NUCLEON IN GEV : 2760.0\n *yheader: RGAMMA\n *xheader: PT IN GEV/c\n *data: x : y";
		ProduceHEPDataFileWithUpperLimitsErrorsSplit(graphDRStat2040, graphDRSysA2040, graphDRSysB2040, graphDRSysC2040, NULL, newBinsComb, nPointsTot, headerRAA2040);
		
		cout << endl;
		TString headerRAA4080 = " *location: Fig. 4\n *dscomment: Double Ratio RGAMMA in 40-80% central inelastic PbPb collisions at center-of-mass energy per nucleon 2.76 TeV\n *reackey: PB PB --> GAMMA X\n *obskey: RGAMMA\n *qual: RE: PB PB --> GAMMA X\n *qual: YRAP : 0.0\n *qual: CENTRALITY : 40.0 TO 80.0\n *qual: SQRT(S)/NUCLEON IN GEV : 2760.0\n *yheader: RGAMMA\n *xheader: PT IN GEV/c\n *data: x : y";
		ProduceHEPDataFileWithUpperLimitsErrorsSplit(graphDRStat4080, graphDRSysA4080, graphDRSysB4080, graphDRSysC4080, NULL, newBinsComb, nPointsTot, headerRAA4080);

		TString headerInc0020 = " *location: not in paper\n *dscomment: Invariant differential yields of inclusive GAMMA produced in 0-20% central inelastic PbPb collisions at center-of-mass energy per nucleon 2.76 TeV\n *reackey: PB PB --> GAMMA X\n *obskey: DN/DPT**2/DYRAP\n *qual: RE: PB PB --> GAMMA X\n *qual: YRAP : 0.0\n *qual: CENTRALITY : 0.0 TO 20.0\n *qual: SQRT(S)/NUCLEON IN GEV : 2760.0\n *yheader: 1/(2*PI*Nevt)*(d**2N)/(PT*dPT*dY) IN (GEV/c)**-2 \n *xheader: PT IN GEV/c\n *data: x : y";
		ProduceHEPDataFileWithUpperLimitsErrorsSplit(graphIncGammaYieldStat0020, graphIncGammaYieldSysA0020, graphIncGammaYieldSysB0020, graphIncGammaYieldSysC0020, NULL, newBinsComb, nPointsTot, headerInc0020);
			
		cout << endl;	
		TString headerInc2040 = " *location: not in paper\n *dscomment: Invariant differential yields of inclusive GAMMA produced in 20-40% central inelastic PbPb collisions at center-of-mass energy per nucleon 2.76 TeV\n *reackey: PB PB --> GAMMA X\n *obskey: DN/DPT**2/DYRAP\n *qual: RE: PB PB --> GAMMA X\n *qual: YRAP : 0.0\n *qual: CENTRALITY : 20.0 TO 40.0\n *qual: SQRT(S)/NUCLEON IN GEV : 2760.0\n *yheader: 1/(2*PI*Nevt)*(d**2N)/(PT*dPT*dY) IN (GEV/c)**-2\n *xheader: PT IN GEV/c\n *data: x : y";
		ProduceHEPDataFileWithUpperLimitsErrorsSplit(graphIncGammaYieldStat2040, graphIncGammaYieldSysA2040, graphIncGammaYieldSysB2040, graphIncGammaYieldSysC2040, NULL, newBinsComb, nPointsTot, headerInc2040);
		
		cout << endl;
		TString headerInc4080 = " *location: not in paper\n *dscomment: Invariant differential yields of inclusive GAMMA produced in 40-80% central inelastic PbPb collisions at center-of-mass energy per nucleon 2.76 TeV\n *reackey: PB PB --> GAMMA X\n *obskey: DN/DPT**2/DYRAP\n *qual: RE: PB PB --> GAMMA X\n *qual: YRAP : 0.0\n *qual: CENTRALITY : 40.0 TO 80.0\n *qual: SQRT(S)/NUCLEON IN GEV : 2760.0\n *yheader: 1/(2*PI*Nevt)*(d**2N)/(PT*dPT*dY) IN (GEV/c)**-2\n *xheader: PT IN GEV/c\n *data: x : y";
		ProduceHEPDataFileWithUpperLimitsErrorsSplit(graphIncGammaYieldStat4080, graphIncGammaYieldSysA4080, graphIncGammaYieldSysB4080, graphIncGammaYieldSysC4080, NULL, newBinsComb, nPointsTot, headerInc4080);
		
	}	
	
}