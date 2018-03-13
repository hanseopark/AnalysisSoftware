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

void ProduceHEPDataFileForPPDirGammaPaper2018(
  TString nameFilepp2760GeV                       = "CombinedGammaResultPP2760GeV_2018_03_08.root",
  TString nameFilepp8TeV                          = "CombinedGammaResultPP8TeV_2018_03_13.root"
){

	gROOT->Reset();
	gROOT->SetStyle("Plain");
	
	StyleSettingsThesis();
	SetPlotStyle();

   TString dateForOutput = ReturnDateStringForOutput();
	
	Double_t newBinsCombpp2760GeV[19] 			        = {	0.4, 0.6,  0.8, 1.0, 1.2,
                                											1.4, 1.6,  1.8, 2.0, 2.2,
                                											2.4, 2.6, 3.0, 3.5, 4.0,
                                											5.0, 6.0, 8.0,  10.0 };
	Double_t newBinsCombpp8TeV[23] 			            = {	0.3, 0.4,  0.6, 0.8, 1.0,
                                											1.2, 1.4,  1.6, 1.8, 2.0,
                                											2.2, 2.4, 2.7, 3.0, 3.5,
                                											4.0, 4.5, 5.0,  6.0, 7.0,
                                                      9.0, 12.0,16.0};
	Int_t nPointsTotpp2760GeV					              = 18;
	Int_t nPointsTotpp8TeV					                = 22;
	
	TFile* fCombResultspp2760GeV= new TFile(nameFilepp2760GeV.Data());
		TDirectory* 		directoryCombGammapp2760GeV 					  = (TDirectory*)fCombResultspp2760GeV->Get("Gamma2.76TeV");
		TGraphAsymmErrors*	graphDirGammaYieldStatpp2760GeV 		= (TGraphAsymmErrors*)directoryCombGammapp2760GeV->Get("graphInvXSectionDirGammaNonFitStatErr");
		TGraphAsymmErrors*	graphDirGammaYieldSyspp2760GeV 			= (TGraphAsymmErrors*)directoryCombGammapp2760GeV->Get("graphInvXSectionDirGammaNonFitSysErr");
    TGraphAsymmErrors*	graphDirGammaYieldUpperpp2760GeV 		= (TGraphAsymmErrors*)directoryCombGammapp2760GeV->Get("graphInvXSectionDirGammaNonFitSumErrAr");
		TGraphAsymmErrors*	graphDRStatpp2760GeV 						    = (TGraphAsymmErrors*)directoryCombGammapp2760GeV->Get("graphRGammaNonFitCombStatErr");
		TGraphAsymmErrors*	graphDRSyspp2760GeV 							  = (TGraphAsymmErrors*)directoryCombGammapp2760GeV->Get("graphRGammaNonFitCombSysErr");
		TGraphAsymmErrors*	graphIncGammaYieldStatpp2760GeV 		= (TGraphAsymmErrors*)directoryCombGammapp2760GeV->Get("graphInvXSectionIncGammaStatErr");
		TGraphAsymmErrors*	graphIncGammaYieldSyspp2760GeV 			= (TGraphAsymmErrors*)directoryCombGammapp2760GeV->Get("graphInvXSectionIncGammaSysErr");
    
    TFile* fCombResultspp8TeV= new TFile(nameFilepp8TeV.Data());

		TDirectory* 		directoryCombGammapp8TeV 					      = (TDirectory*)fCombResultspp8TeV->Get("Gamma8TeV");
		TGraphAsymmErrors*	graphDirGammaYieldStatpp8TeV 				= (TGraphAsymmErrors*)directoryCombGammapp8TeV->Get("graphInvXSectionDirGammaNonFitStatErr");
		TGraphAsymmErrors*	graphDirGammaYieldSyspp8TeV 				= (TGraphAsymmErrors*)directoryCombGammapp8TeV->Get("graphInvXSectionDirGammaNonFitSysErr");;
		TGraphAsymmErrors*	graphDirGammaYieldUpperpp8TeV 			= (TGraphAsymmErrors*)directoryCombGammapp8TeV->Get("graphInvXSectionDirGammaNonFitSumErrAr");
		TGraphAsymmErrors*	graphDRStatpp8TeV 						      = (TGraphAsymmErrors*)directoryCombGammapp8TeV->Get("graphRGammaCombNonFitStatErr");
		TGraphAsymmErrors*	graphDRSyspp8TeV 							      = (TGraphAsymmErrors*)directoryCombGammapp8TeV->Get("graphRGammaCombNonFitSysErr");
		TGraphAsymmErrors*	graphIncGammaYieldStatpp8TeV 				= (TGraphAsymmErrors*)directoryCombGammapp8TeV->Get("graphInvXSectionIncGammaStatErr");
		TGraphAsymmErrors*	graphIncGammaYieldSyspp8TeV 				= (TGraphAsymmErrors*)directoryCombGammapp8TeV->Get("graphInvXSectionIncGammaSysErr");
		
		
	TString headerPaper = "*author: ACHARYA \n*reference: CERN-PH-EP-XXXX-XXX : 2018\n*status: Encoded 13 Mar 2018 by N. Schmidt\n*experiment: CERN-LHC-ALICE\n*detector: ALICE\n\n*title: Direct photon production at low transverse momentum in proton-proton collisions at $\\sqrt{s}=2.76$ and 8 TeV\n*comment: CERN-LHC. Measurements of inclusive and direct photon production at mid-rapidity in pp collisions at $\\sqrt{s}=2.76$ and 8 TeV are presented by the ALICE experiment at the LHC. The results are reported in transverse momentum ranges of $0.4<p_{\\rm T}<10$ GeV/c and $0.3<p_{\\rm T}<16$ GeV/c, respectively. Photons are detected with the electromagnetic calorimeter (EMCal) and via reconstruction of $e^+e^-$ pairs from conversions in the ALICE detector material using the central tracking system. For the final measurement of the inclusive photon spectra the results are combined in the overlapping $p_{\\rm T}$ interval of both methods. Direct photon spectra, or their upper limits at 90% C.L. are extracted using the direct photon excess ratio $R_\\gamma$, which quantifies the ratio of inclusive photons over decay photons generated with a decay-photon simulation. An additional hybrid method, combining photons reconstructed from conversions with those identified in the EMCal, is used for the combination of the direct photon excess ratio $R_\\gamma$, as well as the extraction of direct photon spectra or their upper limits. While no significant signal of direct photons is seen over the full $p_{\\rm T}$ range, $R_\\gamma$ for $p_{\\rm T}>7$ GeV/c is at least one $\\sigma$ above unity and consistent with expectations from next-to-leading order pQCD calculations.\n";
    
		cout << headerPaper.Data() << endl;
		
		cout << endl;
		TString headerRAApp2760GeV = "*location: Fig. 4\n*dscomment: Double Ratio RGAMMA in inelastic pp collisions at center-of-mass energy 2.76 TeV\n*reackey: P P --> GAMMA X\n*obskey: RGAMMA\n*qual: RE: P P --> GAMMA X\n*qual: YRAP : 0.0\n*qual: SQRT(S) IN GEV : 2760.0\n*yheader: RGAMMA \n*xheader: PT IN GEV/c\n*data: x : y";
		ProduceHEPDataFileWithUpperLimits(graphDRStatpp2760GeV, graphDRSyspp2760GeV, NULL, newBinsCombpp2760GeV, nPointsTotpp2760GeV, headerRAApp2760GeV);
    
		cout << endl;
		TString headerRAApp8TeV = "*location: Fig. 4\n*dscomment: Double Ratio RGAMMA in inelastic pp collisions at center-of-mass energy 8 TeV\n*reackey: P P --> GAMMA X\n*obskey: RGAMMA\n*qual: RE: P P --> GAMMA X\n*qual: YRAP : 0.0\n*qual: SQRT(S) IN GEV : 8000.0\n*yheader: RGAMMA \n*xheader: PT IN GEV/c\n*data: x : y";
		ProduceHEPDataFileWithUpperLimits(graphDRStatpp8TeV, graphDRSyspp8TeV, NULL, newBinsCombpp8TeV, nPointsTotpp8TeV, headerRAApp8TeV);
    
		TString headerIncpp2760GeV = "*location: Fig. 5\n*dscomment: Invariant differential cross section of inclusive GAMMA produced in inelastic pp collisions at center-of-mass energy 2.76 TeV, the uncertainty of $\\sigma_{MB}$ of 2.5% is not included in the systematic error.\n*reackey: P P --> GAMMA X\n*obskey: E*D3SIG/D3P\n*qual: RE: P P --> GAMMA X\n*qual: YRAP : 0.0\n*qual: SQRT(S) IN GEV : 2760.0\n*yheader: E d**3sigma/dp**3 IN (pb GEV**-2 c**3) \n*xheader: PT IN GEV/c\n*data: x : y";
		ProduceHEPDataFileWithUpperLimits(graphIncGammaYieldStatpp2760GeV, graphIncGammaYieldSyspp2760GeV, NULL, newBinsCombpp2760GeV, nPointsTotpp2760GeV, headerIncpp2760GeV);
    
		cout << endl;
		TString headerIncpp8TeV = "*location: Fig. 5\n*dscomment: Invariant differential cross section of inclusive GAMMA produced in inelastic pp collisions at center-of-mass energy 8 TeV, the uncertainty of $\\sigma_{MB}$ of 2.6% is not included in the systematic error.\n*reackey: P P --> GAMMA X\n*obskey: E*D3SIG/D3P\n*qual: RE: P P --> GAMMA X\n*qual: YRAP : 0.0\n*qual: SQRT(S) IN GEV : 8000.0\n*yheader: E d**3sigma/dp**3 IN (pb GEV**-2 c**3) \n*xheader: PT IN GEV/c\n*data: x : y";
		ProduceHEPDataFileWithUpperLimits(graphIncGammaYieldStatpp8TeV, graphIncGammaYieldSyspp8TeV, NULL, newBinsCombpp8TeV, nPointsTotpp8TeV, headerIncpp8TeV);
	
    TString headerpp2760GeV = "*location: Fig. 5\n*dscomment: Invariant differential cross section of direct GAMMA produced in inelastic pp collisions at center-of-mass energy 2.76 TeV, the uncertainty of $\\sigma_{MB}$ of 2.5% is not included in the systematic error.\n*reackey: P P --> GAMMA X\n*obskey: E*D3SIG/D3P\n*qual: RE: P P --> GAMMA X\n*qual: YRAP : 0.0\n*qual: SQRT(S) IN GEV : 2760.0\n*yheader: E d**3sigma/dp**3 IN (pb GEV**-2 c**3) \n*xheader: PT IN GEV/c\n*data: x : y";
		ProduceHEPDataFileWithUpperLimits(graphDirGammaYieldStatpp2760GeV, graphDirGammaYieldSyspp2760GeV, graphDirGammaYieldUpperpp2760GeV, newBinsCombpp2760GeV, nPointsTotpp2760GeV, headerpp2760GeV,"(90% CL)");
    
		cout << endl;
		TString headerpp8TeV = "*location: Fig. 5\n*dscomment: Invariant differential cross section of direct GAMMA produced in inelastic pp collisions at center-of-mass energy 8 TeV, the uncertainty of $\\sigma_{MB}$ of 2.6% is not included in the systematic error.\n*reackey: P P --> GAMMA X\n*obskey: E*D3SIG/D3P\n*qual: RE: P P --> GAMMA X\n*qual: YRAP : 0.0\n*qual: SQRT(S) IN GEV : 8000.0\n*yheader: E d**3sigma/dp**3 IN (pb GEV**-2 c**3) \n*xheader: PT IN GEV/c\n*data: x : y";
		ProduceHEPDataFileWithUpperLimits(graphDirGammaYieldStatpp8TeV, graphDirGammaYieldSyspp8TeV, graphDirGammaYieldUpperpp8TeV, newBinsCombpp8TeV, nPointsTotpp8TeV, headerpp8TeV,"(90% CL)");
}