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

void ProduceHEPDataFileForPi0EtaHighPtPaper2017( TString nameFilePP = ""){

    gROOT->Reset();	
    gROOT->SetStyle("Plain");
    
    StyleSettingsThesis();	
    SetPlotStyle();

    TString dateForOutput = ReturnDateStringForOutput();
    

    TString collisionSystemPP = "pp #sqrt{#it{s}} = 2.76 TeV";		
    
    TFile* fileNeutralPionCombDataPP = new TFile(nameFilePP.Data());
    TGraphAsymmErrors* graphInvXSectionPi0Comb2760GeVStatErr= (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("Pi02.76TeV/graphInvCrossSectionPi0Comb2760GeVAStatErr");
    TGraphAsymmErrors* graphInvXSectionPi0Comb2760GeVSysErr = (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("Pi02.76TeV/graphInvCrossSectionPi0Comb2760GeVASysErr");
    TGraphAsymmErrors* graphInvYieldPi0Comb2760GeVStatErr   = (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("Pi02.76TeV/graphInvYieldINELPi0Comb2760GeVAStatErr");
    TGraphAsymmErrors* graphInvYieldPi0Comb2760GeVSysErr    = (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("Pi02.76TeV/graphInvYieldINELPi0Comb2760GeVASysErr");
    TGraphAsymmErrors* graphInvXSectionEtaComb2760GeVStatErr= (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("Eta2.76TeV/graphInvCrossSectionEtaComb2760GeVAStatErr");
    TGraphAsymmErrors* graphInvXSectionEtaComb2760GeVSysErr = (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("Eta2.76TeV/graphInvCrossSectionEtaComb2760GeVASysErr");
    TGraphAsymmErrors* graphInvYieldEtaComb2760GeVStatErr   = (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("Eta2.76TeV/graphInvYieldINELEtaComb2760GeVAStatErr");
    TGraphAsymmErrors* graphInvYieldEtaComb2760GeVSysErr    = (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("Eta2.76TeV/graphInvYieldINELEtaComb2760GeVASysErr");
    TGraphAsymmErrors* graphEtaToPi0Comb2760GeVStatErr      = (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("Eta2.76TeV/graphRatioEtaToPi0Comb2760GeVStatErr");
    TGraphAsymmErrors* graphEtaToPi0Comb2760GeVSysErr       = (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("Eta2.76TeV/graphRatioEtaToPi0Comb2760GeVSysErr");
    
    TString headerPaper = "*author: ACHARYA \n *reference: ARXIV:1702.XXXX : 2017\n *reference:  CERN-EP-2017-019 : 2017\n *reference: \n *status: Encoded 30 Jan 2017 by F. Bock\n *experiment: CERN-LHC-ALICE\n *detector: ALICE\n *inspireId: \n *cdsId: 2243252\n *title: Production of PI0 and ETA mesons up to high transverse momentum in pp collisions at 2.76 TeV\n *comment: CERN-LHC.  The invariant differential cross sections for inclusive PI0 and ETA mesons at midrapidity were measured in pp collisions at $\sqrt{s}=2.76$ TeV for transverse momenta 0.4<pT<40~GeV/c and 0.6<pT<20~GeV/c, respectively, using the ALICE detector. This large range in pT was achieved by combining various analysis techniques and different triggers involving the electromagnetic calorimeter~(EMCal). In particular, a new single-cluster, shower-shape based method was developed for the identification of high-pT neutral pions, which exploits that the showers originating from their decay photons overlap in the EMCal. The measured cross sections are found to exhibit a similar power-law behavior with an exponent of n = 6.2. Next-to-leading order perturbative QCD calculations and generator-level simulations with Pythia 8.2 describe the cross sections to about 30% for the PI0, and between 30-50% for the ETA meson. The new data can therefore be used to further improve the theoretical description of PI0 and ETA meson production\n\n";	
        
    cout << headerPaper.Data() << endl;
    
    TString headerPi0Yield = "*location: Fig. 9\n *dscomment: Invariant differential yields of PI0 produced in inelastic pp collisions at center-of-mass energy 2.76 TeV, the normalization uncertainties of 5.7% are not included in the systematic error \n *reackey: P P --> PI0 X\n *obskey: DN/DPT**2/DYRAP\n *qual: RE: P P --> PI0 X\n *qual: YRAP : 0.0\n *qual: SQRT(S) IN GEV : 2760.0\n *yheader: 1/(2*PI*Nevt)*(d**2N)/(PT*dPT*dY) IN (GEV/c)**-2\n *xheader: PT IN GEV/c\n *data: x : y";
    ProduceHEPDataFile(graphInvYieldPi0Comb2760GeVStatErr, graphInvYieldPi0Comb2760GeVSysErr, headerPi0Yield);
    cout << endl;

    TString headerEtaYield = "*location: Fig. 9\n *dscomment: Invariant differential yields of ETA produced in inelastic pp collisions at center-of-mass energy 2.76 TeV, the normalization uncertainties of 5.7% are not included in the systematic error \n *reackey: P P --> ETA X\n *obskey: DN/DPT**2/DYRAP\n *qual: RE: P P --> ETA X\n *qual: YRAP : 0.0\n *qual: SQRT(S) IN GEV : 2760.0\n *yheader: 1/(2*PI*Nevt)*(d**2N)/(PT*dPT*dY) IN (GEV/c)**-2\n *xheader: PT IN GEV/c\n *data: x : y";
    ProduceHEPDataFile(graphInvYieldEtaComb2760GeVStatErr, graphInvYieldEtaComb2760GeVSysErr, headerEtaYield);
    cout << endl;
    
    TString headerPi0XSec = "*location: Fig. 9\n *dscomment: Invariant differential cross section of PI0 produced in inelastic pp collisions at center-of-mass energy 2.76 TeV, the uncertainty of \sigma_{MB} of 2.5% is not included in the systematic error \n *reackey: P P --> PI0 X\n *obskey: E*D3SIG/D3P\n *qual: RE: P P --> PI0 X\n *qual: YRAP : 0.0\n *qual: SQRT(S) IN GEV : 2760.0\n *yheader: E d**3sigma/dp**3 IN (pb GEV**-2 c**3)\n *xheader: PT IN GEV/c\n *data: x : y";
	ProduceHEPDataFile(graphInvXSectionPi0Comb2760GeVStatErr, graphInvXSectionPi0Comb2760GeVSysErr, headerPi0XSec);
    cout << endl;
    
    TString headerEtaXSec = "*location: Fig. 9\n *dscomment: Invariant differential cross section of ETA produced in inelastic pp collisions at center-of-mass energy 2.76 TeV, the uncertainty of \sigma_{MB} of 2.5% is not included in the systematic error \n *reackey: P P --> ETA X\n *obskey: E*D3SIG/D3P\n *qual: RE: P P --> ETA X\n *qual: YRAP : 0.0\n *qual: SQRT(S) IN GEV : 2760.0\n *yheader: E d**3sigma/dp**3 IN (pb GEV**-2 c**3)\n *xheader: PT IN GEV/c\n *data: x : y";
	ProduceHEPDataFile(graphInvXSectionEtaComb2760GeVStatErr, graphInvXSectionEtaComb2760GeVSysErr, headerEtaXSec);
    cout << endl;
    
    TString headerEtaToPi0 ="*location: Fig.10\n *dscomment: The measured ratio of cross sections for inclusive ETA to PI0 production at a centre-of-mass energy of 2.76 TeV.\n *reackey: P P --> PI0 X\n *reackey: P P --> ETA X\n *obskey: SIG/SIG\n *qual: RE(Q=ETA) : P P --> ETA X\n *qual: RE(Q=PI0) : P P --> PI0 X\n *qual: SQRT(S) IN GEV : 2760.0\n *qual: YRAP : 0.0\n *yheader: SIG(Q=ETA)/SIG(Q=PI0)\n *xheader: PT IN GEV\n *data: x : y";
    ProduceHEPDataFile(graphEtaToPi0Comb2760GeVStatErr, graphEtaToPi0Comb2760GeVSysErr, headerEtaToPi0);
    cout << endl;
   
}