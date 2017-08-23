/****************************************************************************************************************************
****** 		provided by Gamma Conversion Group, PWGGA, 													*****
******		Daniel Mühlheim, d.muehlheim@cern.ch														*****
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

void ProduceHEPDataFileForPi0Eta8TeVPaper2017( TString nameFilePP = ""){

    gROOT->Reset();
    gROOT->SetStyle("Plain");

    StyleSettingsThesis();
    SetPlotStyle();

    TString dateForOutput = ReturnDateStringForOutput();


    TString collisionSystemPP = "pp #sqrt{#it{s}} = 8 TeV";

    TFile* fileNeutralPionCombDataPP = new TFile(nameFilePP.Data());
    TGraphAsymmErrors* graphInvXSectionPi0Comb8TeVStatErr= (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("Pi08TeV/graphInvCrossSectionPi0Comb8TeVAStatErr");
    TGraphAsymmErrors* graphInvXSectionPi0Comb8TeVSysErr = (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("Pi08TeV/graphInvCrossSectionPi0Comb8TeVASysErr");
    TGraphAsymmErrors* graphInvXSectionEtaComb8TeVStatErr= (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("Eta8TeV/graphInvCrossSectionEtaComb8TeVAStatErr");
    TGraphAsymmErrors* graphInvXSectionEtaComb8TeVSysErr = (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("Eta8TeV/graphInvCrossSectionEtaComb8TeVASysErr");
    TGraphAsymmErrors* graphEtaToPi0Comb8TeVStatErr      = (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("Eta8TeV/graphRatioEtaToPi0Comb8TeVStatErr");
    TGraphAsymmErrors* graphEtaToPi0Comb8TeVSysErr       = (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("Eta8TeV/graphRatioEtaToPi0Comb8TeVSysErr");

    TString headerPaper = "*author: ACHARYA \n*reference: XXXXXX \n*reference: ARXIV:XXXX.XXXX : 2017\n*reference:  CERN-EP-2017-216 : 2017\n*status: Encoded 21 Aug 2017 by D. Mühlheim\n*experiment: CERN-LHC-ALICE\n*detector: ALICE\n*inspireId: \n*cdsId: 2280315\n*title: $\pi^0$ and $\eta$ meson production in proton-proton collisions at $\sqrt{s}=8$ TeV \n*comment: CERN-LHC. An invariant differential cross section measurement of inclusive $\pi^0$ and $\eta$ meson production at mid-rapidity in pp collisions at $\sqrt{s}=8$ TeV was carried out by the ALICE experiment at the LHC. The spectra of neutral mesons $\pi^0$ and $\eta$ were measured in transverse momentum ranges of $0.3<p_{\rm T}<35$ GeV/$c$ and $0.5<p_{\rm T}<35$ GeV/$c$, respectively. Next-to-leading order perturbative QCD calculations using fragmentation functions DSS14 for $\pi^0$ and AESSS for $\eta$ overestimate the cross sections of both neutral mesons, but agree with the measured $\eta$/$\pi^0$ ratio within uncertainties. The results are also compared with PYTHIA8.2 predictions for which the Monash2013 tune yields the best agreement with the measured neutral meson spectra. The measurements confirm a universal behavior of the $\eta$/$\pi^0$ ratio observed for NA27, PHENIX and ALICE data for pp collisions from $\sqrt{s}=27.5$ GeV to $\sqrt{s}=8$ TeV within experimental uncertainties. A relation between the $\pi^0$ and $\eta$ production cross sections for pp collisions at $\sqrt{s}=8$ TeV is given by $m_{\rm T}$ scaling for $p_{\rm T}>3.5$ GeV/$c$. However, a deviation from this empirical scaling law is observed for transverse momenta below $p_{\rm T}<3.5$ GeV/$c$ in the $\eta$/$\pi^0$ ratio with a significance of $6.2\sigma$.\n\n";

    cout << headerPaper.Data() << endl;

    TString headerPi0XSec = "*location: Fig. 7\n*dscomment: Invariant differential cross section of $\pi^0$ produced in inelastic pp collisions at center-of-mass energy 8 TeV, the uncertainty of $\sigma_{MB}$ of 2.6% is not included in the systematic error. \n*reackey: P P --> PI0 X\n*obskey: E*D3SIG/D3P\n*qual: RE: P P --> PI0 X\n*qual: YRAP : 0.0\n*qual: SQRT(S) IN GEV : 8000.0\n*yheader: E d**3sigma/dp**3 IN (pb GEV**-2 c**3)\n*xheader: PT IN GEV/c\n*data: x : y";
    ProduceHEPDataFile(graphInvXSectionPi0Comb8TeVStatErr, graphInvXSectionPi0Comb8TeVSysErr, headerPi0XSec);
    cout << endl;

    TString headerEtaXSec = "*location: Fig. 7\n*dscomment: Invariant differential cross section of $\eta$ produced in inelastic pp collisions at center-of-mass energy 8 TeV, the uncertainty of $\sigma_{MB}$ of 2.6% is not included in the systematic error. \n*reackey: P P --> ETA X\n*obskey: E*D3SIG/D3P\n*qual: RE: P P --> ETA X\n*qual: YRAP : 0.0\n*qual: SQRT(S) IN GEV : 8000.0\n*yheader: E d**3sigma/dp**3 IN (pb GEV**-2 c**3)\n*xheader: PT IN GEV/c\n*data: x : y";
    ProduceHEPDataFile(graphInvXSectionEtaComb8TeVStatErr, graphInvXSectionEtaComb8TeVSysErr, headerEtaXSec);
    cout << endl;

    TString headerEtaToPi0 ="*location: Fig. 8\n*dscomment: The measured ratio of cross sections for inclusive $\eta$ to $\pi^0$ production at a centre-of-mass energy of 8 TeV. \n*reackey: P P --> PI0 X\n*reackey: P P --> ETA X\n*obskey: SIG/SIG\n*qual: RE(Q=ETA) : P P --> ETA X\n*qual: RE(Q=PI0) : P P --> PI0 X\n*qual: SQRT(S) IN GEV : 8000.0\n*qual: YRAP : 0.0\n*yheader: SIG(Q=ETA)/SIG(Q=PI0)\n*xheader: PT IN GEV\n*data: x : y";
    ProduceHEPDataFile(graphEtaToPi0Comb8TeVStatErr, graphEtaToPi0Comb8TeVSysErr, headerEtaToPi0);
    cout << endl;

}
