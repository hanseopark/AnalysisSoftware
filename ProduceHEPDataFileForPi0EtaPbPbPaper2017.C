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

void ProduceHEPDataFileForPi0EtaPbPbPaper2017( TString nameFilePbPb = "/home/admin1/leardini/Results/FinalMeson25July2017/pdf/2018_03_06/CombineMesonMeasurementsPbPb2760GeVX/CombinedResultsPaperPbPb2760GeV_2018_03_06.root"
    //"/home/admin1/leardini/Results/FinalMeson25July2017/pdf/2018_02_12/CombineMesonMeasurementsPbPb2760GeVX/CombinedResultsPaperPbPb2760GeV_2018_02_12.root"
){

    TString dateForOutput = ReturnDateStringForOutput();


    TString headerPaper = "*author: ACHARYA \n*reference: XXXX \n*reference: ARXIV: XXXXX : \n*reference: CERN-EP-2018-XXX : 2018 \n*status: Encoded 21 Feb 2018 by L. Leardini \n*experiment: CERN-LHC-ALICE \n*detector: ALICE \n*inspireId: \n*cdsId: 000 \n*title: Neutral pion and ETA meson production at central rapidity in Pb-Pb collisions at $\sqrt{s_{\mbox{NN}}}$~=~2.76~TeV \n*comment: CERN-LHC.  Neutral pion and ETA meson production in the transverse momentum range 1~$< p_{\mbox{T}} <$~20~GeV/$c$ have been measured at midrapidity by the ALICE experiment at the Large Hadron Collider (LHC) in central and semi-central Pb-Pb collisions at $\sqrt{s_{\mbox{NN}}}$~=~2.76~TeV. These results extend the $p_{\mbox{T}}$ reach of the previous ALICE PI0 measurements from 12~GeV/$c$ to 20~GeV/$c$ and present the first measurement of the ETA meson in heavy-ion collisions at the LHC. The ETA to PI0 ratio is similar for the two centralities and reaches a plateau value of 0.457~$\pm$~0.013$^{stat}$~$\pm$~0.018$^{syst}$. A suppression of the same magnitude for PI0 and ETA meson production is observed in Pb-Pb collisions with respect to their production in pp collisions scaled by the number of binary nucleon-nucleon collisions. For both mesons, the nuclear modification factor $R_{\mbox{AA}}$ reaches a minimum at $p_{\mbox{T}} \approx $~7~GeV/$c$ and then increases with transverse momentum. The measurements show a stronger suppression with respect to what was observed at lower center-of-mass energies in the $p_{\mbox{T}}$ range 6~$<$~$p_{\mbox{T}}$$<$~10~GeV/$c$. The results are compared with statistical hadronization models at low-$p_{\mbox{T}}$ and with NLO pQCD jet quenching predictions at high-$p_{\mbox{T}}$.\n\n";
    cout << headerPaper.Data() << endl;

    TFile* fileNeutralMesonsCombDataPbPbLHC11h = new TFile(nameFilePbPb.Data());
    TGraphAsymmErrors* graphInvYieldPi0CombPbPb2760GeVStatErr_0010 = (TGraphAsymmErrors*)fileNeutralMesonsCombDataPbPbLHC11h->Get("graphInvYieldPi0CombPbPb2760GeVStatErr_0010");
    TGraphAsymmErrors* graphInvYieldPi0CombPbPb2760GeVSysErr_0010  = (TGraphAsymmErrors*)fileNeutralMesonsCombDataPbPbLHC11h->Get("graphInvYieldPi0CombPbPb2760GeVSysErr_0010");
    TGraphAsymmErrors* graphInvYieldPi0CombPbPb2760GeVStatErr_2050 = (TGraphAsymmErrors*)fileNeutralMesonsCombDataPbPbLHC11h->Get("graphInvYieldPi0CombPbPb2760GeVStatErr_2050");
    TGraphAsymmErrors* graphInvYieldPi0CombPbPb2760GeVSysErr_2050  = (TGraphAsymmErrors*)fileNeutralMesonsCombDataPbPbLHC11h->Get("graphInvYieldPi0CombPbPb2760GeVSysErr_2050");

    // Pi0 invariant yields
    TString headerPi0InvYield_0010 = "*location: Fig. 2, 6 \n*dscomment: Invariant yields of the PI0 meson in the centrality class 0-10%.\n*reackey: PB PB --> PI0 X\n*obskey: DN/DPT**2/DYRAP\n*qual: RE: Pb-Pb --> PI0 X\n*qual: YRAP : 0.0\n*qual: CENTRALITY [PCT] : 0.0 TO 10.0\n*qual: SQRT(S)/NUCLEON [GEV] : 2760.0\n*yheader: 1/(2*PI*Nevt)*(d**2N)/(PT*dPT*dY) IN (GEV/c)**-2\n*xheader: PT [GEV/c]\n*data: x : y";
    ProduceHEPDataFile(graphInvYieldPi0CombPbPb2760GeVStatErr_0010, graphInvYieldPi0CombPbPb2760GeVSysErr_0010, headerPi0InvYield_0010);
    cout << endl;

    TString headerPi0InvYield_2050 = "*location: Fig. 2, 6 \n*dscomment: Invariant yields of the PI0 meson in the centrality class 20-50%.\n*reackey: PB PB --> PI0 X\n*obskey: DN/DPT**2/DYRAP\n*qual: RE: Pb-Pb --> PI0 X\n*qual: YRAP : 0.0\n*qual: CENTRALITY [PCT] : 20.0 TO 50.0\n*qual: SQRT(S)/NUCLEON [GEV] : 2760.0\n*yheader: 1/(2*PI*Nevt)*(d**2N)/(PT*dPT*dY) IN (GEV/c)**-2\n*xheader: PT [GEV/c]\n*data: x : y";
    ProduceHEPDataFile(graphInvYieldPi0CombPbPb2760GeVStatErr_2050, graphInvYieldPi0CombPbPb2760GeVSysErr_2050, headerPi0InvYield_2050);
    cout << endl;

    // Pi0 nuclear modification factor
    TGraphAsymmErrors* graphRAAPi0CombPbPb2760GeVStatErr_0010        = (TGraphAsymmErrors*)fileNeutralMesonsCombDataPbPbLHC11h->Get("graphRAAPi0CombPbPb2760GeVStatErr_0010");
    TGraphAsymmErrors* graphRAAPi0CombPbPb2760GeVSysErr_0010         = (TGraphAsymmErrors*)fileNeutralMesonsCombDataPbPbLHC11h->Get("graphRAAPi0CombPbPb2760GeVSysErr_0010");
    TGraphAsymmErrors* graphRAAPi0CombPbPb2760GeVStatErr_2050        = (TGraphAsymmErrors*)fileNeutralMesonsCombDataPbPbLHC11h->Get("graphRAAPi0CombPbPb2760GeVStatErr_2050");
    TGraphAsymmErrors* graphRAAPi0CombPbPb2760GeVSysErr_2050         = (TGraphAsymmErrors*)fileNeutralMesonsCombDataPbPbLHC11h->Get("graphRAAPi0CombPbPb2760GeVSysErr_2050");

	TString headerPi0RAA_0010 = "*location: Fig. 4, 5, 8 \n *dscomment: Nuclear modification factor RAA of PI0 produced in 0-10% central inelastic Pb-Pb collisions at center-of-mass energy per nucleon 2.76 TeV, TAA uncertainty of 3.25% and uncertainty of \sigma_{inel} of 3.9% included in systematics\n *reackey: PB PB --> PI0 X\n *obskey: RAA\n *qual: RE: PB PB --> PI0 X\n *qual: YRAP : 0.0\n *qual: CENTRALITY [PCT] : 0.0 TO 10.0\n *qual: SQRT(S)/NUCLEON [GEV] : 2760.0\n *yheader:  RAA\n *xheader: PT [GEV/c]\n *data: x : y";
	ProduceHEPDataFile(graphRAAPi0CombPbPb2760GeVStatErr_0010, graphRAAPi0CombPbPb2760GeVSysErr_0010, headerPi0RAA_0010);
    cout << endl;
    TString headerPi0RAA_2050 = "*location: Fig. 4, 8 \n *dscomment: Nuclear modification factor RAA of PI0 produced in 20-50% central inelastic PbPb collisions at center-of-mass energy per nucleon 2.76 TeV, TAA uncertainty of 3.3% and uncertainty of \sigma_{inel} of 3.9% included in systematics\n *reackey: PB PB --> PI0 X\n *obskey: RAA\n *qual: RE: PB PB --> PI0 X\n *qual: YRAP : 0.0\n *qual: CENTRALITY [PCT] : 20.0 TO 50.0\n *qual: SQRT(S)/NUCLEON [GEV] : 2760.0\n *yheader:  RAA\n *xheader: PT [GEV/c]\n *data: x : y";
	ProduceHEPDataFile(graphRAAPi0CombPbPb2760GeVStatErr_2050, graphRAAPi0CombPbPb2760GeVSysErr_2050, headerPi0RAA_2050);
    cout << endl;


    // Eta invariant yields
    TGraphAsymmErrors* graphInvYieldEtaCombPbPb2760GeVStatErr_0010   = (TGraphAsymmErrors*)fileNeutralMesonsCombDataPbPbLHC11h->Get("graphInvYieldEtaCombPbPb2760GeVStatErr_0010");
    TGraphAsymmErrors* graphInvYieldEtaCombPbPb2760GeVSysErr_0010    = (TGraphAsymmErrors*)fileNeutralMesonsCombDataPbPbLHC11h->Get("graphInvYieldEtaCombPbPb2760GeVSysErr_0010");
    TGraphAsymmErrors* graphInvYieldEtaCombPbPb2760GeVStatErr_2050   = (TGraphAsymmErrors*)fileNeutralMesonsCombDataPbPbLHC11h->Get("graphInvYieldEtaCombPbPb2760GeVStatErr_2050");
    TGraphAsymmErrors* graphInvYieldEtaCombPbPb2760GeVSysErr_2050    = (TGraphAsymmErrors*)fileNeutralMesonsCombDataPbPbLHC11h->Get("graphInvYieldEtaCombPbPb2760GeVSysErr_2050");

    TString headerEtaInvYield_0010 = "*location: Fig. 2, 6 \n*dscomment: Invariant yields of the ETA meson in the centrality class 0-10%.\n*reackey: PB PB --> ETA X\n*obskey: DN/DPT**2/DYRAP\n*qual: RE: Pb-Pb --> ETA X\n*qual: YRAP : 0.0\n*qual: CENTRALITY [PCT] : 0.0 TO 10.0\n*qual: SQRT(S)/NUCLEON [GEV] : 2760.0\n*yheader: 1/(2*PI*Nevt)*(d**2N)/(PT*dPT*dY) IN (GEV/c)**-2\n*xheader: PT [GEV/c]\n*data: x : y";
    ProduceHEPDataFile(graphInvYieldEtaCombPbPb2760GeVStatErr_0010, graphInvYieldEtaCombPbPb2760GeVSysErr_0010, headerEtaInvYield_0010);
    cout << endl;
    TString headerEtaInvYield_2050 = "*location: Fig. 2, 6 \n*dscomment: Invariant yields of the ETA meson in the centrality class 20-50%.\n*reackey: PB PB --> ETA X\n*obskey: DN/DPT**2/DYRAP\n*qual: RE: Pb-Pb --> ETA X\n*qual: YRAP : 0.0\n*qual: CENTRALITY [PCT] : 20.0 TO 50.0\n*qual: SQRT(S)/NUCLEON [GEV] : 2760.0\n*yheader: 1/(2*PI*Nevt)*(d**2N)/(PT*dPT*dY) IN (GEV/c)**-2\n*xheader: PT [GEV/c]\n*data: x : y";
    ProduceHEPDataFile(graphInvYieldEtaCombPbPb2760GeVStatErr_2050, graphInvYieldEtaCombPbPb2760GeVSysErr_2050, headerEtaInvYield_2050);
    cout << endl;


    // Eta nuclear modification factor
    TGraphAsymmErrors* graphRAAEtaCombPbPb2760GeVStatErr_0010        = (TGraphAsymmErrors*)fileNeutralMesonsCombDataPbPbLHC11h->Get("graphRAAEtaCombPbPb2760GeVStatErr_0010");
    TGraphAsymmErrors* graphRAAEtaCombPbPb2760GeVSysErr_0010         = (TGraphAsymmErrors*)fileNeutralMesonsCombDataPbPbLHC11h->Get("graphRAAEtaCombPbPb2760GeVSysErr_0010");
    TGraphAsymmErrors* graphRAAEtaCombPbPb2760GeVStatErr_2050        = (TGraphAsymmErrors*)fileNeutralMesonsCombDataPbPbLHC11h->Get("graphRAAEtaCombPbPb2760GeVStatErr_2050");
    TGraphAsymmErrors* graphRAAEtaCombPbPb2760GeVSysErr_2050         = (TGraphAsymmErrors*)fileNeutralMesonsCombDataPbPbLHC11h->Get("graphRAAEtaCombPbPb2760GeVSysErr_2050");

	TString headerEtaRAA_0010 = "*location: Fig. 4, 5, 8 \n *dscomment: Nuclear modification factor RAA of ETA produced in 0-10% central inelastic Pb-Pb collisions at center-of-mass energy per nucleon 2.76 TeV, TAA uncertainty of 3.25% and uncertainty of \sigma_{inel} of 3.9% included in systematics\n *reackey: PB PB --> ETA X\n *obskey: RAA\n *qual: RE: PB PB --> ETA X\n *qual: YRAP : 0.0\n *qual: CENTRALITY [PCT] : 0.0 TO 10.0\n *qual: SQRT(S)/NUCLEON [GEV] : 2760.0\n *yheader:  RAA\n *xheader: PT [GEV/c]\n *data: x : y";
	ProduceHEPDataFile(graphRAAEtaCombPbPb2760GeVStatErr_0010, graphRAAEtaCombPbPb2760GeVSysErr_0010, headerEtaRAA_0010);
    cout << endl;
    TString headerEtaRAA_2050 = "*location: Fig. 4, 8 \n *dscomment: Nuclear modification factor RAA of ETA produced in 20-50% central inelastic PbPb collisions at center-of-mass energy per nucleon 2.76 TeV, TAA uncertainty of 3.3% and uncertainty of \sigma_{inel} of 3.9% included in systematics\n *reackey: PB PB --> ETA X\n *obskey: RAA\n *qual: RE: PB PB --> ETA X\n *qual: YRAP : 0.0\n *qual: CENTRALITY [PCT] : 20.0 TO 50.0\n *qual: SQRT(S)/NUCLEON [GEV] : 2760.0\n *yheader:  RAA\n *xheader: PT [GEV/c]\n *data: x : y";
	ProduceHEPDataFile(graphRAAEtaCombPbPb2760GeVStatErr_2050, graphRAAEtaCombPbPb2760GeVSysErr_2050, headerEtaRAA_2050);
    cout << endl;


    // Eta/Pi0 ratio
    TGraphAsymmErrors* graphEtaToPi0CombPbPb2760GeVStatErr_0010      = (TGraphAsymmErrors*)fileNeutralMesonsCombDataPbPbLHC11h->Get("graphEtaToPi0RatioCombPbPb2760GeVStatErr_0010");
    TGraphAsymmErrors* graphEtaToPi0CombPbPb2760GeVSysErr_0010       = (TGraphAsymmErrors*)fileNeutralMesonsCombDataPbPbLHC11h->Get("graphEtaToPi0RatioCombPbPb2760GeVSysErr_0010");
    TGraphAsymmErrors* graphEtaToPi0CombPbPb2760GeVStatErr_2050      = (TGraphAsymmErrors*)fileNeutralMesonsCombDataPbPbLHC11h->Get("graphEtaToPi0RatioCombPbPb2760GeVStatErr_2050");
    TGraphAsymmErrors* graphEtaToPi0CombPbPb2760GeVSysErr_2050       = (TGraphAsymmErrors*)fileNeutralMesonsCombDataPbPbLHC11h->Get("graphEtaToPi0RatioCombPbPb2760GeVSysErr_2050");

    graphEtaToPi0CombPbPb2760GeVStatErr_0010->RemovePoint(0);
    graphEtaToPi0CombPbPb2760GeVSysErr_0010->RemovePoint(0);
    TString headerEtaToPi0_0010 ="*location: Fig.3, 7 \n*dscomment: Ratio of the ETA to PI0 invariant yields in the centrality class 0-10%.\n*reackey: PB PB --> PI0 X\n*reackey: PB PB --> ETA X\n*obskey: SIG/SIG\n*qual: RE(Q=ETA) : PB PB --> ETA X\n*qual: RE(Q=PI0) : PB PB --> PI0 X\n*qual: SQRT(S)/NUCLEON [GEV] : 2760.0\n*qual: YRAP : 0.0\n*qual: CENTRALITY [PCT] : 0.0 TO 10.0\n*yheader: SIG(Q=ETA)/SIG(Q=PI0)\n*xheader: PT [GEV/c]\n*data: x : y";
    ProduceHEPDataFile(graphEtaToPi0CombPbPb2760GeVStatErr_0010, graphEtaToPi0CombPbPb2760GeVSysErr_0010, headerEtaToPi0_0010);
    cout << endl;
    graphEtaToPi0CombPbPb2760GeVStatErr_2050->RemovePoint(0);
    graphEtaToPi0CombPbPb2760GeVSysErr_2050->RemovePoint(0);
    TString headerEtaToPi0_2050 ="*location: Fig.3, 7 \n*dscomment: Ratio of the ETA to PI0 invariant yields in the centrality class 20-50%.\n*reackey: PB PB --> PI0 X\n*reackey: PB PB --> ETA X\n*obskey: SIG/SIG\n*qual: RE(Q=ETA) : PB PB --> ETA X\n*qual: RE(Q=PI0) : PB PB --> PI0 X\n*qual: SQRT(S)/NUCLEON [GEV] : 2760.0\n*qual: YRAP : 0.0\n*qual: CENTRALITY [PCT] : 20.0 TO 50.0\n*yheader: SIG(Q=ETA)/SIG(Q=PI0)\n*xheader: PT [GEV/c]\n*data: x : y";
    ProduceHEPDataFile(graphEtaToPi0CombPbPb2760GeVStatErr_2050, graphEtaToPi0CombPbPb2760GeVSysErr_2050, headerEtaToPi0_2050);
    cout << endl;

}