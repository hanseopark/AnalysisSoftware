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
#include "TFitResultPtr.h"
#include "TFitResult.h"


void WorldDataPi0EtapPb( TString fileNameAlice5TeV = ""){

    // Create 1 graph containing all point above 4 GeV
    vector<Double_t> *valuesEtaToPi0            = new vector<Double_t>[4];     // iCent x iParticle x nMeasurements matrix with theory curves
    Int_t currentNumberOfPoint                  = 0;

    Double_t pt[100];
    Double_t value[100];
    Double_t totErr[100];
    Double_t xErr[100];
    Double_t statErr[100];
    Double_t sysErr[100];
    Double_t sysErrB[100];
    Double_t a = 0;
    Double_t b = 0;
    Double_t c = 0;
    Double_t d = 0;

//     @article{Agakishiev:1998mw,
//       author         = "Agakishiev, G. and others",
//       title          = "{Neutral meson production in p Be and p Au collisions at
//                         450-GeV beam energy}",
//       journal        = "Eur. Phys. J.",
//       volume         = "C4",
//       year           = "1998",
//       pages          = "249-257",
//       doi            = "10.1007/s100529800804",
//       SLACcitation   = "%%CITATION = EPHJA,C4,249;%%"
//     }

    ifstream AgakishievpAu29100MeV;
    AgakishievpAu29100MeV.open("OtherExperiments/agakichiev98_pAu_plab450GeV_tapsceres_eta_pi0_ratio.txt");
    cout << "agakichiev98_pAu_plab450GeV_tapsceres_eta_pi0_ratio.txt" << endl;
    Int_t lines = 0;

    while(!AgakishievpAu29100MeV.eof()){
        AgakishievpAu29100MeV >> pt[lines] >> value[lines] >> totErr[lines] ;
        xErr[lines] = 0.;
        lines++;
    }
    AgakishievpAu29100MeV.close();

    TGraph *AgakishievpAu29100MeVGraph = new TGraphErrors(lines-1,pt,value,xErr,totErr);
    AgakishievpAu29100MeVGraph->SetName("AgakishievpAu29100MeV");
    AgakishievpAu29100MeVGraph->SetTitle("p-Au (#sqrt{#it{s}_{_{NN}}}= 29.1 GeV)");
    AgakishievpAu29100MeVGraph->Print();
    for (Int_t i= 0; i< AgakishievpAu29100MeVGraph->GetN(); i++){
        if (AgakishievpAu29100MeVGraph->GetX()[i] > 4){
            valuesEtaToPi0[0].push_back(AgakishievpAu29100MeVGraph->GetX()[i]);
            valuesEtaToPi0[1].push_back(AgakishievpAu29100MeVGraph->GetY()[i]);
            valuesEtaToPi0[2].push_back(AgakishievpAu29100MeVGraph->GetEY()[i]);
            valuesEtaToPi0[3].push_back(AgakishievpAu29100MeVGraph->GetEY()[i]);
            currentNumberOfPoint++;
        }
    }
    cout << "Number of pT bins: "<< currentNumberOfPoint <<  endl;

    ifstream AgakishievpBe29100MeV;
    AgakishievpBe29100MeV.open("OtherExperiments/agakichiev98_pBe_plab450GeV_tapsceres_eta_pi0_ratio.txt");
    cout << "agakichiev98_pBe_plab450GeV_tapsceres_eta_pi0_ratio.txt" << endl;
    lines = 0;

    while(!AgakishievpBe29100MeV.eof()){
        AgakishievpBe29100MeV >> pt[lines] >> value[lines] >> totErr[lines] ;
        xErr[lines] = 0.;
        lines++;
    }
    AgakishievpBe29100MeV.close();

    TGraph *AgakishievpBe29100MeVGraph = new TGraphErrors(lines-1,pt,value,xErr,totErr);
    AgakishievpBe29100MeVGraph->SetName("AgakishievpBe29100MeV");
    AgakishievpBe29100MeVGraph->SetTitle("p-Be (#sqrt{#it{s}_{_{NN}}}= 29.1 GeV)");
    AgakishievpBe29100MeVGraph->Print();
    for (Int_t i= 0; i< AgakishievpBe29100MeVGraph->GetN(); i++){
        if (AgakishievpBe29100MeVGraph->GetX()[i] > 4){
            valuesEtaToPi0[0].push_back(AgakishievpBe29100MeVGraph->GetX()[i]);
            valuesEtaToPi0[1].push_back(AgakishievpBe29100MeVGraph->GetY()[i]);
            valuesEtaToPi0[2].push_back(AgakishievpBe29100MeVGraph->GetEY()[i]);
            valuesEtaToPi0[3].push_back(AgakishievpBe29100MeVGraph->GetEY()[i]);
            currentNumberOfPoint++;
        }
    }
    cout << "Number of pT bins: "<< currentNumberOfPoint <<  endl;

//     @article{Adler:2006bv,
//       author         = "Adler, S. S. and others",
//       title          = "{High transverse momentum $\eta$ meson production in $p^+
//                         p$, $d^+$ Au and Au+Au collisions at $S(NN) ^{(1/2)}$ =
//                         200-GeV}",
//       collaboration  = "PHENIX",
//       journal        = "Phys. Rev.",
//       volume         = "C75",
//       year           = "2007",
//       pages          = "024909",
//       doi            = "10.1103/PhysRevC.75.024909",
//       eprint         = "nucl-ex/0611006",
//       archivePrefix  = "arXiv",
//       primaryClass   = "nucl-ex",
//       SLACcitation   = "%%CITATION = NUCL-EX/0611006;%%"
//     }

    ifstream PhenixdAu200GeV;
    PhenixdAu200GeV.open("OtherExperiments/phenix_dAu_200GeV_eta_pi0_ratio.txt");
    cout << "phenix_dAu_200GeV_eta_pi0_ratio.txt" << endl;
    lines = 0;


    while(!PhenixdAu200GeV.eof()){
        PhenixdAu200GeV >> pt[lines] >> value[lines] >> statErr[lines] >> sysErr[lines] >> sysErrB[lines];
        totErr[lines] = TMath::Sqrt(statErr[lines]*statErr[lines]+ sysErr[lines]*sysErr[lines] + sysErrB[lines]*sysErrB[lines]);
        xErr[lines] = 0.;
        lines++;
    }
    PhenixdAu200GeV.close();
    TGraph *PhenixdAu200GeVGraph = new TGraphErrors(lines-1,pt,value,xErr,totErr);
    PhenixdAu200GeVGraph->SetName("PhenixdAu200GeV");
    PhenixdAu200GeVGraph->SetTitle("d-AU (#sqrt{#it{s_{_{NN}}}}= 200 GeV)");
    PhenixdAu200GeVGraph->Print();
    for (Int_t i= 0; i< PhenixdAu200GeVGraph->GetN(); i++){
        if (PhenixdAu200GeVGraph->GetX()[i] > 4){
            valuesEtaToPi0[0].push_back(PhenixdAu200GeVGraph->GetX()[i]);
            valuesEtaToPi0[1].push_back(PhenixdAu200GeVGraph->GetY()[i]);
            valuesEtaToPi0[2].push_back(PhenixdAu200GeVGraph->GetEY()[i]);
            valuesEtaToPi0[3].push_back(PhenixdAu200GeVGraph->GetEY()[i]);
            currentNumberOfPoint++;
        }
    }
    cout << "Number of pT bins: "<< currentNumberOfPoint <<  endl;


//     @article{Povlis:1983bb,
//         author         = "Povlis, J. and others",
//         title          = "{Nuclear Enhancement of pi0 and eta mesons Produced at
//         Large Transverse Momenta}",
//         booktitle      = "{11th International Symposium on Lepton and Photon
//         Interactions at High Energies Ithaca, New York, August
//         4-9, 1983}",
//         journal        = "Phys. Rev. Lett.",
//         volume         = "51",
//         year           = "1983",
//         pages          = "967",
//         doi            = "10.1103/PhysRevLett.51.967",
//         reportNumber   = "FERMILAB-PUB-83-116-E, PRINT-83-0769",
//         SLACcitation   = "%%CITATION = PRLTA,51,967;%%"
//     }
//      sNN = 19.4 GeV
//     OtherExperiments/povlis83_pC_plab200GeV_fnal629_eta_pi0_ratio.txt
//     OtherExperiments/povlis83_piC_plab200GeV_fnal629_eta_pi0_ratio.txt
//     OtherExperiments/povlis83_pAl_plab200GeV_fnal629_eta_pi0_ratio.txt
//     OtherExperiments/povlis83_pBe_plab200GeV_fnal629_eta_pi0_ratio.txt
    ifstream PovlispBe19400MeV;
    PovlispBe19400MeV.open("OtherExperiments/povlis83_pBe_plab200GeV_fnal629_eta_pi0_ratio.txt");
    cout << "povlis83_pBe_plab200GeV_fnal629_eta_pi0_ratio.txt" << endl;
    lines = 0;

    while(!PovlispBe19400MeV.eof()){
        PovlispBe19400MeV >> pt[lines] >> value[lines] >> totErr[lines] ;
        xErr[lines] = 0.;
        lines++;
    }
    PovlispBe19400MeV.close();

    TGraph *PovlispBe19400MeVGraph = new TGraphErrors(lines-1,pt,value,xErr,totErr);
    PovlispBe19400MeVGraph->SetName("PovlispBe19400MeV");
    PovlispBe19400MeVGraph->SetTitle("p-Be (#sqrt{#it{s}_{_{NN}}}= 19.4 GeV)");
    PovlispBe19400MeVGraph->Print();
    for (Int_t i= 0; i< PovlispBe19400MeVGraph->GetN(); i++){
        if (PovlispBe19400MeVGraph->GetX()[i] > 4){
            valuesEtaToPi0[0].push_back(PovlispBe19400MeVGraph->GetX()[i]);
            valuesEtaToPi0[1].push_back(PovlispBe19400MeVGraph->GetY()[i]);
            valuesEtaToPi0[2].push_back(PovlispBe19400MeVGraph->GetEY()[i]);
            valuesEtaToPi0[3].push_back(PovlispBe19400MeVGraph->GetEY()[i]);
            currentNumberOfPoint++;
        }
    }
    cout << "Number of pT bins: "<< currentNumberOfPoint <<  endl;

    ifstream PovlispAl19400MeV;
    PovlispAl19400MeV.open("OtherExperiments/povlis83_pAl_plab200GeV_fnal629_eta_pi0_ratio.txt");
    cout << "povlis83_pAl_plab200GeV_fnal629_eta_pi0_ratio.txt" << endl;
    lines = 0;

    while(!PovlispAl19400MeV.eof()){
        PovlispAl19400MeV >> pt[lines] >> value[lines] >> totErr[lines] ;
        xErr[lines] = 0.;
        lines++;
    }
    PovlispAl19400MeV.close();

    TGraph *PovlispAl19400MeVGraph = new TGraphErrors(lines-1,pt,value,xErr,totErr);
    PovlispAl19400MeVGraph->SetName("PovlispAl19400MeV");
    PovlispAl19400MeVGraph->SetTitle("p-Al (#sqrt{#it{s}_{_{NN}}}= 19.4 GeV)");
    PovlispAl19400MeVGraph->Print();
    for (Int_t i= 0; i< PovlispAl19400MeVGraph->GetN(); i++){
        if (PovlispAl19400MeVGraph->GetX()[i] > 4){
            valuesEtaToPi0[0].push_back(PovlispAl19400MeVGraph->GetX()[i]);
            valuesEtaToPi0[1].push_back(PovlispAl19400MeVGraph->GetY()[i]);
            valuesEtaToPi0[2].push_back(PovlispAl19400MeVGraph->GetEY()[i]);
            valuesEtaToPi0[3].push_back(PovlispAl19400MeVGraph->GetEY()[i]);
            currentNumberOfPoint++;
        }
    }
    cout << "Number of pT bins: "<< currentNumberOfPoint <<  endl;

    ifstream PovlispC19400MeV;
    PovlispC19400MeV.open("OtherExperiments/povlis83_pC_plab200GeV_fnal629_eta_pi0_ratio.txt");
    cout << "povlis83_pC_plab200GeV_fnal629_eta_pi0_ratio.txt" << endl;
    lines = 0;

    while(!PovlispC19400MeV.eof()){
        PovlispC19400MeV >> pt[lines] >> value[lines] >> totErr[lines] ;
        xErr[lines] = 0.;
        lines++;
    }
    PovlispC19400MeV.close();

    TGraph *PovlispC19400MeVGraph = new TGraphErrors(lines-1,pt,value,xErr,totErr);
    PovlispC19400MeVGraph->SetName("PovlispC19400MeV");
    PovlispC19400MeVGraph->SetTitle("p-C (#sqrt{#it{s}_{_{NN}}}= 19.4 GeV)");
    PovlispC19400MeVGraph->Print();
    for (Int_t i= 0; i< PovlispC19400MeVGraph->GetN(); i++){
        if (PovlispC19400MeVGraph->GetX()[i] > 4){
            valuesEtaToPi0[0].push_back(PovlispC19400MeVGraph->GetX()[i]);
            valuesEtaToPi0[1].push_back(PovlispC19400MeVGraph->GetY()[i]);
            valuesEtaToPi0[2].push_back(PovlispC19400MeVGraph->GetEY()[i]);
            valuesEtaToPi0[3].push_back(PovlispC19400MeVGraph->GetEY()[i]);
            currentNumberOfPoint++;
        }
    }
    cout << "Number of pT bins: "<< currentNumberOfPoint <<  endl;


//     @article{Alverson:1993da,
//         author         = "Alverson, G. and others",
//         title          = "{Production of direct photons and neutral mesons at large
//         transverse momenta by $\pi^{-}$ and $p$ beams at
//         500-GeV/c}",
//         collaboration  = "FERMILAB-E706",
//         journal        = "Phys. Rev.",
//         volume         = "D48",
//         year           = "1993",
//         pages          = "5-28",
//         doi            = "10.1103/PhysRevD.48.5",
//         reportNumber   = "FERMILAB-PUB-93-007-E",
//         SLACcitation   = "%%CITATION = PHRVA,D48,5;%%"
//     }
//        sNN = 30.7 GeV
//          alverson93_pBe_plab500GeV_eta_pi0_ratio.txt
//          alverson93_piminusBe_plab500GeV_eta_pi0_ratio.txt
        ifstream Alverson30100MeV;
        Alverson30100MeV.open("OtherExperiments/alverson93_pBe_plab500GeV_eta_pi0_ratio.txt");
        cout << "alverson93_pBe_plab500GeV_eta_pi0_ratio.txt" << endl;
        lines = 0;

        while(!Alverson30100MeV.eof()){
            Alverson30100MeV >> pt[lines] >> value[lines] >> totErr[lines] ;
            xErr[lines] = 0.;
            lines++;
        }
        Alverson30100MeV.close();

        TGraph *Alverson30100MeVGraph = new TGraphErrors(lines-1,pt,value,xErr,totErr);
        Alverson30100MeVGraph->SetName("Alverson30100MeV");
        Alverson30100MeVGraph->SetTitle("p-Be (#sqrt{#it{s}_{_{NN}}}= 30.7 GeV)");
        Alverson30100MeVGraph->Print();
        for (Int_t i= 0; i< Alverson30100MeVGraph->GetN(); i++){
            if (Alverson30100MeVGraph->GetX()[i] > 4){
                valuesEtaToPi0[0].push_back(Alverson30100MeVGraph->GetX()[i]);
                valuesEtaToPi0[1].push_back(Alverson30100MeVGraph->GetY()[i]);
                valuesEtaToPi0[2].push_back(Alverson30100MeVGraph->GetEY()[i]);
                valuesEtaToPi0[3].push_back(Alverson30100MeVGraph->GetEY()[i]);
                currentNumberOfPoint++;
            }
        }
        cout << "Number of pT bins: "<< currentNumberOfPoint <<  endl;


//         @article{Apanasevich:2003zt,
//             author         = "Apanasevich, L. and others",
//             title          = "{Production of pi0 and eta mesons at large transverse
//             momenta in pi- p and pi- Be interactions at 515-GeV/c}",
//             collaboration  = "Fermilab E706",
//             journal        = "Phys. Rev.",
//             volume         = "D69",
//             year           = "2004",
//             pages          = "032003",
//             doi            = "10.1103/PhysRevD.69.032003",
//             eprint         = "hep-ex/0308022",
//             archivePrefix  = "arXiv",
//             primaryClass   = "hep-ex",
//             reportNumber   = "FERMILAB-PUB-02-230-E",
//             SLACcitation   = "%%CITATION = HEP-EX/0308022;%%"
//     }
//
//          apana03_pBe_plab530GeV_fnal706_eta_pi0_ratio.txt   // SNN 31.6 GeV
//          apana03_pBe_plab800GeV_fnal706_eta_pi0_ratio.txt   // SNN 38.8 GeV

    ifstream ApanapBe31600MeV;
    ApanapBe31600MeV.open("OtherExperiments/apana03_pBe_plab530GeV_fnal706_eta_pi0_ratio_mod.txt");
    cout << "apana03_pBe_plab530GeV_fnal706_eta_pi0_ratio_mod.txt" << endl;
    lines = 0;


    while(!ApanapBe31600MeV.eof()){
        ApanapBe31600MeV >> pt[lines] >> value[lines] >> totErr[lines];
        xErr[lines] = 0.;
        lines++;
    }
    ApanapBe31600MeV.close();
    TGraph *ApanapBe31600MeVGraph = new TGraphErrors(lines-1,pt,value,xErr,totErr);
    ApanapBe31600MeVGraph->SetName("ApanapBe31600MeV");
    ApanapBe31600MeVGraph->SetTitle("p-Be (#sqrt{#it{s_{_{NN}}}}= 31.6 GeV)");
    ApanapBe31600MeVGraph->Print();
    for (Int_t i= 0; i< ApanapBe31600MeVGraph->GetN(); i++){
        if (ApanapBe31600MeVGraph->GetX()[i] > 4){
            valuesEtaToPi0[0].push_back(ApanapBe31600MeVGraph->GetX()[i]);
            valuesEtaToPi0[1].push_back(ApanapBe31600MeVGraph->GetY()[i]);
            valuesEtaToPi0[2].push_back(ApanapBe31600MeVGraph->GetEY()[i]);
            valuesEtaToPi0[3].push_back(ApanapBe31600MeVGraph->GetEY()[i]);
            currentNumberOfPoint++;
        }
    }
    cout << "Number of pT bins: "<< currentNumberOfPoint <<  endl;

    ifstream ApanapBe38800MeV;
    ApanapBe38800MeV.open("OtherExperiments/apana03_pBe_plab800GeV_fnal706_eta_pi0_ratio_mod.txt");
    cout << "apana03_pBe_plab800GeV_fnal706_eta_pi0_ratio_mod.txt" << endl;
    lines = 0;


    while(!ApanapBe38800MeV.eof()){
        ApanapBe38800MeV >> pt[lines] >> value[lines] >> totErr[lines];
        xErr[lines] = 0.;
        lines++;
    }
    ApanapBe38800MeV.close();
    TGraph *ApanapBe38800MeVGraph = new TGraphErrors(lines-1,pt,value,xErr,totErr);
    ApanapBe38800MeVGraph->SetName("ApanapBe38800MeV");
    ApanapBe38800MeVGraph->SetTitle("p-Be (#sqrt{#it{s_{_{NN}}}}= 38.8 GeV)");
    ApanapBe38800MeVGraph->Print();
    for (Int_t i= 0; i< ApanapBe38800MeVGraph->GetN(); i++){
        if (ApanapBe38800MeVGraph->GetX()[i] > 4){
            valuesEtaToPi0[0].push_back(ApanapBe38800MeVGraph->GetX()[i]);
            valuesEtaToPi0[1].push_back(ApanapBe38800MeVGraph->GetY()[i]);
            valuesEtaToPi0[2].push_back(ApanapBe38800MeVGraph->GetEY()[i]);
            valuesEtaToPi0[3].push_back(ApanapBe38800MeVGraph->GetEY()[i]);
            currentNumberOfPoint++;
        }
    }
    cout << "Number of pT bins: "<< currentNumberOfPoint <<  endl;

// NO REFERNZ
//      tikhomirov95_pBe_plab450GeV_Helios_eta_pi0_ratio.txt

    TGraphAsymmErrors* Alice5TeVGraph        = NULL;
    TGraphAsymmErrors* Alice5TeVGraphStat    = NULL;
    TGraphAsymmErrors* Alice5TeVGraphSys     = NULL;
    if (fileNameAlice5TeV.CompareTo("") ){
        TFile* file5TeVAlice                 = new TFile(fileNameAlice5TeV.Data());
        Alice5TeVGraph                       = (TGraphAsymmErrors*)file5TeVAlice->Get("EtapPb_5.023TeV/graphRatioEtaToPi0CombpPb5023GeVTotErr");
        Alice5TeVGraphStat                   = (TGraphAsymmErrors*)file5TeVAlice->Get("EtapPb_5.023TeV/graphRatioEtaToPi0CombpPb5023GeVStatErr");
        Alice5TeVGraphSys                    = (TGraphAsymmErrors*)file5TeVAlice->Get("EtapPb_5.023TeV/graphRatioEtaToPi0CombpPb5023GeVSysErr");
        Alice5TeVGraph->SetName("Alice5TeV");
        Alice5TeVGraphStat->SetName("Alice5TeV_Stat");
        Alice5TeVGraphSys->SetName("Alice5TeV_Sys");
        for (Int_t i= 0; i< Alice5TeVGraph->GetN(); i++){
            if (Alice5TeVGraph->GetX()[i] > 4){
                valuesEtaToPi0[0].push_back(Alice5TeVGraph->GetX()[i]);
                valuesEtaToPi0[1].push_back(Alice5TeVGraph->GetY()[i]);
                valuesEtaToPi0[2].push_back(Alice5TeVGraph->GetEYlow()[i]);
                valuesEtaToPi0[3].push_back(Alice5TeVGraph->GetEYhigh()[i]);
                currentNumberOfPoint++;
            }
        }
        cout << "Number of pT bins: "<< currentNumberOfPoint <<  endl;
    }


    TGraphAsymmErrors* graphWorldHigh               = new TGraphAsymmErrors(currentNumberOfPoint);
    graphWorldHigh->SetName("AllWorldDataAbove4GeV");
    for (Int_t k = 0; k< currentNumberOfPoint; k++){
        graphWorldHigh->SetPoint(k, valuesEtaToPi0[0].at(k), valuesEtaToPi0[1].at(k));
        graphWorldHigh->SetPointError(k, 0.00, 0.00, valuesEtaToPi0[2].at(k), valuesEtaToPi0[3].at(k));
    }
    graphWorldHigh->Print();
    graphWorldHigh->Sort();
    cout << "after sorting" << endl;
    graphWorldHigh->Print();

    TF1* etaToPi0ConstData  = new TF1("etaToPi0ConstData","[0]",4,20);
    graphWorldHigh->Fit(etaToPi0ConstData,"QRME0","",4,20);

    cout  << "***********************************************************************************************************" << endl;
    cout  << "high pt eta/pi0 - data, tot: " << etaToPi0ConstData->GetParameter(0) << "+-"<< etaToPi0ConstData->GetParError(0) << endl;
    cout  << "***********************************************************************************************************" << endl;
    cout  << "***********************************************************************************************************" << endl;

    const char* OutputNameWorld ="WorldDataPi0EtapPb.root";
    TFile* WorldData = new TFile(OutputNameWorld,"RECREATE");
        AgakishievpAu29100MeVGraph->Write();
        AgakishievpBe29100MeVGraph->Write();
        PovlispBe19400MeVGraph->Write();
        PovlispAl19400MeVGraph->Write();
        PovlispC19400MeVGraph->Write();
        Alverson30100MeVGraph->Write();
        PhenixdAu200GeVGraph->Write();
        ApanapBe31600MeVGraph->Write();
        ApanapBe38800MeVGraph->Write();
        graphWorldHigh->Write();
        if (Alice5TeVGraph)Alice5TeVGraph->Write();
        if (Alice5TeVGraphStat)Alice5TeVGraphStat->Write();
        if (Alice5TeVGraphSys)Alice5TeVGraphSys->Write();
        etaToPi0ConstData->Write();
    WorldData->Write();
    WorldData->Close();



}