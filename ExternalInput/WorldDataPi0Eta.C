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


void WorldDataPi0Eta( TString nameFile8TeV      = "",
                      TString nameFile2760GeV   = "CombNeutralMesons/CombinedResultsPaperPP2760GeV_2017_07_10_Pub2017.root"
                    ){

// 	 	const int N_hh = 18*2;
// 	TString files_hh[N_hh] = {
// 		"donaldson78_pp_plab100GeV_eta_pi0_ratio.txt","p+p (#sqrt{#it{s}}= 13.8 GeV)",,
// 		"donaldson78_pp_plab200GeV_eta_pi0_ratio.txt","p+p (#sqrt{#it{s}}= 19.4 GeV)",
// 		"bonesi89_pp_plab280GeV_wa70_eta_pi0_ratio.txt","p+p (#sqrt{#it{s}}= 23 GeV)",
// 		"antille87_pp_24.3GeV_eta_pi0_ratio.txt","p+p (#sqrt{#it{s}}= 24.3 GeV)",
// 		"antille87_ppbar_24.3GeV_eta_pi0_ratio.txt","#bar{p}+p (#sqrt{#it{s}}= 24.3 GeV)",
// 		"aguilar91_pp_plab400GeV_na27_eta_pi0_ratio.txt","p+p (#sqrt{#it{s}}= 27.5 GeV)",
// 		"amaldi79_pp_30.6GeV_isr_eta_pi0_ratio.txt","p+p (#sqrt{#it{s}}= 30.6 GeV)",
// 		"kourkou79_pp_30.6GeV_eta_pi0_ratio.txt","p+p (#sqrt{#it{s}}= 30.6 GeV)",
// 		"apana02_pp_plab530GeV_fnal706_eta_pi0_ratio.txt","p+p (#sqrt{#it{s}}= 31.6 GeV)",
// 		"apana02_pp_plab800GeV_fnal706_eta_pi0_ratio.txt","p+p (#sqrt{#it{s}}= 38.8 GeV)",
// 		"kourkou79_pp_52.7GeV_eta_pi0_ratio.txt","p+p (#sqrt{#it{s}}= 52.7 GeV)",
// 		"akesson85_pbarp_53GeV_eta_pi0_ratio.txt","#bar{p}+p (#sqrt{#it{s}}= 53 GeV)",
// 		"akesson85_pp_53GeV_eta_pi0_ratio.txt","p+p (#sqrt{#it{s}}= 53 GeV)",
// 		"amaldi79_pp_53GeV_isr_eta_pi0_ratio.txt","p+p (#sqrt{#it{s}}= 53.2 GeV)",
// 		"kourkou79_pp_62.4GeV_eta_pi0_ratio.txt","p+p (#sqrt{#it{s}}= 62.4 GeV)",
// 		"akesson83_pp_63GeV_eta_pi0_ratio.txt","p+p (#sqrt{#it{s}}= 63 GeV)",
// 		"phenix_pp_200GeV_eta_pi0_ratio.txt","p+p (#sqrt{#it{s}}= 200 GeV)",
// 		"banner85_ppbar_540GeV_UA2_eta_pi0_ratio.txt","#bar{p}+p (#sqrt{#it{s}}= 540 GeV)",
// 	};

    Double_t pt[400];
    Double_t value[400];
    Double_t totErr[400];
    Double_t xErr[400];
    Double_t statErr[400];
    Double_t sysErr[400];
    Double_t sysErrB[400];
    Double_t a = 0;
    Double_t b = 0;
    Double_t c = 0;
    Double_t d = 0;

    // Create 1 graph containing all point above 4 GeV
    vector<Double_t> *valuesEtaToPi0            = new vector<Double_t>[4];     // iCent x iParticle x nMeasurements matrix with theory curves
    Int_t currentNumberOfPoint                  = 0;
    vector<Double_t> *valuesEtaToPi0ALICE       = new vector<Double_t>[4];     // iCent x iParticle x nMeasurements matrix with theory curves
    Int_t currentNumberOfPointALICE             = 0;

    //Donaldson100GeV
//         @article{Donaldson:1977jt,
//             author         = "Donaldson, G. J. and others",
//             title          = "{Inclusive eta Production at Large Transverse Momenta}",
//             journal        = "Phys. Rev. Lett.",
//             volume         = "40",
//             year           = "1978",
//             pages          = "684",
//             doi            = "10.1103/PhysRevLett.40.684",
//             reportNumber   = "BNL-23702",
//             SLACcitation   = "%%CITATION = PRLTA,40,684;%%"
//         }

    ifstream Donaldson100GeV;
    Donaldson100GeV.open("OtherExperiments/donaldson78_pp_plab100GeV_eta_pi0_ratio.txt");
    cout << "donaldson78_pp_plab100GeV_eta_pi0_ratio.txt" << endl;
    Int_t lines = 0;

    while(!Donaldson100GeV.eof()){
        Donaldson100GeV >> a>> b >> value[lines] >> totErr[lines] >> c >> d;
        pt[lines] = (a+b)/2;
        xErr[lines] = 0.;
        lines++;
    }
    Donaldson100GeV.close();

    TGraph *Donaldson100GeVGraph = new TGraphErrors(lines-1,pt,value,xErr,totErr);
    Donaldson100GeVGraph->SetName("donaldson100GeV");
    Donaldson100GeVGraph->SetTitle("p+p (#sqrt{#it{s}}= 13.8 GeV)");
    Donaldson100GeVGraph->Print();
    for (Int_t i= 0; i< Donaldson100GeVGraph->GetN(); i++){
        if (Donaldson100GeVGraph->GetX()[i] > 4){
            valuesEtaToPi0[0].push_back(Donaldson100GeVGraph->GetX()[i]);
            valuesEtaToPi0[1].push_back(Donaldson100GeVGraph->GetY()[i]);
            valuesEtaToPi0[2].push_back(Donaldson100GeVGraph->GetEY()[i]);
            valuesEtaToPi0[3].push_back(Donaldson100GeVGraph->GetEY()[i]);
            currentNumberOfPoint++;
        }
    }
    cout << "Number of pT bins: "<< currentNumberOfPoint <<  endl;

    //Donaldson200GeV
    ifstream Donaldson200GeV;
    Donaldson200GeV.open("OtherExperiments/donaldson78_pp_plab200GeV_eta_pi0_ratio.txt");
    cout << "donaldson78_pp_plab200GeV_eta_pi0_ratio.txt" << endl;
    lines = 0;

    while(!Donaldson200GeV.eof()){
        Donaldson200GeV >> a>> b >> value[lines] >> totErr[lines] >> c >> d;
        pt[lines] = (a+b)/2;
        xErr[lines] = 0.;
        lines++;
    }
    Donaldson200GeV.close();
    TGraph *Donaldson200GeVGraph = new TGraphErrors(lines-1,pt,value,xErr,totErr);
    Donaldson200GeVGraph->SetName("donaldson200GeV");
    Donaldson200GeVGraph->SetTitle("p+p (#sqrt{#it{s}}= 19.4 GeV)");
    Donaldson200GeVGraph->Print();
    for (Int_t i= 0; i< Donaldson200GeVGraph->GetN(); i++){
        if (Donaldson200GeVGraph->GetX()[i] > 4){
            valuesEtaToPi0[0].push_back(Donaldson200GeVGraph->GetX()[i]);
            valuesEtaToPi0[1].push_back(Donaldson200GeVGraph->GetY()[i]);
            valuesEtaToPi0[2].push_back(Donaldson200GeVGraph->GetEY()[i]);
            valuesEtaToPi0[3].push_back(Donaldson200GeVGraph->GetEY()[i]);
            currentNumberOfPoint++;
        }
    }
    cout << "Number of pT bins: "<< currentNumberOfPoint <<  endl;


    // 		"bonesi89_pp_plab280GeV_wa70_eta_pi0_ratio.txt","p+p (#sqrt{#it{s}}= 23 GeV)",
//     @article{Bonesini:1988rd,
//         author         = "Bonesini, M. and others",
//         title          = "{High Transverse Momentum $\eta$ Production in $\pi^- p$,
//         $\pi^+ p$ and $p p$ Interactions at 280-{GeV}/$c$}",
//         collaboration  = "WA70",
//         journal        = "Z. Phys.",
//         volume         = "C42",
//         year           = "1989",
//         pages          = "527",
//         doi            = "10.1007/BF01557657",
//         reportNumber   = "CERN-EP/88-178",
//         SLACcitation   = "%%CITATION = ZEPYA,C42,527;%%"
//     }
    ifstream Bonesi280GeV;
    Bonesi280GeV.open("OtherExperiments/bonesi89_pp_plab280GeV_wa70_eta_pi0_ratio.txt");
    cout << "bonesi89_pp_plab280GeV_wa70_eta_pi0_ratio.txt" << endl;
    lines = 0;

    while(!Bonesi280GeV.eof()){
        Bonesi280GeV >> pt[lines] >> value[lines] >> totErr[lines];
        xErr[lines] = 0.;
        lines++;
    }
    Bonesi280GeV.close();
    TGraph *Bonesi280GeVGraph = new TGraphErrors(lines-1,pt,value,xErr,totErr);
    Bonesi280GeVGraph->SetName("Bonesi280GeV");
    Bonesi280GeVGraph->SetTitle("p+p (#sqrt{#it{s}}= 23 GeV)");
    Bonesi280GeVGraph->Print();
    for (Int_t i= 0; i< Bonesi280GeVGraph->GetN(); i++){
        if (Bonesi280GeVGraph->GetX()[i] > 4){
            valuesEtaToPi0[0].push_back(Bonesi280GeVGraph->GetX()[i]);
            valuesEtaToPi0[1].push_back(Bonesi280GeVGraph->GetY()[i]);
            valuesEtaToPi0[2].push_back(Bonesi280GeVGraph->GetEY()[i]);
            valuesEtaToPi0[3].push_back(Bonesi280GeVGraph->GetEY()[i]);
            currentNumberOfPoint++;
        }
    }
    cout << "Number of pT bins: "<< currentNumberOfPoint <<  endl;


// 		"antille87_pp_24.3GeV_eta_pi0_ratio.txt","p+p (#sqrt{#it{s}}= 24.3 GeV)",
//     @article{Antille:1987kr,
//         author         = "Antille, J. and others",
//         title          = "{A Measurement of the Inclusive $\pi^0$ and $\eta$
//         Production Cross-sections at High $p_T$ in $p \bar{p}$ and
//         $p p$ Collisions at $\sqrt{s}=24$.3-{GeV}}",
//         collaboration  = "CERN-Lausanne-Michigan-Rockefeller, UA6",
//         journal        = "Phys. Lett.",
//         volume         = "B194",
//         year           = "1987",
//         pages          = "568",
//         doi            = "10.1016/0370-2693(87)90236-X",
//         reportNumber   = "CERN-EP-87-121",
//         SLACcitation   = "%%CITATION = PHLTA,B194,568;%%"
//     }

    ifstream Antille87pp;
    Antille87pp.open("OtherExperiments/antille87_pp_24.3GeV_eta_pi0_ratio.txt");
    cout << "antille87_pp_24.3GeV_eta_pi0_ratio.txt" << endl;
    lines = 0;

    while(!Antille87pp.eof()){
        Antille87pp >> pt[lines] >> value[lines] >> totErr[lines];
        xErr[lines] = 0.;
        lines++;
    }
    Antille87pp.close();
    TGraph *Antille87ppGraph = new TGraphErrors(lines-1,pt,value,xErr,totErr);
    Antille87ppGraph->SetName("Antille24.3GeVpp");
    Antille87ppGraph->SetTitle("p+p (#sqrt{#it{s}}= 24.3 GeV)");
    Antille87ppGraph->Print();
    for (Int_t i= 0; i< Antille87ppGraph->GetN(); i++){
        if (Antille87ppGraph->GetX()[i] > 4){
            valuesEtaToPi0[0].push_back(Antille87ppGraph->GetX()[i]);
            valuesEtaToPi0[1].push_back(Antille87ppGraph->GetY()[i]);
            valuesEtaToPi0[2].push_back(Antille87ppGraph->GetEY()[i]);
            valuesEtaToPi0[3].push_back(Antille87ppGraph->GetEY()[i]);
            currentNumberOfPoint++;
        }
    }
    cout << "Number of pT bins: "<< currentNumberOfPoint <<  endl;

// 		"antille87_ppbar_24.3GeV_eta_pi0_ratio.txt","#bar{p}+p (#sqrt{#it{s}}= 24.3 GeV)",
    ifstream Antille87ppbar;
    Antille87ppbar.open("OtherExperiments/antille87_ppbar_24.3GeV_eta_pi0_ratio.txt");
    cout << "antille87_ppbar_24.3GeV_eta_pi0_ratio.txt" << endl;
    lines = 0;

    while(!Antille87ppbar.eof()){
        Antille87ppbar >> pt[lines] >> value[lines] >> totErr[lines];
        xErr[lines] = 0.;
        lines++;
    }
    Antille87ppbar.close();
    TGraph *Antille87ppbarGraph = new TGraphErrors(lines-1,pt,value,xErr,totErr);
    Antille87ppbarGraph->SetName("Antille24.3GeVppbar");
    Antille87ppbarGraph->SetTitle("#bar{p}+p (#sqrt{#it{s}}= 24.3 GeV)");
    Antille87ppbarGraph->Print();
    for (Int_t i= 0; i< Antille87ppbarGraph->GetN(); i++){
        if (Antille87ppbarGraph->GetX()[i] > 4){
            valuesEtaToPi0[0].push_back(Antille87ppbarGraph->GetX()[i]);
            valuesEtaToPi0[1].push_back(Antille87ppbarGraph->GetY()[i]);
            valuesEtaToPi0[2].push_back(Antille87ppbarGraph->GetEY()[i]);
            valuesEtaToPi0[3].push_back(Antille87ppbarGraph->GetEY()[i]);
            currentNumberOfPoint++;
        }
    }
    cout << "Number of pT bins: "<< currentNumberOfPoint <<  endl;

    // 		"aguilar91_pp_plab400GeV_na27_eta_pi0_ratio.txt","p+p (#sqrt{#it{s}}= 27.5 GeV)",
//     @article{AguilarBenitez:1991yy,
//         author         = "Aguilar-Benitez, M. and others",
//         title          = "{Inclusive particle production in 400-GeV/c p p
//         interactions}",
//         journal        = "Z. Phys.",
//         volume         = "C50",
//         year           = "1991",
//         pages          = "405-426",
//         doi            = "10.1007/BF01551452",
//         reportNumber   = "CERN-PPE-91-21, EXP-CERN-NA27",
//         SLACcitation   = "%%CITATION = ZEPYA,C50,405;%%"
//     }
    ifstream Aguilar400GeV;
    Aguilar400GeV.open("OtherExperiments/aguilar91_pp_plab400GeV_na27_eta_pi0_ratio.txt");
    cout << "aguilar91_pp_plab400GeV_na27_eta_pi0_ratio.txt" << endl;
    lines = 0;

    while(!Aguilar400GeV.eof()){
        Aguilar400GeV >> pt[lines] >> value[lines] >> totErr[lines];
        xErr[lines] = 0.;
        lines++;
    }
    Aguilar400GeV.close();
    TGraph *Aguilar400GeVGraph = new TGraphErrors(lines-1,pt,value,xErr,totErr);
    Aguilar400GeVGraph->SetName("Aguilar400GeV");
    Aguilar400GeVGraph->SetTitle("p+p (#sqrt{#it{s}}= 27.5 GeV)");
    Aguilar400GeVGraph->Print();
    for (Int_t i= 0; i< Aguilar400GeVGraph->GetN(); i++){
        if (Aguilar400GeVGraph->GetX()[i] > 4){
            valuesEtaToPi0[0].push_back(Aguilar400GeVGraph->GetX()[i]);
            valuesEtaToPi0[1].push_back(Aguilar400GeVGraph->GetY()[i]);
            valuesEtaToPi0[2].push_back(Aguilar400GeVGraph->GetEY()[i]);
            valuesEtaToPi0[3].push_back(Aguilar400GeVGraph->GetEY()[i]);
            currentNumberOfPoint++;
        }
    }
    cout << "Number of pT bins: "<< currentNumberOfPoint <<  endl;

// 		"amaldi79_pp_30.6GeV_isr_eta_pi0_ratio.txt","p+p (#sqrt{#it{s}}= 30.6 GeV)",
//     @article{Amaldi:1979zk,
//         author         = "Amaldi, E. and others",
//         title          = "{INCLUSIVE ETA PRODUCTION IN P P COLLISION AT ISR
//         ENERGIES}",
//         journal        = "Nucl. Phys.",
//         volume         = "B158",
//         year           = "1979",
//         pages          = "1-10",
//         doi            = "10.1016/0550-3213(79)90183-4",
//         SLACcitation   = "%%CITATION = NUPHA,B158,1;%%"
//     }
    ifstream Amaldi79pp;
    Amaldi79pp.open("OtherExperiments/amaldi79_pp_30.6GeV_isr_eta_pi0_ratio.txt");
    cout << "amaldi79_pp_30.6GeV_isr_eta_pi0_ratio.txt" << endl;
    lines = 0;

    while(!Amaldi79pp.eof()){
        Amaldi79pp >> pt[lines] >> value[lines] >> totErr[lines];
        xErr[lines] = 0.;
        lines++;
    }
    Amaldi79pp.close();
    TGraph *Amaldi79ppGraph = new TGraphErrors(lines-1,pt,value,xErr,totErr);
    Amaldi79ppGraph->SetName("Amaldi30.6GeVpp");
    Amaldi79ppGraph->SetTitle("p+p (#sqrt{#it{s}}= 30.6 GeV)");
    Amaldi79ppGraph->Print();
    for (Int_t i= 0; i< Amaldi79ppGraph->GetN(); i++){
        if (Amaldi79ppGraph->GetX()[i] > 4){
            valuesEtaToPi0[0].push_back(Amaldi79ppGraph->GetX()[i]);
            valuesEtaToPi0[1].push_back(Amaldi79ppGraph->GetY()[i]);
            valuesEtaToPi0[2].push_back(Amaldi79ppGraph->GetEY()[i]);
            valuesEtaToPi0[3].push_back(Amaldi79ppGraph->GetEY()[i]);
            currentNumberOfPoint++;
        }
    }
    cout << "Number of pT bins: "<< currentNumberOfPoint <<  endl;


// 		"kourkou79_pp_30.6GeV_eta_pi0_ratio.txt","p+p (#sqrt{#it{s}}= 30.6 GeV)",

//     @article{Kourkoumelis:1979tr,
//         author         = "Kourkoumelis, C. and others",
//         title          = "{Inclusive $\pi^0$ Production at Very Large $p_T$ at the
//         {ISR}}",
//         journal        = "Phys. Lett.",
//         volume         = "84B",
//         year           = "1979",
//         pages          = "271-276",
//         doi            = "10.1016/0370-2693(79)90301-0",
//         reportNumber   = "CERN-EP-79-29, BNL-26126",
//         SLACcitation   = "%%CITATION = PHLTA,84B,271;%%"
//     }
//     @article{Kourkoumelis:1979ts,
//         author         = "Kourkoumelis, C. and others",
//         title          = "{Inclusive $\eta$ Production at High $p_T$ at the {ISR}}",
//         journal        = "Phys. Lett.",
//         volume         = "84B",
//         year           = "1979",
//         pages          = "277-280",
//         doi            = "10.1016/0370-2693(79)90302-2",
//         reportNumber   = "CERN-EP-79-30",
//         SLACcitation   = "%%CITATION = PHLTA,84B,277;%%"
//     }
    ifstream Kourkou79pp;
    Kourkou79pp.open("OtherExperiments/kourkou79_pp_30.6GeV_eta_pi0_ratio.txt");
    cout << "kourkou79_pp_30.6GeV_eta_pi0_ratio.txt" << endl;
    lines = 0;

    while(!Kourkou79pp.eof()){
        Kourkou79pp >> pt[lines] >> value[lines] >> totErr[lines];
        xErr[lines] = 0.;
        lines++;
    }
    Kourkou79pp.close();
    TGraph *Kourkou79ppGraph = new TGraphErrors(lines-1,pt,value,xErr,totErr);
    Kourkou79ppGraph->SetName("Kourkou30.6GeVpp");
    Kourkou79ppGraph->SetTitle("p+p (#sqrt{#it{s}}= 30.6 GeV)");
    Kourkou79ppGraph->Print();
    for (Int_t i= 0; i< Kourkou79ppGraph->GetN(); i++){
        if (Kourkou79ppGraph->GetX()[i] > 4){
            valuesEtaToPi0[0].push_back(Kourkou79ppGraph->GetX()[i]);
            valuesEtaToPi0[1].push_back(Kourkou79ppGraph->GetY()[i]);
            valuesEtaToPi0[2].push_back(Kourkou79ppGraph->GetEY()[i]);
            valuesEtaToPi0[3].push_back(Kourkou79ppGraph->GetEY()[i]);
            currentNumberOfPoint++;
        }
    }
    cout << "Number of pT bins: "<< currentNumberOfPoint <<  endl;


// 		"apana02_pp_plab530GeV_fnal706_eta_pi0_ratio.txt","p+p (#sqrt{#it{s}}= 31.6 GeV)",
    ifstream Apana530GeV;
    Apana530GeV.open("OtherExperiments/apana02_pp_plab530GeV_fnal706_eta_pi0_ratio.txt");
    cout << "apana02_pp_plab530GeV_fnal706_eta_pi0_ratio.txt" << endl;
    lines = 0;
//     @article{Apanasevich:2002wt,
//         author         = "Apanasevich, L. and others",
//         title          = "{Production of pi0 and eta mesons at large transverse
//         momenta in pp and p Be interactions at 530-GeV/c and
//         800-GeV/c}",
//         collaboration  = "Fermilab E706",
//         journal        = "Phys. Rev.",
//         volume         = "D68",
//         year           = "2003",
//         pages          = "052001",
//         doi            = "10.1103/PhysRevD.68.052001",
//         eprint         = "hep-ex/0204031",
//         archivePrefix  = "arXiv",
//         primaryClass   = "hep-ex",
//         reportNumber   = "FERMILAB-PUB-01-021-E",
//         SLACcitation   = "%%CITATION = HEP-EX/0204031;%%"
//     }
    while(!Apana530GeV.eof()){
        Apana530GeV >> pt[lines] >> value[lines] >> totErr[lines];
        xErr[lines] = 0.;
        lines++;
    }
    Apana530GeV.close();
    TGraph *Apana530GeVGraph = new TGraphErrors(lines-1,pt,value,xErr,totErr);
    Apana530GeVGraph->SetName("Apana530GeV");
    Apana530GeVGraph->SetTitle("p+p (#sqrt{#it{s}}= 31.6 GeV)");
    Apana530GeVGraph->Print();
    for (Int_t i= 0; i< Apana530GeVGraph->GetN(); i++){
        if (Apana530GeVGraph->GetX()[i] > 4){
            valuesEtaToPi0[0].push_back(Apana530GeVGraph->GetX()[i]);
            valuesEtaToPi0[1].push_back(Apana530GeVGraph->GetY()[i]);
            valuesEtaToPi0[2].push_back(Apana530GeVGraph->GetEY()[i]);
            valuesEtaToPi0[3].push_back(Apana530GeVGraph->GetEY()[i]);
            currentNumberOfPoint++;
        }
    }
    cout << "Number of pT bins: "<< currentNumberOfPoint <<  endl;

    // 		"apana02_pp_plab800GeV_fnal706_eta_pi0_ratio.txt","p+p (#sqrt{#it{s}}= 38.8 GeV)",

    ifstream Apana800GeV;
    Apana800GeV.open("OtherExperiments/apana02_pp_plab800GeV_fnal706_eta_pi0_ratio.txt");
    cout << "apana02_pp_plab800GeV_fnal706_eta_pi0_ratio.txt" << endl;
    lines = 0;

    while(!Apana800GeV.eof()){
        Apana800GeV >> pt[lines] >> value[lines] >> totErr[lines];
        xErr[lines] = 0.;
        lines++;
    }
    Apana800GeV.close();
    TGraph *Apana800GeVGraph = new TGraphErrors(lines-1,pt,value,xErr,totErr);
    Apana800GeVGraph->SetName("Apana800GeV");
    Apana800GeVGraph->SetTitle("p+p (#sqrt{#it{s}}= 38.8 GeV)");
    Apana800GeVGraph->Print();
    for (Int_t i= 0; i< Apana800GeVGraph->GetN(); i++){
        if (Apana800GeVGraph->GetX()[i] > 4){
            valuesEtaToPi0[0].push_back(Apana800GeVGraph->GetX()[i]);
            valuesEtaToPi0[1].push_back(Apana800GeVGraph->GetY()[i]);
            valuesEtaToPi0[2].push_back(Apana800GeVGraph->GetEY()[i]);
            valuesEtaToPi0[3].push_back(Apana800GeVGraph->GetEY()[i]);
            currentNumberOfPoint++;
        }
    }
    cout << "Number of pT bins: "<< currentNumberOfPoint <<  endl;

// 		"kourkou79_pp_52.7GeV_eta_pi0_ratio.txt","p+p (#sqrt{#it{s}}= 52.7 GeV)",
    ifstream Kourkou79pp52;
    Kourkou79pp52.open("OtherExperiments/kourkou79_pp_52.7GeV_eta_pi0_ratio.txt");
    cout << "kourkou79_pp_52.7GeV_eta_pi0_ratio.txt" << endl;
    lines = 0;

    while(!Kourkou79pp52.eof()){
        Kourkou79pp52 >> pt[lines] >> value[lines] >> totErr[lines];
        xErr[lines] = 0.;
        lines++;
    }
    Kourkou79pp52.close();
    TGraph *Kourkou79pp52Graph = new TGraphErrors(lines-1,pt,value,xErr,totErr);
    Kourkou79pp52Graph->SetName("Kourkou52.7GeVpp");
    Kourkou79pp52Graph->SetTitle("p+p (#sqrt{#it{s}}= 52.7 GeV)");
    Kourkou79pp52Graph->Print();
    for (Int_t i= 0; i< Kourkou79pp52Graph->GetN(); i++){
        if (Kourkou79pp52Graph->GetX()[i] > 4){
            valuesEtaToPi0[0].push_back(Kourkou79pp52Graph->GetX()[i]);
            valuesEtaToPi0[1].push_back(Kourkou79pp52Graph->GetY()[i]);
            valuesEtaToPi0[2].push_back(Kourkou79pp52Graph->GetEY()[i]);
            valuesEtaToPi0[3].push_back(Kourkou79pp52Graph->GetEY()[i]);
            currentNumberOfPoint++;
        }
    }
    cout << "Number of pT bins: "<< currentNumberOfPoint <<  endl;


// 		"akesson85_pbarp_53GeV_eta_pi0_ratio.txt","#bar{p}+p (#sqrt{#it{s}}= 53 GeV)",
//     @article{Akesson:1985za,
//         author         = "Akesson, T. and others",
//         title          = "{A comparison of direct photon, pi0, and eta production
//         in p anti-p and pp interactions at the CERN ISR}",
//         collaboration  = "Axial Field Spectrometer",
//         journal        = "Phys. Lett.",
//         volume         = "158B",
//         year           = "1985",
//         pages          = "282-288",
//         doi            = "10.1016/0370-2693(85)90971-2",
//         reportNumber   = "CERN-EP-85-55",
//         SLACcitation   = "%%CITATION = PHLTA,158B,282;%%"
//     }
    ifstream Akesson53GeVppbar;
    Akesson53GeVppbar.open("OtherExperiments/akesson85_pbarp_53GeV_eta_pi0_ratio.txt");
    cout << "akesson85_pbarp_53GeV_eta_pi0_ratio.txt" << endl;
    lines = 0;

    while(!Akesson53GeVppbar.eof()){
        Akesson53GeVppbar >> pt[lines] >> value[lines] >> totErr[lines];
        xErr[lines] = 0.;
        lines++;
    }
    Akesson53GeVppbar.close();
    TGraph *Akesson53GeVppbarGraph = new TGraphErrors(lines-1,pt,value,xErr,totErr);
    Akesson53GeVppbarGraph->SetName("Akesson53GeVppbar");
    Akesson53GeVppbarGraph->SetTitle("#bar{p}+p (#sqrt{#it{s}}= 53 GeV)");
    Akesson53GeVppbarGraph->Print();
    for (Int_t i= 0; i< Akesson53GeVppbarGraph->GetN(); i++){
        if (Akesson53GeVppbarGraph->GetX()[i] > 4){
            valuesEtaToPi0[0].push_back(Akesson53GeVppbarGraph->GetX()[i]);
            valuesEtaToPi0[1].push_back(Akesson53GeVppbarGraph->GetY()[i]);
            valuesEtaToPi0[2].push_back(Akesson53GeVppbarGraph->GetEY()[i]);
            valuesEtaToPi0[3].push_back(Akesson53GeVppbarGraph->GetEY()[i]);
            currentNumberOfPoint++;
        }
    }
    cout << "Number of pT bins: "<< currentNumberOfPoint <<  endl;


    // 		"akesson85_pp_53GeV_eta_pi0_ratio.txt","p+p (#sqrt{#it{s}}= 53 GeV)",
//     @article{Akesson:1986nx,
//         author         = "Akesson, T. and others",
//         title          = "{Inclusive $\eta$ Production at Low Transverse Momentum
//         in 63-{GeV} $p p$ Collisions at the {CERN} Intersecting
//         Storage Rings}",
//         booktitle      = "{Proceedings, 23RD International Conference on High
//         Energy Physics, JULY 16-23, 1986, Berkeley, CA}",
//         collaboration  = "Axial Field Spectrometer",
//         journal        = "Phys. Lett.",
//         volume         = "B178",
//         year           = "1986",
//         pages          = "447",
//         doi            = "10.1016/0370-2693(86)91409-7",
//         reportNumber   = "CERN-EP-86-84",
//         SLACcitation   = "%%CITATION = PHLTA,B178,447;%%"
//     }
    ifstream Akesson53GeVpp;
    Akesson53GeVpp.open("OtherExperiments/akesson85_pp_53GeV_eta_pi0_ratio.txt");
    cout << "akesson85_pp_53GeV_eta_pi0_ratio.txt" << endl;
    lines = 0;

    while(!Akesson53GeVpp.eof()){
        Akesson53GeVpp >> pt[lines] >> value[lines] >> totErr[lines];
        xErr[lines] = 0.;
        lines++;
    }
    Akesson53GeVpp.close();
    TGraph *Akesson53GeVppGraph = new TGraphErrors(lines-1,pt,value,xErr,totErr);
    Akesson53GeVppGraph->SetName("Akesson53GeVpp");
    Akesson53GeVppGraph->SetTitle("p+p (#sqrt{#it{s}}= 53 GeV)");
    Akesson53GeVppGraph->Print();
    for (Int_t i= 0; i< Akesson53GeVppGraph->GetN(); i++){
        if (Akesson53GeVppGraph->GetX()[i] > 4){
            valuesEtaToPi0[0].push_back(Akesson53GeVppGraph->GetX()[i]);
            valuesEtaToPi0[1].push_back(Akesson53GeVppGraph->GetY()[i]);
            valuesEtaToPi0[2].push_back(Akesson53GeVppGraph->GetEY()[i]);
            valuesEtaToPi0[3].push_back(Akesson53GeVppGraph->GetEY()[i]);
            currentNumberOfPoint++;
        }
    }
    cout << "Number of pT bins: "<< currentNumberOfPoint <<  endl;

    // 		"amaldi79_pp_53GeV_isr_eta_pi0_ratio.txt","p+p (#sqrt{#it{s}}= 53.2 GeV)",
    ifstream Amaldi53GeVpp;
    Amaldi53GeVpp.open("OtherExperiments/amaldi79_pp_53GeV_isr_eta_pi0_ratio.txt");
    cout << "amaldi79_pp_53GeV_isr_eta_pi0_ratio.txt" << endl;
    lines = 0;

    while(!Amaldi53GeVpp.eof()){
        Amaldi53GeVpp >> pt[lines] >> value[lines] >> totErr[lines];
        xErr[lines] = 0.;
        lines++;
    }
    Amaldi53GeVpp.close();
    TGraph *Amaldi53GeVppGraph = new TGraphErrors(lines-1,pt,value,xErr,totErr);
    Amaldi53GeVppGraph->SetName("Amaldi53GeVpp");
    Amaldi53GeVppGraph->SetTitle("p+p (#sqrt{#it{s}}= 53.2 GeV)");
    Amaldi53GeVppGraph->Print();
    cout << "Number of pT bins: "<< currentNumberOfPoint <<  endl;
    for (Int_t i= 0; i< Amaldi53GeVppGraph->GetN(); i++){
        if (Amaldi53GeVppGraph->GetX()[i] > 4){
            valuesEtaToPi0[0].push_back(Amaldi53GeVppGraph->GetX()[i]);
            valuesEtaToPi0[1].push_back(Amaldi53GeVppGraph->GetY()[i]);
            valuesEtaToPi0[2].push_back(Amaldi53GeVppGraph->GetEY()[i]);
            valuesEtaToPi0[3].push_back(Amaldi53GeVppGraph->GetEY()[i]);
            currentNumberOfPoint++;
        }
    }
    cout << "Number of pT bins: "<< currentNumberOfPoint <<  endl;


    // 		"kourkou79_pp_62.4GeV_eta_pi0_ratio.txt","p+p (#sqrt{#it{s}}= 62.4 GeV)",
    ifstream Kourkou79pp62;
    Kourkou79pp62.open("OtherExperiments/kourkou79_pp_62.4GeV_eta_pi0_ratio.txt");
    cout << "kourkou79_pp_62.4GeV_eta_pi0_ratio.txt" << endl;
    lines = 0;

    while(!Kourkou79pp62.eof()){
        Kourkou79pp62 >> pt[lines] >> value[lines] >> totErr[lines];
        xErr[lines] = 0.;
        lines++;
    }
    Kourkou79pp62.close();
    TGraph *Kourkou79pp62Graph = new TGraphErrors(lines-1,pt,value,xErr,totErr);
    Kourkou79pp62Graph->SetName("Kourkou79pp62");
    Kourkou79pp62Graph->SetTitle("p+p (#sqrt{#it{s}}= 62.4 GeV)");
    Kourkou79pp62Graph->Print();
    for (Int_t i= 0; i< Kourkou79pp62Graph->GetN(); i++){
        if (Kourkou79pp62Graph->GetX()[i] > 4){
            valuesEtaToPi0[0].push_back(Kourkou79pp62Graph->GetX()[i]);
            valuesEtaToPi0[1].push_back(Kourkou79pp62Graph->GetY()[i]);
            valuesEtaToPi0[2].push_back(Kourkou79pp62Graph->GetEY()[i]);
            valuesEtaToPi0[3].push_back(Kourkou79pp62Graph->GetEY()[i]);
            currentNumberOfPoint++;
        }
    }
    cout << "Number of pT bins: "<< currentNumberOfPoint <<  endl;

    // 		"akesson83_pp_63GeV_eta_pi0_ratio.txt","p+p (#sqrt{#it{s}}= 63 GeV)",
    ifstream Akesson63GeVpp;
    Akesson63GeVpp.open("OtherExperiments/akesson83_pp_63GeV_eta_pi0_ratio.txt");
    cout << "akesson83_pp_63GeV_eta_pi0_ratio.txt" << endl;
    lines = 0;
//     @article{Akesson:1983xq,
//         author         = "Akesson, T. and others",
//         title          = "{Production of $\pi^0$ and Eta0 at 11-degrees in $p p$
//         Collisions at $\sqrt{s}=63$-{GeV}}",
//         collaboration  = "Axial Field Spectrometer",
//         journal        = "Z. Phys.",
//         volume         = "C18",
//         year           = "1983",
//         pages          = "5",
//         doi            = "10.1007/BF01571698",
//         reportNumber   = "CERN-EP-83-17",
//         SLACcitation   = "%%CITATION = ZEPYA,C18,5;%%"
//     }
    while(!Akesson63GeVpp.eof()){
        Akesson63GeVpp >> a >> b >> value[lines] >> totErr[lines];
        pt[lines]   = (a+b)/2;
        xErr[lines] = 0.;
        lines++;
    }
    Akesson63GeVpp.close();
    TGraph *Akesson63GeVppGraph = new TGraphErrors(lines-1,pt,value,xErr,totErr);
    Akesson63GeVppGraph->SetName("Akesson63GeVpp");
    Akesson63GeVppGraph->SetTitle("p+p (#sqrt{#it{s}}= 63 GeV)");
    Akesson63GeVppGraph->Print();
    for (Int_t i= 0; i< Akesson63GeVppGraph->GetN(); i++){
        if (Akesson63GeVppGraph->GetX()[i] > 4){
            valuesEtaToPi0[0].push_back(Akesson63GeVppGraph->GetX()[i]);
            valuesEtaToPi0[1].push_back(Akesson63GeVppGraph->GetY()[i]);
            valuesEtaToPi0[2].push_back(Akesson63GeVppGraph->GetEY()[i]);
            valuesEtaToPi0[3].push_back(Akesson63GeVppGraph->GetEY()[i]);
            currentNumberOfPoint++;
        }
    }
    cout << "Number of pT bins: "<< currentNumberOfPoint <<  endl;



    // 		"phenix_pp_200GeV_eta_pi0_ratio.txt","p+p (#sqrt{#it{s}}= 200 GeV)",
//     @Article{Adler:2006bv,
//         Title                    = {{High transverse momentum $\eta$ meson production in $p^+ p$, $d^+$ Au and Au+Au collisions at $\sqrt{s_{\mbox{\tiny NN}}}$ = 200 GeV}},
//         Author                   = {Adler, S.S. and others},
//         Journal                  = {Phys.Rev.},
//         Year                     = {2007},
//         Pages                    = {024909},
//         Volume                   = {C75},
//
//         Archiveprefix            = {arXiv},
//         Collaboration            = {PHENIX Collaboration},
//         Doi                      = {10.1103/PhysRevC.75.024909},
//         Eprint                   = {nucl-ex/0611006},
//         Primaryclass             = {nucl-ex},
//         Slaccitation             = {%%CITATION = NUCL-EX/0611006;%%}
//     }
//
    ifstream Phenix200GeV;
    Phenix200GeV.open("OtherExperiments/phenix_pp_200GeV_eta_pi0_ratio.txt");
    cout << "phenix_pp_200GeV_eta_pi0_ratio.txt" << endl;
    lines = 0;


    while(!Phenix200GeV.eof()){
        Phenix200GeV >> pt[lines] >> value[lines] >> statErr[lines] >> sysErr[lines] >> sysErrB[lines];
        totErr[lines] = TMath::Sqrt(statErr[lines]*statErr[lines]+ sysErr[lines]*sysErr[lines] + sysErrB[lines]*sysErrB[lines]);
        xErr[lines] = 0.;
        lines++;
    }
    Phenix200GeV.close();
    TGraph *Phenix200GeVGraph = new TGraphErrors(lines-1,pt,value,xErr,totErr);
    Phenix200GeVGraph->SetName("Phenix200GeV");
    Phenix200GeVGraph->SetTitle("p+p (#sqrt{#it{s}}= 200 GeV)");
    Phenix200GeVGraph->Print();
    for (Int_t i= 0; i< Phenix200GeVGraph->GetN(); i++){
        if (Phenix200GeVGraph->GetX()[i] > 4){
            valuesEtaToPi0[0].push_back(Phenix200GeVGraph->GetX()[i]);
            valuesEtaToPi0[1].push_back(Phenix200GeVGraph->GetY()[i]);
            valuesEtaToPi0[2].push_back(Phenix200GeVGraph->GetEY()[i]);
            valuesEtaToPi0[3].push_back(Phenix200GeVGraph->GetEY()[i]);
            currentNumberOfPoint++;
        }
    }

    // 		"banner85_ppbar_540GeV_UA2_eta_pi0_ratio.txt","#bar{p}+p (#sqrt{#it{s}}= 540 GeV)",
//     @article{Banner:1984wh,
//         author         = "Banner, M. and others",
//         title          = "{Inclusive Particle Production in the Transverse Momentum
//         Range Between 0.25-{GeV}/$c$ and 40-{GeV}/$c$ at the
//         {CERN} S $p \bar{p}$ S Collider}",
//         collaboration  = "UA2, Bern-CERN-Copenhagen-Orsay-Pavia-Saclay",
//         journal        = "Z. Phys.",
//         volume         = "C27",
//         year           = "1985",
//         pages          = "329",
//         doi            = "10.1007/BF01548636",
//         reportNumber   = "SACLAY-DPHPE-84-07",
//         SLACcitation   = "%%CITATION = ZEPYA,C27,329;%%"
//     }
    ifstream Banner540GeV;
    Banner540GeV.open("OtherExperiments/banner85_ppbar_540GeV_UA2_eta_pi0_ratio.txt");
    cout << "banner85_ppbar_540GeV_UA2_eta_pi0_ratio.txt" << endl;
    lines = 0;

    while(!Banner540GeV.eof()){
        Banner540GeV >> pt[lines] >> value[lines] >> totErr[lines];
        xErr[lines] = 0.;
        lines++;
    }
    Banner540GeV.close();
    TGraph *Banner540GeVGraph = new TGraphErrors(lines-1,pt,value,xErr,totErr);
    Banner540GeVGraph->SetName("Banner540GeV");
    Banner540GeVGraph->SetTitle("#bar{p}+p (#sqrt{#it{s}}= 540 GeV)");
    Banner540GeVGraph->Print();
    for (Int_t i= 0; i< Banner540GeVGraph->GetN(); i++){
        if (Banner540GeVGraph->GetX()[i] > 4){
            valuesEtaToPi0[0].push_back(Banner540GeVGraph->GetX()[i]);
            valuesEtaToPi0[1].push_back(Banner540GeVGraph->GetY()[i]);
            valuesEtaToPi0[2].push_back(Banner540GeVGraph->GetEY()[i]);
            valuesEtaToPi0[3].push_back(Banner540GeVGraph->GetEY()[i]);
            currentNumberOfPoint++;
        }
    }
    cout << "Number of pT bins: "<< currentNumberOfPoint <<  endl;


    // eta/pi0 7TeV ALICE
    // Experiment: CERN-LHC-ALICE (ALICE)
    // Published in PL B717,162 (DOI:10.1016/j.physletb.2012.09.015)
    // Preprinted as CERN-PH-EP-2012-001
    // Archived as: ARXIV:1205.5724
//     @Article{Abelev:2012cn,
//         Title                    = {{Neutral pion and $\eta$ meson production in proton-proton collisions at $\sqrt{s}=0.9$ TeV and $\sqrt{s}=7$ TeV}},
//         Author                   = {Abelev, B. and others},
//         Journal                  = {Phys.Lett.},
//         Year                     = {2012},
//         Pages                    = {162-172},
//         Volume                   = {B717},
//
//         Archiveprefix            = {arXiv},
//         Collaboration            = {ALICE Collaboration},
//         Doi                      = {10.1016/j.physletb.2012.09.015},
//         Eprint                   = {1205.5724},
//         Primaryclass             = {hep-ex},
//         Reportnumber             = {CERN-PH-EP-2012-001},
//         Slaccitation             = {%%CITATION = ARXIV:1205.5724;%%}
//     }
    cout << "7TeV measurements" << endl;
    double Alice7TeV_xval[] = { 0.55, 0.85, 1.2, 1.6, 2.0, 2.4, 2.8, 3.25, 3.75,
        5.0, 7.0, 9.0, 12.5 };
    double Alice7TeV_xerrminus[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0 };
    double Alice7TeV_xerrplus[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0 };
    double Alice7TeV_yval[] = { 0.1067, 0.1655, 0.2501, 0.2866, 0.3514, 0.3547, 0.3757, 0.3987, 0.4655,
        0.4686, 0.5014, 0.7082, 0.5205 };
    double Alice7TeV_yerrminus[] = { 0.03347013594235912, 0.02529268669002959, 0.02603113520382851, 0.026801865606707307, 0.030696090956341654, 0.02972574641619618, 0.035979021665409415, 0.03927250946909301, 0.04585193561890272,
        0.04081041533726409, 0.08140718641495971, 0.18510386273657284, 0.27690873947927325 };
    double Alice7TeV_yerrplus[] = { 0.03347013594235912, 0.02529268669002959, 0.02603113520382851, 0.026801865606707307, 0.030696090956341654, 0.02972574641619618, 0.035979021665409415, 0.03927250946909301, 0.04585193561890272,
        0.04081041533726409, 0.08140718641495971, 0.18510386273657284, 0.27690873947927325 };
    double Alice7TeV_ystatminus[] = { 0.0259, 0.0234, 0.0179, 0.0175, 0.0201, 0.0219, 0.0232, 0.0267, 0.0336,
        0.0282, 0.0588, 0.1538, 0.1679 };
    double Alice7TeV_ystatplus[] = { 0.0259, 0.0234, 0.0179, 0.0175, 0.0201, 0.0219, 0.0232, 0.0267, 0.0336,
        0.0282, 0.0588, 0.1538, 0.1679 };
    int Alice7TeV_numpoints = 13;
    TGraphAsymmErrors* Alice7TeVGraph = new TGraphAsymmErrors(Alice7TeV_numpoints, Alice7TeV_xval, Alice7TeV_yval, Alice7TeV_xerrminus, Alice7TeV_xerrplus, Alice7TeV_yerrminus, Alice7TeV_yerrplus);
    Alice7TeVGraph->SetName("Alice7TeV");
    TFile* file7TeVAlice                                    = new TFile("CombNeutralMesons/CombinedResultsPP_ShiftedX_PaperRAA_16_May_2014_including7TeVand900GeVpublished.root");
    TGraphAsymmErrors* Alice7TeVGraphStat                   = (TGraphAsymmErrors*)file7TeVAlice->Get("graphEtaToPi0Comb7TeVStat");
    TGraphAsymmErrors* Alice7TeVGraphSys                    = (TGraphAsymmErrors*)file7TeVAlice->Get("graphEtaToPi0Comb7TeVSys");
    Alice7TeVGraphStat->SetName("Alice7TeV_Stat");
    Alice7TeVGraphSys->SetName("Alice7TeV_Sys");
    for (Int_t i= 0; i< Alice7TeVGraph->GetN(); i++){
        if (Alice7TeVGraph->GetX()[i] > 4){
            valuesEtaToPi0[0].push_back(Alice7TeVGraph->GetX()[i]);
            valuesEtaToPi0[1].push_back(Alice7TeVGraph->GetY()[i]);
            valuesEtaToPi0[2].push_back(Alice7TeVGraph->GetEYlow()[i]);
            valuesEtaToPi0[3].push_back(Alice7TeVGraph->GetEYhigh()[i]);
            currentNumberOfPoint++;
            valuesEtaToPi0ALICE[0].push_back(Alice7TeVGraph->GetX()[i]);
            valuesEtaToPi0ALICE[1].push_back(Alice7TeVGraph->GetY()[i]);
            valuesEtaToPi0ALICE[2].push_back(Alice7TeVGraph->GetEYlow()[i]);
            valuesEtaToPi0ALICE[3].push_back(Alice7TeVGraph->GetEYhigh()[i]);
            currentNumberOfPointALICE++;
            cout << "adding points" << endl;
        }
    }
    cout << "Number of pT bins: "<< currentNumberOfPoint <<  endl;
    cout << "Number of pT bins ALICE: "<< currentNumberOfPointALICE <<  endl;

    // eta/pi0 2.76TeV ALICE
    // Experiment: CERN-LHC-ALICE (ALICE)
    // @Article{Acharya:2017hyu,
    //         Title                    = {{Production of ${\pi ^0}$ and $\eta $ mesons up to high
    //         transverse momentum in pp collisions at 2.76 TeV}},
    //         Author                   = {Acharya, Shreyasi and others},
    //         Journal                  = {Eur. Phys. J.},
    //         Year                     = {2017},
    //         Number                   = {5},
    //         Pages                    = {339},
    //         Volume                   = {C77},
    //
    //         Archiveprefix            = {arXiv},
    //         Collaboration            = {ALICE},
    //         Doi                      = {10.1140/epjc/s10052-017-4890-x},
    //         Eprint                   = {1702.00917},
    //         Primaryclass             = {hep-ex},
    //         Reportnumber             = {CERN-EP-2017-019},
    //         Slaccitation             = {%%CITATION = ARXIV:1702.00917;%%}
    //     }
    TGraphAsymmErrors* Alice2760GeVGraph        = NULL;
    TGraphAsymmErrors* Alice2760GeVGraphStat    = NULL;
    TGraphAsymmErrors* Alice2760GeVGraphSys     = NULL;
    if (nameFile2760GeV.CompareTo("") ){
        cout << "2.76TeV measurements" << endl;
        TFile* file2760GeVAlice                 = new TFile(nameFile2760GeV.Data());
        Alice2760GeVGraph                       = (TGraphAsymmErrors*)file2760GeVAlice->Get("Eta2.76TeV/graphRatioEtaToPi0Comb2760GeVTotErr");
        Alice2760GeVGraphStat                   = (TGraphAsymmErrors*)file2760GeVAlice->Get("Eta2.76TeV/graphRatioEtaToPi0Comb2760GeVStatErr");
        Alice2760GeVGraphSys                    = (TGraphAsymmErrors*)file2760GeVAlice->Get("Eta2.76TeV/graphRatioEtaToPi0Comb2760GeVSysErr");
        Alice2760GeVGraph->SetName("Alice2760GeV");
        Alice2760GeVGraphStat->SetName("Alice2760GeV_Stat");
        Alice2760GeVGraphSys->SetName("Alice2760GeV_Sys");
        for (Int_t i= 0; i< Alice2760GeVGraph->GetN(); i++){
            if (Alice2760GeVGraph->GetX()[i] > 4){
                valuesEtaToPi0[0].push_back(Alice2760GeVGraph->GetX()[i]);
                valuesEtaToPi0[1].push_back(Alice2760GeVGraph->GetY()[i]);
                valuesEtaToPi0[2].push_back(Alice2760GeVGraph->GetEYlow()[i]);
                valuesEtaToPi0[3].push_back(Alice2760GeVGraph->GetEYhigh()[i]);
                currentNumberOfPoint++;
                valuesEtaToPi0ALICE[0].push_back(Alice2760GeVGraph->GetX()[i]);
                valuesEtaToPi0ALICE[1].push_back(Alice2760GeVGraph->GetY()[i]);
                valuesEtaToPi0ALICE[2].push_back(Alice2760GeVGraph->GetEYlow()[i]);
                valuesEtaToPi0ALICE[3].push_back(Alice2760GeVGraph->GetEYhigh()[i]);
                currentNumberOfPointALICE++;

            }
        }
        cout << "Number of pT bins: "<< currentNumberOfPoint <<  endl;
        cout << "Number of pT bins ALICE: "<< currentNumberOfPointALICE <<  endl;
    }


    TGraphAsymmErrors* Alice8TeVGraph        = NULL;
    TGraphAsymmErrors* Alice8TeVGraphStat    = NULL;
    TGraphAsymmErrors* Alice8TeVGraphSys     = NULL;
    if (nameFile8TeV.CompareTo("") ){
        cout << "8TeV measurements" << endl;
        TFile* file8TeVAlice                 = new TFile(nameFile8TeV.Data());
        Alice8TeVGraph                       = (TGraphAsymmErrors*)file8TeVAlice->Get("Eta8TeV/graphRatioEtaToPi0Comb8TeVTotErr");
        Alice8TeVGraphStat                   = (TGraphAsymmErrors*)file8TeVAlice->Get("Eta8TeV/graphRatioEtaToPi0Comb8TeVStatErr");
        Alice8TeVGraphSys                    = (TGraphAsymmErrors*)file8TeVAlice->Get("Eta8TeV/graphRatioEtaToPi0Comb8TeVSysErr");
        Alice8TeVGraph->SetName("Alice8TeV");
        Alice8TeVGraphStat->SetName("Alice8TeV_Stat");
        Alice8TeVGraphSys->SetName("Alice8TeV_Sys");
        for (Int_t i= 0; i< Alice8TeVGraph->GetN(); i++){
            if (Alice8TeVGraph->GetX()[i] > 4){
                valuesEtaToPi0[0].push_back(Alice8TeVGraph->GetX()[i]);
                valuesEtaToPi0[1].push_back(Alice8TeVGraph->GetY()[i]);
                valuesEtaToPi0[2].push_back(Alice8TeVGraph->GetEYlow()[i]);
                valuesEtaToPi0[3].push_back(Alice8TeVGraph->GetEYhigh()[i]);
                currentNumberOfPoint++;
                valuesEtaToPi0ALICE[0].push_back(Alice8TeVGraph->GetX()[i]);
                valuesEtaToPi0ALICE[1].push_back(Alice8TeVGraph->GetY()[i]);
                valuesEtaToPi0ALICE[2].push_back(Alice8TeVGraph->GetEYlow()[i]);
                valuesEtaToPi0ALICE[3].push_back(Alice8TeVGraph->GetEYhigh()[i]);
                currentNumberOfPointALICE++;
            }
        }
        cout << "Number of pT bins: "<< currentNumberOfPoint <<  endl;
        cout << "Number of pT bins ALICE: "<< currentNumberOfPointALICE <<  endl;
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

    TGraphAsymmErrors* graphWorldHighALICE               = new TGraphAsymmErrors(currentNumberOfPointALICE);
    graphWorldHighALICE->SetName("AllWorldDataAbove4GeVALICE");
    for (Int_t k = 0; k< currentNumberOfPointALICE; k++){
        graphWorldHighALICE->SetPoint(k, valuesEtaToPi0ALICE[0].at(k), valuesEtaToPi0ALICE[1].at(k));
        graphWorldHighALICE->SetPointError(k, 0.00, 0.00, valuesEtaToPi0ALICE[2].at(k), valuesEtaToPi0ALICE[3].at(k));
    }
    graphWorldHighALICE->Print();
    graphWorldHighALICE->Sort();
    cout << "after sorting" << endl;
    graphWorldHighALICE->Print();

    TF1* etaToPi0ConstDataALICE  = new TF1("etaToPi0ConstDataALICE","[0]",4,20);
    graphWorldHighALICE->Fit(etaToPi0ConstDataALICE,"QRME0","",4,20);

    cout  << "***********************************************************************************************************" << endl;
    cout  << "ALICE high pt eta/pi0 - data, tot: " << etaToPi0ConstDataALICE->GetParameter(0) << "+-"<< etaToPi0ConstDataALICE->GetParError(0) << endl;
    cout  << "***********************************************************************************************************" << endl;
    cout  << "***********************************************************************************************************" << endl;

    const char* OutputNameWorld ="WorldDataPi0Eta.root";
    TFile* WorldData = new TFile(OutputNameWorld,"RECREATE");

        Donaldson100GeVGraph->Write();
        Donaldson200GeVGraph->Write();;
        Bonesi280GeVGraph->Write();
        Antille87ppGraph->Write();
        Antille87ppbarGraph->Write();
        Aguilar400GeVGraph->Write();
        Amaldi79ppGraph->Write();
        Kourkou79ppGraph->Write();
        Apana530GeVGraph->Write();
        Apana800GeVGraph->Write();
        Kourkou79pp52Graph->Write();
        Akesson53GeVppbarGraph->Write();
        Akesson53GeVppGraph->Write();
        Amaldi53GeVppGraph->Write();
        Kourkou79pp62Graph->Write();
        Akesson63GeVppGraph->Write();
        Phenix200GeVGraph->Write();
        Banner540GeVGraph->Write();
        Alice7TeVGraph->Write();
        Alice7TeVGraphStat->Write();
        Alice7TeVGraphSys->Write();
        if (Alice2760GeVGraph) Alice2760GeVGraph->Write();
        if (Alice2760GeVGraphStat) Alice2760GeVGraphStat->Write();
        if (Alice2760GeVGraphSys) Alice2760GeVGraphSys->Write();
        if (Alice8TeVGraph) Alice8TeVGraph->Write();
        if (Alice8TeVGraphStat) Alice8TeVGraphStat->Write();
        if (Alice8TeVGraphSys) Alice8TeVGraphSys->Write();
        // high pt graphs
        graphWorldHigh->Write();
        etaToPi0ConstData->Write();
        // high pt graphs ALICE
        graphWorldHighALICE->Write();
        etaToPi0ConstDataALICE->Write();

    WorldData->Write();
    WorldData->Close();
}