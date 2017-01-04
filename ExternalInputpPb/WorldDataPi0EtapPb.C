void WorldDataPi0EtapPb(){

// apana03_pBe_plab530GeV_fnal706_eta_pi0_ratio.txt  povlis83_pAl_plab200GeV_fnal629_eta_pi0_ratio.txt  tikhomirov95_pBe_plab450GeV_Helios_eta_pi0_ratio.txt
// apana03_pBe_plab800GeV_fnal706_eta_pi0_ratio.txt  povlis83_pBe_plab200GeV_fnal629_eta_pi0_ratio.txt
// alverson93_pBe_plab500GeV_eta_pi0_ratio.txt              eta_pi0_all_pTsystematics.txt                     povlis83_pC_plab200GeV_fnal629_eta_pi0_ratio.txt
// alverson93_piminusBe_plab500GeV_eta_pi0_ratio.txt        phenix_dAu_200GeV_eta_pi0_ratio.txt               povlis83_piC_plab200GeV_fnal629_eta_pi0_ratio.txt
    
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

    ifstream AgakishievpBe29100MeV;
    AgakishievpBe29100MeV.open("OtherExperiments/agakichiev98_pBe_plab450GeV_tapsceres_eta_pi0_ratio.txt");
    cout << "agakichiev98_pBe_plab450GeV_tapsceres_eta_pi0_ratio.txt" << endl;
    Int_t lines = 0;		
    
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
    
//     //Donaldson200GeV
//     ifstream Donaldson200GeV;
//     Donaldson200GeV.open("OtherExperiments/donaldson78_pp_plab200GeV_eta_pi0_ratio.txt");
//     cout << "donaldson78_pp_plab200GeV_eta_pi0_ratio.txt" << endl;
//     Int_t lines = 0;		
//     
//     while(!Donaldson200GeV.eof()){
//         Donaldson200GeV >> a>> b >> value[lines] >> totErr[lines] >> c >> d;
//         pt[lines] = (a+b)/2;
//         xErr[lines] = 0.;	
//         lines++;
//     }
//     Donaldson200GeV.close();
//     TGraph *Donaldson200GeVGraph = new TGraphErrors(lines-1,pt,value,xErr,totErr);	
//     Donaldson200GeVGraph->SetName("donaldson200GeV");
//     Donaldson200GeVGraph->SetTitle("p+p (#sqrt{#it{s}}= 19.4 GeV)");
//     Donaldson200GeVGraph->Print();
//     
// 
//     // 		"bonesi89_pp_plab280GeV_wa70_eta_pi0_ratio.txt","p+p (#sqrt{#it{s}}= 23 GeV)",
//     ifstream Bonesi280GeV;
//     Bonesi280GeV.open("OtherExperiments/bonesi89_pp_plab280GeV_wa70_eta_pi0_ratio.txt");
//     cout << "bonesi89_pp_plab280GeV_wa70_eta_pi0_ratio.txt" << endl;
//     Int_t lines = 0;		
// 
//     while(!Bonesi280GeV.eof()){
//         Bonesi280GeV >> pt[lines] >> value[lines] >> totErr[lines];
//         xErr[lines] = 0.;	
//         lines++;
//     }
//     Bonesi280GeV.close();
//     TGraph *Bonesi280GeVGraph = new TGraphErrors(lines-1,pt,value,xErr,totErr);	
//     Bonesi280GeVGraph->SetName("Bonesi280GeV");
//     Bonesi280GeVGraph->SetTitle("p+p (#sqrt{#it{s}}= 23 GeV)");
//     Bonesi280GeVGraph->Print();
// 
// 
// // 		"antille87_pp_24.3GeV_eta_pi0_ratio.txt","p+p (#sqrt{#it{s}}= 24.3 GeV)",
//     ifstream Antille87pp;
//     Antille87pp.open("OtherExperiments/antille87_pp_24.3GeV_eta_pi0_ratio.txt");
//     cout << "antille87_pp_24.3GeV_eta_pi0_ratio.txt" << endl;
//     Int_t lines = 0;		
// 
//     while(!Antille87pp.eof()){
//         Antille87pp >> pt[lines] >> value[lines] >> totErr[lines];
//         xErr[lines] = 0.;	
//         lines++;
//     }
//     Antille87pp.close();
//     TGraph *Antille87ppGraph = new TGraphErrors(lines-1,pt,value,xErr,totErr);	
//     Antille87ppGraph->SetName("Antille24.3GeVpp");
//     Antille87ppGraph->SetTitle("p+p (#sqrt{#it{s}}= 24.3 GeV)");
//     Antille87ppGraph->Print();
// 
// // 		"antille87_ppbar_24.3GeV_eta_pi0_ratio.txt","#bar{p}+p (#sqrt{#it{s}}= 24.3 GeV)",
//     ifstream Antille87ppbar;
//     Antille87ppbar.open("OtherExperiments/antille87_ppbar_24.3GeV_eta_pi0_ratio.txt");
//     cout << "antille87_ppbar_24.3GeV_eta_pi0_ratio.txt" << endl;
//     Int_t lines = 0;		
// 
//     while(!Antille87ppbar.eof()){
//         Antille87ppbar >> pt[lines] >> value[lines] >> totErr[lines];
//         xErr[lines] = 0.;	
//         lines++;
//     }
//     Antille87ppbar.close();
//     TGraph *Antille87ppbarGraph = new TGraphErrors(lines-1,pt,value,xErr,totErr);	
//     Antille87ppbarGraph->SetName("Antille24.3GeVppbar");
//     Antille87ppbarGraph->SetTitle("#bar{p}+p (#sqrt{#it{s}}= 24.3 GeV)");
//     Antille87ppbarGraph->Print();
// 
//     // 		"aguilar91_pp_plab400GeV_na27_eta_pi0_ratio.txt","p+p (#sqrt{#it{s}}= 27.5 GeV)",
//     ifstream Aguilar400GeV;
//     Aguilar400GeV.open("OtherExperiments/aguilar91_pp_plab400GeV_na27_eta_pi0_ratio.txt");
//     cout << "aguilar91_pp_plab400GeV_na27_eta_pi0_ratio.txt" << endl;
//     Int_t lines = 0;		
//     
//     while(!Aguilar400GeV.eof()){
//         Aguilar400GeV >> pt[lines] >> value[lines] >> totErr[lines];
//         xErr[lines] = 0.;	
//         lines++;
//     }
//     Aguilar400GeV.close();
//     TGraph *Aguilar400GeVGraph = new TGraphErrors(lines-1,pt,value,xErr,totErr);	
//     Aguilar400GeVGraph->SetName("Aguilar400GeV");
//     Aguilar400GeVGraph->SetTitle("p+p (#sqrt{#it{s}}= 27.5 GeV)");
//     Aguilar400GeVGraph->Print();
//     
// // 		"amaldi79_pp_30.6GeV_isr_eta_pi0_ratio.txt","p+p (#sqrt{#it{s}}= 30.6 GeV)",
//     ifstream Amaldi79pp;
//     Amaldi79pp.open("OtherExperiments/amaldi79_pp_30.6GeV_isr_eta_pi0_ratio.txt");
//     cout << "amaldi79_pp_30.6GeV_isr_eta_pi0_ratio.txt" << endl;
//     Int_t lines = 0;		
// 
//     while(!Amaldi79pp.eof()){
//         Amaldi79pp >> pt[lines] >> value[lines] >> totErr[lines];
//         xErr[lines] = 0.;	
//         lines++;
//     }
//     Amaldi79pp.close();
//     TGraph *Amaldi79ppGraph = new TGraphErrors(lines-1,pt,value,xErr,totErr);	
//     Amaldi79ppGraph->SetName("Amaldi30.6GeVpp");
//     Amaldi79ppGraph->SetTitle("p+p (#sqrt{#it{s}}= 30.6 GeV)");
//     Amaldi79ppGraph->Print();
// 
// // 		"kourkou79_pp_30.6GeV_eta_pi0_ratio.txt","p+p (#sqrt{#it{s}}= 30.6 GeV)",
//     ifstream Kourkou79pp;
//     Kourkou79pp.open("OtherExperiments/kourkou79_pp_30.6GeV_eta_pi0_ratio.txt");
//     cout << "kourkou79_pp_30.6GeV_eta_pi0_ratio.txt" << endl;
//     Int_t lines = 0;		
// 
//     while(!Kourkou79pp.eof()){
//         Kourkou79pp >> pt[lines] >> value[lines] >> totErr[lines];
//         xErr[lines] = 0.;	
//         lines++;
//     }
//     Kourkou79pp.close();
//     TGraph *Kourkou79ppGraph = new TGraphErrors(lines-1,pt,value,xErr,totErr);	
//     Kourkou79ppGraph->SetName("Kourkou30.6GeVpp");
//     Kourkou79ppGraph->SetTitle("p+p (#sqrt{#it{s}}= 30.6 GeV)");
//     Kourkou79ppGraph->Print();
// 
// // 		"apana02_pp_plab530GeV_fnal706_eta_pi0_ratio.txt","p+p (#sqrt{#it{s}}= 31.6 GeV)",
//     ifstream Apana530GeV;
//     Apana530GeV.open("OtherExperiments/apana02_pp_plab530GeV_fnal706_eta_pi0_ratio.txt");
//     cout << "apana02_pp_plab530GeV_fnal706_eta_pi0_ratio.txt" << endl;
//     Int_t lines = 0;		
// 
//     while(!Apana530GeV.eof()){
//         Apana530GeV >> pt[lines] >> value[lines] >> totErr[lines];
//         xErr[lines] = 0.;	
//         lines++;
//     }
//     Apana530GeV.close();
//     TGraph *Apana530GeVGraph = new TGraphErrors(lines-1,pt,value,xErr,totErr);	
//     Apana530GeVGraph->SetName("Apana530GeV");
//     Apana530GeVGraph->SetTitle("p+p (#sqrt{#it{s}}= 31.6 GeV)");
//     Apana530GeVGraph->Print();
// 
//     // 		"apana02_pp_plab800GeV_fnal706_eta_pi0_ratio.txt","p+p (#sqrt{#it{s}}= 38.8 GeV)",
// 
//     ifstream Apana800GeV;
//     Apana800GeV.open("OtherExperiments/apana02_pp_plab800GeV_fnal706_eta_pi0_ratio.txt");
//     cout << "apana02_pp_plab800GeV_fnal706_eta_pi0_ratio.txt" << endl;
//     Int_t lines = 0;		
// 
//     while(!Apana800GeV.eof()){
//         Apana800GeV >> pt[lines] >> value[lines] >> totErr[lines];
//         xErr[lines] = 0.;	
//         lines++;
//     }
//     Apana800GeV.close();
//     TGraph *Apana800GeVGraph = new TGraphErrors(lines-1,pt,value,xErr,totErr);	
//     Apana800GeVGraph->SetName("Apana800GeV");
//     Apana800GeVGraph->SetTitle("p+p (#sqrt{#it{s}}= 38.8 GeV)");
//     Apana800GeVGraph->Print();
// 
// // 		"kourkou79_pp_52.7GeV_eta_pi0_ratio.txt","p+p (#sqrt{#it{s}}= 52.7 GeV)",
//     ifstream Kourkou79pp52;
//     Kourkou79pp52.open("OtherExperiments/kourkou79_pp_52.7GeV_eta_pi0_ratio.txt");
//     cout << "kourkou79_pp_52.7GeV_eta_pi0_ratio.txt" << endl;
//     Int_t lines = 0;		
// 
//     while(!Kourkou79pp52.eof()){
//         Kourkou79pp52 >> pt[lines] >> value[lines] >> totErr[lines];
//         xErr[lines] = 0.;	
//         lines++;
//     }
//     Kourkou79pp52.close();
//     TGraph *Kourkou79pp52Graph = new TGraphErrors(lines-1,pt,value,xErr,totErr);	
//     Kourkou79pp52Graph->SetName("Kourkou52.7GeVpp");
//     Kourkou79pp52Graph->SetTitle("p+p (#sqrt{#it{s}}= 52.7 GeV)");
//     Kourkou79pp52Graph->Print();
// 
// 
// // 		"akesson85_pbarp_53GeV_eta_pi0_ratio.txt","#bar{p}+p (#sqrt{#it{s}}= 53 GeV)",
//     ifstream Akesson53GeVppbar;
//     Akesson53GeVppbar.open("OtherExperiments/akesson85_pbarp_53GeV_eta_pi0_ratio.txt");
//     cout << "akesson85_pbarp_53GeV_eta_pi0_ratio.txt" << endl;
//     Int_t lines = 0;		
//     
//     while(!Akesson53GeVppbar.eof()){
//         Akesson53GeVppbar >> pt[lines] >> value[lines] >> totErr[lines];
//         xErr[lines] = 0.;	
//         lines++;
//     }
//     Akesson53GeVppbar.close();
//     TGraph *Akesson53GeVppbarGraph = new TGraphErrors(lines-1,pt,value,xErr,totErr);	
//     Akesson53GeVppbarGraph->SetName("Akesson53GeVppbar");
//     Akesson53GeVppbarGraph->SetTitle("#bar{p}+p (#sqrt{#it{s}}= 53 GeV)");
//     Akesson53GeVppbarGraph->Print();
//                                                     
//     // 		"akesson85_pp_53GeV_eta_pi0_ratio.txt","p+p (#sqrt{#it{s}}= 53 GeV)",
//     ifstream Akesson53GeVpp;
//     Akesson53GeVpp.open("OtherExperiments/akesson85_pp_53GeV_eta_pi0_ratio.txt");
//     cout << "akesson85_pp_53GeV_eta_pi0_ratio.txt" << endl;
//     Int_t lines = 0;		
//     
//     while(!Akesson53GeVpp.eof()){
//         Akesson53GeVpp >> pt[lines] >> value[lines] >> totErr[lines];
//         xErr[lines] = 0.;	
//         lines++;
//     }
//     Akesson53GeVpp.close();
//     TGraph *Akesson53GeVppGraph = new TGraphErrors(lines-1,pt,value,xErr,totErr);	
//     Akesson53GeVppGraph->SetName("Akesson53GeVpp");
//     Akesson53GeVppGraph->SetTitle("p+p (#sqrt{#it{s}}= 53 GeV)");
//     Akesson53GeVppGraph->Print();
//     
//     // 		"amaldi79_pp_53GeV_isr_eta_pi0_ratio.txt","p+p (#sqrt{#it{s}}= 53.2 GeV)",
//     ifstream Amaldi53GeVpp;
//     Amaldi53GeVpp.open("OtherExperiments/amaldi79_pp_53GeV_isr_eta_pi0_ratio.txt");
//     cout << "amaldi79_pp_53GeV_isr_eta_pi0_ratio.txt" << endl;
//     Int_t lines = 0;		
//     
//     while(!Amaldi53GeVpp.eof()){
//         Amaldi53GeVpp >> pt[lines] >> value[lines] >> totErr[lines];
//         xErr[lines] = 0.;	
//         lines++;
//     }
//     Amaldi53GeVpp.close();
//     TGraph *Amaldi53GeVppGraph = new TGraphErrors(lines-1,pt,value,xErr,totErr);	
//     Amaldi53GeVppGraph->SetName("Amaldi53GeVpp");
//     Amaldi53GeVppGraph->SetTitle("p+p (#sqrt{#it{s}}= 53.2 GeV)");
//     Amaldi53GeVppGraph->Print();
//     
//     
//     // 		"kourkou79_pp_62.4GeV_eta_pi0_ratio.txt","p+p (#sqrt{#it{s}}= 62.4 GeV)",
//     ifstream Kourkou79pp62;
//     Kourkou79pp62.open("OtherExperiments/kourkou79_pp_62.4GeV_eta_pi0_ratio.txt");
//     cout << "kourkou79_pp_62.4GeV_eta_pi0_ratio.txt" << endl;
//     Int_t lines = 0;		
//     
//     while(!Kourkou79pp62.eof()){
//         Kourkou79pp62 >> pt[lines] >> value[lines] >> totErr[lines];
//         xErr[lines] = 0.;	
//         lines++;
//     }
//     Kourkou79pp62.close();
//     TGraph *Kourkou79pp62Graph = new TGraphErrors(lines-1,pt,value,xErr,totErr);	
//     Kourkou79pp62Graph->SetName("Kourkou79pp62");
//     Kourkou79pp62Graph->SetTitle("p+p (#sqrt{#it{s}}= 62.4 GeV)");
//     Kourkou79pp62Graph->Print();
//     
//     // 		"akesson83_pp_63GeV_eta_pi0_ratio.txt","p+p (#sqrt{#it{s}}= 63 GeV)",
//     ifstream Akesson63GeVpp;
//     Akesson63GeVpp.open("OtherExperiments/akesson83_pp_63GeV_eta_pi0_ratio.txt");
//     cout << "akesson83_pp_63GeV_eta_pi0_ratio.txt" << endl;
//     Int_t lines = 0;		
//     while(!Akesson63GeVpp.eof()){
//         Akesson63GeVpp >> a >> b >> value[lines] >> totErr[lines];
//         pt[lines]   = (a+b)/2;
//         xErr[lines] = 0.;	
//         lines++;
//     }
//     Akesson63GeVpp.close();
//     TGraph *Akesson63GeVppGraph = new TGraphErrors(lines-1,pt,value,xErr,totErr);	
//     Akesson63GeVppGraph->SetName("Akesson63GeVpp");
//     Akesson63GeVppGraph->SetTitle("p+p (#sqrt{#it{s}}= 63 GeV)");
//     Akesson63GeVppGraph->Print();
//     
    
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
    Int_t lines = 0;		
    
    
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
//     
//     // 		"banner85_ppbar_540GeV_UA2_eta_pi0_ratio.txt","#bar{p}+p (#sqrt{#it{s}}= 540 GeV)",
//     ifstream Banner540GeV;
//     Banner540GeV.open("OtherExperiments/banner85_ppbar_540GeV_UA2_eta_pi0_ratio.txt");
//     cout << "banner85_ppbar_540GeV_UA2_eta_pi0_ratio.txt" << endl;
//     Int_t lines = 0;		
//     
//     while(!Banner540GeV.eof()){
//         Banner540GeV >> pt[lines] >> value[lines] >> totErr[lines];
//         xErr[lines] = 0.;	
//         lines++;
//     }
//     Banner540GeV.close();
//     TGraph *Banner540GeVGraph = new TGraphErrors(lines-1,pt,value,xErr,totErr);	
//     Banner540GeVGraph->SetName("Banner540GeV");
//     Banner540GeVGraph->SetTitle("#bar{p}+p (#sqrt{#it{s}}= 540 GeV)");
//     Banner540GeVGraph->Print();
    

    const char* OutputNameWorld ="WorldDataPi0EtapPb.root";
    WorldData = new TFile(OutputNameWorld,"RECREATE");		

        AgakishievpAu29100MeVGraph->Write();
        AgakishievpBe29100MeVGraph->Write();
//         Donaldson200GeVGraph->Write();;
//         Bonesi280GeVGraph->Write();
//         Antille87ppGraph->Write();
//         Antille87ppbarGraph->Write();
//         Aguilar400GeVGraph->Write();
//         Amaldi79ppGraph->Write();
//         Kourkou79ppGraph->Write();
//         Apana530GeVGraph->Write();
//         Apana800GeVGraph->Write();
//         Kourkou79pp52Graph->Write();
//         Akesson53GeVppbarGraph->Write();
//         Akesson53GeVppGraph->Write();
//         Amaldi53GeVppGraph->Write();
//         Kourkou79pp62Graph->Write();
//         Akesson63GeVppGraph->Write();
        PhenixdAu200GeVGraph->Write();
//         Banner540GeVGraph->Write();
    WorldData->Write();
    WorldData->Close();
    
    
	
}