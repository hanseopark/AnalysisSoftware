void WorldDataPi0Eta(void){

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

	Double_t pt[100];
	Double_t value[100];
	Double_t error[100];
	Double_t xerror[100];
	
	
	//Donaldson100GeV
	ifstream Donaldson100GeV;
	Donaldson100GeV.open("donaldson78_pp_plab100GeV_eta_pi0_ratio.txt");
	cout << "donaldson78_pp_plab100GeV_eta_pi0_ratio.txt" << endl;
	Int_t lines = 0;		
	
	while(!Donaldson100GeV.eof()){
		Donaldson100GeV >> pt[lines] >> value[lines] >> error[lines];
		xerror[lines] = 0.;	
		lines++;
	}
	Donaldson100GeV.close();
	
	TGraph *Donaldson100GeVGraph = new TGraphErrors(lines-1,pt,value,xerror,error);	
	Donaldson100GeVGraph->SetName("donaldson100GeV");
	Donaldson100GeVGraph->SetTitle("p+p (#sqrt{#it{s}}= 13.8 GeV)");
	Donaldson100GeVGraph->Print();

	//Donaldson200GeV
	ifstream Donaldson200GeV;
	Donaldson200GeV.open("donaldson78_pp_plab200GeV_eta_pi0_ratio.txt");
	cout << "donaldson78_pp_plab200GeV_eta_pi0_ratio.txt" << endl;
	Int_t lines = 0;		
	
	while(!Donaldson200GeV.eof()){
		Donaldson200GeV >> pt[lines] >> value[lines] >> error[lines];
		xerror[lines] = 0.;	
		lines++;
	}
	Donaldson200GeV.close();
	TGraph *Donaldson200GeVGraph = new TGraphErrors(lines-1,pt,value,xerror,error);	
	Donaldson200GeVGraph->SetName("donaldson200GeV");
	Donaldson200GeVGraph->SetTitle("p+p (#sqrt{#it{s}}= 19.4 GeV)");
	Donaldson200GeVGraph->Print();
	

	// 		"bonesi89_pp_plab280GeV_wa70_eta_pi0_ratio.txt","p+p (#sqrt{#it{s}}= 23 GeV)",
	ifstream Bonesi280GeV;
	Bonesi280GeV.open("bonesi89_pp_plab280GeV_wa70_eta_pi0_ratio.txt");
	cout << "bonesi89_pp_plab280GeV_wa70_eta_pi0_ratio.txt" << endl;
	Int_t lines = 0;		

	while(!Bonesi280GeV.eof()){
		Bonesi280GeV >> pt[lines] >> value[lines] >> error[lines];
		xerror[lines] = 0.;	
		lines++;
	}
	Bonesi280GeV.close();
	TGraph *Bonesi280GeVGraph = new TGraphErrors(lines-1,pt,value,xerror,error);	
	Bonesi280GeVGraph->SetName("Bonesi280GeV");
	Bonesi280GeVGraph->SetTitle("p+p (#sqrt{#it{s}}= 23 GeV)");
	Bonesi280GeVGraph->Print();


// 		"antille87_pp_24.3GeV_eta_pi0_ratio.txt","p+p (#sqrt{#it{s}}= 24.3 GeV)",
	ifstream Antille87pp;
	Antille87pp.open("antille87_pp_24.3GeV_eta_pi0_ratio.txt");
	cout << "antille87_pp_24.3GeV_eta_pi0_ratio.txt" << endl;
	Int_t lines = 0;		

	while(!Antille87pp.eof()){
		Antille87pp >> pt[lines] >> value[lines] >> error[lines];
		xerror[lines] = 0.;	
		lines++;
	}
	Antille87pp.close();
	TGraph *Antille87ppGraph = new TGraphErrors(lines-1,pt,value,xerror,error);	
	Antille87ppGraph->SetName("Antille24.3GeVpp");
	Antille87ppGraph->SetTitle("p+p (#sqrt{#it{s}}= 24.3 GeV)");
	Antille87ppGraph->Print();

// 		"antille87_ppbar_24.3GeV_eta_pi0_ratio.txt","#bar{p}+p (#sqrt{#it{s}}= 24.3 GeV)",
	ifstream Antille87ppbar;
	Antille87ppbar.open("antille87_ppbar_24.3GeV_eta_pi0_ratio.txt");
	cout << "antille87_ppbar_24.3GeV_eta_pi0_ratio.txt" << endl;
	Int_t lines = 0;		

	while(!Antille87ppbar.eof()){
		Antille87ppbar >> pt[lines] >> value[lines] >> error[lines];
		xerror[lines] = 0.;	
		lines++;
	}
	Antille87ppbar.close();
	TGraph *Antille87ppbarGraph = new TGraphErrors(lines-1,pt,value,xerror,error);	
	Antille87ppbarGraph->SetName("Antille24.3GeVppbar");
	Antille87ppbarGraph->SetTitle("#bar{p}+p (#sqrt{#it{s}}= 24.3 GeV)");
	Antille87ppbarGraph->Print();

	// 		"aguilar91_pp_plab400GeV_na27_eta_pi0_ratio.txt","p+p (#sqrt{#it{s}}= 27.5 GeV)",
	ifstream Aguilar400GeV;
	Aguilar400GeV.open("aguilar91_pp_plab400GeV_na27_eta_pi0_ratio.txt");
	cout << "aguilar91_pp_plab400GeV_na27_eta_pi0_ratio.txt" << endl;
	Int_t lines = 0;		
	
	while(!Aguilar400GeV.eof()){
		Aguilar400GeV >> pt[lines] >> value[lines] >> error[lines];
		xerror[lines] = 0.;	
		lines++;
	}
	Aguilar400GeV.close();
	TGraph *Aguilar400GeVGraph = new TGraphErrors(lines-1,pt,value,xerror,error);	
	Aguilar400GeVGraph->SetName("Aguilar400GeV");
	Aguilar400GeVGraph->SetTitle("p+p (#sqrt{#it{s}}= 27.5 GeV)");
	Aguilar400GeVGraph->Print();
	
// 		"amaldi79_pp_30.6GeV_isr_eta_pi0_ratio.txt","p+p (#sqrt{#it{s}}= 30.6 GeV)",
	ifstream Amaldi79pp;
	Amaldi79pp.open("amaldi79_pp_30.6GeV_isr_eta_pi0_ratio.txt");
	cout << "amaldi79_pp_30.6GeV_isr_eta_pi0_ratio.txt" << endl;
	Int_t lines = 0;		

	while(!Amaldi79pp.eof()){
		Amaldi79pp >> pt[lines] >> value[lines] >> error[lines];
		xerror[lines] = 0.;	
		lines++;
	}
	Amaldi79pp.close();
	TGraph *Amaldi79ppGraph = new TGraphErrors(lines-1,pt,value,xerror,error);	
	Amaldi79ppGraph->SetName("Amaldi30.6GeVpp");
	Amaldi79ppGraph->SetTitle("p+p (#sqrt{#it{s}}= 30.6 GeV)");
	Amaldi79ppGraph->Print();

// 		"kourkou79_pp_30.6GeV_eta_pi0_ratio.txt","p+p (#sqrt{#it{s}}= 30.6 GeV)",
	ifstream Kourkou79pp;
	Kourkou79pp.open("kourkou79_pp_30.6GeV_eta_pi0_ratio.txt");
	cout << "kourkou79_pp_30.6GeV_eta_pi0_ratio.txt" << endl;
	Int_t lines = 0;		

	while(!Kourkou79pp.eof()){
		Kourkou79pp >> pt[lines] >> value[lines] >> error[lines];
		xerror[lines] = 0.;	
		lines++;
	}
	Kourkou79pp.close();
	TGraph *Kourkou79ppGraph = new TGraphErrors(lines-1,pt,value,xerror,error);	
	Kourkou79ppGraph->SetName("Kourkou30.6GeVpp");
	Kourkou79ppGraph->SetTitle("p+p (#sqrt{#it{s}}= 30.6 GeV)");
	Kourkou79ppGraph->Print();

// 		"apana02_pp_plab530GeV_fnal706_eta_pi0_ratio.txt","p+p (#sqrt{#it{s}}= 31.6 GeV)",
	ifstream Apana530GeV;
	Apana530GeV.open("apana02_pp_plab530GeV_fnal706_eta_pi0_ratio.txt");
	cout << "apana02_pp_plab530GeV_fnal706_eta_pi0_ratio.txt" << endl;
	Int_t lines = 0;		

	while(!Apana530GeV.eof()){
		Apana530GeV >> pt[lines] >> value[lines] >> error[lines];
		xerror[lines] = 0.;	
		lines++;
	}
	Apana530GeV.close();
	TGraph *Apana530GeVGraph = new TGraphErrors(lines-1,pt,value,xerror,error);	
	Apana530GeVGraph->SetName("Apana530GeV");
	Apana530GeVGraph->SetTitle("p+p (#sqrt{#it{s}}= 31.6 GeV)");
	Apana530GeVGraph->Print();

	// 		"apana02_pp_plab800GeV_fnal706_eta_pi0_ratio.txt","p+p (#sqrt{#it{s}}= 38.8 GeV)",

	ifstream Apana800GeV;
	Apana800GeV.open("apana02_pp_plab800GeV_fnal706_eta_pi0_ratio.txt");
	cout << "apana02_pp_plab800GeV_fnal706_eta_pi0_ratio.txt" << endl;
	Int_t lines = 0;		

	while(!Apana800GeV.eof()){
		Apana800GeV >> pt[lines] >> value[lines] >> error[lines];
		xerror[lines] = 0.;	
		lines++;
	}
	Apana800GeV.close();
	TGraph *Apana800GeVGraph = new TGraphErrors(lines-1,pt,value,xerror,error);	
	Apana800GeVGraph->SetName("Apana800GeV");
	Apana800GeVGraph->SetTitle("p+p (#sqrt{#it{s}}= 38.8 GeV)");
	Apana800GeVGraph->Print();

// 		"kourkou79_pp_52.7GeV_eta_pi0_ratio.txt","p+p (#sqrt{#it{s}}= 52.7 GeV)",
	ifstream Kourkou79pp52;
	Kourkou79pp52.open("kourkou79_pp_52.7GeV_eta_pi0_ratio.txt");
	cout << "kourkou79_pp_52.7GeV_eta_pi0_ratio.txt" << endl;
	Int_t lines = 0;		

	while(!Kourkou79pp52.eof()){
		Kourkou79pp52 >> pt[lines] >> value[lines] >> error[lines];
		xerror[lines] = 0.;	
		lines++;
	}
	Kourkou79pp52.close();
	TGraph *Kourkou79pp52Graph = new TGraphErrors(lines-1,pt,value,xerror,error);	
	Kourkou79pp52Graph->SetName("Kourkou52.7GeVpp");
	Kourkou79pp52Graph->SetTitle("p+p (#sqrt{#it{s}}= 52.7 GeV)");
	Kourkou79pp52Graph->Print();


// 		"akesson85_pbarp_53GeV_eta_pi0_ratio.txt","#bar{p}+p (#sqrt{#it{s}}= 53 GeV)",
	ifstream Akesson53GeVppbar;
	Akesson53GeVppbar.open("akesson85_pbarp_53GeV_eta_pi0_ratio.txt");
	cout << "akesson85_pbarp_53GeV_eta_pi0_ratio.txt" << endl;
	Int_t lines = 0;		
	
	while(!Akesson53GeVppbar.eof()){
		Akesson53GeVppbar >> pt[lines] >> value[lines] >> error[lines];
		xerror[lines] = 0.;	
		lines++;
	}
	Akesson53GeVppbar.close();
	TGraph *Akesson53GeVppbarGraph = new TGraphErrors(lines-1,pt,value,xerror,error);	
	Akesson53GeVppbarGraph->SetName("Akesson53GeVppbar");
	Akesson53GeVppbarGraph->SetTitle("#bar{p}+p (#sqrt{#it{s}}= 53 GeV)");
	Akesson53GeVppbarGraph->Print();
						    					    
	// 		"akesson85_pp_53GeV_eta_pi0_ratio.txt","p+p (#sqrt{#it{s}}= 53 GeV)",
	ifstream Akesson53GeVpp;
	Akesson53GeVpp.open("akesson85_pp_53GeV_eta_pi0_ratio.txt");
	cout << "akesson85_pp_53GeV_eta_pi0_ratio.txt" << endl;
	Int_t lines = 0;		
	
	while(!Akesson53GeVpp.eof()){
		Akesson53GeVpp >> pt[lines] >> value[lines] >> error[lines];
		xerror[lines] = 0.;	
		lines++;
	}
	Akesson53GeVpp.close();
	TGraph *Akesson53GeVppGraph = new TGraphErrors(lines-1,pt,value,xerror,error);	
	Akesson53GeVppGraph->SetName("Akesson53GeVpp");
	Akesson53GeVppGraph->SetTitle("p+p (#sqrt{#it{s}}= 53 GeV)");
	Akesson53GeVppGraph->Print();
	
	// 		"amaldi79_pp_53GeV_isr_eta_pi0_ratio.txt","p+p (#sqrt{#it{s}}= 53.2 GeV)",
	ifstream Amaldi53GeVpp;
	Amaldi53GeVpp.open("amaldi79_pp_53GeV_isr_eta_pi0_ratio.txt");
	cout << "amaldi79_pp_53GeV_isr_eta_pi0_ratio.txt" << endl;
	Int_t lines = 0;		
	
	while(!Amaldi53GeVpp.eof()){
		Amaldi53GeVpp >> pt[lines] >> value[lines] >> error[lines];
		xerror[lines] = 0.;	
		lines++;
	}
	Amaldi53GeVpp.close();
	TGraph *Amaldi53GeVppGraph = new TGraphErrors(lines-1,pt,value,xerror,error);	
	Amaldi53GeVppGraph->SetName("Amaldi53GeVpp");
	Amaldi53GeVppGraph->SetTitle("p+p (#sqrt{#it{s}}= 53.2 GeV)");
	Amaldi53GeVppGraph->Print();
	
	
	// 		"kourkou79_pp_62.4GeV_eta_pi0_ratio.txt","p+p (#sqrt{#it{s}}= 62.4 GeV)",
	ifstream Kourkou79pp62;
	Kourkou79pp62.open("kourkou79_pp_62.4GeV_eta_pi0_ratio.txt");
	cout << "kourkou79_pp_62.4GeV_eta_pi0_ratio.txt" << endl;
	Int_t lines = 0;		
	
	while(!Kourkou79pp62.eof()){
		Kourkou79pp62 >> pt[lines] >> value[lines] >> error[lines];
		xerror[lines] = 0.;	
		lines++;
	}
	Kourkou79pp62.close();
	TGraph *Kourkou79pp62Graph = new TGraphErrors(lines-1,pt,value,xerror,error);	
	Kourkou79pp62Graph->SetName("Kourkou79pp62");
	Kourkou79pp62Graph->SetTitle("p+p (#sqrt{#it{s}}= 62.4 GeV)");
	Kourkou79pp62Graph->Print();
	
	// 		"akesson83_pp_63GeV_eta_pi0_ratio.txt","p+p (#sqrt{#it{s}}= 63 GeV)",
	ifstream Akesson63GeVpp;
	Akesson63GeVpp.open("akesson83_pp_63GeV_eta_pi0_ratio.txt");
	cout << "akesson83_pp_63GeV_eta_pi0_ratio.txt" << endl;
	Int_t lines = 0;		
	
	while(!Akesson63GeVpp.eof()){
		Akesson63GeVpp >> pt[lines] >> value[lines] >> error[lines];
		xerror[lines] = 0.;	
		lines++;
	}
	Akesson63GeVpp.close();
	TGraph *Akesson63GeVppGraph = new TGraphErrors(lines-1,pt,value,xerror,error);	
	Akesson63GeVppGraph->SetName("Akesson63GeVpp");
	Akesson63GeVppGraph->SetTitle("p+p (#sqrt{#it{s}}= 63 GeV)");
	Akesson63GeVppGraph->Print();
	
	// 		"phenix_pp_200GeV_eta_pi0_ratio.txt","p+p (#sqrt{#it{s}}= 200 GeV)",
	ifstream Phenix200GeV;
	Phenix200GeV.open("phenix_pp_200GeV_eta_pi0_ratio.txt");
	cout << "phenix_pp_200GeV_eta_pi0_ratio.txt" << endl;
	Int_t lines = 0;		
	
	while(!Phenix200GeV.eof()){
		Phenix200GeV >> pt[lines] >> value[lines] >> error[lines];
		xerror[lines] = 0.;	
		lines++;
	}
	Phenix200GeV.close();
	TGraph *Phenix200GeVGraph = new TGraphErrors(lines-1,pt,value,xerror,error);	
	Phenix200GeVGraph->SetName("Phenix200GeV");
	Phenix200GeVGraph->SetTitle("p+p (#sqrt{#it{s}}= 200 GeV)");
	Phenix200GeVGraph->Print();
	
	// 		"banner85_ppbar_540GeV_UA2_eta_pi0_ratio.txt","#bar{p}+p (#sqrt{#it{s}}= 540 GeV)",
	ifstream Banner540GeV;
	Banner540GeV.open("banner85_ppbar_540GeV_UA2_eta_pi0_ratio.txt");
	cout << "banner85_ppbar_540GeV_UA2_eta_pi0_ratio.txt" << endl;
	Int_t lines = 0;		
	
	while(!Banner540GeV.eof()){
		Banner540GeV >> pt[lines] >> value[lines] >> error[lines];
		xerror[lines] = 0.;	
		lines++;
	}
	Banner540GeV.close();
	TGraph *Banner540GeVGraph = new TGraphErrors(lines-1,pt,value,xerror,error);	
	Banner540GeVGraph->SetName("Banner540GeV");
	Banner540GeVGraph->SetTitle("#bar{p}+p (#sqrt{#it{s}}= 540 GeV)");
	Banner540GeVGraph->Print();
	

	
	
	const char* OutputNameWorld ="WorldDataPi0Eta.root";
	WorldData = new TFile(OutputNameWorld,"RECREATE");		

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

	WorldData->Write();
	WorldData->Close();
	
	
	
}