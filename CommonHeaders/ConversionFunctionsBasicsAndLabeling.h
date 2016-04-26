//**************************************************************************************************************************
//******************************** Conversion functions header for basics and plot labeling ********************************
//**************************************************************************************************************************
//*********** this header contains all basic functions for the analysis as well as the correct labeling of plots ***********
//**************************************************************************************************************************
#include "TObjString.h"

//********************* global Variables *********************************************************************
Double_t ncoll0005      = 1684.4;
Double_t ncoll0510      = 1316;
Double_t ncoll0010      = 1500;
Double_t ncoll1020      = 921.2;
Double_t ncoll0020      = 1210.6;
Double_t ncoll2040      = 438.4;
Double_t ncoll2050      = 349.1;
Double_t ncoll4060      = 127.7;
Double_t ncoll4050      = 171.25;
Double_t ncoll5060      = 84.28;
Double_t ncoll6080      = 26.71;
Double_t ncoll4080      = 77.205;
Double_t ncoll6070      = 37.855;
Double_t ncoll7080      = 15.575;
Double_t ncoll8090      = 6.293;
Double_t ncoll7590      = 8.219;

Double_t nCollErr0005   = 190;
Double_t nCollErr0510   = 140;
Double_t nCollErr0010   = 165;
Double_t nCollErr1020   = 96;
Double_t nCollErr0020   = 130.5;
Double_t nCollErr2040   = 42.;
Double_t nCollErr2050   = 51.;
Double_t nCollErr4060   = 11;
Double_t nCollErr4080   = 6.216;
Double_t nCollErr4050   = 16;
Double_t nCollErr5060   = 6.95;
Double_t nCollErr6080   = 2;
Double_t nCollErr6070   = 2.85;
Double_t nCollErr7080   = 1.035;
Double_t nCollErr8090   = 0.325;
Double_t nCollErr7590   = 0.473;

Double_t tAA0005        = 26.32;
Double_t tAA0510        = 20.56;
Double_t tAA0010        = 23.44;
Double_t tAA1020        = 14.39;
Double_t tAA0020        = 18.915;
Double_t tAA2040        = 6.85;
Double_t tAA2050        = 5.46;
Double_t tAA4060        = 1.996;
Double_t tAA6080        = 0.4174;

Double_t tAAErr0005     = 0.84224;
Double_t tAAErr0510     = 0.65792;
Double_t tAAErr0010     = 0.75008;
Double_t tAAErr1020     = 0.44609;
Double_t tAAErr0020     = 0.5958225;
Double_t tAAErr2040     = 0.22605;
Double_t tAAErr2050     = 0.195;
Double_t tAAErr4060     = 0.097804;
Double_t tAAErr6080     = 0.026296;


//************************************************************************************
//********************* Separate cut numbers, old version ****************************
//************************************************************************************
void ReturnSeparatedCutNumber(  TString cutSel, 
                                TString& gammaCutNumber,
                                TString& electronCutNumber,
                                TString& mesonCutNumber, 
                                Bool_t kDalitz=kFALSE){

    TObjArray *arr;
    arr = cutSel.Tokenize("_");
    TObjString* objstrGamma;
    TObjString* objstrElectron;
    TObjString* objstrMeson;

    if (kDalitz){
        objstrGamma = (TObjString*)arr->At(0);
        objstrElectron = (TObjString*)arr->At(1);
        objstrMeson = (TObjString*)arr->At(2);
        
        gammaCutNumber= objstrGamma->GetString();
        electronCutNumber = objstrElectron->GetString();
        mesonCutNumber = objstrMeson->GetString();
    }
    else {
        objstrGamma = (TObjString*)arr->At(0);
        objstrMeson = (TObjString*)arr->At(1);
        
        gammaCutNumber= objstrGamma->GetString();
        mesonCutNumber = objstrMeson->GetString();
    }
    cout << cutSel.Data() << "\t" << gammaCutNumber.Data() << "\t" << electronCutNumber.Data() << "\t" << mesonCutNumber.Data() << endl;
    return;
}

//************************************************************************************
//********************* Separate cut numbers, new version ****************************
//** separates the possible combinations according to the mode into the             ** 
//** respective eventCutNumber, gammaCutNumber (PCM), clusterCutNumber (EMCal/PHOS),**
//** electronCutNumber (Dalitz) and mesonCutNumber                                  **
//************************************************************************************
void ReturnSeparatedCutNumberAdvanced(  TString cutSel, 
                                        TString& eventCutNumber, 
                                        TString& gammaCutNumber, 
                                        TString& clusterCutNumber, 
                                        TString& electronCutNumber,
                                        TString& mesonCutNumber,
                                        Int_t type=0){
    
    TObjArray *arr;
    arr = cutSel.Tokenize("_");
    TObjString* objstrEvent;
    TObjString* objstrGamma;
    TObjString* objstrCluster;
    TObjString* objstrElectron;
    TObjString* objstrMeson;
    
    if (type == 0){ // PCM-PCM
        objstrEvent         = (TObjString*)arr->At(0);
        objstrGamma         = (TObjString*)arr->At(1);
        objstrMeson         = (TObjString*)arr->At(2);
        
        eventCutNumber      = objstrEvent->GetString();
        gammaCutNumber      = objstrGamma->GetString();
        mesonCutNumber      = objstrMeson->GetString();

    } else if (type == 1){ //PCM dalitz
        objstrEvent         = (TObjString*)arr->At(0);
        objstrGamma         = (TObjString*)arr->At(1);
        objstrElectron      = (TObjString*)arr->At(2);
        objstrMeson         = (TObjString*)arr->At(3);
        
        eventCutNumber      = objstrEvent->GetString();
        gammaCutNumber      = objstrGamma->GetString();
        electronCutNumber   = objstrElectron->GetString();
        mesonCutNumber      = objstrMeson->GetString();
    } else if (type == 2){ //PCM-EMCal
        objstrEvent         = (TObjString*)arr->At(0);
        objstrGamma         = (TObjString*)arr->At(1);
        objstrCluster       = (TObjString*)arr->At(2);
        objstrMeson         = (TObjString*)arr->At(3);
        
        eventCutNumber      = objstrEvent->GetString();
        gammaCutNumber      = objstrGamma->GetString();
        clusterCutNumber    = objstrCluster->GetString();
        mesonCutNumber      = objstrMeson->GetString();
    } else if (type == 3){ //PCM-PHOS
        objstrEvent         = (TObjString*)arr->At(0);
        objstrGamma         = (TObjString*)arr->At(1);
        objstrCluster       = (TObjString*)arr->At(2);
        objstrMeson         = (TObjString*)arr->At(3);
        
        eventCutNumber      = objstrEvent->GetString();
        gammaCutNumber      = objstrGamma->GetString();
        clusterCutNumber    = objstrCluster->GetString();
        mesonCutNumber      = objstrMeson->GetString();
    }  else if (type == 4){ //EMCal-EMCal
        objstrEvent         = (TObjString*)arr->At(0);
        objstrCluster       = (TObjString*)arr->At(1);
        objstrMeson         = (TObjString*)arr->At(2);
        
        eventCutNumber      = objstrEvent->GetString();
        clusterCutNumber    = objstrCluster->GetString();
        mesonCutNumber      = objstrMeson->GetString();
    }  else if (type == 5){ //PHOS-PHOS
        objstrEvent         = (TObjString*)arr->At(0);
        objstrCluster       = (TObjString*)arr->At(1);
        objstrMeson         = (TObjString*)arr->At(2);
        
        eventCutNumber      = objstrEvent->GetString();
        clusterCutNumber    = objstrCluster->GetString();
        mesonCutNumber      = objstrMeson->GetString();
    } else if (type == 6 ){ //Dalitz-EMCal
        objstrEvent         = (TObjString*)arr->At(0);
        objstrGamma         = (TObjString*)arr->At(1);
        objstrElectron      = (TObjString*)arr->At(2);
        objstrCluster       = (TObjString*)arr->At(3);
        objstrMeson         = (TObjString*)arr->At(4);
        
        eventCutNumber      = objstrEvent->GetString();
        gammaCutNumber      = objstrGamma->GetString();
        clusterCutNumber    = objstrCluster->GetString();
        electronCutNumber   = objstrElectron->GetString();
        mesonCutNumber      = objstrMeson->GetString();
    } else if ( type == 7 ){ //Dalitz-PHOS
        objstrEvent         = (TObjString*)arr->At(0);
        objstrGamma         = (TObjString*)arr->At(1);
        objstrElectron      = (TObjString*)arr->At(2);
        objstrCluster       = (TObjString*)arr->At(3);
        objstrMeson         = (TObjString*)arr->At(4);
        
        eventCutNumber      = objstrEvent->GetString();
        gammaCutNumber      = objstrGamma->GetString();
        clusterCutNumber    = objstrCluster->GetString();
        electronCutNumber   = objstrElectron->GetString();
        mesonCutNumber      = objstrMeson->GetString();
    } else if ( type == 10 ){ // EMCal merged
        objstrEvent         = (TObjString*)arr->At(0);
        objstrGamma         = (TObjString*)arr->At(1);
        objstrCluster       = (TObjString*)arr->At(2);
        objstrMeson         = (TObjString*)arr->At(3);

        eventCutNumber      = objstrEvent->GetString();
        gammaCutNumber      = objstrGamma->GetString();
        clusterCutNumber    = objstrCluster->GetString();
        mesonCutNumber      = objstrMeson->GetString();
    } else if ( type == 11 ){ // PHOS merged
        objstrEvent         = (TObjString*)arr->At(0);
        objstrGamma         = (TObjString*)arr->At(1);
        objstrCluster       = (TObjString*)arr->At(2);
        objstrMeson         = (TObjString*)arr->At(3);

        eventCutNumber      = objstrEvent->GetString();
        gammaCutNumber      = objstrGamma->GetString();
        clusterCutNumber    = objstrCluster->GetString();
        mesonCutNumber      = objstrMeson->GetString();        
    }    
    
    cout << cutSel.Data() << "\t" << eventCutNumber.Data() << "\t" << gammaCutNumber.Data() << "\t" <<  clusterCutNumber.Data() << "\t" <<electronCutNumber.Data() << "\t" << mesonCutNumber.Data() << endl;
    return;
}

//************************************************************************************
//********************* get number of events for PCM/calo analysis *******************
//************************************************************************************
Int_t GetNEvents (  TH1* histo, 
                    Bool_t doCout=kTRUE){
    if (!histo) cout << "NO EVENT HISTO" << endl;
    if(histo->GetNbinsX()==11){
      if(histo->GetEntries()-histo->GetBinContent(5)-histo->GetBinContent(7)-histo->GetBinContent(4)==0) return 0;
      Int_t nEvents = histo->GetBinContent(1)+(histo->GetBinContent(1)/(histo->GetBinContent(1)+histo->GetBinContent(5)))*histo->GetBinContent(6);
      Int_t nEventsMB = histo->GetEntries()-histo->GetBinContent(4) -histo->GetBinContent(8)-histo->GetBinContent(9);
      for (Int_t i = 1; i<12; i++ ){
            cout << histo->GetBinContent(i) << "\t";
      }    
      cout << nEventsMB  << endl;
      for (Int_t i = 1; i<12; i++ ){
            cout << histo->GetXaxis()->GetBinLabel(i) << "\t" << histo->GetBinContent(i)/nEventsMB << "\n";
      }
      cout  << endl;
      // BinContent 1 - good events
      // BinContent 2 - centrality not selected
      // BinContent 3 - MC corrupt
      // BinContent 4 - no Trigger Bit
      // BinContent 5 - Zvertex-position,
      // BinContent 6 - no Contributors to vtx
      // BinContent 7 - PileUp
      // BinContent 8 - no SDD
      // BinContent 9 - no V0AND
      if(doCout)cout <<"nEvents new: "<< nEvents << "\t nEvents old: "<< histo->GetEntries()-histo->GetBinContent(5)-histo->GetBinContent(7)-histo->GetBinContent(4)<< endl;
      return nEvents;
    }else if(histo->GetNbinsX()==12){
      if(histo->GetEntries()-histo->GetBinContent(5)-histo->GetBinContent(7)-histo->GetBinContent(12)-histo->GetBinContent(4)==0) return 0;
      Int_t nEvents = histo->GetBinContent(1)+(histo->GetBinContent(1)/(histo->GetBinContent(1)+histo->GetBinContent(5)))*histo->GetBinContent(6);
      Int_t nEventsMB = histo->GetEntries()-histo->GetBinContent(4) -histo->GetBinContent(8)-histo->GetBinContent(9);
      for (Int_t i = 1; i<13; i++ ){
            cout << histo->GetBinContent(i) << "\t";
      }    
      cout << nEventsMB  << endl;
      for (Int_t i = 1; i<13; i++ ){
            cout << histo->GetXaxis()->GetBinLabel(i) << "\t" << histo->GetBinContent(i)/nEventsMB << "\n";
      }
      cout << "accepted \t" << (Float_t)nEvents/nEventsMB << endl;
      cout << endl;
      // BinContent 1 - good events
      // BinContent 2 - centrality not selected
      // BinContent 3 - MC corrupt
      // BinContent 4 - no Trigger Bit
      // BinContent 5 - Zvertex-position,
      // BinContent 6 - no Contributors to vtx
      // BinContent 7 - PileUp
      // BinContent 8 - no SDD
      // BinContent 9 - no V0AND
      // BinContent 12 - SPD cluster vs tracklets
      if(doCout)cout <<"nEvents new: "<< nEvents <<  endl;
      return nEvents;
    }else{
      cout << "ERROR: GetNEvents, dimension of histogram not known! Returning 0...!" << endl;
      return 0;
    }
}

//************************************************************************************
//**** Decodes from the mode the respective reco process and return correct label ****
//************************************************************************************
TString ReturnTextReconstructionProcess(Int_t mode){
    switch (mode){
        case 0:
            return "PCM";
        case 1: 
            return "PCM, Dalitz";
        case 2:
            return "PCM, EMCal";
        case 3: 
            return "PCM, PHOS";
        case 4:
            return "EMCal, EMCal";
        case 5: 
            return "PHOS, PHOS";
        case 6:
            return "Dalitz, EMCal";
        case 7:
            return "Dalitz, PHOS";
        default:
            return "not known";
    }
}


//************************************************************************************
//************************ return latex writing of meson name ************************
//************************************************************************************
TString ReturnMesonString ( TString mesonName){
    if(mesonName.CompareTo("Pi0") == 0 || mesonName.CompareTo("Pi0EtaBinning") == 0){
        return "#pi^{0}";
    } else if(mesonName.CompareTo("Eta") == 0)  {
        return "#eta";
    } else if(mesonName.CompareTo("EtaPrim") == 0)  {
        return "#eta'";
    } else {
        cout << "No correct meson has been selected" << endl;
        return "";
    }
}

//************************************************************************************
//************************ return bool for respective meson name *********************
//************************************************************************************
Bool_t ReturnMesonOption ( TString mesonName){
    if(mesonName.CompareTo("Pi0") == 0 || mesonName.CompareTo("Pi0EtaBinning") == 0){
        return kTRUE;
    } else if(mesonName.CompareTo("Eta") == 0)  {
        return kFALSE;
    } else {
        return kFALSE;
    }
}


//************************************************************************************
//***** return detailed labels for different processes depending on energy ***********
//************************************************************************************
TString ReturnFullTextMeson(TString fEnergyFlagOpt, 
                            TString textProcessOpt){ 

    if(fEnergyFlagOpt.CompareTo("7TeV") == 0){
        return Form("pp #rightarrow %s (#rightarrow #gamma#gamma #rightarrow e^{+}e^{-}e^{+}e^{-}) + X @ 7 TeV ",textProcessOpt.Data());
    } else if(fEnergyFlagOpt.CompareTo("8TeV") == 0){
        return Form("pp #rightarrow %s (#rightarrow #gamma#gamma #rightarrow e^{+}e^{-}e^{+}e^{-}) + X @ 8 TeV ",textProcessOpt.Data());
    } else if(fEnergyFlagOpt.CompareTo("13TeV") == 0){
        return Form("pp #rightarrow %s (#rightarrow #gamma#gamma #rightarrow e^{+}e^{-}e^{+}e^{-}) + X @ 13 TeV ",textProcessOpt.Data());
    } else if( fEnergyFlagOpt.CompareTo("900GeV") == 0) {
        return  Form("pp #rightarrow %s (#rightarrow #gamma#gamma #rightarrow e^{+}e^{-}e^{+}e^{-}) + X @ 900 GeV ",textProcessOpt.Data());
    } else if( fEnergyFlagOpt.CompareTo("2.76TeV") == 0) {
        return Form("pp #rightarrow %s (#rightarrow #gamma#gamma #rightarrow e^{+}e^{-}e^{+}e^{-}) + X @ 2.76 TeV ",textProcessOpt.Data());
    } else if( fEnergyFlagOpt.CompareTo("PbPb_2.76TeV") == 0) {
        return  Form("PbPb #rightarrow %s (#rightarrow #gamma#gamma #rightarrow e^{+}e^{-}e^{+}e^{-}) + X @ 2.76 TeV ",textProcessOpt.Data());
    } else if( fEnergyFlagOpt.CompareTo("pPb_5.023TeV") == 0) {
        return  Form("pPb #rightarrow %s (#rightarrow #gamma#gamma #rightarrow e^{+}e^{-}e^{+}e^{-}) + X @ 5.02 TeV ",textProcessOpt.Data());
    } else {
        cout << "No correct collision system specification, has been given" << endl;
        return "";
    }
}

//************************************************************************************
//return detailed labels Pi+pi-gamma channel for different processes depending on energy 
//************************************************************************************
TString ReturnFullTextMesonPiPlPiMiGamma(   TString fEnergyFlagOpt,
                                            TString textProcessOpt){ 

    if(fEnergyFlagOpt.CompareTo("7TeV") == 0){
        return Form("pp #rightarrow %s (#rightarrow #pi^{+}#pi^{-}#gamma) + X @ 7 TeV ",textProcessOpt.Data());
    } else if( fEnergyFlagOpt.CompareTo("8TeV") == 0) {
        return  Form("pp #rightarrow %s (#rightarrow #pi^{+}#pi^{-}#gamma) + X @ 8 TeV ",textProcessOpt.Data());
    } else if( fEnergyFlagOpt.CompareTo("900GeV") == 0) {
        return  Form("pp #rightarrow %s (#rightarrow #pi^{+}#pi^{-}#gamma) + X @ 900 GeV ",textProcessOpt.Data());
    } else if( fEnergyFlagOpt.CompareTo("2.76TeV") == 0) {
        return Form("pp #rightarrow %s (#rightarrow #pi^{+}#pi^{-}#gamma) + X @ 2.76 TeV ",textProcessOpt.Data());
    } else if( fEnergyFlagOpt.CompareTo("PbPb_2.76TeV") == 0) {
        return  Form("PbPb #rightarrow %s (#rightarrow #pi^{+}#pi^{-}#gamma) + X @ 2.76 TeV ",textProcessOpt.Data());
    } else if( fEnergyFlagOpt.CompareTo("pPb_5.023TeV") == 0) {
        return  Form("pPb #rightarrow %s (#rightarrow #pi^{+}#pi^{-}#gamma) + X @ 5.02 TeV ",textProcessOpt.Data());
    } else {
        cout << "No correct collision system specification, has been given" << endl;
        return "";
    }
}

//************************************************************************************
//***************** return proper labeling for collision system **********************
//************************************************************************************
TString ReturnFullCollisionsSystem( TString fEnergyFlagOpt){ 
    if(fEnergyFlagOpt.CompareTo("7TeV") == 0){
        return  "pp, #sqrt{#it{s}} = 7 TeV";
    } else if( fEnergyFlagOpt.CompareTo("8TeV") == 0) {
        return  "pp, #sqrt{#it{s}} = 8 TeV";
    } else if( fEnergyFlagOpt.CompareTo("13TeV") == 0) {
        return  "pp, #sqrt{#it{s}} = 13TeV";
    } else if( fEnergyFlagOpt.CompareTo("900GeV") == 0) {
        return  "pp, #sqrt{#it{s}} = 900 GeV";
    } else if( fEnergyFlagOpt.CompareTo("2.76TeV") == 0) {
        return  "pp, #sqrt{#it{s}} = 2.76 TeV";
    } else if( (fEnergyFlagOpt.CompareTo("PbPb_2.76TeV") == 0) || (fEnergyFlagOpt.CompareTo("HI") == 0) ) {
        return "Pb-Pb, #sqrt{#it{s}_{_{NN}}} = 2.76 TeV";
    } else if( fEnergyFlagOpt.CompareTo("pPb_5.023TeV") == 0) {
        return "p-Pb, #sqrt{#it{s}_{_{NN}}} = 5.02 TeV";
    } else {
        cout << "No correct collision system specification, has been given" << endl;
        return "";
    }
}

//************************************************************************************
//***************** return proper cms-energy for collision system ********************
//************************************************************************************
Double_t ReturnCollisionEnergy( TString fEnergyFlagOpt){ 
    if(fEnergyFlagOpt.CompareTo("7TeV") == 0){
        return  7000;
    } else if( fEnergyFlagOpt.CompareTo("8TeV") == 0) {
        return 8000; 	
    } else if( fEnergyFlagOpt.CompareTo("13TeV") == 0) {
        return 13000; 
    } else if( fEnergyFlagOpt.CompareTo("2.76TeV") == 0) {
        return 2760; 
    } else if( fEnergyFlagOpt.CompareTo("900GeV") == 0) {
        return 900;
    } else if( (fEnergyFlagOpt.CompareTo("PbPb_2.76TeV") == 0) || (fEnergyFlagOpt.CompareTo("HI") == 0) ) {
        return 2760;
    } else if( fEnergyFlagOpt.CompareTo("pPb_5.023TeV") == 0) {
        return 5023;
    } else {
        cout << "No correct collision system energy specification, has been given" << endl;
        return 1;     
    }
}


//************************************************************************************
//********************* EventCuts definition *****************************************
//************************************************************************************
Int_t GetEventSystemCutPosition ()                      {return 0;}
Int_t GetEventCentralityMinCutPosition ()               {return 1;}
Int_t GetEventCentralityMaxCutPosition ()               {return 2;}
Int_t GetEventSelectSpecialTriggerCutPosition ()        {return 3;}
Int_t GetEventSelectSpecialSubTriggerCutPosition ()     {return 4;}
Int_t GetEventRemovePileUpCutPosition ()                {return 5;}
Int_t GetEventRejectExtraSignalsCutPosition ()          {return 6;}
Int_t GetEventVertexCutPosition ()                      {return 7;}

//************************************************************************************
//************************** Conversion photonCut definition *************************
//************************************************************************************
Int_t GetPhotonV0FinderCutPosition (TString gammaCutNumber){
    if(gammaCutNumber.Length() > 23) return 0;
    else return 0;
}
Int_t GetPhotonEtaCutPosition (TString gammaCutNumber){
    if(gammaCutNumber.Length() > 23) return 1;
    else return 1;
}
Int_t GetPhotonMinRCutPosition (TString gammaCutNumber){
    if(gammaCutNumber.Length() > 23) return 2;
    else return 2;
}
Int_t GetPhotonEtaForPhiCutPosition (TString gammaCutNumber){	 
    if(gammaCutNumber.Length() > 23) return 3;
    else return -1;
}
Int_t GetPhotonMinPhiCutPosition (TString gammaCutNumber){	
    if(gammaCutNumber.Length() > 23) return 4;
    else return -1;
}
Int_t GetPhotonMaxPhiCutPosition (TString gammaCutNumber){	 
    if(gammaCutNumber.Length() > 23) return 5;
    else return -1;
}  
Int_t GetPhotonSinglePtCutPosition (TString gammaCutNumber){	 
    if(gammaCutNumber.Length() > 23) return 6;
    else return 3;
}
Int_t GetPhotonClsTPCCutPosition (TString gammaCutNumber){
    if(gammaCutNumber.Length() > 23) return 7;
    else return 4;
}
Int_t GetPhotonEDedxSigmaCutPosition (TString gammaCutNumber){	
    if(gammaCutNumber.Length() > 23) return 8;
    else return 5;
}
Int_t GetPhotonPiDedxSigmaCutPosition (TString gammaCutNumber){
    if(gammaCutNumber.Length() > 23) return 9;
    else return 6;
} 
Int_t GetPhotonPiMomDedxSigmaCutPosition (TString gammaCutNumber){	
    if(gammaCutNumber.Length() > 23) return 10;
    else return 7;
} 
Int_t GetPhotonPiMaxMomDedxSigmaCutPosition (TString gammaCutNumber){		
    if(gammaCutNumber.Length() > 23) return 11;
    else return 8;
} 
Int_t GetPhotonLowPRejectionSigmaCutPosition (TString gammaCutNumber){		
    if(gammaCutNumber.Length() > 23) return 12;
    else return 9;
} 
Int_t GetPhotonTOFelectronPIDCutPosition (TString gammaCutNumber){		
    if(gammaCutNumber.Length() > 23) return 13;
    else return 10;
} 
Int_t GetPhotonQtMaxCutPosition (TString gammaCutNumber){		
    if(gammaCutNumber.Length() == 26) return 16;
    else if(gammaCutNumber.Length() == 24) return 14;
    else return 11;
} 
Int_t GetPhotonChi2GammaCutPosition (TString gammaCutNumber){		
    if(gammaCutNumber.Length() == 26) return 17;
    else if(gammaCutNumber.Length() == 24) return 15;
    else return 12;
} 
Int_t GetPhotonPsiPairCutPosition (TString gammaCutNumber){	
    if(gammaCutNumber.Length() == 26) return 18;
    else if(gammaCutNumber.Length() == 24) return 16;
    else return 13;
} 
Int_t GetPhotonDoPhotonAsymmetryCutPosition (TString gammaCutNumber){	
    if(gammaCutNumber.Length() == 26) return 19;
    else if(gammaCutNumber.Length() == 24) return 17;
    else return 14;
} 
Int_t GetPhotonCosinePointingAngleCutPosition (TString gammaCutNumber){	
    if(gammaCutNumber.Length() == 26) return 20;
    else if(gammaCutNumber.Length() == 24) return 18;
    else return 15;
}
Int_t GetPhotonSharedElectronCutPosition (TString gammaCutNumber){	
    if(gammaCutNumber.Length() == 26) return 21;
    else if(gammaCutNumber.Length() == 24) return 19;
    else return 16;
} 
Int_t GetPhotonRejectToCloseV0sCutPosition (TString gammaCutNumber){		
    if(gammaCutNumber.Length() == 26) return 22;
    else if(gammaCutNumber.Length() == 24) return 20;
    else return 17;
}
Int_t GetPhotonDcaRPrimVtxCutPosition (TString gammaCutNumber){	
    if(gammaCutNumber.Length() == 26) return 23;
    else if(gammaCutNumber.Length() == 24) return 21;
    else return 18;
}
Int_t GetPhotonDcaZPrimVtxCutPosition (TString gammaCutNumber){	
    if(gammaCutNumber.Length() == 26) return 24;
    else if(gammaCutNumber.Length() == 24) return 22;
    else return 19;
}
Int_t GetPhotonEventPlaneCutPosition (TString gammaCutNumber){		
    if(gammaCutNumber.Length() == 26) return 25;
    else if(gammaCutNumber.Length() == 24) return 23;
    else return 20;
}

//************************************************************************************
//************************** ClusterCuts definition **********************************
//************************************************************************************
Int_t GetClusterTypeCutPosition ( TString clusterCutNumber){
    if (clusterCutNumber.Length() == 17) 	return 0;
    else return 0;
}
Int_t GetClusterEtaMinCutPosition (TString clusterCutNumber){
    if (clusterCutNumber.Length() == 17) 	return 1;
    else return 1;
}
Int_t GetClusterEtaMaxCutPosition (TString clusterCutNumber){
    if (clusterCutNumber.Length() == 17) 	return 2;
    else return 2;
}
Int_t GetClusterPhiMinCutPosition (TString clusterCutNumber){
    if (clusterCutNumber.Length() == 17) 	return 3;
    else return 3;
}
Int_t GetClusterPhiMaxCutPosition (TString clusterCutNumber){
    if (clusterCutNumber.Length() == 17) 	return 4;
    else return 4;
}
Int_t GetClusterNonLinearityCutPosition (TString clusterCutNumber){
    if (clusterCutNumber.Length() == 17) 	return -1;
    else return 5;
}
Int_t GetClusterDistanceToBadChannelCutPosition (TString clusterCutNumber){
    if (clusterCutNumber.Length() == 17) 	return 5;
    else return 7;
}
Int_t GetClusterTimingCutPosition (TString clusterCutNumber){
    if (clusterCutNumber.Length() == 17) 	return 6;
    else return 8;
}
Int_t GetClusterTrackMatchingCutPosition (TString clusterCutNumber){
    if (clusterCutNumber.Length() == 17) 	return 7;
    else return 9;
}
Int_t GetClusterExoticCellCutPosition (TString clusterCutNumber){
    if (clusterCutNumber.Length() == 17) 	return 8;
    else return 10;
}
Int_t GetClusterMinEnergyCutPosition (TString clusterCutNumber){
    if (clusterCutNumber.Length() == 17) 	return 9;
    else return 11;
}
Int_t GetClusterMinNCellsCutPosition (TString clusterCutNumber){
    if (clusterCutNumber.Length() == 17) 	return 10;
    else return 12;
}
Int_t GetClusterMinM02CutPosition (TString clusterCutNumber){
    if (clusterCutNumber.Length() == 17) 	return 11;
    else return 13;
}
Int_t GetClusterMaxM02CutPosition (TString clusterCutNumber){
    if (clusterCutNumber.Length() == 17) 	return 12;
    else return 14;
}
Int_t GetClusterMinM20CutPosition (TString clusterCutNumber){
    if (clusterCutNumber.Length() == 17) 	return 13;
    else return 15;
}
Int_t GetClusterMaxM20CutPosition (TString clusterCutNumber){
    if (clusterCutNumber.Length() == 17) 	return 14;
    else return 16;
}
Int_t GetClusterMaximumDispersionCutPosition (TString clusterCutNumber){
    if (clusterCutNumber.Length() == 17) 	return 15;
    else return 17;
}
Int_t GetClusterNLMCutPosition (TString clusterCutNumber){
    if (clusterCutNumber.Length() == 17) 	return 16;
    else return 18;
}

//************************************************************************************
//***************************** MesonCuts defintion **********************************
//************************************************************************************
Int_t GetMesonKindCutPosition ()                        {return 0;}
Int_t GetMesonBGSchemeCutPosition ()                    {return 1;}
Int_t GetMesonNumberOfBGEventsCutPosition ()            {return 2;}
Int_t GetMesonDegreesForRotationMethodCutPosition ()    {return 3;}
Int_t GetMesonRapidityCutPosition ()                    {return 4;}
Int_t GetMesonRCutPosition ()                           {return 5;}
Int_t GetMesonAlphaCutPosition ()                       {return 6;}
Int_t GetMesonSelectionWindowCutPosition ()             {return 7;}
Int_t GetMesonSharedElectronCutPosition ()              {return 8;}
Int_t GetMesonRejectToCloseV0sCutPosition ()            {return 9;}
Int_t GetMesonUseMCPSmearingCutPosition ()              {return 10;}
Int_t GetMesonDcaGammaGammaCutPosition ()               {return 11;}
Int_t GetMesonDcaRCutPosition ()                        {return 12;}
Int_t GetMesonDcaZCutPosition ()                        {return 13;}
Int_t GetMesonOpeningAngleCutPosition ()                {return 14;}

//************************************************************************************
//*********************** date generation for plot labeling **************************
//************************************************************************************
TString ReturnDateString(){
    TDatime today;
    int iDate           = today.GetDate();
    int iYear           = iDate/10000;
    int iMonth          = (iDate%10000)/100;
    int iDay            = iDate%100;
    TString cMonth[12]  = {"Jan","Feb","Mar","Apr","May","Jun",
                           "Jul","Aug","Sep","Oct","Nov","Dec"};
    TString textDayth;
    if (iDay== 11){
        textDayth       = "th";
    } else if  (iDay== 12){
        textDayth       = "th";
    } else if  (iDay== 13){
        textDayth       = "th";
    } else if  (iDay%10 == 1){
        textDayth       = "st";
    } else if (iDay%10 == 2){
        textDayth       = "nd";
    } else if (iDay%10 == 3){
        textDayth       = "rd";
    } else {
        textDayth       = "th";
    }
    return Form("%i^{%s} %s %i",iDay, textDayth.Data(),cMonth[iMonth-1].Data(), iYear);
}

//************************************************************************************
//*********************** date generation for output names ***************************
//************************************************************************************
TString ReturnDateStringForOutput(){
    TDatime today;
    int iDate           = today.GetDate();
    int iYear           = iDate/10000;
    int iMonth          = (iDate%10000)/100;
    int iDay            = iDate%100;
    TString cMonth[12]  = {"Jan","Feb","Mar","Apr","May","Jun",
                           "Jul","Aug","Sep","Oct","Nov","Dec"};
    TString textDayth;
    if (iDay== 11){
        textDayth       = "th";
    } else if  (iDay== 12){
        textDayth       = "th";
    } else if  (iDay== 13){
        textDayth       = "th";
    } else if  (iDay%10 == 1){
        textDayth       = "st";
    } else if (iDay%10 == 2){
        textDayth       = "nd";
    } else if (iDay%10 == 3){
        textDayth       = "rd";
    } else {
        textDayth       = "th";
    }
    return Form("%i_%02d_%02d",iYear, iMonth, iDay);
}

//************************************************************************************
//*********************** time generation for output names ***************************
//************************************************************************************
TString ReturnTimeStringForOutput(){
  TDatime 	today;
  int 		iHour = today.GetHour();
  int     iMinute = today.GetMinute();
  
  if(iHour < 10){
    if(iMinute < 10){
      return Form("_0%i0%i",iHour, iMinute);
    }
    else{
      return Form("_0%i%i",iHour, iMinute);
    }
  }
  else{
    if(iMinute < 10){
      return Form("_%i0%i",iHour, iMinute);
    }
    else{
      return Form("_%i%i",iHour, iMinute);
    }
  }
}

//************************************************************************************
//** Analyzes meson rapidity cut, returns labeling string + double for normalization *
//************************************************************************************
Double_t ReturnRapidityStringAndDouble( TString cutSel, 
                                        TString& rapidityRangeDummy){
    
    TString rapitdityCutNumberDummy     = cutSel(GetMesonRapidityCutPosition(),1);
    if (rapitdityCutNumberDummy.CompareTo("0") == 0){
        cout << "using rapidity of 0.9" << endl;
        rapidityRangeDummy              = "0.9";
        return 1.8;
    } else if (rapitdityCutNumberDummy.CompareTo("1") == 0){
        cout << "using rapidity of 0.8" << endl;
        rapidityRangeDummy              = "0.8";
        return 1.6;
    } else if (rapitdityCutNumberDummy.CompareTo("2") == 0){
        cout << "using rapidity of 0.7" << endl;
        rapidityRangeDummy              = "0.7";
        return 1.4;
    } else if (rapitdityCutNumberDummy.CompareTo("3") == 0){
        cout << "using rapidity of 0.6" << endl;
        rapidityRangeDummy              = "0.6";
        return 1.2;
    } else if (rapitdityCutNumberDummy.CompareTo("4") == 0){
        cout << "using rapidity of 0.5" << endl;
        rapidityRangeDummy              = "0.5";
        return 1.;
    } else if (rapitdityCutNumberDummy.CompareTo("5") == 0){
        cout << "using rapidity of 0.85" << endl;
        rapidityRangeDummy              = "0.85";
        return 1.7;
    } else if (rapitdityCutNumberDummy.CompareTo("6") == 0){
        cout << "using rapidity of 0.75" << endl;
        rapidityRangeDummy              = "0.75";
        return 1.5;
    } else if (rapitdityCutNumberDummy.CompareTo("7") == 0){
        cout << "using rapidity of 0.3" << endl;
        rapidityRangeDummy              = "0.3";
        return 0.6;
    } else if (rapitdityCutNumberDummy.CompareTo("8") == 0){
        cout << "using rapidity of 0.35" << endl;
        rapidityRangeDummy              = "0.35";
        return 0.7;
    } else if (rapitdityCutNumberDummy.CompareTo("9") == 0){
        cout << "using rapidity of 0.4" << endl;
        rapidityRangeDummy              = "0.4";
        return 0.8;   
    } else {
        cout <<  " no rapidity Range selected" << endl;
        return 1.;
    }
}

//************************************************************************************
//************* Analyzes photon eta cut, returns double for normalization ************
//************************************************************************************
Double_t ReturnDeltaEta(TString gammaCutNumber){
    
    TString etaCutNumber(gammaCutNumber(GetPhotonEtaCutPosition(gammaCutNumber),1));
    if (etaCutNumber.CompareTo("0")==0){
        cout << "using eta for gammas of 0.9" << endl;
        return  1.8;
    } else if (etaCutNumber.CompareTo("1")==0){
        cout << "using eta for gammas of 0.6" << endl;
        return  1.2;
    } else if (etaCutNumber.CompareTo("2")==0){
        cout << "using eta for gammas of 1.4" << endl;
        return  2.8;
    } else if (etaCutNumber.CompareTo("3")==0){
        cout << "using eta for gammas of 0.8" << endl;
        return 1.6;
    } else if (etaCutNumber.CompareTo("4")==0){
        cout << "using eta for gammas of 0.75" << endl;
        return  1.5;
    } else if (etaCutNumber.CompareTo("5")==0){
        cout << "using eta for gammas of 0.5" << endl;
        return  1.0;
    } else if (etaCutNumber.CompareTo("6")==0){
        cout << "using eta for gammas of 5.0" << endl;
        return  10.0;
    } else if (etaCutNumber.CompareTo("7")==0){
        cout << "using eta for gammas of 0.3" << endl;
        return  0.6;
    } else if (etaCutNumber.CompareTo("8")==0){
        cout << "using eta for gammas of 0.4" << endl;
        return 0.8;
    }	
    cout << "Eta Value NOT found!!! using eta for gammas of 0.9" << endl;
    return 1.8;
}

//************************************************************************************
//******************* Analyse minimum eta cut for clusters ***************************
//************************************************************************************
Double_t AnalyseClusterMinEtaCut (Int_t etaMin){
    switch (etaMin){
        case 0: 
            return -10.;
        case 1:
            return -0.6687;
        case 2:
            return -0.5;
        case 3:
            return -2;
        case 4:
            return -0.13;
        default:
            return 0;   
    }        
}

//************************************************************************************
//******************* Analyse maximum eta cut for clusters ***************************
//************************************************************************************
Double_t AnalyseClusterMaxEtaCut (Int_t etaMin){
    switch (etaMin){
        case 0: 
            return 10.;
        case 1:
            return 0.6687;
        case 2:
            return 0.5;
        case 3:
            return 2;
        case 4:
            return 0.13;
        default:
            return 0;   
    }        
}

//************************************************************************************
//******************* Analyse minimum phi cut for clusters ***************************
//************************************************************************************
Double_t AnalyseClusterMinPhiCut (Int_t etaMin){
    switch (etaMin){
        case 0: 
            return -10000;
        case 1:
            return 1.39626;
        case 2:
            return 2.1;
        case 3:
            return 2.45;
        case 4:
            return 4.54;
        default:
            return 0;   
    }        
}

//************************************************************************************
//******************* Analyse maximum phi cut for clusters ***************************
//************************************************************************************
Double_t AnalyseClusterMaxPhiCut (Int_t etaMin){
    switch (etaMin){
        case 0: 
            return 10000;
        case 1:
            return 3.15;
        case 2:
            return 2.45;
        case 3:
            return 2.10;
        case 4:
            return 5.59;
        default:
            return 0;   
    }        
}

//************************************************************************************
//****** Analyzes photon (cluster) eta cut, returns double for normalization *********
//************************************************************************************
Double_t ReturnDeltaEtaCalo(TString caloCutNumber){
    
    TString etaMinCutNumber(caloCutNumber(GetClusterEtaMinCutPosition(caloCutNumber),1));
    TString etaMaxCutNumber(caloCutNumber(GetClusterEtaMaxCutPosition(caloCutNumber),1));
   
    Float_t minEtaCut   = AnalyseClusterMinEtaCut(etaMinCutNumber.Atoi()); 
    Float_t maxEtaCut   = AnalyseClusterMaxEtaCut(etaMaxCutNumber.Atoi());
    Float_t deltaEtaCut = TMath::Abs(minEtaCut) + TMath::Abs(maxEtaCut);
    
    return deltaEtaCut;
}

//************************************************************************************
//****** Analyzes photon (cluster) phi cut, returns double for normalization *********
//************************************************************************************
Double_t ReturnDeltaPhiCalo(TString caloCutNumber){
    
    TString phiMinCutNumber(caloCutNumber(GetClusterPhiMinCutPosition(caloCutNumber),1));
    TString phiMaxCutNumber(caloCutNumber(GetClusterPhiMaxCutPosition(caloCutNumber),1));
   
    Float_t minPhiCut   = AnalyseClusterMinPhiCut(phiMinCutNumber.Atoi()); 
    Float_t maxPhiCut   = AnalyseClusterMaxPhiCut(phiMaxCutNumber.Atoi());
    Float_t deltaPhiCut = maxPhiCut - minPhiCut; 
    if ( minPhiCut == -10000 && maxPhiCut == 10000 )
        deltaPhiCut     = 2*TMath::Pi();
    
    return deltaPhiCut;
}


//************************************************************************************
//************* Analyzes meson BG mult cut, returns float for expected BG number *****
//************************************************************************************
Float_t ReturnBackgroundMult(TString cutSel){
    TString fBackgroundMultCutNumberDummy = cutSel(GetMesonNumberOfBGEventsCutPosition(),1);
    if (fBackgroundMultCutNumberDummy.CompareTo("0") == 0){
        cout << "using number of events for BG 5" << endl;
        return 5;
    } else if (fBackgroundMultCutNumberDummy.CompareTo("1") == 0){
        cout << "using number of events for BG 10" << endl;
        return 10;
    } else if (fBackgroundMultCutNumberDummy.CompareTo("2") == 0){
        cout << "using number of events for BG 15" << endl;
        return 15;
    } else if (fBackgroundMultCutNumberDummy.CompareTo("3") == 0){
        cout << "using number of events for BG 20" << endl;
        return 20;
    } else if (fBackgroundMultCutNumberDummy.CompareTo("4") == 0){
        cout << "using number of events for BG 2" << endl;
        return 2;
    } else if (fBackgroundMultCutNumberDummy.CompareTo("5") == 0){
        cout << "using number of events for BG 50" << endl;
        return 50;
    } else if (fBackgroundMultCutNumberDummy.CompareTo("6") == 0){
        cout << "using number of events for BG 80" << endl;
        return 80;
    } else if (fBackgroundMultCutNumberDummy.CompareTo("7") == 0){
        cout << "using number of events for BG 100" << endl;
        return 100;
    }
    return 0;
}   

//************************************************************************************
//****** Analyzes the eventCutnumber and returns the NColl for the respective cent ***
//************************************************************************************
Double_t GetNCollFromCutNumber (TString cutNumber){
    TString systemCutNumber = cutNumber(GetEventSystemCutPosition(),1);
    TString centralityCutNumber = cutNumber(GetEventCentralityMinCutPosition(),2);
    if (systemCutNumber.CompareTo("1") == 0 || systemCutNumber.CompareTo("2") == 0 || systemCutNumber.CompareTo("5") == 0){
        if (centralityCutNumber.CompareTo("02") == 0){ //0-20%
            return ncoll0020;
        } else if (centralityCutNumber.CompareTo("01") == 0){ //0-10%
            return ncoll0010;
        } else if (centralityCutNumber.CompareTo("12") == 0){ //10-20%
            return ncoll1020;
        } else if (centralityCutNumber.CompareTo("24") == 0){ //20-40%
            return ncoll2040;
        } else if (centralityCutNumber.CompareTo("25") == 0){ //20-50%
            return ncoll2050;
        } else if (centralityCutNumber.CompareTo("46") == 0){ //40-60%
            return ncoll4060;
        } else if (centralityCutNumber.CompareTo("48") == 0){ //40-80%
            return ncoll4080;		
        } else if (centralityCutNumber.CompareTo("45") == 0){ //40-50%
            return ncoll4050;   
        } else if (centralityCutNumber.CompareTo("56") == 0){ //50-60%
            return ncoll5060;   
        } else if (centralityCutNumber.CompareTo("68") == 0){ //60-80%
            return ncoll6080;
        } else if (centralityCutNumber.CompareTo("67") == 0){ //60-70%
            return ncoll6070;   
        } else if (centralityCutNumber.CompareTo("78") == 0){ //70-80%
            return ncoll7080;   
        } else if (centralityCutNumber.CompareTo("89") == 0){ //80-90%
            return ncoll8090;   
        } else if (centralityCutNumber.CompareTo("48") == 0){ //40-80%
            return 77.1;   
        } else if (centralityCutNumber.CompareTo("04") == 0){ //0-40%
            return 740.;   
        }      
    } else if (systemCutNumber.CompareTo("3") == 0 || systemCutNumber.CompareTo("6") == 0){   
        if (centralityCutNumber.CompareTo("01") == 0){ //0-5%
            return ncoll0005;
        } else if (centralityCutNumber.CompareTo("12") == 0){ //0-5%
            return ncoll0510;
        } 
    } else if (systemCutNumber.CompareTo("4") == 0 || systemCutNumber.CompareTo("7") == 0){
        if (centralityCutNumber.CompareTo("69") == 0){ //75-90%
            return ncoll7590;
        }
    } else return 1.;
    return 1.;
}

//************************************************************************************
//** Analyzes the eventCutnumber and returns the NColl error for the respective cent *
//************************************************************************************
Double_t GetNCollErrFromCutNumber (TString cutNumber){
    TString systemCutNumber = cutNumber(GetEventSystemCutPosition(),1);
    TString centralityCutNumber = cutNumber(GetEventCentralityMinCutPosition(),2);
    if (systemCutNumber.CompareTo("1") == 0 || systemCutNumber.CompareTo("2") == 0 || systemCutNumber.CompareTo("5") == 0){
        if (centralityCutNumber.CompareTo("02") == 0){ //0-20%
            return nCollErr0020;
        } else if (centralityCutNumber.CompareTo("01") == 0){ //0-10%
            return nCollErr0010;
        } else if (centralityCutNumber.CompareTo("12") == 0){ //10-20%
            return nCollErr1020;
        } else if (centralityCutNumber.CompareTo("24") == 0){ //20-40%
            return nCollErr2040;
        } else if (centralityCutNumber.CompareTo("25") == 0){ //20-40%
            return nCollErr2050;
        } else if (centralityCutNumber.CompareTo("46") == 0){ //40-60%
            return nCollErr4060;
        } else if (centralityCutNumber.CompareTo("48") == 0){ //40-80%
            return nCollErr4080;
        } else if (centralityCutNumber.CompareTo("45") == 0){ //40-50%
            return nCollErr4050;   
        } else if (centralityCutNumber.CompareTo("56") == 0){ //50-60%
            return nCollErr5060;   
        } else if (centralityCutNumber.CompareTo("68") == 0){ //60-80%
            return nCollErr6080;
        } else if (centralityCutNumber.CompareTo("67") == 0){ //60-70%
            return nCollErr6070;   
        } else if (centralityCutNumber.CompareTo("78") == 0){ //70-80%
            return nCollErr7080;   
        } else if (centralityCutNumber.CompareTo("89") == 0){ //80-90%
            return nCollErr8090;   
        }      
    } else if (systemCutNumber.CompareTo("3") == 0 || systemCutNumber.CompareTo("6") == 0){   
        if (centralityCutNumber.CompareTo("01") == 0){ //0-5%
            return nCollErr0005;
        } else if (centralityCutNumber.CompareTo("12") == 0){ //0-5%
            return nCollErr0510;
        } 
    } else if (systemCutNumber.CompareTo("4") == 0 || systemCutNumber.CompareTo("7") == 0){
        if (centralityCutNumber.CompareTo("69") == 0){ //75-90%
            return nCollErr7590;
        }
    } else return 1.;
    return 1;
}


//************************************************************************************
//***** Analyzes the name of the cent and return the NColl for the respective cent ***
//************************************************************************************
Double_t GetNCollFromName (TString name){
    if (name.CompareTo("0020") == 0){ //0-20%
        return ncoll0020;
    } else if (name.CompareTo("0005") == 0){ //0-5%
        return ncoll0005;
    } else if (name.CompareTo("0510") == 0){ //0-5%
        return ncoll0510;
    } else if (name.CompareTo("0010") == 0){ //0-10%
        return ncoll0010;
    } else if (name.CompareTo("1020") == 0){ //10-20%
        return ncoll1020;
    } else if (name.CompareTo("2040") == 0){ //20-40%
        return ncoll2040;
    } else if (name.CompareTo("2050") == 0){ //20-50%
        return ncoll2050;
    } else if (name.CompareTo("4060") == 0){ //40-60%
        return ncoll4060;
    } else if (name.CompareTo("4080") == 0){ //40-80%
        return ncoll4080;	
    } else if (name.CompareTo("4050") == 0){ //40-50%
        return ncoll4050;
    } else if (name.CompareTo("5060") == 0){ //40-60%
        return ncoll5060;   
    } else if (name.CompareTo("6080") == 0){ //60-80%
        return ncoll6080;
    } else if (name.CompareTo("6070") == 0){ //60-80%
        return ncoll6070;   
    } else if (name.CompareTo("7080") == 0){ //60-80%
        return ncoll7080;      
    } else if (name.CompareTo("8090") == 0){ //60-80%
        return ncoll8090;
    } else if (name.CompareTo("7590") == 0){ //60-80%
        return ncoll7590;            
    } 
    else return 1.;
}

//************************************************************************************
//* Analyzes the name of the cent and return the NColl error for the respective cent *
//************************************************************************************
Double_t GetNCollErrFromName (TString name){
    if (name.CompareTo("0020") == 0){ //0-20%
        return nCollErr0020;
    } else if (name.CompareTo("0005") == 0){ //0-5%
        return nCollErr0005;
    } else if (name.CompareTo("0510") == 0){ //0-5%
        return nCollErr0510;
    } else if (name.CompareTo("0010") == 0){ //0-10%
        return nCollErr0010;
    } else if (name.CompareTo("1020") == 0){ //10-20%
        return nCollErr1020;
    } else if (name.CompareTo("2040") == 0){ //20-40%
        return nCollErr2040;
    } else if (name.CompareTo("2050") == 0){ //20-40%
        return nCollErr2050;
    } else if (name.CompareTo("4060") == 0){ //40-60%
        return nCollErr4060;
    } else if (name.CompareTo("4080") == 0){ //40-80%
        return nCollErr4080;	
    } else if (name.CompareTo("4050") == 0){ //40-50%
        return nCollErr4050;
    } else if (name.CompareTo("5060") == 0){ //50-60%
        return nCollErr5060;   
    } else if (name.CompareTo("6080") == 0){ //60-80%
        return nCollErr6080;
    } else if (name.CompareTo("6070") == 0){ //60-70%
        return nCollErr6070;
    } else if (name.CompareTo("7080") == 0){ //70-80%
        return nCollErr7080;   
    } else if (name.CompareTo("8090") == 0){ //80-90%
        return nCollErr8090;      
    } else if (name.CompareTo("7590") == 0){ //75-90%
        return nCollErr7590;      
    } 
    else return 1.;
    return 1.;
}


//************************************************************************************
//******* Analyzes the eventCutnumber and returns the secondary scaling factor for ***
//******* K0s for the respective cent                                              ***
//************************************************************************************
Double_t GetScalingFactorSecCorrection(TString cutNumber){

    TString systemCutNumber = cutNumber(GetEventSystemCutPosition(),1);
    TString centralityCutNumber = cutNumber(GetEventCentralityMinCutPosition(),2);
    if (systemCutNumber.CompareTo("1") == 0 || systemCutNumber.CompareTo("2") == 0 || systemCutNumber.CompareTo("5") == 0){
        if (centralityCutNumber.CompareTo("02") == 0){ //0-20%
            return 1./0.302 -1.; 
        } else if (centralityCutNumber.CompareTo("01") == 0){ //0-10%
            return 1./0.303 -1.;
        } else if (centralityCutNumber.CompareTo("12") == 0){ //10-20%
            return 1./0.2989 -1.;
        } else if (centralityCutNumber.CompareTo("24") == 0){ //20-40%
            return 1./0.308 -1.;
        } else if (centralityCutNumber.CompareTo("25") == 0){ //20-50%
            return 1./0.308 -1.;
        } else if (centralityCutNumber.CompareTo("46") == 0){ //40-60%
            return 1./0.344 -1.;
        } else if (centralityCutNumber.CompareTo("45") == 0){ //40-50%
            return 1./0.342 -1.;   
        } else if (centralityCutNumber.CompareTo("56") == 0){ //50-60%
            return 1./0.346 -1.;   
        } else if (centralityCutNumber.CompareTo("68") == 0){ //60-80%
            return 1./0.3979 -1.;
        } else if (centralityCutNumber.CompareTo("67") == 0){ //60-70%
            return 1./0.39 -1.; 
        } else if (centralityCutNumber.CompareTo("78") == 0){ //70-80%
            return 1./0.41 -1.;
        } else if (centralityCutNumber.CompareTo("89") == 0){ //80-90%
            return 1./0.42 -1.;
        } else if (centralityCutNumber.CompareTo("48") == 0){ //40-80%
            return 1./0.36 -1.;
        } else if (centralityCutNumber.CompareTo("04") == 0){ //0-40%
            return 1./0.303 -1.;
        } else if (centralityCutNumber.CompareTo("08") == 0){ //0-80%
            return 1./0.3 -1.;
        }
    } else if (systemCutNumber.CompareTo("3") == 0 || systemCutNumber.CompareTo("6") == 0){   
        if (centralityCutNumber.CompareTo("01") == 0){ //0-5%
            return 1./0.306  -1.;
        } else if (centralityCutNumber.CompareTo("12") == 0){ //0-5%
            return 1./0.301  -1.;
        } 
    } else if (systemCutNumber.CompareTo("4") == 0 || systemCutNumber.CompareTo("7") == 0){
        if (centralityCutNumber.CompareTo("69") == 0){ //75-90%
        return 1./0.42 -1.;
        }
    } else if (systemCutNumber.CompareTo("8") == 0 || systemCutNumber.CompareTo("9") == 0){
        if (centralityCutNumber.CompareTo("00") == 0){ //MB
            return 1./0.62 -1.; //value for HIJING MC, for DPMJet: 1./0.52 -1.;
        } else if (centralityCutNumber.CompareTo("02") == 0){ //00-20%
            return 1./0.24 -1.;
        } else if (centralityCutNumber.CompareTo("24") == 0){ //20-40%
            return 1./0.39 -1.;
        } else if (centralityCutNumber.CompareTo("46") == 0){ //40-60%
            return 1./0.61 -1.;
        } else if (centralityCutNumber.CompareTo("68") == 0){ //60-80%
            return 1./1.11 -1.;
        } else if (centralityCutNumber.CompareTo("60") == 0){ //60-100%
            return 1./1.62 -1.; 
        } else if (centralityCutNumber.CompareTo("80") == 0){ //80-100%
            return 1./3.00 -1.;
        } 
    } else return 0.;
    return 0.;
}

//************************************************************************************
//******* Analyzes the eventCutnumber and returns the secondary scaling factor for ***
//******* K0s for the respective collision system                                  ***
//************************************************************************************
Double_t ReturnCorrectK0ScalingFactor(TString fEnergyFlagOpt, TString cutNr){
    if(fEnergyFlagOpt.CompareTo("7TeV") == 0){
        return 1./0.75 -1.;
    } else if( fEnergyFlagOpt.CompareTo("8TeV") == 0) {
        return  1./0.75 -1.;	
    } else if( fEnergyFlagOpt.CompareTo("13TeV") == 0) {
        cout << "Caution: no correct K0 Scaling factor for 13TeV available yet" << endl;
        return  1./0.75 -1.;
    } else if( fEnergyFlagOpt.CompareTo("2.76TeV") == 0) {
        return  1./0.685 -1.;
    } else if( fEnergyFlagOpt.CompareTo("900GeV") == 0) {
        return 1./0.6 -1.;
    } else if( (fEnergyFlagOpt.CompareTo("PbPb_2.76TeV") == 0) || (fEnergyFlagOpt.CompareTo("HI") == 0) ) {
        return GetScalingFactorSecCorrection(cutNr.Data());
    } else if( fEnergyFlagOpt.CompareTo("pPb_5.023TeV") == 0) {
        // return 0.;
        return  GetScalingFactorSecCorrection(cutNr.Data());
    } else {
        cout << "No correct collision system specification, has been given" << endl;
        return 0.;     
    }
}

//************************************************************************************
//***** Analyzes the name of the cent and return the TAA for the respective cent *****
//************************************************************************************
Double_t GetTAAFromName (TString name){
    if (name.CompareTo("0020") == 0){ //0-20%
        return tAA0020;
    } else if (name.CompareTo("0005") == 0){ //0-5%
        return tAA0005;
    } else if (name.CompareTo("0510") == 0){ //0-5%
        return tAA0510;
    } else if (name.CompareTo("0010") == 0){ //0-10%
        return tAA0010;
    } else if (name.CompareTo("1020") == 0){ //10-20%
        return tAA1020;
    } else if (name.CompareTo("2040") == 0){ //20-40%
        return tAA2040;
    } else if (name.CompareTo("4060") == 0){ //40-60%
        return tAA4060;
    } else if (name.CompareTo("6080") == 0){ //60-80%
        return tAA6080;
    } 
    else return 1.;
}

//************************************************************************************
//** Analyzes the name of the cent and return the TAA error for the respective cent **
//************************************************************************************
Double_t GetTAAErrFromName (TString name){
    if (name.CompareTo("0020") == 0){ //0-20%
        return tAAErr0020;
    } else if (name.CompareTo("0005") == 0){ //0-5%
        return tAAErr0005;
    } else if (name.CompareTo("0510") == 0){ //5-10%
        return tAAErr0510;
    } else if (name.CompareTo("0010") == 0){ //0-10%
        return tAAErr0010;
    } else if (name.CompareTo("1020") == 0){ //10-20%
        return tAAErr1020;
    } else if (name.CompareTo("2040") == 0){ //20-40%
        return tAAErr2040;
    } else if (name.CompareTo("4060") == 0){ //40-60%
        return tAAErr4060;
    } else if (name.CompareTo("6080") == 0){ //60-80%
        return tAAErr6080;
    } 
    else return 1.;
    return 1.;
}

//************************************************************************************
//** Analyzes the eventCutNumber and returns the staring parameters for the fit ******
//************************************************************************************
void ReturnParameterSetFittingPbPb(TString cutsel, Double_t* parameters){
    TString centralityCutNumber = cutsel(GetEventSystemCutPosition(),3);
    cout << "bla here" << endl;
    if ( centralityCutNumber.CompareTo("102") == 0 || 
         centralityCutNumber.CompareTo("101") == 0 || 
         centralityCutNumber.CompareTo("112") == 0 ||
         centralityCutNumber.CompareTo("502") == 0 || 
         centralityCutNumber.CompareTo("501") == 0 || 
         centralityCutNumber.CompareTo("512") == 0 ){ //0-20%
        parameters[0] = 100.;
        parameters[1] = 1500.;
        parameters[2] = 4.;
        parameters[3] = 40.;
        parameters[4] = 0.07;
        parameters[5] = 0.7;
        parameters[6] = 3.;
        parameters[7] = 4.5;
        parameters[8] = 2.;
        parameters[9] = 18.;
    } else if (centralityCutNumber.CompareTo("301") == 0 ||
               centralityCutNumber.CompareTo("601") == 0){ //0-5%
        parameters[0] = 50.;
        parameters[1] = 1500.;
        parameters[2] = 4.;
        parameters[3] = 40.;
        parameters[4] = 0.07;
        parameters[5] = 0.7;
        parameters[6] = 3.;
        parameters[7] = 6.;
        parameters[8] = 2.;
        parameters[9] = 18.;
    } else if (centralityCutNumber.CompareTo("312") == 0 ||
               centralityCutNumber.CompareTo("612") == 0 ){ //5-10%
        parameters[0] = 50.;
        parameters[1] = 1500.;
        parameters[2] = 4.;
        parameters[3] = 40.;
        parameters[4] = 0.07;
        parameters[5] = 0.7;
        parameters[6] = 3.;
        parameters[7] = 6.;
        parameters[8] = 2.;
        parameters[9] = 18.;
    } else if (centralityCutNumber.CompareTo("124") == 0 || 
               centralityCutNumber.CompareTo("524") == 0){ //20-40%
        parameters[0] = 50.;
        parameters[1] = 900.;
        parameters[2] = 3.;
        parameters[3] = 30.;
        parameters[4] = 0.06;
        parameters[5] = 0.5;
        parameters[6] = 3.;
        parameters[7] = 4.;
        parameters[8] = 2.;
        parameters[9] = 18.;
    } else if (centralityCutNumber.CompareTo("525") == 0){ //20-50%
        parameters[0] = 50.;
        parameters[1] = 900.;
        parameters[2] = 3.;
        parameters[3] = 30.;
        parameters[4] = 0.06;
        parameters[5] = 0.5;
        parameters[6] = 3.;
        parameters[7] = 4.;
        parameters[8] = 2.;
        parameters[9] = 18.;
    } else if (centralityCutNumber.CompareTo("146") == 0 || 
               centralityCutNumber.CompareTo("546") == 0){ //40-60%
        parameters[0] = 10.;
        parameters[1] = 200.;
        parameters[2] = 4.;
        parameters[3] = 40.;
        parameters[4] = 0.07;
        parameters[5] = 0.7;
        parameters[6] = 4.2;
        parameters[7] = 4.5;
        parameters[8] = 2.;
        parameters[9] = 18.;
    } else if (centralityCutNumber.CompareTo("168") == 0 || 
               centralityCutNumber.CompareTo("568") == 0){ //60-80%
        parameters[0] = 1.;
        parameters[1] = 80.;
        parameters[2] = 2.5;
        parameters[3] = 30.;
        parameters[4] = 0.05;
        parameters[5] = 0.5;
        parameters[6] = 2.;
        parameters[7] = 3.3;
        parameters[8] = 2.;
        parameters[9] = 18.;
    } else {
        parameters[0] = 1.;
        parameters[1] = 80.;
        parameters[2] = 2.5;
        parameters[3] = 30.;
        parameters[4] = 0.05;
        parameters[5] = 0.5;
        parameters[6] = 2.;
        parameters[7] = 3.3;
        parameters[8] = 2.;
        parameters[9] = 18.;
    } 
    return;
}

//************************************************************************************
//** Analyzes the name of the cent and returns the staring parameters for the fit ****
//************************************************************************************
void ReturnParameterSetFittingPbPbFromString(TString centralityCutNumber, Double_t* parameters){
    if (centralityCutNumber.CompareTo("0020") == 0 || 
        centralityCutNumber.CompareTo("0010") == 0 || 
        centralityCutNumber.CompareTo("1020") == 0  ){ //0-20%
        parameters[0] = 100.;
        parameters[1] = 1500.;
        parameters[2] = 4.;
        parameters[3] = 40.;
        parameters[4] = 0.07;
        parameters[5] = 0.7;
        parameters[6] = 3.;
        parameters[7] = 4.5;
        parameters[8] = 2.;
        parameters[9] = 18.;
        parameters[10] = 150.;
        parameters[11] = 11.5;
        parameters[12] = 0.135;
        parameters[13] = 4.3;
        parameters[14] = 7.6;
    } else if (centralityCutNumber.CompareTo("0005") == 0){ //0-5%
        parameters[0] = 50.;
        parameters[1] = 1500.;
        parameters[2] = 4.;
        parameters[3] = 40.;
        parameters[4] = 0.07;
        parameters[5] = 0.7;
        parameters[6] = 3.;
        parameters[7] = 6.;
        parameters[8] = 2.;
        parameters[9] = 18.;
        parameters[10] = 150.;
        parameters[11] = 11.5;
        parameters[12] = 0.135;
        parameters[13] = 4.3;
        parameters[14] = 7.6;
    } else if (centralityCutNumber.CompareTo("0510") == 0 ){ //5-10%
        parameters[0] = 50.;
        parameters[1] = 1500.;
        parameters[2] = 4.;
        parameters[3] = 40.;
        parameters[4] = 0.07;
        parameters[5] = 0.7;
        parameters[6] = 3.;
        parameters[7] = 6.;
        parameters[8] = 2.;
        parameters[9] = 18.;
        parameters[10] = 150.;
        parameters[11] = 11.5;
        parameters[12] = 0.135;
        parameters[13] = 4.3;
        parameters[14] = 7.6;
    } else if (centralityCutNumber.CompareTo("2040") == 0){ //20-40%
        parameters[0] = 50.;
        parameters[1] = 900.;
        parameters[2] = 3.;
        parameters[3] = 30.;
        parameters[4] = 0.06;
        parameters[5] = 0.5;
        parameters[6] = 3.;
        parameters[7] = 4.;
        parameters[8] = 2.;
        parameters[9] = 18.;
        parameters[10] = 90.;
        parameters[11] = 10.;
        parameters[12] = 0.135;
        parameters[13] = 3.5;
        parameters[14] = 7.7;
    } else if (centralityCutNumber.CompareTo("4060") == 0){ //40-60%
        parameters[0] = 10.;
        parameters[1] = 200.;
        parameters[2] = 4.;
        parameters[3] = 40.;
        parameters[4] = 0.07;
        parameters[5] = 0.7;
        parameters[6] = 4.2;
        parameters[7] = 4.5;
        parameters[8] = 2.;
        parameters[9] = 18.;
        parameters[10] = 27.;
        parameters[11] = 10.;
        parameters[12] = 0.135;
        parameters[13] = 4.4;
        parameters[14] = 7.5;
    } else if (centralityCutNumber.CompareTo("6080") == 0){ //60-80%
        parameters[0] = 1.;
        parameters[1] = 80.;
        parameters[2] = 2.5;
        parameters[3] = 30.;
        parameters[4] = 0.05;
        parameters[5] = 0.5;
        parameters[6] = 2.;
        parameters[7] = 3.3;
        parameters[8] = 2.;
        parameters[9] = 18.;
        parameters[10] = 10.;
        parameters[11] = 10.;
        parameters[12] = 0.135;
        parameters[13] = 3.;
        parameters[14] = 7.;
        cout << "bla here" << endl;
    } 	else {
    // 		Double_t parameter6080[10] = {10.,80.,2.5,30.,0.05,.5,3.,3.3,2.,18.};
        parameters[0] = 1.;
        parameters[1] = 80.;
        parameters[2] = 2.5;
        parameters[3] = 30.;
        parameters[4] = 0.05;
        parameters[5] = 0.5;
        parameters[6] = 2.;
        parameters[7] = 3.3;
        parameters[8] = 2.;
        parameters[9] = 18.;
        parameters[10] = 10.;
        parameters[11] = 10.;
        parameters[12] = 0.135;
        parameters[13] = 3.;
        parameters[14] = 7.;
    } 
    return;
}


//************************************************************************************
//** Analyzes the eventCutNumber and returns the centrality name with % for plotting *
//************************************************************************************
TString GetCentralityString(TString cutNumber){
    TString centralityCutNumberStart    = cutNumber(GetEventCentralityMinCutPosition(),1);
    TString centralityCutNumberEnd      = cutNumber(GetEventCentralityMaxCutPosition(),1);
    TString ppCutNumber                 = cutNumber(GetEventSystemCutPosition(),1);
    if (ppCutNumber.CompareTo("0") ==0){
        return "pp"; 
    } else if ( ppCutNumber.CompareTo("1") ==0 || ppCutNumber.CompareTo("2") ==0 || ppCutNumber.CompareTo("5") ==0 || ppCutNumber.CompareTo("8") ==0 || ppCutNumber.CompareTo("9") ==0){       
        if (centralityCutNumberStart.CompareTo("0") == 0 && centralityCutNumberEnd.CompareTo("0") == 0  ){
            return "0-100%"; 
        } else if( centralityCutNumberEnd.CompareTo("0")!=0){
            return Form("%i-%i%s", centralityCutNumberStart.Atoi()*10,centralityCutNumberEnd.Atoi()*10,"%");
        } else {
            return Form("%i-%i%s", centralityCutNumberStart.Atoi()*10,centralityCutNumberEnd.Atoi()*10,"%");
        }
    } else if (ppCutNumber.CompareTo("3") ==0 || ppCutNumber.CompareTo("6") ==0){
        if (centralityCutNumberStart.CompareTo("0") == 0 && centralityCutNumberEnd.CompareTo("0") == 0  ){
            return "0-45%"; 
        } else {
            return Form("%i-%i%s", centralityCutNumberStart.Atoi()*5,centralityCutNumberEnd.Atoi()*5,"%");
        }
    } else if (ppCutNumber.CompareTo("4") ==0 || ppCutNumber.CompareTo("7") ==0){
        if (centralityCutNumberStart.CompareTo("0") == 0 && centralityCutNumberEnd.CompareTo("0") == 0  ){
            return "45-95%"; 
        } else {
            return Form("%i-%i%s",45+centralityCutNumberStart.Atoi()*5,45+centralityCutNumberEnd.Atoi()*5,"%");
        }
    } else return "";
}	

//************************************************************************************
//** Analyzes the eventCutNumber and returns the centrality name for labeling w/o % **
//************************************************************************************
TString GetCentralityStringWoPer(TString cutNumber){
    TString centralityCutNumberStart    = cutNumber(GetEventCentralityMinCutPosition(),1);
    TString centralityCutNumberEnd      = cutNumber(GetEventCentralityMaxCutPosition(),1);
    TString ppCutNumber                 = cutNumber(GetEventSystemCutPosition(),1);
    if (ppCutNumber.CompareTo("0") ==0){
        return "pp"; 
    } else if ( ppCutNumber.CompareTo("1") ==0 || ppCutNumber.CompareTo("2") ==0 || ppCutNumber.CompareTo("5") ==0 || ppCutNumber.CompareTo("8") ==0 || ppCutNumber.CompareTo("9") ==0){       
        if (centralityCutNumberStart.CompareTo("0") == 0 && centralityCutNumberEnd.CompareTo("0") == 0  ){
            return "0-100"; 
        } else {
            return Form("%i-%i", centralityCutNumberStart.Atoi()*10,centralityCutNumberEnd.Atoi()*10);
        }
    } else if (ppCutNumber.CompareTo("3") ==0 || ppCutNumber.CompareTo("6") ==0){
        if (centralityCutNumberStart.CompareTo("0") == 0 && centralityCutNumberEnd.CompareTo("0") == 0  ){
            return "0-45"; 
        } else {
            return Form("%i-%i", centralityCutNumberStart.Atoi()*5,centralityCutNumberEnd.Atoi()*5);
        }
    } else if (ppCutNumber.CompareTo("4") ==0 || ppCutNumber.CompareTo("7") ==0){
        if (centralityCutNumberStart.CompareTo("0") == 0 && centralityCutNumberEnd.CompareTo("0") == 0  ){
            return "45-95"; 
        } else {
            return Form("%i-%i",45+centralityCutNumberStart.Atoi()*5,45+centralityCutNumberEnd.Atoi()*5);
        }
    } else return ""; 
}  

//************************************************************************************
//** Analyzes the eventCutNumber and returns the centrality name as joing numbers ****
//************************************************************************************
TString GetCentralityStringOutput(TString cutNumber){
    TString centralityCutNumberStart    = cutNumber(GetEventCentralityMinCutPosition(),1);
    TString centralityCutNumberEnd      = cutNumber(GetEventCentralityMaxCutPosition(),1);
    TString ppCutNumber                 = cutNumber(GetEventSystemCutPosition(),1);
    if (ppCutNumber.CompareTo("0") ==0){
        return "pp"; 
    } else if ( ppCutNumber.CompareTo("1") ==0 || ppCutNumber.CompareTo("2") ==0 || ppCutNumber.CompareTo("5") ==0 || ppCutNumber.CompareTo("8") ==0 || ppCutNumber.CompareTo("9") ==0){       
        if (centralityCutNumberStart.CompareTo("0") == 0 && centralityCutNumberEnd.CompareTo("0") == 0  ){
            return "00100"; 
        } else {
            if (centralityCutNumberStart.CompareTo("0") == 0){
                return Form("00%i", centralityCutNumberEnd.Atoi()*10);
            } else {	
                return Form("%i%i", centralityCutNumberStart.Atoi()*10,centralityCutNumberEnd.Atoi()*10);
            }	
        }
    } else if (ppCutNumber.CompareTo("3") ==0 || ppCutNumber.CompareTo("6") ==0){
        if (centralityCutNumberStart.CompareTo("0") == 0 && centralityCutNumberEnd.CompareTo("0") == 0  ){
            return "0045"; 
        } else {
            if (centralityCutNumberStart.CompareTo("0") == 0){
                return Form("00%i", centralityCutNumberEnd.Atoi()*5);
            } else {
                return Form("%i%i", centralityCutNumberStart.Atoi()*5,centralityCutNumberEnd.Atoi()*5);
            }	
        }
    } else if (ppCutNumber.CompareTo("4") ==0 || ppCutNumber.CompareTo("7") ==0){
        if (centralityCutNumberStart.CompareTo("0") == 0 && centralityCutNumberEnd.CompareTo("0") == 0  ){
            return "4595"; 
        } else {
            return Form("%i%i",45+centralityCutNumberStart.Atoi()*5,45+centralityCutNumberEnd.Atoi()*5);
        }
    } else return ""; 
}  

//************************************************************************************
//** Analyzes the TPC cluster cut, return correct cut label **************************
//************************************************************************************
TString AnalyseTPCClusterCut(Int_t clsTPCCut){   
    switch(clsTPCCut){
        case 0: // 0
            return "min. TPC cluster: 0";
        case 1:  // 60
            return "min. TPC cluster: 60";
        case 2:  // 80
            return "min. TPC cluster: 80";
        case 3:  // 100
            return "min. TPC cluster: 100";     
        case 4:  // 95% of findable clusters
            return "TPC cluster/uncorr findable Clusters: 0.95";     
        case 5:  // 0% of findable clusters
            return "TPC cluster/corr findable Clusters: 0.";     
        case 6:  // 80% of findable clusters
            return "TPC cluster/corr findable Clusters: 0.7";     
        case 7:  // 0% of findable clusters
            return "TPC cluster/uncorr findable Clusters: 0.35";          
        case 8:
            return "TPC cluster/corr findable Clusters: 0.35";          
        case 9:
            return "TPC cluster/corr findable Clusters: 0.6";          
        default:
            return "no cluster cut defined";
    }  
}

//************************************************************************************
//** Analyzes the TPC dEdx electron cut, return correct cut label ********************
//************************************************************************************
TString AnalyseTPCdEdxCutElectronLine(Int_t ededxSigmaCut){
    switch(ededxSigmaCut){
        case 0: // -10,10
            return "-10 < #sigma_{e} < 10";
        case 1: // -5,5
            return "-5 < #sigma_{e} < 5";
        case 2: // -3,5
            return "-3 < #sigma_{e} < 5";
        case 3: // -4,5
            return "-4 < #sigma_{e} < 5";
        case 4: // -6,7
            return "-6 < #sigma_{e} < 7";
        case 5: // -4,4
            return "-4 < #sigma_{e} < 4";
        case 6: // -2.5,4
            return "-2.5 < #sigma_{e} < 4";
        case 7: // -2,3.5
            return "-2 < #sigma_{e} < 3.5";
        default:
            return "no dEdx cut defined";
    }
    return kTRUE;
}

//************************************************************************************
//** Analyzes the TPC dEdx pion cuts, return correct cut label ***********************
//************************************************************************************
TString AnalyseTPCdEdxCutPionLine(TString sPionCut){
    cout << sPionCut << endl;
    TString sPidedxSigmaCut                 = sPionCut(0,1);
    Int_t pidedxSigmaCut                    = sPidedxSigmaCut.Atoi();
    cout << "pidedxSigmaCut: " << pidedxSigmaCut << endl;
    Double_t fPIDnSigmaAbovePionLine        = 0;
    Double_t fPIDnSigmaAbovePionLineHighPt  = 0;
    switch(pidedxSigmaCut){
        case 0:  // -10
            fPIDnSigmaAbovePionLine         = -10;
            fPIDnSigmaAbovePionLineHighPt   = -10;
            break;
        case 1:   // 0
            fPIDnSigmaAbovePionLine         = 0;
            fPIDnSigmaAbovePionLineHighPt   = -10;
            break;
        case 2:  // 1
            fPIDnSigmaAbovePionLine         = 1;
            fPIDnSigmaAbovePionLineHighPt   = -10;
            break;
        case 3:  // 1
            fPIDnSigmaAbovePionLine         = 2.5;
            fPIDnSigmaAbovePionLineHighPt   = -10;
            break;
        case 4:  // 1
            fPIDnSigmaAbovePionLine         = 3.;
            fPIDnSigmaAbovePionLineHighPt   = 1.;
            break;
        case 5:  // 1
            fPIDnSigmaAbovePionLine         = 2.;
            fPIDnSigmaAbovePionLineHighPt   = -10;
            break;
        case 6:  // 1
            fPIDnSigmaAbovePionLine         = 2.;
            fPIDnSigmaAbovePionLineHighPt   = 0.5;
            break;
        case 7:  // 1
            fPIDnSigmaAbovePionLine         = 3.5;
            fPIDnSigmaAbovePionLineHighPt   = -10;
            break;
        case 8:  // 1
            fPIDnSigmaAbovePionLine         = 2.;
            fPIDnSigmaAbovePionLineHighPt   = 1.;
            break;
        case 9:
            fPIDnSigmaAbovePionLine         = 3.0; // We need a bit less tight cut on dE/dx
            fPIDnSigmaAbovePionLineHighPt   = -10;
            break;
        default:
            cout << "pion line cut unknown" << endl;
    }

    TString sPiMomdedxSigmaCut              = sPionCut(1,1);
    Int_t piMomdedxSigmaCut                 = sPiMomdedxSigmaCut.Atoi();
    Double_t fPIDMinPnSigmaAbovePionLine    = 0;
    cout << "piMomdedxSigmaCut: " << piMomdedxSigmaCut << endl;
    switch(piMomdedxSigmaCut){
        case 0:  // 0.5 GeV
            fPIDMinPnSigmaAbovePionLine     = 0.5;
            break;
        case 1:  // 1. GeV
            fPIDMinPnSigmaAbovePionLine     = 1.;
            break;
        case 2:  // 1.5 GeV
            fPIDMinPnSigmaAbovePionLine     = 1.5;
            break;
        case 3:  // 20.0 GeV
            fPIDMinPnSigmaAbovePionLine     = 20.;
            break;
        case 4:  // 50.0 GeV
            fPIDMinPnSigmaAbovePionLine     = 50.;
            break;
        case 5:  // 0.3 GeV
            fPIDMinPnSigmaAbovePionLine     = 0.3;
            break;
        case 6:  // 0.25 GeV
            fPIDMinPnSigmaAbovePionLine     = 0.25;
            break;
        case 7:  // 0.4 GeV
            fPIDMinPnSigmaAbovePionLine     = 0.4;
            break;
        case 8:  // 0.2 GeV
            fPIDMinPnSigmaAbovePionLine     = 0.2;
            break;
        default:
            cout << "pion line minimum pt cut unknown" << endl;            
    }

    TString sPiMaxMomdedxSigmaCut           = sPionCut(2,1);
    Int_t piMaxMomdedxSigmaCut              = sPiMaxMomdedxSigmaCut.Atoi();
    Double_t fPIDMaxPnSigmaAbovePionLine    = 0;
    cout << "piMaxMomdedxSigmaCut: " << piMaxMomdedxSigmaCut << endl;
    switch(piMaxMomdedxSigmaCut){
        case 0:  // 100. GeV
            fPIDMaxPnSigmaAbovePionLine     = 100.;
            break;
        case 1:  // 5. GeV
            fPIDMaxPnSigmaAbovePionLine     = 5.;
            break;
        case 2:  // 4. GeV
            fPIDMaxPnSigmaAbovePionLine     = 4.;
            break;
        case 3:  // 3.5 GeV
            fPIDMaxPnSigmaAbovePionLine     = 3.5;
            break;
        case 4:  // 3. GeV
            fPIDMaxPnSigmaAbovePionLine     = 3.;
            break;
        case 5:  // 7. GeV
            fPIDMaxPnSigmaAbovePionLine     = 7.;
            break;
        case 6:  // 2. GeV
            fPIDMaxPnSigmaAbovePionLine     = 2.;
            break;
        default:
            cout << "pion line minimum pt cut unknown" << endl;
    }
    return Form("rejected #pi for #sigma_{#pi} < %.2f (%.2f GeV/c < p_{T}_{#pi} < %.2f GeV/c), #sigma_{#pi} < %.2f (p_{T}_{#pi} > %.2f GeV/c)",fPIDnSigmaAbovePionLine, fPIDMinPnSigmaAbovePionLine,fPIDMaxPnSigmaAbovePionLine,fPIDnSigmaAbovePionLineHighPt,fPIDMaxPnSigmaAbovePionLine);
}


//************************************************************************************
//** Analyzes the TOF electron PID cut, return correct cut label *********************
//************************************************************************************
TString AnalyseTOFelectronPIDCut(Int_t TOFpidCut){   
    switch(TOFpidCut){
        case 0: // 0
            return "-100 < #sigma^{TOF}_{e} < 100";
            break;
        case 1: 
            return "-7 < #sigma^{TOF}_{e} < 7";
            break;
        case 2: 
            return "-5 < #sigma^{TOF}_{e} < 5";
            break;
        case 3: 
            return "-3 < #sigma^{TOF}_{e} < 5";
            break;
        case 4: 
            return "-2 < #sigma^{TOF}_{e} < 3";
            break;
        case 5: 
            return "-3 < #sigma^{TOF}_{e} < 3";
            break;
        default:
            return "no TOF cut defined";
    }  
}


//************************************************************************************
//********* Analyzes the chi2 gamma cut, return correct cut label ********************
//************************************************************************************
TString AnalyseChi2GammaCut(    Int_t chi2GammaCut, 
                                Int_t psiPairCut){   // Set Cut

    TString psiPairCutString            = "";
    Bool_t k2DPsiPairChi2               = kFALSE;
    switch(psiPairCut) {
        case 0:
            psiPairCutString = "|#Psi_{Pair}| < 10000";
            break;
        case 1:
            psiPairCutString = "|#Psi_{Pair}| < 0.1";
            break;
        case 2:
            psiPairCutString = "|#Psi_{Pair}| < 0.05";
            break;
        case 3:
            psiPairCutString = "|#Psi_{Pair}| < 0.035";
            break;
        case 4:
            psiPairCutString = "|#Psi_{Pair}| < 0.2";
            break;
        case 5:
            psiPairCutString = "|#Psi_{Pair}| < 0.1";
            k2DPsiPairChi2= kTRUE;
            break;
        case 6:
            psiPairCutString = "|#Psi_{Pair}| < 0.05";
            k2DPsiPairChi2= kTRUE;
            break;
        case 7:
            psiPairCutString = "|#Psi_{Pair}| < 0.035";
            k2DPsiPairChi2= kTRUE;
            break;
        case 8:
            psiPairCutString =  "|#Psi_{Pair}| < 0.2";
            k2DPsiPairChi2= kTRUE;
            break;
        case 9:
            psiPairCutString =  "|#Psi_{Pair}| < 0.5";
            break;
        default:
            psiPairCutString =  "#Psi_{Pair} cut not defined";
            break;
    }

    TString chi2CutString               = "";
    switch(chi2GammaCut){
    case 0: // 100
        chi2CutString = "#chi_{#gamma}^{2} < 100";
        break;
    case 1:  // 50
        chi2CutString = "#chi_{#gamma}^{2} < 50";
        break;
    case 2:  // 30
        chi2CutString = "#chi_{#gamma}^{2} < 30";
        break;
    case 3:
        chi2CutString = "#chi_{#gamma}^{2} < 200";
        break;
    case 4:
        chi2CutString = "#chi_{#gamma}^{2} < 500";
        break;
    case 5:
        chi2CutString = "#chi_{#gamma}^{2} < 1000000";
        break;
    case 6:
        chi2CutString = "#chi_{#gamma}^{2} < 5";
        break;
    case 7:
        chi2CutString = "#chi_{#gamma}^{2} < 10";
        break;
    case 8:
        chi2CutString = "#chi_{#gamma}^{2} < 20";
        break;
    case 9:
        chi2CutString = "#chi_{#gamma}^{2} < 15";
        break;
    default:
        chi2CutString = "#chi_{#gamma}^{2} cut unknown";
        break;
    }
    if (k2DPsiPairChi2){
        return Form("2D cut: %s, %s", chi2CutString.Data(), psiPairCutString.Data());
    } else {
        return Form("1D cut: %s", chi2CutString.Data());
    }
    return "";
}

//************************************************************************************
//********* Analyzes the qt gamma cut, return correct cut label **********************
//************************************************************************************
TString AnalyseQtMaxCut(Int_t QtMaxCut){   
		
    switch(QtMaxCut){
        case 0: //
            return "no q_{T}_{#gamma} cut applied";
        case 1:
            return "q_{T}_{#gamma} < 0.1 GeV/c";
        case 2:
            return "2D ellipse q_{T}_{#gamma} < 0.06 GeV/c, #alpha < 0.95";
        case 3:
            return "q_{T}_{#gamma} < 0.05 GeV/c";
        case 4:
            return "q_{T}_{#gamma} < 0.03 GeV/c";
        case 5:
            return "q_{T}_{#gamma} < 0.02 GeV/c";
        case 6:
            return "2D ellipse q_{T}_{#gamma} < 0.02 GeV/c, #alpha < 0.95";
        case 7:
            return "q_{T}_{#gamma} < 0.15 GeV/c";  
        case 8:
            return "2D ellipse q_{T}_{#gamma} < 0.05 GeV/c, #alpha < 0.95";
        case 9:
            return "2D ellipse q_{T}_{#gamma} < 0.03 GeV/c, #alpha < 0.95";
        default:
            return "no q_{T}_{#gamma} cut defined";
    }
}

//************************************************************************************
//********* Analyzes the single leg pt cut, return correct cut label *****************
//************************************************************************************
TString AnalyseSinglePtCut(Int_t singlePtCut){

    switch(singlePtCut){
        case 0: // 0.050 GeV
            return "p_{T}_{e^{#pm}} > 0.050 GeV/c";
        case 1:  // 0.100 GeV
            return "p_{T}_{e^{#pm}} > 0.100 GeV/c";
        case 2:  // 0.150 GeV
            return "p_{T}_{e^{#pm}} > 0.150 GeV/c";
        case 3:  // 0.200 GeV
            return "p_{T}_{e^{#pm}} > 0.200 GeV/c";
        case 4:  // 0.075 GeV
            return "p_{T}_{e^{#pm}} > 0.075 GeV/c";
        case 5:  // 0.125 GeV
            return "p_{T}_{e^{#pm}} > 0.125 GeV/c";
        case 6:  // 0.04 GeV
            return "p_{T}_{e^{#pm}} > 0.040 GeV/c";
        case 7:  // 0.0 GeV
            return "p_{T}_{e^{#pm}} > 0.0 GeV/c";
        default:
            return "p_{T}_{e^{#pm}} cut not defined";
    }
    
}

//************************************************************************************
//********* Analyzes the dca z photon cu, return correct cut label *******************
//************************************************************************************
TString AnalyseDCAZPhotonCut(Int_t dcaZPhoton){
    switch(dcaZPhoton){
        case 0:  //
            return "|dca_{Z}| < 1000 cm"; 
        case 1:  //
            return "|dca_{Z}| < 10 cm"; 
        case 2:  //
            return "|dca_{Z}| < 5 cm"; 
        case 3:  //
            return "|dca_{Z}| < 4 cm"; 
        case 4:  //
            return "|dca_{Z}| < 3 cm"; 
        case 5:  //
            return "|dca_{Z}| < 2.5 cm"; 
        case 6:  //
            return "|dca_{Z}| < 2 cm"; 
        case 7:  //
            return "|dca_{Z}| < 1.5 cm"; 
        case 8:  //
            return "|dca_{Z}| < 1 cm"; 
        case 9:  //
            return "|dca_{Z}| < 0.5 cm"; 
        default:
            return "|dca_{Z}| cut not defined";
    }
}

//************************************************************************************
//********* Analyzes the dca z gamma cut, return max dca value ***********************
//************************************************************************************
Double_t AnalyseDCAZPhotonCutValue(Int_t dcaZPhoton){
    switch(dcaZPhoton){   
        case 0:  //
            return 1000.; 
        case 1:  //
            return 10.; 
        case 2:  //
            return 5.; 
        case 3:  //
            return 4.; 
        case 4:  //
            return 3.; 
        case 5:  //
            return 2.5; 
        case 6:  //
            return 2.; 
        case 7:  //
            return 1.5; 
        case 8:  //
            return 1.; 
        case 9:  //
            return 0.5; 
        default:
            return 1000;
    }
}

//************************************************************************************
//********* Analyzes the cos(theta_point) gamma cut, return correct cut label ********
//************************************************************************************
TString AnalyseCosPointCut(Int_t cosPoint){
    switch(cosPoint){
        case 0:  //
            return "cos(#Theta_{point}) > -1"; 
        case 1:  //
            return "cos(#Theta_{point}) > 0"; 
        case 2:  //
            return "cos(#Theta_{point}) > 0.5"; 
        case 3:  //
            return "cos(#Theta_{point}) > 0.75"; 
        case 4:  //
            return "cos(#Theta_{point}) > 0.85"; 
        case 5:  //
            return "cos(#Theta_{point}) > 0.88"; 
        case 6:  //
            return "cos(#Theta_{point}) > 0.9"; 
        case 7:  //
            return "cos(#Theta_{point}) > 0.95"; 
        default:
            return "cos(#Theta_{point}) cut not defined";
    }
}

//************************************************************************************
//********* Analyzes the eta gamma/electron cut, return correct cut label ************
//************************************************************************************
TString AnalyseEtaCut(Int_t etaCut){ 

    switch(etaCut){
        case 0: // 0.9
            return "#eta_{#gamma,e^{#pm}} < 0.9";
        case 1:  // 0.6
            return "#eta_{#gamma,e^{#pm}} < 0.6";
        case 2:  // 1.4
            return "#eta_{#gamma,e^{#pm}} < 1.4";
        case 3: // 0.65
            return "#eta_{#gamma,e^{#pm}} < 0.65";
        case 4: // 0.75
            return "#eta_{#gamma,e^{#pm}} < 0.75";
        case 5: // 0.5
            return "#eta_{#gamma,e^{#pm}} < 0.5";
        case 6: // 5.
            return "#eta_{#gamma,e^{#pm}} < 5";
        case 7: // 0.1 - 0.8
            return "#eta_{#gamma,e^{#pm}} < 0.7";
        case 8: // 0.1 - 0.8
            return "#eta_{#gamma,e^{#pm}} < 0.4";
        case 9: // 10
            return "#eta_{#gamma,e^{#pm}} < 10";
        default:
            return "no #eta_{#gamma,e^{#pm}} cut defined";
    }
}

//************************************************************************************
//****** Analyzes the eta gamma/electron cut for pPb, return correct cut label *******
//************************************************************************************
TString AnalyseEtaCutpPb(Int_t etaCut){ 

    switch(etaCut){
        case 0: // 0.9
            return "#eta_{#gamma,e^{#pm}} < 0.9";
        case 1:  // 1.2
            return "#eta_{#gamma,e^{#pm}} < 0.6";
        case 2:  // 1.4
            return "#eta_{#gamma,e^{#pm}} < 1.4";
        case 3: // 0.8
            return "#eta_{#gamma,e^{#pm}} < 0.8";
        case 4: // 0.75
            return "#eta_{#gamma,e^{#pm}} < 0.75";
        case 5: // 0.9 - 1.4
            return "#eta_{#gamma,e^{#pm}} < 0.5";
        case 6: // 5.
            return "#eta_{#gamma,e^{#pm}} < 5";
        case 7: // 0.1 - 0.8
            return "#eta_{#gamma,e^{#pm}} < 0.3";
        case 8: // 0.1 - 0.8
            return "#eta_{#gamma,e^{#pm}} < 0.4";
        case 9: // 10
            return "#eta_{#gamma,e^{#pm}} < 10";
        default:
            return "no #eta_{#gamma,e^{#pm}} cut defined";
    }
}

//************************************************************************************
//**************** Analyzes the R gamma cut, return correct cut label ****************
//************************************************************************************
TString AnalyseRCut(Int_t RCut){
    // Set Cut
    switch(RCut){
        case 0:
            return "0 cm < R_{conv, #gamma} < 180 cm";
        case 1:
            return "2.8 cm < R_{conv, #gamma} < 180 cm";
        case 2:
            return "5 cm < R_{conv, #gamma} < 180 cm";
        case 3:
            return "10 cm < R_{conv, #gamma} < 70 cm";
        case 4:
            return "5 cm < R_{conv, #gamma} < 70 cm";
        case 5:
            return "10 cm < R_{conv, #gamma} < 180 cm";
        case 6:
            return "20 cm < R_{conv, #gamma} < 180 cm";
        case 7:
            return "35 cm < R_{conv, #gamma} < 180 cm"; // before 26 (19.4.2014)qdel 
        case 8:
            return "12.5 cm < R_{conv, #gamma} < 180 cm";
        case 9:
            return "7.5 cm < R_{conv, #gamma} < 180 cm";
        default:
            return "R cut not defined";
    }
}

//************************************************************************************
// Analyzes the R gamma cut, together with photonQuality cut, return correct cut label
//************************************************************************************
TString AnalyseRCutAndQuality(Int_t RCut, Int_t photonQualitCut){
    // Set Cut
    TString stringRCut = "";
    switch(RCut){
        case 0:
            stringRCut= "0 cm < R_{conv, #gamma} < 180 cm";
            break;
        case 1:
            stringRCut= "2.8 cm < R_{conv, #gamma} < 180 cm";
            break;
        case 2:
            stringRCut=  "5 cm < R_{conv, #gamma} < 180 cm";
            break;
        case 3:
            stringRCut= "10 cm < R_{conv, #gamma} < 70 cm";
            break;
        case 4:
            stringRCut= "5 cm < R_{conv, #gamma} < 70 cm";
            break;
        case 5:
            stringRCut= "10 cm < R_{conv, #gamma} < 180 cm";
            break;
        case 6:
            stringRCut= "20 cm < R_{conv, #gamma} < 180 cm";
            break;
        case 7:
            stringRCut= "35 cm < R_{conv, #gamma} < 180 cm"; // before 26 (19.4.2014)qdel 
            break;
        case 8:
            stringRCut= "12.5 cm < R_{conv, #gamma} < 180 cm";
            break;
        case 9:
            stringRCut= "7.5 cm < R_{conv, #gamma} < 180 cm";
            break;
        default:
            stringRCut= "R cut not defined";
            break;
    }
    
    TString stringPhotonQuality = "";
    switch(photonQualitCut) {
        case 0:
            stringPhotonQuality =  "photon Quality: 1,2,3";
            break;
        case 2:
            stringPhotonQuality =  "photon Quality: 1";
            break;
        case 3:
            stringPhotonQuality =  "photon Quality: 2";
            break;
        case 4:
            stringPhotonQuality =  "photon Quality: 3";
            break;
        default:
            stringPhotonQuality =  "photon Quality cut not defined";
            break;
    }
    return Form("%s, %s", stringRCut.Data(), stringPhotonQuality.Data());
    
}

//************************************************************************************
//************* Analyzes the psi_pair gamma cut, return correct cut label ************
//************************************************************************************
TString AnalysePsiPair(Int_t PsiPairCut, Int_t chi2GammaCut){

    TString psiPairCutString = "";
    Bool_t k2DPsiPairChi2 = kFALSE;
    switch(PsiPairCut) {
        case 0:
            psiPairCutString = "|#Psi_{Pair}| < 10000";
            break;
        case 1:
            psiPairCutString = "|#Psi_{Pair}| < 0.1";
            break;
        case 2:
            psiPairCutString = "|#Psi_{Pair}| < 0.05";
            break;
        case 3:
            psiPairCutString = "|#Psi_{Pair}| < 0.035";
            break;
        case 4:
            psiPairCutString = "|#Psi_{Pair}| < 0.2";
            break;
        case 5:
            psiPairCutString = "|#Psi_{Pair}| < 0.1";
            k2DPsiPairChi2= kTRUE;
            break;
        case 6:
            psiPairCutString = "|#Psi_{Pair}| < 0.05";
            k2DPsiPairChi2= kTRUE;
            break;
        case 7:
            psiPairCutString = "|#Psi_{Pair}| < 0.035";
            k2DPsiPairChi2= kTRUE;
            break;
        case 8:
            psiPairCutString =  "|#Psi_{Pair}| < 0.2";
            k2DPsiPairChi2= kTRUE;
            break;
        case 9:
            psiPairCutString =  "|#Psi_{Pair}| < 0.5";
            break;
        default:
            psiPairCutString =  "#Psi_{Pair} cut not defined";
            break;
    }
    
    TString chi2CutString = "";
    switch(chi2GammaCut){
        case 0: // 100
            chi2CutString = "#chi_{#gamma}^{2} < 100";
            break;
        case 1:  // 50
            chi2CutString = "#chi_{#gamma}^{2} < 50";
            break;
        case 2:  // 30
            chi2CutString = "#chi_{#gamma}^{2} < 30";
            break;
        case 3:
            chi2CutString = "#chi_{#gamma}^{2} < 200";
            break;
        case 4:
            chi2CutString = "#chi_{#gamma}^{2} < 500";
            break;
        case 5:
            chi2CutString = "#chi_{#gamma}^{2} < 1000000";
            break;
        case 6:
            chi2CutString = "#chi_{#gamma}^{2} < 5";
            break;
        case 7:
            chi2CutString = "#chi_{#gamma}^{2} < 10";
            break;
        case 8:
            chi2CutString = "#chi_{#gamma}^{2} < 20";
            break;
        case 9:
            chi2CutString = "#chi_{#gamma}^{2} < 15";
            break;
        default:
            chi2CutString = "#chi_{#gamma}^{2} cut unknown";
            break;
    }
    if (k2DPsiPairChi2){
        return Form("2D cut: %s, %s", psiPairCutString.Data(), chi2CutString.Data());
    } else {
        return Form("1D cut: %s", psiPairCutString.Data());
    }
    return "";
}   

//************************************************************************************
//******** Analyzes the psi_pair and R gamma cut, return correct cut label ***********
//************************************************************************************
TString AnalysePsiPairAndR(Int_t PsiPairCut, Int_t RCut ){
    TString psiPairCut = "";
    switch(PsiPairCut) {
        case 0:
            psiPairCut = "|#Psi_{Pair}| < 10000";
            break;
        case 1:
            psiPairCut = "|#Psi_{Pair}| < 0.1";
            break;
        case 2:
            psiPairCut = "|#Psi_{Pair}| < 0.05";
            break;
        case 3:
            psiPairCut = "|#Psi_{Pair}| < 0.035";
            break;
        case 4:
            psiPairCut = "|#Psi_{Pair}| < 0.15";
            break;
        case 5:
            psiPairCut = "#Psi_{Pair} < (0.1 - 0.1/1 * #Delta #Phi )";
            break;
        case 6:
            psiPairCut = "#Psi_{Pair} < (0.05 - 0.05/1 * #Delta #Phi )";
            break;
        case 7:
            psiPairCut = "#Psi_{Pair} < (0.035 - 0.035/1 * #Delta #Phi )";
            break;
        case 8:
            psiPairCut = "#Psi_{Pair} < (0.2 - 0.2/1 * #Delta #Phi )";
            break;
        case 9:
            psiPairCut = "#Psi_{Pair} < 0.5";
            break;
        default:
            psiPairCut = "#Psi_{Pair} cut not defined";
            break;
    }
    TString RCutString = "";
    switch(RCut){
        case 0:
            RCutString =  "0 cm < R_{conv, #gamma} < 180 cm";
            break;
        case 1:
            RCutString =  "2.8 cm < R_{conv, #gamma} < 180 cm";
            break;
        case 2:
            RCutString =  "5 cm < R_{conv, #gamma} < 180 cm";
            break;
        case 3:
            RCutString =  "10 cm < R_{conv, #gamma} < 70 cm";
            break;
        case 4:
            RCutString =  "5 cm < R_{conv, #gamma} < 70 cm";
            break;
        case 5:
            RCutString =  "10 cm < R_{conv, #gamma} < 180 cm";
            break;
        case 6:
            RCutString =  "20 cm < R_{conv, #gamma} < 180 cm";
            break;
        case 7:
            RCutString =  "26 cm < R_{conv, #gamma} < 180 cm";
            break;
        case 8:
            RCutString =  "35 cm < R_{conv, #gamma} < 180 cm";
            break;
        case 9:
            RCutString =  "5 cm < R_{conv, #gamma} < 35 cm";
            break;
        default:
            RCutString =  "R cut not defined";
            break;
    }
    return Form("%s, %s", psiPairCut.Data(), RCutString.Data());
}   

//************************************************************************************
//********** Analyzes the V0reader gamma cut, return correct cut label ***************
//************************************************************************************
TString AnalyseV0ReaderCut(Int_t V0ReaderCut){ 
    // Set Cut
    switch(V0ReaderCut){
        case 0:  //
            return "Onfly V0finder";
        case 1:  //
            return "Offline V0finder";
        default:
            return "V0finder cut not defined";
    }
}


//************************************************************************************
//***** Analyzes the 2D psi pair & chi2 gamma cut, return correct cut label **********
//************************************************************************************
TString AnalyseChi2PsiPair(Int_t PsiPairCut ){
    switch(PsiPairCut) {
        case 26:
            return "2D: #chi_{#gamma}^{2} = 30 |#Psi_{Pair}| < 0.05";
        case 22:
            return "1D: #chi_{#gamma}^{2} = 30 |#Psi_{Pair}| < 0.05";
        case 16:
            return "2D: #chi_{#gamma}^{2} = 50 |#Psi_{Pair}| < 0.05";
        case 86:
            return "2D: #chi_{#gamma}^{2} = 20|#Psi_{Pair}| < 0.05";
        case 25:
            return "2D: #chi_{#gamma}^{2} = 30 |#Psi_{Pair}| < 0.1";
        case 27:
            return "2D: #chi_{#gamma}^{2} = 30 |#Psi_{Pair}| < 0.035";
        default:
            return " #chi_{#gamma}^{2} #Psi_{Pair} cut not defined";
    }
}  

//************************************************************************************
//***** Analyzes the photon quality gamma cut, return correct cut label **************
//************************************************************************************
TString AnalysePhotonQuality(Int_t photonQualitCut ){
    switch(photonQualitCut) {
        case 0:
            return "photon Quality: 1,2,3";
        case 2:
            return "photon Quality: 1 (TPC only photons)";
        case 3:
            return "photon Quality: 2 (1 leg with #geq 2 ITS hits)";
        case 4:
            return "photon Quality: 3 (both legs with #geq 2 ITS hits)";
        default:
            return "photon Quality cut not defined";
    }
}  

//************************************************************************************
//******** Analyzes the phi exclusion photon cuts, return correct cut label **********
//************************************************************************************
TString AnalyseConvPhiExclusionCut(TString gammaCutNumber){

    Int_t minCutNumberPhi   = GetPhotonMinPhiCutPosition(gammaCutNumber);
    if (minCutNumberPhi == -1) return "Full TPC acceptance";
    
    TString minPhiCutNumber(gammaCutNumber(minCutNumberPhi,1));
    TString maxPhiCutNumber(gammaCutNumber(GetPhotonMaxPhiCutPosition(gammaCutNumber),1));
    TString etaExclusion(gammaCutNumber(GetPhotonEtaForPhiCutPosition(gammaCutNumber),1));
    
    TString etaAcc          = "";
    if (etaExclusion.CompareTo("1") == 0) etaAcc = "#eta < 0" ;
        else if (etaExclusion.CompareTo("2") == 0) etaAcc = "#eta > 0" ;
        
    if ( etaExclusion.CompareTo("0") == 0 && minPhiCutNumber.CompareTo("0")==0 && maxPhiCutNumber.CompareTo("0")==0 ){
        return "Full TPC acceptance";
    } else if (minPhiCutNumber.CompareTo("0")==0 && maxPhiCutNumber.CompareTo("0")==0){
        return Form("Full $varphi acceptance, %s",etaAcc.Data()) ;
    }	
        
    Double_t fMinPhiCut     = 0;
    Double_t fMaxPhiCut     = 2* TMath::Pi();
    switch(minPhiCutNumber.Atoi()) {
        case 0:
            fMinPhiCut = 0; //no cut on phi
            break;
        case 1:
            fMinPhiCut = 1.7; //OROC C08
            break;
        case 2:
            fMinPhiCut = 4.4; //EMCal
            break;
        case 3:
            fMinPhiCut = 1.0; //PHOS
            break;
        case 4:
            fMinPhiCut = 3.4; //EMCal tight
            break;
        case 5:
            fMinPhiCut = 2.0; //OROC C08 medium cut 
            break;
        case 6:
            fMinPhiCut = 2.2; //OROC C08 small cut
            break;
        case 7:
            fMinPhiCut = 2.4; //OROC C08 tightest cut
            break;

        default:
            fMinPhiCut = 0.;
            break;
    }	

    switch (maxPhiCutNumber.Atoi()) {
        case 0:
            fMaxPhiCut = 2* TMath::Pi(); //no cut
            break;
        case 1:
            fMaxPhiCut = 4.3; //OROC C08
            break;
        case 2:
            fMaxPhiCut = 5.8; //EMCal
            break;
        case 3:
            fMaxPhiCut = 3.0; //PHOS
            break;
        case 4:
            fMaxPhiCut = 1.; //EMCal
            break;
        case 5:
            fMaxPhiCut = 4.0; //OROC C08 medium cut 
            break;
        case 6:
            fMaxPhiCut = 3.8; //OROC C08 small cut
            break;
        case 7:
            fMaxPhiCut = 3.6; //OROC C08 tighest cut
            break;

        default:
            fMaxPhiCut = 2* TMath::Pi(); //no cut
            break;
    }

    if (fMinPhiCut < fMaxPhiCut){
        if (etaAcc.CompareTo("") == 0)
            return Form("exclude #gamma_{conv} in %1.3f < #varphi < %1.3f", fMinPhiCut, fMaxPhiCut);
        else 
            return Form("exclude #gamma_{conv} in %1.3f < #varphi < %1.3f, %s", fMinPhiCut, fMaxPhiCut, etaAcc.Data());
    } else {
        return Form("restrict #gamma_{conv} to %1.3f < #varphi < %1.3f",  fMaxPhiCut, fMinPhiCut );
    }	
    
    return "";
}


//************************************************************************************
//******** Analyzes the trigger event cut, return correct cut label ******************
//************************************************************************************
TString AnalyseSpecialTriggerCut(Int_t SpecialTrigger, TString periodName = "No"){ 
    // Set Cut
    switch(SpecialTrigger){
        case 0:
            if (periodName.CompareTo("LHC11a")==0||periodName.CompareTo("LHC11a_p4_woSDD")==0)
                return "without SDD, V0OR";
            else if (periodName.Contains("LHC12") || periodName.Contains("LHC13") )
                return "V0AND";
            else 	
                return "V0OR";
        case 1:
            if (periodName.CompareTo("LHC11a")==0||periodName.CompareTo("LHC11a_p4_woSDD")==0)
                return "without SDD, V0OR";
            else 
                return "V0OR";
        case 3:
            if (periodName.CompareTo("LHC11a")==0||periodName.CompareTo("LHC11a_p4_wSDD")==0)
                return "with SDD, V0OR";
            else 
                return "V0OR";
        case 10:
            return "V0AND";
        case 11:
            return "T0AND";
        case 12:
            if (periodName.CompareTo("LHC11a")==0||periodName.CompareTo("LHC11a_p4_wSDD")==0)
                return "V0AND";
        case 13:
            if (periodName.CompareTo("LHC11a")==0||periodName.CompareTo("LHC11a_p4_wSDD")==0)
                return "with SDD, V0AND";
            else 	
                return "V0AND";
        case 51: 	
            return "EMC L0, INT1";
        case 52: 	
            return "EMC L0, INT7";
        case 53: 	
            return "EMC L0, INT8";
        case 81: 	
            return "EMC L1-GA, INT7";
        case 83: 	
            return "EMC L1-G1, INT7";
        case 85: 	
            return "EMC L1-G2, INT7";
        case 91: 	
            return "EMC L1-JE, INT7";
        case 93: 	
            return "EMC L1-J1, INT7";
        case 95: 	
            return "EMC L1-J2, INT7";

        default:
            return "special Trigger cut not defined";
    }
}

//************************************************************************************
//******** Analyzes the multiplicity cut in pp events ********************************
//************************************************************************************
TString AnalysePPMultiplicityCut(Int_t minMult, Int_t maxMult){ 
    // Set Cut
    Int_t multBins[9] =  {0,   2,   5,    10,   15,  30,  50,  100,  1000 };
    
    TString multiplicityString = "";
    if ( (minMult == 0 && maxMult == 0) || maxMult<minMult){
        multiplicityString = "No selection";
    } else {
        multiplicityString =  Form("%d #leq TPC track mult < %d", multBins[minMult], multBins[maxMult]);
    }    
    return multiplicityString;
}


//************************************************************************************
//******** Analyzes the cluster track matching cuts, return correct cut label ********
//************************************************************************************
TString AnalyseTrackMatchingCut(Int_t trackmatching ){
    switch(trackmatching) {
        case 0:
            return "no track missmatch cut for V0s";
        case 1:
            return "TrMatch excl. #Delta#eta < 0.008, -0.03 < #Delta#varphi_{+} < 0.03, -0.03 < #Delta#varphi_{-} < 0.03 for V0s";
        case 2:
            return "TrMatch excl. #Delta#eta < 0.012, -0.05 < #Delta#varphi_{+} < 0.04, -0.04 < #Delta#varphi_{-} < 0.05 for V0s";
        case 3:
            return "TrMatch excl. #Delta#eta < 0.016, -0.09 < #Delta#varphi_{+} < 0.06, -0.06 < #Delta#varphi_{-} < 0.09 for V0s";
        case 4:
            return "TrMatch excl. #Delta#eta < 0.018, -0.11 < #Delta#varphi_{+} < 0.07, -0.07 < #Delta#varphi_{-} < 0.11 for V0s";
        case 5:
            return "TrMatch excl. #Delta#eta < 0.020, -0.13 < #Delta#varphi_{+} < 0.08, -0.08 < #Delta#varphi_{-} < 0.13 for V0s";
        case 6:
            return "TrMatch excl. #Delta#eta < 0.022, -0.15 < #Delta#varphi_{+} < 0.10, -0.10 < #Delta#varphi_{-} < 0.15 for V0s";
        case 7:
            return "TrMatch excl. #Delta#eta < 0.005, -0.03 < #Delta#varphi_{+} < 0.03, -0.03 < #Delta#varphi_{-} < 0.03 for V0s";
        case 8:
            return "TrMatch excl. #Delta#eta < 0.010, -0.09 < #Delta#varphi_{+} < 0.07, -0.07 < #Delta#varphi_{-} < 0.09 for V0s";
        case 9:
            return "TrMatch excl. #Delta#eta < 0.015, -0.15 < #Delta#varphi_{+} < 0.11, -0.11 < #Delta#varphi_{-} < 0.15 for V0s";
        default:
            return "track missmatch not defined ";
    }
}  

//************************************************************************************
//**** Analyzes the cluster track matching cuts for calo, return correct cut label ***
//************************************************************************************
TString AnalyseTrackMatchingCaloCut(Int_t trackmatching ){
    switch(trackmatching) {
        case 0:
            return "no track matching exclusion";
        case 1:
            return "TrMatch excl. #Delta#eta < 0.008, -0.03 < #Delta#varphi_{+} < 0.03, -0.03 < #Delta#varphi_{-} < 0.03";
        case 2:
            return "TrMatch excl. #Delta#eta < 0.012, -0.05 < #Delta#varphi_{+} < 0.04, -0.04 < #Delta#varphi_{-} < 0.05";
        case 3:
            return "TrMatch excl. #Delta#eta < 0.016, -0.09 < #Delta#varphi_{+} < 0.06, -0.06 < #Delta#varphi_{-} < 0.09";
        case 4:
            return "TrMatch excl. #Delta#eta < 0.018, -0.11 < #Delta#varphi_{+} < 0.07, -0.07 < #Delta#varphi_{-} < 0.11";
        case 5:
            return "TrMatch excl. #Delta#eta < 0.020, -0.13 < #Delta#varphi_{+} < 0.08, -0.08 < #Delta#varphi_{-} < 0.13";
        case 6:
            return "TrMatch excl. #Delta#eta < 0.022, -0.15 < #Delta#varphi_{+} < 0.10, -0.10 < #Delta#varphi_{-} < 0.15";
        case 7:
            return "TrMatch excl. #Delta#eta < 0.005, -0.03 < #Delta#varphi_{+} < 0.03, -0.03 < #Delta#varphi_{-} < 0.03";
        case 8:
            return "TrMatch excl. #Delta#eta < 0.010, -0.09 < #Delta#varphi_{+} < 0.07, -0.07 < #Delta#varphi_{-} < 0.09";
        case 9:
            return "TrMatch excl. #Delta#eta < 0.015, -0.15 < #Delta#varphi_{+} < 0.11, -0.11 < #Delta#varphi_{-} < 0.15";
        default:
            return "track missmatch not defined ";
    }
}  

//************************************************************************************
//***************** Analyzes the min energy, return correct cut value ****************
//************************************************************************************
Double_t ReturnMinClusterEnergy(TString clusterCutSelection){
    TString minEnergyCutNumber(clusterCutSelection(GetClusterMinEnergyCutPosition(clusterCutSelection),1));
    Int_t minEnergyCut          = minEnergyCutNumber.Atoi();
    switch(minEnergyCut){
        case 0: 
            return 0;
            break;
        case 1: 
            return 0.5; 
            break;
        case 2: 
            return 0.6; 
            break;
        case 3: 
            return 0.7; 
            break;
        case 4: 
            return 0.8; 
            break;
        case 5: 
            return 0.9; 
            break;
        case 6: 
            return 4.5; 
            break;
        case 7: 
            return 5.0; 
            break;
        case 8: 
            return 5.5; 
            break;
        case 9: 
            return 6.0; 
            break;
        default:
            return 0;
    }
}

//************************************************************************************
//***************** Analyzes the min energy, return correct cut value ****************
//************************************************************************************
Double_t ReturnMinNCells(TString clusterCutSelection){
    TString nCellsCutNumber(clusterCutSelection(GetClusterMinNCellsCutPosition(clusterCutSelection),1));
    return nCellsCutNumber.Atoi();
}

//************************************************************************************
//***************** Analyzes the N cells cut, return correct cut label ****************
//************************************************************************************
TString AnalyseNCellsCut (Int_t nCellsCut){
    switch(nCellsCut){
        case 0:
            return "No Ncells cut";
        case 1:
            return "n Cells > 1";
        case 2:
            return "n Cells > 2";
        case 3:
            return "n Cells > 3";
        case 4:
            return "n Cells > 4";
        case 5:
            return "n Cells > 5";
        case 6:
            return "n Cells > 6";
        default:
            return "n Cells cut not defined";
    }
}

//************************************************************************************
//********************* Analyzes the M02 cut, return correct cut label ***************
//************************************************************************************
TString AnalyseM02Cut(Int_t minM02, Int_t maxM02){
    Double_t fMinM02 = 0.;
    cout << minM02 << "\t" << maxM02 << endl;
    switch(minM02){
        case 0: 
            fMinM02=0;
            break;
        case 1: 
            fMinM02=0.002; 
            break;
        case 2: 
            fMinM02=0.1; 
            break;
        case 3: 
            fMinM02=0.2; 
            break;
        default:
            fMinM02 = -10;
            break;
      
    }
    
    Double_t fMaxM02 = 1000.;
    switch(maxM02){
        case 0: 
            fMaxM02=1000;
            break;
        case 1: 
            fMaxM02=1.; 
            break;
        case 2: 
            fMaxM02=0.7; 
            break;
        case 3: 
            fMaxM02=0.5; 
            break;
        case 4: 
            fMaxM02=0.4; 
            break;
        default:
            fMaxM02 = -10;
            break;
    }
    return Form("%1.3f < M_{02} < %3.3f", fMinM02, fMaxM02);
  
    
}	

//************************************************************************************
//********************* Analyzes the phi cut, return correct cut value ***************
//************************************************************************************
TString AnalyseAcceptanceCutPhiCluster(Int_t minPhi, Int_t maxPhi){
    if (minPhi == 0 && maxPhi == 0){
        return "Full acceptance";
    } else if (minPhi == 1 && maxPhi ==1){
        return "restricted to EMCal acceptance";
    } else if (minPhi == 2 && maxPhi ==1){
        return "EMCal 2012/13 geometry with TRD (2.1 < #varphi < 3.15)";
    } else if (minPhi == 1 && maxPhi ==3){
        return "EMCal 2012/13 geometry no TRD (1.39626 < #varphi < 2.1)";
    } else if (minPhi == 3 && maxPhi ==1){
        return "EMCal 2011 geometry with TRD (2.45 < #varphi < 3.15)";
    } else if (minPhi == 1 && maxPhi ==2){
        return "EMCal 2011 geometry no TRD (1.39636 < #varphi < 2.45)";
    }	
    return "";
}

//************************************************************************************
//***************** Analyzes the min energy cut, return correct cut label ************
//************************************************************************************
TString AnalyseMinEnergyCut(Int_t minEnergyCut){
    Double_t fMinEnergy = 0.;
    switch(minEnergyCut){
        case 0: 
            fMinEnergy = 0;
            break;
        case 1: 
            fMinEnergy = 0.5; 
            break;
        case 2: 
            fMinEnergy = 0.6; 
            break;
        case 3: 
            fMinEnergy = 0.7; 
            break;
        case 4: 
            fMinEnergy = 0.8; 
            break;
        case 5: 
            fMinEnergy = 0.9; 
            break;
        case 6: 
            fMinEnergy = 4.5; 
            break;
        case 7: 
            fMinEnergy = 5.0; 
            break;
        case 8: 
            fMinEnergy = 5.5; 
            break;
        case 9: 
            fMinEnergy = 6.0; 
            break;
        default:
            fMinEnergy = 0;
            break;
    }
    
    return Form("E_{clus} > %3.3f GeV/c", fMinEnergy);
}

//************************************************************************************
//***************** Analyzes the min energy cut, return correct cut label ************
//************************************************************************************
TString AnalyseMinEnergyMergedCut(Int_t minEnergyCut){
    Double_t fMinEnergy = 0.;
    switch(minEnergyCut){
        case 0: 
            fMinEnergy = 0.1;
            break;
        case 1: 
            fMinEnergy = 4.0; 
            break;
        case 2: 
            fMinEnergy = 5.0; 
            break;
        case 3: 
            fMinEnergy = 6.0; 
            break;
        case 4: 
            fMinEnergy = 7.0; 
            break;
        case 5: 
            fMinEnergy = 7.5; 
            break;
        case 6: 
            fMinEnergy = 8.0; 
            break;
        case 7: 
            fMinEnergy = 8.5; 
            break;
        case 8: 
            fMinEnergy = 9.0; 
            break;
        case 9: 
            fMinEnergy = 9.5; 
            break;
        default:
            fMinEnergy = 0;
            break;
    }
    
    return Form("E_{clus} > %3.3f GeV/c", fMinEnergy);
}


//************************************************************************************
//** Analyzes the cluster timing cut, return correct cut label **************************
//************************************************************************************
TString AnalyseClusterTimingCut(Int_t timing){
    switch(timing){
        case 0: 
            return "-500 s < t_{clus} < 500 s";
        case 1: 
            return "-1 #mu s < t_{clus} < 1 #mu s";
        case 2: 
            return "-500 ns < t_{clus} < 500 ns";
        case 3: 
            return "-200 ns < t_{clus} < 200 ns";
        case 4: 
            return "-100 ns < t_{clus} < 100 ns";
        case 5: 
            return "-50 ns < t_{clus} < 50 ns";
        case 6:
            return "-30 ns < t_{clus} < 35 ns";
        case 7:
            return "-30 ns < t_{clus} < 30 ns";
        case 8:
            return "-20 ns < t_{clus} < 30 ns";
        case 9:
            return "-20 ns < t_{clus} < 25 ns";
        default:
            return "no timing cut defined";
    }
}

//************************************************************************************
//** Analyzes the cluster timing cut, return correct cut label **************************
//************************************************************************************
TString AnalyseClusterNonLinearityCut(Int_t nonlinearity){

  
//   return Form("non linearity cut %d",nonlinearity);
  switch(nonlinearity){
      case 0: 
          return "none";
      case 1: 
          return "SDM (Jason)";
      case 11:
          return "ConvCalo, tight timing";
      case 12:
          return "SDM w/o data corr (new), tight timing";
      case 13:
          return "TBv3 w/ ConvCalo-mass corr";
      case 14:
          return "TBv3 w/ Calo-mass corr";
      case 15: 
          return "ConvCalo with data shifted";
      case 16: 
          return "SDM w/ data corr (new)";
      case 21:
          return "ConvCalo, mass param, open timing";
      case 22:
          return "SDM w/o data corr (new), mass param, open timing";
      case 31:
          return "ConvCalo, open timing";
      case 32:
          return "SDM w/o data corr (new), open timing";
      default:
          return Form("non linearity cut %d",nonlinearity);
  }    
  //  switch(nonlinearity){
//    case 0:
//      return "from tender?";
//    default:
//      return "no nonlinearity defined";
//  }
}

//************************************************************************************
//***************************** Returns the cluster NLM ******************************
//************************************************************************************
Int_t ReturnClusterNLM( TString clusterCutString){
    TString nCellsCutNumber(clusterCutString(GetClusterNLMCutPosition(clusterCutString),1));
    return nCellsCutNumber.Atoi();
}
    
//************************************************************************************
//** Analyzes the alpha meson cut, return correct cut label **************************
//************************************************************************************
TString AnalyseAlphaMesonCut(Int_t alphaMesonCut){ 
    switch(alphaMesonCut){
        case 0:	// 0- 0.7
            return "|#alpha_{meson}| < 0.7";
        case 1:	// 0-0.5
            return "0.65*tanh(1.8*x)"; //"|#alpha_{meson}| < 0.5";
        case 2:	// 0.5-1
            return "0.5 < |#alpha_{meson}| < 1";
        case 3:	// 0.0-1
            return "|#alpha_{meson}| < 1";
        case 4:	// 0-0.65
            return "|#alpha_{meson}| < 0.65";
        case 5:	// 0-0.75
            return "|#alpha_{meson}| < 0.75";
        case 6:	// 0-0.8
            return "|#alpha_{meson}| < 0.8";
        case 7:	// 0.0-0.85
            return "|#alpha_{meson}| < 0.85";
        case 8:	// 0.0-0.6
            return "|#alpha_{meson}| < 0.6";
        case 9:	// 0.0-0.3
            return "|#alpha_{meson}| < 0.3";
        default:
            return "no alpha cut defined";
    }
}

//************************************************************************************
//** Analyzes the MC smearing cut for mesons, return correct cut label ***************
//************************************************************************************
TString AnalyseMCSmearingCut(Int_t mcSmearingCut){ 
    switch(mcSmearingCut){
        case 0:
            return "no additional smearing" ;
        case 1:
            return "#sqrt{0.011^{2} + 0.007^{2} #times P_{#gamma}^{2}} #times Rand.Gaus(0,1)" ;
        case 2:
            return "#sqrt{0.022^{2} + 0.007^{2} #times P_{#gamma}^{2}} #times Rand.Gaus(0,1)" ;
        case 3:
            return "#sqrt{0.044^{2} + 0.007^{2} #times P_{#gamma}^{2}} #times Rand.Gaus(0,1)" ;
        case 7:
            return "#sqrt{0.011^{2} + 0.014^{2} #times P_{#gamma}^{2}} #times Rand.Gaus(0,1)" ;
        case 8:
            return "#sqrt{0.011^{2} + 0.0035^{2} #times P_{#gamma}^{2}} #times Rand.Gaus(0,1)" ;
        case 9:
            return "#sqrt{0.011^{2} + 0.028^{2} #times P_{#gamma}^{2}} #times Rand.Gaus(0,1)" ;
        default:
            return "smearing cut not defined";
    }
}

//************************************************************************************
//********* Analyzes the meson BG scheme, return correct cut label *******************
//************************************************************************************
TString AnalyseBackgroundScheme(TString sBackgroundScheme){
    TString sBackgroundSchemeB      = sBackgroundScheme(0,1);
    Int_t BackgroundScheme          = sBackgroundSchemeB.Atoi();
    TString bGScheme                = "";
    cout << "BackgroundScheme: " << BackgroundScheme << endl;
    
    switch(BackgroundScheme){
        case 0: //Rotation
            bGScheme="Rotation";
            break;
        case 1: // mixed event with V0 multiplicity
            bGScheme="mixed event, V0 mult";
            break;
        case 2: // mixed event with track multiplicity
            bGScheme="mixed event, track mult";
            break;
        case 3: //Rotation
            bGScheme="Rotation, with prob.";
            break;
        case 4: //No BG calculation
            bGScheme="no BG calculated";
            break;
        case 5: //Rotation
            bGScheme="Rotation, new BG handler";
            break;
        case 6: // mixed event with V0 multiplicity
            bGScheme="mixed event, new BG handler, V0 mult";
            break;
        case 7: // mixed event with track multiplicity
            bGScheme="mixed event, new BG handler, track mult";
            break;
        case 8: //Rotation
            bGScheme="Rotation, new BG handler, with prob. ";
            break;
        default:
            bGScheme="no BG method selected";
        
    }
    
    TString sNumberOfBGEvents       = sBackgroundScheme(1,1);
    Int_t NumberOfBGEvents          = sNumberOfBGEvents.Atoi();
    Int_t fNumberOfBGEvents         = 0;
    cout << "NumberOfBGEvents: " << NumberOfBGEvents << endl;

    switch(NumberOfBGEvents){
        case 0:
            fNumberOfBGEvents = 5;
            break;
        case 1:
            fNumberOfBGEvents = 10;
            break;
        case 2:
            fNumberOfBGEvents = 15;
            break;
        case 3:
            fNumberOfBGEvents = 20;
            break;
        case 4:
            fNumberOfBGEvents = 2;
            break;
        case 5:
            fNumberOfBGEvents = 50;
            break;
        case 6:
            fNumberOfBGEvents = 80;
            break;
        case 7:
            fNumberOfBGEvents = 100;
            break;
        default:
            cout<<"Warning: NumberOfBGEvents not defined "<<NumberOfBGEvents<<endl;
    }

    TString sDegreesForRotationMethod   = sBackgroundScheme(1,2);
    Int_t DegreesForRotationMethod      = sDegreesForRotationMethod.Atoi();
    Int_t fnDegreeRotationPMForBG       = 0;
    cout << "DegreesForRotationMethod: " << DegreesForRotationMethod << endl;
    
    switch(DegreesForRotationMethod){
        case 0:
            fnDegreeRotationPMForBG = 5;
            break;
        case 1:
            fnDegreeRotationPMForBG = 10;
            break;
        case 2:
            fnDegreeRotationPMForBG = 15;
            break;
        case 3:
            fnDegreeRotationPMForBG = 20;
            break;
        default:
            cout<<"Warning: DegreesForRotationMethod not defined "<<DegreesForRotationMethod<<endl;
    }

    if (BackgroundScheme == 0 || BackgroundScheme == 3 || BackgroundScheme == 5 || BackgroundScheme == 8){
        return Form("%s within #pm %i°, %i photons per pool",bGScheme.Data(), fnDegreeRotationPMForBG, fNumberOfBGEvents)   ;
    } else {
        return Form("%s with %i photons per pool",bGScheme.Data(), fNumberOfBGEvents)   ;
    }
    return "";
}

//************************************************************************************
//*********** Analyzes the rapidity meson cut, return correct cut label **************
//************************************************************************************
TString AnalyseRapidityMesonCut(Int_t RapidityMesonCut){ 
    // Set Cut
    switch(RapidityMesonCut){
        case 0:  //
            return "|y_{meson}| < 0.9";
        case 1:  //
            return "|y_{meson}| < 0.8";
        case 2:  //
            return "|y_{meson}| < 0.7";
        case 3:  //
            return "|y_{meson}| < 0.6";
        case 4:  //
            return "|y_{meson}| < 0.5";
        case 5:  //
            return "|y_{meson}| < 0.65";
        case 6:  //
            return "|y_{meson}| < 0.75";
        case 7:  //
            return "|y_{meson}| < 0.3";
        case 8:  //
            return "|y_{meson}| < 0.35";
        case 9:  //
            return "|y_{meson}| < 0.4";        
        default:
            return "rapidity cut not defined";
    }
}

//************************************************************************************
//*********** Analyzes the rapidity meson cut for pPv, return correct cut label ******
//************************************************************************************
TString AnalyseRapidityMesonCutpPb(Int_t RapidityMesonCut){ 
    // Set Cut
    switch(RapidityMesonCut){
        case 0:  //
            return "|y_{meson}| < 0.9";
        case 1:  //
            return "|y_{meson}| < 0.8";
        case 2:  //
            return "|y_{meson}| < 0.7";
        case 3:  //
            return "|y_{meson}| < 0.6";
        case 4:  //
            return "|y_{meson}| < 0.5";
        case 5:  //
            return "|y_{meson}| < 0.65";
        case 6:  //
            return "|y_{meson}| < 0.75";
        case 7:  //
            return "0.165 < y_{c.m.s, meson} < 0.765";
        case 8:  //
            return "|y_{meson}| < 0.35";
        case 9:  //
            return "|y_{meson}| < 0.4";        
        default:
            return "rapidity cut not defined";
    }
}

//************************************************************************************
//******** Analyzes the chi2 meson cut, return correct cut label *********************
//************************************************************************************
TString AnalyseChi2MesonCut(Int_t chi2GammaCut){   // Set Cut
    switch(chi2GammaCut){
        case 0: // 100
            return "#chi_{meson}^{2} < 100";
        case 1:  // 50
            return "#chi_{meson}^{2} < 50";
        case 2:  // 30
            return "#chi_{meson}^{2} < 30";
        case 3:
            return "#chi_{meson}^{2} < 200";
        case 4:
            return "#chi_{meson}^{2} < 500";
        case 5:
            return "#chi_{meson}^{2} < 1000";  
        default:
            return "#chi_{meson}^{2} cut unknown";
    }
}

//************************************************************************************
//******** Analyzes the openign angle meson cut, return correct cut label *********************
//************************************************************************************
 TString AnalyseMesonOpeningAngleCut(Int_t OpeningAngleCut){   // Set Cut
   switch(OpeningAngleCut){
     case 0:
       return "#theta_{meson} > 0";
     case 1:
       return "#theta_{meson} > 0.005";
     case 2:
       return "#theta_{meson} > f(#it{p}_{T})";
     case 3:
       return "#theta_{meson} > 0.01";
     case 4:
       return "#theta_{meson} > 0.0152";
     case 5:
       return "#theta_{meson} > 0.0202";
     case 6:
       return "#theta_{meson} > 0.0404";
     default:
       return "#theta_{meson} cut unknown";
   }
}

// ******************************************************************************************
// ************************ Set correct xSection for pp *************************************
// ******************************************************************************************
Double_t xSection8TeVV0AND      = 55.74*1e-3;   // from https://indico.cern.ch/event/276276/contribution/3/material/slides/0.pdf
Double_t xSection8TeVErrUp      = 0.46;         // from https://indico.cern.ch/event/276276/contribution/3/material/slides/0.pdf
Double_t xSection8TeVErrDown    = 0.46;         // from https://indico.cern.ch/event/276276/contribution/3/material/slides/0.pdf
Double_t xSection8TeVT0AND      = 55.74*1e-3;   // from https://indico.cern.ch/event/276276/contribution/3/material/slides/0.pdf
Double_t xSection8TeVT0ErrUp    = 0.46;         // from https://indico.cern.ch/event/276276/contribution/3/material/slides/0.pdf
Double_t xSection8TeVT0ErrDown  = 0.46;         // from https://indico.cern.ch/event/276276/contribution/3/material/slides/0.pdf
Double_t xSection7TeV           = 62.22*1e-3;
Double_t xSection7TeVV0AND      = 54.31*1e-3;
Double_t xSection7TeVErrUp      = 2.18;
Double_t xSection7TeVErrDown    = 2.18;
Double_t xSection900GeV         = 47.78*1e-3;
Double_t xSection900GeVV0AND    = 40.06*1e-3;
Double_t xSection900GeVErrUp    = 2.39;
Double_t xSection900GeVErrDown  = 1.86;
Double_t xSection2760GeV        = 55.416*1e-3;
Double_t xSection2760GeVV0AND   = 47.73*1e-3;
Double_t xSection2760GeVINEL    = 62.8*1e9;
Double_t xSection2760GeVErr     = 3.9;
Double_t recalcBarn             = 1e12;         //NLO in pbarn!!!!
Double_t factorToInel           = 1/1.12;       // this factor is multiplied with Raa and comes from trigger inelastic effiency for pp

Double_t ReturnCorrectXSection ( TString energy, 
                                 Int_t isV0AND){
    
Double_t xSectionInt = 0;    
    if(energy.CompareTo("7TeV") == 0){
        if (isV0AND == 1){
            xSectionInt = xSection7TeVV0AND;
            cout << "V0AND xSection taken: \t" << xSectionInt << endl;
        } else {
            xSectionInt = xSection7TeV;
            cout << "V0OR xSection taken: \t" << xSectionInt << endl;
        }    
    } else if(energy.CompareTo("2.76TeV") == 0){
        if (isV0AND == 1){
            xSectionInt = xSection2760GeVV0AND;
            cout << "V0AND xSection taken: \t" << xSectionInt << endl;
        } else if (isV0AND == 3){
            xSectionInt = xSection2760GeVINEL;
            cout << "V0AND xSection taken: \t" << xSectionInt << endl;            
        } else {    
            xSectionInt = xSection2760GeV;
            cout << "V0OR xSection taken: \t" << xSectionInt << endl;  
        }
    } else if(energy.CompareTo("900GeV") == 0){
        if (isV0AND == 1){
            xSectionInt = xSection900GeVV0AND;
            cout << "V0AND xSection taken: \t" << xSectionInt << endl;
        } else {
            xSectionInt = xSection900GeV;
            cout << "V0OR xSection taken: \t" << xSectionInt << endl;
        }    
    } else if(energy.CompareTo("8TeV") == 0){
        if (isV0AND == 1){
            xSectionInt = xSection8TeVV0AND;
            cout << "V0AND xSection taken: \t" << xSectionInt << endl;
        } else if (isV0AND == 2){
            xSectionInt = xSection8TeVT0AND;
            cout << "T0AND xSection taken: \t" << xSectionInt << endl;
        } else {
            cout << "ERROR: V0OR xSection not deterimined, set to \t" << xSectionInt << endl;
        }
    } else {
        cout << "ERROR: energy not deterimined, xsection set to \t" << xSectionInt << endl;
    }    
    return xSectionInt;
}


//************************************************************************************
//* Decodes from the mode the respective reco process and return correct label + details 
//************************************************************************************
TString ReturnFullTextReconstructionProcess( Int_t mode, Int_t separate = 0, TString meson = "", TString clusterCutNumber = "" ){ 
    if (separate == 0){
        switch (mode){
            case 0:
                return "#gamma's rec. with PCM";
            case 1: 
                return "#gamma's rec. with PCM, Dalitz";
            case 2:
                return "#gamma's rec. with PCM, EMCal";
            case 3: 
                return "#gamma's rec. with PCM, PHOS";
            case 4:
                return "#gamma's rec. with EMCal";
            case 5:
                return "#gamma's rec. with PHOS";
            case 6: 
                return "#gamma's rec. with EMCal, Dalitz";
            case 7:
                return "#gamma's rec. with PHOS, Dalitz";
            case 10:
                if (clusterCutNumber.CompareTo("") != 0){
                    Int_t nlm = ReturnClusterNLM(clusterCutNumber);
                    if (meson.CompareTo("") == 0) return Form("rec. w/ EMCal m. cl., %d lm", nlm );
                    else return Form("%s rec. w/ EMCal m. cl., %d lm", meson.Data(), nlm);                    
                } else {    
                    if (meson.CompareTo("") == 0) return "rec. w/ EMCal m. cl.";
                    else return Form("%s rec. w/ EMCal m. cl.", meson.Data());
                }    
            case 11:
                if (clusterCutNumber.CompareTo("") != 0){
                    Int_t nlm = ReturnClusterNLM(clusterCutNumber);
                    if (meson.CompareTo("") == 0) return Form("rec. w/ PHOS m. cl., %d lm", nlm );
                    else return Form("%s rec. w/ PHOS m. cl., %d lm", meson.Data(), nlm);                    
                } else {    
                    if (meson.CompareTo("") == 0) return "rec. w/ PHOS m. cl.";
                    else return Form("%s rec. w/ PHOS m. cl.", meson.Data());
                }    
            default:
                return "not known";
        }
    } else if (separate == 1){
        switch (mode){
            case 0:
            case 1:
            case 2:
            case 3:    
                return "#gamma's rec. with PCM";
            case 6:
            case 7:    
                return "#gamma*'s rec. with Dalitz";
            default:
                return "not known";
        }
    } else if (separate == 2){
        switch (mode){
            case 1:
                return "#gamma*'s rec. with Dalitz";
            case 0:
            case 2:
            case 4:
            case 6:
                return "#gamma's rec. with EMCal";
            case 3:
            case 5:
            case 7:
                return "#gamma's rec. with PHOS";

            default:
                return "not known";
        }    
    }    
    
}
