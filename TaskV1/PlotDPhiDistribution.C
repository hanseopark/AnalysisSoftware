/**************************************************************************************************************************
 ****** 		provided by Gamma Conversion Group, PWG4, 																						*****
 ******		Friederike Bock, friederike.bock@cern.ch																							*****
 *****************************************************************************************************************************
 *** This macro can be used to display the Photon Characteristics of the conversion method in ALICE, it can be operated  *****
 *** on the output of the GammaConversionTask. It can take 2 input files, the second one should be MC, if this is not 	 *****
 *** the case all histograms including MC need to be commented out otherwise the running will crash.							 *****
 ****************************************************************************************************************************/

#include <Riostream.h>
#include <fstream>
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
#include "TH3F.h"
#include "TF1.h"
#include "TVirtualFitter.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TDatabasePDG.h"
#include "TMinuit.h"
#include "TLatex.h"
#include "TASImage.h"
#include "TPostScript.h"
#include "TGraphErrors.h"
#include "TArrow.h"
#include "TMarker.h"
#include "TGraphAsymmErrors.h" 
#include "../CommonHeaders/PlottingGammaConversionHistos.h"
#include "../CommonHeaders/PlottingGammaConversionAdditional.h"
#include "../CommonHeaders/FittingGammaConversion.h"
#include "../CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "../CommonHeaders/ConversionFunctions.h"

TString collisionSystem;
TString textDate;

void PlotStandard2D( TH2* histo2D, TString nameOutput, TString title, TString xTitle, TString yTitle, Bool_t kRangeY, Double_t startY, Double_t endY, Bool_t kRangeX, Double_t startX, Double_t endX, Int_t logX, Int_t logZ, Float_t* floatLogo, Int_t canvasSizeX = 500, Int_t canvasSizeY = 500, TString generator ="" , TString textPeriod ="", TString drawingAdditional = ""){
    TCanvas * canvasStandard = new TCanvas("canvasStandard","",10,10,canvasSizeX,canvasSizeY);  // gives the page size		
    canvasStandard->SetLogx(logX);
    canvasStandard->SetLogz(logZ);
    canvasStandard->SetRightMargin(0.12); 		
    canvasStandard->SetLeftMargin(0.12); 		
    canvasStandard->SetBottomMargin(0.1); 		
    canvasStandard->SetTopMargin(0.04); 		
    canvasStandard->cd();
    histo2D->SetTitle("");
    DrawAutoGammaHisto2D(	histo2D,
                                    title.Data(), xTitle.Data(), yTitle.Data(),"",kRangeY, startY, endY, kRangeX, startX, endX);
    cout << generator.Data() << endl;
    cout << textPeriod.Data() << endl;
    TString labeltext = "";
    if (!drawingAdditional.CompareTo("dPhiEPi")){
        labeltext = "#gamma_{fake} rec from comb #pi + e";
    } else if (!drawingAdditional.CompareTo("dPhiEK")){
        labeltext = "#gamma_{fake} rec from comb K + e";
    } else if (!drawingAdditional.CompareTo("dPhiEP")){
        labeltext = "#gamma_{fake} rec from comb P + e";
    } else if (!drawingAdditional.CompareTo("dPhiPiK")){
        labeltext = "#gamma_{fake} rec from comb K + #pi";
    } else if (!drawingAdditional.CompareTo("dPhiEK")){
        labeltext = "#gamma_{fake} rec from comb P + #pi";
    }		
    
    if (labeltext.CompareTo("")){
        TLatex *add = new TLatex(floatLogo[0],(floatLogo[1]-3*floatLogo[2]*1.15),Form("%s",labeltext.Data())); // Bo: this was modified
        add->SetNDC();
        add->SetTextColor(1);
        add->SetTextFont(62);
        add->SetTextSize(floatLogo[2]);
        add->SetLineWidth(2);
        add->Draw("same");	

    }	
    DrawLabelsEvents(floatLogo[0],floatLogo[1],floatLogo[2], 0.00, collisionSystem, generator, textPeriod);
        
    canvasStandard->Update();
    canvasStandard->SaveAs(nameOutput.Data());
    delete canvasStandard;
}

					
void  PlotDPhiDistribution( TString optMC = "GammaConvV1_138.root", TString cutSel = "5010001_00216609297002008250400000_01525065000000", const char *suffix = "eps", TString optEnergy="PbPb_2.76TeV", TString optMCGenerator="Hijing", TString optPeriod="LHC11h", Int_t mode = 0){	
    cout << optMC.Data() << "\t" << cutSel.Data() << endl;
    gROOT->Reset();	
    gSystem->Load("libCore.so");
    gSystem->Load("libTree.so");
    gSystem->Load("libGeom.so");
    gSystem->Load("libVMC.so");
    gSystem->Load("libPhysics.so");
    gSystem->Load("libSTEERBase");
    gSystem->Load("libESD");
    gSystem->Load("libAOD");
    gSystem->Load("libANALYSIS");
    gSystem->Load("libANALYSISalice");
    gSystem->Load("libCORRFW.so");
    gROOT->SetStyle("Plain");
    
    StyleSettingsThesis();	
    //StyleSettings();
    
    SetPlotStyle();
    
    collisionSystem = ReturnFullCollisionsSystem(optEnergy);
    TString centralityString = GetCentralityString(cutSel);
    if (centralityString.CompareTo("pp")!=0 && !centralityString.Contains("0-100%") ){
        collisionSystem = Form("%s %s", centralityString.Data(), collisionSystem.Data());	
    }	

    
    TString outputDirectory;
    TString outputDirectory1;
    if(optPeriod.CompareTo("") ==0 || optPeriod.CompareTo("All") ==0){
        if (optMCGenerator.CompareTo("")!=0){
            outputDirectory     = Form("%s/%s/%s/%s/DPhiFakePhoton",optMCGenerator.Data(),cutSel.Data(),optEnergy.Data(),suffix);
            outputDirectory1    = Form("%s/%s",optMCGenerator.Data(),cutSel.Data());
        } else {
            outputDirectory     = Form("%s/%s/%s/DPhiFakePhoton",cutSel.Data(),optEnergy.Data(),suffix);
            outputDirectory1    = Form("%s",cutSel.Data());
        }
        gSystem->Exec("mkdir -p "+outputDirectory);
    } else if (optMCGenerator.CompareTo("")!=0){
        outputDirectory         = Form("%s/%s/%s/%s/%s/DPhiFakePhoton",optMCGenerator.Data(),cutSel.Data(),optEnergy.Data(),suffix,optPeriod.Data());
        outputDirectory1        = Form("%s/%s",optMCGenerator.Data(),cutSel.Data());
        gSystem->Exec("mkdir -p "+outputDirectory);
    } else {
        outputDirectory         = Form("%s/%s/%s/%s/DPhiFakePhoton",cutSel.Data(),optEnergy.Data(),suffix,optPeriod.Data());
        outputDirectory1        = Form("%s/%s",optMCGenerator.Data(),cutSel.Data());
        gSystem->Exec("mkdir -p "+outputDirectory);
    }

    Float_t floatLocationRightUp2D[4] = {0.42,0.95,0.035, 0.02};
    
// 	TString textStandardYAxis = "#frac{dN_{#gamma}}{N_{ch}}";
// 	TString textYAxisEtaElec = "#frac{dN_{e^{-}}}{N_{e^{-}} d#eta}";
// 	TString textYAxisEtaPos = "#frac{dN_{e^{+}}}{N_{e^{+}} d#eta}";
// 	TString textYAxisEtaGamma = "#frac{dN_{#gamma}}{N_{#gamma} d#eta}";
// 	TString textYAxisPtDist = "#frac{dN_{#gamma}}{N_{#gamma} dp_{t}} (GeV/c)^{-1}";
// 	TString textYAxisPtDistE = "#frac{dN_{e^{-}}}{N_{e^{-}} dp_{t}} (GeV/c)^{-1}";
// 	TString textYAxisPtDistP = "#frac{dN_{e^{+}}}{N_{e^{+}} dp_{t}} (GeV/c)^{-1}";
// 	TString textYAxisChi2 = "#frac{dN_{#gamma}}{N_{#gamma} d#chi^{2}}";
// 	TString textYAxisPsiPair = "#frac{dN_{#gamma}}{N_{#gamma} d#psi_{pair}}";
// 	TString textYAxisPhi = "#frac{dN_{#gamma}}{N_{#gamma} d#phi}";
// 	TString textYAxisMEP = "#frac{dN_{#gamma}}{N_{ch} dM_{e^{+}e^{-}}} (Gev/c^{2})^{-1}";
// 	TString textYAxisMGG = "#frac{dN_{#gamma}}{N_{ch} dM_{#gamma#gamma} } (Gev/c^{2})^{-1}";
// 	TString textYAxisAlpha = "#frac{dN_{#gamma}}{N_{ch} d#alpha}";
// 	TString textYAxisPoint = "#frac{dN_{#gamma}}{N_{ch} d cos(#theta_{pointing})}";
// 	TString textYAxisDCA = "#frac{dN_{#gamma}}{N_{ch} d DCA} (cm^{-1})";
// 	TString textYAxisNDCA = "#frac{1}{N_{ch}} #frac{dN_{#gamma}}{d DCA/#sigma}";
// 	TString textYAxisLike = "#frac{dN_{#gamma}}{N_{ch} d Likelyhood}";

    textDate                        = ReturnDateString();

// 	Float_t floatLineWidth = 1;
    
    // ---------------------------- LOAD ALL FILES ---------------------------------------------------
    TFile* fileMC                   = new TFile(fileMonteCarloInput);  
    TString autoDetectedMainDir     = AutoDetectMainTList(mode , fileMC);
    if (autoDetectedMainDir.CompareTo("") == 0){
        cout << "ERROR: trying to read file, which is incompatible with mode selected" << endl;;
        return;
    }

    //************************** Container Loading ********************************************************************
    TList *TopDirMonteCarlo         = (TList*)fileMC->Get(autoDetectedMainDir.Data());
    if(TopDirMonteCarlo == NULL){
        cout<<"ERROR: TopDirMonteCarlo not Found"<<endl;
        return;
    }

    TList *listHistosGammaConversion = (TList*)TopDirMonteCarlo->FindObject(Form("Cut Number %s",cutSel.Data()));
    TList *TrueConversionContainer   = (TList*)listHistosGammaConversion->FindObject(Form("%s True histograms",cutSel.Data()));
    
    TH2F* histdPhiTrueEK             = (TH2F*)TrueConversionContainer->FindObject("ESD_TrueCombinatorial_Pt_DeltaPhi_ek");
    TH2F* histdPhiTrueEP             = (TH2F*)TrueConversionContainer->FindObject("ESD_TrueCombinatorial_Pt_DeltaPhi_ep");
    TH2F* histdPhiTrueEPi            = (TH2F*)TrueConversionContainer->FindObject("ESD_TrueCombinatorial_Pt_DeltaPhi_epi");
    TH2F* histdPhiTruePiK            = (TH2F*)TrueConversionContainer->FindObject("ESD_TrueCombinatorial_Pt_DeltaPhi_pik");
    TH2F* histdPhiTruePiP            = (TH2F*)TrueConversionContainer->FindObject("ESD_TrueCombinatorial_Pt_DeltaPhi_pip");
        
// 	histdPhiTrueEK->GetZaxis()->SetRangeUser(1e-8,5e-4);
    if (histdPhiTrueEK->GetEntries() > 0){
        PlotStandard2D( histdPhiTrueEK , Form("%s/dPhi_fakePhoton_ek.%s",outputDirectory.Data(),suffix), "", "#it{p}_{#gamma_{cand},T} (GeV/#it{c})", "|#varphi_{K}-#varphi_{#gamma_{cand}}| (rad)",
                        kTRUE, 0., 1.5, kTRUE, 0.0, 4.,
                        0, 1, floatLocationRightUp2D,500,500,"MC", optPeriod, "dPhiEK");
    }
    if (histdPhiTrueEP->GetEntries() > 0){
        PlotStandard2D( histdPhiTrueEP , Form("%s/dPhi_fakePhoton_ep.%s",outputDirectory.Data(),suffix), "", "#it{p}_{#gamma_{cand},T} (GeV/#it{c})", "|#varphi_{p}-#varphi_{#gamma_{cand}}| (rad)",   
                        kTRUE, 0., 1.5, kTRUE, 0.0, 4.,
                        0, 1, floatLocationRightUp2D,500,500,"MC", optPeriod, "dPhiEP");
    }
    if (histdPhiTrueEPi->GetEntries() > 0){
        PlotStandard2D( histdPhiTrueEPi , Form("%s/dPhi_fakePhoton_epi.%s",outputDirectory.Data(),suffix), "", "#it{p}_{#gamma_{cand},T} (GeV/#it{c})", "|#varphi_{#pi}-#varphi_{#gamma_{cand}}| (rad)",  
                        kTRUE, 0., 1.5, kTRUE, 0.0, 4., 
                        0, 1, floatLocationRightUp2D,500,500,"MC", optPeriod, "dPhiEPi");
    }
    if (histdPhiTruePiK->GetEntries() > 0){
        PlotStandard2D( histdPhiTruePiK , Form("%s/dPhi_fakePhoton_pik.%s",outputDirectory.Data(),suffix), "", "#it{p}_{#gamma_{cand},T} (GeV/#it{c})", "|#varphi_{K}-#varphi_{#gamma_{cand}}| (rad)",  
                        kTRUE, 0., 1.5, kTRUE, 0.0, 4.,
                        0, 1, floatLocationRightUp2D,500,500,"MC", optPeriod, "dPhiPiK");
    }
    if (histdPhiTruePiP->GetEntries() > 0){
        PlotStandard2D( histdPhiTruePiP , Form("%s/dPhi_fakePhoton_pip.%s",outputDirectory.Data(),suffix), "", "#it{p}_{#gamma_{cand},T} (GeV/#it{c})", "|#varphi_{P}-#varphi_{#gamma_{cand}}| (rad)",  
                        kTRUE, 0., 1.5, kTRUE, 0.0, 4., 
                        0, 1, floatLocationRightUp2D,500,500,"MC", optPeriod, "dPhiPiP");
    }

}

