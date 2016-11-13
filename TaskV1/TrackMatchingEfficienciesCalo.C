/***********************************************************************************************
*** provided by Gamma Conversion Group, PWGGA,                                            ******
***     Friederike Bock, fbock@cern.ch                                                    ******
************************************************************************************************
************************************************************************************************/

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
#include "../CommonHeaders/ExtractSignalBinning.h"
// #include "../CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "../CommonHeaders/ConversionFunctions.h"

TString textGenerator;
TString collisionSystem;
TString textPeriod;
TString textDate;
TString mesonLatex;
TString detectionProcess;

//**********************************************************************************
//******************** Helper function for 2D histo plotting ***********************
//**********************************************************************************
void DrawAutoGammaHistoPaper2D( TH2* histo1,
                    TString Title, TString XTitle, TString YTitle,
                    Bool_t YRangeMax, Double_t YMaxFactor, Double_t YMinimum,
                    Bool_t YRange, Double_t YMin ,Double_t YMax,
                    Bool_t XRange, Double_t XMin, Double_t XMax, Double_t xOffset, Double_t yOffset) {
    if (YRangeMax && !XRange){
        YRange = kFALSE;
        Double_t maxRangeR = histo1->GetMaximum();
        Double_t minRangeR = histo1->GetMinimum();
        if(YMinimum > minRangeR){minRangeR = YMinimum;}
        histo1->GetYaxis()->SetRangeUser(minRangeR, maxRangeR*YMaxFactor);
    }
    if (YRangeMax && XRange){
        YRange = kFALSE;
        Double_t maxRangeR = histo1->GetMaximum();
        Double_t minRangeR = histo1->GetMinimum();
        if(YMinimum > minRangeR){minRangeR = YMinimum;}
        histo1->GetYaxis()->SetRangeUser(minRangeR, maxRangeR*YMaxFactor);
        histo1->GetXaxis()->SetRangeUser(XMin, XMax);
    }
    if (YRange && XRange){
        histo1->GetYaxis()->SetRangeUser(YMin, YMax);
        histo1->GetXaxis()->SetRangeUser(XMin, XMax);
    }
    if (!YRangeMax && !YRange && XRange){
        histo1->GetXaxis()->SetRangeUser(XMin, XMax);
    }

    if (YRange && !XRange){
        histo1->GetYaxis()->SetRangeUser(YMin, YMax);
    }

    histo1->SetTitle(Title.Data());

    if(XTitle.CompareTo("") != 0){
        histo1->SetXTitle(XTitle.Data());
    }
    if(YTitle.CompareTo("") != 0){
        histo1->SetYTitle(YTitle.Data());
    }

    histo1->GetYaxis()->SetLabelFont(42);
    histo1->GetXaxis()->SetLabelFont(42);
    histo1->GetYaxis()->SetTitleFont(62);
    histo1->GetXaxis()->SetTitleFont(62);
    histo1->GetYaxis()->SetLabelSize(0.045);
    histo1->GetYaxis()->SetTitleSize(0.05);
    histo1->GetYaxis()->SetDecimals();
    histo1->GetXaxis()->SetTitleOffset(xOffset);
    histo1->GetYaxis()->SetTitleOffset(yOffset);
    histo1->GetXaxis()->SetTitleSize(0.05);
    histo1->GetXaxis()->SetLabelSize(0.045);
}

//**********************************************************************************
//******************** Plotting 2D resolution matrix *******************************
//**********************************************************************************
void PlotStandard2D( TH2* histo2D, 
                     TString nameOutput,
                     TString title, 
                     TString xTitle, 
                     TString yTitle, 
                     Bool_t kRangeYAutoMax,
                     Double_t maxFacY,
                     Double_t startYAuto,
                     Bool_t kRangeY, 
                     Double_t startY, 
                     Double_t endY, 
                     Bool_t kRangeX, 
                     Double_t startX, 
                     Double_t endX,                
                     Double_t startZ,
                     Int_t logX, 
                     Int_t logY, 
                     Int_t logZ, 
                     Float_t* floatLogo, 
                     Double_t offsetX,
                     Double_t offsetY,
                     TString additionalLabel  = "",
                     Int_t canvasSizeX        = 500, 
                     Int_t canvasSizeY        = 500
                   ){
  
    TCanvas *canvasStandard = new TCanvas("canvasStandard","canvasStandard",canvasSizeX,canvasSizeY);
    DrawGammaCanvasSettings( canvasStandard, 0.115, 0.11, 0.02, 0.11); 
    canvasStandard->SetLogx(logX);
    canvasStandard->SetLogy(logY);
    canvasStandard->SetLogz(logZ);
    canvasStandard->cd();
    if (logX) histo2D->GetXaxis()->SetLabelOffset(-0.01);
    
    DrawAutoGammaHistoPaper2D(  histo2D,
                                title.Data(), xTitle.Data(), yTitle.Data(),
                                kRangeYAutoMax, maxFacY, startYAuto, 
                                kRangeY, startY, endY, 
                                kRangeX, startX, endX,
                                offsetX, offsetY                               
                             );
    histo2D->GetZaxis()->SetRangeUser(startZ, histo2D->GetMaximum());

    histo2D->SetTitle(title.Data());
    histo2D->Draw("colz");
    
    PutProcessLabelAndEnergyOnPlot(floatLogo[0], floatLogo[1], 28, collisionSystem.Data(), detectionProcess.Data(), additionalLabel.Data(), 63, 0.03);

    canvasStandard->Update();
    canvasStandard->SaveAs(nameOutput.Data());
    delete canvasStandard;
}


//**********************************************************************************
//******************* return minimum for 2 D histo  ********************************
//**********************************************************************************
Double_t FindSmallestEntryIn2D(TH2* histo){
    Double_t minimum = 1;
    for (Int_t i = 1; i<histo->GetNbinsX(); i++){
        for (Int_t j = 1; j<histo->GetNbinsY(); j++){
            if (histo->GetBinContent(i,j) < minimum && histo->GetBinContent(i,j) > 0){
                minimum = histo->GetBinContent(i,j);
            }
        }
    }
    return minimum;
}

//**********************************************************************************
//************* Main function for track matching efficiencies **********************
//**********************************************************************************
void TrackMatchingEfficienciesCalo(   TString fileMonteCarloInput          = "", 
                                      TString optionCutSelection           = "", 
                                      TString suffix                       = "eps", 
                                      TString optEnergy                    = "", 
                                      TString optMCGenerator               = "",
                                      TString optPeriod                    = "",
                                      Int_t mode                           = 0){ 
  
    gROOT->Reset(); 
    gROOT->SetStyle("Plain");
    
    StyleSettings();  
    SetPlotStyle();
    
    
    
    if (mode == 0 || mode == 1 || mode == 9) {
      cout << "ERROR: you are running in the PCM or PCM-Dalitz mode, this macro can't run in these modes, aborting." << endl;
      return;
    }
    InitializeClusterBinning(optEnergy, mode);
    
    collisionSystem                 = ReturnFullCollisionsSystem(optEnergy);
    detectionProcess                = ReturnFullTextReconstructionProcess(mode);
    if (collisionSystem.CompareTo("") == 0){
        cout << "No correct collision system specification, has been given" << endl;
        return;     
    }
    
    TString centrality              = "";
    TString firstCutnumber          = optionCutSelection(GetEventSystemCutPosition(),1);
    if (firstCutnumber.CompareTo("0") != 0){
        centrality                  = GetCentralityString(optionCutSelection);
        collisionSystem             = Form("%s %s", centrality.Data(), collisionSystem.Data());
    }   
    cout << centrality.Data() << endl;
    
    
    if(optMCGenerator.CompareTo("") ==0){
        textGenerator               = "";
    } else {
        textGenerator               = optMCGenerator;
    }
    textDate                        = ReturnDateString();

    //************************************ Separate cutstrings ***********************************
    TString fEventCutNumber         = "";
    TString fGammaCutNumber         = "";
    TString fClusterCutNumber       = "";
    TString fClusterMergedCutNumber = "";
    TString fMesonCutNumber         = "";
    TString dummyString             = "";
    if (mode != 10 && mode != 11 ){
        ReturnSeparatedCutNumberAdvanced( optionCutSelection,fEventCutNumber, fGammaCutNumber, fClusterCutNumber, dummyString, fMesonCutNumber, mode);
    } else {
        ReturnSeparatedCutNumberAdvanced( optionCutSelection, fEventCutNumber, fClusterCutNumber, fClusterMergedCutNumber, dummyString, fMesonCutNumber, mode);
    }
    TString trigger                             = fEventCutNumber(GetEventSelectSpecialTriggerCutPosition(),2);
    TString nameTrigger                         = ReturnTriggerName(trigger.Atoi());

    
    TString outputDirectory;
    if(optPeriod.CompareTo("") ==0 || optPeriod.CompareTo("All") ==0){
        textPeriod                  = "";
        outputDirectory             = Form("%s/%s/%s/TrackMatchingEfficienciesCalo",optionCutSelection.Data(),optEnergy.Data(),suffix.Data());
        gSystem->Exec("mkdir -p "+outputDirectory);
    } else {
        textPeriod                  = optPeriod;
        outputDirectory             = Form("%s/%s/%s/%s/TrackMatchingEfficienciesCalo",optionCutSelection.Data(),optEnergy.Data(),suffix.Data(),optPeriod.Data());
        gSystem->Exec("mkdir -p "+outputDirectory);
    }
    
    
    //**********************************************************************************************************************
    //******************************** Definition of some plotting variables ***********************************************
    //**********************************************************************************************************************

    //Array defintion for printing Logo 
    Float_t floatLocationRightDown2D[4]     = {0.15,0.25,0.11, 0.02};
    Float_t floatLocationUpDown2D[4]        = {0.15,0.95,0.11, 0.02};
    //**********************************************************************************************************************
    //****************************************** Loading of Histograms *****************************************************
    //**********************************************************************************************************************
    
    TFile fileMC(fileMonteCarloInput.Data());  
    
    //************************** Container Loading ********************************************************************
    TString nameOutputContainer     = "";
    if (mode == 2 || mode == 3) 
        nameOutputContainer         = "GammaConvCalo";
    else if (mode == 4 || mode == 5) 
        nameOutputContainer         = "GammaCalo";
    else if (mode == 10 || mode == 11) 
        nameOutputContainer         = "GammaCaloMerged";
      
    TList* TopDir                   = (TList*)fileMC.Get(nameOutputContainer.Data());
    if(TopDir == NULL){
        cout<<"ERROR: TopDir not Found"<<endl;
        return;
    }
    //***************************** Histo Loading ********************************************************************
    TList* HistosMainCut                = (TList*)TopDir->FindObject(Form("Cut Number %s",optionCutSelection.Data()));
    TList* CaloCutContainer             = (TList*)HistosMainCut->FindObject(Form("CaloCuts_%s",fClusterCutNumber.Data()));
    TH2F* histoTMEffiInput              = (TH2F*)CaloCutContainer->FindObject(Form("TMEffiInputHisto %s", fClusterCutNumber.Data()));
    TH2F* histoClE_TrE_ChCl             = (TH2F*)CaloCutContainer->FindObject(Form("ClusterE_TrackE_ChargedCluster %s", fClusterCutNumber.Data()));
    TH2F* histoClE_TrE_ChClLeadMatch    = (TH2F*)CaloCutContainer->FindObject(Form("ClusterE_TrackE_ChargedCluster_LeadMatched %s", fClusterCutNumber.Data()));
    TH2F* histoClE_TrE_NeCl             = (TH2F*)CaloCutContainer->FindObject(Form("ClusterE_TrackE_NeutralCluster %s", fClusterCutNumber.Data()));
    TH2F* histoClE_TrE_NeClSubCh        = (TH2F*)CaloCutContainer->FindObject(Form("ClusterE_TrackE_NeutralClusterSubCharged %s", fClusterCutNumber.Data()));
    TH2F* histoClE_TrE_GaCl             = (TH2F*)CaloCutContainer->FindObject(Form("ClusterE_TrackE_GammaCluster %s", fClusterCutNumber.Data()));
    TH2F* histoClE_TrE_GaClSubCh        = (TH2F*)CaloCutContainer->FindObject(Form("ClusterE_TrackE_GammaClusterSubCharged %s", fClusterCutNumber.Data()));
    TH2F* histoClE_TrE_ConvCl           = (TH2F*)CaloCutContainer->FindObject(Form("ClusterE_TrackE_ConvCluster %s", fClusterCutNumber.Data()));
    TH2F* histoClE_NMatch_NeCl          = (TH2F*)CaloCutContainer->FindObject(Form("ClusterE_NMatches_NeutralCluster %s", fClusterCutNumber.Data()));
    TH2F* histoClE_NMatch_ChCl          = (TH2F*)CaloCutContainer->FindObject(Form("ClusterE_NMatches_ChargedCluster %s", fClusterCutNumber.Data()));
    
    // ********************************************************************************************************************
    // ****************************** Plotting 2D distributions ***********************************************************
    // ********************************************************************************************************************
    PlotStandard2D( histoClE_TrE_ChCl, Form("%s/ClusterE_TrackE_ChargedClusters.%s",outputDirectory.Data(),suffix.Data()), "", 
                    "#it{E}_{cl} (GeV)", "#it{E}_{tr} (GeV)", 
                    0, 0, 0,
                    0, -10, 10, 
                    0, 0, 50, FindSmallestEntryIn2D(histoClE_TrE_ChCl),
                    1, 1, 1, floatLocationUpDown2D, 1, 1., "lead. part. in cl.: charged ", 1000, 800);
     PlotStandard2D( histoClE_TrE_ChClLeadMatch, Form("%s/ClusterE_TrackE_ChargedClusters_MatchedLeadParticle.%s",outputDirectory.Data(),suffix.Data()), "", 
                    "#it{E}_{cl} (GeV)", "#it{E}_{tr} (GeV)", 
                    0, 0, 0,
                    0, -10, 10, 
                    0, 0, 50, FindSmallestEntryIn2D(histoClE_TrE_ChClLeadMatch),
                    1, 1, 1, floatLocationUpDown2D, 1, 1., "lead. part. in cl.: charged, matched", 1000, 800);
     PlotStandard2D( histoClE_TrE_NeCl, Form("%s/ClusterE_TrackE_NeutralClusters.%s",outputDirectory.Data(),suffix.Data()), "", 
                    "#it{E}_{cl} (GeV)", "#it{E}_{tr} (GeV)", 
                    0, 0, 0,
                    0, -10, 10, 
                    0, 0, 50, FindSmallestEntryIn2D(histoClE_TrE_NeCl),
                    1, 1, 1, floatLocationUpDown2D, 1, 1., "lead. part. in cl.: neutral", 1000, 800);
     PlotStandard2D( histoClE_TrE_NeClSubCh, Form("%s/ClusterE_TrackE_NeutralClustersSubCharged.%s",outputDirectory.Data(),suffix.Data()), "", 
                    "#it{E}_{cl} (GeV)", "#it{E}_{tr} (GeV)", 
                    0, 0, 0,
                    0, -10, 10, 
                    0, 0, 50, FindSmallestEntryIn2D(histoClE_TrE_NeClSubCh),
                    1, 1, 1, floatLocationUpDown2D, 1, 1., "lead. part. in cl.: neutral, sub charged", 1000, 800);
     PlotStandard2D( histoClE_TrE_GaCl, Form("%s/ClusterE_TrackE_GammaClusters.%s",outputDirectory.Data(),suffix.Data()), "", 
                    "#it{E}_{cl} (GeV)", "#it{E}_{tr} (GeV)", 
                    0, 0, 0,
                    0, -10, 10, 
                    0, 0, 50, FindSmallestEntryIn2D(histoClE_TrE_GaCl),
                    1, 1, 1, floatLocationUpDown2D, 1, 1., "lead. part. in cl.: #gamma", 1000, 800);
     PlotStandard2D( histoClE_TrE_GaClSubCh, Form("%s/ClusterE_TrackE_GammaClustersSubCharged.%s",outputDirectory.Data(),suffix.Data()), "", 
                    "#it{E}_{cl} (GeV)", "#it{E}_{tr} (GeV)", 
                    0, 0, 0,
                    0, -10, 10, 
                    0, 0, 50, FindSmallestEntryIn2D(histoClE_TrE_GaClSubCh),
                    1, 1, 1, floatLocationUpDown2D, 1, 1., "lead. part. in cl.: #gamma, sub charged", 1000, 800);
     PlotStandard2D( histoClE_TrE_ConvCl, Form("%s/ClusterE_TrackE_ConvClusters.%s",outputDirectory.Data(),suffix.Data()), "", 
                    "#it{E}_{cl} (GeV)", "#it{E}_{tr} (GeV)", 
                    0, 0, 0,
                    0, -10, 10, 
                    0, 0, 50, FindSmallestEntryIn2D(histoClE_TrE_ConvCl),
                    1, 1, 1, floatLocationUpDown2D, 1, 1., "lead. part. in cl.: #gamma_{conv}", 1000, 800);
     PlotStandard2D( histoClE_NMatch_ChCl, Form("%s/ClusterE_NMatches_ChargedClusters.%s",outputDirectory.Data(),suffix.Data()), "", 
                    "#it{E}_{cl} (GeV)", "#it{N}_{matches}", 
                    0, 0, 0,
                    0, -10, 10, 
                    0, 0, 50, FindSmallestEntryIn2D(histoClE_NMatch_ChCl),
                    1, 0, 1, floatLocationUpDown2D, 1, 1., "lead. part. in cl.: charged", 1000, 800);
     PlotStandard2D( histoClE_NMatch_NeCl, Form("%s/ClusterE_NMatches_NeutralClusters.%s",outputDirectory.Data(),suffix.Data()), "", 
                    "#it{E}_{cl} (GeV)", "#it{N}_{matches}", 
                    0, 0, 0,
                    0, -10, 10, 
                    0, 0, 50, FindSmallestEntryIn2D(histoClE_NMatch_NeCl),
                    1, 0, 1, floatLocationUpDown2D, 1, 1., "lead. part. in cl.: neutral", 1000, 800);
    
    // ********************************************************************************************************************
    // ********************************* projections for efficiency calculations ******************************************
    // ********************************************************************************************************************
    histoTMEffiInput->Sumw2();

    TH1D*   fDeltaPtCluster       = new TH1D("fDeltaPtCluster","",fNBinsClusterPt,fBinsClusterPt);
    for(Int_t iPt=1;iPt<fNBinsClusterPt+1;iPt++){
        fDeltaPtCluster->SetBinContent(iPt,fBinsClusterPt[iPt]-fBinsClusterPt[iPt-1]);
        fDeltaPtCluster->SetBinError(iPt,0);
    }
    
    // *********** all inputs split in classes
    TH1D* histoClE_AllCl            = (TH1D*) histoTMEffiInput->ProjectionX("histoClE_AllCl",1,1,"e");
    TH1D* histoClE_AllClReb         = (TH1D*)histoClE_AllCl->Rebin(fNBinsClusterPt,"histoClE_AllClReb",fBinsClusterPt);
    histoClE_AllClReb->Divide(fDeltaPtCluster);
    TH1D* histoClE_ChCl             = (TH1D*) histoTMEffiInput->ProjectionX("histoClE_ChCl",2,2,"e");
    TH1D* histoClE_ChClReb          = (TH1D*)histoClE_ChCl->Rebin(fNBinsClusterPt,"histoClE_ChClReb",fBinsClusterPt);
    histoClE_ChClReb->Divide(fDeltaPtCluster);
    TH1D* histoClE_NeCl             = (TH1D*) histoTMEffiInput->ProjectionX("histoClE_NeCl",3,3,"e");
    TH1D* histoClE_NeClReb          = (TH1D*)histoClE_NeCl->Rebin(fNBinsClusterPt,"histoClE_NeClReb",fBinsClusterPt);
    histoClE_NeClReb->Divide(fDeltaPtCluster);
    TH1D* histoClE_NeClSubCh        = (TH1D*) histoTMEffiInput->ProjectionX("histoClE_NeClSubCh",4,4,"e");
    TH1D* histoClE_NeClSubChReb     = (TH1D*)histoClE_NeClSubCh->Rebin(fNBinsClusterPt,"histoClE_NeClSubChReb",fBinsClusterPt);
    histoClE_NeClSubChReb->Divide(fDeltaPtCluster);
    TH1D* histoClE_GaCl             = (TH1D*) histoTMEffiInput->ProjectionX("histoClE_GaCl",5,5,"e");
    TH1D* histoClE_GaClReb          = (TH1D*)histoClE_GaCl->Rebin(fNBinsClusterPt,"histoClE_GaClReb",fBinsClusterPt);
    histoClE_GaClReb->Divide(fDeltaPtCluster);
    TH1D* histoClE_GaClSubCh        = (TH1D*) histoTMEffiInput->ProjectionX("histoClE_GaClSubCh",6,6,"e");
    TH1D* histoClE_GaClSubChReb     = (TH1D*)histoClE_GaClSubCh->Rebin(fNBinsClusterPt,"histoClE_GaClSubChReb",fBinsClusterPt);
    histoClE_GaClSubChReb->Divide(fDeltaPtCluster);
    TH1D* histoClE_ConvCl           = (TH1D*) histoTMEffiInput->ProjectionX("histoClE_ConvCl",7,7,"e");
    TH1D* histoClE_ConvClReb        = (TH1D*)histoClE_ConvCl->Rebin(fNBinsClusterPt,"histoClE_ConvClReb",fBinsClusterPt);
    histoClE_ConvClReb->Divide(fDeltaPtCluster);
    TH1D* histoClE_ChPrimCl         = (TH1D*) histoTMEffiInput->ProjectionX("histoClE_ChPrimCl",8,8,"e");
    TH1D* histoClE_ChPrimClReb      = (TH1D*)histoClE_ChPrimCl->Rebin(fNBinsClusterPt,"histoClE_ChPrimClReb",fBinsClusterPt);
    histoClE_ChPrimClReb->Divide(fDeltaPtCluster);
    TH1D* histoClE_ElCl             = (TH1D*) histoTMEffiInput->ProjectionX("histoClE_ElCl",9,9,"e");
    TH1D* histoClE_ElClReb          = (TH1D*)histoClE_ElCl->Rebin(fNBinsClusterPt,"histoClE_ElClReb",fBinsClusterPt);
    histoClE_ElClReb->Divide(fDeltaPtCluster);
    
    // *********** all inputs split in classes with matches
    TH1D* histoClE_AllCl_mat        = (TH1D*) histoTMEffiInput->ProjectionX("histoClE_AllCl_mat",10,10,"e");
    TH1D* histoClE_AllCl_matReb     = (TH1D*)histoClE_AllCl_mat->Rebin(fNBinsClusterPt,"histoClE_AllCl_matReb",fBinsClusterPt);
    histoClE_AllCl_matReb->Divide(fDeltaPtCluster);
    TH1D* histoClE_ChCl_mat         = (TH1D*) histoTMEffiInput->ProjectionX("histoClE_ChCl_mat",11,11,"e");
    TH1D* histoClE_ChCl_matReb      = (TH1D*)histoClE_ChCl_mat->Rebin(fNBinsClusterPt,"histoClE_ChCl_matReb",fBinsClusterPt);
    histoClE_ChCl_matReb->Divide(fDeltaPtCluster);
    TH1D* histoClE_ChCl_leadmat     = (TH1D*) histoTMEffiInput->ProjectionX("histoClE_ChCl_leadmat",12,12,"e");
    TH1D* histoClE_ChCl_leadmatReb  = (TH1D*)histoClE_ChCl_leadmat->Rebin(fNBinsClusterPt,"histoClE_ChCl_leadmatReb",fBinsClusterPt);
    histoClE_ChCl_leadmatReb->Divide(fDeltaPtCluster);
    TH1D* histoClE_NeCl_mat         = (TH1D*) histoTMEffiInput->ProjectionX("histoClE_NeCl_mat",13,13,"e");
    TH1D* histoClE_NeCl_matReb      = (TH1D*)histoClE_NeCl_mat->Rebin(fNBinsClusterPt,"histoClE_NeCl_matReb",fBinsClusterPt);
    histoClE_NeCl_matReb->Divide(fDeltaPtCluster);
    TH1D* histoClE_NeClSubCh_mat    = (TH1D*) histoTMEffiInput->ProjectionX("histoClE_NeClSubCh_mat",14,14,"e");
    TH1D* histoClE_NeClSubCh_matReb = (TH1D*)histoClE_NeClSubCh_mat->Rebin(fNBinsClusterPt,"histoClE_NeClSubCh_matReb",fBinsClusterPt);
    histoClE_NeClSubCh_matReb->Divide(fDeltaPtCluster);
    TH1D* histoClE_GaCl_mat         = (TH1D*) histoTMEffiInput->ProjectionX("histoClE_GaCl_mat",15,15,"e");
    TH1D* histoClE_GaCl_matReb      = (TH1D*)histoClE_GaCl_mat->Rebin(fNBinsClusterPt,"histoClE_GaCl_matReb",fBinsClusterPt);
    histoClE_GaCl_matReb->Divide(fDeltaPtCluster);
    TH1D* histoClE_GaClSubCh_mat    = (TH1D*) histoTMEffiInput->ProjectionX("histoClE_GaClSubCh_mat",16,16,"e");
    TH1D* histoClE_GaClSubCh_matReb = (TH1D*)histoClE_GaClSubCh_mat->Rebin(fNBinsClusterPt,"histoClE_GaClSubCh_matReb",fBinsClusterPt);
    histoClE_GaClSubCh_matReb->Divide(fDeltaPtCluster);
    TH1D* histoClE_ConvCl_mat       = (TH1D*) histoTMEffiInput->ProjectionX("histoClE_ConvCl_mat",17,17,"e");
    TH1D* histoClE_ConvCl_matReb    = (TH1D*)histoClE_ConvCl_mat->Rebin(fNBinsClusterPt,"histoClE_ConvCl_matReb",fBinsClusterPt);
    histoClE_ConvCl_matReb->Divide(fDeltaPtCluster);
    TH1D* histoClE_ConvCl_leadmat   = (TH1D*) histoTMEffiInput->ProjectionX("histoClE_ConvCl_leadmat",18,18,"e");
    TH1D* histoClE_ConvCl_leadmatReb= (TH1D*)histoClE_ConvCl_leadmat->Rebin(fNBinsClusterPt,"histoClE_ConvCl_leadmatReb",fBinsClusterPt);
    histoClE_ConvCl_leadmatReb->Divide(fDeltaPtCluster);
    TH1D* histoClE_ChPrimCl_mat     = (TH1D*) histoTMEffiInput->ProjectionX("histoClE_ChPrimCl_mat",19,19,"e");
    TH1D* histoClE_ChPrimCl_matReb  = (TH1D*)histoClE_ChPrimCl_mat->Rebin(fNBinsClusterPt,"histoClE_ChPrimCl_matReb",fBinsClusterPt);
    histoClE_ChPrimCl_matReb->Divide(fDeltaPtCluster);
    TH1D* histoClE_ChPrimCl_leadmat = (TH1D*) histoTMEffiInput->ProjectionX("histoClE_ChPrimCl_leadmat",20,20,"e");
    TH1D* histoClE_ChPrimCl_leadmatReb  = (TH1D*)histoClE_ChPrimCl_leadmat->Rebin(fNBinsClusterPt,"histoClE_ChPrimCl_leadmatReb",fBinsClusterPt);
    histoClE_ChPrimCl_leadmatReb->Divide(fDeltaPtCluster);
    TH1D* histoClE_ElCl_mat         = (TH1D*) histoTMEffiInput->ProjectionX("histoClE_ElCl_mat",21,21,"e");
    TH1D* histoClE_ElCl_matReb      = (TH1D*)histoClE_ElCl_mat->Rebin(fNBinsClusterPt,"histoClE_ElCl_matReb",fBinsClusterPt);
    histoClE_ElCl_matReb->Divide(fDeltaPtCluster);
    TH1D* histoClE_ElCl_leadmat     = (TH1D*) histoTMEffiInput->ProjectionX("histoClE_ElCl_leadmat",22,22,"e");
    TH1D* histoClE_ElCl_leadmatReb  = (TH1D*)histoClE_ElCl_leadmat->Rebin(fNBinsClusterPt,"histoClE_ElCl_leadmatReb",fBinsClusterPt);
    histoClE_ElCl_leadmatReb->Divide(fDeltaPtCluster);
    
    TH1D* histoEffi_AllCl           = (TH1D*)histoClE_AllCl_mat->Clone("histoEffi_AllCl");
    histoEffi_AllCl->Divide(histoEffi_AllCl, histoClE_AllCl, 1, 1, "B");
    TH1D* histoEffi_ChCl            = (TH1D*)histoClE_ChCl_mat->Clone("histoEffi_ChCl");
    histoEffi_ChCl->Divide(histoEffi_ChCl, histoClE_ChCl, 1, 1, "B");
    TH1D* histoEffi_ChPrimCl        = (TH1D*)histoClE_ChPrimCl_mat->Clone("histoEffi_ChPrimCl");
    histoEffi_ChPrimCl->Divide(histoEffi_ChPrimCl, histoClE_ChPrimCl, 1, 1, "B");
    TH1D* histoEffi_NeCl            = (TH1D*)histoClE_NeCl_mat->Clone("histoEffi_NeCl");
    histoEffi_NeCl->Divide(histoEffi_NeCl, histoClE_NeCl, 1, 1, "B");
    TH1D* histoEffi_NeSubChCl       = (TH1D*)histoClE_NeClSubCh_mat->Clone("histoEffi_NeSubChCl");
    histoEffi_NeSubChCl->Divide(histoEffi_NeSubChCl, histoClE_NeClSubCh, 1, 1, "B");
    TH1D* histoEffi_GaCl            = (TH1D*)histoClE_GaCl_mat->Clone("histoEffi_GaCl");
    histoEffi_GaCl->Divide(histoEffi_GaCl, histoClE_GaCl, 1, 1, "B");
    TH1D* histoEffi_GaSubChCl       = (TH1D*)histoClE_GaClSubCh_mat->Clone("histoEffi_GaSubChCl");
    histoEffi_GaSubChCl->Divide(histoEffi_GaSubChCl, histoClE_GaClSubCh, 1, 1, "B");
    TH1D* histoEffi_ConvCl          = (TH1D*)histoClE_ConvCl_mat->Clone("histoEffi_ConvCl");
    histoEffi_ConvCl->Divide(histoEffi_ConvCl, histoClE_ConvCl, 1, 1, "B");
    TH1D* histoEffi_ElCl            = (TH1D*)histoClE_ElCl_mat->Clone("histoEffi_ConvCl");
    histoEffi_ElCl->Divide(histoEffi_ElCl, histoClE_ElCl, 1, 1, "B");
    
    TH1D* histoEffi_AllClReb        = (TH1D*)histoClE_AllCl_matReb->Clone("histoEffi_AllClReb");
    histoEffi_AllClReb->Divide(histoEffi_AllClReb, histoClE_AllClReb, 1, 1, "B");
    TH1D* histoEffi_ChClReb         = (TH1D*)histoClE_ChCl_matReb->Clone("histoEffi_ChClReb");
    histoEffi_ChClReb->Divide(histoEffi_ChClReb, histoClE_ChClReb, 1, 1, "B");
    TH1D* histoEffi_ChPrimClReb     = (TH1D*)histoClE_ChPrimCl_matReb->Clone("histoEffi_ChPrimClReb");
    histoEffi_ChPrimClReb->Divide(histoEffi_ChPrimClReb, histoClE_ChPrimClReb, 1, 1, "B");
    TH1D* histoEffi_ChCl_leadReb    = (TH1D*)histoClE_ChCl_leadmatReb->Clone("histoEffi_ChClLeadReb");
    histoEffi_ChCl_leadReb->Divide(histoEffi_ChCl_leadReb, histoClE_ChClReb, 1, 1, "B");
    TH1D* histoEffi_ChPrimCl_leadReb= (TH1D*)histoClE_ChPrimCl_leadmatReb->Clone("histoEffi_ChPrimClLeadReb");
    histoEffi_ChPrimCl_leadReb->Divide(histoEffi_ChPrimCl_leadReb, histoClE_ChPrimClReb, 1, 1, "B");

    TH1D* histoEffi_NeClReb         = (TH1D*)histoClE_NeCl_matReb->Clone("histoEffi_NeClReb");
    histoEffi_NeClReb->Divide(histoEffi_NeClReb, histoClE_NeClReb, 1, 1, "B");
    TH1D* histoEffi_NeSubChClReb    = (TH1D*)histoClE_NeClSubCh_matReb->Clone("histoEffi_NeSubChClReb");
    histoEffi_NeSubChClReb->Divide(histoEffi_NeSubChClReb, histoClE_NeClSubChReb, 1, 1, "B");
    TH1D* histoEffi_GaClReb         = (TH1D*)histoClE_GaCl_matReb->Clone("histoEffi_GaClReb");
    histoEffi_GaClReb->Divide(histoEffi_GaClReb, histoClE_GaClReb, 1, 1, "B");
    TH1D* histoEffi_GaSubChClReb    = (TH1D*)histoClE_GaClSubCh_matReb->Clone("histoEffi_GaSubChClReb");
    histoEffi_GaSubChClReb->Divide(histoEffi_GaSubChClReb, histoClE_GaClSubChReb, 1, 1, "B");
    TH1D* histoEffi_ConvClReb       = (TH1D*)histoClE_ConvCl_matReb->Clone("histoEffi_ConvCl_leadReb");
    histoEffi_ConvClReb->Divide(histoEffi_ConvClReb, histoClE_ConvClReb, 1, 1, "B");
    TH1D* histoEffi_ConvCl_leadReb  = (TH1D*)histoClE_ConvCl_leadmatReb->Clone("histoEffi_ConvClReb");
    histoEffi_ConvCl_leadReb->Divide(histoEffi_ConvCl_leadReb, histoClE_ConvClReb, 1, 1, "B");
    TH1D* histoEffi_ElClReb         = (TH1D*)histoClE_ElCl_matReb->Clone("histoEffi_ConvClReb");
    histoEffi_ElClReb->Divide(histoEffi_ElClReb, histoClE_ElClReb, 1, 1, "B");
    TH1D* histoEffi_ElCl_leadReb    = (TH1D*)histoClE_ElCl_leadmatReb->Clone("histoEffi_ElCl_leadReb");
    histoEffi_ElCl_leadReb->Divide(histoEffi_ElCl_leadReb, histoClE_ElClReb, 1, 1, "B");
    
    // Plot TM efficiency for different cluster type
    TCanvas* canvasTMEffi = new TCanvas("canvasTMEffi","",200,10,1100,900);  // gives the page size
    DrawGammaCanvasSettings( canvasTMEffi, 0.075, 0.02, 0.02, 0.09);

        Double_t rangeTMEffi[2]     = {0., 1.02};

        Double_t maxPt              = fBinsClusterPt[fNBinsClusterPt];
        if (nameTrigger.Contains("INT") || nameTrigger.Contains("MB")){
            maxPt                   = 15;
            rangeTMEffi[1]          = 1.3;
        }    
        TH2F * histo2DTMEff         = new TH2F("histo2DTMEff", "histo2DTMEff",1000, 0., maxPt , 1000, rangeTMEffi[0], rangeTMEffi[1] );
        SetStyleHistoTH2ForGraphs(  histo2DTMEff, "#it{p}_{T} (GeV/#it{c})", "#varepsilon_{TM}", 
                                    0.85*0.04, 0.04, 0.85*0.04, 0.04, 0.95, 0.94);
        histo2DTMEff->GetYaxis()->SetLabelOffset(0.001);
        histo2DTMEff->DrawCopy(); 
        
        DrawGammaSetMarker(histoEffi_ChClReb, 24, 1.5, kAzure+2, kAzure+2);
        histoEffi_ChClReb->DrawCopy("same,e1"); 
        DrawGammaSetMarker(histoEffi_ElClReb, 21, 1.5, kGreen-2, kGreen-2);
        histoEffi_ElClReb->DrawCopy("same,e1"); 
        DrawGammaSetMarker(histoEffi_ConvClReb, 25, 1.5, kGreen+2, kGreen+2);
        histoEffi_ConvClReb->DrawCopy("same,e1"); 
        DrawGammaSetMarker(histoEffi_ChPrimClReb, 24, 1.5, kBlue+1, kBlue+1);
        histoEffi_ChPrimClReb->DrawCopy("same,e1"); 

        TLegend* legendTMEffiCh           = GetAndSetLegend2(0.42, 0.95, 0.7, 0.95-(4*0.04),32);
        legendTMEffiCh->AddEntry(histoEffi_ConvClReb,"conv #gamma cl","p");
        legendTMEffiCh->AddEntry(histoEffi_ElClReb,"e^{#pm} cl","p");
        legendTMEffiCh->AddEntry(histoEffi_ChPrimClReb,"other charged cl., prod. vtx R #leq 5 cm","p");
        legendTMEffiCh->AddEntry(histoEffi_ChClReb,"other charged cl., prod. vtx R > 5 cm","p");
        legendTMEffiCh->Draw();

        PutProcessLabelAndEnergyOnPlot(0.65, 0.20, 32, collisionSystem.Data(), detectionProcess.Data(),"",  63, 0.029);

    canvasTMEffi->Update();
    canvasTMEffi->SaveAs(Form("%s/TMEffiDifferentChargedClusters_%s.%s",outputDirectory.Data(),optionCutSelection.Data(),suffix.Data()));

        histo2DTMEff->DrawCopy(); 
        
        DrawGammaSetMarker(histoEffi_ChCl_leadReb, 24, 1.5, kAzure+2, kAzure+2);
        histoEffi_ChCl_leadReb->DrawCopy("same,e1"); 
        DrawGammaSetMarker(histoEffi_ElCl_leadReb, 21, 1.5, kGreen-2, kGreen-2);
        histoEffi_ElCl_leadReb->DrawCopy("same,e1"); 
        DrawGammaSetMarker(histoEffi_ConvCl_leadReb, 25, 1.5, kGreen+2, kGreen+2);
        histoEffi_ConvCl_leadReb->DrawCopy("same,e1"); 
        DrawGammaSetMarker(histoEffi_ChPrimCl_leadReb, 24, 1.5, kBlue+1, kBlue+1);
        histoEffi_ChPrimCl_leadReb->DrawCopy("same,e1"); 

        TLegend* legendTMEffiChlead           = GetAndSetLegend2(0.42, 0.95, 0.7, 0.95-(4*0.04),32);
        legendTMEffiChlead->AddEntry(histoEffi_ConvCl_leadReb,"conv #gamma cl","p");
        legendTMEffiChlead->AddEntry(histoEffi_ElCl_leadReb,"e^{#pm} cl","p");
        legendTMEffiChlead->AddEntry(histoEffi_ChPrimCl_leadReb,"other charged cl., prod. vtx R #leq 5 cm","p");
        legendTMEffiChlead->AddEntry(histoEffi_ChCl_leadReb,"other charged cl., prod. vtx R > 5 cm","p");
        legendTMEffiChlead->Draw();

        PutProcessLabelAndEnergyOnPlot(0.105, 0.95, 32, collisionSystem.Data(), detectionProcess.Data(),"matched w/ lead h^{#pm} of cl.",  63, 0.029);

    canvasTMEffi->Update();
    canvasTMEffi->SaveAs(Form("%s/TMEffiDifferentChargedClustersMatchedWithLead_%s.%s",outputDirectory.Data(),optionCutSelection.Data(),suffix.Data()));
    
    
        TH2F * histo2DTMEffNe        = new TH2F("histo2DTMEffNe", "histo2DTMEffNe",1000, 0., maxPt , 1000, rangeTMEffi[0], rangeTMEffi[1] );
        SetStyleHistoTH2ForGraphs(  histo2DTMEffNe, "#it{p}_{T} (GeV/#it{c})", "#varepsilon_{miss,TM}", 
                                    0.85*0.04, 0.04, 0.85*0.04, 0.04, 0.95, 0.94);
        histo2DTMEffNe->GetYaxis()->SetLabelOffset(0.001);
        histo2DTMEffNe->DrawCopy(); 
//         DrawGammaSetMarker(histoEffi_NeSubChClReb, 20, 1.5, kRed-6, kRed-6);
//         histoEffi_NeSubChClReb->DrawCopy("same,e1"); 
        DrawGammaSetMarker(histoEffi_NeClReb, 24, 1.5, kRed+2, kRed+2);
        histoEffi_NeClReb->DrawCopy("same,e1"); 
        DrawGammaSetMarker(histoEffi_GaClReb, 25, 1.5, kOrange, kOrange);
        histoEffi_GaClReb->DrawCopy("same,e1"); 
        DrawGammaSetMarker(histoEffi_GaSubChClReb, 21, 1.5, 807, 807);
        histoEffi_GaSubChClReb->DrawCopy("same,e1"); 

        TLegend* legendTMEffiNe           = GetAndSetLegend2(0.42, 0.95, 0.7, 0.95-(3*0.04),32);
        legendTMEffiNe->AddEntry(histoEffi_GaClReb,"#gamma cl.","p");
        legendTMEffiNe->AddEntry(histoEffi_GaSubChClReb,"#gamma cl. w/ ch. cont","p");
        legendTMEffiNe->AddEntry(histoEffi_NeClReb,"other neutral cl.","p");
//         legendTMEffiNe->AddEntry(histoEffi_NeSubChClReb,"other neutral cl. w/ ch. cont.","p");
        legendTMEffiNe->Draw();

        PutProcessLabelAndEnergyOnPlot(0.105, 0.95, 32, collisionSystem.Data(), detectionProcess.Data(),"",  63, 0.029);

    canvasTMEffi->Update();
    canvasTMEffi->SaveAs(Form("%s/TMEffiDifferentNeutralClusters_%s.%s",outputDirectory.Data(),optionCutSelection.Data(),suffix.Data()));
    
    return;

}
