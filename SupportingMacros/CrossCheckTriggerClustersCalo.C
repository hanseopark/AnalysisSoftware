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
    DrawGammaCanvasSettings( canvasStandard, 0.115, 0.11, 0.05, 0.11); 
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
//************* Main function for cross check of trigger clusters ******************
//**********************************************************************************
void CrossCheckTriggerClustersCalo(   TString fileDataInput                = "", 
                                      TString optionCutSelection           = "", 
                                      TString suffix                       = "eps", 
                                      TString optEnergy                    = "", 
                                      TString optPeriod                    = "",
                                      Int_t mode                           = 0){ 
  
    //************************************* general  plotting settings ****************************
    gROOT->Reset(); 
    gROOT->SetStyle("Plain");
    
    StyleSettings();  
    SetPlotStyle();
    
    if (mode != 4 ) {
      cout << "ERROR: you not running with the EMC mode" << endl;
      return;
    }
    InitializeClusterBinning(optEnergy, mode);
    
    textDate                        = ReturnDateString();
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
    TString trigger                 = fEventCutNumber(GetEventSelectSpecialTriggerCutPosition(),2);
    TString nameTrigger             = ReturnTriggerName(trigger.Atoi());

    Double_t triggerEMin            = 0;
    if (optEnergy.CompareTo("2.76TeV") == 0){
        if (trigger.Atoi() == 51){
            triggerEMin             = 4.5;
        } else if (trigger.Atoi() == 52){
            triggerEMin             = 2.7;
        } else if (trigger.Atoi() == 85){
            triggerEMin             = 4.7;
        } else if (trigger.Atoi() == 83){    
            triggerEMin             = 6.5;
        }    
    }
    
    TString outputDirectory;
    if(optPeriod.CompareTo("") ==0 || optPeriod.CompareTo("All") ==0){
        textPeriod                  = "";
        outputDirectory             = Form("%s/%s/%s/CrossCheckTriggerClusterCalo",optionCutSelection.Data(),optEnergy.Data(),suffix.Data());
        gSystem->Exec("mkdir -p "+outputDirectory);
    } else {
        textPeriod                  = optPeriod;
        outputDirectory             = Form("%s/%s/%s/%s/CrossCheckTriggerClusterCalo",optionCutSelection.Data(),optEnergy.Data(),suffix.Data(),optPeriod.Data());
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
    
    TFile* fileData                     = new TFile(fileDataInput.Data());  
    TString autoDetectedMainDir         = AutoDetectMainTList(mode , fileData);
    if (autoDetectedMainDir.CompareTo("") == 0){
        cout << "ERROR: trying to read file, which is incompatible with mode selected" << endl;;
        return;
    }
    
    TList* TopDir                       = (TList*)fileData->Get(autoDetectedMainDir.Data());
    if(TopDir == NULL){
        cout<<"ERROR: TopDir not Found"<<endl;
        return;
    }
    //***************************** Histo Loading ********************************************************************
    TList* HistosMainCut                = (TList*)TopDir->FindObject(Form("Cut Number %s",optionCutSelection.Data()));
    
    TList* MainCutContainer             = (TList*)HistosMainCut->FindObject(Form("%s ESD histograms",optionCutSelection.Data()));
    TH2F* histoCellIDvsClusterEnergy    = (TH2F*)MainCutContainer->FindObject("CellIDvsClusterEnergy");
    
    TH2F* histoCellIDvsClusterEnergyMax = (TH2F*)MainCutContainer->FindObject("CellIDvsClusterEnergyMax");
    if (!histoCellIDvsClusterEnergy){
        cout << "INFO: current contained: " << optionCutSelection.Data() << " does not contain necessary inputs, abort." << endl;
        return;
    }
    histoCellIDvsClusterEnergy->Sumw2();
    histoCellIDvsClusterEnergyMax->Sumw2();
    TH1D* fEventQuality                 = (TH1D*)MainCutContainer->FindObject("NEvents");
    Double_t fNEvents                   = 0;
    if (optEnergy.Contains("PbPb") || optEnergy.Contains("pPb") ){
        fNEvents                        = fEventQuality->GetBinContent(1);
    } else {
        fNEvents                        =  GetNEvents(fEventQuality);
    }
    
    
    // ********************************************************************************************************************
    // ****************************** Plotting 2D distributions ***********************************************************
    // ********************************************************************************************************************
    PlotStandard2D( histoCellIDvsClusterEnergy, Form("%s/ClusterE_CellID_%s.png",outputDirectory.Data(),optionCutSelection.Data()), "", 
                    "#it{E}_{cl} (GeV)", "ID_{cell} (GeV)", 
                    0, 0, 0,
                    0, 0, 12000, 
                    0, 0, 50, FindSmallestEntryIn2D(histoCellIDvsClusterEnergy),
                    1, 0, 1, floatLocationUpDown2D, 1, 1., "cluster #it{E}_{cl} vs. lead cell ID", 1000, 800);
    PlotStandard2D( histoCellIDvsClusterEnergyMax, Form("%s/ClusterEMaxevent_CellID_%s.png",outputDirectory.Data(),optionCutSelection.Data()), "", 
                    "#it{E}_{cl} (GeV)", "ID_{cell} (GeV)", 
                    0, 0, 0,
                    0, 0, 12000, 
                    0, 0, 50, FindSmallestEntryIn2D(histoCellIDvsClusterEnergyMax),
                    1, 0, 1, floatLocationUpDown2D, 1, 1., "cluster #it{E}_{max,cl} in event vs. lead cell ID", 1000, 800);
    
    // ********************************************************************************************************************
    // ********************************* projections for efficiency calculations ******************************************
    // ********************************************************************************************************************
    histoCellIDvsClusterEnergy->Sumw2();

    Int_t maxChannel        = 12672;
    if ( optEnergy.CompareTo("2.76TeV") == 0 || optEnergy.CompareTo("8TeV") == 0 || optEnergy.CompareTo("pPb_5.023TeV") == 0 || optEnergy.CompareTo("PbPb_2.76TeV") == 0 )
        maxChannel = 288*4*10;
    
    TH1D*   fDeltaPtCluster       = new TH1D("fDeltaPtCluster","",fNBinsClusterPt,fBinsClusterPt);
    for(Int_t iPt=1;iPt<fNBinsClusterPt+1;iPt++){
        fDeltaPtCluster->SetBinContent(iPt,fBinsClusterPt[iPt]-fBinsClusterPt[iPt-1]);
        fDeltaPtCluster->SetBinError(iPt,0);
    }
    
    // *********** all inputs vs E
    TH1D* histoClE_AllCell          = (TH1D*)histoCellIDvsClusterEnergy->ProjectionX("histoClE_AllCell",0,14000,"e");
    TH1D* histoClE_AllCellReb       = (TH1D*)histoClE_AllCell->Rebin(fNBinsClusterPt,"histoClE_AllCellReb",fBinsClusterPt);
    histoClE_AllCellReb->Divide(fDeltaPtCluster);
    TH1D* histoClEMax_AllCell       = (TH1D*)histoCellIDvsClusterEnergyMax->ProjectionX("histoClEMax_AllCell",0,14000,"e");
    TH1D* histoClEMax_AllCellReb    = (TH1D*)histoClEMax_AllCell->Rebin(fNBinsClusterPt,"histoClEMax_AllCellReb",fBinsClusterPt);
    histoClEMax_AllCellReb->Divide(fDeltaPtCluster);
    // *********** Cell ID filled
    
    cout << nameTrigger.Data() << "\t" << fNEvents << "\t event in max E hist: "<< histoClEMax_AllCell->Integral(0,50) << endl;
    
    
    TH1D* histoAllCells             = (TH1D*)histoCellIDvsClusterEnergy->ProjectionY("histoAllCells",0,histoCellIDvsClusterEnergy->GetNbinsX(),"e");
    histoAllCells->Sumw2();
    histoAllCells->GetXaxis()->SetRange(1,maxChannel+1);
    Double_t maxFreq                = histoAllCells->GetMaximum()/fNEvents;
    TH1D* histoAllCells_Erange      = (TH1D*)histoCellIDvsClusterEnergy->ProjectionY("histoAllCells_Erange",histoCellIDvsClusterEnergy->GetXaxis()->FindBin(triggerEMin), histoCellIDvsClusterEnergy->GetNbinsX(),"e");
    histoAllCells_Erange->Sumw2();
    histoAllCells_Erange->GetXaxis()->SetRange(1,maxChannel+1);
    TH1D* histoBadChannel           = (TH1D*)histoAllCells->Clone("histoBadChannels");
    TH1D* histoBadChannel_Erange    = (TH1D*)histoAllCells_Erange->Clone("histoBadChannel_Erange");
    histoAllCells->Scale(1/fNEvents);
    histoAllCells_Erange->Scale(1/fNEvents);
    
    TH1D* histoAllCells_maxE        = (TH1D*)histoCellIDvsClusterEnergyMax->ProjectionY("histoAllCells_maxE",0,histoCellIDvsClusterEnergyMax->GetNbinsX(),"e");
    histoAllCells_maxE->Sumw2();
    histoAllCells_maxE->GetXaxis()->SetRange(1,maxChannel+1);
    TH1D* histoAllCells_maxE_Erange = (TH1D*)histoCellIDvsClusterEnergyMax->ProjectionY("histoAllCells_maxE_Erange",histoCellIDvsClusterEnergyMax->GetXaxis()->FindBin(triggerEMin),histoCellIDvsClusterEnergyMax->GetNbinsX(),"e");
    histoAllCells_maxE_Erange->Sumw2();
    histoAllCells_maxE_Erange->GetXaxis()->SetRange(1,maxChannel+1);
    Double_t minFreq                = 1/fNEvents;
    TH1D* histoBadChannelTrigg      = (TH1D*)histoAllCells_maxE->Clone("histoBadChannelsTrigg");
    TH1D* histoBadChannelTrigg_Erange   = (TH1D*)histoAllCells_maxE_Erange->Clone("histoBadChannelTrigg_Erange");
    histoAllCells_maxE->Scale(1/fNEvents);
    histoAllCells_maxE_Erange->Scale(1/fNEvents);

    Int_t nBadChanE         = 0;
    Int_t nBadChanEhigh     = 0;
    Int_t nBadChanTrigg     = 0;
    Int_t nBadChanTrigghigh = 0;    
    Int_t nBadNotOverlap    = 0;    
    Int_t nBadNotOverlaphigh= 0;       
    for (Int_t i = 1; i < maxChannel+1; i++){
        if (histoBadChannel->GetBinContent(i) > 0) 
            histoBadChannel->SetBinContent(i,1);
        else 
            nBadChanE++;
        if (histoBadChannel_Erange->GetBinContent(i) > 0) 
            histoBadChannel_Erange->SetBinContent(i,1);
        else 
            nBadChanEhigh++;        
        if (histoBadChannelTrigg->GetBinContent(i) > 0) 
            histoBadChannelTrigg->SetBinContent(i,1);
        else 
            nBadChanTrigg++;
        if (histoBadChannelTrigg_Erange->GetBinContent(i) > 0) 
            histoBadChannelTrigg_Erange->SetBinContent(i,1);
        else 
            nBadChanTrigghigh++;
    }
    TH1D* histoBadChannelOverlap    = (TH1D*)histoBadChannel->Clone("histoBadChannelOverlap");
    histoBadChannelOverlap->Add(histoBadChannelTrigg, -1);
    TH1D* histoBadChannelOverlapHigh= (TH1D*)histoBadChannel_Erange->Clone("histoBadChannelOverlapHigh");
    histoBadChannelOverlapHigh->Add(histoBadChannelTrigg_Erange, -1);
    
    for (Int_t i = 1; i < maxChannel+1; i++){
        if (histoBadChannelOverlap->GetBinContent(i) != 0) nBadNotOverlap++;
        if (histoBadChannelOverlapHigh->GetBinContent(i) != 0) nBadNotOverlaphigh++;
    }    
    cout << "n Bad channels without min E cut: " << nBadChanE << "\t only max E clus per event: " << nBadChanTrigg << endl;
    cout << "n Bad channels with min E cut of "<< triggerEMin << "GeV : " << nBadChanEhigh << "\t only max E clus per event: " << nBadChanTrigghigh << endl;
        
    cout << "bad channels diff all clus vs max clus per event: " << nBadNotOverlap << "\t" << (Double_t)nBadNotOverlap/maxChannel*100 << "%" << endl;
    cout << "bad channels diff all clus vs max clus per event with E>"<< triggerEMin<< "GeV : " << nBadNotOverlaphigh << "\t" << (Double_t)nBadNotOverlaphigh/maxChannel*100 << "%" << endl;
    
    // Plot TM efficiency for different cluster type
    TCanvas* canvasClusterE = new TCanvas("canvasClusterE","",200,10,1100,900);  // gives the page size
    DrawGammaCanvasSettings( canvasClusterE, 0.075, 0.02, 0.02, 0.09);
    canvasClusterE->SetLogy(1);
//     canvasClusterE->SetLogx(1);
    
        TH1F * histo1DEAllCell      = new TH1F("histo1DEAllCell", "histo1DEAllCell",1000, 0., fBinsClusterPt[fNBinsClusterPt] );
        SetStyleHistoTH1ForGraphs(  histo1DEAllCell, "#it{E}_{clus} (GeV)", "counts", 
                                    0.85*0.04, 0.04, 0.85*0.04, 0.04, 0.95, 0.94);
        histo1DEAllCell->GetYaxis()->SetRangeUser(histoClEMax_AllCell->GetMinimum(0.1)/10, histoClE_AllCell->GetMaximum()*10 );
        histo1DEAllCell->GetYaxis()->SetLabelOffset(0.001);
        histo1DEAllCell->DrawCopy(); 
        
        DrawGammaSetMarker(histoClEMax_AllCell, 21, 1.5, kGreen-2, kGreen-2);
        histoClEMax_AllCell->DrawCopy("same,e1"); 
        DrawGammaSetMarker(histoClE_AllCell, 24, 1.5, kAzure+2, kAzure+2);
        histoClE_AllCell->DrawCopy("same,e1"); 

        TLegend* legendClEAllCell   = GetAndSetLegend2(0.42, 0.95, 0.7, 0.95-(2*0.04),32);
        legendClEAllCell->AddEntry(histoClE_AllCell,"all clusters","p");
        legendClEAllCell->AddEntry(histoClEMax_AllCell,"only maximum cluster per event","p");
        legendClEAllCell->Draw();

        PutProcessLabelAndEnergyOnPlot(0.65, 0.95-(2*0.04), 32, collisionSystem.Data(), detectionProcess.Data(),"",  63, 0.029);

    canvasClusterE->Update();
    canvasClusterE->SaveAs(Form("%s/SpectraAllCells_%s.%s",outputDirectory.Data(),optionCutSelection.Data(),suffix.Data()));

    TCanvas* canvasClusterF = new TCanvas("canvasClusterF","",200,10,1600,900);  // gives the page size
    DrawGammaCanvasSettings( canvasClusterF, 0.06, 0.035, 0.02, 0.09);    
    canvasClusterF->SetLogy(1);
    canvasClusterF->cd();
    
    TH1F * histo1DFAllCell      = new TH1F("histo1DFAllCell", "histo1DFAllCell",maxChannel, 0., maxChannel );
    SetStyleHistoTH1ForGraphs(  histo1DFAllCell, "ID_{cell}", "#it{f}", 
                                0.85*0.04, 0.04, 0.85*0.04, 0.04, 0.95, 0.9);
    histo1DFAllCell->GetYaxis()->SetRangeUser(minFreq/10, maxFreq*10 );
    histo1DFAllCell->GetYaxis()->SetLabelOffset(0.001);
    
    cout << minFreq << "\t" << maxFreq << endl;
    
    histo1DFAllCell->DrawCopy(); 
        DrawGammaSetMarker(histoAllCells, 24, 1.5, kAzure+2, kAzure+2);
        histoAllCells->DrawCopy("same,e1"); 
//         DrawGammaSetMarker(histoAllCells_Erange, 21, 1.5, kGreen-2, kGreen-2);
//         histoAllCells_Erange->DrawCopy("same,e1"); 
//         DrawGammaSetMarker(histoAllCells_maxE, 25, 1.5, kGreen+2, kGreen+2);
//         histoAllCells_maxE->DrawCopy("same,e1"); 
//         DrawGammaSetMarker(histoAllCells_maxE_Erange, 24, 1.5, kBlue+1, kBlue+1);
//         histoAllCells_maxE_Erange->DrawCopy("same,e1"); 

//         TLegend* legendTMEffiChlead           = GetAndSetLegend2(0.42, 0.95, 0.7, 0.95-(4*0.04),32);
//         legendTMEffiChlead->AddEntry(histoEffi_ConvCl_leadReb,"conv #gamma cl","p");
//         legendTMEffiChlead->AddEntry(histoEffi_ElCl_leadReb,"e^{#pm} cl","p");
//         legendTMEffiChlead->AddEntry(histoEffi_ChPrimCl_leadReb,"other charged cl., prod. vtx R #leq 5 cm","p");
//         legendTMEffiChlead->AddEntry(histoEffi_ChCl_leadReb,"other charged cl., prod. vtx R > 5 cm","p");
//         legendTMEffiChlead->Draw();

        PutProcessLabelAndEnergyOnPlot(0.105, 0.19, 32, collisionSystem.Data(), detectionProcess.Data(),"",  63, 0.029);

    canvasClusterF->Update();
    canvasClusterF->SaveAs(Form("%s/CellFrequency_%s.png",outputDirectory.Data(),optionCutSelection.Data()));
//         
//         
//             TH2F * histo2DTMEffNe        = new TH2F("histo2DTMEffNe", "histo2DTMEffNe",1000, 0., maxPt , 1000, rangeTMEffi[0], rangeTMEffi[1] );
//             SetStyleHistoTH2ForGraphs(  histo2DTMEffNe, "#it{E}_{clus} (GeV)", "#varepsilon_{miss,TM}", 
//                                         0.85*0.04, 0.04, 0.85*0.04, 0.04, 0.95, 0.94);
//             histo2DTMEffNe->GetYaxis()->SetLabelOffset(0.001);
//             histo2DTMEffNe->DrawCopy(); 
//             DrawGammaSetMarker(histoEffi_NeSubChClReb, 20, 1.5, kRed-6, kRed-6);
//             histoEffi_NeSubChClReb->DrawCopy("same,e1"); 
//             DrawGammaSetMarker(histoEffi_NeClReb, 24, 1.5, kRed+2, kRed+2);
//             histoEffi_NeClReb->DrawCopy("same,e1"); 
//             DrawGammaSetMarker(histoEffi_GaClReb, 25, 1.5, kOrange, kOrange);
//             histoEffi_GaClReb->DrawCopy("same,e1"); 
//             DrawGammaSetMarker(histoEffi_GaSubChClReb, 21, 1.5, 807, 807);
//             histoEffi_GaSubChClReb->DrawCopy("same,e1"); 
// 
//             TLegend* legendTMEffiNe           = GetAndSetLegend2(0.42, 0.95, 0.7, 0.95-(4*0.04),32);
//             legendTMEffiNe->AddEntry(histoEffi_GaClReb,"#gamma cl.","p");
//             legendTMEffiNe->AddEntry(histoEffi_GaSubChClReb,"#gamma cl. w/ ch. cont","p");
//             legendTMEffiNe->AddEntry(histoEffi_NeClReb,"other neutral cl.","p");
//             legendTMEffiNe->AddEntry(histoEffi_NeSubChClReb,"other neutral cl. w/ ch. cont.","p");
//             legendTMEffiNe->Draw();
// 
//             PutProcessLabelAndEnergyOnPlot(0.105, 0.95, 32, collisionSystem.Data(), detectionProcess.Data(),"",  63, 0.029);
// 
//         canvasTMEffi->Update();
//         canvasTMEffi->SaveAs(Form("%s/TMEffiDifferentNeutralClusters_%s.%s",outputDirectory.Data(),optionCutSelection.Data(),suffix.Data()));
// 
//             TH2F * histFrac             = new TH2F("histFrac", "histFrac",1000, 0., maxPt , 1000, 0, rangeTMEffi[1] );
//             SetStyleHistoTH2ForGraphs(  histFrac, "#it{E}_{clus} (GeV)", "N_{id}/N_{all}", 
//                                         0.85*0.04, 0.04, 0.85*0.04, 0.04, 0.95, 0.94);
//             histFrac->GetYaxis()->SetLabelOffset(0.001);
//             histFrac->DrawCopy(); 
//             
//             DrawGammaSetMarker(histoRatioClE_ChtoAllReb, 24, 1.5, kAzure+2, kAzure+2);
//             histoRatioClE_ChtoAllReb->DrawCopy("same,e1"); 
//             DrawGammaSetMarker(histoRatioClE_EltoAllReb, 21, 1.5, kGreen-2, kGreen-2);
//             histoRatioClE_EltoAllReb->DrawCopy("same,e1"); 
//             DrawGammaSetMarker(histoRatioClE_ConvtoAllReb, 25, 1.5, kGreen+2, kGreen+2);
//             histoRatioClE_ConvtoAllReb->DrawCopy("same,e1"); 
//             DrawGammaSetMarker(histoRatioClE_ChPrimtoAllReb, 24, 1.5, kBlue+1, kBlue+1);
//             histoRatioClE_ChPrimtoAllReb->DrawCopy("same,e1"); 
// 
//             DrawGammaSetMarker(histoRatioClE_NeSubChtoAllReb, 20, 1.5, kRed-6, kRed-6);
//             histoRatioClE_NeSubChtoAllReb->DrawCopy("same,e1"); 
//             DrawGammaSetMarker(histoRatioClE_NetoAllReb, 24, 1.5, kRed+2, kRed+2);
//             histoRatioClE_NetoAllReb->DrawCopy("same,e1"); 
//             DrawGammaSetMarker(histoRatioClE_GatoAllReb, 25, 1.5, kOrange, kOrange);
//             histoRatioClE_GatoAllReb->DrawCopy("same,e1"); 
//             DrawGammaSetMarker(histoRatioClE_GaSubChtoAllReb, 21, 1.5, 807, 807);
//             histoRatioClE_GaSubChtoAllReb->DrawCopy("same,e1"); 
// 
//             TLegend* legendRatio           = GetAndSetLegend2(0.42, 0.95, 0.7, 0.95-(8*0.04),32);
//             legendRatio->AddEntry(histoEffi_ConvClReb,"conv #gamma cl","p");
//             legendRatio->AddEntry(histoEffi_ElClReb,"e^{#pm} cl","p");
//             legendRatio->AddEntry(histoEffi_ChPrimClReb,"other charged cl., prod. vtx R #leq 5 cm","p");
//             legendRatio->AddEntry(histoEffi_ChClReb,"other charged cl., prod. vtx R > 5 cm","p");
//             legendRatio->AddEntry(histoEffi_GaClReb,"#gamma cl.","p");
//             legendRatio->AddEntry(histoEffi_GaSubChClReb,"#gamma cl. w/ ch. cont","p");
//             legendRatio->AddEntry(histoEffi_NeClReb,"other neutral cl.","p");
//             legendRatio->AddEntry(histoEffi_NeSubChClReb,"other neutral cl. w/ ch. cont.","p");
//             legendRatio->Draw();
// 
//             PutProcessLabelAndEnergyOnPlot(0.105, 0.95, 32, collisionSystem.Data(), detectionProcess.Data(),"",  63, 0.029);
// 
//         canvasTMEffi->Update();
//         canvasTMEffi->SaveAs(Form("%s/TMFractioOfClusters_%s.%s",outputDirectory.Data(),optionCutSelection.Data(),suffix.Data()));
//         
    
    TFile* fileOutput = new TFile("TriggerCrossCheck.root","UPDATE");  
        fEventQuality->Write(Form("NEvents_%s",nameTrigger.Data()));
        histoCellIDvsClusterEnergy->Write(Form("CellIDvsClusterEnergy_%s", nameTrigger.Data()));
        histoCellIDvsClusterEnergyMax->Write(Form("CellIDvsClusterEnergyMax_%s", nameTrigger.Data())); 
        histoClE_AllCell->Write(Form("ClusterE_AllCells_%s", nameTrigger.Data()));
        histoClEMax_AllCell->Write(Form("ClusterEmaxClPerEvent_AllCells_%s", nameTrigger.Data()));
        histoAllCells->Write(Form("Frequency_CellID_AllCells_%s", nameTrigger.Data()));
        histoAllCells_Erange->Write(Form("Frequency_CellID_AllCells_aboveThreshold_%s", nameTrigger.Data()));
        histoBadChannel->Write(Form("BadChannel_CellID_AllCells_%s", nameTrigger.Data()));
        histoBadChannel_Erange->Write(Form("BadChannel_CellID_AllCells_aboveThreshold_%s", nameTrigger.Data()));
    
        histoAllCells_maxE->Write(Form("Frequency_CellIDEmaxCl_AllCells_%s", nameTrigger.Data()));
        histoAllCells_maxE_Erange->Write(Form("Frequency_CellIDEmaxCl_AllCells_aboveThreshold_%s", nameTrigger.Data()));
        histoBadChannelTrigg->Write(Form("BadChannel_CellIDEmaxCl_AllCells_%s", nameTrigger.Data()));
        histoBadChannelTrigg_Erange->Write(Form("BadChannel_CellIDEmaxCl_AllCells_aboveThreshold_%s", nameTrigger.Data()));
    fileOutput->Write();
    fileOutput->Close();
    
    return;

}
