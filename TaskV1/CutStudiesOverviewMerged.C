//***********************************************************************************************
//**************************** CutStudiesOverview ***********************************************
//***********************************************************************************************
/************************************************************************************************
 ******     provided by Gamma Conversion Group, PWG4,                                         *****
 ******        Ana Marin, marin@physi.uni-heidelberg.de                                        *****
 ******        Friederike Bock, friederike.bock@cern.ch                                        *****
 ***********************************************************************************************/

#include <Riostream.h>
#include <fstream>
#include "TMath.h"
#include <stdio.h>
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
#include "../CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "../CommonHeaders/ConversionFunctions.h"

struct SysErrorConversion {
   Double_t value;
   Double_t error;
   //    TString name;
};


void CutStudiesOverviewMerged(  TString CombineCutsName     = "CombineCuts.dat", 
                                TString suffix              = "gif", 
                                TString meson               = "", 
                                TString isMC                = "", 
                                TString optionEnergy        = "", 
                                TString cutVariationName    = "", 
                                Int_t numberOfCuts          = 1, 
                                TString optionPeriod        = "No", 
                                Int_t mode                  = 10,
                                Bool_t doBarlow             = kFALSE ){

    
    if (!(mode == 10 || mode == 11 )){
        cout << "incorrect mode: " << mode << endl;
        return ;
    }    
    
    // Define global arrays
    TString     cutNumber               [50];
    TString     cutNumberAdv            [50];
    Double_t    nColls                  [50];
    TString     prefix2                                         = "";
    Bool_t      correctionFilesAvail                            = kTRUE;
    
    // Set common default plot style
    StyleSettingsThesis();
    SetPlotStyle();

    // Set cutvariation-name to "" for no explicit name
    if (cutVariationName.CompareTo("None")==0){
        cutVariationName                                        = "";
    }
        
    // Define Output Directory    
    TString outputDir                                           = Form("CutStudies/%s",optionEnergy.Data());
    if (cutVariationName.CompareTo("None")!=0) outputDir        = Form("CutStudies/%s/%s",optionEnergy.Data(),cutVariationName.Data());
    TString outputDirRootFile                                   = Form("CutStudies/%s",optionEnergy.Data());
    gSystem->Exec("mkdir -p "+outputDir);
    
    // Define meson names for plots
    TString textMeson;
    Bool_t isEta                                                = ReturnMesonString ( meson );
    if ( meson.CompareTo("Eta")==0 ){
        isEta                                                   = kTRUE;
    } 

    // Define input and output MC/data
    if ( isMC.CompareTo("kTRUE") ==0){
        prefix2                                                 = "MC";
    } else {
        prefix2                                                 = "data";
    }
    
    // Set collisions system
    TString collisionSystem                                     = ReturnFullCollisionsSystem(optionEnergy);   
    if (collisionSystem.CompareTo("") == 0){
        cout << "No correct collision system specification, has been given" << endl;
        return;     
    }
    TString detectionProcess                                    = ReturnFullTextReconstructionProcess(mode, 0, textMeson);
    TString fNLMString                                          = "";
    Int_t fNLMmin                                               = 0;
    
    // Define colors for differnt cuts
    Color_t color[20]                                           = { kBlack, kAzure, kGreen+2,kOrange+2,kRed, kViolet,  kBlue-9, kMagenta+4,
                                                                    kCyan+3, kCyan-10, kCyan, kGreen+4, kGreen-9, 
                                                                    kGreen,  kYellow+4, kYellow+3, kSpring+10,
                                                                    kMagenta-8, kGray, kGray+3};

    // Read cuts from CutSelection file                    
    ifstream in(CombineCutsName.Data());
    cout<<"Available Cuts:"<<endl;
    string TempCutNumber;
    Int_t Number = 0;
    while(getline(in, TempCutNumber)){
        TString tempCutNumberAdv                                = TempCutNumber;
        cutNumberAdv[Number]                                    = tempCutNumberAdv;
        cutNumber[Number]                                       = tempCutNumberAdv;
        cout<< cutNumber[Number]<<endl;
        Number++;
    }
    cout<<"=========================="<<endl;

    cout << "analysing " << cutVariationName << " cut variations" << endl;

    cout << " " << endl;

    // Determine number of collisions for PbPb case
    if( optionEnergy.Contains("Pb") ) {
        for (Int_t i=0; i< numberOfCuts; i++){
            nColls[i]                                             = GetNCollFromCutNumber(cutNumber[i]);
        }
    } else {
        for (Int_t i=0; i< numberOfCuts; i++) nColls[i]         = 1.;
    }

    // Define necessary histogram/file/string arrays 
    const Int_t constNumberOfCuts                               = numberOfCuts;
    const Int_t maxNumberOfCuts                                 = 20;
    if(constNumberOfCuts > maxNumberOfCuts){
        cout << "Too many cuts, beware!" << endl;
        return;
    }
    TString nameCorrectedFile           [maxNumberOfCuts];
    TString nameUnCorrectedFile         [maxNumberOfCuts];
    TFile*  fileCutCorr                 [constNumberOfCuts];
    TFile*  fileCutUncorr               [constNumberOfCuts];
    TH1D*   histoCorrectedYieldCut      [constNumberOfCuts];
    TH1D*   histoRawClusterPtCut        [constNumberOfCuts];
    TH1D*   histoEffiCut                [constNumberOfCuts];
    TH1D*   histoPurityCut              [constNumberOfCuts];
    TH1D*   histoAcceptanceCut          [constNumberOfCuts];
    TH1D*   histoRawYieldCut            [constNumberOfCuts];
    TH1D*   histoRatioCorrectedYieldCut [constNumberOfCuts];
    TH1D*   histoRatioEffiCut           [constNumberOfCuts];
    TH1D*   histoRatioPurityCut         [constNumberOfCuts];
    TH1D*   histoRatioAcceptanceCut     [constNumberOfCuts];
    TH1D*   histoRatioRawYieldCut       [constNumberOfCuts];
    TH1D*   histoRatioRawClusterPtCut   [constNumberOfCuts];
    TString cutStringsName              [maxNumberOfCuts];

    Double_t           maxPt                                    = 0;
    Bool_t             kSpecialTrigger                          = kFALSE;
    
    if (cutVariationName.CompareTo("SpecialTrigg") == 0){
        kSpecialTrigger                                         = kTRUE;
    }        
    for (Int_t i=0; i< numberOfCuts; i++){
        // Define CutSelections
        TString fEventCutSelection                  = "";
        TString fClusterCutSelection                = "";
        TString fClusterMergedCutSelection          = "";
        TString dummyString                         = "";
        TString fMesonCutSelection                  = "";
        ReturnSeparatedCutNumberAdvanced( cutNumberAdv[i].Data(), fEventCutSelection, fClusterCutSelection, fClusterMergedCutSelection, dummyString, fMesonCutSelection, mode);
        if (i == 0){
            fNLMmin = ReturnClusterNLM(fClusterMergedCutSelection);
            if (fNLMmin)
                fNLMString                          = Form("%i local maximum", fNLMmin);        
            else 
                fNLMString                          = Form("%i local maxima", fNLMmin);
        } else {
            if ( fNLMmin != ReturnClusterNLM(fClusterMergedCutSelection))
                fNLMString                          = "";
        }
        // Check if there was a special trigger among the cuts
        TString fTrigger = fEventCutSelection(GetEventSelectSpecialTriggerCutPosition(),1);
//         if (fTrigger.Atoi()>3 && cutVariationName.CompareTo("") == 0){
//             cutVariationName                                  ="SpecialTrigg";
//             kSpecialTrigger                                   = kTRUE;
//         }    
        // only read corrected file if "Special trigger was used"

        nameCorrectedFile[i]                                = Form( "%s/%s/%s_%s_GammaMergedCorrection_%s.root", cutNumberAdv[i].Data(), optionEnergy.Data(), meson.Data(), prefix2.Data(), 
                                                                    cutNumber[i].Data());
        cout<< nameCorrectedFile[i] << endl;
        fileCutCorr[i]                                      = new TFile(nameCorrectedFile[i]);    
        if (fileCutCorr[i]->IsZombie()){
            correctionFilesAvail                            = kFALSE;
        }

        // read uncorrected file
        nameUnCorrectedFile[i]                                  = Form("%s/%s/%s_%s_GammaMergedWithoutCorrection_%s.root",cutNumberAdv[i].Data(), optionEnergy.Data(), meson.Data(), prefix2.Data(), 
                                                                        cutNumber[i].Data());
        cout<< nameUnCorrectedFile[i] << endl;
        fileCutUncorr[i]                                        = new TFile(nameUnCorrectedFile[i]);
        if (fileCutUncorr[i]->IsZombie()) return;
        
        // put proper cutvariation labeling for plots
        if (cutVariationName.Contains("SpecialTrigg")){
            fTrigger                                            = fEventCutSelection(GetEventSelectSpecialTriggerCutPosition(),2);
            cutStringsName[i]                                   = AnalyseSpecialTriggerCut(fTrigger.Atoi(), optionPeriod);      
        } else if (cutVariationName.Contains("Rapidity")){
            TString fRapidityCut                                = fMesonCutSelection(GetMesonRapidityCutPosition(),1);
            cutStringsName[i]                                   = AnalyseRapidityMesonCut(fRapidityCut.Atoi());      
        } else if (cutVariationName.Contains("MesonAlpha")){
            TString fAlphaCut                                   = fMesonCutSelection(GetMesonAlphaCutPosition(),1);
            cutStringsName[i]                                   = fAlphaCut;
        } else if (cutVariationName.Contains("MesonOpeningAngle")){
          TString fMesonOpeningAngleCut                         = fMesonCutSelection(GetMesonOpeningAngleCutPosition(),1);
          cutStringsName[i]                                     = fMesonOpeningAngleCut;
        } else if (cutVariationName.Contains("MesonMass")){
          TString fMassMassCut                                  = fMesonCutSelection(GetMesonSelectionWindowCutPosition(),1);
          cutStringsName[i]                                     = fMassMassCut;
        } else if (cutVariationName.Contains("DiffRapWindow")){
            TString fRapidityCut                                = fMesonCutSelection(GetMesonRapidityCutPosition(),1);
            cutStringsName[i]                                   = AnalyseRapidityMesonCutpPb(fRapidityCut.Atoi());      
        } else if (cutVariationName.Contains("ClusterTrackMatchingCalo")){    
            TString fTrackMatching                              = fClusterCutSelection(GetClusterTrackMatchingCutPosition(fClusterCutSelection),1);
            cutStringsName[i]                                   = AnalyseTrackMatchingCaloCut(fTrackMatching.Atoi());
        } else if (cutVariationName.Contains("ClusterMaterialTRD")){    
            TString fMinPhi                                     = fClusterCutSelection(GetClusterPhiMinCutPosition(fClusterCutSelection),1);
            TString fMaxPhi                                     = fClusterCutSelection(GetClusterPhiMaxCutPosition(fClusterCutSelection),1);
            cutStringsName[i]                                   = AnalyseAcceptanceCutPhiCluster(fMinPhi.Atoi(), fMaxPhi.Atoi());
        } else if (cutVariationName.Contains("ClusterM02")){    
            TString fMinM02Cut                                  = fClusterMergedCutSelection(GetClusterMinM02CutPosition(fClusterMergedCutSelection),1);
            TString fMaxM02Cut                                  = fClusterMergedCutSelection(GetClusterMaxM02CutPosition(fClusterMergedCutSelection),1);
            cutStringsName[i]                                   = Form("%s%s",fMinM02Cut.Data(),fMaxM02Cut.Data());
//             AnalyseM02Cut(fMinM02Cut.Atoi(), fMaxM02Cut.Atoi());
        } else if (cutVariationName.Contains("ClusterNCells")){    
            TString fNCellsCut                                  = fClusterCutSelection(GetClusterMinNCellsCutPosition(fClusterCutSelection),1);
            cutStringsName[i]                                   = AnalyseNCellsCut(fNCellsCut.Atoi());
        } else if (cutVariationName.Contains("ClusterMinEnergy")){    
            TString fMinEnergyCut                               = fClusterMergedCutSelection(GetClusterMinEnergyCutPosition(fClusterMergedCutSelection),1);
            cout << fMinEnergyCut << "\t" << GetClusterMinEnergyCutPosition(fClusterMergedCutSelection) << "\t"<< fClusterCutSelection.Length()<<endl;
            cutStringsName[i]                                   = AnalyseMinEnergyMergedCut(fMinEnergyCut.Atoi());
        } else if (cutVariationName.Contains("ClusterTiming")){    
            TString fTimingCut                                  = fClusterCutSelection(GetClusterTimingCutPosition(fClusterCutSelection),1);
            cutStringsName[i]                                   = AnalyseClusterTimingCut(fTimingCut.Atoi());
        } else if (cutVariationName.Contains("ClusterNonLinearity")){
          TString fClusterNonLinearity                          = fClusterCutSelection(GetClusterNonLinearityCutPosition(fClusterCutSelection),2);
          cutStringsName[i]                                     = AnalyseClusterNonLinearityCut(fClusterNonLinearity.Atoi());
        } else {
            cutStringsName[i]                                   = cutNumberAdv[i].Data();
        }
        
        // read histograms from corrected file
        if (correctionFilesAvail){
            // for first cut read yield extraction errors as well            
            TString nameCorrectedYield                          = "CorrectedYieldTrueEff";
            TString nameEfficiency                              = "TrueMesonEffiPrimPt";
            TString namePurity                                  = "TruePi0PurityMergedPt";
            TString nameAcceptance                              = "fHistoMCAcceptancePt";
            
            cout << "line " << __LINE__ << endl;
            histoCorrectedYieldCut[i]                           = (TH1D*)fileCutCorr[i]->Get(nameCorrectedYield.Data());
            histoCorrectedYieldCut[i]->SetName(Form("%s_%s", nameCorrectedYield.Data(),cutNumber[i].Data()));
            cout << "line " << __LINE__ << endl;
            histoEffiCut[i]                                     = (TH1D*)fileCutCorr[i]->Get(nameEfficiency.Data());
            histoEffiCut[i]->SetName(Form("%s_%s", nameEfficiency.Data(), cutNumber[i].Data()));
            cout << "line " << __LINE__ << endl;
            histoPurityCut[i]                                   = (TH1D*)fileCutCorr[i]->Get(namePurity.Data());
            histoPurityCut[i]->SetName(Form("%s_%s", namePurity.Data(), cutNumber[i].Data()));
            cout << "line " << __LINE__ << endl;
            histoAcceptanceCut[i]                               = (TH1D*)fileCutCorr[i]->Get(nameAcceptance.Data());
            histoAcceptanceCut[i]->SetName(Form("%s_%s", nameAcceptance.Data(), cutNumber[i].Data()));
            cout << "line " << __LINE__ << endl;
            histoRawYieldCut[i]                                 = (TH1D*)fileCutCorr[i]->Get("histoYieldMesonPerEvent");
            histoRawYieldCut[i]->SetName(Form("histoYieldMesonPerEvent_%s",cutNumber[i].Data()));
            histoRawYieldCut[i]->Scale(1./nColls[i]);

            
        }
        // read histograms from uncorrected file
        cout << "line " << __LINE__ << endl;
        histoRawClusterPtCut[i]                                 = (TH1D*)fileCutUncorr[i]->Get("ClusterPtPerEvent");
        
        // calculate ratios for meson measurements
        if (correctionFilesAvail){
            histoRatioCorrectedYieldCut[i]                      = (TH1D*) histoCorrectedYieldCut[i]->Clone(Form("histoRatioCorrectedYieldCut_%s", cutNumber[i].Data()));
            histoRatioCorrectedYieldCut[i]->Divide(histoRatioCorrectedYieldCut[i],histoCorrectedYieldCut[0],1.,1.,"B");
            if (i > 0){
                maxPt= histoCorrectedYieldCut[i]->GetBinCenter(histoCorrectedYieldCut[i]->GetNbinsX()) + 0.5* histoCorrectedYieldCut[i]->GetBinWidth(histoCorrectedYieldCut[i]->GetNbinsX());
            }
            cout << "line " << __LINE__ << endl;
            histoRatioEffiCut[i]                            = (TH1D*) histoEffiCut[i]->Clone(Form("histoRatioEffiCut_%s", cutNumber[i].Data()));
            histoRatioEffiCut[i]->Divide(histoRatioEffiCut[i],histoEffiCut[0],1.,1.,"B");
            cout << "line " << __LINE__ << endl;
            histoRatioPurityCut[i]                            = (TH1D*) histoPurityCut[i]->Clone(Form("histoRatioPurityCut_%s", cutNumber[i].Data()));
            histoRatioPurityCut[i]->Divide(histoRatioPurityCut[i],histoPurityCut[0],1.,1.,"B");
            cout << "line " << __LINE__ << endl;
            histoRatioAcceptanceCut[i]                          = (TH1D*) histoAcceptanceCut[i]->Clone(Form("histoRatioAcceptanceCut_%s", cutNumber[i].Data()));
            histoRatioAcceptanceCut[i]->Divide(histoRatioAcceptanceCut[i],histoAcceptanceCut[0],1.,1.,"B");
            cout << "line " << __LINE__ << endl;
        }    
        histoRatioRawYieldCut[i]                                = (TH1D*) histoRawYieldCut[i]->Clone(Form("histoRatioRawYieldCut_%s", cutNumber[i].Data()));
        if (!kSpecialTrigger){
            histoRatioRawYieldCut[i]->Divide(histoRatioRawYieldCut[i],histoRawYieldCut[0],1.,1.,"B");
        } else {
            histoRatioRawYieldCut[i]->Divide(histoRatioRawYieldCut[i],histoRawYieldCut[0],1.,1.,"");
        }    
        cout << "line " << __LINE__ << endl;
        histoRatioRawClusterPtCut[i]                        = (TH1D*) histoRawClusterPtCut[i]->Clone(Form("histoRatioRawClusterPtCut%s", cutNumber[i].Data()));
        if (!kSpecialTrigger){
            histoRatioRawClusterPtCut[i]->Divide(histoRatioRawClusterPtCut[i],histoRawClusterPtCut[0],1.,1.,"B");
        } else {
            histoRatioRawClusterPtCut[i]->Divide(histoRatioRawClusterPtCut[i],histoRawClusterPtCut[0],1.,1.,"");
        }    
        cout << "line " << __LINE__ << endl;
        
    }

    cout<<"=========================="<<endl;



    if (cutVariationName.Contains("SpecialTrigg")){
        //**************************************************************************************
        //********************* Plotting RAW-Yield for special triggers  ***********************
        //**************************************************************************************        

        TCanvas* canvasRawYieldsTrigg = new TCanvas("canvasRawYieldsTrigg","",1000,1000);  // gives the page size
        DrawGammaCanvasSettings( canvasRawYieldsTrigg,  0.1, 0.02, 0.04, 0.08);
        canvasRawYieldsTrigg->SetLogy(1);

        // Set legend
        
        TLegend* legendSpecRawTrigger = GetAndSetLegend2(0.15,0.15,0.3,0.15+1.15*0.032*numberOfCuts, 1000*0.032); 
        // Draw Raw yield for different triggers
        for(Int_t i = 0; i< numberOfCuts; i++){
            if(i == 0){
                Double_t scaleFactorRaw = 10.;
                DrawGammaSetMarker(histoRawYieldCut[i], 20, 1., color[0], color[0]);
                DrawAutoGammaMesonHistos( histoRawYieldCut[i],
                                        "", "#it{p}_{T} (GeV/#it{c})", Form("%s RAW Yield/(#it{N}_{ev})",textMeson.Data()),
                                        kTRUE, scaleFactorRaw, 10e-10, kTRUE,
                                        kFALSE, 0.0, 0.030,
                                        kFALSE, 0., 10.);
                legendSpecRawTrigger->AddEntry(histoRawYieldCut[i],Form("standard: %s",cutStringsName[i].Data()));
            } else {
                if(i<20){
                    DrawGammaSetMarker(histoRawYieldCut[i], 20+i, 1.,color[i],color[i]);
                } else {
                    DrawGammaSetMarker(histoRawYieldCut[i], 20+i, 1.,color[i-20],color[i-20]);
                }
                histoRawYieldCut[i]->DrawCopy("same,e1,p");
                legendSpecRawTrigger->AddEntry(histoRawYieldCut[i],cutStringsName[i].Data());
            }
        }
        legendSpecRawTrigger->Draw();
        // labeling of the plot
        PutProcessLabelAndEnergyOnPlot( 0.65, 0.95, 0.032, collisionSystem, fNLMString, detectionProcess, 42, 0.03, optionPeriod);
        
        canvasRawYieldsTrigg->Update();
        canvasRawYieldsTrigg->SaveAs(Form("%s/%s_%s_TriggerYieldSpectra.%s",outputDir.Data(),meson.Data(),prefix2.Data(),suffix.Data()));
        delete canvasRawYieldsTrigg;

        //**************************************************************************************
        //***************** Plotting RAW-Yield ratios for special triggers  ********************
        //**************************************************************************************                
        TCanvas* canvasRatioRawYields = new TCanvas("canvasRatioRawYields","",1000,1000);  // gives the page size
        DrawGammaCanvasSettings( canvasRatioRawYields,  0.1, 0.02, 0.04, 0.08);
        canvasRatioRawYields->SetLogy(1); 
        // create legend
        TLegend* legendRatioRaw = GetAndSetLegend2(0.55,0.15,0.7,0.15+1.15*0.032*numberOfCuts, 1000*0.032);
        // find min und max
        Double_t maxRatio = 2e3;
        Double_t minRatio = 0;
        
        // plot ratios in canvas
        for(Int_t i = 1; i< numberOfCuts; i++){
            if(i == 1){
                DrawGammaSetMarker(histoRatioRawYieldCut[i], 20, 1., color[1], color[1]);
                DrawAutoGammaMesonHistos( histoRatioRawYieldCut[i],
                                        "", "#it{p}_{T} (GeV/#it{c})", Form("triggered Yield/standard"),
                                        kFALSE, 5., 10e-10, kTRUE,
                                        kTRUE, minRatio, maxRatio,
                                        kFALSE, 0., 10.);
                legendRatioRaw->AddEntry(histoRatioRawYieldCut[i],Form("%s",cutStringsName[i].Data()));
            } else {
                if(i<20){
                    DrawGammaSetMarker(histoRatioRawYieldCut[i], 20+i, 1.,color[i],color[i]);
                } else {
                    DrawGammaSetMarker(histoRatioRawYieldCut[i], 20+i, 1.,color[i-20],color[i-20]);
                }
                histoRatioRawYieldCut[i]->DrawCopy("same,e1,p");
                legendRatioRaw->AddEntry(histoRatioRawYieldCut[i],cutStringsName[i].Data());
            }
        }
        // labeling
        legendRatioRaw->Draw();
        PutProcessLabelAndEnergyOnPlot( 0.65, 0.95, 0.032, collisionSystem, fNLMString, detectionProcess, 42, 0.03, optionPeriod);
        
        canvasRatioRawYields->Update();
        canvasRatioRawYields->SaveAs(Form("%s/%s_%s_TriggerYieldRatio.%s",outputDir.Data(),meson.Data(),prefix2.Data(),suffix.Data()));
        delete canvasRatioRawYields;
    }
    
        
    //**************************************************************************************
    //********************* Plotting RAW-Cluster Yield *********************************************
    //**************************************************************************************

    TCanvas* canvasRawClusterPt = new TCanvas("canvasRawClusterPt","",1350,1500);  
    DrawGammaCanvasSettings( canvasRawClusterPt,  0.13, 0.02, 0.02, 0.09);
    // Upper pad definition
    TPad* padRawClusterPt = new TPad("padRawClusterPt", "", 0., 0.33, 1., 1.,-1, -1, -2);
    DrawGammaPadSettings( padRawClusterPt, 0.12, 0.02, 0.02, 0.);
    padRawClusterPt->SetLogy();
    padRawClusterPt->Draw();
    // lower pad definition
    TPad* padRawClusterPtRatios = new TPad("padRawClusterPtRatios", "", 0., 0., 1., 0.33,-1, -1, -2);
    DrawGammaPadSettings( padRawClusterPtRatios, 0.12, 0.02, 0.0, 0.2);
    padRawClusterPtRatios->Draw();

    padRawClusterPt->cd();
    padRawClusterPt->SetTickx();
    padRawClusterPt->SetTicky();

    // Plot raw yield in uppper panel
    padRawClusterPt->cd();
    TLegend* legendRawCluster = GetAndSetLegend2(0.15,0.02,0.3,0.02+1.15*0.032*numberOfCuts, 1500*0.75*0.032);
    for(Int_t i = 0; i< numberOfCuts; i++){
        if(i == 0){
            Double_t scaleFactorRaw = 5.;
            if (kSpecialTrigger) scaleFactorRaw = 5.;
            
            DrawGammaSetMarker(histoRawClusterPtCut[i], 20, 1., color[0], color[0]);
            DrawAutoGammaMesonHistos( histoRawClusterPtCut[i],
                                    "", "#it{p}_{T} (GeV/#it{c})", "cluster RAW Yield/(#it{N}_{ev} #it{N}_{coll})",
                                    kTRUE, scaleFactorRaw, 10e-10, kTRUE,
                                    kFALSE, 0.0, 0.030,
                                    kFALSE, 0., 10.);
            legendRawCluster->AddEntry(histoRawClusterPtCut[i],Form("standard: %s",cutStringsName[i].Data()));
        }
        else {
            if(i<20){
                DrawGammaSetMarker(histoRawClusterPtCut[i], 20+i, 1.,color[i],color[i]);
            } else {
                DrawGammaSetMarker(histoRawClusterPtCut[i], 20+i, 1.,color[i-20],color[i-20]);
            }
            histoRawClusterPtCut[i]->DrawCopy("same,e1,p");
            legendRawCluster->AddEntry(histoRawClusterPtCut[i],cutStringsName[i].Data());
        }
        
    }
    legendRawCluster->Draw();
    // Labeling of plot
    PutProcessLabelAndEnergyOnPlot( 0.65, 0.95, 0.032, collisionSystem, "#gamma candidates", detectionProcess, 42, 0.03, optionPeriod); 
        
    padRawClusterPtRatios->cd();
    for(Int_t i = 0; i< numberOfCuts; i++){
        if(i==0){
            // Set ratio min and max
            Double_t minYRatio = 0.45;
            Double_t maxYRatio = 2.05; //qui
            SetStyleHistoTH1ForGraphs(histoRatioRawClusterPtCut[i], "#it{p}_{T} (GeV/#it{c})", "#frac{modified}{standard}", 0.08, 0.11, 0.07, 0.1, 0.75, 0.5, 510,505);
            DrawGammaSetMarker(histoRatioRawClusterPtCut[i], 20, 1.,color[0],color[0]);
            histoRatioRawClusterPtCut[i]->GetYaxis()->SetRangeUser(minYRatio,maxYRatio);
            histoRatioRawClusterPtCut[i]->DrawCopy("p,e1");
        } else{
            if(i<20){
                DrawGammaSetMarker(histoRatioRawClusterPtCut[i], 20+i, 1.,color[i],color[i]);
            } else {
                DrawGammaSetMarker(histoRatioRawClusterPtCut[i], 20+i, 1.,color[i-20],color[i-20]);
            }
            histoRatioRawClusterPtCut[i]->DrawCopy("same,e1,p");
        }
        DrawGammaLines(0., 50,1., 1.,0.1);
    }

    canvasRawClusterPt->Update();
    canvasRawClusterPt->SaveAs(Form("%s/%s_%s_RAWYieldCluster.%s",outputDir.Data(),meson.Data(),prefix2.Data(),suffix.Data()));
    delete canvasRawClusterPt;

    if (cutVariationName.Contains("SpecialTrigg") && meson.CompareTo("Pi0") == 0 ){
 
        //**************************************************************************************
        //********************* Plotting RAW-Yield for special triggers  ***********************
        //**************************************************************************************        

        TCanvas* canvasRawClusterPtTrigg = new TCanvas("canvasRawClusterPtTrigg","",1000,1000);  // gives the page size
        DrawGammaCanvasSettings( canvasRawClusterPtTrigg,  0.1, 0.02, 0.04, 0.08);
        canvasRawClusterPtTrigg->SetLogy(1);

        // Set legend
        TLegend* legendSpecRawClusterTrigger =GetAndSetLegend2(0.15,0.11,0.3,0.11+1.15*0.032*numberOfCuts, 1000*0.032);
        // Draw Raw yield for different triggers
        for(Int_t i = 0; i< numberOfCuts; i++){
            if(i == 0){
                Double_t scaleFactorRaw = 10.;
                DrawGammaSetMarker(histoRawClusterPtCut[i], 20, 1., color[0], color[0]);
                DrawAutoGammaMesonHistos( histoRawClusterPtCut[i],
                                        "", "#it{p}_{T,cluster} (GeV/#it{c})", "cluster RAW Yield/(#it{N}_{ev})",
                                        kTRUE, scaleFactorRaw, 10e-10, kTRUE,
                                        kFALSE, 0.0, 0.030,
                                        kFALSE, 0., 10.);
                legendSpecRawClusterTrigger->AddEntry(histoRawClusterPtCut[i],Form("standard: %s",cutStringsName[i].Data()));
            } else {
                if(i<20){
                    DrawGammaSetMarker(histoRawClusterPtCut[i], 20+i, 1.,color[i],color[i]);
                } else {
                    DrawGammaSetMarker(histoRawClusterPtCut[i], 20+i, 1.,color[i-20],color[i-20]);
                }
                histoRawClusterPtCut[i]->DrawCopy("same,e1,p");
                legendSpecRawClusterTrigger->AddEntry(histoRawClusterPtCut[i],cutStringsName[i].Data());
            }
        }
        legendSpecRawClusterTrigger->Draw();
        // labeling of the plot
        PutProcessLabelAndEnergyOnPlot( 0.65, 0.95, 0.032, collisionSystem, "#gamma candidates", detectionProcess, 42, 0.03, optionPeriod);
        
        canvasRawClusterPtTrigg->Update();
        canvasRawClusterPtTrigg->SaveAs(Form("%s/%s_TriggerYieldCluster.%s",outputDir.Data(),prefix2.Data(),suffix.Data()));
        delete canvasRawClusterPtTrigg;

        //**************************************************************************************
        //***************** Plotting RAW-Yield ratios for special triggers  ********************
        //**************************************************************************************                
        TCanvas* canvasRatioRawClusterPt = new TCanvas("canvasRatioRawClusterPt","",1000,1000);  // gives the page size
        DrawGammaCanvasSettings( canvasRatioRawClusterPt,  0.1, 0.02, 0.04, 0.08);
        canvasRatioRawClusterPt->SetLogy(1); 
        // create legend
        TLegend* legendRatioRawClusterTrigger = GetAndSetLegend2(0.55,0.15,0.7,0.15+1.15*0.032*numberOfCuts, 1000*0.032);
        // find min und max
        Double_t maxRatio = 2e3;
        Double_t minRatio = 0;
        // plot ratios in canvas
        for(Int_t i = 1; i< numberOfCuts; i++){
            if(i == 1){
                DrawGammaSetMarker(histoRatioRawClusterPtCut[i], 20, 1., color[1], color[1]);
                DrawAutoGammaMesonHistos( histoRatioRawClusterPtCut[i],
                                        "", "#it{p}_{T,cluster} (GeV/#it{c})", Form("triggered Yield/standard"),
                                        kFALSE, 5., 10e-10, kTRUE,
                                        kTRUE, minRatio, maxRatio,
                                        kFALSE, 0., 10.);
                legendRatioRawClusterTrigger->AddEntry(histoRatioRawClusterPtCut[i],Form("%s",cutStringsName[i].Data()));
            } else {
                if(i<20){
                    DrawGammaSetMarker(histoRatioRawClusterPtCut[i], 20+i, 1.,color[i],color[i]);
                } else {
                    DrawGammaSetMarker(histoRatioRawClusterPtCut[i], 20+i, 1.,color[i-20],color[i-20]);
                }
                histoRatioRawClusterPtCut[i]->DrawCopy("same,e1,p");
                legendRatioRawClusterTrigger->AddEntry(histoRatioRawClusterPtCut[i],cutStringsName[i].Data());
            }
        }
        // labeling
        legendRatioRawClusterTrigger->Draw();
        PutProcessLabelAndEnergyOnPlot( 0.65, 0.95, 0.032, collisionSystem, "#gamma candidates", detectionProcess, 42, 0.03, optionPeriod);

        canvasRatioRawClusterPt->Update();
        canvasRatioRawClusterPt->SaveAs(Form("%s/%s_TriggerYieldClusterRatio.%s",outputDir.Data(),prefix2.Data(),suffix.Data()));
        delete canvasRatioRawClusterPt;
    }
        
    if (correctionFilesAvail){
        //**************************************************************************************
        //********************* Plotting RAW-Yield *********************************************
        //**************************************************************************************

        TCanvas* canvasRawYieldMeson = new TCanvas("canvasRawYieldMeson","",1350,1500);  
        DrawGammaCanvasSettings( canvasRawYieldMeson,  0.13, 0.02, 0.02, 0.09);
        // Upper pad definition
        TPad* padRawYield = new TPad("padRawYield", "", 0., 0.33, 1., 1.,-1, -1, -2);
        DrawGammaPadSettings( padRawYield, 0.12, 0.02, 0.02, 0.);
        padRawYield->SetLogy();
        padRawYield->Draw();
        // lower pad definition
        TPad* padRawYieldRatios = new TPad("padRawYieldRatios", "", 0., 0., 1., 0.33,-1, -1, -2);
        DrawGammaPadSettings( padRawYieldRatios, 0.12, 0.02, 0.0, 0.2);
        padRawYieldRatios->Draw();

        padRawYield->cd();
        padRawYield->SetTickx();
        padRawYield->SetTicky();

        // Plot raw yield in uppper panel
        padRawYield->cd();
        TLegend* legendRawMeson = GetAndSetLegend2(0.15,0.02,0.3,0.02+1.15*0.032*numberOfCuts, 1500*0.75*0.032);
        for(Int_t i = 0; i< numberOfCuts; i++){
            if(i == 0){
                Double_t scaleFactorRaw = 5.;
                if (kSpecialTrigger) scaleFactorRaw = 100.;
                
                DrawGammaSetMarker(histoRawYieldCut[i], 20, 1., color[0], color[0]);
                DrawAutoGammaMesonHistos( histoRawYieldCut[i],
                                        "", "#it{p}_{T} (GeV/#it{c})", Form("%s RAW Yield/(#it{N}_{ev} #it{N}_{coll})",textMeson.Data()),
                                        kTRUE, scaleFactorRaw, 10e-10, kTRUE,
                                        kFALSE, 0.0, 0.030,
                                        kFALSE, 0., 10.);
                legendRawMeson->AddEntry(histoRawYieldCut[i],Form("standard: %s",cutStringsName[i].Data()));
            }
            else {
                if(i<20){
                    DrawGammaSetMarker(histoRawYieldCut[i], 20+i, 1.,color[i],color[i]);
                } else {
                    DrawGammaSetMarker(histoRawYieldCut[i], 20+i, 1.,color[i-20],color[i-20]);
                }
                histoRawYieldCut[i]->DrawCopy("same,e1,p");
                legendRawMeson->AddEntry(histoRawYieldCut[i],cutStringsName[i].Data());
            }
            
        }
        legendRawMeson->Draw();
        // Labeling of plot
        PutProcessLabelAndEnergyOnPlot( 0.65, 0.95, 0.032, collisionSystem, fNLMString, detectionProcess, 42, 0.03, optionPeriod);
                                                
        padRawYieldRatios->cd();
        for(Int_t i = 0; i< numberOfCuts; i++){
            if(i==0){
                // Set ratio min and max
                Double_t minYRatio = 0.01;
                Double_t maxYRatio = 2.1; 
                SetStyleHistoTH1ForGraphs(histoRatioRawYieldCut[i], "#it{p}_{T} (GeV/#it{c})", "#frac{modified}{standard}", 0.08, 0.11, 0.07, 0.1, 0.75, 0.5, 510,505);
                DrawGammaSetMarker(histoRatioRawYieldCut[i], 20, 1.,color[0],color[0]);
                histoRatioRawYieldCut[i]->GetYaxis()->SetRangeUser(minYRatio,maxYRatio);
                
                for(Int_t b = 0; b< histoRatioRawYieldCut[i]->GetNbinsX(); b++){
                    histoRatioRawYieldCut[i]->SetBinError(b+1,histoRawYieldCut[i]->GetBinError(b+1)/histoRawYieldCut[i]->GetBinContent(b+1));
                }
                histoRatioRawYieldCut[i]->SetFillColor(kGray+2);
                histoRatioRawYieldCut[i]->SetFillStyle(0);
                histoRatioRawYieldCut[i]->DrawCopy("p,e2");  

            } else{
                if(i<20){
                    DrawGammaSetMarker(histoRatioRawYieldCut[i], 20+i, 1.,color[i],color[i]);
                } else {
                    DrawGammaSetMarker(histoRatioRawYieldCut[i], 20+i, 1.,color[i-20],color[i-20]);
                }
                histoRatioRawYieldCut[i]->DrawCopy("same,e1,p");
            }
            DrawGammaLines(0., maxPt,1., 1.,0.1);
        }

        canvasRawYieldMeson->Update();
        canvasRawYieldMeson->SaveAs(Form("%s/%s_%s_RAWYield.%s",outputDir.Data(),meson.Data(),prefix2.Data(),suffix.Data()));
        delete canvasRawYieldMeson;

        //*****************************************************************************************
        //******************* Compare Corrected Yields ********************************************
        //*****************************************************************************************
        // Define canvas
        TCanvas* canvasCorrectedYieldMeson = new TCanvas("canvasCorrectedYieldMeson","",1350,1500);  
        DrawGammaCanvasSettings( canvasCorrectedYieldMeson,  0.13, 0.02, 0.02, 0.09);
        // Define upper panel
        TPad* padCorrectedYield = new TPad("padCorrectedYield", "", 0., 0.33, 1., 1.,-1, -1, -2);
        DrawGammaPadSettings( padCorrectedYield, 0.12, 0.02, 0.02, 0.);
        padCorrectedYield->SetLogy(1);
        padCorrectedYield->Draw();
        // Define lower panel
        TPad* padCorrectedYieldRatios = new TPad("padCorrectedYieldRatios", "", 0., 0., 1., 0.33,-1, -1, -2);
        DrawGammaPadSettings( padCorrectedYieldRatios, 0.12, 0.02, 0.0, 0.2);
        padCorrectedYieldRatios->SetLogy(0);
        padCorrectedYieldRatios->Draw();

        // Plot corrected yield in upper panel
        padCorrectedYield->cd();            
        TLegend* legendCorrectedYieldMeson = GetAndSetLegend2(0.15,0.02,0.3,0.02+1.15*0.032*numberOfCuts, 1500*0.75*0.032); 
        for(Int_t i = 0; i< numberOfCuts; i++){
            if(i == 0){
                DrawAutoGammaMesonHistos( histoCorrectedYieldCut[i],
                                        "", "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)",
                                        kTRUE, 5., 5e-10,kTRUE,
                                        kFALSE, 0.0, 0.030,
                                        kFALSE, 0., 10.);
                DrawGammaSetMarker(histoCorrectedYieldCut[i], 20, 1., color[0], color[0]);
                histoCorrectedYieldCut[i]->DrawCopy("e1,p");
                legendCorrectedYieldMeson->AddEntry(histoCorrectedYieldCut[i], Form("standard: %s",cutStringsName[i].Data()));
            }
            else{
                if(i<20){
                    DrawGammaSetMarker(histoCorrectedYieldCut[i], 20+i, 1.,color[i],color[i]);
                } else {
                    DrawGammaSetMarker(histoCorrectedYieldCut[i], 20+i, 1.,color[i-20],color[i-20]);
                }
                histoCorrectedYieldCut[i]->DrawCopy("same,e1,p");
                legendCorrectedYieldMeson->AddEntry(histoCorrectedYieldCut[i], cutStringsName[i].Data());
            }
        }
        legendCorrectedYieldMeson->Draw();

        // labeling the plot
        PutProcessLabelAndEnergyOnPlot( 0.65, 0.95, 0.032, collisionSystem, fNLMString, detectionProcess, 42, 0.03, optionPeriod);
        
        // plot ratio of corrected yields in lower panel    
        padCorrectedYieldRatios->cd();
        for(Int_t i = 0; i< numberOfCuts; i++){
            if(i==0){
                // Set ratio min and max
                Double_t minYRatio = 0.45;
                Double_t maxYRatio = 1.55;
                SetStyleHistoTH1ForGraphs(histoRatioCorrectedYieldCut[i], "#it{p}_{T} (GeV/#it{c})", "#frac{modified}{standard}", 0.08, 0.11, 0.07, 0.1, 0.75, 0.5, 510,505);
                DrawGammaSetMarker(histoRatioCorrectedYieldCut[i], 20, 1.,color[0],color[0]);
                histoRatioCorrectedYieldCut[i]->GetYaxis()->SetRangeUser(minYRatio,maxYRatio);
                
                for(Int_t b = 0; b< histoRatioCorrectedYieldCut[i]->GetNbinsX(); b++){
                    histoRatioCorrectedYieldCut[i]->SetBinError(b+1,histoCorrectedYieldCut[i]->GetBinError(b+1)/histoCorrectedYieldCut[i]->GetBinContent(b+1));
                }
                histoRatioCorrectedYieldCut[i]->SetFillColor(kGray+2);
                histoRatioCorrectedYieldCut[i]->SetFillStyle(0);
                histoRatioCorrectedYieldCut[i]->DrawCopy("p,e2");  
    
//                 histoRatioCorrectedYieldCut[i]->DrawCopy("p,e1");                
            } else{
                if(i<20){
                    DrawGammaSetMarker(histoRatioCorrectedYieldCut[i], 20+i, 1.,color[i],color[i]);
                } else {
                    DrawGammaSetMarker(histoRatioCorrectedYieldCut[i], 20+i, 1.,color[i-20],color[i-20]);
                }
                histoRatioCorrectedYieldCut[i]->DrawCopy("same,e1,p");
            }
        }
        


        DrawGammaLines(0., maxPt,1., 1.,0.1);

        canvasCorrectedYieldMeson->Update();
        canvasCorrectedYieldMeson->SaveAs(Form("%s/%s_%s_CorrectedYield.%s",outputDir.Data(), meson.Data(),prefix2.Data(),suffix.Data()));
        delete canvasCorrectedYieldMeson;


        //**************************************************************************************
        //********************* Plotting efficiencies *********************************************
        //**************************************************************************************
        // Define canvas    
        TCanvas* canvasEffiMeson = new TCanvas("canvasEffiMeson","",1350,1500);  // gives the page size
        DrawGammaCanvasSettings( canvasEffiMeson,  0.13, 0.02, 0.02, 0.09);
        // Define upper panel
        TPad* padEffi = new TPad("padEffi", "", 0., 0.25, 1., 1.,-1, -1, -2);
        DrawGammaPadSettings( padEffi, 0.12, 0.02, 0.04, 0.);
        padEffi->Draw();
        // Define lower panel
        TPad* padEffiRatios = new TPad("padEffiRatios", "", 0., 0., 1., 0.25,-1, -1, -2);
        DrawGammaPadSettings( padEffiRatios, 0.12, 0.02, 0.0, 0.2);
        padEffiRatios->Draw();

        // draw efficiency in upper panel
        padEffi->cd();
        padEffi->SetLogy(1);

        TLegend* legendEffiMeson = GetAndSetLegend2(0.15,0.92-1.15*0.032*numberOfCuts,0.3,0.92, 1500*0.75*0.032); 
        for(Int_t i = 0; i< numberOfCuts; i++){
            if(i == 0){
                DrawGammaSetMarker(histoEffiCut[i], 20, 1., color[0], color[0]);
                DrawAutoGammaMesonHistos( histoEffiCut[i],
                                        "", "#it{p}_{T} (GeV/#it{c})", Form("#epsilon_{eff,%s}",textMeson.Data()),
                                        kTRUE, 20., 10e-5,kFALSE,
                                        kFALSE, 0., 0.4,
                                        kFALSE, 0., 10.);
                histoEffiCut[i]->DrawCopy("e1,p");
                legendEffiMeson->AddEntry(histoEffiCut[i],Form("standard: %s",cutStringsName[i].Data()));
            } else {
                if(i<20){
                    DrawGammaSetMarker(histoEffiCut[i], 20+i, 1.,color[i],color[i]);
                } else {
                    DrawGammaSetMarker(histoEffiCut[i], 20+i, 1.,color[i-20],color[i-20]);
                }
                histoEffiCut[i]->DrawCopy("same,e1,p");
                legendEffiMeson->AddEntry(histoEffiCut[i],cutStringsName[i].Data());
            }
            
        }
        legendEffiMeson->Draw();
    
        // Efficiency plot labeling
        PutProcessLabelAndEnergyOnPlot( 0.65, 0.2, 0.032, collisionSystem, fNLMString, detectionProcess, 42, 0.03, optionPeriod);
                
        // Draw ratio of efficiencies in lower panel
        padEffiRatios->cd();
        if( optionEnergy.Contains("Pb") ) padEffiRatios->SetLogy(0);
        else padEffiRatios->SetLogy(0);
        for(Int_t i = 0; i< numberOfCuts; i++){
            if(i==0){      
                Double_t minYRatio = 0.45;
                Double_t maxYRatio = 1.55;
                if (cutVariationName.Contains("MultiplicityPP")){
                    minYRatio = 0.5;
                    maxYRatio = 2.7;
                    
                }    
                SetStyleHistoTH1ForGraphs(histoRatioEffiCut[i], "#it{p}_{T} (GeV/#it{c})", "#frac{modified}{standard}", 0.08, 0.11, 0.07, 0.1, 0.75, 0.5, 510,505);
                DrawGammaSetMarker(histoRatioEffiCut[i], 20, 1.,color[0],color[0]);
                histoRatioEffiCut[i]->GetYaxis()->SetRangeUser(minYRatio,maxYRatio);
                
                for(Int_t b = 0; b< histoRatioEffiCut[i]->GetNbinsX(); b++){
                    histoRatioEffiCut[i]->SetBinError(b+1,histoEffiCut[i]->GetBinError(b+1)/histoEffiCut[i]->GetBinContent(b+1));
                }
                histoRatioEffiCut[i]->SetFillColor(kGray+2);
                histoRatioEffiCut[i]->SetFillStyle(0);
                histoRatioEffiCut[i]->DrawCopy("p,e2");  
            } else{
                if(i<20){
                    DrawGammaSetMarker(histoRatioEffiCut[i], 20+i, 1.,color[i],color[i]);
                } else {
                    DrawGammaSetMarker(histoRatioEffiCut[i], 20+i, 1.,color[i-20],color[i-20]);
                }
                histoRatioEffiCut[i]->DrawCopy("same,e1,p");
            }

            DrawGammaLines(0., maxPt,1., 1.,0.1);
        }

        canvasEffiMeson->Update();
        canvasEffiMeson->SaveAs(Form("%s/%s_%s_Efficiencies.%s",outputDir.Data(),meson.Data(),prefix2.Data(),suffix.Data()));
        delete canvasEffiMeson;

        //**************************************************************************************
        //********************* Plotting purities *********************************************
        //**************************************************************************************
        // Define canvas    
        TCanvas* canvasPurityMeson = new TCanvas("canvasPurityMeson","",1350,1500);  // gives the page size
        DrawGammaCanvasSettings( canvasPurityMeson,  0.13, 0.02, 0.02, 0.09);
        // Define upper panel
        TPad* padPurity = new TPad("padPurity", "", 0., 0.25, 1., 1.,-1, -1, -2);
        DrawGammaPadSettings( padPurity, 0.12, 0.02, 0.04, 0.);
        padPurity->Draw();
        // Define lower panel
        TPad* padPurityRatios = new TPad("padPurityRatios", "", 0., 0., 1., 0.25,-1, -1, -2);
        DrawGammaPadSettings( padPurityRatios, 0.12, 0.02, 0.0, 0.2);
        padPurityRatios->Draw();

        // draw efficiency in upper panel
        padPurity->cd();

        TLegend* legendPurityMeson = GetAndSetLegend2(0.15,0.92-1.15*0.032*numberOfCuts,0.3,0.92, 1500*0.75*0.032); 
        for(Int_t i = 0; i< numberOfCuts; i++){
            if(i == 0){
                DrawGammaSetMarker(histoPurityCut[i], 20, 1., color[0], color[0]);
                DrawAutoGammaMesonHistos( histoPurityCut[i],
                                        "", "#it{p}_{T} (GeV/#it{c})", Form("#epsilon_{pur,%s}",textMeson.Data()),
                                        kFALSE, 5., 10e-10,kFALSE,
                                        kTRUE, 0.4, 1.05,
                                        kFALSE, 0., 10.);
                histoPurityCut[i]->DrawCopy("e1,p");
                legendPurityMeson->AddEntry(histoPurityCut[i],Form("standard: %s",cutStringsName[i].Data()));
            } else {
                if(i<20){
                    DrawGammaSetMarker(histoPurityCut[i], 20+i, 1.,color[i],color[i]);
                } else {
                    DrawGammaSetMarker(histoPurityCut[i], 20+i, 1.,color[i-20],color[i-20]);
                }
                histoPurityCut[i]->DrawCopy("same,e1,p");
                legendPurityMeson->AddEntry(histoPurityCut[i],cutStringsName[i].Data());
            }
            
        }
        legendPurityMeson->Draw();
    
        // Purity plot labeling
        PutProcessLabelAndEnergyOnPlot( 0.65, 0.2, 0.032, collisionSystem, fNLMString, detectionProcess, 42, 0.03, optionPeriod);
                
        // Draw ratio of purities in lower panel
        padPurityRatios->cd();
        if( optionEnergy.Contains("Pb") ) padPurityRatios->SetLogy(0);
        else padPurityRatios->SetLogy(0);
        for(Int_t i = 0; i< numberOfCuts; i++){
            if(i==0){      
                Double_t minYRatio = 0.45;
                Double_t maxYRatio = 1.55;
                if (cutVariationName.Contains("MultiplicityPP")){
                    minYRatio = 0.5;
                    maxYRatio = 2.7;
                    
                }    
                SetStyleHistoTH1ForGraphs(histoRatioPurityCut[i], "#it{p}_{T} (GeV/#it{c})", "#frac{modified}{standard}", 0.08, 0.11, 0.07, 0.1, 0.75, 0.5, 510,505);
                DrawGammaSetMarker(histoRatioPurityCut[i], 20, 1.,color[0],color[0]);
                histoRatioPurityCut[i]->GetYaxis()->SetRangeUser(minYRatio,maxYRatio);
                
                for(Int_t b = 0; b< histoRatioPurityCut[i]->GetNbinsX(); b++){
                    histoRatioPurityCut[i]->SetBinError(b+1,histoPurityCut[i]->GetBinError(b+1)/histoPurityCut[i]->GetBinContent(b+1));
                }
                histoRatioPurityCut[i]->SetFillColor(kGray+2);
                histoRatioPurityCut[i]->SetFillStyle(0);
                histoRatioPurityCut[i]->DrawCopy("p,e2");  
            } else{
                if(i<20){
                    DrawGammaSetMarker(histoRatioPurityCut[i], 20+i, 1.,color[i],color[i]);
                } else {
                    DrawGammaSetMarker(histoRatioPurityCut[i], 20+i, 1.,color[i-20],color[i-20]);
                }
                histoRatioPurityCut[i]->DrawCopy("same,e1,p");
            }

            DrawGammaLines(0., maxPt,1., 1.,0.1);
        }

        canvasPurityMeson->Update();
        canvasPurityMeson->SaveAs(Form("%s/%s_%s_Purity.%s",outputDir.Data(),meson.Data(),prefix2.Data(),suffix.Data()));
        delete canvasPurityMeson;
        

        //**************************************************************************************
        //********************* Plotting Acceptance *********************************************
        //**************************************************************************************
        // Define canvas    
        TCanvas* canvasAcceptanceMeson = new TCanvas("canvasAcceptanceMeson","",1350,1500);  // gives the page size
        DrawGammaCanvasSettings( canvasAcceptanceMeson,  0.13, 0.02, 0.02, 0.09);
        // Define upper panel
        TPad* padAcceptance = new TPad("padAcceptance", "", 0., 0.25, 1., 1.,-1, -1, -2);
        DrawGammaPadSettings( padAcceptance, 0.12, 0.02, 0.04, 0.);
        padAcceptance->Draw();
        // Define lower panel
        TPad* padAcceptanceRatios = new TPad("padAcceptanceRatios", "", 0., 0., 1., 0.25,-1, -1, -2);
        DrawGammaPadSettings( padAcceptanceRatios, 0.12, 0.02, 0.0, 0.2);
        padAcceptanceRatios->Draw();

        // draw acceptance in upper panel
        padAcceptance->cd();

        TLegend* legendAcceptMeson = GetAndSetLegend2(0.15,0.92-1.15*0.032*numberOfCuts,0.3,0.92, 1500*0.75*0.032);
        for(Int_t i = 0; i< numberOfCuts; i++){
            if(i == 0){
                DrawGammaSetMarker(histoAcceptanceCut[i], 20, 1., color[0], color[0]);
                DrawAutoGammaMesonHistos( histoAcceptanceCut[i],
                                        "", "#it{p}_{T} (GeV/#it{c})", Form("A_{%s}",textMeson.Data()),
                                        kFALSE, 5., 10e-10,kFALSE,
                                        kFALSE, 0.15, 0.30,
                                        kFALSE, 0., 10.);
                histoAcceptanceCut[i]->GetYaxis()->SetRangeUser(0.15, 0.3);
                histoAcceptanceCut[i]->Draw("e1,p");
                legendAcceptMeson->AddEntry(histoAcceptanceCut[i],Form("standard: %s",cutStringsName[i].Data()));
                cout << "I plotted correctly" << endl;
            } else {
                if(i<20){
                    DrawGammaSetMarker(histoAcceptanceCut[i], 20+i, 1.,color[i],color[i]);
                } else {
                    DrawGammaSetMarker(histoAcceptanceCut[i], 20+i, 1.,color[i-20],color[i-20]);
                }
                histoAcceptanceCut[i]->DrawCopy("same,e1,p");
                legendAcceptMeson->AddEntry(histoAcceptanceCut[i],cutStringsName[i].Data());
            }
            
        }
        legendAcceptMeson->Draw();
    
        // Acceptance plot labeling
        PutProcessLabelAndEnergyOnPlot( 0.65, 0.20, 0.032, collisionSystem, fNLMString, detectionProcess, 42, 0.03, optionPeriod);
                
        // Draw ratio of acceptance in lower panel
        padAcceptanceRatios->cd();
        if( optionEnergy.Contains("Pb") ) padAcceptanceRatios->SetLogy(0);
        else padAcceptanceRatios->SetLogy(0);
        for(Int_t i = 0; i< numberOfCuts; i++){
            if(i==0){      
                Double_t minYRatio = 0.85;
                Double_t maxYRatio = 1.15;
                if( optionEnergy.Contains("Pb")){
                    minYRatio = 0.05;        maxYRatio = 2.4;
                } 
                SetStyleHistoTH1ForGraphs(histoRatioAcceptanceCut[i], "#it{p}_{T} (GeV/#it{c})", "#frac{modified}{standard}", 0.08, 0.11, 0.07, 0.1, 0.75, 0.5, 510,505);
                DrawGammaSetMarker(histoRatioAcceptanceCut[i], 20, 1.,color[0],color[0]);
                histoRatioAcceptanceCut[i]->GetYaxis()->SetRangeUser(minYRatio,maxYRatio);
                histoRatioAcceptanceCut[i]->DrawCopy("p,e1");
            } else{
                if(i<20){
                    DrawGammaSetMarker(histoRatioAcceptanceCut[i], 20+i, 1.,color[i],color[i]);
                } else {
                    DrawGammaSetMarker(histoRatioAcceptanceCut[i], 20+i, 1.,color[i-20],color[i-20]);
                }
                histoRatioAcceptanceCut[i]->DrawCopy("same,e1,p");
            }

            DrawGammaLines(0., maxPt,1., 1.,0.1);
        }

        canvasAcceptanceMeson->Update();
        canvasAcceptanceMeson->SaveAs(Form("%s/%s_%s_Acceptance.%s",outputDir.Data(),meson.Data(),prefix2.Data(),suffix.Data()));
        delete canvasAcceptanceMeson;
        
        //*************************************************************************************************
        //******************** Output of the systematic Error due to Signal extraction for Meson ************
        //*************************************************************************************************
        // Determine number of bins
        Int_t NBinsPt = histoCorrectedYieldCut[0]->GetNbinsX();
        const Int_t NBinstPtConst = NBinsPt+1;
        
        // Create array of bin boundaries
        Double_t  BinsXCenter[NBinstPtConst];
        Double_t  BinsXWidth[NBinstPtConst];
        BinsXCenter[0]                      = 0;
        BinsXWidth[0]                       = 0.;
        for (Int_t i = 1; i < NBinsPt +1; i++){
            BinsXCenter[i]                  = histoCorrectedYieldCut[0]->GetBinCenter(i);
            BinsXWidth[i]                   = histoCorrectedYieldCut[0]->GetBinWidth(i)/2.;
        }

        // Create array of Sys Err Objects and fill them
        SysErrorConversion SysErrCut    [constNumberOfCuts][NBinstPtConst];
        SysErrorConversion SysErrCutRaw [constNumberOfCuts][NBinstPtConst];
        for (Int_t j = 0; j < numberOfCuts; j++){
            for (Int_t i = 1; i < NBinsPt +1; i++){
                SysErrCut[j][i].value       = histoCorrectedYieldCut[j]->GetBinContent(i);
                SysErrCut[j][i].error       = histoCorrectedYieldCut[j]->GetBinError(i);
                SysErrCutRaw[j][i].value    = histoRawYieldCut[j]->GetBinContent(i);
                SysErrCutRaw[j][i].error    = histoRawYieldCut[j]->GetBinError(i);
            }
        }

        // Create Difference arrays
        Double_t DifferenceCut          [constNumberOfCuts][NBinstPtConst];
        Double_t DifferenceErrorCut     [constNumberOfCuts][NBinstPtConst];
        Double_t RelDifferenceCut       [constNumberOfCuts][NBinstPtConst];
        Double_t RelDifferenceErrorCut  [constNumberOfCuts][NBinstPtConst];
        Double_t RelDifferenceRawCut    [constNumberOfCuts][NBinstPtConst];
            
        // Create largest difference array
        Double_t LargestDiffNeg         [NBinstPtConst];
        Double_t LargestDiffPos         [NBinstPtConst];
        Double_t LargestDiffErrorNeg    [NBinstPtConst];
        Double_t LargestDiffErrorPos    [NBinstPtConst];
        Double_t LargestDiffRelNeg      [NBinstPtConst];
        Double_t LargestDiffRelPos      [NBinstPtConst];
        Double_t LargestDiffRelErrorNeg [NBinstPtConst];
        Double_t LargestDiffRelErrorPos [NBinstPtConst];

        // Initialize all differences with 0
        for (Int_t j = 0; j < numberOfCuts; j++){
            for ( Int_t i = 0; i < NBinstPtConst; i++) {
                DifferenceCut[j][i]         = 0.;
                DifferenceErrorCut[j][i]    = 0.;
                LargestDiffNeg[i]           = 0.;
                LargestDiffPos[i]           = 0.;
                LargestDiffErrorNeg[i]      = 0.;
                LargestDiffErrorPos[i]      = 0.;
                RelDifferenceCut[j][i]      = 0.;
                RelDifferenceRawCut[j][i]   = 0.;
                RelDifferenceErrorCut[j][i] = 0.;
            }
        }


        // Calculate largest difference among cut variation 
        for(Int_t j = 1; j < numberOfCuts; j++){
            for (Int_t i = 0; i < NBinsPt +1; i++){
                // Calculate difference (rel/abs) and error for corrected yield
                DifferenceCut[j][i]         = SysErrCut[j][i].value - SysErrCut[0][i].value;
                DifferenceErrorCut[j][i]    = TMath::Sqrt(TMath::Abs(TMath::Power(SysErrCut[j][i].error,2)-TMath::Power(SysErrCut[0][i].error,2)));
                if(SysErrCut[0][i].value != 0){
                    RelDifferenceCut[j][i]      = DifferenceCut[j][i]/SysErrCut[0][i].value*100. ;
                    RelDifferenceErrorCut[j][i] = DifferenceErrorCut[j][i]/SysErrCut[0][i].value*100. ;
                } else {
                    RelDifferenceCut[j][i]      = -10000.;
                    RelDifferenceErrorCut[j][i] = 100. ;
                }
                // Calculate relativ difference for raw yield
                if(SysErrCutRaw[0][i].value != 0){
                    RelDifferenceRawCut[j][i]   = (SysErrCutRaw[j][i].value - SysErrCutRaw[0][i].value)/SysErrCutRaw[0][i].value*100. ;
                } else {
                    RelDifferenceRawCut[j][i]   = -10000.;
                }
                if (i == 0){
                    RelDifferenceRawCut[j][i]   = 0.;
                    RelDifferenceCut[j][i]      = 0. ;
                    RelDifferenceErrorCut[j][i] = 0. ;
                    DifferenceCut[j][i]         = 0.;
                    DifferenceErrorCut[j][i]    = 0.;
                }  
                    
                if(doBarlow){ 
                    // !!! => Careful, this is meant to be a cross check. If it has to be used for the syst errors
                    // the syste error macros has to be changed accordingly (the mean of pos and  neg dev cannot be used)
                    
                    // Calculate largest differences in positiv and negative direction
                    if(DifferenceCut[j][i] < 0){ // largest negativ deviation
                        // Take deviation if larger than previous largest deviation 
                        // and relative raw yield loss less than 75%
                        if (TMath::Abs(LargestDiffNeg[i]) < TMath::Abs(DifferenceCut[j][i]) && RelDifferenceRawCut[j][i] > -75.){
                            if( TMath::Abs(DifferenceCut[j][i])/TMath::Abs(DifferenceErrorCut[j][i]) > 1. ){
                                LargestDiffNeg[i]       = DifferenceCut[j][i];
                                LargestDiffErrorNeg[i]  = DifferenceErrorCut[j][i];
                            } else if( TMath::Abs(DifferenceCut[j][i])/TMath::Abs(DifferenceErrorCut[j][i]) < 1.){
                                cout << "Largest negative difference not updated" << endl;
                            }
                        }
                    } else { // largest positive deviation
                        // Take deviation if larger than previous largest deviation 
                        // and relative raw yield loss less than 75%
                        if (TMath::Abs(LargestDiffPos[i]) < TMath::Abs(DifferenceCut[j][i]) && RelDifferenceRawCut[j][i] > -75.){
                            if( TMath::Abs(DifferenceCut[j][i])/TMath::Abs(DifferenceErrorCut[j][i]) > 1. ){
                                LargestDiffPos[i]       = DifferenceCut[j][i];
                                LargestDiffErrorPos[i]  = DifferenceErrorCut[j][i];
                            } else if( TMath::Abs(DifferenceCut[j][i])/TMath::Abs(DifferenceErrorCut[j][i]) < 1.){
                                cout << "Largest positive difference not updated" << endl;
                            }
                        }
                    }
                } else {
                    
                    // Calculate largest differences in positiv and negative direction
                    if(DifferenceCut[j][i] < 0){ // largest negativ deviation
                        // Take deviation if larger than previous largest deviation 
                        // and relative raw yield loss less than 75%
                        if (TMath::Abs(LargestDiffNeg[i]) < TMath::Abs(DifferenceCut[j][i]) && RelDifferenceRawCut[j][i] > -75.){
                            LargestDiffNeg[i]       = DifferenceCut[j][i];
                            LargestDiffErrorNeg[i]  = DifferenceErrorCut[j][i];
                        }
                    } else { // largest positive deviation
                        // Take deviation if larger than previous largest deviation 
                        // and relative raw yield loss less than 75%
                        if (TMath::Abs(LargestDiffPos[i]) < TMath::Abs(DifferenceCut[j][i]) && RelDifferenceRawCut[j][i] > -75.){
                            LargestDiffPos[i]       = DifferenceCut[j][i];
                            LargestDiffErrorPos[i]  = DifferenceErrorCut[j][i];
                        }
                    }
                    if (i == 0){
                        LargestDiffPos[i]       = 0;
                        LargestDiffErrorPos[i]  = 0;
                        LargestDiffNeg[i]       = 0;
                        LargestDiffErrorNeg[i]  = 0;                      
                    }
                }
            }
        }
        
        if(doBarlow){ 
            TString SysErrCheckname = Form("%s/%s_%s_SystematicErrorBarlowCheck.dat",outputDir.Data(),meson.Data(),prefix2.Data());
            fstream SysErrDatCheck;
            SysErrDatCheck.open(SysErrCheckname.Data(), ios::out);
            SysErrDatCheck << "Barlow check for the systematic error" << endl;
            for (Int_t l=0; l< numberOfCuts; l++){
                if (l == 0) {
                    SysErrDatCheck << endl <<"Bin" << "\t" << cutNumber[l] << "\t" <<endl;
                    for(Int_t i = 1; i < (NBinsPt +1); i++){
                        SysErrDatCheck << BinsXCenter[i] << "\t" << SysErrCut[l][i].value << "\t" << SysErrCut[l][i].error << endl;    
                    }
                } else{
                    for(Int_t i = 1; i < (NBinsPt +1); i++){
                        SysErrDatCheck << endl <<"Cut" << "\t" << cutNumber[l] << "\t" <<endl;
                        SysErrDatCheck << "\t Barlow check for " << BinsXCenter[i] << endl;
                        SysErrDatCheck << "Delta = |a_1 - a_2| \t" << TMath::Abs(DifferenceCut[l][i]) << endl;
                        SysErrDatCheck << "sigma_Delta = sqrt( |sigma_2^2 - sigma_1^2| ) \t" << TMath::Abs(DifferenceErrorCut[l][i]) << endl;
                        SysErrDatCheck << "Check: Delta/sigma_Delta < 1? " << (TMath::Abs(DifferenceCut[l][i]))/(TMath::Abs(DifferenceErrorCut[l][i])) << endl;
                    }
                }
            }            
            SysErrDatCheck.close();
        }                        
                    
        // Write systematic error input to log file
        TString SysErrDatname   = Form("%s/%s_%s_SystematicErrorCutStudies.dat",outputDir.Data(),meson.Data(),prefix2.Data());
        fstream SysErrDat;
        SysErrDat.open(SysErrDatname.Data(), ios::out);
        SysErrDat << "Calculation of the systematic error due to the yield cuts" << endl;

        cout << "works" << endl;
        for (Int_t l=0; l< numberOfCuts; l++){
            if (l == 0) {
                SysErrDat << endl <<"Bin" << "\t" << cutNumber[l] << "\t" <<endl;
                for(Int_t i = 1; i < (NBinsPt +1); i++){
                SysErrDat << BinsXCenter[i] << "\t" << SysErrCut[l][i].value << "\t" << SysErrCut[l][i].error << endl;    
                }
            } else{
                SysErrDat << endl <<"Bin" << "\t" << cutNumber[l] << "\t" << "Error " << "\t Dif to Cut1" << endl;
                for(Int_t i = 1; i < (NBinsPt +1); i++){
                    if (RelDifferenceRawCut[l][i] > -75.){
                        SysErrDat << BinsXCenter[i] << "\t" << SysErrCut[l][i].value << "\t" << SysErrCut[l][i].error << "\t" <<  DifferenceCut[l][i] << "\t"<< DifferenceErrorCut[l][i] << "\t" <<
                        RelDifferenceCut[l][i] <<  "\t" << RelDifferenceErrorCut[l][i] <<"\t" << RelDifferenceRawCut[l][i]<< endl;
                    } else {
                        SysErrDat << BinsXCenter[i] << "\t" << SysErrCut[l][i].value << "\t" << SysErrCut[l][i].error << "\t" <<  DifferenceCut[l][i] << "\t"<< DifferenceErrorCut[l][i] << "\t" << 
                        RelDifferenceCut[l][i] <<  "\t" << RelDifferenceErrorCut[l][i] <<"\t" << RelDifferenceRawCut[l][i]  <<"\t not considered in largest dev" <<endl;
                    }
                }
            }
        }
        SysErrDat << endl;
        SysErrDat << endl;
        SysErrDat << "Bin" << "\t" << "Largest Dev Neg" << "\t" << "Largest Dev Pos"  << endl;
        for(Int_t i = 1; i < (NBinsPt +1); i++){
            SysErrDat << BinsXCenter[i]  << "\t" << LargestDiffNeg[i] << "\t" <<LargestDiffErrorNeg[i]<< "\t" << LargestDiffPos[i] << "\t" << LargestDiffErrorPos[i]<<endl;
        }
        SysErrDat << endl << endl <<"Bin" << "\t" << "Largest Dev Neg rel" << "\t" << "Largest Dev Pos rel"  << endl;
        // Calculate largest relative deviations
        for(Int_t i = 0; i < (NBinsPt +1); i++){
            if ( SysErrCut[0][i].value != 0.){
                LargestDiffRelNeg[i]        = - LargestDiffNeg[i]/SysErrCut[0][i].value*100.;
                LargestDiffRelPos[i]        = LargestDiffPos[i]/SysErrCut[0][i].value*100.;
                LargestDiffRelErrorNeg[i]   = - LargestDiffErrorNeg[i]/SysErrCut[0][i].value*100.;
                LargestDiffRelErrorPos[i]   = LargestDiffErrorPos[i]/SysErrCut[0][i].value*100.;
                if (i > 0){
                    SysErrDat << BinsXCenter[i] << "\t" << LargestDiffNeg[i]/SysErrCut[0][i].value*100. << "\t" << LargestDiffErrorNeg[i]/SysErrCut[0][i].value*100. << "\t" <<     
                    LargestDiffPos[i]/SysErrCut[0][i].value*100. << "\t" << LargestDiffErrorPos[i]/SysErrCut[0][i].value*100.<<endl;
                } else {
                    LargestDiffRelNeg[i]        = 0.;
                    LargestDiffRelPos[i]        = 0.;
                    LargestDiffRelErrorNeg[i]   = 0.;
                    LargestDiffRelErrorPos[i]   = 0.;
                }  
            } else {
                LargestDiffRelNeg[i]        = 0.;
                LargestDiffRelPos[i]        = 0.;
                LargestDiffRelErrorNeg[i]   = 0.;
                LargestDiffRelErrorPos[i]   = 0.;
            }
        }
        SysErrDat.close();

        // Create sys-err graphs
        TGraphAsymmErrors* SystErrGraphNeg  = new TGraphAsymmErrors(NBinsPt+1, BinsXCenter, LargestDiffRelNeg, BinsXWidth, BinsXWidth, LargestDiffRelErrorNeg, LargestDiffRelErrorNeg);
        SystErrGraphNeg->SetName(Form("%s_SystErrorRelNeg_%s",meson.Data(),cutVariationName.Data()));
        cout << "Negative error graph" << endl; 
        SystErrGraphNeg->Print();
        TGraphAsymmErrors* SystErrGraphPos  = new TGraphAsymmErrors(NBinsPt+1, BinsXCenter, LargestDiffRelPos, BinsXWidth, BinsXWidth, LargestDiffRelErrorPos, LargestDiffRelErrorPos);
        SystErrGraphPos->SetName(Form("%s_SystErrorRelPos_%s",meson.Data(),cutVariationName.Data()));
        cout << "positive error graph" << endl;
        SystErrGraphPos->Print();
        
        // Write sys-err graph to root output file
        TString Outputname                  = Form("%s/%s_%s_SystematicErrorCuts.root",outputDirRootFile.Data(),meson.Data(),prefix2.Data());
        TFile* SystematicErrorFile          = new TFile(Outputname.Data(),"UPDATE");
            SystErrGraphPos->Write(Form("%s_SystErrorRelPos_%s",meson.Data(),cutVariationName.Data()),TObject::kOverwrite);
            SystErrGraphNeg->Write(Form("%s_SystErrorRelNeg_%s",meson.Data(),cutVariationName.Data()),TObject::kOverwrite);
        SystematicErrorFile->Write();
        SystematicErrorFile->Close();
    }
    return;
}
