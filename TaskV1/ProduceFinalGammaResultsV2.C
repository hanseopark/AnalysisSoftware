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
#include "TBenchmark.h"
#include "TRandom.h"
#include "TLatex.h"
#include "TASImage.h"
#include "TPostScript.h"
#include "TGraphErrors.h"
#include "TArrow.h"
#include "TGaxis.h"
#include "TMarker.h"
#include "TPaveText.h"
#include "TGraphAsymmErrors.h"
#include "../CommonHeaders/PlottingGammaConversionHistos.h"
#include "../CommonHeaders/PlottingGammaConversionAdditional.h"
#include "../CommonHeaders/FittingGammaConversion.h"
#include "../CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "../CommonHeaders/ConversionFunctions.h"

void ProduceFinalGammaResultsV2(    TString cutSel          = "50100013_00200009847005008750404000_0152501500000000",
                                    TString optionEnergy    = "PbPb_2.76TeV",
                                    TString suffix          = "eps",
                                    Int_t mode              = 0
                                ){

    StyleSettingsThesis();
    SetPlotStyle();

    TString dateForOutput           = ReturnDateStringForOutput();
    TString collisionSystem         = ReturnFullCollisionsSystem(optionEnergy);
    TString collisionSystemOutput   = ReturnCollisionEnergyOutputString(optionEnergy);
    TString centrality              = GetCentralityString(cutSel);
    TString centralityW0Per         = GetCentralityStringWoPer(cutSel);
    if (centrality.CompareTo("pp") == 0 || centrality.CompareTo("0-100%") == 0 ) {
        centrality                  = "";
        centralityW0Per             = "";
    }    
    TString detectionProcess        = ReturnFullTextReconstructionProcess(mode);

    TString system                  = "PCM";
    if (mode == 2) system           = "PCM-EMCAL";
    if (mode == 3) system           = "PCM-PHOS";
    if (mode == 4) system           = "EMCAL-EMCAL";
    if (mode == 5) system           = "PHOS-PHOS";
    if (mode == 10) system          = "EMC-merged";
    if (mode == 11) system          = "PHOS-merged";
    
    // defining output directory
    TString outputDir       = Form("%s/%s/FinalGammaResults_%s_%s%s", suffix.Data(), dateForOutput.Data(), system.Data(),centralityW0Per.Data(), collisionSystemOutput.Data() );
    gSystem->Exec("mkdir -p "+outputDir);

    //*************************************************************************************************
    //******************** read from input files ******************************************************
    //*************************************************************************************************
    TString inputFileName               = Form("%s/%s/Gamma_Pi0_data_GammaConvV1_InclusiveRatio.root", cutSel.Data(), optionEnergy.Data());
    cout << "trying to read: " << inputFileName.Data() << endl;
    TFile *fileInput                    = new TFile(inputFileName.Data());
    if (fileInput->IsZombie()) {
        cout << "file couldn't be read, aborting....";
        return;
    }    
    
    // read stat error hists
    TH1D* histoIncGamma                     = (TH1D*) fileInput->Get("histoGammaSpecCorrPurity");
    TH1D* histoPi0Spectrum                  = (TH1D*) fileInput->Get("CorrectedYieldTrueEff");
    TH1D* histoIncRatio                     = (TH1D*) fileInput->Get("IncRatioPurity_trueEff");
    TH1D* histoIncRatioPi0Fit               = (TH1D*) fileInput->Get("histoIncRatioFitPurity");
    TH1D* histoDR                           = (TH1D*) fileInput->Get("DoubleRatioTrueEffPurity");
    TH1D* histoDRFit                        = (TH1D*) fileInput->Get("DoubleRatioFitPurity");
    TH1D* histoDirGamma                     = (TH1D*) fileInput->Get("histoDirectPhotonSpectrum");
    // read sys error graphs
    TGraphAsymmErrors* graphIncGammaSysErr      = (TGraphAsymmErrors*) fileInput->Get("histoGammaSpecCorrPurity_SystErr");
    TGraphAsymmErrors* graphIncRatioSysErr      = (TGraphAsymmErrors*) fileInput->Get("IncRatioPurity_trueEff_SystErr");
    TGraphAsymmErrors* graphIncRatioPi0FitSysErr= (TGraphAsymmErrors*) fileInput->Get("histoIncRatioFitPurity_SystErr");
    TGraphAsymmErrors* graphDRSysErr            = (TGraphAsymmErrors*) fileInput->Get("DoubleRatioTrueEffPurity_SystErr");
    TGraphAsymmErrors* graphDRPi0FitSysErr      = (TGraphAsymmErrors*) fileInput->Get("DoubleRatioFitPurity_SystErr");
    // read theory graphs
    TGraphAsymmErrors* graphNLODR               = (TGraphAsymmErrors*) fileInput->Get("graphNLODoubleRatio");
    TGraphAsymmErrors* graphNLOGammaDir         = (TGraphAsymmErrors*) fileInput->Get("graphNLODirGamma");
    TGraphAsymmErrors* graphNLOGammaPrompt      = (TGraphAsymmErrors*) fileInput->Get("graphPromptPhotonNLO");
    TGraphAsymmErrors* graphNLOGammaFrag        = (TGraphAsymmErrors*) fileInput->Get("graphFragmentationPhotonNLO");
    
    
    //***************************************************************************************************
    //***************************************************************************************************
    //***************************************************************************************************
    Int_t nLinesNLOLegends  = 2;
    if (optionEnergy.CompareTo("PbPb_2.76TeV") == 0 || optionEnergy.CompareTo("pPb_5.023TeV") == 0) 
        nLinesNLOLegends    = 3;
    Double_t minPt                              = 0.2;
    Double_t maxPtLog                        = 40.;
    Double_t minY                            = 0.85;
    Double_t maxY                            = 1.65;
    if(mode == 0 && optionEnergy.CompareTo("PbPb_2.76TeV") == 0){
        maxY                                 = 2.;
        minY                                 = 0.75;        
    } else if (mode==4) {
        minPt                                = 1.5;
        minY                                 = 0.5;
    }
    
    TCanvas *canvasDoublRatioNLO = new TCanvas("canvasDoublRatioNLO","",0.095,0.09,1000,815);
    DrawGammaCanvasSettings( canvasDoublRatioNLO, 0.09, 0.02, 0.02, 0.09);
    canvasDoublRatioNLO->cd();
    canvasDoublRatioNLO->SetLogx();
            
    TH2F * histo2DDoubleRatioPlotting       = new TH2F("histo2DDoubleRatioPlotting","histo2DDoubleRatioPlotting",1000,minPt,maxPtLog,1000,minY,maxY);
    SetStyleHistoTH2ForGraphs(histo2DDoubleRatioPlotting, "#it{p}_{T} (GeV/#it{c})","(#gamma_{inc}/#pi^{0})/(#gamma_{decay}/#pi^{0})", 0.04,0.04, 0.04,0.04, 1.,1.);
    if (optionEnergy.CompareTo("PbPb_2.76TeV") == 0)
        histo2DDoubleRatioPlotting->GetXaxis()->SetRangeUser(0.4,histoDR->GetXaxis()->GetBinUpEdge(histoDR->GetNbinsX())*1.5);
    else
        histo2DDoubleRatioPlotting->GetXaxis()->SetRangeUser(0.,histoDR->GetXaxis()->GetBinUpEdge(histoDR->GetNbinsX())*1.5);
    histo2DDoubleRatioPlotting->GetXaxis()->SetTitleFont(62);
    histo2DDoubleRatioPlotting->GetYaxis()->SetTitleFont(62);
    histo2DDoubleRatioPlotting->GetXaxis()->SetLabelOffset(-1e-2);
    histo2DDoubleRatioPlotting->DrawCopy();
 
        histo2DDoubleRatioPlotting->DrawCopy();
            DrawGammaLines(0., histoDR->GetXaxis()->GetBinUpEdge(histoDR->GetNbinsX())*1.5, 1.0, 1.0,2.0, kGray+2 ,7);
            graphNLODR->Draw("lp3");

            DrawGammaSetMarker(histoDR, 20, 2.0, kBlack, kBlack);
            DrawGammaSetMarker(histoDRFit, 20, 2.0, kBlue+2, kBlue+2);
            histoDR->DrawCopy("same");
            histoDRFit->DrawCopy("same");
            if (graphDRSysErr)     graphDRSysErr->Draw("p,2,same");
            if (graphDRPi0FitSysErr)  graphDRPi0FitSysErr->Draw("p,2,same");
            
            TLegend* legendDoubleConversionFitNLO = GetAndSetLegend2(0.14,0.93-(nLinesNLOLegends+1+0.5)*0.045,0.5,0.93,0.045,1,"",42,0.2); 
            legendDoubleConversionFitNLO->AddEntry(histoDR,"Data","p");
            legendDoubleConversionFitNLO->AddEntry(histoDRFit,"Data, fitted #pi^{0}","p");
            if(optionEnergy.CompareTo("PbPb_2.76TeV") == 0)       legendDoubleConversionFitNLO->AddEntry(graphNLODR,"pp #gamma_{dir} NLO #sqrt{s} = 2.76 TeV","l");
            else if(optionEnergy.CompareTo("pPb_5.023TeV") == 0)  legendDoubleConversionFitNLO->AddEntry(graphNLODR,"pp #gamma_{dir} NLO #sqrt{s} = 5.02 TeV ","l");
            else                                             legendDoubleConversionFitNLO->AddEntry(graphNLODR,"pp NLO Direct Photon","l");
            if (optionEnergy.CompareTo("PbPb_2.76TeV") == 0 || optionEnergy.CompareTo("pPb_5.023TeV") == 0) legendDoubleConversionFitNLO->AddEntry((TObject*)0,"scaled by N_{coll}","");

            legendDoubleConversionFitNLO->Draw();

//             PlotLatexLegend(0.93, 0.96-3*0.045, 0.045,collisionSystem,detectionProcess,3,31);
        
        canvasDoublRatioNLO->Print(Form("%s/DoubleRatioComparison_NLO.%s",outputDir.Data(),suffix.Data()));

    
    //*************************************************************************************************
    // put everything in common output per system
    //*************************************************************************************************
    TString optionOutput                = "pp";
    TString fileNameOutputComp          = Form("%s_%sResultsFullCorrection_PP.root","data",system.Data());
    if (optionEnergy.Contains("pPb")){
        fileNameOutputComp              = Form("%s_%sResultsFullCorrection_pPb.root","data",system.Data());
        optionOutput                    = "pPb";    
    } else if (optionEnergy.Contains("PbPb")){
        fileNameOutputComp              = Form("%s_%sResultsFullCorrection_PbPb.root","data",system.Data());
        optionOutput                    = "PbPb";    
    }
    TFile* fileGammaFinal               = new TFile(fileNameOutputComp,"UPDATE");
    
        // create subdirectory for respective energy
        TDirectoryFile* directoryGamma      = (TDirectoryFile*)fileGammaFinal->Get(Form("Gamma_%s%s",collisionSystemOutput.Data(), centrality.Data() )); 
        if (!directoryGamma){
            fileGammaFinal->mkdir(Form("Gamma_%s%s",collisionSystemOutput.Data(), centrality.Data() ));
            directoryGamma      = (TDirectoryFile*)fileGammaFinal->Get(Form("Gamma_%s%s",collisionSystemOutput.Data(), centrality.Data() )); 
            fileGammaFinal->cd(Form("Gamma_%s%s",collisionSystemOutput.Data(), centrality.Data() ));
        } else {
            fileGammaFinal->cd(Form("Gamma_%s%s",collisionSystemOutput.Data(), centrality.Data() ));
        }
    
            // writing double ratio quantities
            if (histoDR){
                SetHistogramm(histoDR,"#it{p}_{T} (GeV/#it{c})", "(#it{N}_{#gamma_{inc}}/#it{N}_{#pi^{0}})/(#it{N}_{#gamma_{decay}}/#it{N}_{#pi^{0}})");
                histoDR->Write("DoubleRatioStatError",TObject::kOverwrite);
            }    
            if (graphDRSysErr) graphDRSysErr->Write("DoubleRatioSystError",TObject::kOverwrite);

            if (histoDRFit){
                SetHistogramm(histoDRFit,"#it{p}_{T} (GeV/#it{c})", "(#it{N}_{#gamma_{inc}}/#it{N}_{#pi^{0}})/(#it{N}_{#gamma_{decay}}/#it{N}_{#pi^{0}})");
                histoDRFit->Write("DoubleRatioPi0FitStatError",TObject::kOverwrite);
            }
            if(graphDRPi0FitSysErr) graphDRPi0FitSysErr->Write("DoubleRatioPi0FitSystError",TObject::kOverwrite);

            // writing inclusive ratio quantities
            if (histoIncRatio){
                SetHistogramm(histoIncRatio,"#it{p}_{T} (GeV/#it{c})", "#gamma_{inc}/#pi^{0}");
                histoIncRatio->Write("IncRatioStatError",TObject::kOverwrite);
            }    
            if (graphIncRatioSysErr) graphIncRatioSysErr->Write("IncRatioSystError",TObject::kOverwrite);

            if (histoIncRatioPi0Fit){
                SetHistogramm(histoIncRatioPi0Fit,"#it{p}_{T} (GeV/#it{c})", "#gamma_{inc}/#pi^{0}");
                histoIncRatioPi0Fit->Write("IncRatioPi0FitStatError",TObject::kOverwrite);
            }
            if(graphIncRatioPi0FitSysErr) graphIncRatioPi0FitSysErr->Write("IncRatioPi0FitSystError",TObject::kOverwrite);

            // writing inclusive gamma spectrum
            if (histoIncGamma){
                SetHistogramm(histoIncGamma,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}N_{#gamma}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c}^{2})");
                histoIncGamma->Write("IncGammaStatError",TObject::kOverwrite);
            }    
            if (graphIncGammaSysErr) graphIncGammaSysErr->Write("IncGammaSystError",TObject::kOverwrite);

            // writing pi0 used pi0 spectrum
            if (histoPi0Spectrum){
                SetHistogramm(histoPi0Spectrum,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}N_{#pi^{0}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c}^{2})");
                histoPi0Spectrum->Write("Pi0StatError",TObject::kOverwrite);
            }
            // writing direct photon spectrum
            if(histoDirGamma)            histoDirGamma->Write("DirGammaUpperLimits",TObject::kOverwrite);
                
//             if (histoCocktailAllGamma)histoCocktailAllGamma->Write("CocktailSumGamma",TObject::kOverwrite);
//             if (histoCocktailPi0Gamma)histoCocktailPi0Gamma->Write("CocktailPi0Gamma",TObject::kOverwrite);
//             if (histoCocktailEtaGamma)histoCocktailEtaGamma->Write("CocktailEtaGamma",TObject::kOverwrite);
//             if (histoCocktailOmegaGamma)histoCocktailOmegaGamma->Write("CocktailOmegaGamma",TObject::kOverwrite);
//             if (histoCocktailEtapGamma)histoCocktailEtapGamma->Write("CocktailEtapGamma",TObject::kOverwrite);
//             if (histoCocktailPhiGamma)histoCocktailPhiGamma->Write("CocktailPhiGamma",TObject::kOverwrite);
//             if (histoCocktailRhoGamma)histoCocktailRhoGamma->Write("CocktailRhoGamma",TObject::kOverwrite);
//             if (histoCocktailSigmaGamma)histoCocktailSigmaGamma->Write("CocktailSigmaGamma",TObject::kOverwrite);
                
        fileGammaFinal->mkdir("Theory");
        fileGammaFinal->cd("Theory");
            if (graphNLODR) graphNLODR->Write(Form("NLODoubleRatio_%s%s",centrality.Data(),collisionSystemOutput.Data()),TObject::kOverwrite);
            if (graphNLOGammaDir) graphNLOGammaDir->Write(Form("NLOGammaDir_%s%s",centrality.Data(),collisionSystemOutput.Data()),TObject::kOverwrite);
            if (graphNLOGammaPrompt) graphNLOGammaPrompt->Write(Form("NLOGammaPrompt_%s%s",centrality.Data(),collisionSystemOutput.Data()),TObject::kOverwrite);
            if (graphNLOGammaFrag) graphNLOGammaFrag->Write(Form("NLOGammaFrag_%s%s",centrality.Data(),collisionSystemOutput.Data()),TObject::kOverwrite);
                    

    
    fileGammaFinal->Write();
    fileGammaFinal->Close();
 
}

