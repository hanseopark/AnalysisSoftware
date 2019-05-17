/****************************************************************************************************************************
******      provided by Gamma Conversion Group, PWGGA,                                                   *****
******      Friederike Bock, friederike.bock@cern.ch                                                    *****
*****************************************************************************************************************************/

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
#include "CommonHeaders/PlottingGammaConversionHistos.h"
#include "CommonHeaders/PlottingGammaConversionAdditional.h"
#include "CommonHeaders/FittingGammaConversion.h"
#include "CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "CommonHeaders/ConversionFunctions.h"

extern TRandom* gRandom;
extern TBenchmark*  gBenchmark;
extern TSystem* gSystem;
extern TMinuit*     gMinuit;
//________________________________________________________________________________________________________________________
TH1D* ConvertYieldHisto(TH1D* input, Bool_t DivideBy2pi, Bool_t DivideByPt, Bool_t MultiplyBy2pi, Bool_t MultiplyByPt){
    if (!input) {
        cout << "Error: Histogram is NULL" << endl;
        return NULL;
    }

    Int_t nBins                 = input->GetNbinsX();
    Double_t newValue           = 0;
    Double_t newErrorValue      = 0;
    Double_t correctionValue    = 1;

    //correct by 2pi if specified
    if (DivideBy2pi) input->Scale(1/(2*TMath::Pi()));
    if (MultiplyBy2pi) input->Scale(2*TMath::Pi());

    for(Int_t i=0;i<nBins;i++){

        //correct by 1/Pt if specified
        if(DivideByPt)    correctionValue  = 1/(input->GetBinCenter(i+1));
        if(MultiplyByPt)  correctionValue  = input->GetBinCenter(i+1);

        //set the value and error of the bin
        input->SetBinContent(i+1,input->GetBinContent(i+1)*correctionValue);
        input->SetBinError(i+1,input->GetBinError(i+1)*correctionValue);
    }

    return input;
}

TH1D* CombineDiffCentsYields  (TH1D* hist1, TH1D* hist2, TString name, Bool_t isSys){
    if (hist1 && hist2){
        TH1D* histOut = (TH1D*)hist1->Clone(name.Data());
        histOut->Add(hist2);
        histOut->Scale(0.5);

        if (isSys){
            for (Int_t i = 1; i < histOut->GetNbinsX()+1; i++){
                Double_t relErrLowerCent = 0;
                if (hist1->GetBinContent(i) != 0){
                    relErrLowerCent= hist1->GetBinError(i)/hist1->GetBinContent(i)*100 ;
                }
                Double_t relErrHigherCent = 0;
                if (hist2->GetBinContent(i) != 0){
                    relErrHigherCent = hist2->GetBinError(i)/hist2->GetBinContent(i)*100 ;
                }

                if (relErrHigherCent > relErrLowerCent){
                    histOut->SetBinError(i, histOut->GetBinContent(i)*relErrHigherCent/100);
                } else {
                    histOut->SetBinError(i, histOut->GetBinContent(i)*relErrLowerCent/100);
                }
            }
        }
        return histOut;
    } else {
        return NULL;
    }
}

TGraphErrors* ConvertHistoToGraphAndRemoveZeros(TH1D* hist, TString name){
    if (hist){
        TGraphErrors* graphOut = new TGraphErrors(hist);
        graphOut->SetName(name.Data());
        while (graphOut->GetY()[0] < 0 || graphOut->GetY()[0] == 0 )
            graphOut->RemovePoint(0);
        while (graphOut->GetY()[graphOut->GetN()-1] < 0 || graphOut->GetY()[graphOut->GetN()-1] == 0)
            graphOut->RemovePoint(graphOut->GetN()-1);
        if (graphOut->GetN() > 0)
            return graphOut;
        else
            return NULL;
    } else {
        return NULL;
    }
}

void RemoveZerosFromGraph(TGraph* graph){
    if (graph){
        while (graph->GetY()[0] < 0 || graph->GetY()[0] == 0 )
            graph->RemovePoint(0);
        while (graph->GetY()[graph->GetN()-1] < 0 || graph->GetY()[graph->GetN()-1] == 0)
            graph->RemovePoint(graph->GetN()-1);
    }
}

void PrepareChargedPionDataALICE_PbPb(TString energy = "PbPb_5.02TeV"){

    TString dateForOutput                       = ReturnDateStringForOutput();

    if (energy.CompareTo("PbPb_2.76TeV") == 0){

        TString centName[15]            = { "0005", "0510", "1020", "2040", "4060", "6080", "2030", "3040", "4050", "3050",
                                            "5080", "0010", "0020", "6070", "7080"};
        TString centNameOutput[15]      = { "0-5%", "5-10%", "10-20%", "20-40%", "40-60%", "60-80%", "20-30%", "30-40%", "40-50%", "30-50%",
                                            "50-80%", "0-10%", "0-20%", "60-70%", "70-80%"};
        TString centNameReadChId[15]    = {"0_5", "5_10", "10_20", "20_40", "40_60", "60_80", "20_30", "30_40", "40_50", "30_50",
                                            "50_80", "0_10", "0_20", "60_70", "70_80" };
        Int_t tobeMerged[2][2]          = {{0,1}, {11,2}};

        //*********************************** Final pi, K, p - PbPb 2.76TeV spectra ***************************************************************
        // Documentation: JIRA_PWGLF-258
        TFile* fileChargedIdentifiedSpectraFinal    = new TFile("ExternalInputPbPb/IdentifiedCharged_2.76TeV/PbPb276.fullpT.INEL.20140329.root");
        TFile* fileChargedIdentifiedSpectraAndRAAAdd= new TFile("ExternalInputPbPb/IdentifiedCharged_2.76TeV/PionsSpectraAndRAAForDMesonAnalysis24062014.root");
        TFile* fileChargedIdentifiedRatiosFinal     = new TFile("ExternalInputPbPb/IdentifiedCharged_2.76TeV/PbPb276.fullpT.RATIOS.20140329.root");
        TH1D* histoChargedPionSpecStat[15]           = {NULL};
        TH1D* histoChargedPionSpecSyst[15]           = {NULL};
        TH1D* histoChargedKaonSpecStat[15]           = {NULL};
        TH1D* histoChargedKaonSpecSyst[15]           = {NULL};
        TH1D* histoChargedProtonSpecStat[15]         = {NULL};
        TH1D* histoChargedProtonSpecSyst[15]         = {NULL};
        TH1D* histoKaonPionRatioStat[15]             = {NULL};
        TH1D* histoKaonPionRatioSyst[15]             = {NULL};
        TH1D* histoProtonPionRatioStat[15]           = {NULL};
        TH1D* histoProtonPionRatioSyst[15]           = {NULL};
        TGraphErrors* graphChargedPionSpecStat[15]   = {NULL};
        TGraphErrors* graphChargedPionSpecSyst[15]   = {NULL};
        TGraphErrors* graphChargedKaonSpecStat[15]   = {NULL};
        TGraphErrors* graphChargedKaonSpecSyst[15]   = {NULL};
        TGraphErrors* graphChargedProtonSpecStat[15] = {NULL};
        TGraphErrors* graphChargedProtonSpecSyst[15] = {NULL};
        TGraphErrors* graphKaonPionRatioStat[15]     = {NULL};
        TGraphErrors* graphKaonPionRatioSyst[15]     = {NULL};
        TGraphErrors* graphProtonPionRatioStat[15]   = {NULL};
        TGraphErrors* graphProtonPionRatioSyst[15]   = {NULL};

        for (Int_t cent = 0; cent < 6; cent++){
            // read pi+ + pi- spectrum
            histoChargedPionSpecStat[cent]      = (TH1D*)fileChargedIdentifiedSpectraFinal->Get(Form("hstat_PbPb276_%s_pion_sum",centName[cent].Data() ));
            histoChargedPionSpecSyst[cent]      = (TH1D*)fileChargedIdentifiedSpectraFinal->Get(Form("hsys_PbPb276_%s_pion_sum",centName[cent].Data() ));
            // calculate average charged pions spectrum
            if (histoChargedPionSpecStat[cent] && histoChargedPionSpecSyst[cent] ){
                histoChargedPionSpecStat[cent]->Sumw2();
                histoChargedPionSpecSyst[cent]->Sumw2();
                histoChargedPionSpecStat[cent]->Scale(0.5);
                histoChargedPionSpecSyst[cent]->Scale(0.5);
                // set correct name
                histoChargedPionSpecStat[cent]->SetName(Form("histoChargedPionSpecStat%s",centName[cent].Data()));
                histoChargedPionSpecSyst[cent]->SetName(Form("histoChargedPionSpecSyst%s",centName[cent].Data()));
            } else {
                cout << "couldn't find charged pion inputs for " << centName[cent].Data() << endl;
            }
            // read K+ + K- spectrum
            histoChargedKaonSpecStat[cent]      = (TH1D*)fileChargedIdentifiedSpectraFinal->Get(Form("hstat_PbPb276_%s_kaon_sum",centName[cent].Data() ));
            histoChargedKaonSpecSyst[cent]      = (TH1D*)fileChargedIdentifiedSpectraFinal->Get(Form("hsys_PbPb276_%s_kaon_sum",centName[cent].Data() ));
            // calculate average charged pions spectrum
            if (histoChargedKaonSpecStat[cent] && histoChargedKaonSpecSyst[cent] ){
                histoChargedKaonSpecStat[cent]->Sumw2();
                histoChargedKaonSpecSyst[cent]->Sumw2();
                histoChargedKaonSpecStat[cent]->Scale(0.5);
                histoChargedKaonSpecSyst[cent]->Scale(0.5);
                // set correct name
                histoChargedKaonSpecStat[cent]->SetName(Form("histoChargedKaonSpecStat%s",centName[cent].Data()));
                histoChargedKaonSpecSyst[cent]->SetName(Form("histoChargedKaonSpecSyst%s",centName[cent].Data()));
            } else {
                cout << "couldn't find charged kaon inputs for " << centName[cent].Data() << endl;
            }

            // read p + anti-p spectrum
            histoChargedProtonSpecStat[cent]      = (TH1D*)fileChargedIdentifiedSpectraFinal->Get(Form("hstat_PbPb276_%s_proton_sum",centName[cent].Data() ));
            histoChargedProtonSpecSyst[cent]      = (TH1D*)fileChargedIdentifiedSpectraFinal->Get(Form("hsys_PbPb276_%s_proton_sum",centName[cent].Data() ));
            // calculate average charged pions spectrum
            if (histoChargedProtonSpecStat[cent] && histoChargedProtonSpecSyst[cent] ){
                histoChargedProtonSpecStat[cent]->Sumw2();
                histoChargedProtonSpecSyst[cent]->Sumw2();
                histoChargedProtonSpecStat[cent]->Scale(0.5);
                histoChargedProtonSpecSyst[cent]->Scale(0.5);
                // set correct name
                histoChargedProtonSpecStat[cent]->SetName(Form("histoChargedProtonSpecStat%s",centName[cent].Data()));
                histoChargedProtonSpecSyst[cent]->SetName(Form("histoChargedProtonSpecSyst%s",centName[cent].Data()));
            } else {
                cout << "couldn't find proton inputs for " << centName[cent].Data() << endl;
            }
            // Kaon over pion ratio
            histoKaonPionRatioStat[cent]  = (TH1D*)fileChargedIdentifiedRatiosFinal->Get(Form("hstat_PbPb276_%s_kaon_to_pion_sum", centName[cent].Data() ));
            histoKaonPionRatioSyst[cent]  = (TH1D*)fileChargedIdentifiedRatiosFinal->Get(Form("hsys_PbPb276_%s_kaon_to_pion_sum", centName[cent].Data() ));
            if (histoKaonPionRatioStat[cent] && histoKaonPionRatioSyst[cent] ){
                histoKaonPionRatioStat[cent]->Sumw2();
                histoKaonPionRatioSyst[cent]->Sumw2();
                // set correct name
                histoKaonPionRatioStat[cent]->SetName(Form("histoKaonPionRatioStat%s",centName[cent].Data()));
                histoKaonPionRatioSyst[cent]->SetName(Form("histoKaonPionRatioSyst%s",centName[cent].Data()));
            }

            // Proton over pion ratio
            histoProtonPionRatioStat[cent]  = (TH1D*)fileChargedIdentifiedRatiosFinal->Get(Form("hstat_PbPb276_%s_proton_to_pion_sum", centName[cent].Data() ));
            histoProtonPionRatioSyst[cent]  = (TH1D*)fileChargedIdentifiedRatiosFinal->Get(Form("hsys_PbPb276_%s_proton_to_pion_sum", centName[cent].Data() ));
            if (histoProtonPionRatioStat[cent] && histoProtonPionRatioSyst[cent] ){
                histoProtonPionRatioStat[cent]->Sumw2();
                histoProtonPionRatioSyst[cent]->Sumw2();
                // set correct name
                histoProtonPionRatioStat[cent]->SetName(Form("histoProtonPionRatioStat%s",centName[cent].Data()));
                histoProtonPionRatioSyst[cent]->SetName(Form("histoProtonPionRatioSyst%s",centName[cent].Data()));
            }
        }

        for (Int_t centb = 0; centb < 5; centb++){
            // read pi+ + pi- spectrum
            histoChargedPionSpecStat[centb+6]      = (TH1D*)fileChargedIdentifiedSpectraAndRAAAdd->Get(Form("hPionSpectrum_%s",centNameReadChId[centb+6].Data() ));
            histoChargedPionSpecSyst[centb+6]      = (TH1D*)fileChargedIdentifiedSpectraAndRAAAdd->Get(Form("hPionSpectrumSyst_%s",centNameReadChId[centb+6].Data() ));
            // calculate average charged pions spectrum
            if (histoChargedPionSpecStat[centb+6] && histoChargedPionSpecSyst[centb] ){
                histoChargedPionSpecStat[centb+6]->Sumw2();
                histoChargedPionSpecSyst[centb+6]->Sumw2();
                histoChargedPionSpecStat[centb+6]->Scale(0.5);
                histoChargedPionSpecSyst[centb+6]->Scale(0.5);
                // set correct name
                histoChargedPionSpecStat[centb+6]->SetName(Form("histoChargedPionSpecStat%s",centName[centb+6].Data()));
                histoChargedPionSpecSyst[centb+6]->SetName(Form("histoChargedPionSpecSyst%s",centName[centb+6].Data()));
            } else {
                cout << "couldn't find charged pion inputs for " << centName[centb+6].Data() << endl;
            }
        }

        for (Int_t centb = 0; centb < 2; centb++){
            if (histoChargedPionSpecStat[centb+11]==NULL && histoChargedPionSpecStat[centb+11] == NULL){
                histoChargedPionSpecStat[centb+11]      = CombineDiffCentsYields  ( histoChargedPionSpecStat[tobeMerged[centb][0]], histoChargedPionSpecStat[tobeMerged[centb][1]],
                                                                                Form("histoChargedPionSpecStat%s",centName[centb+11].Data()) , kFALSE);
                histoChargedPionSpecSyst[centb+11]      = CombineDiffCentsYields  ( histoChargedPionSpecSyst[tobeMerged[centb][0]], histoChargedPionSpecSyst[tobeMerged[centb][1]],
                                                                                Form("histoChargedPionSpecSyst%s",centName[centb+11].Data()) , kTRUE);
            }
            histoChargedKaonSpecStat[centb+11]      = CombineDiffCentsYields  ( histoChargedKaonSpecStat[tobeMerged[centb][0]], histoChargedKaonSpecStat[tobeMerged[centb][1]],
                                                                               Form("histoChargedKaonSpecStat%s",centName[centb+11].Data()) , kFALSE);
            histoChargedKaonSpecSyst[centb+11]      = CombineDiffCentsYields  ( histoChargedKaonSpecSyst[tobeMerged[centb][0]], histoChargedKaonSpecSyst[tobeMerged[centb][1]],
                                                                               Form("histoChargedKaonSpecSyst%s",centName[centb+11].Data()) , kTRUE);
            histoChargedProtonSpecStat[centb+11]    = CombineDiffCentsYields  ( histoChargedProtonSpecStat[tobeMerged[centb][0]], histoChargedProtonSpecStat[tobeMerged[centb][1]],
                                                                               Form("histoChargedProtonSpecStat%s",centName[centb+11].Data()) , kFALSE);
            histoChargedProtonSpecSyst[centb+11]    = CombineDiffCentsYields  ( histoChargedProtonSpecSyst[tobeMerged[centb][0]], histoChargedProtonSpecSyst[tobeMerged[centb][1]],
                                                                               Form("histoChargedProtonSpecSyst%s",centName[centb+11].Data()) , kTRUE);

            histoKaonPionRatioStat[centb+11]        = CombineDiffCentsYields  ( histoKaonPionRatioStat[tobeMerged[centb][0]], histoKaonPionRatioStat[tobeMerged[centb][1]],
                                                                               Form("histoKaonPionRatioStat%s",centName[centb+11].Data()) , kFALSE);
            histoKaonPionRatioSyst[centb+11]        = CombineDiffCentsYields  ( histoKaonPionRatioSyst[tobeMerged[centb][0]], histoKaonPionRatioSyst[tobeMerged[centb][1]],
                                                                               Form("histoKaonPionRatioSyst%s",centName[centb+11].Data()) , kTRUE);
            histoProtonPionRatioStat[centb+11]        = CombineDiffCentsYields  ( histoProtonPionRatioStat[tobeMerged[centb][0]], histoProtonPionRatioStat[tobeMerged[centb][1]],
                                                                                 Form("histoProtonPionRatioStat%s",centName[centb+11].Data()) , kFALSE);
            histoProtonPionRatioSyst[centb+11]        = CombineDiffCentsYields  ( histoProtonPionRatioSyst[tobeMerged[centb][0]], histoProtonPionRatioSyst[tobeMerged[centb][1]],
                                                                                 Form("histoProtonPionRatioSyst%s",centName[centb+11].Data()) , kTRUE);
        }

        //*********************************** Final K0s, Lambda - PbPb 2.76TeV spectra ***************************************************************
        TFile* fileK0sFinalPbPb = new TFile("ExternalInputPbPb/NeutralKaon/k0s_lambda_final_spectra_12112013.root");
        TH1D* histoNeutralKaonSpecStat[15]       = {NULL};
        TH1D* histoNeutralKaonSpecSyst[15]       = {NULL};
        TH1D* histoNeutralLambdaSpecStat[15]     = {NULL};
        TH1D* histoNeutralLambdaSpecSyst[15]     = {NULL};
        TGraphErrors* graphNeutralKaonSpecStat[15]       = {NULL};
        TGraphErrors* graphNeutralKaonSpecSyst[15]       = {NULL};
        TGraphErrors* graphNeutralLambdaSpecStat[15]     = {NULL};
        TGraphErrors* graphNeutralLambdaSpecSyst[15]     = {NULL};

        for (Int_t cent = 0; cent < 6; cent++){
            // read K0s spectrum
            histoNeutralKaonSpecStat[cent]      = (TH1D*)fileK0sFinalPbPb->Get(Form("statonly_cent%s_K0s",centName[cent].Data() ));
            histoNeutralKaonSpecSyst[cent]      = (TH1D*)fileK0sFinalPbPb->Get(Form("systonly_cent%s_K0s",centName[cent].Data() ));
            // calculate average charged pions spectrum
            if (histoNeutralKaonSpecStat[cent] && histoNeutralKaonSpecSyst[cent] ){
                histoNeutralKaonSpecStat[cent]->Sumw2();
                histoNeutralKaonSpecSyst[cent]->Sumw2();
                ConvertYieldHisto(histoNeutralKaonSpecStat[cent], kTRUE, kTRUE, kFALSE, kFALSE);
                ConvertYieldHisto(histoNeutralKaonSpecSyst[cent], kTRUE, kTRUE, kFALSE, kFALSE);
                // set correct name
                histoNeutralKaonSpecStat[cent]->SetName(Form("histoNeutralKaonSpecStat%s",centName[cent].Data()));
                histoNeutralKaonSpecSyst[cent]->SetName(Form("histoNeutralKaonSpecSyst%s",centName[cent].Data()));
            } else {
                cout << "couldn't find neural kaon inputs for " << centName[cent].Data() << endl;
            }
            // read Lambda spectrum
            histoNeutralLambdaSpecStat[cent]      = (TH1D*)fileK0sFinalPbPb->Get(Form("statonly_cent%s_Lambda",centName[cent].Data() ));
            histoNeutralLambdaSpecSyst[cent]      = (TH1D*)fileK0sFinalPbPb->Get(Form("systonly_cent%s_Lambda",centName[cent].Data() ));
            // calculate average charged pions spectrum
            if (histoNeutralLambdaSpecStat[cent] && histoNeutralLambdaSpecSyst[cent] ){
                histoNeutralLambdaSpecStat[cent]->Sumw2();
                histoNeutralLambdaSpecSyst[cent]->Sumw2();
                ConvertYieldHisto(histoNeutralLambdaSpecStat[cent], kTRUE, kTRUE, kFALSE, kFALSE);
                ConvertYieldHisto(histoNeutralLambdaSpecSyst[cent], kTRUE, kTRUE, kFALSE, kFALSE);
                // set correct name
                histoNeutralLambdaSpecStat[cent]->SetName(Form("histoNeutralLambdaSpecStat%s",centName[cent].Data()));
                histoNeutralLambdaSpecSyst[cent]->SetName(Form("histoNeutralLambdaSpecSyst%s",centName[cent].Data()));
            } else {
                cout << "couldn't find neural Lambda inputs for " << centName[cent].Data() << endl;
            }
        }

        for (Int_t centb = 0; centb < 2; centb++){
            histoNeutralKaonSpecStat[centb+11]      = CombineDiffCentsYields  ( histoNeutralKaonSpecStat[tobeMerged[centb][0]], histoNeutralKaonSpecStat[tobeMerged[centb][1]],
                                                                               Form("histoNeutralKaonSpecStat%s",centName[centb+11].Data()) , kFALSE);
            histoNeutralKaonSpecSyst[centb+11]      = CombineDiffCentsYields  ( histoNeutralKaonSpecSyst[tobeMerged[centb][0]], histoNeutralKaonSpecSyst[tobeMerged[centb][1]],
                                                                               Form("histoNeutralKaonSpecSyst%s",centName[centb+11].Data()) , kTRUE);
            histoNeutralLambdaSpecStat[centb+11]      = CombineDiffCentsYields  ( histoNeutralLambdaSpecStat[tobeMerged[centb][0]], histoNeutralLambdaSpecStat[tobeMerged[centb][1]],
                                                                                 Form("histoNeutralLambdaSpecStat%s",centName[centb+11].Data()) , kFALSE);
            histoNeutralLambdaSpecSyst[centb+11]      = CombineDiffCentsYields  ( histoNeutralLambdaSpecSyst[tobeMerged[centb][0]], histoNeutralLambdaSpecSyst[tobeMerged[centb][1]],
                                                                                 Form("histoNeutralLambdaSpecSyst%s",centName[centb+11].Data()) , kTRUE);
        }


        //*********************************** Final pi, K, p - PbPb 2.76TeV RAA ***************************************************************
        // Documentation: JIRA_PWGLF-258
        TFile* fileChargedPionRAAFinal2014      = new TFile("ExternalInputPbPb/IdentifiedCharged_2.76TeV/RAA_Pion_08052014.root");
        TFile* fileChargedKaonRAAFinal2014      = new TFile("ExternalInputPbPb/IdentifiedCharged_2.76TeV/RAA_Kaon_08052014.root");
        TFile* fileChargedProtonRAAFinal2014    = new TFile("ExternalInputPbPb/IdentifiedCharged_2.76TeV/RAA_Proton_08052014.root");
        TH1D* histoChargedPionRAAStat[15]        = {NULL};
        TH1D* histoChargedPionRAASyst[15]        = {NULL};
        TH1D* histoChargedKaonRAAStat[15]        = {NULL};
        TH1D* histoChargedKaonRAASyst[15]        = {NULL};
        TH1D* histoChargedProtonRAAStat[15]      = {NULL};
        TH1D* histoChargedProtonRAASyst[15]      = {NULL};
        TGraphErrors* graphChargedPionRAAStat[15]       = {NULL};
        TGraphErrors* graphChargedPionRAASyst[15]       = {NULL};
        TGraphErrors* graphChargedKaonRAAStat[15]       = {NULL};
        TGraphErrors* graphChargedKaonRAASyst[15]       = {NULL};
        TGraphErrors* graphChargedProtonRAAStat[15]     = {NULL};
        TGraphErrors* graphChargedProtonRAASyst[15]     = {NULL};

        for (Int_t cent = 0; cent < 6; cent++){
            // read charged pion RAA
            histoChargedPionRAAStat[cent]      = (TH1D*)fileChargedPionRAAFinal2014->Get(Form("RAAPion_Stat_%s",centNameReadChId[cent].Data() ));
            histoChargedPionRAASyst[cent]      = (TH1D*)fileChargedPionRAAFinal2014->Get(Form("RAAPion_Syst_%s",centNameReadChId[cent].Data() ));
            // read charged kaon RAA
            histoChargedKaonRAAStat[cent]      = (TH1D*)fileChargedKaonRAAFinal2014->Get(Form("RAAKaon_Stat_%s",centNameReadChId[cent].Data() ));
            histoChargedKaonRAASyst[cent]      = (TH1D*)fileChargedKaonRAAFinal2014->Get(Form("RAAKaon_Syst_%s",centNameReadChId[cent].Data() ));
            // read proton RAA
            histoChargedProtonRAAStat[cent]    = (TH1D*)fileChargedProtonRAAFinal2014->Get(Form("RAAProton_Stat_%s",centNameReadChId[cent].Data() ));
            histoChargedProtonRAASyst[cent]    = (TH1D*)fileChargedProtonRAAFinal2014->Get(Form("RAAProton_Syst_%s",centNameReadChId[cent].Data() ));
        }
        for (Int_t centb = 0; centb < 6; centb++){
            // read pi+ + pi- RAA
            histoChargedPionRAAStat[centb+6]   = (TH1D*)fileChargedIdentifiedSpectraAndRAAAdd->Get(Form("hPionRAA_%s",centNameReadChId[centb+6].Data() ));
            histoChargedPionRAASyst[centb+6]   = (TH1D*)fileChargedIdentifiedSpectraAndRAAAdd->Get(Form("hPionRAASyst_%s",centNameReadChId[centb+6].Data() ));
        }

        // convert histos also to graph and remove 0s
        for (Int_t cent = 0; cent < 15; cent++){
            graphChargedPionSpecStat[cent]          = ConvertHistoToGraphAndRemoveZeros(histoChargedPionSpecStat[cent],Form("graphChargedPionSpecStat%s",centName[cent].Data()) );
            graphChargedPionSpecSyst[cent]          = ConvertHistoToGraphAndRemoveZeros(histoChargedPionSpecSyst[cent],Form("graphChargedPionSpecSyst%s",centName[cent].Data()) );
            graphChargedKaonSpecStat[cent]          = ConvertHistoToGraphAndRemoveZeros(histoChargedKaonSpecStat[cent],Form("graphChargedKaonSpecStat%s",centName[cent].Data()) );
            graphChargedKaonSpecSyst[cent]          = ConvertHistoToGraphAndRemoveZeros(histoChargedKaonSpecSyst[cent],Form("graphChargedKaonSpecSyst%s",centName[cent].Data()) );
            graphChargedProtonSpecStat[cent]        = ConvertHistoToGraphAndRemoveZeros(histoChargedProtonSpecStat[cent],Form("graphChargedProtonSpecStat%s",centName[cent].Data()) );
            graphChargedProtonSpecSyst[cent]        = ConvertHistoToGraphAndRemoveZeros(histoChargedProtonSpecSyst[cent],Form("graphChargedProtonSpecSyst%s",centName[cent].Data()) );
            graphKaonPionRatioStat[cent]            = ConvertHistoToGraphAndRemoveZeros(histoKaonPionRatioStat[cent],Form("graphKaonPionRatioStat%s",centName[cent].Data()) );
            graphKaonPionRatioSyst[cent]            = ConvertHistoToGraphAndRemoveZeros(histoKaonPionRatioSyst[cent],Form("graphKaonPionRatioSyst%s",centName[cent].Data()) );
            graphProtonPionRatioStat[cent]          = ConvertHistoToGraphAndRemoveZeros(histoProtonPionRatioStat[cent],Form("graphProtonPionRatioStat%s",centName[cent].Data()) );
            graphProtonPionRatioSyst[cent]          = ConvertHistoToGraphAndRemoveZeros(histoProtonPionRatioSyst[cent],Form("graphProtonPionRatioSyst%s",centName[cent].Data()) );

            graphNeutralKaonSpecStat[cent]          = ConvertHistoToGraphAndRemoveZeros(histoNeutralKaonSpecStat[cent],Form("graphNeutralKaonSpecStat%s",centName[cent].Data()) );
            graphNeutralKaonSpecSyst[cent]          = ConvertHistoToGraphAndRemoveZeros(histoNeutralKaonSpecSyst[cent],Form("graphNeutralKaonSpecSyst%s",centName[cent].Data()) );
            graphNeutralLambdaSpecStat[cent]        = ConvertHistoToGraphAndRemoveZeros(histoNeutralLambdaSpecStat[cent],Form("graphNeutralLambdaSpecStat%s",centName[cent].Data()) );
            graphNeutralLambdaSpecSyst[cent]        = ConvertHistoToGraphAndRemoveZeros(histoNeutralLambdaSpecSyst[cent],Form("graphNeutralLambdaSpecSyst%s",centName[cent].Data()) );

            graphChargedPionRAAStat[cent]          = ConvertHistoToGraphAndRemoveZeros(histoChargedPionRAAStat[cent],Form("graphChargedPionRAAStat%s",centName[cent].Data()) );
            graphChargedPionRAASyst[cent]          = ConvertHistoToGraphAndRemoveZeros(histoChargedPionRAASyst[cent],Form("graphChargedPionRAASyst%s",centName[cent].Data()) );
            graphChargedKaonRAAStat[cent]          = ConvertHistoToGraphAndRemoveZeros(histoChargedKaonRAAStat[cent],Form("graphChargedKaonRAAStat%s",centName[cent].Data()) );
            graphChargedKaonRAASyst[cent]          = ConvertHistoToGraphAndRemoveZeros(histoChargedKaonRAASyst[cent],Form("graphChargedKaonRAASyst%s",centName[cent].Data()) );
            graphChargedProtonRAAStat[cent]        = ConvertHistoToGraphAndRemoveZeros(histoChargedProtonRAAStat[cent],Form("graphChargedProtonRAAStat%s",centName[cent].Data()) );
            graphChargedProtonRAASyst[cent]        = ConvertHistoToGraphAndRemoveZeros(histoChargedProtonRAASyst[cent],Form("graphChargedProtonRAASyst%s",centName[cent].Data()) );
        }



        //*********************************** Final ch hadron - PbPb 2.76TeV RAA & spectra **********************************************************
        // Documentation: ???
        TFile* fileChargedSpectraAndRAAFinal            = new TFile("ExternalInputPbPb/IdentifiedCharged_2.76TeV/PbPb_RAA_sigma_2760GeV_20120809.root");
        TGraphErrors* graphChargedHadronSpecStat[15]    = {NULL};
        TGraphErrors* graphChargedHadronSpecSyst[15]    = {NULL};
        TGraphErrors* graphChargedHadronRAAStat[15]     = {NULL};
        TGraphErrors* graphChargedHadronRAASyst[15]     = {NULL};

        for (Int_t cent = 0; cent < 15; cent++){
            graphChargedHadronSpecStat[cent]        = (TGraphErrors*)fileChargedSpectraAndRAAFinal->Get(Form("pt_c%s_stat",centNameReadChId[cent].Data() ));
            graphChargedHadronSpecSyst[cent]        = (TGraphErrors*)fileChargedSpectraAndRAAFinal->Get(Form("pt_c%s_syst",centNameReadChId[cent].Data() ));
            graphChargedHadronRAAStat[cent]         = (TGraphErrors*)fileChargedSpectraAndRAAFinal->Get(Form("raa_c%s_stat",centNameReadChId[cent].Data() ));
            graphChargedHadronRAASyst[cent]         = (TGraphErrors*)fileChargedSpectraAndRAAFinal->Get(Form("raa_c%s_syst",centNameReadChId[cent].Data() ));
            RemoveZerosFromGraph(graphChargedHadronSpecStat[cent]);
            RemoveZerosFromGraph(graphChargedHadronSpecSyst[cent]);
            RemoveZerosFromGraph(graphChargedHadronRAAStat[cent]);
            RemoveZerosFromGraph(graphChargedHadronRAASyst[cent]);
        }

        //--------------------- Write Files--------------------------------------------------------------

        TString outputFileName = Form("ExternalInputPbPb/IdentifiedParticleCollection_ALICE_%s.root",dateForOutput.Data());
        TFile* fileOutput = new TFile(outputFileName,"UPDATE");

        for (Int_t cent = 0; cent < 15; cent++){
            TString directoryName       = Form("%sPbPb_2.76TeV", centNameOutput[cent].Data());
            TDirectoryFile* directory   = (TDirectoryFile*)fileOutput->Get(directoryName.Data());
            if (!directory)
                fileOutput->mkdir(directoryName.Data());
            fileOutput->cd(directoryName.Data());
            if(histoChargedPionSpecStat[cent])  histoChargedPionSpecStat[cent]->Write("histoChargedPionSpecStat", TObject::kOverwrite);
            if(histoChargedPionSpecSyst[cent])  histoChargedPionSpecSyst[cent]->Write("histoChargedPionSpecSyst", TObject::kOverwrite);
            if(histoChargedKaonSpecStat[cent])  histoChargedKaonSpecStat[cent]->Write("histoChargedKaonSpecStat", TObject::kOverwrite);
            if(histoChargedKaonSpecSyst[cent])  histoChargedKaonSpecSyst[cent]->Write("histoChargedKaonSpecSyst", TObject::kOverwrite);
            if(histoChargedProtonSpecStat[cent])  histoChargedProtonSpecStat[cent]->Write("histoProtonSpecStat", TObject::kOverwrite);
            if(histoChargedProtonSpecSyst[cent])  histoChargedProtonSpecSyst[cent]->Write("histoProtonSpecSyst", TObject::kOverwrite);
            if(histoNeutralKaonSpecStat[cent])  histoNeutralKaonSpecStat[cent]->Write("histoNeutralKaonSpecStat", TObject::kOverwrite);
            if(histoNeutralKaonSpecSyst[cent])  histoNeutralKaonSpecSyst[cent]->Write("histoNeutralKaonSpecSyst", TObject::kOverwrite);
            if(histoNeutralLambdaSpecStat[cent])  histoNeutralLambdaSpecStat[cent]->Write("histoLambdaSpecStat", TObject::kOverwrite);
            if(histoNeutralLambdaSpecSyst[cent])  histoNeutralLambdaSpecSyst[cent]->Write("histoLambdaSpecSyst", TObject::kOverwrite);

            if(histoChargedPionRAAStat[cent])  histoChargedPionRAAStat[cent]->Write("histoChargedPionRAAStat", TObject::kOverwrite);
            if(histoChargedPionRAASyst[cent])  histoChargedPionRAASyst[cent]->Write("histoChargedPionRAASyst", TObject::kOverwrite);
            if(histoChargedKaonRAAStat[cent])  histoChargedKaonRAAStat[cent]->Write("histoChargedKaonRAAStat", TObject::kOverwrite);
            if(histoChargedKaonRAASyst[cent])  histoChargedKaonRAASyst[cent]->Write("histoChargedKaonRAASyst", TObject::kOverwrite);
            if(histoKaonPionRatioStat[cent])  histoKaonPionRatioStat[cent]->Write("histoKaonPionRatioStat", TObject::kOverwrite);
            if(histoKaonPionRatioSyst[cent])  histoKaonPionRatioSyst[cent]->Write("histoKaonPionRatioSyst", TObject::kOverwrite);
            if(histoProtonPionRatioStat[cent])  histoProtonPionRatioStat[cent]->Write("histoProtonPionRatioStat", TObject::kOverwrite);
            if(histoProtonPionRatioSyst[cent])  histoProtonPionRatioSyst[cent]->Write("histoProtonPionRatioSyst", TObject::kOverwrite);

            if(graphChargedPionSpecStat[cent])  graphChargedPionSpecStat[cent]->Write("graphChargedPionSpecStat", TObject::kOverwrite);
            if(graphChargedPionSpecSyst[cent])  graphChargedPionSpecSyst[cent]->Write("graphChargedPionSpecSyst", TObject::kOverwrite);
            if(graphChargedKaonSpecStat[cent])  graphChargedKaonSpecStat[cent]->Write("graphChargedKaonSpecStat", TObject::kOverwrite);
            if(graphChargedKaonSpecSyst[cent])  graphChargedKaonSpecSyst[cent]->Write("graphChargedKaonSpecSyst", TObject::kOverwrite);
            if(graphChargedProtonSpecStat[cent])  graphChargedProtonSpecStat[cent]->Write("graphProtonSpecStat", TObject::kOverwrite);
            if(graphChargedProtonSpecSyst[cent])  graphChargedProtonSpecSyst[cent]->Write("graphProtonSpecSyst", TObject::kOverwrite);
            if(graphNeutralKaonSpecStat[cent])  graphNeutralKaonSpecStat[cent]->Write("graphNeutralKaonSpecStat", TObject::kOverwrite);
            if(graphNeutralKaonSpecSyst[cent])  graphNeutralKaonSpecSyst[cent]->Write("graphNeutralKaonSpecSyst", TObject::kOverwrite);
            if(graphNeutralLambdaSpecStat[cent])  graphNeutralLambdaSpecStat[cent]->Write("graphLambdaSpecStat", TObject::kOverwrite);
            if(graphNeutralLambdaSpecSyst[cent])  graphNeutralLambdaSpecSyst[cent]->Write("graphLambdaSpecSyst", TObject::kOverwrite);

            if(graphChargedPionRAAStat[cent])  graphChargedPionRAAStat[cent]->Write("graphChargedPionRAAStat", TObject::kOverwrite);
            if(graphChargedPionRAASyst[cent])  graphChargedPionRAASyst[cent]->Write("graphChargedPionRAASyst", TObject::kOverwrite);
            if(graphChargedKaonRAAStat[cent])  graphChargedKaonRAAStat[cent]->Write("graphChargedKaonRAAStat", TObject::kOverwrite);
            if(graphChargedKaonRAASyst[cent])  graphChargedKaonRAASyst[cent]->Write("graphChargedKaonRAASyst", TObject::kOverwrite);
            if(graphChargedProtonRAAStat[cent])  graphChargedProtonRAAStat[cent]->Write("graphChargedProtonRAAStat", TObject::kOverwrite);
            if(graphChargedProtonRAASyst[cent])  graphChargedProtonRAASyst[cent]->Write("graphChargedProtonRAASyst", TObject::kOverwrite);
            if(graphKaonPionRatioStat[cent])  graphKaonPionRatioStat[cent]->Write("graphKaonPionRatioStat", TObject::kOverwrite);
            if(graphKaonPionRatioSyst[cent])  graphKaonPionRatioSyst[cent]->Write("graphKaonPionRatioSyst", TObject::kOverwrite);
            if(graphProtonPionRatioStat[cent])  graphProtonPionRatioStat[cent]->Write("graphProtonPionRatioStat", TObject::kOverwrite);
            if(graphProtonPionRatioSyst[cent])  graphProtonPionRatioSyst[cent]->Write("graphProtonPionRatioSyst", TObject::kOverwrite);

            if(graphChargedHadronSpecStat[cent])  graphChargedHadronSpecStat[cent]->Write("graphChargedHadronSpecStat", TObject::kOverwrite);
            if(graphChargedHadronSpecSyst[cent])  graphChargedHadronSpecSyst[cent]->Write("graphChargedHadronSpecSyst", TObject::kOverwrite);
            if(graphChargedHadronRAAStat[cent])  graphChargedHadronRAAStat[cent]->Write("graphChargedHadronRAAStat", TObject::kOverwrite);
            if(graphChargedHadronRAASyst[cent])  graphChargedHadronRAASyst[cent]->Write("graphChargedHadronRAASyst", TObject::kOverwrite);

        }

        fileOutput->Write();
        fileOutput->Close();
    } else if(energy.CompareTo("PbPb_5.02TeV") == 0){

        TString centName[15]            = {"0005", "0510", "1020", "2030", "3040", "4050", "5060", "6070", "7080", "8090", "0010", "2040", "4060", "6080", "0020"};
        TString centNameOutput[15]      = {"0-5%", "5-10%", "10-20%", "20-30%", "30-40%", "40-50%", "50-60%", "60-70%", "70-80%", "80-90%", "0-10%", "20-40%", "40-60%", "60-80%", "0-20%"};
        TString centNameReadChId[10]    = {"0.00to5.00", "5.00to10.00", "10.00to20.00", "20.00to30.00", "30.00to40.00", "40.00to50.00", "50.00to60.00", "60.00to70.00", "70.00to80.00", "80.00to90.00"};

        Int_t tobeMerged[5][2]          = {{0,1}, {3,4}, {5,6}, {7,8}, {10,2}};

        TFile* filePbPb5TeV             = new TFile("ExternalInputPbPb/IdentifiedCharged_5.02TeV/Spectra_PbPbLHC15o_Combined_23-11-18_Histograms.root");
        TList* listSummedPionStat       = (TList*)filePbPb5TeV->Get("Summed_Pion");
        TList* listSummedPionSyst       = (TList*)filePbPb5TeV->Get("Summed_Pion_Sys");
        TList* listSummedKaonStat       = (TList*)filePbPb5TeV->Get("Summed_Kaon");
        TList* listSummedKaonSyst       = (TList*)filePbPb5TeV->Get("Summed_Kaon_Sys");
        TList* listSummedProtonStat     = (TList*)filePbPb5TeV->Get("Summed_Proton");
        TList* listSummedProtonSyst     = (TList*)filePbPb5TeV->Get("Summed_Proton_Sys");
        TList* listSummedKaonPionStat   = (TList*)filePbPb5TeV->Get("Summed_Kaon_Over_Summed_Pion");
        TList* listSummedKaonPionSyst   = (TList*)filePbPb5TeV->Get("Summed_Kaon_Over_Summed_Pion_Sys");
        TList* listSummedProtonPionStat = (TList*)filePbPb5TeV->Get("Summed_Proton_Over_Summed_Pion");
        TList* listSummedProtonPionSyst = (TList*)filePbPb5TeV->Get("Summed_Proton_Over_Summed_Pion_Sys");

        TH1D* histoChargedPionSpecStat[15]      = {NULL};
        TH1D* histoChargedPionSpecSyst[15]      = {NULL};
        TH1D* histoChargedKaonSpecStat[15]      = {NULL};
        TH1D* histoChargedKaonSpecSyst[15]      = {NULL};
        TH1D* histoChargedProtonSpecStat[15]    = {NULL};
        TH1D* histoChargedProtonSpecSyst[15]    = {NULL};
        TH1D* histoKaonPionRatioStat[15]        = {NULL};
        TH1D* histoKaonPionRatioSyst[15]        = {NULL};
        TH1D* histoProtonPionRatioStat[15]      = {NULL};
        TH1D* histoProtonPionRatioSyst[15]      = {NULL};
        TGraphErrors* graphChargedPionSpecStat[15]      = {NULL};
        TGraphErrors* graphChargedPionSpecSyst[15]      = {NULL};
        TGraphErrors* graphChargedKaonSpecStat[15]      = {NULL};
        TGraphErrors* graphChargedKaonSpecSyst[15]      = {NULL};
        TGraphErrors* graphChargedProtonSpecStat[15]    = {NULL};
        TGraphErrors* graphChargedProtonSpecSyst[15]    = {NULL};
        TGraphErrors* graphKaonPionRatioStat[15]        = {NULL};
        TGraphErrors* graphKaonPionRatioSyst[15]        = {NULL};
        TGraphErrors* graphProtonPionRatioStat[15]      = {NULL};
        TGraphErrors* graphProtonPionRatioSyst[15]      = {NULL};

        for (Int_t cent = 0; cent < 10; cent++){
            // read pi+ + pi- spectrum
            histoChargedPionSpecStat[cent]      = (TH1D*)listSummedPionStat->FindObject(Form("hSpectraSummedPion_PbPb_Combined_%s",centNameReadChId[cent].Data() ));
            histoChargedPionSpecSyst[cent]      = (TH1D*)listSummedPionSyst->FindObject(Form("hSpectraSummedPion_PbPb_Combined_%s",centNameReadChId[cent].Data() ));
            // calculate average charged pions spectrum
            if (histoChargedPionSpecStat[cent] && histoChargedPionSpecSyst[cent] ){
                histoChargedPionSpecStat[cent]->Sumw2();
                histoChargedPionSpecSyst[cent]->Sumw2();
                histoChargedPionSpecStat[cent]->Scale(0.5);
                histoChargedPionSpecSyst[cent]->Scale(0.5);
                ConvertYieldHisto(histoChargedPionSpecStat[cent], kTRUE, kTRUE, kFALSE, kFALSE);
                ConvertYieldHisto(histoChargedPionSpecSyst[cent], kTRUE, kTRUE, kFALSE, kFALSE);
                // set correct name
                histoChargedPionSpecStat[cent]->SetName(Form("histoChargedPionSpecStat%s",centName[cent].Data()));
                histoChargedPionSpecSyst[cent]->SetName(Form("histoChargedPionSpecSyst%s",centName[cent].Data()));
            } else {
                cout << "couldn't find charged pion inputs for " << centName[cent].Data() << endl;
            }

            // read K+ + K- spectrum
            histoChargedKaonSpecStat[cent]      = (TH1D*)listSummedKaonStat->FindObject(Form("hSpectraSummedKaon_PbPb_Combined_%s",centNameReadChId[cent].Data() ));
            histoChargedKaonSpecSyst[cent]      = (TH1D*)listSummedKaonSyst->FindObject(Form("hSpectraSummedKaon_PbPb_Combined_%s",centNameReadChId[cent].Data() ));
            // calculate average charged kaons spectrum
            if (histoChargedKaonSpecStat[cent] && histoChargedKaonSpecSyst[cent] ){
                histoChargedKaonSpecStat[cent]->Sumw2();
                histoChargedKaonSpecSyst[cent]->Sumw2();
                histoChargedKaonSpecStat[cent]->Scale(0.5);
                histoChargedKaonSpecSyst[cent]->Scale(0.5);
                ConvertYieldHisto(histoChargedKaonSpecStat[cent], kTRUE, kTRUE, kFALSE, kFALSE);
                ConvertYieldHisto(histoChargedKaonSpecSyst[cent], kTRUE, kTRUE, kFALSE, kFALSE);
                // set correct name
                histoChargedKaonSpecStat[cent]->SetName(Form("histoChargedKaonSpecStat%s",centName[cent].Data()));
                histoChargedKaonSpecSyst[cent]->SetName(Form("histoChargedKaonSpecSyst%s",centName[cent].Data()));
            } else {
                cout << "couldn't find charged kaon inputs for " << centName[cent].Data() << endl;
            }

            // read p + anti-p spectrum
            histoChargedProtonSpecStat[cent]      = (TH1D*)listSummedProtonStat->FindObject(Form("hSpectraSummedProton_PbPb_Combined_%s",centNameReadChId[cent].Data() ));
            histoChargedProtonSpecSyst[cent]      = (TH1D*)listSummedProtonSyst->FindObject(Form("hSpectraSummedProton_PbPb_Combined_%s",centNameReadChId[cent].Data() ));
            // calculate average charged proton spectrum
            if (histoChargedProtonSpecStat[cent] && histoChargedProtonSpecSyst[cent] ){
                histoChargedProtonSpecStat[cent]->Sumw2();
                histoChargedProtonSpecSyst[cent]->Sumw2();
                histoChargedProtonSpecStat[cent]->Scale(0.5);
                histoChargedProtonSpecSyst[cent]->Scale(0.5);
                ConvertYieldHisto(histoChargedProtonSpecStat[cent], kTRUE, kTRUE, kFALSE, kFALSE);
                ConvertYieldHisto(histoChargedProtonSpecSyst[cent], kTRUE, kTRUE, kFALSE, kFALSE);
                // set correct name
                histoChargedProtonSpecStat[cent]->SetName(Form("histoChargedProtonSpecStat%s",centName[cent].Data()));
                histoChargedProtonSpecSyst[cent]->SetName(Form("histoChargedProtonSpecSyst%s",centName[cent].Data()));
            } else {
                cout << "couldn't find charged proton inputs for " << centName[cent].Data() << endl;
            }
            // Kaon over pion ratio
            histoKaonPionRatioStat[cent]  = (TH1D*)listSummedKaonPionStat->FindObject(Form("hSpectraRatioSummedKaon_over_SummedPion_PbPb_Combined_%s", centNameReadChId[cent].Data() ));
            histoKaonPionRatioSyst[cent]  = (TH1D*)listSummedKaonPionSyst->FindObject(Form("hSpectraRatioSummedKaon_over_SummedPion_PbPb_Combined_%s", centNameReadChId[cent].Data() ));
            if (histoKaonPionRatioStat[cent] && histoKaonPionRatioSyst[cent] ){
                histoKaonPionRatioStat[cent]->Sumw2();
                histoKaonPionRatioSyst[cent]->Sumw2();
                // set correct name
                histoKaonPionRatioStat[cent]->SetName(Form("histoKaonPionRatioStat%s",centName[cent].Data()));
                histoKaonPionRatioSyst[cent]->SetName(Form("histoKaonPionRatioSyst%s",centName[cent].Data()));
            }

            // Proton over pion ratio
            histoProtonPionRatioStat[cent]  = (TH1D*)listSummedProtonPionStat->FindObject(Form("hSpectraRatioSummedProton_over_SummedPion_PbPb_Combined_%s", centNameReadChId[cent].Data() ));
            histoProtonPionRatioSyst[cent]  = (TH1D*)listSummedProtonPionSyst->FindObject(Form("hSpectraRatioSummedProton_over_SummedPion_PbPb_Combined_%s", centNameReadChId[cent].Data() ));
            if (histoProtonPionRatioStat[cent] && histoProtonPionRatioSyst[cent] ){
                histoProtonPionRatioStat[cent]->Sumw2();
                histoProtonPionRatioSyst[cent]->Sumw2();
                // set correct name
                histoProtonPionRatioStat[cent]->SetName(Form("histoProtonPionRatioStat%s",centName[cent].Data()));
                histoProtonPionRatioSyst[cent]->SetName(Form("histoProtonPionRatioSyst%s",centName[cent].Data()));
            }
        }

        for (Int_t centb = 0; centb < 5; centb++){
            histoChargedPionSpecStat[centb+10]      = CombineDiffCentsYields  ( histoChargedPionSpecStat[tobeMerged[centb][0]], histoChargedPionSpecStat[tobeMerged[centb][1]],
                                                                                Form("histoChargedPionSpecStat%s",centName[centb+10].Data()) , kFALSE);
            histoChargedPionSpecSyst[centb+10]      = CombineDiffCentsYields  ( histoChargedPionSpecSyst[tobeMerged[centb][0]], histoChargedPionSpecSyst[tobeMerged[centb][1]],
                                                                                Form("histoChargedPionSpecSyst%s",centName[centb+10].Data()) , kTRUE);
            histoChargedKaonSpecStat[centb+10]      = CombineDiffCentsYields  ( histoChargedKaonSpecStat[tobeMerged[centb][0]], histoChargedKaonSpecStat[tobeMerged[centb][1]],
                                                                                Form("histoChargedKaonSpecStat%s",centName[centb+10].Data()) , kFALSE);
            histoChargedKaonSpecSyst[centb+10]      = CombineDiffCentsYields  ( histoChargedKaonSpecSyst[tobeMerged[centb][0]], histoChargedKaonSpecSyst[tobeMerged[centb][1]],
                                                                                Form("histoChargedKaonSpecSyst%s",centName[centb+10].Data()) , kTRUE);
            histoChargedProtonSpecStat[centb+10]    = CombineDiffCentsYields  ( histoChargedProtonSpecStat[tobeMerged[centb][0]], histoChargedProtonSpecStat[tobeMerged[centb][1]],
                                                                                Form("histoChargedProtonSpecStat%s",centName[centb+10].Data()) , kFALSE);
            histoChargedProtonSpecSyst[centb+10]    = CombineDiffCentsYields  ( histoChargedProtonSpecSyst[tobeMerged[centb][0]], histoChargedProtonSpecSyst[tobeMerged[centb][1]],
                                                                                Form("histoChargedProtonSpecSyst%s",centName[centb+10].Data()) , kTRUE);

            histoKaonPionRatioStat[centb+10]        = CombineDiffCentsYields  ( histoKaonPionRatioStat[tobeMerged[centb][0]], histoKaonPionRatioStat[tobeMerged[centb][1]],
                                                                                Form("histoKaonPionRatioStat%s",centName[centb+10].Data()) , kFALSE);
            histoKaonPionRatioSyst[centb+10]        = CombineDiffCentsYields  ( histoKaonPionRatioSyst[tobeMerged[centb][0]], histoKaonPionRatioSyst[tobeMerged[centb][1]],
                                                                                Form("histoKaonPionRatioSyst%s",centName[centb+10].Data()) , kTRUE);
            histoProtonPionRatioStat[centb+10]      = CombineDiffCentsYields  ( histoProtonPionRatioStat[tobeMerged[centb][0]], histoProtonPionRatioStat[tobeMerged[centb][1]],
                                                                                Form("histoProtonPionRatioStat%s",centName[centb+10].Data()) , kFALSE);
            histoProtonPionRatioSyst[centb+10]      = CombineDiffCentsYields  ( histoProtonPionRatioSyst[tobeMerged[centb][0]], histoProtonPionRatioSyst[tobeMerged[centb][1]],
                                                                                Form("histoProtonPionRatioSyst%s",centName[centb+10].Data()) , kTRUE);
        }

        // convert histos also to graph and remove 0s
        for (Int_t cent = 0; cent < 15; cent++){
            graphChargedPionSpecStat[cent]          = ConvertHistoToGraphAndRemoveZeros(histoChargedPionSpecStat[cent],Form("graphChargedPionSpecStat%s",centName[cent].Data()) );
            graphChargedPionSpecSyst[cent]          = ConvertHistoToGraphAndRemoveZeros(histoChargedPionSpecSyst[cent],Form("graphChargedPionSpecSyst%s",centName[cent].Data()) );
            graphChargedKaonSpecStat[cent]          = ConvertHistoToGraphAndRemoveZeros(histoChargedKaonSpecStat[cent],Form("graphChargedKaonSpecStat%s",centName[cent].Data()) );
            graphChargedKaonSpecSyst[cent]          = ConvertHistoToGraphAndRemoveZeros(histoChargedKaonSpecSyst[cent],Form("graphChargedKaonSpecSyst%s",centName[cent].Data()) );
            graphChargedProtonSpecStat[cent]        = ConvertHistoToGraphAndRemoveZeros(histoChargedProtonSpecStat[cent],Form("graphChargedProtonSpecStat%s",centName[cent].Data()) );
            graphChargedProtonSpecSyst[cent]        = ConvertHistoToGraphAndRemoveZeros(histoChargedProtonSpecSyst[cent],Form("graphChargedProtonSpecSyst%s",centName[cent].Data()) );
            graphKaonPionRatioStat[cent]            = ConvertHistoToGraphAndRemoveZeros(histoKaonPionRatioStat[cent],Form("graphKaonPionRatioStat%s",centName[cent].Data()) );
            graphKaonPionRatioSyst[cent]            = ConvertHistoToGraphAndRemoveZeros(histoKaonPionRatioSyst[cent],Form("graphKaonPionRatioSyst%s",centName[cent].Data()) );
            graphProtonPionRatioStat[cent]          = ConvertHistoToGraphAndRemoveZeros(histoProtonPionRatioStat[cent],Form("graphProtonPionRatioStat%s",centName[cent].Data()) );
            graphProtonPionRatioSyst[cent]          = ConvertHistoToGraphAndRemoveZeros(histoProtonPionRatioSyst[cent],Form("graphProtonPionRatioSyst%s",centName[cent].Data()) );
        }

        TString outputFileName = Form("ExternalInputPbPb/IdentifiedParticleCollection_ALICE_%s.root",dateForOutput.Data());
        TFile* fileOutput = new TFile(outputFileName,"UPDATE");

        for (Int_t cent = 0; cent < 15; cent++){
            TString directoryName       = Form("%sPbPb_5.02TeV", centNameOutput[cent].Data());
            TDirectoryFile* directory   = (TDirectoryFile*)fileOutput->Get(directoryName.Data());
            if (!directory)
                fileOutput->mkdir(directoryName.Data());
            fileOutput->cd(directoryName.Data());
            if(histoChargedPionSpecStat[cent])  histoChargedPionSpecStat[cent]->Write("histoChargedPionSpecStat", TObject::kOverwrite);
            if(histoChargedPionSpecSyst[cent])  histoChargedPionSpecSyst[cent]->Write("histoChargedPionSpecSyst", TObject::kOverwrite);
            if(histoChargedKaonSpecStat[cent])  histoChargedKaonSpecStat[cent]->Write("histoChargedKaonSpecStat", TObject::kOverwrite);
            if(histoChargedKaonSpecSyst[cent])  histoChargedKaonSpecSyst[cent]->Write("histoChargedKaonSpecSyst", TObject::kOverwrite);
            if(histoChargedProtonSpecStat[cent])  histoChargedProtonSpecStat[cent]->Write("histoProtonSpecStat", TObject::kOverwrite);
            if(histoChargedProtonSpecSyst[cent])  histoChargedProtonSpecSyst[cent]->Write("histoProtonSpecSyst", TObject::kOverwrite);
            if(histoKaonPionRatioStat[cent])  histoKaonPionRatioStat[cent]->Write("histoKaonPionRatioStat", TObject::kOverwrite);
            if(histoKaonPionRatioSyst[cent])  histoKaonPionRatioSyst[cent]->Write("histoKaonPionRatioSyst", TObject::kOverwrite);
            if(histoProtonPionRatioStat[cent])  histoProtonPionRatioStat[cent]->Write("histoProtonPionRatioStat", TObject::kOverwrite);
            if(histoProtonPionRatioSyst[cent])  histoProtonPionRatioSyst[cent]->Write("histoProtonPionRatioSyst", TObject::kOverwrite);

            if(graphChargedPionSpecStat[cent])  graphChargedPionSpecStat[cent]->Write("graphChargedPionSpecStat", TObject::kOverwrite);
            if(graphChargedPionSpecSyst[cent])  graphChargedPionSpecSyst[cent]->Write("graphChargedPionSpecSyst", TObject::kOverwrite);
            if(graphChargedKaonSpecStat[cent])  graphChargedKaonSpecStat[cent]->Write("graphChargedKaonSpecStat", TObject::kOverwrite);
            if(graphChargedKaonSpecSyst[cent])  graphChargedKaonSpecSyst[cent]->Write("graphChargedKaonSpecSyst", TObject::kOverwrite);
            if(graphChargedProtonSpecStat[cent])  graphChargedProtonSpecStat[cent]->Write("graphProtonSpecStat", TObject::kOverwrite);
            if(graphChargedProtonSpecSyst[cent])  graphChargedProtonSpecSyst[cent]->Write("graphProtonSpecSyst", TObject::kOverwrite);
            if(graphKaonPionRatioStat[cent])  graphKaonPionRatioStat[cent]->Write("graphKaonPionRatioStat", TObject::kOverwrite);
            if(graphKaonPionRatioSyst[cent])  graphKaonPionRatioSyst[cent]->Write("graphKaonPionRatioSyst", TObject::kOverwrite);
            if(graphProtonPionRatioStat[cent])  graphProtonPionRatioStat[cent]->Write("graphProtonPionRatioStat", TObject::kOverwrite);
            if(graphProtonPionRatioSyst[cent])  graphProtonPionRatioSyst[cent]->Write("graphProtonPionRatioSyst", TObject::kOverwrite);

        }

        fileOutput->Write();
        fileOutput->Close();

    }else{
        cout << "ERROR: no external files know for this system" << endl;
    }
}

