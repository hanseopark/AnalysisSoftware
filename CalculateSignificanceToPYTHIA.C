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
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "CommonHeaders/PlottingGammaConversionHistos.h"
#include "CommonHeaders/PlottingGammaConversionAdditional.h"
#include "CommonHeaders/FittingGammaConversion.h"
#include "CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "CommonHeaders/ConversionFunctions.h"
#include "CommonHeaders/CombinationFunctions.h"

extern TRandom*    gRandom;
extern TBenchmark* gBenchmark;
extern TSystem*    gSystem;
extern TMinuit*    gMinuit;

struct SysErrorConversion {
    Double_t value;
    Double_t error;
    //    TString name;
};


//*************************************************************************************************************
// Fill chi2 histogram for null hypothesis
//*************************************************************************************************************
void FillChi2HistForNullHypoPValue    (  Int_t    n_pseudo_exp,
                                TGraphAsymmErrors*    graphTrue,
                                TH1F*    &histo,
                                TGraph*  &g_rel_stat_plus_type_a_error,
                                TGraph*  &g_rel_type_b_error,
                                Double_t &rel_type_c_error,
                                Bool_t   anti_corr_type_b,
                                Bool_t useFixedValue
                                      ) {

    // null hypothesis:
    Double_t R_true                             = 1.073;
    if(useFixedValue) cout << "expected value is " << R_true << endl;

    // random numbers
    TRandom rndm;

    for (Int_t i_pseudo_exp=0; i_pseudo_exp<n_pseudo_exp; i_pseudo_exp++) {
        Double_t sumsq                          = 0;
        Double_t eps_b                          = rndm.Gaus(0,1);
        Double_t eps_c                          = rndm.Gaus(0,1);

        // create pseudo data set
        for (Int_t ip=0; ip<g_rel_type_b_error->GetN(); ip++) {
            Double_t rel_type_b_error           = g_rel_type_b_error->GetY()[ip];
            Double_t R_mod                      = R_true;// * (1. + eps_b * rel_type_b_error) * (1. + eps_c * rel_type_c_error); // ignore error correlation
            if(graphTrue && !useFixedValue){
                R_true                          = graphTrue->Eval(g_rel_stat_plus_type_a_error->GetX()[ip]);
                R_mod                           = R_true;
            }
            Double_t rel_stat_plus_type_a_err   = g_rel_stat_plus_type_a_error->GetY()[ip];
            Double_t abs_stat_plus_type_a_err_scaled    = R_mod * rel_stat_plus_type_a_err;

            Double_t y                          = rndm.Gaus(R_mod, abs_stat_plus_type_a_err_scaled);
            Double_t nsig                       = (y - R_true)/(R_true * rel_stat_plus_type_a_err);
            sumsq                               += nsig*nsig;
        }
        histo->Fill(sumsq);
    }
}

Double_t Chi2ForNullHypoPValue(TGraphErrors* g, TGraphAsymmErrors*    graphTrue,  Bool_t useFixedValue) {

    Double_t R_true                             = 1.073;
    if(useFixedValue)
        cout << "expected value is " << R_true << endl;
    Double_t sumsq                              = 0;

    for (Int_t i=0; i<g->GetN(); i++) {
        if(graphTrue && !useFixedValue){
                R_true                          = graphTrue->Eval(g->GetX()[i]);
            }
        Double_t y                              = g->GetY()[i];
        Double_t y_sig                          = g->GetEY()[i];
        Double_t y_sig_rel                      = y_sig/y;
//         Double_t nsig                   = (y-R_true)/(R_true * y_sig_rel);
        Double_t nsig                           = (y-R_true)/(y_sig);
        sumsq                                   += nsig*nsig;
        cout << "ytrue: " << R_true << " ymeas: " << y << " y_sig: " << y_sig << " y_sig_rel: " << y_sig_rel << " nsig: " << nsig << endl;
    }
    return sumsq;
}

void drawLatexAdd(TString latextext, Double_t textcolumn, Double_t textrow, Double_t textSizePixel,Bool_t setFont = kFALSE, Bool_t setFont2 = kFALSE, Bool_t alignRight = kFALSE, Color_t textcolor = kBlack){
    TLatex *latexDummy                          = new TLatex(textcolumn ,textrow,latextext);
    SetStyleTLatex( latexDummy, textSizePixel,4);
    if(setFont)
        latexDummy->SetTextFont(62);
    if(setFont2)
        latexDummy->SetTextFont(43);
    if(alignRight)
        latexDummy->SetTextAlign(31);
    latexDummy->SetTextColor(textcolor);
    latexDummy->Draw();
}

void CalculateSignificanceToPYTHIA(TString suffix = "pdf",Int_t minBinSig = 0,
        Int_t maxBinSig = 2, Double_t systematiccorrelation = 1,Bool_t useFixedValue = kFALSE

){

    // set global variables
    TString date                                = ReturnDateString();
    gROOT->Reset();
    gROOT->SetStyle("Plain");

    StyleSettingsThesis();
    SetPlotStyle();

    TString dateForOutput                       = ReturnDateStringForOutput();
    cout << dateForOutput.Data() << endl;
    //___________________________________ Declaration of files _____________________________________________
    TString collisionSystem8TeV                 = "pp, #sqrt{#it{s}} = 8 TeV";

    TString outputDir                           = Form("%s/%s/calculatePValues",suffix.Data(),dateForOutput.Data());
    gSystem->Exec("mkdir -p "+outputDir);


    Double_t mesonMassExpectPi0                 = TDatabasePDG::Instance()->GetParticle(111)->Mass();
    Double_t mesonMassExpectEta                 = TDatabasePDG::Instance()->GetParticle(221)->Mass();

    // load input files
    TFile* inputFileSpectrum1                   = new TFile("CalculatePValues/7TeV_data_PCMResultsFullCorrection_PP.root");
    TFile* inputFileSpectrum2                   = new TFile("CalculatePValues/data_PCMResultsFullCorrection_PP_8TeV_pi0andEtaNewPileup.root");
    if(!inputFileSpectrum1 || !inputFileSpectrum1)
        return;

    TDirectory* directoryPi01                   = (TDirectory*)inputFileSpectrum1->Get("Pi07TeV");
    TDirectory* directoryPi02                   = (TDirectory*)inputFileSpectrum2->Get("Pi08TeV");
    if(!directoryPi01 || !directoryPi02)
        return;

    TGraphAsymmErrors* graphCorrectedYieldPi01          = (TGraphAsymmErrors*)directoryPi01->Get("graphInvCrossSectionPi0");
    TGraphAsymmErrors* graphCorrectedYieldPi0SysTot1    = (TGraphAsymmErrors*)directoryPi01->Get("InvCrossSectionPi0Sys");
    TGraphAsymmErrors* graphCorrectedYieldPi02          = (TGraphAsymmErrors*)directoryPi02->Get("graphInvCrossSectionPi0");
    TGraphAsymmErrors* graphCorrectedYieldPi0SysTot2    = (TGraphAsymmErrors*)directoryPi02->Get("InvCrossSectionPi0Sys");

    TGraphAsymmErrors* graphCorrectedYieldRatio12       = CalculateAsymGraphRatioToGraph(graphCorrectedYieldPi02,graphCorrectedYieldPi01);
    TGraphAsymmErrors* graphCorrectedYieldRatioSysTot12 = CalculateAsymGraphRatioToGraph(graphCorrectedYieldPi0SysTot2,graphCorrectedYieldPi0SysTot1);

    //___________________________________ Calculation of ratio w/o material __________________________________
    for (Int_t i=0; i<graphCorrectedYieldRatioSysTot12->GetN(); i++) {
        Double_t tempValue                              = systematiccorrelation*(TMath::Sqrt(TMath::Power(((graphCorrectedYieldRatioSysTot12->GetEYlow()[i])/graphCorrectedYieldRatioSysTot12->GetY()[i]),2)-2*TMath::Power(0.09,2)));
        graphCorrectedYieldRatioSysTot12->GetEYlow()[i] = tempValue*graphCorrectedYieldRatioSysTot12->GetY()[i];
        graphCorrectedYieldRatioSysTot12->GetEYhigh()[i]= tempValue*graphCorrectedYieldRatioSysTot12->GetY()[i];
    }

    TString fileNameTheory                      = "ExternalInput/Theory/TheoryCompilationPP.root";
    TFile* fileTheoryCompilation                = new TFile(fileNameTheory.Data());

        // 8TeV Pythia8 Monash2013:
        TH1F* histoPythia8InvXSectionPi08TeV    = (TH1F*) fileTheoryCompilation->Get("histoInvSecPythia8Monash2013LegoPi08TeV");
        // 7TeV Pythia8 Monash2013:
        TH1F* histoPythia8InvXSectionPi07TeV    = (TH1F*) fileTheoryCompilation->Get("histoInvSecPythia8Monash2013LegoPi07TeV");

        // plotting data 8TeV/7TeV
    TH1F *ratioOfEnergiesToPythia               = (TH1F*)histoPythia8InvXSectionPi08TeV->Clone("ratioOfEnergiesToPythia");
    ratioOfEnergiesToPythia->SetTitle("");
    ratioOfEnergiesToPythia->SetLineColor(kRed+2);
    ratioOfEnergiesToPythia->SetMarkerColor(kRed+2);
    ratioOfEnergiesToPythia->SetMarkerStyle(24);
    ratioOfEnergiesToPythia->SetLineWidth(2.);

    for (Int_t i=1; i<=ratioOfEnergiesToPythia->GetNbinsX(); i++) {
      Double_t temp                             = histoPythia8InvXSectionPi08TeV->GetBinContent(i)/histoPythia8InvXSectionPi07TeV->GetBinContent(i);
      ratioOfEnergiesToPythia->SetBinContent(i,temp);
      ratioOfEnergiesToPythia->SetBinError(i,0.);
    }
    TGraphAsymmErrors* graphratioPythia         = new TGraphAsymmErrors(ratioOfEnergiesToPythia);


    // ***************************************************************************************************************
    // ******************************* Plotting eta/pi0 ratio for single measurements ********************************
    // ***************************************************************************************************************
    Style_t markerstyles[4]                     ={34,20,24,29};
    Size_t markersize[4]                        ={2.3,2.3,2.3,3.4};
    Width_t  widthLinesBoxes                    = 1.4;
    Color_t colorData[4]                        ={kRed+2,kBlue+2,kGreen+2,kMagenta+2};
    Double_t textSizeLabelsPixel                = 54;
    TCanvas* canvas87ratiocombo                 = new TCanvas("canvas87ratiocombo","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvas87ratiocombo, 0.1, 0.01, 0.01, 0.125);
    canvas87ratiocombo->SetLogx();

    Double_t textsizeLabels87ratio              = 0;
    Double_t textsizeFac87ratio                 = 0;
    if (canvas87ratiocombo->XtoPixel(canvas87ratiocombo->GetX2()) <canvas87ratiocombo->YtoPixel(canvas87ratiocombo->GetY1()) ){
        textsizeLabels87ratio                   = (Double_t)textSizeLabelsPixel/canvas87ratiocombo->XtoPixel(canvas87ratiocombo->GetX2()) ;
        textsizeFac87ratio                      = (Double_t)1./canvas87ratiocombo->XtoPixel(canvas87ratiocombo->GetX2()) ;
    } else {
        textsizeLabels87ratio                   = (Double_t)textSizeLabelsPixel/canvas87ratiocombo->YtoPixel(canvas87ratiocombo->GetY1());
        textsizeFac87ratio                      = (Double_t)1./canvas87ratiocombo->YtoPixel(canvas87ratiocombo->GetY1());
    }
    Double_t yMin                               = 0.71;
    Double_t yMax                               = 2.9;
    Double_t xMin                               = 0.23;
    Double_t xMax                               = 19;

    TH2F * histo2D87ratiocombo;
    histo2D87ratiocombo                         = new TH2F("histo2D87ratiocombo","histo2D87ratiocombo",1000,xMin,xMax,1000,yMin,yMax    );
    SetStyleHistoTH2ForGraphs(histo2D87ratiocombo, "#it{p}_{T} (GeV/#it{c})","Ratio 8 to 7 TeV", 0.85*textsizeLabels87ratio, textsizeLabels87ratio,0.85*textsizeLabels87ratio,1.1*textsizeLabels87ratio, 0.9, 0.65, 510, 510);
    histo2D87ratiocombo->GetXaxis()->SetMoreLogLabels();
    histo2D87ratiocombo->GetXaxis()->SetLabelOffset(-0.01);
    histo2D87ratiocombo->Draw();
    Double_t markerScale                        = 1.3;

    TGraphAsymmErrors* graphCorrectedYieldRatio12_WOXErr = (TGraphAsymmErrors*) graphCorrectedYieldRatio12->Clone("graphCorrectedYieldRatio12_WOXErr");
    ProduceGraphAsymmWithoutXErrors(graphCorrectedYieldRatio12_WOXErr);

    DrawGammaSetMarkerTGraphAsym(graphCorrectedYieldRatio12_WOXErr, markerstyles[1], markerScale*markersize[1], colorData[1], colorData[1], widthLinesBoxes, kFALSE);
    graphCorrectedYieldRatio12_WOXErr->SetLineWidth(widthLinesBoxes);

    DrawGammaSetMarkerTGraphAsym(graphCorrectedYieldRatioSysTot12, markerstyles[1], markerScale*markersize[1], colorData[1] , colorData[1], widthLinesBoxes, kTRUE);
    graphCorrectedYieldRatioSysTot12->SetLineWidth(0);
    graphCorrectedYieldRatioSysTot12->Draw("2,same");

    graphCorrectedYieldRatio12_WOXErr->Draw("p,same");

    DrawGammaSetMarkerTGraphAsym(graphratioPythia, markerstyles[2], markerScale*markersize[2], colorData[2] , colorData[2]);
    graphratioPythia->Draw("p,same,e");

    TLegend* legend87ratio;
    legend87ratio                               = GetAndSetLegend2(0.25, 0.85, 0.48, 0.85+(textsizeLabels87ratio*2*0.9), textSizeLabelsPixel);
    legend87ratio->AddEntry(graphCorrectedYieldRatioSysTot12,"8 to 7 TeV ratio","pf");
    legend87ratio->AddEntry(graphratioPythia,"pythia ratio","p");
    legend87ratio->Draw();

    drawLatexAdd("ALICE this work",0.92,0.92,0.85*textsizeLabels87ratio,kFALSE,kFALSE,kTRUE);
    drawLatexAdd("pp, PCM",0.92, 0.92-(1*textsizeLabels87ratio*0.9),0.85*textsizeLabels87ratio,kFALSE,kFALSE,kTRUE);

    canvas87ratiocombo->Update();
    canvas87ratiocombo->SaveAs(Form("%s/%s.%s",outputDir.Data(),"Ratios", suffix.Data()));


        //*********************************************************************************
        // Calculating Significance of 8 to 7 TeV pi0 x-section ratio to pythia prediction
        //*********************************************************************************

        const Int_t nPointSig                   = maxBinSig- minBinSig +1;

        // filling of graphs
        TGraph* graphAbsTypeAPlusStatErr        = new TGraph(nPointSig);
        TGraph* graphRelTypeAPlusStatErr        = new TGraph(nPointSig);
        TGraph* graphRelTypeBErr                = new TGraph(nPointSig);
        TGraphErrors* graphTypeAPlusStatErr     = new TGraphErrors(nPointSig);
        for (Int_t i = 0; i < nPointSig; i++){
            Double_t pt                         = graphCorrectedYieldRatioSysTot12->GetX()[minBinSig+i];
            Double_t R                          = graphCorrectedYieldRatioSysTot12->GetY()[minBinSig+i];
            Double_t statErr                    = graphCorrectedYieldRatio12->GetEYlow()[minBinSig+i];
            Double_t sysAErr                    = graphCorrectedYieldRatioSysTot12->GetEYlow()[minBinSig+i];
            Double_t statSysAErr                = TMath::Sqrt(statErr*statErr+sysAErr*sysAErr);
            Double_t relErrStatSysA             = statSysAErr/ R;
            Double_t relErrSysB                 = graphCorrectedYieldRatioSysTot12->GetEYlow()[minBinSig+i]/ graphCorrectedYieldRatioSysTot12->GetY()[minBinSig+i];
            cout << "Using point with pT: " << Form("%1.2f", pt) << " and value: " << Form("%1.2f", R) << "+-" << Form("%1.2f stat", statErr) << "+-" << Form("%1.2f sys", sysAErr) << endl;
            graphAbsTypeAPlusStatErr->SetPoint(i, pt, statSysAErr);
            graphRelTypeAPlusStatErr->SetPoint(i, pt, relErrStatSysA);
            graphRelTypeBErr->SetPoint(i, pt, relErrSysB);
            graphTypeAPlusStatErr->SetPoint(i, pt, R);
            graphTypeAPlusStatErr->SetPointError(i,0,statSysAErr);
        }
        Double_t relErrSysC                     = graphCorrectedYieldRatioSysTot12->GetEYlow()[minBinSig]/ graphCorrectedYieldRatioSysTot12->GetY()[minBinSig];

        // Generating pseudo data
        TH1F* histoChi2NullHypo                 = new TH1F("histoChi2NullHypo", "histoChi2NullHypo", 3000, 0., 1500);
        FillChi2HistForNullHypoPValue(100000000, graphratioPythia, histoChi2NullHypo, graphRelTypeAPlusStatErr, graphRelTypeBErr, relErrSysC, kFALSE,useFixedValue);

        // calculating significance
        Double_t chi2Data                       = Chi2ForNullHypoPValue(graphTypeAPlusStatErr,graphratioPythia,useFixedValue);
        Int_t binChi2Data                       = histoChi2NullHypo->FindBin(chi2Data);
        Int_t binFirst                          = 1;
        Int_t binLast                           = histoChi2NullHypo->GetNbinsX();
        Double_t intTot                         = histoChi2NullHypo->Integral(binFirst, binLast);
        Double_t intData                        = histoChi2NullHypo->Integral(binChi2Data, binLast);
        Double_t pValue                         = intData/ intTot;
        Double_t nSigma                         = PValueToNSigma(pValue);

        cout << "pp 8TeV: chi2data = " << chi2Data << "\t, pVal = " << pValue << "\t, nSigma = " << nSigma << endl;

        // preparing histo for drawing
        TH1F* histoChi2NullData                 = (TH1F*)histoChi2NullHypo->Clone("histoChi2NullData");
        for (Int_t i = 1; i < binChi2Data; i++ ){
            histoChi2NullData->SetBinContent(i,-1);
        }
        for (Int_t i = 1; i < binLast; i++ ){
            histoChi2NullData->SetBinError(i,0);
        }

        //*********************************************************************************
        // Plotting Significance of 8 to 7 TeV ratio
        //*********************************************************************************


        TCanvas* canvasSignificance             = new TCanvas("canvasSignificance","",200,10,1400,1100);  // gives the page size
        DrawGammaCanvasSettings( canvasSignificance,  0.09, 0.08, 0.015, 0.1);
        canvasSignificance->SetLogy(1);

        Int_t textSizeLabelsPixelSignificance   = 48;
        Double_t textsizeLabelsSignificance     = 0;
        if (canvasSignificance->XtoPixel(canvasSignificance->GetX2()) < canvasSignificance->YtoPixel(canvasSignificance->GetY1())){
            textsizeLabelsSignificance          = (Double_t)textSizeLabelsPixelSignificance/canvasSignificance->XtoPixel(canvasSignificance->GetX2()) ;
        } else {
            textsizeLabelsSignificance          = (Double_t)textSizeLabelsPixelSignificance/canvasSignificance->YtoPixel(canvasSignificance->GetY1());
        }

        SetStyleHistoTH1ForGraphs(histoChi2NullHypo, "Test statistics #it{t}","Counts",
                                0.85*textsizeLabelsSignificance, textsizeLabelsSignificance,
                                0.85*textsizeLabelsSignificance, textsizeLabelsSignificance,
                                0.95, 1., 510, 510);
        histoChi2NullHypo->GetYaxis()->SetLabelOffset(0.005);

        Int_t firstAbove                        = histoChi2NullHypo->FindFirstBinAbove();
        Int_t lastAbove                         = histoChi2NullHypo->FindLastBinAbove();

        histoChi2NullHypo->GetXaxis()->SetRange(1,lastAbove+1);
        DrawGammaSetMarker(histoChi2NullHypo, 1, 0.1, kGray+2 , kGray+2);
        histoChi2NullHypo->DrawCopy("");

        Color_t colorHypo                       = kRed+2;

        DrawGammaSetMarker(histoChi2NullData, 1, 0.1, colorHypo , colorHypo);
        histoChi2NullData->SetFillColor(colorHypo);
        histoChi2NullData->SetFillStyle(3356);
        histoChi2NullData->Draw("same,lf");
        histoChi2NullHypo->DrawCopy("same");
        histoChi2NullHypo->DrawCopy("same,axis");

            TLatex *labelSignificanceEnergy     = new TLatex(0.58,0.92,"pp #sqrt{s} = 7 and 8 TeV");
            SetStyleTLatex( labelSignificanceEnergy, 0.85*textsizeLabelsSignificance,4);
            labelSignificanceEnergy->Draw();

            TLatex *labelSignificanceALICE      = new TLatex(0.58,0.87,"ALICE pseudo data");
            SetStyleTLatex( labelSignificanceALICE, 0.85*textsizeLabelsSignificance,4);
            labelSignificanceALICE->Draw();

            DrawGammaLines(chi2Data, chi2Data, 0, histoChi2NullHypo->GetMaximum()*0.1, 3, colorHypo, 7);

            TLatex *labelTData                  = new TLatex(chi2Data+10,histoChi2NullHypo->GetMaximum()*0.1,"#it{t}_{data}");
            SetStyleTLatex( labelTData, 0.85*textsizeLabelsSignificance,4,colorHypo,42,kFALSE);
            labelTData->Draw();

            TLatex *labelSignificancePValue     = new TLatex(0.48,0.57,Form("#it{p}-value = %1.4f",pValue));
            SetStyleTLatex( labelSignificancePValue, 0.85*textsizeLabelsSignificance,4);
            labelSignificancePValue->Draw();
            TLatex *labelSignificanceNSigma     = new TLatex(0.48,0.53,Form("%1.2f #sigma",nSigma));
            SetStyleTLatex( labelSignificanceNSigma, 0.85*textsizeLabelsSignificance,4);
            labelSignificanceNSigma->Draw();


        canvasSignificance->Update();
        canvasSignificance->Print(Form("%s/lowpTpointssignificanceTest.%s",outputDir.Data(),suffix.Data()));

}