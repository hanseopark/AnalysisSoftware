/****************************************************************************************************************************
******         provided by Gamma Conversion Group, PWG4,                                                     *****
******        Daniel Mühlheim, d.muehlheim@cern.ch                                                        *****
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
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "CommonHeaders/PlottingGammaConversionHistos.h"
#include "CommonHeaders/PlottingGammaConversionAdditional.h"
#include "CommonHeaders/FittingGammaConversion.h"
#include "CommonHeaders/ExtractSignalBinning.h"
// #include "CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "CommonHeaders/ConversionFunctions.h"
//#include "CommonHeaders/CombinationFunctions.h"

//____________________________________________________________________________________________________________________________________________
void CombineMesonMeasurementsYieldPt(TString suffix="eps"){

    TString date = ReturnDateString(kTRUE);

    TString ALICEperfor = "ALICE performance";
    
    gROOT->Reset();
    gROOT->SetStyle("Plain");
    
    StyleSettingsThesis();
    SetPlotStyle();
    
    gStyle->SetEndErrorSize(0);
    
    TString dateForOutput                       = ReturnDateStringForOutput();
    cout << dateForOutput.Data() << endl;

    TString fileNameIntYield8TeV                = "ExternalInput/CocktailInputPP_IntegYieldComp8TeV.root";
    TString fileNameIntYield7TeV                = "ExternalInput/CocktailInputPP_IntegYieldComp7TeV.root";
    TString fileNameIntYield2760GeV             = "ExternalInput/CocktailInputPP_IntegYieldComp2_76TeV.root";

    TString outputDir                           = Form("%s/%s/CombineMesonMeasurementsYieldPt",suffix.Data(),dateForOutput.Data());
    gSystem->Exec("mkdir -p "+outputDir);

      // **********************************************************************************************************************
      // ************************* intYield+mean pT histos ********************************************************************
      // **********************************************************************************************************************

      TFile* fileIntYield8TeV                              = new TFile(fileNameIntYield8TeV.Data());
      if(fileIntYield8TeV->IsZombie()) return;
      TList* directoryIntYield8TeV                    = (TList*)fileIntYield8TeV->Get("pp_8TeV");

      TFile* fileIntYield2760GeV                           = new TFile(fileNameIntYield2760GeV.Data());
      if(fileIntYield2760GeV->IsZombie()) return;
      TList* directoryIntYield2760GeV                 = (TList*)fileIntYield2760GeV->Get("pp_2.76TeV");

      TGraphAsymmErrors* yieldStat8TeV = (TGraphAsymmErrors*) directoryIntYield8TeV->FindObject("IntegratedYieldStatVsMassMeson");
      TGraphAsymmErrors* yieldSys8TeV = (TGraphAsymmErrors*) directoryIntYield8TeV->FindObject("IntegratedYieldSysVsMassMeson");
      TGraphAsymmErrors* yieldSysFunc8TeV = (TGraphAsymmErrors*) directoryIntYield8TeV->FindObject("IntegratedYieldSysFuncVsMassMeson");

      TGraphAsymmErrors* yieldStat2760GeVtemp = (TGraphAsymmErrors*) directoryIntYield2760GeV->FindObject("IntegratedYieldStatVsMassMeson");
      TGraphAsymmErrors* yieldSys2760GeVtemp = (TGraphAsymmErrors*) directoryIntYield2760GeV->FindObject("IntegratedYieldSysVsMassMeson");
      TGraphAsymmErrors* yieldSysFunc2760GeVtemp = (TGraphAsymmErrors*) directoryIntYield2760GeV->FindObject("IntegratedYieldSysFuncVsMassMeson");

      TGraphAsymmErrors* yieldStat2760GeV = (TGraphAsymmErrors*) yieldStat2760GeVtemp->Clone("stat2760GeV");
      TGraphAsymmErrors* yieldSys2760GeV = (TGraphAsymmErrors*) yieldSys2760GeVtemp->Clone("stat2760GeV");
      TGraphAsymmErrors* yieldSysFunc2760GeV = (TGraphAsymmErrors*) yieldSysFunc2760GeVtemp->Clone("stat2760GeV");

      yieldStat2760GeV->RemovePoint(1);
      yieldSys2760GeV->RemovePoint(1);
      yieldSysFunc2760GeV->RemovePoint(1);
      yieldStat2760GeV->RemovePoint(1);
      yieldSys2760GeV->RemovePoint(1);
      yieldSysFunc2760GeV->RemovePoint(1);
      yieldStat2760GeV->RemovePoint(1);
      yieldSys2760GeV->RemovePoint(1);
      yieldSysFunc2760GeV->RemovePoint(1);

      TGraphAsymmErrors* yieldStat2760GeV_charged = (TGraphAsymmErrors*) yieldStat2760GeVtemp->Clone("stat2760GeV_charged");
      TGraphAsymmErrors* yieldSys2760GeV_charged = (TGraphAsymmErrors*) yieldSys2760GeVtemp->Clone("sys2760GeV_charged");
      TGraphAsymmErrors* yieldSysFunc2760GeV_charged = (TGraphAsymmErrors*) yieldSysFunc2760GeVtemp->Clone("sysFunc2760GeV_charged");

      yieldStat2760GeV_charged->RemovePoint(0);
      yieldSys2760GeV_charged->RemovePoint(0);
      yieldSysFunc2760GeV_charged->RemovePoint(0);
      yieldStat2760GeV_charged->RemovePoint(1);
      yieldSys2760GeV_charged->RemovePoint(1);
      yieldSysFunc2760GeV_charged->RemovePoint(1);


      TCanvas* canvasY                         = GetAndSetCanvas("canvasY", 0, 0, 1700, 1800);
      DrawCanvasSettings(canvasY, 0.07, 0.015, 0.015, 0.07);
      canvasY->SetLogy();

      TH1D* dummyHisto                          = new TH1D("dummyHisto", "", 100000, 0., 0.69);
      SetHistogramm(dummyHisto,"#it{M} (GeV/#it{c}^{2})","d#it{N}/d#it{y}", 0.05, 6, 0.73, 0.7);
      dummyHisto->GetXaxis()->SetLabelOffset(-0.001);
      dummyHisto->GetXaxis()->SetRangeUser(0,2);
      dummyHisto->GetXaxis()->SetTickLength(0.02);
      dummyHisto->Draw();

      DrawGammaSetMarkerTGraphAsym(yieldSys2760GeV, GetDefaultMarkerStyle("2.76TeV","",""), 2, kMagenta-6, kMagenta-6);
      gStyle->SetEndErrorSize(10);
      yieldSys2760GeV->DrawClone("e||");
      DrawGammaSetMarkerTGraphAsym(yieldSysFunc2760GeV, GetDefaultMarkerStyle("2.76TeV","",""), 2, kMagenta+4, kMagenta+4);
      gStyle->SetEndErrorSize(15);
      yieldSysFunc2760GeV->DrawClone("e[]");
      DrawGammaSetMarkerTGraphAsym(yieldStat2760GeV, GetDefaultMarkerStyle("2.76TeV","",""), 4, kMagenta+2, kMagenta+2);
      yieldStat2760GeV->SetLineWidth(2);
      yieldStat2760GeV->DrawClone("pZ");

      DrawGammaSetMarkerTGraphAsym(yieldSys8TeV, GetDefaultMarkerStyle("8TeV","",""), 2, kGreen-6, kGreen-6);
      gStyle->SetEndErrorSize(10);
      yieldSys8TeV->DrawClone("e||");
      DrawGammaSetMarkerTGraphAsym(yieldSysFunc8TeV, GetDefaultMarkerStyle("8TeV","",""), 2, kGreen+4, kGreen+4);
      gStyle->SetEndErrorSize(15);
      yieldSysFunc8TeV->DrawClone("e[]");
      DrawGammaSetMarkerTGraphAsym(yieldStat8TeV, GetDefaultMarkerStyle("8TeV","",""), 4, kGreen+2, kGreen+2);
      yieldStat8TeV->SetLineWidth(2);
      yieldStat8TeV->DrawClone("pZ");

//      TLatex *labelALICE                     = new TLatex(0.94, 0.925, "ALICE");
//      SetStyleTLatex( labelALICE, 0.04,4, 1, 42, kTRUE, 31);
//      labelALICE->Draw();

//      TLatex *labelpp                     = new TLatex(0.94, 0.885, "pp");
//      SetStyleTLatex( labelpp, 0.04,4, 1, 42, kTRUE, 31);
//      labelpp->Draw();

      TLegend* legend                         = NULL;
      TLegend* legend2                        = NULL;

      Int_t nColumns                          = 2;
      Double_t widthLegend                    = 0.35;
      Double_t widthLegend2                   = 0.2;
      Double_t yLegend                        = 0.95;
      Double_t xLegend                        = 0.18;
      Double_t xLegend2                       = 0.52;
      Double_t margin                         = 0.27;

      legend                              = GetAndSetLegend2(xLegend, 0.10, xLegend+widthLegend*nColumns, 0.10+(1.05*0.042*3),75,nColumns, "" ,43,margin);
      legend2                             = GetAndSetLegend2(xLegend2, 0.10, xLegend2+widthLegend2, 0.10+(1.05*0.042*3),75,1, "" ,43,0);

      TGraphAsymmErrors* yieldStat2760GeV_c = (TGraphAsymmErrors*) yieldStat2760GeV->Clone();
      TGraphAsymmErrors* yieldStat8TeV_c = (TGraphAsymmErrors*) yieldStat8TeV->Clone();
      yieldStat2760GeV_c->SetMarkerSize(5);
      yieldStat8TeV_c->SetMarkerSize(5);

      legend->AddEntry(yieldStat2760GeV_c, "   ", "pe");
      legend->AddEntry(yieldStat8TeV_c, "        ", "pe");
      legend2->AddEntry((TObject*)0, "stat. unc.", "");

      yieldSys2760GeV->SetTitle("");
      legend->AddEntry(yieldSys2760GeV, "", "l");
      yieldSys8TeV->SetTitle("");
      legend->AddEntry(yieldSys8TeV, "", "l");

      legend2->AddEntry((TObject*)0, "sys. unc.", "");
      legend->AddEntry(yieldSysFunc2760GeV, " ", "l");
      legend->AddEntry(yieldSysFunc8TeV, " ", "l");
      legend2->AddEntry((TObject*)0, "sys. unc. func. form", "");

      legend->Draw();
      legend2->Draw();


      TLatex *labelpp                     = new TLatex(0.31, 0.3, "ALICE, pp");
      SetStyleTLatex( labelpp, 0.04,4, 1, 42, kTRUE, 31);
      labelpp->Draw();

      TLatex *labelpp2                     = new TLatex(0.36, 0.25, "#sqrt{#it{s}} = 2.76TeV");
      SetStyleTLatex( labelpp2, 0.04,4, 1, 42, kTRUE, 31);
      labelpp2->Draw();

      TLatex *labelpp3                     = new TLatex(0.55, 0.25, "#sqrt{#it{s}} = 8TeV");
      SetStyleTLatex( labelpp3, 0.04,4, 1, 42, kTRUE, 31);
      labelpp3->Draw();

      yieldSys8TeV->DrawClone("e||");

      canvasY->SaveAs(Form("%s/ParticleYield.%s", outputDir.Data(), suffix.Data()));

      //--------------
      //----mean pT
      //--------------


      TGraphAsymmErrors* meanPtStat8TeV = (TGraphAsymmErrors*) directoryIntYield8TeV->FindObject("meanPtStatVsMassMeson");
      TGraphAsymmErrors* meanPtSys8TeV = (TGraphAsymmErrors*) directoryIntYield8TeV->FindObject("meanPtSysVsMassMeson");
      TGraphAsymmErrors* meanPtSysFunc8TeV = (TGraphAsymmErrors*) directoryIntYield8TeV->FindObject("meanPtSysFuncVsMassMeson");

      TGraphAsymmErrors* meanPtStat2760GeVtemp = (TGraphAsymmErrors*) directoryIntYield2760GeV->FindObject("meanPtStatVsMassMeson");
      TGraphAsymmErrors* meanPtSys2760GeVtemp = (TGraphAsymmErrors*) directoryIntYield2760GeV->FindObject("meanPtSysVsMassMeson");
      TGraphAsymmErrors* meanPtSysFunc2760GeVtemp = (TGraphAsymmErrors*) directoryIntYield2760GeV->FindObject("meanPtSysFuncVsMassMeson");

      TGraphAsymmErrors* meanPtStat2760GeV = (TGraphAsymmErrors*) meanPtStat2760GeVtemp->Clone("meanPtStat2760GeV");
      TGraphAsymmErrors* meanPtSys2760GeV = (TGraphAsymmErrors*) meanPtSys2760GeVtemp->Clone("meanPtSys2760GeV");
      TGraphAsymmErrors* meanPtSysFunc2760GeV = (TGraphAsymmErrors*) meanPtSysFunc2760GeVtemp->Clone("meanPtSysFunc2760GeV");

      meanPtStat2760GeV->RemovePoint(1);
      meanPtSys2760GeV->RemovePoint(1);
      meanPtSysFunc2760GeV->RemovePoint(1);
      meanPtStat2760GeV->RemovePoint(1);
      meanPtSys2760GeV->RemovePoint(1);
      meanPtSysFunc2760GeV->RemovePoint(1);
      meanPtStat2760GeV->RemovePoint(1);
      meanPtSys2760GeV->RemovePoint(1);
      meanPtSysFunc2760GeV->RemovePoint(1);


      dummyHisto->GetYaxis()->SetTitle("#LT #it{p}_{T} #GT (GeV/#it{c})");
      dummyHisto->GetYaxis()->SetRangeUser(0.095,1.2);
      //dummyHisto->GetYaxis()->SetNoExponent();
      //dummyHisto->GetYaxis()->SetMoreLogLabels();
      dummyHisto->Draw();

      DrawGammaSetMarkerTGraphAsym(meanPtSys2760GeV, GetDefaultMarkerStyle("2.76TeV","",""), 2, kMagenta-6, kMagenta-6);
      gStyle->SetEndErrorSize(10);
      meanPtSys2760GeV->DrawClone("e||");
      DrawGammaSetMarkerTGraphAsym(meanPtSysFunc2760GeV, GetDefaultMarkerStyle("2.76TeV","",""), 2, kMagenta+4, kMagenta+4);
      gStyle->SetEndErrorSize(15);
      meanPtSysFunc2760GeV->DrawClone("e[]");
      DrawGammaSetMarkerTGraphAsym(meanPtStat2760GeV, GetDefaultMarkerStyle("2.76TeV","",""), 4, kMagenta+2, kMagenta+2);
      meanPtStat2760GeV->SetLineWidth(2);
      meanPtStat2760GeV->DrawClone("pZ");

      DrawGammaSetMarkerTGraphAsym(meanPtSys8TeV, GetDefaultMarkerStyle("8TeV","",""), 2, kGreen-6, kGreen-6);
      gStyle->SetEndErrorSize(10);
      meanPtSys8TeV->DrawClone("e||");
      DrawGammaSetMarkerTGraphAsym(meanPtSysFunc8TeV, GetDefaultMarkerStyle("8TeV","",""), 2, kGreen+4, kGreen+4);
      gStyle->SetEndErrorSize(15);
      meanPtSysFunc8TeV->DrawClone("e[]");
      DrawGammaSetMarkerTGraphAsym(meanPtStat8TeV, GetDefaultMarkerStyle("8TeV","",""), 4, kGreen+2, kGreen+2);
      meanPtStat8TeV->SetLineWidth(2);
      meanPtStat8TeV->DrawClone("pZ");

//      TLatex *labelALICE                     = new TLatex(0.94, 0.925, "ALICE");
//      SetStyleTLatex( labelALICE, 0.04,4, 1, 42, kTRUE, 31);
//      labelALICE->Draw();

//      TLatex *labelpp                     = new TLatex(0.94, 0.885, "pp");
//      SetStyleTLatex( labelpp, 0.04,4, 1, 42, kTRUE, 31);
//      labelpp->Draw();

//      TLegend* legend                         = NULL;
//      TLegend* legend2                        = NULL;

//      Int_t nColumns                          = 2;
//      Double_t widthLegend                    = 0.35;
//      Double_t yLegend                        = 0.95;
//      Double_t xLegend                        = 0.18;
//      Double_t xLegend2                       = 0.47;
//      Double_t margin                         = 0.27;

//      legend                              = GetAndSetLegend2(xLegend, 0.10, xLegend+widthLegend*nColumns, 0.10+(1.05*0.042*3),75,nColumns, "" ,43,margin);
//      legend2                             = GetAndSetLegend2(xLegend2, 0.10, xLegend2+widthLegend, 0.10+(1.05*0.042*3),75,1, "" ,43,0);

//      TGraphAsymmErrors* meanPtStat2760GeV_c = (TGraphAsymmErrors*) meanPtStat2760GeV->Clone();
//      TGraphAsymmErrors* meanPtStat8TeV_c = (TGraphAsymmErrors*) meanPtStat8TeV->Clone();
//      meanPtStat2760GeV_c->SetMarkerSize(5);
//      meanPtStat8TeV_c->SetMarkerSize(5);

//      legend->AddEntry(meanPtStat2760GeV_c, "   ", "pe");
//      legend->AddEntry(meanPtStat8TeV_c, "        ", "pe");
//      legend2->AddEntry((TObject*)0, "stat. err.", "");

//      meanPtSys2760GeV->SetTitle("");
//      legend->AddEntry(meanPtSys2760GeV, "", "l");
//      meanPtSys8TeV->SetTitle("");
//      legend->AddEntry(meanPtSys8TeV, "", "l");

//      legend2->AddEntry((TObject*)0, "sys. err.", "");
//      legend->AddEntry(meanPtSysFunc2760GeV, " ", "l");
//      legend->AddEntry(meanPtSysFunc8TeV, " ", "l");
//      legend2->AddEntry((TObject*)0, "sys. err. func. form", "");

      legend->Draw();
      legend2->Draw();


      labelpp->Draw();
      labelpp2->Draw();
      labelpp3->Draw();

      meanPtSys8TeV->DrawClone("e||");

      canvasY->SaveAs(Form("%s/ParticleMeanPt.%s", outputDir.Data(), suffix.Data()));


      TFile* fileIntYield7TeV                           = new TFile(fileNameIntYield7TeV.Data());
      if(fileIntYield7TeV->IsZombie()) return;
      TList* directoryIntYield7TeV                 = (TList*)fileIntYield7TeV->Get("pp_7TeV");

      TGraphAsymmErrors* yieldStat7TeVtemp = (TGraphAsymmErrors*) directoryIntYield7TeV->FindObject("IntegratedYieldStatVsMassMeson");
      TGraphAsymmErrors* yieldSys7TeVtemp = (TGraphAsymmErrors*) directoryIntYield7TeV->FindObject("IntegratedYieldSysVsMassMeson");
      TGraphAsymmErrors* yieldSysFunc7TeVtemp = (TGraphAsymmErrors*) directoryIntYield7TeV->FindObject("IntegratedYieldSysFuncVsMassMeson");

      TGraphAsymmErrors* yieldStat7TeV = (TGraphAsymmErrors*) yieldStat7TeVtemp->Clone("stat7TeV");
      TGraphAsymmErrors* yieldSys7TeV = (TGraphAsymmErrors*) yieldSys7TeVtemp->Clone("stat7TeV");
      TGraphAsymmErrors* yieldSysFunc7TeV = (TGraphAsymmErrors*) yieldSysFunc7TeVtemp->Clone("stat7TeV");

      yieldStat7TeV->RemovePoint(1);
      yieldSys7TeV->RemovePoint(1);
      yieldSysFunc7TeV->RemovePoint(1);
      yieldStat7TeV->RemovePoint(1);
      yieldSys7TeV->RemovePoint(1);
      yieldSysFunc7TeV->RemovePoint(1);
      yieldStat7TeV->RemovePoint(1);
      yieldSys7TeV->RemovePoint(1);
      yieldSysFunc7TeV->RemovePoint(1);

      TGraphAsymmErrors* yieldStat7TeV_charged = (TGraphAsymmErrors*) yieldStat7TeVtemp->Clone("stat7TeV_charged");
      TGraphAsymmErrors* yieldSys7TeV_charged = (TGraphAsymmErrors*) yieldSys7TeVtemp->Clone("sys7TeV_charged");
      TGraphAsymmErrors* yieldSysFunc7TeV_charged = (TGraphAsymmErrors*) yieldSysFunc7TeVtemp->Clone("sysFunc7TeV_charged");

      yieldStat7TeV_charged->RemovePoint(0);
      yieldSys7TeV_charged->RemovePoint(0);
      yieldSysFunc7TeV_charged->RemovePoint(0);
      yieldStat7TeV_charged->RemovePoint(1);
      yieldSys7TeV_charged->RemovePoint(1);
      yieldSysFunc7TeV_charged->RemovePoint(1);

      dummyHisto->GetYaxis()->SetTitle("d#it{N}/d#it{y}");
      dummyHisto->GetXaxis()->SetRangeUser(0.1,0.65);
      dummyHisto->GetXaxis()->SetNoExponent();
      dummyHisto->GetXaxis()->SetMoreLogLabels();
      dummyHisto->GetYaxis()->SetRangeUser(0.05,6);
      dummyHisto->Draw();


      yieldSys2760GeV->DrawClone("e||");
      yieldSysFunc2760GeV->DrawClone("e[]");
      yieldStat2760GeV->DrawClone("pZ");

      yieldSys8TeV->DrawClone("e||");
      yieldSysFunc8TeV->DrawClone("e[]");
      yieldStat8TeV->DrawClone("pZ");

      DrawGammaSetMarkerTGraphAsym(yieldSys7TeV, GetDefaultMarkerStyle("7TeV","",""), 2, kBlue-6, kBlue-6);
      gStyle->SetEndErrorSize(10);
      yieldSys7TeV->DrawClone("e||");
      DrawGammaSetMarkerTGraphAsym(yieldSysFunc7TeV, GetDefaultMarkerStyle("7TeV","",""), 2, kBlue+4, kBlue+4);
      gStyle->SetEndErrorSize(15);
      yieldSysFunc7TeV->DrawClone("e[]");
      DrawGammaSetMarkerTGraphAsym(yieldStat7TeV, GetDefaultMarkerStyle("7TeV","",""), 4, kBlue+2, kBlue+2);
      yieldStat7TeV->SetLineWidth(2);
      yieldStat7TeV->DrawClone("pZ");


      DrawGammaSetMarkerTGraphAsym(yieldSys2760GeV_charged, 30, 2, kMagenta-6, kMagenta-6);
      gStyle->SetEndErrorSize(10);
      yieldSys2760GeV_charged->DrawClone("e||");
      DrawGammaSetMarkerTGraphAsym(yieldSysFunc2760GeV_charged, 30, 2, kMagenta+4, kMagenta+4);
      gStyle->SetEndErrorSize(15);
      yieldSysFunc2760GeV_charged->DrawClone("e[]");
      DrawGammaSetMarkerTGraphAsym(yieldStat2760GeV_charged, 30, 4, kMagenta+2, kMagenta+2);
      yieldStat2760GeV_charged->SetLineWidth(2);
      yieldStat2760GeV_charged->DrawClone("pZ");

      DrawGammaSetMarkerTGraphAsym(yieldSys7TeV_charged, 24, 2, kBlue-6, kBlue-6);
      gStyle->SetEndErrorSize(10);
      yieldSys7TeV_charged->DrawClone("e||");
      DrawGammaSetMarkerTGraphAsym(yieldSysFunc7TeV_charged, 24, 2, kBlue+4, kBlue+4);
      gStyle->SetEndErrorSize(15);
      yieldSysFunc7TeV_charged->DrawClone("e[]");
      DrawGammaSetMarkerTGraphAsym(yieldStat7TeV_charged, 24, 4, kBlue+2, kBlue+2);
      yieldStat7TeV_charged->SetLineWidth(2);
      yieldStat7TeV_charged->DrawClone("pZ");

      nColumns = 3;
      legend                              = GetAndSetLegend2(0.18, 0.10, 0.48, 0.10+(1.05*0.042*3),75,nColumns, "" ,43,margin);
      legend2                             = GetAndSetLegend2(0.45, 0.10, 0.45+widthLegend2, 0.10+(1.05*0.042*3),75,1, "" ,43,0);

      TGraphAsymmErrors* yieldStat2760GeV_c2 = (TGraphAsymmErrors*) yieldStat2760GeV->Clone();
      TGraphAsymmErrors* yieldStat7TeV_c2 = (TGraphAsymmErrors*) yieldStat7TeV->Clone();
      TGraphAsymmErrors* yieldStat8TeV_c2 = (TGraphAsymmErrors*) yieldStat8TeV->Clone();
      yieldStat2760GeV_c2->SetMarkerSize(5);
      yieldStat7TeV_c2->SetMarkerSize(5);
      yieldStat8TeV_c2->SetMarkerSize(5);

      legend->AddEntry(yieldStat2760GeV_c2, "   ", "pe");
      legend->AddEntry(yieldStat7TeV_c2, "   ", "pe");
      legend->AddEntry(yieldStat8TeV_c2, "   ", "pe");
      legend2->AddEntry((TObject*)0, "stat. unc.", "");

      yieldSys2760GeV->SetTitle("");
      legend->AddEntry(yieldSys2760GeV, "", "l");
      yieldSys7TeV->SetTitle("");
      legend->AddEntry(yieldSys7TeV, "", "l");
      yieldSys8TeV->SetTitle("");
      legend->AddEntry(yieldSys8TeV, "", "l");

      legend2->AddEntry((TObject*)0, "sys. unc.", "");
      legend->AddEntry(yieldSysFunc2760GeV, " ", "l");
      legend->AddEntry(yieldSysFunc7TeV, " ", "l");
      legend->AddEntry(yieldSysFunc8TeV, " ", "l");
      legend2->AddEntry((TObject*)0, "sys. unc. func. form", "");

      legend->Draw();
      legend2->Draw();

      TLatex *label_                     = new TLatex(0.9, 0.9, "charged particles");
      SetStyleTLatex( label_, 0.04,4, 1, 42, kTRUE, 31);
      label_->Draw();
      TLatex *label_2                     = new TLatex(0.9, 0.85, "with open markers");
      SetStyleTLatex( label_2, 0.04,4, 1, 42, kTRUE, 31);
      label_2->Draw();

      TLatex *labelpp_                     = new TLatex(0.36, 0.3, "ALICE, pp@#sqrt{#it{s}}");
      SetStyleTLatex( labelpp_, 0.04,4, 1, 42, kTRUE, 31);
      labelpp_->Draw();

      TLatex *labelpp2_                     = new TLatex(0.25, 0.25, "2.76TeV");
      SetStyleTLatex( labelpp2_, 0.04,4, 1, 42, kTRUE, 31);
      labelpp2_->Draw();

      TLatex *labelpp2_1                     = new TLatex(0.35, 0.25, "7TeV");
      SetStyleTLatex( labelpp2_1, 0.04,4, 1, 42, kTRUE, 31);
      labelpp2_1->Draw();

      TLatex *labelpp3_                     = new TLatex(0.45, 0.25, "8TeV");
      SetStyleTLatex( labelpp3_, 0.04,4, 1, 42, kTRUE, 31);
      labelpp3_->Draw();
      labelpp_->Draw();

      yieldSys8TeV->DrawClone("e||");

      canvasY->SetLogx();
      canvasY->SaveAs(Form("%s/ParticleYield_with7TeV_pass4.%s", outputDir.Data(), suffix.Data()));

      //--------------
      //----mean pT
      //--------------

      TGraphAsymmErrors* meanPtStat7TeVtemp = (TGraphAsymmErrors*) directoryIntYield7TeV->FindObject("meanPtStatVsMassMeson");
      TGraphAsymmErrors* meanPtSys7TeVtemp = (TGraphAsymmErrors*) directoryIntYield7TeV->FindObject("meanPtSysVsMassMeson");
      TGraphAsymmErrors* meanPtSysFunc7TeVtemp = (TGraphAsymmErrors*) directoryIntYield7TeV->FindObject("meanPtSysFuncVsMassMeson");

      TGraphAsymmErrors* meanPtStat7TeV = (TGraphAsymmErrors*) meanPtStat7TeVtemp->Clone("meanPtStat");
      TGraphAsymmErrors* meanPtSys7TeV = (TGraphAsymmErrors*) meanPtSys7TeVtemp->Clone("meanPtSys");
      TGraphAsymmErrors* meanPtSysFunc7TeV = (TGraphAsymmErrors*) meanPtSysFunc7TeVtemp->Clone("meanPtSysFunc");

      meanPtStat7TeV->RemovePoint(1);
      meanPtSys7TeV->RemovePoint(1);
      meanPtSysFunc7TeV->RemovePoint(1);
      meanPtStat7TeV->RemovePoint(1);
      meanPtSys7TeV->RemovePoint(1);
      meanPtSysFunc7TeV->RemovePoint(1);
      meanPtStat7TeV->RemovePoint(1);
      meanPtSys7TeV->RemovePoint(1);
      meanPtSysFunc7TeV->RemovePoint(1);

      TGraphAsymmErrors* meanPtStat7TeV_charged = (TGraphAsymmErrors*) meanPtStat7TeVtemp->Clone("meanPtstat7TeV_charged");
      TGraphAsymmErrors* meanPtSys7TeV_charged = (TGraphAsymmErrors*) meanPtSys7TeVtemp->Clone("meanPtsys7TeV_charged");
      TGraphAsymmErrors* meanPtSysFunc7TeV_charged = (TGraphAsymmErrors*) meanPtSysFunc7TeVtemp->Clone("meanPtsysFunc7TeV_charged");

      meanPtStat7TeV_charged->RemovePoint(0);
      meanPtSys7TeV_charged->RemovePoint(0);
      meanPtSysFunc7TeV_charged->RemovePoint(0);
      meanPtStat7TeV_charged->RemovePoint(1);
      meanPtSys7TeV_charged->RemovePoint(1);
      meanPtSysFunc7TeV_charged->RemovePoint(1);


      TGraphAsymmErrors* meanPtStat2760GeV_charged = (TGraphAsymmErrors*) meanPtStat2760GeVtemp->Clone("meanPtstat2760GeV_charged");
      TGraphAsymmErrors* meanPtSys2760GeV_charged = (TGraphAsymmErrors*) meanPtSys2760GeVtemp->Clone("meanPtsys2760GeV_charged");
      TGraphAsymmErrors* meanPtSysFunc2760GeV_charged = (TGraphAsymmErrors*) meanPtSysFunc2760GeVtemp->Clone("meanPtsysFunc2760GeV_charged");

      meanPtStat2760GeV_charged->RemovePoint(0);
      meanPtSys2760GeV_charged->RemovePoint(0);
      meanPtSysFunc2760GeV_charged->RemovePoint(0);
      meanPtStat2760GeV_charged->RemovePoint(1);
      meanPtSys2760GeV_charged->RemovePoint(1);
      meanPtSysFunc2760GeV_charged->RemovePoint(1);

      dummyHisto->GetYaxis()->SetTitle("#LT #it{p}_{T} #GT (GeV/#it{c})");
      dummyHisto->GetYaxis()->SetRangeUser(0.095,1.2);
      //dummyHisto->GetYaxis()->SetNoExponent();
      //dummyHisto->GetYaxis()->SetMoreLogLabels();
      dummyHisto->Draw();

      meanPtSys2760GeV->DrawClone("e||");
      meanPtSysFunc2760GeV->DrawClone("e[]");
      meanPtStat2760GeV->DrawClone("pZ");

      meanPtSys8TeV->DrawClone("e||");
      meanPtSysFunc8TeV->DrawClone("e[]");
      meanPtStat8TeV->DrawClone("pZ");

      DrawGammaSetMarkerTGraphAsym(meanPtSys7TeV, GetDefaultMarkerStyle("7TeV","",""), 2, kBlue-6, kBlue-6);
      gStyle->SetEndErrorSize(10);
      meanPtSys7TeV->DrawClone("e||");
      DrawGammaSetMarkerTGraphAsym(meanPtSysFunc7TeV, GetDefaultMarkerStyle("7TeV","",""), 2, kBlue+4, kBlue+4);
      gStyle->SetEndErrorSize(15);
      meanPtSysFunc7TeV->DrawClone("e[]");
      DrawGammaSetMarkerTGraphAsym(meanPtStat7TeV, GetDefaultMarkerStyle("7TeV","",""), 4, kBlue+2, kBlue+2);
      meanPtStat7TeV->SetLineWidth(2);
      meanPtStat7TeV->DrawClone("pZ");

      DrawGammaSetMarkerTGraphAsym(meanPtSys7TeV_charged, 24, 2, kBlue-6, kBlue-6);
      gStyle->SetEndErrorSize(10);
      meanPtSys7TeV_charged->DrawClone("e||");
      DrawGammaSetMarkerTGraphAsym(meanPtSysFunc7TeV_charged, 24, 2, kBlue+4, kBlue+4);
      gStyle->SetEndErrorSize(15);
      meanPtSysFunc7TeV_charged->DrawClone("e[]");
      DrawGammaSetMarkerTGraphAsym(meanPtStat7TeV_charged, 24, 4, kBlue+2, kBlue+2);
      meanPtStat7TeV_charged->SetLineWidth(2);
      meanPtStat7TeV_charged->DrawClone("pZ");

      DrawGammaSetMarkerTGraphAsym(meanPtSys2760GeV_charged, 30, 2, kMagenta-6, kMagenta-6);
      gStyle->SetEndErrorSize(10);
      meanPtSys2760GeV_charged->DrawClone("e||");
      DrawGammaSetMarkerTGraphAsym(meanPtSysFunc2760GeV_charged, 30, 2, kMagenta+4, kMagenta+4);
      gStyle->SetEndErrorSize(15);
      meanPtSysFunc2760GeV_charged->DrawClone("e[]");
      DrawGammaSetMarkerTGraphAsym(meanPtStat2760GeV_charged, 30, 4, kMagenta+2, kMagenta+2);
      meanPtStat2760GeV_charged->SetLineWidth(2);
      meanPtStat2760GeV_charged->DrawClone("pZ");

      legend->Draw();
      legend2->Draw();

      labelpp_->Draw();
      labelpp2_->Draw();
      labelpp2_1->Draw();
      labelpp3_->Draw();
      labelpp_->Draw();


      meanPtSys8TeV->DrawClone("e||");

      canvasY->SaveAs(Form("%s/ParticleMeanPt_with7TeV_pass4.%s", outputDir.Data(), suffix.Data()));
      return;
}
