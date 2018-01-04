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
#include "CommonHeaders/PlottingGammaConversionHistos.h"
#include "CommonHeaders/PlottingGammaConversionAdditional.h"
#include "CommonHeaders/ExtractSignalBinning.h"
#include "CommonHeaders/FittingGammaConversion.h"
#include "CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "CommonHeaders/ConversionFunctions.h"
#include "CommonHeaders/CombinationFunctions.h"

void PlotM02PtDependentFunctions (){

    StyleSettingsThesis();
    SetPlotStyle();

    gStyle->SetLineStyleString(11,"40 20 4 20 4 20");
    gStyle->SetLineStyleString(12,"70 20 4 20");
    gStyle->SetLineStyleString(13,"30 20 4 20 4 20");
    gStyle->SetLineStyleString(14,"15 20 4 20");
    gStyle->SetLineStyleString(15,"15 20 4 10");
    gStyle->SetLineStyleString(16,"10 20 2 20 10 20 2 20");

    TF1* M02CutsPtDep[16];
    Bool_t enableM02[16]	= {1,		1,		1,		1,		1,		1,		1,		1,		1,		1,		1,		1,		1,		1,		1,		1};
    Double_t param0[16]     = {0.5,     0.4,    0.7,    0.5,    0.5,    0.7,    0.39,   0.5,    0.7,    0.5,    0.5,    0.7, 	0.5,	0.7, 	0.7,	0.7};
    Double_t param1[16]     = {0.32,    0.27,   0.36,   0.31,   0.30,   0.35,   0.25,   0.25,   0.37,   0.32,   0.32,   0.32,	0.27,  	0.27,	0.32, 	0.34};
    Double_t param2[16]     = {0.0072,  0.0072, 0.0072, 0.0072, 0.0072, 0.0072, 0.0072, 0.0072, 0.0072, 0.0152, 0.0238, 0.0238,	0.0092,	0.0092,	0.0072, 0.0072};
    Color_t colors[16]      = {kBlack, kRed+1, kBlue+1, kOrange, kMagenta+1, kGreen+2, kViolet+1, kYellow+2, kCyan+2, kTeal+2, kPink+2, kRed-6, kGray+1, kGreen-6, kBlue-6, kMagenta-6};
    Style_t styles[16]      = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
    Double_t maxPtPart1[16];
	Int_t nParamsTot		= 0;
    for (Int_t i = 0; i < 16; i++){
        M02CutsPtDep[i]     =  new TF1(Form("M02Cut%i",i),"(x<=TMath::Sqrt(([0]-[1])/[2]))*([1]+[2]*TMath::Power(x,2)) + (x>TMath::Sqrt(([0]-[1])/[2])) * [0]");
        M02CutsPtDep[i]->SetParameter(0, param0[i]);
        M02CutsPtDep[i]->SetParameter(1, param1[i]);
        M02CutsPtDep[i]->SetParameter(2, param2[i]);
        M02CutsPtDep[i]->SetLineColor(colors[i]);
        M02CutsPtDep[i]->SetLineStyle(styles[i]);
        M02CutsPtDep[i]->SetRange(0,10);
        maxPtPart1[i]       = TMath::Sqrt((param0[i]-param1[i])/param2[i]);
		if (enableM02[i])nParamsTot;
    }

    TCanvas* canvasM02 = new TCanvas("canvasM02","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasM02, 0.08, 0.01, 0.01, 0.09);
//     canvasM02->SetLogx();

        TH2F * histo2DM02;
        histo2DM02 = new TH2F("histo2DM02","histo2DM02",11000,0, 10,1000,0,1.41);
        SetStyleHistoTH2ForGraphs(histo2DM02, "#it{p}_{T} (GeV/#it{c})","#sigma_{long}^{2}",0.035,0.04, 0.035,0.04, 1.,1.);
    //     histo2DM02->GetXaxis()->SetMoreLogLabels();
    //     histo2DM02->GetXaxis()->SetNoExponent();
        histo2DM02->Draw("copy");

        TLegend* legendM02   = GetAndSetLegend2(0.11, 0.97-(0.04*(nParamsTot)/2), 0.95, 0.97, 30, 2, "", 43, 0.08);
        for (Int_t i = 0; i < 16; i++){
			if (enableM02[i]){
	            M02CutsPtDep[i]->Draw("same");
    	        legendM02->AddEntry(M02CutsPtDep[i],Form("#it{E} #leq %1.1f: %1.2f + %1.4f#it{E}^{2}, #it{E} > %1.1f: %1.2f", maxPtPart1[i],param1[i],param2[i],maxPtPart1[i],param0[i] ),"l");
			}
        }
        legendM02->Draw();

    canvasM02->SaveAs("M02PtDependentCuts.eps");
    canvasM02->SaveAs("M02PtDependentCuts.pdf");
}

