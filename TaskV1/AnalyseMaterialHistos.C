// provided by Gamma Conversion Group, $ALICE_ROOT/PWG4/GammaConv ;https://twiki.cern.ch/twiki/bin/view/ALICE/PWG4GammaConversion

//This file is not supposed to be run on outputfiles of the GammaConv-Software before the 30th Sept 2010.

#include <stdlib.h>
#include <iostream>
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
#include "TProfile2D.h"
#include "TF1.h"
#include "TVirtualFitter.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TDatabasePDG.h"
#include "TMinuit.h"
#include "../CommonHeaders/PlottingMeson.h"
#include "TLatex.h"
#include "TASImage.h"
#include "TMath.h"
#include "TPostScript.h"
#include "TGraphErrors.h"
#include "TArrow.h"
#include "TGraphAsymmErrors.h"
#include "TGaxis.h"
#include "TMarker.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "../CommonHeaders/PlottingGammaConversionHistos.h"
#include "../CommonHeaders/PlottingGammaConversionAdditional.h"
#include "../CommonHeaders/FittingGammaConversion.h"
#include "../CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "../CommonHeaders/ConversionFunctions.h"
#include "AnalyseMaterialHistos.h"


TH1F* ScaleByIntegralWithinLimits(const TH1F* hist, Double_t rmin, Double_t rmax, Double_t mcGasCorrectionFactor )
{
    Double_t rminTPC = rmin; //cm
    Double_t rmaxTPC = rmax;
    TH1F * histClone = (TH1F*)hist->Clone();
    TH1F * histSca = (TH1F*)hist->Clone();
    histClone->SetAxisRange(rmin,rmax);
    Double_t nconvInGas =  histClone->Integral();
    GammaScalingHistogramm(histSca,1./(mcGasCorrectionFactor*nconvInGas));
    cout << "nconvInGas::"<< " " <<nconvInGas<< "  " << TMath::Sqrt(nconvInGas)<< " "<<  TMath::Sqrt(nconvInGas)/nconvInGas<<endl;

    return histSca;
}




void AnalyseMaterialHistos(
			   TString fileName               = "", 
			   TString cutSelection           = "",
			   TString Suffix                 = "", 
			   TString fileNameMC               = "", 
			   //			   Int_t mode                     = 0 , 
			   TString optionPeriod           = "", TString optionGrid="Plots", Double_t fMinPt=0.){
    
    gROOT->Reset();
    gROOT->SetStyle("Plain");
    StyleSettingsThesis();
    SetPlotStyle();

 
    TObjArray *arr;
    arr = cutSelection.Tokenize("_");
    TObjString* objstrEvent;
    TObjString* objstrGamma;
    TString eventCutNumber; 
    TString  gammaCutNumber; 

    Double_t rMinGas=95.;
    Double_t rMaxGas=145.;

    Int_t rebinRPlots=2;
    //    Int_t rebinRPlots=10.;
    Int_t rebinZPlots=4;
    Int_t rebinPtPlots=4;

    TLine * 	lineRLimits[13];
    TLine * 	lineRLimitsPurity[13];
    // =	 		new TLine (0,1,2*TMath::Pi(),1);
 
    const int 	nBinsR = 				12;  
    //   Float_t         arrayRBins[14] =                {0.,1.5,3.5,5.7,8.6,13.,21.,33.5,41.,55.,72.,90.,150.,180};

    Double_t         arrayZBins[13] =                {10.,15.,20.,40.,40.,60.,60.,80.,120.,120.,200.,200.,200};
    // {0.,1.5,5.,8.6,13.,21.,33.5,41.,55.,72.,90.,150.,180};
    Double_t         arrayRBins[13] =                {0.,1.5,5.,8.6,13.,21.,33.5,41.,55.,72.,95.,145.,180};
    TString         arrayNamesRangesRBins[12]=      {"0 cm < R < 1.5 cm",       //0
					       "1.5 cm < R < 5. cm",    //1
					       "5. cm < R < 8.6 cm",    //2
					       "8.6 cm < R < 13 cm",    //3
					       "13 cm < R < 21 cm",     //4
					       "21 cm < R < 33.5 cm",   //5
					       "33.5 cm < R < 41 cm",   //6   
					       "41 cm < R < 55 cm",     //7
					       "55 cm < R < 72 cm",     //8
					       "72 cm < R < 90 cm",     //9
					       "90 cm < R < 150 cm",    //10
					       "150 cm < R < 180 cm"};  //11



    TString         arrayNamesRBins[12]=    {"Vertex",          //0
					     "BeamPipe+SPD 1",           //1
					     "SPD 2",           //2
					     "Thermal shield/Support between SPD/SDD",       //3
					     "SDD 1 +Thermal shield",                         //4
					     "SDD 2 +Thermal shield",                        //5 
					     "SSD 1",                                        //6
					     "SSD 2",                                        //7
					     "Air + TPC in. cont. vessel + CO_{2}",          //8
					     "CO_{2} + TPC in. field cage vessel+TPC rods ", //9
					     "Ne: CO_{2}: N_{2}",                            //10
					     "Ne: CO_{2}: N_{2}"};                           //11

    // const int 	nBinsR = 				10;  
    // Float_t         arrayRBins[11] =                {0.,2.,3.5,13.,21.,33.5,55.,72.,90.,150.,180};
    // TString         arrayNamesRangesRBins[11]=      {"0 cm < R < 2 cm",       //0
    // 						     "2. cm < R < 3.5 cm",    //1
    // 						     "3.5 cm < R < 13 cm",    //2
    // 						     "13 cm < R < 21. cm",   //3
    // 						     "21 cm < R < 33.5 cm",   //4
    // 						     "33.5 cm < R < 55 cm",   //5   
    // 						     "55 cm < R < 72 cm",     //6
    // 						     "72 cm < R < 90 cm",     //7
    // 						     "90 cm < R < 150 cm",    //8
    // 						     "150 cm < R < 180 cm"};  //9

    cout<<"Hola -1"<< endl;
    Double_t doubleLatexNamingBinsX = 0.16;
    Double_t doubleLatexNamingBinsRatioX = 0.11;
    Double_t doubleLatexNamingBinsX2 = 0.24;
    Double_t doubleLatexNamingBinsY = 0.9;
    Double_t doubleLatexNamingBinsY2 = 0.86;
    Double_t doubleLatexNamingBinsY3 = 0.82;
    Size_t sizeTextNameBins = 0.05;
    
    Double_t doubleLatexNamingCutX=0.16;
    Double_t doubleLatexNamingCutY=0.20;

    TString 	nameHistoRatioPhiInRMC;
    TString 	nameHistoRatioZInRMC;

    TH1F*	histoMappingRatioPhiInR[nBinsR];
    TH1F*	histoMappingRatioZInR[nBinsR];
    TH1F*	histoMappingRatioR;
    TH1F*	histoMappingMidPtRatioR;
    TH1F*	histoMappingPurityR;
    TH1F*	histoMappingPurityPt;
    TH1F*	histoMappingPurityPt5cm;

    objstrEvent         = (TObjString*)arr->At(0);
    objstrGamma         = (TObjString*)arr->At(1);
    fEventCutSelection      = objstrEvent->GetString();
    fGammaCutSelection      = objstrGamma->GetString();
    cout<<"Hola -2"<< endl;
    fCutSelectionRead = Form("%s_%s",fEventCutSelection.Data(), fGammaCutSelection.Data());
    cout << "Cuts::" <<fEventCutSelection.Data() <<"  "<< fGammaCutSelection.Data()<< endl;

    TString outputDirectoryBase= Form("%s",optionGrid.Data());
    gSystem->Exec("mkdir -p "+outputDirectoryBase);

    TString outputDirectory= Form("%s/%s/%s",optionGrid.Data(),optionPeriod.Data(),fGammaCutSelection.Data());
    //    outputDirectory="Plots";
    gSystem->Exec("mkdir -p "+outputDirectory);


// factor needed for Run1 simulations. LHC10b
//static const Float_t mcGasCorrectionFactor = 0.960035693454;
    Double_t mcGasCorrectionFactor = 1.;

    if (optionPeriod.CompareTo("LHC10b") == 0 || optionPeriod.CompareTo("LHC10bc") == 0 ){
      mcGasCorrectionFactor = 0.960035693454;
    }


    TFile f(fileName.Data());
    TString nameMainDir = "";
    nameMainDir = "GammaConvMaterial";
    TList *TopDir =(TList*)f.Get(nameMainDir.Data());
    if(TopDir == NULL){
        cout<<"ERROR: TopDir not Found"<<endl;
        return;
    }
  
    TList *HistosGammaConversion       = (TList*)TopDir->FindObject(Form("Cut Number %s",fCutSelectionRead.Data()));
    if(HistosGammaConversion == NULL){
        cout<<"ERROR: " << Form("Cut Number %s",fCutSelectionRead.Data()) << " not Found in File"<<endl;
        return;
    }
    TList *ESDContainer                = (TList*)HistosGammaConversion->FindObject(Form("%s ESD histograms",fCutSelectionRead.Data()));



    TFile fMC(fileNameMC.Data());
    TList *TopDirMC =(TList*)fMC.Get(nameMainDir.Data());
    if(TopDirMC == NULL){
        cout<<"ERROR: TopDirMC not Found"<<endl;
        return;
    }
  
    TList *HistosGammaConversionMC       = (TList*)TopDirMC->FindObject(Form("Cut Number %s",fCutSelectionRead.Data()));
    if(HistosGammaConversionMC == NULL){
        cout<<"ERROR: " << Form("Cut Number %s",fCutSelectionRead.Data()) << " not Found in FileMC"<<endl;
        return;
    }

    TList *ESDContainerMC             = (TList*)HistosGammaConversionMC->FindObject(Form("%s ESD histograms",fCutSelectionRead.Data()));
    TList *MCContainer                = (TList*)HistosGammaConversionMC->FindObject(Form("%s MC histograms",fCutSelectionRead.Data()));
    TList *TrueConversionContainerMC  = (TList*)HistosGammaConversionMC->FindObject(Form("%s True histograms",fCutSelectionRead.Data()));
   
    TString nameNTrackHisto = "GoodESDTracks09";
    // TString nameNTrackHisto = "GoodESDTracks14";
    // TString nameNTrackHisto = "GoodESDTracks09_14";

    TH1F* histoEventQualityData=     (TH1F*)ESDContainer->FindObject("NEvents");
    TH2F *histoMappingRPhiData =     (TH2F*)ESDContainer->FindObject("ESD_ConversionMapping_RPhi");
    TH2F *histoMappingRZData =       (TH2F*)ESDContainer->FindObject("ESD_ConversionMapping_RZ");
    TH1F *histoMappingRData    =     (TH1F*)ESDContainer->FindObject("ESD_ConversionMapping_R");
    TH1F *histoMappingMidPtRData =   (TH1F*)ESDContainer->FindObject("ESD_ConversionMappingMidPt_R");
    TH1F *histoMappingPtData =       (TH1F*)ESDContainer->FindObject("ESD_ConversionMapping_Pt");
    TH1F *histoMappingPt5cmData =    (TH1F*)ESDContainer->FindObject("ESD_ConversionMapping_Pt5cm");
    TH2F *histoMappingPtRData =      (TH2F*)ESDContainer->FindObject("ESD_ConversionMapping_Pt_R");
 
    TH1F *histoMappingDCAData    =   (TH1F*)ESDContainer->FindObject("ESD_ConversionMapping_DCA");
    TH1F *histoMappingPsiPairData =  (TH1F*)ESDContainer->FindObject("ESD_ConversionMapping_PsiPair");
    TH1F *histoMappingChi2Data    =  (TH1F*)ESDContainer->FindObject("ESD_ConversionMapping_Chi2");
    TH1F *histoMappingMassData    =  (TH1F*)ESDContainer->FindObject("ESD_ConversionMapping_Mass"); 
    TH1F *histoNumberOfGoodESDTracksData =     (TH1F*)ESDContainer->FindObject(nameNTrackHisto.Data());
    TH1F *histoMappingHighPtRData =   (TH1F*)ESDContainer->FindObject("ESD_ConversionMappingHighPt_R");   

    TH1F* histoEventQualityMC =     (TH1F*)ESDContainerMC->FindObject("NEvents");
    TH2F *histoMappingRPhiMC  =     (TH2F*)ESDContainerMC->FindObject("ESD_ConversionMapping_RPhi");
    TH2F *histoMappingRZMC  =       (TH2F*)ESDContainerMC->FindObject("ESD_ConversionMapping_RZ");
    TH1F *histoMappingRMC     =     (TH1F*)ESDContainerMC->FindObject("ESD_ConversionMapping_R");
    TH1F *histoMappingMidPtRMC =     (TH1F*)ESDContainerMC->FindObject("ESD_ConversionMappingMidPt_R");
    TH1F *histoMappingPtMC =         (TH1F*)ESDContainerMC->FindObject("ESD_ConversionMapping_Pt");
    TH1F *histoMappingPt5cmMC =      (TH1F*)ESDContainerMC->FindObject("ESD_ConversionMapping_Pt5cm");
    TH2F *histoMappingPtRMC =         (TH2F*)ESDContainerMC->FindObject("ESD_ConversionMapping_Pt_R");
    TH1F *histoMappingDCAMC    =      (TH1F*)ESDContainerMC->FindObject("ESD_ConversionMapping_DCA");
    TH1F *histoMappingPsiPairMC =     (TH1F*)ESDContainerMC->FindObject("ESD_ConversionMapping_PsiPair");
    TH1F *histoMappingChi2MC    =     (TH1F*)ESDContainerMC->FindObject("ESD_ConversionMapping_Chi2");
    TH1F *histoMappingMassMC    =     (TH1F*)ESDContainerMC->FindObject("ESD_ConversionMapping_Mass"); 

    TH2F *histoMappingPtTrueMC =         (TH2F*)TrueConversionContainerMC->FindObject("ESD_TrueConversionMapping_Pt");
    TH2F *histoMappingPtTrue5cmMC =      (TH2F*)TrueConversionContainerMC->FindObject("ESD_TrueConversionMapping_Pt5cm");
    TH1F *histoMappingRTrueMC     =      (TH1F*)TrueConversionContainerMC->FindObject("ESD_TrueConversionMapping_R");
    TH1F *histoMappingMidPtRTrueMC     =     (TH1F*)TrueConversionContainerMC->FindObject("ESD_TrueConversionMappingMidPt_R");
    TH1F *histoMappingRTruePi0DalMC     =     (TH1F*)TrueConversionContainerMC->FindObject("ESD_TruePi0DalConversionMapping_R");
    TH1F *histoMappingRTrueEtaDalMC     =     (TH1F*)TrueConversionContainerMC->FindObject("ESD_TrueEtaDalConversionMapping_R");
    TH1F *histoMappingRTrueCombMC     =     (TH1F*)TrueConversionContainerMC->FindObject("ESD_TrueCombConversionMapping_R");
    TH1F *histoNumberOfGoodESDTracksMC  =     (TH1F*)ESDContainerMC->FindObject(nameNTrackHisto.Data());
    TH1F *histoMappingHighPtRMC =     (TH1F*)ESDContainerMC->FindObject("ESD_ConversionMappingHighPt_R");
    
    Double_t ffMinPt=fMinPt;
    Double_t fMaxPt=20.;
    Int_t startBin  = histoMappingPtRData ->GetXaxis()->FindBin(ffMinPt+0.001);
    Int_t endBin    = histoMappingPtRData->GetXaxis()->FindBin(fMaxPt-0.001);
    cout<< "start bin::"<< startBin<< " "<< endBin<<endl;

    histoMappingPtRData ->ProjectionY("histoMappingRPt300Data",startBin,endBin,"e");
    TH1F * histoMappingRPt300Data =(TH1F*)gDirectory->Get("histoMappingRPt300Data");

    histoMappingPtRMC ->ProjectionY("histoMappingRPt300MC",startBin,endBin,"e");
    TH1F * histoMappingRPt300MC =(TH1F*)gDirectory->Get("histoMappingRPt300MC");


    Double_t fMinPt1=0.8;
    Int_t startBin1  = histoMappingPtRData ->GetXaxis()->FindBin(fMinPt1+0.001);
    histoMappingPtRData ->ProjectionY("histoMappingRPt800",startBin1,endBin,"e");
    TH1F * histoMappingRPt800 =(TH1F*)gDirectory->Get("histoMappingRPt800");


    Float_t numberGoodEventsData =          histoEventQualityData->GetBinContent(1);
    Double_t meanMultiplitcityData =        histoNumberOfGoodESDTracksData->GetMean();
    Float_t normFactorReconstData=          1./numberGoodEventsData*1./meanMultiplitcityData;
    // normFactorReconstData=          1./numberGoodEventsData;
    //normFactorReconstData=1;

    Float_t numberGoodEventsMC =          histoEventQualityMC->GetBinContent(1);
    Double_t meanMultiplitcityMC =        histoNumberOfGoodESDTracksMC->GetMean();
    Float_t normFactorReconstMC=          1./numberGoodEventsMC*1./meanMultiplitcityMC;
    // normFactorReconstMC=          1./numberGoodEventsMC;
    //    normFactorReconstMC=1;
    cout<< "Normalization factor Data/ MC::"<< meanMultiplitcityData<< " " << meanMultiplitcityMC<< " " << meanMultiplitcityData/meanMultiplitcityMC<< endl;



    GammaScalingHistogramm(histoNumberOfGoodESDTracksData,normFactorReconstData);
    GammaScalingHistogramm(histoNumberOfGoodESDTracksMC,normFactorReconstMC);

    GammaScalingHistogramm(histoMappingRData,normFactorReconstData);
    GammaScalingHistogramm(histoMappingRPt300Data,normFactorReconstData);
    GammaScalingHistogramm(histoMappingRPt800,normFactorReconstData);

    GammaScalingHistogramm(histoMappingRMC,normFactorReconstMC);
    GammaScalingHistogramm(histoMappingRPt300MC,normFactorReconstData);
    GammaScalingHistogramm(histoMappingRTrueMC,normFactorReconstMC);
    GammaScalingHistogramm(histoMappingRTruePi0DalMC,normFactorReconstMC);
    GammaScalingHistogramm(histoMappingRTrueEtaDalMC,normFactorReconstMC);
    GammaScalingHistogramm(histoMappingRTrueCombMC,normFactorReconstMC);
    GammaScalingHistogramm(histoMappingMidPtRData,normFactorReconstData);
    GammaScalingHistogramm(histoMappingMidPtRMC,normFactorReconstMC);
    GammaScalingHistogramm(histoMappingMidPtRTrueMC,normFactorReconstMC);

    GammaScalingHistogramm(histoMappingPtData,normFactorReconstData);
    GammaScalingHistogramm(histoMappingPtMC,normFactorReconstMC);
    GammaScalingHistogramm(histoMappingPtTrueMC,normFactorReconstMC);
 
    GammaScalingHistogramm(histoMappingPt5cmData,normFactorReconstData);
    GammaScalingHistogramm(histoMappingPt5cmMC,normFactorReconstMC);
    GammaScalingHistogramm(histoMappingPtTrue5cmMC,normFactorReconstMC);


    GammaScalingHistogramm(histoMappingDCAData,normFactorReconstData);
    GammaScalingHistogramm(histoMappingDCAMC,normFactorReconstMC);


    GammaScalingHistogramm(histoMappingPsiPairData,normFactorReconstData);
    GammaScalingHistogramm(histoMappingPsiPairMC,normFactorReconstMC);


    GammaScalingHistogramm(histoMappingChi2Data,normFactorReconstData);
    GammaScalingHistogramm(histoMappingChi2MC,normFactorReconstMC);


    GammaScalingHistogramm(histoMappingMassData,normFactorReconstData);
    GammaScalingHistogramm(histoMappingMassMC,normFactorReconstMC);


    ConvGammaRebinWithBinCorrection(histoMappingRData,rebinRPlots);
    ConvGammaRebinWithBinCorrection(histoMappingRMC,rebinRPlots);
    ConvGammaRebinWithBinCorrection(histoMappingRTrueMC,rebinRPlots);
    ConvGammaRebinWithBinCorrection(histoMappingRTruePi0DalMC,rebinRPlots);
    ConvGammaRebinWithBinCorrection(histoMappingRTrueEtaDalMC,rebinRPlots);
    ConvGammaRebinWithBinCorrection(histoMappingRTrueCombMC,rebinRPlots);
    ConvGammaRebinWithBinCorrection(histoMappingMidPtRData,rebinRPlots);
    ConvGammaRebinWithBinCorrection(histoMappingMidPtRMC,rebinRPlots);
    ConvGammaRebinWithBinCorrection(histoMappingMidPtRTrueMC,rebinRPlots);


    ConvGammaRebinWithBinCorrection(histoMappingPtData,rebinPtPlots);
    ConvGammaRebinWithBinCorrection(histoMappingPtMC,rebinPtPlots);
    ConvGammaRebinWithBinCorrection(histoMappingPtTrueMC,rebinPtPlots);



    // Create set of histograms normalized with conversions in the gas Data and MC


    TH1F * histoMappingRDataScaledToGas = ScaleByIntegralWithinLimits(histoMappingRData, rMinGas, rMaxGas, 1. );
    TH1F * histoMappingRPt300DataScaledToGas  = ScaleByIntegralWithinLimits(histoMappingRPt300Data, rMinGas, rMaxGas, 1. );

    TH1F * histoMappingRMCScaledToGas = ScaleByIntegralWithinLimits(histoMappingRMC, rMinGas, rMaxGas, mcGasCorrectionFactor );
    TH1F * histoMappingRPt300MCScaledToGas  = ScaleByIntegralWithinLimits(histoMappingRPt300MC, rMinGas, rMaxGas, mcGasCorrectionFactor );


    TH1F * histoMappingRDataScaledToGasRebin = (TH1F*)histoMappingRDataScaledToGas->Rebin(nBinsR,"histoMappingRDataScaledToGasRebin", arrayRBins);
    TH1F * histoMappingRMCScaledToGasRebin = (TH1F*)histoMappingRMCScaledToGas->Rebin(nBinsR,"histoMappingRMCScaledToGasRebin", arrayRBins);

    TH1F * histoMappingRPt300DataScaledToGasRebin = (TH1F*)histoMappingRPt300DataScaledToGas->Rebin(nBinsR,"histoMappingRPt300DataScaledToGasRebin", arrayRBins);
    TH1F * histoMappingRPt300MCScaledToGasRebin = (TH1F*)histoMappingRPt300MCScaledToGas->Rebin(nBinsR,"histoMappingRPt300MCScaledToGasRebin", arrayRBins);






    TH1F* histoMappingDoubleRatioR= (TH1F*)histoMappingRDataScaledToGasRebin->Clone();
    histoMappingDoubleRatioR->SetName("histoMappingDoubleRatioR"); 
    histoMappingDoubleRatioR->Divide(histoMappingRDataScaledToGasRebin,histoMappingRMCScaledToGasRebin,1.,1.,"B");

    TH1F* histoMappingDoubleRatioRPt300= (TH1F*)histoMappingRPt300DataScaledToGasRebin->Clone();
    histoMappingDoubleRatioRPt300->SetName("histoMappingDoubleRatioRPt300"); 
    histoMappingDoubleRatioRPt300->Divide(histoMappingRPt300DataScaledToGasRebin,histoMappingRPt300MCScaledToGasRebin,1.,1.,"B");
  


    // Start plotting

    TLine * 	linePhi =	 		new TLine (0,1,2*TMath::Pi(),1);
    linePhi->SetLineColor(1);
 

    TLine * 	lineZ =	 		new TLine (-100,1,100.,1.);
    lineZ->SetLineColor(1);
    TLine * 	lineR =	 		new TLine (0,1,190.,1);
    lineR->SetLineColor(1);
    TLine * 	linePt =	 		new TLine (0,1,20.,1);
    lineR->SetLineColor(1);
 

 
  
    // nameHistoRatioR=Form("histoMappingRatioR");
    histoMappingRatioR= (TH1F*)histoMappingRData->Clone();
    histoMappingRatioR->SetName("histoMappingRatioR"); 
    histoMappingRatioR->Divide(histoMappingRatioR,histoMappingRMC,1.,1.,"B");


    histoMappingMidPtRatioR= (TH1F*)histoMappingMidPtRData->Clone();
    histoMappingMidPtRatioR->SetName("histoMappingRatioR"); 
    histoMappingMidPtRatioR->Divide(histoMappingMidPtRatioR,histoMappingMidPtRMC,1.,1.,"B");


    histoMappingPurityR = (TH1F*)histoMappingRTrueMC->Clone();
    histoMappingPurityR->Sumw2();
    histoMappingPurityR->SetName("histoMappingPurityR"); 
    histoMappingPurityR->Divide(histoMappingPurityR,histoMappingRMC,1.,1.,"B");

    histoMappingPurityPt = (TH1F*)histoMappingPtTrueMC->Clone();
    histoMappingPurityPt->Sumw2();
    histoMappingPurityPt->SetName("histoMappingPurityPt"); 
    histoMappingPurityPt->Divide(histoMappingPurityPt,histoMappingPtMC,1.,1.,"B");

    histoMappingPurityPt5cm = (TH1F*)histoMappingPtTrue5cmMC->Clone();
    histoMappingPurityPt5cm->Sumw2();
    histoMappingPurityPt5cm->SetName("histoMappingPurityPt5cm"); 
    histoMappingPurityPt5cm->Divide(histoMappingPurityPt5cm,histoMappingPt5cmMC,1.,1.,"B");


    // -Multiplicity----
    TCanvas * canvasNTracks = new TCanvas("canvasNTracks","",1200,1000);  // gives the page size
    DrawGammaCanvasSettings( canvasNTracks, 0.12, 0.015, 0.095, 0.09);
    //   canvasNTracks->Divide(2,2);
    canvasNTracks->SetLogy(1);
    


    DrawAutoGammaHistosMaterial( histoNumberOfGoodESDTracksData, 
				 histoNumberOfGoodESDTracksMC, 
			 "","Good ESD Tracks ","Counts",
			 kTRUE, 1.2,2.e-8,
			 kFALSE, 0.,0.,
			 kTRUE, 0.,100.);
    histoNumberOfGoodESDTracksMC->DrawCopy("same,hist");
    canvasNTracks->Print(Form("%s/NumberOfGoodESDTracks_%s_%s.pdf",outputDirectory.Data(),optionPeriod.Data(),fGammaCutSelection.Data() ));



    TCanvas * canvasScaledToGas= new TCanvas("canvasNTracks","",1200,1000);
    // histoMappingRDataScaledToGas->DrawCopy();
    //histoMappingRMCScaledToGas->DrawCopy("same");
    canvasScaledToGas->Divide(1,2);
    canvasScaledToGas->cd(1);
    DrawAutoGammaHistosMaterial( histoMappingRDataScaledToGas, 
    				 histoMappingRMCScaledToGas, 
    				 "","R (cm) ","Counts",
    				 kFALSE, 1.2,2.e-8,
    				 kFALSE, 0.,0.,
    				 kTRUE, 0.,180.);
    canvasScaledToGas->cd(2);
    histoMappingDoubleRatioR->SetAxisRange(0.,179.);
    histoMappingDoubleRatioR->SetMinimum(0.75);
    histoMappingDoubleRatioR->SetMaximum(1.25);
    histoMappingDoubleRatioR->DrawCopy();
    lineR->Draw("same");

    histoMappingDoubleRatioRPt300->SetLineColor(3);
    histoMappingDoubleRatioRPt300->DrawCopy("same");


    // canvasScaledToGas->cd(3);
    // fProfileContainingMaterialBudgetWeights->Draw();
    // fProfileContainingMaterialBudgetWeights->Print();



    canvasScaledToGas->Print(Form("%s/RScaledtoGas_%s_%s.pdf",outputDirectory.Data(),optionPeriod.Data(),fGammaCutSelection.Data() ));

    // TCanvas * canvasRScaledToGas= new TCanvas("canvasRPlot","",1200,1000);
    // canvasRScaledToGas->SetLogy(1);
    // histoMappingRData->SetXTitle("R(cm)");
    // histoMappingRData->SetMaximum(2e-3);
    // histoMappingRData->SetMinimum(1e-6);
    // histoMappingRData->SetLineColor(4);
    // histoMappingRData->DrawCopy("hist");
    // histoMappingRPt300Data->SetLineColor(2);
    // histoMappingRPt300Data->DrawCopy("hist,same");
    // histoMappingRPt800->SetLineColor(6);
    // histoMappingRPt800->DrawCopy("hist,same");
    // canvasRScaledToGas->Print(Form("%s/RScaledtoGasData_%s.pdf",outputDirectory.Data(),fGammaCutSelection.Data() ));




    //----------------------------------------------

    TLatex* latexCutName= new TLatex(doubleLatexNamingCutX, doubleLatexNamingCutY,fGammaCutSelection.Data()); 
    SetStyleTLatex( latexCutName, sizeTextNameBins,2);    
  
    TCanvas * canvasPhotonR = new TCanvas("canvasPhotonR","",1200,1000);  // gives the page size
    DrawGammaCanvasSettings( canvasPhotonR, 0.12, 0.015, 0.095, 0.09);
    canvasPhotonR->Divide(2,2);
    canvasPhotonR->SetLogy(1);
    
    canvasPhotonR->cd(1);
    canvasPhotonR->cd(1)->SetLogx(1);

    DrawAutoGammaHistosMaterial( histoMappingRData, 
				 histoMappingRMC, 
			 "","R (cm) ","Counts",
			 kTRUE, 1.2,2.e-8,
			 kTRUE, 0.,0.,
			 kTRUE, 0.,190.);
    histoMappingRTrueMC->DrawCopy("same,hist");
    histoMappingRTruePi0DalMC->SetLineColor(kBlue-9);
    histoMappingRTruePi0DalMC->SetFillColor(kBlue-9);
    histoMappingRTruePi0DalMC->SetFillStyle(3244);
    histoMappingRTruePi0DalMC->DrawCopy("same,hist");

    histoMappingRTrueEtaDalMC->SetLineColor(kBlue-3);
    histoMappingRTrueEtaDalMC->SetFillColor(kBlue-3);
    histoMappingRTrueEtaDalMC->SetFillStyle(3002);
    histoMappingRTrueEtaDalMC->DrawCopy("same,hist");
    histoMappingRTrueCombMC->DrawCopy("same,hist");

    histoMappingMidPtRData->SetLineStyle(2);
    histoMappingMidPtRData->DrawCopy("same,histo");
    
    histoMappingMidPtRMC->SetLineStyle(2);
    histoMappingMidPtRMC->SetLineColor(kRed);
    histoMappingMidPtRMC->DrawCopy("same,histo");

    histoMappingMidPtRTrueMC->SetLineColor(kBlue-9);
    histoMappingMidPtRTrueMC->SetLineStyle(2);
    histoMappingMidPtRTrueMC->DrawCopy("same,histo");


    latexCutName->Draw("same");

    for(Int_t iR = 1; iR < nBinsR; iR++){
      lineRLimits[iR] = new TLine (arrayRBins[iR],2e-6,arrayRBins[iR],1.03e-3);
      lineRLimits[iR] ->SetLineColor(1);
      lineRLimits[iR] ->SetLineWidth(0.1);
      lineRLimits[iR] ->Draw("same");
      lineRLimitsPurity[iR] = new TLine (arrayRBins[iR],0.5,arrayRBins[iR],1.5);
    }
    
    canvasPhotonR->cd(2);
    
    DrawAutoGammaHisto( histoMappingRatioR, 
			"","R (cm) ","Ratio",
			kTRUE, 1.2,0.5,
			kFALSE,0. ,0.,
			kTRUE, 0.,190.);
    lineR->Draw("same");
    
    for(Int_t iR = 1; iR < nBinsR; iR++){
      lineRLimitsPurity[iR] ->Draw("same");
    }
    
    latexCutName->Draw("same");
    histoMappingMidPtRatioR->SetLineColor(kMagenta);
    //histoMappingMidPtRatioR->SetLineStyle(2);
    histoMappingMidPtRatioR->DrawCopy("same,histo");
    
    canvasPhotonR->cd(3);
    canvasPhotonR->cd(3)->SetLogy(1);
    DrawAutoGammaHistosMaterial( histoMappingRData, 
			 histoMappingRMC, 
			 "","R (cm) ","Counts",
			 kTRUE, 1.2,2.e-6,
			 kTRUE,1.e-8 ,1e-3,
			 kTRUE, 0.,190.);

    histoMappingRTrueMC->DrawCopy("same,hist");
    histoMappingRTruePi0DalMC->DrawCopy("same,hist");
    histoMappingRTrueEtaDalMC->DrawCopy("same,hist");
    histoMappingRTrueCombMC->DrawCopy("same,hist");

    latexCutName->Draw("same");

    for(Int_t iR = 1; iR < nBinsR; iR++){
      lineRLimits[iR] ->Draw("same");
    }

    histoMappingMidPtRData->DrawCopy("same,histo");
    histoMappingMidPtRMC->DrawCopy("same,histo");
    histoMappingMidPtRTrueMC->DrawCopy("same,histo");

                
    canvasPhotonR->cd(4);
    canvasPhotonR->cd(4)->SetLogy(1);
    DrawAutoGammaHisto( histoMappingRatioR, 
			"","R (cm) ","Ratio",
			 kTRUE, 1.2,0.15,
			 kFALSE,0. ,0.,
			 kTRUE, 0.,190.);
  
    lineR->Draw("same");
    for(Int_t iR = 1; iR < nBinsR; iR++){
      lineRLimitsPurity[iR] ->Draw("same");
    }

    latexCutName->Draw("same");
    canvasPhotonR->Print(Form("%s/PhotonRadius_%s_%s.pdf",outputDirectory.Data(),optionPeriod.Data(),fGammaCutSelection.Data() ));






  TCanvas * canvasPhotonPurity = new TCanvas("canvasPhotonPurity","",1200,1000);  // gives the page size
    DrawGammaCanvasSettings( canvasPhotonPurity, 0.12, 0.015, 0.095, 0.09);
    canvasPhotonPurity->Divide(1,2);
    canvasPhotonPurity->cd(1);
    DrawAutoGammaHisto( histoMappingPurityR,
                         "","R (cm) ","Purity",
			 kFALSE, 0., 0.,
			 kTRUE,0.8,1.05,
			 kTRUE,0.,190.);
    lineR->Draw("same");

    for(Int_t iR = 1; iR < nBinsR; iR++){
      lineRLimitsPurity[iR] ->Draw("same");
    }
    latexCutName->Draw("same");

    //    histoMappingPurityR->DrawCopy();
    canvasPhotonPurity->cd(2);

    //    histoMappingPurityPt->DrawCopy();
    DrawAutoGammaHistosMaterial( histoMappingPurityPt,
				 histoMappingPurityPt5cm,
				 "","Pt (GeV/c) ","Purity",
			 kFALSE, 0., 0.,
			 kTRUE,0.8 ,1.05,
			 kTRUE,0.,10.);

    linePt->Draw("same");
    latexCutName->Draw("same");
    canvasPhotonPurity->Print(Form("%s/PhotonPurity_%s_%s.pdf",outputDirectory.Data(),optionPeriod.Data(),fGammaCutSelection.Data() ));

   

    Int_t rebinPhiPlots=3;



    TCanvas * canvasPhotonChar = new TCanvas("canvasPhotonChar","",1200,1000);  // gives the page size
    DrawGammaCanvasSettings( canvasPhotonChar, 0.12, 0.015, 0.075, 0.09);
    canvasPhotonChar->Divide(2,2);
    canvasPhotonChar->SetLogy(1);
    
    canvasPhotonChar->cd(1);
    canvasPhotonChar->cd(1)->SetLogy(1);
  
    
    DrawAutoGammaHistosMaterial( histoMappingDCAData, 
			 histoMappingDCAMC, 
			 "","DCA ","Counts",
			 kTRUE, 1.2,2.e-8,
			 kFALSE,0. ,0.,
			 kTRUE, 0.,190.);

    latexCutName->Draw("same");

    canvasPhotonChar->cd(2);
    canvasPhotonChar->cd(2)->SetLogy(1);
    DrawAutoGammaHistosMaterial( histoMappingPsiPairData, 
			 histoMappingPsiPairMC, 
			 "","PsiPair ","Counts",
			 kTRUE, 1.2,2.e-8,
			 kFALSE,0. ,0.,
			 kTRUE, 0.,190.);
    latexCutName->Draw("same");                

    canvasPhotonChar->cd(3);
    canvasPhotonChar->cd(3)->SetLogy(1);
    DrawAutoGammaHistosMaterial( histoMappingChi2Data, 
			 histoMappingChi2MC, 
			 "","Chi2/NDF ","Counts",
			 kTRUE, 1.2,2.e-8,
			 kFALSE,0. ,0.,
			 kTRUE, 0.,190.);
                
    latexCutName->Draw("same");
    canvasPhotonChar->cd(4);
    canvasPhotonChar->cd(4)->SetLogy(1);
    DrawAutoGammaHistosMaterial( histoMappingMassData, 
			 histoMappingMassMC, 
			 "","Mass (GeV) ","Counts",
			 kTRUE, 1.2,2.e-8,
			 kFALSE,0. ,0.,
			 kTRUE, 0.,190.);

    latexCutName->Draw("same");
    canvasPhotonChar->Print(Form("%s/PhotonChar_%s_%s.pdf",outputDirectory.Data(),optionPeriod.Data(),fGammaCutSelection.Data() ));


    //    GammaScalingHistogramm(histoMappingPtData,normFactorReconstData);
    //   GammaScalingHistogramm(histoMappingPtMC,normFactorReconstMC);
    //   GammaScalingHistogramm(histoMappingPtTrueMC,normFactorReconstMC);



    // histoMappingRData->SetLineColor(1);
    // histoMappingRData->DrawCopy("hist");
    // histoMappingRMC->SetLineColor(2);
    // histoMappingRMC->DrawCopy("hist,same");
 

    TH1D*   histoMappingPhiInRData[nBinsR];
    TH1D*   histoMappingPhiInRMC[nBinsR];

    for(Int_t iR = 0; iR < nBinsR; iR++){
      histoMappingPhiInRData[iR] =  histoMappingRPhiData->ProjectionX( Form("histoMappingPhiInRData_%i",iR), histoMappingRPhiData->GetYaxis()->FindBin(arrayRBins[iR]), histoMappingRPhiData->GetYaxis()->FindBin(arrayRBins[iR+1]));
      GammaScalingHistogramm(histoMappingPhiInRData[iR],normFactorReconstData);
      ConvGammaRebinWithBinCorrection(histoMappingPhiInRData[iR],rebinPhiPlots);
      // ConvGammaRebinWithBinCorrection(histoMappingPhiInRData[iR],rebinPhiPlots);
      //GammaScalingHistogramm(histoMappingPhiInRData[iR],normFactorReconstData*1/600*1/(TMath::Abs(histoMappingRPhiData->GetXaxis()->GetBinUpEdge(histoMappingRPhiData->GetXaxis()->FindBin(arrayRBins[iR+1]))-histoMappingRPhiData->GetXaxis()->GetBinLowEdge(histoMappingRPhiData->GetXaxis()->FindBin(arrayRBins[iR])))));





      histoMappingPhiInRMC[iR] =  histoMappingRPhiMC->ProjectionX( Form("histoMappingPhiInRMC_%i",iR), histoMappingRPhiMC->GetYaxis()->FindBin(arrayRBins[iR]), histoMappingRPhiMC->GetYaxis()->FindBin(arrayRBins[iR+1]));
      ConvGammaRebinWithBinCorrection(histoMappingPhiInRMC[iR],rebinPhiPlots);
      GammaScalingHistogramm(histoMappingPhiInRMC[iR],normFactorReconstMC);

      
      nameHistoRatioPhiInRMC=Form("histoMappingRatioPhiInR_%02d",iR);
      histoMappingRatioPhiInR[iR]= (TH1F*)histoMappingPhiInRData[iR]->Clone();
      histoMappingRatioPhiInR[iR]->SetName(nameHistoRatioPhiInRMC); 
      histoMappingRatioPhiInR[iR]->Divide(histoMappingRatioPhiInR[iR],histoMappingPhiInRMC[iR]);
      
      TCanvas * canvasSingleBin = new TCanvas("canvasSingleBin","",500,350);  // gives the page size
      //if (iR == 11)   canvasSingleBin->SetLogy(1);
      canvasSingleBin->Divide(1,2);
      canvasSingleBin->SetTopMargin(0.07);
      canvasSingleBin->cd(1);
      TString nameHistoPhiInRData=Form("#varphi in R:  %s",arrayNamesRangesRBins[iR].Data());
      TLatex *latexBinning = new TLatex(doubleLatexNamingBinsX, doubleLatexNamingBinsY,nameHistoPhiInRData); 
      SetStyleTLatex( latexBinning, sizeTextNameBins,2);
      TLatex *latexBinning2 = new TLatex(doubleLatexNamingBinsX2, doubleLatexNamingBinsY2,arrayNamesRBins[iR]); 
      SetStyleTLatex( latexBinning2, sizeTextNameBins,2);
      


      DrawAutoGammaHistosMaterial( histoMappingPhiInRData[iR], 
			    histoMappingPhiInRMC[iR],
			   //	    histoTrueMappingPhiInRMonteCarlo[iR],
			    "","#varphi","Counts",
			    kTRUE,1.2 ,1e-10,
			    kFALSE,0. ,0.,
			    kFALSE, 0.,200.);

      latexBinning->Draw("same");  
      latexBinning2->Draw("same");  
      latexCutName->Draw("same");

      canvasSingleBin->cd(2);

      DrawAutoGammaHistosMaterial( histoMappingRatioPhiInR[iR], 
				   histoMappingRatioPhiInR[iR],
				   //	    histoTrueMappingPhiInRMonteCarlo[iR],
				   "","#varphi","Ratio",
				   kTRUE,1.2 ,1e-10,
				   kFALSE,0. ,0.,
				   kFALSE, 0.,200.);
      linePhi->Draw("same");
      latexCutName->Draw("same");

      // if (outer) {
      // 	TLatex *latexBinning3 = new TLatex(doubleLatexNamingBinsX2, doubleLatexNamingBinsY3,arrayNamesAddRBins[iR]);
      // 	SetStyleTLatex( latexBinning3, sizeTextNameBins,2);
      // 	latexBinning3->Draw("same"); 
      // }
      
      //                      if(!thesis)DrawAliceText(floatLocationRightUpText[0],floatLocationRightUpText[1], floatLocationRightUpText[3]);             
      
      canvasSingleBin->Update();
      canvasSingleBin->SaveAs(Form("%s/Phi_in_R_%i_%s_%s.pdf",outputDirectory.Data(), iR,optionPeriod.Data(),fGammaCutSelection.Data()));
      delete canvasSingleBin;
      
    }

    TH1D*   histoMappingZInRData[nBinsR];
    TH1D*   histoMappingZInRMC[nBinsR];

    Float_t materialMarkerSize=0.015;

    for(Int_t iR = 0; iR < nBinsR; iR++){
      if (iR>9) rebinZPlots=4;
      histoMappingZInRData[iR] =  histoMappingRZData->ProjectionX( Form("histoMappingZInRData_%i",iR), histoMappingRZData->GetYaxis()->FindBin(arrayRBins[iR]), histoMappingRZData->GetYaxis()->FindBin(arrayRBins[iR+1]));
      GammaScalingHistogramm(histoMappingZInRData[iR],normFactorReconstData);
      ConvGammaRebinWithBinCorrection(histoMappingZInRData[iR],rebinZPlots);
      
      histoMappingZInRMC[iR] =  histoMappingRZMC->ProjectionX( Form("histoMappingZInRMC_%i",iR), histoMappingRZMC->GetYaxis()->FindBin(arrayRBins[iR]), histoMappingRZMC->GetYaxis()->FindBin(arrayRBins[iR+1]));
      ConvGammaRebinWithBinCorrection(histoMappingZInRMC[iR],rebinZPlots);
      GammaScalingHistogramm(histoMappingZInRMC[iR],normFactorReconstMC);

      
      nameHistoRatioZInRMC=Form("histoMappingRatioZInR_%02d",iR);
      histoMappingRatioZInR[iR]= (TH1F*)histoMappingZInRData[iR]->Clone();
      histoMappingRatioZInR[iR]->SetName(nameHistoRatioZInRMC); 
      histoMappingRatioZInR[iR]->Divide(histoMappingRatioZInR[iR],histoMappingZInRMC[iR]);
     
      TCanvas * canvasSingleBin = new TCanvas("canvasSingleBin","",500,350);  // gives the page size
      //if (iR == 11)   canvasSingleBin->SetLogy(1);
      canvasSingleBin->Divide(1,2);
      canvasSingleBin->SetTopMargin(0.07);
      canvasSingleBin->cd(1);
      //   canvasSingleBin->cd(1)->SetLogy(1);
      TString nameHistoZInRData=Form("#varphi in R:  %s",arrayNamesRangesRBins[iR].Data());
      TLatex *latexBinning = new TLatex(doubleLatexNamingBinsX, doubleLatexNamingBinsY,nameHistoZInRData); 
      SetStyleTLatex( latexBinning, sizeTextNameBins,2);
      TLatex *latexBinning2 = new TLatex(doubleLatexNamingBinsX2, doubleLatexNamingBinsY2,arrayNamesRBins[iR]); 
      SetStyleTLatex( latexBinning2, sizeTextNameBins,2);
      
      StylingSliceHistos(histoMappingZInRData[iR], materialMarkerSize);
      //      StylingSliceHistos(histoMappingZInRMC[iR], materialMarkerSize);

      DrawAutoGammaHistosMaterial( histoMappingZInRData[iR], 
			    histoMappingZInRMC[iR],
			   //	    histoTrueMappingPhiInRMonteCarlo[iR],
			    "","#varphi","counts",
			    kTRUE,1.2 ,1e-10,
			    kFALSE,0. ,0.,
			    kTRUE, -arrayZBins[iR],arrayZBins[iR]);

      latexBinning->Draw("same");  
      latexBinning2->Draw("same");  
      latexCutName->Draw("same");

      canvasSingleBin->cd(2);

      DrawAutoGammaHistosMaterial( histoMappingRatioZInR[iR], 
				   histoMappingRatioZInR[iR],
				   //	    histoTrueMappingPhiInRMonteCarlo[iR],
				   "","#varphi","Ratio",
				   kFALSE,1.2 ,1e-10,
				   kTRUE,0.4 ,2.,
				   kTRUE,-arrayZBins[iR],arrayZBins[iR] );
      lineZ->Draw("same");
      latexCutName->Draw("same");
      canvasSingleBin->Update();
      canvasSingleBin->SaveAs(Form("%s/Z_in_R_%i_%s_%s.pdf",outputDirectory.Data(), iR,optionPeriod.Data(),fGammaCutSelection.Data()));
      delete canvasSingleBin;
    }

    // for(Int_t iR = 0; iR < nBinsR; iR++){
    //   new TCanvas;
    //   histoMappingPhiInRMC[iR]->SetLineColor(2);
    //   histoMappingPhiInRMC[iR]->DrawCopy("hist");
    // }




   TProfile* fProfileContainingMaterialBudgetWeights = new TProfile("profileContainingMaterialBudgetWeights_manyRadialBins","profileContainingMaterialBudgetWeights_manyRadialBins",nBinsR,arrayRBins);
    fProfileContainingMaterialBudgetWeights->GetXaxis()->SetTitle("R (cm)");
    fProfileContainingMaterialBudgetWeights->GetYaxis()->SetTitle("weight_per_gamma");
   
    for(Int_t i=0; i<nBinsR; i++){
      cout<< arrayRBins[i] << " " << histoMappingDoubleRatioR->GetBinContent(i+1) <<endl;
      fProfileContainingMaterialBudgetWeights->Fill(arrayRBins[i],histoMappingDoubleRatioR->GetBinContent(i+1));
    }

    cout<< " pT 300"<< endl; 
   TProfile* fProfileContainingMaterialBudgetWeightsPt300 = new TProfile("profileContainingMaterialBudgetWeightsPt300_manyRadialBins","profileContainingMaterialBudgetWeightsPt300_manyRadialBins",nBinsR,arrayRBins);
    fProfileContainingMaterialBudgetWeightsPt300->GetXaxis()->SetTitle("R (cm)");
    fProfileContainingMaterialBudgetWeightsPt300->GetYaxis()->SetTitle("weight_per_gamma");
   
    for(Int_t i=0; i<nBinsR; i++){
      cout<< arrayRBins[i] << " " << histoMappingDoubleRatioRPt300->GetBinContent(i+1) <<endl;
      fProfileContainingMaterialBudgetWeightsPt300->Fill(arrayRBins[i],histoMappingDoubleRatioRPt300->GetBinContent(i+1));
    }
    //Form("%s/PhotonChar_%s.pdf",optionPeriod.Data(),fGammaCutSelection.Data() ))
    cout<< Form("MCInputFileMaterialBudgetWeights_%s_%s.root",optionPeriod.Data(),fGammaCutSelection.Data())<< endl;
    TFile outFile(Form("MCInputFileMaterialBudgetWeights_%s_%s.root",optionPeriod.Data(),fGammaCutSelection.Data()) ,"RECREATE");
    fProfileContainingMaterialBudgetWeights->Write();
    fProfileContainingMaterialBudgetWeightsPt300->Write();

    histoMappingRData->Write("Data");
    histoMappingRMC->Write("MC");
    histoMappingMidPtRData->Write("DataMidtPt");
    histoMappingMidPtRMC->Write("MCMidPt");
    histoMappingRDataScaledToGas->Write("DataScaledToGas");
    histoMappingRMCScaledToGas ->Write("MCScaledToGas");
    outFile.Close();



} 
