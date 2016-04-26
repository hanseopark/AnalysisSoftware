
/****************************************************************************************************************************
******      provided by Gamma Conversion Group, PWG4,                                        *****
******      Ana Marin, marin@physi.uni-heidelberg.de                                      *****
******         Kathrin Koch, kkoch@physi.uni-heidelberg.de                                      *****
******      Friederike Bock, friederike.bock@cern.ch                                      *****
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
#include "CommonHeaders/CombinationFunctions.h"

extern TRandom*   gRandom;
extern TBenchmark*   gBenchmark;
extern TSystem*   gSystem;
extern TMinuit*   gMinuit;

void CompareChargedAndNeutralPionDataALICEpPb_FullpTMB(TString outputDir = "pdf/CombineMesonMeasurementspPb", TString suffix = "pdf"){

   gROOT->Reset();   
   gROOT->SetStyle("Plain");
   
   StyleSettingsThesis();  
   SetPlotStyle();
   
   Double_t xSection2760GeVpp =     55.416*1e-3;
   Double_t xSection2760GeVErrpp =  3.9;
   Double_t xSection2760GeVppINEL = 62.8*1e9;
   Double_t xSection900GeVppINEL = 52.5*1e9;
   Double_t xSection7TeVppINEL = 73.2*1e9;   
   Double_t recalcBarn =         1e12; //NLO in pbarn!!!!

   gSystem->Exec("mkdir -p "+outputDir);
   
   Width_t  widthLinesBoxes;
   
   TString collisionSystempPb = "p-Pb #sqrt{s_{NN}} = 5.02 TeV";     
   
   Size_t markerSizeComparison = 1.5;
   TString nameHistoPCM = "CorrectedYieldPi0";
   TString nameHistoDalitz = "Pi0SystError";
   TString nameGraphPCM = "Pi0SystError";
   

   // TFile* fileNeutralPionPCMDatapPb = new TFile("data_PCMResults_pPb_20141023.root");
   TFile* fileNeutralPionPCMDatapPb = new TFile("data_PCMResults_pPb_20150624_dc4.root");
   //  TFile* fileNeutralPionPCMDatapPb = new TFile("data_PCMResults_pPb_RCut_20140424MB_CatAPileup.root");
   TFile* fileNeutralPionDalitzpPb = new TFile("ExternalInputpPb/data_PCMResults_Dalitz_pPb_2014-09-12.root");
   //   TFile*   filePHOSpPb =       new TFile("ExternalInputpPb/data_PHOSResultsFullCorrection_pPb_140422.root");
   TFile*   filePHOSpPb =       new TFile("ExternalInputpPb/PHOS/data_PHOSResultsFullCorrection_pPb-23102014.root");
   //  TFile*   filePHOSpPbSys =       new TFile("ExternalInputpPb/sys_err_PHOS_20140130.root");
   TString fileNameEMCalpPb = "ExternalInputpPb/EMCAL/EMCALResults_11Sept.root"; 

	TString dateDalitz = "2014-09-12";
	TString datePCM = "2015-06-24";
	TString datePHOS = "2014-10-23";
	TString dateEMCal = "2014-09-11";

    //****************************************************************************************************
   //************************** Read data for PCM *******************************************************
   //****************************************************************************************************
   TDirectory* directoryPCMPi0pPb =             (TDirectory*)fileNeutralPionPCMDatapPb->Get("Pi0_pPb_5.023TeV_0-100%"); 
   TH1D* histoPCMYieldPi0pPb =            (TH1D*)directoryPCMPi0pPb->Get(nameHistoPCM.Data());
   TGraphAsymmErrors* graphPCMYieldPi0SysErrpPb=    (TGraphAsymmErrors*)directoryPCMPi0pPb->Get(nameGraphPCM.Data()); 
   //****************************************************************************************************
   //************************** Read data for Dalitz *******************************************************
   //****************************************************************************************************  

  TDirectory* directoryDalitzPi0pPb =             (TDirectory*)fileNeutralPionDalitzpPb->Get("Pi0_pPb_5.023TeV_0-100%"); 
  //   TH1D* histoDalitzYieldPi0pPb =            (TH1D*)directoryDalitzPi0pPb->Get(nameHistoPCM.Data());
   TGraphAsymmErrors* histoDalitzYieldPi0pPb =            (TGraphAsymmErrors*)directoryDalitzPi0pPb->Get(nameHistoDalitz.Data());
   TGraphAsymmErrors* graphDalitzYieldPi0SysErrpPb=    (TGraphAsymmErrors*)directoryDalitzPi0pPb->Get(nameGraphPCM.Data()); 
   histoDalitzYieldPi0pPb->RemovePoint(0);
   graphDalitzYieldPi0SysErrpPb->RemovePoint(0);
   //****************************************************************************************************
   //************************** Read data for PHOS *******************************************************
   //****************************************************************************************************
 
   TDirectory *directoryPHOSPi0pPb = 	(TDirectory*)filePHOSpPb->Get("Pi0_pPb_5.023TeV_0-100%"); 
   TH1D* histoPHOSYieldPi0pPb = 	(TH1D*)directoryPHOSPi0pPb->Get(nameHistoPCM.Data());      
    
   //  TH1D* histoPHOSYieldPi0pPbSys = 	(TH1D*)filePHOSpPbSys->Get("yeild1_GS_All_cen0");      
   TH1D* histoPHOSYieldPi0pPbSys = 	(TH1D*)directoryPHOSPi0pPb->Get("Pi0SystError");      
   Int_t PHOS_start=histoPHOSYieldPi0pPbSys->FindBin(2.3);
   Int_t PHOS_stopp=histoPHOSYieldPi0pPbSys->FindBin(19.);
   Int_t PHOS_Bin=histoPHOSYieldPi0pPbSys->GetNbinsX();
   cout <<"PHOS range "  <<  PHOS_start << "   "<< PHOS_stopp ;
   for (int i=1; i<PHOS_start;i++ ){
     histoPHOSYieldPi0pPbSys->SetBinError(i,0.);
     histoPHOSYieldPi0pPbSys->SetBinContent(i,0.);
    histoPHOSYieldPi0pPb->SetBinError(i,0.);
     histoPHOSYieldPi0pPb->SetBinContent(i,0.);
   }
   for (int i=PHOS_stopp; i<=PHOS_Bin;i++ ){
     histoPHOSYieldPi0pPbSys->SetBinError(i,0.);
     histoPHOSYieldPi0pPbSys->SetBinContent(i,0.);
    histoPHOSYieldPi0pPb->SetBinError(i,0.);
     histoPHOSYieldPi0pPb->SetBinContent(i,0.);
   }  	
			
   //****************************************************************************************************
   //************************** Read data for EMCal *******************************************************
   //****************************************************************************************************
	TString nameHistoEMCal = "CorrectedYieldPi0";
	TString nameHistoEMCalSysErrors = "Pi0SystError";
   TFile*	fileEMCalpPb = 					new TFile(fileNameEMCalpPb);
   TDirectory*	directoryEMCalPi0pPb = 				(TDirectory*)fileEMCalpPb->Get("Pi05.02TeV_pPb"); 
   TH1D*	histoEMCalYieldPi0pPb = 				(TH1D*)directoryEMCalPi0pPb->Get(nameHistoEMCal.Data());
   TH1D*   	histoEMCalYieldPi0pPbSys= 	(TH1D*)directoryEMCalPi0pPb->Get(nameHistoEMCalSysErrors.Data());



  
   cout << "*************************************************************************"<< endl;  
   cout << "****************** charged pions Full pT ********************************"<< endl;
   cout << "*************************************************************************"<< endl;



 TFile* fileChargedPionsFullpT = new TFile("ExternalInputpPb/FullPreliminaryPionSpectra.root");
   TH1D* histoChargedPionSys0_5= (TH1D*)fileChargedPionsFullpT->Get("hPionSpectrum_Syst_0_5"); //0-5% 
   TH1D* histoChargedPionStat0_5= (TH1D*)fileChargedPionsFullpT->Get("hPionSpectrum_Stat_0_5");
   TH1D* histoChargedPionSys5_10= (TH1D*)fileChargedPionsFullpT->Get("hPionSpectrum_Syst_5_10"); //0-5% 
   TH1D* histoChargedPionStat5_10= (TH1D*)fileChargedPionsFullpT->Get("hPionSpectrum_Stat_5_10");
   TH1D* histoChargedPionSys10_20= (TH1D*)fileChargedPionsFullpT->Get("hPionSpectrum_Syst_10_20"); //0-5% 
   TH1D* histoChargedPionStat10_20= (TH1D*)fileChargedPionsFullpT->Get("hPionSpectrum_Stat_10_20");
   TH1D* histoChargedPionSys20_40= (TH1D*)fileChargedPionsFullpT->Get("hPionSpectrum_Syst_20_40"); //0-5% 
   TH1D* histoChargedPionStat20_40= (TH1D*)fileChargedPionsFullpT->Get("hPionSpectrum_Stat_20_40");
   TH1D* histoChargedPionSys40_60= (TH1D*)fileChargedPionsFullpT->Get("hPionSpectrum_Syst_40_60"); //0-5% 
   TH1D* histoChargedPionStat40_60= (TH1D*)fileChargedPionsFullpT->Get("hPionSpectrum_Stat_40_60");
   TH1D* histoChargedPionSys60_80= (TH1D*)fileChargedPionsFullpT->Get("hPionSpectrum_Syst_60_80"); //0-5% 
   TH1D* histoChargedPionStat60_80= (TH1D*)fileChargedPionsFullpT->Get("hPionSpectrum_Stat_60_80");
   TH1D* histoChargedPionSys80_100= (TH1D*)fileChargedPionsFullpT->Get("hPionSpectrum_Syst_80_100"); //0-5% 
   TH1D* histoChargedPionStat80_100= (TH1D*)fileChargedPionsFullpT->Get("hPionSpectrum_Stat_80_100");

   histoChargedPionSys0_5->Scale(0.5);
   histoChargedPionStat0_5->Scale(0.5);
   histoChargedPionSys5_10->Scale(0.5);
   histoChargedPionStat5_10->Scale(0.5);
   histoChargedPionSys10_20->Scale(0.5);
   histoChargedPionStat10_20->Scale(0.5);
   histoChargedPionSys20_40->Scale(0.5);
   histoChargedPionStat20_40->Scale(0.5);
   histoChargedPionSys40_60->Scale(0.5);
   histoChargedPionStat40_60->Scale(0.5);
   histoChargedPionSys60_80->Scale(0.5);
   histoChargedPionStat60_80->Scale(0.5);
   histoChargedPionSys80_100->Scale(0.5);
   histoChargedPionStat80_100->Scale(0.5);


   TH1D *histoChargedPionSpecSys=(TH1D*)histoChargedPionSys0_5->Clone();
   TH1D *histoChargedPionSpecStat=(TH1D*)histoChargedPionStat0_5->Clone();

   histoChargedPionSpecSys->Add(histoChargedPionSys0_5,histoChargedPionSys5_10,0.05,0.05); //weight with the centrality bin width
   histoChargedPionSpecSys->Add(histoChargedPionSys10_20,0.1);
   histoChargedPionSpecSys->Add(histoChargedPionSys20_40,0.2);
   histoChargedPionSpecSys->Add(histoChargedPionSys40_60,0.2);
   histoChargedPionSpecSys->Add(histoChargedPionSys60_80,0.2);
   histoChargedPionSpecSys->Add(histoChargedPionSys80_100,0.2);

   histoChargedPionSpecStat->Add(histoChargedPionStat0_5,histoChargedPionStat5_10,0.05,0.05); //weight with the centrality bin width
   histoChargedPionSpecStat->Add(histoChargedPionStat10_20,0.1);
   histoChargedPionSpecStat->Add(histoChargedPionStat20_40,0.2);
   histoChargedPionSpecStat->Add(histoChargedPionStat40_60,0.2);
   histoChargedPionSpecStat->Add(histoChargedPionStat60_80,0.2);
   histoChargedPionSpecStat->Add(histoChargedPionStat80_100,0.2);

   TH1D *histoChargedPionSpecSys0020=(TH1D*)histoChargedPionSys0_5->Clone();
   TH1D *histoChargedPionSpecStat0020=(TH1D*)histoChargedPionStat0_5->Clone();
   histoChargedPionSpecSys0020->Add(histoChargedPionSys0_5,histoChargedPionSys5_10,0.25,0.25); //weight with the centrality bin width
   histoChargedPionSpecSys0020->Add(histoChargedPionSys10_20,0.5);
   histoChargedPionSpecStat0020->Add(histoChargedPionStat0_5,histoChargedPionStat5_10,0.25,0.25); //weight with the centrality bin width
   histoChargedPionSpecStat0020->Add(histoChargedPionStat10_20,0.5);

   TH1D *histoChargedPionSpecSys2040=(TH1D*)histoChargedPionSys20_40->Clone();
   TH1D *histoChargedPionSpecStat2040=(TH1D*)histoChargedPionStat20_40->Clone();

   TH1D *histoChargedPionSpecSys4060=(TH1D*)histoChargedPionSys40_60->Clone();
   TH1D *histoChargedPionSpecStat4060=(TH1D*)histoChargedPionStat40_60->Clone();

   TH1D *histoChargedPionSpecSys6080=(TH1D*)histoChargedPionSys60_80->Clone();
   TH1D *histoChargedPionSpecStat6080=(TH1D*)histoChargedPionStat60_80->Clone();
   TH1D *histoChargedPionSpecSys60100=(TH1D*)histoChargedPionSys60_80->Clone();
   TH1D *histoChargedPionSpecStat60100=(TH1D*)histoChargedPionStat60_80->Clone();
   histoChargedPionSpecSys60100->Add(histoChargedPionSys60_80,histoChargedPionSys80_100,0.5,0.5);
   histoChargedPionSpecStat60100->Add(histoChargedPionStat60_80,histoChargedPionStat80_100,0.5,0.5);


   Double_t SysErrorRel[7]={0};
   Double_t RelErrorMB=0;
   Double_t RelError0020=0;
   Double_t RelError60100=0;
   for(Int_t i = 1; i < histoChargedPionSys0_5->GetNbinsX()+1; i++){

     if (histoChargedPionSys0_5->GetBinContent(i) != 0){
         SysErrorRel[0]= histoChargedPionSys0_5->GetBinError(i)/histoChargedPionSys0_5->GetBinContent(i)*100 ;
      } else    SysErrorRel[0] = 0;
     if (histoChargedPionSys5_10->GetBinContent(i) != 0){
         SysErrorRel[1]= histoChargedPionSys5_10->GetBinError(i)/histoChargedPionSys5_10->GetBinContent(i)*100 ;
      } else    SysErrorRel[1] = 0;
     if (histoChargedPionSys10_20->GetBinContent(i) != 0){
         SysErrorRel[2]= histoChargedPionSys10_20->GetBinError(i)/histoChargedPionSys10_20->GetBinContent(i)*100 ;
      } else    SysErrorRel[2] = 0;
     if (histoChargedPionSys20_40->GetBinContent(i) != 0){
         SysErrorRel[3]= histoChargedPionSys20_40->GetBinError(i)/histoChargedPionSys20_40->GetBinContent(i)*100 ;
      } else    SysErrorRel[3] = 0;
     if (histoChargedPionSys40_60->GetBinContent(i) != 0){
         SysErrorRel[4]= histoChargedPionSys40_60->GetBinError(i)/histoChargedPionSys40_60->GetBinContent(i)*100 ;
      } else    SysErrorRel[4] = 0;
     if (histoChargedPionSys60_80->GetBinContent(i) != 0){
         SysErrorRel[5]= histoChargedPionSys60_80->GetBinError(i)/histoChargedPionSys60_80->GetBinContent(i)*100 ;
      } else    SysErrorRel[5] = 0;
     if (histoChargedPionSys80_100->GetBinContent(i) != 0){
         SysErrorRel[6]= histoChargedPionSys80_100->GetBinError(i)/histoChargedPionSys80_100->GetBinContent(i)*100 ;
      } else    SysErrorRel[6] = 0;
     RelErrorMB=SysErrorRel[0];
     for(Int_t j= 1; j < 7; j++){
     
       if(SysErrorRel[j-1]<SysErrorRel[j])      RelErrorMB=SysErrorRel[j];
     }
     RelError0020=SysErrorRel[0];
     for(Int_t j= 1; j < 3; j++){
     
       if(SysErrorRel[j-1]<SysErrorRel[j])      RelError0020=SysErrorRel[j];
     }
     RelError60100=SysErrorRel[5];
         if(SysErrorRel[5]<SysErrorRel[6])      RelError0020=SysErrorRel[6];

	 histoChargedPionSpecSys->SetBinError(i, histoChargedPionSpecSys->GetBinContent(i)*RelErrorMB/100);
	 histoChargedPionSpecSys0020->SetBinError(i, histoChargedPionSpecSys0020->GetBinContent(i)*RelError0020/100);
	 histoChargedPionSpecSys60100->SetBinError(i, histoChargedPionSpecSys60100->GetBinContent(i)*RelError60100/100);

   }






 

   cout << "*************************************************************************"<< endl;  
   cout << "******************************  pPb *************************************"<< endl;
   cout << "*************************************************************************"<< endl;
   
   TGraphAsymmErrors* graphPCMYieldPi0SysErrpPbCopy = (TGraphAsymmErrors*) graphPCMYieldPi0SysErrpPb->Clone("graphPCMYieldPi0SysErrpPbCopy");
   TGraphAsymmErrors* graphDalitzYieldPi0SysErrpPbCopy = (TGraphAsymmErrors*) graphDalitzYieldPi0SysErrpPb->Clone("graphDalitzYieldPi0SysErrpPbCopy");   
   TGraphAsymmErrors* graphDalitzYieldPi0StatErrpPbCopy = (TGraphAsymmErrors*) histoDalitzYieldPi0pPb->Clone("graphDalitzYieldPi0StatErrpPbCopy");   
   cout << "*************************************************************************"<< endl;  
   cout << "******************************  pPb MinBias *****************************"<< endl;
   cout << "*************************************************************************"<< endl;
 TGraphErrors* bla1 = NULL;
 TGraphErrors* bla2 = NULL;
 TGraphErrors* bla3 = NULL;
 TGraphErrors* bla4 = NULL;
   cout << "PCM " << endl;
   TGraphErrors* graphRatioLowPtChargedPionsPCMpPb = CalculateRatioBetweenSpectraWithDifferentBinning(histoPCMYieldPi0pPb, graphPCMYieldPi0SysErrpPbCopy, histoChargedPionSpecStat, histoChargedPionSpecSys,  kTRUE,  kTRUE,&bla1,&bla2,&bla3,&bla4)  ;
   graphRatioLowPtChargedPionsPCMpPb->Print();
   cout << "Dalitz" << endl;
   TGraphErrors* graphRatioLowPtChargedPionsDalitzpPb = CalculateRatioBetweenSpectraWithDifferentBinning(graphDalitzYieldPi0StatErrpPbCopy, graphDalitzYieldPi0SysErrpPbCopy, histoChargedPionSpecStat, histoChargedPionSpecSys,  kTRUE,  kTRUE,&bla1,&bla2,&bla3,&bla4)  ;
   graphRatioLowPtChargedPionsDalitzpPb->Print();
   
   cout << "PHOS " << endl;
   TGraphErrors* graphRatioLowPtChargedPionsPHOSpPb = CalculateRatioBetweenSpectraWithDifferentBinning(histoPHOSYieldPi0pPb, histoPHOSYieldPi0pPbSys, histoChargedPionSpecStat, histoChargedPionSpecSys,  kTRUE,  kTRUE,&bla1,&bla2,&bla3,&bla4)  ;
   graphRatioLowPtChargedPionsPHOSpPb->Print();
   
  cout << "EMCal" << endl;
   TGraphErrors* graphRatioLowPtChargedPionsEMCalpPb = CalculateRatioBetweenSpectraWithDifferentBinning(histoEMCalYieldPi0pPb, histoEMCalYieldPi0pPbSys, histoChargedPionSpecStat, histoChargedPionSpecSys,  kTRUE,  kTRUE,&bla1,&bla2,&bla3,&bla4)  ;
   graphRatioLowPtChargedPionsEMCalpPb->Print();
   

   //************************************************************************************************************
   //******************************  plotting just minBias individual measurements ******************************
   //************************************************************************************************************

   TCanvas* canvasCompYieldpPbInd = new TCanvas("canvasCompYieldpPbInd","",200,10,700,500);  // gives the page size
   DrawGammaCanvasSettings( canvasCompYieldpPbInd,  0.12, 0.02, 0.02, 0.12);
   
   canvasCompYieldpPbInd->SetLogx();
   TH2F * histo2DCompCombinedRatio2;
   histo2DCompCombinedRatio2 = new TH2F("histo2DCompCombinedRatio2","histo2DCompCombinedRatio2",1000,0.23,30.,1000,0.2,4.   );
   SetStyleHistoTH2ForGraphs(histo2DCompCombinedRatio2, "#it{p}_{T} (GeV/#it{c})","#pi^{0}/#pi^{#pm}", 0.05,0.064, 0.05,0.06, 0.8,0.9, 512, 505);
   histo2DCompCombinedRatio2->GetXaxis()->SetRangeUser(0.23,30.);
   histo2DCompCombinedRatio2->GetYaxis()->SetRangeUser(0.1,2.1);
   histo2DCompCombinedRatio2->DrawCopy();

      DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPCMpPb,20,markerSizeComparison, kBlue+2, kBlue+2);
         graphRatioLowPtChargedPionsPCMpPb->DrawClone("E1psame");
      DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsDalitzpPb,25,markerSizeComparison, kBlue, kBlue);
      graphRatioLowPtChargedPionsDalitzpPb->Draw("E1psame");

      DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPHOSpPb,20,markerSizeComparison,  kMagenta+1, kMagenta+1);
         graphRatioLowPtChargedPionsPHOSpPb->Draw("E1psame");
	 DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsEMCalpPb,28,markerSizeComparison,  kGreen+2, kGreen+2);
	 graphRatioLowPtChargedPionsEMCalpPb->Draw("E1psame");
	    graphRatioLowPtChargedPionsPCMpPb->Draw("E1psame");
      TLatex *labelRatioPi0pPb = new TLatex(0.16,0.9,collisionSystempPb.Data());
      SetStyleTLatex( labelRatioPi0pPb, 0.06,4);
      labelRatioPi0pPb->Draw();

      TLegend* legendPi0CompIndChargedPionspPb = new TLegend(0.13,0.15,0.76,0.35);
      legendPi0CompIndChargedPionspPb->SetFillColor(0);
      legendPi0CompIndChargedPionspPb->SetLineColor(0);
      legendPi0CompIndChargedPionspPb->SetNColumns(1);
      //   legendPi0CompIndChargedPionspPb->SetNColumns(2);
      legendPi0CompIndChargedPionspPb->SetTextSize(0.038);
      legendPi0CompIndChargedPionspPb->SetMargin(0.14);
      legendPi0CompIndChargedPionspPb->AddEntry(graphRatioLowPtChargedPionsPCMpPb,Form("#pi^{0}/#pi^{#pm}  (PCM) -0.8 < y_{lab} < 0.8, %s",datePCM.Data()),"p"); 
      legendPi0CompIndChargedPionspPb->AddEntry(graphRatioLowPtChargedPionsDalitzpPb,Form("#pi^{0}/#pi^{#pm}  (Dalitz) -0.8 < y_{lab} < 0.8, %s",dateDalitz.Data()),"p");
      legendPi0CompIndChargedPionspPb->AddEntry(graphRatioLowPtChargedPionsPHOSpPb,Form("#pi^{0}/#pi^{#pm}  (PHOS)  -0.135 < y_{lab} < 0.135, %s",datePHOS.Data()),"p");
      legendPi0CompIndChargedPionspPb->AddEntry(graphRatioLowPtChargedPionsEMCalpPb,Form("#pi^{0}/#pi^{#pm}  (EMCal) -0.7 < y_{lab} < 0.7, %s",dateEMCal.Data()),"p");  

 //  legendPi0CompIndChargedPionspPb->AddEntry(graphRatioLowPtChargedPionsPCMpPb,"#pi^{0}/#pi^{#pm}  (PCM) -0.8 < y_{lab} < 0.8","p"); 
//       legendPi0CompIndChargedPionspPb->AddEntry(graphRatioLowPtChargedPionsDalitzpPb,"#pi^{0}/#pi^{#pm}  (Dalitz) -0.8 < y_{lab} < 0.8","p");
//       legendPi0CompIndChargedPionspPb->AddEntry(graphRatioLowPtChargedPionsPHOSpPb,"#pi^{0}/#pi^{#pm}  (PHOS)  -0.135 < y_{lab} < 0.135","p");
//       legendPi0CompIndChargedPionspPb->AddEntry(graphRatioLowPtChargedPionsEMCalpPb,"#pi^{0}/#pi^{#pm}  (EMCal) -0.7 < y_{lab} < 0.7","p");
//       legendPi0CompIndChargedPionspPb->Draw();

      legendPi0CompIndChargedPionspPb->Draw();
      DrawGammaLines(0., 30. , 1, 1 ,1,kGray);
   
  
   canvasCompYieldpPbInd->Update();
   canvasCompYieldpPbInd->Print(Form("%s/ComparisonChargedToNeutralPionsAll_pPb.%s",outputDir.Data(),suffix.Data()));

//    TCanvas * c1 = new TCanvas("c1","c1",600,800);
//       graphRatioLowPtChargedPionsPCMpPb->DrawClone("E1psame");

//    c1->Update();
}
