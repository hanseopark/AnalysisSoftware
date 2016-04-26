/****************************************************************************************************************************
 ******  provided by Gamma Conversion Group, PWG4,                                           *****
 ******     Ana Marin, marin@physi.uni-heidelberg.de                                      *****
 ******        Kathrin Koch, kkoch@physi.uni-heidelberg.de                                      *****
 ******     Friederike Bock, friederike.bock@cern.ch                                      *****
 *****************************************************************************************************************************/

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
#include "../CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "../CommonHeaders/ConversionFunctions.h"

void MergePileUpCorrectionWithProperWeighting(TString fCutSelection, TString mesonType, TString fSuffix, TString fEnergyFlag, TString fMCFlag){
   
   TString nameFilesDCAAna[5] = {  
      Form("%s/%s/%s_%s_GammaConvV1DCATestAnalysedLHC10b.root",fCutSelection.Data(), fEnergyFlag.Data(), mesonType.Data(), fMCFlag.Data()),
      Form("%s/%s/%s_%s_GammaConvV1DCATestAnalysedLHC10c.root",fCutSelection.Data(), fEnergyFlag.Data(), mesonType.Data(), fMCFlag.Data()),
      Form("%s/%s/%s_%s_GammaConvV1DCATestAnalysedLHC10d.root",fCutSelection.Data(), fEnergyFlag.Data(), mesonType.Data(), fMCFlag.Data()),
      Form("%s/%s/%s_%s_GammaConvV1DCATestAnalysedLHC10e.root",fCutSelection.Data(), fEnergyFlag.Data(), mesonType.Data(), fMCFlag.Data()),
      Form("%s/%s/%s_%s_GammaConvV1DCATestAnalysed.root",fCutSelection.Data(), fEnergyFlag.Data(), mesonType.Data(), fMCFlag.Data())
   };
   
   TString outputDir = Form("%s/%s/%s/AnalyseDCATests",fCutSelection.Data(),fEnergyFlag.Data(),fSuffix.Data());
   gSystem->Exec("mkdir -p "+outputDir);
   
   Double_t numberEvents[5];
   numberEvents[4] = 0;
   
   for (Int_t nFiles = 0; nFiles < 4; nFiles ++){
      TFile* fileInput = new TFile(nameFilesDCAAna[nFiles]);
      TH1F* fEventQuality = (TH1F*)fileInput->Get("NEvents");
      if (fEnergyFlag.CompareTo("PbPb_2.76TeV")==0){
         numberEvents[nFiles] = fEventQuality->GetBinContent(1);
      } else {
         numberEvents[nFiles]=  GetNEvents(fEventQuality);
      }
      numberEvents[4] = numberEvents[4]+numberEvents[nFiles];
      cout << "Period " << nFiles << "\t" <<numberEvents[nFiles] << endl;
   }  
   cout << "full data set " << numberEvents[4] << endl;
   
   Double_t ratioData[4];
   ratioData[0] = numberEvents[0]/numberEvents[4];
   ratioData[1] = numberEvents[1]/numberEvents[4];
   ratioData[2] = numberEvents[2]/numberEvents[4];
   ratioData[3] = numberEvents[3]/numberEvents[4];

   cout<< "Ratios Data:  B/Full =" << ratioData[0] << endl;
   cout<< "Ratios Data:  C/Full =" << ratioData[1] << endl;
   cout<< "Ratios Data:  D/Full =" << ratioData[2] << endl;
   cout<< "Ratios Data:  E/Full =" << ratioData[3] << endl;

   TFile* fileDCAAnalysisData[4];
   TH1D *histoFracCatvsPt[5][6];
   TH1D *histoFracIntHistBGvsPt[5][6];
   TH1D *histoCorrectionFactorsHistCat[5][6];
   TH1D* histoCorrectionFactorsHistvsPt[5];
   TH1D* histoCorrectionFactorsFitvsPt[5];
   TH1D* histoCorrectionFactorsHistvsPtCatA[5];
   TH1D* histoCorrectionFactorsHistvsPtCatC[5];
   TH1D* histoCorrectionFactorsHistvsPtCatD[5];
   TH1D* histoDCAZUnderMesonAllCat_AllPt[5];
   for ( Int_t i = 0; i < 4; i++){
      fileDCAAnalysisData[i] =         new TFile(nameFilesDCAAna[i]);
      if (fileDCAAnalysisData[i]->IsZombie()){
         cout << "ERROR:: file " << nameFilesDCAAna[i].Data() << " does not exist, cant do the weighting" << endl;
         return;
      }   
   
      histoCorrectionFactorsHistvsPt[i] =             (TH1D*)fileDCAAnalysisData[i]->Get("fHistCorrectionFactorsHistAllCat_vsPt");
      histoCorrectionFactorsFitvsPt[i] =             (TH1D*)fileDCAAnalysisData[i]->Get("fHistCorrectionFactorsFitAllCat_vsPt");
      histoCorrectionFactorsHistvsPtCatA[i] =             (TH1D*)fileDCAAnalysisData[i]->Get("fHistCorrectionFactorsHistvsPt_0");
      histoCorrectionFactorsHistvsPtCatC[i] =             (TH1D*)fileDCAAnalysisData[i]->Get("fHistCorrectionFactorsHistvsPt_1");
      histoCorrectionFactorsHistvsPtCatD[i] =             (TH1D*)fileDCAAnalysisData[i]->Get("fHistCorrectionFactorsHistvsPt_2");
      histoDCAZUnderMesonAllCat_AllPt[i] =  (TH1D*)fileDCAAnalysisData[i]->Get("HistDCAZUnderMesonAllCat_AllPt");
      for (Int_t j = 0; j < 6 ; j++){
         histoFracCatvsPt[i][j] = (TH1D*)fileDCAAnalysisData[i]->Get(Form("fHistFracCat_%i_vsPt",j+1));
         histoFracIntHistBGvsPt[i][j] = (TH1D*)fileDCAAnalysisData[i]->Get(Form("fHistFracIntHistBGvsPt_Cat_%i_Variant_1",j+1));
         histoCorrectionFactorsHistCat[i][j] = (TH1D*)fileDCAAnalysisData[i]->Get(Form("fHistCorrectionFactorsHistCat_%i_Variant_1_vsPt",j+1));
      }
      
   }

   Color_t colorMethod[5] = {kBlack, kCyan+2, kRed+2, kGreen+2, kBlue+2};
   Style_t styleMethod[5] = {20,24,21,29,33};
   TString namePeriod[5] = {"LHC10B","LHC10C","LHC10D","LHC10E","Full"};

   TCanvas* canvasCorrFrac = new TCanvas("canvasCorrFrac","",200,10,1350,900);  // gives the page size
   DrawGammaCanvasSettings( canvasCorrFrac, 0.08, 0.02, 0.02, 0.09);

   canvasCorrFrac->cd();
   DrawAutoGammaMesonHistos( histoCorrectionFactorsHistvsPt[0], 
                           "", "p_{T,#pi^{0}} (GeV/c)", "Contamination from Pileup (%)", 
                           kFALSE, 2.,1e-8, kFALSE,
                           kTRUE, 0, 30, 
                           kFALSE, 0., 10.);
   TLegend* legendFractionCat = new TLegend(0.5,0.8,0.95,0.95);
   legendFractionCat->SetTextSize(0.03);         
   legendFractionCat->SetFillColor(0);
   legendFractionCat->SetLineColor(0);

   for ( Int_t i = 0; i < 4; i++){
      DrawGammaSetMarker(histoCorrectionFactorsHistvsPt[i], styleMethod[i], 1.2, colorMethod[i], colorMethod[i]);      
      histoCorrectionFactorsHistvsPt[i]->GetYaxis()->SetTitleOffset(0.8);
      if (i == 0) histoCorrectionFactorsHistvsPt[i]->DrawCopy("p,e1"); 
         else  histoCorrectionFactorsHistvsPt[i]->DrawCopy("same,p,e1"); 
      legendFractionCat->AddEntry(histoCorrectionFactorsHistvsPt[i],Form("%s",namePeriod[i].Data()),"p");
   }
   
   
   legendFractionCat->Draw();
   canvasCorrFrac->Update(); 
   canvasCorrFrac->SaveAs(Form("%s/%s_CorrectionFactorForDifferentPeriods.%s",outputDir.Data(),mesonType.Data(),fSuffix.Data()));

   histoCorrectionFactorsFitvsPt[4] = NULL;
   for ( Int_t i = 0; i < 4; i++){
      if ( i == 0){
         histoCorrectionFactorsHistvsPt[4] = (TH1D*)histoCorrectionFactorsHistvsPt[i]->Clone();
         histoCorrectionFactorsHistvsPt[4]->Scale(ratioData[i])  ;
         if (histoCorrectionFactorsFitvsPt[i]){
            histoCorrectionFactorsFitvsPt[4] = (TH1D*)histoCorrectionFactorsFitvsPt[i]->Clone();
            histoCorrectionFactorsFitvsPt[4]->Scale(ratioData[i])  ;
         }   
         histoCorrectionFactorsHistvsPtCatA[4] =  (TH1D*)histoCorrectionFactorsHistvsPtCatA[i]->Clone();
         histoCorrectionFactorsHistvsPtCatA[4]->Scale(ratioData[i])  ;
         histoCorrectionFactorsHistvsPtCatC[4] =  (TH1D*)histoCorrectionFactorsHistvsPtCatC[i]->Clone()  ;
         histoCorrectionFactorsHistvsPtCatC[4]->Scale(ratioData[i])  ;
         histoCorrectionFactorsHistvsPtCatD[4] =  (TH1D*)histoCorrectionFactorsHistvsPtCatD[i]->Clone()  ;
         histoCorrectionFactorsHistvsPtCatD[4]->Scale(ratioData[i])  ;
         histoDCAZUnderMesonAllCat_AllPt[4] =  (TH1D*)histoDCAZUnderMesonAllCat_AllPt[i]->Clone()  ;
         histoDCAZUnderMesonAllCat_AllPt[4]->Scale(ratioData[i])  ;
         for (Int_t j = 0; j < 6 ; j++){
            histoFracCatvsPt[4][j] =  (TH1D*)histoFracCatvsPt[i][j]->Clone() ;
            histoFracCatvsPt[i][j]->Scale(ratioData[i])  ;
            histoFracIntHistBGvsPt[4][j] =  (TH1D*)histoFracIntHistBGvsPt[i][j]->Clone() ;
            histoFracIntHistBGvsPt[4][j]->Scale(ratioData[i])  ;
            histoCorrectionFactorsHistCat[4][j] =  (TH1D*)histoCorrectionFactorsHistCat[i][j]->Clone();
            histoCorrectionFactorsHistCat[4][j]->Scale(ratioData[i])  ;
         }  
      } else {
         histoCorrectionFactorsHistvsPt[4]->Add(histoCorrectionFactorsHistvsPt[i],ratioData[i] ) ;
         if (histoCorrectionFactorsFitvsPt[i]) histoCorrectionFactorsFitvsPt[4]->Add(histoCorrectionFactorsFitvsPt[i],ratioData[i] ) ;
         histoCorrectionFactorsHistvsPtCatA[4]->Add(histoCorrectionFactorsHistvsPtCatA[i],ratioData[i] ) ;
         histoCorrectionFactorsHistvsPtCatC[4]->Add(histoCorrectionFactorsHistvsPtCatC[i],ratioData[i] ) ;
         histoCorrectionFactorsHistvsPtCatD[4]->Add(histoCorrectionFactorsHistvsPtCatD[i],ratioData[i] ) ;
         histoDCAZUnderMesonAllCat_AllPt[4]->Add(histoDCAZUnderMesonAllCat_AllPt[i],ratioData[i] ) ;
         for (Int_t j = 0; j < 6 ; j++){
            histoFracCatvsPt[4][i]->Add(histoFracCatvsPt[i][j],ratioData[i] ) ;
            histoFracIntHistBGvsPt[4][j]->Add(histoFracIntHistBGvsPt[i][j],ratioData[i] ) ;
            histoCorrectionFactorsHistCat[4][j]->Add(histoCorrectionFactorsHistCat[i][j],ratioData[i] ) ;
         }
      }
   }
   
   canvasCorrFrac->cd();
   DrawAutoGammaMesonHistos( histoCorrectionFactorsHistvsPt[0], 
                           "", "p_{T,#pi^{0}} (GeV/c)", "Contamination from Pileup (%)", 
                           kFALSE, 2.,1e-8, kFALSE,
                           kTRUE, 0, 30, 
                           kFALSE, 0., 10.);
   TLegend* legendFractionCatFinal = new TLegend(0.5,0.8,0.95,0.95);
   legendFractionCatFinal->SetTextSize(0.03);         
   legendFractionCatFinal->SetFillColor(0);
   legendFractionCatFinal->SetLineColor(0);

   for ( Int_t i = 0; i < 5; i++){
      DrawGammaSetMarker(histoCorrectionFactorsHistvsPt[i], styleMethod[i], 1.2, colorMethod[i], colorMethod[i]);      
      histoCorrectionFactorsHistvsPt[i]->GetYaxis()->SetTitleOffset(0.8);
      if (i == 0) histoCorrectionFactorsHistvsPt[i]->DrawCopy("p,e1"); 
         else  histoCorrectionFactorsHistvsPt[i]->DrawCopy("same,p,e1"); 
      legendFractionCatFinal->AddEntry(histoCorrectionFactorsHistvsPt[i],Form("%s",namePeriod[i].Data()),"p");
   }
   
   
   legendFractionCatFinal->Draw();
   canvasCorrFrac->Update(); 
   canvasCorrFrac->SaveAs(Form("%s/%s_CorrectionFactorWeighted.%s",outputDir.Data(),mesonType.Data(),fSuffix.Data()));

   
   cout << nameFilesDCAAna[4].Data() << endl;
   TFile* fileSum = new TFile (nameFilesDCAAna[4],"RECREATE");
         if (histoCorrectionFactorsHistvsPt[4])histoCorrectionFactorsHistvsPt[4]->Write();
         if (histoCorrectionFactorsFitvsPt[4]) histoCorrectionFactorsFitvsPt[4]->Write();
         if (histoCorrectionFactorsHistvsPtCatA[4]) histoCorrectionFactorsHistvsPtCatA[4]->Write();
         if (histoCorrectionFactorsHistvsPtCatC[4])histoCorrectionFactorsHistvsPtCatC[4]->Write();
         if (histoCorrectionFactorsHistvsPtCatD[4])histoCorrectionFactorsHistvsPtCatD[4]->Write();
         if (histoDCAZUnderMesonAllCat_AllPt[4])histoDCAZUnderMesonAllCat_AllPt[4]->Write();
         for (Int_t j = 0; j < 6 ; j++){
            if (histoFracCatvsPt[4][j]) histoFracCatvsPt[4][j]->Write();
            if (histoFracIntHistBGvsPt[4][j]) histoFracIntHistBGvsPt[4][j]->Write();
            if (histoCorrectionFactorsHistCat[4][j]) histoCorrectionFactorsHistCat[4][j]->Write();
         }
   
   fileSum->Write();
   fileSum->Close();
   
   
}
