#include "../CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "../CommonHeaders/PlottingGammaConversionHistos.h"

void SetStyle(TH1D* histo, Int_t i){

    if(histo==NULL) return;
    //cout << "Set style for " << histo->GetName() << endl;
    
    if(i==3){
        histo->SetMarkerStyle(20);
        histo->SetMarkerColor(kBlack);
        histo->SetLineColor(kBlack);
    } else if(i==2){
        histo->SetMarkerStyle(25);
        histo->SetMarkerColor(kBlue);
        histo->SetLineColor(kBlue);
    } else if(i==1){
        histo->SetMarkerStyle(24);
        histo->SetMarkerColor(kRed);
        histo->SetLineColor(kRed);
    } else if(i==4){
        histo->SetMarkerStyle(24);
        histo->SetMarkerColor(kGreen+2);
        histo->SetLineColor(kGreen+2);
    } else if(i==5){
        histo->SetMarkerStyle(2);
        histo->SetMarkerColor(kOrange-3);
        histo->SetLineColor(kOrange-3);
    } else if(i==6){
        histo->SetMarkerStyle(2);
        histo->SetMarkerColor(kGray+2);
        histo->SetLineColor(kGray+2);
    }
    
}

void CompareFitParameters(TString fileData1 = "/home/meike/analysis/results/photonconvResults/PbPb/remainingBGAddChecks/10110a13_00200009247600008250404000_0152501500000000/PbPb_5.02TeV/afterburner_AOD_MCmergedWithAddSigSep_linBG/Pi0_data_GammaConvV1WithoutCorrection_10110a13_00200009247600008250404000_0152501500000000.root",
                          TString fileData2 = "/home/meike/analysis/results/photonconvResults/PbPb/remainingBGAddChecks/10110a13_00200009247600008250404000_0152501500000000/PbPb_5.02TeV/afterburner_AOD_MCmergedWithAddSigSep_pol2BG/Pi0_data_GammaConvV1WithoutCorrection_10110a13_00200009247600008250404000_0152501500000000.root",
                          TString fileMC1 = "/home/meike/analysis/results/photonconvResults/PbPb/remainingBGAddChecks/10110a13_00200009247600008250404000_0152501500000000/PbPb_5.02TeV/afterburner_AOD_MCmergedWithAddSigSep_linBG/Pi0_MC_GammaConvV1WithoutCorrection_10110a13_00200009247600008250404000_0152501500000000.root",
                          TString fileMC2 = "/home/meike/analysis/results/photonconvResults/PbPb/remainingBGAddChecks/10110a13_00200009247600008250404000_0152501500000000/PbPb_5.02TeV/afterburner_AOD_MCmergedWithAddSigSep_pol2BG/Pi0_MC_GammaConvV1WithoutCorrection_10110a13_00200009247600008250404000_0152501500000000.root",
                          TString fileMCTrue1 = "/home/meike/analysis/results/photonconvResults/PbPb/remainingBGAddChecks/10110a13_00200009247600008250404000_0152501500000000/PbPb_5.02TeV/afterburner_AOD_MCmergedWithAddSigSep_linBG/Pi0_MC_GammaConvV1CorrectionHistos_10110a13_00200009247600008250404000_0152501500000000.root",
                          TString fileMCTrue2 = "/home/meike/analysis/results/photonconvResults/PbPb/remainingBGAddChecks/10110a13_00200009247600008250404000_0152501500000000/PbPb_5.02TeV/afterburner_AOD_MCmergedWithAddSigSep_pol2BG/Pi0_MC_GammaConvV1CorrectionHistos_10110a13_00200009247600008250404000_0152501500000000.root",
                          TString meson = "Pi0",
                          TString cut = "10110a13_00200009247600008250404000_0152501500000000",
                          TString label1 = "lin",
                          TString label2 = "pol2",
                          TString suffix = "pdf"){

    //StyleSettingsThesis(suffix);
    SetPlotStyle();

    TString dateForOutput = ReturnDateStringForOutput();
    TString outputDir     = Form("%s/%s",suffix.Data(), dateForOutput.Data());
    gSystem->Exec("mkdir -p "+outputDir); 
    
   TFile* fData1 = new TFile(fileData1.Data());
   TH1D* histoMassMesonData1  = (TH1D*) fData1->Get("histoMassMeson");
   TH1D* histoFWHMMesonData1  = (TH1D*) fData1->Get("histoFWHMMeson");
   TH1D* histoLambdaTailData1 = (TH1D*) fData1->Get("histoLambdaTail");
   TH1D* histoAmplitudeData1  = (TH1D*) fData1->Get("histoAmplitude");
   TH1D* histoSigmaData1      = (TH1D*) fData1->Get("histoSigma");
   
   TFile* fData2 = new TFile(fileData2.Data());
   TH1D* histoMassMesonData2  = (TH1D*) fData2->Get("histoMassMeson");
   TH1D* histoFWHMMesonData2  = (TH1D*) fData2->Get("histoFWHMMeson");
   TH1D* histoLambdaTailData2 = (TH1D*) fData2->Get("histoLambdaTail");
   TH1D* histoAmplitudeData2  = (TH1D*) fData2->Get("histoAmplitude");
   TH1D* histoSigmaData2      = (TH1D*) fData2->Get("histoSigma");
   
   TFile* fMC1     = new TFile(fileMC1.Data());
   TFile* fMCTrue1 = new TFile(fileMCTrue1.Data());
   TH1D* histoMassMesonMC1      = (TH1D*) fMC1->Get("histoMassMeson");
   TH1D* histoFWHMMesonMC1      = (TH1D*) fMC1->Get("histoFWHMMeson");
   TH1D* histoLambdaTailMC1     = (TH1D*) fMC1->Get("histoLambdaTail");
   TH1D* histoAmplitudeMC1      = (TH1D*) fMC1->Get("histoAmplitude");
   TH1D* histoSigmaMC1          = (TH1D*) fMC1->Get("histoSigma");
   
   TH1D* histoTrueMassMesonMC1  = (TH1D*) fMCTrue1->Get("histoTrueMassMeson");
   TH1D* histoTrueFWHMMesonMC1  = (TH1D*) fMCTrue1->Get("histoTrueFWHMMeson");
   TH1D* histoTrueLambdaTailMC1 = (TH1D*) fMCTrue1->Get("histoTrueLambdaTail");
   TH1D* histoTrueAmplitudeMC1  = (TH1D*) fMCTrue1->Get("histoTrueAmplitude");
   TH1D* histoTrueSigmaMC1      = (TH1D*) fMCTrue1->Get("histoTrueSigma");
     
   TFile* fMC2     = new TFile(fileMC2.Data());
   TFile* fMCTrue2 = new TFile(fileMCTrue2.Data());
   TH1D* histoMassMesonMC2      = (TH1D*) fMC2->Get("histoMassMeson");
   TH1D* histoFWHMMesonMC2      = (TH1D*) fMC2->Get("histoFWHMMeson");
   TH1D* histoLambdaTailMC2     = (TH1D*) fMC2->Get("histoLambdaTail");
   TH1D* histoAmplitudeMC2      = (TH1D*) fMC2->Get("histoAmplitude");
   TH1D* histoSigmaMC2          = (TH1D*) fMC2->Get("histoSigma");
   
   TH1D* histoTrueMassMesonMC2  = (TH1D*) fMCTrue2->Get("histoTrueMassMeson");
   TH1D* histoTrueFWHMMesonMC2  = (TH1D*) fMCTrue2->Get("histoTrueFWHMMeson");
   TH1D* histoTrueLambdaTailMC2 = (TH1D*) fMCTrue2->Get("histoTrueLambdaTail");
   TH1D* histoTrueAmplitudeMC2  = (TH1D*) fMCTrue2->Get("histoTrueAmplitude");
   TH1D* histoTrueSigmaMC2      = (TH1D*) fMCTrue2->Get("histoTrueSigma");

   
   TLegend* legend = new TLegend(0.8,0.8,0.9,0.9);

   // data 1
   SetStyle(histoMassMesonData1, 1);
   SetStyle(histoFWHMMesonData1, 1);
   SetStyle(histoLambdaTailData1, 1);
   SetStyle(histoAmplitudeData1, 1);
   SetStyle(histoSigmaData1, 1);
   
   // data 2
   SetStyle(histoMassMesonData2, 2);
   SetStyle(histoFWHMMesonData2, 2);
   SetStyle(histoLambdaTailData2, 2);
   SetStyle(histoAmplitudeData2, 2);
   SetStyle(histoSigmaData2, 2);
   
   // MC 1
   SetStyle(histoMassMesonMC1, 3);
   SetStyle(histoFWHMMesonMC1, 3);
   SetStyle(histoLambdaTailMC1, 3);
   SetStyle(histoAmplitudeMC1, 3);
   SetStyle(histoSigmaMC1, 3);
   
   SetStyle(histoTrueMassMesonMC1, 4);
   SetStyle(histoTrueFWHMMesonMC1, 4);
   SetStyle(histoTrueLambdaTailMC1, 4);
   SetStyle(histoTrueAmplitudeMC1, 4);
   SetStyle(histoTrueSigmaMC1, 4);
   
   // MC 2
   SetStyle(histoMassMesonMC2, 5);
   SetStyle(histoFWHMMesonMC2, 5);
   SetStyle(histoLambdaTailMC2, 5);
   SetStyle(histoAmplitudeMC2, 5);
   SetStyle(histoSigmaMC2, 5);
   
   SetStyle(histoTrueMassMesonMC2, 6);
   SetStyle(histoTrueFWHMMesonMC2, 6);
   SetStyle(histoTrueLambdaTailMC2, 6);
   SetStyle(histoTrueAmplitudeMC2, 6);
   SetStyle(histoTrueSigmaMC2, 6);

   
   TString plotName = "";
   Int_t logx = 1;

   Double_t MinMass   = 0.131;
   Double_t MaxMass   = 0.139;
   Double_t MinWidth  = 0.004;
   Double_t MaxWidth  = 0.018; 
   Double_t MinSigma  = 0.000;
   Double_t MaxSigma  = 0.011;
   Double_t MinLambda = 0.000;
   Double_t MaxLambda = 0.018;
   Double_t MinAmplitude = 10;
   Double_t MaxAmplitude = 1000000;
   
   // Mass Data
   
   plotName = "CompareMassData";
   TCanvas *canvasData = new TCanvas("canvasData","",1550,1200);
   canvasData->SetLogx(logx);
   histoMassMesonData1->Draw("same");
   histoMassMesonData2->Draw("same");
   if(meson.CompareTo("Pi0")==0){
       histoMassMesonData1->GetYaxis()->SetRangeUser(MinMass, MaxMass);
   }
   legend->AddEntry(histoMassMesonData1, Form("Data %s",label1.Data()), "p");
   legend->AddEntry(histoMassMesonData2, Form("Data %s",label2.Data()), "p");
   legend->Draw();
   canvasData->SaveAs(Form("%s/%s_%s_%s.%s", outputDir.Data(), meson.Data(), plotName.Data(), cut.Data(), suffix.Data()));
   
   // Mass rec and true MC
   
   plotName = "CompareMassMC";
   TCanvas *canvasMC = new TCanvas("canvasMC","",1550,1200);
   canvasMC->SetLogx(logx);
   histoMassMesonMC1->Draw("same");
   histoMassMesonMC2->Draw("same");
   histoTrueMassMesonMC1->Draw("same");
   histoTrueMassMesonMC2->Draw("same");
   if(meson.CompareTo("Pi0")==0){
       histoMassMesonData1->GetYaxis()->SetRangeUser(MinMass, MaxMass);
   }
   legend->Clear();
   legend->AddEntry(histoMassMesonMC1, Form("rec MC %s",label1.Data()), "p");
   legend->AddEntry(histoMassMesonMC2, Form("rec MC %s",label2.Data()), "p");
   legend->AddEntry(histoTrueMassMesonMC1, Form("true MC %s",label1.Data()), "p");
   legend->AddEntry(histoTrueMassMesonMC2, Form("true MC %s",label2.Data()), "p");
   legend->Draw();
   canvasMC->SaveAs(Form("%s/%s_%s_%s.%s", outputDir.Data(), meson.Data(), plotName.Data(), cut.Data(), suffix.Data()));
      
   // Mass Data and rec MC

   plotName = "CompareMassDataAndRecMC";
   TCanvas *canvasDataAndRecMC = new TCanvas("canvasDataAndRecMC","",1550,1200);
   canvasDataAndRecMC->SetLogx(logx);
   histoMassMesonData1->Draw("same");
   histoMassMesonData2->Draw("same");
   histoMassMesonMC1->Draw("same");
   histoMassMesonMC2->Draw("same");
   if(meson.CompareTo("Pi0")==0){
       histoMassMesonData1->GetYaxis()->SetRangeUser(MinMass, MaxMass);
   }
   legend->Clear();
   legend->AddEntry(histoMassMesonData1, Form("Data %s",label1.Data()), "p");
   legend->AddEntry(histoMassMesonData2, Form("Data %s",label2.Data()), "p");
   legend->AddEntry(histoMassMesonMC1, Form("rec MC %s",label1.Data()), "p");
   legend->AddEntry(histoMassMesonMC2, Form("rec MC %s",label2.Data()), "p");
   legend->Draw();
   canvasDataAndRecMC->SaveAs(Form("%s/%s_%s_%s.%s", outputDir.Data(), meson.Data(), plotName.Data(), cut.Data(), suffix.Data()));

   // Mass Data and true MC

   plotName = "CompareMassDataAndTrueMC";
   TCanvas *canvasDataAndTrueMC = new TCanvas("canvasDataAndTrueMC","",1550,1200);
   canvasDataAndTrueMC->SetLogx(logx);
   histoMassMesonData1->Draw("same");
   histoMassMesonData2->Draw("same");
   histoTrueMassMesonMC1->Draw("same");
   histoTrueMassMesonMC2->Draw("same");
   if(meson.CompareTo("Pi0")==0){
       histoMassMesonData1->GetYaxis()->SetRangeUser(MinMass, MaxMass);
   }
   legend->Clear();
   legend->AddEntry(histoMassMesonData1, Form("Data %s",label1.Data()), "p");
   legend->AddEntry(histoMassMesonData2, Form("Data %s",label2.Data()), "p");
   legend->AddEntry(histoTrueMassMesonMC1, Form("true MC %s",label1.Data()), "p");
   legend->AddEntry(histoTrueMassMesonMC2, Form("true MC %s",label2.Data()), "p");
   legend->Draw();
   canvasDataAndTrueMC->SaveAs(Form("%s/%s_%s_%s.%s", outputDir.Data(), meson.Data(), plotName.Data(), cut.Data(), suffix.Data()));

   delete canvasData, canvasMC, canvasDataAndRecMC, canvasDataAndTrueMC;

  // FWHM Data
   
   plotName = "CompareFWHMData";
   canvasData = new TCanvas("canvasData","",1550,1200);
   canvasData->SetLogx(logx);
   histoFWHMMesonData1->Draw("same");
   histoFWHMMesonData2->Draw("same");
   if(meson.CompareTo("Pi0")==0){
       histoFWHMMesonData1->GetYaxis()->SetRangeUser(MinWidth, MaxWidth);
   }
   legend->Clear();
   legend->AddEntry(histoFWHMMesonData1, Form("Data %s",label1.Data()), "p");
   legend->AddEntry(histoFWHMMesonData2, Form("Data %s",label2.Data()), "p");
   legend->Draw();
   canvasData->SaveAs(Form("%s/%s_%s_%s.%s", outputDir.Data(), meson.Data(), plotName.Data(), cut.Data(), suffix.Data()));
   
   // FWHM rec and true MC
   
   plotName = "CompareFWHMMC";
   canvasMC = new TCanvas("canvasMC","",1550,1200);
   canvasMC->SetLogx(logx);
   histoFWHMMesonMC1->Draw("same");
   histoFWHMMesonMC2->Draw("same");
   histoTrueFWHMMesonMC1->Draw("same");
   histoTrueFWHMMesonMC2->Draw("same");
   if(meson.CompareTo("Pi0")==0){
       histoFWHMMesonData1->GetYaxis()->SetRangeUser(MinWidth, MaxWidth);
   }
   legend->Clear();
   legend->AddEntry(histoFWHMMesonMC1, Form("rec MC %s",label1.Data()), "p");
   legend->AddEntry(histoFWHMMesonMC2, Form("rec MC %s",label2.Data()), "p");
   legend->AddEntry(histoTrueFWHMMesonMC1, Form("true MC %s",label1.Data()), "p");
   legend->AddEntry(histoTrueFWHMMesonMC2, Form("true MC %s",label2.Data()), "p");
   legend->Draw();
   canvasMC->SaveAs(Form("%s/%s_%s_%s.%s", outputDir.Data(), meson.Data(), plotName.Data(), cut.Data(), suffix.Data()));
      
   // FWHM Data and rec MC

   plotName = "CompareFWHMDataAndRecMC";
   canvasDataAndRecMC = new TCanvas("canvasDataAndRecMC","",1550,1200);
   canvasDataAndRecMC->SetLogx(logx);
   histoFWHMMesonData1->Draw("same");
   histoFWHMMesonData2->Draw("same");
   histoFWHMMesonMC1->Draw("same");
   histoFWHMMesonMC2->Draw("same");
   if(meson.CompareTo("Pi0")==0){
       histoFWHMMesonData1->GetYaxis()->SetRangeUser(MinWidth, MaxWidth);
   }
   legend->Clear();
   legend->AddEntry(histoFWHMMesonData1, Form("Data %s",label1.Data()), "p");
   legend->AddEntry(histoFWHMMesonData2, Form("Data %s",label2.Data()), "p");
   legend->AddEntry(histoFWHMMesonMC1, Form("rec MC %s",label1.Data()), "p");
   legend->AddEntry(histoFWHMMesonMC2, Form("rec MC %s",label2.Data()), "p");
   legend->Draw();
   canvasDataAndRecMC->SaveAs(Form("%s/%s_%s_%s.%s", outputDir.Data(), meson.Data(), plotName.Data(), cut.Data(), suffix.Data()));

   // FWHM Data and true MC

   plotName = "CompareFWHMDataAndTrueMC";
   canvasDataAndTrueMC = new TCanvas("canvasDataAndTrueMC","",1550,1200);
   canvasDataAndTrueMC->SetLogx(logx);
   histoFWHMMesonData1->Draw("same");
   histoFWHMMesonData2->Draw("same");
   histoTrueFWHMMesonMC1->Draw("same");
   histoTrueFWHMMesonMC2->Draw("same");
   if(meson.CompareTo("Pi0")==0){
       histoFWHMMesonData1->GetYaxis()->SetRangeUser(MinWidth, MaxWidth);
   }
   legend->Clear();
   legend->AddEntry(histoFWHMMesonData1, Form("Data %s",label1.Data()), "p");
   legend->AddEntry(histoFWHMMesonData2, Form("Data %s",label2.Data()), "p");
   legend->AddEntry(histoTrueFWHMMesonMC1, Form("true MC %s",label1.Data()), "p");
   legend->AddEntry(histoTrueFWHMMesonMC2, Form("true MC %s",label2.Data()), "p");
   legend->Draw();
   canvasDataAndTrueMC->SaveAs(Form("%s/%s_%s_%s.%s", outputDir.Data(), meson.Data(), plotName.Data(), cut.Data(), suffix.Data()));

   delete canvasData, canvasMC, canvasDataAndRecMC, canvasDataAndTrueMC;


  // Lambda Tail Data
   
   plotName = "CompareLambdaTailData";
   canvasData = new TCanvas("canvasData","",1550,1200);
   canvasData->SetLogx(logx);
   histoLambdaTailData1->Draw("same");
   histoLambdaTailData2->Draw("same");
   if(meson.CompareTo("Pi0")==0){
       histoLambdaTailData1->GetYaxis()->SetRangeUser(MinLambda, MaxLambda);
   }
   legend->Clear();
   legend->AddEntry(histoLambdaTailData1, Form("Data %s",label1.Data()), "p");
   legend->AddEntry(histoLambdaTailData2, Form("Data %s",label2.Data()), "p");
   legend->Draw();
   canvasData->SaveAs(Form("%s/%s_%s_%s.%s", outputDir.Data(), meson.Data(), plotName.Data(), cut.Data(), suffix.Data()));
   
   // Lambda Tail rec and true MC
   
   plotName = "CompareLambdaTailMC";
   canvasMC = new TCanvas("canvasMC","",1550,1200);
   canvasMC->SetLogx(logx);
   histoLambdaTailMC1->Draw("same");
   histoLambdaTailMC2->Draw("same");
   histoTrueLambdaTailMC1->Draw("same");
   histoTrueLambdaTailMC2->Draw("same");
   if(meson.CompareTo("Pi0")==0){
       histoLambdaTailData1->GetYaxis()->SetRangeUser(MinLambda, MaxLambda);
   }
   legend->Clear();
   legend->AddEntry(histoLambdaTailMC1, Form("rec MC %s",label1.Data()), "p");
   legend->AddEntry(histoLambdaTailMC2, Form("rec MC %s",label2.Data()), "p");
   legend->AddEntry(histoTrueLambdaTailMC1, Form("true MC %s",label1.Data()), "p");
   legend->AddEntry(histoTrueLambdaTailMC2, Form("true MC %s",label2.Data()), "p");
   legend->Draw();
   canvasMC->SaveAs(Form("%s/%s_%s_%s.%s", outputDir.Data(), meson.Data(), plotName.Data(), cut.Data(), suffix.Data()));
      
   // Lambda Tail Data and rec MC

   plotName = "CompareLambdaTailDataAndRecMC";
   canvasDataAndRecMC = new TCanvas("canvasDataAndRecMC","",1550,1200);
   canvasDataAndRecMC->SetLogx(logx);
   histoLambdaTailData1->Draw("same");
   histoLambdaTailData2->Draw("same");
   histoLambdaTailMC1->Draw("same");
   histoLambdaTailMC2->Draw("same");
   if(meson.CompareTo("Pi0")==0){
       histoLambdaTailData1->GetYaxis()->SetRangeUser(MinLambda, MaxLambda);
   }
   legend->Clear();
   legend->AddEntry(histoLambdaTailData1, Form("Data %s",label1.Data()), "p");
   legend->AddEntry(histoLambdaTailData2, Form("Data %s",label2.Data()), "p");
   legend->AddEntry(histoLambdaTailMC1, Form("rec MC %s",label1.Data()), "p");
   legend->AddEntry(histoLambdaTailMC2, Form("rec MC %s",label2.Data()), "p");
   legend->Draw();
   canvasDataAndRecMC->SaveAs(Form("%s/%s_%s_%s.%s", outputDir.Data(), meson.Data(), plotName.Data(), cut.Data(), suffix.Data()));

   // Lambda Tail Data and true MC

   plotName = "CompareLambdaTailDataAndTrueMC";
   canvasDataAndTrueMC = new TCanvas("canvasDataAndTrueMC","",1550,1200);
   canvasDataAndTrueMC->SetLogx(logx);
   histoLambdaTailData1->Draw("same");
   histoLambdaTailData2->Draw("same");
   histoTrueLambdaTailMC1->Draw("same");
   histoTrueLambdaTailMC2->Draw("same");
   if(meson.CompareTo("Pi0")==0){
       histoLambdaTailData1->GetYaxis()->SetRangeUser(MinLambda, MaxLambda);
   }
   legend->Clear();
   legend->AddEntry(histoLambdaTailData1, Form("Data %s",label1.Data()), "p");
   legend->AddEntry(histoLambdaTailData2, Form("Data %s",label2.Data()), "p");
   legend->AddEntry(histoTrueLambdaTailMC1, Form("true MC %s",label1.Data()), "p");
   legend->AddEntry(histoTrueLambdaTailMC2, Form("true MC %s",label2.Data()), "p");
   legend->Draw();
   canvasDataAndTrueMC->SaveAs(Form("%s/%s_%s_%s.%s", outputDir.Data(), meson.Data(), plotName.Data(), cut.Data(), suffix.Data()));

   delete canvasData, canvasMC, canvasDataAndRecMC, canvasDataAndTrueMC;

   
  // Sigma Data
   
   plotName = "CompareSigmaData";
   canvasData = new TCanvas("canvasData","",1550,1200);
   canvasData->SetLogx(logx);
   histoSigmaData1->Draw("same");
   histoSigmaData2->Draw("same");
   if(meson.CompareTo("Pi0")==0){
       histoSigmaData1->GetYaxis()->SetRangeUser(MinSigma, MaxSigma);
   }
   legend->Clear();
   legend->AddEntry(histoSigmaData1, Form("Data %s",label1.Data()), "p");
   legend->AddEntry(histoSigmaData2, Form("Data %s",label2.Data()), "p");
   legend->Draw();
   canvasData->SaveAs(Form("%s/%s_%s_%s.%s", outputDir.Data(), meson.Data(), plotName.Data(), cut.Data(), suffix.Data()));
   
   // Sigma rec and true MC
   
   plotName = "CompareSigmaMC";
   canvasMC = new TCanvas("canvasMC","",1550,1200);
   canvasMC->SetLogx(logx);
   histoSigmaMC1->Draw("same");
   histoSigmaMC2->Draw("same");
   histoTrueSigmaMC1->Draw("same");
   histoTrueSigmaMC2->Draw("same");
   if(meson.CompareTo("Pi0")==0){
       histoSigmaData1->GetYaxis()->SetRangeUser(MinSigma, MaxSigma);
   }
   legend->Clear();
   legend->AddEntry(histoSigmaMC1, Form("rec MC %s",label1.Data()), "p");
   legend->AddEntry(histoSigmaMC2, Form("rec MC %s",label2.Data()), "p");
   legend->AddEntry(histoTrueSigmaMC1, Form("true MC %s",label1.Data()), "p");
   legend->AddEntry(histoTrueSigmaMC2, Form("true MC %s",label2.Data()), "p");
   legend->Draw();
   canvasMC->SaveAs(Form("%s/%s_%s_%s.%s", outputDir.Data(), meson.Data(), plotName.Data(), cut.Data(), suffix.Data()));
      
   // Sigma Data and rec MC

   plotName = "CompareSigmaDataAndRecMC";
   canvasDataAndRecMC = new TCanvas("canvasDataAndRecMC","",1550,1200);
   canvasDataAndRecMC->SetLogx(logx);
   histoSigmaData1->Draw("same");
   histoSigmaData2->Draw("same");
   histoSigmaMC1->Draw("same");
   histoSigmaMC2->Draw("same");
   if(meson.CompareTo("Pi0")==0){
       histoSigmaData1->GetYaxis()->SetRangeUser(MinSigma, MaxSigma);
   }
   legend->Clear();
   legend->AddEntry(histoSigmaData1, Form("Data %s",label1.Data()), "p");
   legend->AddEntry(histoSigmaData2, Form("Data %s",label2.Data()), "p");
   legend->AddEntry(histoSigmaMC1, Form("rec MC %s",label1.Data()), "p");
   legend->AddEntry(histoSigmaMC2, Form("rec MC %s",label2.Data()), "p");
   legend->Draw();
   canvasDataAndRecMC->SaveAs(Form("%s/%s_%s_%s.%s", outputDir.Data(), meson.Data(), plotName.Data(), cut.Data(), suffix.Data()));

   // Sigma Data and true MC

   plotName = "CompareSigmaDataAndTrueMC";
   canvasDataAndTrueMC = new TCanvas("canvasDataAndTrueMC","",1550,1200);
   canvasDataAndTrueMC->SetLogx(logx);
   histoSigmaData1->Draw("same");
   histoSigmaData2->Draw("same");
   histoTrueSigmaMC1->Draw("same");
   histoTrueSigmaMC2->Draw("same");
   if(meson.CompareTo("Pi0")==0){
       histoSigmaData1->GetYaxis()->SetRangeUser(MinSigma, MaxSigma);
   }
   legend->Clear();
   legend->AddEntry(histoSigmaData1, Form("Data %s",label1.Data()), "p");
   legend->AddEntry(histoSigmaData2, Form("Data %s",label2.Data()), "p");
   legend->AddEntry(histoTrueSigmaMC1, Form("true MC %s",label1.Data()), "p");
   legend->AddEntry(histoTrueSigmaMC2, Form("true MC %s",label2.Data()), "p");
   legend->Draw();
   canvasDataAndTrueMC->SaveAs(Form("%s/%s_%s_%s.%s", outputDir.Data(), meson.Data(), plotName.Data(), cut.Data(), suffix.Data()));

   delete canvasData, canvasMC, canvasDataAndRecMC, canvasDataAndTrueMC;
      
  // Amplitude Data
   
   plotName = "CompareAmplitudeData";
   canvasData = new TCanvas("canvasData","",1550,1200);
   canvasData->SetLogx(logx);
   canvasData->SetLogy();
   histoAmplitudeData1->Draw("same");
   histoAmplitudeData2->Draw("same");
   if(meson.CompareTo("Pi0")==0){
       histoAmplitudeData1->GetYaxis()->SetRangeUser(MinAmplitude, MaxAmplitude);
   }
   legend->Clear();
   legend->AddEntry(histoAmplitudeData1, Form("Data %s",label1.Data()), "p");
   legend->AddEntry(histoAmplitudeData2, Form("Data %s",label2.Data()), "p");
   legend->Draw();
   canvasData->SaveAs(Form("%s/%s_%s_%s.%s", outputDir.Data(), meson.Data(), plotName.Data(), cut.Data(), suffix.Data()));
   
   // Amplitude rec and true MC
   
   plotName = "CompareAmplitudeMC";
   canvasMC = new TCanvas("canvasMC","",1550,1200);
   canvasMC->SetLogx(logx);
   canvasMC->SetLogy();
   histoAmplitudeMC1->Draw("same");
   histoAmplitudeMC2->Draw("same");
   histoTrueAmplitudeMC1->Draw("same");
   histoTrueAmplitudeMC2->Draw("same");
   if(meson.CompareTo("Pi0")==0){
       histoAmplitudeData1->GetYaxis()->SetRangeUser(MinAmplitude, MaxAmplitude);
   }
   legend->Clear();
   legend->AddEntry(histoAmplitudeMC1, Form("rec MC %s",label1.Data()), "p");
   legend->AddEntry(histoAmplitudeMC2, Form("rec MC %s",label2.Data()), "p");
   legend->AddEntry(histoTrueAmplitudeMC1, Form("true MC %s",label1.Data()), "p");
   legend->AddEntry(histoTrueAmplitudeMC2, Form("true MC %s",label2.Data()), "p");
   legend->Draw();
   canvasMC->SaveAs(Form("%s/%s_%s_%s.%s", outputDir.Data(), meson.Data(), plotName.Data(), cut.Data(), suffix.Data()));
      
   // Amplitude Data and rec MC

   plotName = "CompareAmplitudeDataAndRecMC";
   canvasDataAndRecMC = new TCanvas("canvasDataAndRecMC","",1550,1200);
   canvasDataAndRecMC->SetLogx(logx);
   canvasDataAndRecMC->SetLogy();
   histoAmplitudeData1->Draw("same");
   histoAmplitudeData2->Draw("same");
   histoAmplitudeMC1->Draw("same");
   histoAmplitudeMC2->Draw("same");
   if(meson.CompareTo("Pi0")==0){
       histoAmplitudeData1->GetYaxis()->SetRangeUser(MinAmplitude, MaxAmplitude);
   }
   legend->Clear();
   legend->AddEntry(histoAmplitudeData1, Form("Data %s",label1.Data()), "p");
   legend->AddEntry(histoAmplitudeData2, Form("Data %s",label2.Data()), "p");
   legend->AddEntry(histoAmplitudeMC1, Form("rec MC %s",label1.Data()), "p");
   legend->AddEntry(histoAmplitudeMC2, Form("rec MC %s",label2.Data()), "p");
   legend->Draw();
   canvasDataAndRecMC->SaveAs(Form("%s/%s_%s_%s.%s", outputDir.Data(), meson.Data(), plotName.Data(), cut.Data(), suffix.Data()));

   // Amplitude Data and true MC

   plotName = "CompareAmplitudeDataAndTrueMC";
   canvasDataAndTrueMC = new TCanvas("canvasDataAndTrueMC","",1550,1200);
   canvasDataAndTrueMC->SetLogx(logx);
   canvasDataAndTrueMC->SetLogy();
   histoAmplitudeData1->Draw("same");
   histoAmplitudeData2->Draw("same");
   histoTrueAmplitudeMC1->Draw("same");
   histoTrueAmplitudeMC2->Draw("same");
   if(meson.CompareTo("Pi0")==0){
       histoAmplitudeData1->GetYaxis()->SetRangeUser(MinAmplitude, MaxAmplitude);
   }
   legend->Clear();
   legend->AddEntry(histoAmplitudeData1, Form("Data %s",label1.Data()), "p");
   legend->AddEntry(histoAmplitudeData2, Form("Data %s",label2.Data()), "p");
   legend->AddEntry(histoTrueAmplitudeMC1, Form("true MC %s",label1.Data()), "p");
   legend->AddEntry(histoTrueAmplitudeMC2, Form("true MC %s",label2.Data()), "p");
   legend->Draw();
   canvasDataAndTrueMC->SaveAs(Form("%s/%s_%s_%s.%s", outputDir.Data(), meson.Data(), plotName.Data(), cut.Data(), suffix.Data()));

   delete canvasData, canvasMC, canvasDataAndRecMC, canvasDataAndTrueMC;
   // end new

   
   delete legend;
   delete fData1, fData2, fMC1, fMC2;

}



