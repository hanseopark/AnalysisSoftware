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
#include <string>
#include "TGaxis.h"
#include "TFractionFitter.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TF1.h"
#include "THStack.h"
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
#include "TEllipse.h"
#include "DirectPhotonFlowFunctions.h"

void  projection_LTM(   TString CentralityLow = "0",
                        TString CentralityHigh = "20",
                        TString Cutnumber = "50200013_00200009007000008250400000",
                        TString Cutnumbereff = "50200013_00200009300000008250400000_0152304500900000"
                     ){
  
  
    //==================================================================================
    //DEFINING SOME STUFF
    //==================================================================================
    
    StyleSettingsThesis();  
    SetPlotStyle();
    
    Int_t Trainconfig = 59;
    
    TString InfoSystem = Form("%s-%s %% PbPb, #sqrt{s}=2.76TeV",CentralityLow.Data(),CentralityHigh.Data());
    TString InfoData = "Data LHC10h";
    TString InfoMC = "MC LHC13d2";

    TFile* fileMC = new TFile("/home/mike/3_PbPb_dirg/1_data/161103_PbPb_v2_MC/GammaConvFlow_59.root");
    TList* listMC_1   = new TList();
    listMC_1     = (TList*)fileMC->Get(Form("GammaConvV1_%i_v2",Trainconfig));
    TList* listMC_2   = new TList();
    listMC_2     = (TList*)listMC_1->FindObject(Form("Cut Number %s",Cutnumber.Data()));
    TList* listMC_3   = new TList();
    listMC_3     = (TList*)listMC_2->FindObject(Form("%s ESD histograms",Cutnumber.Data()));
    
    TFile* fileData  = new TFile("/home/mike/3_PbPb_dirg/1_data/161103_PbPb_v2_data/GammaConvFlow_59.root");
    TList* listData_1   = new TList();
    listData_1     = (TList*)fileData->Get(Form("GammaConvV1_%i_v2",Trainconfig));
    TList* listData_2   = new TList();
    listData_2     = (TList*)listData_1->FindObject(Form("Cut Number %s",Cutnumber.Data()));
    TList* listData_3   = new TList();
    listData_3     = (TList*)listData_2->FindObject(Form("%s ESD histograms",Cutnumber.Data()));
    
    TH2F* h2d_LTM_data = (TH2F*)listData_3->FindObject("LTM_Pt");
    h2d_LTM_data->SetName("LTM_Pt_data");
    TH2F* h2d_LTM_MC = (TH2F*)listMC_3->FindObject("LTM_Pt");
    h2d_LTM_MC->SetName("LTM_Pt_MC");
    TH2F* h2d_LTM_MCgen = (TH2F*)listMC_3->FindObject("LTM_Pt_MCgen");
    h2d_LTM_MCgen->SetName("LTM_Pt_MCgen");
    
    TFile* filepurMC = new TFile("/home/mike/3_PbPb_dirg/1_data/161214_PbPb_v2_MC/GammaConvV1_312_RGamma.root");
    TList* listpurMC_1   = new TList();
    listpurMC_1     = (TList*)filepurMC->Get("GammaConvV1");
    cout << "listpurMC_1 " << listpurMC_1 << endl;
    TList* listpurMC_2   = new TList();
    listpurMC_2     = (TList*)listpurMC_1->FindObject(Form("Cut Number %s",Cutnumbereff.Data()));
    cout << "listpurMC_2 " << listpurMC_2 << endl;
    TList* listpurMC_3   = new TList();
    listpurMC_3     = (TList*)listpurMC_2->FindObject(Form("%s MC histograms",Cutnumbereff.Data()));
    cout << "listpurMC_3 " << listpurMC_3 << endl;
    TList* listpurMC_4   = new TList();
    listpurMC_4     = (TList*)listpurMC_2->FindObject(Form("%s True histograms",Cutnumbereff.Data()));
    cout << "listpurMC_4 " << listpurMC_4 << endl;
    
    TH1F* MC_ConvGamma_efficiency = (TH1F*)listpurMC_3->FindObject("MC_ConvGamma_Pt");
    TH1F* ESD_ConvGamma_efficiency = (TH1F*)listpurMC_4->FindObject("ESD_TrueConvGamma_Pt");
    ESD_ConvGamma_efficiency->Divide(ESD_ConvGamma_efficiency,MC_ConvGamma_efficiency,1,1,"B");
    ESD_ConvGamma_efficiency->SetTitle("");
    ESD_ConvGamma_efficiency->SetMarkerStyle(20);
    ESD_ConvGamma_efficiency->SetMarkerSize(1.2);
    ESD_ConvGamma_efficiency->SetMarkerColor(kGreen+2);
    ESD_ConvGamma_efficiency->SetLineColor(kGreen+2);
    ESD_ConvGamma_efficiency->GetYaxis()->SetTitle("#epsilon");
    
    TH1F*  histoInclusiveEmpty = (TH1F*)Getv2HistStyle();
    histoInclusiveEmpty->GetYaxis()->SetRangeUser(0,0.49);
    
    TH1F*  histoInclusiveEmptyRatio = (TH1F*)Getv2HistStyle();
    histoInclusiveEmptyRatio->GetYaxis()->SetRangeUser(0.91,1.99);
    histoInclusiveEmptyRatio->GetYaxis()->SetTitle("Ratio");
    
    TFile* fileIncl  = new TFile(Form("/home/mike/Dropbox/msas_ALICE/PbPb_v2dirg/PCM_InclusivePhotonFlow_Syst_%s%s.root",CentralityLow.Data(),CentralityHigh.Data()));
    TH1F*  histoInclusive = (TH1F*)fileIncl->Get(Form("histoInclusive_%s%s",CentralityLow.Data(),CentralityHigh.Data()));
    
    TH1D* histoInclusive_corrected_LTM = (TH1D*)histoInclusive->Clone();
    histoInclusive_corrected_LTM->SetMarkerStyle(20);
    histoInclusive_corrected_LTM->SetMarkerSize(1.2);
    histoInclusive_corrected_LTM->SetMarkerColor(kGreen+2);
    histoInclusive_corrected_LTM->SetLineColor(kGreen+2);
    
    TH1F* histoInclusive_corrected_LTM_ratio = (TH1F*)histoInclusive->Clone();
    
    TCanvas* c2 = new TCanvas("c2","",800,400);
    c2->Divide(2,1,0.001,0.001);
    c2->cd(1);
    h2d_LTM_data->Draw();
    c2->cd(2);
    h2d_LTM_MC->Draw();
    
    Double_t newBinsComb[19] = {0.0, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 3.0, 3.3, 3.7, 4.1, 4.6, 5.4, 6.0, 7.0};
    TString newBinsName[19] = {"00", "09", "11", "13", "15", "17", "19", "21", "23", "25", "27", "30", "33", "37", "41", "46", "54", "62", "70"};
    Double_t centralbinvalues[19];
    Double_t binwidths[19];
    Double_t effvariations[19];
    
    for(int i=0;i<18;i++){
      centralbinvalues[i]=newBinsComb[i]+(newBinsComb[i+1]-newBinsComb[i])/2;
      binwidths[i]=(newBinsComb[i+1]-newBinsComb[i])/2;
    }
    
    TH1D* hEff_pt = new TH1D("hEff_pt","hEff_pt",18,newBinsComb);
    hEff_pt->SetTitle("");
    hEff_pt->GetYaxis()->SetTitle("#epsilon  ");
    hEff_pt->GetXaxis()->SetTitle("#it{p}_{T}(GeV/c)");
    hEff_pt->GetXaxis()->SetTitleSize(0.06);
    hEff_pt->GetYaxis()->SetTitleSize(0.06);
    hEff_pt->GetXaxis()->SetTitleOffset(0.8);
    hEff_pt->GetYaxis()->SetTitleOffset(0.85);
    hEff_pt->GetYaxis()->SetRangeUser(0.0,1.0);
    hEff_pt->SetMarkerStyle(20);
    hEff_pt->SetMarkerSize(1.2);
    hEff_pt->SetMarkerColor(kRed+2);
    hEff_pt->SetLineColor(kRed+2);
    
    
    TH1D* hEff_pt_1 = (TH1D*)hEff_pt->Clone();
    hEff_pt_1->SetMarkerStyle(21);
    hEff_pt_1->SetMarkerSize(1.2);
    hEff_pt_1->SetMarkerColor(kGreen+2);
    hEff_pt_1->SetLineColor(kGreen+2);
    TH1D* hEff_pt_2 = (TH1D*)hEff_pt->Clone();
    hEff_pt_2->SetMarkerStyle(33);
    hEff_pt_2->SetMarkerSize(1.2);
    hEff_pt_2->SetMarkerColor(kBlue+2);
    hEff_pt_2->SetLineColor(kBlue+2);
    TH1D* hEff_pt_3 = (TH1D*)hEff_pt->Clone();
    hEff_pt_3->SetMarkerStyle(34);
    hEff_pt_3->SetMarkerSize(1.2);
    hEff_pt_3->SetMarkerColor(kMagenta+2);
    hEff_pt_3->SetLineColor(kMagenta+2);
    
    TH1D* hLTM_Eff_pt = new TH1D("hLTM_Eff_pt","hLTM_Eff_pt",18,newBinsComb);
    hLTM_Eff_pt->SetTitle("");
    hLTM_Eff_pt->GetYaxis()->SetTitle("#epsilon ");
    hLTM_Eff_pt->GetXaxis()->SetTitle("#it{p}_{T}(GeV/c)");
    hLTM_Eff_pt->GetXaxis()->SetTitleSize(0.06);
    hLTM_Eff_pt->GetYaxis()->SetTitleSize(0.06);
    hLTM_Eff_pt->GetXaxis()->SetTitleOffset(0.8);
    hLTM_Eff_pt->GetYaxis()->SetTitleOffset(0.85);
    hLTM_Eff_pt->GetYaxis()->SetRangeUser(0.0,0.49);
    hLTM_Eff_pt->SetMarkerStyle(20);
    hLTM_Eff_pt->SetMarkerSize(1.2);
    hLTM_Eff_pt->SetMarkerColor(kRed+2);
    hLTM_Eff_pt->SetLineColor(kRed+2);
    
    TH1D* hLTM_w2_pt = new TH1D("hLTM_w2_pt","hLTM_w2_pt",18,newBinsComb);
    hLTM_w2_pt->SetTitle("");
    hLTM_w2_pt->GetYaxis()->SetTitle("w_{2} ");
    hLTM_w2_pt->GetXaxis()->SetTitle("#it{p}_{T}(GeV/c)");
    hLTM_w2_pt->GetXaxis()->SetTitleSize(0.06);
    hLTM_w2_pt->GetYaxis()->SetTitleSize(0.06);
    hLTM_w2_pt->GetXaxis()->SetTitleOffset(0.8);
    hLTM_w2_pt->GetYaxis()->SetTitleOffset(0.85);
    hLTM_w2_pt->GetYaxis()->SetRangeUser(0.0,0.1);
    hLTM_w2_pt->SetMarkerStyle(20);
    hLTM_w2_pt->SetMarkerSize(1.2);
    hLTM_w2_pt->SetMarkerColor(kRed+2);
    hLTM_w2_pt->SetLineColor(kRed+2);
    
    Double_t LTMval1 = 40;Double_t LTMval2 = 120;Double_t LTMval3 = 40;Double_t LTMval4 = 60;Double_t LTMval5 = 60;Double_t LTMval6 = 80;Double_t LTMval7 = 80;Double_t LTMval8 = 120;
    if(CentralityLow.CompareTo("20")==0){ LTMval1 = 20; LTMval2 = 80; LTMval3 = 20; LTMval4 = 40; LTMval5 = 40; LTMval6 = 60; LTMval7 = 60; LTMval8 = 80; }
    if(CentralityLow.CompareTo("40")==0){ LTMval1 = 5; LTMval2 = 40; LTMval3 = 5; LTMval4 = 20; LTMval5 = 20; LTMval6 = 30; LTMval7 = 30; LTMval8 = 40; }
    
    TLatex T1;
    T1.SetTextSize(0.06);
    T1.SetTextAlign(12);
    T1.SetNDC();
    
    Double_t xlow, xhigh;
    if(CentralityLow.CompareTo("0")==0){ xlow = 40; xhigh = 120;}
    if(CentralityLow.CompareTo("20")==0){ xlow = 20; xhigh = 80;}
    if(CentralityLow.CompareTo("40")==0){ xlow = 5; xhigh = 40;}
    TF1* fit = new TF1("fit","pol0",xlow,xhigh);
    
    gSystem->mkdir("Results");
    gSystem->mkdir(Form("Results/%s%sLTM",CentralityLow.Data(),CentralityHigh.Data()));
    
    //==================================================================================
    //MAIN LOOP OVER P_{T}
    //==================================================================================

     for(int i=0; i<18; i++){
   
      Double_t PtLow 	= newBinsComb[i];
      Double_t PtHigh 	= newBinsComb[i+1];
      
      h2d_LTM_data->GetYaxis()->SetRangeUser(PtLow,PtHigh);
      h2d_LTM_MC->GetYaxis()->SetRangeUser(PtLow,PtHigh);
      h2d_LTM_MCgen->GetYaxis()->SetRangeUser(PtLow,PtHigh);

      TH1D* h1d_LTM_data = (TH1D*)h2d_LTM_data->ProjectionX();
      TH1D* h1d_LTM_MC = (TH1D*)h2d_LTM_MC->ProjectionX();
      TH1D* h1d_LTM_MCgen = (TH1D*)h2d_LTM_MCgen->ProjectionX();
      
      h1d_LTM_data->SetTitle("");
      h1d_LTM_data->GetYaxis()->SetTitle("N  ");
      h1d_LTM_data->GetXaxis()->SetTitle("LTM");
      
      h1d_LTM_MC->SetTitle("");
      h1d_LTM_MC->GetYaxis()->SetTitle("N  ");
      h1d_LTM_MC->GetXaxis()->SetTitle("LTM");
      
      h1d_LTM_MCgen->SetTitle("");
      h1d_LTM_MCgen->GetYaxis()->SetTitle("N  ");
      h1d_LTM_MCgen->GetXaxis()->SetTitle("LTM");
      
      TH1D* LTM_eff = (TH1D*)h1d_LTM_MC->Clone();
      LTM_eff->Divide(h1d_LTM_MCgen);
      LTM_eff->GetYaxis()->SetTitle("#epsilon  ");
      LTM_eff->GetXaxis()->SetTitle("LTM");
      LTM_eff->Fit("fit","RN");
      LTM_eff->GetYaxis()->SetRangeUser(0,1);
      cout << "fit par0= " << fit->GetParameter(0) << endl;
      
      TH1D* LTM_eff_dev = (TH1D*)LTM_eff->Clone();
      LTM_eff_dev->GetYaxis()->SetTitle("#epsilon-fit  ");
      LTM_eff_dev->GetXaxis()->SetTitle("LTM");
      
      for(Int_t j=0; j<LTM_eff->GetNbinsX(); j++){
        if(LTM_eff->FindBin(LTM_eff_dev->GetBinCenter(j))>xlow && LTM_eff->FindBin(LTM_eff_dev->GetBinCenter(j))<xhigh){
          LTM_eff_dev->SetBinContent(LTM_eff->FindBin(LTM_eff_dev->GetBinCenter(j)),(LTM_eff->GetBinContent(LTM_eff->FindBin(LTM_eff_dev->GetBinCenter(j)))-fit->GetParameter(0)));
        }else{
          LTM_eff_dev->SetBinContent(LTM_eff->FindBin(LTM_eff_dev->GetBinCenter(j)),0);
        }
      }
      
      cout << "ptlow= " << PtLow << " bincenter= " << hEff_pt->GetBinCenter(i+1) << " pthigh= " << PtHigh << endl;
      hEff_pt->SetBinContent(i+1,h1d_LTM_MC->Integral(h1d_LTM_MC->FindBin(LTMval1),h1d_LTM_MC->FindBin(LTMval2))/h1d_LTM_MCgen->Integral(h1d_LTM_MCgen->FindBin(LTMval1),h1d_LTM_MCgen->FindBin(LTMval2)));
      hEff_pt_1->SetBinContent(i+1,h1d_LTM_MC->Integral(h1d_LTM_MC->FindBin(LTMval3),h1d_LTM_MC->FindBin(LTMval4))/h1d_LTM_MCgen->Integral(h1d_LTM_MCgen->FindBin(LTMval3),h1d_LTM_MCgen->FindBin(LTMval4)));
      hEff_pt_2->SetBinContent(i+1,h1d_LTM_MC->Integral(h1d_LTM_MC->FindBin(LTMval5),h1d_LTM_MC->FindBin(LTMval6))/h1d_LTM_MCgen->Integral(h1d_LTM_MCgen->FindBin(LTMval5),h1d_LTM_MCgen->FindBin(LTMval6)));
      hEff_pt_3->SetBinContent(i+1,h1d_LTM_MC->Integral(h1d_LTM_MC->FindBin(LTMval7),h1d_LTM_MC->FindBin(LTMval8))/h1d_LTM_MCgen->Integral(h1d_LTM_MCgen->FindBin(LTMval7),h1d_LTM_MCgen->FindBin(LTMval8)));
      
      LTM_eff_dev->GetYaxis()->SetRangeUser(LTM_eff_dev->GetMinimum()*1.2,LTM_eff_dev->GetMaximum()*1.2);
      
      Double_t LTM_Mean_Data = h1d_LTM_data->GetBinCenter(h1d_LTM_data->GetMaximumBin());
      Double_t InclusiveValue = histoInclusive->GetBinContent(histoInclusive->FindBin((newBinsComb[i]+newBinsComb[i+1])/2));
      cout << "InclusiveValue =  " << InclusiveValue << endl;
      Double_t LTM_Low_Data = LTM_Mean_Data*(1-2*InclusiveValue);
      Double_t LTM_High_Data = LTM_Mean_Data*(1+2*InclusiveValue);
      cout << LTM_Low_Data << " " << LTM_Mean_Data << " " << LTM_High_Data << endl;
      
      Double_t LTM_var = TMath::Sqrt( pow((LTM_eff->GetBinContent(LTM_eff->GetBin(LTM_Mean_Data))-LTM_eff->GetBinContent(LTM_eff->GetBin(LTM_Low_Data))),2) + pow((LTM_eff->GetBinContent(LTM_eff->GetBin(LTM_Mean_Data))-LTM_eff->GetBinContent(LTM_eff->GetBin(LTM_High_Data))),2) ) / TMath::Sqrt(2);
      
      Double_t LTM_w2 = LTM_var / (4 * ESD_ConvGamma_efficiency->GetBinContent(ESD_ConvGamma_efficiency->FindBin((newBinsComb[i]+newBinsComb[i+1])/2)));
      cout << "ptbin = " << hEff_pt->GetBinCenter(i+1) << " eff = " << ESD_ConvGamma_efficiency->GetBinContent(ESD_ConvGamma_efficiency->FindBin((newBinsComb[i]+newBinsComb[i+1])/2)) << " w2 = " << LTM_w2 << endl;
      
      histoInclusive_corrected_LTM->SetBinContent(histoInclusive->FindBin((newBinsComb[i]+newBinsComb[i+1])/2),InclusiveValue+LTM_w2);
      
      cout << LTM_Low_Data << " " << LTM_Mean_Data << " " << LTM_High_Data << " variation =  " << LTM_var << endl;
      
      effvariations[i] = LTM_var;
      if(i>0 && i<17){
        hLTM_Eff_pt->SetBinContent(i+1,LTM_var);
        hLTM_Eff_pt->SetBinError(i+1,LTM_var / 2);
        hLTM_w2_pt->SetBinContent(i+1,LTM_w2);
        hLTM_w2_pt->SetBinError(i+1,TMath::Sqrt(pow(1/(4*ESD_ConvGamma_efficiency->GetBinContent(ESD_ConvGamma_efficiency->FindBin((newBinsComb[i]+newBinsComb[i+1])/2))),2)*pow(LTM_var / TMath::Sqrt(2),2)+pow(LTM_var/(4*ESD_ConvGamma_efficiency->GetBinContent(ESD_ConvGamma_efficiency->FindBin((newBinsComb[i]+newBinsComb[i+1])/2))*ESD_ConvGamma_efficiency->GetBinContent(ESD_ConvGamma_efficiency->FindBin((newBinsComb[i]+newBinsComb[i+1])/2))),2)*pow(ESD_ConvGamma_efficiency->GetBinError(ESD_ConvGamma_efficiency->FindBin((newBinsComb[i]+newBinsComb[i+1])/2)),2)));
      }
      
      TCanvas* c1 = new TCanvas("c1","",1200,800);
      gStyle->SetOptStat(0);
      c1->Divide(3,2,0.001,0.001);
      
      c1->cd(1);
      h1d_LTM_data->Draw("HIST");
      T1.DrawLatex(0.75, 0.15, "data");
      T1.DrawLatex(0.4, 0.9, Form("%2.1f < P_{t}(GeV/c) < %2.1f",newBinsComb[i],newBinsComb[i+1]));
      
      c1->cd(2);
      h1d_LTM_MC->Draw("HIST");
      T1.DrawLatex(0.75, 0.15, "MC");
      
      c1->cd(3);
      h1d_LTM_MCgen->Draw("HIST");
      T1.DrawLatex(0.75, 0.15, "MCgen");
      
      c1->cd(4);
      LTM_eff->Draw("HIST");
      fit->Draw("SAME");
      T1.DrawLatex(0.65, 0.15, "MC/MCgen");
      T1.DrawLatex(0.65, 0.25, Form("%2.2f",h1d_LTM_MC->GetEntries()/h1d_LTM_MCgen->GetEntries()));
      
      c1->cd(5);
      LTM_eff_dev->Draw("HIST");
      
      c1->SaveAs(Form("Results/%s%sLTM/LTM_%s-%s.eps",CentralityLow.Data(),CentralityHigh.Data(),newBinsName[i].Data(),newBinsName[i+1].Data()));
      
    }
    
    Double_t legxlow = 0.16; Double_t legxhigh = 0.4; Double_t legylow = 0.65; Double_t legyhigh = 0.85;
    if(CentralityLow.CompareTo("40")==0){legxlow = 0.5; legxhigh = 0.75; legylow = 0.2; legyhigh = 0.4;}
    TLegend* leg = new TLegend(legxlow,legylow,legxhigh,legyhigh);
    leg->SetTextSize(0.04);
    leg->AddEntry(hEff_pt,Form("%2.1f<LTM<%2.1f",LTMval1,LTMval2),"lp");
    leg->AddEntry(hEff_pt_1,Form("%2.1f<LTM<%2.1f",LTMval3,LTMval4),"lp");
    leg->AddEntry(hEff_pt_2,Form("%2.1f<LTM<%2.1f",LTMval5,LTMval6),"lp");
    leg->AddEntry(hEff_pt_3,Form("%2.1f<LTM<%2.1f",LTMval7,LTMval8),"lp");
    leg->SetBorderSize(0);
    
    TCanvas* c3 = new TCanvas("c3","",800,800);
    gStyle->SetOptStat(0);
    SetProperMargins();
    hEff_pt->Draw("p");
    hEff_pt_1->Draw("pSAME");
    hEff_pt_2->Draw("pSAME");
    hEff_pt_3->Draw("pSAME");
    leg->Draw();
    DrawInfoLabelOnPlot(CentralityLow.Data(),CentralityHigh.Data(),2);
    
    TCanvas* c4 = new TCanvas("c4","",800,800);
    gStyle->SetOptStat(0);
    SetProperMargins();
    hLTM_Eff_pt->Draw("p");
    DrawInfoLabelOnPlot(CentralityLow.Data(),CentralityHigh.Data(),2);
    
    TCanvas* c5 = new TCanvas("c5","",800,800);
    gStyle->SetOptStat(0);
    SetProperMargins();
    hLTM_w2_pt->Draw("p");
    DrawInfoLabelOnPlot(CentralityLow.Data(),CentralityHigh.Data(),2);
    
    TLegend* leg1 = new TLegend(0.16,0.75,0.6,0.86);
    leg1->SetTextSize(0.05);
    leg1->AddEntry(histoInclusive,"v_{2,meas}","lp");
    leg1->AddEntry(histoInclusive_corrected_LTM,"v_{2,cor} = v_{2,meas}+w_{2}","lp");
    leg1->SetBorderSize(0);
    
    TCanvas* c6 = new TCanvas("c6","",800,800);
    gStyle->SetOptStat(0);
    SetProperMargins();
    histoInclusiveEmpty->Draw();
    histoInclusive->Draw("SAME");
    histoInclusive_corrected_LTM->Draw("SAME");
    leg1->Draw();
    DrawInfoLabelOnPlot(CentralityLow.Data(),CentralityHigh.Data(),2);
    
    TCanvas* c7 = new TCanvas("c7","",800,800);
    gStyle->SetOptStat(0);
    SetProperMargins();
    ESD_ConvGamma_efficiency->GetXaxis()->SetRangeUser(0,6.2);
    ESD_ConvGamma_efficiency->GetYaxis()->SetRangeUser(0,1);
    ESD_ConvGamma_efficiency->Draw("p");
    DrawInfoLabelOnPlot(CentralityLow.Data(),CentralityHigh.Data(),2);
    
    TCanvas* c8 = new TCanvas("c8","",800,800);
    gStyle->SetOptStat(0);
    SetProperMargins();
    histoInclusiveEmptyRatio->Draw();
    histoInclusive_corrected_LTM_ratio->Divide(histoInclusive_corrected_LTM,histoInclusive,1,1,"B");
    MaskPoints(histoInclusive_corrected_LTM_ratio,0,10);
    histoInclusive_corrected_LTM_ratio->Draw("SAME");
    DrawInfoLabelOnPlot(CentralityLow.Data(),CentralityHigh.Data(),2);
    
    c3->SaveAs(Form("Results/%s%sLTM/%s%seff_LTM.eps",CentralityLow.Data(),CentralityHigh.Data(),CentralityLow.Data(),CentralityHigh.Data()));
    c4->SaveAs(Form("Results/%s%sLTM/%s%sIO_plane_eff_LTM.eps",CentralityLow.Data(),CentralityHigh.Data(),CentralityLow.Data(),CentralityHigh.Data()));
    c5->SaveAs(Form("Results/%s%sLTM/%s%sIO_plane_w2_LTM.eps",CentralityLow.Data(),CentralityHigh.Data(),CentralityLow.Data(),CentralityHigh.Data()));
    c6->SaveAs(Form("Results/%s%sLTM/%s%s_v2_corrected_LTM.eps",CentralityLow.Data(),CentralityHigh.Data(),CentralityLow.Data(),CentralityHigh.Data()));
    c7->SaveAs(Form("Results/%s%sLTM/%s%s_v2_efficiency_from_spectrum.eps",CentralityLow.Data(),CentralityHigh.Data(),CentralityLow.Data(),CentralityHigh.Data()));
    c8->SaveAs(Form("Results/%s%sLTM/%s%s_v2_corrected_LTM_ratio.eps",CentralityLow.Data(),CentralityHigh.Data(),CentralityLow.Data(),CentralityHigh.Data()));
    
    TFile *output_File = new TFile(Form("Results/%s%sLTM_efficiencies.root",CentralityLow.Data(),CentralityHigh.Data()),"RECREATE");
    hLTM_Eff_pt->Write("hLTM_Eff_pt");
    ESD_ConvGamma_efficiency->Write("ConvGamma_efficiency_spectrum");
    output_File->Close();
    
}
