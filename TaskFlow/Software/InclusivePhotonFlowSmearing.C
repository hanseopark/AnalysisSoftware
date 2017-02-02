#include "DirectPhotonFlowFunctions.h"
#include "TMatrixD.h"

void  InclusivePhotonFlowSmearing(
                                      Int_t Trainconfig = 55,
                                      TString CentralityLow = "20",
                                      TString CentralityHigh = "40",
                                      TString Cutnumber = "52400013_00200009007000008250400000",
                                      TString Cutnumbereff = "52400013_00200009300000008250400000_0152304500900000"
                                 ){
  
  //========================================================
  //Opening files and creating histograms
  //========================================================
  
  TFile* fileInclusive  = new TFile(Form("/home/mike/git_afterburner/AnalysisSoftware/TaskFlow/Results/PCM_InclusivePhotonFlow_Syst_%s%s.root",CentralityLow.Data(),CentralityHigh.Data()));
  TH1F*  histoInclusive = (TH1F*)fileInclusive->Get(Form("histoInclusive_%s%s",CentralityLow.Data(),CentralityHigh.Data()));
  if(!histoInclusive) cout << "histoInclusive not found in fileInclusive!!" << endl;
  
  TH1F*  histoInclusiveEmptyRatio = (TH1F*)Getv2HistStyle();
  histoInclusiveEmptyRatio->GetYaxis()->SetRangeUser(0.851,1.149);
  histoInclusiveEmptyRatio->GetYaxis()->SetTitle("Ratio");
  
  TFile* fileMC = new TFile(Form("/home/mike/3_PbPb_dirg/3_GridJobs/161215_PbPb_v2_MC_test/GammaConvFlow_%i.root",Trainconfig));
  TList* listMC_1   = new TList();
  listMC_1     = (TList*)fileMC->Get(Form("GammaConvV1_%i_v2",Trainconfig));
  cout << listMC_1 << endl;
  TList* listMC_2   = new TList();
  listMC_2     = (TList*)listMC_1->FindObject(Form("Cut Number %s",Cutnumber.Data()));
  cout << listMC_2 << endl;
  TList* listMC_3   = new TList();
  listMC_3     = (TList*)listMC_2->FindObject(Form("%s ESD histograms",Cutnumber.Data()));
  cout << listMC_3 << endl;
  
  TH2F* hPt_TruePt = (TH2F*)listMC_3->FindObject("hPt_TruePt");
  cout << hPt_TruePt << endl;
  
  TH2F*  hPt_TruePt_copy = (TH2F*)hPt_TruePt->Clone("hPt_TruePt_copy");
  TH1D* hPt_TruePt_1D_X = (TH1D*)hPt_TruePt_copy->ProjectionX();
  TH1D* hPt_TruePt_1D_Y = (TH1D*)hPt_TruePt_copy->ProjectionY();
  hPt_TruePt_1D_Y->SetLineColor(kRed);
  hPt_TruePt_1D_Y->SetMarkerColor(kRed);
  hPt_TruePt_1D_X->SetTitle("");
  
  hPt_TruePt->SetTitle("");
  hPt_TruePt->GetXaxis()->SetTitle("rec. P_{T}(GeV/c)");
  hPt_TruePt->GetYaxis()->SetTitle("true P_{T}(GeV/c)");
  hPt_TruePt->GetXaxis()->SetTitleSize(0.06);
  hPt_TruePt->GetYaxis()->SetTitleSize(0.07);
  hPt_TruePt->GetXaxis()->SetTitleOffset(0.8);
  hPt_TruePt->GetYaxis()->SetTitleOffset(0.80);
  
  TString InfoSystem = Form("%s-%s %% PbPb, #sqrt{s}=2.76TeV",CentralityLow.Data(),CentralityHigh.Data());
  TString InfoMC = "MC LHC13d2";
  
  //==================================================================================
  //CALCULATING THE SMEARING FROM MC HISTOGRAM
  //==================================================================================
  
  Double_t newBinsComb[17] = {0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 3.0, 3.3, 3.7, 4.1, 4.6, 5.4, 6.2};
  TString newBinsName[17] = {"09", "11", "13", "15", "17", "19", "21", "23", "25", "27", "30", "33", "37", "41", "46", "54", "62"};
  Double_t centralbinvalues[16];
  Double_t binwidths[16];
  
  for(int i=0;i<16;i++){
    centralbinvalues[i]=newBinsComb[i]+(newBinsComb[i+1]-newBinsComb[i])/2;
    binwidths[i]=(newBinsComb[i+1]-newBinsComb[i])/2;
  }
  
  TLatex T1;
  T1.SetTextSize(0.04);
  T1.SetTextAlign(12);
  T1.SetNDC();
  gStyle->SetOptStat(0);
  
//   TH1F*  histopTcorrections = (TH1F*)histoInclusive->Clone();
  
  TCanvas* cProj = new TCanvas("cProj","",400,400);
  
  Double_t hMeans[16];
  TMatrixD m(29,29);
  Double_t m_pt[29];
//   Double_t pTCorrection[16];
  
  for(int i=0; i<16; i++){
    
    Double_t PtLow 	= newBinsComb[i];
    Double_t PtHigh 	= newBinsComb[i+1];
        
    hPt_TruePt->GetYaxis()->SetRangeUser(PtLow,PtHigh);

    TH1D* hPt_TruePt_1D = (TH1D*)hPt_TruePt->ProjectionX();
    hPt_TruePt_1D->SetTitle("");
    hPt_TruePt_1D->GetYaxis()->SetTitle("N  ");
    hPt_TruePt_1D->GetXaxis()->SetTitle("rec p_{T}(GeV/c)");
    hPt_TruePt_1D->GetXaxis()->SetTitleSize(0.06);
    hPt_TruePt_1D->GetYaxis()->SetTitleSize(0.07);
    hPt_TruePt_1D->GetXaxis()->SetRangeUser(0,7.0);
    hPt_TruePt_1D->GetXaxis()->SetTitleOffset(0.8);
    hPt_TruePt_1D->GetYaxis()->SetTitleOffset(0.85);
    hPt_TruePt_1D->SetFillColor(12);
    hPt_TruePt_1D->SetLineColor(12);
    hPt_TruePt_1D->SetFillStyle(3004);
//     pTCorrection[i] = hMeans[i]/centralbinvalues[i];
//     histopTcorrections->SetBinContent(histopTcorrections->FindBin(centralbinvalues[i]),pTCorrection[i]);
    
    TH1F*  hPt_TruePt_1D_rebinned = (TH1F*)histoInclusive->Clone();
    hPt_TruePt_1D_rebinned->SetTitle("");
    hPt_TruePt_1D_rebinned->GetYaxis()->SetTitle("N  ");
    hPt_TruePt_1D_rebinned->GetXaxis()->SetTitle("rec p_{T}(GeV/c)");
    hPt_TruePt_1D_rebinned->GetXaxis()->SetTitleSize(0.06);
    hPt_TruePt_1D_rebinned->GetYaxis()->SetTitleSize(0.07);
    hPt_TruePt_1D_rebinned->GetXaxis()->SetRangeUser(0,7.0);
    hPt_TruePt_1D_rebinned->GetXaxis()->SetTitleOffset(0.8);
    hPt_TruePt_1D_rebinned->GetYaxis()->SetTitleOffset(0.85);
    hPt_TruePt_1D_rebinned->SetFillColor(12);
    hPt_TruePt_1D_rebinned->SetLineColor(12);
    hPt_TruePt_1D_rebinned->SetFillStyle(3004);
    for(Int_t i=0;i<hPt_TruePt_1D_rebinned->GetNbinsX();i++){
      hPt_TruePt_1D_rebinned->SetBinContent(i,0);
    }
    
    for(int j=0; j<hPt_TruePt_1D->GetNbinsX(); j++){
//       cout << "filling pt " << hPt_TruePt_1D->GetBinCenter(j) << " with amount " << hPt_TruePt_1D->GetBinContent(j) << endl;
//       cout << " pt rebinned value " << hPt_TruePt_1D_rebinned->GetBinCenter(hPt_TruePt_1D_rebinned->FindBin(hPt_TruePt_1D->GetBinCenter(j))) << endl;
      hPt_TruePt_1D_rebinned->SetBinContent(hPt_TruePt_1D_rebinned->FindBin(hPt_TruePt_1D->GetBinCenter(j)),hPt_TruePt_1D_rebinned->GetBinContent(hPt_TruePt_1D_rebinned->FindBin(hPt_TruePt_1D->GetBinCenter(j)))+hPt_TruePt_1D->GetBinContent(j));
    }
    
    hMeans[i] = hPt_TruePt_1D_rebinned->GetMean();
    hPt_TruePt_1D_rebinned->Scale(1./hPt_TruePt_1D_rebinned->Integral());
    
    for(int j=0; j<hPt_TruePt_1D_rebinned->GetNbinsX(); j++){
      m[j][i] = hPt_TruePt_1D_rebinned->GetBinContent(j);
      m_pt[j] = hPt_TruePt_1D_rebinned->GetBinCenter(j);
    }
    
    cProj->cd();
    SetProperMargins();
    cProj->SetTopMargin(0.05);
    hPt_TruePt_1D_rebinned->Draw("HIST");
    T1.DrawLatex(0.15, 0.97, InfoSystem.Data());
    T1.DrawLatex(0.65, 0.97, InfoMC.Data());
    if(i>11){
      T1.DrawLatex(0.18, 0.85, Form("%2.1f < P_{T}(GeV/c) < %2.1f",newBinsComb[i],newBinsComb[i+1]));
      T1.DrawLatex(0.18, 0.80, Form("Bincenter = %2.1f",centralbinvalues[i]));
      T1.DrawLatex(0.18, 0.75, Form("Mean = %2.3f",hMeans[i]));
      T1.DrawLatex(0.18, 0.70, Form("Max = %2.3f",hPt_TruePt_1D_rebinned->GetXaxis()->GetBinCenter(hPt_TruePt_1D_rebinned->GetMaximumBin())));
      T1.DrawLatex(0.18, 0.66, Form("leak = %2.3f",1.-hPt_TruePt_1D_rebinned->GetBinContent(hPt_TruePt_1D_rebinned->GetMaximumBin())));
    }else{
      T1.DrawLatex(0.58, 0.85, Form("%2.1f < P_{T}(GeV/c) < %2.1f",newBinsComb[i],newBinsComb[i+1]));
      T1.DrawLatex(0.58, 0.80, Form("Bincenter = %2.1f",centralbinvalues[i]));
      T1.DrawLatex(0.58, 0.75, Form("Mean = %2.3f",hMeans[i]));
      T1.DrawLatex(0.58, 0.70, Form("Max = %2.3f",hPt_TruePt_1D_rebinned->GetXaxis()->GetBinCenter(hPt_TruePt_1D_rebinned->GetMaximumBin())));
      T1.DrawLatex(0.58, 0.66, Form("leak = %2.3f",1.-hPt_TruePt_1D_rebinned->GetBinContent(hPt_TruePt_1D_rebinned->GetMaximumBin())));
    }
    
    gSystem->mkdir(Form("smearing_studies_%s",Cutnumber.Data()));
    cProj->SaveAs(Form("smearing_studies_%s/smearing_%s-%s.eps",Cutnumber.Data(),newBinsName[i].Data(),newBinsName[i+1].Data()));
  
  }
  
  for(int pt=0;pt<29;pt++){
    cout << "pt = " << histoInclusive->GetBinCenter(pt) << endl;
  }
  for(int row=0;row<29;row++){
    for(int col=0;col<29;col++){
      cout << fixed;
      cout << setprecision(4);
      cout << setw(6) << m[row][col] << " ";
    }
    cout << endl;
  }
  for(int k=0;k<29;k++){
    cout << k << " " << m_pt[k] << endl;
  }
  
  
  //========================================================
  //TOY
  //========================================================
  
  TFile* fileData = new TFile("/home/mike/3_PbPb_dirg/1_data/161214_PbPb_v2_MC/GammaConvV1_312_RGamma.root");
  TList* listData_1   = new TList();
  listData_1     = (TList*)fileData->Get("GammaConvV1");
  cout << "listData_1 " << listData_1 << endl;
  TList* listData_2   = new TList();
  listData_2     = (TList*)listData_1->FindObject(Form("Cut Number %s",Cutnumbereff.Data()));
  cout << "listData_2 " << listData_2 << endl;
  TList* listData_3   = new TList();
  listData_3     = (TList*)listData_2->FindObject(Form("%s True histograms",Cutnumbereff.Data()));
  cout << "listData_3 " << listData_3 << endl;
  
  TH1F* hGammaYield = (TH1F*)listData_3->FindObject("ESD_TrueConvGamma_Pt");
  hGammaYield->SetTitle("");
  
  TF1 *fitYield = new TF1("fitYield","[0]*x^(-1.*[1])",0.9,5.4);
  fitYield->SetParameter(0,10e6);
  fitYield->SetParameter(1,6);
  hGammaYield->Fit("fitYield","R");
  
  TRandom3* rand = new TRandom3();
  rand->SetSeed(0);
  double_t random_value;
  
  Double_t y[10000], v2[10000];
  TH1F* hRandYield = new TH1F("hRandYield","",16,newBinsComb);
  hRandYield->SetLineColor(kBlack);
  hRandYield->SetMarkerColor(kBlack);
  hRandYield->SetMarkerStyle(25);
  hRandYield->GetXaxis()->SetTitle("p_{T}(GeV/c)");
  TH1F* hRandYieldRec = new TH1F("hRandYieldRec","",16,newBinsComb);
  hRandYieldRec->SetLineColor(kRed);
  hRandYieldRec->SetMarkerColor(kRed);
  hRandYieldRec->SetMarkerStyle(24);
  TH1F* hRandv2 = new TH1F("hRandv2","",16,newBinsComb);
  hRandv2->SetLineColor(kBlack);
  hRandv2->SetMarkerColor(kBlack);
  hRandv2->SetMarkerStyle(25);
  TH1F* hRandv2Rec = new TH1F("hRandv2Rec","",16,newBinsComb);
  hRandv2Rec->SetLineColor(kRed);
  hRandv2Rec->SetMarkerColor(kRed);
  hRandv2Rec->SetMarkerStyle(24);
  for(Int_t i=0;i<10000;i++){
    y[i] = fitYield->GetRandom();
    hRandYield->Fill(y[i]);
    v2[i] = histoInclusive->GetBinContent(histoInclusive->FindBin(y[i]));
    hRandv2->Fill(y[i],v2[i]);
    random_value = rand->Uniform(0,1);
    
    Double_t newBinsComb[17] = {0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 3.0, 3.3, 3.7, 4.1, 4.6, 5.4, 6.2};
    
    Double_t integral = 0;
    Double_t pt_new;
    if(      y[i] > 0.9 && y[i] < 1.1){
      integral = 0;
      for(Int_t c=0;c<29;c++){
        integral += m[c][0];
        if(integral>random_value){
          pt_new = m_pt[c];
          break;
        }
      }
//       cout << "energy was: " << y[i] << " and it now becomes: " << pt_new << endl;
    }else if(y[i] > 1.1 && y[i] < 1.3){
      integral = 0;
      for(Int_t c=0;c<29;c++){
        integral += m[c][1];
        if(integral>random_value){
          pt_new = m_pt[c];
          break;
        }
      }
//       cout << "energy was: " << y[i] << " and it now becomes: " << pt_new << endl;
    }else if(y[i] > 1.3 && y[i] < 1.5){
      integral = 0;
      for(Int_t c=0;c<29;c++){
        integral += m[c][2];
        if(integral>random_value){
          pt_new = m_pt[c];
          break;
        }
      }
//       cout << "energy was: " << y[i] << " and it now becomes: " << pt_new << endl;
    }else if(y[i] > 1.5 && y[i] < 1.7){
      integral = 0;
      for(Int_t c=0;c<29;c++){
        integral += m[c][3];
        if(integral>random_value){
          pt_new = m_pt[c];
          break;
        }
      }
//       cout << "energy was: " << y[i] << " and it now becomes: " << pt_new << endl;
    }else if(y[i] > 1.7 && y[i] < 1.9){
      integral = 0;
      for(Int_t c=0;c<29;c++){
        integral += m[c][4];
        if(integral>random_value){
          pt_new = m_pt[c];
          break;
        }
      }
//       cout << "energy was: " << y[i] << " and it now becomes: " << pt_new << endl;
    }else if(y[i] > 1.9 && y[i] < 2.1){
      integral = 0;
      for(Int_t c=0;c<29;c++){
        integral += m[c][5];
        if(integral>random_value){
          pt_new = m_pt[c];
          break;
        }
      }
//       cout << "energy was: " << y[i] << " and it now becomes: " << pt_new << endl;
    }else if(y[i] > 2.1 && y[i] < 2.3){
      integral = 0;
      for(Int_t c=0;c<29;c++){
        integral += m[c][6];
        if(integral>random_value){
          pt_new = m_pt[c];
          break;
        }
      }
//       cout << "energy was: " << y[i] << " and it now becomes: " << pt_new << endl;
    }else if(y[i] > 2.3 && y[i] < 2.5){
      integral = 0;
      for(Int_t c=0;c<29;c++){
        integral += m[c][7];
        if(integral>random_value){
          pt_new = m_pt[c];
          break;
        }
      }
//       cout << "energy was: " << y[i] << " and it now becomes: " << pt_new << endl;
    }else if(y[i] > 2.5 && y[i] < 2.7){
      integral = 0;
      for(Int_t c=0;c<29;c++){
        integral += m[c][8];
        if(integral>random_value){
          pt_new = m_pt[c];
          break;
        }
      }
//       cout << "energy was: " << y[i] << " and it now becomes: " << pt_new << endl;
    }else if(y[i] > 2.7 && y[i] < 3.0){
      integral = 0;
      for(Int_t c=0;c<29;c++){
        integral += m[c][9];
        if(integral>random_value){
          pt_new = m_pt[c];
          break;
        }
      }
//       cout << "energy was: " << y[i] << " and it now becomes: " << pt_new << endl;
    }else if(y[i] > 3.0 && y[i] < 3.3){
      integral = 0;
      for(Int_t c=0;c<29;c++){
        integral += m[c][10];
        if(integral>random_value){
          pt_new = m_pt[c];
          break;
        }
      }
//       cout << "energy was: " << y[i] << " and it now becomes: " << pt_new << endl;
    }else if(y[i] > 3.3 && y[i] < 3.7){
      integral = 0;
      for(Int_t c=0;c<29;c++){
        integral += m[c][11];
        if(integral>random_value){
          pt_new = m_pt[c];
          break;
        }
      }
//       cout << "energy was: " << y[i] << " and it now becomes: " << pt_new << endl;
    }else if(y[i] > 3.7 && y[i] < 4.1){
      integral = 0;
      for(Int_t c=0;c<29;c++){
        integral += m[c][12];
        if(integral>random_value){
          pt_new = m_pt[c];
          break;
        }
      }
//       cout << "energy was: " << y[i] << " and it now becomes: " << pt_new << endl;
    }else if(y[i] > 4.1 && y[i] < 4.6){
      integral = 0;
      for(Int_t c=0;c<29;c++){
        integral += m[c][13];
        if(integral>random_value){
          pt_new = m_pt[c];
          break;
        }
      }
//       cout << "energy was: " << y[i] << " and it now becomes: " << pt_new << endl;
    }else if(y[i] > 4.6 && y[i] < 5.4){
      integral = 0;
      for(Int_t c=0;c<29;c++){
        integral += m[c][14];
        if(integral>random_value){
          pt_new = m_pt[c];
          break;
        }
      }
//       cout << "energy was: " << y[i] << " and it now becomes: " << pt_new << endl;
    }else if(y[i] > 5.4 && y[i] < 6.2){
      integral = 0;
      for(Int_t c=0;c<29;c++){
        integral += m[c][15];
        if(integral>random_value){
          pt_new = m_pt[c];
          break;
        }
      }
//       cout << "energy was: " << y[i] << " and it now becomes: " << pt_new << endl;
    }
    if(pt_new < 5.4){
      hRandYieldRec->Fill(pt_new);
      hRandv2Rec->Fill(pt_new,v2[i]);
    }
  }
  hRandv2->Divide(hRandYield);
  hRandv2Rec->Divide(hRandYieldRec);
  
  for(Int_t i=0;i<hRandv2->GetNbinsX();i++){
    hRandv2->SetBinError(i,histoInclusive->GetBinError(histoInclusive->FindBin(hRandv2->GetBinCenter(i))));
    hRandv2Rec->SetBinError(i,histoInclusive->GetBinError(histoInclusive->FindBin(hRandv2->GetBinCenter(i))));
  }
  
  TCanvas* cToy1 = new TCanvas("cToy1","",800,800);
  gStyle->SetOptStat(0);
  gPad->SetLogy();
  SetProperMargins();
  hPt_TruePt_1D_X->Draw();
  hPt_TruePt_1D_Y->Draw("SAME");
  DrawInfoLabelOnPlot(CentralityLow.Data(),CentralityHigh.Data(),2);
  
  TH1F*  histoInclusiveEmpty = (TH1F*)Getv2HistStyle();
  histoInclusiveEmpty->GetYaxis()->SetRangeUser(0,0.49);
  
  TLegend* legToy2 = new TLegend(0.26,0.75,0.6,0.86);
  legToy2->SetTextSize(0.05);
  legToy2->AddEntry(hRandYield,"yield, true","lp");
  legToy2->AddEntry(hRandYieldRec,"yield, rec","lp");
  legToy2->SetBorderSize(0);
  
  TCanvas* cToy2 = new TCanvas("cToy2","",800,800);
  gStyle->SetOptStat(0);
//   gPad->SetLogy();
  SetProperMargins();
  hRandYield->Draw("p");
  hRandYieldRec->Draw("SAMEp");
  legToy2->Draw();
  DrawInfoLabelOnPlot(CentralityLow.Data(),CentralityHigh.Data(),2);
  
  TLegend* legToy3 = new TLegend(0.16,0.75,0.6,0.86);
  legToy3->SetTextSize(0.05);
  legToy3->AddEntry(hRandv2,"v_{2,true}","lp");
  legToy3->AddEntry(hRandv2Rec,"v_{2,rec}","lp");
  legToy3->SetBorderSize(0);
  
  TCanvas* cToy3 = new TCanvas("cToy3","",800,800);
  gStyle->SetOptStat(0);
  SetProperMargins();
  histoInclusiveEmpty->Draw();
  hRandv2->Draw("SAMEp");
  hRandv2Rec->Draw("SAMEp");
  legToy3->Draw();
  DrawInfoLabelOnPlot(CentralityLow.Data(),CentralityHigh.Data(),2);
  
  TLine *line1 = new TLine(0,1,6.2,1);
  TLine *line2 = new TLine(0,1.01,6.2,1.01);
  line2->SetLineStyle(2);
  TLine *line3 = new TLine(0,0.99,6.2,0.99);
  line3->SetLineStyle(2);
  
  TCanvas* cToy4 = new TCanvas("cToy4","",800,800);
  gStyle->SetOptStat(0);
  SetProperMargins();
  histoInclusiveEmptyRatio->Draw();
  line1->Draw(); line2->Draw(); line3->Draw();
  T1.DrawLatex(0.18, 0.90, "v_{2,true}/v_{2,rec}");
  TH1F*  hRandv2_ratio = (TH1F*)hRandv2->Clone();
  hRandv2_ratio->Divide(hRandv2,hRandv2Rec,1,1,"B");
  hRandv2_ratio->Draw("SAMEp");
  DrawInfoLabelOnPlot(CentralityLow.Data(),CentralityHigh.Data(),2);
  
  cToy1->SaveAs(Form("smearing_studies_%s/pt_spectrum_%s-%s.eps",Cutnumber.Data(),CentralityLow.Data(),CentralityHigh.Data()));
  cToy2->SaveAs(Form("smearing_studies_%s/rand_pt_spectrum_%s-%s.eps",Cutnumber.Data(),CentralityLow.Data(),CentralityHigh.Data()));
  cToy3->SaveAs(Form("smearing_studies_%s/rand_v2_spectrum_%s-%s.eps",Cutnumber.Data(),CentralityLow.Data(),CentralityHigh.Data()));
  cToy4->SaveAs(Form("smearing_studies_%s/rand_v2_ratio_%s-%s.eps",Cutnumber.Data(),CentralityLow.Data(),CentralityHigh.Data()));
  
  
  
//   //========================================================
//   //Fitting
//   //========================================================
//   TF1* fit_pTShift = new TF1("fit_pTShift","pol0",1.0,5.8);
//   histopTcorrections->Fit("fit_pTShift","R");
//   
//   TF1* fit_uncorrected = new TF1("fit_uncorrected","pol3",1.0,5.8);
//   fit_uncorrected->SetLineColor(kRed+2);
//   histoInclusive->Fit("fit_uncorrected","R");
//   
//   Double_t par0 = fit_uncorrected->GetParameter(0);
//   Double_t par1 = fit_uncorrected->GetParameter(1);
//   Double_t par2 = fit_uncorrected->GetParameter(2);
//   Double_t par3 = fit_uncorrected->GetParameter(3);
//   Double_t par4 = 1/fit_pTShift->GetParameter(0);
//   
//   TF1* fit_shift = new TF1("fit_shift","[0]+[1]*(x*[4])+[2]*(x*[4])*(x*[4])+[3]*(x*[4])*(x*[4])*(x*[4])",1.0,5.8);
//   fit_shift->SetParameter(0,par0);
//   fit_shift->SetParameter(1,par1);
//   fit_shift->SetParameter(2,par2);
//   fit_shift->SetParameter(3,par3);
//   fit_shift->SetParameter(4,par4);
// //   TF1* fit_shift = new TF1("fit_shift","0.98",1.0,5.8);
//   fit_shift->SetLineColor(kGreen+2);
//   
// //   TF1* fit_corrected = new TF1("fit_corrected","fit_uncorrected*fit_shift",1.0,5.8);
// //   fit_corrected->SetLineColor(kGreen+2);
//   
//   TF1* fit_ratio = new TF1("fit_ratio","fit_shift/fit_uncorrected",1.0,5.8);
//   fit_ratio->SetLineColor(kBlack);
//   
//   TLine* line = new TLine(0,1,6.2,1);
//   line->SetLineStyle(2);
//   line->SetLineWidth(2);
//   line->SetLineColor(12);
//   
//   //========================================================
//   //Correcting
//   //========================================================
//   
//   TH1F*  histoInclusive_corrected = (TH1F*)histoInclusive->Clone();
//   histoInclusive_corrected->SetMarkerColor(kGreen+2);
//   histoInclusive_corrected->SetLineColor(kGreen+2);
//   
//   for(Int_t i=0;i<16;i++){
//     histoInclusive_corrected->SetBinContent(histoInclusive_corrected->FindBin(centralbinvalues[i]),histoInclusive->GetBinContent(histoInclusive->FindBin(centralbinvalues[i]))*fit_ratio->Eval(centralbinvalues[i]));
//   }
//   
//   TH1F*  histoInclusive_divided = (TH1F*)histoInclusive_corrected->Clone();
//   histoInclusive_divided->SetMarkerColor(kBlack);
//   histoInclusive_divided->SetLineColor(kBlack);
//   histoInclusive_divided->Divide(histoInclusive_divided,histoInclusive,1,1,"b");
//   MaskPoints(histoInclusive_divided,0,10);
  
//   //========================================================
//   //Plotting and Saving
//   //========================================================
//   
//   gStyle->SetOptStat(0);
//   gStyle->SetPadTickY(1);
//   gStyle->SetPadTickX(1);
//   
//   TLegend* leg = new TLegend(0.16,0.15,0.6,0.26);
//   leg->SetTextSize(0.04);
//   leg->AddEntry(histoInclusive,"v_{2}^{#gamma,incl}","lp");
//   leg->SetBorderSize(0);
//   
//   TCanvas* c1 = new TCanvas("c1","",400,800);
//   c1->Divide(1,2,0.001,0.001);
//   c1->cd(1);
//   SetProperMargins();
//   histoInclusiveEmpty->Draw();
//   histoInclusive->Draw("SAME");
//   fit_shift->Draw("SAME");
//   histoInclusive_corrected->Draw("SAME");
//   DrawInfoLabelOnPlot(CentralityLow.Data(),CentralityHigh.Data(),2);
//   T1.DrawLatex(0.17, 0.94, "Smearing");
//   leg->Draw();
//   
//   c1->cd(2);
//   SetProperMargins();
//   histoInclusiveEmptyRatio->Draw();
//   fit_ratio->Draw("SAME");
//   histoInclusive_divided->Draw("SAME");
//   T1.DrawLatex(0.17, 0.94, "Correction to v_{2}");
//   line->Draw();
//   
//   TCanvas* c2 = new TCanvas("c2","",400,400);
//   SetProperMargins();
//   histoInclusiveEmptyRatio->Draw();
//   histopTcorrections->Draw("SAME");
//   T1.DrawLatex(0.17, 0.94, "Correction to p_{T}");
//   DrawInfoLabelOnPlot(CentralityLow.Data(),CentralityHigh.Data(),2);
//   line->Draw();
//   
//   TCanvas* c3 = new TCanvas("c3","",400,400);
//   gPad->SetLogz(1);
//   SetProperMargins();
//   c3->SetTopMargin(0.05);
//   c3->SetRightMargin(0.05);
//   hPt_TruePt->GetYaxis()->SetRangeUser(0,25);
//   hPt_TruePt->Draw("COLZ");
//   T1.DrawLatex(0.15, 0.97, InfoSystem.Data());
//   T1.DrawLatex(0.65, 0.97, InfoMC.Data());
//   
//   c1->SaveAs(Form("smearing_studies_%s/v2_incl_smearing.eps",Cutnumber.Data()));
//   c2->SaveAs(Form("smearing_studies_%s/pT_correction_smearing.eps",Cutnumber.Data()));
//   c3->SaveAs(Form("smearing_studies_%s/pT_correction_smearing.eps",Cutnumber.Data()));

  
}