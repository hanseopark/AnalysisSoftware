#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TFile.h"
#include <random>
#include "TStopwatch.h"
#include "TLatex.h"
// #include <omp.h>
#include <iostream>


double myfunction(double *x, double *par)
{
   float xx =x[0];
   TF1 *f = new TF1("f","[0]*exp(-0.5*((x-[1])/[2])**2)",0,100);
   f->SetParameters(par[0], par[1], par[2]);
   double value = f->Integral(0,xx);
   return value;
}


void NiceHisto(TH1 *h, double WStyle = 43 , double WSize = 45 , const char *xtitle="",const char *ytitle="", int Color = 1, double MarkerSize = 1.2 ){
  h->SetStats(0);
  h->SetTitle("");
  h->SetMarkerStyle(20);
  h->SetMarkerSize(MarkerSize);
  h->SetLineWidth(MarkerSize);
  h->SetMarkerColor(Color);
  h->SetLineColor(Color);
  h->GetXaxis()->SetTitle(xtitle);
  h->GetYaxis()->SetTitle(ytitle);
  h->GetYaxis()->SetTitleOffset(1.2);
  h->GetXaxis()->SetTitleOffset(1.);
  h->GetYaxis()->SetTitleFont(WStyle);
  h->GetYaxis()->SetTitleSize(WSize);
  h->GetYaxis()->SetLabelFont(WStyle);
  h->GetYaxis()->SetLabelSize(WSize);
  h->GetXaxis()->SetTitleFont(WStyle);
  h->GetXaxis()->SetTitleSize(WSize);
  h->GetXaxis()->SetLabelFont(WStyle);
  h->GetXaxis()->SetLabelSize(WSize);
}

void NiceLatex(Float_t TextSize = 35, Float_t PosX = 0.5, Float_t PosY = 0.5, TString string1 = "", TString string2 = "", TString string3 = "", TString string4 = "", TString string5 = "",  Double_t dDist = 0.065){

    Float_t textFont = 42;
    if(TextSize > 1) textFont = 43;
  Double_t ddistance = dDist;
  if(string1.Sizeof() > 1.){
  TLatex *l = 			new TLatex(PosX, PosY,string1);
						l->SetNDC();
						l->SetTextFont(textFont);
						l->SetTextColor(1);
						l->SetTextSize(TextSize);
						l->Draw("same");
  }
  if(string2.Sizeof() > 1.){
	TLatex *l2 = 		new TLatex(PosX, PosY-dDist,string2);
						l2->SetNDC();
						l2->SetTextFont(textFont);
						l2->SetTextColor(1);
						l2->SetTextSize(TextSize);
						l2->Draw("same");

	dDist += ddistance;
  }
  if(string3.Sizeof() > 1.){
	TLatex *l3 = 		new TLatex(PosX, PosY-dDist,string3);
						l3->SetNDC();
						l3->SetTextFont(textFont);
						l3->SetTextColor(1);
						l3->SetTextSize(TextSize);
						l3->Draw("same");

  dDist += ddistance;
  }

  if(string4.Sizeof() > 1.){
	TLatex *l4 = 		new TLatex(PosX, PosY-dDist,string4);
						l4->SetNDC();
						l4->SetTextFont(textFont);
						l4->SetTextColor(1);
						l4->SetTextSize(TextSize);
						l4->Draw("same");
  dDist += ddistance;
  }

  if(string5.Sizeof() > 1.){
	TLatex *l5 = 		new TLatex(PosX, PosY-dDist,string5);
						l5->SetNDC();
						l5->SetTextFont(textFont);
						l5->SetTextColor(1);
						l5->SetTextSize(TextSize);
						l5->Draw("same");
  dDist += ddistance;
  }
}

int main(){

    TString Sbase = "/media/gustav/external_drive/Data/pp_13TeV_1617_date1703";

    std::vector<TString> SFile = {"ConvCaloCalibration_2_4.root","ConvCaloCalibration_2_5.root","ConvCaloCalibration_2_6.root"};
    std::vector<TString> STopDir = {"ConvCaloCalibration_2_0_4","ConvCaloCalibration_2_0_5","ConvCaloCalibration_2_0_6"};
    std::vector<TString> SFolder = {"Cut Number 2_00010113_4117900067032230000_01631031000000d0","Cut Number 2_0008e113_4117900067032230000_01631031000000d0","Cut Number 2_0008d113_4117911067032230000_01631031000000d0"};
    std::vector<TString> SFolderESD = {"2_00010113_4117900067032230000_01631031000000d0 ESD histograms","2_0008e113_4117900067032230000_01631031000000d0 ESD histograms","2_0008d113_4117911067032230000_01631031000000d0 ESD histograms"};
    TString SHistName = "ClusGamma_E_SM";

    TF1 myfunc("myfunc",myfunction,0,20,3);

    TH1F* hMB = nullptr;
    TH1F* hEG2 = nullptr;
    TH1F* hEG1 = nullptr;
    TH1F* hTriggTurnOnEG2 = nullptr;
    TH1F* hTriggTurnOnEG1 = nullptr;

    TFile fout("f_out.root","RECREATE");

    TH1S h_EG2("EMCalL1G2", "EMCalL1G2", 20, 0, 20);
    TH1S h_EG1("EMCalL1G1", "EMCalL1G1", 20, 0, 20);
    TH1S h_L0("EMCalL0","EMCalL0",1,0,1);

    for(int iSM = 0; iSM <= 20; ++iSM){
        TFile fileINT7(Form("%s/%s",Sbase.Data(),SFile[0].Data()));
        TFile fileEG2(Form("%s/%s",Sbase.Data(),SFile[1].Data()));
        TFile fileEG1(Form("%s/%s",Sbase.Data(),SFile[2].Data()));
        if(iSM != 20){
            hMB = (TH1F*) fileINT7.Get(Form("%s",STopDir[0].Data()))->FindObject(Form("%s",SFolder[0].Data()))->FindObject(Form("%s",SFolderESD[0].Data()))->FindObject(Form("%s%i",SHistName.Data(), iSM));
            hEG2 = (TH1F*) fileEG2.Get(Form("%s",STopDir[1].Data()))->FindObject(Form("%s",SFolder[1].Data()))->FindObject(Form("%s",SFolderESD[1].Data()))->FindObject(Form("%s%i",SHistName.Data(), iSM));
            hEG1 = (TH1F*) fileEG1.Get(Form("%s",STopDir[2].Data()))->FindObject(Form("%s",SFolder[2].Data()))->FindObject(Form("%s",SFolderESD[2].Data()))->FindObject(Form("%s%i",SHistName.Data(), iSM));
        } else {
            hMB = (TH1F*) fileINT7.Get(Form("%s",STopDir[0].Data()))->FindObject(Form("%s",SFolder[0].Data()))->FindObject(Form("%s",SFolderESD[0].Data()))->FindObject("ClusGamma_E");
            hEG2 = (TH1F*) fileEG2.Get(Form("%s",STopDir[1].Data()))->FindObject(Form("%s",SFolder[1].Data()))->FindObject(Form("%s",SFolderESD[1].Data()))->FindObject("ClusGamma_E");
            hEG1 = (TH1F*) fileEG1.Get(Form("%s",STopDir[2].Data()))->FindObject(Form("%s",SFolder[2].Data()))->FindObject(Form("%s",SFolderESD[2].Data()))->FindObject("ClusGamma_E");
        }

        hTriggTurnOnEG2 = (TH1F*) hEG2->Clone("hTriggTurnOnEG2");
        hTriggTurnOnEG2->Divide(hMB);
        hTriggTurnOnEG1 = (TH1F*) hEG1->Clone("hTriggTurnOnEG1");
        hTriggTurnOnEG1->Divide(hMB);

        TCanvas Can("Can","Can",1200,1000);
        Can.cd();
        gPad->SetLeftMargin(0.13);
        myfunc.SetRange(3.8,10);
        myfunc.SetParameters(70,2,0.3);
        myfunc.SetParLimits(0,50,700);
        myfunc.SetParLimits(1,3,4);
        hTriggTurnOnEG2->Fit(&myfunc,"MR0");
        hTriggTurnOnEG2->GetYaxis()->SetRangeUser(0.001,100);
        hTriggTurnOnEG2->GetXaxis()->SetRangeUser(0.,15);
        NiceHisto(hTriggTurnOnEG2, 43,45,"#it{E}_{cluster}", "EG2/INT7",kBlue,3.5);
        hTriggTurnOnEG2->Draw();
        myfunc.DrawCopy("same");
        NiceLatex(43,0.4,0.22,Form("corected mean value: %.02f",myfunc.GetParameter(1) *0.95), Form("corected #sigma value: %.03f",myfunc.GetParameter(2)*0.8));
        if(iSM != 20)h_EG2.SetBinContent(iSM + 1,abs(myfunc.GetParameter(1))*0.95* 100);
        if(iSM != 20)h_EG2.SetBinError(iSM + 1,abs(myfunc.GetParameter(2))*0.8* 100);
        if(iSM == 20) Can.SaveAs("TurnOnEG2.png");
        else Can.SaveAs(Form("TurnOnEG2_SM%i.png",iSM));

        myfunc.SetRange(7.5,16);
        myfunc.SetParLimits(1,7.5,9.5);
        hTriggTurnOnEG1->Fit(&myfunc,"MR0");
        hTriggTurnOnEG1->GetXaxis()->SetRangeUser(5.,20);
        hTriggTurnOnEG1->GetYaxis()->SetRangeUser(0.001,700);
        NiceHisto(hTriggTurnOnEG1, 43,45,"#it{E}_{cluster}", "EG1/INT7",kBlue,4.5);
        hTriggTurnOnEG1->Draw();
        myfunc.DrawCopy("same");
        NiceLatex(43,0.4,0.22,Form("corected mean value: %.02f",myfunc.GetParameter(1)*0.95), Form("corected #sigma value: %.03f",myfunc.GetParameter(2)*0.8));
        if(iSM != 20)h_EG1.SetBinContent(iSM + 1,abs(myfunc.GetParameter(1)) *0.95 * 100);
        if(iSM != 20)h_EG1.SetBinError(iSM + 1,abs(myfunc.GetParameter(2))*0.8* 100);
        if(iSM == 20) Can.SaveAs("TurnOnEG1.png");
        else Can.SaveAs(Form("TurnOnEG1_SM%i.png",iSM));


        // std::cout<<"mean: "<<myfunc.GetParameter(1)<<"\n";
        // std::cout<<"sigma: "<<myfunc.GetParameter(2)<<"\n";
        //
        // std::cout<<"corected mean value: "<<myfunc.GetParameter(1)*3.75/3.9<<"\n";
        // std::cout<<"corected sigma value: "<<myfunc.GetParameter(2)*0.3/0.42<<"\n";


    }
    h_L0.SetBinContent(1,220);
    h_L0.SetBinError(1,20);

    fout.cd();
    h_EG1.Write();
    h_EG2.Write();
    h_L0.Write();
    fout.Close();
    return 1;

}


void GetSMWiseTriggTurnOn(){
    main();
}
