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
//#include "AliHEPDataParser.h"
#include "../CommonHeaders/PlottingGammaConversionHistos.h"
#include "../CommonHeaders/PlottingGammaConversionAdditional.h"
#include "../CommonHeaders/FittingGammaConversion.h"
#include "../CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "../CommonHeaders/ConversionFunctions.h"
#include "CalculateGammaToPi0.h"


TString GetCentralityStringA (TString cutNumber){
   TString centralityCutNumber = cutNumber(0,3);

   if (centralityCutNumber.CompareTo("504") == 0) { // 0-40%
      return "0-20%";
   }else if (centralityCutNumber.CompareTo("601") == 0) { // 0-40%
      return "0-10%";
   }else if (centralityCutNumber.CompareTo("612") == 0) { // 0-40%
      return "0-10%";
   }else if (centralityCutNumber.CompareTo("502") == 0) { // 0-40%
      return "0-20%";
   }else if (centralityCutNumber.CompareTo("501") == 0) { // 0-40%
      return "0-10%";
   }else if (centralityCutNumber.CompareTo("512") == 0) { // 0-40%
      return "10-20%";
   } else if (centralityCutNumber.CompareTo("548") == 0){ //40-80%
      return "40-80%";
   } else if (centralityCutNumber.CompareTo("502") == 0){ //0-20%
      return "0-20%";
   } else if (centralityCutNumber.CompareTo("524") == 0){ //20-40%
      return "20-40%";
   } else if (centralityCutNumber.CompareTo("546") == 0){ //40-60%
      return "40-60%";
   } else if (centralityCutNumber.CompareTo("568") == 0){ //60-80%
      return "60-80%";
   } else if (centralityCutNumber.CompareTo("508") == 0){ //00-80%
      return "0-80%";

   } 	else return "pp";

}

TString GetCentralityStringB (TString cutNumber){
   TString centralityCutNumber = cutNumber(0,3);

   if (centralityCutNumber.CompareTo("504") == 0) { // 0-40%
      return "0-20%";
   }else if (centralityCutNumber.CompareTo("601") == 0) { // 0-40%
      return "0-10%";
   }else if (centralityCutNumber.CompareTo("612") == 0) { // 0-40%
      return "0-10%";
   }else if (centralityCutNumber.CompareTo("502") == 0) { // 0-40%
      return "0-10%";
   }else if (centralityCutNumber.CompareTo("501") == 0) { // 0-40%
      return "0-10%";
   }else if (centralityCutNumber.CompareTo("512") == 0) { // 0-40%
      return "10-20%";
   } else if (centralityCutNumber.CompareTo("548") == 0){ //40-80%
      return "40-80%";
   } else if (centralityCutNumber.CompareTo("502") == 0){ //0-20%
      return "0-10%";
   } else if (centralityCutNumber.CompareTo("524") == 0){ //20-40%
      return "20-40%";
   } else if (centralityCutNumber.CompareTo("546") == 0){ //40-60%
      return "40-60%";
   } else if (centralityCutNumber.CompareTo("568") == 0){ //60-80%
      return "60-80%";
   } else if (centralityCutNumber.CompareTo("508") == 0){ //00-80%
      return "0-80%";
   } 	else return "pp";
}

TString GetCentralityStringC (TString cutNumber){
   TString centralityCutNumber = cutNumber(0,3);

   if (centralityCutNumber.CompareTo("504") == 0) { // 0-40%
      return "0-40%";
   }else if (centralityCutNumber.CompareTo("601") == 0) { // 0-40%
      return "0-5%";
   }else if (centralityCutNumber.CompareTo("612") == 0) { // 0-40%
      return "5-10%";
   }else if (centralityCutNumber.CompareTo("502") == 0) { // 0-40%
      return "0-20%";
   }else if (centralityCutNumber.CompareTo("501") == 0) { // 0-40%
      return "0-10%";
   }else if (centralityCutNumber.CompareTo("512") == 0) { // 0-40%
      return "10-20%";
   } else if (centralityCutNumber.CompareTo("548") == 0){ //40-80%
      return "40-80%";
   } else if (centralityCutNumber.CompareTo("502") == 0){ //0-20%
      return "0-20%";
   } else if (centralityCutNumber.CompareTo("524") == 0){ //20-40%
      return "20-40%";
   } else if (centralityCutNumber.CompareTo("546") == 0){ //40-60%
      return "40-60%";
   } else if (centralityCutNumber.CompareTo("568") == 0){ //60-80%
      return "60-80%";
   } else if (centralityCutNumber.CompareTo("508") == 0){ //00-80%
      return "0-80%";
   } 	else return "pp";
}

Double_t BWdndpt(Double_t *x, Double_t *par) {

   Double_t pt = x[0]*1000.; //GeV->MeV
   Double_t t=par[1];
   Double_t beta=abs(par[2]);
   Double_t yt=0.5*TMath::Log((1+beta)/(1-beta)); // atanh(beta)
   Double_t m0=0.; //pi0
   Double_t mt=TMath::Sqrt(m0*m0+pt*pt);

   Double_t f = par[0]*pt*mt*TMath::BesselI(0,pt*sinh(yt)/t)*TMath::BesselK(1,TMath::Abs(mt*cosh(yt)/t));

   return f;
}
//-----------------------------------------------------------------------------
Double_t FitSpectrum(Double_t *x, Double_t *p)
{
   // Double_t hagd = p[0]*(p[2]-1)*(p[2]-2)/p[1]/p[2]/TMath::TwoPi()/TMath::Power((1.+x[0]/p[1]/p[2]),p[2]);
   Double_t hagd = p[0]/TMath::TwoPi()/TMath::Power((1.+x[0]/p[1]/p[2]),p[2]);
   Double_t expo = p[3]*TMath::Exp(-x[0]/p[4]);
   return hagd+expo;
}
//-----------------------------------------------------------------------------
Double_t Hagedorn(Double_t *x, Double_t *p)
{
   Double_t hagd = p[0]/TMath::TwoPi()/TMath::Power((1.+x[0]/p[1]/p[2]),p[2]);
   return hagd;
}
//-----------------------------------------------------------------------------
Double_t Exponent(Double_t *x, Double_t *p)
{
   Double_t expo = p[0]*TMath::Exp(-x[0]/p[1]);
   return expo;
}

TGraphAsymmErrors *getPhotonPoints(TH1* error, TH1* spectrum, float iFlag, float ScaleArrowlength){

   const int kMaxPoints = 50;
   int nPoints = 0;
   float x[kMaxPoints],xel[kMaxPoints],xeh[kMaxPoints];
   float y[kMaxPoints],yel[kMaxPoints],yeh[kMaxPoints];

   ScaleArrowlength = ScaleArrowlength*1.3;
   //  const float  ScaleArrowlength = 0.5;
   //  const float  ScaleArrowlength = 0.5;

   for(int ib = 1;ib <= error->GetNbinsX();ib++){
      float val =  spectrum->GetBinContent(ib);
      float pt =   spectrum->GetBinCenter(ib);
      float err =  val - error->GetBinContent(ib);    float errt =  error->GetBinContent(ib);
      float errX =  spectrum->GetBinWidth(ib);

      float lowbound = val - 1.28*abs(err); // 90% CL

      if(iFlag==0){ // use only data points with error not consistent with zerp
         if(val>0&&errt>0&&lowbound>=0){
            x[nPoints] = pt;
            xel[nPoints] = errX/2;
            xeh[nPoints] = errX/2;
            y[nPoints] = val;
            yel[nPoints] =  abs(err);
            yeh[nPoints] =  abs(err);
            nPoints++;
         }
      }
      else if(iFlag==1){ // points and errors > 0 not confidencelevel not ok
         if(val>0&&errt>0&&lowbound<=0){
            x[nPoints] = pt;
            xel[nPoints] = errX/2;
            xeh[nPoints] = errX/2;
            y[nPoints] = val;
            yel[nPoints] =  abs(err);
            yeh[nPoints] =  abs(err);
            nPoints++;
         }
      }
      else if(iFlag==10){ // points and errors > 0 not confidencelevel not ok
         if(val>0&&errt>0&&lowbound<=0){
            x[nPoints] = pt;
            xel[nPoints] = 0;
            xeh[nPoints] = 0;
            y[nPoints] = val;
            yel[nPoints] = 0;
            yeh[nPoints] = 1.28*abs(err);
            nPoints++;
         }
      }
      else if(iFlag==0.5){ // points and errors > 0 not confidencelevel not ok
         if(val>0&&errt>0){
            x[nPoints] = pt;
            xel[nPoints] = errX/2;
            xeh[nPoints] = errX/2;
            y[nPoints] = val;
            yel[nPoints] =  abs(err);
            yeh[nPoints] =  abs(err);
            nPoints++;
         }
      }
      else if(iFlag==0.75){ // points and errors > 0 not confidencelevel not ok
         if(val>0&&errt>0){
            x[nPoints] = pt;
            xel[nPoints] = 0;
            xeh[nPoints] = 0;
            y[nPoints] = val;
            yel[nPoints] =  abs(err);
            yeh[nPoints] =  abs(err);
            nPoints++;
         }
      }
      else if(iFlag==2){ // upperlimit for points with  errors consistent with zero
         if(val>0&&errt<=0&&lowbound<=0){
            x[nPoints] = pt;
            xel[nPoints] = 0;
            xeh[nPoints] = 0;
            y[nPoints] = val;
            yel[nPoints] = 0;
            yeh[nPoints] = 1.28*abs(err);
            nPoints++;
         }
      }
      else if(iFlag==3){ // arrow for points with errors consistent with zero
         if(val>0&&errt>0&&lowbound<=0){
            x[nPoints] = pt;
            xel[nPoints] = 0;
            xeh[nPoints] = 0;
            y[nPoints] = val+1.28*abs(err);
            yel[nPoints] = (1.28 * abs(err) + ScaleArrowlength * val);
            yeh[nPoints] = 0;
            nPoints++;
         }
      }
      else if(iFlag==4){ // upperlimit for points with zero
         if(val<=0&&errt<=0){
            x[nPoints] = pt;
            xel[nPoints] = 0;
            xeh[nPoints] = 0;
            y[nPoints] = 1.28*abs(err);
            yel[nPoints] = 1.28*abs(err);
            yeh[nPoints] = 0;
            nPoints++;
         }
      }
      else if(iFlag==5){ // arrow for points with zero
         if(val>0&&errt<0){
            x[nPoints] = pt;
            xel[nPoints] = 0;
            xeh[nPoints] = 0;
            y[nPoints] = 1.28*abs(err)+val;
            //yel[nPoints] = (1.575 * abs(err) + ScaleArrowlength * val);
            yel[nPoints] = 1.28*abs(err)+(val*ScaleArrowlength);
            yeh[nPoints] = 0;
            nPoints++;
         }
      }
      else if(iFlag==6){ // arrow for point consist with zero
         if(val<=0&&errt<=0&&lowbound<=0){
            x[nPoints] = pt;
            xel[nPoints] = 0;
            xeh[nPoints] = 0;
            y[nPoints] = 1.28*abs(err)+val*1.6;
            yel[nPoints] =  (1.48 * abs(err) + ScaleArrowlength * abs(val));
            yeh[nPoints] =  0;
            nPoints++;
         }
      }
   }
   
   if(nPoints > 0){
      TGraphAsymmErrors *gr = new TGraphAsymmErrors(nPoints,x,y,xel,xeh,yel,yeh);
      return gr;
   }

   return 0;
}



void ProduceFinalGammaResults(TString inputFileName = "",TString systematicFileName = "",TString cutSel = "", TString option = ""){

   gSystem->Load("libCore.so");
   gSystem->Load("libTree.so");
   gSystem->Load("libGeom.so");
   gSystem->Load("libVMC.so");
   gSystem->Load("libPhysics.so");
   gSystem->Load("libMinuit");
   gSystem->Load("libSTEERBase");
   gSystem->Load("libESD");
   gSystem->Load("libAOD");
   gSystem->Load("libANALYSIS");
   gSystem->Load("libANALYSISalice");

   StyleSettingsThesis();
   SetPlotStyle();
   TColor *white = gROOT->GetColor(0);
   white->SetAlpha(0.);



   TString outputDir = Form("FinalGammaResults/%s",option.Data());
   gSystem->Exec("mkdir -p "+outputDir);


   TDatime now;
   int iDate = now.GetDate();
   int iYear=iDate/10000;
   int iMonth=(iDate%10000)/100;
   int iDay=iDate%100;
   char* cMonth[12]={"Jan","Feb","Mar","Apr","May","Jun",
                     "Jul","Aug","Sep","Oct","Nov","Dec"};
   char cStamp1[25],cStamp2[25];
   sprintf(cStamp1,"%i_%s_%i",iDay, cMonth[iMonth-1], iYear);
   sprintf(cStamp2,"%i/%.2d/%i",iDay, iMonth, iYear);


   const Int_t ErrorDoubleRatio = 13;
   const Int_t ErrorIncRatio = 10;
   const Int_t ErrorPi0 = 10;
   const Int_t ErrorGammaSpec = 8;

   TString nameCutVariationsDoubleRatio[ErrorDoubleRatio] = {"Chi2","PsiPair","TPCClst","edEdx","PidEdx","qT","SinglePt","TOF","Alpha","IntRange","CocktailEta","CocktailParam","CocktailEtaNorm"};//,"Fit"};
   TString nameCutVariationsIncRatio[ErrorIncRatio] = {"Chi2","PsiPair","TPCClst","edEdx","PidEdx","qT","SinglePt","TOF","Alpha","IntRange"};//,"Fit"};
   TString nameCutVariationsPi0[ErrorPi0] = {"Chi2","PsiPair","TPCClst","edEdx","PidEdx","qT","SinglePt","TOF","Alpha","IntRange"};//,"Fit"};
   TString nameCutVariationsGamma[ErrorGammaSpec] = {"Chi2","PsiPair","TPCClst","edEdx","PidEdx","qT","SinglePt","TOF"};
   TString nameCutVariationsGammaLegend[ErrorDoubleRatio] = {"#chi^{2} Cut","#psi_{pair} Cut","Findable Cls. TPC","e^{#pm} dEdx","#pi^{#pm} dEdx","q_{T} Cut","min. e^{#pm} #it{p}_{T}","TOF Cut","#pi^{0} #alpha_{asym} Cut","#pi^{0} Int. Range","Cocktail #eta Shape","Cocktail #pi^{0} Shape","Cocktail #eta Norm" };

   TFile *fileSystematics = new TFile(systematicFileName);
   TFile *fileInput = new TFile(inputFileName);



   TH1D *DoubleRatio = (TH1D*) fileInput->Get("DoubleRatioConversionTrueEffPurity");
   TH1D *DoubleRatioPi0Fit = (TH1D*) fileInput->Get("DoubleRatioConversionFitPurityA");
   TH1D *IncRatio = (TH1D*) fileInput->Get("IncRatioPurity_trueEff");
   TH1D *IncRatioPi0Fit = (TH1D*) fileInput->Get("histoIncRatioFitPurityA");
   TH1D *Gamma = (TH1D*) fileInput->Get("histoGammaSpecCorrPurity");
   TH1D *Pi0 = (TH1D*) fileInput->Get("CorrectedYieldTrueEff");
   TH1D *Pi0Fit = (TH1D*) fileInput->Get("CorrectedYieldTrueEffPi0Fit");
   TH1D *DoubleRatioChargedPions = (TH1D*) fileInput->Get("DoubleRatioConversionFitPurityPionsAverage");
   TH1D *CocktailChargeToNeutral = (TH1D*) fileInput->Get("cocktailChargeToNeutral");
   TH1D *ChargedPionsSyst = (TH1D*) fileInput->Get("histoChargedPionsSystAverage");
   TH1D *ChargedPionsStat = (TH1D*) fileInput->Get("histoChargedPionsStatAverage");
   TH1D *IncRatioChargedPions = (TH1D*) fileInput->Get("histoIncRatioFitPurityPionsStat");
   TH1D *CocktailRatioChargedPions = (TH1D*) fileInput->Get("cocktailAllGammaPi0Pions");
   CocktailRatioChargedPions->SetName("cocktailAllGammaPi0Pions");
   TH1D *CocktailRatioPi0 = (TH1D*) fileInput->Get("cocktailAllGammaPi0");
   CocktailRatioPi0->SetName("cocktailAllGammaPi0");


   TGraphAsymmErrors *PbPb276CTEQ61EPS09BFG2_sum_pdferr_InvYield = (TGraphAsymmErrors*)fileInput->Get("PbPb276CTEQ61EPS09BFG2_sum_pdferr_InvYield");
   TGraphAsymmErrors *pp276CT10BFG2_sum_pdferr_InvYield = (TGraphAsymmErrors*)fileInput->Get("pp276CT10BFG2_sum_pdferr_InvYield");

   TGraphAsymmErrors *PbPb276EPS09BFG2_sum_scale_InvYield = (TGraphAsymmErrors*)fileInput->Get("PbPb276EPS09BFG2_sum_scale_InvYield");
   TGraphAsymmErrors *pp276CT10BFG2_sum_scale_InvYield = (TGraphAsymmErrors*)fileInput->Get("pp276CT10BFG2_sum_scale_InvYield");

   PbPb276CTEQ61EPS09BFG2_sum_pdferr_InvYield->RemovePoint(PbPb276CTEQ61EPS09BFG2_sum_pdferr_InvYield->GetN()-1);
   pp276CT10BFG2_sum_pdferr_InvYield->RemovePoint(pp276CT10BFG2_sum_pdferr_InvYield->GetN()-1);
   
   PbPb276EPS09BFG2_sum_scale_InvYield->RemovePoint(PbPb276EPS09BFG2_sum_scale_InvYield->GetN()-1);
   pp276CT10BFG2_sum_scale_InvYield->RemovePoint(pp276CT10BFG2_sum_scale_InvYield->GetN()-1);

   for(Int_t i = 0;i<PbPb276CTEQ61EPS09BFG2_sum_pdferr_InvYield->GetN();i++){
      Double_t yerrlow1 = PbPb276CTEQ61EPS09BFG2_sum_pdferr_InvYield->GetErrorYlow(i);
      Double_t yerrlow2 = 1*PbPb276EPS09BFG2_sum_scale_InvYield->GetErrorYlow(i);
      Double_t yerrhigh1 = PbPb276CTEQ61EPS09BFG2_sum_pdferr_InvYield->GetErrorYhigh(i);
      Double_t yerrhigh2 = 1*PbPb276EPS09BFG2_sum_scale_InvYield->GetErrorYhigh(i);
      Double_t xerrlow = PbPb276EPS09BFG2_sum_scale_InvYield->GetErrorXhigh(i);

      PbPb276CTEQ61EPS09BFG2_sum_pdferr_InvYield->SetPointError(i,xerrlow,xerrlow,
                                                                sqrt(yerrlow1*yerrlow1+ yerrlow2*yerrlow2),
                                                                sqrt(yerrhigh1*yerrhigh1+ yerrhigh2*yerrhigh2));
   }

   for(Int_t i = 0;i<pp276CT10BFG2_sum_scale_InvYield->GetN();i++){
      Double_t yerrlow1 = pp276CT10BFG2_sum_pdferr_InvYield->GetErrorYlow(i);
      Double_t yerrlow2 = pp276CT10BFG2_sum_scale_InvYield->GetErrorYlow(i);
      Double_t yerrhigh1 = pp276CT10BFG2_sum_pdferr_InvYield->GetErrorYhigh(i);
      Double_t yerrhigh2 = pp276CT10BFG2_sum_scale_InvYield->GetErrorYhigh(i);
      Double_t xerrlow = pp276CT10BFG2_sum_scale_InvYield->GetErrorXhigh(i);

      pp276CT10BFG2_sum_pdferr_InvYield->SetPointError(i,xerrlow,xerrlow,
                                                       sqrt(yerrlow1*yerrlow1+ yerrlow2*yerrlow2),
                                                       sqrt(yerrhigh1*yerrhigh1+ yerrhigh2*yerrhigh2));
   }



   TFile *PhosGammaFile = new TFile("PHOS_gamma_PbPb_PCM_06022014.root");
   TFile *PhosPi0File = new TFile("PHOS_pi0_PbPb.root");

   Int_t binsX = DoubleRatioPi0Fit->GetNbinsX();

   Double_t *binning = new Double_t[binsX+1];

   for(Int_t i = 0; i<binsX+1;i++){
      binning[i] = DoubleRatioPi0Fit->GetBinLowEdge(i+1);
      //cout<<binning[i]<<endl;
      cout<<DoubleRatioPi0Fit->GetBinCenter(i+1)<<endl;

   }
   Int_t icent = 0;
   if(GetCentralityStringA(cutSel).CompareTo("0-10%") == 0)
      icent = 7;
   if(GetCentralityStringA(cutSel).CompareTo("10-20%") == 0)
      icent = 2;
   if(GetCentralityStringA(cutSel).CompareTo("0-20%") == 0)
      icent = 6;
   if(GetCentralityStringA(cutSel).CompareTo("0-40%") == 0)
      icent = 9;
   if(GetCentralityStringA(cutSel).CompareTo("20-40%") == 0)
      icent = 3;
   if(GetCentralityStringA(cutSel).CompareTo("40-80%") == 0)
      icent = 8;



   TFile *PHOSDoubleRatios = new TFile("PHOS_doubleratio_07022014.root");
   TH1D *DoubleRatioPHOSSyst = (TH1D*)PHOSDoubleRatios->Get(Form("hGamma_PbPb_cen%i_SystRatio",icent));
   TH1D *DoubleRatioPHOSSystA = (TH1D*)PHOSDoubleRatios->Get(Form("hGamma_PbPb_cen%i_SystA",icent));
   TH1D *DoubleRatioPHOSSystB = (TH1D*)DoubleRatioPHOSSystA->Clone(Form("hGamma_PbPb_cen%i_SystB",icent));
   TH1D *DoubleRatioPHOSSystC = (TH1D*)PHOSDoubleRatios->Get(Form("hGamma_PbPb_cen%i_SystC",icent));
   TH1D *DoubleRatioPHOSStat = (TH1D*)PHOSDoubleRatios->Get(Form("hGamma_PbPb_cen%i_Stat",icent));
   // for(Int_t i = 0;i<DoubleRatioPHOSSyst->GetNbinsX();i++){
   //    DoubleRatioPHOSSyst->SetBinError(i+1,sqrt(pow(DoubleRatioPHOSSyst->GetBinError(i+1),2)+pow(DoubleRatioPHOSSystC->GetBinError(i+1),2)));
   // }
   TGraphErrors *GraphDoubleRatioPHOSSyst = new TGraphErrors(DoubleRatioPHOSSyst);
   TGraphErrors *GraphDoubleRatioPHOSSystA = new TGraphErrors(DoubleRatioPHOSSystA);

   TH1D* PhosStatGamma = (TH1D*)PhosGammaFile->Get(Form("hGamma_PbPb_cen%i_Stat",icent));
   TH1D* PhosSystGamma = (TH1D*)PhosGammaFile->Get(Form("hGamma_PbPb_cen%i_Syst",icent));
   TH1D* PhosSystGammaA = (TH1D*)PhosGammaFile->Get(Form("hGamma_PbPb_cen%i_SystA",icent));
   TH1D* PhosSystGammaBnl = (TH1D*)PhosGammaFile->Get(Form("hGamma_PbPb_cen%i_SystBnl",icent));
   TH1D* PhosSystGammaBglobalE = (TH1D*)PhosGammaFile->Get(Form("hGamma_PbPb_cen%i_SystBglobalE",icent));
   TH1D *PhosSystGammaAB = (TH1D*)PhosSystGammaA->Clone("hGamma_PbPb_cen%i_SystSumAB");
   for(Int_t i = 0;i<PhosSystGammaAB->GetNbinsX();i++){
      PhosSystGammaAB->SetBinError(i+1,sqrt( pow(PhosSystGammaAB->GetBinError(i+1),2)+pow(PhosSystGammaBnl->GetBinError(i+1),2)+pow(PhosSystGammaBglobalE->GetBinError(i+1),2)));
   }

   TGraphErrors *GraphGammaPHOSSyst = new TGraphErrors(PhosSystGamma);

   TH1D *histoBinningRatio = new TH1D("","",binsX,binning);

   TF1 *qcdFit = FitObject("qmpt","QCDbinShiftFitForGamma","Pi0");
   TGraphAsymmErrors *gammaAssymError = new TGraphAsymmErrors(Gamma);
   gammaAssymError->RemovePoint(0);
   gammaAssymError = ApplyXshift(gammaAssymError,qcdFit);

   gStyle->SetErrorX(0);
   gStyle->SetMarkerStyle(20);
   gStyle->SetMarkerSize(2);
   gStyle->SetLineWidth(1.);
   gStyle->SetEndErrorSize(5);
   gStyle->SetFillColor(kGray);



   TH1D *cocktailAllGammaPi0 = (TH1D*) fileInput->Get("sumgammapi0");
   // cocktailAllGammaPi0->Rebin(2);
   // cocktailAllGammaPi0->Scale(0.5);
   TH1D *cocktailPi0GammaPi0 = (TH1D*)fileInput->Get("pi0gammapi0");
   // cocktailPi0GammaPi0->Rebin(2);
   // cocktailPi0GammaPi0->Scale(0.5);
   TH1D *cocktailEtaGammaPi0 = (TH1D*)fileInput->Get("etagammapi0");
   // cocktailEtaGammaPi0->Rebin(2);
   // cocktailEtaGammaPi0->Scale(0.5);
   TH1D *cocktailOmegaGammaPi0 = (TH1D*)fileInput->Get("omegagammapi0");
   // cocktailOmegaGammaPi0->Rebin(2);
   // cocktailOmegaGammaPi0->Scale(0.5);
   TH1D *cocktailEtapGammaPi0 = (TH1D*)fileInput->Get("etapgammapi0");
   // cocktailEtapGammaPi0->Rebin(2);
   // cocktailEtapGammaPi0->Scale(0.5);
   TH1D *cocktailPhiGammaPi0 = (TH1D*)fileInput->Get("phigammapi0");
   // cocktailPhiGammaPi0->Rebin(2);
   // cocktailPhiGammaPi0->Scale(0.5);
   TH1D *cocktailRhoGammaPi0 = (TH1D*)fileInput->Get("rhogammapi0");
   // cocktailRhoGammaPi0->Rebin(2);
   // cocktailRhoGammaPi0->Scale(0.5);

   TH1D *cocktailAllGamma = (TH1D* )fileInput->Get("ptg");
   // cocktailAllGamma->Rebin(2);
   // cocktailAllGamma->Scale(0.5);
   TH1D *cocktailPi0Gamma = (TH1D* )fileInput->Get("ptgPi0");
   // cocktailPi0Gamma->Rebin(2);
   // cocktailPi0Gamma->Scale(0.5);
   TH1D *cocktailEtaGamma = (TH1D* )fileInput->Get("ptgEta");
   // cocktailEtaGamma->Rebin(2);
   // cocktailEtaGamma->Scale(0.5);
   TH1D *cocktailEtapGamma = (TH1D* )fileInput->Get("ptgEtaprime");
   // cocktailEtapGamma->Rebin(2);
   // cocktailEtapGamma->Scale(0.5);
   TH1D *cocktailOmegaGamma = (TH1D* )fileInput->Get("ptgOmega");
   // cocktailOmegaGamma->Rebin(2);
   // cocktailOmegaGamma->Scale(0.5);
   TH1D *cocktailPhiGamma = (TH1D* )fileInput->Get("ptgPhi");
   // cocktailPhiGamma->Rebin(2);
   // cocktailPhiGamma->Scale(0.5);
   TH1D *cocktailRhoGamma = (TH1D* )fileInput->Get("ptgRho");
   TH1D *cocktailSigmaGamma = (TH1D* )fileInput->Get("ptgSigma");
   // cocktailRhoGamma->Rebin(2);
   // cocktailRhoGamma->Scale(0.5);
   TH1D *cocktailPion = (TH1D* )fileInput->Get("ptPi0");
   TH1D *cocktailEta = (TH1D* )fileInput->Get("ptEta");
   TH1D *cocktailOmega = (TH1D* )fileInput->Get("ptOmega");
   TH1D *cocktailEtap = (TH1D* )fileInput->Get("ptEtaprime");
   TH1D *cocktailPhi = (TH1D* )fileInput->Get("ptPhi");
   TH1D *cocktailSigma = (TH1D* )fileInput->Get("ptSigma");
   TH1D *K0s = (TH1D* )fileInput->Get("K0s");
   TH1D *cocktailEtaLow = (TH1D*)fileInput->Get("ptEtaLow");
   TH1D *cocktailEtaHigh = (TH1D*)fileInput->Get("ptEtaHigh");
   TH1D *cocktailEtaGammaLow = (TH1D*)fileInput->Get("ptgEtaLow");
   TH1D *cocktailEtaGammaHigh = (TH1D*)fileInput->Get("ptgEtaHigh");


   cocktailEtaLow->DrawCopy();
   cocktailEtaHigh->DrawCopy("same");


   for(Int_t i = 0;i<cocktailEtaLow->GetNbinsX();i++){
      if(cocktailEtaLow->GetBinContent(i+1)>cocktailEtaHigh->GetBinContent(i+1)){
         cocktailEtaLow->SetBinContent(i+1,(cocktailEtaLow->GetBinContent(i+1)/cocktailEtaHigh->GetBinContent(i+1)-1)*0.5);
      }
      else cocktailEtaLow->SetBinContent(i+1,(cocktailEtaHigh->GetBinContent(i+1)/cocktailEtaLow->GetBinContent(i+1)-1)*0.5);
   }

   for(Int_t i = 0;i<cocktailEtaGammaLow->GetNbinsX();i++){
      if(cocktailEtaGammaLow->GetBinContent(i+1)>cocktailEtaGammaHigh->GetBinContent(i+1)){
         cocktailEtaGammaLow->SetBinContent(i+1,(cocktailEtaGammaLow->GetBinContent(i+1)/cocktailEtaGammaHigh->GetBinContent(i+1)-1)*0.5);
      }
      else cocktailEtaGammaLow->SetBinContent(i+1,(cocktailEtaGammaHigh->GetBinContent(i+1)/cocktailEtaGammaLow->GetBinContent(i+1)-1)*0.5);
   }

   
   // cocktailPion->Rebin(2);
   // cocktailPion->Scale(0.5);

   Int_t nBins = cocktailAllGammaPi0->GetNbinsX();
   Int_t deletedBins = 0;
   for(Int_t i = 0; i< cocktailAllGammaPi0->GetNbinsX();i++){
      if(cocktailAllGammaPi0->GetBinLowEdge(i+1)<1)
         deletedBins++;
   }

   Double_t *binsCocktail = new Double_t[nBins-deletedBins];
   for(Int_t i = deletedBins+1; i< cocktailAllGammaPi0->GetNbinsX()+1;i++){
      binsCocktail[i-deletedBins] = cocktailAllGammaPi0->GetBinLowEdge(i);
      cout<<binsCocktail[i-deletedBins]<<endl;
   }
   TH1D *RebinCocktail = new TH1D("RebinCocktail","RebinCocktail",nBins-deletedBins,binsCocktail);
   // cocktailAllGammaPi0 = RebinTH1D(cocktailAllGammaPi0,RebinCocktail);
   // cocktailPi0GammaPi0 = RebinTH1D(cocktailPi0GammaPi0,RebinCocktail);
   // cocktailEtaGammaPi0 = RebinTH1D(cocktailEtaGammaPi0,RebinCocktail);
   // cocktailOmegaGammaPi0 = RebinTH1D(cocktailOmegaGammaPi0,RebinCocktail);
   // cocktailEtapGammaPi0 = RebinTH1D(cocktailEtapGammaPi0,RebinCocktail);
   // cocktailPhiGammaPi0 = RebinTH1D(cocktailPhiGammaPi0,RebinCocktail);
   // cocktailRhoGammaPi0 = RebinTH1D(cocktailRhoGammaPi0,RebinCocktail);
   // cocktailAllGamma = RebinTH1D(cocktailAllGamma,RebinCocktail);
   // cocktailPi0Gamma = RebinTH1D(cocktailPi0Gamma,RebinCocktail);
   // cocktailEtaGamma = RebinTH1D(cocktailEtaGamma,RebinCocktail);
   // cocktailEtapGamma = RebinTH1D(cocktailEtapGamma,RebinCocktail);
   // cocktailOmegaGamma = RebinTH1D(cocktailOmegaGamma,RebinCocktail);
   // cocktailPhiGamma = RebinTH1D(cocktailPhiGamma,RebinCocktail);
   // cocktailRhoGamma = RebinTH1D(cocktailRhoGamma,RebinCocktail);
   // cocktailPion = RebinTH1D(cocktailPion,RebinCocktail);




   TH1D *cocktailGammaPionsAllGamma = (TH1D*)cocktailPi0Gamma->Clone("cocktailGammaPionsAllGamma");
   cocktailGammaPionsAllGamma->Divide(cocktailAllGamma,cocktailGammaPionsAllGamma,1,1,"B");
   TH1D *cocktailGammaEtaAllGamma = (TH1D*)cocktailPi0Gamma->Clone("cocktailGammaEtaAllGamma");
   cocktailGammaEtaAllGamma->Divide(cocktailEtaGamma,cocktailGammaEtaAllGamma,1,1,"B");
   TH1D *cocktailGammaOmegaAllGamma = (TH1D*)cocktailPi0Gamma->Clone("cocktailGammaOmegaAllGamma");
   cocktailGammaOmegaAllGamma->Divide(cocktailOmegaGamma,cocktailGammaOmegaAllGamma,1,1,"B");
   TH1D *cocktailGammaSigmaAllGamma = (TH1D*)cocktailPi0Gamma->Clone("cocktailGammaSigmaAllGamma");
   cocktailGammaSigmaAllGamma->Divide(cocktailSigmaGamma,cocktailGammaSigmaAllGamma,1,1,"B");
   TH1D *cocktailGammaEtapAllGamma = (TH1D*)cocktailPi0Gamma->Clone("cocktailGammaEtapAllGamma");
   cocktailGammaEtapAllGamma->Divide(cocktailEtapGamma,cocktailGammaEtapAllGamma,1,1,"B");
   TH1D *cocktailGammaPhiAllGamma = (TH1D*)cocktailPi0Gamma->Clone("cocktailGammaPhiAllGamma");
   cocktailGammaPhiAllGamma->Divide(cocktailPhiGamma,cocktailGammaPhiAllGamma,1,1,"B");


   // return;
   // cocktailAllGamma = RebinTH1D(cocktailAllGamma,cocktailPion);
   // cocktailAllGammaPi0->Divide(cocktailAllGamma,cocktailPion);
   
   // // cocktailAllGammaPi0->Draw();
   // // cocktailPion->Draw("same");
   
   // return;
   // cocktailPi0GammaPi0->Divide(cocktailPi0Gamma,cocktailPion);
   // cocktailEtaGammaPi0->Divide(cocktailEtaGamma,cocktailPion);
   // cocktailOmegaGammaPi0->Divide(cocktailOmegaGamma,cocktailPion);
   // cocktailEtapGammaPi0->Divide(cocktailEtapGamma,cocktailPion);
   // cocktailPhiGammaPi0->Divide(cocktailPhiGamma,cocktailPion);
   // cocktailRhoGammaPi0->Divide(cocktailRhoGamma,cocktailPion);



   TGraphErrors *DirectPhotonDoubleNLOhalf =(TGraphErrors*) fileInput->Get("doubleRatioNLOhalf");
   TGraphErrors *DirectPhotonDoubleNLOone =(TGraphErrors*) fileInput->Get("doubleRatioNLOone");
   TGraphErrors *DirectPhotonDoubleNLOtwo  =(TGraphErrors*) fileInput->Get("doubleRatioNLOtwo");

   DirectPhotonDoubleNLOhalf->RemovePoint(0);
   DirectPhotonDoubleNLOone->RemovePoint(0);
   DirectPhotonDoubleNLOtwo->RemovePoint(0);

   TGraphErrors *DirectPhotonNLOhalf =(TGraphErrors*) fileInput->Get("graphNLOCalcMuHalf");
   TGraphErrors *DirectPhotonNLOone =(TGraphErrors*) fileInput->Get("graphNLOCalcMuOne");
   TGraphErrors *DirectPhotonNLOtwo  =(TGraphErrors*) fileInput->Get("graphNLOCalcMuTwo");

   DirectPhotonNLOhalf->RemovePoint(0);
   DirectPhotonNLOone->RemovePoint(0);
   DirectPhotonNLOtwo->RemovePoint(0);


   Double_t *errorup = new Double_t[DirectPhotonDoubleNLOhalf->GetN()];
   Double_t *errorlow = new Double_t[DirectPhotonDoubleNLOtwo->GetN()];

   Double_t *errorSpecup = new Double_t[DirectPhotonDoubleNLOhalf->GetN()];
   Double_t *errorSpeclow = new Double_t[DirectPhotonDoubleNLOtwo->GetN()];


   Double_t *yHalf = DirectPhotonDoubleNLOhalf->GetY();
   Double_t *yOne = DirectPhotonDoubleNLOone->GetY();
   Double_t *yTwo = DirectPhotonDoubleNLOtwo->GetY();

   Double_t *ySpecHalf = DirectPhotonNLOhalf->GetY();
   Double_t *ySpecOne = DirectPhotonNLOone->GetY();
   Double_t *ySpecTwo = DirectPhotonNLOtwo->GetY();

   for(Int_t i = 0;i<DirectPhotonDoubleNLOhalf->GetN(); i++){
      errorup[i] = yHalf[i]-yOne[i];
      errorlow[i] = -yTwo[i]+yOne[i];
      errorSpecup[i] = ySpecHalf[i]-ySpecOne[i];
      errorSpeclow[i] = -ySpecTwo[i]+ySpecOne[i];
   }

   Int_t reduceBins = 0;
   if(GetCentralityStringA(cutSel).CompareTo("0-20%") == 0) reduceBins = 34;
   if(GetCentralityStringA(cutSel).CompareTo("20-40%") == 0) reduceBins = 34;
   if(GetCentralityStringA(cutSel).CompareTo("40-80%") == 0) reduceBins = 34;
   if(GetCentralityStringA(cutSel).CompareTo("0-40%") == 0) reduceBins = 34;
   if(GetCentralityStringA(cutSel).CompareTo("0-80%") == 0) reduceBins = 34;
   reduceBins = 40;

   TGraphAsymmErrors *NLODoubleRatio = new TGraphAsymmErrors( DirectPhotonDoubleNLOhalf->GetN()-reduceBins, DirectPhotonDoubleNLOone->GetX(), DirectPhotonDoubleNLOone->GetY(), DirectPhotonDoubleNLOone->GetEX(), DirectPhotonDoubleNLOone->GetEX(), errorlow,errorup );
   TGraphAsymmErrors *NLO = new TGraphAsymmErrors( DirectPhotonDoubleNLOhalf->GetN()-reduceBins, DirectPhotonNLOone->GetX(), DirectPhotonNLOone->GetY(), DirectPhotonNLOone->GetEX(), DirectPhotonNLOone->GetEX(),errorSpeclow ,errorSpecup );


   NLODoubleRatio->SetLineColor(kBlue-2);
   NLODoubleRatio->SetFillColor(kBlue-2);
   NLODoubleRatio->SetLineWidth(3.0);
   NLODoubleRatio->SetMarkerSize(0);
   NLO->SetLineColor(kBlue-2);
   NLO->SetFillColor(kBlue-2);
   NLO->SetLineWidth(3.0);
   NLO->SetMarkerSize(0);
   Int_t nBinsDoubleRatio = DoubleRatioPi0Fit->GetNbinsX();
   Int_t nBinsIncRatio = IncRatio->GetNbinsX();
   Int_t nBinsPi0 = Pi0->GetNbinsX();
   Int_t nBinsGamma = Gamma->GetNbinsX();
   Double_t* ptBinsDoubleRatio;
   Double_t* ptBinsDoubleRatioErr;
   Double_t* ptBinsIncRatio;
   Double_t* ptBinsIncRatioErr;
   Double_t* ptBinsPi0;
   Double_t* ptBinsPi0Err;
   Double_t* ptBinsGamma;
   Double_t* ptBinsGammaErr;

   TGraphAsymmErrors  **graphPosErrorsDoubleRatio = new TGraphAsymmErrors*[ErrorDoubleRatio];
   TGraphAsymmErrors  **graphNegErrorsDoubleRatio = new TGraphAsymmErrors*[ErrorDoubleRatio];
   TGraphAsymmErrors  **graphPosErrorsDoubleRatioA = new TGraphAsymmErrors*[ErrorDoubleRatio];
   TGraphAsymmErrors  **graphNegErrorsDoubleRatioA = new TGraphAsymmErrors*[ErrorDoubleRatio];
   TGraphAsymmErrors  **graphPosErrorsDoubleRatioB = new TGraphAsymmErrors*[ErrorDoubleRatio];
   TGraphAsymmErrors  **graphNegErrorsDoubleRatioB = new TGraphAsymmErrors*[ErrorDoubleRatio];
   TGraphAsymmErrors  **graphPosErrorsDoubleRatioC = new TGraphAsymmErrors*[ErrorDoubleRatio];
   TGraphAsymmErrors  **graphNegErrorsDoubleRatioC = new TGraphAsymmErrors*[ErrorDoubleRatio];

   TGraphAsymmErrors  **graphPosErrorsDoubleRatioPi0Fit = new TGraphAsymmErrors*[ErrorDoubleRatio];
   TGraphAsymmErrors  **graphNegErrorsDoubleRatioPi0Fit = new TGraphAsymmErrors*[ErrorDoubleRatio];
   TGraphAsymmErrors  **graphPosErrorsDoubleRatioPi0FitA = new TGraphAsymmErrors*[ErrorDoubleRatio];
   TGraphAsymmErrors  **graphNegErrorsDoubleRatioPi0FitA = new TGraphAsymmErrors*[ErrorDoubleRatio];
   TGraphAsymmErrors  **graphPosErrorsDoubleRatioPi0FitB = new TGraphAsymmErrors*[ErrorDoubleRatio];
   TGraphAsymmErrors  **graphNegErrorsDoubleRatioPi0FitB = new TGraphAsymmErrors*[ErrorDoubleRatio];
   TGraphAsymmErrors  **graphPosErrorsDoubleRatioPi0FitC = new TGraphAsymmErrors*[ErrorDoubleRatio];
   TGraphAsymmErrors  **graphNegErrorsDoubleRatioPi0FitC = new TGraphAsymmErrors*[ErrorDoubleRatio];

   TGraphAsymmErrors  **graphPosErrorsIncRatio = new TGraphAsymmErrors*[ErrorIncRatio];
   TGraphAsymmErrors  **graphNegErrorsIncRatio = new TGraphAsymmErrors*[ErrorIncRatio];
   TGraphAsymmErrors  **graphPosErrorsIncRatioA = new TGraphAsymmErrors*[ErrorIncRatio];
   TGraphAsymmErrors  **graphNegErrorsIncRatioA = new TGraphAsymmErrors*[ErrorIncRatio];
   TGraphAsymmErrors  **graphPosErrorsIncRatioB = new TGraphAsymmErrors*[ErrorIncRatio];
   TGraphAsymmErrors  **graphNegErrorsIncRatioB = new TGraphAsymmErrors*[ErrorIncRatio];
   TGraphAsymmErrors  **graphPosErrorsIncRatioC = new TGraphAsymmErrors*[ErrorIncRatio];
   TGraphAsymmErrors  **graphNegErrorsIncRatioC = new TGraphAsymmErrors*[ErrorIncRatio];

   TGraphAsymmErrors  **graphPosErrorsIncRatioPi0Fit = new TGraphAsymmErrors*[ErrorIncRatio];
   TGraphAsymmErrors  **graphNegErrorsIncRatioPi0Fit = new TGraphAsymmErrors*[ErrorIncRatio];
   TGraphAsymmErrors  **graphPosErrorsIncRatioPi0FitA = new TGraphAsymmErrors*[ErrorIncRatio];
   TGraphAsymmErrors  **graphNegErrorsIncRatioPi0FitA = new TGraphAsymmErrors*[ErrorIncRatio];
   TGraphAsymmErrors  **graphPosErrorsIncRatioPi0FitB = new TGraphAsymmErrors*[ErrorIncRatio];
   TGraphAsymmErrors  **graphNegErrorsIncRatioPi0FitB = new TGraphAsymmErrors*[ErrorIncRatio];
   TGraphAsymmErrors  **graphPosErrorsIncRatioPi0FitC = new TGraphAsymmErrors*[ErrorIncRatio];
   TGraphAsymmErrors  **graphNegErrorsIncRatioPi0FitC = new TGraphAsymmErrors*[ErrorIncRatio];

   TGraphAsymmErrors  **graphPosErrorsPi0 = new TGraphAsymmErrors*[ErrorPi0];
   TGraphAsymmErrors  **graphNegErrorsPi0 = new TGraphAsymmErrors*[ErrorPi0];

   TGraphAsymmErrors  **graphPosErrorsPi0Fit = new TGraphAsymmErrors*[ErrorPi0];
   TGraphAsymmErrors  **graphNegErrorsPi0Fit = new TGraphAsymmErrors*[ErrorPi0];

   TGraphAsymmErrors  **graphPosErrorsGamma = new TGraphAsymmErrors*[ErrorGammaSpec];
   TGraphAsymmErrors  **graphNegErrorsGamma = new TGraphAsymmErrors*[ErrorGammaSpec];
   TGraphAsymmErrors  **graphPosErrorsGammaA = new TGraphAsymmErrors*[ErrorGammaSpec];
   TGraphAsymmErrors  **graphNegErrorsGammaA = new TGraphAsymmErrors*[ErrorGammaSpec];
   TGraphAsymmErrors  **graphPosErrorsGammaB = new TGraphAsymmErrors*[ErrorGammaSpec];
   TGraphAsymmErrors  **graphNegErrorsGammaB = new TGraphAsymmErrors*[ErrorGammaSpec];
   TGraphAsymmErrors  **graphPosErrorsGammaC = new TGraphAsymmErrors*[ErrorGammaSpec];
   TGraphAsymmErrors  **graphNegErrorsGammaC = new TGraphAsymmErrors*[ErrorGammaSpec];

   const Int_t ConstnBinsDoubleRatio = nBinsDoubleRatio;
   const Int_t ConstnBinsIncRatio = nBinsIncRatio;
   const Int_t ConstnBinsPi0 = nBinsPi0;
   const Int_t ConstnBinsGamma = nBinsGamma;

   Double_t* errorsNegDoubleRatio[ErrorDoubleRatio];
   Double_t errorsNegCorrDoubleRatio[ErrorDoubleRatio][ConstnBinsDoubleRatio];
   Double_t errorsNegSummedDoubleRatio[ConstnBinsDoubleRatio];
   Double_t errorsNegCorrSummedDoubleRatio[ConstnBinsDoubleRatio];
   Double_t errorsNegCorrMatSummedDoubleRatio[ConstnBinsDoubleRatio];

   Double_t* errorsNegErrDoubleRatio[ErrorDoubleRatio];
   Double_t errorsNegErrCorrDoubleRatio[ErrorDoubleRatio][ConstnBinsDoubleRatio];
   Double_t errorsNegErrSummedDoubleRatio[ConstnBinsDoubleRatio];
   Double_t errorsNegErrCorrSummedDoubleRatio[ConstnBinsDoubleRatio];

   Double_t* errorsPosDoubleRatio[ErrorDoubleRatio];
   Double_t errorsPosCorrDoubleRatio[ErrorDoubleRatio][ConstnBinsDoubleRatio];
   Double_t errorsPosSummedDoubleRatio[ConstnBinsDoubleRatio];
   Double_t errorsPosCorrSummedDoubleRatio[ConstnBinsDoubleRatio];
   Double_t errorsPosCorrMatSummedDoubleRatio[ConstnBinsDoubleRatio];

   Double_t* errorsPosErrDoubleRatio[ErrorDoubleRatio];
   Double_t errorsPosErrSummedDoubleRatio[ConstnBinsDoubleRatio];
   Double_t errorsPosErrCorrDoubleRatio[ErrorDoubleRatio][ConstnBinsDoubleRatio];
   Double_t errorsPosErrCorrSummedDoubleRatio[ConstnBinsDoubleRatio];

   Double_t* errorsNegDoubleRatioPi0Fit[ErrorDoubleRatio];
   Double_t errorsNegCorrDoubleRatioPi0Fit[ErrorDoubleRatio][ConstnBinsDoubleRatio];
   Double_t errorsNegSummedDoubleRatioPi0Fit[ConstnBinsDoubleRatio];
   Double_t errorsNegCorrSummedDoubleRatioPi0Fit[ConstnBinsDoubleRatio];
   Double_t errorsNegCorrMatSummedDoubleRatioPi0Fit[ConstnBinsDoubleRatio];

   Double_t* errorsNegErrDoubleRatioPi0Fit[ErrorDoubleRatio];
   Double_t errorsNegErrCorrDoubleRatioPi0Fit[ErrorDoubleRatio][ConstnBinsDoubleRatio];
   Double_t errorsNegErrSummedDoubleRatioPi0Fit[ConstnBinsDoubleRatio];
   Double_t errorsNegErrCorrSummedDoubleRatioPi0Fit[ConstnBinsDoubleRatio];

   Double_t* errorsPosDoubleRatioPi0Fit[ErrorDoubleRatio];
   Double_t errorsPosCorrDoubleRatioPi0Fit[ErrorDoubleRatio][ConstnBinsDoubleRatio];
   Double_t errorsPosSummedDoubleRatioPi0Fit[ConstnBinsDoubleRatio];
   Double_t errorsPosCorrSummedDoubleRatioPi0Fit[ConstnBinsDoubleRatio];
   Double_t errorsPosCorrMatSummedDoubleRatioPi0Fit[ConstnBinsDoubleRatio];

   Double_t* errorsPosErrDoubleRatioPi0Fit[ErrorDoubleRatio];
   Double_t errorsPosErrSummedDoubleRatioPi0Fit[ConstnBinsDoubleRatio];
   Double_t errorsPosErrCorrDoubleRatioPi0Fit[ErrorDoubleRatio][ConstnBinsDoubleRatio];
   Double_t errorsPosErrCorrSummedDoubleRatioPi0Fit[ConstnBinsDoubleRatio];

   Double_t* errorsNegDoubleRatioA[ErrorDoubleRatio];
   Double_t errorsNegCorrDoubleRatioA[ErrorDoubleRatio][ConstnBinsDoubleRatio];
   Double_t errorsNegSummedDoubleRatioA[ConstnBinsDoubleRatio];
   Double_t errorsNegCorrSummedDoubleRatioA[ConstnBinsDoubleRatio];
   Double_t errorsNegCorrMatSummedDoubleRatioA[ConstnBinsDoubleRatio];

   Double_t* errorsNegErrDoubleRatioA[ErrorDoubleRatio];
   Double_t errorsNegErrCorrDoubleRatioA[ErrorDoubleRatio][ConstnBinsDoubleRatio];
   Double_t errorsNegErrSummedDoubleRatioA[ConstnBinsDoubleRatio];
   Double_t errorsNegErrCorrSummedDoubleRatioA[ConstnBinsDoubleRatio];

   Double_t* errorsPosDoubleRatioA[ErrorDoubleRatio];
   Double_t errorsPosCorrDoubleRatioA[ErrorDoubleRatio][ConstnBinsDoubleRatio];
   Double_t errorsPosSummedDoubleRatioA[ConstnBinsDoubleRatio];
   Double_t errorsPosCorrSummedDoubleRatioA[ConstnBinsDoubleRatio];
   Double_t errorsPosCorrMatSummedDoubleRatioA[ConstnBinsDoubleRatio];

   Double_t* errorsPosErrDoubleRatioA[ErrorDoubleRatio];
   Double_t errorsPosErrSummedDoubleRatioA[ConstnBinsDoubleRatio];
   Double_t errorsPosErrCorrDoubleRatioA[ErrorDoubleRatio][ConstnBinsDoubleRatio];
   Double_t errorsPosErrCorrSummedDoubleRatioA[ConstnBinsDoubleRatio];

   Double_t* errorsNegDoubleRatioPi0FitA[ErrorDoubleRatio];
   Double_t errorsNegCorrDoubleRatioPi0FitA[ErrorDoubleRatio][ConstnBinsDoubleRatio];
   Double_t errorsNegSummedDoubleRatioPi0FitA[ConstnBinsDoubleRatio];
   Double_t errorsNegCorrSummedDoubleRatioPi0FitA[ConstnBinsDoubleRatio];
   Double_t errorsNegCorrMatSummedDoubleRatioPi0FitA[ConstnBinsDoubleRatio];

   Double_t* errorsNegErrDoubleRatioPi0FitA[ErrorDoubleRatio];
   Double_t errorsNegErrCorrDoubleRatioPi0FitA[ErrorDoubleRatio][ConstnBinsDoubleRatio];
   Double_t errorsNegErrSummedDoubleRatioPi0FitA[ConstnBinsDoubleRatio];
   Double_t errorsNegErrCorrSummedDoubleRatioPi0FitA[ConstnBinsDoubleRatio];

   Double_t* errorsPosDoubleRatioPi0FitA[ErrorDoubleRatio];
   Double_t errorsPosCorrDoubleRatioPi0FitA[ErrorDoubleRatio][ConstnBinsDoubleRatio];
   Double_t errorsPosSummedDoubleRatioPi0FitA[ConstnBinsDoubleRatio];
   Double_t errorsPosCorrSummedDoubleRatioPi0FitA[ConstnBinsDoubleRatio];
   Double_t errorsPosCorrMatSummedDoubleRatioPi0FitA[ConstnBinsDoubleRatio];

   Double_t* errorsPosErrDoubleRatioPi0FitA[ErrorDoubleRatio];
   Double_t errorsPosErrSummedDoubleRatioPi0FitA[ConstnBinsDoubleRatio];
   Double_t errorsPosErrCorrDoubleRatioPi0FitA[ErrorDoubleRatio][ConstnBinsDoubleRatio];
   Double_t errorsPosErrCorrSummedDoubleRatioPi0FitA[ConstnBinsDoubleRatio];

   Double_t* errorsNegDoubleRatioB[ErrorDoubleRatio];
   Double_t errorsNegCorrDoubleRatioB[ErrorDoubleRatio][ConstnBinsDoubleRatio];
   Double_t errorsNegSummedDoubleRatioB[ConstnBinsDoubleRatio];
   Double_t errorsNegCorrSummedDoubleRatioB[ConstnBinsDoubleRatio];
   Double_t errorsNegCorrMatSummedDoubleRatioB[ConstnBinsDoubleRatio];

   Double_t* errorsNegErrDoubleRatioB[ErrorDoubleRatio];
   Double_t errorsNegErrCorrDoubleRatioB[ErrorDoubleRatio][ConstnBinsDoubleRatio];
   Double_t errorsNegErrSummedDoubleRatioB[ConstnBinsDoubleRatio];
   Double_t errorsNegErrCorrSummedDoubleRatioB[ConstnBinsDoubleRatio];

   Double_t* errorsPosDoubleRatioB[ErrorDoubleRatio];
   Double_t errorsPosCorrDoubleRatioB[ErrorDoubleRatio][ConstnBinsDoubleRatio];
   Double_t errorsPosSummedDoubleRatioB[ConstnBinsDoubleRatio];
   Double_t errorsPosCorrSummedDoubleRatioB[ConstnBinsDoubleRatio];
   Double_t errorsPosCorrMatSummedDoubleRatioB[ConstnBinsDoubleRatio];

   Double_t* errorsPosErrDoubleRatioB[ErrorDoubleRatio];
   Double_t errorsPosErrSummedDoubleRatioB[ConstnBinsDoubleRatio];
   Double_t errorsPosErrCorrDoubleRatioB[ErrorDoubleRatio][ConstnBinsDoubleRatio];
   Double_t errorsPosErrCorrSummedDoubleRatioB[ConstnBinsDoubleRatio];

   Double_t* errorsNegDoubleRatioPi0FitB[ErrorDoubleRatio];
   Double_t errorsNegCorrDoubleRatioPi0FitB[ErrorDoubleRatio][ConstnBinsDoubleRatio];
   Double_t errorsNegSummedDoubleRatioPi0FitB[ConstnBinsDoubleRatio];
   Double_t errorsNegCorrSummedDoubleRatioPi0FitB[ConstnBinsDoubleRatio];
   Double_t errorsNegCorrMatSummedDoubleRatioPi0FitB[ConstnBinsDoubleRatio];

   Double_t* errorsNegErrDoubleRatioPi0FitB[ErrorDoubleRatio];
   Double_t errorsNegErrCorrDoubleRatioPi0FitB[ErrorDoubleRatio][ConstnBinsDoubleRatio];
   Double_t errorsNegErrSummedDoubleRatioPi0FitB[ConstnBinsDoubleRatio];
   Double_t errorsNegErrCorrSummedDoubleRatioPi0FitB[ConstnBinsDoubleRatio];

   Double_t* errorsPosDoubleRatioPi0FitB[ErrorDoubleRatio];
   Double_t errorsPosCorrDoubleRatioPi0FitB[ErrorDoubleRatio][ConstnBinsDoubleRatio];
   Double_t errorsPosSummedDoubleRatioPi0FitB[ConstnBinsDoubleRatio];
   Double_t errorsPosCorrSummedDoubleRatioPi0FitB[ConstnBinsDoubleRatio];
   Double_t errorsPosCorrMatSummedDoubleRatioPi0FitB[ConstnBinsDoubleRatio];

   Double_t* errorsPosErrDoubleRatioPi0FitB[ErrorDoubleRatio];
   Double_t errorsPosErrSummedDoubleRatioPi0FitB[ConstnBinsDoubleRatio];
   Double_t errorsPosErrCorrDoubleRatioPi0FitB[ErrorDoubleRatio][ConstnBinsDoubleRatio];
   Double_t errorsPosErrCorrSummedDoubleRatioPi0FitB[ConstnBinsDoubleRatio];

   Double_t* errorsNegDoubleRatioC[ErrorDoubleRatio];
   Double_t errorsNegCorrDoubleRatioC[ErrorDoubleRatio][ConstnBinsDoubleRatio];
   Double_t errorsNegSummedDoubleRatioC[ConstnBinsDoubleRatio];
   Double_t errorsNegCorrSummedDoubleRatioC[ConstnBinsDoubleRatio];
   Double_t errorsNegCorrMatSummedDoubleRatioC[ConstnBinsDoubleRatio];

   Double_t* errorsNegErrDoubleRatioC[ErrorDoubleRatio];
   Double_t errorsNegErrCorrDoubleRatioC[ErrorDoubleRatio][ConstnBinsDoubleRatio];
   Double_t errorsNegErrSummedDoubleRatioC[ConstnBinsDoubleRatio];
   Double_t errorsNegErrCorrSummedDoubleRatioC[ConstnBinsDoubleRatio];

   Double_t* errorsPosDoubleRatioC[ErrorDoubleRatio];
   Double_t errorsPosCorrDoubleRatioC[ErrorDoubleRatio][ConstnBinsDoubleRatio];
   Double_t errorsPosSummedDoubleRatioC[ConstnBinsDoubleRatio];
   Double_t errorsPosCorrSummedDoubleRatioC[ConstnBinsDoubleRatio];
   Double_t errorsPosCorrMatSummedDoubleRatioC[ConstnBinsDoubleRatio];

   Double_t* errorsPosErrDoubleRatioC[ErrorDoubleRatio];
   Double_t errorsPosErrSummedDoubleRatioC[ConstnBinsDoubleRatio];
   Double_t errorsPosErrCorrDoubleRatioC[ErrorDoubleRatio][ConstnBinsDoubleRatio];
   Double_t errorsPosErrCorrSummedDoubleRatioC[ConstnBinsDoubleRatio];

   Double_t* errorsNegDoubleRatioPi0FitC[ErrorDoubleRatio];
   Double_t errorsNegCorrDoubleRatioPi0FitC[ErrorDoubleRatio][ConstnBinsDoubleRatio];
   Double_t errorsNegSummedDoubleRatioPi0FitC[ConstnBinsDoubleRatio];
   Double_t errorsNegCorrSummedDoubleRatioPi0FitC[ConstnBinsDoubleRatio];
   Double_t errorsNegCorrMatSummedDoubleRatioPi0FitC[ConstnBinsDoubleRatio];

   Double_t* errorsNegErrDoubleRatioPi0FitC[ErrorDoubleRatio];
   Double_t errorsNegErrCorrDoubleRatioPi0FitC[ErrorDoubleRatio][ConstnBinsDoubleRatio];
   Double_t errorsNegErrSummedDoubleRatioPi0FitC[ConstnBinsDoubleRatio];
   Double_t errorsNegErrCorrSummedDoubleRatioPi0FitC[ConstnBinsDoubleRatio];

   Double_t* errorsPosDoubleRatioPi0FitC[ErrorDoubleRatio];
   Double_t errorsPosCorrDoubleRatioPi0FitC[ErrorDoubleRatio][ConstnBinsDoubleRatio];
   Double_t errorsPosSummedDoubleRatioPi0FitC[ConstnBinsDoubleRatio];
   Double_t errorsPosCorrSummedDoubleRatioPi0FitC[ConstnBinsDoubleRatio];
   Double_t errorsPosCorrMatSummedDoubleRatioPi0FitC[ConstnBinsDoubleRatio];

   Double_t* errorsPosErrDoubleRatioPi0FitC[ErrorDoubleRatio];
   Double_t errorsPosErrSummedDoubleRatioPi0FitC[ConstnBinsDoubleRatio];
   Double_t errorsPosErrCorrDoubleRatioPi0FitC[ErrorDoubleRatio][ConstnBinsDoubleRatio];
   Double_t errorsPosErrCorrSummedDoubleRatioPi0FitC[ConstnBinsDoubleRatio];

   Double_t* errorsNegIncRatio[ErrorIncRatio];
   Double_t errorsNegCorrIncRatio[ErrorIncRatio][ConstnBinsIncRatio];
   Double_t errorsNegSummedIncRatio[ConstnBinsIncRatio];
   Double_t errorsNegCorrSummedIncRatio[ConstnBinsIncRatio];
   Double_t errorsNegCorrMatSummedIncRatio[ConstnBinsIncRatio];

   Double_t* errorsNegErrIncRatio[ErrorIncRatio];
   Double_t errorsNegErrCorrIncRatio[ErrorIncRatio][ConstnBinsIncRatio];
   Double_t errorsNegErrSummedIncRatio[ConstnBinsIncRatio];
   Double_t errorsNegErrCorrSummedIncRatio[ConstnBinsIncRatio];

   Double_t* errorsPosIncRatio[ErrorIncRatio];
   Double_t errorsPosCorrIncRatio[ErrorIncRatio][ConstnBinsIncRatio];
   Double_t errorsPosSummedIncRatio[ConstnBinsIncRatio];
   Double_t errorsPosCorrSummedIncRatio[ConstnBinsIncRatio];
   Double_t errorsPosCorrMatSummedIncRatio[ConstnBinsIncRatio];

   Double_t* errorsPosErrIncRatio[ErrorIncRatio];
   Double_t errorsPosErrSummedIncRatio[ConstnBinsIncRatio];
   Double_t errorsPosErrCorrIncRatio[ErrorIncRatio][ConstnBinsIncRatio];
   Double_t errorsPosErrCorrSummedIncRatio[ConstnBinsIncRatio];

   Double_t* errorsNegIncRatioPi0Fit[ErrorIncRatio];
   Double_t errorsNegCorrIncRatioPi0Fit[ErrorIncRatio][ConstnBinsIncRatio];
   Double_t errorsNegSummedIncRatioPi0Fit[ConstnBinsIncRatio];
   Double_t errorsNegCorrSummedIncRatioPi0Fit[ConstnBinsIncRatio];
   Double_t errorsNegCorrMatSummedIncRatioPi0Fit[ConstnBinsIncRatio];

   Double_t* errorsNegErrIncRatioPi0Fit[ErrorIncRatio];
   Double_t errorsNegErrCorrIncRatioPi0Fit[ErrorIncRatio][ConstnBinsIncRatio];
   Double_t errorsNegErrSummedIncRatioPi0Fit[ConstnBinsIncRatio];
   Double_t errorsNegErrCorrSummedIncRatioPi0Fit[ConstnBinsIncRatio];

   Double_t* errorsPosIncRatioPi0Fit[ErrorIncRatio];
   Double_t errorsPosCorrIncRatioPi0Fit[ErrorIncRatio][ConstnBinsIncRatio];
   Double_t errorsPosSummedIncRatioPi0Fit[ConstnBinsIncRatio];
   Double_t errorsPosCorrSummedIncRatioPi0Fit[ConstnBinsIncRatio];
   Double_t errorsPosCorrMatSummedIncRatioPi0Fit[ConstnBinsIncRatio];

   Double_t* errorsPosErrIncRatioPi0Fit[ErrorIncRatio];
   Double_t errorsPosErrSummedIncRatioPi0Fit[ConstnBinsIncRatio];
   Double_t errorsPosErrCorrIncRatioPi0Fit[ErrorIncRatio][ConstnBinsIncRatio];
   Double_t errorsPosErrCorrSummedIncRatioPi0Fit[ConstnBinsIncRatio];

   Double_t* errorsNegIncRatioA[ErrorIncRatio];
   Double_t errorsNegCorrIncRatioA[ErrorIncRatio][ConstnBinsIncRatio];
   Double_t errorsNegSummedIncRatioA[ConstnBinsIncRatio];
   Double_t errorsNegCorrSummedIncRatioA[ConstnBinsIncRatio];
   Double_t errorsNegCorrMatSummedIncRatioA[ConstnBinsIncRatio];

   Double_t* errorsNegErrIncRatioA[ErrorIncRatio];
   Double_t errorsNegErrCorrIncRatioA[ErrorIncRatio][ConstnBinsIncRatio];
   Double_t errorsNegErrSummedIncRatioA[ConstnBinsIncRatio];
   Double_t errorsNegErrCorrSummedIncRatioA[ConstnBinsIncRatio];

   Double_t* errorsPosIncRatioA[ErrorIncRatio];
   Double_t errorsPosCorrIncRatioA[ErrorIncRatio][ConstnBinsIncRatio];
   Double_t errorsPosSummedIncRatioA[ConstnBinsIncRatio];
   Double_t errorsPosCorrSummedIncRatioA[ConstnBinsIncRatio];
   Double_t errorsPosCorrMatSummedIncRatioA[ConstnBinsIncRatio];

   Double_t* errorsPosErrIncRatioA[ErrorIncRatio];
   Double_t errorsPosErrSummedIncRatioA[ConstnBinsIncRatio];
   Double_t errorsPosErrCorrIncRatioA[ErrorIncRatio][ConstnBinsIncRatio];
   Double_t errorsPosErrCorrSummedIncRatioA[ConstnBinsIncRatio];

   Double_t* errorsNegIncRatioPi0FitA[ErrorIncRatio];
   Double_t errorsNegCorrIncRatioPi0FitA[ErrorIncRatio][ConstnBinsIncRatio];
   Double_t errorsNegSummedIncRatioPi0FitA[ConstnBinsIncRatio];
   Double_t errorsNegCorrSummedIncRatioPi0FitA[ConstnBinsIncRatio];
   Double_t errorsNegCorrMatSummedIncRatioPi0FitA[ConstnBinsIncRatio];

   Double_t* errorsNegErrIncRatioPi0FitA[ErrorIncRatio];
   Double_t errorsNegErrCorrIncRatioPi0FitA[ErrorIncRatio][ConstnBinsIncRatio];
   Double_t errorsNegErrSummedIncRatioPi0FitA[ConstnBinsIncRatio];
   Double_t errorsNegErrCorrSummedIncRatioPi0FitA[ConstnBinsIncRatio];

   Double_t* errorsPosIncRatioPi0FitA[ErrorIncRatio];
   Double_t errorsPosCorrIncRatioPi0FitA[ErrorIncRatio][ConstnBinsIncRatio];
   Double_t errorsPosSummedIncRatioPi0FitA[ConstnBinsIncRatio];
   Double_t errorsPosCorrSummedIncRatioPi0FitA[ConstnBinsIncRatio];
   Double_t errorsPosCorrMatSummedIncRatioPi0FitA[ConstnBinsIncRatio];

   Double_t* errorsPosErrIncRatioPi0FitA[ErrorIncRatio];
   Double_t errorsPosErrSummedIncRatioPi0FitA[ConstnBinsIncRatio];
   Double_t errorsPosErrCorrIncRatioPi0FitA[ErrorIncRatio][ConstnBinsIncRatio];
   Double_t errorsPosErrCorrSummedIncRatioPi0FitA[ConstnBinsIncRatio];

   Double_t* errorsNegIncRatioB[ErrorIncRatio];
   Double_t errorsNegCorrIncRatioB[ErrorIncRatio][ConstnBinsIncRatio];
   Double_t errorsNegSummedIncRatioB[ConstnBinsIncRatio];
   Double_t errorsNegCorrSummedIncRatioB[ConstnBinsIncRatio];
   Double_t errorsNegCorrMatSummedIncRatioB[ConstnBinsIncRatio];

   Double_t* errorsNegErrIncRatioB[ErrorIncRatio];
   Double_t errorsNegErrCorrIncRatioB[ErrorIncRatio][ConstnBinsIncRatio];
   Double_t errorsNegErrSummedIncRatioB[ConstnBinsIncRatio];
   Double_t errorsNegErrCorrSummedIncRatioB[ConstnBinsIncRatio];

   Double_t* errorsPosIncRatioB[ErrorIncRatio];
   Double_t errorsPosCorrIncRatioB[ErrorIncRatio][ConstnBinsIncRatio];
   Double_t errorsPosSummedIncRatioB[ConstnBinsIncRatio];
   Double_t errorsPosCorrSummedIncRatioB[ConstnBinsIncRatio];
   Double_t errorsPosCorrMatSummedIncRatioB[ConstnBinsIncRatio];

   Double_t* errorsPosErrIncRatioB[ErrorIncRatio];
   Double_t errorsPosErrSummedIncRatioB[ConstnBinsIncRatio];
   Double_t errorsPosErrCorrIncRatioB[ErrorIncRatio][ConstnBinsIncRatio];
   Double_t errorsPosErrCorrSummedIncRatioB[ConstnBinsIncRatio];

   Double_t* errorsNegIncRatioPi0FitB[ErrorIncRatio];
   Double_t errorsNegCorrIncRatioPi0FitB[ErrorIncRatio][ConstnBinsIncRatio];
   Double_t errorsNegSummedIncRatioPi0FitB[ConstnBinsIncRatio];
   Double_t errorsNegCorrSummedIncRatioPi0FitB[ConstnBinsIncRatio];
   Double_t errorsNegCorrMatSummedIncRatioPi0FitB[ConstnBinsIncRatio];

   Double_t* errorsNegErrIncRatioPi0FitB[ErrorIncRatio];
   Double_t errorsNegErrCorrIncRatioPi0FitB[ErrorIncRatio][ConstnBinsIncRatio];
   Double_t errorsNegErrSummedIncRatioPi0FitB[ConstnBinsIncRatio];
   Double_t errorsNegErrCorrSummedIncRatioPi0FitB[ConstnBinsIncRatio];

   Double_t* errorsPosIncRatioPi0FitB[ErrorIncRatio];
   Double_t errorsPosCorrIncRatioPi0FitB[ErrorIncRatio][ConstnBinsIncRatio];
   Double_t errorsPosSummedIncRatioPi0FitB[ConstnBinsIncRatio];
   Double_t errorsPosCorrSummedIncRatioPi0FitB[ConstnBinsIncRatio];
   Double_t errorsPosCorrMatSummedIncRatioPi0FitB[ConstnBinsIncRatio];

   Double_t* errorsPosErrIncRatioPi0FitB[ErrorIncRatio];
   Double_t errorsPosErrSummedIncRatioPi0FitB[ConstnBinsIncRatio];
   Double_t errorsPosErrCorrIncRatioPi0FitB[ErrorIncRatio][ConstnBinsIncRatio];
   Double_t errorsPosErrCorrSummedIncRatioPi0FitB[ConstnBinsIncRatio];

   Double_t* errorsNegIncRatioC[ErrorIncRatio];
   Double_t errorsNegCorrIncRatioC[ErrorIncRatio][ConstnBinsIncRatio];
   Double_t errorsNegSummedIncRatioC[ConstnBinsIncRatio];
   Double_t errorsNegCorrSummedIncRatioC[ConstnBinsIncRatio];
   Double_t errorsNegCorrMatSummedIncRatioC[ConstnBinsIncRatio];

   Double_t* errorsNegErrIncRatioC[ErrorIncRatio];
   Double_t errorsNegErrCorrIncRatioC[ErrorIncRatio][ConstnBinsIncRatio];
   Double_t errorsNegErrSummedIncRatioC[ConstnBinsIncRatio];
   Double_t errorsNegErrCorrSummedIncRatioC[ConstnBinsIncRatio];

   Double_t* errorsPosIncRatioC[ErrorIncRatio];
   Double_t errorsPosCorrIncRatioC[ErrorIncRatio][ConstnBinsIncRatio];
   Double_t errorsPosSummedIncRatioC[ConstnBinsIncRatio];
   Double_t errorsPosCorrSummedIncRatioC[ConstnBinsIncRatio];
   Double_t errorsPosCorrMatSummedIncRatioC[ConstnBinsIncRatio];

   Double_t* errorsPosErrIncRatioC[ErrorIncRatio];
   Double_t errorsPosErrSummedIncRatioC[ConstnBinsIncRatio];
   Double_t errorsPosErrCorrIncRatioC[ErrorIncRatio][ConstnBinsIncRatio];
   Double_t errorsPosErrCorrSummedIncRatioC[ConstnBinsIncRatio];

   Double_t* errorsNegIncRatioPi0FitC[ErrorIncRatio];
   Double_t errorsNegCorrIncRatioPi0FitC[ErrorIncRatio][ConstnBinsIncRatio];
   Double_t errorsNegSummedIncRatioPi0FitC[ConstnBinsIncRatio];
   Double_t errorsNegCorrSummedIncRatioPi0FitC[ConstnBinsIncRatio];
   Double_t errorsNegCorrMatSummedIncRatioPi0FitC[ConstnBinsIncRatio];

   Double_t* errorsNegErrIncRatioPi0FitC[ErrorIncRatio];
   Double_t errorsNegErrCorrIncRatioPi0FitC[ErrorIncRatio][ConstnBinsIncRatio];
   Double_t errorsNegErrSummedIncRatioPi0FitC[ConstnBinsIncRatio];
   Double_t errorsNegErrCorrSummedIncRatioPi0FitC[ConstnBinsIncRatio];

   Double_t* errorsPosIncRatioPi0FitC[ErrorIncRatio];
   Double_t errorsPosCorrIncRatioPi0FitC[ErrorIncRatio][ConstnBinsIncRatio];
   Double_t errorsPosSummedIncRatioPi0FitC[ConstnBinsIncRatio];
   Double_t errorsPosCorrSummedIncRatioPi0FitC[ConstnBinsIncRatio];
   Double_t errorsPosCorrMatSummedIncRatioPi0FitC[ConstnBinsIncRatio];

   Double_t* errorsPosErrIncRatioPi0FitC[ErrorIncRatio];
   Double_t errorsPosErrSummedIncRatioPi0FitC[ConstnBinsIncRatio];
   Double_t errorsPosErrCorrIncRatioPi0FitC[ErrorIncRatio][ConstnBinsIncRatio];
   Double_t errorsPosErrCorrSummedIncRatioPi0FitC[ConstnBinsIncRatio];

   Double_t* errorsNegPi0[ErrorPi0];
   Double_t errorsNegCorrPi0[ErrorPi0][ConstnBinsPi0];
   Double_t errorsNegSummedPi0[ConstnBinsPi0];
   Double_t errorsNegCorrSummedPi0[ConstnBinsPi0];
   Double_t errorsNegCorrMatSummedPi0[ConstnBinsPi0];

   Double_t* errorsNegErrPi0[ErrorPi0];
   Double_t errorsNegErrCorrPi0[ErrorPi0][ConstnBinsPi0];
   Double_t errorsNegErrSummedPi0[ConstnBinsPi0];
   Double_t errorsNegErrCorrSummedPi0[ConstnBinsPi0];

   Double_t* errorsPosPi0[ErrorPi0];
   Double_t errorsPosCorrPi0[ErrorPi0][ConstnBinsPi0];
   Double_t errorsPosSummedPi0[ConstnBinsPi0];
   Double_t errorsPosCorrSummedPi0[ConstnBinsPi0];
   Double_t errorsPosCorrMatSummedPi0[ConstnBinsPi0];

   Double_t* errorsPosErrPi0[ErrorPi0];
   Double_t errorsPosErrSummedPi0[ConstnBinsPi0];
   Double_t errorsPosErrCorrPi0[ErrorPi0][ConstnBinsPi0];
   Double_t errorsPosErrCorrSummedPi0[ConstnBinsPi0];

   Double_t* errorsNegPi0Fit[ErrorPi0];
   Double_t errorsNegCorrPi0Fit[ErrorPi0][ConstnBinsPi0];
   Double_t errorsNegSummedPi0Fit[ConstnBinsPi0];
   Double_t errorsNegCorrSummedPi0Fit[ConstnBinsPi0];
   Double_t errorsNegCorrMatSummedPi0Fit[ConstnBinsPi0];

   Double_t* errorsNegErrPi0Fit[ErrorPi0];
   Double_t errorsNegErrCorrPi0Fit[ErrorPi0][ConstnBinsPi0];
   Double_t errorsNegErrSummedPi0Fit[ConstnBinsPi0];
   Double_t errorsNegErrCorrSummedPi0Fit[ConstnBinsPi0];

   Double_t* errorsPosPi0Fit[ErrorPi0];
   Double_t errorsPosCorrPi0Fit[ErrorPi0][ConstnBinsPi0];
   Double_t errorsPosSummedPi0Fit[ConstnBinsPi0];
   Double_t errorsPosCorrSummedPi0Fit[ConstnBinsPi0];
   Double_t errorsPosCorrMatSummedPi0Fit[ConstnBinsPi0];

   Double_t* errorsPosErrPi0Fit[ErrorPi0];
   Double_t errorsPosErrSummedPi0Fit[ConstnBinsPi0];
   Double_t errorsPosErrCorrPi0Fit[ErrorPi0][ConstnBinsPi0];
   Double_t errorsPosErrCorrSummedPi0Fit[ConstnBinsPi0];


   Double_t* errorsNegGamma[ErrorGammaSpec];
   Double_t errorsNegCorrGamma[ErrorGammaSpec][ConstnBinsGamma];
   Double_t errorsNegSummedGamma[ConstnBinsGamma];
   Double_t errorsNegCorrSummedGamma[ConstnBinsGamma];
   Double_t errorsNegCorrMatSummedGamma[ConstnBinsGamma];

   Double_t* errorsNegErrGamma[ErrorGammaSpec];
   Double_t errorsNegErrCorrGamma[ErrorGammaSpec][ConstnBinsGamma];
   Double_t errorsNegErrSummedGamma[ConstnBinsGamma];
   Double_t errorsNegErrCorrSummedGamma[ConstnBinsGamma];

   Double_t* errorsPosGamma[ErrorGammaSpec];
   Double_t errorsPosCorrGamma[ErrorGammaSpec][ConstnBinsGamma];
   Double_t errorsPosSummedGamma[ConstnBinsGamma];
   Double_t errorsPosCorrSummedGamma[ConstnBinsGamma];
   Double_t errorsPosCorrMatSummedGamma[ConstnBinsGamma];

   Double_t* errorsPosErrGamma[ErrorGammaSpec];
   Double_t errorsPosErrSummedGamma[ConstnBinsGamma];
   Double_t errorsPosErrCorrGamma[ErrorGammaSpec][ConstnBinsGamma];
   Double_t errorsPosErrCorrSummedGamma[ConstnBinsGamma];

   Double_t* errorsNegGammaA[ErrorGammaSpec];
   Double_t errorsNegCorrGammaA[ErrorGammaSpec][ConstnBinsGamma];
   Double_t errorsNegSummedGammaA[ConstnBinsGamma];
   Double_t errorsNegCorrSummedGammaA[ConstnBinsGamma];
   Double_t errorsNegCorrMatSummedGammaA[ConstnBinsGamma];

   Double_t* errorsNegErrGammaA[ErrorGammaSpec];
   Double_t errorsNegErrCorrGammaA[ErrorGammaSpec][ConstnBinsGamma];
   Double_t errorsNegErrSummedGammaA[ConstnBinsGamma];
   Double_t errorsNegErrCorrSummedGammaA[ConstnBinsGamma];

   Double_t* errorsPosGammaA[ErrorGammaSpec];
   Double_t errorsPosCorrGammaA[ErrorGammaSpec][ConstnBinsGamma];
   Double_t errorsPosSummedGammaA[ConstnBinsGamma];
   Double_t errorsPosCorrSummedGammaA[ConstnBinsGamma];
   Double_t errorsPosCorrMatSummedGammaA[ConstnBinsGamma];

   Double_t* errorsPosErrGammaA[ErrorGammaSpec];
   Double_t errorsPosErrSummedGammaA[ConstnBinsGamma];
   Double_t errorsPosErrCorrGammaA[ErrorGammaSpec][ConstnBinsGamma];
   Double_t errorsPosErrCorrSummedGammaA[ConstnBinsGamma];

   Double_t* errorsNegGammaB[ErrorGammaSpec];
   Double_t errorsNegCorrGammaB[ErrorGammaSpec][ConstnBinsGamma];
   Double_t errorsNegSummedGammaB[ConstnBinsGamma];
   Double_t errorsNegCorrSummedGammaB[ConstnBinsGamma];
   Double_t errorsNegCorrMatSummedGammaB[ConstnBinsGamma];

   Double_t* errorsNegErrGammaB[ErrorGammaSpec];
   Double_t errorsNegErrCorrGammaB[ErrorGammaSpec][ConstnBinsGamma];
   Double_t errorsNegErrSummedGammaB[ConstnBinsGamma];
   Double_t errorsNegErrCorrSummedGammaB[ConstnBinsGamma];

   Double_t* errorsPosGammaB[ErrorGammaSpec];
   Double_t errorsPosCorrGammaB[ErrorGammaSpec][ConstnBinsGamma];
   Double_t errorsPosSummedGammaB[ConstnBinsGamma];
   Double_t errorsPosCorrSummedGammaB[ConstnBinsGamma];
   Double_t errorsPosCorrMatSummedGammaB[ConstnBinsGamma];

   Double_t* errorsPosErrGammaB[ErrorGammaSpec];
   Double_t errorsPosErrSummedGammaB[ConstnBinsGamma];
   Double_t errorsPosErrCorrGammaB[ErrorGammaSpec][ConstnBinsGamma];
   Double_t errorsPosErrCorrSummedGammaB[ConstnBinsGamma];

   Double_t* errorsNegGammaC[ErrorGammaSpec];
   Double_t errorsNegCorrGammaC[ErrorGammaSpec][ConstnBinsGamma];
   Double_t errorsNegSummedGammaC[ConstnBinsGamma];
   Double_t errorsNegCorrSummedGammaC[ConstnBinsGamma];
   Double_t errorsNegCorrMatSummedGammaC[ConstnBinsGamma];

   Double_t* errorsNegErrGammaC[ErrorGammaSpec];
   Double_t errorsNegErrCorrGammaC[ErrorGammaSpec][ConstnBinsGamma];
   Double_t errorsNegErrSummedGammaC[ConstnBinsGamma];
   Double_t errorsNegErrCorrSummedGammaC[ConstnBinsGamma];

   Double_t* errorsPosGammaC[ErrorGammaSpec];
   Double_t errorsPosCorrGammaC[ErrorGammaSpec][ConstnBinsGamma];
   Double_t errorsPosSummedGammaC[ConstnBinsGamma];
   Double_t errorsPosCorrSummedGammaC[ConstnBinsGamma];
   Double_t errorsPosCorrMatSummedGammaC[ConstnBinsGamma];

   Double_t* errorsPosErrGammaC[ErrorGammaSpec];
   Double_t errorsPosErrSummedGammaC[ConstnBinsGamma];
   Double_t errorsPosErrCorrGammaC[ErrorGammaSpec][ConstnBinsGamma];
   Double_t errorsPosErrCorrSummedGammaC[ConstnBinsGamma];

   Double_t errorsMeanDoubleRatio[ErrorDoubleRatio][ConstnBinsDoubleRatio];
   Double_t errorsMeanCorrDoubleRatio[ErrorDoubleRatio][ConstnBinsDoubleRatio];
   Double_t errorsMeanSummedDoubleRatio[ConstnBinsDoubleRatio];
   Double_t errorsMeanCorrSummedDoubleRatio[ConstnBinsDoubleRatio];
   Double_t errorsMeanCorrMatSummedDoubleRatio[ConstnBinsDoubleRatio];

   Double_t errorsMeanErrDoubleRatio[ErrorDoubleRatio][ConstnBinsDoubleRatio];
   Double_t errorsMeanErrCorrDoubleRatio[ErrorDoubleRatio][ConstnBinsDoubleRatio];
   Double_t errorsMeanErrSummedDoubleRatio[ConstnBinsDoubleRatio];
   Double_t errorsMeanErrCorrSummedDoubleRatio[ConstnBinsDoubleRatio];
   Double_t errorsMeanErrCorrMatSummedDoubleRatio[ConstnBinsDoubleRatio];

   Double_t errorsMeanDoubleRatioPi0Fit[ErrorDoubleRatio][ConstnBinsDoubleRatio];
   Double_t errorsMeanCorrDoubleRatioPi0Fit[ErrorDoubleRatio][ConstnBinsDoubleRatio];
   Double_t errorsMeanSummedDoubleRatioPi0Fit[ConstnBinsDoubleRatio];
   Double_t errorsMeanCorrSummedDoubleRatioPi0Fit[ConstnBinsDoubleRatio];
   Double_t errorsMeanCorrMatSummedDoubleRatioPi0Fit[ConstnBinsDoubleRatio];

   Double_t errorsMeanErrDoubleRatioPi0Fit[ErrorDoubleRatio][ConstnBinsDoubleRatio];
   Double_t errorsMeanErrCorrDoubleRatioPi0Fit[ErrorDoubleRatio][ConstnBinsDoubleRatio];
   Double_t errorsMeanErrSummedDoubleRatioPi0Fit[ConstnBinsDoubleRatio];
   Double_t errorsMeanErrCorrSummedDoubleRatioPi0Fit[ConstnBinsDoubleRatio];
   Double_t errorsMeanErrCorrMatSummedDoubleRatioPi0Fit[ConstnBinsDoubleRatio];

   Double_t errorsMeanDoubleRatioA[ErrorDoubleRatio][ConstnBinsDoubleRatio];
   Double_t errorsMeanCorrDoubleRatioA[ErrorDoubleRatio][ConstnBinsDoubleRatio];
   Double_t errorsMeanSummedDoubleRatioA[ConstnBinsDoubleRatio];
   Double_t errorsMeanCorrSummedDoubleRatioA[ConstnBinsDoubleRatio];
   Double_t errorsMeanCorrMatSummedDoubleRatioA[ConstnBinsDoubleRatio];

   Double_t errorsMeanErrDoubleRatioA[ErrorDoubleRatio][ConstnBinsDoubleRatio];
   Double_t errorsMeanErrCorrDoubleRatioA[ErrorDoubleRatio][ConstnBinsDoubleRatio];
   Double_t errorsMeanErrSummedDoubleRatioA[ConstnBinsDoubleRatio];
   Double_t errorsMeanErrCorrSummedDoubleRatioA[ConstnBinsDoubleRatio];
   Double_t errorsMeanErrCorrMatSummedDoubleRatioA[ConstnBinsDoubleRatio];

   Double_t errorsMeanDoubleRatioB[ErrorDoubleRatio][ConstnBinsDoubleRatio];
   Double_t errorsMeanCorrDoubleRatioB[ErrorDoubleRatio][ConstnBinsDoubleRatio];
   Double_t errorsMeanSummedDoubleRatioB[ConstnBinsDoubleRatio];
   Double_t errorsMeanCorrSummedDoubleRatioB[ConstnBinsDoubleRatio];
   Double_t errorsMeanCorrMatSummedDoubleRatioB[ConstnBinsDoubleRatio];

   Double_t errorsMeanErrDoubleRatioB[ErrorDoubleRatio][ConstnBinsDoubleRatio];
   Double_t errorsMeanErrCorrDoubleRatioB[ErrorDoubleRatio][ConstnBinsDoubleRatio];
   Double_t errorsMeanErrSummedDoubleRatioB[ConstnBinsDoubleRatio];
   Double_t errorsMeanErrCorrSummedDoubleRatioB[ConstnBinsDoubleRatio];
   Double_t errorsMeanErrCorrMatSummedDoubleRatioB[ConstnBinsDoubleRatio];

   Double_t errorsMeanDoubleRatioC[ErrorDoubleRatio][ConstnBinsDoubleRatio];
   Double_t errorsMeanCorrDoubleRatioC[ErrorDoubleRatio][ConstnBinsDoubleRatio];
   Double_t errorsMeanSummedDoubleRatioC[ConstnBinsDoubleRatio];
   Double_t errorsMeanCorrSummedDoubleRatioC[ConstnBinsDoubleRatio];
   Double_t errorsMeanCorrMatSummedDoubleRatioC[ConstnBinsDoubleRatio];

   Double_t errorsMeanErrDoubleRatioC[ErrorDoubleRatio][ConstnBinsDoubleRatio];
   Double_t errorsMeanErrCorrDoubleRatioC[ErrorDoubleRatio][ConstnBinsDoubleRatio];
   Double_t errorsMeanErrSummedDoubleRatioC[ConstnBinsDoubleRatio];
   Double_t errorsMeanErrCorrSummedDoubleRatioC[ConstnBinsDoubleRatio];
   Double_t errorsMeanErrCorrMatSummedDoubleRatioC[ConstnBinsDoubleRatio];


   Double_t errorsMeanDoubleRatioPi0FitA[ErrorDoubleRatio][ConstnBinsDoubleRatio];
   Double_t errorsMeanCorrDoubleRatioPi0FitA[ErrorDoubleRatio][ConstnBinsDoubleRatio];
   Double_t errorsMeanSummedDoubleRatioPi0FitA[ConstnBinsDoubleRatio];
   Double_t errorsMeanCorrSummedDoubleRatioPi0FitA[ConstnBinsDoubleRatio];
   Double_t errorsMeanCorrMatSummedDoubleRatioPi0FitA[ConstnBinsDoubleRatio];

   Double_t errorsMeanErrDoubleRatioPi0FitA[ErrorDoubleRatio][ConstnBinsDoubleRatio];
   Double_t errorsMeanErrCorrDoubleRatioPi0FitA[ErrorDoubleRatio][ConstnBinsDoubleRatio];
   Double_t errorsMeanErrSummedDoubleRatioPi0FitA[ConstnBinsDoubleRatio];
   Double_t errorsMeanErrCorrSummedDoubleRatioPi0FitA[ConstnBinsDoubleRatio];
   Double_t errorsMeanErrCorrMatSummedDoubleRatioPi0FitA[ConstnBinsDoubleRatio];

   Double_t errorsMeanDoubleRatioPi0FitB[ErrorDoubleRatio][ConstnBinsDoubleRatio];
   Double_t errorsMeanCorrDoubleRatioPi0FitB[ErrorDoubleRatio][ConstnBinsDoubleRatio];
   Double_t errorsMeanSummedDoubleRatioPi0FitB[ConstnBinsDoubleRatio];
   Double_t errorsMeanCorrSummedDoubleRatioPi0FitB[ConstnBinsDoubleRatio];
   Double_t errorsMeanCorrMatSummedDoubleRatioPi0FitB[ConstnBinsDoubleRatio];

   Double_t errorsMeanErrDoubleRatioPi0FitB[ErrorDoubleRatio][ConstnBinsDoubleRatio];
   Double_t errorsMeanErrCorrDoubleRatioPi0FitB[ErrorDoubleRatio][ConstnBinsDoubleRatio];
   Double_t errorsMeanErrSummedDoubleRatioPi0FitB[ConstnBinsDoubleRatio];
   Double_t errorsMeanErrCorrSummedDoubleRatioPi0FitB[ConstnBinsDoubleRatio];
   Double_t errorsMeanErrCorrMatSummedDoubleRatioPi0FitB[ConstnBinsDoubleRatio];

   Double_t errorsMeanDoubleRatioPi0FitC[ErrorDoubleRatio][ConstnBinsDoubleRatio];
   Double_t errorsMeanCorrDoubleRatioPi0FitC[ErrorDoubleRatio][ConstnBinsDoubleRatio];
   Double_t errorsMeanSummedDoubleRatioPi0FitC[ConstnBinsDoubleRatio];
   Double_t errorsMeanCorrSummedDoubleRatioPi0FitC[ConstnBinsDoubleRatio];
   Double_t errorsMeanCorrMatSummedDoubleRatioPi0FitC[ConstnBinsDoubleRatio];

   Double_t errorsMeanErrDoubleRatioPi0FitC[ErrorDoubleRatio][ConstnBinsDoubleRatio];
   Double_t errorsMeanErrCorrDoubleRatioPi0FitC[ErrorDoubleRatio][ConstnBinsDoubleRatio];
   Double_t errorsMeanErrSummedDoubleRatioPi0FitC[ConstnBinsDoubleRatio];
   Double_t errorsMeanErrCorrSummedDoubleRatioPi0FitC[ConstnBinsDoubleRatio];
   Double_t errorsMeanErrCorrMatSummedDoubleRatioPi0FitC[ConstnBinsDoubleRatio];

   Double_t errorsMeanCorrMatSummedDoubleRatioPi0FitWithoutMat[ConstnBinsDoubleRatio];
   Double_t errorsMaterialBudgetDoubleRatioPi0Fit[ConstnBinsDoubleRatio];
   Double_t errorsMeanErrCorrMatSummedDoubleRatioPi0FitWithoutMat[ConstnBinsDoubleRatio];
   Double_t errorsErrMaterialBudgetDoubleRatioPi0Fit[ConstnBinsDoubleRatio];

   Double_t errorsMeanIncRatio[ErrorIncRatio][ConstnBinsIncRatio];
   Double_t errorsMeanCorrIncRatio[ErrorIncRatio][ConstnBinsIncRatio];
   Double_t errorsMeanSummedIncRatio[ConstnBinsIncRatio];
   Double_t errorsMeanCorrSummedIncRatio[ConstnBinsIncRatio];
   Double_t errorsMeanCorrMatSummedIncRatio[ConstnBinsIncRatio];

   Double_t errorsMeanErrIncRatio[ErrorIncRatio][ConstnBinsIncRatio];
   Double_t errorsMeanErrCorrIncRatio[ErrorIncRatio][ConstnBinsIncRatio];
   Double_t errorsMeanErrSummedIncRatio[ConstnBinsIncRatio];
   Double_t errorsMeanErrCorrSummedIncRatio[ConstnBinsIncRatio];
   Double_t errorsMeanErrCorrMatSummedIncRatio[ConstnBinsIncRatio];

   Double_t errorsMeanIncRatioPi0Fit[ErrorIncRatio][ConstnBinsIncRatio];
   Double_t errorsMeanCorrIncRatioPi0Fit[ErrorIncRatio][ConstnBinsIncRatio];
   Double_t errorsMeanSummedIncRatioPi0Fit[ConstnBinsIncRatio];
   Double_t errorsMeanCorrSummedIncRatioPi0Fit[ConstnBinsIncRatio];
   Double_t errorsMeanCorrMatSummedIncRatioPi0Fit[ConstnBinsIncRatio];

   Double_t errorsMeanErrIncRatioPi0Fit[ErrorIncRatio][ConstnBinsIncRatio];
   Double_t errorsMeanErrCorrIncRatioPi0Fit[ErrorIncRatio][ConstnBinsIncRatio];
   Double_t errorsMeanErrSummedIncRatioPi0Fit[ConstnBinsIncRatio];
   Double_t errorsMeanErrCorrSummedIncRatioPi0Fit[ConstnBinsIncRatio];
   Double_t errorsMeanErrCorrMatSummedIncRatioPi0Fit[ConstnBinsIncRatio];

   Double_t errorsMeanIncRatioA[ErrorIncRatio][ConstnBinsIncRatio];
   Double_t errorsMeanCorrIncRatioA[ErrorIncRatio][ConstnBinsIncRatio];
   Double_t errorsMeanSummedIncRatioA[ConstnBinsIncRatio];
   Double_t errorsMeanCorrSummedIncRatioA[ConstnBinsIncRatio];
   Double_t errorsMeanCorrMatSummedIncRatioA[ConstnBinsIncRatio];

   Double_t errorsMeanErrIncRatioA[ErrorIncRatio][ConstnBinsIncRatio];
   Double_t errorsMeanErrCorrIncRatioA[ErrorIncRatio][ConstnBinsIncRatio];
   Double_t errorsMeanErrSummedIncRatioA[ConstnBinsIncRatio];
   Double_t errorsMeanErrCorrSummedIncRatioA[ConstnBinsIncRatio];
   Double_t errorsMeanErrCorrMatSummedIncRatioA[ConstnBinsIncRatio];

   Double_t errorsMeanIncRatioB[ErrorIncRatio][ConstnBinsIncRatio];
   Double_t errorsMeanCorrIncRatioB[ErrorIncRatio][ConstnBinsIncRatio];
   Double_t errorsMeanSummedIncRatioB[ConstnBinsIncRatio];
   Double_t errorsMeanCorrSummedIncRatioB[ConstnBinsIncRatio];
   Double_t errorsMeanCorrMatSummedIncRatioB[ConstnBinsIncRatio];

   Double_t errorsMeanErrIncRatioB[ErrorIncRatio][ConstnBinsIncRatio];
   Double_t errorsMeanErrCorrIncRatioB[ErrorIncRatio][ConstnBinsIncRatio];
   Double_t errorsMeanErrSummedIncRatioB[ConstnBinsIncRatio];
   Double_t errorsMeanErrCorrSummedIncRatioB[ConstnBinsIncRatio];
   Double_t errorsMeanErrCorrMatSummedIncRatioB[ConstnBinsIncRatio];

   Double_t errorsMeanIncRatioC[ErrorIncRatio][ConstnBinsIncRatio];
   Double_t errorsMeanCorrIncRatioC[ErrorIncRatio][ConstnBinsIncRatio];
   Double_t errorsMeanSummedIncRatioC[ConstnBinsIncRatio];
   Double_t errorsMeanCorrSummedIncRatioC[ConstnBinsIncRatio];
   Double_t errorsMeanCorrMatSummedIncRatioC[ConstnBinsIncRatio];

   Double_t errorsMeanErrIncRatioC[ErrorIncRatio][ConstnBinsIncRatio];
   Double_t errorsMeanErrCorrIncRatioC[ErrorIncRatio][ConstnBinsIncRatio];
   Double_t errorsMeanErrSummedIncRatioC[ConstnBinsIncRatio];
   Double_t errorsMeanErrCorrSummedIncRatioC[ConstnBinsIncRatio];
   Double_t errorsMeanErrCorrMatSummedIncRatioC[ConstnBinsIncRatio];

   Double_t errorsMeanIncRatioPi0FitA[ErrorIncRatio][ConstnBinsIncRatio];
   Double_t errorsMeanCorrIncRatioPi0FitA[ErrorIncRatio][ConstnBinsIncRatio];
   Double_t errorsMeanSummedIncRatioPi0FitA[ConstnBinsIncRatio];
   Double_t errorsMeanCorrSummedIncRatioPi0FitA[ConstnBinsIncRatio];
   Double_t errorsMeanCorrMatSummedIncRatioPi0FitA[ConstnBinsIncRatio];

   Double_t errorsMeanErrIncRatioPi0FitA[ErrorIncRatio][ConstnBinsIncRatio];
   Double_t errorsMeanErrCorrIncRatioPi0FitA[ErrorIncRatio][ConstnBinsIncRatio];
   Double_t errorsMeanErrSummedIncRatioPi0FitA[ConstnBinsIncRatio];
   Double_t errorsMeanErrCorrSummedIncRatioPi0FitA[ConstnBinsIncRatio];
   Double_t errorsMeanErrCorrMatSummedIncRatioPi0FitA[ConstnBinsIncRatio];

   Double_t errorsMeanIncRatioPi0FitB[ErrorIncRatio][ConstnBinsIncRatio];
   Double_t errorsMeanCorrIncRatioPi0FitB[ErrorIncRatio][ConstnBinsIncRatio];
   Double_t errorsMeanSummedIncRatioPi0FitB[ConstnBinsIncRatio];
   Double_t errorsMeanCorrSummedIncRatioPi0FitB[ConstnBinsIncRatio];
   Double_t errorsMeanCorrMatSummedIncRatioPi0FitB[ConstnBinsIncRatio];

   Double_t errorsMeanErrIncRatioPi0FitB[ErrorIncRatio][ConstnBinsIncRatio];
   Double_t errorsMeanErrCorrIncRatioPi0FitB[ErrorIncRatio][ConstnBinsIncRatio];
   Double_t errorsMeanErrSummedIncRatioPi0FitB[ConstnBinsIncRatio];
   Double_t errorsMeanErrCorrSummedIncRatioPi0FitB[ConstnBinsIncRatio];
   Double_t errorsMeanErrCorrMatSummedIncRatioPi0FitB[ConstnBinsIncRatio];

   Double_t errorsMeanIncRatioPi0FitC[ErrorIncRatio][ConstnBinsIncRatio];
   Double_t errorsMeanCorrIncRatioPi0FitC[ErrorIncRatio][ConstnBinsIncRatio];
   Double_t errorsMeanSummedIncRatioPi0FitC[ConstnBinsIncRatio];
   Double_t errorsMeanCorrSummedIncRatioPi0FitC[ConstnBinsIncRatio];
   Double_t errorsMeanCorrMatSummedIncRatioPi0FitC[ConstnBinsIncRatio];

   Double_t errorsMeanErrIncRatioPi0FitC[ErrorIncRatio][ConstnBinsIncRatio];
   Double_t errorsMeanErrCorrIncRatioPi0FitC[ErrorIncRatio][ConstnBinsIncRatio];
   Double_t errorsMeanErrSummedIncRatioPi0FitC[ConstnBinsIncRatio];
   Double_t errorsMeanErrCorrSummedIncRatioPi0FitC[ConstnBinsIncRatio];
   Double_t errorsMeanErrCorrMatSummedIncRatioPi0FitC[ConstnBinsIncRatio];

   Double_t errorsMeanPi0[ErrorPi0][ConstnBinsPi0];
   Double_t errorsMeanCorrPi0[ErrorPi0][ConstnBinsPi0];
   Double_t errorsMeanSummedPi0[ConstnBinsPi0];
   Double_t errorsMeanCorrSummedPi0[ConstnBinsPi0];
   Double_t errorsMeanCorrMatSummedPi0[ConstnBinsPi0];
   Double_t errorsMeanCorrMatSummedPi0WithoutMaterial[ConstnBinsPi0];

   Double_t errorsMeanErrPi0[ErrorPi0][ConstnBinsPi0];
   Double_t errorsMeanErrCorrPi0[ErrorPi0][ConstnBinsPi0];
   Double_t errorsMeanErrSummedPi0[ConstnBinsPi0];
   Double_t errorsMeanErrCorrSummedPi0[ConstnBinsPi0];
   Double_t errorsMeanErrCorrMatSummedPi0[ConstnBinsPi0];

   Double_t errorsMeanPi0Fit[ErrorPi0][ConstnBinsPi0];
   Double_t errorsMeanCorrPi0Fit[ErrorPi0][ConstnBinsPi0];
   Double_t errorsMeanSummedPi0Fit[ConstnBinsPi0];
   Double_t errorsMeanCorrSummedPi0Fit[ConstnBinsPi0];
   Double_t errorsMeanCorrMatSummedPi0Fit[ConstnBinsPi0];
   Double_t errorsMeanCorrMatSummedPi0FitWithoutMaterial[ConstnBinsPi0];

   Double_t errorsMeanErrPi0Fit[ErrorPi0][ConstnBinsPi0];
   Double_t errorsMeanErrCorrPi0Fit[ErrorPi0][ConstnBinsPi0];
   Double_t errorsMeanErrSummedPi0Fit[ConstnBinsPi0];
   Double_t errorsMeanErrCorrSummedPi0Fit[ConstnBinsPi0];
   Double_t errorsMeanErrCorrMatSummedPi0Fit[ConstnBinsPi0];

   Double_t errorsMeanGamma[ErrorGammaSpec][ConstnBinsGamma];
   Double_t errorsMeanCorrGamma[ErrorGammaSpec][ConstnBinsGamma];
   Double_t errorsMeanSummedGamma[ConstnBinsGamma];
   Double_t errorsMeanCorrSummedGamma[ConstnBinsGamma];
   Double_t errorsMeanCorrMatSummedGamma[ConstnBinsGamma];
   Double_t errorsMeanCorrMatSummedGammaWithoutMat[ConstnBinsGamma];
   Double_t errorsMaterialBudget[ConstnBinsGamma];

   Double_t errorsMeanErrGamma[ErrorGammaSpec][ConstnBinsGamma];
   Double_t errorsMeanErrCorrGamma[ErrorGammaSpec][ConstnBinsGamma];
   Double_t errorsMeanErrSummedGamma[ConstnBinsGamma];
   Double_t errorsMeanErrCorrSummedGamma[ConstnBinsGamma];
   Double_t errorsMeanErrCorrMatSummedGamma[ConstnBinsGamma];
   Double_t errorsMeanErrCorrMatSummedGammaWithoutMat[ConstnBinsGamma];
   Double_t errorsErrMaterialBudget[ConstnBinsGamma];

   Double_t errorsMeanGammaA[ErrorGammaSpec][ConstnBinsGamma];
   Double_t errorsMeanCorrGammaA[ErrorGammaSpec][ConstnBinsGamma];
   Double_t errorsMeanSummedGammaA[ConstnBinsGamma];
   Double_t errorsMeanCorrSummedGammaA[ConstnBinsGamma];
   Double_t errorsMeanCorrMatSummedGammaA[ConstnBinsGamma];

   Double_t errorsMeanErrGammaA[ErrorGammaSpec][ConstnBinsGamma];
   Double_t errorsMeanErrCorrGammaA[ErrorGammaSpec][ConstnBinsGamma];
   Double_t errorsMeanErrSummedGammaA[ConstnBinsGamma];
   Double_t errorsMeanErrCorrSummedGammaA[ConstnBinsGamma];
   Double_t errorsMeanErrCorrMatSummedGammaA[ConstnBinsGamma];

   Double_t errorsMeanGammaB[ErrorGammaSpec][ConstnBinsGamma];
   Double_t errorsMeanCorrGammaB[ErrorGammaSpec][ConstnBinsGamma];
   Double_t errorsMeanSummedGammaB[ConstnBinsGamma];
   Double_t errorsMeanCorrSummedGammaB[ConstnBinsGamma];
   Double_t errorsMeanCorrMatSummedGammaB[ConstnBinsGamma];

   Double_t errorsMeanErrGammaB[ErrorGammaSpec][ConstnBinsGamma];
   Double_t errorsMeanErrCorrGammaB[ErrorGammaSpec][ConstnBinsGamma];
   Double_t errorsMeanErrSummedGammaB[ConstnBinsGamma];
   Double_t errorsMeanErrCorrSummedGammaB[ConstnBinsGamma];
   Double_t errorsMeanErrCorrMatSummedGammaB[ConstnBinsGamma];

   Double_t errorsMeanGammaC[ErrorGammaSpec][ConstnBinsGamma];
   Double_t errorsMeanCorrGammaC[ErrorGammaSpec][ConstnBinsGamma];
   Double_t errorsMeanSummedGammaC[ConstnBinsGamma];
   Double_t errorsMeanCorrSummedGammaC[ConstnBinsGamma];
   Double_t errorsMeanCorrMatSummedGammaC[ConstnBinsGamma];

   Double_t errorsMeanErrGammaC[ErrorGammaSpec][ConstnBinsGamma];
   Double_t errorsMeanErrCorrGammaC[ErrorGammaSpec][ConstnBinsGamma];
   Double_t errorsMeanErrSummedGammaC[ConstnBinsGamma];
   Double_t errorsMeanErrCorrSummedGammaC[ConstnBinsGamma];
   Double_t errorsMeanErrCorrMatSummedGammaC[ConstnBinsGamma];


   for (Int_t l = 0; l < ConstnBinsDoubleRatio; l++){
      errorsPosSummedDoubleRatio[l] = 0.;
      errorsNegSummedDoubleRatio[l] = 0.;
      errorsMeanSummedDoubleRatio[l] = 0.;
      errorsPosCorrSummedDoubleRatio[l] = 0.;
      errorsNegCorrSummedDoubleRatio[l] = 0.;
      errorsMeanCorrSummedDoubleRatio[l] = 0.;

      errorsPosSummedDoubleRatioPi0Fit[l] = 0.;
      errorsNegSummedDoubleRatioPi0Fit[l] = 0.;
      errorsMeanSummedDoubleRatioPi0Fit[l] = 0.;
      errorsPosCorrSummedDoubleRatioPi0Fit[l] = 0.;
      errorsNegCorrSummedDoubleRatioPi0Fit[l] = 0.;
      errorsMeanCorrSummedDoubleRatioPi0Fit[l] = 0.;

      errorsPosSummedDoubleRatioA[l] = 0.;
      errorsNegSummedDoubleRatioA[l] = 0.;
      errorsMeanSummedDoubleRatioA[l] = 0.;
      errorsPosCorrSummedDoubleRatioA[l] = 0.;
      errorsNegCorrSummedDoubleRatioA[l] = 0.;
      errorsMeanCorrSummedDoubleRatioA[l] = 0.;

      errorsPosSummedDoubleRatioPi0FitA[l] = 0.;
      errorsNegSummedDoubleRatioPi0FitA[l] = 0.;
      errorsMeanSummedDoubleRatioPi0FitA[l] = 0.;
      errorsPosCorrSummedDoubleRatioPi0FitA[l] = 0.;
      errorsNegCorrSummedDoubleRatioPi0FitA[l] = 0.;
      errorsMeanCorrSummedDoubleRatioPi0FitA[l] = 0.;

      errorsPosSummedDoubleRatioB[l] = 0.;
      errorsNegSummedDoubleRatioB[l] = 0.;
      errorsMeanSummedDoubleRatioB[l] = 0.;
      errorsPosCorrSummedDoubleRatioB[l] = 0.;
      errorsNegCorrSummedDoubleRatioB[l] = 0.;
      errorsMeanCorrSummedDoubleRatioB[l] = 0.;

      errorsPosSummedDoubleRatioPi0FitB[l] = 0.;
      errorsNegSummedDoubleRatioPi0FitB[l] = 0.;
      errorsMeanSummedDoubleRatioPi0FitB[l] = 0.;
      errorsPosCorrSummedDoubleRatioPi0FitB[l] = 0.;
      errorsNegCorrSummedDoubleRatioPi0FitB[l] = 0.;
      errorsMeanCorrSummedDoubleRatioPi0FitB[l] = 0.;

      errorsPosSummedDoubleRatioC[l] = 0.;
      errorsNegSummedDoubleRatioC[l] = 0.;
      errorsMeanSummedDoubleRatioC[l] = 0.;
      errorsPosCorrSummedDoubleRatioC[l] = 0.;
      errorsNegCorrSummedDoubleRatioC[l] = 0.;
      errorsMeanCorrSummedDoubleRatioC[l] = 0.;

      errorsPosSummedDoubleRatioPi0FitC[l] = 0.;
      errorsNegSummedDoubleRatioPi0FitC[l] = 0.;
      errorsMeanSummedDoubleRatioPi0FitC[l] = 0.;
      errorsPosCorrSummedDoubleRatioPi0FitC[l] = 0.;
      errorsNegCorrSummedDoubleRatioPi0FitC[l] = 0.;
      errorsMeanCorrSummedDoubleRatioPi0FitC[l] = 0.;
   }


   for (Int_t l = 0; l < ConstnBinsIncRatio; l++){
      errorsPosSummedIncRatio[l] = 0.;
      errorsNegSummedIncRatio[l] = 0.;
      errorsMeanSummedIncRatio[l] = 0.;
      errorsPosCorrSummedIncRatio[l] = 0.;
      errorsNegCorrSummedIncRatio[l] = 0.;
      errorsMeanCorrSummedIncRatio[l] = 0.;

      errorsPosSummedIncRatioPi0Fit[l] = 0.;
      errorsNegSummedIncRatioPi0Fit[l] = 0.;
      errorsMeanSummedIncRatioPi0Fit[l] = 0.;
      errorsPosCorrSummedIncRatioPi0Fit[l] = 0.;
      errorsNegCorrSummedIncRatioPi0Fit[l] = 0.;
      errorsMeanCorrSummedIncRatioPi0Fit[l] = 0.;

      errorsPosSummedIncRatioA[l] = 0.;
      errorsNegSummedIncRatioA[l] = 0.;
      errorsMeanSummedIncRatioA[l] = 0.;
      errorsPosCorrSummedIncRatioA[l] = 0.;
      errorsNegCorrSummedIncRatioA[l] = 0.;
      errorsMeanCorrSummedIncRatioA[l] = 0.;

      errorsPosSummedIncRatioPi0FitA[l] = 0.;
      errorsNegSummedIncRatioPi0FitA[l] = 0.;
      errorsMeanSummedIncRatioPi0FitA[l] = 0.;
      errorsPosCorrSummedIncRatioPi0FitA[l] = 0.;
      errorsNegCorrSummedIncRatioPi0FitA[l] = 0.;
      errorsMeanCorrSummedIncRatioPi0FitA[l] = 0.;

      errorsPosSummedIncRatioB[l] = 0.;
      errorsNegSummedIncRatioB[l] = 0.;
      errorsMeanSummedIncRatioB[l] = 0.;
      errorsPosCorrSummedIncRatioB[l] = 0.;
      errorsNegCorrSummedIncRatioB[l] = 0.;
      errorsMeanCorrSummedIncRatioB[l] = 0.;

      errorsPosSummedIncRatioPi0FitB[l] = 0.;
      errorsNegSummedIncRatioPi0FitB[l] = 0.;
      errorsMeanSummedIncRatioPi0FitB[l] = 0.;
      errorsPosCorrSummedIncRatioPi0FitB[l] = 0.;
      errorsNegCorrSummedIncRatioPi0FitB[l] = 0.;
      errorsMeanCorrSummedIncRatioPi0FitB[l] = 0.;

      errorsPosSummedIncRatioC[l] = 0.;
      errorsNegSummedIncRatioC[l] = 0.;
      errorsMeanSummedIncRatioC[l] = 0.;
      errorsPosCorrSummedIncRatioC[l] = 0.;
      errorsNegCorrSummedIncRatioC[l] = 0.;
      errorsMeanCorrSummedIncRatioC[l] = 0.;

      errorsPosSummedIncRatioPi0FitC[l] = 0.;
      errorsNegSummedIncRatioPi0FitC[l] = 0.;
      errorsMeanSummedIncRatioPi0FitC[l] = 0.;
      errorsPosCorrSummedIncRatioPi0FitC[l] = 0.;
      errorsNegCorrSummedIncRatioPi0FitC[l] = 0.;
      errorsMeanCorrSummedIncRatioPi0FitC[l] = 0.;

   }

   for (Int_t l = 0; l < ConstnBinsPi0; l++){
      errorsPosSummedPi0[l] = 0.;
      errorsNegSummedPi0[l] = 0.;
      errorsMeanSummedPi0[l] = 0.;
      errorsPosCorrSummedPi0[l] = 0.;
      errorsNegCorrSummedPi0[l] = 0.;
      errorsMeanCorrSummedPi0[l] = 0.;

      errorsPosSummedPi0Fit[l] = 0.;
      errorsNegSummedPi0Fit[l] = 0.;
      errorsMeanSummedPi0Fit[l] = 0.;
      errorsPosCorrSummedPi0Fit[l] = 0.;
      errorsNegCorrSummedPi0Fit[l] = 0.;
      errorsMeanCorrSummedPi0Fit[l] = 0.;
   }


   for (Int_t l = 0; l < ConstnBinsGamma; l++){
      errorsPosSummedGamma[l] = 0.;
      errorsNegSummedGamma[l] = 0.;
      errorsMeanSummedGamma[l] = 0.;
      errorsPosCorrSummedGamma[l] = 0.;
      errorsNegCorrSummedGamma[l] = 0.;
      errorsMeanCorrSummedGamma[l] = 0.;

      errorsPosSummedGammaA[l] = 0.;
      errorsNegSummedGammaA[l] = 0.;
      errorsMeanSummedGammaA[l] = 0.;
      errorsPosCorrSummedGammaA[l] = 0.;
      errorsNegCorrSummedGammaA[l] = 0.;
      errorsMeanCorrSummedGammaA[l] = 0.;

      errorsPosSummedGammaB[l] = 0.;
      errorsNegSummedGammaB[l] = 0.;
      errorsMeanSummedGammaB[l] = 0.;
      errorsPosCorrSummedGammaB[l] = 0.;
      errorsNegCorrSummedGammaB[l] = 0.;
      errorsMeanCorrSummedGammaB[l] = 0.;

      errorsPosSummedGammaC[l] = 0.;
      errorsNegSummedGammaC[l] = 0.;
      errorsMeanSummedGammaC[l] = 0.;
      errorsPosCorrSummedGammaC[l] = 0.;
      errorsNegCorrSummedGammaC[l] = 0.;
      errorsMeanCorrSummedGammaC[l] = 0.;
   }

   TGraphErrors* negativeErrorsDoubleRatio[ErrorDoubleRatio];
   TGraphErrors* negativeErrorsSummedDoubleRatio;
   TGraphErrors* positiveErrorsDoubleRatio[ErrorDoubleRatio];
   TGraphErrors* positiveErrorsSummedDoubleRatio;
   TGraphErrors* negativeErrorsCorrDoubleRatio[ErrorDoubleRatio];
   TGraphErrors* negativeErrorsCorrSummedDoubleRatio;
   TGraphErrors* positiveErrorsCorrDoubleRatio[ErrorDoubleRatio];
   TGraphErrors* positiveErrorsCorrSummedDoubleRatio;
   TGraphErrors* meanErrorsDoubleRatio[ErrorDoubleRatio];
   TGraphErrors* meanErrorsSummedDoubleRatio;
   TGraphErrors* meanErrorsCorrDoubleRatio[ErrorDoubleRatio];
   TGraphErrors* meanErrorsCorrSummedDoubleRatio;
   TGraphErrors* meanErrorsCorrSummedIncMatDoubleRatio;

   TGraphErrors* negativeErrorsDoubleRatioPi0Fit[ErrorDoubleRatio];
   TGraphErrors* negativeErrorsSummedDoubleRatioPi0Fit;
   TGraphErrors* positiveErrorsDoubleRatioPi0Fit[ErrorDoubleRatio];
   TGraphErrors* positiveErrorsSummedDoubleRatioPi0Fit;
   TGraphErrors* negativeErrorsCorrDoubleRatioPi0Fit[ErrorDoubleRatio];
   TGraphErrors* negativeErrorsCorrSummedDoubleRatioPi0Fit;
   TGraphErrors* positiveErrorsCorrDoubleRatioPi0Fit[ErrorDoubleRatio];
   TGraphErrors* positiveErrorsCorrSummedDoubleRatioPi0Fit;
   TGraphErrors* meanErrorsDoubleRatioPi0Fit[ErrorDoubleRatio];
   TGraphErrors* meanErrorsSummedDoubleRatioPi0Fit;
   TGraphErrors* meanErrorsCorrDoubleRatioPi0Fit[ErrorDoubleRatio];
   TGraphErrors* meanErrorsCorrSummedDoubleRatioPi0Fit;
   TGraphErrors* meanErrorsCorrSummedIncMatDoubleRatioPi0Fit;
   TGraphErrors* meanErrorsCorrSummedWithoutMatDoubleRatioPi0Fit;
   TGraphErrors* meanErrorsCorrSummedMaterialDoubleRatioPi0Fit;
   TGraphErrors* meanErrorsCorrSummedMaterialDoubleRatioPi0FitCocktail;

   TGraphErrors* negativeErrorsDoubleRatioA[ErrorDoubleRatio];
   TGraphErrors* negativeErrorsSummedDoubleRatioA;
   TGraphErrors* positiveErrorsDoubleRatioA[ErrorDoubleRatio];
   TGraphErrors* positiveErrorsSummedDoubleRatioA;
   TGraphErrors* negativeErrorsCorrDoubleRatioA[ErrorDoubleRatio];
   TGraphErrors* negativeErrorsCorrSummedDoubleRatioA;
   TGraphErrors* positiveErrorsCorrDoubleRatioA[ErrorDoubleRatio];
   TGraphErrors* positiveErrorsCorrSummedDoubleRatioA;
   TGraphErrors* meanErrorsDoubleRatioA[ErrorDoubleRatio];
   TGraphErrors* meanErrorsSummedDoubleRatioA;
   TGraphErrors* meanErrorsCorrDoubleRatioA[ErrorDoubleRatio];
   TGraphErrors* meanErrorsCorrSummedDoubleRatioA;
   TGraphErrors* meanErrorsCorrSummedIncMatDoubleRatioA;

   TGraphErrors* negativeErrorsDoubleRatioPi0FitA[ErrorDoubleRatio];
   TGraphErrors* negativeErrorsSummedDoubleRatioPi0FitA;
   TGraphErrors* positiveErrorsDoubleRatioPi0FitA[ErrorDoubleRatio];
   TGraphErrors* positiveErrorsSummedDoubleRatioPi0FitA;
   TGraphErrors* negativeErrorsCorrDoubleRatioPi0FitA[ErrorDoubleRatio];
   TGraphErrors* negativeErrorsCorrSummedDoubleRatioPi0FitA;
   TGraphErrors* positiveErrorsCorrDoubleRatioPi0FitA[ErrorDoubleRatio];
   TGraphErrors* positiveErrorsCorrSummedDoubleRatioPi0FitA;
   TGraphErrors* meanErrorsDoubleRatioPi0FitA[ErrorDoubleRatio];
   TGraphErrors* meanErrorsSummedDoubleRatioPi0FitA;
   TGraphErrors* meanErrorsCorrDoubleRatioPi0FitA[ErrorDoubleRatio];
   TGraphErrors* meanErrorsCorrSummedDoubleRatioPi0FitA;
   TGraphErrors* meanErrorsCorrSummedIncMatDoubleRatioPi0FitA;
   TGraphErrors* meanErrorsCorrSummedWithoutMatDoubleRatioPi0FitA;
   TGraphErrors* meanErrorsCorrSummedMaterialDoubleRatioPi0FitA;

   TGraphErrors* negativeErrorsDoubleRatioB[ErrorDoubleRatio];
   TGraphErrors* negativeErrorsSummedDoubleRatioB;
   TGraphErrors* positiveErrorsDoubleRatioB[ErrorDoubleRatio];
   TGraphErrors* positiveErrorsSummedDoubleRatioB;
   TGraphErrors* negativeErrorsCorrDoubleRatioB[ErrorDoubleRatio];
   TGraphErrors* negativeErrorsCorrSummedDoubleRatioB;
   TGraphErrors* positiveErrorsCorrDoubleRatioB[ErrorDoubleRatio];
   TGraphErrors* positiveErrorsCorrSummedDoubleRatioB;
   TGraphErrors* meanErrorsDoubleRatioB[ErrorDoubleRatio];
   TGraphErrors* meanErrorsSummedDoubleRatioB;
   TGraphErrors* meanErrorsCorrDoubleRatioB[ErrorDoubleRatio];
   TGraphErrors* meanErrorsCorrSummedDoubleRatioB;
   TGraphErrors* meanErrorsCorrSummedIncMatDoubleRatioB;

   TGraphErrors* negativeErrorsDoubleRatioPi0FitB[ErrorDoubleRatio];
   TGraphErrors* negativeErrorsSummedDoubleRatioPi0FitB;
   TGraphErrors* positiveErrorsDoubleRatioPi0FitB[ErrorDoubleRatio];
   TGraphErrors* positiveErrorsSummedDoubleRatioPi0FitB;
   TGraphErrors* negativeErrorsCorrDoubleRatioPi0FitB[ErrorDoubleRatio];
   TGraphErrors* negativeErrorsCorrSummedDoubleRatioPi0FitB;
   TGraphErrors* positiveErrorsCorrDoubleRatioPi0FitB[ErrorDoubleRatio];
   TGraphErrors* positiveErrorsCorrSummedDoubleRatioPi0FitB;
   TGraphErrors* meanErrorsDoubleRatioPi0FitB[ErrorDoubleRatio];
   TGraphErrors* meanErrorsSummedDoubleRatioPi0FitB;
   TGraphErrors* meanErrorsCorrDoubleRatioPi0FitB[ErrorDoubleRatio];
   TGraphErrors* meanErrorsCorrSummedDoubleRatioPi0FitB;
   TGraphErrors* meanErrorsCorrSummedIncMatDoubleRatioPi0FitB;
   TGraphErrors* meanErrorsCorrSummedWithoutMatDoubleRatioPi0FitB;
   TGraphErrors* meanErrorsCorrSummedMaterialDoubleRatioPi0FitB;

   TGraphErrors* negativeErrorsDoubleRatioC[ErrorDoubleRatio];
   TGraphErrors* negativeErrorsSummedDoubleRatioC;
   TGraphErrors* positiveErrorsDoubleRatioC[ErrorDoubleRatio];
   TGraphErrors* positiveErrorsSummedDoubleRatioC;
   TGraphErrors* negativeErrorsCorrDoubleRatioC[ErrorDoubleRatio];
   TGraphErrors* negativeErrorsCorrSummedDoubleRatioC;
   TGraphErrors* positiveErrorsCorrDoubleRatioC[ErrorDoubleRatio];
   TGraphErrors* positiveErrorsCorrSummedDoubleRatioC;
   TGraphErrors* meanErrorsDoubleRatioC[ErrorDoubleRatio];
   TGraphErrors* meanErrorsSummedDoubleRatioC;
   TGraphErrors* meanErrorsCorrDoubleRatioC[ErrorDoubleRatio];
   TGraphErrors* meanErrorsCorrSummedDoubleRatioC;
   TGraphErrors* meanErrorsCorrSummedIncMatDoubleRatioC;

   TGraphErrors* negativeErrorsDoubleRatioPi0FitC[ErrorDoubleRatio];
   TGraphErrors* negativeErrorsSummedDoubleRatioPi0FitC;
   TGraphErrors* positiveErrorsDoubleRatioPi0FitC[ErrorDoubleRatio];
   TGraphErrors* positiveErrorsSummedDoubleRatioPi0FitC;
   TGraphErrors* negativeErrorsCorrDoubleRatioPi0FitC[ErrorDoubleRatio];
   TGraphErrors* negativeErrorsCorrSummedDoubleRatioPi0FitC;
   TGraphErrors* positiveErrorsCorrDoubleRatioPi0FitC[ErrorDoubleRatio];
   TGraphErrors* positiveErrorsCorrSummedDoubleRatioPi0FitC;
   TGraphErrors* meanErrorsDoubleRatioPi0FitC[ErrorDoubleRatio];
   TGraphErrors* meanErrorsSummedDoubleRatioPi0FitC;
   TGraphErrors* meanErrorsCorrDoubleRatioPi0FitC[ErrorDoubleRatio];
   TGraphErrors* meanErrorsCorrSummedDoubleRatioPi0FitC;
   TGraphErrors* meanErrorsCorrSummedIncMatDoubleRatioPi0FitC;
   TGraphErrors* meanErrorsCorrSummedWithoutMatDoubleRatioPi0FitC;
   TGraphErrors* meanErrorsCorrSummedMaterialDoubleRatioPi0FitC;


   TGraphErrors* negativeErrorsIncRatio[ErrorIncRatio];
   TGraphErrors* negativeErrorsSummedIncRatio;
   TGraphErrors* positiveErrorsIncRatio[ErrorIncRatio];
   TGraphErrors* positiveErrorsSummedIncRatio;
   TGraphErrors* negativeErrorsCorrIncRatio[ErrorIncRatio];
   TGraphErrors* negativeErrorsCorrSummedIncRatio;
   TGraphErrors* positiveErrorsCorrIncRatio[ErrorIncRatio];
   TGraphErrors* positiveErrorsCorrSummedIncRatio;
   TGraphErrors* meanErrorsIncRatio[ErrorIncRatio];
   TGraphErrors* meanErrorsSummedIncRatio;
   TGraphErrors* meanErrorsCorrIncRatio[ErrorIncRatio];
   TGraphErrors* meanErrorsCorrSummedIncRatio;
   TGraphErrors* meanErrorsCorrSummedIncMatIncRatio;

   TGraphErrors* negativeErrorsIncRatioPi0Fit[ErrorIncRatio];
   TGraphErrors* negativeErrorsSummedIncRatioPi0Fit;
   TGraphErrors* positiveErrorsIncRatioPi0Fit[ErrorIncRatio];
   TGraphErrors* positiveErrorsSummedIncRatioPi0Fit;
   TGraphErrors* negativeErrorsCorrIncRatioPi0Fit[ErrorIncRatio];
   TGraphErrors* negativeErrorsCorrSummedIncRatioPi0Fit;
   TGraphErrors* positiveErrorsCorrIncRatioPi0Fit[ErrorIncRatio];
   TGraphErrors* positiveErrorsCorrSummedIncRatioPi0Fit;
   TGraphErrors* meanErrorsIncRatioPi0Fit[ErrorIncRatio];
   TGraphErrors* meanErrorsSummedIncRatioPi0Fit;
   TGraphErrors* meanErrorsCorrIncRatioPi0Fit[ErrorIncRatio];
   TGraphErrors* meanErrorsCorrSummedIncRatioPi0Fit;
   TGraphErrors* meanErrorsCorrSummedIncMatIncRatioPi0Fit;

   TGraphErrors* negativeErrorsIncRatioA[ErrorIncRatio];
   TGraphErrors* negativeErrorsSummedIncRatioA;
   TGraphErrors* positiveErrorsIncRatioA[ErrorIncRatio];
   TGraphErrors* positiveErrorsSummedIncRatioA;
   TGraphErrors* negativeErrorsCorrIncRatioA[ErrorIncRatio];
   TGraphErrors* negativeErrorsCorrSummedIncRatioA;
   TGraphErrors* positiveErrorsCorrIncRatioA[ErrorIncRatio];
   TGraphErrors* positiveErrorsCorrSummedIncRatioA;
   TGraphErrors* meanErrorsIncRatioA[ErrorIncRatio];
   TGraphErrors* meanErrorsSummedIncRatioA;
   TGraphErrors* meanErrorsCorrIncRatioA[ErrorIncRatio];
   TGraphErrors* meanErrorsCorrSummedIncRatioA;
   TGraphErrors* meanErrorsCorrSummedIncMatIncRatioA;

   TGraphErrors* negativeErrorsIncRatioPi0FitA[ErrorIncRatio];
   TGraphErrors* negativeErrorsSummedIncRatioPi0FitA;
   TGraphErrors* positiveErrorsIncRatioPi0FitA[ErrorIncRatio];
   TGraphErrors* positiveErrorsSummedIncRatioPi0FitA;
   TGraphErrors* negativeErrorsCorrIncRatioPi0FitA[ErrorIncRatio];
   TGraphErrors* negativeErrorsCorrSummedIncRatioPi0FitA;
   TGraphErrors* positiveErrorsCorrIncRatioPi0FitA[ErrorIncRatio];
   TGraphErrors* positiveErrorsCorrSummedIncRatioPi0FitA;
   TGraphErrors* meanErrorsIncRatioPi0FitA[ErrorIncRatio];
   TGraphErrors* meanErrorsSummedIncRatioPi0FitA;
   TGraphErrors* meanErrorsCorrIncRatioPi0FitA[ErrorIncRatio];
   TGraphErrors* meanErrorsCorrSummedIncRatioPi0FitA;
   TGraphErrors* meanErrorsCorrSummedIncMatIncRatioPi0FitA;
   TGraphErrors* meanErrorsCorrSummedWithoutMatIncRatioPi0FitA;
   TGraphErrors* meanErrorsCorrSummedMaterialIncRatioPi0FitA;

   TGraphErrors* negativeErrorsIncRatioB[ErrorIncRatio];
   TGraphErrors* negativeErrorsSummedIncRatioB;
   TGraphErrors* positiveErrorsIncRatioB[ErrorIncRatio];
   TGraphErrors* positiveErrorsSummedIncRatioB;
   TGraphErrors* negativeErrorsCorrIncRatioB[ErrorIncRatio];
   TGraphErrors* negativeErrorsCorrSummedIncRatioB;
   TGraphErrors* positiveErrorsCorrIncRatioB[ErrorIncRatio];
   TGraphErrors* positiveErrorsCorrSummedIncRatioB;
   TGraphErrors* meanErrorsIncRatioB[ErrorIncRatio];
   TGraphErrors* meanErrorsSummedIncRatioB;
   TGraphErrors* meanErrorsCorrIncRatioB[ErrorIncRatio];
   TGraphErrors* meanErrorsCorrSummedIncRatioB;
   TGraphErrors* meanErrorsCorrSummedIncMatIncRatioB;

   TGraphErrors* negativeErrorsIncRatioPi0FitB[ErrorIncRatio];
   TGraphErrors* negativeErrorsSummedIncRatioPi0FitB;
   TGraphErrors* positiveErrorsIncRatioPi0FitB[ErrorIncRatio];
   TGraphErrors* positiveErrorsSummedIncRatioPi0FitB;
   TGraphErrors* negativeErrorsCorrIncRatioPi0FitB[ErrorIncRatio];
   TGraphErrors* negativeErrorsCorrSummedIncRatioPi0FitB;
   TGraphErrors* positiveErrorsCorrIncRatioPi0FitB[ErrorIncRatio];
   TGraphErrors* positiveErrorsCorrSummedIncRatioPi0FitB;
   TGraphErrors* meanErrorsIncRatioPi0FitB[ErrorIncRatio];
   TGraphErrors* meanErrorsSummedIncRatioPi0FitB;
   TGraphErrors* meanErrorsCorrIncRatioPi0FitB[ErrorIncRatio];
   TGraphErrors* meanErrorsCorrSummedIncRatioPi0FitB;
   TGraphErrors* meanErrorsCorrSummedIncMatIncRatioPi0FitB;
   TGraphErrors* meanErrorsCorrSummedWithoutMatIncRatioPi0FitB;
   TGraphErrors* meanErrorsCorrSummedMaterialIncRatioPi0FitB;

   TGraphErrors* negativeErrorsIncRatioC[ErrorIncRatio];
   TGraphErrors* negativeErrorsSummedIncRatioC;
   TGraphErrors* positiveErrorsIncRatioC[ErrorIncRatio];
   TGraphErrors* positiveErrorsSummedIncRatioC;
   TGraphErrors* negativeErrorsCorrIncRatioC[ErrorIncRatio];
   TGraphErrors* negativeErrorsCorrSummedIncRatioC;
   TGraphErrors* positiveErrorsCorrIncRatioC[ErrorIncRatio];
   TGraphErrors* positiveErrorsCorrSummedIncRatioC;
   TGraphErrors* meanErrorsIncRatioC[ErrorIncRatio];
   TGraphErrors* meanErrorsSummedIncRatioC;
   TGraphErrors* meanErrorsCorrIncRatioC[ErrorIncRatio];
   TGraphErrors* meanErrorsCorrSummedIncRatioC;
   TGraphErrors* meanErrorsCorrSummedIncMatIncRatioC;

   TGraphErrors* negativeErrorsIncRatioPi0FitC[ErrorIncRatio];
   TGraphErrors* negativeErrorsSummedIncRatioPi0FitC;
   TGraphErrors* positiveErrorsIncRatioPi0FitC[ErrorIncRatio];
   TGraphErrors* positiveErrorsSummedIncRatioPi0FitC;
   TGraphErrors* negativeErrorsCorrIncRatioPi0FitC[ErrorIncRatio];
   TGraphErrors* negativeErrorsCorrSummedIncRatioPi0FitC;
   TGraphErrors* positiveErrorsCorrIncRatioPi0FitC[ErrorIncRatio];
   TGraphErrors* positiveErrorsCorrSummedIncRatioPi0FitC;
   TGraphErrors* meanErrorsIncRatioPi0FitC[ErrorIncRatio];
   TGraphErrors* meanErrorsSummedIncRatioPi0FitC;
   TGraphErrors* meanErrorsCorrIncRatioPi0FitC[ErrorIncRatio];
   TGraphErrors* meanErrorsCorrSummedIncRatioPi0FitC;
   TGraphErrors* meanErrorsCorrSummedIncMatIncRatioPi0FitC;
   TGraphErrors* meanErrorsCorrSummedWithoutMatIncRatioPi0FitC;
   TGraphErrors* meanErrorsCorrSummedMaterialIncRatioPi0FitC;

   TGraphErrors* negativeErrorsPi0[ErrorPi0];
   TGraphErrors* negativeErrorsSummedPi0;
   TGraphErrors* positiveErrorsPi0[ErrorPi0];
   TGraphErrors* positiveErrorsSummedPi0;
   TGraphErrors* negativeErrorsCorrPi0[ErrorPi0];
   TGraphErrors* negativeErrorsCorrSummedPi0;
   TGraphErrors* positiveErrorsCorrPi0[ErrorPi0];
   TGraphErrors* positiveErrorsCorrSummedPi0;
   TGraphErrors* meanErrorsPi0[ErrorPi0];
   TGraphErrors* meanErrorsSummedPi0;
   TGraphErrors* meanErrorsCorrPi0[ErrorPi0];
   TGraphErrors* meanErrorsCorrSummedPi0;
   TGraphErrors* meanErrorsCorrSummedIncMatPi0;
   TGraphErrors* meanErrorsCorrSummedIncMatPi0WithoutMat;

   TGraphErrors* negativeErrorsPi0Fit[ErrorPi0];
   TGraphErrors* negativeErrorsSummedPi0Fit;
   TGraphErrors* positiveErrorsPi0Fit[ErrorPi0];
   TGraphErrors* positiveErrorsSummedPi0Fit;
   TGraphErrors* negativeErrorsCorrPi0Fit[ErrorPi0];
   TGraphErrors* negativeErrorsCorrSummedPi0Fit;
   TGraphErrors* positiveErrorsCorrPi0Fit[ErrorPi0];
   TGraphErrors* positiveErrorsCorrSummedPi0Fit;
   TGraphErrors* meanErrorsPi0Fit[ErrorPi0];
   TGraphErrors* meanErrorsSummedPi0Fit;
   TGraphErrors* meanErrorsCorrPi0Fit[ErrorPi0];
   TGraphErrors* meanErrorsCorrSummedPi0Fit;
   TGraphErrors* meanErrorsCorrSummedIncMatPi0Fit;
   TGraphErrors* meanErrorsCorrSummedIncMatPi0FitWithoutMat;

   TGraphErrors* negativeErrorsGamma[ErrorGammaSpec];
   TGraphErrors* negativeErrorsSummedGamma;
   TGraphErrors* positiveErrorsGamma[ErrorGammaSpec];
   TGraphErrors* positiveErrorsSummedGamma;
   TGraphErrors* negativeErrorsCorrGamma[ErrorGammaSpec];
   TGraphErrors* negativeErrorsCorrSummedGamma;
   TGraphErrors* positiveErrorsCorrGamma[ErrorGammaSpec];
   TGraphErrors* positiveErrorsCorrSummedGamma;
   TGraphErrors* meanErrorsGamma[ErrorGammaSpec];
   TGraphErrors* meanErrorsSummedGamma;
   TGraphErrors* meanErrorsCorrGamma[ErrorGammaSpec];
   TGraphErrors* meanErrorsCorrSummedGamma;
   TGraphErrors* meanErrorsCorrSummedIncMatGamma;
   TGraphErrors* meanErrorsCorrSummedWithoutMatGamma;
   TGraphErrors* meanErrorsCorrSummedMaterialGamma;

   TGraphErrors* negativeErrorsGammaA[ErrorGammaSpec];
   TGraphErrors* negativeErrorsSummedGammaA;
   TGraphErrors* positiveErrorsGammaA[ErrorGammaSpec];
   TGraphErrors* positiveErrorsSummedGammaA;
   TGraphErrors* negativeErrorsCorrGammaA[ErrorGammaSpec];
   TGraphErrors* negativeErrorsCorrSummedGammaA;
   TGraphErrors* positiveErrorsCorrGammaA[ErrorGammaSpec];
   TGraphErrors* positiveErrorsCorrSummedGammaA;
   TGraphErrors* meanErrorsGammaA[ErrorGammaSpec];
   TGraphErrors* meanErrorsSummedGammaA;
   TGraphErrors* meanErrorsCorrGammaA[ErrorGammaSpec];
   TGraphErrors* meanErrorsCorrSummedGammaA;
   TGraphErrors* meanErrorsCorrSummedIncMatGammaA;

   TGraphErrors* negativeErrorsGammaB[ErrorGammaSpec];
   TGraphErrors* negativeErrorsSummedGammaB;
   TGraphErrors* positiveErrorsGammaB[ErrorGammaSpec];
   TGraphErrors* positiveErrorsSummedGammaB;
   TGraphErrors* negativeErrorsCorrGammaB[ErrorGammaSpec];
   TGraphErrors* negativeErrorsCorrSummedGammaB;
   TGraphErrors* positiveErrorsCorrGammaB[ErrorGammaSpec];
   TGraphErrors* positiveErrorsCorrSummedGammaB;
   TGraphErrors* meanErrorsGammaB[ErrorGammaSpec];
   TGraphErrors* meanErrorsSummedGammaB;
   TGraphErrors* meanErrorsCorrGammaB[ErrorGammaSpec];
   TGraphErrors* meanErrorsCorrSummedGammaB;
   TGraphErrors* meanErrorsCorrSummedIncMatGammaB;

   TGraphErrors* negativeErrorsGammaC[ErrorGammaSpec];
   TGraphErrors* negativeErrorsSummedGammaC;
   TGraphErrors* positiveErrorsGammaC[ErrorGammaSpec];
   TGraphErrors* positiveErrorsSummedGammaC;
   TGraphErrors* negativeErrorsCorrGammaC[ErrorGammaSpec];
   TGraphErrors* negativeErrorsCorrSummedGammaC;
   TGraphErrors* positiveErrorsCorrGammaC[ErrorGammaSpec];
   TGraphErrors* positiveErrorsCorrSummedGammaC;
   TGraphErrors* meanErrorsGammaC[ErrorGammaSpec];
   TGraphErrors* meanErrorsSummedGammaC;
   TGraphErrors* meanErrorsCorrGammaC[ErrorGammaSpec];
   TGraphErrors* meanErrorsCorrSummedGammaC;
   TGraphErrors* meanErrorsCorrSummedIncMatGammaC;



   Int_t nBinsGraphDoubeRatio = 0;
   Int_t nBinsGraphGamma = 0;
   for(Int_t i = 0; i<ErrorDoubleRatio; i++){
      cout<<"------------------------------------------------>  "<<nameCutVariationsDoubleRatio[i]<<endl;

      graphPosErrorsDoubleRatio[i] = (TGraphAsymmErrors*)fileSystematics->Get(Form("DoubleRatio_SystErrorRelPos_%s_%s",nameCutVariationsDoubleRatio[i].Data(),(GetCentralityStringA(cutSel)).Data()));
      graphNegErrorsDoubleRatio[i] = (TGraphAsymmErrors*)fileSystematics->Get(Form("DoubleRatio_SystErrorRelNeg_%s_%s",nameCutVariationsDoubleRatio[i].Data(),(GetCentralityStringA(cutSel)).Data()));
      if(!nameCutVariationsDoubleRatio[i].CompareTo("CocktailEtaNorm")){
         graphPosErrorsDoubleRatio[i] = (TGraphAsymmErrors*)fileSystematics->Get(Form("DoubleRatio_SystErrorRelPos_%s_%s",nameCutVariationsDoubleRatio[i].Data(),(GetCentralityStringB(cutSel)).Data()));
         graphNegErrorsDoubleRatio[i] = (TGraphAsymmErrors*)fileSystematics->Get(Form("DoubleRatio_SystErrorRelNeg_%s_%s",nameCutVariationsDoubleRatio[i].Data(),(GetCentralityStringB(cutSel)).Data()));
         for(Int_t jj = 0;jj<graphPosErrorsDoubleRatio[i]->GetN();jj++){
            graphPosErrorsDoubleRatio[i]->SetPoint(jj, graphPosErrorsDoubleRatio[i]->GetX()[jj], 1.4);
            graphNegErrorsDoubleRatio[i]->SetPoint(jj, graphNegErrorsDoubleRatio[i]->GetX()[jj], 1.4);
         }
      }

      graphPosErrorsDoubleRatio[i]->RemovePoint(0);
      graphNegErrorsDoubleRatio[i]->RemovePoint(0);
      graphPosErrorsDoubleRatio[i]->RemovePoint(0);
      graphNegErrorsDoubleRatio[i]->RemovePoint(0);

      graphPosErrorsDoubleRatioPi0Fit[i] = (TGraphAsymmErrors*)fileSystematics->Get(Form("DoubleRatioFit_SystErrorRelPos_%s_%s",nameCutVariationsDoubleRatio[i].Data(),(GetCentralityStringA(cutSel)).Data()));
      graphNegErrorsDoubleRatioPi0Fit[i] = (TGraphAsymmErrors*)fileSystematics->Get(Form("DoubleRatioFit_SystErrorRelNeg_%s_%s",nameCutVariationsDoubleRatio[i].Data(),(GetCentralityStringA(cutSel)).Data()));
      if(!nameCutVariationsDoubleRatio[i].CompareTo("CocktailEtaNorm")){
         graphPosErrorsDoubleRatioPi0Fit[i] = (TGraphAsymmErrors*)fileSystematics->Get(Form("DoubleRatioFit_SystErrorRelPos_%s_%s",nameCutVariationsDoubleRatio[i].Data(),(GetCentralityStringB(cutSel)).Data()));
         graphNegErrorsDoubleRatioPi0Fit[i] = (TGraphAsymmErrors*)fileSystematics->Get(Form("DoubleRatioFit_SystErrorRelNeg_%s_%s",nameCutVariationsDoubleRatio[i].Data(),(GetCentralityStringB(cutSel)).Data()));
         for(Int_t jj = 0;jj<graphPosErrorsDoubleRatioPi0Fit[i]->GetN();jj++){
            graphPosErrorsDoubleRatioPi0Fit[i]->SetPoint(jj, graphPosErrorsDoubleRatioPi0Fit[i]->GetX()[jj], 1.4);
            graphNegErrorsDoubleRatioPi0Fit[i]->SetPoint(jj, graphNegErrorsDoubleRatioPi0Fit[i]->GetX()[jj], 1.4);
         }
      }

      graphPosErrorsDoubleRatioPi0Fit[i]->RemovePoint(0);
      graphNegErrorsDoubleRatioPi0Fit[i]->RemovePoint(0);
      graphPosErrorsDoubleRatioPi0Fit[i]->RemovePoint(0);
      graphNegErrorsDoubleRatioPi0Fit[i]->RemovePoint(0);

      Int_t type = 0;
      if(!nameCutVariationsDoubleRatio[i].CompareTo("Chi2") || !nameCutVariationsDoubleRatio[i].CompareTo("edEdx") || !nameCutVariationsDoubleRatio[i].CompareTo("PsiPair")  || !nameCutVariationsDoubleRatio[i].CompareTo("TPCClst") || !nameCutVariationsDoubleRatio[i].CompareTo("CocktailParam") || !nameCutVariationsDoubleRatio[i].CompareTo("CocktailEta") || !nameCutVariationsDoubleRatio[i].CompareTo("TOF") || !nameCutVariationsDoubleRatio[i].CompareTo("SinglePt") || !nameCutVariationsDoubleRatio[i].CompareTo("IntRange") || !nameCutVariationsDoubleRatio[i].CompareTo("Alpha") ){ // A Errors
         type = 1;
      }
      if(!nameCutVariationsDoubleRatio[i].CompareTo("qT") || !nameCutVariationsDoubleRatio[i].CompareTo("PidEdx")){ // B Errors
         type = 2;
      }
      if(!nameCutVariationsDoubleRatio[i].CompareTo("CocktailEtaNorm")){ // C Errors
         type = 3;
      }

      if(type == 1){ // A Errors
         graphPosErrorsDoubleRatioA[i] = (TGraphAsymmErrors*)fileSystematics->Get(Form("DoubleRatio_SystErrorRelPos_%s_%s",nameCutVariationsDoubleRatio[i].Data(),(GetCentralityStringA(cutSel)).Data()));
         graphNegErrorsDoubleRatioA[i] = (TGraphAsymmErrors*)fileSystematics->Get(Form("DoubleRatio_SystErrorRelNeg_%s_%s",nameCutVariationsDoubleRatio[i].Data(),(GetCentralityStringA(cutSel)).Data()));

         graphPosErrorsDoubleRatioA[i]->RemovePoint(0);
         graphNegErrorsDoubleRatioA[i]->RemovePoint(0);
         graphPosErrorsDoubleRatioA[i]->RemovePoint(0);
         graphNegErrorsDoubleRatioA[i]->RemovePoint(0);

      }
      if(type == 2){ // B Errors
         graphPosErrorsDoubleRatioB[i] = (TGraphAsymmErrors*)fileSystematics->Get(Form("DoubleRatio_SystErrorRelPos_%s_%s",nameCutVariationsDoubleRatio[i].Data(),(GetCentralityStringA(cutSel)).Data()));
         graphNegErrorsDoubleRatioB[i] = (TGraphAsymmErrors*)fileSystematics->Get(Form("DoubleRatio_SystErrorRelNeg_%s_%s",nameCutVariationsDoubleRatio[i].Data(),(GetCentralityStringA(cutSel)).Data()));

         graphPosErrorsDoubleRatioB[i]->RemovePoint(0);
         graphNegErrorsDoubleRatioB[i]->RemovePoint(0);
         graphPosErrorsDoubleRatioB[i]->RemovePoint(0);
         graphNegErrorsDoubleRatioB[i]->RemovePoint(0);
      }
      if(type == 3){ // C Errors
         graphPosErrorsDoubleRatioC[i] = (TGraphAsymmErrors*)fileSystematics->Get(Form("DoubleRatio_SystErrorRelPos_%s_%s",nameCutVariationsDoubleRatio[i].Data(),(GetCentralityStringA(cutSel)).Data()));
         graphNegErrorsDoubleRatioC[i] = (TGraphAsymmErrors*)fileSystematics->Get(Form("DoubleRatio_SystErrorRelNeg_%s_%s",nameCutVariationsDoubleRatio[i].Data(),(GetCentralityStringA(cutSel)).Data()));
         if(!nameCutVariationsDoubleRatio[i].CompareTo("CocktailEtaNorm")){
            for(Int_t jj = 0;jj<graphPosErrorsDoubleRatioC[i]->GetN();jj++){
               graphPosErrorsDoubleRatio[i]->SetPoint(jj, graphPosErrorsDoubleRatio[i]->GetX()[jj], 1.4);
               graphNegErrorsDoubleRatio[i]->SetPoint(jj, graphNegErrorsDoubleRatio[i]->GetX()[jj], 1.4);
            }
         }
         graphPosErrorsDoubleRatioC[i]->RemovePoint(0);
         graphNegErrorsDoubleRatioC[i]->RemovePoint(0);
         graphPosErrorsDoubleRatioC[i]->RemovePoint(0);
         graphNegErrorsDoubleRatioC[i]->RemovePoint(0);
      }


      if(type == 1){
         graphPosErrorsDoubleRatioPi0FitA[i] = (TGraphAsymmErrors*)fileSystematics->Get(Form("DoubleRatioFit_SystErrorRelPos_%s_%s",nameCutVariationsDoubleRatio[i].Data(),(GetCentralityStringA(cutSel)).Data()));
         graphNegErrorsDoubleRatioPi0FitA[i] = (TGraphAsymmErrors*)fileSystematics->Get(Form("DoubleRatioFit_SystErrorRelNeg_%s_%s",nameCutVariationsDoubleRatio[i].Data(),(GetCentralityStringA(cutSel)).Data()));

         graphPosErrorsDoubleRatioPi0FitA[i]->RemovePoint(0);
         graphNegErrorsDoubleRatioPi0FitA[i]->RemovePoint(0);
         graphPosErrorsDoubleRatioPi0FitA[i]->RemovePoint(0);
         graphNegErrorsDoubleRatioPi0FitA[i]->RemovePoint(0);
      }
      if(type == 2){
         graphPosErrorsDoubleRatioPi0FitB[i] = (TGraphAsymmErrors*)fileSystematics->Get(Form("DoubleRatioFit_SystErrorRelPos_%s_%s",nameCutVariationsDoubleRatio[i].Data(),(GetCentralityStringA(cutSel)).Data()));
         graphNegErrorsDoubleRatioPi0FitB[i] = (TGraphAsymmErrors*)fileSystematics->Get(Form("DoubleRatioFit_SystErrorRelNeg_%s_%s",nameCutVariationsDoubleRatio[i].Data(),(GetCentralityStringA(cutSel)).Data()));

         graphPosErrorsDoubleRatioPi0FitB[i]->RemovePoint(0);
         graphNegErrorsDoubleRatioPi0FitB[i]->RemovePoint(0);
         graphPosErrorsDoubleRatioPi0FitB[i]->RemovePoint(0);
         graphNegErrorsDoubleRatioPi0FitB[i]->RemovePoint(0);
      }
      if(type == 3){
         graphPosErrorsDoubleRatioPi0FitC[i] = (TGraphAsymmErrors*)fileSystematics->Get(Form("DoubleRatioFit_SystErrorRelPos_%s_%s",nameCutVariationsDoubleRatio[i].Data(),(GetCentralityStringB(cutSel)).Data()));
         graphNegErrorsDoubleRatioPi0FitC[i] = (TGraphAsymmErrors*)fileSystematics->Get(Form("DoubleRatioFit_SystErrorRelNeg_%s_%s",nameCutVariationsDoubleRatio[i].Data(),(GetCentralityStringB(cutSel)).Data()));
         if(!nameCutVariationsDoubleRatio[i].CompareTo("CocktailEtaNorm")){
            for(Int_t jj = 0;jj<graphPosErrorsDoubleRatioPi0FitC[i]->GetN();jj++){
               graphPosErrorsDoubleRatioPi0Fit[i]->SetPoint(jj, graphPosErrorsDoubleRatioPi0Fit[i]->GetX()[jj], 1.4);
               graphNegErrorsDoubleRatioPi0Fit[i]->SetPoint(jj, graphNegErrorsDoubleRatioPi0Fit[i]->GetX()[jj], 1.4);
            }
         }
         graphPosErrorsDoubleRatioPi0FitC[i]->RemovePoint(0);
         graphNegErrorsDoubleRatioPi0FitC[i]->RemovePoint(0);
         graphPosErrorsDoubleRatioPi0FitC[i]->RemovePoint(0);
         graphNegErrorsDoubleRatioPi0FitC[i]->RemovePoint(0);
      }

      if(i==0){
         ptBinsDoubleRatio = graphPosErrorsDoubleRatio[i]->GetX();
         nBinsGraphDoubeRatio = graphPosErrorsDoubleRatio[i]->GetN();
         ptBinsDoubleRatioErr = graphPosErrorsDoubleRatio[i]->GetEXhigh();
      }

      errorsNegDoubleRatio[i] = graphNegErrorsDoubleRatio[i]->GetY();
      errorsNegErrDoubleRatio[i] = graphNegErrorsDoubleRatio[i]->GetEYhigh();
      errorsPosDoubleRatio[i] = graphPosErrorsDoubleRatio[i]->GetY();
      errorsPosErrDoubleRatio[i] = graphPosErrorsDoubleRatio[i]->GetEYhigh();

      errorsPosDoubleRatioPi0Fit[i] = graphPosErrorsDoubleRatioPi0Fit[i]->GetY();
      errorsPosErrDoubleRatioPi0Fit[i] = graphPosErrorsDoubleRatioPi0Fit[i]->GetEYhigh();
      errorsNegDoubleRatioPi0Fit[i] = graphNegErrorsDoubleRatioPi0Fit[i]->GetY();
      errorsNegErrDoubleRatioPi0Fit[i] = graphNegErrorsDoubleRatioPi0Fit[i]->GetEYhigh();

      CalculateMeanSysErr(errorsMeanDoubleRatio[i], errorsMeanErrDoubleRatio[i], errorsPosDoubleRatio[i], errorsNegDoubleRatio[i], ConstnBinsDoubleRatio);
      CalculateMeanSysErr(errorsMeanDoubleRatioPi0Fit[i], errorsMeanErrDoubleRatioPi0Fit[i], errorsPosDoubleRatioPi0Fit[i], errorsNegDoubleRatioPi0Fit[i], ConstnBinsDoubleRatio);

      CorrectSystematicErrorsWithMean(errorsPosDoubleRatio[i],errorsPosErrDoubleRatio[i], errorsPosCorrDoubleRatio[i], errorsPosErrCorrDoubleRatio[i], ConstnBinsDoubleRatio);
      CorrectSystematicErrorsWithMean(errorsNegDoubleRatio[i],errorsNegErrDoubleRatio[i], errorsNegCorrDoubleRatio[i], errorsNegErrCorrDoubleRatio[i], ConstnBinsDoubleRatio);
      CorrectSystematicErrorsWithMean(errorsMeanDoubleRatio[i], errorsMeanErrDoubleRatio[i], errorsMeanCorrDoubleRatio[i], errorsMeanErrCorrDoubleRatio[i], ConstnBinsDoubleRatio);

      CorrectSystematicErrorsWithMean(errorsPosDoubleRatioPi0Fit[i],errorsPosErrDoubleRatioPi0Fit[i], errorsPosCorrDoubleRatioPi0Fit[i], errorsPosErrCorrDoubleRatioPi0Fit[i], ConstnBinsDoubleRatio);
      CorrectSystematicErrorsWithMean(errorsNegDoubleRatioPi0Fit[i],errorsNegErrDoubleRatioPi0Fit[i], errorsNegCorrDoubleRatioPi0Fit[i], errorsNegErrCorrDoubleRatioPi0Fit[i], ConstnBinsDoubleRatio);
      CorrectSystematicErrorsWithMean(errorsMeanDoubleRatioPi0Fit[i], errorsMeanErrDoubleRatioPi0Fit[i], errorsMeanCorrDoubleRatioPi0Fit[i], errorsMeanErrCorrDoubleRatioPi0Fit[i], ConstnBinsDoubleRatio);


      negativeErrorsDoubleRatio[i] = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsNegDoubleRatio[i] ,ptBinsDoubleRatioErr ,errorsNegErrDoubleRatio[i] );
      meanErrorsDoubleRatio[i] = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanDoubleRatio[i] ,ptBinsDoubleRatioErr ,errorsMeanErrDoubleRatio[i] );
      positiveErrorsDoubleRatio[i] = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsPosDoubleRatio[i] ,ptBinsDoubleRatioErr ,errorsPosErrDoubleRatio[i] );
      negativeErrorsCorrDoubleRatio[i] = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsNegCorrDoubleRatio[i] ,ptBinsDoubleRatioErr ,errorsNegErrCorrDoubleRatio[i] );
      meanErrorsCorrDoubleRatio[i] = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanCorrDoubleRatio[i] ,ptBinsDoubleRatioErr ,errorsMeanErrCorrDoubleRatio[i] );
      positiveErrorsCorrDoubleRatio[i] = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsPosCorrDoubleRatio[i] ,ptBinsDoubleRatioErr ,errorsPosErrCorrDoubleRatio[i] );

      negativeErrorsDoubleRatioPi0Fit[i] = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsNegDoubleRatioPi0Fit[i] ,ptBinsDoubleRatioErr ,errorsNegErrDoubleRatioPi0Fit[i] );
      meanErrorsDoubleRatioPi0Fit[i] = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanDoubleRatioPi0Fit[i] ,ptBinsDoubleRatioErr ,errorsMeanErrDoubleRatioPi0Fit[i] );
      positiveErrorsDoubleRatioPi0Fit[i] = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsPosDoubleRatioPi0Fit[i] ,ptBinsDoubleRatioErr ,errorsPosErrDoubleRatioPi0Fit[i] );
      negativeErrorsCorrDoubleRatioPi0Fit[i] = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsNegCorrDoubleRatioPi0Fit[i] ,ptBinsDoubleRatioErr ,errorsNegErrCorrDoubleRatioPi0Fit[i] );
      meanErrorsCorrDoubleRatioPi0Fit[i] = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanCorrDoubleRatioPi0Fit[i] ,ptBinsDoubleRatioErr ,errorsMeanErrCorrDoubleRatioPi0Fit[i] );
      positiveErrorsCorrDoubleRatioPi0Fit[i] = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsPosCorrDoubleRatioPi0Fit[i] ,ptBinsDoubleRatioErr ,errorsPosErrCorrDoubleRatioPi0Fit[i] );


      if(type==1){
         errorsNegDoubleRatioA[i] = graphNegErrorsDoubleRatioA[i]->GetY();
         errorsNegErrDoubleRatioA[i] = graphNegErrorsDoubleRatioA[i]->GetEYhigh();
         errorsPosDoubleRatioA[i] = graphPosErrorsDoubleRatioA[i]->GetY();
         errorsPosErrDoubleRatioA[i] = graphPosErrorsDoubleRatioA[i]->GetEYhigh();

         errorsPosDoubleRatioPi0FitA[i] = graphPosErrorsDoubleRatioPi0FitA[i]->GetY();
         errorsPosErrDoubleRatioPi0FitA[i] = graphPosErrorsDoubleRatioPi0FitA[i]->GetEYhigh();
         errorsNegDoubleRatioPi0FitA[i] = graphNegErrorsDoubleRatioPi0FitA[i]->GetY();
         errorsNegErrDoubleRatioPi0FitA[i] = graphNegErrorsDoubleRatioPi0FitA[i]->GetEYhigh();

         CalculateMeanSysErr(errorsMeanDoubleRatioA[i], errorsMeanErrDoubleRatioA[i], errorsPosDoubleRatioA[i], errorsNegDoubleRatioA[i], ConstnBinsDoubleRatio);
         CalculateMeanSysErr(errorsMeanDoubleRatioPi0FitA[i], errorsMeanErrDoubleRatioPi0FitA[i], errorsPosDoubleRatioPi0FitA[i], errorsNegDoubleRatioPi0FitA[i], ConstnBinsDoubleRatio);

         CorrectSystematicErrorsWithMean(errorsPosDoubleRatioA[i],errorsPosErrDoubleRatioA[i], errorsPosCorrDoubleRatioA[i], errorsPosErrCorrDoubleRatioA[i], ConstnBinsDoubleRatio);
         CorrectSystematicErrorsWithMean(errorsNegDoubleRatioA[i],errorsNegErrDoubleRatioA[i], errorsNegCorrDoubleRatioA[i], errorsNegErrCorrDoubleRatioA[i], ConstnBinsDoubleRatio);
         CorrectSystematicErrorsWithMean(errorsMeanDoubleRatioA[i], errorsMeanErrDoubleRatioA[i], errorsMeanCorrDoubleRatioA[i], errorsMeanErrCorrDoubleRatioA[i], ConstnBinsDoubleRatio);

         CorrectSystematicErrorsWithMean(errorsPosDoubleRatioPi0FitA[i],errorsPosErrDoubleRatioPi0FitA[i], errorsPosCorrDoubleRatioPi0FitA[i], errorsPosErrCorrDoubleRatioPi0FitA[i], ConstnBinsDoubleRatio);
         CorrectSystematicErrorsWithMean(errorsNegDoubleRatioPi0FitA[i],errorsNegErrDoubleRatioPi0FitA[i], errorsNegCorrDoubleRatioPi0FitA[i], errorsNegErrCorrDoubleRatioPi0FitA[i], ConstnBinsDoubleRatio);
         CorrectSystematicErrorsWithMean(errorsMeanDoubleRatioPi0FitA[i], errorsMeanErrDoubleRatioPi0FitA[i], errorsMeanCorrDoubleRatioPi0FitA[i], errorsMeanErrCorrDoubleRatioPi0FitA[i], ConstnBinsDoubleRatio);

         negativeErrorsDoubleRatioA[i] = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsNegDoubleRatioA[i] ,ptBinsDoubleRatioErr ,errorsNegErrDoubleRatioA[i] );
         meanErrorsDoubleRatioA[i] = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanDoubleRatioA[i] ,ptBinsDoubleRatioErr ,errorsMeanErrDoubleRatioA[i] );
         positiveErrorsDoubleRatioA[i] = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsPosDoubleRatioA[i] ,ptBinsDoubleRatioErr ,errorsPosErrDoubleRatioA[i] );
         negativeErrorsCorrDoubleRatioA[i] = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsNegCorrDoubleRatioA[i] ,ptBinsDoubleRatioErr ,errorsNegErrCorrDoubleRatioA[i] );
         meanErrorsCorrDoubleRatioA[i] = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanCorrDoubleRatioA[i] ,ptBinsDoubleRatioErr ,errorsMeanErrCorrDoubleRatioA[i] );
         positiveErrorsCorrDoubleRatioA[i] = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsPosCorrDoubleRatioA[i] ,ptBinsDoubleRatioErr ,errorsPosErrCorrDoubleRatioA[i] );

         negativeErrorsDoubleRatioPi0FitA[i] = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsNegDoubleRatioPi0FitA[i] ,ptBinsDoubleRatioErr ,errorsNegErrDoubleRatioPi0FitA[i] );
         meanErrorsDoubleRatioPi0FitA[i] = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanDoubleRatioPi0FitA[i] ,ptBinsDoubleRatioErr ,errorsMeanErrDoubleRatioPi0FitA[i] );
         positiveErrorsDoubleRatioPi0FitA[i] = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsPosDoubleRatioPi0FitA[i] ,ptBinsDoubleRatioErr ,errorsPosErrDoubleRatioPi0FitA[i] );
         negativeErrorsCorrDoubleRatioPi0FitA[i] = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsNegCorrDoubleRatioPi0FitA[i] ,ptBinsDoubleRatioErr ,errorsNegErrCorrDoubleRatioPi0FitA[i] );
         meanErrorsCorrDoubleRatioPi0FitA[i] = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanCorrDoubleRatioPi0FitA[i] ,ptBinsDoubleRatioErr ,errorsMeanErrCorrDoubleRatioPi0FitA[i] );
         positiveErrorsCorrDoubleRatioPi0FitA[i] = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsPosCorrDoubleRatioPi0FitA[i] ,ptBinsDoubleRatioErr ,errorsPosErrCorrDoubleRatioPi0FitA[i] );
      }
      if(type==2){
         errorsNegDoubleRatioB[i] = graphNegErrorsDoubleRatioB[i]->GetY();
         errorsNegErrDoubleRatioB[i] = graphNegErrorsDoubleRatioB[i]->GetEYhigh();
         errorsPosDoubleRatioB[i] = graphPosErrorsDoubleRatioB[i]->GetY();
         errorsPosErrDoubleRatioB[i] = graphPosErrorsDoubleRatioB[i]->GetEYhigh();

         errorsPosDoubleRatioPi0FitB[i] = graphPosErrorsDoubleRatioPi0FitB[i]->GetY();
         errorsPosErrDoubleRatioPi0FitB[i] = graphPosErrorsDoubleRatioPi0FitB[i]->GetEYhigh();
         errorsNegDoubleRatioPi0FitB[i] = graphNegErrorsDoubleRatioPi0FitB[i]->GetY();
         errorsNegErrDoubleRatioPi0FitB[i] = graphNegErrorsDoubleRatioPi0FitB[i]->GetEYhigh();

         CalculateMeanSysErr(errorsMeanDoubleRatioB[i], errorsMeanErrDoubleRatioB[i], errorsPosDoubleRatioB[i], errorsNegDoubleRatioB[i], ConstnBinsDoubleRatio);
         CalculateMeanSysErr(errorsMeanDoubleRatioPi0FitB[i], errorsMeanErrDoubleRatioPi0FitB[i], errorsPosDoubleRatioPi0FitB[i], errorsNegDoubleRatioPi0FitB[i], ConstnBinsDoubleRatio);

         CorrectSystematicErrorsWithMean(errorsPosDoubleRatioB[i],errorsPosErrDoubleRatioB[i], errorsPosCorrDoubleRatioB[i], errorsPosErrCorrDoubleRatioB[i], ConstnBinsDoubleRatio);
         CorrectSystematicErrorsWithMean(errorsNegDoubleRatioB[i],errorsNegErrDoubleRatioB[i], errorsNegCorrDoubleRatioB[i], errorsNegErrCorrDoubleRatioB[i], ConstnBinsDoubleRatio);
         CorrectSystematicErrorsWithMean(errorsMeanDoubleRatioB[i], errorsMeanErrDoubleRatioB[i], errorsMeanCorrDoubleRatioB[i], errorsMeanErrCorrDoubleRatioB[i], ConstnBinsDoubleRatio);

         CorrectSystematicErrorsWithMean(errorsPosDoubleRatioPi0FitB[i],errorsPosErrDoubleRatioPi0FitB[i], errorsPosCorrDoubleRatioPi0FitB[i], errorsPosErrCorrDoubleRatioPi0FitB[i], ConstnBinsDoubleRatio);
         CorrectSystematicErrorsWithMean(errorsNegDoubleRatioPi0FitB[i],errorsNegErrDoubleRatioPi0FitB[i], errorsNegCorrDoubleRatioPi0FitB[i], errorsNegErrCorrDoubleRatioPi0FitB[i], ConstnBinsDoubleRatio);
         CorrectSystematicErrorsWithMean(errorsMeanDoubleRatioPi0FitB[i], errorsMeanErrDoubleRatioPi0FitB[i], errorsMeanCorrDoubleRatioPi0FitB[i], errorsMeanErrCorrDoubleRatioPi0FitB[i], ConstnBinsDoubleRatio);

         negativeErrorsDoubleRatioB[i] = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsNegDoubleRatioB[i] ,ptBinsDoubleRatioErr ,errorsNegErrDoubleRatioB[i] );
         meanErrorsDoubleRatioB[i] = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanDoubleRatioB[i] ,ptBinsDoubleRatioErr ,errorsMeanErrDoubleRatioB[i] );
         positiveErrorsDoubleRatioB[i] = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsPosDoubleRatioB[i] ,ptBinsDoubleRatioErr ,errorsPosErrDoubleRatioB[i] );
         negativeErrorsCorrDoubleRatioB[i] = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsNegCorrDoubleRatioB[i] ,ptBinsDoubleRatioErr ,errorsNegErrCorrDoubleRatioB[i] );
         meanErrorsCorrDoubleRatioB[i] = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanCorrDoubleRatioB[i] ,ptBinsDoubleRatioErr ,errorsMeanErrCorrDoubleRatioB[i] );
         positiveErrorsCorrDoubleRatioB[i] = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsPosCorrDoubleRatioB[i] ,ptBinsDoubleRatioErr ,errorsPosErrCorrDoubleRatioB[i] );

         negativeErrorsDoubleRatioPi0FitB[i] = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsNegDoubleRatioPi0FitB[i] ,ptBinsDoubleRatioErr ,errorsNegErrDoubleRatioPi0FitB[i] );
         meanErrorsDoubleRatioPi0FitB[i] = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanDoubleRatioPi0FitB[i] ,ptBinsDoubleRatioErr ,errorsMeanErrDoubleRatioPi0FitB[i] );
         positiveErrorsDoubleRatioPi0FitB[i] = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsPosDoubleRatioPi0FitB[i] ,ptBinsDoubleRatioErr ,errorsPosErrDoubleRatioPi0FitB[i] );
         negativeErrorsCorrDoubleRatioPi0FitB[i] = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsNegCorrDoubleRatioPi0FitB[i] ,ptBinsDoubleRatioErr ,errorsNegErrCorrDoubleRatioPi0FitB[i] );
         meanErrorsCorrDoubleRatioPi0FitB[i] = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanCorrDoubleRatioPi0FitB[i] ,ptBinsDoubleRatioErr ,errorsMeanErrCorrDoubleRatioPi0FitB[i] );
         positiveErrorsCorrDoubleRatioPi0FitB[i] = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsPosCorrDoubleRatioPi0FitB[i] ,ptBinsDoubleRatioErr ,errorsPosErrCorrDoubleRatioPi0FitB[i] );
      }
      if(type==3){
         errorsNegDoubleRatioC[i] = graphNegErrorsDoubleRatioC[i]->GetY();
         errorsNegErrDoubleRatioC[i] = graphNegErrorsDoubleRatioC[i]->GetEYhigh();
         errorsPosDoubleRatioC[i] = graphPosErrorsDoubleRatioC[i]->GetY();
         errorsPosErrDoubleRatioC[i] = graphPosErrorsDoubleRatioC[i]->GetEYhigh();

         errorsPosDoubleRatioPi0FitC[i] = graphPosErrorsDoubleRatioPi0FitC[i]->GetY();
         errorsPosErrDoubleRatioPi0FitC[i] = graphPosErrorsDoubleRatioPi0FitC[i]->GetEYhigh();
         errorsNegDoubleRatioPi0FitC[i] = graphNegErrorsDoubleRatioPi0FitC[i]->GetY();
         errorsNegErrDoubleRatioPi0FitC[i] = graphNegErrorsDoubleRatioPi0FitC[i]->GetEYhigh();

         CalculateMeanSysErr(errorsMeanDoubleRatioC[i], errorsMeanErrDoubleRatioC[i], errorsPosDoubleRatioC[i], errorsNegDoubleRatioC[i], ConstnBinsDoubleRatio);
         CalculateMeanSysErr(errorsMeanDoubleRatioPi0FitC[i], errorsMeanErrDoubleRatioPi0FitC[i], errorsPosDoubleRatioPi0FitC[i], errorsNegDoubleRatioPi0FitC[i], ConstnBinsDoubleRatio);

         CorrectSystematicErrorsWithMean(errorsPosDoubleRatioC[i],errorsPosErrDoubleRatioC[i], errorsPosCorrDoubleRatioC[i], errorsPosErrCorrDoubleRatioC[i], ConstnBinsDoubleRatio);
         CorrectSystematicErrorsWithMean(errorsNegDoubleRatioC[i],errorsNegErrDoubleRatioC[i], errorsNegCorrDoubleRatioC[i], errorsNegErrCorrDoubleRatioC[i], ConstnBinsDoubleRatio);
         CorrectSystematicErrorsWithMean(errorsMeanDoubleRatioC[i], errorsMeanErrDoubleRatioC[i], errorsMeanCorrDoubleRatioC[i], errorsMeanErrCorrDoubleRatioC[i], ConstnBinsDoubleRatio);

         CorrectSystematicErrorsWithMean(errorsPosDoubleRatioPi0FitC[i],errorsPosErrDoubleRatioPi0FitC[i], errorsPosCorrDoubleRatioPi0FitC[i], errorsPosErrCorrDoubleRatioPi0FitC[i], ConstnBinsDoubleRatio);
         CorrectSystematicErrorsWithMean(errorsNegDoubleRatioPi0FitC[i],errorsNegErrDoubleRatioPi0FitC[i], errorsNegCorrDoubleRatioPi0FitC[i], errorsNegErrCorrDoubleRatioPi0FitC[i], ConstnBinsDoubleRatio);
         CorrectSystematicErrorsWithMean(errorsMeanDoubleRatioPi0FitC[i], errorsMeanErrDoubleRatioPi0FitC[i], errorsMeanCorrDoubleRatioPi0FitC[i], errorsMeanErrCorrDoubleRatioPi0FitC[i], ConstnBinsDoubleRatio);

         negativeErrorsDoubleRatioC[i] = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsNegDoubleRatioC[i] ,ptBinsDoubleRatioErr ,errorsNegErrDoubleRatioC[i] );
         meanErrorsDoubleRatioC[i] = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanDoubleRatioC[i] ,ptBinsDoubleRatioErr ,errorsMeanErrDoubleRatioC[i] );
         positiveErrorsDoubleRatioC[i] = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsPosDoubleRatioC[i] ,ptBinsDoubleRatioErr ,errorsPosErrDoubleRatioC[i] );
         negativeErrorsCorrDoubleRatioC[i] = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsNegCorrDoubleRatioC[i] ,ptBinsDoubleRatioErr ,errorsNegErrCorrDoubleRatioC[i] );
         meanErrorsCorrDoubleRatioC[i] = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanCorrDoubleRatioC[i] ,ptBinsDoubleRatioErr ,errorsMeanErrCorrDoubleRatioC[i] );
         positiveErrorsCorrDoubleRatioC[i] = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsPosCorrDoubleRatioC[i] ,ptBinsDoubleRatioErr ,errorsPosErrCorrDoubleRatioC[i] );

         negativeErrorsDoubleRatioPi0FitC[i] = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsNegDoubleRatioPi0FitC[i] ,ptBinsDoubleRatioErr ,errorsNegErrDoubleRatioPi0FitC[i] );
         meanErrorsDoubleRatioPi0FitC[i] = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanDoubleRatioPi0FitC[i] ,ptBinsDoubleRatioErr ,errorsMeanErrDoubleRatioPi0FitC[i] );
         positiveErrorsDoubleRatioPi0FitC[i] = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsPosDoubleRatioPi0FitC[i] ,ptBinsDoubleRatioErr ,errorsPosErrDoubleRatioPi0FitC[i] );
         negativeErrorsCorrDoubleRatioPi0FitC[i] = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsNegCorrDoubleRatioPi0FitC[i] ,ptBinsDoubleRatioErr ,errorsNegErrCorrDoubleRatioPi0FitC[i] );
         meanErrorsCorrDoubleRatioPi0FitC[i] = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanCorrDoubleRatioPi0FitC[i] ,ptBinsDoubleRatioErr ,errorsMeanErrCorrDoubleRatioPi0FitC[i] );
         positiveErrorsCorrDoubleRatioPi0FitC[i] = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsPosCorrDoubleRatioPi0FitC[i] ,ptBinsDoubleRatioErr ,errorsPosErrCorrDoubleRatioPi0FitC[i] );
      }

      if(i<ErrorIncRatio){
         graphPosErrorsIncRatio[i] = (TGraphAsymmErrors*)fileSystematics->Get(Form("IncRatio_SystErrorRelPos_%s_%s",nameCutVariationsDoubleRatio[i].Data(),(GetCentralityStringA(cutSel)).Data()));
         graphNegErrorsIncRatio[i] = (TGraphAsymmErrors*)fileSystematics->Get(Form("IncRatio_SystErrorRelNeg_%s_%s",nameCutVariationsDoubleRatio[i].Data(),(GetCentralityStringA(cutSel)).Data()));

         graphPosErrorsIncRatio[i]->RemovePoint(0);
         graphNegErrorsIncRatio[i]->RemovePoint(0);
         graphPosErrorsIncRatio[i]->RemovePoint(0);
         graphNegErrorsIncRatio[i]->RemovePoint(0);

         graphPosErrorsIncRatioPi0Fit[i] = (TGraphAsymmErrors*)fileSystematics->Get(Form("IncRatioFit_SystErrorRelPos_%s_%s",nameCutVariationsDoubleRatio[i].Data(),(GetCentralityStringA(cutSel)).Data()));
         graphNegErrorsIncRatioPi0Fit[i] = (TGraphAsymmErrors*)fileSystematics->Get(Form("IncRatioFit_SystErrorRelNeg_%s_%s",nameCutVariationsDoubleRatio[i].Data(),(GetCentralityStringA(cutSel)).Data()));

         graphPosErrorsIncRatioPi0Fit[i]->RemovePoint(0);
         graphNegErrorsIncRatioPi0Fit[i]->RemovePoint(0);
         graphPosErrorsIncRatioPi0Fit[i]->RemovePoint(0);
         graphNegErrorsIncRatioPi0Fit[i]->RemovePoint(0);

         if(type == 1){ // A Errors
            graphPosErrorsIncRatioA[i] = (TGraphAsymmErrors*)fileSystematics->Get(Form("IncRatio_SystErrorRelPos_%s_%s",nameCutVariationsIncRatio[i].Data(),(GetCentralityStringA(cutSel)).Data()));
            graphNegErrorsIncRatioA[i] = (TGraphAsymmErrors*)fileSystematics->Get(Form("IncRatio_SystErrorRelNeg_%s_%s",nameCutVariationsIncRatio[i].Data(),(GetCentralityStringA(cutSel)).Data()));

            graphPosErrorsIncRatioA[i]->RemovePoint(0);
            graphNegErrorsIncRatioA[i]->RemovePoint(0);
            graphPosErrorsIncRatioA[i]->RemovePoint(0);
            graphNegErrorsIncRatioA[i]->RemovePoint(0);
         }
         if(type == 2){ // B Errors
            graphPosErrorsIncRatioB[i] = (TGraphAsymmErrors*)fileSystematics->Get(Form("IncRatio_SystErrorRelPos_%s_%s",nameCutVariationsIncRatio[i].Data(),(GetCentralityStringA(cutSel)).Data()));
            graphNegErrorsIncRatioB[i] = (TGraphAsymmErrors*)fileSystematics->Get(Form("IncRatio_SystErrorRelNeg_%s_%s",nameCutVariationsIncRatio[i].Data(),(GetCentralityStringA(cutSel)).Data()));

            graphPosErrorsIncRatioB[i]->RemovePoint(0);
            graphNegErrorsIncRatioB[i]->RemovePoint(0);
            graphPosErrorsIncRatioB[i]->RemovePoint(0);
            graphNegErrorsIncRatioB[i]->RemovePoint(0);
         }
         if(type == 3){ // C Errors
            graphPosErrorsIncRatioC[i] = (TGraphAsymmErrors*)fileSystematics->Get(Form("IncRatio_SystErrorRelPos_%s_%s",nameCutVariationsIncRatio[i].Data(),(GetCentralityStringB(cutSel)).Data()));
            graphNegErrorsIncRatioC[i] = (TGraphAsymmErrors*)fileSystematics->Get(Form("IncRatio_SystErrorRelNeg_%s_%s",nameCutVariationsIncRatio[i].Data(),(GetCentralityStringB(cutSel)).Data()));

            graphPosErrorsIncRatioC[i]->RemovePoint(0);
            graphNegErrorsIncRatioC[i]->RemovePoint(0);
            graphPosErrorsIncRatioC[i]->RemovePoint(0);
            graphNegErrorsIncRatioC[i]->RemovePoint(0);
         }
         if(type == 1){
            graphPosErrorsIncRatioPi0FitA[i] = (TGraphAsymmErrors*)fileSystematics->Get(Form("IncRatioFit_SystErrorRelPos_%s_%s",nameCutVariationsIncRatio[i].Data(),(GetCentralityStringA(cutSel)).Data()));
            graphNegErrorsIncRatioPi0FitA[i] = (TGraphAsymmErrors*)fileSystematics->Get(Form("IncRatioFit_SystErrorRelNeg_%s_%s",nameCutVariationsIncRatio[i].Data(),(GetCentralityStringA(cutSel)).Data()));

            graphPosErrorsIncRatioPi0FitA[i]->RemovePoint(0);
            graphNegErrorsIncRatioPi0FitA[i]->RemovePoint(0);
            graphPosErrorsIncRatioPi0FitA[i]->RemovePoint(0);
            graphNegErrorsIncRatioPi0FitA[i]->RemovePoint(0);
         }
         if(type == 2){
            graphPosErrorsIncRatioPi0FitB[i] = (TGraphAsymmErrors*)fileSystematics->Get(Form("IncRatioFit_SystErrorRelPos_%s_%s",nameCutVariationsIncRatio[i].Data(),(GetCentralityStringA(cutSel)).Data()));
            graphNegErrorsIncRatioPi0FitB[i] = (TGraphAsymmErrors*)fileSystematics->Get(Form("IncRatioFit_SystErrorRelNeg_%s_%s",nameCutVariationsIncRatio[i].Data(),(GetCentralityStringA(cutSel)).Data()));

            graphPosErrorsIncRatioPi0FitB[i]->RemovePoint(0);
            graphNegErrorsIncRatioPi0FitB[i]->RemovePoint(0);
            graphPosErrorsIncRatioPi0FitB[i]->RemovePoint(0);
            graphNegErrorsIncRatioPi0FitB[i]->RemovePoint(0);
         }
         if(type == 3){
            graphPosErrorsIncRatioPi0FitC[i] = (TGraphAsymmErrors*)fileSystematics->Get(Form("IncRatioFit_SystErrorRelPos_%s_%s",nameCutVariationsIncRatio[i].Data(),(GetCentralityStringB(cutSel)).Data()));
            graphNegErrorsIncRatioPi0FitC[i] = (TGraphAsymmErrors*)fileSystematics->Get(Form("IncRatioFit_SystErrorRelNeg_%s_%s",nameCutVariationsIncRatio[i].Data(),(GetCentralityStringB(cutSel)).Data()));

            graphPosErrorsIncRatioPi0FitC[i]->RemovePoint(0);
            graphNegErrorsIncRatioPi0FitC[i]->RemovePoint(0);
            graphPosErrorsIncRatioPi0FitC[i]->RemovePoint(0);
            graphNegErrorsIncRatioPi0FitC[i]->RemovePoint(0);
         }


         if(i==0){
            ptBinsIncRatio = graphPosErrorsIncRatio[i]->GetX();
            nBinsGraphDoubeRatio = graphPosErrorsIncRatio[i]->GetN();
            ptBinsIncRatioErr = graphPosErrorsIncRatio[i]->GetEXhigh();
         }

         errorsNegIncRatio[i] = graphNegErrorsIncRatio[i]->GetY();
         errorsNegErrIncRatio[i] = graphNegErrorsIncRatio[i]->GetEYhigh();
         errorsPosIncRatio[i] = graphPosErrorsIncRatio[i]->GetY();
         errorsPosErrIncRatio[i] = graphPosErrorsIncRatio[i]->GetEYhigh();

         errorsPosIncRatioPi0Fit[i] = graphPosErrorsIncRatioPi0Fit[i]->GetY();
         errorsPosErrIncRatioPi0Fit[i] = graphPosErrorsIncRatioPi0Fit[i]->GetEYhigh();
         errorsNegIncRatioPi0Fit[i] = graphNegErrorsIncRatioPi0Fit[i]->GetY();
         errorsNegErrIncRatioPi0Fit[i] = graphNegErrorsIncRatioPi0Fit[i]->GetEYhigh();


         CalculateMeanSysErr(errorsMeanIncRatio[i], errorsMeanErrIncRatio[i], errorsPosIncRatio[i], errorsNegIncRatio[i], ConstnBinsIncRatio);
         CalculateMeanSysErr(errorsMeanIncRatioPi0Fit[i], errorsMeanErrIncRatioPi0Fit[i], errorsPosIncRatioPi0Fit[i], errorsNegIncRatioPi0Fit[i], ConstnBinsIncRatio);

         CorrectSystematicErrorsWithMean(errorsPosIncRatio[i],errorsPosErrIncRatio[i], errorsPosCorrIncRatio[i], errorsPosErrCorrIncRatio[i], ConstnBinsIncRatio);
         CorrectSystematicErrorsWithMean(errorsNegIncRatio[i],errorsNegErrIncRatio[i], errorsNegCorrIncRatio[i], errorsNegErrCorrIncRatio[i], ConstnBinsIncRatio);
         CorrectSystematicErrorsWithMean(errorsMeanIncRatio[i], errorsMeanErrIncRatio[i], errorsMeanCorrIncRatio[i], errorsMeanErrCorrIncRatio[i], ConstnBinsIncRatio);

         CorrectSystematicErrorsWithMean(errorsPosIncRatioPi0Fit[i],errorsPosErrIncRatioPi0Fit[i], errorsPosCorrIncRatioPi0Fit[i], errorsPosErrCorrIncRatioPi0Fit[i], ConstnBinsIncRatio);
         CorrectSystematicErrorsWithMean(errorsNegIncRatioPi0Fit[i],errorsNegErrIncRatioPi0Fit[i], errorsNegCorrIncRatioPi0Fit[i], errorsNegErrCorrIncRatioPi0Fit[i], ConstnBinsIncRatio);
         CorrectSystematicErrorsWithMean(errorsMeanIncRatioPi0Fit[i], errorsMeanErrIncRatioPi0Fit[i], errorsMeanCorrIncRatioPi0Fit[i], errorsMeanErrCorrIncRatioPi0Fit[i], ConstnBinsIncRatio);

         negativeErrorsIncRatio[i] = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsNegIncRatio[i] ,ptBinsIncRatioErr ,errorsNegErrIncRatio[i] );
         meanErrorsIncRatio[i] = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanIncRatio[i] ,ptBinsIncRatioErr ,errorsMeanErrIncRatio[i] );
         positiveErrorsIncRatio[i] = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsPosIncRatio[i] ,ptBinsIncRatioErr ,errorsPosErrIncRatio[i] );
         negativeErrorsCorrIncRatio[i] = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsNegCorrIncRatio[i] ,ptBinsIncRatioErr ,errorsNegErrCorrIncRatio[i] );
         meanErrorsCorrIncRatio[i] = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanCorrIncRatio[i] ,ptBinsIncRatioErr ,errorsMeanErrCorrIncRatio[i] );
         positiveErrorsCorrIncRatio[i] = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsPosCorrIncRatio[i] ,ptBinsIncRatioErr ,errorsPosErrCorrIncRatio[i] );

         negativeErrorsIncRatioPi0Fit[i] = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsNegIncRatioPi0Fit[i] ,ptBinsIncRatioErr ,errorsNegErrIncRatioPi0Fit[i] );
         meanErrorsIncRatioPi0Fit[i] = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanIncRatioPi0Fit[i] ,ptBinsIncRatioErr ,errorsMeanErrIncRatioPi0Fit[i] );
         positiveErrorsIncRatioPi0Fit[i] = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsPosIncRatioPi0Fit[i] ,ptBinsIncRatioErr ,errorsPosErrIncRatioPi0Fit[i] );
         negativeErrorsCorrIncRatioPi0Fit[i] = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsNegCorrIncRatioPi0Fit[i] ,ptBinsIncRatioErr ,errorsNegErrCorrIncRatioPi0Fit[i] );
         meanErrorsCorrIncRatioPi0Fit[i] = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanCorrIncRatioPi0Fit[i] ,ptBinsIncRatioErr ,errorsMeanErrCorrIncRatioPi0Fit[i] );
         positiveErrorsCorrIncRatioPi0Fit[i] = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsPosCorrIncRatioPi0Fit[i] ,ptBinsIncRatioErr ,errorsPosErrCorrIncRatioPi0Fit[i] );

         if(type==1){
            errorsNegIncRatioA[i] = graphNegErrorsIncRatioA[i]->GetY();
            errorsNegErrIncRatioA[i] = graphNegErrorsIncRatioA[i]->GetEYhigh();
            errorsPosIncRatioA[i] = graphPosErrorsIncRatioA[i]->GetY();
            errorsPosErrIncRatioA[i] = graphPosErrorsIncRatioA[i]->GetEYhigh();

            errorsPosIncRatioPi0FitA[i] = graphPosErrorsIncRatioPi0FitA[i]->GetY();
            errorsPosErrIncRatioPi0FitA[i] = graphPosErrorsIncRatioPi0FitA[i]->GetEYhigh();
            errorsNegIncRatioPi0FitA[i] = graphNegErrorsIncRatioPi0FitA[i]->GetY();
            errorsNegErrIncRatioPi0FitA[i] = graphNegErrorsIncRatioPi0FitA[i]->GetEYhigh();

            CalculateMeanSysErr(errorsMeanIncRatioA[i], errorsMeanErrIncRatioA[i], errorsPosIncRatioA[i], errorsNegIncRatioA[i], ConstnBinsIncRatio);
            CalculateMeanSysErr(errorsMeanIncRatioPi0FitA[i], errorsMeanErrIncRatioPi0FitA[i], errorsPosIncRatioPi0FitA[i], errorsNegIncRatioPi0FitA[i], ConstnBinsIncRatio);

            CorrectSystematicErrorsWithMean(errorsPosIncRatioA[i],errorsPosErrIncRatioA[i], errorsPosCorrIncRatioA[i], errorsPosErrCorrIncRatioA[i], ConstnBinsIncRatio);
            CorrectSystematicErrorsWithMean(errorsNegIncRatioA[i],errorsNegErrIncRatioA[i], errorsNegCorrIncRatioA[i], errorsNegErrCorrIncRatioA[i], ConstnBinsIncRatio);
            CorrectSystematicErrorsWithMean(errorsMeanIncRatioA[i], errorsMeanErrIncRatioA[i], errorsMeanCorrIncRatioA[i], errorsMeanErrCorrIncRatioA[i], ConstnBinsIncRatio);

            CorrectSystematicErrorsWithMean(errorsPosIncRatioPi0FitA[i],errorsPosErrIncRatioPi0FitA[i], errorsPosCorrIncRatioPi0FitA[i], errorsPosErrCorrIncRatioPi0FitA[i], ConstnBinsIncRatio);
            CorrectSystematicErrorsWithMean(errorsNegIncRatioPi0FitA[i],errorsNegErrIncRatioPi0FitA[i], errorsNegCorrIncRatioPi0FitA[i], errorsNegErrCorrIncRatioPi0FitA[i], ConstnBinsIncRatio);
            CorrectSystematicErrorsWithMean(errorsMeanIncRatioPi0FitA[i], errorsMeanErrIncRatioPi0FitA[i], errorsMeanCorrIncRatioPi0FitA[i], errorsMeanErrCorrIncRatioPi0FitA[i], ConstnBinsIncRatio);

            negativeErrorsIncRatioA[i] = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsNegIncRatioA[i] ,ptBinsIncRatioErr ,errorsNegErrIncRatioA[i] );
            meanErrorsIncRatioA[i] = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanIncRatioA[i] ,ptBinsIncRatioErr ,errorsMeanErrIncRatioA[i] );
            positiveErrorsIncRatioA[i] = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsPosIncRatioA[i] ,ptBinsIncRatioErr ,errorsPosErrIncRatioA[i] );
            negativeErrorsCorrIncRatioA[i] = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsNegCorrIncRatioA[i] ,ptBinsIncRatioErr ,errorsNegErrCorrIncRatioA[i] );
            meanErrorsCorrIncRatioA[i] = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanCorrIncRatioA[i] ,ptBinsIncRatioErr ,errorsMeanErrCorrIncRatioA[i] );
            positiveErrorsCorrIncRatioA[i] = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsPosCorrIncRatioA[i] ,ptBinsIncRatioErr ,errorsPosErrCorrIncRatioA[i] );

            negativeErrorsIncRatioPi0FitA[i] = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsNegIncRatioPi0FitA[i] ,ptBinsIncRatioErr ,errorsNegErrIncRatioPi0FitA[i] );
            meanErrorsIncRatioPi0FitA[i] = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanIncRatioPi0FitA[i] ,ptBinsIncRatioErr ,errorsMeanErrIncRatioPi0FitA[i] );
            positiveErrorsIncRatioPi0FitA[i] = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsPosIncRatioPi0FitA[i] ,ptBinsIncRatioErr ,errorsPosErrIncRatioPi0FitA[i] );
            negativeErrorsCorrIncRatioPi0FitA[i] = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsNegCorrIncRatioPi0FitA[i] ,ptBinsIncRatioErr ,errorsNegErrCorrIncRatioPi0FitA[i] );
            meanErrorsCorrIncRatioPi0FitA[i] = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanCorrIncRatioPi0FitA[i] ,ptBinsIncRatioErr ,errorsMeanErrCorrIncRatioPi0FitA[i] );
            positiveErrorsCorrIncRatioPi0FitA[i] = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsPosCorrIncRatioPi0FitA[i] ,ptBinsIncRatioErr ,errorsPosErrCorrIncRatioPi0FitA[i] );
         }
         if(type==2){
            errorsNegIncRatioB[i] = graphNegErrorsIncRatioB[i]->GetY();
            errorsNegErrIncRatioB[i] = graphNegErrorsIncRatioB[i]->GetEYhigh();
            errorsPosIncRatioB[i] = graphPosErrorsIncRatioB[i]->GetY();
            errorsPosErrIncRatioB[i] = graphPosErrorsIncRatioB[i]->GetEYhigh();

            errorsPosIncRatioPi0FitB[i] = graphPosErrorsIncRatioPi0FitB[i]->GetY();
            errorsPosErrIncRatioPi0FitB[i] = graphPosErrorsIncRatioPi0FitB[i]->GetEYhigh();
            errorsNegIncRatioPi0FitB[i] = graphNegErrorsIncRatioPi0FitB[i]->GetY();
            errorsNegErrIncRatioPi0FitB[i] = graphNegErrorsIncRatioPi0FitB[i]->GetEYhigh();

            CalculateMeanSysErr(errorsMeanIncRatioB[i], errorsMeanErrIncRatioB[i], errorsPosIncRatioB[i], errorsNegIncRatioB[i], ConstnBinsIncRatio);
            CalculateMeanSysErr(errorsMeanIncRatioPi0FitB[i], errorsMeanErrIncRatioPi0FitB[i], errorsPosIncRatioPi0FitB[i], errorsNegIncRatioPi0FitB[i], ConstnBinsIncRatio);

            CorrectSystematicErrorsWithMean(errorsPosIncRatioB[i],errorsPosErrIncRatioB[i], errorsPosCorrIncRatioB[i], errorsPosErrCorrIncRatioB[i], ConstnBinsIncRatio);
            CorrectSystematicErrorsWithMean(errorsNegIncRatioB[i],errorsNegErrIncRatioB[i], errorsNegCorrIncRatioB[i], errorsNegErrCorrIncRatioB[i], ConstnBinsIncRatio);
            CorrectSystematicErrorsWithMean(errorsMeanIncRatioB[i], errorsMeanErrIncRatioB[i], errorsMeanCorrIncRatioB[i], errorsMeanErrCorrIncRatioB[i], ConstnBinsIncRatio);

            CorrectSystematicErrorsWithMean(errorsPosIncRatioPi0FitB[i],errorsPosErrIncRatioPi0FitB[i], errorsPosCorrIncRatioPi0FitB[i], errorsPosErrCorrIncRatioPi0FitB[i], ConstnBinsIncRatio);
            CorrectSystematicErrorsWithMean(errorsNegIncRatioPi0FitB[i],errorsNegErrIncRatioPi0FitB[i], errorsNegCorrIncRatioPi0FitB[i], errorsNegErrCorrIncRatioPi0FitB[i], ConstnBinsIncRatio);
            CorrectSystematicErrorsWithMean(errorsMeanIncRatioPi0FitB[i], errorsMeanErrIncRatioPi0FitB[i], errorsMeanCorrIncRatioPi0FitB[i], errorsMeanErrCorrIncRatioPi0FitB[i], ConstnBinsIncRatio);

            negativeErrorsIncRatioB[i] = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsNegIncRatioB[i] ,ptBinsIncRatioErr ,errorsNegErrIncRatioB[i] );
            meanErrorsIncRatioB[i] = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanIncRatioB[i] ,ptBinsIncRatioErr ,errorsMeanErrIncRatioB[i] );
            positiveErrorsIncRatioB[i] = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsPosIncRatioB[i] ,ptBinsIncRatioErr ,errorsPosErrIncRatioB[i] );
            negativeErrorsCorrIncRatioB[i] = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsNegCorrIncRatioB[i] ,ptBinsIncRatioErr ,errorsNegErrCorrIncRatioB[i] );
            meanErrorsCorrIncRatioB[i] = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanCorrIncRatioB[i] ,ptBinsIncRatioErr ,errorsMeanErrCorrIncRatioB[i] );
            positiveErrorsCorrIncRatioB[i] = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsPosCorrIncRatioB[i] ,ptBinsIncRatioErr ,errorsPosErrCorrIncRatioB[i] );

            negativeErrorsIncRatioPi0FitB[i] = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsNegIncRatioPi0FitB[i] ,ptBinsIncRatioErr ,errorsNegErrIncRatioPi0FitB[i] );
            meanErrorsIncRatioPi0FitB[i] = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanIncRatioPi0FitB[i] ,ptBinsIncRatioErr ,errorsMeanErrIncRatioPi0FitB[i] );
            positiveErrorsIncRatioPi0FitB[i] = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsPosIncRatioPi0FitB[i] ,ptBinsIncRatioErr ,errorsPosErrIncRatioPi0FitB[i] );
            negativeErrorsCorrIncRatioPi0FitB[i] = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsNegCorrIncRatioPi0FitB[i] ,ptBinsIncRatioErr ,errorsNegErrCorrIncRatioPi0FitB[i] );
            meanErrorsCorrIncRatioPi0FitB[i] = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanCorrIncRatioPi0FitB[i] ,ptBinsIncRatioErr ,errorsMeanErrCorrIncRatioPi0FitB[i] );
            positiveErrorsCorrIncRatioPi0FitB[i] = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsPosCorrIncRatioPi0FitB[i] ,ptBinsIncRatioErr ,errorsPosErrCorrIncRatioPi0FitB[i] );
         }
         if(type==3){
            errorsNegIncRatioC[i] = graphNegErrorsIncRatioC[i]->GetY();
            errorsNegErrIncRatioC[i] = graphNegErrorsIncRatioC[i]->GetEYhigh();
            errorsPosIncRatioC[i] = graphPosErrorsIncRatioC[i]->GetY();
            errorsPosErrIncRatioC[i] = graphPosErrorsIncRatioC[i]->GetEYhigh();

            errorsPosIncRatioPi0FitC[i] = graphPosErrorsIncRatioPi0FitC[i]->GetY();
            errorsPosErrIncRatioPi0FitC[i] = graphPosErrorsIncRatioPi0FitC[i]->GetEYhigh();
            errorsNegIncRatioPi0FitC[i] = graphNegErrorsIncRatioPi0FitC[i]->GetY();
            errorsNegErrIncRatioPi0FitC[i] = graphNegErrorsIncRatioPi0FitC[i]->GetEYhigh();

            CalculateMeanSysErr(errorsMeanIncRatioC[i], errorsMeanErrIncRatioC[i], errorsPosIncRatioC[i], errorsNegIncRatioC[i], ConstnBinsIncRatio);
            CalculateMeanSysErr(errorsMeanIncRatioPi0FitC[i], errorsMeanErrIncRatioPi0FitC[i], errorsPosIncRatioPi0FitC[i], errorsNegIncRatioPi0FitC[i], ConstnBinsIncRatio);

            CorrectSystematicErrorsWithMean(errorsPosIncRatioC[i],errorsPosErrIncRatioC[i], errorsPosCorrIncRatioC[i], errorsPosErrCorrIncRatioC[i], ConstnBinsIncRatio);
            CorrectSystematicErrorsWithMean(errorsNegIncRatioC[i],errorsNegErrIncRatioC[i], errorsNegCorrIncRatioC[i], errorsNegErrCorrIncRatioC[i], ConstnBinsIncRatio);
            CorrectSystematicErrorsWithMean(errorsMeanIncRatioC[i], errorsMeanErrIncRatioC[i], errorsMeanCorrIncRatioC[i], errorsMeanErrCorrIncRatioC[i], ConstnBinsIncRatio);

            CorrectSystematicErrorsWithMean(errorsPosIncRatioPi0FitC[i],errorsPosErrIncRatioPi0FitC[i], errorsPosCorrIncRatioPi0FitC[i], errorsPosErrCorrIncRatioPi0FitC[i], ConstnBinsIncRatio);
            CorrectSystematicErrorsWithMean(errorsNegIncRatioPi0FitC[i],errorsNegErrIncRatioPi0FitC[i], errorsNegCorrIncRatioPi0FitC[i], errorsNegErrCorrIncRatioPi0FitC[i], ConstnBinsIncRatio);
            CorrectSystematicErrorsWithMean(errorsMeanIncRatioPi0FitC[i], errorsMeanErrIncRatioPi0FitC[i], errorsMeanCorrIncRatioPi0FitC[i], errorsMeanErrCorrIncRatioPi0FitC[i], ConstnBinsIncRatio);

            negativeErrorsIncRatioC[i] = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsNegIncRatioC[i] ,ptBinsIncRatioErr ,errorsNegErrIncRatioC[i] );
            meanErrorsIncRatioC[i] = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanIncRatioC[i] ,ptBinsIncRatioErr ,errorsMeanErrIncRatioC[i] );
            positiveErrorsIncRatioC[i] = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsPosIncRatioC[i] ,ptBinsIncRatioErr ,errorsPosErrIncRatioC[i] );
            negativeErrorsCorrIncRatioC[i] = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsNegCorrIncRatioC[i] ,ptBinsIncRatioErr ,errorsNegErrCorrIncRatioC[i] );
            meanErrorsCorrIncRatioC[i] = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanCorrIncRatioC[i] ,ptBinsIncRatioErr ,errorsMeanErrCorrIncRatioC[i] );
            positiveErrorsCorrIncRatioC[i] = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsPosCorrIncRatioC[i] ,ptBinsIncRatioErr ,errorsPosErrCorrIncRatioC[i] );

            negativeErrorsIncRatioPi0FitC[i] = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsNegIncRatioPi0FitC[i] ,ptBinsIncRatioErr ,errorsNegErrIncRatioPi0FitC[i] );
            meanErrorsIncRatioPi0FitC[i] = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanIncRatioPi0FitC[i] ,ptBinsIncRatioErr ,errorsMeanErrIncRatioPi0FitC[i] );
            positiveErrorsIncRatioPi0FitC[i] = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsPosIncRatioPi0FitC[i] ,ptBinsIncRatioErr ,errorsPosErrIncRatioPi0FitC[i] );
            negativeErrorsCorrIncRatioPi0FitC[i] = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsNegCorrIncRatioPi0FitC[i] ,ptBinsIncRatioErr ,errorsNegErrCorrIncRatioPi0FitC[i] );
            meanErrorsCorrIncRatioPi0FitC[i] = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanCorrIncRatioPi0FitC[i] ,ptBinsIncRatioErr ,errorsMeanErrCorrIncRatioPi0FitC[i] );
            positiveErrorsCorrIncRatioPi0FitC[i] = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsPosCorrIncRatioPi0FitC[i] ,ptBinsIncRatioErr ,errorsPosErrCorrIncRatioPi0FitC[i] );
         }

         for (Int_t l = 0; l < ConstnBinsIncRatio; l++){
            errorsPosSummedIncRatio[l] = errorsPosSummedIncRatio[l]+pow(errorsPosIncRatio[i][l],2);
            errorsNegSummedIncRatio[l] = errorsNegSummedIncRatio[l]+ pow(errorsNegIncRatio[i][l],2);
            errorsMeanSummedIncRatio[l] = errorsMeanSummedIncRatio[l]+ pow(errorsMeanIncRatio[i][l],2);
            errorsPosCorrSummedIncRatio[l] = errorsPosCorrSummedIncRatio[l]+pow(errorsPosCorrIncRatio[i][l],2);
            errorsNegCorrSummedIncRatio[l] = errorsNegCorrSummedIncRatio[l] +pow(errorsNegCorrIncRatio[i][l],2);
            errorsMeanCorrSummedIncRatio[l] =errorsMeanCorrSummedIncRatio[l]+ pow(errorsMeanCorrIncRatio[i][l],2);

            errorsPosSummedIncRatioPi0Fit[l] = errorsPosSummedIncRatioPi0Fit[l]+pow(errorsPosIncRatioPi0Fit[i][l],2);
            errorsNegSummedIncRatioPi0Fit[l] = errorsNegSummedIncRatioPi0Fit[l]+ pow(errorsNegIncRatioPi0Fit[i][l],2);
            errorsMeanSummedIncRatioPi0Fit[l] = errorsMeanSummedIncRatioPi0Fit[l]+ pow(errorsMeanIncRatioPi0Fit[i][l],2);
            errorsPosCorrSummedIncRatioPi0Fit[l] = errorsPosCorrSummedIncRatioPi0Fit[l]+pow(errorsPosCorrIncRatioPi0Fit[i][l],2);
            errorsNegCorrSummedIncRatioPi0Fit[l] = errorsNegCorrSummedIncRatioPi0Fit[l] +pow(errorsNegCorrIncRatioPi0Fit[i][l],2);
            errorsMeanCorrSummedIncRatioPi0Fit[l] =errorsMeanCorrSummedIncRatioPi0Fit[l]+ pow(errorsMeanCorrIncRatioPi0Fit[i][l],2);
            if(type == 1){
               errorsPosSummedIncRatioA[l] = errorsPosSummedIncRatioA[l]+pow(errorsPosIncRatioA[i][l],2);
               errorsNegSummedIncRatioA[l] = errorsNegSummedIncRatioA[l]+ pow(errorsNegIncRatioA[i][l],2);
               errorsMeanSummedIncRatioA[l] = errorsMeanSummedIncRatioA[l]+ pow(errorsMeanIncRatioA[i][l],2);
               errorsPosCorrSummedIncRatioA[l] = errorsPosCorrSummedIncRatioA[l]+pow(errorsPosCorrIncRatioA[i][l],2);
               errorsNegCorrSummedIncRatioA[l] = errorsNegCorrSummedIncRatioA[l] +pow(errorsNegCorrIncRatioA[i][l],2);
               errorsMeanCorrSummedIncRatioA[l] =errorsMeanCorrSummedIncRatioA[l]+ pow(errorsMeanCorrIncRatioA[i][l],2);

               errorsPosSummedIncRatioPi0FitA[l] = errorsPosSummedIncRatioPi0FitA[l]+pow(errorsPosIncRatioPi0FitA[i][l],2);
               errorsNegSummedIncRatioPi0FitA[l] = errorsNegSummedIncRatioPi0FitA[l]+ pow(errorsNegIncRatioPi0FitA[i][l],2);
               errorsMeanSummedIncRatioPi0FitA[l] = errorsMeanSummedIncRatioPi0FitA[l]+ pow(errorsMeanIncRatioPi0FitA[i][l],2);
               errorsPosCorrSummedIncRatioPi0FitA[l] = errorsPosCorrSummedIncRatioPi0FitA[l]+pow(errorsPosCorrIncRatioPi0FitA[i][l],2);
               errorsNegCorrSummedIncRatioPi0FitA[l] = errorsNegCorrSummedIncRatioPi0FitA[l] +pow(errorsNegCorrIncRatioPi0FitA[i][l],2);
               errorsMeanCorrSummedIncRatioPi0FitA[l] =errorsMeanCorrSummedIncRatioPi0FitA[l]+ pow(errorsMeanCorrIncRatioPi0FitA[i][l],2);
            }
            if(type == 2){
               errorsPosSummedIncRatioB[l] = errorsPosSummedIncRatioB[l]+pow(errorsPosIncRatioB[i][l],2);
               errorsNegSummedIncRatioB[l] = errorsNegSummedIncRatioB[l]+ pow(errorsNegIncRatioB[i][l],2);
               errorsMeanSummedIncRatioB[l] = errorsMeanSummedIncRatioB[l]+ pow(errorsMeanIncRatioB[i][l],2);
               errorsPosCorrSummedIncRatioB[l] = errorsPosCorrSummedIncRatioB[l]+pow(errorsPosCorrIncRatioB[i][l],2);
               errorsNegCorrSummedIncRatioB[l] = errorsNegCorrSummedIncRatioB[l] +pow(errorsNegCorrIncRatioB[i][l],2);
               errorsMeanCorrSummedIncRatioB[l] =errorsMeanCorrSummedIncRatioB[l]+ pow(errorsMeanCorrIncRatioB[i][l],2);

               errorsPosSummedIncRatioPi0FitB[l] = errorsPosSummedIncRatioPi0FitB[l]+pow(errorsPosIncRatioPi0FitB[i][l],2);
               errorsNegSummedIncRatioPi0FitB[l] = errorsNegSummedIncRatioPi0FitB[l]+ pow(errorsNegIncRatioPi0FitB[i][l],2);
               errorsMeanSummedIncRatioPi0FitB[l] = errorsMeanSummedIncRatioPi0FitB[l]+ pow(errorsMeanIncRatioPi0FitB[i][l],2);
               errorsPosCorrSummedIncRatioPi0FitB[l] = errorsPosCorrSummedIncRatioPi0FitB[l]+pow(errorsPosCorrIncRatioPi0FitB[i][l],2);
               errorsNegCorrSummedIncRatioPi0FitB[l] = errorsNegCorrSummedIncRatioPi0FitB[l] +pow(errorsNegCorrIncRatioPi0FitB[i][l],2);
               errorsMeanCorrSummedIncRatioPi0FitB[l] =errorsMeanCorrSummedIncRatioPi0FitB[l]+ pow(errorsMeanCorrIncRatioPi0FitB[i][l],2);
            }
            if(type == 3){
               errorsPosSummedIncRatioC[l] = errorsPosSummedIncRatioC[l]+pow(errorsPosIncRatioC[i][l],2);
               errorsNegSummedIncRatioC[l] = errorsNegSummedIncRatioC[l]+ pow(errorsNegIncRatioC[i][l],2);
               errorsMeanSummedIncRatioC[l] = errorsMeanSummedIncRatioC[l]+ pow(errorsMeanIncRatioC[i][l],2);
               errorsPosCorrSummedIncRatioC[l] = errorsPosCorrSummedIncRatioC[l]+pow(errorsPosCorrIncRatioC[i][l],2);
               errorsNegCorrSummedIncRatioC[l] = errorsNegCorrSummedIncRatioC[l] +pow(errorsNegCorrIncRatioC[i][l],2);
               errorsMeanCorrSummedIncRatioC[l] =errorsMeanCorrSummedIncRatioC[l]+ pow(errorsMeanCorrIncRatioC[i][l],2);

               errorsPosSummedIncRatioPi0FitC[l] = errorsPosSummedIncRatioPi0FitC[l]+pow(errorsPosIncRatioPi0FitC[i][l],2);
               errorsNegSummedIncRatioPi0FitC[l] = errorsNegSummedIncRatioPi0FitC[l]+ pow(errorsNegIncRatioPi0FitC[i][l],2);
               errorsMeanSummedIncRatioPi0FitC[l] = errorsMeanSummedIncRatioPi0FitC[l]+ pow(errorsMeanIncRatioPi0FitC[i][l],2);
               errorsPosCorrSummedIncRatioPi0FitC[l] = errorsPosCorrSummedIncRatioPi0FitC[l]+pow(errorsPosCorrIncRatioPi0FitC[i][l],2);
               errorsNegCorrSummedIncRatioPi0FitC[l] = errorsNegCorrSummedIncRatioPi0FitC[l] +pow(errorsNegCorrIncRatioPi0FitC[i][l],2);
               errorsMeanCorrSummedIncRatioPi0FitC[l] =errorsMeanCorrSummedIncRatioPi0FitC[l]+ pow(errorsMeanCorrIncRatioPi0FitC[i][l],2);
            }
         }
      }

      if(i<ErrorPi0){

         graphPosErrorsPi0[i] = (TGraphAsymmErrors*)fileSystematics->Get(Form("Pi0_SystErrorRelPos_%s_%s",nameCutVariationsDoubleRatio[i].Data(),(GetCentralityStringA(cutSel)).Data()));
         graphNegErrorsPi0[i] = (TGraphAsymmErrors*)fileSystematics->Get(Form("Pi0_SystErrorRelNeg_%s_%s",nameCutVariationsDoubleRatio[i].Data(),(GetCentralityStringA(cutSel)).Data()));

         graphPosErrorsPi0Fit[i] = (TGraphAsymmErrors*)fileSystematics->Get(Form("Pi0Fit_SystErrorRelPos_%s_%s",nameCutVariationsDoubleRatio[i].Data(),(GetCentralityStringA(cutSel)).Data()));
         graphNegErrorsPi0Fit[i] = (TGraphAsymmErrors*)fileSystematics->Get(Form("Pi0Fit_SystErrorRelNeg_%s_%s",nameCutVariationsDoubleRatio[i].Data(),(GetCentralityStringA(cutSel)).Data()));

         graphPosErrorsPi0[i]->RemovePoint(0);
         graphNegErrorsPi0[i]->RemovePoint(0);
         graphPosErrorsPi0[i]->RemovePoint(0);
         graphNegErrorsPi0[i]->RemovePoint(0);

         graphPosErrorsPi0Fit[i]->RemovePoint(0);
         graphNegErrorsPi0Fit[i]->RemovePoint(0);
         graphPosErrorsPi0Fit[i]->RemovePoint(0);
         graphNegErrorsPi0Fit[i]->RemovePoint(0);

         if(i==0){
            ptBinsPi0 = graphPosErrorsPi0[i]->GetX();
            nBinsGraphDoubeRatio = graphPosErrorsPi0[i]->GetN();
            ptBinsPi0Err = graphPosErrorsPi0[i]->GetEXhigh();
         }



         errorsNegPi0[i] = graphNegErrorsPi0[i]->GetY();
         errorsNegErrPi0[i] = graphNegErrorsPi0[i]->GetEYhigh();
         errorsPosPi0[i] = graphPosErrorsPi0[i]->GetY();
         errorsPosErrPi0[i] = graphPosErrorsPi0[i]->GetEYhigh();


         errorsPosPi0Fit[i] = graphPosErrorsPi0Fit[i]->GetY();
         errorsPosErrPi0Fit[i] = graphPosErrorsPi0Fit[i]->GetEYhigh();
         errorsNegPi0Fit[i] = graphNegErrorsPi0Fit[i]->GetY();
         errorsNegErrPi0Fit[i] = graphNegErrorsPi0Fit[i]->GetEYhigh();


         CalculateMeanSysErr(errorsMeanPi0[i], errorsMeanErrPi0[i], errorsPosPi0[i], errorsNegPi0[i], ConstnBinsPi0);
         CalculateMeanSysErr(errorsMeanPi0Fit[i], errorsMeanErrPi0Fit[i], errorsPosPi0Fit[i], errorsNegPi0Fit[i], ConstnBinsPi0);

         CorrectSystematicErrorsWithMean(errorsPosPi0[i],errorsPosErrPi0[i], errorsPosCorrPi0[i], errorsPosErrCorrPi0[i], ConstnBinsPi0);
         CorrectSystematicErrorsWithMean(errorsNegPi0[i],errorsNegErrPi0[i], errorsNegCorrPi0[i], errorsNegErrCorrPi0[i], ConstnBinsPi0);
         CorrectSystematicErrorsWithMean(errorsMeanPi0[i], errorsMeanErrPi0[i], errorsMeanCorrPi0[i], errorsMeanErrCorrPi0[i], ConstnBinsPi0);

         CorrectSystematicErrorsWithMean(errorsPosPi0Fit[i],errorsPosErrPi0Fit[i], errorsPosCorrPi0Fit[i], errorsPosErrCorrPi0Fit[i], ConstnBinsPi0);
         CorrectSystematicErrorsWithMean(errorsNegPi0Fit[i],errorsNegErrPi0Fit[i], errorsNegCorrPi0Fit[i], errorsNegErrCorrPi0Fit[i], ConstnBinsPi0);
         CorrectSystematicErrorsWithMean(errorsMeanPi0Fit[i], errorsMeanErrPi0Fit[i], errorsMeanCorrPi0Fit[i], errorsMeanErrCorrPi0Fit[i], ConstnBinsPi0);

         negativeErrorsPi0[i] = new TGraphErrors(ConstnBinsPi0,ptBinsPi0 ,errorsNegPi0[i] ,ptBinsPi0Err ,errorsNegErrPi0[i] );
         meanErrorsPi0[i] = new TGraphErrors(ConstnBinsPi0,ptBinsPi0 ,errorsMeanPi0[i] ,ptBinsPi0Err ,errorsMeanErrPi0[i] );
         positiveErrorsPi0[i] = new TGraphErrors(ConstnBinsPi0,ptBinsPi0 ,errorsPosPi0[i] ,ptBinsPi0Err ,errorsPosErrPi0[i] );
         negativeErrorsCorrPi0[i] = new TGraphErrors(ConstnBinsPi0,ptBinsPi0 ,errorsNegCorrPi0[i] ,ptBinsPi0Err ,errorsNegErrCorrPi0[i] );
         meanErrorsCorrPi0[i] = new TGraphErrors(ConstnBinsPi0,ptBinsPi0 ,errorsMeanCorrPi0[i] ,ptBinsPi0Err ,errorsMeanErrCorrPi0[i] );
         positiveErrorsCorrPi0[i] = new TGraphErrors(ConstnBinsPi0,ptBinsPi0 ,errorsPosCorrPi0[i] ,ptBinsPi0Err ,errorsPosErrCorrPi0[i] );

         negativeErrorsPi0Fit[i] = new TGraphErrors(ConstnBinsPi0,ptBinsPi0 ,errorsNegPi0Fit[i] ,ptBinsPi0Err ,errorsNegErrPi0Fit[i] );
         meanErrorsPi0Fit[i] = new TGraphErrors(ConstnBinsPi0,ptBinsPi0 ,errorsMeanPi0Fit[i] ,ptBinsPi0Err ,errorsMeanErrPi0Fit[i] );
         positiveErrorsPi0Fit[i] = new TGraphErrors(ConstnBinsPi0,ptBinsPi0 ,errorsPosPi0Fit[i] ,ptBinsPi0Err ,errorsPosErrPi0Fit[i] );
         negativeErrorsCorrPi0Fit[i] = new TGraphErrors(ConstnBinsPi0,ptBinsPi0 ,errorsNegCorrPi0Fit[i] ,ptBinsPi0Err ,errorsNegErrCorrPi0Fit[i] );
         meanErrorsCorrPi0Fit[i] = new TGraphErrors(ConstnBinsPi0,ptBinsPi0 ,errorsMeanCorrPi0Fit[i] ,ptBinsPi0Err ,errorsMeanErrCorrPi0Fit[i] );
         positiveErrorsCorrPi0Fit[i] = new TGraphErrors(ConstnBinsPi0,ptBinsPi0 ,errorsPosCorrPi0Fit[i] ,ptBinsPi0Err ,errorsPosErrCorrPi0Fit[i] );

         for (Int_t l = 0; l < ConstnBinsPi0; l++){

            errorsPosSummedPi0[l] = errorsPosSummedPi0[l]+pow(errorsPosPi0[i][l],2);
            errorsNegSummedPi0[l] = errorsNegSummedPi0[l]+ pow(errorsNegPi0[i][l],2);
            errorsMeanSummedPi0[l] = errorsMeanSummedPi0[l]+ pow(errorsMeanPi0[i][l],2);
            errorsPosCorrSummedPi0[l] = errorsPosCorrSummedPi0[l]+pow(errorsPosCorrPi0[i][l],2);
            errorsNegCorrSummedPi0[l] = errorsNegCorrSummedPi0[l] +pow(errorsNegCorrPi0[i][l],2);
            errorsMeanCorrSummedPi0[l] =errorsMeanCorrSummedPi0[l]+ pow(errorsMeanCorrPi0[i][l],2);

            errorsPosSummedPi0Fit[l] = errorsPosSummedPi0Fit[l]+pow(errorsPosPi0Fit[i][l],2);
            errorsNegSummedPi0Fit[l] = errorsNegSummedPi0Fit[l]+ pow(errorsNegPi0Fit[i][l],2);
            errorsMeanSummedPi0Fit[l] = errorsMeanSummedPi0Fit[l]+ pow(errorsMeanPi0Fit[i][l],2);
            errorsPosCorrSummedPi0Fit[l] = errorsPosCorrSummedPi0Fit[l]+pow(errorsPosCorrPi0Fit[i][l],2);
            errorsNegCorrSummedPi0Fit[l] = errorsNegCorrSummedPi0Fit[l] +pow(errorsNegCorrPi0Fit[i][l],2);
            errorsMeanCorrSummedPi0Fit[l] =errorsMeanCorrSummedPi0Fit[l]+ pow(errorsMeanCorrPi0Fit[i][l],2);
         }
      }


      if(i<ErrorGammaSpec){
         graphPosErrorsGamma[i] = (TGraphAsymmErrors*)fileSystematics->Get(Form("Gamma_SystErrorRelPos_%s_%s",nameCutVariationsDoubleRatio[i].Data(),(GetCentralityStringA(cutSel)).Data()));
         graphNegErrorsGamma[i] = (TGraphAsymmErrors*)fileSystematics->Get(Form("Gamma_SystErrorRelNeg_%s_%s",nameCutVariationsDoubleRatio[i].Data(),(GetCentralityStringA(cutSel)).Data()));

         graphPosErrorsGamma[i]->RemovePoint(0);
         graphNegErrorsGamma[i]->RemovePoint(0);
         graphPosErrorsGamma[i]->RemovePoint(0);
         graphNegErrorsGamma[i]->RemovePoint(0);

         if(type == 1){ // A Errors
            graphPosErrorsGammaA[i] = (TGraphAsymmErrors*)fileSystematics->Get(Form("Gamma_SystErrorRelPos_%s_%s",nameCutVariationsGamma[i].Data(),(GetCentralityStringA(cutSel)).Data()));
            graphNegErrorsGammaA[i] = (TGraphAsymmErrors*)fileSystematics->Get(Form("Gamma_SystErrorRelNeg_%s_%s",nameCutVariationsGamma[i].Data(),(GetCentralityStringA(cutSel)).Data()));

            graphPosErrorsGammaA[i]->RemovePoint(0);
            graphNegErrorsGammaA[i]->RemovePoint(0);
            graphPosErrorsGammaA[i]->RemovePoint(0);
            graphNegErrorsGammaA[i]->RemovePoint(0);
         }
         if(type == 2){ // B Errors
            graphPosErrorsGammaB[i] = (TGraphAsymmErrors*)fileSystematics->Get(Form("Gamma_SystErrorRelPos_%s_%s",nameCutVariationsGamma[i].Data(),(GetCentralityStringA(cutSel)).Data()));
            graphNegErrorsGammaB[i] = (TGraphAsymmErrors*)fileSystematics->Get(Form("Gamma_SystErrorRelNeg_%s_%s",nameCutVariationsGamma[i].Data(),(GetCentralityStringA(cutSel)).Data()));

            graphPosErrorsGammaB[i]->RemovePoint(0);
            graphNegErrorsGammaB[i]->RemovePoint(0);
            graphPosErrorsGammaB[i]->RemovePoint(0);
            graphNegErrorsGammaB[i]->RemovePoint(0);
         }
         if(type == 3){ // C Errors
            graphPosErrorsGammaC[i] = (TGraphAsymmErrors*)fileSystematics->Get(Form("Gamma_SystErrorRelPos_%s_%s",nameCutVariationsGamma[i].Data(),(GetCentralityStringB(cutSel)).Data()));
            graphNegErrorsGammaC[i] = (TGraphAsymmErrors*)fileSystematics->Get(Form("Gamma_SystErrorRelNeg_%s_%s",nameCutVariationsGamma[i].Data(),(GetCentralityStringB(cutSel)).Data()));

            graphPosErrorsGammaC[i]->RemovePoint(0);
            graphNegErrorsGammaC[i]->RemovePoint(0);
            graphPosErrorsGammaC[i]->RemovePoint(0);
            graphNegErrorsGammaC[i]->RemovePoint(0);
         }
         if(i==0){
            ptBinsGamma = graphPosErrorsGamma[i]->GetX();
            nBinsGraphGamma = graphPosErrorsGamma[i]->GetN();
            ptBinsGammaErr = graphPosErrorsGamma[i]->GetEXhigh();
         }

         errorsPosGamma[i] = graphPosErrorsGamma[i]->GetY();
         errorsPosErrGamma[i] = graphPosErrorsGamma[i]->GetEYhigh();
         errorsNegGamma[i] = graphNegErrorsGamma[i]->GetY();
         errorsNegErrGamma[i] = graphNegErrorsGamma[i]->GetEYhigh();

         CalculateMeanSysErr(errorsMeanGamma[i], errorsMeanErrGamma[i], errorsPosGamma[i], errorsNegGamma[i], ConstnBinsGamma);
         CorrectSystematicErrorsWithMean(errorsPosGamma[i],errorsPosErrGamma[i], errorsPosCorrGamma[i], errorsPosErrCorrGamma[i], ConstnBinsGamma);
         CorrectSystematicErrorsWithMean(errorsNegGamma[i],errorsNegErrGamma[i], errorsNegCorrGamma[i], errorsNegErrCorrGamma[i], ConstnBinsGamma);
         CorrectSystematicErrorsWithMean(errorsMeanGamma[i], errorsMeanErrGamma[i], errorsMeanCorrGamma[i], errorsMeanErrCorrGamma[i], ConstnBinsGamma);

         negativeErrorsGamma[i] = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsNegGamma[i] ,ptBinsGammaErr ,errorsNegErrGamma[i] );
         meanErrorsGamma[i] = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsMeanGamma[i] ,ptBinsGammaErr ,errorsMeanErrGamma[i] );
         positiveErrorsGamma[i] = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsPosGamma[i] ,ptBinsGammaErr ,errorsPosErrGamma[i] );
         negativeErrorsCorrGamma[i] = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsNegCorrGamma[i] ,ptBinsGammaErr ,errorsNegErrCorrGamma[i] );
         meanErrorsCorrGamma[i] = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsMeanCorrGamma[i] ,ptBinsGammaErr ,errorsMeanErrCorrGamma[i] );
         positiveErrorsCorrGamma[i] = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsPosCorrGamma[i] ,ptBinsGammaErr ,errorsPosErrCorrGamma[i] );


         if(type==1){
            errorsNegGammaA[i] = graphNegErrorsGammaA[i]->GetY();
            errorsNegErrGammaA[i] = graphNegErrorsGammaA[i]->GetEYhigh();
            errorsPosGammaA[i] = graphPosErrorsGammaA[i]->GetY();
            errorsPosErrGammaA[i] = graphPosErrorsGammaA[i]->GetEYhigh();

            CalculateMeanSysErr(errorsMeanGammaA[i], errorsMeanErrGammaA[i], errorsPosGammaA[i], errorsNegGammaA[i], ConstnBinsGamma);

            CorrectSystematicErrorsWithMean(errorsPosGammaA[i],errorsPosErrGammaA[i], errorsPosCorrGammaA[i], errorsPosErrCorrGammaA[i], ConstnBinsGamma);
            CorrectSystematicErrorsWithMean(errorsNegGammaA[i],errorsNegErrGammaA[i], errorsNegCorrGammaA[i], errorsNegErrCorrGammaA[i], ConstnBinsGamma);
            CorrectSystematicErrorsWithMean(errorsMeanGammaA[i], errorsMeanErrGammaA[i], errorsMeanCorrGammaA[i], errorsMeanErrCorrGammaA[i], ConstnBinsGamma);

            negativeErrorsGammaA[i] = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsNegGammaA[i] ,ptBinsGammaErr ,errorsNegErrGammaA[i] );
            meanErrorsGammaA[i] = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsMeanGammaA[i] ,ptBinsGammaErr ,errorsMeanErrGammaA[i] );
            positiveErrorsGammaA[i] = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsPosGammaA[i] ,ptBinsGammaErr ,errorsPosErrGammaA[i] );
            negativeErrorsCorrGammaA[i] = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsNegCorrGammaA[i] ,ptBinsGammaErr ,errorsNegErrCorrGammaA[i] );
            meanErrorsCorrGammaA[i] = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsMeanCorrGammaA[i] ,ptBinsGammaErr ,errorsMeanErrCorrGammaA[i] );
            positiveErrorsCorrGammaA[i] = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsPosCorrGammaA[i] ,ptBinsGammaErr ,errorsPosErrCorrGammaA[i] );
         }
         if(type==2){
            errorsNegGammaB[i] = graphNegErrorsGammaB[i]->GetY();
            errorsNegErrGammaB[i] = graphNegErrorsGammaB[i]->GetEYhigh();
            errorsPosGammaB[i] = graphPosErrorsGammaB[i]->GetY();
            errorsPosErrGammaB[i] = graphPosErrorsGammaB[i]->GetEYhigh();

            CalculateMeanSysErr(errorsMeanGammaB[i], errorsMeanErrGammaB[i], errorsPosGammaB[i], errorsNegGammaB[i], ConstnBinsGamma);

            CorrectSystematicErrorsWithMean(errorsPosGammaB[i],errorsPosErrGammaB[i], errorsPosCorrGammaB[i], errorsPosErrCorrGammaB[i], ConstnBinsGamma);
            CorrectSystematicErrorsWithMean(errorsNegGammaB[i],errorsNegErrGammaB[i], errorsNegCorrGammaB[i], errorsNegErrCorrGammaB[i], ConstnBinsGamma);
            CorrectSystematicErrorsWithMean(errorsMeanGammaB[i], errorsMeanErrGammaB[i], errorsMeanCorrGammaB[i], errorsMeanErrCorrGammaB[i], ConstnBinsGamma);

            negativeErrorsGammaB[i] = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsNegGammaB[i] ,ptBinsGammaErr ,errorsNegErrGammaB[i] );
            meanErrorsGammaB[i] = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsMeanGammaB[i] ,ptBinsGammaErr ,errorsMeanErrGammaB[i] );
            positiveErrorsGammaB[i] = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsPosGammaB[i] ,ptBinsGammaErr ,errorsPosErrGammaB[i] );
            negativeErrorsCorrGammaB[i] = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsNegCorrGammaB[i] ,ptBinsGammaErr ,errorsNegErrCorrGammaB[i] );
            meanErrorsCorrGammaB[i] = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsMeanCorrGammaB[i] ,ptBinsGammaErr ,errorsMeanErrCorrGammaB[i] );
            positiveErrorsCorrGammaB[i] = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsPosCorrGammaB[i] ,ptBinsGammaErr ,errorsPosErrCorrGammaB[i] );
         }
         if(type==3){
            errorsNegGammaC[i] = graphNegErrorsGammaC[i]->GetY();
            errorsNegErrGammaC[i] = graphNegErrorsGammaC[i]->GetEYhigh();
            errorsPosGammaC[i] = graphPosErrorsGammaC[i]->GetY();
            errorsPosErrGammaC[i] = graphPosErrorsGammaC[i]->GetEYhigh();

            CalculateMeanSysErr(errorsMeanGammaC[i], errorsMeanErrGammaC[i], errorsPosGammaC[i], errorsNegGammaC[i], ConstnBinsGamma);

            CorrectSystematicErrorsWithMean(errorsPosGammaC[i],errorsPosErrGammaC[i], errorsPosCorrGammaC[i], errorsPosErrCorrGammaC[i], ConstnBinsGamma);
            CorrectSystematicErrorsWithMean(errorsNegGammaC[i],errorsNegErrGammaC[i], errorsNegCorrGammaC[i], errorsNegErrCorrGammaC[i], ConstnBinsGamma);
            CorrectSystematicErrorsWithMean(errorsMeanGammaC[i], errorsMeanErrGammaC[i], errorsMeanCorrGammaC[i], errorsMeanErrCorrGammaC[i], ConstnBinsGamma);

            negativeErrorsGammaC[i] = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsNegGammaC[i] ,ptBinsGammaErr ,errorsNegErrGammaC[i] );
            meanErrorsGammaC[i] = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsMeanGammaC[i] ,ptBinsGammaErr ,errorsMeanErrGammaC[i] );
            positiveErrorsGammaC[i] = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsPosGammaC[i] ,ptBinsGammaErr ,errorsPosErrGammaC[i] );
            negativeErrorsCorrGammaC[i] = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsNegCorrGammaC[i] ,ptBinsGammaErr ,errorsNegErrCorrGammaC[i] );
            meanErrorsCorrGammaC[i] = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsMeanCorrGammaC[i] ,ptBinsGammaErr ,errorsMeanErrCorrGammaC[i] );
            positiveErrorsCorrGammaC[i] = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsPosCorrGammaC[i] ,ptBinsGammaErr ,errorsPosErrCorrGammaC[i] );
         }

         for (Int_t l = 0; l < ConstnBinsGamma; l++){
            errorsPosSummedGamma[l] = errorsPosSummedGamma[l]+pow(errorsPosGamma[i][l],2);
            errorsNegSummedGamma[l] = errorsNegSummedGamma[l]+ pow(errorsNegGamma[i][l],2);
            errorsMeanSummedGamma[l] = errorsMeanSummedGamma[l]+ pow(errorsMeanGamma[i][l],2);
            errorsPosCorrSummedGamma[l] = errorsPosCorrSummedGamma[l]+pow(errorsPosCorrGamma[i][l],2);
            errorsNegCorrSummedGamma[l] = errorsNegCorrSummedGamma[l] +pow(errorsNegCorrGamma[i][l],2);
            errorsMeanCorrSummedGamma[l] =errorsMeanCorrSummedGamma[l]+ pow(errorsMeanCorrGamma[i][l],2);
            if(type == 1){
               errorsPosSummedGammaA[l] = errorsPosSummedGammaA[l]+pow(errorsPosGammaA[i][l],2);
               errorsNegSummedGammaA[l] = errorsNegSummedGammaA[l]+ pow(errorsNegGammaA[i][l],2);
               errorsMeanSummedGammaA[l] = errorsMeanSummedGammaA[l]+ pow(errorsMeanGammaA[i][l],2);
               errorsPosCorrSummedGammaA[l] = errorsPosCorrSummedGammaA[l]+pow(errorsPosCorrGammaA[i][l],2);
               errorsNegCorrSummedGammaA[l] = errorsNegCorrSummedGammaA[l] +pow(errorsNegCorrGammaA[i][l],2);
               errorsMeanCorrSummedGammaA[l] =errorsMeanCorrSummedGammaA[l]+ pow(errorsMeanCorrGammaA[i][l],2);
            }
            if(type == 2){
               errorsPosSummedGammaB[l] = errorsPosSummedGammaB[l]+pow(errorsPosGammaB[i][l],2);
               errorsNegSummedGammaB[l] = errorsNegSummedGammaB[l]+ pow(errorsNegGammaB[i][l],2);
               errorsMeanSummedGammaB[l] = errorsMeanSummedGammaB[l]+ pow(errorsMeanGammaB[i][l],2);
               errorsPosCorrSummedGammaB[l] = errorsPosCorrSummedGammaB[l]+pow(errorsPosCorrGammaB[i][l],2);
               errorsNegCorrSummedGammaB[l] = errorsNegCorrSummedGammaB[l] +pow(errorsNegCorrGammaB[i][l],2);
               errorsMeanCorrSummedGammaB[l] =errorsMeanCorrSummedGammaB[l]+ pow(errorsMeanCorrGammaB[i][l],2);
            }
            if(type == 3){
               errorsPosSummedGammaC[l] = errorsPosSummedGammaC[l]+pow(errorsPosGammaC[i][l],2);
               errorsNegSummedGammaC[l] = errorsNegSummedGammaC[l]+ pow(errorsNegGammaC[i][l],2);
               errorsMeanSummedGammaC[l] = errorsMeanSummedGammaC[l]+ pow(errorsMeanGammaC[i][l],2);
               errorsPosCorrSummedGammaC[l] = errorsPosCorrSummedGammaC[l]+pow(errorsPosCorrGammaC[i][l],2);
               errorsNegCorrSummedGammaC[l] = errorsNegCorrSummedGammaC[l] +pow(errorsNegCorrGammaC[i][l],2);
               errorsMeanCorrSummedGammaC[l] =errorsMeanCorrSummedGammaC[l]+ pow(errorsMeanCorrGammaC[i][l],2);
            }
         }
      }

      for (Int_t l = 0; l < ConstnBinsDoubleRatio; l++){

         errorsPosSummedDoubleRatio[l] = errorsPosSummedDoubleRatio[l]+pow(errorsPosDoubleRatio[i][l],2);
         errorsNegSummedDoubleRatio[l] = errorsNegSummedDoubleRatio[l]+ pow(errorsNegDoubleRatio[i][l],2);
         errorsMeanSummedDoubleRatio[l] = errorsMeanSummedDoubleRatio[l]+ pow(errorsMeanDoubleRatio[i][l],2);
         errorsPosCorrSummedDoubleRatio[l] = errorsPosCorrSummedDoubleRatio[l]+pow(errorsPosCorrDoubleRatio[i][l],2);
         errorsNegCorrSummedDoubleRatio[l] = errorsNegCorrSummedDoubleRatio[l] +pow(errorsNegCorrDoubleRatio[i][l],2);
         errorsMeanCorrSummedDoubleRatio[l] =errorsMeanCorrSummedDoubleRatio[l]+ pow(errorsMeanCorrDoubleRatio[i][l],2);

         errorsPosSummedDoubleRatioPi0Fit[l] = errorsPosSummedDoubleRatioPi0Fit[l]+pow(errorsPosDoubleRatioPi0Fit[i][l],2);
         errorsNegSummedDoubleRatioPi0Fit[l] = errorsNegSummedDoubleRatioPi0Fit[l]+ pow(errorsNegDoubleRatioPi0Fit[i][l],2);
         errorsMeanSummedDoubleRatioPi0Fit[l] = errorsMeanSummedDoubleRatioPi0Fit[l]+ pow(errorsMeanDoubleRatioPi0Fit[i][l],2);
         errorsPosCorrSummedDoubleRatioPi0Fit[l] = errorsPosCorrSummedDoubleRatioPi0Fit[l]+pow(errorsPosCorrDoubleRatioPi0Fit[i][l],2);
         errorsNegCorrSummedDoubleRatioPi0Fit[l] = errorsNegCorrSummedDoubleRatioPi0Fit[l] +pow(errorsNegCorrDoubleRatioPi0Fit[i][l],2);
         errorsMeanCorrSummedDoubleRatioPi0Fit[l] = errorsMeanCorrSummedDoubleRatioPi0Fit[l]+ pow(errorsMeanCorrDoubleRatioPi0Fit[i][l],2);

         if(type == 1){
            errorsPosSummedDoubleRatioA[l] = errorsPosSummedDoubleRatioA[l]+pow(errorsPosDoubleRatioA[i][l],2);
            errorsNegSummedDoubleRatioA[l] = errorsNegSummedDoubleRatioA[l]+ pow(errorsNegDoubleRatioA[i][l],2);
            errorsMeanSummedDoubleRatioA[l] = errorsMeanSummedDoubleRatioA[l]+ pow(errorsMeanDoubleRatioA[i][l],2);
            errorsPosCorrSummedDoubleRatioA[l] = errorsPosCorrSummedDoubleRatioA[l]+pow(errorsPosCorrDoubleRatioA[i][l],2);
            errorsNegCorrSummedDoubleRatioA[l] = errorsNegCorrSummedDoubleRatioA[l] +pow(errorsNegCorrDoubleRatioA[i][l],2);
            errorsMeanCorrSummedDoubleRatioA[l] =errorsMeanCorrSummedDoubleRatioA[l]+ pow(errorsMeanCorrDoubleRatioA[i][l],2);

            errorsPosSummedDoubleRatioPi0FitA[l] = errorsPosSummedDoubleRatioPi0FitA[l]+pow(errorsPosDoubleRatioPi0FitA[i][l],2);
            errorsNegSummedDoubleRatioPi0FitA[l] = errorsNegSummedDoubleRatioPi0FitA[l]+ pow(errorsNegDoubleRatioPi0FitA[i][l],2);
            errorsMeanSummedDoubleRatioPi0FitA[l] = errorsMeanSummedDoubleRatioPi0FitA[l]+ pow(errorsMeanDoubleRatioPi0FitA[i][l],2);
            errorsPosCorrSummedDoubleRatioPi0FitA[l] = errorsPosCorrSummedDoubleRatioPi0FitA[l]+pow(errorsPosCorrDoubleRatioPi0FitA[i][l],2);
            errorsNegCorrSummedDoubleRatioPi0FitA[l] = errorsNegCorrSummedDoubleRatioPi0FitA[l] +pow(errorsNegCorrDoubleRatioPi0FitA[i][l],2);
            errorsMeanCorrSummedDoubleRatioPi0FitA[l] = errorsMeanCorrSummedDoubleRatioPi0FitA[l]+ pow(errorsMeanCorrDoubleRatioPi0FitA[i][l],2);
         }
         if(type == 2){
            errorsPosSummedDoubleRatioB[l] = errorsPosSummedDoubleRatioB[l]+pow(errorsPosDoubleRatioB[i][l],2);
            errorsNegSummedDoubleRatioB[l] = errorsNegSummedDoubleRatioB[l]+ pow(errorsNegDoubleRatioB[i][l],2);
            errorsMeanSummedDoubleRatioB[l] = errorsMeanSummedDoubleRatioB[l]+ pow(errorsMeanDoubleRatioB[i][l],2);
            errorsPosCorrSummedDoubleRatioB[l] = errorsPosCorrSummedDoubleRatioB[l]+pow(errorsPosCorrDoubleRatioB[i][l],2);
            errorsNegCorrSummedDoubleRatioB[l] = errorsNegCorrSummedDoubleRatioB[l] +pow(errorsNegCorrDoubleRatioB[i][l],2);
            errorsMeanCorrSummedDoubleRatioB[l] =errorsMeanCorrSummedDoubleRatioB[l]+ pow(errorsMeanCorrDoubleRatioB[i][l],2);

            errorsPosSummedDoubleRatioPi0FitB[l] = errorsPosSummedDoubleRatioPi0FitB[l]+pow(errorsPosDoubleRatioPi0FitB[i][l],2);
            errorsNegSummedDoubleRatioPi0FitB[l] = errorsNegSummedDoubleRatioPi0FitB[l]+ pow(errorsNegDoubleRatioPi0FitB[i][l],2);
            errorsMeanSummedDoubleRatioPi0FitB[l] = errorsMeanSummedDoubleRatioPi0FitB[l]+ pow(errorsMeanDoubleRatioPi0FitB[i][l],2);
            errorsPosCorrSummedDoubleRatioPi0FitB[l] = errorsPosCorrSummedDoubleRatioPi0FitB[l]+pow(errorsPosCorrDoubleRatioPi0FitB[i][l],2);
            errorsNegCorrSummedDoubleRatioPi0FitB[l] = errorsNegCorrSummedDoubleRatioPi0FitB[l] +pow(errorsNegCorrDoubleRatioPi0FitB[i][l],2);
            errorsMeanCorrSummedDoubleRatioPi0FitB[l] = errorsMeanCorrSummedDoubleRatioPi0FitB[l]+ pow(errorsMeanCorrDoubleRatioPi0FitB[i][l],2);
         }
         if(type == 3){
            errorsPosSummedDoubleRatioC[l] = errorsPosSummedDoubleRatioC[l]+pow(errorsPosDoubleRatioC[i][l],2);
            errorsNegSummedDoubleRatioC[l] = errorsNegSummedDoubleRatioC[l]+ pow(errorsNegDoubleRatioC[i][l],2);
            errorsMeanSummedDoubleRatioC[l] = errorsMeanSummedDoubleRatioC[l]+ pow(errorsMeanDoubleRatioC[i][l],2);
            errorsPosCorrSummedDoubleRatioC[l] = errorsPosCorrSummedDoubleRatioC[l]+pow(errorsPosCorrDoubleRatioC[i][l],2);
            errorsNegCorrSummedDoubleRatioC[l] = errorsNegCorrSummedDoubleRatioC[l] +pow(errorsNegCorrDoubleRatioC[i][l],2);
            errorsMeanCorrSummedDoubleRatioC[l] =errorsMeanCorrSummedDoubleRatioC[l]+ pow(errorsMeanCorrDoubleRatioC[i][l],2);

            errorsPosSummedDoubleRatioPi0FitC[l] = errorsPosSummedDoubleRatioPi0FitC[l]+pow(errorsPosDoubleRatioPi0FitC[i][l],2);
            errorsNegSummedDoubleRatioPi0FitC[l] = errorsNegSummedDoubleRatioPi0FitC[l]+ pow(errorsNegDoubleRatioPi0FitC[i][l],2);
            errorsMeanSummedDoubleRatioPi0FitC[l] = errorsMeanSummedDoubleRatioPi0FitC[l]+ pow(errorsMeanDoubleRatioPi0FitC[i][l],2);
            errorsPosCorrSummedDoubleRatioPi0FitC[l] = errorsPosCorrSummedDoubleRatioPi0FitC[l]+pow(errorsPosCorrDoubleRatioPi0FitC[i][l],2);
            errorsNegCorrSummedDoubleRatioPi0FitC[l] = errorsNegCorrSummedDoubleRatioPi0FitC[l] +pow(errorsNegCorrDoubleRatioPi0FitC[i][l],2);
            errorsMeanCorrSummedDoubleRatioPi0FitC[l] = errorsMeanCorrSummedDoubleRatioPi0FitC[l]+ pow(errorsMeanCorrDoubleRatioPi0FitC[i][l],2);
         }
      }
   }


   for(Int_t l = 0; l <ConstnBinsDoubleRatio; l++){


      errorsPosSummedDoubleRatio[l] = pow(errorsPosSummedDoubleRatio[l],0.5);
      errorsMeanSummedDoubleRatio[l] = pow(errorsMeanSummedDoubleRatio[l],0.5);
      errorsPosErrSummedDoubleRatio[l] = errorsPosSummedDoubleRatio[l]*0.001;
      errorsMeanErrSummedDoubleRatio[l] = errorsMeanSummedDoubleRatio[l]*0.001;
      errorsNegSummedDoubleRatio[l] = -pow(errorsNegSummedDoubleRatio[l],0.5);
      errorsNegErrSummedDoubleRatio[l] = errorsNegSummedDoubleRatio[l]*0.001;

      errorsPosCorrMatSummedDoubleRatio[l] = pow(errorsPosCorrSummedDoubleRatio[l]+ pow(4.50 ,2.),0.5);
      errorsMeanCorrMatSummedDoubleRatio[l] = pow(errorsMeanCorrSummedDoubleRatio[l]+ pow(4.50 ,2.),0.5);
      errorsNegCorrMatSummedDoubleRatio[l] = -pow(errorsNegCorrSummedDoubleRatio[l]+ pow(4.50 ,2.),0.5);

      errorsPosCorrSummedDoubleRatio[l] = pow(errorsPosCorrSummedDoubleRatio[l],0.5);
      errorsMeanCorrSummedDoubleRatio[l] = pow(errorsMeanCorrSummedDoubleRatio[l],0.5);
      errorsPosErrCorrSummedDoubleRatio[l] = errorsPosCorrSummedDoubleRatio[l]*0.001;
      errorsMeanErrCorrSummedDoubleRatio[l] = errorsMeanCorrSummedDoubleRatio[l]*0.001;
      errorsMeanErrCorrMatSummedDoubleRatio[l] = errorsMeanCorrMatSummedDoubleRatio[l]*0.001;
      errorsNegCorrSummedDoubleRatio[l] = -pow(errorsNegCorrSummedDoubleRatio[l],0.5);
      errorsNegErrCorrSummedDoubleRatio[l] = errorsNegCorrSummedDoubleRatio[l]*0.001;

      errorsPosSummedDoubleRatioA[l] = pow(errorsPosSummedDoubleRatioA[l],0.5);
      errorsMeanSummedDoubleRatioA[l] = pow(errorsMeanSummedDoubleRatioA[l],0.5);
      errorsPosErrSummedDoubleRatioA[l] = errorsPosSummedDoubleRatioA[l]*0.001;
      errorsMeanErrSummedDoubleRatioA[l] = errorsMeanSummedDoubleRatioA[l]*0.001;
      errorsNegSummedDoubleRatioA[l] = -pow(errorsNegSummedDoubleRatioA[l],0.5);
      errorsNegErrSummedDoubleRatioA[l] = errorsNegSummedDoubleRatioA[l]*0.001;

      errorsPosCorrMatSummedDoubleRatioA[l] = pow(errorsPosCorrSummedDoubleRatioA[l],0.5);
      errorsMeanCorrMatSummedDoubleRatioA[l] = pow(errorsMeanCorrSummedDoubleRatioA[l],0.5);
      errorsNegCorrMatSummedDoubleRatioA[l] = -pow(errorsNegCorrSummedDoubleRatioA[l],0.5);

      errorsPosCorrSummedDoubleRatioA[l] = pow(errorsPosCorrSummedDoubleRatioA[l],0.5);
      errorsMeanCorrSummedDoubleRatioA[l] = pow(errorsMeanCorrSummedDoubleRatioA[l],0.5);
      errorsPosErrCorrSummedDoubleRatioA[l] = errorsPosCorrSummedDoubleRatioA[l]*0.001;
      errorsMeanErrCorrSummedDoubleRatioA[l] = errorsMeanCorrSummedDoubleRatioA[l]*0.001;
      errorsMeanErrCorrMatSummedDoubleRatioA[l] = errorsMeanCorrMatSummedDoubleRatioA[l]*0.001;
      errorsNegCorrSummedDoubleRatioA[l] = -pow(errorsNegCorrSummedDoubleRatioA[l],0.5);
      errorsNegErrCorrSummedDoubleRatioA[l] = errorsNegCorrSummedDoubleRatioA[l]*0.001;

      errorsPosSummedDoubleRatioB[l] = pow(errorsPosSummedDoubleRatioB[l],0.5);
      errorsMeanSummedDoubleRatioB[l] = pow(errorsMeanSummedDoubleRatioB[l],0.5);
      errorsPosErrSummedDoubleRatioB[l] = errorsPosSummedDoubleRatioB[l]*0.001;
      errorsMeanErrSummedDoubleRatioB[l] = errorsMeanSummedDoubleRatioB[l]*0.001;
      errorsNegSummedDoubleRatioB[l] = -pow(errorsNegSummedDoubleRatioB[l],0.5);
      errorsNegErrSummedDoubleRatioB[l] = errorsNegSummedDoubleRatioB[l]*0.001;

      errorsPosCorrMatSummedDoubleRatioB[l] = pow(errorsPosCorrSummedDoubleRatioB[l],0.5);
      errorsMeanCorrMatSummedDoubleRatioB[l] = pow(errorsMeanCorrSummedDoubleRatioB[l],0.5);
      errorsNegCorrMatSummedDoubleRatioB[l] = -pow(errorsNegCorrSummedDoubleRatioB[l],0.5);

      errorsPosCorrSummedDoubleRatioB[l] = pow(errorsPosCorrSummedDoubleRatioB[l],0.5);
      errorsMeanCorrSummedDoubleRatioB[l] = pow(errorsMeanCorrSummedDoubleRatioB[l],0.5);
      errorsPosErrCorrSummedDoubleRatioB[l] = errorsPosCorrSummedDoubleRatioB[l]*0.001;
      errorsMeanErrCorrSummedDoubleRatioB[l] = errorsMeanCorrSummedDoubleRatioB[l]*0.001;
      errorsMeanErrCorrMatSummedDoubleRatioB[l] = errorsMeanCorrMatSummedDoubleRatioB[l]*0.001;
      errorsNegCorrSummedDoubleRatioB[l] = -pow(errorsNegCorrSummedDoubleRatioB[l],0.5);
      errorsNegErrCorrSummedDoubleRatioB[l] = errorsNegCorrSummedDoubleRatioB[l]*0.001;

      errorsPosSummedDoubleRatioC[l] = pow(errorsPosSummedDoubleRatioC[l],0.5);
      errorsMeanSummedDoubleRatioC[l] = pow(errorsMeanSummedDoubleRatioC[l],0.5);
      errorsPosErrSummedDoubleRatioC[l] = errorsPosSummedDoubleRatioC[l]*0.001;
      errorsMeanErrSummedDoubleRatioC[l] = errorsMeanSummedDoubleRatioC[l]*0.001;
      errorsNegSummedDoubleRatioC[l] = -pow(errorsNegSummedDoubleRatioC[l],0.5);
      errorsNegErrSummedDoubleRatioC[l] = errorsNegSummedDoubleRatioC[l]*0.001;

      errorsPosCorrMatSummedDoubleRatioC[l] = pow(errorsPosCorrSummedDoubleRatioC[l]+ pow(4.50 ,2.),0.5);
      errorsMeanCorrMatSummedDoubleRatioC[l] = pow(errorsMeanCorrSummedDoubleRatioC[l]+ pow(4.50 ,2.),0.5);
      errorsNegCorrMatSummedDoubleRatioC[l] = -pow(errorsNegCorrSummedDoubleRatioC[l]+ pow(4.50 ,2.),0.5);

      errorsPosCorrSummedDoubleRatioC[l] = pow(errorsPosCorrSummedDoubleRatioC[l],0.5);
      errorsMeanCorrSummedDoubleRatioC[l] = pow(errorsMeanCorrSummedDoubleRatioC[l],0.5);
      errorsPosErrCorrSummedDoubleRatioC[l] = errorsPosCorrSummedDoubleRatioC[l]*0.001;
      errorsMeanErrCorrSummedDoubleRatioC[l] = errorsMeanCorrSummedDoubleRatioC[l]*0.001;
      errorsMeanErrCorrMatSummedDoubleRatioC[l] = errorsMeanCorrMatSummedDoubleRatioC[l]*0.001;
      errorsNegCorrSummedDoubleRatioC[l] = -pow(errorsNegCorrSummedDoubleRatioC[l],0.5);
      errorsNegErrCorrSummedDoubleRatioC[l] = errorsNegCorrSummedDoubleRatioC[l]*0.001;


      errorsPosSummedDoubleRatioPi0Fit[l] = pow(errorsPosSummedDoubleRatioPi0Fit[l],0.5);
      errorsMeanSummedDoubleRatioPi0Fit[l] = pow(errorsMeanSummedDoubleRatioPi0Fit[l],0.5);
      errorsPosErrSummedDoubleRatioPi0Fit[l] = errorsPosSummedDoubleRatioPi0Fit[l]*0.001;
      errorsMeanErrSummedDoubleRatioPi0Fit[l] = errorsMeanSummedDoubleRatioPi0Fit[l]*0.001;
      errorsNegSummedDoubleRatioPi0Fit[l] = -pow(errorsNegSummedDoubleRatioPi0Fit[l],0.5);
      errorsNegErrSummedDoubleRatioPi0Fit[l] = errorsNegSummedDoubleRatioPi0Fit[l]*0.001;

      errorsPosCorrMatSummedDoubleRatioPi0Fit[l] = pow(errorsPosCorrSummedDoubleRatioPi0Fit[l]+ pow(4.50 ,2.),0.5);
      errorsMeanCorrMatSummedDoubleRatioPi0Fit[l] = pow(errorsMeanCorrSummedDoubleRatioPi0Fit[l]+ pow(4.50 ,2.),0.5);
      errorsNegCorrMatSummedDoubleRatioPi0Fit[l] = -pow(errorsNegCorrSummedDoubleRatioPi0Fit[l]+ pow(4.50 ,2.),0.5);

      errorsPosCorrSummedDoubleRatioPi0Fit[l] = pow(errorsPosCorrSummedDoubleRatioPi0Fit[l],0.5);
      errorsMeanCorrSummedDoubleRatioPi0Fit[l] = pow(errorsMeanCorrSummedDoubleRatioPi0Fit[l],0.5);
      errorsPosErrCorrSummedDoubleRatioPi0Fit[l] = errorsPosCorrSummedDoubleRatioPi0Fit[l]*0.001;
      errorsMeanErrCorrSummedDoubleRatioPi0Fit[l] = errorsMeanCorrSummedDoubleRatioPi0Fit[l]*0.001;
      errorsMeanErrCorrMatSummedDoubleRatioPi0Fit[l] = errorsMeanCorrMatSummedDoubleRatioPi0Fit[l]*0.001;
      errorsNegCorrSummedDoubleRatioPi0Fit[l] = -pow(errorsNegCorrSummedDoubleRatioPi0Fit[l],0.5);
      errorsNegErrCorrSummedDoubleRatioPi0Fit[l] = errorsNegCorrSummedDoubleRatioPi0Fit[l]*0.001;

      errorsMeanErrCorrMatSummedDoubleRatioPi0FitWithoutMat[l] = errorsMeanCorrMatSummedDoubleRatioPi0FitWithoutMat[l]*0.001;
      errorsErrMaterialBudgetDoubleRatioPi0Fit[l] = errorsMaterialBudgetDoubleRatioPi0Fit[l]*0.001;


      errorsPosSummedDoubleRatioPi0FitA[l] = pow(errorsPosSummedDoubleRatioPi0FitA[l],0.5);
      errorsMeanSummedDoubleRatioPi0FitA[l] = pow(errorsMeanSummedDoubleRatioPi0FitA[l],0.5);
      errorsPosErrSummedDoubleRatioPi0FitA[l] = errorsPosSummedDoubleRatioPi0FitA[l]*0.001;
      errorsMeanErrSummedDoubleRatioPi0FitA[l] = errorsMeanSummedDoubleRatioPi0FitA[l]*0.001;
      errorsNegSummedDoubleRatioPi0FitA[l] = -pow(errorsNegSummedDoubleRatioPi0FitA[l],0.5);
      errorsNegErrSummedDoubleRatioPi0FitA[l] = errorsNegSummedDoubleRatioPi0FitA[l]*0.001;

      errorsPosCorrMatSummedDoubleRatioPi0FitA[l] = pow(errorsPosCorrSummedDoubleRatioPi0FitA[l],0.5);
      errorsMeanCorrMatSummedDoubleRatioPi0FitA[l] = pow(errorsMeanCorrSummedDoubleRatioPi0FitA[l],0.5);
      errorsNegCorrMatSummedDoubleRatioPi0FitA[l] = -pow(errorsNegCorrSummedDoubleRatioPi0FitA[l],0.5);

      errorsPosCorrSummedDoubleRatioPi0FitA[l] = pow(errorsPosCorrSummedDoubleRatioPi0FitA[l],0.5);
      errorsMeanCorrSummedDoubleRatioPi0FitA[l] = pow(errorsMeanCorrSummedDoubleRatioPi0FitA[l],0.5);
      errorsPosErrCorrSummedDoubleRatioPi0FitA[l] = errorsPosCorrSummedDoubleRatioPi0FitA[l]*0.001;
      errorsMeanErrCorrSummedDoubleRatioPi0FitA[l] = errorsMeanCorrSummedDoubleRatioPi0FitA[l]*0.001;
      errorsMeanErrCorrMatSummedDoubleRatioPi0FitA[l] = errorsMeanCorrMatSummedDoubleRatioPi0FitA[l]*0.001;
      errorsNegCorrSummedDoubleRatioPi0FitA[l] = -pow(errorsNegCorrSummedDoubleRatioPi0FitA[l],0.5);
      errorsNegErrCorrSummedDoubleRatioPi0FitA[l] = errorsNegCorrSummedDoubleRatioPi0FitA[l]*0.001;

      errorsPosSummedDoubleRatioPi0FitB[l] = pow(errorsPosSummedDoubleRatioPi0FitB[l],0.5);
      errorsMeanSummedDoubleRatioPi0FitB[l] = pow(errorsMeanSummedDoubleRatioPi0FitB[l],0.5);
      errorsPosErrSummedDoubleRatioPi0FitB[l] = errorsPosSummedDoubleRatioPi0FitB[l]*0.001;
      errorsMeanErrSummedDoubleRatioPi0FitB[l] = errorsMeanSummedDoubleRatioPi0FitB[l]*0.001;
      errorsNegSummedDoubleRatioPi0FitB[l] = -pow(errorsNegSummedDoubleRatioPi0FitB[l],0.5);
      errorsNegErrSummedDoubleRatioPi0FitB[l] = errorsNegSummedDoubleRatioPi0FitB[l]*0.001;

      errorsPosCorrMatSummedDoubleRatioPi0FitB[l] = pow(errorsPosCorrSummedDoubleRatioPi0FitB[l],0.5);
      errorsMeanCorrMatSummedDoubleRatioPi0FitB[l] = pow(errorsMeanCorrSummedDoubleRatioPi0FitB[l],0.5);
      errorsNegCorrMatSummedDoubleRatioPi0FitB[l] = -pow(errorsNegCorrSummedDoubleRatioPi0FitB[l],0.5);

      errorsPosCorrSummedDoubleRatioPi0FitB[l] = pow(errorsPosCorrSummedDoubleRatioPi0FitB[l],0.5);
      errorsMeanCorrSummedDoubleRatioPi0FitB[l] = pow(errorsMeanCorrSummedDoubleRatioPi0FitB[l],0.5);
      errorsPosErrCorrSummedDoubleRatioPi0FitB[l] = errorsPosCorrSummedDoubleRatioPi0FitB[l]*0.001;
      errorsMeanErrCorrSummedDoubleRatioPi0FitB[l] = errorsMeanCorrSummedDoubleRatioPi0FitB[l]*0.001;
      errorsMeanErrCorrMatSummedDoubleRatioPi0FitB[l] = errorsMeanCorrMatSummedDoubleRatioPi0FitB[l]*0.001;
      errorsNegCorrSummedDoubleRatioPi0FitB[l] = -pow(errorsNegCorrSummedDoubleRatioPi0FitB[l],0.5);
      errorsNegErrCorrSummedDoubleRatioPi0FitB[l] = errorsNegCorrSummedDoubleRatioPi0FitB[l]*0.001;

      errorsPosSummedDoubleRatioPi0FitC[l] = pow(errorsPosSummedDoubleRatioPi0FitC[l],0.5);
      errorsMeanSummedDoubleRatioPi0FitC[l] = pow(errorsMeanSummedDoubleRatioPi0FitC[l],0.5);
      errorsPosErrSummedDoubleRatioPi0FitC[l] = errorsPosSummedDoubleRatioPi0FitC[l]*0.001;
      errorsMeanErrSummedDoubleRatioPi0FitC[l] = errorsMeanSummedDoubleRatioPi0FitC[l]*0.001;
      errorsNegSummedDoubleRatioPi0FitC[l] = -pow(errorsNegSummedDoubleRatioPi0FitC[l],0.5);
      errorsNegErrSummedDoubleRatioPi0FitC[l] = errorsNegSummedDoubleRatioPi0FitC[l]*0.001;

      errorsPosCorrMatSummedDoubleRatioPi0FitC[l] = pow(errorsPosCorrSummedDoubleRatioPi0FitC[l]+ pow(4.50 ,2.),0.5);
      errorsMeanCorrMatSummedDoubleRatioPi0FitC[l] = pow(errorsMeanCorrSummedDoubleRatioPi0FitC[l]+ pow(4.50 ,2.),0.5);
      errorsNegCorrMatSummedDoubleRatioPi0FitC[l] = -pow(errorsNegCorrSummedDoubleRatioPi0FitC[l]+ pow(4.50 ,2.),0.5);

      errorsPosCorrSummedDoubleRatioPi0FitC[l] = pow(errorsPosCorrSummedDoubleRatioPi0FitC[l],0.5);
      errorsMeanCorrSummedDoubleRatioPi0FitC[l] = pow(errorsMeanCorrSummedDoubleRatioPi0FitC[l],0.5);
      errorsPosErrCorrSummedDoubleRatioPi0FitC[l] = errorsPosCorrSummedDoubleRatioPi0FitC[l]*0.001;
      errorsMeanErrCorrSummedDoubleRatioPi0FitC[l] = errorsMeanCorrSummedDoubleRatioPi0FitC[l]*0.001;
      errorsMeanErrCorrMatSummedDoubleRatioPi0FitC[l] = errorsMeanCorrMatSummedDoubleRatioPi0FitC[l]*0.001;
      errorsNegCorrSummedDoubleRatioPi0FitC[l] = -pow(errorsNegCorrSummedDoubleRatioPi0FitC[l],0.5);
      errorsNegErrCorrSummedDoubleRatioPi0FitC[l] = errorsNegCorrSummedDoubleRatioPi0FitC[l]*0.001;

   }

   for (Int_t l = 0; l <ConstnBinsIncRatio; l++){

      errorsPosSummedIncRatio[l] = pow(errorsPosSummedIncRatio[l],0.5);
      errorsMeanSummedIncRatio[l] = pow(errorsMeanSummedIncRatio[l],0.5);
      errorsPosErrSummedIncRatio[l] = errorsPosSummedIncRatio[l]*0.001;
      errorsMeanErrSummedIncRatio[l] = errorsMeanSummedIncRatio[l]*0.001;
      errorsNegSummedIncRatio[l] = -pow(errorsNegSummedIncRatio[l],0.5);
      errorsNegErrSummedIncRatio[l] = errorsNegSummedIncRatio[l]*0.001;
      errorsPosCorrMatSummedIncRatio[l] = pow(errorsPosCorrSummedIncRatio[l]+ pow(4.50 ,2.),0.5);
      errorsMeanCorrMatSummedIncRatio[l] = pow(errorsMeanCorrSummedIncRatio[l]+ pow(4.50 ,2.),0.5);
      errorsNegCorrMatSummedIncRatio[l] = -pow(errorsNegCorrSummedIncRatio[l]+ pow(4.50 ,2.),0.5);

      errorsPosCorrSummedIncRatio[l] = pow(errorsPosCorrSummedIncRatio[l],0.5);
      errorsMeanCorrSummedIncRatio[l] = pow(errorsMeanCorrSummedIncRatio[l],0.5);
      errorsPosErrCorrSummedIncRatio[l] = errorsPosCorrSummedIncRatio[l]*0.001;
      errorsMeanErrCorrSummedIncRatio[l] = errorsMeanCorrSummedIncRatio[l]*0.001;
      errorsMeanErrCorrMatSummedIncRatio[l] = errorsMeanCorrMatSummedIncRatio[l]*0.001;
      errorsNegCorrSummedIncRatio[l] = -pow(errorsNegCorrSummedIncRatio[l],0.5);
      errorsNegErrCorrSummedIncRatio[l] = errorsNegCorrSummedIncRatio[l]*0.001;


      errorsPosSummedIncRatioA[l] = pow(errorsPosSummedIncRatioA[l],0.5);
      errorsMeanSummedIncRatioA[l] = pow(errorsMeanSummedIncRatioA[l],0.5);
      errorsPosErrSummedIncRatioA[l] = errorsPosSummedIncRatioA[l]*0.001;
      errorsMeanErrSummedIncRatioA[l] = errorsMeanSummedIncRatioA[l]*0.001;
      errorsNegSummedIncRatioA[l] = -pow(errorsNegSummedIncRatioA[l],0.5);
      errorsNegErrSummedIncRatioA[l] = errorsNegSummedIncRatioA[l]*0.001;
      errorsPosCorrMatSummedIncRatioA[l] = pow(errorsPosCorrSummedIncRatioA[l],0.5);
      errorsMeanCorrMatSummedIncRatioA[l] = pow(errorsMeanCorrSummedIncRatioA[l],0.5);
      errorsNegCorrMatSummedIncRatioA[l] = -pow(errorsNegCorrSummedIncRatioA[l],0.5);

      errorsPosCorrSummedIncRatioA[l] = pow(errorsPosCorrSummedIncRatioA[l],0.5);
      errorsMeanCorrSummedIncRatioA[l] = pow(errorsMeanCorrSummedIncRatioA[l],0.5);
      errorsPosErrCorrSummedIncRatioA[l] = errorsPosCorrSummedIncRatioA[l]*0.001;
      errorsMeanErrCorrSummedIncRatioA[l] = errorsMeanCorrSummedIncRatioA[l]*0.001;
      errorsMeanErrCorrMatSummedIncRatioA[l] = errorsMeanCorrMatSummedIncRatioA[l]*0.001;
      errorsNegCorrSummedIncRatioA[l] = -pow(errorsNegCorrSummedIncRatioA[l],0.5);
      errorsNegErrCorrSummedIncRatioA[l] = errorsNegCorrSummedIncRatioA[l]*0.001;

      errorsPosSummedIncRatioB[l] = pow(errorsPosSummedIncRatioB[l],0.5);
      errorsMeanSummedIncRatioB[l] = pow(errorsMeanSummedIncRatioB[l],0.5);
      errorsPosErrSummedIncRatioB[l] = errorsPosSummedIncRatioB[l]*0.001;
      errorsMeanErrSummedIncRatioB[l] = errorsMeanSummedIncRatioB[l]*0.001;
      errorsNegSummedIncRatioB[l] = -pow(errorsNegSummedIncRatioB[l],0.5);
      errorsNegErrSummedIncRatioB[l] = errorsNegSummedIncRatioB[l]*0.001;
      errorsPosCorrMatSummedIncRatioB[l] = pow(errorsPosCorrSummedIncRatioB[l],0.5);
      errorsMeanCorrMatSummedIncRatioB[l] = pow(errorsMeanCorrSummedIncRatioB[l],0.5);
      errorsNegCorrMatSummedIncRatioB[l] = -pow(errorsNegCorrSummedIncRatioB[l],0.5);

      errorsPosCorrSummedIncRatioB[l] = pow(errorsPosCorrSummedIncRatioB[l],0.5);
      errorsMeanCorrSummedIncRatioB[l] = pow(errorsMeanCorrSummedIncRatioB[l],0.5);
      errorsPosErrCorrSummedIncRatioB[l] = errorsPosCorrSummedIncRatioB[l]*0.001;
      errorsMeanErrCorrSummedIncRatioB[l] = errorsMeanCorrSummedIncRatioB[l]*0.001;
      errorsMeanErrCorrMatSummedIncRatioB[l] = errorsMeanCorrMatSummedIncRatioB[l]*0.001;
      errorsNegCorrSummedIncRatioB[l] = -pow(errorsNegCorrSummedIncRatioB[l],0.5);
      errorsNegErrCorrSummedIncRatioB[l] = errorsNegCorrSummedIncRatioB[l]*0.001;

      errorsPosSummedIncRatioC[l] = pow(errorsPosSummedIncRatioC[l],0.5);
      errorsMeanSummedIncRatioC[l] = pow(errorsMeanSummedIncRatioC[l],0.5);
      errorsPosErrSummedIncRatioC[l] = errorsPosSummedIncRatioC[l]*0.001;
      errorsMeanErrSummedIncRatioC[l] = errorsMeanSummedIncRatioC[l]*0.001;
      errorsNegSummedIncRatioC[l] = -pow(errorsNegSummedIncRatioC[l],0.5);
      errorsNegErrSummedIncRatioC[l] = errorsNegSummedIncRatioC[l]*0.001;
      errorsPosCorrMatSummedIncRatioC[l] = pow(errorsPosCorrSummedIncRatioC[l]+ pow(4.50 ,2.),0.5);
      errorsMeanCorrMatSummedIncRatioC[l] = pow(errorsMeanCorrSummedIncRatioC[l]+ pow(4.50 ,2.),0.5);
      errorsNegCorrMatSummedIncRatioC[l] = -pow(errorsNegCorrSummedIncRatioC[l]+ pow(4.50 ,2.),0.5);

      errorsPosCorrSummedIncRatioC[l] = pow(errorsPosCorrSummedIncRatioC[l],0.5);
      errorsMeanCorrSummedIncRatioC[l] = pow(errorsMeanCorrSummedIncRatioC[l],0.5);
      errorsPosErrCorrSummedIncRatioC[l] = errorsPosCorrSummedIncRatioC[l]*0.001;
      errorsMeanErrCorrSummedIncRatioC[l] = errorsMeanCorrSummedIncRatioC[l]*0.001;
      errorsMeanErrCorrMatSummedIncRatioC[l] = errorsMeanCorrMatSummedIncRatioC[l]*0.001;
      errorsNegCorrSummedIncRatioC[l] = -pow(errorsNegCorrSummedIncRatioC[l],0.5);
      errorsNegErrCorrSummedIncRatioC[l] = errorsNegCorrSummedIncRatioC[l]*0.001;

      errorsPosSummedIncRatioPi0Fit[l] = pow(errorsPosSummedIncRatioPi0Fit[l],0.5);
      errorsMeanSummedIncRatioPi0Fit[l] = pow(errorsMeanSummedIncRatioPi0Fit[l],0.5);
      errorsPosErrSummedIncRatioPi0Fit[l] = errorsPosSummedIncRatioPi0Fit[l]*0.001;
      errorsMeanErrSummedIncRatioPi0Fit[l] = errorsMeanSummedIncRatioPi0Fit[l]*0.001;
      errorsNegSummedIncRatioPi0Fit[l] = -pow(errorsNegSummedIncRatioPi0Fit[l],0.5);
      errorsNegErrSummedIncRatioPi0Fit[l] = errorsNegSummedIncRatioPi0Fit[l]*0.001;
      errorsPosCorrMatSummedIncRatioPi0Fit[l] = pow(errorsPosCorrSummedIncRatioPi0Fit[l]+ pow(4.50 ,2.),0.5);
      errorsMeanCorrMatSummedIncRatioPi0Fit[l] = pow(errorsMeanCorrSummedIncRatioPi0Fit[l]+ pow(4.50 ,2.),0.5);
      errorsNegCorrMatSummedIncRatioPi0Fit[l] = -pow(errorsNegCorrSummedIncRatioPi0Fit[l]+ pow(4.50 ,2.),0.5);

      errorsPosCorrSummedIncRatioPi0Fit[l] = pow(errorsPosCorrSummedIncRatioPi0Fit[l],0.5);
      errorsMeanCorrSummedIncRatioPi0Fit[l] = pow(errorsMeanCorrSummedIncRatioPi0Fit[l],0.5);
      errorsPosErrCorrSummedIncRatioPi0Fit[l] = errorsPosCorrSummedIncRatioPi0Fit[l]*0.001;
      errorsMeanErrCorrSummedIncRatioPi0Fit[l] = errorsMeanCorrSummedIncRatioPi0Fit[l]*0.001;
      errorsMeanErrCorrMatSummedIncRatioPi0Fit[l] = errorsMeanCorrMatSummedIncRatioPi0Fit[l]*0.001;
      errorsNegCorrSummedIncRatioPi0Fit[l] = -pow(errorsNegCorrSummedIncRatioPi0Fit[l],0.5);
      errorsNegErrCorrSummedIncRatioPi0Fit[l] = errorsNegCorrSummedIncRatioPi0Fit[l]*0.001;


      errorsPosSummedIncRatioPi0FitA[l] = pow(errorsPosSummedIncRatioPi0FitA[l],0.5);
      errorsMeanSummedIncRatioPi0FitA[l] = pow(errorsMeanSummedIncRatioPi0FitA[l],0.5);
      errorsPosErrSummedIncRatioPi0FitA[l] = errorsPosSummedIncRatioPi0FitA[l]*0.001;
      errorsMeanErrSummedIncRatioPi0FitA[l] = errorsMeanSummedIncRatioPi0FitA[l]*0.001;
      errorsNegSummedIncRatioPi0FitA[l] = -pow(errorsNegSummedIncRatioPi0FitA[l],0.5);
      errorsNegErrSummedIncRatioPi0FitA[l] = errorsNegSummedIncRatioPi0FitA[l]*0.001;
      errorsPosCorrMatSummedIncRatioPi0FitA[l] = pow(errorsPosCorrSummedIncRatioPi0FitA[l],0.5);
      errorsMeanCorrMatSummedIncRatioPi0FitA[l] = pow(errorsMeanCorrSummedIncRatioPi0FitA[l],0.5);
      errorsNegCorrMatSummedIncRatioPi0FitA[l] = -pow(errorsNegCorrSummedIncRatioPi0FitA[l],0.5);

      errorsPosCorrSummedIncRatioPi0FitA[l] = pow(errorsPosCorrSummedIncRatioPi0FitA[l],0.5);
      errorsMeanCorrSummedIncRatioPi0FitA[l] = pow(errorsMeanCorrSummedIncRatioPi0FitA[l],0.5);
      errorsPosErrCorrSummedIncRatioPi0FitA[l] = errorsPosCorrSummedIncRatioPi0FitA[l]*0.001;
      errorsMeanErrCorrSummedIncRatioPi0FitA[l] = errorsMeanCorrSummedIncRatioPi0FitA[l]*0.001;
      errorsMeanErrCorrMatSummedIncRatioPi0FitA[l] = errorsMeanCorrMatSummedIncRatioPi0FitA[l]*0.001;
      errorsNegCorrSummedIncRatioPi0FitA[l] = -pow(errorsNegCorrSummedIncRatioPi0FitA[l],0.5);
      errorsNegErrCorrSummedIncRatioPi0FitA[l] = errorsNegCorrSummedIncRatioPi0FitA[l]*0.001;

      errorsPosSummedIncRatioPi0FitB[l] = pow(errorsPosSummedIncRatioPi0FitB[l],0.5);
      errorsMeanSummedIncRatioPi0FitB[l] = pow(errorsMeanSummedIncRatioPi0FitB[l],0.5);
      errorsPosErrSummedIncRatioPi0FitB[l] = errorsPosSummedIncRatioPi0FitB[l]*0.001;
      errorsMeanErrSummedIncRatioPi0FitB[l] = errorsMeanSummedIncRatioPi0FitB[l]*0.001;
      errorsNegSummedIncRatioPi0FitB[l] = -pow(errorsNegSummedIncRatioPi0FitB[l],0.5);
      errorsNegErrSummedIncRatioPi0FitB[l] = errorsNegSummedIncRatioPi0FitB[l]*0.001;
      errorsPosCorrMatSummedIncRatioPi0FitB[l] = pow(errorsPosCorrSummedIncRatioPi0FitB[l],0.5);
      errorsMeanCorrMatSummedIncRatioPi0FitB[l] = pow(errorsMeanCorrSummedIncRatioPi0FitB[l],0.5);
      errorsNegCorrMatSummedIncRatioPi0FitB[l] = -pow(errorsNegCorrSummedIncRatioPi0FitB[l],0.5);

      errorsPosCorrSummedIncRatioPi0FitB[l] = pow(errorsPosCorrSummedIncRatioPi0FitB[l],0.5);
      errorsMeanCorrSummedIncRatioPi0FitB[l] = pow(errorsMeanCorrSummedIncRatioPi0FitB[l],0.5);
      errorsPosErrCorrSummedIncRatioPi0FitB[l] = errorsPosCorrSummedIncRatioPi0FitB[l]*0.001;
      errorsMeanErrCorrSummedIncRatioPi0FitB[l] = errorsMeanCorrSummedIncRatioPi0FitB[l]*0.001;
      errorsMeanErrCorrMatSummedIncRatioPi0FitB[l] = errorsMeanCorrMatSummedIncRatioPi0FitB[l]*0.001;
      errorsNegCorrSummedIncRatioPi0FitB[l] = -pow(errorsNegCorrSummedIncRatioPi0FitB[l],0.5);
      errorsNegErrCorrSummedIncRatioPi0FitB[l] = errorsNegCorrSummedIncRatioPi0FitB[l]*0.001;

      errorsPosSummedIncRatioPi0FitC[l] = pow(errorsPosSummedIncRatioPi0FitC[l],0.5);
      errorsMeanSummedIncRatioPi0FitC[l] = pow(errorsMeanSummedIncRatioPi0FitC[l],0.5);
      errorsPosErrSummedIncRatioPi0FitC[l] = errorsPosSummedIncRatioPi0FitC[l]*0.001;
      errorsMeanErrSummedIncRatioPi0FitC[l] = errorsMeanSummedIncRatioPi0FitC[l]*0.001;
      errorsNegSummedIncRatioPi0FitC[l] = -pow(errorsNegSummedIncRatioPi0FitC[l],0.5);
      errorsNegErrSummedIncRatioPi0FitC[l] = errorsNegSummedIncRatioPi0FitC[l]*0.001;
      errorsPosCorrMatSummedIncRatioPi0FitC[l] = pow(errorsPosCorrSummedIncRatioPi0FitC[l]+ pow(4.50 ,2.),0.5);
      errorsMeanCorrMatSummedIncRatioPi0FitC[l] = pow(errorsMeanCorrSummedIncRatioPi0FitC[l]+ pow(4.50 ,2.),0.5);
      errorsNegCorrMatSummedIncRatioPi0FitC[l] = -pow(errorsNegCorrSummedIncRatioPi0FitC[l]+ pow(4.50 ,2.),0.5);

      errorsPosCorrSummedIncRatioPi0FitC[l] = pow(errorsPosCorrSummedIncRatioPi0FitC[l],0.5);
      errorsMeanCorrSummedIncRatioPi0FitC[l] = pow(errorsMeanCorrSummedIncRatioPi0FitC[l],0.5);
      errorsPosErrCorrSummedIncRatioPi0FitC[l] = errorsPosCorrSummedIncRatioPi0FitC[l]*0.001;
      errorsMeanErrCorrSummedIncRatioPi0FitC[l] = errorsMeanCorrSummedIncRatioPi0FitC[l]*0.001;
      errorsMeanErrCorrMatSummedIncRatioPi0FitC[l] = errorsMeanCorrMatSummedIncRatioPi0FitC[l]*0.001;
      errorsNegCorrSummedIncRatioPi0FitC[l] = -pow(errorsNegCorrSummedIncRatioPi0FitC[l],0.5);
      errorsNegErrCorrSummedIncRatioPi0FitC[l] = errorsNegCorrSummedIncRatioPi0FitC[l]*0.001;



   }

   for (Int_t l = 0; l <ConstnBinsPi0; l++){

      errorsPosSummedPi0[l] = pow(errorsPosSummedPi0[l],0.5);
      errorsMeanSummedPi0[l] = pow(errorsMeanSummedPi0[l],0.5);
      errorsPosErrSummedPi0[l] = errorsPosSummedPi0[l]*0.001;
      errorsMeanErrSummedPi0[l] = errorsMeanSummedPi0[l]*0.001;
      errorsNegSummedPi0[l] = -pow(errorsNegSummedPi0[l],0.5);
      errorsNegErrSummedPi0[l] = errorsNegSummedPi0[l]*0.001;
      errorsPosCorrMatSummedPi0[l] = pow(errorsPosCorrSummedPi0[l]+ pow(9.0 ,2.),0.5);
      errorsMeanCorrMatSummedPi0[l] = pow(errorsMeanCorrSummedPi0[l]+ pow(9.0 ,2.),0.5);
      errorsNegCorrMatSummedPi0[l] = -pow(errorsNegCorrSummedPi0[l]+ pow(9.0 ,2.),0.5);
      errorsMeanCorrMatSummedPi0WithoutMaterial[l] = pow(errorsMeanCorrSummedPi0[l],0.5);


      errorsPosCorrSummedPi0[l] = pow(errorsPosCorrSummedPi0[l],0.5);
      errorsMeanCorrSummedPi0[l] = pow(errorsMeanCorrSummedPi0[l],0.5);
      errorsPosErrCorrSummedPi0[l] = errorsPosCorrSummedPi0[l]*0.001;
      errorsMeanErrCorrSummedPi0[l] = errorsMeanCorrSummedPi0[l]*0.001;
      errorsMeanErrCorrMatSummedPi0[l] = errorsMeanCorrMatSummedPi0[l]*0.001;
      errorsNegCorrSummedPi0[l] = -pow(errorsNegCorrSummedPi0[l],0.5);
      errorsNegErrCorrSummedPi0[l] = errorsNegCorrSummedPi0[l]*0.001;

      errorsPosSummedPi0Fit[l] = pow(errorsPosSummedPi0Fit[l],0.5);
      errorsMeanSummedPi0Fit[l] = pow(errorsMeanSummedPi0Fit[l],0.5);
      errorsPosErrSummedPi0Fit[l] = errorsPosSummedPi0Fit[l]*0.001;
      errorsMeanErrSummedPi0Fit[l] = errorsMeanSummedPi0Fit[l]*0.001;
      errorsNegSummedPi0Fit[l] = -pow(errorsNegSummedPi0Fit[l],0.5);
      errorsNegErrSummedPi0Fit[l] = errorsNegSummedPi0Fit[l]*0.001;
      errorsPosCorrMatSummedPi0Fit[l] = pow(errorsPosCorrSummedPi0Fit[l]+ pow(9.0 ,2.),0.5);
      errorsMeanCorrMatSummedPi0Fit[l] = pow(errorsMeanCorrSummedPi0Fit[l]+ pow(9.0 ,2.),0.5);
      errorsNegCorrMatSummedPi0Fit[l] = -pow(errorsNegCorrSummedPi0Fit[l]+ pow(9.0 ,2.),0.5);
      errorsMeanCorrMatSummedPi0FitWithoutMaterial[l] = pow(errorsMeanCorrSummedPi0Fit[l],0.5);


      errorsPosCorrSummedPi0Fit[l] = pow(errorsPosCorrSummedPi0Fit[l],0.5);
      errorsMeanCorrSummedPi0Fit[l] = pow(errorsMeanCorrSummedPi0Fit[l],0.5);
      errorsPosErrCorrSummedPi0Fit[l] = errorsPosCorrSummedPi0Fit[l]*0.001;
      errorsMeanErrCorrSummedPi0Fit[l] = errorsMeanCorrSummedPi0Fit[l]*0.001;
      errorsMeanErrCorrMatSummedPi0Fit[l] = errorsMeanCorrMatSummedPi0Fit[l]*0.001;
      errorsNegCorrSummedPi0Fit[l] = -pow(errorsNegCorrSummedPi0Fit[l],0.5);
      errorsNegErrCorrSummedPi0Fit[l] = errorsNegCorrSummedPi0Fit[l]*0.001;
   }


   for (Int_t l = 0; l <ConstnBinsGamma; l++){
      errorsPosSummedGamma[l] = pow(errorsPosSummedGamma[l],0.5);
      errorsMeanSummedGamma[l] = pow(errorsMeanSummedGamma[l],0.5);
      errorsPosErrSummedGamma[l] = errorsPosSummedGamma[l]*0.001;
      errorsMeanErrSummedGamma[l] = errorsMeanSummedGamma[l]*0.001;
      errorsNegSummedGamma[l] = -pow(errorsNegSummedGamma[l],0.5);
      errorsNegErrSummedGamma[l] = errorsNegSummedGamma[l]*0.001;
      errorsPosCorrMatSummedGamma[l] = pow(errorsPosCorrSummedGamma[l]+ pow(4.50 ,2.),0.5);
      errorsMeanCorrMatSummedGamma[l] = pow(errorsMeanCorrSummedGamma[l]+ pow(4.50 ,2.),0.5);
      errorsNegCorrMatSummedGamma[l] = -pow(errorsNegCorrSummedGamma[l]+ pow(4.50 ,2.),0.5);

      errorsMeanCorrMatSummedGammaWithoutMat[l] =  pow(errorsMeanCorrSummedGamma[l],0.5);
      errorsMaterialBudget[l] = 4.5;

      errorsPosCorrSummedGamma[l] = pow(errorsPosCorrSummedGamma[l],0.5);
      errorsMeanCorrSummedGamma[l] = pow(errorsMeanCorrSummedGamma[l],0.5);
      errorsPosErrCorrSummedGamma[l] = errorsPosCorrSummedGamma[l]*0.001;
      errorsMeanErrCorrSummedGamma[l] = errorsMeanCorrSummedGamma[l]*0.001;
      errorsMeanErrCorrMatSummedGamma[l] = errorsMeanCorrMatSummedGamma[l]*0.001;
      errorsNegCorrSummedGamma[l] = -pow(errorsNegCorrSummedGamma[l],0.5);
      errorsNegErrCorrSummedGamma[l] = errorsNegCorrSummedGamma[l]*0.001;

      errorsMeanErrCorrMatSummedGammaWithoutMat[l] = errorsMeanCorrMatSummedGammaWithoutMat[l]*0.001;
      errorsErrMaterialBudget[l] = errorsMaterialBudget[l]*0.001;

      errorsPosSummedGammaA[l] = pow(errorsPosSummedGammaA[l],0.5);
      errorsMeanSummedGammaA[l] = pow(errorsMeanSummedGammaA[l],0.5);
      errorsPosErrSummedGammaA[l] = errorsPosSummedGammaA[l]*0.001;
      errorsMeanErrSummedGammaA[l] = errorsMeanSummedGammaA[l]*0.001;
      errorsNegSummedGammaA[l] = -pow(errorsNegSummedGammaA[l],0.5);
      errorsNegErrSummedGammaA[l] = errorsNegSummedGammaA[l]*0.001;
      errorsPosCorrMatSummedGammaA[l] = pow(errorsPosCorrSummedGammaA[l],0.5);
      errorsMeanCorrMatSummedGammaA[l] = pow(errorsMeanCorrSummedGammaA[l],0.5);
      errorsNegCorrMatSummedGammaA[l] = -pow(errorsNegCorrSummedGammaA[l],0.5);

      errorsPosCorrSummedGammaA[l] = pow(errorsPosCorrSummedGammaA[l],0.5);
      errorsMeanCorrSummedGammaA[l] = pow(errorsMeanCorrSummedGammaA[l],0.5);
      errorsPosErrCorrSummedGammaA[l] = errorsPosCorrSummedGammaA[l]*0.001;
      errorsMeanErrCorrSummedGammaA[l] = errorsMeanCorrSummedGammaA[l]*0.001;
      errorsMeanErrCorrMatSummedGammaA[l] = errorsMeanCorrMatSummedGammaA[l]*0.001;
      errorsNegCorrSummedGammaA[l] = -pow(errorsNegCorrSummedGammaA[l],0.5);
      errorsNegErrCorrSummedGammaA[l] = errorsNegCorrSummedGammaA[l]*0.001;

      errorsPosSummedGammaB[l] = pow(errorsPosSummedGammaB[l],0.5);
      errorsMeanSummedGammaB[l] = pow(errorsMeanSummedGammaB[l],0.5);
      errorsPosErrSummedGammaB[l] = errorsPosSummedGammaB[l]*0.001;
      errorsMeanErrSummedGammaB[l] = errorsMeanSummedGammaB[l]*0.001;
      errorsNegSummedGammaB[l] = -pow(errorsNegSummedGammaB[l],0.5);
      errorsNegErrSummedGammaB[l] = errorsNegSummedGammaB[l]*0.001;
      errorsPosCorrMatSummedGammaB[l] = pow(errorsPosCorrSummedGammaB[l],0.5);
      errorsMeanCorrMatSummedGammaB[l] = pow(errorsMeanCorrSummedGammaB[l],0.5);
      errorsNegCorrMatSummedGammaB[l] = -pow(errorsNegCorrSummedGammaB[l],0.5);

      errorsPosCorrSummedGammaB[l] = pow(errorsPosCorrSummedGammaB[l],0.5);
      errorsMeanCorrSummedGammaB[l] = pow(errorsMeanCorrSummedGammaB[l],0.5);
      errorsPosErrCorrSummedGammaB[l] = errorsPosCorrSummedGammaB[l]*0.001;
      errorsMeanErrCorrSummedGammaB[l] = errorsMeanCorrSummedGammaB[l]*0.001;
      errorsMeanErrCorrMatSummedGammaB[l] = errorsMeanCorrMatSummedGammaB[l]*0.001;
      errorsNegCorrSummedGammaB[l] = -pow(errorsNegCorrSummedGammaB[l],0.5);
      errorsNegErrCorrSummedGammaB[l] = errorsNegCorrSummedGammaB[l]*0.001;

      errorsPosSummedGammaC[l] = pow(errorsPosSummedGammaC[l],0.5);
      errorsMeanSummedGammaC[l] = pow(errorsMeanSummedGammaC[l],0.5);
      errorsPosErrSummedGammaC[l] = errorsPosSummedGammaC[l]*0.001;
      errorsMeanErrSummedGammaC[l] = errorsMeanSummedGammaC[l]*0.001;
      errorsNegSummedGammaC[l] = -pow(errorsNegSummedGammaC[l],0.5);
      errorsNegErrSummedGammaC[l] = errorsNegSummedGammaC[l]*0.001;
      errorsPosCorrMatSummedGammaC[l] = pow(errorsPosCorrSummedGammaC[l]+ pow(4.50 ,2.),0.5);
      errorsMeanCorrMatSummedGammaC[l] = pow(errorsMeanCorrSummedGammaC[l]+ pow(4.50 ,2.),0.5);
      errorsNegCorrMatSummedGammaC[l] = -pow(errorsNegCorrSummedGammaC[l]+ pow(4.50 ,2.),0.5);

      errorsPosCorrSummedGammaC[l] = pow(errorsPosCorrSummedGammaC[l],0.5);
      errorsMeanCorrSummedGammaC[l] = pow(errorsMeanCorrSummedGammaC[l],0.5);
      errorsPosErrCorrSummedGammaC[l] = errorsPosCorrSummedGammaC[l]*0.001;
      errorsMeanErrCorrSummedGammaC[l] = errorsMeanCorrSummedGammaC[l]*0.001;
      errorsMeanErrCorrMatSummedGammaC[l] = errorsMeanCorrMatSummedGammaC[l]*0.001;
      errorsNegCorrSummedGammaC[l] = -pow(errorsNegCorrSummedGammaC[l],0.5);
      errorsNegErrCorrSummedGammaC[l] = errorsNegCorrSummedGammaC[l]*0.001;
   }

   negativeErrorsSummedDoubleRatio = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsNegSummedDoubleRatio ,ptBinsDoubleRatioErr ,errorsNegErrSummedDoubleRatio );
   negativeErrorsCorrSummedDoubleRatio = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsNegCorrSummedDoubleRatio ,ptBinsDoubleRatioErr ,errorsNegErrCorrSummedDoubleRatio );
   positiveErrorsSummedDoubleRatio = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsPosSummedDoubleRatio ,ptBinsDoubleRatioErr ,errorsPosErrSummedDoubleRatio );
   positiveErrorsCorrSummedDoubleRatio = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsPosCorrSummedDoubleRatio ,ptBinsDoubleRatioErr ,errorsPosErrCorrSummedDoubleRatio );
   meanErrorsSummedDoubleRatio = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanSummedDoubleRatio ,ptBinsDoubleRatioErr ,errorsMeanErrSummedDoubleRatio );
   meanErrorsCorrSummedDoubleRatio = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanCorrSummedDoubleRatio ,ptBinsDoubleRatioErr ,errorsMeanErrCorrSummedDoubleRatio );
   meanErrorsCorrSummedIncMatDoubleRatio = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanCorrMatSummedDoubleRatio ,ptBinsDoubleRatioErr ,errorsMeanErrCorrMatSummedDoubleRatio );

   negativeErrorsSummedDoubleRatioPi0Fit = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsNegSummedDoubleRatioPi0Fit ,ptBinsDoubleRatioErr ,errorsNegErrSummedDoubleRatioPi0Fit );
   negativeErrorsCorrSummedDoubleRatioPi0Fit = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsNegCorrSummedDoubleRatioPi0Fit ,ptBinsDoubleRatioErr ,errorsNegErrCorrSummedDoubleRatioPi0Fit );
   positiveErrorsSummedDoubleRatioPi0Fit = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsPosSummedDoubleRatioPi0Fit ,ptBinsDoubleRatioErr ,errorsPosErrSummedDoubleRatioPi0Fit );
   positiveErrorsCorrSummedDoubleRatioPi0Fit = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsPosCorrSummedDoubleRatioPi0Fit ,ptBinsDoubleRatioErr ,errorsPosErrCorrSummedDoubleRatioPi0Fit );
   meanErrorsSummedDoubleRatioPi0Fit = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanSummedDoubleRatioPi0Fit ,ptBinsDoubleRatioErr ,errorsMeanErrSummedDoubleRatioPi0Fit );
   meanErrorsCorrSummedDoubleRatioPi0Fit = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanCorrSummedDoubleRatioPi0Fit ,ptBinsDoubleRatioErr ,errorsMeanErrCorrSummedDoubleRatioPi0Fit );
   meanErrorsCorrSummedIncMatDoubleRatioPi0Fit = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanCorrMatSummedDoubleRatioPi0Fit ,ptBinsDoubleRatioErr ,errorsMeanErrCorrMatSummedDoubleRatioPi0Fit );

   meanErrorsCorrSummedMaterialDoubleRatioPi0FitCocktail = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanCorrMatSummedDoubleRatioPi0Fit ,ptBinsDoubleRatioErr ,errorsMeanErrCorrMatSummedDoubleRatioPi0Fit );
   for(Int_t i = 0;i<meanErrorsCorrSummedMaterialDoubleRatioPi0FitCocktail->GetN();i++){
      Double_t x,y;
      meanErrorsCorrSummedMaterialDoubleRatioPi0FitCocktail->GetPoint(i,x,y);
      meanErrorsCorrSummedMaterialDoubleRatioPi0FitCocktail->SetPoint(i,x,0);
   }

   negativeErrorsSummedDoubleRatioA = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsNegSummedDoubleRatioA ,ptBinsDoubleRatioErr ,errorsNegErrSummedDoubleRatioA );
   negativeErrorsCorrSummedDoubleRatioA = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsNegCorrSummedDoubleRatioA ,ptBinsDoubleRatioErr ,errorsNegErrCorrSummedDoubleRatioA );
   positiveErrorsSummedDoubleRatioA = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsPosSummedDoubleRatioA ,ptBinsDoubleRatioErr ,errorsPosErrSummedDoubleRatioA );
   positiveErrorsCorrSummedDoubleRatioA = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsPosCorrSummedDoubleRatioA ,ptBinsDoubleRatioErr ,errorsPosErrCorrSummedDoubleRatioA );
   meanErrorsSummedDoubleRatioA = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanSummedDoubleRatioA ,ptBinsDoubleRatioErr ,errorsMeanErrSummedDoubleRatioA );
   meanErrorsCorrSummedDoubleRatioA = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanCorrSummedDoubleRatioA ,ptBinsDoubleRatioErr ,errorsMeanErrCorrSummedDoubleRatioA );
   meanErrorsCorrSummedIncMatDoubleRatioA = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanCorrMatSummedDoubleRatioA ,ptBinsDoubleRatioErr ,errorsMeanErrCorrMatSummedDoubleRatioA );

   negativeErrorsSummedDoubleRatioPi0FitA = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsNegSummedDoubleRatioPi0FitA ,ptBinsDoubleRatioErr ,errorsNegErrSummedDoubleRatioPi0FitA );
   negativeErrorsCorrSummedDoubleRatioPi0FitA = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsNegCorrSummedDoubleRatioPi0FitA ,ptBinsDoubleRatioErr ,errorsNegErrCorrSummedDoubleRatioPi0FitA );
   positiveErrorsSummedDoubleRatioPi0FitA = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsPosSummedDoubleRatioPi0FitA ,ptBinsDoubleRatioErr ,errorsPosErrSummedDoubleRatioPi0FitA );
   positiveErrorsCorrSummedDoubleRatioPi0FitA = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsPosCorrSummedDoubleRatioPi0FitA ,ptBinsDoubleRatioErr ,errorsPosErrCorrSummedDoubleRatioPi0FitA );
   meanErrorsSummedDoubleRatioPi0FitA = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanSummedDoubleRatioPi0FitA ,ptBinsDoubleRatioErr ,errorsMeanErrSummedDoubleRatioPi0FitA );
   meanErrorsCorrSummedDoubleRatioPi0FitA = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanCorrSummedDoubleRatioPi0FitA ,ptBinsDoubleRatioErr ,errorsMeanErrCorrSummedDoubleRatioPi0FitA );
   meanErrorsCorrSummedIncMatDoubleRatioPi0FitA = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanCorrMatSummedDoubleRatioPi0FitA ,ptBinsDoubleRatioErr ,errorsMeanErrCorrMatSummedDoubleRatioPi0FitA );

   negativeErrorsSummedDoubleRatioB = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsNegSummedDoubleRatioB ,ptBinsDoubleRatioErr ,errorsNegErrSummedDoubleRatioB );
   negativeErrorsCorrSummedDoubleRatioB = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsNegCorrSummedDoubleRatioB ,ptBinsDoubleRatioErr ,errorsNegErrCorrSummedDoubleRatioB );
   positiveErrorsSummedDoubleRatioB = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsPosSummedDoubleRatioB ,ptBinsDoubleRatioErr ,errorsPosErrSummedDoubleRatioB );
   positiveErrorsCorrSummedDoubleRatioB = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsPosCorrSummedDoubleRatioB ,ptBinsDoubleRatioErr ,errorsPosErrCorrSummedDoubleRatioB );
   meanErrorsSummedDoubleRatioB = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanSummedDoubleRatioB ,ptBinsDoubleRatioErr ,errorsMeanErrSummedDoubleRatioB );
   meanErrorsCorrSummedDoubleRatioB = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanCorrSummedDoubleRatioB ,ptBinsDoubleRatioErr ,errorsMeanErrCorrSummedDoubleRatioB );
   meanErrorsCorrSummedIncMatDoubleRatioB = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanCorrMatSummedDoubleRatioB ,ptBinsDoubleRatioErr ,errorsMeanErrCorrMatSummedDoubleRatioB );

   negativeErrorsSummedDoubleRatioPi0FitB = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsNegSummedDoubleRatioPi0FitB ,ptBinsDoubleRatioErr ,errorsNegErrSummedDoubleRatioPi0FitB );
   negativeErrorsCorrSummedDoubleRatioPi0FitB = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsNegCorrSummedDoubleRatioPi0FitB ,ptBinsDoubleRatioErr ,errorsNegErrCorrSummedDoubleRatioPi0FitB );
   positiveErrorsSummedDoubleRatioPi0FitB = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsPosSummedDoubleRatioPi0FitB ,ptBinsDoubleRatioErr ,errorsPosErrSummedDoubleRatioPi0FitB );
   positiveErrorsCorrSummedDoubleRatioPi0FitB = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsPosCorrSummedDoubleRatioPi0FitB ,ptBinsDoubleRatioErr ,errorsPosErrCorrSummedDoubleRatioPi0FitB );
   meanErrorsSummedDoubleRatioPi0FitB = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanSummedDoubleRatioPi0FitB ,ptBinsDoubleRatioErr ,errorsMeanErrSummedDoubleRatioPi0FitB );
   meanErrorsCorrSummedDoubleRatioPi0FitB = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanCorrSummedDoubleRatioPi0FitB ,ptBinsDoubleRatioErr ,errorsMeanErrCorrSummedDoubleRatioPi0FitB );
   meanErrorsCorrSummedIncMatDoubleRatioPi0FitB = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanCorrMatSummedDoubleRatioPi0FitB ,ptBinsDoubleRatioErr ,errorsMeanErrCorrMatSummedDoubleRatioPi0FitB );

   negativeErrorsSummedDoubleRatioC = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsNegSummedDoubleRatioC ,ptBinsDoubleRatioErr ,errorsNegErrSummedDoubleRatioC );
   negativeErrorsCorrSummedDoubleRatioC = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsNegCorrSummedDoubleRatioC ,ptBinsDoubleRatioErr ,errorsNegErrCorrSummedDoubleRatioC );
   positiveErrorsSummedDoubleRatioC = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsPosSummedDoubleRatioC ,ptBinsDoubleRatioErr ,errorsPosErrSummedDoubleRatioC );
   positiveErrorsCorrSummedDoubleRatioC = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsPosCorrSummedDoubleRatioC ,ptBinsDoubleRatioErr ,errorsPosErrCorrSummedDoubleRatioC );
   meanErrorsSummedDoubleRatioC = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanSummedDoubleRatioC ,ptBinsDoubleRatioErr ,errorsMeanErrSummedDoubleRatioC );
   meanErrorsCorrSummedDoubleRatioC = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanCorrSummedDoubleRatioC ,ptBinsDoubleRatioErr ,errorsMeanErrCorrSummedDoubleRatioC );
   meanErrorsCorrSummedIncMatDoubleRatioC = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanCorrMatSummedDoubleRatioC ,ptBinsDoubleRatioErr ,errorsMeanErrCorrMatSummedDoubleRatioC );

   negativeErrorsSummedDoubleRatioPi0FitC = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsNegSummedDoubleRatioPi0FitC ,ptBinsDoubleRatioErr ,errorsNegErrSummedDoubleRatioPi0FitC );
   negativeErrorsCorrSummedDoubleRatioPi0FitC = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsNegCorrSummedDoubleRatioPi0FitC ,ptBinsDoubleRatioErr ,errorsNegErrCorrSummedDoubleRatioPi0FitC );
   positiveErrorsSummedDoubleRatioPi0FitC = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsPosSummedDoubleRatioPi0FitC ,ptBinsDoubleRatioErr ,errorsPosErrSummedDoubleRatioPi0FitC );
   positiveErrorsCorrSummedDoubleRatioPi0FitC = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsPosCorrSummedDoubleRatioPi0FitC ,ptBinsDoubleRatioErr ,errorsPosErrCorrSummedDoubleRatioPi0FitC );
   meanErrorsSummedDoubleRatioPi0FitC = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanSummedDoubleRatioPi0FitC ,ptBinsDoubleRatioErr ,errorsMeanErrSummedDoubleRatioPi0FitC );
   meanErrorsCorrSummedDoubleRatioPi0FitC = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanCorrSummedDoubleRatioPi0FitC ,ptBinsDoubleRatioErr ,errorsMeanErrCorrSummedDoubleRatioPi0FitC );
   meanErrorsCorrSummedIncMatDoubleRatioPi0FitC = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanCorrMatSummedDoubleRatioPi0FitC ,ptBinsDoubleRatioErr ,errorsMeanErrCorrMatSummedDoubleRatioPi0FitC );


   meanErrorsCorrSummedWithoutMatDoubleRatioPi0Fit = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio , errorsMeanCorrMatSummedDoubleRatioPi0FitWithoutMat,ptBinsDoubleRatioErr , errorsMeanErrCorrMatSummedDoubleRatioPi0FitWithoutMat);
   meanErrorsCorrSummedMaterialDoubleRatioPi0Fit =   new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio , errorsMaterialBudgetDoubleRatioPi0Fit,ptBinsDoubleRatioErr , errorsErrMaterialBudgetDoubleRatioPi0Fit);

   negativeErrorsSummedIncRatio = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsNegSummedIncRatio ,ptBinsIncRatioErr ,errorsNegErrSummedIncRatio );
   negativeErrorsCorrSummedIncRatio = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsNegCorrSummedIncRatio ,ptBinsIncRatioErr ,errorsNegErrCorrSummedIncRatio );
   positiveErrorsSummedIncRatio = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsPosSummedIncRatio ,ptBinsIncRatioErr ,errorsPosErrSummedIncRatio );
   positiveErrorsCorrSummedIncRatio = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsPosCorrSummedIncRatio ,ptBinsIncRatioErr ,errorsPosErrCorrSummedIncRatio );
   meanErrorsSummedIncRatio = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanSummedIncRatio ,ptBinsIncRatioErr ,errorsMeanErrSummedIncRatio );
   meanErrorsCorrSummedIncRatio = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanCorrSummedIncRatio ,ptBinsIncRatioErr ,errorsMeanErrCorrSummedIncRatio );
   meanErrorsCorrSummedIncMatIncRatio = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanCorrMatSummedIncRatio ,ptBinsIncRatioErr ,errorsMeanErrCorrMatSummedIncRatio );

   negativeErrorsSummedIncRatioPi0Fit = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsNegSummedIncRatioPi0Fit ,ptBinsIncRatioErr ,errorsNegErrSummedIncRatioPi0Fit );
   negativeErrorsCorrSummedIncRatioPi0Fit = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsNegCorrSummedIncRatioPi0Fit ,ptBinsIncRatioErr ,errorsNegErrCorrSummedIncRatioPi0Fit );
   positiveErrorsSummedIncRatioPi0Fit = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsPosSummedIncRatioPi0Fit ,ptBinsIncRatioErr ,errorsPosErrSummedIncRatioPi0Fit );
   positiveErrorsCorrSummedIncRatioPi0Fit = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsPosCorrSummedIncRatioPi0Fit ,ptBinsIncRatioErr ,errorsPosErrCorrSummedIncRatioPi0Fit );
   meanErrorsSummedIncRatioPi0Fit = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanSummedIncRatioPi0Fit ,ptBinsIncRatioErr ,errorsMeanErrSummedIncRatioPi0Fit );
   meanErrorsCorrSummedIncRatioPi0Fit = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanCorrSummedIncRatioPi0Fit ,ptBinsIncRatioErr ,errorsMeanErrCorrSummedIncRatioPi0Fit );
   meanErrorsCorrSummedIncMatIncRatioPi0Fit = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanCorrMatSummedIncRatioPi0Fit ,ptBinsIncRatioErr ,errorsMeanErrCorrMatSummedIncRatioPi0Fit );

   negativeErrorsSummedIncRatioA = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsNegSummedIncRatioA ,ptBinsIncRatioErr ,errorsNegErrSummedIncRatioA );
   negativeErrorsCorrSummedIncRatioA = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsNegCorrSummedIncRatioA ,ptBinsIncRatioErr ,errorsNegErrCorrSummedIncRatioA );
   positiveErrorsSummedIncRatioA = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsPosSummedIncRatioA ,ptBinsIncRatioErr ,errorsPosErrSummedIncRatioA );
   positiveErrorsCorrSummedIncRatioA = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsPosCorrSummedIncRatioA ,ptBinsIncRatioErr ,errorsPosErrCorrSummedIncRatioA );
   meanErrorsSummedIncRatioA = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanSummedIncRatioA ,ptBinsIncRatioErr ,errorsMeanErrSummedIncRatioA );
   meanErrorsCorrSummedIncRatioA = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanCorrSummedIncRatioA ,ptBinsIncRatioErr ,errorsMeanErrCorrSummedIncRatioA );
   meanErrorsCorrSummedIncMatIncRatioA = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanCorrMatSummedIncRatioA ,ptBinsIncRatioErr ,errorsMeanErrCorrMatSummedIncRatioA );

   negativeErrorsSummedIncRatioPi0FitA = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsNegSummedIncRatioPi0FitA ,ptBinsIncRatioErr ,errorsNegErrSummedIncRatioPi0FitA );
   negativeErrorsCorrSummedIncRatioPi0FitA = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsNegCorrSummedIncRatioPi0FitA ,ptBinsIncRatioErr ,errorsNegErrCorrSummedIncRatioPi0FitA );
   positiveErrorsSummedIncRatioPi0FitA = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsPosSummedIncRatioPi0FitA ,ptBinsIncRatioErr ,errorsPosErrSummedIncRatioPi0FitA );
   positiveErrorsCorrSummedIncRatioPi0FitA = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsPosCorrSummedIncRatioPi0FitA ,ptBinsIncRatioErr ,errorsPosErrCorrSummedIncRatioPi0FitA );
   meanErrorsSummedIncRatioPi0FitA = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanSummedIncRatioPi0FitA ,ptBinsIncRatioErr ,errorsMeanErrSummedIncRatioPi0FitA );
   meanErrorsCorrSummedIncRatioPi0FitA = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanCorrSummedIncRatioPi0FitA ,ptBinsIncRatioErr ,errorsMeanErrCorrSummedIncRatioPi0FitA );
   meanErrorsCorrSummedIncMatIncRatioPi0FitA = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanCorrMatSummedIncRatioPi0FitA ,ptBinsIncRatioErr ,errorsMeanErrCorrMatSummedIncRatioPi0FitA );

   negativeErrorsSummedIncRatioB = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsNegSummedIncRatioB ,ptBinsIncRatioErr ,errorsNegErrSummedIncRatioB );
   negativeErrorsCorrSummedIncRatioB = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsNegCorrSummedIncRatioB ,ptBinsIncRatioErr ,errorsNegErrCorrSummedIncRatioB );
   positiveErrorsSummedIncRatioB = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsPosSummedIncRatioB ,ptBinsIncRatioErr ,errorsPosErrSummedIncRatioB );
   positiveErrorsCorrSummedIncRatioB = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsPosCorrSummedIncRatioB ,ptBinsIncRatioErr ,errorsPosErrCorrSummedIncRatioB );
   meanErrorsSummedIncRatioB = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanSummedIncRatioB ,ptBinsIncRatioErr ,errorsMeanErrSummedIncRatioB );
   meanErrorsCorrSummedIncRatioB = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanCorrSummedIncRatioB ,ptBinsIncRatioErr ,errorsMeanErrCorrSummedIncRatioB );
   meanErrorsCorrSummedIncMatIncRatioB = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanCorrMatSummedIncRatioB ,ptBinsIncRatioErr ,errorsMeanErrCorrMatSummedIncRatioB );

   negativeErrorsSummedIncRatioPi0FitB = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsNegSummedIncRatioPi0FitB ,ptBinsIncRatioErr ,errorsNegErrSummedIncRatioPi0FitB );
   negativeErrorsCorrSummedIncRatioPi0FitB = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsNegCorrSummedIncRatioPi0FitB ,ptBinsIncRatioErr ,errorsNegErrCorrSummedIncRatioPi0FitB );
   positiveErrorsSummedIncRatioPi0FitB = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsPosSummedIncRatioPi0FitB ,ptBinsIncRatioErr ,errorsPosErrSummedIncRatioPi0FitB );
   positiveErrorsCorrSummedIncRatioPi0FitB = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsPosCorrSummedIncRatioPi0FitB ,ptBinsIncRatioErr ,errorsPosErrCorrSummedIncRatioPi0FitB );
   meanErrorsSummedIncRatioPi0FitB = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanSummedIncRatioPi0FitB ,ptBinsIncRatioErr ,errorsMeanErrSummedIncRatioPi0FitB );
   meanErrorsCorrSummedIncRatioPi0FitB = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanCorrSummedIncRatioPi0FitB ,ptBinsIncRatioErr ,errorsMeanErrCorrSummedIncRatioPi0FitB );
   meanErrorsCorrSummedIncMatIncRatioPi0FitB = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanCorrMatSummedIncRatioPi0FitB ,ptBinsIncRatioErr ,errorsMeanErrCorrMatSummedIncRatioPi0FitB );

   negativeErrorsSummedIncRatioC = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsNegSummedIncRatioC ,ptBinsIncRatioErr ,errorsNegErrSummedIncRatioC );
   negativeErrorsCorrSummedIncRatioC = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsNegCorrSummedIncRatioC ,ptBinsIncRatioErr ,errorsNegErrCorrSummedIncRatioC );
   positiveErrorsSummedIncRatioC = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsPosSummedIncRatioC ,ptBinsIncRatioErr ,errorsPosErrSummedIncRatioC );
   positiveErrorsCorrSummedIncRatioC = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsPosCorrSummedIncRatioC ,ptBinsIncRatioErr ,errorsPosErrCorrSummedIncRatioC );
   meanErrorsSummedIncRatioC = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanSummedIncRatioC ,ptBinsIncRatioErr ,errorsMeanErrSummedIncRatioC );
   meanErrorsCorrSummedIncRatioC = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanCorrSummedIncRatioC ,ptBinsIncRatioErr ,errorsMeanErrCorrSummedIncRatioC );
   meanErrorsCorrSummedIncMatIncRatioC = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanCorrMatSummedIncRatioC ,ptBinsIncRatioErr ,errorsMeanErrCorrMatSummedIncRatioC );

   negativeErrorsSummedIncRatioPi0FitC = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsNegSummedIncRatioPi0FitC ,ptBinsIncRatioErr ,errorsNegErrSummedIncRatioPi0FitC );
   negativeErrorsCorrSummedIncRatioPi0FitC = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsNegCorrSummedIncRatioPi0FitC ,ptBinsIncRatioErr ,errorsNegErrCorrSummedIncRatioPi0FitC );
   positiveErrorsSummedIncRatioPi0FitC = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsPosSummedIncRatioPi0FitC ,ptBinsIncRatioErr ,errorsPosErrSummedIncRatioPi0FitC );
   positiveErrorsCorrSummedIncRatioPi0FitC = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsPosCorrSummedIncRatioPi0FitC ,ptBinsIncRatioErr ,errorsPosErrCorrSummedIncRatioPi0FitC );
   meanErrorsSummedIncRatioPi0FitC = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanSummedIncRatioPi0FitC ,ptBinsIncRatioErr ,errorsMeanErrSummedIncRatioPi0FitC );
   meanErrorsCorrSummedIncRatioPi0FitC = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanCorrSummedIncRatioPi0FitC ,ptBinsIncRatioErr ,errorsMeanErrCorrSummedIncRatioPi0FitC );
   meanErrorsCorrSummedIncMatIncRatioPi0FitC = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanCorrMatSummedIncRatioPi0FitC ,ptBinsIncRatioErr ,errorsMeanErrCorrMatSummedIncRatioPi0FitC );

   negativeErrorsSummedPi0 = new TGraphErrors(ConstnBinsPi0,ptBinsPi0 ,errorsNegSummedPi0 ,ptBinsPi0Err ,errorsNegErrSummedPi0 );
   negativeErrorsCorrSummedPi0 = new TGraphErrors(ConstnBinsPi0,ptBinsPi0 ,errorsNegCorrSummedPi0 ,ptBinsPi0Err ,errorsNegErrCorrSummedPi0 );
   positiveErrorsSummedPi0 = new TGraphErrors(ConstnBinsPi0,ptBinsPi0 ,errorsPosSummedPi0 ,ptBinsPi0Err ,errorsPosErrSummedPi0 );
   positiveErrorsCorrSummedPi0 = new TGraphErrors(ConstnBinsPi0,ptBinsPi0 ,errorsPosCorrSummedPi0 ,ptBinsPi0Err ,errorsPosErrCorrSummedPi0 );
   meanErrorsSummedPi0 = new TGraphErrors(ConstnBinsPi0,ptBinsPi0 ,errorsMeanSummedPi0 ,ptBinsPi0Err ,errorsMeanErrSummedPi0 );
   meanErrorsCorrSummedPi0 = new TGraphErrors(ConstnBinsPi0,ptBinsPi0 ,errorsMeanCorrSummedPi0 ,ptBinsPi0Err ,errorsMeanErrCorrSummedPi0 );
   meanErrorsCorrSummedIncMatPi0 = new TGraphErrors(ConstnBinsPi0,ptBinsPi0 ,errorsMeanCorrMatSummedPi0 ,ptBinsPi0Err ,errorsMeanErrCorrMatSummedPi0 );
   meanErrorsCorrSummedIncMatPi0WithoutMat = new TGraphErrors(ConstnBinsPi0,ptBinsPi0 ,errorsMeanCorrMatSummedPi0WithoutMaterial ,ptBinsPi0Err ,errorsMeanErrCorrMatSummedPi0 );

   negativeErrorsSummedPi0Fit = new TGraphErrors(ConstnBinsPi0,ptBinsPi0 ,errorsNegSummedPi0Fit ,ptBinsPi0Err ,errorsNegErrSummedPi0Fit );
   negativeErrorsCorrSummedPi0Fit = new TGraphErrors(ConstnBinsPi0,ptBinsPi0 ,errorsNegCorrSummedPi0Fit ,ptBinsPi0Err ,errorsNegErrCorrSummedPi0Fit );
   positiveErrorsSummedPi0Fit = new TGraphErrors(ConstnBinsPi0,ptBinsPi0 ,errorsPosSummedPi0Fit ,ptBinsPi0Err ,errorsPosErrSummedPi0Fit );
   positiveErrorsCorrSummedPi0Fit = new TGraphErrors(ConstnBinsPi0,ptBinsPi0 ,errorsPosCorrSummedPi0Fit ,ptBinsPi0Err ,errorsPosErrCorrSummedPi0Fit );
   meanErrorsSummedPi0Fit = new TGraphErrors(ConstnBinsPi0,ptBinsPi0 ,errorsMeanSummedPi0Fit ,ptBinsPi0Err ,errorsMeanErrSummedPi0Fit );
   meanErrorsCorrSummedPi0Fit = new TGraphErrors(ConstnBinsPi0,ptBinsPi0 ,errorsMeanCorrSummedPi0Fit ,ptBinsPi0Err ,errorsMeanErrCorrSummedPi0Fit );
   meanErrorsCorrSummedIncMatPi0Fit = new TGraphErrors(ConstnBinsPi0,ptBinsPi0 ,errorsMeanCorrMatSummedPi0Fit ,ptBinsPi0Err ,errorsMeanErrCorrMatSummedPi0Fit );
   meanErrorsCorrSummedIncMatPi0FitWithoutMat = new TGraphErrors(ConstnBinsPi0,ptBinsPi0 ,errorsMeanCorrMatSummedPi0FitWithoutMaterial ,ptBinsPi0Err ,errorsMeanErrCorrMatSummedPi0Fit );

   negativeErrorsSummedGamma = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsNegSummedGamma ,ptBinsGammaErr ,errorsNegErrSummedGamma );
   negativeErrorsCorrSummedGamma = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsNegCorrSummedGamma ,ptBinsGammaErr ,errorsNegErrCorrSummedGamma );
   positiveErrorsSummedGamma = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsPosSummedGamma ,ptBinsGammaErr ,errorsPosErrSummedGamma );
   positiveErrorsCorrSummedGamma = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsPosCorrSummedGamma ,ptBinsGammaErr ,errorsPosErrCorrSummedGamma );
   meanErrorsSummedGamma = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsMeanSummedGamma ,ptBinsGammaErr ,errorsMeanErrSummedGamma );
   meanErrorsCorrSummedGamma = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsMeanCorrSummedGamma ,ptBinsGammaErr ,errorsMeanErrCorrSummedGamma );
   meanErrorsCorrSummedIncMatGamma = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsMeanCorrMatSummedGamma ,ptBinsGammaErr ,errorsMeanErrCorrMatSummedGamma );

   meanErrorsCorrSummedWithoutMatGamma = new TGraphErrors(ConstnBinsGamma,ptBinsGamma , errorsMeanCorrMatSummedGammaWithoutMat,ptBinsGammaErr , errorsMeanErrCorrMatSummedGammaWithoutMat);
   meanErrorsCorrSummedMaterialGamma = new TGraphErrors(ConstnBinsGamma,ptBinsGamma , errorsMaterialBudget,ptBinsGammaErr , errorsErrMaterialBudget);

   negativeErrorsSummedGammaA = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsNegSummedGammaA ,ptBinsGammaErr ,errorsNegErrSummedGammaA );
   negativeErrorsCorrSummedGammaA = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsNegCorrSummedGammaA ,ptBinsGammaErr ,errorsNegErrCorrSummedGammaA );
   positiveErrorsSummedGammaA = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsPosSummedGammaA ,ptBinsGammaErr ,errorsPosErrSummedGammaA );
   positiveErrorsCorrSummedGammaA = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsPosCorrSummedGammaA ,ptBinsGammaErr ,errorsPosErrCorrSummedGammaA );
   meanErrorsSummedGammaA = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsMeanSummedGammaA ,ptBinsGammaErr ,errorsMeanErrSummedGammaA );
   meanErrorsCorrSummedGammaA = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsMeanCorrSummedGammaA ,ptBinsGammaErr ,errorsMeanErrCorrSummedGammaA );
   meanErrorsCorrSummedIncMatGammaA = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsMeanCorrMatSummedGammaA ,ptBinsGammaErr ,errorsMeanErrCorrMatSummedGammaA );

   negativeErrorsSummedGammaB = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsNegSummedGammaB ,ptBinsGammaErr ,errorsNegErrSummedGammaB );
   negativeErrorsCorrSummedGammaB = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsNegCorrSummedGammaB ,ptBinsGammaErr ,errorsNegErrCorrSummedGammaB );
   positiveErrorsSummedGammaB = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsPosSummedGammaB ,ptBinsGammaErr ,errorsPosErrSummedGammaB );
   positiveErrorsCorrSummedGammaB = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsPosCorrSummedGammaB ,ptBinsGammaErr ,errorsPosErrCorrSummedGammaB );
   meanErrorsSummedGammaB = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsMeanSummedGammaB ,ptBinsGammaErr ,errorsMeanErrSummedGammaB );
   meanErrorsCorrSummedGammaB = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsMeanCorrSummedGammaB ,ptBinsGammaErr ,errorsMeanErrCorrSummedGammaB );
   meanErrorsCorrSummedIncMatGammaB = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsMeanCorrMatSummedGammaB ,ptBinsGammaErr ,errorsMeanErrCorrMatSummedGammaB );

   negativeErrorsSummedGammaC = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsNegSummedGammaC ,ptBinsGammaErr ,errorsNegErrSummedGammaC );
   negativeErrorsCorrSummedGammaC = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsNegCorrSummedGammaC ,ptBinsGammaErr ,errorsNegErrCorrSummedGammaC );
   positiveErrorsSummedGammaC = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsPosSummedGammaC ,ptBinsGammaErr ,errorsPosErrSummedGammaC );
   positiveErrorsCorrSummedGammaC = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsPosCorrSummedGammaC ,ptBinsGammaErr ,errorsPosErrCorrSummedGammaC );
   meanErrorsSummedGammaC = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsMeanSummedGammaC ,ptBinsGammaErr ,errorsMeanErrSummedGammaC );
   meanErrorsCorrSummedGammaC = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsMeanCorrSummedGammaC ,ptBinsGammaErr ,errorsMeanErrCorrSummedGammaC );
   meanErrorsCorrSummedIncMatGammaC = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsMeanCorrMatSummedGammaC ,ptBinsGammaErr ,errorsMeanErrCorrMatSummedGammaC );



   // Axis Range
   // 0-40%
   Double_t gammaSpec[2];
   Double_t gammaErrors[2];
   Double_t doubleRatioErrors[2];
   Double_t doubleRatio[2];
   Double_t incRatio[2];
   Double_t doubleRatioX[2];
   Double_t directPhoton[2];
   TString cent ="";


   Color_t	colorComb0005				= kOrange+10;
   Color_t	colorComb0510				= kOrange+7;
   Color_t	colorComb0010				= kRed-7;
   Color_t	colorComb0020				= kRed+1;
   Color_t	colorComb0040				= kMagenta+2;
   Color_t	colorComb1020				= 800;
   Color_t	colorComb2040				= kGreen+2;
   Color_t	colorComb4060				= kCyan+2;
   Color_t	colorComb6080				= kBlue+1;

   Color_t colorcent;
   TString centralityCutNumber;
   if(option.CompareTo("2.76TeV") == 0){
      gammaSpec[0] = 1e-5; gammaSpec[1] = 10;
      gammaErrors[0] = 0; gammaErrors[1] = 35;
      doubleRatioErrors[0] = 0.0; doubleRatioErrors[1] = 35;
      doubleRatio[0] = 0.55; doubleRatio[1] = 3.0;
      incRatio[0] = 0.0; incRatio[1] = 1.7;
      doubleRatioX[0] = 0.8; doubleRatioX[1] = 7.2;
      directPhoton[0] = 1e-5; directPhoton[1] = 10;
      cent = "2760GeV";
   }
   if(option.CompareTo("PbPb_2.76TeV") == 0){
      centralityCutNumber = cutSel(GetEventSystemCutPosition(),3);
      if(centralityCutNumber.CompareTo("502") == 0 ||
         centralityCutNumber.CompareTo("501") == 0 ||
         centralityCutNumber.CompareTo("512") == 0 ||
         centralityCutNumber.CompareTo("601") == 0 ||
         centralityCutNumber.CompareTo("612") == 0 ||
         centralityCutNumber.CompareTo("524") == 0 ||
         centralityCutNumber.CompareTo("548") == 0 ||
         centralityCutNumber.CompareTo("504") == 0 ||
         centralityCutNumber.CompareTo("508") == 0){
         gammaSpec[0] = 1e-5; gammaSpec[1] = 1000;
         gammaErrors[0] = 0; gammaErrors[1] = 35;
         doubleRatioErrors[0] = 0.0; doubleRatioErrors[1] = 32;
         doubleRatio[0] = 0.75; doubleRatio[1] = 2.0;
         incRatio[0] = 0.0; incRatio[1] = 1.7;
         doubleRatioX[0] = 0.8; doubleRatioX[1] = 8;
         directPhoton[0] = 1e-5; directPhoton[1] = 8;
         cent = centralityCutNumber;
      }
   }
   if(centralityCutNumber.CompareTo("504") == 0)  colorcent=colorComb0040;
   if(centralityCutNumber.CompareTo("502") == 0)  colorcent=colorComb0020;
   if(centralityCutNumber.CompareTo("501") == 0)  colorcent=colorComb0010;
   if(centralityCutNumber.CompareTo("612") == 0)  colorcent=colorComb0510;
   if(centralityCutNumber.CompareTo("601") == 0)  colorcent=colorComb0005;
   if(centralityCutNumber.CompareTo("512") == 0)  colorcent=colorComb1020;
   if(centralityCutNumber.CompareTo("524") == 0)  colorcent=colorComb2040;
   if(centralityCutNumber.CompareTo("548") == 0)  colorcent=colorComb4060;


   // Draw Final Gamma Pictures
   TF1 *One = new TF1("One","1",0,16);
   One->SetLineWidth(1.2);
   One->SetLineColor(kGray+2);
   //Int_t color[20] = {860,894,620,880,591,403,802,923,634,432,422,435,624,407,2,830,404,608,920,1};
   Int_t color[13] = {kBlue-6,kPink-2,kRed+2,kOrange+7,kCyan-3,kMagenta+2,kRed-9,kAzure-1,kAzure+5,kGreen-3,kSpring+9,kMagenta-5,kOrange-8};

   TH1D *dummy = new TH1D("dummy","",160,0.0,16.0);//,1500,-10,10);
   dummy->GetXaxis()->SetLabelOffset(-0.015);
   //dummy->GetXaxis()->SetNdivisions(510);

   Double_t *xErrorDummy = new Double_t[1];
   xErrorDummy[0] = 5;
   TGraphErrors *GraphErrorsDummy = new TGraphErrors(1,0,0,xErrorDummy,0);
   GraphErrorsDummy->SetMarkerSize(2);
   GraphErrorsDummy->SetMarkerStyle(20);
   GraphErrorsDummy->SetMarkerColor(colorcent);
   GraphErrorsDummy->SetLineColor(1);
   GraphErrorsDummy->SetFillColor(kGray);

   TGraphErrors *GraphErrorsDummyB = new TGraphErrors(1,0,0,xErrorDummy,0);
   GraphErrorsDummyB->SetMarkerSize(2);
   GraphErrorsDummyB->SetMarkerStyle(20);
   GraphErrorsDummyB->SetMarkerColor(colorcent);
   GraphErrorsDummyB->SetLineColor(1);
   GraphErrorsDummyB->SetFillColor(kGray);


   TCanvas *cocktailCanvasSpec = GetAndSetCanvas("cocktailCanvasSpec");
   cocktailCanvasSpec->SetLogy();
   cocktailAllGamma->SetTitle("");
   cocktailAllGamma ->GetYaxis()->SetTitle("#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c}^{2})");
   cocktailAllGamma->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
   TLatex *gammaL= new TLatex(2.0,5.0,"#gamma from ALICE Data");
   TLatex *cocktail= new TLatex(9.0,1000,"all decay #gamma");
   TLatex *tpi = new TLatex(9.0,100.,"#pi^{0} #rightarrow #gamma#gamma (e^{+}e^{-}#gamma)");
   TLatex *teta = new TLatex(9.0,10.,"#eta #rightarrow #gamma#gamma (#pi^{+}#pi^{-}#gamma,e^{+}e^{-}#gamma,#pi^{0}#gamma#gamma)");
   TLatex *tomega = new TLatex(9.0,0.8,"#omega #rightarrow #pi^{0}#gamma (#eta#gamma)");
   TLatex *tetaprime = new TLatex(9.0,0.09,"#eta' #rightarrow #rho#gamma (#omega#gamma, #gamma#gamma)");
   TLatex *tsigma = new TLatex(9.0,0.09,"#eta' #rightarrow #rho#gamma (#omega#gamma, #gamma#gamma)");
   TLatex *tphi = new TLatex(9.0,0.007,"#phi #rightarrow #eta#gamma (#pi^{0}#gamma, #omega#gamma)");
   TLatex *trho = new TLatex(9.0,0.0007,"#rho #rightarrow #pi^{+}#pi^{-}#gamma (#pi^{0}#gamma, #eta#gamma)");

   cocktailAllGamma->GetYaxis()->SetRangeUser(1e-7,1e3);
   cocktailAllGamma->GetYaxis()->SetTitleSize(0.035);
   cocktailAllGamma->GetYaxis()->SetTitleOffset(1.3);
   cocktailAllGamma->GetXaxis()->SetRangeUser(0.4,5);
   cocktailAllGamma->SetMarkerSize(0);
   cocktailAllGamma->SetLineColor(kBlack);
   cocktailAllGamma->SetLineWidth(2);
   cocktailAllGamma->Draw("chist9");
   cocktailPi0Gamma->Draw("csamehist9");
   cocktailPi0Gamma->SetLineColor(kRed);
   cocktailPi0Gamma->SetLineWidth(2);
   cocktailEtaGamma->Draw("csamehist9");
   cocktailEtaGamma->SetLineColor(kBlue);
   cocktailEtaGamma->SetLineWidth(2);
   cocktailEtapGamma->Draw("csamehist");
   cocktailEtapGamma->SetLineColor(kGray+2);
   cocktailEtapGamma->SetLineWidth(2);
   cocktailOmegaGamma->Draw("csamehist");
   cocktailOmegaGamma->SetLineColor(kTeal+3);
   cocktailOmegaGamma->SetLineWidth(2);
   cocktailPhiGamma->Draw("csamehist");
   cocktailPhiGamma->SetLineColor(kViolet-1);
   cocktailPhiGamma->SetLineWidth(2);
   cocktailRhoGamma->Draw("csamehist");
   cocktailRhoGamma->SetLineColor(kOrange+4);
   cocktailRhoGamma->SetLineWidth(2);

   cocktail->SetTextColor(kBlack);
   cocktail->SetTextSize(0.04);
   cocktail->SetTextFont(42);
   cocktail->Draw();
   tpi->SetTextColor(kRed);
   tpi->SetTextSize(0.04);
   tpi->SetTextFont(42);
   tpi->Draw();
   teta->SetTextColor(kBlue);
   teta->SetTextSize(0.04);
   teta->SetTextFont(42);
   teta->Draw();
   tomega->SetTextColor(kTeal+3);
   tomega->SetTextSize(0.04);
   tomega->SetTextFont(42);
   tomega->Draw();
   tetaprime->SetTextColor(kGray+2);
   tetaprime->SetTextSize(0.04);
   tetaprime->SetTextFont(42);
   tetaprime->Draw();
   tphi->SetTextColor(kViolet-1);
   tphi->SetTextSize(0.04);
   tphi->SetTextFont(42);
   tphi->Draw();
   trho->SetTextColor(kOrange+4);
   trho->SetTextSize(0.04);
   trho->SetTextFont(42);
   trho->Draw();


   TLegend* legendCocktailGammaSpec = GetAndSetLegend(0.53,0.67,6.2,2,"",0.5);
   legendCocktailGammaSpec->SetTextSize(0.05);
   legendCocktailGammaSpec->AddEntry(cocktailAllGamma,"all #gamma_{decay}","l");
   legendCocktailGammaSpec->AddEntry(cocktailPi0Gamma,"#gamma_{decay}^{#pi^{0}}","l");
   legendCocktailGammaSpec->AddEntry(cocktailEtaGamma,"#gamma_{decay}^{#eta}","l");
   legendCocktailGammaSpec->AddEntry(cocktailOmegaGamma,"#gamma_{decay}^{#omega}","l");
   legendCocktailGammaSpec->AddEntry(cocktailEtapGamma,"#gamma_{decay}^{#eta'}","l");
   legendCocktailGammaSpec->AddEntry(cocktailPhiGamma,"#gamma_{decay}^{#phi}","l");
   legendCocktailGammaSpec->Draw();


   DrawSystem(0.27,0.16,option.Data(),(GetCentralityStringC(cutSel)).Data());
   cocktailCanvasSpec->Print(Form("%s/CocktailDecaySpectra_%s.eps",outputDir.Data(),cent.Data()));
   cocktailCanvasSpec->Print(Form("%s/CocktailDecaySpectra_%s.C",outputDir.Data(),cent.Data()));




   TCanvas *cocktailCanvasMesonSpec = GetAndSetCanvas("cocktailCanvasMesonSpec");
   cocktailCanvasMesonSpec->SetLogy();
   cocktailAllGamma->SetTitle("");
   cocktailAllGamma ->GetYaxis()->SetTitle("#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c}^{2})");
   cocktailAllGamma->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
   TLatex *gammaL= new TLatex(2.0,5.0,"#gamma from ALICE Data");
   TLatex *cocktail= new TLatex(9.0,1000,"all decay #gamma");
   TLatex *tpi = new TLatex(9.0,100.,"#pi^{0}");
   TLatex *teta = new TLatex(9.0,10.,"#eta");
   TLatex *tomega = new TLatex(9.0,0.8,"#omega");
   TLatex *tetaprime = new TLatex(9.0,0.09,"#eta'");
   TLatex *tsigma = new TLatex(9.0,0.09,"#Sigma^{0}");
   TLatex *tphi = new TLatex(9.0,0.007,"#phi");
   TLatex *tsigma = new TLatex(9.0,0.0007,"#Sigma^{0}");
   TLatex *trho = new TLatex(9.0,0.0007,"#rho #rightarrow #pi^{+}#pi^{-}#gamma (#pi^{0}#gamma, #eta#gamma)");


   cocktailPion->SetTitle("");
   cocktailPion->GetYaxis()->SetRangeUser(5e-4,1e3);
   cocktailPion->GetYaxis()->SetTitleSize(0.035);
   cocktailPion->GetYaxis()->SetTitleOffset(1.3);
   cocktailPion->GetXaxis()->SetRangeUser(0.4,5);
   cocktailPion->SetMarkerSize(0);
   cocktailPion->SetLineColor(kBlack);
   cocktailPion->SetLineWidth(2);
   cocktailPion->GetYaxis()->SetTitle("#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c}^{2})");
   cocktailPion->GetXaxis()->SetTitle("p_{T} (Gev/c)");
   cocktailPion->Draw("chist9");
   cocktailPion->Draw("csamehist9");
   cocktailPion->SetLineColor(kRed);
   cocktailPion->SetLineWidth(2);
   cocktailEta->DrawCopy("csamehist9");
   for(Int_t i = 0;i<cocktailEta->GetNbinsX();i++){
      cocktailEtaGamma->SetBinError(i+1,cocktailEtaGamma->GetBinContent(i+1)*sqrt( pow(0.2,2) + pow(cocktailEtaLow->GetBinContent(i+1),2) ));
      cocktailEta->SetBinError(i+1,cocktailEta->GetBinContent(i+1)*cocktailEtaGamma->GetBinError(i+1)/cocktailEtaGamma->GetBinContent(i+1));
   }
   cocktailEta->SetFillStyle(3003);
   cocktailEta->SetLineColor(kBlue);
   cocktailEta->SetMarkerSize(0);
   cocktailEta->SetFillColor(kBlue);
   cocktailEta->DrawCopy("sameE4");

   cocktailEta->SetLineColor(kBlue);
   cocktailEta->SetLineWidth(2);
   cocktailEtap->Draw("csamehist");
   cocktailEtap->SetLineColor(kGray+2);
   cocktailEtap->SetLineWidth(2);
   cocktailOmega->Draw("csamehist");
   cocktailOmega->SetLineColor(kTeal+3);
   cocktailOmega->SetLineWidth(2);
   cocktailPhi->Draw("csamehist");
   cocktailPhi->SetLineColor(kViolet-1);
   cocktailPhi->SetLineWidth(2);
   cocktailSigma->Draw("csamehist");
   cocktailSigma->SetLineColor(kRed-2);
   cocktailSigma->SetLineWidth(2);

   // cocktailRho->Draw("csamehist");
   // cocktailRho->SetLineColor(kOrange+4);
   // cocktailRho->SetLineWidth(2);

   Pi0->SetMarkerSize(2);
   Pi0->Draw("E1same");
   K0s->SetMarkerSize(2);
   K0s->Draw("E1same");

   tpi->SetTextColor(kRed);
   tpi->SetTextSize(0.04);
   tpi->SetTextFont(42);
   tpi->Draw();
   teta->SetTextColor(kBlue);
   teta->SetTextSize(0.04);
   teta->SetTextFont(42);
   teta->Draw();
   tomega->SetTextColor(kTeal+3);
   tomega->SetTextSize(0.04);
   tomega->SetTextFont(42);
   tomega->Draw();
   tetaprime->SetTextColor(kGray+2);
   tetaprime->SetTextSize(0.04);
   tetaprime->SetTextFont(42);
   tetaprime->Draw();
   tphi->SetTextColor(kViolet-1);
   tphi->SetTextSize(0.04);
   tphi->SetTextFont(42);
   tphi->Draw();
   // trho->SetTextColor(kOrange+4);
   // trho->SetTextSize(0.04);
   // trho->SetTextFont(42);
   // trho->Draw();
   tsigma->SetTextColor(kRed-2);
   tsigma->SetTextSize(0.04);
   tsigma->Draw();

   TLegend* legendCocktailMesonSpec = GetAndSetLegend(0.43,0.67,6.2,2,"",1.);
   legendCocktailMesonSpec->SetTextSize(0.05);
   legendCocktailMesonSpec->AddEntry(cocktailPi0,"#pi^{0}","l");
   legendCocktailMesonSpec->AddEntry(Pi0,"#pi^{0}","p");
   legendCocktailMesonSpec->AddEntry(cocktailEta,"#eta","l");
   legendCocktailMesonSpec->AddEntry(K0s,"scaled K_{s}^{0} x 0.799","p");
   legendCocktailMesonSpec->AddEntry(cocktailOmega,"#omega","l");
   legendCocktailMesonSpec->AddEntry(cocktailEtap,"#eta'","l");
   legendCocktailMesonSpec->AddEntry(cocktailPhi,"#phi","l");
   legendCocktailMesonSpec->AddEntry(cocktailSigma,"#Sigma^{0}","l");
   legendCocktailMesonSpec->Draw();



   DrawSystem(0.27,0.16,option.Data(),(GetCentralityStringC(cutSel)).Data());
   cocktailCanvasMesonSpec->Print(Form("%s/CocktailMesonSpectra_%s.eps",outputDir.Data(),cent.Data()));
   cocktailCanvasMesonSpec->Print(Form("%s/CocktailMesonSpectra_%s.C",outputDir.Data(),cent.Data()));





   TCanvas *cocktailCanvasRatio = GetAndSetCanvas("cocktailCanvasRatio");
   cocktailCanvasRatio->SetLogy();
   cocktailAllGammaPi0->SetTitle("");
   cocktailAllGammaPi0->GetXaxis()->SetRangeUser(0.0,14);
   cocktailAllGammaPi0->GetYaxis()->SetRangeUser(1e-5,5000);
   cocktailAllGammaPi0->SetMarkerSize(0);

   cocktailAllGammaPi0->SetLineColor(kBlack);
   cocktailAllGammaPi0->SetLineWidth(2);
   cocktailAllGammaPi0->GetYaxis()->SetTitle("#gamma_{decay} / #pi^{0}");
   cocktailAllGammaPi0->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
   cocktailAllGammaPi0->Draw("chist");
   cocktailPi0GammaPi0->Draw("csamehist");
   cocktailPi0GammaPi0->SetLineColor(kRed);
   cocktailPi0GammaPi0->SetLineWidth(2);
   cocktailEtaGammaPi0->Draw("csamehist");
   cocktailEtaGammaPi0->SetLineColor(kBlue);
   cocktailEtaGammaPi0->SetLineWidth(2);
   cocktailEtapGammaPi0->Draw("csamehist");
   cocktailEtapGammaPi0->SetLineColor(kGray+2);
   cocktailEtapGammaPi0->SetLineWidth(2);
   cocktailOmegaGammaPi0->Draw("csamehist");
   cocktailOmegaGammaPi0->SetLineColor(kTeal+3);
   cocktailOmegaGammaPi0->SetLineWidth(2);
   cocktailPhiGammaPi0->Draw("csamehist");
   cocktailPhiGammaPi0->SetLineColor(kViolet-1);
   cocktailPhiGammaPi0->SetLineWidth(2);
   cocktailRhoGammaPi0->Draw("csamehist");
   cocktailRhoGammaPi0->SetLineColor(kOrange+4);
   cocktailRhoGammaPi0->SetLineWidth(2);


   cocktail= new TLatex(9.0,1200.,"all decay #gamma");
   gammaL= new TLatex(4.7,0.1,"#gamma_{inc}/#pi^{0} from ALICE Data");
   tpi = new TLatex(9.0,400.,"#pi^{0} #rightarrow #gamma#gamma (e^{+}e^{-}#gamma)");
   teta = new TLatex(9.0,140.,"#eta #rightarrow #gamma#gamma (#pi^{+}#pi^{-}#gamma,e^{+}e^{-}#gamma,#pi^{0}#gamma#gamma)");
   tomega = new TLatex(9.0,40.,"#omega #rightarrow #pi^{0}#gamma (#eta#gamma)");
   tetaprime = new TLatex(9.0,12,"#eta' #rightarrow #rho#gamma (#omega#gamma, #gamma#gamma)");
   tphi = new TLatex(9.0,4.,"#phi #rightarrow #eta#gamma (#pi^{0}#gamma, #omega#gamma)");
   trho = new TLatex(9.0,1,"#rho #rightarrow #pi^{+}#pi^{-}#gamma (#pi^{0}#gamma, #eta#gamma)");


   cocktail->SetTextColor(kBlack);
   cocktail->SetTextSize(0.04);
   cocktail->SetTextFont(42);
   cocktail->Draw();
   tpi->SetTextColor(kRed);
   tpi->SetTextSize(0.04);
   tpi->SetTextFont(42);
   tpi->Draw();
   teta->SetTextColor(kBlue);
   teta->SetTextSize(0.04);
   teta->SetTextFont(42);
   teta->Draw();
   tomega->SetTextColor(kTeal+3);
   tomega->SetTextSize(0.04);
   tomega->SetTextFont(42);
   tomega->Draw();
   tetaprime->SetTextColor(kGray+2);
   tetaprime->SetTextSize(0.04);
   tetaprime->SetTextFont(42);
   tetaprime->Draw();
   tphi->SetTextColor(kViolet-1);
   tphi->SetTextSize(0.04);
   tphi->SetTextFont(42);
   tphi->Draw();
   trho->SetTextColor(kOrange+4);
   trho->SetTextSize(0.04);
   trho->SetTextFont(42);
   trho->Draw();


   DrawSystem(0.24,0.14,option.Data(),(GetCentralityStringC(cutSel)).Data());


   cocktailCanvasRatio->Print(Form("%s/CocktailDecayRatios_%s.eps",outputDir.Data(),cent.Data()));
   cocktailCanvasRatio->Print(Form("%s/CocktailDecayRatios_%s.C",outputDir.Data(),cent.Data()));







   TCanvas *cocktailCanvasGammaRatio = GetAndSetCanvas("cocktailCanvasGammaRatio");

   cocktailCanvasGammaRatio->SetLogy(1);
   cocktailGammaPionsAllGamma->SetTitle("");
   cocktailGammaPionsAllGamma->GetXaxis()->SetRangeUser(0.4,5);
   cocktailGammaPionsAllGamma->GetYaxis()->SetRangeUser(0.00001,0.25);
   //cocktailGammaPionsAllGamma->GetYaxis()->SetRangeUser(0.00001,0.5);
   cocktailGammaPionsAllGamma->SetMarkerSize(0);


   cocktailGammaPionsAllGamma->SetLineColor(kRed);
   cocktailGammaPionsAllGamma->SetLineWidth(2);
   cocktailGammaPionsAllGamma->GetYaxis()->SetTitle("#gamma_{decay}^{hadron} / #gamma_{decay}^{#pi^{0}}");
   cocktailGammaPionsAllGamma->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
   cocktailGammaPionsAllGamma->Draw("chist");
   cocktailGammaEtaAllGamma->SetLineColor(kBlue);
   cocktailGammaEtaAllGamma->SetLineWidth(2);
   cocktailGammaEtaAllGamma->DrawCopy("csamehist");
   

   for(Int_t i = 0;i<cocktailGammaEtaAllGamma->GetNbinsX();i++){
      cocktailEtaGamma->SetBinError(i+1,cocktailEtaGamma->GetBinContent(i+1)*sqrt( pow(0.2,2) + pow(cocktailEtaGammaLow->GetBinContent(i+1),2) ));
      cocktailGammaEtaAllGamma->SetBinError(i+1,cocktailGammaEtaAllGamma->GetBinContent(i+1)*cocktailEtaGamma->GetBinError(i+1)/cocktailEtaGamma->GetBinContent(i+1));
   }

   cocktailGammaEtaAllGamma->SetFillStyle(3003);
   cocktailGammaEtaAllGamma->SetLineColor(kBlue);
   cocktailGammaEtaAllGamma->SetMarkerSize(0);
   cocktailGammaEtaAllGamma->SetFillColor(kBlue);
   cocktailGammaEtaAllGamma->DrawCopy("sameE4");
   cocktailGammaOmegaAllGamma->SetLineColor(kTeal+3);
   cocktailGammaOmegaAllGamma->SetLineWidth(2);
   cocktailGammaOmegaAllGamma->Draw("csamehist");
   cocktailGammaSigmaAllGamma->SetLineColor(kRed-2);
   cocktailGammaSigmaAllGamma->SetLineWidth(2);
   cocktailGammaSigmaAllGamma->Draw("csamehist");
   cocktailGammaEtapAllGamma->SetLineColor(kGray+2);
   cocktailGammaEtapAllGamma->SetLineWidth(2);
   cocktailGammaEtapAllGamma->Draw("csamehist");
   cocktailGammaPhiAllGamma->SetLineColor(kViolet-1);
   cocktailGammaPhiAllGamma->SetLineWidth(2);
   cocktailGammaPhiAllGamma->Draw("csamehist");

   cocktail= new TLatex(9.0,1.,"all decay #gamma");
   // gammaL= new TLatex(4.7,0.1,"#gamma_{inc}/#pi^{0} from ALICE Data");
   tpi = new TLatex(9.0,.7,"#gamma_{decay}^{sum}/#gamma_{decay}^{#pi^{0}}");
   teta = new TLatex(9.0,.1,"#gamma_{decay}^{#eta}/#gamma_{decay}^{#pi^{0}}");
   tomega = new TLatex(9.0,0.02,"#gamma_{decay}^{#omega}/#gamma_{decay}^{#pi^{0}}");
   tetaprime = new TLatex(9.0,0.004,"#gamma_{decay}^{#eta'}/#gamma_{decay}^{#pi^{0}}");
   tphi = new TLatex(9.0,0.0005,"#gamma_{decay}^{#phi}/#gamma_{decay}^{#pi^{0}}");
   trho = new TLatex(9.0,1,"#gamma_{decay}^{#pi^{0}}/#gamma_{decay}^{#pi^{0}}");
   tsigma = new TLatex(9.0,0.00005,"#gamma_{decay}^{#Sigma^{0}}/#gamma_{decay}^{#pi^{0}}");
   
   // //cocktail->Draw();
   // tpi->SetTextColor(kRed);
   // tpi->SetTextSize(0.04);
   // tpi->SetTextFont(42);
   // tpi->Draw();
   // teta->SetTextColor(kBlue);
   // teta->SetTextSize(0.04);
   // teta->SetTextFont(42);
   // teta->Draw();
   // tomega->SetTextColor(kTeal+3);
   // tomega->SetTextSize(0.04);
   // tomega->SetTextFont(42);
   // tomega->Draw();
   // tphi->SetTextColor(kViolet-1);
   // tphi->SetTextSize(0.04);
   // tphi->SetTextFont(42);
   // tphi->Draw();
   // tetaprime->SetTextColor(kGray+2);
   // tetaprime->SetTextSize(0.04);
   // tetaprime->Draw();
   // tsigma->SetTextColor(kRed-2);
   // tsigma->SetTextSize(0.04);
   // tsigma->Draw();
   
   //TString ttt = "measured #it{p}_{T} > 0.9 GeV/c   "+GetCentralityStringC(cutSel);

   //DrawSystem(0.5,0.3,option.Data(),GetCentralityStringC(cutSel));

   TLegend* legendCocktailGammaRatios = GetAndSetLegend(0.13,0.10,6.2,5,"",.5);
   legendCocktailGammaRatios->SetTextSize(0.05);
   legendCocktailGammaRatios->AddEntry(cocktailEta,"#gamma_{#eta}/#gamma_{#pi^{0}}","l");
   legendCocktailGammaRatios->AddEntry(cocktailOmega,"#gamma_{#omega}/#gamma_{#pi^{0}}","l");
   legendCocktailGammaRatios->AddEntry(cocktailEtap,"#gamma_{#eta'}/#gamma_{#pi^{0}}","l");
   legendCocktailGammaRatios->AddEntry(cocktailGammaPhiAllGamma,"#gamma_{#phi}/#gamma_{#pi^{0}}","l");
   legendCocktailGammaRatios->AddEntry(cocktailGammaSigmaAllGamma,"#gamma_{#Sigma^{0}}/#gamma_{#pi^{0}}","l");
   legendCocktailGammaRatios->Draw();


   cocktailCanvasGammaRatio->Print(Form("%s/CocktailDecayGammaRatios_%s.eps",outputDir.Data(),cent.Data()));
   cocktailCanvasGammaRatio->Print(Form("%s/CocktailDecayGammaRatios_%s.C",outputDir.Data(),cent.Data()));



   // Shift Uncert
   //Double_t fBinsPi0HIDirectPhotonPt[22] = {0.0,0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,3.0,3.3,3.7,4.1,4.6,5.4,6.2,7.0,8.0,11.0,14.0};
   //Double_t fBinsPi0HIPt[22] =             {0.0,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,3.0,3.5,4.0,5.0,6.0,8.0,10.0,12.0,14.0};



   TGraphErrors *meanPidEdxTypeGammaB;
   TGraphErrors *meanQtTypeGammaB;
   TGraphErrors *meanChi2TypeGammaB;
   TGraphErrors *meanedEdxTypeGammaB;

   Bool_t createdGammaCuts = kFALSE;
   TGraphErrors *meanErrorsGammaGammaCuts;
   Bool_t createdPID = kFALSE;
   TGraphErrors *meanErrorsGammaPID;
   Bool_t createdTrackCuts = kFALSE;
   TGraphErrors *meanErrorsGammaTrackCuts;

   for(Int_t i = 0; i< ErrorGammaSpec; i++){
      DrawGammaSetMarkerTGraphErr(meanErrorsGamma[i], 20+i, 1.,color[i],color[i]);
      if(!nameCutVariationsDoubleRatio[i].CompareTo("PidEdx"))
         meanPidEdxTypeGammaB = (TGraphErrors*) meanErrorsGamma[i]->Clone("GammaError_Syst_TypeB_PidEdx");
      if(!nameCutVariationsDoubleRatio[i].CompareTo("qT"))
         meanQtTypeGammaB = (TGraphErrors*)meanErrorsGamma[i]->Clone("GammaError_Syst_TypeB_qT");
      if(!nameCutVariationsDoubleRatio[i].CompareTo("Chi2"))
         meanChi2TypeGammaB = (TGraphErrors*)meanErrorsGamma[i]->Clone("GammaError_Syst_TypeB_Chi2");
      if(!nameCutVariationsDoubleRatio[i].CompareTo("edEdx"))
         meanedEdxTypeGammaB = (TGraphErrors*)meanErrorsGamma[i]->Clone("GammaError_Syst_TypeB_edEdx");


      if(!nameCutVariationsGamma[i].CompareTo("edEdx") || !nameCutVariationsGamma[i].CompareTo("PidEdx") || !nameCutVariationsGamma[i].CompareTo("TOF")){
         if(!createdPID){
            meanErrorsGammaPID = (TGraphErrors*)meanErrorsGamma[i]->Clone("GammaError_Syst_PID");
            createdPID = kTRUE;
         }
         else {
            for(Int_t j = 0;j<meanErrorsGamma[i]->GetN();j++){
               Double_t x1,x2,y1,y2;
               meanErrorsGammaPID->GetPoint(j,x1,y1);
               meanErrorsGamma[i]->GetPoint(j,x2,y2);
               meanErrorsGammaPID->SetPoint(j,x1,sqrt(y1*y1+y2*y2));
            }
         }
      }
      if(!nameCutVariationsGamma[i].CompareTo("Chi2") || !nameCutVariationsGamma[i].CompareTo("PsiPair") || !nameCutVariationsGamma[i].CompareTo("qT")){
         if(!createdGammaCuts){
            meanErrorsGammaGammaCuts = (TGraphErrors*)meanErrorsGamma[i]->Clone("GammaError_Syst_GammaCuts");
            createdGammaCuts = kTRUE;
         }
         else {
            for(Int_t j = 0;j<meanErrorsGamma[i]->GetN();j++){
               Double_t x1,x2,y1,y2;
               meanErrorsGammaGammaCuts->GetPoint(j,x1,y1);
               meanErrorsGamma[i]->GetPoint(j,x2,y2);
               meanErrorsGammaGammaCuts->SetPoint(j,x1,sqrt(y1*y1+y2*y2));

            }
         }
      }
      if(!nameCutVariationsGamma[i].CompareTo("TPCClst") || !nameCutVariationsGamma[i].CompareTo("SinglePt")){
         if(!createdTrackCuts){
            meanErrorsGammaTrackCuts = (TGraphErrors*)meanErrorsGamma[i]->Clone("GammaError_Syst_TrackCuts");
            createdTrackCuts = kTRUE;
         }
         else {
            for(Int_t j = 0;j<meanErrorsGamma[i]->GetN();j++){
               Double_t x1,x2,y1,y2;
               meanErrorsGammaTrackCuts->GetPoint(j,x1,y1);
               meanErrorsGamma[i]->GetPoint(j,x2,y2);
               meanErrorsGammaTrackCuts->SetPoint(j,x1,sqrt(y1*y1+y2*y2));

            }
         }
      }
   }


   // ------------------------- Gamma Errors --------------------------------
   TCanvas *canvasGammaErrors  = GetAndSetCanvas("canvasGammaErrorsMean");
   TLegend* legendGammaErrors = GetAndSetLegend(0.13,0.63,5.2,2,Form("Inclusive #gamma Systematic Uncertainties in %s",(GetCentralityStringC(cutSel)).Data()));

   //dummy->GetYaxis()->SetRangeUser(gammaErrors[0],gammaErrors[1]);
   dummy->GetYaxis()->SetRangeUser(0,30);
   dummy->GetXaxis()->SetRangeUser(0.8,8);
   dummy->GetXaxis()->SetLabelOffset(0);
   SetHistogramm(dummy, "#it{p}_{T} (GeV/#it{c})", "Systematic Uncertanties %");
   dummy->DrawCopy();

   DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummedIncMatGamma, 20, 1.,1,1);
   DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummedWithoutMatGamma, 20, 1.,kGray+2,kGray+2);
   DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummedMaterialGamma, 20, 1.,kGray+2,kGray+2);


   for(Int_t i = 0; i< ErrorGammaSpec; i++){
      DrawGammaSetMarkerTGraphErr(meanErrorsGamma[i], 20+i, 1.,color[i],color[i]);
      meanErrorsGamma[i]->Draw("p,csame");
      cout << "trying to create legend" << endl;
      legendGammaErrors->AddEntry(meanErrorsGamma[i],nameCutVariationsGammaLegend[i].Data(),"p");
   }
   legendGammaErrors->AddEntry(meanErrorsCorrSummedIncMatGamma,"Quadratically Summed","p");

   meanErrorsCorrSummedIncMatGamma->Draw("p,csame");

   //meanErrorsCorrSummedWithoutMatGamma->Draw("p,csame");
   meanErrorsCorrSummedMaterialGamma->Draw("p,csame");
   legendGammaErrors->AddEntry(meanErrorsCorrSummedMaterialGamma,"Material Budget","p");

   legendGammaErrors->Draw();
   canvasGammaErrors->Update();

   canvasGammaErrors->Print(Form("%s/GammaErrors_%s.eps",outputDir.Data(),cent.Data()));

   // ------------------------- Gamma Errors Grouped --------------------------------
   TCanvas *canvasGammaErrorsGrouped  = GetAndSetCanvas("canvasGammaErrorsGroupedMean");
   TLegend* legendGammaErrorsGrouped = GetAndSetLegend(0.13,0.65,4,1,Form("Inclusive #gamma Systematic Uncertainties in %s",(GetCentralityStringC(cutSel)).Data()));

   //dummy->GetYaxis()->SetRangeUser(gammaErrors[0],gammaErrors[1]);
   dummy->GetYaxis()->SetRangeUser(0,30);
   dummy->GetXaxis()->SetRangeUser(0.8,8);
   dummy->GetXaxis()->SetLabelOffset(0);
   SetHistogramm(dummy, "#it{p}_{T} (GeV/#it{c})", "Systematic Uncertanties %");
   dummy->DrawCopy();


   DrawGammaSetMarkerTGraphErr(meanErrorsGammaPID, 20,1, kGreen-8, kGreen-8, 1, kTRUE);
   DrawGammaSetMarkerTGraphErr(meanErrorsGammaTrackCuts , 20,1, kRed-8, kRed-8, 1, kTRUE);
   DrawGammaSetMarkerTGraphErr(meanErrorsGammaGammaCuts , 20,1, kBlue-8, kBlue-8, 1, kTRUE);


   meanErrorsGammaPID->Draw("p,csame");
   meanErrorsGammaTrackCuts->Draw("p,csame");
   meanErrorsGammaGammaCuts->Draw("p,csame");
   meanErrorsCorrSummedMaterialGamma->Draw("p,csame");

   Double_t *SumGammaGrouped = new Double_t[meanErrorsCorrSummedIncMatGammaA->GetN()];

   for(Int_t i = 0;i<meanErrorsCorrSummedIncMatGammaA->GetN();i++)
      {
         Double_t x1,x2,x3,y1,y2,y3;
         meanErrorsGammaPID->GetPoint(i,x1,y1);
         meanErrorsGammaTrackCuts->GetPoint(i,x2,y2);
         meanErrorsGammaGammaCuts->GetPoint(i,x3,y3);
         SumGammaGrouped[i] = sqrt( y1*y1+y2*y2+y3*y3+4.5*4.5 );
         //cout<<x1<<"  "<<y1<<" "<<y2<<" "<<y3<<" "<<sqrt( y1*y1+y2*y2+y3*y3+4.5*4.5 )<<endl;
      }

   TGraphErrors *meanErrorsCorrSummedIncMatGammaGrouped = new TGraphErrors(ConstnBinsGamma,ptBinsGamma , SumGammaGrouped,ptBinsGammaErr ,errorsMeanErrCorrMatSummedGamma );
   DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummedIncMatGammaGrouped, 20, 1.,1,1);
   meanErrorsCorrSummedIncMatGammaGrouped->Draw("p,csame");

   legendGammaErrorsGrouped->AddEntry(meanErrorsGammaPID,"PID Uncertainties","pl");
   legendGammaErrorsGrouped->AddEntry(meanErrorsGammaTrackCuts,"Track Selection Uncertainties","pl");
   legendGammaErrorsGrouped->AddEntry(meanErrorsGammaGammaCuts,"Gamma Cut Uncertainties","pl");
   legendGammaErrorsGrouped->AddEntry(meanErrorsCorrSummedMaterialGamma,"Material Budget Uncertainties","pl");
   legendGammaErrorsGrouped->AddEntry(meanErrorsCorrSummedIncMatGammaGrouped,"Summed Uncertainty","pl");
   legendGammaErrorsGrouped->Draw();

   canvasGammaErrorsGrouped->Update();

   canvasGammaErrorsGrouped->Print(Form("%s/GammaErrorsGrouped_%s.eps",outputDir.Data(),cent.Data()));


   // ------------------------- Gamma Errors ABC --------------------------------
   TCanvas *canvasGammaErrorsABC  = GetAndSetCanvas("canvasGammaErrorsABCMean");
   TLegend* legendGammaErrorsABC = GetAndSetLegend(0.13,0.65,4,1,Form("Inclusive #gamma Systematic Uncertainties in %s",(GetCentralityStringC(cutSel)).Data()));

   //dummy->GetYaxis()->SetRangeUser(gammaErrors[0],gammaErrors[1]);
   dummy->GetYaxis()->SetRangeUser(0,30);
   dummy->GetXaxis()->SetRangeUser(0.8,8);
   dummy->GetXaxis()->SetLabelOffset(0);
   SetHistogramm(dummy, "#it{p}_{T} (GeV/#it{c})", "Systematic Uncertanties %");
   dummy->DrawCopy();


   DrawGammaSetMarkerTGraphErr( meanErrorsCorrSummedIncMatGammaA, 20,1, kGreen-8, kGreen-8, 1, kTRUE);
   DrawGammaSetMarkerTGraphErr( meanErrorsCorrSummedIncMatGammaB, 20,1, kRed-8, kRed-8, 1, kTRUE);
   DrawGammaSetMarkerTGraphErr( meanErrorsCorrSummedIncMatGammaC, 20,1, kBlue-8, kBlue-8, 1, kTRUE);


   meanErrorsCorrSummedIncMatGammaA->Draw("p,csame");
   meanErrorsCorrSummedIncMatGammaB->Draw("p,csame");
   meanErrorsCorrSummedIncMatGammaC->Draw("p,csame");

   Double_t *SumGammaABC = new Double_t[meanErrorsCorrSummedIncMatGammaA->GetN()];

   for(Int_t i = 0;i<meanErrorsCorrSummedIncMatGammaA->GetN();i++)
      {
         Double_t x1,x2,x3,y1,y2,y3;
         meanErrorsCorrSummedIncMatGammaA->GetPoint(i,x1,y1);
         meanErrorsCorrSummedIncMatGammaB->GetPoint(i,x2,y2);
         meanErrorsCorrSummedIncMatGammaC->GetPoint(i,x3,y3);
         SumGammaABC[i] = sqrt( y1*y1+y2*y2+y3*y3 );
      }

   TGraphErrors *meanErrorsCorrSummedIncMatGammaABC = new TGraphErrors(ConstnBinsGamma,ptBinsGamma , SumGammaABC,ptBinsGammaErr ,errorsMeanErrCorrMatSummedGamma );
   DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummedIncMatGammaABC, 20, 1.,1,1);
   meanErrorsCorrSummedIncMatGammaABC->Draw("p,csame");

   legendGammaErrorsABC->AddEntry(meanErrorsCorrSummedIncMatGammaA,"Errors Type A","pl");
   legendGammaErrorsABC->AddEntry(meanErrorsCorrSummedIncMatGammaB,"Errors Type B","pl");
   legendGammaErrorsABC->AddEntry(meanErrorsCorrSummedIncMatGammaC,"Errors Type C","pl");
   legendGammaErrorsABC->AddEntry(meanErrorsCorrSummedIncMatGammaABC,"Quad. Sum of Errors Type A,B,C","pl");
   legendGammaErrorsABC->Draw();

   canvasGammaErrorsABC->Update();

   canvasGammaErrorsABC->Print(Form("%s/GammaErrorsABC_%s.eps",outputDir.Data(),cent.Data()));





   TGraphErrors *meanPidEdxTypeDoubleRatioPi0FitB;
   TGraphErrors *meanQtTypeDoubleRatioPi0FitB;
   TGraphErrors *meanedEdxTypeDoubleRatioPi0FitB;
   TGraphErrors *meanChi2TypeDoubleRatioPi0FitB;

   createdPID = kFALSE;
   TGraphErrors *meanErrorsDoubleRatioPi0FitPID;
   createdGammaCuts = kFALSE;
   TGraphErrors *meanErrorsDoubleRatioPi0FitGammaCuts;
   createdTrackCuts = kFALSE;
   TGraphErrors *meanErrorsDoubleRatioPi0FitTrackCuts;
   Double_t createdPi0Extraction = kFALSE;
   TGraphErrors *meanErrorsDoubleRatioPi0FitPi0Extraction;
   Bool_t createdCocktail = kFALSE;
   TGraphErrors *meanErrorsDoubleRatioPi0FitCocktail;
   TGraphErrors *meanErrorsDoubleRatioPi0FitCocktailEtaNorm;

   for(Int_t i = 0; i< ErrorDoubleRatio; i++){
      DrawGammaSetMarkerTGraphErr(meanErrorsDoubleRatioPi0Fit[i], 20+i, 1.,color[i],color[i]);
      if(!nameCutVariationsDoubleRatio[i].CompareTo("PidEdx"))
         meanPidEdxTypeDoubleRatioPi0FitB = (TGraphErrors*) meanErrorsDoubleRatioPi0Fit[i]->Clone("DoubleRatioFitError_Syst_TypeB_PidEdx");
      if(!nameCutVariationsDoubleRatio[i].CompareTo("qT"))
         meanQtTypeDoubleRatioPi0FitB = (TGraphErrors*) meanErrorsDoubleRatioPi0Fit[i]->Clone("DoubleRatioFitError_Syst_TypeB_qT");
      if(!nameCutVariationsDoubleRatio[i].CompareTo("edEdx"))
         meanedEdxTypeDoubleRatioPi0FitB = (TGraphErrors*) meanErrorsDoubleRatioPi0Fit[i]->Clone("DoubleRatioFitError_Syst_TypeB_edEdx");
      if(!nameCutVariationsDoubleRatio[i].CompareTo("Chi2"))
         meanChi2TypeDoubleRatioPi0FitB = (TGraphErrors*) meanErrorsDoubleRatioPi0Fit[i]->Clone("DoubleRatioFitError_Syst_TypeB_Chi2");

      if(!nameCutVariationsDoubleRatio[i].CompareTo("edEdx") || !nameCutVariationsDoubleRatio[i].CompareTo("PidEdx") || !nameCutVariationsDoubleRatio[i].CompareTo("TOF")){
         if(!createdPID){
            meanErrorsDoubleRatioPi0FitPID = (TGraphErrors*)meanErrorsDoubleRatioPi0Fit[i]->Clone("DoubleRatioPi0FitError_Syst_PID");
            createdPID = kTRUE;
         }
         else {
            for(Int_t j = 0;j<meanErrorsDoubleRatioPi0Fit[i]->GetN();j++){
               Double_t x1,x2,y1,y2;
               meanErrorsDoubleRatioPi0FitPID->GetPoint(j,x1,y1);
               meanErrorsDoubleRatioPi0Fit[i]->GetPoint(j,x2,y2);
               meanErrorsDoubleRatioPi0FitPID->SetPoint(j,x1,sqrt(y1*y1+y2*y2));
            }
         }
      }
      if(!nameCutVariationsDoubleRatio[i].CompareTo("Chi2") || !nameCutVariationsDoubleRatio[i].CompareTo("PsiPair") || !nameCutVariationsDoubleRatio[i].CompareTo("qT")){
         if(!createdGammaCuts){
            meanErrorsDoubleRatioPi0FitGammaCuts = (TGraphErrors*)meanErrorsDoubleRatioPi0Fit[i]->Clone("DoubleRatioPi0FitError_Syst_DoubleRatioPi0FitCuts");
            createdGammaCuts = kTRUE;
         }
         else {
            for(Int_t j = 0;j<meanErrorsDoubleRatioPi0Fit[i]->GetN();j++){
               Double_t x1,x2,y1,y2;
               meanErrorsDoubleRatioPi0FitGammaCuts->GetPoint(j,x1,y1);
               meanErrorsDoubleRatioPi0Fit[i]->GetPoint(j,x2,y2);
               meanErrorsDoubleRatioPi0FitGammaCuts->SetPoint(j,x1,sqrt(y1*y1+y2*y2));

            }
         }
      }
      if(!nameCutVariationsDoubleRatio[i].CompareTo("TPCClst") || !nameCutVariationsDoubleRatio[i].CompareTo("SinglePt")){
         if(!createdTrackCuts){
            meanErrorsDoubleRatioPi0FitTrackCuts = (TGraphErrors*)meanErrorsDoubleRatioPi0Fit[i]->Clone("DoubleRatioPi0FitError_Syst_TrackCuts");
            createdTrackCuts = kTRUE;
         }
         else {
            for(Int_t j = 0;j<meanErrorsDoubleRatioPi0Fit[i]->GetN();j++){
               Double_t x1,x2,y1,y2;
               meanErrorsDoubleRatioPi0FitTrackCuts->GetPoint(j,x1,y1);
               meanErrorsDoubleRatioPi0Fit[i]->GetPoint(j,x2,y2);
               meanErrorsDoubleRatioPi0FitTrackCuts->SetPoint(j,x1,sqrt(y1*y1+y2*y2));

            }
         }
      }
      if(!nameCutVariationsDoubleRatio[i].CompareTo("Alpha") || !nameCutVariationsDoubleRatio[i].CompareTo("IntRange")){
         if(!createdPi0Extraction){
            meanErrorsDoubleRatioPi0FitPi0Extraction = (TGraphErrors*)meanErrorsDoubleRatioPi0Fit[i]->Clone("DoubleRatioPi0FitError_Syst_Pi0Extraction");
            createdPi0Extraction = kTRUE;
         }
         else {
            for(Int_t j = 0;j<meanErrorsDoubleRatioPi0Fit[i]->GetN();j++){
               Double_t x1,x2,y1,y2;
               meanErrorsDoubleRatioPi0FitPi0Extraction->GetPoint(j,x1,y1);
               meanErrorsDoubleRatioPi0Fit[i]->GetPoint(j,x2,y2);
               meanErrorsDoubleRatioPi0FitPi0Extraction->SetPoint(j,x1,sqrt(y1*y1+y2*y2));

            }
         }
      }
      if(!nameCutVariationsDoubleRatio[i].CompareTo("CocktailEta") || !nameCutVariationsDoubleRatio[i].CompareTo("CocktailParam") || !nameCutVariationsDoubleRatio[i].CompareTo("CocktailEtaNorm")){
         if(!createdCocktail){
            meanErrorsDoubleRatioPi0FitCocktail = (TGraphErrors*)meanErrorsDoubleRatioPi0Fit[i]->Clone("DoubleRatioPi0FitError_Syst_Cocktail");
            createdCocktail = kTRUE;
         }
         else {
            for(Int_t j = 0;j<meanErrorsDoubleRatioPi0Fit[i]->GetN();j++){
               Double_t x1,x2,y1,y2;
               meanErrorsDoubleRatioPi0FitCocktail->GetPoint(j,x1,y1);
               meanErrorsDoubleRatioPi0Fit[i]->GetPoint(j,x2,y2);
               meanErrorsDoubleRatioPi0FitCocktail->SetPoint(j,x1,sqrt(y1*y1+y2*y2));

            }
         }
         if(!nameCutVariationsDoubleRatio[i].CompareTo("CocktailEtaNorm")){
            meanErrorsDoubleRatioPi0FitCocktailEtaNorm = (TGraphErrors*)meanErrorsDoubleRatioPi0Fit[i]->Clone("DoubleRatioPi0FitError_Syst_CocktailEtaNorm");
         }
      }
   }



   // ------------------------- Double Ratio Errors --------------------------------
   TCanvas *canvasDoubleRatioErrors = GetAndSetCanvas("canvasDoubleRatioErrorsMean");
   TLegend* legendDoubleRatioErrors = GetAndSetLegend(0.2,0.6,6.2,2,Form("Systematics Double Ratio in %s",(GetCentralityStringC(cutSel)).Data()));
   dummy->GetXaxis()->SetLabelOffset(0);
   SetHistogramm(dummy, "#it{p}_{T} (GeV/#it{c})", "Systematic Uncertanties %");
   dummy->GetYaxis()->SetRangeUser(doubleRatioErrors[0],doubleRatioErrors[1]);
   dummy->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
   dummy->DrawCopy();
   DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummedIncMatDoubleRatio, 20, 1.,1,1);
   for(Int_t i = 0; i< ErrorDoubleRatio ; i++){
      if(nameCutVariationsDoubleRatio[i].CompareTo("Fit") != 0){
         DrawGammaSetMarkerTGraphErr(meanErrorsDoubleRatio[i], 20+i, 1.,color[i],color[i]);
         //meanErrorsDoubleRatio[i]->RemovePoint(0);
         meanErrorsDoubleRatio[i]->Draw("p,csame");

         cout << "trying to create legend" << endl;
         legendDoubleRatioErrors->AddEntry(meanErrorsDoubleRatio[i],nameCutVariationsDoubleRatio[i].Data(),"p");
      }
   }

   //meanErrorsCorrSummedIncMatDoubleRatio->RemovePoint(0);
   meanErrorsCorrSummedIncMatDoubleRatio->Draw("p,csame");
   legendDoubleRatioErrors->AddEntry(meanErrorsCorrSummedIncMatDoubleRatio,"quadratically summed plus material","p");
   legendDoubleRatioErrors->Draw();
   canvasDoubleRatioErrors->Update();

   canvasDoubleRatioErrors->Print(Form("%s/DoubleRatioErrors_%s.eps",outputDir.Data(),cent.Data()));


   // ------------------------- Double Ratio ErrorsABC --------------------------------
   TCanvas *canvasDoubleRatioErrorsABC = GetAndSetCanvas("canvasDoubleRatioErrorsABCMean");
   TLegend* legendDoubleRatioErrorsABC = GetAndSetLegend(0.2,0.65,4,1,Form("Systematics Double Ratio in %s",(GetCentralityStringC(cutSel)).Data()));

   dummy->GetXaxis()->SetLabelOffset(0);
   SetHistogramm(dummy, "#it{p}_{T} (GeV/#it{c})", "Systematic Uncertanties %");
   dummy->GetYaxis()->SetRangeUser(doubleRatioErrors[0],doubleRatioErrors[1]);
   dummy->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
   dummy->DrawCopy();



   DrawGammaSetMarkerTGraphErr( meanErrorsCorrSummedIncMatDoubleRatioA, 20,1, kGreen-8, kGreen-8, 1, kTRUE);
   DrawGammaSetMarkerTGraphErr( meanErrorsCorrSummedIncMatDoubleRatioB, 20,1, kRed-8, kRed-8, 1, kTRUE);
   DrawGammaSetMarkerTGraphErr( meanErrorsCorrSummedIncMatDoubleRatioC, 20,1, kBlue-8, kBlue-8, 1, kTRUE);


   meanErrorsCorrSummedIncMatDoubleRatioA->Draw("p,csame");
   meanErrorsCorrSummedIncMatDoubleRatioB->Draw("p,csame");
   meanErrorsCorrSummedIncMatDoubleRatioC->Draw("p,csame");

   Double_t *SumDoubleRatioABC = new Double_t[meanErrorsCorrSummedIncMatDoubleRatioA->GetN()];

   for(Int_t i = 0;i<meanErrorsCorrSummedIncMatDoubleRatioA->GetN();i++)
      {
         Double_t x1,x2,x3,y1,y2,y3;
         meanErrorsCorrSummedIncMatDoubleRatioA->GetPoint(i,x1,y1);
         meanErrorsCorrSummedIncMatDoubleRatioB->GetPoint(i,x2,y2);
         meanErrorsCorrSummedIncMatDoubleRatioC->GetPoint(i,x3,y3);
         SumDoubleRatioABC[i] = sqrt( y1*y1+y2*y2+y3*y3 );
      }

   TGraphErrors *meanErrorsCorrSummedIncMatDoubleRatioABC = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio , SumDoubleRatioABC,ptBinsDoubleRatioErr ,errorsMeanErrCorrMatSummedDoubleRatio );
   DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummedIncMatDoubleRatioABC, 20, 1.,1,1);
   meanErrorsCorrSummedIncMatDoubleRatioABC->Draw("p,csame");


   legendDoubleRatioErrorsABC->AddEntry(meanErrorsCorrSummedIncMatDoubleRatioA,"Errors Type A","pl");
   legendDoubleRatioErrorsABC->AddEntry(meanErrorsCorrSummedIncMatDoubleRatioB,"Errors Type B","pl");
   legendDoubleRatioErrorsABC->AddEntry(meanErrorsCorrSummedIncMatDoubleRatioC,"Errors Type C","pl");
   legendDoubleRatioErrorsABC->AddEntry(meanErrorsCorrSummedIncMatDoubleRatioABC,"Quad. Sum of Errors Type A,B,C","pl");
   legendDoubleRatioErrorsABC->Draw();

   canvasDoubleRatioErrorsABC->Update();

   canvasDoubleRatioErrorsABC->Print(Form("%s/DoubleRatioErrorsABC_%s.eps",outputDir.Data(),cent.Data()));

   // ------------------------- Double Ratio Fit Errors --------------------------------


   TCanvas *canvasDoubleRatioPi0FitErrors  = GetAndSetCanvas("canvasDoubleRatioPi0FitErrorsMean");
   TLegend* legendDoubleRatioPi0FitErrors = GetAndSetLegend(0.13,0.60,6.2,2,Form("Systematics Double Ratio in %s",(GetCentralityStringC(cutSel)).Data()));

   dummy->GetYaxis()->SetRangeUser(0,35);
   SetHistogramm(dummy, "#it{p}_{T} (GeV/#it{c})", "Systematic Uncertanties %");
   dummy->GetYaxis()->SetRangeUser(doubleRatioErrors[0],doubleRatioErrors[1]);
   dummy->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
   dummy->DrawCopy();
   DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummedIncMatDoubleRatioPi0Fit, 20, 1.,1,1);
   for(Int_t i = 0; i< ErrorDoubleRatio ; i++){
      if(i==ErrorDoubleRatio)color[i]=2;
      if(nameCutVariationsDoubleRatio[i].CompareTo("Fit") == 0)color[i]=3;
      DrawGammaSetMarkerTGraphErr(meanErrorsDoubleRatioPi0Fit[i], 20+i, 1.,color[i],color[i]);

      if(!nameCutVariationsDoubleRatio[i].CompareTo("CocktailEta") ||
         !nameCutVariationsDoubleRatio[i].CompareTo("CocktailParam") ||
         !nameCutVariationsDoubleRatio[i].CompareTo("CocktailEtaNorm") ){
         for(Int_t j = 0;j<meanErrorsCorrSummedMaterialDoubleRatioPi0FitCocktail->GetN();j++){
            Double_t x1,y1,x2,y2;
            meanErrorsCorrSummedMaterialDoubleRatioPi0FitCocktail->GetPoint(j,x1,y1);
            meanErrorsDoubleRatioPi0Fit[i]->GetPoint(j,x2,y2);
            meanErrorsCorrSummedMaterialDoubleRatioPi0FitCocktail->SetPoint(j,x1,sqrt(y1*y1+y2*y2));
         }
      }
      meanErrorsDoubleRatioPi0Fit[i]->Draw("p,csame");

      cout << "trying to create legend" << endl;
      legendDoubleRatioPi0FitErrors->AddEntry(meanErrorsDoubleRatioPi0Fit[i],nameCutVariationsGammaLegend[i].Data(),"p");
   }


   meanErrorsCorrSummedIncMatDoubleRatioPi0Fit->Draw("p,csame");
   meanErrorsCorrSummedMaterialGamma->Draw("p,csame");
   legendDoubleRatioPi0FitErrors->AddEntry(meanErrorsCorrSummedMaterialGamma,"Material Budget","p");


   legendDoubleRatioPi0FitErrors->AddEntry(meanErrorsCorrSummedIncMatDoubleRatioPi0Fit,"Quadratically Summed","p");
   legendDoubleRatioPi0FitErrors->Draw();
   canvasDoubleRatioPi0FitErrors->Update();


   canvasDoubleRatioPi0FitErrors->Print(Form("%s/DoubleRatioFitErrors_%s.eps",outputDir.Data(),cent.Data()));


   // ------------------------- Double Ratio Fit ErrorsGrouped --------------------------------
   TCanvas *canvasDoubleRatioPi0FitErrorsGrouped = GetAndSetCanvas("canvasDoubleRatioPi0FitErrorsGroupedMean");
   TLegend* legendDoubleRatioPi0FitErrorsGrouped = GetAndSetLegend(0.17,0.55,6,1,Form("Systematics Double Ratio in %s",(GetCentralityStringC(cutSel)).Data()));

   dummy->GetXaxis()->SetLabelOffset(0);
   SetHistogramm(dummy, "#it{p}_{T} (GeV/#it{c})", "Systematic Uncertanties %");
   dummy->GetYaxis()->SetRangeUser(doubleRatioErrors[0],doubleRatioErrors[1]);
   dummy->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
   dummy->DrawCopy();


   DrawGammaSetMarkerTGraphErr(meanErrorsDoubleRatioPi0FitPID , 20,1, kGreen-8, kGreen-8, 1, kTRUE);
   DrawGammaSetMarkerTGraphErr(meanErrorsDoubleRatioPi0FitGammaCuts , 20,1, kRed-8, kRed-8, 1, kTRUE);
   DrawGammaSetMarkerTGraphErr(meanErrorsDoubleRatioPi0FitTrackCuts , 20,1, kBlue-8, kBlue-8, 1, kTRUE);
   DrawGammaSetMarkerTGraphErr(meanErrorsDoubleRatioPi0FitPi0Extraction , 20,1, kCyan-8, kCyan-8, 1, kTRUE);
   DrawGammaSetMarkerTGraphErr(meanErrorsDoubleRatioPi0FitCocktail , 20,1, kYellow-8, kYellow-8, 1, kTRUE);

   meanErrorsDoubleRatioPi0FitPID->Draw("p,csame");
   meanErrorsDoubleRatioPi0FitGammaCuts->Draw("p,csame");
   meanErrorsDoubleRatioPi0FitTrackCuts->Draw("p,csame");
   meanErrorsDoubleRatioPi0FitPi0Extraction->Draw("p,csame");
   meanErrorsDoubleRatioPi0FitCocktail->Draw("p,csame");
   meanErrorsCorrSummedMaterialGamma->Draw("p,csame");

   Double_t *SumDoubleRatioPi0FitGrouped = new Double_t[meanErrorsCorrSummedIncMatDoubleRatioPi0FitA->GetN()];

   for(Int_t i = 0;i<meanErrorsCorrSummedIncMatDoubleRatioPi0FitA->GetN();i++)
      {
         Double_t x1,x2,x3,x4,x5,y1,y2,y3,y4,y5;
         meanErrorsDoubleRatioPi0FitPID->GetPoint(i,x1,y1);
         meanErrorsDoubleRatioPi0FitGammaCuts->GetPoint(i,x2,y2);
         meanErrorsDoubleRatioPi0FitTrackCuts->GetPoint(i,x3,y3);
         meanErrorsDoubleRatioPi0FitPi0Extraction->GetPoint(i,x4,y4);
         meanErrorsDoubleRatioPi0FitCocktail->GetPoint(i,x5,y5);
         SumDoubleRatioPi0FitGrouped[i] = sqrt( y1*y1+y2*y2+y3*y3+y4*y4+y5*y5+4.5*4.5 );
         cout<<x1<<" Gamma "<<y2<<" PID  "<<y1<<" Track "<<y3<<" Pi0 "<<y4<<" Cocktail "<<y5<<" Sum "<<sqrt( y1*y1+y2*y2+y3*y3+y4*y4+y5*y5+4.5*4.5 )<<endl;
      }
   //return;

   TGraphErrors *meanErrorsCorrSummedIncMatDoubleRatioPi0FitGrouped = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio , SumDoubleRatioPi0FitGrouped,ptBinsDoubleRatioErr ,errorsMeanErrCorrMatSummedDoubleRatioPi0Fit );
   DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummedIncMatDoubleRatioPi0FitGrouped, 20, 1.,1,1);
   meanErrorsCorrSummedIncMatDoubleRatioPi0FitGrouped->Draw("p,csame");

   legendDoubleRatioPi0FitErrorsGrouped->AddEntry(meanErrorsDoubleRatioPi0FitPID,"PID Uncertainties","pl");
   legendDoubleRatioPi0FitErrorsGrouped->AddEntry(meanErrorsDoubleRatioPi0FitTrackCuts,"Track Selection Uncertainties","pl");
   legendDoubleRatioPi0FitErrorsGrouped->AddEntry(meanErrorsDoubleRatioPi0FitGammaCuts,"#gamma Cut Uncertainties","pl");
   legendDoubleRatioPi0FitErrorsGrouped->AddEntry(meanErrorsDoubleRatioPi0FitPi0Extraction,"#pi^{0} Signal Extraction","pl");
   legendDoubleRatioPi0FitErrorsGrouped->AddEntry(meanErrorsDoubleRatioPi0FitCocktail,"Decay Photon Cocktail","pl");
   legendDoubleRatioPi0FitErrorsGrouped->AddEntry(meanErrorsCorrSummedMaterialGamma,"Material Budget Uncertainties","pl");
   legendDoubleRatioPi0FitErrorsGrouped->AddEntry(meanErrorsCorrSummedIncMatDoubleRatioPi0FitGrouped,"Summed Uncertainty","pl");
   legendDoubleRatioPi0FitErrorsGrouped->Draw();

   canvasDoubleRatioPi0FitErrorsGrouped->Update();

   canvasDoubleRatioPi0FitErrorsGrouped->Print(Form("%s/DoubleRatioPi0FitErrorsGrouped_%s.eps",outputDir.Data(),cent.Data()));



   // ------------------------- Double Ratio Fit ErrorsABC --------------------------------
   TCanvas *canvasDoubleRatioPi0FitErrorsABC = GetAndSetCanvas("canvasDoubleRatioPi0FitErrorsABCMean");
   TLegend* legendDoubleRatioPi0FitErrorsABC = GetAndSetLegend(0.2,0.65,4,1,Form("Systematics Double Ratio in %s",(GetCentralityStringC(cutSel)).Data()));

   dummy->GetXaxis()->SetLabelOffset(0);
   SetHistogramm(dummy, "#it{p}_{T} (GeV/#it{c})", "Systematic Uncertanties %");
   dummy->GetYaxis()->SetRangeUser(doubleRatioErrors[0],doubleRatioErrors[1]);
   dummy->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
   dummy->DrawCopy();

   DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummedIncMatDoubleRatioPi0FitA , 20,1, kGreen-8, kGreen-8, 1, kTRUE);
   DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummedIncMatDoubleRatioPi0FitB , 20,1, kRed-8, kRed-8, 1, kTRUE);
   DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummedIncMatDoubleRatioPi0FitC , 20,1, kBlue-8, kBlue-8, 1, kTRUE);


   meanErrorsCorrSummedIncMatDoubleRatioPi0FitA->Draw("p,csame");
   meanErrorsCorrSummedIncMatDoubleRatioPi0FitB->Draw("p,csame");
   meanErrorsCorrSummedIncMatDoubleRatioPi0FitC->Draw("p,csame");

   Double_t *SumDoubleRatioPi0FitABC = new Double_t[meanErrorsCorrSummedIncMatDoubleRatioPi0FitA->GetN()];

   for(Int_t i = 0;i<meanErrorsCorrSummedIncMatDoubleRatioPi0FitA->GetN();i++)
      {
         Double_t x1,x2,x3,y1,y2,y3;
         meanErrorsCorrSummedIncMatDoubleRatioPi0FitA->GetPoint(i,x1,y1);
         meanErrorsCorrSummedIncMatDoubleRatioPi0FitB->GetPoint(i,x2,y2);
         meanErrorsCorrSummedIncMatDoubleRatioPi0FitC->GetPoint(i,x3,y3);
         SumDoubleRatioPi0FitABC[i] = sqrt( y1*y1+y2*y2+y3*y3 );
         
         //cout<<x1<<" DR ABC  Type A "<<y1<<" Type B  "<<y2<<" Type C "<<y3<<" Sum "<<sqrt( y1*y1+y2*y2+y3*y3)<<endl;
         
      }


   TGraphErrors *meanErrorsCorrSummedIncMatDoubleRatioPi0FitABC = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio , SumDoubleRatioPi0FitABC,ptBinsDoubleRatioErr ,errorsMeanErrCorrMatSummedDoubleRatioPi0Fit );
   DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummedIncMatDoubleRatioPi0FitABC, 20, 1.,1,1);
   meanErrorsCorrSummedIncMatDoubleRatioPi0FitABC->Draw("p,csame");

   legendDoubleRatioPi0FitErrorsABC->AddEntry(meanErrorsCorrSummedIncMatDoubleRatioPi0FitA,"Errors Type A","pl");
   legendDoubleRatioPi0FitErrorsABC->AddEntry(meanErrorsCorrSummedIncMatDoubleRatioPi0FitB,"Errors Type B","pl");
   legendDoubleRatioPi0FitErrorsABC->AddEntry(meanErrorsCorrSummedIncMatDoubleRatioPi0FitC,"Errors Type C","pl");
   legendDoubleRatioPi0FitErrorsABC->AddEntry(meanErrorsCorrSummedIncMatDoubleRatioPi0FitABC,"Quad. Sum of Errors Type A,B,C","pl");
   legendDoubleRatioPi0FitErrorsABC->Draw();

   canvasDoubleRatioPi0FitErrorsABC->Update();

   canvasDoubleRatioPi0FitErrorsABC->Print(Form("%s/DoubleRatioFitErrorsABC_%s.eps",outputDir.Data(),cent.Data()));


   // ------------------------- Double Ratio  --------------------------------

   Int_t nbinsDoubleRatio = meanErrorsCorrSummedIncMatDoubleRatioPi0Fit->GetN();
   Double_t *xDoubleRatio = meanErrorsCorrSummedIncMatDoubleRatioPi0Fit->GetX();
   Double_t *xErrorDoubleRatio = meanErrorsCorrSummedIncMatDoubleRatio->GetEX();

   Double_t *yGraphDoubleRatio = meanErrorsCorrSummedIncMatDoubleRatio->GetY();
   Double_t *yGraphDoubleRatioA = meanErrorsCorrSummedIncMatDoubleRatioA->GetY();
   Double_t *yGraphDoubleRatioB = meanErrorsCorrSummedIncMatDoubleRatioB->GetY();
   Double_t *yGraphDoubleRatioC = meanErrorsCorrSummedIncMatDoubleRatioC->GetY();
   Double_t *yGraphDoubleRatioPi0Fit = meanErrorsCorrSummedIncMatDoubleRatioPi0Fit->GetY();
   Double_t *yGraphDoubleRatioCharged = meanErrorsCorrSummedIncMatGamma->GetY();
   Double_t *yGraphDoubleRatioChargedLow = meanErrorsCorrSummedIncMatGamma->GetY();
   Double_t *yGraphDoubleRatioChargedHigh = meanErrorsCorrSummedIncMatGamma->GetY();

   Double_t *yGraphDoubleRatioPi0FitPtCorr = meanErrorsCorrSummedWithoutMatDoubleRatioPi0Fit->GetY();
   Double_t *yGraphDoubleRatioPi0FitPtUnCorr = meanErrorsCorrSummedMaterialDoubleRatioPi0Fit->GetY();

   Double_t *yGraphDoubleRatioPi0FitA = meanErrorsCorrSummedIncMatDoubleRatioPi0FitA->GetY();
   Double_t *yGraphDoubleRatioPi0FitB = meanErrorsCorrSummedIncMatDoubleRatioPi0FitB->GetY();
   Double_t *yGraphDoubleRatioPi0FitC = meanErrorsCorrSummedIncMatDoubleRatioPi0FitC->GetY();
   
   Double_t *yGraphDoubleRatioPi0FitPidEdxB = meanPidEdxTypeDoubleRatioPi0FitB->GetY();
   Double_t *yGraphDoubleRatioPi0FitqTB = meanQtTypeDoubleRatioPi0FitB->GetY();
   Double_t *yGraphDoubleRatioPi0FitChi2B = meanChi2TypeDoubleRatioPi0FitB->GetY();
   Double_t *yGraphDoubleRatioPi0FitedEdxB = meanedEdxTypeDoubleRatioPi0FitB->GetY();

   Double_t *yGraphDoubleRatioPi0FitMaterialC = meanErrorsCorrSummedMaterialGamma->GetY();
   Double_t *yGraphDoubleRatioPi0FitEtaNormC= meanErrorsDoubleRatioPi0FitCocktailEtaNorm->GetY();

   Double_t *yDoubleRatio = new Double_t[nbinsDoubleRatio];
   Double_t *yDoubleRatioA = new Double_t[nbinsDoubleRatio];
   Double_t *yDoubleRatioB = new Double_t[nbinsDoubleRatio];
   Double_t *yDoubleRatioC = new Double_t[nbinsDoubleRatio];
   Double_t *yDoubleRatioPi0Fit = new Double_t[nbinsDoubleRatio];
   Double_t *yDoubleRatioCharged = new Double_t[nbinsDoubleRatio];
   Double_t *yDoubleRatioChargedLow = new Double_t[nbinsDoubleRatio];
   Double_t *yDoubleRatioChargedHigh = new Double_t[nbinsDoubleRatio];
   Double_t *errorDoubleRatio = new Double_t[nbinsDoubleRatio];
   Double_t *errorDoubleRatioA = new Double_t[nbinsDoubleRatio];
   Double_t *errorDoubleRatioB = new Double_t[nbinsDoubleRatio];
   Double_t *errorDoubleRatioC = new Double_t[nbinsDoubleRatio];
   Double_t *errorDoubleRatioPi0Fit = new Double_t[nbinsDoubleRatio];
   Double_t *errorDoubleRatioPi0WithFit = new Double_t[nbinsDoubleRatio];
   Double_t *errorDoubleRatioCharged = new Double_t[nbinsDoubleRatio];
   Double_t *errorDoubleRatioChargedLow = new Double_t[nbinsDoubleRatio];
   Double_t *errorDoubleRatioChargedHigh = new Double_t[nbinsDoubleRatio];

   Double_t *errorDoubleRatioPi0FitPtCorr = new Double_t[nbinsDoubleRatio];
   Double_t *errorDoubleRatioPi0FitPtUnCorr = new Double_t[nbinsDoubleRatio];
   Double_t *errorDoubleRatioPi0FitA = new Double_t[nbinsDoubleRatio];
   Double_t *errorDoubleRatioPi0FitB = new Double_t[nbinsDoubleRatio];
   Double_t *errorDoubleRatioPi0FitC = new Double_t[nbinsDoubleRatio];

   Double_t *errorDoubleRatioPi0FitPidEdxB = new Double_t[nbinsDoubleRatio];
   Double_t *errorDoubleRatioPi0FitqTB = new Double_t[nbinsDoubleRatio];
   Double_t *errorDoubleRatioPi0FitedEdxB = new Double_t[nbinsDoubleRatio];
   Double_t *errorDoubleRatioPi0FitChi2B = new Double_t[nbinsDoubleRatio];

   Double_t *errorDoubleRatioPi0FitMaterialC = new Double_t[nbinsDoubleRatio];
   Double_t *errorDoubleRatioPi0FitEtaNormC = new Double_t[nbinsDoubleRatio];

   // Cocktail Error

   Double_t *yGraphCocktailCharged = meanErrorsCorrSummedMaterialDoubleRatioPi0FitCocktail->GetY();
   Double_t *yCocktailCharged = new Double_t[nbinsDoubleRatio];
   Double_t *errorCocktailCharged = new Double_t[nbinsDoubleRatio];


   for(Int_t i = 0; i<nbinsDoubleRatio; i++){

      yDoubleRatio[i] = DoubleRatio->GetBinContent(i+2);
      yDoubleRatioA[i] = DoubleRatio->GetBinContent(i+2);
      yDoubleRatioB[i] = DoubleRatio->GetBinContent(i+2);
      yDoubleRatioB[i] = DoubleRatio->GetBinContent(i+2);

      yDoubleRatioPi0Fit[i] = DoubleRatioPi0Fit->GetBinContent(i+2);
      yDoubleRatioCharged[i] = DoubleRatioChargedPions->GetBinContent(i+2);
      errorDoubleRatio[i] = yGraphDoubleRatio[i]*(DoubleRatio->GetBinContent(i+2)/100);
      errorDoubleRatioA[i] = yGraphDoubleRatioA[i]*(DoubleRatio->GetBinContent(i+2)/100);
      errorDoubleRatioB[i] = yGraphDoubleRatioB[i]*(DoubleRatio->GetBinContent(i+2)/100);
      errorDoubleRatioC[i] = yGraphDoubleRatioC[i]*(DoubleRatio->GetBinContent(i+2)/100);
      errorDoubleRatioPi0Fit[i] = yGraphDoubleRatioPi0Fit[i]*(DoubleRatioPi0Fit->GetBinContent(i+2)/100);
      errorDoubleRatioPi0WithFit[i] = yGraphDoubleRatioPi0Fit[i]*(DoubleRatio->GetBinContent(i+2)/100);
      errorDoubleRatioCharged[i] = yGraphDoubleRatioCharged[i]*(DoubleRatioChargedPions->GetBinContent(i+2)/100);

      errorDoubleRatioPi0FitPtCorr[i] = yGraphDoubleRatioPi0FitPtCorr[i]*(DoubleRatioPi0Fit->GetBinContent(i+2)/100);
      errorDoubleRatioPi0FitPtUnCorr[i] = yGraphDoubleRatioPi0FitPtUnCorr[i]*(DoubleRatioPi0Fit->GetBinContent(i+2)/100);
      errorDoubleRatioPi0FitA[i] = yGraphDoubleRatioPi0FitA[i]*(DoubleRatioPi0Fit->GetBinContent(i+2)/100);
      errorDoubleRatioPi0FitB[i] = yGraphDoubleRatioPi0FitB[i]*(DoubleRatioPi0Fit->GetBinContent(i+2)/100);
      errorDoubleRatioPi0FitC[i] = yGraphDoubleRatioPi0FitC[i]*(DoubleRatioPi0Fit->GetBinContent(i+2)/100);

      errorDoubleRatioPi0FitPidEdxB[i] = yGraphDoubleRatioPi0FitPidEdxB[i]*(DoubleRatioPi0Fit->GetBinContent(i+2)/100);
      errorDoubleRatioPi0FitqTB[i] = yGraphDoubleRatioPi0FitqTB[i]*(DoubleRatioPi0Fit->GetBinContent(i+2)/100);
      errorDoubleRatioPi0FitChi2B[i] = yGraphDoubleRatioPi0FitChi2B[i]*(DoubleRatioPi0Fit->GetBinContent(i+2)/100);
      errorDoubleRatioPi0FitedEdxB[i] = yGraphDoubleRatioPi0FitedEdxB[i]*(DoubleRatioPi0Fit->GetBinContent(i+2)/100);

      yCocktailCharged[i] = CocktailRatioChargedPions->GetBinContent(i+2);
      errorCocktailCharged[i] = yGraphCocktailCharged[i]*(CocktailRatioChargedPions->GetBinContent(i+2)/100);

      errorDoubleRatioPi0FitMaterialC[i] = yGraphDoubleRatioPi0FitMaterialC[i]*(DoubleRatioPi0Fit->GetBinContent(i+2)/100);
      errorDoubleRatioPi0FitEtaNormC[i] = yGraphDoubleRatioPi0FitEtaNormC[i]*(DoubleRatioPi0Fit->GetBinContent(i+2)/100);

   }



   TGraphErrors *DoubleRatioErrors = new TGraphErrors(nbinsDoubleRatio,xDoubleRatio,yDoubleRatio,xErrorDoubleRatio,errorDoubleRatio);
   TGraphErrors *DoubleRatioErrorsA = new TGraphErrors(nbinsDoubleRatio,xDoubleRatio,yDoubleRatio,xErrorDoubleRatio,errorDoubleRatioA);
   TGraphErrors *DoubleRatioErrorsB = new TGraphErrors(nbinsDoubleRatio,xDoubleRatio,yDoubleRatio,xErrorDoubleRatio,errorDoubleRatioB);
   TGraphErrors *DoubleRatioErrorsC = new TGraphErrors(nbinsDoubleRatio,xDoubleRatio,yDoubleRatio,xErrorDoubleRatio,errorDoubleRatioC);

   TGraphErrors *DoubleRatioPi0FitErrors = new TGraphErrors(nbinsDoubleRatio,xDoubleRatio,yDoubleRatioPi0Fit,xErrorDoubleRatio,errorDoubleRatioPi0Fit);
   TGraphErrors *DoubleRatioPi0FitErrorsA = new TGraphErrors(nbinsDoubleRatio,xDoubleRatio,yDoubleRatioPi0Fit,xErrorDoubleRatio,errorDoubleRatioPi0FitA);
   TGraphErrors *DoubleRatioPi0FitErrorsB = new TGraphErrors(nbinsDoubleRatio,xDoubleRatio,yDoubleRatioPi0Fit,xErrorDoubleRatio,errorDoubleRatioPi0FitB);
   TGraphErrors *DoubleRatioPi0FitErrorsC = new TGraphErrors(nbinsDoubleRatio,xDoubleRatio,yDoubleRatioPi0Fit,xErrorDoubleRatio,errorDoubleRatioPi0FitC);

   TGraphErrors *DoubleRatioPi0FitErrorsPidEdxB = new TGraphErrors(nbinsDoubleRatio,xDoubleRatio,yDoubleRatioPi0Fit,xErrorDoubleRatio,errorDoubleRatioPi0FitPidEdxB);
   TGraphErrors *DoubleRatioPi0FitErrorsqTB = new TGraphErrors(nbinsDoubleRatio,xDoubleRatio,yDoubleRatioPi0Fit,xErrorDoubleRatio,errorDoubleRatioPi0FitqTB);
   TGraphErrors *DoubleRatioPi0FitErrorsedEdxB = new TGraphErrors(nbinsDoubleRatio,xDoubleRatio,yDoubleRatioPi0Fit,xErrorDoubleRatio,errorDoubleRatioPi0FitedEdxB);
   TGraphErrors *DoubleRatioPi0FitErrorsChi2B = new TGraphErrors(nbinsDoubleRatio,xDoubleRatio,yDoubleRatioPi0Fit,xErrorDoubleRatio,errorDoubleRatioPi0FitChi2B);
   TGraphErrors *DoubleRatioPi0FitErrorsMaterialC = new TGraphErrors(nbinsDoubleRatio,xDoubleRatio,yDoubleRatioPi0Fit,xErrorDoubleRatio,errorDoubleRatioPi0FitMaterialC);
   TGraphErrors *DoubleRatioPi0FitErrorsEtaNormC = new TGraphErrors(nbinsDoubleRatio,xDoubleRatio,yDoubleRatioPi0Fit,xErrorDoubleRatio,errorDoubleRatioPi0FitEtaNormC);

   TGraphAsymmErrors *DoubleRatioPi0FitErrorsAAsym= new TGraphAsymmErrors(nbinsDoubleRatio,xDoubleRatio,yDoubleRatioPi0Fit,xErrorDoubleRatio,xErrorDoubleRatio,errorDoubleRatioPi0FitA,errorDoubleRatioPi0FitA);
   TGraphAsymmErrors *DoubleRatioPi0FitErrorsBAsym= new TGraphAsymmErrors(nbinsDoubleRatio,xDoubleRatio,yDoubleRatioPi0Fit,xErrorDoubleRatio,xErrorDoubleRatio,errorDoubleRatioPi0FitB,errorDoubleRatioPi0FitB);
   TGraphAsymmErrors *DoubleRatioPi0FitErrorsCAsym= new TGraphAsymmErrors(nbinsDoubleRatio,xDoubleRatio,yDoubleRatioPi0Fit,xErrorDoubleRatio,xErrorDoubleRatio,errorDoubleRatioPi0FitC,errorDoubleRatioPi0FitC);
   TGraphAsymmErrors *DoubleRatioPi0FitErrorsAsym = new TGraphAsymmErrors(nbinsDoubleRatio,xDoubleRatio,yDoubleRatioPi0Fit,xErrorDoubleRatio,xErrorDoubleRatio,errorDoubleRatioPi0Fit,errorDoubleRatioPi0Fit);
   TGraphErrors *DoubleRatioWithFitErrors = new TGraphErrors(nbinsDoubleRatio,xDoubleRatio,yDoubleRatio,xErrorDoubleRatio,errorDoubleRatioPi0WithFit);

   TGraphErrors *DoubleRatioPi0FitErrorsPtCorr = new TGraphErrors(nbinsDoubleRatio,meanErrorsCorrSummedWithoutMatDoubleRatioPi0Fit->GetX(),yDoubleRatioPi0Fit,meanErrorsCorrSummedWithoutMatDoubleRatioPi0Fit->GetEX(),errorDoubleRatioPi0FitPtCorr);
   TGraphErrors *DoubleRatioPi0FitErrorsPtUnCorr = new TGraphErrors(nbinsDoubleRatio,meanErrorsCorrSummedMaterialDoubleRatioPi0Fit->GetX(),yDoubleRatioPi0Fit,meanErrorsCorrSummedMaterialDoubleRatioPi0Fit->GetEX(),errorDoubleRatioPi0FitPtUnCorr);

   TGraphErrors *CocktailChargedErrors = new TGraphErrors(nbinsDoubleRatio,xDoubleRatio,yCocktailCharged,xErrorDoubleRatio,errorCocktailCharged);


   TCanvas *canvasDoubleRatio = GetAndSetCanvas("canvasDoubleRatioFinal");
   canvasDoubleRatio->SetLogx();
   DrawGammaSetMarker(DoubleRatio, 20, 2., 2, 1);
   DoubleRatio->SetLineColor(colorcent);
   DoubleRatio->SetFillColor(0);
   DoubleRatio->SetMarkerColor(colorcent);


   SetHistogramm(dummy, "#it{p}_{T} (GeV/#it{c})", "(#gamma_{inc}/#pi^{0})/(#gamma_{decay}/#pi^{0})");
   dummy->GetYaxis()->SetTitleSize(0.05);
   dummy->GetYaxis()->SetTitleOffset(.85);
   dummy->GetXaxis()->SetTitleSize(0.05);
   dummy->GetXaxis()->SetTitleOffset(.85);
   dummy->GetYaxis()->SetLabelSize(0.045);
   dummy->GetXaxis()->SetLabelSize(0.045);
   dummy->GetXaxis()->SetLabelOffset(-0.015);
   dummy->GetYaxis()->SetRangeUser(doubleRatio[0], doubleRatio[1]);
   dummy->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
   dummy->DrawCopy();

   TH1D *DoubleRatioRebinedNorm = (TH1D*) DoubleRatio->Clone("DoubleRatioRebinedNorm");
   DoubleRatioRebinedNorm = RebinTH1D(DoubleRatioRebinedNorm,histoBinningRatio);
   TH1D *DoubleRatioRebinedFit = (TH1D*) DoubleRatioPi0Fit->Clone("DoubleRatioRebined");
   DoubleRatioRebinedFit = RebinTH1D(DoubleRatioRebinedFit,histoBinningRatio);
   TH1D *DoubleRatioRebinedPions =  (TH1D*) DoubleRatioChargedPions->Clone("DoubleRatioRebinedPions");
   DoubleRatioChargedPions = RebinTH1D(DoubleRatioRebinedPions,histoBinningRatio);

   DrawGammaSetMarkerTGraphErr(DoubleRatioErrors , 20,1, colorcent, colorcent, 1, kTRUE);
   //DoubleRatioErrors->Draw("e2,same");
   One->Draw("same");
   NLODoubleRatio->Draw("p3lsame");
   DrawGammaSetMarkerTGraphErr(DoubleRatioPi0FitErrors , 20,1, colorcent, colorcent, 1, kTRUE);
   DoubleRatioPi0FitErrors->Draw("E2same");
   //DoubleRatio->DrawCopy("same");
   DoubleRatioPi0Fit->DrawCopy("same");


   TLegend* legendDoubleRatio = GetAndSetLegend(0.15,0.75,4);
   legendDoubleRatio->AddEntry(DoubleRatioErrors,Form("ALICE, %s Pb-Pb #sqrt{s_{_{NN}}} = 2.76 TeV","pf",GetCentralityStringC(cutSel).Data()),"pf");
   legendDoubleRatio->AddEntry(NLODoubleRatio,"NLO prediction: 1 + (N_{coll}#gamma_{direct,pp,NLO}/#gamma_{decay})","l");
   //legendDoubleRatio->AddEntry((TObject*)0, "for #mu = 0.5,1.0,2.0 #it{p}_{T}", "");
   legendDoubleRatio->AddEntry((TObject*)0, "for #mu = 0.5 to 2.0 #it{p}_{T}", "");
   legendDoubleRatio->Draw();


   canvasDoubleRatio->Print(Form("%s/DoubleRatio_%s.eps",outputDir.Data(),cent.Data()));
   canvasDoubleRatio->Print(Form("%s/DoubleRatio_%s.C",outputDir.Data(),cent.Data()));


   // ------------------------- Double Ratio ABC --------------------------------

   DrawGammaSetMarker(DoubleRatio, 20, 2., 1, 1);

   TCanvas *canvasDoubleRatioABC = GetAndSetCanvas("canvasDoubleRatioFinal");
   canvasDoubleRatioABC->SetLogx();

   dummy->GetYaxis()->SetRangeUser(doubleRatio[0], doubleRatio[1]);
   dummy->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
   dummy->DrawCopy();

   DoubleRatio->SetMarkerColor(colorcent);
   DoubleRatio->SetLineColor(colorcent);
   DoubleRatio->SetFillColor(0);


   DoubleRatio = RebinTH1D(DoubleRatio,histoBinningRatio);
   One->Draw("same");
   NLODoubleRatio->Draw("p3lsame");
   DoubleRatio->DrawCopy("same");
   DrawGammaSetMarkerTGraphErr(DoubleRatioErrorsA , 20,1, kGreen-8, kGreen-8, 1, kTRUE);
   DrawGammaSetMarkerTGraphErr(DoubleRatioErrorsB , 20,1, kRed-8, kRed-8, 1, kTRUE);
   DrawGammaSetMarkerTGraphErr(DoubleRatioErrorsC , 20,1, kBlue-8, kBlue-8, 1, kTRUE);


   DoubleRatioErrorsA->Draw("E2same");
   DoubleRatioErrorsB->Draw("E2same");
   DoubleRatioErrorsC->Draw("E2same");
   DoubleRatioErrors->Draw("E2same");

   DrawGammaSetMarkerTGraphErr(GraphErrorsDummy , 20,2, colorcent, colorcent, 1, kTRUE);


   TLegend* legendDoubleRatioABC = GetAndSetLegend(0.15,0.55,7);
   legendDoubleRatioABC->AddEntry(GraphErrorsDummy,Form("ALICE, %s Pb-Pb #sqrt{s_{_{NN}}} = 2.76 TeV",GetCentralityStringC(cutSel).Data()),"pf");
   legendDoubleRatioABC->AddEntry(DoubleRatioErrorsA,"Error Type A","f");
   legendDoubleRatioABC->AddEntry(DoubleRatioErrorsB,"Error Type B","f");
   legendDoubleRatioABC->AddEntry(DoubleRatioErrorsC,"Error Type C","f");
   legendDoubleRatioABC->AddEntry(NLODoubleRatio,"NLO prediction: 1 + (N_{coll}#gamma_{direct,pp,NLO}/#gamma_{decay})","l");
   legendDoubleRatioABC->AddEntry((TObject*)0, "for #mu = 0.5 to 2.0 #it{p}_{T}", "");
   DoubleRatioErrors->SetFillColor(0);

   legendDoubleRatioABC->Draw();

   canvasDoubleRatioABC->Print(Form("%s/DoubleRatio_ABC_%s.eps",outputDir.Data(),cent.Data()));
   canvasDoubleRatioABC->Print(Form("%s/DoubleRatio_ABC_%s.C",outputDir.Data(),cent.Data()));





   // ------------------------- Double Ratio Fit --------------------------------

   DrawGammaSetMarker(DoubleRatioPi0Fit, 20, 2., 1, 1);

   TCanvas *canvasDoubleRatioPi0Fit = GetAndSetCanvas("canvasDoubleRatioPi0FitFinal");
   canvasDoubleRatioPi0Fit->SetLogx();

   dummy->GetYaxis()->SetRangeUser(doubleRatio[0], doubleRatio[1]);
   dummy->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
   dummy->DrawCopy();

   DoubleRatioPi0Fit->SetMarkerColor(colorcent);
   DoubleRatioPi0Fit->SetLineColor(colorcent);
   DoubleRatioPi0Fit->SetFillColor(0);


   DoubleRatioPi0Fit = RebinTH1D(DoubleRatioPi0Fit,histoBinningRatio);
   One->Draw("same");
   NLODoubleRatio->Draw("p3lsame");
   DoubleRatioPi0Fit->DrawCopy("same");
   DrawGammaSetMarkerTGraphErr(DoubleRatioPi0FitErrors , 20,1, colorcent, colorcent, 1, kTRUE);
   DoubleRatioPi0FitErrors->Draw("E2same");

   DrawGammaSetMarkerTGraphErr(GraphErrorsDummy , 20,2, colorcent, colorcent, 1, kTRUE);


   TLegend* legendDoubleRatioFit = GetAndSetLegend(0.15,0.75,4);
   legendDoubleRatioFit->AddEntry(GraphErrorsDummy,Form("ALICE, %s Pb-Pb #sqrt{s_{_{NN}}} = 2.76 TeV",GetCentralityStringC(cutSel).Data()),"pf");
   legendDoubleRatioFit->AddEntry(NLODoubleRatio,"NLO prediction: 1 + (N_{coll}#gamma_{direct,pp,NLO}/#gamma_{decay})","l");
   legendDoubleRatioFit->AddEntry((TObject*)0, "for #mu = 0.5 to 2.0 #it{p}_{T}", "");
   DoubleRatioPi0FitErrors->SetFillColor(0);

   legendDoubleRatioFit->Draw();



      // TF1 *constFitPrelim = new TF1("constFitPrelim","[0]",1.,2.4);
   TF1 *constFitNew = new TF1("constFitNew","[0]",1.,2.4);
      // PrelimDoubleRatio->Fit(constFitPrelim,"NRE+","",1.,2.4);
   DoubleRatioPi0Fit->Fit(constFitNew,"NRE+","",1.,2.4);
   


   // constFitPrelim->Draw("same");
      // constFitNew->Draw("same");

   canvasDoubleRatioPi0Fit->Print(Form("%s/DoubleRatioFit_%s.eps",outputDir.Data(),cent.Data()));
   canvasDoubleRatioPi0Fit->Print(Form("%s/DoubleRatioFit_%s.C",outputDir.Data(),cent.Data()));


   // ------------------------- Double Ratio Fit ABC --------------------------------

   DrawGammaSetMarker(DoubleRatioPi0Fit, 20, 2., 1, 1);

   TCanvas *canvasDoubleRatioPi0FitABC = GetAndSetCanvas("canvasDoubleRatioPi0FitFinal");
   canvasDoubleRatioPi0FitABC->SetLogx();

   dummy->GetYaxis()->SetRangeUser(doubleRatio[0], doubleRatio[1]);
   dummy->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
   dummy->DrawCopy();

   DoubleRatioPi0Fit->SetMarkerColor(colorcent);
   DoubleRatioPi0Fit->SetLineColor(colorcent);
   DoubleRatioPi0Fit->SetFillColor(0);


   DoubleRatioPi0Fit = RebinTH1D(DoubleRatioPi0Fit,histoBinningRatio);
   One->Draw("same");
   NLODoubleRatio->Draw("p3lsame");
   DoubleRatioPi0Fit->DrawCopy("same");
   DrawGammaSetMarkerTGraphErr(DoubleRatioPi0FitErrorsA , 20,1, kGreen-8, kGreen-8, 1, kTRUE);
   DrawGammaSetMarkerTGraphErr(DoubleRatioPi0FitErrorsB , 20,1, kRed-8, kRed-8, 1, kTRUE);
   DrawGammaSetMarkerTGraphErr(DoubleRatioPi0FitErrorsC , 20,1, kBlue-8, kBlue-8, 1, kTRUE);


   DoubleRatioPi0FitErrorsA->Draw("E2same");
   DoubleRatioPi0FitErrorsB->Draw("E2same");
   DoubleRatioPi0FitErrorsC->Draw("E2same");
   DoubleRatioPi0FitErrors->Draw("E2same");

   DrawGammaSetMarkerTGraphErr(GraphErrorsDummy , 20,2, colorcent, colorcent, 1, kTRUE);


   TLegend* legendDoubleRatioFitABC = GetAndSetLegend(0.15,0.55,7);
   legendDoubleRatioFitABC->AddEntry(GraphErrorsDummy,Form("ALICE, %s Pb-Pb #sqrt{s_{_{NN}}} = 2.76 TeV",GetCentralityStringC(cutSel).Data()),"pf");
   legendDoubleRatioFitABC->AddEntry(DoubleRatioPi0FitErrorsA,"Error Type A","f");
   legendDoubleRatioFitABC->AddEntry(DoubleRatioPi0FitErrorsB,"Error Type B","f");
   legendDoubleRatioFitABC->AddEntry(DoubleRatioPi0FitErrorsC,"Error Type C","f");
   legendDoubleRatioFitABC->AddEntry(NLODoubleRatio,"NLO prediction: 1 + (N_{coll}#gamma_{direct,pp,NLO}/#gamma_{decay})","l");
   legendDoubleRatioFitABC->AddEntry((TObject*)0, "for #mu = 0.5 to 2.0 #it{p}_{T}", "");
   DoubleRatioPi0FitErrors->SetFillColor(0);

   legendDoubleRatioFitABC->Draw();

   canvasDoubleRatioPi0FitABC->Print(Form("%s/DoubleRatioFit_ABC_%s.eps",outputDir.Data(),cent.Data()));
   canvasDoubleRatioPi0FitABC->Print(Form("%s/DoubleRatioFit_ABC_%s.C",outputDir.Data(),cent.Data()));




   TCanvas *canvasDoubleRatioPi0FitPHOS = GetAndSetCanvas("canvasDoubleRatioPi0FitPHOSFinal");

   canvasDoubleRatioPi0FitPHOS->SetLogx();

   dummy->GetYaxis()->SetRangeUser(doubleRatio[0], doubleRatio[1]);
   dummy->GetXaxis()->SetRangeUser(doubleRatioX[0],20);
   dummy->DrawCopy();

   DoubleRatioPi0Fit->SetMarkerColor(colorcent);
   DoubleRatioPi0Fit->SetLineColor(colorcent);
   DoubleRatioPi0Fit->SetFillColor(0);


   DoubleRatioChargedPions->SetMarkerStyle(24);
   DoubleRatioChargedPions->SetMarkerColor(kBlack);
   DoubleRatioChargedPions->SetLineColor(kBlack);
   DoubleRatioChargedPions->SetFillColor(0);

   TF1* ONE2 = new TF1("ooo","1",0,16);
   ONE2->SetLineWidth(1.2);
   ONE2->SetLineColor(kGray+2);

   TH1D *DoubleRatioPi0FitRebined = (TH1D*)DoubleRatioPi0Fit->Clone();
   DoubleRatioPi0FitRebined = RebinTH1D(DoubleRatioPi0FitRebined,histoBinningRatio);
   ONE2->Draw("same");
   //NLODoubleRatio->Draw("p3lsame");

   DrawGammaSetMarkerTGraphErr(DoubleRatioPi0FitErrors , 20,1, kBlack, kBlack, 1, kTRUE);
   DoubleRatioPi0FitRebined->SetMarkerColor(kBlack);
   DoubleRatioPi0FitRebined->SetLineColor(kBlack);
   DoubleRatioPi0FitRebined->DrawCopy("e1same");

   DoubleRatioPHOSSyst->SetLineColor(kRed+2);
   DoubleRatioPHOSStat->SetLineColor(kRed+2);
   DoubleRatioPHOSStat->SetMarkerColor(kRed+2);
   DoubleRatioPHOSStat->SetMarkerSize(2);
   DoubleRatioPHOSSyst->SetMarkerColor(kRed+2);
   DoubleRatioPHOSSyst->SetFillColor(kRed+2);
   DoubleRatioPHOSSyst->SetMarkerSize(2);
   DoubleRatioPi0FitErrors->Draw("E2same");

   DrawGammaSetMarkerTGraphErr(GraphDoubleRatioPHOSSyst , 20,1, kRed+2, kRed+2, 1, kTRUE);
   GraphDoubleRatioPHOSSyst->DrawClone("E2same");
   DoubleRatioPHOSStat->DrawCopy("e1same");

   DrawGammaSetMarkerTGraphErr(GraphErrorsDummy , 20,2, kBlack, kBlack, 1, kTRUE);
   DrawGammaSetMarkerTGraphErr(GraphErrorsDummyB , 20,2, kRed+2, kRed+2, 1, kTRUE);


   TLegend* legendDoubleRatioFitPhos = GetAndSetLegend(0.15,0.68,3);
   legendDoubleRatioFitPhos->AddEntry(GraphErrorsDummy,"PCM" ,"pf");
   legendDoubleRatioFitPhos->AddEntry(GraphErrorsDummyB,"PHOS","pf");

   DoubleRatioPi0FitErrors->SetFillColor(0);

   legendDoubleRatioFitPhos->Draw();

   DrawSystem(0.2,0.86,option.Data(),(GetCentralityString(cutSel)).Data());

   canvasDoubleRatioPi0FitPHOS->Print(Form("%s/DoubleRatioFitPlusPHOS_%s.eps",outputDir.Data(),cent.Data()));

   // ------------------------- Combined Double Ratio --------------------------------

   TCanvas *canvasDoubleRatioPi0FitCombined = GetAndSetCanvas("canvasDoubleRatioPi0FitCombinedFinal");

   canvasDoubleRatioPi0FitCombined->SetLogx();

   dummy->GetYaxis()->SetRangeUser(doubleRatio[0], doubleRatio[1]);
   dummy->GetXaxis()->SetRangeUser(doubleRatioX[0],20);
   dummy->DrawCopy();

   DoubleRatioPi0Fit->SetMarkerColor(colorcent);
   DoubleRatioPi0Fit->SetLineColor(colorcent);
   DoubleRatioPi0Fit->SetFillColor(0);


   DoubleRatioChargedPions->SetMarkerStyle(24);
   DoubleRatioChargedPions->SetMarkerColor(kBlack);
   DoubleRatioChargedPions->SetLineColor(kBlack);
   DoubleRatioChargedPions->SetFillColor(0);



   Double_t newBins[18] = {0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,3.0,3.3,3.7,4.1,4.6,5.4,6.2,7.0};
   Double_t newBinsPHOS[19] = {1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,3.0,3.3,3.7,4.1,4.6,5.4,6.2,7.0,8.0,11.0};
   Double_t newBinsComb[20] = {0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,3.0,3.3,3.7,4.1,4.6,5.4,6.2,7.0,8.0,11.0};
   TH1D *newBinning = new TH1D("","",17,newBins);
   TH1D *newBinningPHOS = new TH1D("","",18,newBinsPHOS);
   TH1D *newBinningComb = new TH1D("","",19,newBinsComb);

   DoubleRatioPi0FitRebined = RebinTH1D(DoubleRatioPi0FitRebined,newBinning);
   ONE2->Draw("same");

   DoubleRatioPHOSSyst = RebinTH1D(DoubleRatioPHOSSyst,newBinningPHOS);
   DoubleRatioPHOSSystA = RebinTH1D(DoubleRatioPHOSSystA,newBinningPHOS);
   DoubleRatioPHOSSystB = RebinTH1D(DoubleRatioPHOSSystB,newBinningPHOS);
   DoubleRatioPHOSSystC = RebinTH1D(DoubleRatioPHOSSystC,newBinningPHOS);
   DoubleRatioPHOSStat = RebinTH1D(DoubleRatioPHOSStat,newBinningPHOS);
   TGraphAsymmErrors *GraphDoubleRatioPHOSSystRebinned = new TGraphAsymmErrors(DoubleRatioPHOSSyst);
   DoubleRatioPi0FitErrorsAsym->RemovePoint(17);
   DoubleRatioPi0FitErrorsAsym->RemovePoint(17);
   DoubleRatioPi0FitErrorsAsym->RemovePoint(17);
   TGraphAsymmErrors *GraphDoubleRatioPHOSSystARebinned = new TGraphAsymmErrors(DoubleRatioPHOSSystA);
   DoubleRatioPi0FitErrorsAAsym->RemovePoint(17);
   DoubleRatioPi0FitErrorsAAsym->RemovePoint(17);
   DoubleRatioPi0FitErrorsAAsym->RemovePoint(17);
   for(Int_t i = 0;i<DoubleRatioPHOSSystB->GetNbinsX();i++){
      DoubleRatioPHOSSystB->SetBinError(i+1,0);
   }
   TGraphAsymmErrors *GraphDoubleRatioPHOSSystBRebinned = new TGraphAsymmErrors(DoubleRatioPHOSSystB);
   DoubleRatioPi0FitErrorsBAsym->RemovePoint(17);
   DoubleRatioPi0FitErrorsBAsym->RemovePoint(17);
   DoubleRatioPi0FitErrorsBAsym->RemovePoint(17);
   TGraphAsymmErrors *GraphDoubleRatioPHOSSystCRebinned = new TGraphAsymmErrors(DoubleRatioPHOSSystC);
   DoubleRatioPi0FitErrorsCAsym->RemovePoint(17);
   DoubleRatioPi0FitErrorsCAsym->RemovePoint(17);
   DoubleRatioPi0FitErrorsCAsym->RemovePoint(17);


   TGraphAsymmErrors *DoubleRatioCombinedSyst;
   TGraphAsymmErrors *DoubleRatioCombinedStat;
   TGraphAsymmErrors *DoubleRatioCombinedTest;

   TGraphAsymmErrors *DoubleRatioCombinedSystA;
   TGraphAsymmErrors *DoubleRatioCombinedStatA;
   TGraphAsymmErrors *DoubleRatioCombinedTestA;

   TGraphAsymmErrors *DoubleRatioCombinedSystB;
   TGraphAsymmErrors *DoubleRatioCombinedStatB;
   TGraphAsymmErrors *DoubleRatioCombinedTestB;

   TGraphAsymmErrors *DoubleRatioCombinedSystC;
   TGraphAsymmErrors *DoubleRatioCombinedStatC;
   TGraphAsymmErrors *DoubleRatioCombinedTestC;


   for(Int_t i = 0; i<GraphDoubleRatioPHOSSystRebinned->GetN();i++){
      GraphDoubleRatioPHOSSystRebinned->SetPointEXhigh(i,newBinningPHOS->GetBinWidth(i+1)*0.5);
      GraphDoubleRatioPHOSSystRebinned->SetPointEXlow(i,newBinningPHOS->GetBinWidth(i+1)*0.5);
   }

   //TCanvas *jjj = new TCanvas();

   //DoubleRatioPi0FitErrorsAsym->DrawClone("apE1");
   // for(Int_t i = 0; i< ErrorDoubleRatio ; i++){
   //    if(!nameCutVariationsDoubleRatio[i].CompareTo("CocktailEta") ||
   //       !nameCutVariationsDoubleRatio[i].CompareTo("CocktailParam") ||
   //       !nameCutVariationsDoubleRatio[i].CompareTo("CocktailEtaNorm") ){
   //       for(Int_t ii = 0; ii<meanErrorsDoubleRatioPi0Fit[i]->GetN()  ; ii++){
   //          Double_t x1,y1,x2,y2,x3,y3;
   //          meanErrorsDoubleRatioPi0Fit[i]->GetPoint(ii,x1,y1);
   //          for(Int_t j = 0;j<DoubleRatioPi0FitErrorsAsym->GetN();j++){
               
   //             DoubleRatioPi0FitErrorsAsym->GetPoint(j,x2,y2);
   //             Double_t xerrlow = DoubleRatioPi0FitErrorsAsym->GetErrorXlow(j);
   //             Double_t xerrhigh = DoubleRatioPi0FitErrorsAsym->GetErrorXhigh(j);
   //             Double_t yerrlow = DoubleRatioPi0FitErrorsAsym->GetErrorYlow(j);
   //             Double_t yerrhigh = DoubleRatioPi0FitErrorsAsym->GetErrorYhigh(j);
   //             if(x1>x2-xerrlow && x1<x2+xerrhigh){
   //                cout<<yerrlow<<endl;
   //                DoubleRatioPi0FitErrorsAsym->SetPointError(j,xerrlow,xerrhigh,y2*sqrt( pow(yerrlow/y2,2) - (y1/100)*(y1/100)),y2*sqrt( pow(yerrhigh/y2,2) - (y1/100)*(y1/100)));
   //                cout<<" B "<<y2*sqrt( pow(yerrlow/y2,2) - (y1/100)*(y1/100))<<endl;
   //             }
   //          }
   //          for(Int_t j = 0;j<GraphDoubleRatioPHOSSystRebinned->GetN();j++){
   //             Double_t x1,y1,x2,y2,x3,y3;
   //             GraphDoubleRatioPHOSSystRebinned->GetPoint(j,x3,y3);
   //             Double_t xerrlowP = GraphDoubleRatioPHOSSystRebinned->GetErrorXlow(j);
   //             Double_t xerrhighP = GraphDoubleRatioPHOSSystRebinned->GetErrorXhigh(j);
   //             Double_t yerrlowP = GraphDoubleRatioPHOSSystRebinned->GetErrorYlow(j);
   //             Double_t yerrhighP = GraphDoubleRatioPHOSSystRebinned->GetErrorYhigh(j);
               
   //             if(x1>x3-xerrlowP && x1<x3+xerrhighP){
   //                GraphDoubleRatioPHOSSystRebinned->SetPointError(j,xerrlowP,xerrhighP,y3*sqrt( pow(yerrlowP/y3,2) - (y1/100)*(y1/100)),y3*sqrt( pow(yerrhighP/y3,2) - (y1/100)*(y1/100)));
   //             }
   //          }
   //       }
   //    }
   // }
   // DoubleRatioPi0FitErrorsAsym->DrawClone("pE1same");

   // return;

   DoubleRatioCombinedTest = CombinePtPointsSpectra(DoubleRatioPi0FitRebined,(TGraphAsymmErrors*)DoubleRatioPi0FitErrorsAsym,
                                                    DoubleRatioPHOSStat,(TGraphAsymmErrors*)GraphDoubleRatioPHOSSystRebinned,
                                                    DoubleRatioCombinedStat,DoubleRatioCombinedSyst,newBinsComb,20,0,0,1);

   DoubleRatioCombinedTestA = CombinePtPointsSpectra(DoubleRatioPi0FitRebined,(TGraphAsymmErrors*)DoubleRatioPi0FitErrorsAAsym,
                                                     DoubleRatioPHOSStat,(TGraphAsymmErrors*)GraphDoubleRatioPHOSSystARebinned,
                                                     DoubleRatioCombinedStatA,DoubleRatioCombinedSystA,newBinsComb,20,0,0,1);

   DoubleRatioCombinedTestB = CombinePtPointsSpectra(DoubleRatioPi0FitRebined,(TGraphAsymmErrors*)DoubleRatioPi0FitErrorsBAsym,
                                                     DoubleRatioPHOSStat,(TGraphAsymmErrors*)GraphDoubleRatioPHOSSystBRebinned,
                                                     DoubleRatioCombinedStatB,DoubleRatioCombinedSystB,newBinsComb,20,0,0,1);

   DoubleRatioCombinedTestC = CombinePtPointsSpectra(DoubleRatioPi0FitRebined,(TGraphAsymmErrors*)DoubleRatioPi0FitErrorsCAsym,
                                                     DoubleRatioPHOSStat,(TGraphAsymmErrors*)GraphDoubleRatioPHOSSystCRebinned,
                                                     DoubleRatioCombinedStatC,DoubleRatioCombinedSystC,newBinsComb,20,0,0,1);
   

   
   // for(Int_t i = 0; i< ErrorDoubleRatio ; i++){
   //    if(//!nameCutVariationsDoubleRatio[i].CompareTo("CocktailEta") ||
   //       //!nameCutVariationsDoubleRatio[i].CompareTo("CocktailParam") ||
   //       !nameCutVariationsDoubleRatio[i].CompareTo("CocktailEtaNorm") ){
   //       for(Int_t ii = 0; ii<meanErrorsDoubleRatioPi0Fit[i]->GetN()  ; ii++){
   //          Double_t x1,y1,x2,y2,x3,y3;
   //          meanErrorsDoubleRatioPi0Fit[i]->GetPoint(ii,x1,y1);
   //          for(Int_t j = 0;j<DoubleRatioCombinedSyst->GetN();j++){
   //             DoubleRatioCombinedSyst->GetPoint(j,x2,y2);
   //             Double_t xerrlow = DoubleRatioCombinedSyst->GetErrorXlow(j);
   //             Double_t xerrhigh = DoubleRatioCombinedSyst->GetErrorXhigh(j);
   //             Double_t yerrlow = DoubleRatioCombinedSyst->GetErrorYlow(j);
   //             Double_t yerrhigh = DoubleRatioCombinedSyst->GetErrorYhigh(j);
   //             if(x1>x2-xerrlow && x1<x2+xerrhigh){
   //                cout<<yerrlow<<endl;
   //                DoubleRatioCombinedSyst->SetPointError(j,xerrlow,xerrhigh,y2*sqrt( pow(yerrlow/y2,2) + (y1/100)*(y1/100)),y2*sqrt( pow(yerrhigh/y2,2) + (y1/100)*(y1/100)));
   //                cout<<" B "<<y2*sqrt( pow(yerrlow/y2,2) - (y1/100)*(y1/100))<<endl;
   //             }
   //          }
   //       }
   //    }
   // }





   DrawGammaSetMarkerTGraphAsym(DoubleRatioPi0FitErrorsAsym , 20,1, kBlack, kBlack,1, kTRUE);
   DrawGammaSetMarkerTGraphAsym(GraphDoubleRatioPHOSSystRebinned , 20,2, kRed+2, kRed+2,1, kTRUE);
   DoubleRatioPi0FitRebined->SetMarkerColor(kBlack);

   DoubleRatioCombinedStat->SetMarkerColor(1);
   DoubleRatioCombinedStat->SetMarkerStyle(20);
   DoubleRatioCombinedStat->SetMarkerSize(2);
   DoubleRatioCombinedSyst->SetMarkerColor(1);

   DoubleRatioCombinedStatA->SetMarkerColor(2);
   DoubleRatioCombinedStatA->SetMarkerStyle(20);
   DoubleRatioCombinedStatA->SetMarkerSize(2);
   DoubleRatioCombinedSystA->SetMarkerColor(2);

    for(Int_t i = 0;i<GraphDoubleRatioPHOSSystRebinned->GetN();i++){
       GraphDoubleRatioPHOSSystRebinned->SetPointError(i,DoubleRatioPHOSStat->GetBinWidth(i+1)*0.5,DoubleRatioPHOSStat->GetBinWidth(i+1)*0.5,GraphDoubleRatioPHOSSystRebinned->GetErrorYlow(i),GraphDoubleRatioPHOSSystRebinned->GetErrorYhigh(i));
    }


   for(Int_t i = 0;i<DoubleRatioCombinedStat->GetN();i++){
      DoubleRatioCombinedSyst->SetPointError(i,newBinningComb->GetBinWidth(i+1)*0.5,newBinningComb->GetBinWidth(i+1)*0.5,DoubleRatioCombinedSyst->GetErrorYlow(i),DoubleRatioCombinedSyst->GetErrorYhigh(i));
      DoubleRatioCombinedStat->SetPointError(i,0,0,DoubleRatioCombinedStat->GetErrorYlow(i),DoubleRatioCombinedStat->GetErrorYhigh(i));
      DoubleRatioCombinedStatA->SetPointError(i,0,0,DoubleRatioCombinedStat->GetErrorYlow(i),DoubleRatioCombinedStat->GetErrorYhigh(i));
      DoubleRatioCombinedSystA->SetPointError(i,newBinningComb->GetBinWidth(i+1)*0.5,newBinningComb->GetBinWidth(i+1)*0.5,DoubleRatioCombinedSystA->GetErrorYlow(i),DoubleRatioCombinedSystA->GetErrorYhigh(i));
      DoubleRatioCombinedStatB->SetPointError(i,0,0,DoubleRatioCombinedStat->GetErrorYlow(i),DoubleRatioCombinedStat->GetErrorYhigh(i));
      DoubleRatioCombinedSystB->SetPointError(i,newBinningComb->GetBinWidth(i+1)*0.5,newBinningComb->GetBinWidth(i+1)*0.5,DoubleRatioCombinedSystB->GetErrorYlow(i),DoubleRatioCombinedSystB->GetErrorYhigh(i));
      DoubleRatioCombinedStatC->SetPointError(i,0,0,DoubleRatioCombinedStat->GetErrorYlow(i),DoubleRatioCombinedStat->GetErrorYhigh(i));
      DoubleRatioCombinedSystC->SetPointError(i,newBinningComb->GetBinWidth(i+1)*0.5,newBinningComb->GetBinWidth(i+1)*0.5,DoubleRatioCombinedSystC->GetErrorYlow(i),DoubleRatioCombinedSystC->GetErrorYhigh(i));
      Double_t x,y;
      DoubleRatioCombinedStat->GetPoint(i,x,y);
      DoubleRatioCombinedStatA->SetPoint(i,x,y);
      DoubleRatioCombinedSystA->SetPoint(i,x,y);
      DoubleRatioCombinedStatB->SetPoint(i,x,y);
      DoubleRatioCombinedSystB->SetPoint(i,x,y);
      DoubleRatioCombinedStatC->SetPoint(i,x,y);
      DoubleRatioCombinedSystC->SetPoint(i,x,y);
   }

   DrawGammaSetMarkerTGraphAsym(DoubleRatioCombinedSyst , 20,1, colorcent, colorcent,1, kTRUE);
   DrawGammaSetMarkerTGraphAsym(DoubleRatioCombinedStat , 20,2, colorcent, colorcent,1, kTRUE);

   DrawGammaSetMarkerTGraphAsym(DoubleRatioCombinedSystA , 20,1, 1, 1, 1, kTRUE);
   DrawGammaSetMarkerTGraphAsym(DoubleRatioCombinedStatA , 20,2, 1, 1, 1, kTRUE);


   DoubleRatioCombinedSyst->Draw("e2psame");
   DoubleRatioCombinedStat->Draw("e1psame");

   // DoubleRatioCombinedSystA->Draw("e2psame");
   // DoubleRatioCombinedStatA->Draw("e1psame");

   //DrawGammaSetMarkerTGraphErr(GraphErrorsDummy , 20,1, 1,colorcent, colorcent, 1, kTRUE);



   TLegend* legendDoubleRatioFitCombined = GetAndSetLegend(0.15,0.78,3);
   legendDoubleRatioFitCombined->AddEntry(GraphErrorsDummy,Form("ALICE, %s Pb-Pb #sqrt{s_{_{NN}}} = 2.76 TeV",GetCentralityStringC(cutSel).Data()),"pf");
   DoubleRatioPi0FitErrors->SetFillColor(0);

   legendDoubleRatioFitCombined->Draw();

   canvasDoubleRatioPi0FitCombined->Print(Form("%s/DoubleRatioFitCombined_%s.eps",outputDir.Data(),cent.Data()));

   // ------------------------- Double Ratio Fit Error --------------------------------

   TCanvas *canvasDoubleRatioPi0WithFitError = GetAndSetCanvas("canvasDoubleRatioPi0WithFitErrorFinal");
   canvasDoubleRatioPi0WithFitError->SetLogx();

   dummy->GetYaxis()->SetRangeUser(doubleRatio[0], doubleRatio[1]);
   dummy->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
   dummy->GetXaxis()->SetRangeUser(doubleRatioX[0],15);
   dummy->DrawCopy();

   One->Draw("same");
   NLODoubleRatio->Draw("p3lsame");
   DoubleRatio->DrawCopy("same");
   DrawGammaSetMarkerTGraphErr(DoubleRatioWithFitErrors , 20,1, colorcent, colorcent, 1, kTRUE);
   DoubleRatioWithFitErrors->Draw("E2same");

   DoubleRatioRebinedFit->SetFillColor(kGray);
   TLegend* legendDoubleRatioWithFitError = GetAndSetLegend(0.15,0.75,4);

   legendDoubleRatioWithFitError->AddEntry(GraphErrorsDummy,Form("ALICE, %s Pb-Pb #sqrt{s_{_{NN}}} = 2.76 TeV",GetCentralityStringC(cutSel).Data()),"pf");
   //legendDoubleRatioWithFitError->AddEntry(DoubleRatioPHOSsyst,Form("PHOS"),"pf");
   //legendDoubleRatioWithFitError->AddEntry(PrelimDoubleRatioSyst,Form("Prelim. ALICE, 0-40% Pb-Pb #sqrt{s_{_{NN}}} = 2.76 TeV"),"pf");

   // legendDoubleRatioWithFitError->AddEntry(DoubleRatioRebinedNorm,Form("Point to Point Ratio"));
   // legendDoubleRatioWithFitError->AddEntry(DoubleRatioChargedPions,Form("Charged Pions"),"pe");
   // legendDoubleRatioWithFitError->AddEntry(NLODoubleRatio,"NLO prediction: 1 + (N_{coll}#gamma_{direct,pp,NLO}/#gamma_{decay})","l");

   legendDoubleRatioWithFitError->AddEntry((TObject*)0, "for #mu = 0.5 to 2.0 #it{p}_{T}", "");


   legendDoubleRatioWithFitError->Draw();

   //DrawSystem(0.3,0.86,option.Data(),(GetCentralityString(cutSel)).Data());

   canvasDoubleRatioPi0WithFitError->Print(Form("%s/DoubleRatioWithFitError_%s.eps",outputDir.Data(),cent.Data()));
   canvasDoubleRatioPi0WithFitError->Print(Form("%s/DoubleRatioWithFitError_%s.C",outputDir.Data(),cent.Data()));



   TFile *PrelimFile = new TFile("DirectPhotons_00-40_QM.root");
   TH1D *PrelimDoubleRatio = (TH1D*)PrelimFile->Get("DoubleRatio0-40_stat");
   TGraphErrors *PrelimDoubleRatioSyst = (TGraphErrors*)PrelimFile->Get("DoubleRatio0-40_syst");



   if(centralityCutNumber.CompareTo("504") == 0){


      TCanvas *canvasDoubleRatioPi0WithFitErrorPlusPrelim = GetAndSetCanvas("canvasDoubleRatioPi0WithFitErrorPlusPrelimFinal");
      canvasDoubleRatioPi0WithFitErrorPlusPrelim->SetLogx();

      dummy->GetYaxis()->SetRangeUser(doubleRatio[0], doubleRatio[1]);
      dummy->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
      //dummy->GetXaxis()->SetRangeUser(doubleRatioX[0],15);
      dummy->GetXaxis()->SetRangeUser(0.8,11.4);
      dummy->DrawCopy();
      One->SetRange(0,16);
      One->Draw("same");
      PrelimDoubleRatioSyst->Draw("samee2");
      NLODoubleRatio->Draw("p3lsame");
      //PrelimDoubleRatio->DrawCopy("e1same");
      TGraphErrors *GraphPrelimDoubleRatio = new TGraphErrors(PrelimDoubleRatio);
      GraphPrelimDoubleRatio->RemovePoint(0);
      GraphPrelimDoubleRatio->Draw("pe1same");
      //DoubleRatio->DrawCopy("same");
      //DrawGammaSetMarkerTGraphErr(DoubleRatioWithFitErrors , 20,1, colorcent, colorcent, 1, kTRUE);
      DrawGammaSetMarkerTGraphErr(DoubleRatioPi0FitErrors , 20,1, colorcent, colorcent, 1, kTRUE);
      // DoubleRatioWithFitErrors->Draw("E2same");
      DoubleRatioPi0Fit->DrawCopy("e1same");
      DoubleRatioPi0FitErrors->Draw("E2same");
      DoubleRatioPi0FitErrors->SetMarkerSize(2);

      // TF1 *constFitPrelim = new TF1("constFitPrelim","[0]",1.,2.4);
      // TF1 *constFitNew = new TF1("constFitNew","[0]",1.,2.4);
      // PrelimDoubleRatio->Fit(constFitPrelim,"NRE+","",1.,2.4);
      // DoubleRatioPi0Fit->Fit(constFitNew,"NRE+","",1.,2.4);
      // constFitPrelim->Draw("same");
      // constFitNew->Draw("same");

      TLegend* legendDoubleRatioWithFitError = GetAndSetLegend(0.15,0.72,4.7);

      legendDoubleRatioWithFitError->AddEntry(DoubleRatioPi0FitErrors,Form("ALICE PCM, %s Pb-Pb #sqrt{s_{_{NN}}} = 2.76 TeV",GetCentralityStringC(cutSel).Data()),"pf");
      legendDoubleRatioWithFitError->AddEntry(PrelimDoubleRatioSyst,Form("Prelim. ALICE, 0-40% Pb-Pb #sqrt{s_{_{NN}}} = 2.76 TeV"),"pfe");
      legendDoubleRatioWithFitError->AddEntry(NLODoubleRatio,"NLO prediction: 1 + (N_{coll}#gamma_{direct,pp,NLO}/#gamma_{decay})","l");
      legendDoubleRatioWithFitError->AddEntry((TObject*)0, "for #mu = 0.5 to 2.0 #it{p}_{T}", "");


      legendDoubleRatioWithFitError->Draw();
      gPad->RedrawAxis();   

      canvasDoubleRatioPi0WithFitErrorPlusPrelim->Print(Form("%s/DoubleRatioWithFitErrorPlusPrelim_%s.eps",outputDir.Data(),cent.Data()));
      canvasDoubleRatioPi0WithFitErrorPlusPrelim->Print(Form("%s/DoubleRatioWithFitErrorPlusPrelim_%s.C",outputDir.Data(),cent.Data()));



   }
   // ------------------------- Double Ratio Fit Error No NLO --------------------------------


   TCanvas *canvasDoubleRatioPi0WithFitErrorNoNLO = GetAndSetCanvas("canvasDoubleRatioPi0WithFitErrorNoNLOFinal");


   dummy->DrawCopy();
   DoubleRatioWithFitErrors->Draw("pe2,same");
   One->Draw("same");
   DoubleRatioRebinedFit->DrawCopy("same");
   //DoubleRatioChargedPions->SetMarkerStyle(20);
   //DoubleRatioChargedPions->SetMarkerSize(2);

   //DoubleRatioChargedPions->DrawCopy("same");

   TLegend* legendDoubleRatioWithFitErrorNoNLO = GetAndSetLegend(0.15,0.65,1.5);
   legendDoubleRatioWithFitErrorNoNLO->AddEntry(GraphErrorsDummy,Form("ALICE, %s Pb-Pb #sqrt{s_{_{NN}}} = 2.76 TeV",GetCentralityStringC(cutSel).Data()),"pf");
   //legendDoubleRatioWithFitErrorNoNLO->AddEntry(DoubleRatioChargedPions,"Direct photon double ratio from charged pions","pf");
   legendDoubleRatioWithFitErrorNoNLO->Draw();

   DrawSystem(0.3,0.86,option.Data(),(GetCentralityStringC(cutSel)).Data());

   canvasDoubleRatioPi0WithFitErrorNoNLO->Print(Form("%s/DoubleRatioWithFitErrorNoNLO_%s.eps",outputDir.Data(),cent.Data()));
   canvasDoubleRatioPi0WithFitErrorNoNLO->Print(Form("%s/DoubleRatioWithFitErrorNoNLO_%s.C",outputDir.Data(),cent.Data()));


   // ------------------------- Inc Ratio Errors --------------------------------

   TCanvas *canvasIncRatioErrors = GetAndSetCanvas("canvasIncRatioErrorsMean");
   TLegend* legendIncRatioErrors = GetAndSetLegend(0.2,0.45,10);

   SetHistogramm(dummy, "#it{p}_{T} (GeV/#it{c})", "Systematic Uncertanties %");
   dummy->GetYaxis()->SetRangeUser(doubleRatioErrors[0],doubleRatioErrors[1]);
   dummy->DrawCopy();
   DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummedIncMatIncRatio, 20, 1.,1,1);
   for(Int_t i = 0; i< ErrorIncRatio; i++){
      if(nameCutVariationsDoubleRatio[i].CompareTo("Fit") != 0){
         cout<<meanErrorsIncRatio[i]<<endl;
         DrawGammaSetMarkerTGraphErr(meanErrorsIncRatio[i], 20+i, 1.,color[i],color[i]);
         //meanErrorsIncRatio[i]->RemovePoint(0);
         meanErrorsIncRatio[i]->Draw("p,csame");
         cout << "trying to create legend" << endl;

         legendIncRatioErrors->AddEntry(meanErrorsIncRatio[i],nameCutVariationsDoubleRatio[i].Data(),"p");
      }
   }
   //meanErrorsCorrSummedIncMatIncRatio->RemovePoint(0);
   meanErrorsCorrSummedIncMatIncRatio->Draw("p,csame");
   legendIncRatioErrors->AddEntry(meanErrorsSummedIncRatio,"quadratically summed plus material","p");
   legendIncRatioErrors->Draw();
   canvasIncRatioErrors->Update();

   canvasIncRatioErrors->Print(Form("%s/IncRatioErrors_%s.eps",outputDir.Data(),cent.Data()));


   // ------------------------- Inc Ratio Fit ErrorsABC --------------------------------
   TCanvas *canvasIncRatioErrorsABC = GetAndSetCanvas("canvasIncRatioErrorsABCMean");
   TLegend* legendIncRatioErrorsABC = GetAndSetLegend(0.2,0.55,4,1,Form("Systematics Double Ratio in %s",(GetCentralityStringC(cutSel)).Data()));

   dummy->GetXaxis()->SetLabelOffset(0);
   SetHistogramm(dummy, "#it{p}_{T} (GeV/#it{c})", "Systematic Uncertanties %");
   dummy->GetYaxis()->SetRangeUser(doubleRatioErrors[0],doubleRatioErrors[1]);
   dummy->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
   dummy->DrawCopy();

   DrawGammaSetMarkerTGraphErr( meanErrorsCorrSummedIncMatIncRatioA, 20,1, kGreen-8, kGreen-8, 1, kTRUE);
   DrawGammaSetMarkerTGraphErr( meanErrorsCorrSummedIncMatIncRatioB, 20,1, kRed-8, kRed-8, 1, kTRUE);
   DrawGammaSetMarkerTGraphErr( meanErrorsCorrSummedIncMatIncRatioC, 20,1, kBlue-8, kBlue-8, 1, kTRUE);


   meanErrorsCorrSummedIncMatIncRatioA->Draw("p,csame");
   meanErrorsCorrSummedIncMatIncRatioB->Draw("p,csame");
   meanErrorsCorrSummedIncMatIncRatioC->Draw("p,csame");

   Double_t *SumIncRatioABC = new Double_t[meanErrorsCorrSummedIncMatIncRatioA->GetN()];

   for(Int_t i = 0;i<meanErrorsCorrSummedIncMatIncRatioA->GetN();i++)
      {
         Double_t x1,x2,x3,y1,y2,y3;
         meanErrorsCorrSummedIncMatIncRatioA->GetPoint(i,x1,y1);
         meanErrorsCorrSummedIncMatIncRatioB->GetPoint(i,x2,y2);
         meanErrorsCorrSummedIncMatIncRatioC->GetPoint(i,x3,y3);
         SumIncRatioABC[i] = sqrt( y1*y1+y2*y2+y3*y3 );
      }

   TGraphErrors *meanErrorsCorrSummedIncMatIncRatioABC = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio , SumIncRatioABC,ptBinsDoubleRatioErr ,errorsMeanErrCorrMatSummedIncRatio );
   DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummedIncMatIncRatioABC, 20, 1.,1,1);
   meanErrorsCorrSummedIncMatIncRatioABC->Draw("p,csame");

   legendIncRatioErrorsABC->AddEntry(meanErrorsCorrSummedIncMatIncRatioA,"Errors Type A","pl");
   legendIncRatioErrorsABC->AddEntry(meanErrorsCorrSummedIncMatIncRatioB,"Errors Type B","pl");
   legendIncRatioErrorsABC->AddEntry(meanErrorsCorrSummedIncMatIncRatioC,"Errors Type C","pl");
   legendIncRatioErrorsABC->AddEntry(meanErrorsCorrSummedIncMatIncRatioABC,"Quad. Sum of Errors Type A,B,C","pl");
   legendIncRatioErrorsABC->Draw();

   canvasIncRatioErrorsABC->Update();

   canvasIncRatioErrorsABC->Print(Form("%s/IncRatioErrorsABC_%s.eps",outputDir.Data(),cent.Data()));


   // ------------------------- Inc Ratio Fit Errors --------------------------------

   TGraphErrors *meanPidEdxTypeIncRatioPi0FitB;
   TGraphErrors *meanQtTypeIncRatioPi0FitB;
   TGraphErrors *meanChi2TypeIncRatioPi0FitB;
   TGraphErrors *meanedEdxTypeIncRatioPi0FitB;

   createdPID = kFALSE;
   TGraphErrors *meanErrorsIncRatioPi0FitPID;
   createdGammaCuts = kFALSE;
   TGraphErrors *meanErrorsIncRatioPi0FitGammaCuts;
   createdTrackCuts = kFALSE;
   TGraphErrors *meanErrorsIncRatioPi0FitTrackCuts;
   createdPi0Extraction = kFALSE;
   TGraphErrors *meanErrorsIncRatioPi0FitPi0Extraction;

   for(Int_t i = 0; i< ErrorIncRatio; i++){
      DrawGammaSetMarkerTGraphErr(meanErrorsIncRatioPi0Fit[i], 20+i, 1.,color[i],color[i]);
      if(!nameCutVariationsIncRatio[i].CompareTo("PidEdx"))
         meanPidEdxTypeIncRatioPi0FitB = (TGraphErrors*) meanErrorsIncRatioPi0Fit[i]->Clone("IncRatioFitError_Syst_TypeB_PidEdx");
      if(!nameCutVariationsIncRatio[i].CompareTo("qT"))
         meanQtTypeIncRatioPi0FitB = (TGraphErrors*) meanErrorsIncRatioPi0Fit[i]->Clone("IncRatioFitError_Syst_TypeB_qT");
      if(!nameCutVariationsIncRatio[i].CompareTo("edEdx"))
         meanedEdxTypeIncRatioPi0FitB = (TGraphErrors*) meanErrorsIncRatioPi0Fit[i]->Clone("IncRatioFitError_Syst_TypeB_edEdx");
      if(!nameCutVariationsIncRatio[i].CompareTo("Chi2"))
         meanChi2TypeIncRatioPi0FitB = (TGraphErrors*) meanErrorsIncRatioPi0Fit[i]->Clone("IncRatioFitError_Syst_TypeB_Chi2");


      if(!nameCutVariationsIncRatio[i].CompareTo("edEdx") || !nameCutVariationsIncRatio[i].CompareTo("PidEdx") || !nameCutVariationsIncRatio[i].CompareTo("TOF")){
         if(!createdPID){
            meanErrorsIncRatioPi0FitPID = (TGraphErrors*)meanErrorsIncRatioPi0Fit[i]->Clone("IncRatioPi0FitError_Syst_PID");
            createdPID = kTRUE;
         }
         else {
            for(Int_t j = 0;j<meanErrorsIncRatioPi0Fit[i]->GetN();j++){
               Double_t x1,x2,y1,y2;
               meanErrorsIncRatioPi0FitPID->GetPoint(j,x1,y1);
               meanErrorsIncRatioPi0Fit[i]->GetPoint(j,x2,y2);
               meanErrorsIncRatioPi0FitPID->SetPoint(j,x1,sqrt(y1*y1+y2*y2));
            }
         }
      }
      if(!nameCutVariationsIncRatio[i].CompareTo("Chi2") || !nameCutVariationsIncRatio[i].CompareTo("PsiPair") || !nameCutVariationsIncRatio[i].CompareTo("qT")){
         if(!createdGammaCuts){
            meanErrorsIncRatioPi0FitGammaCuts = (TGraphErrors*)meanErrorsIncRatioPi0Fit[i]->Clone("IncRatioPi0FitError_Syst_IncRatioPi0FitCuts");
            createdGammaCuts = kTRUE;
         }
         else {
            for(Int_t j = 0;j<meanErrorsIncRatioPi0Fit[i]->GetN();j++){
               Double_t x1,x2,y1,y2;
               meanErrorsIncRatioPi0FitGammaCuts->GetPoint(j,x1,y1);
               meanErrorsIncRatioPi0Fit[i]->GetPoint(j,x2,y2);
               meanErrorsIncRatioPi0FitGammaCuts->SetPoint(j,x1,sqrt(y1*y1+y2*y2));

            }
         }
      }
      if(!nameCutVariationsIncRatio[i].CompareTo("TPCClst") || !nameCutVariationsIncRatio[i].CompareTo("SinglePt")){
         if(!createdTrackCuts){
            meanErrorsIncRatioPi0FitTrackCuts = (TGraphErrors*)meanErrorsIncRatioPi0Fit[i]->Clone("IncRatioPi0FitError_Syst_TrackCuts");
            createdTrackCuts = kTRUE;
         }
         else {
            for(Int_t j = 0;j<meanErrorsIncRatioPi0Fit[i]->GetN();j++){
               Double_t x1,x2,y1,y2;
               meanErrorsIncRatioPi0FitTrackCuts->GetPoint(j,x1,y1);
               meanErrorsIncRatioPi0Fit[i]->GetPoint(j,x2,y2);
               meanErrorsIncRatioPi0FitTrackCuts->SetPoint(j,x1,sqrt(y1*y1+y2*y2));

            }
         }
      }
      if(!nameCutVariationsIncRatio[i].CompareTo("Alpha") || !nameCutVariationsIncRatio[i].CompareTo("IntRange")){
         if(!createdPi0Extraction){
            meanErrorsIncRatioPi0FitPi0Extraction = (TGraphErrors*)meanErrorsIncRatioPi0Fit[i]->Clone("IncRatioPi0FitError_Syst_Pi0Extraction");
            createdPi0Extraction = kTRUE;
         }
         else {
            for(Int_t j = 0;j<meanErrorsIncRatioPi0Fit[i]->GetN();j++){
               Double_t x1,x2,y1,y2;
               meanErrorsIncRatioPi0FitPi0Extraction->GetPoint(j,x1,y1);
               meanErrorsIncRatioPi0Fit[i]->GetPoint(j,x2,y2);
               meanErrorsIncRatioPi0FitPi0Extraction->SetPoint(j,x1,sqrt(y1*y1+y2*y2));

            }
         }
      }
   }



   TCanvas *canvasIncRatioPi0FitErrors  = GetAndSetCanvas("canvasIncRatioPi0FitErrorsMean");
   TLegend* legendIncRatioPi0FitErrors = GetAndSetLegend(0.2,0.45,10);


   dummy->GetYaxis()->SetRangeUser(0,35);
   SetHistogramm(dummy, "#it{p}_{T} (GeV/#it{c})", "Systematic Uncertanties %");
   dummy->GetYaxis()->SetRangeUser(doubleRatioErrors[0],doubleRatioErrors[1]);
   dummy->DrawCopy();
   DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummedIncMatIncRatioPi0Fit, 20, 1.,1,1);
   for(Int_t i = 0; i< ErrorIncRatio ; i++){
      DrawGammaSetMarkerTGraphErr(meanErrorsIncRatioPi0Fit[i], 20+i, 1.,color[i],color[i]);
      //meanErrorsIncRatioPi0Fit[i]->RemovePoint(0);
      meanErrorsIncRatioPi0Fit[i]->Draw("p,csame");
      cout << "trying to create legend" << endl;
      legendIncRatioPi0FitErrors->AddEntry(meanErrorsIncRatioPi0Fit[i],nameCutVariationsDoubleRatio[i].Data(),"p");
   }

   //meanErrorsCorrSummedIncMatIncRatioPi0Fit->RemovePoint(0);
   meanErrorsCorrSummedIncMatIncRatioPi0Fit->Draw("p,csame");
   legendIncRatioPi0FitErrors->AddEntry(meanErrorsSummedIncRatioPi0Fit,"quadratically summed plus material","p");
   legendIncRatioPi0FitErrors->Draw();
   canvasIncRatioPi0FitErrors->Update();


   canvasIncRatioPi0FitErrors->Print(Form("%s/IncRatioFitErrors_%s.eps",outputDir.Data(),cent.Data()));


   // ------------------------- Inc Ratio Fit ErrorsGrouped --------------------------------
   TCanvas *canvasIncRatioPi0FitErrorsGrouped = GetAndSetCanvas("canvasIncRatioPi0FitErrorsGroupedMean");
   TLegend* legendIncRatioPi0FitErrorsGrouped = GetAndSetLegend(0.17,0.58,5,1,Form("Systematics Double Ratio in %s",(GetCentralityStringC(cutSel)).Data()));

   dummy->GetXaxis()->SetLabelOffset(0);
   SetHistogramm(dummy, "#it{p}_{T} (GeV/#it{c})", "Systematic Uncertanties %");
   dummy->GetYaxis()->SetRangeUser(doubleRatioErrors[0],doubleRatioErrors[1]);
   dummy->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
   dummy->DrawCopy();



   DrawGammaSetMarkerTGraphErr(meanErrorsIncRatioPi0FitPID , 20,1, kGreen-8, kGreen-8, 1, kTRUE);
   DrawGammaSetMarkerTGraphErr(meanErrorsIncRatioPi0FitGammaCuts , 20,1, kRed-8, kRed-8, 1, kTRUE);
   DrawGammaSetMarkerTGraphErr(meanErrorsIncRatioPi0FitTrackCuts , 20,1, kBlue-8, kBlue-8, 1, kTRUE);
   DrawGammaSetMarkerTGraphErr(meanErrorsIncRatioPi0FitPi0Extraction , 20,1, kCyan-8, kCyan-8, 1, kTRUE);

   meanErrorsIncRatioPi0FitPID->Draw("p,csame");
   meanErrorsIncRatioPi0FitGammaCuts->Draw("p,csame");
   meanErrorsIncRatioPi0FitTrackCuts->Draw("p,csame");
   meanErrorsIncRatioPi0FitPi0Extraction->Draw("p,csame");
   meanErrorsCorrSummedMaterialGamma->Draw("p,csame");

   Double_t *SumIncRatioPi0FitGrouped = new Double_t[meanErrorsCorrSummedIncMatIncRatioPi0FitA->GetN()];

   for(Int_t i = 0;i<meanErrorsCorrSummedIncMatIncRatioPi0FitA->GetN();i++)
      {
         Double_t x1,x2,x3,x4,y1,y2,y3,y4;
         meanErrorsIncRatioPi0FitPID->GetPoint(i,x1,y1);
         meanErrorsIncRatioPi0FitGammaCuts->GetPoint(i,x2,y2);
         meanErrorsIncRatioPi0FitTrackCuts->GetPoint(i,x3,y3);
         meanErrorsIncRatioPi0FitPi0Extraction->GetPoint(i,x4,y4);
         SumIncRatioPi0FitGrouped[i] = sqrt( y1*y1+y2*y2+y3*y3+y4*y4+4.5*4.5 );

         cout<<x1<<" Gamma "<<y2<<" PID  "<<y1<<" Track "<<y3<<" Pi0 "<<y4<<" Sum "<<sqrt( y1*y1+y2*y2+y3*y3+y4*y4+4.5*4.5 )<<endl;
      }
   //return;


   TGraphErrors *meanErrorsCorrSummedIncMatIncRatioPi0FitGrouped = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio , SumIncRatioPi0FitGrouped,ptBinsDoubleRatioErr ,errorsMeanErrCorrMatSummedIncRatioPi0Fit );
   DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummedIncMatIncRatioPi0FitGrouped, 20, 1.,1,1);
   meanErrorsCorrSummedIncMatIncRatioPi0FitGrouped->Draw("p,csame");

   legendIncRatioPi0FitErrorsGrouped->AddEntry(meanErrorsIncRatioPi0FitPID,"PID Uncertainties","pl");
   legendIncRatioPi0FitErrorsGrouped->AddEntry(meanErrorsIncRatioPi0FitTrackCuts,"Track Selection Uncertainties","pl");
   legendIncRatioPi0FitErrorsGrouped->AddEntry(meanErrorsIncRatioPi0FitGammaCuts,"#gamma Cut Uncertainties","pl");
   legendIncRatioPi0FitErrorsGrouped->AddEntry(meanErrorsIncRatioPi0FitPi0Extraction,"#pi^{0} Signal Extraction","pl");
   legendIncRatioPi0FitErrorsGrouped->AddEntry(meanErrorsCorrSummedMaterialGamma,"Material Budget Uncertainties","pl");
   legendIncRatioPi0FitErrorsGrouped->AddEntry(meanErrorsCorrSummedIncMatIncRatioPi0FitGrouped,"Summed Uncertainty","pl");
   legendIncRatioPi0FitErrorsGrouped->Draw();

   canvasIncRatioPi0FitErrorsGrouped->Update();

   canvasIncRatioPi0FitErrorsGrouped->Print(Form("%s/IncRatioPi0FitErrorsGrouped_%s.eps",outputDir.Data(),cent.Data()));



   // ------------------------- Inc Ratio Fit ErrorsABC --------------------------------
   TCanvas *canvasIncRatioPi0FitErrorsABC = GetAndSetCanvas("canvasIncRatioPi0FitErrorsABCMean");
   TLegend* legendIncRatioPi0FitErrorsABC = GetAndSetLegend(0.2,0.55,4,1,Form("Systematics Double Ratio in %s",(GetCentralityStringC(cutSel)).Data()));

   dummy->GetXaxis()->SetLabelOffset(0);
   SetHistogramm(dummy, "#it{p}_{T} (GeV/#it{c})", "Systematic Uncertanties %");
   dummy->GetYaxis()->SetRangeUser(doubleRatioErrors[0],doubleRatioErrors[1]);
   dummy->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
   dummy->DrawCopy();



   DrawGammaSetMarkerTGraphErr( meanErrorsCorrSummedIncMatIncRatioPi0FitA, 20,1, kGreen-8, kGreen-8, 1, kTRUE);
   DrawGammaSetMarkerTGraphErr( meanErrorsCorrSummedIncMatIncRatioPi0FitB, 20,1, kRed-8, kRed-8, 1, kTRUE);
   DrawGammaSetMarkerTGraphErr( meanErrorsCorrSummedIncMatIncRatioPi0FitC, 20,1, kBlue-8, kBlue-8, 1, kTRUE);

   meanErrorsCorrSummedIncMatIncRatioPi0FitA->Draw("p,csame");
   meanErrorsCorrSummedIncMatIncRatioPi0FitB->Draw("p,csame");
   meanErrorsCorrSummedIncMatIncRatioPi0FitC->Draw("p,csame");

   Double_t *SumIncRatioPi0FitABC = new Double_t[meanErrorsCorrSummedIncMatIncRatioPi0FitA->GetN()];

   for(Int_t i = 0;i<meanErrorsCorrSummedIncMatIncRatioPi0FitA->GetN();i++)
      {
         Double_t x1,x2,x3,y1,y2,y3;
         meanErrorsCorrSummedIncMatIncRatioPi0FitA->GetPoint(i,x1,y1);
         meanErrorsCorrSummedIncMatIncRatioPi0FitB->GetPoint(i,x2,y2);
         meanErrorsCorrSummedIncMatIncRatioPi0FitC->GetPoint(i,x3,y3);
         SumIncRatioPi0FitABC[i] = sqrt( y1*y1+y2*y2+y3*y3 );
      }

   TGraphErrors *meanErrorsCorrSummedIncMatIncRatioPi0FitABC = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio , SumIncRatioPi0FitABC,ptBinsDoubleRatioErr ,errorsMeanErrCorrMatSummedIncRatioPi0Fit );
   DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummedIncMatIncRatioPi0FitABC, 20, 1.,1,1);
   meanErrorsCorrSummedIncMatIncRatioPi0FitABC->Draw("p,csame");

   legendIncRatioPi0FitErrorsABC->AddEntry(meanErrorsCorrSummedIncMatIncRatioPi0FitA,"Errors Type A","pl");
   legendIncRatioPi0FitErrorsABC->AddEntry(meanErrorsCorrSummedIncMatIncRatioPi0FitB,"Errors Type B","pl");
   legendIncRatioPi0FitErrorsABC->AddEntry(meanErrorsCorrSummedIncMatIncRatioPi0FitC,"Errors Type C","pl");
   legendIncRatioPi0FitErrorsABC->AddEntry(meanErrorsCorrSummedIncMatIncRatioPi0FitABC,"Quad. Sum of Errors Type A,B,C","pl");
   legendIncRatioPi0FitErrorsABC->Draw();

   canvasIncRatioPi0FitErrorsABC->Update();

   canvasIncRatioPi0FitErrorsABC->Print(Form("%s/IncRatioPi0FitErrorsABC_%s.eps",outputDir.Data(),cent.Data()));


   // ------------------------- Inc Ratio --------------------------------

   Int_t nbinsIncRatio = meanErrorsCorrSummedIncMatIncRatioPi0Fit->GetN();
   Double_t *xIncRatio = meanErrorsCorrSummedIncMatIncRatioPi0Fit->GetX();
   Double_t *xErrorIncRatio = meanErrorsCorrSummedIncMatIncRatio->GetEX();

   Double_t *yGraphIncRatio = meanErrorsCorrSummedIncMatIncRatio->GetY();
   Double_t *yGraphIncRatioA = meanErrorsCorrSummedIncMatIncRatioA->GetY();
   Double_t *yGraphIncRatioB = meanErrorsCorrSummedIncMatIncRatioB->GetY();
   Double_t *yGraphIncRatioC = meanErrorsCorrSummedIncMatIncRatioC->GetY();
   Double_t *yGraphIncRatioPi0Fit = meanErrorsCorrSummedIncMatIncRatioPi0Fit->GetY();
   Double_t *yGraphIncRatioPi0FitA = meanErrorsCorrSummedIncMatIncRatioPi0FitA->GetY();
   Double_t *yGraphIncRatioPi0FitB = meanErrorsCorrSummedIncMatIncRatioPi0FitB->GetY();
   Double_t *yGraphIncRatioPi0FitPidEdxB = meanPidEdxTypeIncRatioPi0FitB->GetY();
   Double_t *yGraphIncRatioPi0FitqTB = meanQtTypeIncRatioPi0FitB->GetY();
   Double_t *yGraphIncRatioPi0FitChi2B = meanChi2TypeIncRatioPi0FitB->GetY();
   Double_t *yGraphIncRatioPi0FitedEdxB = meanedEdxTypeIncRatioPi0FitB->GetY();
   Double_t *yGraphIncRatioPi0FitC = meanErrorsCorrSummedIncMatIncRatioPi0FitC->GetY();


   Double_t *yIncRatio = new Double_t[nbinsIncRatio];
   Double_t *yIncRatioPi0Fit = new Double_t[nbinsIncRatio];
   Double_t *errorIncRatio = new Double_t[nbinsIncRatio];
   Double_t *errorIncRatioA = new Double_t[nbinsIncRatio];
   Double_t *errorIncRatioB = new Double_t[nbinsIncRatio];
   Double_t *errorIncRatioC = new Double_t[nbinsIncRatio];
   Double_t *errorIncRatioPi0Fit = new Double_t[nbinsIncRatio];
   Double_t *errorIncRatioPi0FitA = new Double_t[nbinsIncRatio];
   Double_t *errorIncRatioPi0FitB = new Double_t[nbinsIncRatio];
   Double_t *errorIncRatioPi0FitPidEdxB = new Double_t[nbinsIncRatio];
   Double_t *errorIncRatioPi0FitqTB = new Double_t[nbinsIncRatio];
   Double_t *errorIncRatioPi0FitedEdxB = new Double_t[nbinsIncRatio];
   Double_t *errorIncRatioPi0FitChi2B = new Double_t[nbinsIncRatio];
   Double_t *errorIncRatioPi0FitC = new Double_t[nbinsIncRatio];

   for(Int_t i = 0; i<nbinsIncRatio; i++){
      yIncRatio[i] = IncRatio->GetBinContent(i+2);
      yIncRatioPi0Fit[i] = IncRatioPi0Fit->GetBinContent(i+2);
      errorIncRatio[i] = yGraphIncRatio[i]*(IncRatio->GetBinContent(i+2)/100);
      errorIncRatioA[i] = yGraphIncRatioA[i]*(IncRatio->GetBinContent(i+2)/100);
      errorIncRatioB[i] = yGraphIncRatioB[i]*(IncRatio->GetBinContent(i+2)/100);
      errorIncRatioC[i] = yGraphIncRatioC[i]*(IncRatio->GetBinContent(i+2)/100);
      errorIncRatioPi0Fit[i] = yGraphIncRatioPi0Fit[i]*(IncRatio->GetBinContent(i+2)/100);
      errorIncRatioPi0FitA[i] = yGraphIncRatioPi0FitA[i]*(IncRatio->GetBinContent(i+2)/100);
      errorIncRatioPi0FitB[i] = yGraphIncRatioPi0FitB[i]*(IncRatio->GetBinContent(i+2)/100);
      errorIncRatioPi0FitPidEdxB[i] = yGraphIncRatioPi0FitPidEdxB[i]*(IncRatio->GetBinContent(i+2)/100);
      errorIncRatioPi0FitqTB[i] = yGraphIncRatioPi0FitqTB[i]*(IncRatio->GetBinContent(i+2)/100);
      errorIncRatioPi0FitChi2B[i] = yGraphIncRatioPi0FitChi2B[i]*(IncRatio->GetBinContent(i+2)/100);
      errorIncRatioPi0FitedEdxB[i] = yGraphIncRatioPi0FitedEdxB[i]*(IncRatio->GetBinContent(i+2)/100);
      errorIncRatioPi0FitC[i] = yGraphIncRatioPi0FitC[i]*(IncRatio->GetBinContent(i+2)/100);


   }


   TGraphErrors *IncRatioErrors = new TGraphErrors(nbinsIncRatio,xIncRatio,yIncRatio,xErrorIncRatio,errorIncRatio);
   TGraphErrors *IncRatioErrorsA = new TGraphErrors(nbinsIncRatio,xIncRatio,yIncRatio,xErrorIncRatio,errorIncRatioA);
   TGraphErrors *IncRatioErrorsB = new TGraphErrors(nbinsIncRatio,xIncRatio,yIncRatio,xErrorIncRatio,errorIncRatioB);
   TGraphErrors *IncRatioErrorsC = new TGraphErrors(nbinsIncRatio,xIncRatio,yIncRatio,xErrorIncRatio,errorIncRatioC);
   TGraphErrors *IncRatioPi0FitErrors = new TGraphErrors(nbinsIncRatio,xIncRatio,yIncRatioPi0Fit,xErrorIncRatio,errorIncRatioPi0Fit);
   TGraphErrors *IncRatioPi0FitErrorsA = new TGraphErrors(nbinsIncRatio,xIncRatio,yIncRatioPi0Fit,xErrorIncRatio,errorIncRatioPi0FitA);
   TGraphErrors *IncRatioPi0FitErrorsB = new TGraphErrors(nbinsIncRatio,xIncRatio,yIncRatioPi0Fit,xErrorIncRatio,errorIncRatioPi0FitB);
   TGraphErrors *IncRatioPi0FitErrorsPidEdxB = new TGraphErrors(nbinsIncRatio,xIncRatio,yIncRatioPi0Fit,xErrorIncRatio,errorIncRatioPi0FitPidEdxB);
   TGraphErrors *IncRatioPi0FitErrorsqTB = new TGraphErrors(nbinsIncRatio,xIncRatio,yIncRatioPi0Fit,xErrorIncRatio,errorIncRatioPi0FitqTB);
   TGraphErrors *IncRatioPi0FitErrorsedEdxB = new TGraphErrors(nbinsIncRatio,xIncRatio,yIncRatioPi0Fit,xErrorIncRatio,errorIncRatioPi0FitedEdxB);
   TGraphErrors *IncRatioPi0FitErrorsChi2B = new TGraphErrors(nbinsIncRatio,xIncRatio,yIncRatioPi0Fit,xErrorIncRatio,errorIncRatioPi0FitChi2B);
   TGraphErrors *IncRatioPi0FitErrorsC = new TGraphErrors(nbinsIncRatio,xIncRatio,yIncRatioPi0Fit,xErrorIncRatio,errorIncRatioPi0FitC);
   TGraphErrors *IncRatioPi0WithFitErrors = new TGraphErrors(nbinsIncRatio,xIncRatio,yIncRatio,xErrorIncRatio,errorIncRatioPi0Fit);

   IncRatioErrors->RemovePoint(IncRatioErrors->GetN()-1);
   IncRatioPi0FitErrors->RemovePoint(IncRatioPi0FitErrors->GetN()-1);
   IncRatioPi0WithFitErrors->RemovePoint(IncRatioPi0WithFitErrors->GetN()-1);
   IncRatioErrors->SetTitle("");
   IncRatioPi0FitErrors->SetTitle("");
   IncRatioPi0WithFitErrors->SetTitle("");
   TCanvas *canvasIncRatio = GetAndSetCanvas("canvasIncRatioFinal");

   DrawGammaSetMarker(IncRatio, 20, 2., 1, 1);

   SetHistogramm(dummy, "#it{p}_{T} (GeV/#it{c})", "#gamma_{inc}/#pi^{0}");
   dummy->GetYaxis()->SetRangeUser(incRatio[0], incRatio[1]);
   dummy->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
   dummy->GetXaxis()->SetLabelOffset(0.015);
   dummy->DrawCopy();

   IncRatio->SetMarkerColor(colorcent);
   IncRatio->SetLineColor(colorcent);
   IncRatio->SetFillColor(0);

   IncRatio->SetTitle("");
   IncRatio = RebinTH1D(IncRatio,histoBinningRatio);
   IncRatio->DrawCopy("same");

   DrawGammaSetMarkerTGraphErr(IncRatioErrors , 20,1, colorcent, colorcent, 1, kTRUE);
   IncRatioErrors->Draw("E2same");

   IncRatio->DrawCopy("same");

   IncRatio->SetFillColor(kGray);

   TLegend* legendIncRatio = GetAndSetLegend(0.12,0.2,1.5);
   legendIncRatio->AddEntry(GraphErrorsDummy,"Inclusive #gamma to #pi^{0} ratio","pf");
   legendIncRatio->Draw();


   DrawSystem(0.15,0.90,option.Data(),(GetCentralityStringC(cutSel)).Data());


   canvasIncRatio->Print(Form("%s/IncRatio_%s.eps",outputDir.Data(),cent.Data()));
   canvasIncRatio->Print(Form("%s/IncRatio_%s.C",outputDir.Data(),cent.Data()));

   // ------------------------- Inc Ratio Fit --------------------------------

   TCanvas *canvasIncRatioPi0Fit = GetAndSetCanvas("canvasIncRatioPi0FitFinal");

   DrawGammaSetMarker(IncRatioPi0Fit, 20, 2., 1, 1);



   IncRatioPi0FitErrors->SetFillColor(kGray);
   dummy->GetYaxis()->SetRangeUser(incRatio[0], incRatio[1]);
   dummy->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
   dummy->DrawCopy();
   IncRatioPi0Fit = RebinTH1D(IncRatioPi0Fit,histoBinningRatio);

   IncRatioPi0Fit->SetMarkerColor(colorcent);
   IncRatioPi0Fit->SetLineColor(colorcent);
   IncRatioPi0Fit->SetFillColor(0);

   DrawGammaSetMarkerTGraphErr(IncRatioPi0FitErrors , 20,1, colorcent, colorcent, 1, kTRUE);
   IncRatioPi0FitErrors->Draw("E2same");
   IncRatioPi0Fit->DrawCopy("same");

   TLegend* legendIncRatioFit = GetAndSetLegend(0.12,0.2,1.5);
   legendIncRatioFit->AddEntry(GraphErrorsDummy,"Inclusive #gamma to #pi^{0} (fit on #pi^{0})","pf");

   legendIncRatioFit->Draw();

   DrawSystem(0.15,0.90,option.Data(),(GetCentralityStringC(cutSel)).Data());

   canvasIncRatioPi0Fit->Print(Form("%s/IncRatioFit_%s.eps",outputDir.Data(),cent.Data()));


   // ------------------------- Inc Ratio with Fit Errors --------------------------------
   TCanvas *canvasIncRatioWithFitErrors = GetAndSetCanvas("canvasIncRatioWithFitErrorsFinal");
   IncRatioPi0WithFitErrors->SetFillColor(kGray);
   dummy->GetYaxis()->SetRangeUser(incRatio[0], incRatio[1]);
   dummy->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
   dummy->DrawCopy();

   DrawGammaSetMarkerTGraphErr( IncRatioPi0WithFitErrors, 20,1, colorcent, colorcent, 1, kTRUE);
   IncRatioPi0WithFitErrors->Draw("E2same");
   IncRatio->DrawCopy("same");

   TLegend* legendIncRatioWithFitError = GetAndSetLegend(0.12,0.2,1.5);
   legendIncRatioWithFitError->AddEntry(GraphErrorsDummy,"Inclusive #gamma to #pi^{0} ratio","pf");

   legendIncRatioWithFitError->Draw();

   DrawSystem(0.15,0.90,option.Data(),(GetCentralityStringC(cutSel)).Data());

   canvasIncRatioWithFitErrors->Print(Form("%s/IncRatioWithFitError_%s.eps",outputDir.Data(),cent.Data()));
   canvasIncRatioWithFitErrors->Print(Form("%s/IncRatioWithFitError_%s.C",outputDir.Data(),cent.Data()));


   // ------------------------- Pi0 Errors --------------------------------

   TCanvas *canvasPi0Errors = GetAndSetCanvas("canvasPi0ErrorsMean");
   TLegend* legendPi0Errors = GetAndSetLegend(0.2,0.45,10);

   SetHistogramm(dummy, "#it{p}_{T} (GeV/#it{c})", "Systematic Uncertanties %");
   dummy->GetYaxis()->SetRangeUser(gammaErrors[0],gammaErrors[1]);
   dummy->DrawCopy();
   DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummedIncMatPi0, 20, 1.,1,1);
   for(Int_t i = 0; i< ErrorPi0; i++){
      if(nameCutVariationsDoubleRatio[i].CompareTo("Fit") != 0){
         cout<<meanErrorsPi0[i]<<endl;
         DrawGammaSetMarkerTGraphErr(meanErrorsPi0[i], 20+i, 1.,color[i],color[i]);
         //meanErrorsPi0[i]->RemovePoint(0);
         meanErrorsPi0[i]->Draw("p,csame");
         cout << "trying to create legend" << endl;
         legendPi0Errors->AddEntry(meanErrorsPi0[i],nameCutVariationsDoubleRatio[i].Data(),"p");
      }
   }
   //meanErrorsCorrSummedIncMatPi0->RemovePoint(0);
   meanErrorsCorrSummedIncMatPi0->Draw("p,csame");
   legendPi0Errors->AddEntry(meanErrorsSummedPi0,"quadratically summed plus material","p");
   legendPi0Errors->Draw();
   canvasPi0Errors->Update();

   canvasPi0Errors->Print(Form("%s/Pi0Errors_%s.eps",outputDir.Data(),cent.Data()));


   // ------------------------- Pi0 Fit Errors --------------------------------

   TCanvas *canvasPi0FitErrors  = GetAndSetCanvas("canvasPi0FitErrorsMean");
   TLegend* legendPi0FitErrors = GetAndSetLegend(0.2,0.45,10);

   dummy->GetYaxis()->SetRangeUser(0,35);

   SetHistogramm(dummy, "#it{p}_{T} (GeV/#it{c})", "Systematic Uncertanties %");
   dummy->GetYaxis()->SetRangeUser(gammaErrors[0],gammaErrors[1]);
   dummy->DrawCopy();
   DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummedIncMatPi0Fit, 20, 1.,1,1);
   for(Int_t i = 0; i< ErrorPi0 ; i++){
      DrawGammaSetMarkerTGraphErr(meanErrorsPi0Fit[i], 20+i, 1.,color[i],color[i]);
      //meanErrorsPi0Fit[i]->RemovePoint(0);
      meanErrorsPi0Fit[i]->Draw("p,csame");
      cout << "trying to create legend" << endl;
      legendPi0FitErrors->AddEntry(meanErrorsPi0Fit[i],nameCutVariationsDoubleRatio[i].Data(),"p");
   }

   //meanErrorsCorrSummedIncMatPi0Fit->RemovePoint(0);
   meanErrorsCorrSummedIncMatPi0Fit->Draw("p,csame");
   legendPi0FitErrors->AddEntry(meanErrorsSummedPi0Fit,"quadratically summed plus material","p");
   legendPi0FitErrors->Draw();
   canvasPi0FitErrors->Update();


   canvasPi0FitErrors->Print(Form("%s/Pi0FitErrors_%s.eps",outputDir.Data(),cent.Data()));

   // ------------------------- Pi0 --------------------------------

   Int_t nbinsPi0 = meanErrorsCorrSummedIncMatPi0Fit->GetN();
   Double_t *xPi0 = meanErrorsCorrSummedIncMatPi0Fit->GetX();
   Double_t *xErrorPi0 = meanErrorsCorrSummedIncMatPi0->GetEX();

   Double_t *yGraphPi0 = meanErrorsCorrSummedIncMatPi0->GetY();
   Double_t *yGraphPi0Fit = meanErrorsCorrSummedIncMatPi0Fit->GetY();
   Double_t *yGraphPi0WithoutMat = meanErrorsCorrSummedIncMatPi0WithoutMat->GetY();
   Double_t *yGraphPi0FitWithoutMat = meanErrorsCorrSummedIncMatPi0FitWithoutMat->GetY();

   Double_t *yPi0 = new Double_t[nbinsPi0];
   Double_t *yPi0Fit = new Double_t[nbinsPi0];
   Double_t *errorPi0 = new Double_t[nbinsPi0];
   Double_t *errorPi0WithoutMat = new Double_t[nbinsPi0];
   Double_t *errorPi0Fit = new Double_t[nbinsPi0];
   Double_t *errorPi0FitWithoutMat = new Double_t[nbinsPi0];


   for(Int_t i = 0; i<nbinsPi0; i++){
      yPi0[i] = Pi0->GetBinContent(i+2);
      yPi0Fit[i] = Pi0Fit->GetBinContent(i+2);
      errorPi0[i] = yGraphPi0[i]*(Pi0->GetBinContent(i+2)/100);
      errorPi0Fit[i] = yGraphPi0Fit[i]*(Pi0Fit->GetBinContent(i+2)/100);
      errorPi0WithoutMat[i] = yGraphPi0WithoutMat[i]*(Pi0->GetBinContent(i+2)/100);
      errorPi0FitWithoutMat[i] = yGraphPi0FitWithoutMat[i]*(Pi0Fit->GetBinContent(i+2)/100);
   }


   TGraphErrors *Pi0Errors = new TGraphErrors(nbinsPi0,xPi0,yPi0,xErrorPi0,errorPi0);
   TGraphErrors *Pi0FitErrors = new TGraphErrors(nbinsPi0,xPi0,yPi0Fit,xErrorPi0,errorPi0Fit);
   TGraphErrors *Pi0ErrorsWithoutMat = new TGraphErrors(nbinsPi0,xPi0,yPi0,xErrorPi0,errorPi0WithoutMat);
   TGraphErrors *Pi0FitErrorsWithoutMat = new TGraphErrors(nbinsPi0,xPi0,yPi0Fit,xErrorPi0,errorPi0FitWithoutMat);
   TGraphErrors *Pi0WithFitErrors = new TGraphErrors(nbinsPi0,xPi0,yPi0,xErrorPi0,errorPi0Fit);

   Pi0Errors->RemovePoint(Pi0Errors->GetN()-1);
   Pi0FitErrors->RemovePoint(Pi0FitErrors->GetN()-1);
   Pi0WithFitErrors->RemovePoint(Pi0WithFitErrors->GetN()-1);
   Pi0Errors->SetTitle("");
   Pi0FitErrors->SetTitle("");
   Pi0WithFitErrors->SetTitle("");

   TCanvas *canvasPi0 = GetAndSetCanvas("canvasPi0Final");
   canvasPi0->SetLogy();
   DrawGammaSetMarker(Pi0, 20, 2., 1, 1);

   SetHistogramm(dummy, "#it{p}_{T} (GeV/#it{c})", "#gamma_{inc}/#pi^{0}");

   dummy->GetYaxis()->SetRangeUser(Pi0->GetBinContent(Pi0->GetNbinsX())/10,Pi0->GetBinContent(2)*10);
   dummy->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
   dummy->GetXaxis()->SetLabelOffset(0.015);
   dummy->DrawCopy();

   Pi0->SetMarkerColor(colorcent);
   Pi0->SetLineColor(colorcent);
   Pi0->SetFillColor(0);

   Pi0->SetTitle("");
   Pi0 = RebinTH1D(Pi0,histoBinningRatio);
   Pi0->DrawCopy("same");

   DrawGammaSetMarkerTGraphErr(Pi0Errors , 20,1, colorcent, colorcent, 1, kTRUE);
   Pi0Errors->Draw("E2same");

   Pi0->DrawCopy("same");

   Pi0->SetFillColor(kGray);

   TLegend* legendPi0 = GetAndSetLegend(0.12,0.2,1.5);
   legendPi0->AddEntry(GraphErrorsDummy,"#pi^{0} spectrum","pf");
   legendPi0->Draw();

   DrawSystem(0.15,0.90,option.Data(),(GetCentralityStringC(cutSel)).Data());


   canvasPi0->Print(Form("%s/Pi0_%s.eps",outputDir.Data(),cent.Data()));
   canvasPi0->Print(Form("%s/Pi0_%s.C",outputDir.Data(),cent.Data()));

   // ------------------------- Pi0 Fit --------------------------------

   TCanvas *canvasPi0Fit = GetAndSetCanvas("canvasPi0FitFinal");
   canvasPi0Fit->SetLogy();
   DrawGammaSetMarker(Pi0Fit, 20, 2., 1, 1);



   Pi0FitErrors->SetFillColor(kGray);
   dummy->GetYaxis()->SetRangeUser(Pi0->GetBinContent(Pi0->GetNbinsX())/10,Pi0->GetBinContent(2)*10);
   dummy->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
   dummy->DrawCopy();
   Pi0Fit = RebinTH1D(Pi0Fit,histoBinningRatio);

   Pi0Fit->SetMarkerColor(colorcent);
   Pi0Fit->SetLineColor(colorcent);
   Pi0Fit->SetFillColor(0);

   DrawGammaSetMarkerTGraphErr(Pi0FitErrors , 20,1, colorcent, colorcent, 1, kTRUE);
   Pi0FitErrors->Draw("E2same");
   Pi0Fit->DrawCopy("same");

   TLegend* legendPi0Fit = GetAndSetLegend(0.12,0.2,1.5);
   legendPi0Fit->AddEntry(GraphErrorsDummy,"Inclusive #gamma to #pi^{0} (fit on #pi^{0})","pf");

   legendPi0Fit->Draw();

   DrawSystem(0.15,0.90,option.Data(),(GetCentralityStringC(cutSel)).Data());

   canvasPi0Fit->Print(Form("%s/Pi0Fit_%s.eps",outputDir.Data(),cent.Data()));


   // ------------------------- Pi0 with Fit Errors --------------------------------
   TCanvas *canvasPi0WithFitErrors = GetAndSetCanvas("canvasPi0WithFitErrorsFinal");
   canvasPi0WithFitErrors->SetLogy();
   Pi0WithFitErrors->SetFillColor(kGray);
   dummy->GetYaxis()->SetRangeUser(Pi0->GetBinContent(Pi0->GetNbinsX())/10,Pi0->GetBinContent(2)*10);
   dummy->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
   dummy->DrawCopy();

   DrawGammaSetMarkerTGraphErr( Pi0WithFitErrors, 20,1, colorcent, colorcent, 1, kTRUE);
   Pi0WithFitErrors->Draw("E2same");
   Pi0->DrawCopy("same");

   TLegend* legendPi0WithFitError = GetAndSetLegend(0.12,0.2,1.5);
   legendPi0WithFitError->AddEntry(GraphErrorsDummy,"#pi^{0} spectrum","pf");

   legendPi0WithFitError->Draw();

   DrawSystem(0.15,0.90,option.Data(),(GetCentralityStringC(cutSel)).Data());


   canvasPi0WithFitErrors->Print(Form("%s/Pi0WithFitError_%s.eps",outputDir.Data(),cent.Data()));
   canvasPi0WithFitErrors->Print(Form("%s/Pi0WithFitError_%s.C",outputDir.Data(),cent.Data()));

   //---------------------

   //---------------------

   Int_t nbinsGamma = meanErrorsCorrSummedIncMatGamma->GetN();
   Double_t *xGamma = meanErrorsCorrSummedIncMatGamma->GetX();
   Double_t *xErrorGamma = meanErrorsCorrSummedIncMatGamma->GetEX();
   Double_t *yGraphGamma = meanErrorsCorrSummedIncMatGamma->GetY();
   Double_t *yGraphGammaPtCorr = meanErrorsCorrSummedWithoutMatGamma->GetY();
   Double_t *yGraphGammaPtUnCorr = meanErrorsCorrSummedMaterialGamma->GetY();
   Double_t *yGraphGammaA = meanErrorsCorrSummedIncMatGammaA->GetY();
   Double_t *yGraphGammaB = meanErrorsCorrSummedIncMatGammaB->GetY();
   Double_t *yGraphGammaC = meanErrorsCorrSummedIncMatGammaC->GetY();

   Double_t *yGraphGammaPidEdxB = meanPidEdxTypeGammaB->GetY();
   Double_t *yGraphGammaqTB = meanQtTypeGammaB->GetY();
   Double_t *yGraphGammaChi2B = meanChi2TypeGammaB->GetY();
   Double_t *yGraphGammaedEdxB = meanedEdxTypeGammaB->GetY();

   Double_t *yGamma = new Double_t[nbinsGamma];
   Double_t *errorGamma = new Double_t[nbinsGamma];
   Double_t *errorGammaPtCorr = new Double_t[nbinsGamma];
   Double_t *errorGammaPtUnCorr = new Double_t[nbinsGamma];
   Double_t *errorGammaA = new Double_t[nbinsGamma];
   Double_t *errorGammaB = new Double_t[nbinsGamma];
   Double_t *errorGammaC = new Double_t[nbinsGamma];

   Double_t *errorGammaPidEdxB = new Double_t[nbinsGamma];
   Double_t *errorGammaqTB = new Double_t[nbinsGamma];
   Double_t *errorGammaedEdxB = new Double_t[nbinsGamma];
   Double_t *errorGammaChi2B = new Double_t[nbinsGamma];

   for(Int_t i = 0; i<nbinsGamma; i++){
      yGamma[i] = Gamma->GetBinContent(i+2);
      errorGamma[i] = yGraphGamma[i]*(Gamma->GetBinContent(i+2)/100);
      errorGammaPtCorr[i] = yGraphGammaPtCorr[i]*(Gamma->GetBinContent(i+2)/100);
      errorGammaPtUnCorr[i] = yGraphGammaPtUnCorr[i]*(Gamma->GetBinContent(i+2)/100);
      errorGammaA[i] = yGraphGammaA[i]*(Gamma->GetBinContent(i+2)/100);
      errorGammaB[i] = yGraphGammaB[i]*(Gamma->GetBinContent(i+2)/100);
      errorGammaC[i] = yGraphGammaC[i]*(Gamma->GetBinContent(i+2)/100);

      errorGammaPidEdxB[i] = yGraphGammaPidEdxB[i]*(Gamma->GetBinContent(i+2)/100);
      errorGammaqTB[i] = yGraphGammaqTB[i]*(Gamma->GetBinContent(i+2)/100);
      errorGammaChi2B[i] = yGraphGammaChi2B[i]*(Gamma->GetBinContent(i+2)/100);
      errorGammaedEdxB[i] = yGraphGammaedEdxB[i]*(Gamma->GetBinContent(i+2)/100);
   }


   // ------------------------- Gamma Spectrum --------------------------------
   TCanvas *canvasGammaSpectrum = GetAndSetCanvas("canvasGammaSpectrumFinal");
   canvasGammaSpectrum->SetLogy();
   dummy->GetXaxis()->SetLabelOffset(0.015);
   //canvasGammaSpectrum->SetLogx();
   TGraphErrors *GammaErrors = new TGraphErrors(nbinsGamma,xGamma,yGamma,xErrorGamma,errorGamma);
   TGraphAsymmErrors *GammaAsymErrors = new TGraphAsymmErrors(nbinsGamma,xGamma,yGamma,xErrorGamma,xErrorGamma,errorGamma,errorGamma);
   TGraphAsymmErrors *GammaAsymErrorsA = new TGraphAsymmErrors(nbinsGamma,xGamma,yGamma,xErrorGamma,xErrorGamma,errorGammaA,errorGammaA);
   TGraphErrors *GammaErrorsPtCorr = new TGraphErrors(nbinsGamma,meanErrorsCorrSummedWithoutMatGamma->GetX(),yGamma,meanErrorsCorrSummedWithoutMatGamma->GetEX(),errorGammaPtCorr);
   TGraphErrors *GammaErrorsPtUnCorr = new TGraphErrors(nbinsGamma,meanErrorsCorrSummedMaterialGamma->GetX(),yGamma,meanErrorsCorrSummedMaterialGamma->GetEX(),errorGammaPtUnCorr);
   TGraphErrors *GammaErrorsA = new TGraphErrors(nbinsGamma,xGamma,yGamma,xErrorGamma,errorGammaA);
   TGraphErrors *GammaErrorsB = new TGraphErrors(nbinsGamma,xGamma,yGamma,xErrorGamma,errorGammaB);
   TGraphErrors *GammaErrorsC = new TGraphErrors(nbinsGamma,xGamma,yGamma,xErrorGamma,errorGammaC);

   TGraphErrors *GammaErrorsPidEdxB = new TGraphErrors(nbinsGamma,xGamma,yGamma,xErrorGamma,errorGammaPidEdxB);
   TGraphErrors *GammaErrorsqTB = new TGraphErrors(nbinsGamma,xGamma,yGamma,xErrorGamma,errorGammaqTB);
   TGraphErrors *GammaErrorsedEdxB = new TGraphErrors(nbinsGamma,xGamma,yGamma,xErrorGamma,errorGammaedEdxB);
   TGraphErrors *GammaErrorsChi2B = new TGraphErrors(nbinsGamma,xGamma,yGamma,xErrorGamma,errorGammaChi2B);

   GammaErrors->SetTitle("");
   GammaErrors->SetFillColor(0);
   SetHistogramm(dummy, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}N}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c}^{2})");
   dummy->GetYaxis()->SetRangeUser(Gamma->GetBinContent(Gamma->GetNbinsX())/10,Gamma->GetBinContent(2)*10);
   dummy->GetXaxis()->SetRangeUser(0.0,11.5);
   dummy->Draw();

   SetHistogramm(Gamma, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}N}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c}^{2})");
   DrawGammaSetMarker(Gamma, 20, 2.0, 1, 1);

   Double_t *shiftedGammaY = gammaAssymError->GetY();
   Double_t *GammaX = GammaErrors->GetX();


   for(Int_t i = 0; i<gammaAssymError->GetN(); i++){
      GammaErrors->SetPoint(i,GammaX[i],shiftedGammaY[i]);
      gammaAssymError->SetPointEXhigh(i, 0.0);
      gammaAssymError->SetPointEXlow(i, 0.0);
   }
   gammaAssymError->SetMarkerSize(2);
   gammaAssymError->SetMarkerStyle(20);
   gammaAssymError->SetMarkerColor(colorcent);
   gammaAssymError->SetLineColor(colorcent);
   DrawGammaSetMarkerTGraphErr(GammaErrors , 20,2, colorcent, colorcent, 1, kTRUE);
   GammaErrors->SetFillColor(0);
   GammaErrors->SetMarkerSize(0);
   GammaErrors->SetLineColor(colorcent);

   TFile *FileGammaSpectrum = new TFile("DoubleRatios_PCM_Combined_110314.root","UPDATE");
   // gammaAssymError->Write(Form("InclusiveGamma_PCM_stat_%s",(GetCentralityStringC(cutSel)).Data()),TObject::kOverwrite);
   // GammaErrors->Write(Form("InclusiveGamma_PCM_Allsyst_%s",(GetCentralityStringC(cutSel)).Data()),TObject::kOverwrite);
   // GammaErrorsA->Write(Form("InclusiveGamma_PCM_AllsystTypeA_%s",(GetCentralityStringC(cutSel)).Data()),TObject::kOverwrite);
   // GammaErrorsB->Write(Form("InclusiveGamma_PCM_AllsystTypeB_%s",(GetCentralityStringC(cutSel)).Data()),TObject::kOverwrite);
   // GammaErrorsPidEdxB->Write(Form("InclusiveGamma_PCM_systTypeB_PidEdx_%s",(GetCentralityStringC(cutSel)).Data()),TObject::kOverwrite);
   // GammaErrorsqTB->Write(Form("InclusiveGamma_PCM_systTypeB_qT_%s",(GetCentralityStringC(cutSel)).Data()),TObject::kOverwrite);
   // GammaErrorsC->Write(Form("GammaErrorsC_syst_%s",(GetCentralityStringC(cutSel)).Data()),TObject::kOverwrite);
   // meanErrorsCorrSummedIncMatGammaA->Write(Form("InclusveGammaRelErrors_AllsystTypeA_%s",(GetCentralityStringC(cutSel)).Data()),TObject::kOverwrite);
   // meanErrorsCorrSummedIncMatGammaB->Write(Form("InclusveGammaRelErrors_AllsystTypeB_%s",(GetCentralityStringC(cutSel)).Data()),TObject::kOverwrite);
   // meanErrorsCorrSummedIncMatGammaC->Write(Form("InclusveGammaRelErrors_AllsystTypeC_%s",(GetCentralityStringC(cutSel)).Data()),TObject::kOverwrite);
   // meanPidEdxTypeGammaB->Write(Form("InclusveGammaRelErrors_systTypeB_PidEdx_%s",(GetCentralityStringC(cutSel)).Data()),TObject::kOverwrite);
   // meanQtTypeGammaB->Write(Form("InclusveGammaRelErrors_systTypeB_qT_%s",(GetCentralityStringC(cutSel)).Data()),TObject::kOverwrite);

   DoubleRatioRebinedFit->Write(Form("DoubleRatio_PCM_stat_%s",(GetCentralityStringC(cutSel)).Data()),TObject::kOverwrite);
   DoubleRatioPi0FitErrors->Write(Form("DoubleRatio_PCM_Allsyst_%s",(GetCentralityStringC(cutSel)).Data()),TObject::kOverwrite);
   DoubleRatioPi0FitErrorsA->Write(Form("DoubleRatio_PCM_AllsystTypeA_%s",(GetCentralityStringC(cutSel)).Data()),TObject::kOverwrite);
   DoubleRatioPi0FitErrorsB->Write(Form("DoubleRatio_PCM_AllsystTypeB_%s",(GetCentralityStringC(cutSel)).Data()),TObject::kOverwrite);
   DoubleRatioPi0FitErrorsC->Write(Form("DoubleRatio_PCM_AllsystTypeC_%s",(GetCentralityStringC(cutSel)).Data()),TObject::kOverwrite);
   DoubleRatioPi0FitErrorsPidEdxB->Write(Form("DoubleRatio_PCM_systTypeB_PidEdx_%s",(GetCentralityStringC(cutSel)).Data()),TObject::kOverwrite);
   DoubleRatioPi0FitErrorsqTB->Write(Form("DoubleRatio_PCM_systTypeB_qT_%s",(GetCentralityStringC(cutSel)).Data()),TObject::kOverwrite);
   DoubleRatioPi0FitErrorsMaterialC->Write(Form("DoubleRatio_PCM_systTypeC_Material_%s",(GetCentralityStringC(cutSel)).Data()),TObject::kOverwrite);
   DoubleRatioPi0FitErrorsEtaNormC->Write(Form("DoubleRatio_PCM_systTypeC_EtaNorm_%s",(GetCentralityStringC(cutSel)).Data()),TObject::kOverwrite);
   meanErrorsCorrSummedIncMatDoubleRatioPi0FitA->Write(Form("DoubleRatioRelErrors_AllsystTypeA_%s",(GetCentralityStringC(cutSel)).Data()),TObject::kOverwrite);
   meanErrorsCorrSummedIncMatDoubleRatioPi0FitB->Write(Form("DoubleRatioRelErrors_AllsystTypeB_%s",(GetCentralityStringC(cutSel)).Data()),TObject::kOverwrite);
   meanErrorsCorrSummedIncMatDoubleRatioPi0FitC->Write(Form("DoubleRatioRelErrors_AllsystTypeC_%s",(GetCentralityStringC(cutSel)).Data()),TObject::kOverwrite);   
   meanPidEdxTypeDoubleRatioPi0FitB->Write(Form("DoubleRatioRelErrors_systTypeB_PidEdx_%s",(GetCentralityStringC(cutSel)).Data()),TObject::kOverwrite);
   meanQtTypeDoubleRatioPi0FitB->Write(Form("DoubleRatioRelErrors_systTypeB_qT_%s",(GetCentralityStringC(cutSel)).Data()),TObject::kOverwrite);
   meanErrorsDoubleRatioPi0FitCocktailEtaNorm->Write(Form("DoubleRatioRelErrors_systTypeC_EtaNorm_%s",(GetCentralityStringC(cutSel)).Data()),TObject::kOverwrite);
   meanErrorsCorrSummedMaterialGamma->Write(Form("DoubleRatioRelErrors_systTypeC_Material_%s",(GetCentralityStringC(cutSel)).Data()),TObject::kOverwrite);

   NLODoubleRatio->Write(Form("NLODoubleRatio_%s",(GetCentralityStringC(cutSel)).Data()),TObject::kOverwrite);
   DoubleRatioCombinedSyst->Write(Form("DoubleRatioCombined_Allsyst_%s",(GetCentralityStringC(cutSel)).Data()),TObject::kOverwrite);
   DoubleRatioCombinedStat->Write(Form("DoubleRatioCombined_stat_%s",(GetCentralityStringC(cutSel)).Data()),TObject::kOverwrite);
   DoubleRatioCombinedSystA->Write(Form("DoubleRatioCombined_AllsystTypeA_%s",(GetCentralityStringC(cutSel)).Data()),TObject::kOverwrite);
   DoubleRatioCombinedSystB->Write(Form("DoubleRatioCombined_AllsystTypeB_%s",(GetCentralityStringC(cutSel)).Data()),TObject::kOverwrite);
   DoubleRatioCombinedSystC->Write(Form("DoubleRatioCombined_AllsystTypeC_%s",(GetCentralityStringC(cutSel)).Data()),TObject::kOverwrite);


   // IncRatioPi0Fit->Write(Form("IncRatio_stat_%s",(GetCentralityStringC(cutSel)).Data()),TObject::kOverwrite);
   // IncRatioPi0FitErrors->Write(Form("IncRatio_syst_%s",(GetCentralityStringC(cutSel)).Data()),TObject::kOverwrite);
   // IncRatioPi0FitErrorsA->Write(Form("IncRatioErrorsA_syst_%s",(GetCentralityStringC(cutSel)).Data()),TObject::kOverwrite);
   // IncRatioPi0FitErrorsB->Write(Form("IncRatioErrorsB_syst_%s",(GetCentralityStringC(cutSel)).Data()),TObject::kOverwrite);
   // IncRatioPi0FitErrorsPidEdxB->Write(Form("IncRatioErrorB_PidEdx_syst_%s",(GetCentralityStringC(cutSel)).Data()),TObject::kOverwrite);
   // IncRatioPi0FitErrorsqTB->Write(Form("IncRatioErrorB_qT_syst_%s",(GetCentralityStringC(cutSel)).Data()),TObject::kOverwrite);
   // IncRatioPi0FitErrorsChi2B->Write(Form("IncRatioErrorB_Chi2_syst_%s",(GetCentralityStringC(cutSel)).Data()),TObject::kOverwrite);
   // IncRatioPi0FitErrorsedEdxB->Write(Form("IncRatioErrorB_edEdx_syst_%s",(GetCentralityStringC(cutSel)).Data()),TObject::kOverwrite);
   // IncRatioPi0FitErrorsC->Write(Form("IncRatioErrorsC_syst_%s",(GetCentralityStringC(cutSel)).Data()),TObject::kOverwrite);
   // meanErrorsCorrSummedIncMatIncRatioPi0FitA->Write(Form("IncRatioErrorsA_percent_syst_%s",(GetCentralityStringC(cutSel)).Data()),TObject::kOverwrite);
   // meanErrorsCorrSummedIncMatIncRatioPi0FitB->Write(Form("IncRatioErrorsB_percent_syst_%s",(GetCentralityStringC(cutSel)).Data()),TObject::kOverwrite);
   // meanPidEdxTypeIncRatioPi0FitB->Write(Form("IncRatioErrorB_percent_PidEdx_syst_%s",(GetCentralityStringC(cutSel)).Data()),TObject::kOverwrite);
   // meanedEdxTypeIncRatioPi0FitB->Write(Form("IncRatioErrorB_percent_edEdx_syst_%s",(GetCentralityStringC(cutSel)).Data()),TObject::kOverwrite);
   // meanChi2TypeIncRatioPi0FitB->Write(Form("IncRatioErrorB_percent_Chi2_syst_%s",(GetCentralityStringC(cutSel)).Data()),TObject::kOverwrite);
   // meanQtTypeIncRatioPi0FitB->Write(Form("IncRatioErrorB_percent_qT_syst_%s",(GetCentralityStringC(cutSel)).Data()),TObject::kOverwrite);
   // meanErrorsCorrSummedIncMatIncRatioPi0FitC->Write(Form("IncRatioErrorsC_percent_syst_%s",(GetCentralityStringC(cutSel)).Data()),TObject::kOverwrite);
   if(GetCentralityStringC(cutSel).CompareTo("0-40%") == 0){
      PrelimDoubleRatioSyst->Write("Prelimniary_QM_DoubleRatio_Syst",TObject::kOverwrite);
      PrelimDoubleRatio->Write("Prelimniary_QM_DoubleRatio_Stat",TObject::kOverwrite);
   }


   // Gamma->SetLineColor(0);
   // Gamma->SetFillColor(kGray);
   // Gamma->SetLineColor(1);


   // if(GetCentralityString(cutSel).CompareTo("0-40%") == 0){
   // GraphErrorsDummy->SetFillColor(kGray);
   // GraphErrorsDummy->SetMarkerColor(2);
   // GraphErrorsDummy->SetMarkerSize(2);
   // GraphErrorsDummy->SetMarkerStyle(20);
   // legendGammaSpectrum->AddEntry(GraphErrorsDummy,"Inclusive #gamma (0-40%)","pf");
   // perGamma->SetMarkerSize(2);
   // perGamma->SetMarkerColor(1);
   // perGamma->SetFillColor(kGray);
   // legendGammaSpectrum->AddEntry(perGamma,"Inclusive #gamma (40-80%)","pf");
   // }
   // else legendGammaSpectrum->AddEntry(GraphErrorsDummy,Form("Inclusive #gamma Spectrum for %s",(GetCentralityString(cutSel)).Data()),"pf");




   // ------------------------- Gamma Spectrum Single --------------------------------
   TCanvas *canvasGammaSpectrumSingle = GetAndSetCanvas("canvasGammaSpectrumSingleFinal");
   canvasGammaSpectrumSingle->SetLogy();
   //canvasGammaSpectrumSingle->SetLogx();
   dummy->GetXaxis()->SetRangeUser(0.3,11.4);
   dummy->GetYaxis()->SetRangeUser(Gamma->GetBinContent(Gamma->GetNbinsX())/10,Gamma->GetBinContent(2)*10);
   dummy->DrawCopy();



   DrawGammaSetMarkerTGraphErr(GammaErrors , 20,2, colorcent, colorcent, 1, kTRUE);
   GammaErrors->Draw("E2same");
   gammaAssymError->Draw("pZsame");

   Gamma->SetMarkerColor(1);
   Gamma->SetMarkerSize(2);
   Gamma->SetLineColor(1);
   TLegend* legendGammaSpectrumSingle = GetAndSetLegend(0.15,0.2,2.);
   legendGammaSpectrumSingle->AddEntry(GraphErrorsDummy,"Inclusive #gamma","pf");

   legendGammaSpectrumSingle->Draw();

   canvasGammaSpectrumSingle->Print(Form("%s/GammaSpectrumSingle_%s.eps",outputDir.Data(),cent.Data()));

   TCanvas *canvasGammaSpectrumCombined = GetAndSetCanvas("canvasGammaSpectrumCombinedFinal");
   canvasGammaSpectrumCombined->SetLogy();
   //canvasGammaSpectrumCombined->SetLogx();
   dummy->GetXaxis()->SetRangeUser(0.3,15.2);
   dummy->GetYaxis()->SetRangeUser(1e-8,100);
   dummy->DrawCopy();


   TH1D *GammaRebined = (TH1D*)Gamma->Clone();
   GammaRebined= RebinTH1D(GammaRebined,newBinning);

   PhosSystGamma = RebinTH1D(PhosSystGamma,newBinningPHOS);
   PhosSystGammaAB = RebinTH1D(PhosSystGammaAB,newBinningPHOS);
   PhosStatGamma = RebinTH1D(PhosStatGamma,newBinningPHOS);
   TGraphAsymmErrors *GraphGammaPHOSSystRebinned = new TGraphAsymmErrors(PhosSystGamma);
   GammaAsymErrors->RemovePoint(17);
   GammaAsymErrors->RemovePoint(17);
   GammaAsymErrors->RemovePoint(17);

   TGraphAsymmErrors *GammaCombinedSyst;
   TGraphAsymmErrors *GammaCombinedStat;
   TGraphAsymmErrors *GammaCombinedTest;


   GammaCombinedTest = CombinePtPointsSpectra(GammaRebined,(TGraphAsymmErrors*)GammaAsymErrors,
                                              PhosStatGamma,(TGraphAsymmErrors*)GraphGammaPHOSSystRebinned,
                                              GammaCombinedStat,GammaCombinedSyst,newBinsComb,20,0,0,1);

   TGraphAsymmErrors *GraphGammaPHOSSystABRebinned = new TGraphAsymmErrors(PhosSystGammaAB);
   GammaAsymErrorsA->RemovePoint(17);
   GammaAsymErrorsA->RemovePoint(17);
   GammaAsymErrorsA->RemovePoint(17);

   TGraphAsymmErrors *GammaCombinedSystAB;
   TGraphAsymmErrors *GammaCombinedStatAB;
   TGraphAsymmErrors *GammaCombinedTestAB;


   GammaCombinedTestAB = CombinePtPointsSpectra(GammaRebined,(TGraphAsymmErrors*)GammaAsymErrorsA,
                                              PhosStatGamma,(TGraphAsymmErrors*)GraphGammaPHOSSystABRebinned,
                                              GammaCombinedStatAB,GammaCombinedSystAB,newBinsComb,20,0,0,1);

   for(Int_t i = 0;i<GammaCombinedSystAB->GetN();i++){
      GammaCombinedSystAB->SetPoint(i,GammaCombinedTest->GetX()[i],GammaCombinedTest->GetY()[i]);
      GammaCombinedStatAB->SetPoint(i,GammaCombinedTest->GetX()[i],GammaCombinedTest->GetY()[i]);
      GammaCombinedTestAB->SetPoint(i,GammaCombinedTest->GetX()[i],GammaCombinedTest->GetY()[i]);
   }


   DrawGammaSetMarkerTGraphAsym(GammaAsymErrors , 20,2, colorcent, colorcent, 1, kTRUE);
   DrawGammaSetMarkerTGraphAsym(GraphGammaPHOSSystRebinned , 20,2, 1, 1, 1, kTRUE);
   DrawGammaSetMarkerTGraphAsym(GammaCombinedTest , 20,1, 3, 3, 1, kTRUE);


   GammaCombinedTest->Draw("pE1same");
   GammaAsymErrors->Draw("E2same");
   GraphGammaPHOSSystRebinned->Draw("E2same");

   Gamma->SetMarkerColor(1);
   Gamma->SetMarkerSize(2);
   Gamma->SetLineColor(1);
   TLegend* legendGammaSpectrumCombined = GetAndSetLegend(0.15,0.2,2.);
   legendGammaSpectrumCombined->AddEntry(GraphErrorsDummy,"Inclusive #gamma","pf");

   legendGammaSpectrumCombined->Draw();

   canvasGammaSpectrumCombined->Print(Form("%s/GammaSpectrumCombined_%s.eps",outputDir.Data(),cent.Data()));


   // TH1D *DoubleRatioPHOSStatPCMBins = RebinTH1D(DoubleRatioPHOSStat,newBinning);
   // TH1D *PhosStatGammaPCMBins = RebinTH1D(PhosStatGamma,newBinning);

   // DoubleRatioPHOSStatPCMBins->Multiply(PhosStatGammaPCMBins);

   // TH1D *DoubleRatioPi0FitRebinedTimesGamma = (TH1D*)DoubleRatioPi0FitRebined->Clone();
   // DoubleRatioPi0FitRebinedTimesGamma->Multiply(GammaRebined);
   // DoubleRatioPHOSStatPCMBins->Divide(DoubleRatioPi0FitRebinedTimesGamma);
   // DoubleRatioPHOSStatPCMBins->Draw();
   // //PhosStatGammaPCMBins->Scale(0.95);
   // PhosStatGammaPCMBins->Divide(GammaRebined);


   // PhosStatGammaPCMBins->Draw();
   // return;

   //Inc Ratio By Hand
   // Pi0WithFitErrors
   // Pi0FitErrors
   // Pi0Errors
   // GammaErrors


   GammaErrors->RemovePoint(GammaErrors->GetN()-1);

   Double_t *gammaToPi0x = new Double_t[GammaErrors->GetN()];
   Double_t *gammaToPi0y = new Double_t[GammaErrors->GetN()];
   Double_t *gammaToPi0errx = new Double_t[GammaErrors->GetN()];
   Double_t *gammaToPi0erry = new Double_t[GammaErrors->GetN()];

   Double_t *gammaToPi0Fitx = new Double_t[GammaErrors->GetN()];
   Double_t *gammaToPi0Fity = new Double_t[GammaErrors->GetN()];
   Double_t *gammaToPi0Fiterrx = new Double_t[GammaErrors->GetN()];
   Double_t *gammaToPi0Fiterry = new Double_t[GammaErrors->GetN()];
   Double_t *gammaToPi0FitCancellationSize = new Double_t[GammaErrors->GetN()];


   for(Int_t i = 0; i<GammaErrors->GetN();i++){
      Double_t x1,x2,y1,y2,yerr1,yerr2;
      GammaErrors->GetPoint(i,x1,y1);
      Pi0Errors->GetPoint(i,x2,y2);
      if(x1!=x2) return;
      gammaToPi0x[i] = x1;
      gammaToPi0y[i] = y1/y2;
      yerr1 = GammaErrors->GetErrorY(i);
      yerr2 = Pi0Errors->GetErrorY(i);
      gammaToPi0errx[i] = GammaErrors->GetErrorX(i);
      gammaToPi0erry[i] = sqrt( (y1*y1*yerr2*yerr2+y2*y2*yerr1*yerr1)/(y2*y2*y2*y2));

      GammaErrors->GetPoint(i,x1,y1);
      Pi0FitErrors->GetPoint(i,x2,y2);
      if(x1!=x2) return;
      gammaToPi0Fitx[i] = x1;
      gammaToPi0Fity[i] = y1/y2;
      yerr1 = GammaErrors->GetErrorY(i);
      yerr2 = Pi0FitErrors->GetErrorY(i);
      gammaToPi0Fiterrx[i] = GammaErrors->GetErrorX(i);
      gammaToPi0Fiterry[i] = sqrt( (y1*y1*yerr2*yerr2+y2*y2*yerr1*yerr1)/(y2*y2*y2*y2));

      gammaToPi0FitCancellationSize[i] = IncRatioPi0FitErrors->GetErrorY(i)/gammaToPi0Fiterry[i];

   }

   TGraphErrors *GammaToPi0Ratio = new TGraphErrors(GammaErrors->GetN(),gammaToPi0x,yIncRatio,gammaToPi0errx,gammaToPi0erry);
   TGraphErrors *GammaToPi0FitRatio = new TGraphErrors(GammaErrors->GetN(),gammaToPi0Fitx,yIncRatioPi0Fit,gammaToPi0Fiterrx,gammaToPi0Fiterry);
   TGraph *ErrorCancellationSize = new TGraphErrors(GammaErrors->GetN(),gammaToPi0Fitx,gammaToPi0FitCancellationSize);



   TCanvas *canvasErrorCancellation = GetAndSetCanvas("canvasErrorCancellationFinal");

   SetHistogramm(dummy, "#it{p}_{T} (GeV/#it{c})", "Cancellation Size");
   dummy->GetYaxis()->SetRangeUser(incRatio[0], incRatio[1]);
   dummy->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
   dummy->GetXaxis()->SetLabelOffset(0.015);
   dummy->DrawCopy();

   ErrorCancellationSize->SetLineColor(1);
   ErrorCancellationSize->SetMarkerStyle(20);
   ErrorCancellationSize->SetMarkerSize(2);
   ErrorCancellationSize->Draw("pcsame");

   TLegend* legendErrorCancellation = GetAndSetLegend(0.12,0.15,2.5);
   legendErrorCancellation->AddEntry(ErrorCancellationSize,"Uncertainty Ratio of #gamma to #pi^{0}, Error Cancellation","lp");
   legendErrorCancellation->Draw();

   DrawSystem(0.15,0.90,option.Data(),(GetCentralityStringC(cutSel)).Data());

   canvasErrorCancellation->Print(Form("%s/ErrorCancellation_%s.eps",outputDir.Data(),cent.Data()));
   canvasErrorCancellation->Print(Form("%s/ErrorCancellation_%s.C",outputDir.Data(),cent.Data()));




   // ------------------------- Inc Ratio --------------------------------

   TCanvas *canvasGammaToPi0 = GetAndSetCanvas("canvasGammaToPi0Final");

   SetHistogramm(dummy, "#it{p}_{T} (GeV/#it{c})", "#gamma_{inc}/#pi^{0}");
   dummy->GetYaxis()->SetRangeUser(incRatio[0], incRatio[1]);
   dummy->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
   dummy->GetXaxis()->SetLabelOffset(0.015);
   dummy->DrawCopy();

   DrawGammaSetMarkerTGraphErr(GammaToPi0Ratio , 20,1, 1, 1, 1, kTRUE);
   GammaToPi0Ratio->Draw("E2same");
   IncRatioErrors->Draw("E2same");
   IncRatio->DrawCopy("e1same");

   TLegend* legendGammaToPi0 = GetAndSetLegend(0.12,0.15,2.5);
   legendGammaToPi0->AddEntry(GraphErrorsDummy,"Inclusive #gamma to #pi^{0} ratio","pf");
   legendGammaToPi0->AddEntry(GammaToPi0Ratio,"Uncet. w.o. cancellation","f");
   legendGammaToPi0->Draw();

   DrawSystem(0.15,0.90,option.Data(),(GetCentralityStringC(cutSel)).Data());

   canvasGammaToPi0->Print(Form("%s/GammaToPi0_%s.eps",outputDir.Data(),cent.Data()));
   canvasGammaToPi0->Print(Form("%s/GammaToPi0_%s.C",outputDir.Data(),cent.Data()));

   // ------------------------- Inc Ratio Fit --------------------------------

   TCanvas *canvasGammaToPi0Pi0Fit = GetAndSetCanvas("canvasGammaToPi0Pi0FitFinal");

   dummy->GetYaxis()->SetRangeUser(incRatio[0], incRatio[1]);
   dummy->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
   dummy->DrawCopy();

   DrawGammaSetMarkerTGraphErr(GammaToPi0FitRatio , 20,1, 1, 1, 1, kTRUE);
   IncRatioPi0FitErrors->Draw("E2same");
   GammaToPi0FitRatio->Draw("E2same");
   IncRatioPi0Fit->DrawCopy("e1same");

   TLegend* legendGammaToPi0Fit = GetAndSetLegend(0.12,0.15,2.5);
   legendGammaToPi0Fit->AddEntry(GraphErrorsDummy,"Inclusive #gamma to #pi^{0} (fit on #pi^{0})","pf");
   legendGammaToPi0Fit->AddEntry(GammaToPi0FitRatio,"Uncet. w.o. cancellation","f");
   legendGammaToPi0Fit->Draw();

   DrawSystem(0.15,0.90,option.Data(),(GetCentralityStringC(cutSel)).Data());

   canvasGammaToPi0Pi0Fit->Print(Form("%s/GammaToPi0Fit_%s.eps",outputDir.Data(),cent.Data()));

   // ------------------------- Inc Ratio With Charged --------------------------------

   Double_t *xValueIncRatioCharged = new Double_t[GammaErrors->GetN()];
   Double_t *yValueIncRatioCharged = new Double_t[GammaErrors->GetN()];
   Double_t *xErrorIncRatioCharged = new Double_t[GammaErrors->GetN()];
   Double_t *yErrorIncRatioCharged = new Double_t[GammaErrors->GetN()];

   TGraphErrors *GraphChargedPions = new TGraphErrors(ChargedPionsSyst);
   GraphChargedPions->RemovePoint(0);
   for(Int_t i = 0; i<GammaErrors->GetN();i++){
      xValueIncRatioCharged[i] = IncRatioChargedPions->GetBinCenter(i+2);
      xErrorIncRatioCharged[i] = IncRatioChargedPions->GetBinWidth(i+2)*0.5;
      yValueIncRatioCharged[i] = IncRatioChargedPions->GetBinContent(i+2);
      Double_t x1,x2,y1,y2;
      GammaErrors->GetPoint(i,x1,y1);
      GraphChargedPions->GetPoint(i,x2,y2);
      Double_t yerr1 = GammaErrors->GetErrorY(i);
      Double_t yerr2 = GraphChargedPions->GetErrorY(i);
      yErrorIncRatioCharged[i] = sqrt( (yerr1*yerr1*y2*y2 + yerr2*yerr2*y1*y1)/(y2*y2*y2*y2) );
   }

   TGraphErrors *IncRatioChargedErrors = new TGraphErrors(GammaErrors->GetN(),xValueIncRatioCharged,yValueIncRatioCharged,xErrorIncRatioCharged,yErrorIncRatioCharged);

   Double_t *xValueDoubleRatioCharged = new Double_t[GammaErrors->GetN()];
   Double_t *yValueDoubleRatioCharged = new Double_t[GammaErrors->GetN()];
   Double_t *xErrorDoubleRatioCharged = new Double_t[GammaErrors->GetN()];
   Double_t *yErrorDoubleRatioCharged = new Double_t[GammaErrors->GetN()];

   for(Int_t i = 0;i<CocktailChargedErrors->GetN()-1;i++){
      xValueDoubleRatioCharged[i] = DoubleRatioChargedPions->GetBinCenter(i+2);
      xErrorDoubleRatioCharged[i] = DoubleRatioChargedPions->GetBinWidth(i+2)*0.5;
      yValueDoubleRatioCharged[i] = DoubleRatioChargedPions->GetBinContent(i+2);
      Double_t x1,x2,y1,y2;
      IncRatioChargedErrors->GetPoint(i,x1,y1);
      CocktailChargedErrors->GetPoint(i,x2,y2);
      Double_t yerr1 = IncRatioChargedErrors->GetErrorY(i);
      Double_t yerr2 = CocktailChargedErrors->GetErrorY(i);
      yErrorDoubleRatioCharged[i] = sqrt( (yerr1*yerr1*y2*y2 + yerr2*yerr2*y1*y1)/(y2*y2*y2*y2) );
   }
   TGraphErrors *DoubleRatioChargedErrors = new TGraphErrors(GammaErrors->GetN(),xValueDoubleRatioCharged,yValueDoubleRatioCharged,xErrorDoubleRatioCharged,yErrorDoubleRatioCharged);

   // ------------------------- Double Ratio Fit Plus Charged --------------------------------

   DrawGammaSetMarker(DoubleRatioPi0Fit, 20, 2., 1, 1);

   TCanvas *canvasDoubleRatioPi0FitCharged = GetAndSetCanvas("canvasDoubleRatioPi0FitChargedFinal");
   canvasDoubleRatioPi0FitCharged->SetLogx();

   dummy->GetYaxis()->SetRangeUser(doubleRatio[0], doubleRatio[1]);
   dummy->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
   dummy->DrawCopy();

   DoubleRatioPi0Fit->SetMarkerColor(colorcent);
   DoubleRatioPi0Fit->SetLineColor(colorcent);
   DoubleRatioPi0Fit->SetFillColor(0);


   DoubleRatioChargedPions->SetMarkerStyle(24);
   DoubleRatioChargedPions->SetMarkerColor(kBlack);
   DoubleRatioChargedPions->SetLineColor(kBlack);
   DoubleRatioChargedPions->SetFillColor(0);


   DoubleRatioPi0Fit = RebinTH1D(DoubleRatioPi0Fit,histoBinningRatio);
   One->Draw("same");
   NLODoubleRatio->Draw("p3lsame");
   DoubleRatioPi0Fit->DrawCopy("same");
   DrawGammaSetMarkerTGraphErr(DoubleRatioPi0FitErrors , 20,1, colorcent, colorcent, 1, kTRUE);
   DrawGammaSetMarkerTGraphErr(DoubleRatioChargedErrors , 20,1, 1, 1, 1, kTRUE);

   DoubleRatioChargedPions->SetLineColor(kBlack);
   DoubleRatioPi0FitErrors->Draw("E2same");
   DoubleRatioChargedPions->DrawCopy("same");
   DoubleRatioChargedErrors->Draw("E2same");

   // DoubleRatioChargedOnlyGammaLow->Draw("E2same");
   // DoubleRatioChargedOnlyGammaHigh->Draw("E2same");

   DrawGammaSetMarkerTGraphErr(GraphErrorsDummy , 20,2, colorcent, colorcent, 1, kTRUE);


   TLegend* legendDoubleRatioFitCharged = GetAndSetLegend(0.15,0.68,5.7);
   legendDoubleRatioFitCharged->AddEntry(GraphErrorsDummy,Form("ALICE, %s Pb-Pb #sqrt{s_{_{NN}}} = 2.76 TeV",GetCentralityStringC(cutSel).Data()),"pf");
   legendDoubleRatioFitCharged->AddEntry(DoubleRatioChargedPions,"charged Pions (averaged fits to low / high); erros from #gamma only","pe");
   legendDoubleRatioFitCharged->AddEntry(NLODoubleRatio,"NLO prediction: 1 + (N_{coll}#gamma_{direct,pp,NLO}/#gamma_{decay})","l");
   legendDoubleRatioFitCharged->AddEntry((TObject*)0, "for #mu = 0.5 to 2.0 #it{p}_{T}", "");
   DoubleRatioPi0FitErrors->SetFillColor(0);

   legendDoubleRatioFitCharged->Draw();

   canvasDoubleRatioPi0FitCharged->Print(Form("%s/DoubleRatioFitPlusCharged_%s.eps",outputDir.Data(),cent.Data()));
   canvasDoubleRatioPi0FitCharged->Print(Form("%s/DoubleRatioFitPlusCharged_%s.C",outputDir.Data(),cent.Data()));

   //
   // TGraphErrors *Pi0ErrorsWithoutMat = new TGraphErrors(nbinsPi0,xPi0,yPi0,xErrorPi0,errorPi0WithoutMat);
   // TGraphErrors *Pi0FitErrorsWithoutMat = new TGraphErrors(nbinsPi0,xPi0,yPi0Fit,xErrorPi0,errorPi0FitWithoutMat);


   TCanvas *canvasDoubleRatioRatio = GetAndSetCanvas("canvasDoubleRatioRatioFinal");
   canvasDoubleRatioRatio->SetLogx();

   dummy->GetYaxis()->SetRangeUser(doubleRatio[0], doubleRatio[1]);
   dummy->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);

   dummy->GetYaxis()->SetTitle("DR^{#pi^{#pm}}/DR^{#pi^{0}}");
   dummy->DrawCopy();


   TGraphErrors *GammaToPi0ForComp = (TGraphErrors*)GammaErrorsPtUnCorr->Clone("GammaToPi0ForComp");
   TGraphErrors *GammaToChargedForComp = (TGraphErrors*)GammaErrorsPtUnCorr->Clone("GammaToCharged");
   TGraphErrors *ComparisonChargedNeutral = (TGraphErrors*)GammaErrorsPtUnCorr->Clone("GammaToCharged");
   TGraphErrors *ChargedPionsSystGraph = new TGraphErrors(ChargedPionsSyst);
   TGraphErrors *ChargedPionsStatGraph = new TGraphErrors(ChargedPionsStat);
   TGraphErrors *Pi0FitStatGraph = new TGraphErrors(Pi0Fit);
   ChargedPionsSystGraph->RemovePoint(0);
   ChargedPionsStatGraph->RemovePoint(0);
   Pi0FitStatGraph->RemovePoint(0);


   GammaToPi0ForComp->RemovePoint(GammaToPi0ForComp->GetN()-1);
   GammaToChargedForComp->RemovePoint(GammaToChargedForComp->GetN()-1);
   ComparisonChargedNeutral->RemovePoint(ComparisonChargedNeutral->GetN()-1);

   TGraphErrors *NeutralToChargedPions = (TGraphErrors*)Pi0FitErrors->Clone("GammaToPi0ForComp");

   for(Int_t i = 0;i<Pi0FitErrors->GetN();i++){

      Double_t x1,y1,x2,y2,x3,y3;
      Pi0FitErrors->GetPoint(i,x1,y1);
      Double_t yerr1 = Pi0FitErrors->GetErrorY(i);
      ChargedPionsSystGraph->GetPoint(i,x2,y2);
      Double_t yerr2 = ChargedPionsSystGraph->GetErrorY(i);
      NeutralToChargedPions->SetPoint(i,x1,(y1/CocktailRatioChargedPions->GetBinContent(i+2))/(y2/CocktailRatioPi0->GetBinContent(i+2)));
      NeutralToChargedPions->SetPointError(i,DoubleRatioPi0Fit->GetBinWidth(i+2)*0.5,sqrt( (y1*y1*yerr2*yerr2+y2*y2*yerr1*yerr1)/(y2*y2*y2*y2) ) );
      NeutralToChargedPions->GetPoint(i,x3,y3);
      Double_t yerr3 = NeutralToChargedPions->GetErrorY(i);
      NeutralToChargedPions->SetPointError(i,DoubleRatioPi0Fit->GetBinWidth(i+2)*0.5,sqrt(yerr3*yerr3 + pow(y3*(1-CocktailChargeToNeutral->GetBinContent(i+2)),2)));
      cout<<x1<<"  "<<x2<<"  "<<y1/y2<<endl;
   }




   DrawGammaSetMarkerTGraphErr( NeutralToChargedPions, 20,2, 1, 1, 1, kTRUE);
   One->Draw("same");
   NeutralToChargedPions->Draw("pe2same");


   TLegend* legendDoubleRatioRatio = GetAndSetLegend(0.15,0.68,2);
   legendDoubleRatioRatio->AddEntry(NeutralToChargedPions,"Ratio of #pi^{#pm} DR with #pi^{0} DR","pf");

   legendDoubleRatioRatio->Draw();



   // for(Int_t i = 0;i<GammaToPi0ForComp->GetN();i++){
   //    Double_t x1,y1,x2,y2,x3,y3,x4,y4,x5,y5;
   //    GammaErrorsPtUnCorr->GetPoint(i,x1,y1);
   //    Pi0FitErrorsWithoutMat->GetPoint(i,x2,y2);
   //    Double_t yerr1 = GammaErrorsPtUnCorr->GetErrorY(i);
   //    Double_t yerr2 = Pi0FitErrorsWithoutMat->GetErrorY(i);
   //    Double_t yerrStat2 = Pi0FitStatGraph->GetErrorY(i);

   //    cout<<yerrStat2<<"  "<<yerr2<<endl;
   //    yerr2 = sqrt(yerr2*yerr2+yerrStat2*yerrStat2);
   //    GammaToPi0ForComp->SetPoint(i,x1,y1/y2);
   //    GammaToPi0ForComp->SetPointError(i,DoubleRatioPi0Fit->GetBinWidth(i+2)*0.5,sqrt( (y1*y1*yerr2*yerr2+y2*y2*yerr1*yerr1)/(y2*y2*y2*y2) ));

   //    ChargedPionsSystGraph->GetPoint(i,x3,y3);
   //    Double_t yerr3 = ChargedPionsSystGraph->GetErrorY(i);
   //    Double_t yerrStat3 = ChargedPionsStatGraph->GetErrorY(i);
   //    yerr3 = sqrt(yerr3*yerr3+yerrStat3*yerrStat3);
   //    GammaToChargedForComp->SetPoint(i,x1,y1/y3);
   //    GammaToChargedForComp->SetPointError(i,DoubleRatioPi0Fit->GetBinWidth(i+2)*.5,sqrt( (y1*y1*yerr3*yerr3+y3*y3*yerr1*yerr1)/(y3*y3*y3*y3) ));

   //    //cout<<sqrt( (y1*y1*yerr3*yerr3+y3*y3*yerr1*yerr1)/(y3*y3*y3*y3) )<< "  "<< sqrt( (y1*y1*yerr2*yerr2+y2*y2*yerr1*yerr1)/(y2*y2*y2*y2) )<<endl;

   //    GammaToPi0ForComp->GetPoint(i,x4,y4);
   //    Double_t yerr4 = GammaToPi0ForComp->GetErrorY(i);
   //    GammaToChargedForComp->GetPoint(i,x5,y5);
   //    Double_t yerr5 = GammaToChargedForComp->GetErrorY(i);
   //    ComparisonChargedNeutral->SetPoint(i,x1,(y4/CocktailRatioPi0->GetBinContent(i+2))/(y5/CocktailRatioChargedPions->GetBinContent(i+2)));
   //    ComparisonChargedNeutral->SetPointError(i,DoubleRatioPi0Fit->GetBinWidth(i+2)*.5, sqrt( (y4*y4*yerr5*yerr5+y5*y5*yerr4*yerr4)/(y5*y5*y5*y5)));
   //    Double_t yerr6 = ComparisonChargedNeutral->GetErrorY(i);
   //    ComparisonChargedNeutral->SetPointError(i,DoubleRatioPi0Fit->GetBinWidth(i+2)*.5, sqrt( (yerr6*yerr6) + pow(1-CocktailChargeToNeutral->GetBinContent(i+2),2)   ) );




   //    //cout<<x1<<"  "<<x2<<" "<<DoubleRatioPi0Fit->GetBinCenter(i+2)<<"   "<<y1/y2<<endl;

   // }

   // DrawGammaSetMarkerTGraphErr(GammaToPi0ForComp , 20,1, colorcent, colorcent, 1, kTRUE);
   // DrawGammaSetMarkerTGraphErr(GammaToChargedForComp , 20,1, 1, 1, 1, kTRUE);
   // DrawGammaSetMarkerTGraphErr(ComparisonChargedNeutral , 20,1, 1, 1, 1, kTRUE);
   // ComparisonChargedNeutral->Draw("ape1");
   // One->Draw("same");

   DrawSystem(0.15,0.90,option.Data(),(GetCentralityStringC(cutSel)).Data());

   canvasDoubleRatioRatio->Print(Form("%s/DoubleRatioRatioCharged_%s.eps",outputDir.Data(),cent.Data()));
   canvasDoubleRatioRatio->Print(Form("%s/DoubleRatioRatioCharged_%s.C",outputDir.Data(),cent.Data()));



   // ------------------------- Direct Gamma Spectrum --------------------------------

   TH1D *histoDirectPhotonSpectrum = (TH1D*)Gamma->Clone("DirectPhotonSpectrum");
   TH1D *histoDirectPhotonSpectrumSummed = (TH1D*)Gamma->Clone("DirectPhotonSpectrumSummed");
   TH1D *histoDirectPhotonSpectrumSyst = (TH1D*)Gamma->Clone("DirectPhotonSpectrumSyst");
   TH1D *histoDirectPhotonSpectrumSystA = (TH1D*)Gamma->Clone("DirectPhotonSpectrumSystA");
   TH1D *histoDirectPhotonSpectrumSystB = (TH1D*)Gamma->Clone("DirectPhotonSpectrumSystB");
   TH1D *histoDirectPhotonSpectrumSystC = (TH1D*)Gamma->Clone("DirectPhotonSpectrumSystC");
   TH1D *histoDirectPhotonSpectrumSystWithoutMat = (TH1D*)Gamma->Clone("DirectPhotonSpectrumSyst");
   TH1D *histoDirectPhotonSpectrumStat = (TH1D*)Gamma->Clone("DirectPhotonSpectrumStat");
   TH1D *histoDirectPhotonSpectrumSubNLO = (TH1D*)Gamma->Clone("histoDirectPhotonSpectrumSubNLO");
   TH1D *histoDirectPhotonSpectrumSubNLOSummend = (TH1D*)Gamma->Clone("histoDirectPhotonSpectrumSubNLOSummend");
   TH1D *histoDirectPhotonSpectrumSubNLOStat = (TH1D*)Gamma->Clone("histoDirectPhotonSpectrumSubNLOStat");
   TH1D *histoDirectPhotonSpectrumSubNLOSyst = (TH1D*)Gamma->Clone("histoDirectPhotonSpectrumSubNLOSyst");

   Double_t *systErrorsDoubleRatio = new Double_t[DoubleRatioWithFitErrors->GetN()];
   systErrorsDoubleRatio =  meanErrorsCorrSummedIncMatDoubleRatioPi0Fit->GetY();
   Double_t *systErrorsDoubleRatioA = new Double_t[DoubleRatioWithFitErrors->GetN()];
   systErrorsDoubleRatioA =  meanErrorsCorrSummedIncMatDoubleRatioPi0FitA->GetY();
   Double_t *systErrorsDoubleRatioB = new Double_t[DoubleRatioWithFitErrors->GetN()];
   systErrorsDoubleRatioB =  meanErrorsCorrSummedIncMatDoubleRatioPi0FitB->GetY();
   Double_t *systErrorsDoubleRatioC = new Double_t[DoubleRatioWithFitErrors->GetN()];
   systErrorsDoubleRatioC =  meanErrorsCorrSummedIncMatDoubleRatioPi0FitC->GetY();
   Double_t *systErrorsDoubleRatioX = new Double_t[DoubleRatioWithFitErrors->GetN()];
   systErrorsDoubleRatioX =  meanErrorsCorrSummedIncMatDoubleRatioPi0Fit->GetX();
   Double_t *systErrorsWithoutMatDoubleRatio = new Double_t[DoubleRatioWithFitErrors->GetN()];
   systErrorsWithoutMatDoubleRatio =  meanErrorsCorrSummedDoubleRatioPi0Fit->GetY();




   //DoubleRatioChargedPions
   //DoubleRatioRebinedFit
   Bool_t isCharged = kFALSE;
   Bool_t drawCharged = kFALSE;
   TH1D* DoubleRatioWithSummedErrors = (TH1D*) DoubleRatioRebinedFit->Clone("DoubleRatioWithSummedErrors");
   TH1D* DoubleRatioWithStatErrors = (TH1D*) DoubleRatioRebinedFit->Clone("DoubleRatioWithStatErrors");
   TH1D* DoubleRatioWithSystErrors = (TH1D*) DoubleRatioRebinedFit->Clone("DoubleRatioWithSystErrors");
   TH1D* DoubleRatioWithSystErrorsA = (TH1D*) DoubleRatioRebinedFit->Clone("DoubleRatioWithSystErrorsA");
   TH1D* DoubleRatioWithSystErrorsB = (TH1D*) DoubleRatioRebinedFit->Clone("DoubleRatioWithSystErrorsB");
   TH1D* DoubleRatioWithSystErrorsC = (TH1D*) DoubleRatioRebinedFit->Clone("DoubleRatioWithSystErrorsC");
   TH1D* DoubleRatioWithSystErrorsWithoutMat = (TH1D*) DoubleRatioRebinedFit->Clone("DoubleRatioWithSystErrors");

   Gamma = RebinTH1D(Gamma,DoubleRatioRebinedFit);
   histoDirectPhotonSpectrum = RebinTH1D(histoDirectPhotonSpectrum,DoubleRatioRebinedFit);
   histoDirectPhotonSpectrumSummed = RebinTH1D(histoDirectPhotonSpectrumSummed,DoubleRatioRebinedFit);
   histoDirectPhotonSpectrumSyst = RebinTH1D(histoDirectPhotonSpectrumSyst,DoubleRatioRebinedFit);
   histoDirectPhotonSpectrumSystA = RebinTH1D(histoDirectPhotonSpectrumSystA,DoubleRatioRebinedFit);
   histoDirectPhotonSpectrumSystB = RebinTH1D(histoDirectPhotonSpectrumSystB,DoubleRatioRebinedFit);
   histoDirectPhotonSpectrumSystC = RebinTH1D(histoDirectPhotonSpectrumSystC,DoubleRatioRebinedFit);
   histoDirectPhotonSpectrumStat = RebinTH1D(histoDirectPhotonSpectrumStat,DoubleRatioRebinedFit);
   DoubleRatioWithSummedErrors = RebinTH1D(DoubleRatioWithSummedErrors,DoubleRatioRebinedFit);
   DoubleRatioWithStatErrors = RebinTH1D(DoubleRatioWithStatErrors,DoubleRatioRebinedFit);
   DoubleRatioWithSystErrors = RebinTH1D(DoubleRatioWithSystErrors,DoubleRatioRebinedFit);


   for(Int_t i = 1; i<DoubleRatioWithSummedErrors->GetNbinsX();i++){
      cout<<systErrorsDoubleRatioX[i-1]<<"  "<<DoubleRatioWithSummedErrors->GetBinCenter(i+1)<<endl;

      Double_t binErrorSummed = sqrt( pow( (DoubleRatioWithSummedErrors->GetBinError(i+1)/DoubleRatioWithSummedErrors->GetBinContent(i+1))*100,2) + pow(systErrorsDoubleRatio[i-1],2) );
      Double_t binErrorSyst = systErrorsDoubleRatio[i-1];
      Double_t binErrorSystA = systErrorsDoubleRatioA[i-1];
      Double_t binErrorSystB = systErrorsDoubleRatioB[i-1];
      Double_t binErrorSystC = systErrorsDoubleRatioC[i-1];
      Double_t binErrorSystWitoutMat = systErrorsWithoutMatDoubleRatio[i-1];
      Double_t binErrorStat = (DoubleRatioWithSummedErrors->GetBinError(i+1)/DoubleRatioWithSummedErrors->GetBinContent(i+1))*100;

      DoubleRatioWithSummedErrors->SetBinError(i+1,(binErrorSummed/100)*DoubleRatioWithSummedErrors->GetBinContent(i+1));
      DoubleRatioWithSystErrors->SetBinError(i+1,(binErrorSyst/100)*DoubleRatioWithSystErrors->GetBinContent(i+1));
      DoubleRatioWithSystErrorsA->SetBinError(i+1,(binErrorSystA/100)*DoubleRatioWithSystErrorsA->GetBinContent(i+1));
      DoubleRatioWithSystErrorsB->SetBinError(i+1,(binErrorSystB/100)*DoubleRatioWithSystErrorsB->GetBinContent(i+1));
      DoubleRatioWithSystErrorsC->SetBinError(i+1,(binErrorSystC/100)*DoubleRatioWithSystErrorsC->GetBinContent(i+1));
      DoubleRatioWithSystErrorsWithoutMat->SetBinError(i+1,(binErrorSystWitoutMat/100)*DoubleRatioWithSystErrorsWithoutMat->GetBinContent(i+1));
      DoubleRatioWithStatErrors->SetBinError(i+1,(binErrorStat/100)*DoubleRatioWithStatErrors->GetBinContent(i+1));

      // Fehler im sinne von +- GetBinError()
   }


   TH1D *DoubleRatioCombinedStatErr = (TH1D*)newBinningComb->Clone();
   TH1D *DoubleRatioCombinedSystErr = (TH1D*)newBinningComb->Clone();
   TH1D *DoubleRatioCombinedSystErrA = (TH1D*)newBinningComb->Clone();
   TH1D *DoubleRatioCombinedSumErr = (TH1D*)newBinningComb->Clone();

   for(Int_t i = 0;i<DoubleRatioCombinedSyst->GetN();i++){
      Double_t x,y;
      DoubleRatioCombinedSyst->GetPoint(i,x,y);

      DoubleRatioCombinedSystErr->SetBinContent(i+1,y);
      DoubleRatioCombinedSystErr->SetBinError(i+1,DoubleRatioCombinedSyst->GetErrorYhigh(i));

      DoubleRatioCombinedSystErrA->SetBinContent(i+1,y);
      DoubleRatioCombinedSystErrA->SetBinError(i+1,DoubleRatioCombinedSystA->GetErrorYhigh(i));

      DoubleRatioCombinedStatErr->SetBinContent(i+1,y);
      DoubleRatioCombinedStatErr->SetBinError(i+1,DoubleRatioCombinedStat->GetErrorYhigh(i));

      DoubleRatioCombinedSumErr->SetBinContent(i+1,y);
      DoubleRatioCombinedSumErr->SetBinError(i+1,sqrt(pow(DoubleRatioCombinedSyst->GetErrorYhigh(i),2) + pow(DoubleRatioCombinedStat->GetErrorYhigh(i),2) )  );
   }

   TH1D *histoDirectPhotonSpectrumCombinedSyst = (TH1D*)newBinningComb->Clone("DirectPhotonSpectrumCombinedSyst");
   TH1D *histoDirectPhotonSpectrumCombinedSystA = (TH1D*)newBinningComb->Clone("DirectPhotonSpectrumCombinedSystA");
   TH1D *histoDirectPhotonSpectrumCombinedStat = (TH1D*)newBinningComb->Clone("DirectPhotonSpectrumCombinedStat");
   TH1D *histoDirectPhotonSpectrumCombinedSum = (TH1D*)newBinningComb->Clone("DirectPhotonSpectrumCombinedSum");
   TH1D *histoDirectPhotonSpectrumCombined = (TH1D*)newBinningComb->Clone("DirectPhotonSpectrumCombined");
   TH1D *CombinedGamma = (TH1D*)newBinningComb->Clone();
   TH1D *CombinedGammaSyst = (TH1D*)newBinningComb->Clone();

   for(Int_t i = 0;i<newBinningComb->GetNbinsX();i++){
      Double_t x,y;
      GammaCombinedStat->GetPoint(i,x,y);
      CombinedGamma->SetBinContent(i+1,y);
      CombinedGammaSyst->SetBinContent(i+1,y);
      CombinedGamma->SetBinError(i+1,GammaCombinedStat->GetErrorYhigh(i));
      CombinedGammaSyst->SetBinError(i+1,GammaCombinedSyst->GetErrorYhigh(i));
   }

   for(Int_t i = 0; i<newBinningComb->GetNbinsX();i++){
      histoDirectPhotonSpectrumCombinedSyst->SetBinContent(i+1,-1);
      histoDirectPhotonSpectrumCombinedSystA->SetBinContent(i+1,-1);
      histoDirectPhotonSpectrumCombinedStat->SetBinContent(i+1,-1);
      histoDirectPhotonSpectrumCombinedSum->SetBinContent(i+1,-1);
      histoDirectPhotonSpectrumCombined->SetBinContent(i+1,-1);

      Double_t binContent = -1;
      Double_t binError = -1;

      binContent = 1-(1./DoubleRatioCombinedSumErr->GetBinContent(i+1));
      binError = 1-(1./( DoubleRatioCombinedSumErr->GetBinContent(i+1) - DoubleRatioCombinedSumErr->GetBinError(i+1) ));
      binError = binError*CombinedGamma->GetBinContent(i+1);
      histoDirectPhotonSpectrumCombined->SetBinError(i+1,binError);
      binContent = binContent*CombinedGamma->GetBinContent(i+1);
      histoDirectPhotonSpectrumCombined->SetBinContent(i+1,binContent);

      binContent = 1-(1./( DoubleRatioCombinedSystErr->GetBinContent(i+1) - DoubleRatioCombinedSystErr->GetBinError(i+1) ));
      binContent = binContent*CombinedGamma->GetBinContent(i+1);
      histoDirectPhotonSpectrumCombinedSyst->SetBinContent(i+1,binContent);
      histoDirectPhotonSpectrumCombinedSyst->SetBinError(i+1,0);

      binContent = 1-(1./( DoubleRatioCombinedSystErrA->GetBinContent(i+1) - DoubleRatioCombinedSystErrA->GetBinError(i+1) ));
      binContent = binContent*CombinedGamma->GetBinContent(i+1);
      histoDirectPhotonSpectrumCombinedSystA->SetBinContent(i+1,binContent);
      histoDirectPhotonSpectrumCombinedSystA->SetBinError(i+1,0);

      binContent = 1-(1./( DoubleRatioCombinedStatErr->GetBinContent(i+1) - DoubleRatioCombinedStatErr->GetBinError(i+1) ));
      binContent = binContent*CombinedGamma->GetBinContent(i+1);
      histoDirectPhotonSpectrumCombinedStat->SetBinContent(i+1,binContent);
      histoDirectPhotonSpectrumCombinedStat->SetBinError(i+1,0);

      binContent = 1-(1./( DoubleRatioCombinedSumErr->GetBinContent(i+1) - DoubleRatioCombinedSumErr->GetBinError(i+1) ));
      binContent = binContent*CombinedGamma->GetBinContent(i+1);
      histoDirectPhotonSpectrumCombinedSum->SetBinContent(i+1,binContent);
      histoDirectPhotonSpectrumCombinedSum->SetBinError(i+1,0);
   }


   Double_t *yValueSubtractedNLO = PbPb276CTEQ61EPS09BFG2_sum_pdferr_InvYield->GetY();
   Double_t *xValueSubtractedNLO = PbPb276CTEQ61EPS09BFG2_sum_pdferr_InvYield->GetX();


   for(Int_t i = 1; i<Gamma->GetNbinsX();i++){
      histoDirectPhotonSpectrum->SetBinContent(i+1,-1);

      histoDirectPhotonSpectrumSummed->SetBinContent(i+1,-1);
      histoDirectPhotonSpectrumSyst->SetBinContent(i+1,-1);
      histoDirectPhotonSpectrumSyst->SetBinContent(i+1,-1);
      histoDirectPhotonSpectrumSystA->SetBinContent(i+1,-1);
      histoDirectPhotonSpectrumSystB->SetBinContent(i+1,-1);
      histoDirectPhotonSpectrumSystC->SetBinContent(i+1,-1);
      histoDirectPhotonSpectrumStat->SetBinContent(i+1,-1);

      Double_t binContent = -1;
      Double_t binError = -1;
      binContent = 1-(1./DoubleRatioWithSummedErrors->GetBinContent(i+1));
      binError = 1-(1./( DoubleRatioWithSummedErrors->GetBinContent(i+1) - DoubleRatioWithSummedErrors->GetBinError(i+1) ));
      binError = binError*Gamma->GetBinContent(i+1);
      histoDirectPhotonSpectrum->SetBinError(i+1,binError);
      binContent = binContent*Gamma->GetBinContent(i+1);
      histoDirectPhotonSpectrum->SetBinContent(i+1,binContent);
      histoDirectPhotonSpectrum->SetBinError(i+1,binContent*sqrt( pow(binError/binContent,2) + pow(meanErrorsCorrSummedIncMatGamma->GetY()[i-1]/100,2) ));
      if(i>4){
         histoDirectPhotonSpectrumSubNLO->SetBinContent(i+1,binContent-yValueSubtractedNLO[i-4]);
      }
      else histoDirectPhotonSpectrumSubNLO->SetBinContent(i+1,binContent);
      histoDirectPhotonSpectrumSubNLO->SetBinError(i+1,binError);

      binContent = 1-(1./( DoubleRatioWithSummedErrors->GetBinContent(i+1) - DoubleRatioWithSummedErrors->GetBinError(i+1) ));
      binContent = binContent*Gamma->GetBinContent(i+1);
      histoDirectPhotonSpectrumSummed->SetBinContent(i+1,binContent);
      histoDirectPhotonSpectrumSummed->SetBinError(i+1,0);
      if(i>4){
         histoDirectPhotonSpectrumSubNLOSummend->SetBinContent(i+1,binContent-yValueSubtractedNLO[i-4]);
      }
      else histoDirectPhotonSpectrumSubNLOSummend->SetBinContent(i+1,binContent);
      histoDirectPhotonSpectrumSubNLOSummend->SetBinError(i+1,0);

      binContent = 1-(1./( DoubleRatioWithStatErrors->GetBinContent(i+1) - DoubleRatioWithStatErrors->GetBinError(i+1) ));
      binContent = binContent*Gamma->GetBinContent(i+1);
      histoDirectPhotonSpectrumStat->SetBinContent(i+1,binContent);
      histoDirectPhotonSpectrumStat->SetBinError(i+1,0);
      if(i>4){
         histoDirectPhotonSpectrumSubNLOStat->SetBinContent(i+1,binContent-yValueSubtractedNLO[i-4]);
      }
      else histoDirectPhotonSpectrumSubNLOStat->SetBinContent(i+1,binContent);
      histoDirectPhotonSpectrumSubNLOStat->SetBinError(i+1,0);

      // histoDirectPhotonSpectrumStat->SetBinContent(i+1,binContent);
      // binError = 1-(1./( DoubleRatioWithStatErrors->GetBinContent(i+1) + DoubleRatioWithStatErrors->GetBinError(i+1)));
      // binError = binError*Gamma->GetBinContent(i+1);
      // binError = binError-binContent; ??????
      // histoDirectPhotonSpectrumStat->SetBinError(i+1,binError);

      binContent = 1-(1./( DoubleRatioWithSystErrors->GetBinContent(i+1) - DoubleRatioWithSystErrors->GetBinError(i+1) ));
      binContent = binContent*Gamma->GetBinContent(i+1);
      histoDirectPhotonSpectrumSyst->SetBinContent(i+1,binContent);
      histoDirectPhotonSpectrumSyst->SetBinError(i+1,0);
      if(i>4){
         histoDirectPhotonSpectrumSubNLOSyst->SetBinContent(i+1,binContent-yValueSubtractedNLO[i-4]);
      }
      else histoDirectPhotonSpectrumSubNLOSyst->SetBinContent(i+1,binContent);
      histoDirectPhotonSpectrumSubNLOSyst->SetBinError(i+1,0);

      binContent = 1-(1./( DoubleRatioWithSystErrorsA->GetBinContent(i+1) - DoubleRatioWithSystErrorsA->GetBinError(i+1) ));
      binContent = binContent*Gamma->GetBinContent(i+1);
      histoDirectPhotonSpectrumSystA->SetBinContent(i+1,binContent);
      histoDirectPhotonSpectrumSystA->SetBinError(i+1,0);

      binContent = 1-(1./( DoubleRatioWithSystErrorsB->GetBinContent(i+1) - DoubleRatioWithSystErrorsB->GetBinError(i+1) ));
      binContent = binContent*Gamma->GetBinContent(i+1);
      histoDirectPhotonSpectrumSystB->SetBinContent(i+1,binContent);
      histoDirectPhotonSpectrumSystB->SetBinError(i+1,0);

      binContent = 1-(1./( DoubleRatioWithSystErrorsC->GetBinContent(i+1) - DoubleRatioWithSystErrorsC->GetBinError(i+1) ));
      binContent = binContent*Gamma->GetBinContent(i+1);
      histoDirectPhotonSpectrumSystC->SetBinContent(i+1,binContent);
      histoDirectPhotonSpectrumSystC->SetBinError(i+1,0);

      binContent = 1-(1./( DoubleRatioWithSystErrorsWithoutMat->GetBinContent(i+1) - DoubleRatioWithSystErrorsWithoutMat->GetBinError(i+1) ));
      binContent = binContent*Gamma->GetBinContent(i+1);
      histoDirectPhotonSpectrumSystWithoutMat->SetBinContent(i+1,binContent);
      histoDirectPhotonSpectrumSystWithoutMat->SetBinError(i+1,0);
   }

   TGraphAsymmErrors *grStatFull = getPhotonPoints(histoDirectPhotonSpectrumStat,histoDirectPhotonSpectrum,0.75,0.5);
   if(grStatFull) grStatFull->SetName("graphDataPointsStatError");
   TGraphAsymmErrors *grStatPoints1 = getPhotonPoints(histoDirectPhotonSpectrumStat,histoDirectPhotonSpectrum,1,0.5);
   if(grStatPoints1)grStatPoints1->SetName("graphErrorsConsistWithZero");//12
   TGraphAsymmErrors *grStatLimit1 = getPhotonPoints(histoDirectPhotonSpectrumStat,histoDirectPhotonSpectrum,2,0.5);
   if(grStatLimit1) grStatLimit1->SetName("graphUpperLimitsWithPoints");//12
   TGraphAsymmErrors *grStatArrow1 = getPhotonPoints(histoDirectPhotonSpectrumStat,histoDirectPhotonSpectrum,3,0.5);
   if(grStatArrow1) grStatArrow1->SetName("graphArrowForErrorConsitsWithZero");//12
   TGraphAsymmErrors *grStatLimit2 = getPhotonPoints(histoDirectPhotonSpectrumStat,histoDirectPhotonSpectrum,4,0.5);
   if(grStatLimit2) grStatLimit2->SetName("graphUpperLimitWithPointConsitsWithZero");//1
   TGraphAsymmErrors *grStatArrow2 = getPhotonPoints(histoDirectPhotonSpectrumStat,histoDirectPhotonSpectrum,5,0.5);
   if(grStatArrow2) grStatArrow2->SetName("graphArrowForPointsConsitsWithZero");//1
   TGraphAsymmErrors *grStatAll = getPhotonPoints(histoDirectPhotonSpectrumStat,histoDirectPhotonSpectrum,6,0.5);
   if(grStatAll) grStatAll->SetName("graphAll");


   TGraphAsymmErrors *grSystFull = getPhotonPoints(histoDirectPhotonSpectrumSyst,histoDirectPhotonSpectrum,0.5,0.5);
   if(grSystFull) grSystFull->SetName("graphDataPointsSystError");
   TGraphAsymmErrors *grSystFullA = getPhotonPoints(histoDirectPhotonSpectrumSystA,histoDirectPhotonSpectrum,0.5,0.5);
   if(grSystFullA) grSystFullA->SetName("graphDataPointsSystErrorA");
   TGraphAsymmErrors *grSystFullB = getPhotonPoints(histoDirectPhotonSpectrumSystB,histoDirectPhotonSpectrum,0.5,0.5);
   if(grSystFullB) grSystFullB->SetName("graphDataPointsSystErrorB");
   TGraphAsymmErrors *grSystFullC = getPhotonPoints(histoDirectPhotonSpectrumSystC,histoDirectPhotonSpectrum,0.5,0.5);
   if(grSystFullC) grSystFullC->SetName("graphDataPointsSystErrorC");
   TGraphAsymmErrors *grSystPoints1 = getPhotonPoints(histoDirectPhotonSpectrumSyst,histoDirectPhotonSpectrum,1,0.5);
   if(grSystPoints1)grSystPoints1->SetName("graphErrorsConsistWithZero");//12
   TGraphAsymmErrors *grSystLimit1 = getPhotonPoints(histoDirectPhotonSpectrumSyst,histoDirectPhotonSpectrum,2,0.5);
   if(grSystLimit1) grSystLimit1->SetName("graphUpperLimitsWithPoints");//12
   TGraphAsymmErrors *grSystArrow1 = getPhotonPoints(histoDirectPhotonSpectrumSyst,histoDirectPhotonSpectrum,3,0.5);
   if(grSystArrow1) grSystArrow1->SetName("graphArrowForErrorConsitsWithZero");//12
   TGraphAsymmErrors *grSystLimit2 = getPhotonPoints(histoDirectPhotonSpectrumSyst,histoDirectPhotonSpectrum,4,0.5);
   if(grSystLimit2) grSystLimit2->SetName("graphUpperLimitWithPointConsitsWithZero");//1
   TGraphAsymmErrors *grSystArrow2 = getPhotonPoints(histoDirectPhotonSpectrumSyst,histoDirectPhotonSpectrum,5,0.5);
   if(grSystArrow2) grSystArrow2->SetName("graphArrowForPointsConsitsWithZero");//1
   TGraphAsymmErrors *grSystAll = getPhotonPoints(histoDirectPhotonSpectrumSyst,histoDirectPhotonSpectrum,6,0.5);
   if(grSystAll) grSystAll->SetName("graphAll");


   TGraphAsymmErrors *grSystWithoutMatFull = getPhotonPoints(histoDirectPhotonSpectrumSystWithoutMat,histoDirectPhotonSpectrum,0.5,0.5);
   if(grSystWithoutMatFull) grSystWithoutMatFull->SetName("graphDataPointsSystWithoutMatError");
   TGraphAsymmErrors *grSystWithoutMatPoints1 = getPhotonPoints(histoDirectPhotonSpectrumSystWithoutMat,histoDirectPhotonSpectrum,1,0.5);
   if(grSystWithoutMatPoints1)grSystWithoutMatPoints1->SetName("graphErrorsConsistWithZero");//12
   TGraphAsymmErrors *grSystWithoutMatLimit1 = getPhotonPoints(histoDirectPhotonSpectrumSystWithoutMat,histoDirectPhotonSpectrum,2,0.5);
   if(grSystWithoutMatLimit1) grSystWithoutMatLimit1->SetName("graphUpperLimitsWithPoints");//12
   TGraphAsymmErrors *grSystWithoutMatArrow1 = getPhotonPoints(histoDirectPhotonSpectrumSystWithoutMat,histoDirectPhotonSpectrum,3,0.5);
   if(grSystWithoutMatArrow1) grSystWithoutMatArrow1->SetName("graphArrowForErrorConsitsWithZero");//12
   TGraphAsymmErrors *grSystWithoutMatLimit2 = getPhotonPoints(histoDirectPhotonSpectrumSystWithoutMat,histoDirectPhotonSpectrum,4,0.5);
   if(grSystWithoutMatLimit2) grSystWithoutMatLimit2->SetName("graphUpperLimitWithPointConsitsWithZero");//1
   TGraphAsymmErrors *grSystWithoutMatArrow2 = getPhotonPoints(histoDirectPhotonSpectrumSystWithoutMat,histoDirectPhotonSpectrum,5,0.5);
   if(grSystWithoutMatArrow2) grSystWithoutMatArrow2->SetName("graphArrowForPointsConsitsWithZero");//1
   TGraphAsymmErrors *grSystWithoutMatAll = getPhotonPoints(histoDirectPhotonSpectrumSystWithoutMat,histoDirectPhotonSpectrum,6,0.5);
   if(grSystWithoutMatAll) grSystWithoutMatAll->SetName("graphAll");


   TGraphAsymmErrors *graphPointsErrorsAbove0 = getPhotonPoints(histoDirectPhotonSpectrumSummed,histoDirectPhotonSpectrum,0,.5);
   if(graphPointsErrorsAbove0) graphPointsErrorsAbove0->SetName("graphPointsErrorsAbove0");
   TGraphAsymmErrors *graphPointsErrorsAboveButConfidence = getPhotonPoints(histoDirectPhotonSpectrumSummed,histoDirectPhotonSpectrum,1,.5);
   if(graphPointsErrorsAboveButConfidence)graphPointsErrorsAboveButConfidence->SetName("graphPointsErrorsAboveButConfidence");//12
   TGraphAsymmErrors *graphUpperLimitForPointAboveButError = getPhotonPoints(histoDirectPhotonSpectrumSummed,histoDirectPhotonSpectrum,2,.5);
   if(graphUpperLimitForPointAboveButError) graphUpperLimitForPointAboveButError->SetName("graphUpperLimitForPointAboveButError");//12
   TGraphAsymmErrors *graphArrowConfi = getPhotonPoints(histoDirectPhotonSpectrumSummed,histoDirectPhotonSpectrum,3,0.5);
   if(graphArrowConfi) graphArrowConfi->SetName("graphArrowConfi");//12
   TGraphAsymmErrors *graphUpperLimitConfi = getPhotonPoints(histoDirectPhotonSpectrumSummed,histoDirectPhotonSpectrum,10,0.5);
   if(graphUpperLimitConfi) graphUpperLimitConfi->SetName("graphUpperLimitConfi");//12
   TGraphAsymmErrors *graphUpperLimitPointsAndErrorBelow0 = getPhotonPoints(histoDirectPhotonSpectrumSummed,histoDirectPhotonSpectrum,4,.5);
   if(graphUpperLimitPointsAndErrorBelow0) graphUpperLimitPointsAndErrorBelow0->SetName("graphUpperLimitPointsAndErrorBelow0");//1
   TGraphAsymmErrors *graphArrorErrorBelow0 = getPhotonPoints(histoDirectPhotonSpectrumSummed,histoDirectPhotonSpectrum,5,.5);
   if(graphArrorErrorBelow0) graphArrorErrorBelow0->SetName("graphArrorErrorBelow0");//1
   TGraphAsymmErrors *graphArrowPointErrorBelow0 = getPhotonPoints(histoDirectPhotonSpectrumSummed,histoDirectPhotonSpectrum,6,0.5);
   if(graphArrowPointErrorBelow0) graphArrowPointErrorBelow0->SetName("graphArrowPointErrorBelow0");


   TGraphAsymmErrors *grStatFullSubNLO = getPhotonPoints(histoDirectPhotonSpectrumSubNLOStat,histoDirectPhotonSpectrumSubNLO,0.75,0.5);
   if(grStatFullSubNLO) grStatFullSubNLO->SetName("graphDataPointsStatError");
   TGraphAsymmErrors *grSystFullSubNLO = getPhotonPoints(histoDirectPhotonSpectrumSubNLOSyst,histoDirectPhotonSpectrumSubNLO,0.5,0.5);
   if(grSystFullSubNLO) grSystFullSubNLO->SetName("graphDataPointsSystError");
   TGraphAsymmErrors *graphUpperLimitForPointAboveButErrorSubNLO = getPhotonPoints(histoDirectPhotonSpectrumSubNLOSummend,histoDirectPhotonSpectrumSubNLO,2,.5);
   if(graphUpperLimitForPointAboveButErrorSubNLO) graphUpperLimitForPointAboveButErrorSubNLO->SetName("graphUpperLimitForPointAboveButErrorSubNLO");//12
   TGraphAsymmErrors *graphUpperLimitPointsAndErrorBelow0SubNLO = getPhotonPoints(histoDirectPhotonSpectrumSubNLOSummend,histoDirectPhotonSpectrumSubNLO,4,.5);
   if(graphUpperLimitPointsAndErrorBelow0SubNLO) graphUpperLimitPointsAndErrorBelow0SubNLO->SetName("graphUpperLimitPointsAndErrorBelow0SubNLO");//1
   TGraphAsymmErrors *graphArrorErrorBelow0SubNLO = getPhotonPoints(histoDirectPhotonSpectrumSubNLOSummend,histoDirectPhotonSpectrumSubNLO,5,.5);
   if(graphArrorErrorBelow0SubNLO) graphArrorErrorBelow0SubNLO->SetName("graphArrorErrorBelow0SubNLO");//1
   TGraphAsymmErrors *graphArrowPointErrorBelow0SubNLO = getPhotonPoints(histoDirectPhotonSpectrumSubNLOSummend,histoDirectPhotonSpectrumSubNLO,6,0.5);
   if(graphArrowPointErrorBelow0SubNLO) graphArrowPointErrorBelow0SubNLO->SetName("graphArrowPointErrorBelow0SubNLO");

   TGraphAsymmErrors *grFullStatErrorAllPoints = getPhotonPoints(histoDirectPhotonSpectrumStat,histoDirectPhotonSpectrum,0.5,0.5);
   grFullStatErrorAllPoints->SetName("grFullStatErrorAllPoints");
   TGraphAsymmErrors *grFullSystErrorAllPoints = getPhotonPoints(histoDirectPhotonSpectrumSyst,histoDirectPhotonSpectrum,0.5,0.5);
   grFullSystErrorAllPoints->SetName("grFullSystErrorAllPoints");
   TGraphAsymmErrors *grFullSystWithoutMatErrorAllPoints = getPhotonPoints(histoDirectPhotonSpectrumSystWithoutMat,histoDirectPhotonSpectrum,0.5,0.5);
   grFullSystWithoutMatErrorAllPoints->SetName("grFullSystWithoutMatErrorAllPoints");
   TGraphAsymmErrors *grFullSummedErrorAllPoints = getPhotonPoints(histoDirectPhotonSpectrumSummed,histoDirectPhotonSpectrum,0.5,0.5);
   if(grFullSummedErrorAllPoints)grFullSummedErrorAllPoints->SetName("grFullSummedErrorAllPointsNotShifted");

   cout<<"------------------------------------"<<endl;

   for(Int_t i = 0;i<grSystFull->GetN();i++){
      for(Int_t j = 0;j<meanErrorsCorrSummedIncMatGamma->GetN();j++){
         if(grSystFull->GetX()[i]==meanErrorsCorrSummedIncMatGamma->GetX()[j]){
            grSystFull->SetPointError(i,grSystFull->GetErrorXlow(i),grSystFull->GetErrorXhigh(i),
                                      grSystFull->GetY()[i]*sqrt( pow(grSystFull->GetErrorYlow(i)/grSystFull->GetY()[i],2) + pow(meanErrorsCorrSummedIncMatGamma->GetY()[j]/100,2) ),
                                      grSystFull->GetY()[i]*sqrt( pow(grSystFull->GetErrorYhigh(i)/grSystFull->GetY()[i],2) + pow(meanErrorsCorrSummedIncMatGamma->GetY()[j]/100,2)));
         }
      }
   }
   for(Int_t i = 0;i<grSystFullA->GetN();i++){
      for(Int_t j = 0;j<meanErrorsCorrSummedIncMatGammaA->GetN();j++){
         if(grSystFullA->GetX()[i]==meanErrorsCorrSummedIncMatGammaA->GetX()[j]){
            grSystFullA->SetPointError(i,grSystFullA->GetErrorXlow(i),grSystFullA->GetErrorXhigh(i),
                                      grSystFullA->GetY()[i]*sqrt( pow(grSystFullA->GetErrorYlow(i)/grSystFullA->GetY()[i],2) + pow(meanErrorsCorrSummedIncMatGammaA->GetY()[j]/100,2) ),
                                      grSystFullA->GetY()[i]*sqrt( pow(grSystFullA->GetErrorYhigh(i)/grSystFullA->GetY()[i],2) + pow(meanErrorsCorrSummedIncMatGammaA->GetY()[j]/100,2)));
         }
      }
   }
   for(Int_t i = 0;i<grSystFullB->GetN();i++){
      for(Int_t j = 0;j<meanErrorsCorrSummedIncMatGammaB->GetN();j++){
         if(grSystFullB->GetX()[i]==meanErrorsCorrSummedIncMatGammaB->GetX()[j]){
            grSystFullB->SetPointError(i,grSystFullB->GetErrorXlow(i),grSystFullB->GetErrorXhigh(i),
                                      grSystFullB->GetY()[i]*sqrt( pow(grSystFullB->GetErrorYlow(i)/grSystFullB->GetY()[i],2) + pow(meanErrorsCorrSummedIncMatGammaB->GetY()[j]/100,2) ),
                                      grSystFullB->GetY()[i]*sqrt( pow(grSystFullB->GetErrorYhigh(i)/grSystFullB->GetY()[i],2) + pow(meanErrorsCorrSummedIncMatGammaB->GetY()[j]/100,2)));
         }
      }
   }
   for(Int_t i = 0;i<grSystFullC->GetN();i++){
      for(Int_t j = 0;j<meanErrorsCorrSummedIncMatGammaC->GetN();j++){
         if(grSystFullC->GetX()[i]==meanErrorsCorrSummedIncMatGammaC->GetX()[j]){
            grSystFullC->SetPointError(i,grSystFullC->GetErrorXlow(i),grSystFullC->GetErrorXhigh(i),
                                      grSystFullC->GetY()[i]*sqrt( pow(grSystFullC->GetErrorYlow(i)/grSystFullC->GetY()[i],2) + pow(meanErrorsCorrSummedIncMatGammaC->GetY()[j]/100,2) ),
                                      grSystFullC->GetY()[i]*sqrt( pow(grSystFullC->GetErrorYhigh(i)/grSystFullC->GetY()[i],2) + pow(meanErrorsCorrSummedIncMatGammaC->GetY()[j]/100,2)));
         }
      }
   }


   TGraphAsymmErrors *grStatFullCombined = getPhotonPoints(histoDirectPhotonSpectrumCombinedStat,histoDirectPhotonSpectrumCombined,0.75,0.5);
   if(grStatFullCombined) grStatFullCombined->SetName("graphDataPointsStatError");
   TGraphAsymmErrors *grStatPoints1Combined = getPhotonPoints(histoDirectPhotonSpectrumCombinedStat,histoDirectPhotonSpectrumCombined,1,0.5);
   if(grStatPoints1Combined)grStatPoints1Combined->SetName("graphErrorsConsistWithZero");//12
   TGraphAsymmErrors *grStatLimit1Combined = getPhotonPoints(histoDirectPhotonSpectrumCombinedStat,histoDirectPhotonSpectrumCombined,2,0.5);
   if(grStatLimit1Combined) grStatLimit1Combined->SetName("graphUpperLimitsWithPoints");//12
   TGraphAsymmErrors *grStatArrow1Combined = getPhotonPoints(histoDirectPhotonSpectrumCombinedStat,histoDirectPhotonSpectrumCombined,3,0.5);
   if(grStatArrow1Combined) grStatArrow1Combined->SetName("graphArrowForErrorConsitsWithZero");//12
   TGraphAsymmErrors *grStatLimit2Combined = getPhotonPoints(histoDirectPhotonSpectrumCombinedStat,histoDirectPhotonSpectrumCombined,4,0.5);
   if(grStatLimit2Combined) grStatLimit2Combined->SetName("graphUpperLimitWithPointConsitsWithZero");//1
   TGraphAsymmErrors *grStatArrow2Combined = getPhotonPoints(histoDirectPhotonSpectrumCombinedStat,histoDirectPhotonSpectrumCombined,5,0.5);
   if(grStatArrow2Combined) grStatArrow2Combined->SetName("graphArrowForPointsConsitsWithZero");//1
   TGraphAsymmErrors *grStatAllCombined = getPhotonPoints(histoDirectPhotonSpectrumCombinedStat,histoDirectPhotonSpectrumCombined,6,0.5);
   if(grStatAllCombined) grStatAllCombined->SetName("graphAll");

   TGraphAsymmErrors *grSystFullCombined = getPhotonPoints(histoDirectPhotonSpectrumCombinedSyst,histoDirectPhotonSpectrumCombined,0.5,0.5);
   if(grSystFullCombined) grSystFullCombined->SetName("graphDataPointsSystError");
   TGraphAsymmErrors *grSystFullCombinedA = getPhotonPoints(histoDirectPhotonSpectrumCombinedSystA,histoDirectPhotonSpectrumCombined,0.5,0.5);
   if(grSystFullCombinedA) grSystFullCombinedA->SetName("graphDataPointsSystErrorA");
   TGraphAsymmErrors *grSystPoints1Combined = getPhotonPoints(histoDirectPhotonSpectrumCombinedSyst,histoDirectPhotonSpectrumCombined,1,0.5);
   if(grSystPoints1Combined)grSystPoints1Combined->SetName("graphErrorsConsistWithZero");//12
   TGraphAsymmErrors *grSystLimit1Combined = getPhotonPoints(histoDirectPhotonSpectrumCombinedSyst,histoDirectPhotonSpectrumCombined,2,0.5);
   if(grSystLimit1Combined) grSystLimit1Combined->SetName("graphUpperLimitsWithPoints");//12
   TGraphAsymmErrors *grSystArrow1Combined = getPhotonPoints(histoDirectPhotonSpectrumCombinedSyst,histoDirectPhotonSpectrumCombined,3,0.5);
   if(grSystArrow1Combined) grSystArrow1Combined->SetName("graphArrowForErrorConsitsWithZero");//12
   TGraphAsymmErrors *grSystLimit2Combined = getPhotonPoints(histoDirectPhotonSpectrumCombinedSyst,histoDirectPhotonSpectrumCombined,4,0.5);
   if(grSystLimit2Combined) grSystLimit2Combined->SetName("graphUpperLimitWithPointConsitsWithZero");//1
   TGraphAsymmErrors *grSystArrow2Combined = getPhotonPoints(histoDirectPhotonSpectrumCombinedSyst,histoDirectPhotonSpectrumCombined,5,0.5);
   if(grSystArrow2Combined) grSystArrow2Combined->SetName("graphArrowForPointsConsitsWithZero");//1
   TGraphAsymmErrors *grSystAllCombined = getPhotonPoints(histoDirectPhotonSpectrumCombinedSyst,histoDirectPhotonSpectrumCombined,6,0.5);
   if(grSystAllCombined) grSystAllCombined->SetName("graphAll");


   TGraphAsymmErrors *graphPointsErrorsAbove0Combined = getPhotonPoints(histoDirectPhotonSpectrumCombinedSum,histoDirectPhotonSpectrumCombined,0,.5);
   if(graphPointsErrorsAbove0Combined) graphPointsErrorsAbove0Combined->SetName("graphPointsErrorsAbove0");
   TGraphAsymmErrors *graphPointsErrorsAboveButConfidenceCombined = getPhotonPoints(histoDirectPhotonSpectrumCombinedSum,histoDirectPhotonSpectrumCombined,1,.5);
   if(graphPointsErrorsAboveButConfidenceCombined)graphPointsErrorsAboveButConfidenceCombined->SetName("graphPointsErrorsAboveButConfidence");//12
   TGraphAsymmErrors *graphUpperLimitForPointAboveButErrorCombined = getPhotonPoints(histoDirectPhotonSpectrumCombinedSum,histoDirectPhotonSpectrumCombined,2,.5);
   if(graphUpperLimitForPointAboveButErrorCombined) graphUpperLimitForPointAboveButErrorCombined->SetName("graphUpperLimitForPointAboveButError");//12
   TGraphAsymmErrors *graphArrowConfiCombined = getPhotonPoints(histoDirectPhotonSpectrumCombinedSum,histoDirectPhotonSpectrumCombined,3,0.5);
   if(graphArrowConfiCombined) graphArrowConfiCombined->SetName("graphArrowConfi");//12
   TGraphAsymmErrors *graphUpperLimitConfiCombined = getPhotonPoints(histoDirectPhotonSpectrumCombinedSum,histoDirectPhotonSpectrumCombined,10,0.5);
   if(graphUpperLimitConfiCombined) graphUpperLimitConfiCombined->SetName("graphUpperLimitConfi");//12
   TGraphAsymmErrors *graphUpperLimitPointsAndErrorBelow0Combined = getPhotonPoints(histoDirectPhotonSpectrumCombinedSum,histoDirectPhotonSpectrumCombined,4,.5);
   if(graphUpperLimitPointsAndErrorBelow0Combined) graphUpperLimitPointsAndErrorBelow0Combined->SetName("graphUpperLimitPointsAndErrorBelow0");//1
   TGraphAsymmErrors *graphArrorErrorBelow0Combined = getPhotonPoints(histoDirectPhotonSpectrumCombinedSum,histoDirectPhotonSpectrumCombined,5,.5);
   if(graphArrorErrorBelow0Combined) graphArrorErrorBelow0Combined->SetName("graphArrorErrorBelow0");//1
   TGraphAsymmErrors *graphArrowPointErrorBelow0Combined = getPhotonPoints(histoDirectPhotonSpectrumCombinedSum,histoDirectPhotonSpectrumCombined,6,0.5);
   if(graphArrowPointErrorBelow0Combined) graphArrowPointErrorBelow0Combined->SetName("graphArrowPointErrorBelow0");


   TGraphAsymmErrors *grFullStatErrorAllPointsCombined = getPhotonPoints(histoDirectPhotonSpectrumCombinedStat,histoDirectPhotonSpectrumCombined,0.5,0.5);
   grFullStatErrorAllPointsCombined->SetName("grFullStatErrorAllPoints");
   TGraphAsymmErrors *grFullSystErrorAllPointsCombined = getPhotonPoints(histoDirectPhotonSpectrumCombinedSyst,histoDirectPhotonSpectrumCombined,0.5,0.5);
   grFullSystErrorAllPointsCombined->SetName("grFullSystErrorAllPoints");
   TGraphAsymmErrors *grFullSummedErrorAllPointsCombined = getPhotonPoints(histoDirectPhotonSpectrumCombinedSum,histoDirectPhotonSpectrumCombined,0.5,0.5);
   grFullSummedErrorAllPointsCombined->SetName("grFullSummedErrorAllPointsNotShifted");



   for(Int_t i = 0;i<grSystFullCombined->GetN();i++){
      for(Int_t j = 0;j<GammaCombinedSyst->GetN();j++){
         if(grSystFullCombined->GetX()[i]==GammaCombinedSyst->GetX()[j]){
            grSystFullCombined->SetPointError(i,grSystFullCombined->GetErrorXlow(i),grSystFullCombined->GetErrorXhigh(i),
                                              grSystFullCombined->GetY()[i]*sqrt( pow(grSystFullCombined->GetErrorYlow(i)/grSystFullCombined->GetY()[i],2) + pow(GammaCombinedSyst->GetErrorYlow(j)/GammaCombinedSyst->GetY()[j],2)),
                                              grSystFullCombined->GetY()[i]*sqrt( pow(grSystFullCombined->GetErrorYhigh(i)/grSystFullCombined->GetY()[i],2) + pow(GammaCombinedSyst->GetErrorYhigh(j)/GammaCombinedSyst->GetY()[j],2)));
         }
      }
   }
   for(Int_t i = 0;i<grSystFullCombinedA->GetN();i++){
      for(Int_t j = 0;j<GammaCombinedSystAB->GetN();j++){
         if(grSystFullCombinedA->GetX()[i]==GammaCombinedSystAB->GetX()[j]){
            grSystFullCombinedA->SetPointError(i,grSystFullCombinedA->GetErrorXlow(i),grSystFullCombinedA->GetErrorXhigh(i),
                                              grSystFullCombinedA->GetY()[i]*sqrt( pow(grSystFullCombinedA->GetErrorYlow(i)/grSystFullCombinedA->GetY()[i],2) + pow(GammaCombinedSystAB->GetErrorYlow(j)/GammaCombinedSystAB->GetY()[j],2)),
                                              grSystFullCombinedA->GetY()[i]*sqrt( pow(grSystFullCombinedA->GetErrorYhigh(i)/grSystFullCombinedA->GetY()[i],2) + pow(GammaCombinedSystAB->GetErrorYhigh(j)/GammaCombinedSystAB->GetY()[j],2)));

         }
      }
   }

   // histoDirectPhotonSpectrumCombined->Draw("e1");
   // grFullSystErrorAllPointsCombined->Draw("e1p");
   // grFullStatErrorAllPointsCombined->Draw("e1p");
   // return;
   cout<<"------------------------------------"<<endl;



   TGraphAsymmErrors *grFullStatErrorAllPointsSubNLO = getPhotonPoints(histoDirectPhotonSpectrumSubNLOStat,histoDirectPhotonSpectrumSubNLO,0.5,0.5);
   if(grFullStatErrorAllPointsSubNLO)grFullStatErrorAllPointsSubNLO->SetName("grFullStatErrorAllPointsSubNLO");
   TGraphAsymmErrors *grFullSystErrorAllPointsSubNLO = getPhotonPoints(histoDirectPhotonSpectrumSubNLOSyst,histoDirectPhotonSpectrumSubNLO,0.5,0.5);
   if(grFullSystErrorAllPointsSubNLO)grFullSystErrorAllPointsSubNLO->SetName("grFullSystErrorAllPointsSubNLO");
   TGraphAsymmErrors *grFullSummedErrorAllPointsSubNLO = getPhotonPoints(histoDirectPhotonSpectrumSubNLOSummend,histoDirectPhotonSpectrumSubNLO,0.5,0.5);
   if(grFullSummedErrorAllPointsSubNLO)grFullSummedErrorAllPointsSubNLO->SetName("grFullSummedErrorAllPointsNotShiftedSubNLO");

   TGraphAsymmErrors *graAllShifted;
   if(grFullSummedErrorAllPoints) graAllShifted = (TGraphAsymmErrors*)grFullSummedErrorAllPoints->Clone("graAllShifted");
   TGraphAsymmErrors *graAllShiftedStat = (TGraphAsymmErrors*)grFullStatErrorAllPoints->Clone("graAllShiftedStat");
   TGraphAsymmErrors *graAllShiftedSyst = (TGraphAsymmErrors*)grFullSystErrorAllPoints->Clone("graAllShiftedSyst");
   TGraphAsymmErrors *graAllShiftedSystWithoutMat = (TGraphAsymmErrors*)grFullSystWithoutMatErrorAllPoints->Clone("graAllShiftedSystWithoutMat");
   TGraphAsymmErrors *graAllShiftedSubtractedNLO;
   if(grFullSummedErrorAllPointsSubNLO) graAllShiftedSubtractedNLO = (TGraphAsymmErrors*)grFullSummedErrorAllPointsSubNLO->Clone("graAllShiftedSystSubtractedNLO");
   TGraphAsymmErrors *graAllShiftedSystSubtractedNLO;
   if(grFullSystErrorAllPointsSubNLO)graAllShiftedSystSubtractedNLO = (TGraphAsymmErrors*)grFullSystErrorAllPointsSubNLO->Clone("graAllShiftedSystSubtractedNLO");
   TGraphAsymmErrors *graAllShiftedStatSubtractedNLO;
   if(grFullStatErrorAllPointsSubNLO)graAllShiftedStatSubtractedNLO = (TGraphAsymmErrors*)grFullStatErrorAllPointsSubNLO->Clone("graAllShiftedStatSubtractedNLO");
   TGraphAsymmErrors *graAllShiftedCombined = (TGraphAsymmErrors*)grFullSummedErrorAllPointsCombined->Clone("graAllShifted");
   TGraphAsymmErrors *graAllShiftedStatCombined = (TGraphAsymmErrors*)grFullStatErrorAllPointsCombined->Clone("graAllShiftedStatCombined");
   TGraphAsymmErrors *graAllShiftedSystCombined = (TGraphAsymmErrors*)grFullSystErrorAllPointsCombined->Clone("graAllShiftedSystCombined");

   Double_t Params[5];
   GetFitParameter("l",GetCentralityStringC(cutSel),Params);

   TF1 *qcdFit2 = FitObject("tmpt","QCDbinShiftFitForDirectGamma2","Pi0");
   //TF1 *qcdFit2 = FitObject("h","fitNLOMuOne","Pi0",graAllShifted,.9,12.,Params);//,paramGraphNLOMuOne);

   qcdFit2->SetRange(0,25);
   if(graAllShifted)graAllShifted = ApplyXshift(graAllShifted,qcdFit2);
   if(graAllShifted)graAllShifted->SetName("graphAllXshifted");
   qcdFit2 = FitObject("tmpt","QCDbinShiftFitForDirectGamma2","Pi0");
   if(graAllShiftedStat)graAllShiftedStat = ApplyXshift(graAllShiftedStat,qcdFit2);
   if(graAllShiftedStat)graAllShiftedStat->SetName("graphAllXshiftedStat");
   qcdFit2 = FitObject("tmpt","QCDbinShiftFitForDirectGamma2","Pi0");
   if(graAllShiftedSyst)graAllShiftedSyst = ApplyXshift(graAllShiftedSyst,qcdFit2);
   if(graAllShiftedSyst)graAllShiftedSyst->SetName("graphAllXshiftedSyst");
   qcdFit2 = FitObject("tmpt","QCDbinShiftFitForDirectGamma2","Pi0");
   if(graAllShiftedSystWithoutMat)graAllShiftedSystWithoutMat = ApplyXshift(graAllShiftedSystWithoutMat,qcdFit2);
   if(graAllShiftedSystWithoutMat)graAllShiftedSystWithoutMat->SetName("graphAllXshiftedSystWithoutMat");
   if(graAllShiftedCombined)graAllShiftedCombined = ApplyXshift(graAllShiftedCombined,qcdFit2);
   if(graAllShiftedCombined)graAllShiftedCombined->SetName("graphAllXshiftedCombined");
   if(graAllShiftedStatCombined)graAllShiftedStatCombined = ApplyXshift(graAllShiftedStatCombined,qcdFit2);
   if(graAllShiftedStatCombined)graAllShiftedStatCombined->SetName("graphAllXshiftedStatCombined");
   qcdFit2 = FitObject("tmpt","QCDbinShiftFitForDirectGamma2","Pi0");
   if(graAllShiftedSystCombined)graAllShiftedSystCombined = ApplyXshift(graAllShiftedSystCombined,qcdFit2);
   if(graAllShiftedSystCombined)graAllShiftedSystCombined->SetName("graphAllXshiftedSystCombined");


   if(!(!centralityCutNumber.CompareTo("524") || !centralityCutNumber.CompareTo("548"))){
      Double_t x,y;
      graAllShiftedSubtractedNLO->GetPoint(graAllShiftedSubtractedNLO->GetN()-1,x,y);
      if(x>4)graAllShiftedSubtractedNLO->RemovePoint(graAllShiftedSubtractedNLO->GetN()-1);
      graAllShiftedSubtractedNLO->GetPoint(graAllShiftedSubtractedNLO->GetN()-1,x,y);
      if(x>4)graAllShiftedSubtractedNLO->RemovePoint(graAllShiftedSubtractedNLO->GetN()-1);
      graAllShiftedSubtractedNLO->GetPoint(graAllShiftedSubtractedNLO->GetN()-1,x,y);
      if(x>4)graAllShiftedSubtractedNLO->RemovePoint(graAllShiftedSubtractedNLO->GetN()-1);
      graAllShiftedSubtractedNLO->GetPoint(graAllShiftedSubtractedNLO->GetN()-1,x,y);
      if(x>4)graAllShiftedSubtractedNLO->RemovePoint(graAllShiftedSubtractedNLO->GetN()-1);

      qcdFit2 = FitObject("tmpt","QCDbinShiftFitForDirectGamma2","Pi0");
      graAllShiftedSubtractedNLO = ApplyXshift(graAllShiftedSubtractedNLO,qcdFit2);
      graAllShiftedSubtractedNLO->SetName("graphAllXshiftedSubtractedNLO");

      graAllShiftedStatSubtractedNLO->GetPoint(graAllShiftedStatSubtractedNLO->GetN()-1,x,y);
      if(x>4)graAllShiftedStatSubtractedNLO->RemovePoint(graAllShiftedStatSubtractedNLO->GetN()-1);
      graAllShiftedStatSubtractedNLO->GetPoint(graAllShiftedStatSubtractedNLO->GetN()-1,x,y);
      if(x>4)graAllShiftedStatSubtractedNLO->RemovePoint(graAllShiftedStatSubtractedNLO->GetN()-1);
      graAllShiftedStatSubtractedNLO->GetPoint(graAllShiftedStatSubtractedNLO->GetN()-1,x,y);
      if(x>4)graAllShiftedStatSubtractedNLO->RemovePoint(graAllShiftedStatSubtractedNLO->GetN()-1);
      graAllShiftedStatSubtractedNLO->GetPoint(graAllShiftedStatSubtractedNLO->GetN()-1,x,y);
      if(x>4)graAllShiftedStatSubtractedNLO->RemovePoint(graAllShiftedStatSubtractedNLO->GetN()-1);
      qcdFit2 = FitObject("tmpt","QCDbinShiftFitForDirectGamma2","Pi0");
      graAllShiftedStatSubtractedNLO = ApplyXshift(graAllShiftedStatSubtractedNLO,qcdFit2);
      graAllShiftedStatSubtractedNLO->SetName("graphAllXshiftedStatSubtractedNLO");

      graAllShiftedSystSubtractedNLO->GetPoint(graAllShiftedSystSubtractedNLO->GetN()-1,x,y);
      if(x>4)graAllShiftedSystSubtractedNLO->RemovePoint(graAllShiftedSystSubtractedNLO->GetN()-1);
      graAllShiftedSystSubtractedNLO->GetPoint(graAllShiftedSystSubtractedNLO->GetN()-1,x,y);
      if(x>4)graAllShiftedSystSubtractedNLO->RemovePoint(graAllShiftedSystSubtractedNLO->GetN()-1);
      graAllShiftedSystSubtractedNLO->GetPoint(graAllShiftedSystSubtractedNLO->GetN()-1,x,y);
      if(x>4)graAllShiftedSystSubtractedNLO->RemovePoint(graAllShiftedSystSubtractedNLO->GetN()-1);
      graAllShiftedSystSubtractedNLO->GetPoint(graAllShiftedSystSubtractedNLO->GetN()-1,x,y);
      if(x>4)graAllShiftedSystSubtractedNLO->RemovePoint(graAllShiftedSystSubtractedNLO->GetN()-1);

      qcdFit2 = FitObject("tmpt","QCDbinShiftFitForDirectGamma2","Pi0");
      graAllShiftedSystSubtractedNLO = ApplyXshift(graAllShiftedSystSubtractedNLO,qcdFit2);
      graAllShiftedSystSubtractedNLO->SetName("graphAllXshiftedSystSubtractedNLO");
   }



   Double_t *xValusShifted;
   Double_t *xErrowLow;
   Double_t *xErrowHigh;
   Double_t *yValusShifted;
   if(graAllShifted){
      xValusShifted = new Double_t[graAllShifted->GetN()];
      xErrowLow = new Double_t[graAllShifted->GetN()];
      xErrowHigh = new Double_t[graAllShifted->GetN()];
      yValusShifted = new Double_t[graAllShifted->GetN()];


      xValusShifted = graAllShifted->GetX();
      yValusShifted = graAllShifted->GetY();
      xErrowLow = graAllShifted->GetEXlow();
      xErrowHigh = graAllShifted->GetEXhigh();

      for(Int_t i = 0;i<graAllShifted->GetN();i++){
         Double_t lowEdge = xValusShifted[i]-xErrowLow[i];
         Double_t upEdge = xValusShifted[i]+xErrowHigh[i];
         Double_t *xOld = grSystFull->GetX();
         cout<<lowEdge<<"  "<<upEdge<<" "<<xValusShifted[i]<<endl;
         for(Int_t j = 0;j<grSystFull->GetN();j++){
            if(xOld[j]>lowEdge && xOld[j]<upEdge){
               if(grSystFull)grSystFull->SetPoint(j,xValusShifted[i],yValusShifted[i]);
               if(grSystFull)grSystFull->SetPointEXhigh(j,xErrowHigh[i]);
               if(grSystFull)grSystFull->SetPointEXlow(j,xErrowLow[i]);
            }
         }
         xOld = grSystFullA->GetX();
         for(Int_t j = 0;j<grSystFullA->GetN();j++){
            if(xOld[j]>lowEdge && xOld[j]<upEdge){
               if(grSystFullA)grSystFullA->SetPoint(j,xValusShifted[i],yValusShifted[i]);
               if(grSystFullA)grSystFullA->SetPointEXhigh(j,xErrowHigh[i]);
               if(grSystFullA)grSystFullA->SetPointEXlow(j,xErrowLow[i]);
            }
         }
         xOld = grSystFullB->GetX();
         for(Int_t j = 0;j<grSystFullB->GetN();j++){
            if(xOld[j]>lowEdge && xOld[j]<upEdge){
               if(grSystFullB)grSystFullB->SetPoint(j,xValusShifted[i],yValusShifted[i]);
               if(grSystFullB)grSystFullB->SetPointEXhigh(j,xErrowHigh[i]);
               if(grSystFullB)grSystFullB->SetPointEXlow(j,xErrowLow[i]);
            }
         }
         xOld = grSystFullC->GetX();
         for(Int_t j = 0;j<grSystFullC->GetN();j++){
            if(xOld[j]>lowEdge && xOld[j]<upEdge){
               if(grSystFullC)grSystFullC->SetPoint(j,xValusShifted[i],yValusShifted[i]);
               if(grSystFullC)grSystFullC->SetPointEXhigh(j,xErrowHigh[i]);
               if(grSystFullC)grSystFullC->SetPointEXlow(j,xErrowLow[i]);
            }
         }
      }

      for(Int_t i = 0;i<graAllShifted->GetN();i++){
         Double_t lowEdge = xValusShifted[i]-xErrowLow[i];
         Double_t upEdge = xValusShifted[i]+xErrowHigh[i];
         Double_t *xOld = grStatFull->GetX();
         for(Int_t j = 0;j<grStatFull->GetN();j++){
            if(xOld[j]>lowEdge && xOld[j]<upEdge){
               grStatFull->SetPoint(j,xValusShifted[i],yValusShifted[i]);
            }
         }
      }
   }

   Double_t *xValusShiftedCombined = new Double_t[graAllShiftedCombined->GetN()];
   Double_t *xErrowLowCombined = new Double_t[graAllShiftedCombined->GetN()];
   Double_t *xErrowHighCombined = new Double_t[graAllShiftedCombined->GetN()];
   Double_t *yValusShiftedCombined = new Double_t[graAllShiftedCombined->GetN()];
   xValusShiftedCombined = graAllShiftedCombined->GetX();
   yValusShiftedCombined = graAllShiftedCombined->GetY();
   xErrowLowCombined = graAllShiftedCombined->GetEXlow();
   xErrowHighCombined = graAllShiftedCombined->GetEXhigh();

   for(Int_t i = 0;i<graAllShiftedCombined->GetN();i++){
      Double_t lowEdgeCombined = xValusShiftedCombined[i]-xErrowLowCombined[i];
      Double_t upEdgeCombined = xValusShiftedCombined[i]+xErrowHighCombined[i];
      Double_t *xOldCombined = grSystFullCombined->GetX();
      cout<<lowEdgeCombined<<"  "<<upEdgeCombined<<" "<<xValusShiftedCombined[i]<<endl;
      for(Int_t j = 0;j<grSystFullCombined->GetN();j++){
         if(xOldCombined[j]>lowEdgeCombined && xOldCombined[j]<upEdgeCombined){
            if(grSystFullCombined)grSystFullCombined->SetPoint(j,xValusShiftedCombined[i],yValusShiftedCombined[i]);
            if(grSystFullCombined)grSystFullCombined->SetPointEXhigh(j,xErrowHighCombined[i]);
            if(grSystFullCombined)grSystFullCombined->SetPointEXlow(j,xErrowLowCombined[i]);
         }
      }
      Double_t *xOldCombinedA = grSystFullCombinedA->GetX();
      for(Int_t j = 0;j<grSystFullCombinedA->GetN();j++){
         if(xOldCombinedA[j]>lowEdgeCombined && xOldCombinedA[j]<upEdgeCombined){
            if(grSystFullCombinedA)grSystFullCombinedA->SetPoint(j,xValusShiftedCombined[i],yValusShiftedCombined[i]);
            if(grSystFullCombinedA)grSystFullCombinedA->SetPointEXhigh(j,xErrowHighCombined[i]);
            if(grSystFullCombinedA)grSystFullCombinedA->SetPointEXlow(j,xErrowLowCombined[i]);
         }
      }
   }

   for(Int_t i = 0;i<graAllShiftedCombined->GetN();i++){
      Double_t lowEdgeCombined = xValusShiftedCombined[i]-xErrowLowCombined[i];
      Double_t upEdgeCombined = xValusShiftedCombined[i]+xErrowHighCombined[i];
      Double_t *xOldCombined = grStatFullCombined->GetX();
      for(Int_t j = 0;j<grStatFullCombined->GetN();j++){
         if(xOldCombined[j]>lowEdgeCombined && xOldCombined[j]<upEdgeCombined){
            grStatFullCombined->SetPoint(j,xValusShiftedCombined[i],yValusShiftedCombined[i]);
         }
      }
   }


   if(!(!centralityCutNumber.CompareTo("524") || !centralityCutNumber.CompareTo("548"))){
      Double_t *xValusShiftedSubNLO = new Double_t[graAllShiftedSubtractedNLO->GetN()];
      Double_t *xErrowLowSubNLO = new Double_t[graAllShiftedSubtractedNLO->GetN()];
      Double_t *xErrowHighSubNLO = new Double_t[graAllShiftedSubtractedNLO->GetN()];
      Double_t *yValusShiftedSubNLO = new Double_t[graAllShiftedSubtractedNLO->GetN()];

      xValusShiftedSubNLO = graAllShiftedStatSubtractedNLO->GetX();
      yValusShiftedSubNLO = graAllShiftedStatSubtractedNLO->GetY();
      xErrowLowSubNLO = graAllShiftedStatSubtractedNLO->GetEXlow();
      xErrowHighSubNLO = graAllShiftedStatSubtractedNLO->GetEXhigh();

      for(Int_t i = 0;i<graAllShiftedSubtractedNLO->GetN();i++){
         Double_t lowEdge = xValusShiftedSubNLO[i]-xErrowLowSubNLO[i];
         Double_t upEdge = xValusShiftedSubNLO[i]+xErrowHighSubNLO[i];
         Double_t *xOld = grSystFullSubNLO->GetX();
         cout<<lowEdge<<"  "<<upEdge<<" "<<xValusShiftedSubNLO[i]<<endl;
         for(Int_t j = 0;j<grSystFullSubNLO->GetN();j++){
            if(xOld[j]>lowEdge && xOld[j]<upEdge){
               if(grSystFullSubNLO)grSystFullSubNLO->SetPoint(j,xValusShiftedSubNLO[i],yValusShiftedSubNLO[i]);
               if(grSystFullSubNLO)grSystFullSubNLO->SetPointEXhigh(j,xErrowHighSubNLO[i]);
               if(grSystFullSubNLO)grSystFullSubNLO->SetPointEXlow(j,xErrowLowSubNLO[i]);
            }
         }
      }

      for(Int_t i = 0;i<graAllShiftedSubtractedNLO->GetN();i++){
         Double_t lowEdge = xValusShiftedSubNLO[i]-xErrowLowSubNLO[i];
         Double_t upEdge = xValusShiftedSubNLO[i]+xErrowHighSubNLO[i];
         Double_t *xOld = grStatFullSubNLO->GetX();
         for(Int_t j = 0;j<grStatFullSubNLO->GetN();j++){
            if(xOld[j]>lowEdge && xOld[j]<upEdge){
               if(grStatFullSubNLO)grStatFullSubNLO->SetPoint(j,xValusShiftedSubNLO[i],yValusShiftedSubNLO[i]);
            }
         }
      }
   }

   TCanvas *canvasDirectPhotonWithShift = GetAndSetCanvas("canvasDirectPhotonWithShift",0.15,0.1, 1100 ,1250);
   dummy->GetYaxis()->SetTitle("#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c}^{2})");
   dummy->GetXaxis()->SetRangeUser(0.50,11.4);
   dummy->GetYaxis()->SetTitleOffset(1.5);
   dummy->GetXaxis()->SetLabelOffset(0);
   dummy->DrawCopy();
   canvasDirectPhotonWithShift->SetLogy();
   dummy->GetYaxis()->SetRangeUser(1e-7,4e1);
   if(!centralityCutNumber.CompareTo("548")){
      dummy->GetYaxis()->SetRangeUser(2e-10,2e1);
   }
   dummy->DrawCopy();



   pp276CT10BFG2_sum_pdferr_InvYield->SetFillColor(kBlue-8);
   pp276CT10BFG2_sum_pdferr_InvYield->SetFillStyle(3008);

   PbPb276CTEQ61EPS09BFG2_sum_pdferr_InvYield->SetFillColor(kAzure+6);
   PbPb276CTEQ61EPS09BFG2_sum_pdferr_InvYield->SetFillStyle(3008);
   //pp276CT10BFG2_sum_pdferr_InvYield->RemovePoint(pp276CT10BFG2_sum_pdferr_InvYield->GetN()-1);
   //PbPb276CTEQ61EPS09BFG2_sum_pdferr_InvYield->RemovePoint(PbPb276CTEQ61EPS09BFG2_sum_pdferr_InvYield->GetN()-1);

   DrawGammaSetMarkerTGraphAsym(grSystFull , 20,2, colorcent, colorcent, 1, kTRUE);
   DrawGammaSetMarkerTGraphAsym(grSystFullA , 20,2, kGreen-8, kGreen-8, 1, kTRUE);
   DrawGammaSetMarkerTGraphAsym(grSystFullB , 20,2, kRed-8, kRed-8, 1, kTRUE);
   DrawGammaSetMarkerTGraphAsym(grSystFullC , 20,2, kBlue-8, kBlue-8, 1, kTRUE);
   grStatFull->SetMarkerColor(colorcent);
   grStatFull->SetLineColor(colorcent);

   if(!centralityCutNumber.CompareTo("601")){
      grSystFull->RemovePoint(0);
      grSystFull->RemovePoint(10);
      grStatFull->RemovePoint(10);
   }
   if(!centralityCutNumber.CompareTo("612")){
      grSystFull->RemovePoint(0);
      grStatFull->RemovePoint(0);

      grSystFull->RemovePoint(4);
      grStatFull->RemovePoint(4);
      grStatFull->RemovePoint(4);
      grStatFull->RemovePoint(4);
      grStatFull->RemovePoint(4);
      // grStatFull->RemovePoint(5);
      // grStatFull->RemovePoint(5);

      // grSystFull->RemovePoint(7);
      // grSystFull->RemovePoint(7);
      // grStatFull->RemovePoint(7);
      // grStatFull->RemovePoint(7);

   }
   if(!centralityCutNumber.CompareTo("502")){
      //grSystFull->RemovePoint(0);
      //grStatFull->RemovePoint(0);
      grSystFull->RemovePoint(16);
      grStatFull->RemovePoint(16);
   }
   if(!centralityCutNumber.CompareTo("504")){
      grStatFull->RemovePoint(0);
      grStatFull->RemovePoint(16);
   }
   if(!centralityCutNumber.CompareTo("501")){
      grStatFull->RemovePoint(0);
      grStatFull->RemovePoint(5);
      grStatFull->RemovePoint(6);
      grStatFull->RemovePoint(7);
      grStatFull->RemovePoint(7);
      grSystFull->RemovePoint(5);
      grSystFull->RemovePoint(6);
      grSystFull->RemovePoint(7);
      grSystFull->RemovePoint(12);
   }
   if(!centralityCutNumber.CompareTo("512")){

   }
   if(!centralityCutNumber.CompareTo("524")){
      grStatFull->RemovePoint(0);
      grStatFull->RemovePoint(0);
      grStatFull->RemovePoint(0);
      grStatFull->RemovePoint(0);
      grStatFull->RemovePoint(0);
      grStatFull->RemovePoint(0);
      grStatFull->RemovePoint(0);
      grStatFull->RemovePoint(0);
      grStatFull->RemovePoint(0);
      grStatFull->RemovePoint(1);
      grSystFull->RemovePoint(0);
      grSystFull->RemovePoint(1);

      graphArrowPointErrorBelow0->RemovePoint(0);
      graphArrowPointErrorBelow0->RemovePoint(0);
   }

   if(!centralityCutNumber.CompareTo("548")){
      graphArrowPointErrorBelow0->RemovePoint(0);
      graphArrowPointErrorBelow0->RemovePoint(0);
      graphArrowPointErrorBelow0->RemovePoint(0);
      graphArrowPointErrorBelow0->RemovePoint(0);
      graphArrowPointErrorBelow0->RemovePoint(0);

      grStatFull->RemovePoint(0);
      grStatFull->RemovePoint(0);
      grStatFull->RemovePoint(0);
      grStatFull->RemovePoint(0);
      grStatFull->RemovePoint(0);
      grSystFull->RemovePoint(0);


   }

   Double_t fitRangeLow = 1.1;
   Double_t fitRangeHigh = 2.3;

   TF1 *exponetial = new TF1("Exponential Fit","[0]*exp(-x/[1])",fitRangeLow,fitRangeHigh);
   TF1 *exponetialSyst = new TF1("Exponential Fit","[0]*exp(-x/[1])",fitRangeLow,fitRangeHigh);
   TF1 *exponetialStat = new TF1("Exponential Fit","[0]*exp(-x/[1])",fitRangeLow,fitRangeHigh);

   exponetial->SetParameters(100,0.300);
   exponetial->SetLineColor(colorComb4060);
   exponetialStat->SetParameters(100,0.300);
   exponetialStat->SetLineColor(colorComb4060);
   exponetialSyst->SetParameters(100,0.300);
   exponetialSyst->SetLineColor(colorComb4060);
   for(Int_t i = 0;i<grSystFullA->GetN();i++){
      grSystFullA->SetPointEXhigh(i,0);
      grSystFullA->SetPointEXlow(i,0);
   }
   for(Int_t i = 0;i<grStatFull->GetN();i++){
      grStatFull->SetPointEXhigh(i,0);
      grStatFull->SetPointEXlow(i,0);
   }
   if(grSystFullA)grSystFullA->Fit(exponetial,"NRE+","",fitRangeLow,fitRangeHigh);
   if(grStatFull)grStatFull->Fit(exponetialStat,"NRE+","",fitRangeLow,fitRangeHigh);
   //if(graAllShiftedSystWithoutMat)graAllShiftedSystWithoutMat->Fit(exponetialSyst,"NRE+","",fitRangeLow,fitRangeHigh);
   //graphUpperLimitForPointAboveButError->Fit(exponetialSyst,"NRE+","",fitRangeLow,fitRangeHigh);

   if(graphArrorErrorBelow0)graphArrorErrorBelow0->SetLineColor(colorcent);
   if(graphUpperLimitForPointAboveButError)graphUpperLimitForPointAboveButError->SetLineColor(colorcent);
   if(graphUpperLimitPointsAndErrorBelow0)graphUpperLimitPointsAndErrorBelow0->SetLineColor(colorcent);
   if(grSystFull)grSystFull->Draw("Z2");
   // if(grSystFullA)grSystFullA->Draw("Z2");
   // if(grSystFullB)grSystFullB->Draw("Z2");
   // if(grSystFullC)grSystFullC->Draw("Z2");


   if(grStatFull)grStatFull->Draw("PZe1");
   if(graphUpperLimitForPointAboveButError){
      graphUpperLimitForPointAboveButError->SetLineWidth(1);
      graphUpperLimitForPointAboveButError->SetLineColor(colorcent);
      graphUpperLimitForPointAboveButError->SetMarkerStyle(1);
      graphUpperLimitForPointAboveButError->Draw("||");
   }


   // if(graphArrowConfi){
   //    graphArrowConfi->SetLineColor(colorcent);
   //    graphArrowConfi->SetMarkerStyle(1);
   //    graphArrowConfi->Draw(">");
   // }
   // if(graphUpperLimitConfi){
   //    graphUpperLimitConfi->SetLineColor(colorcent);
   //    graphUpperLimitConfi->SetMarkerStyle(1);
   //    graphUpperLimitConfi->Draw("||");
   // }
   if(graphUpperLimitPointsAndErrorBelow0){
      graphUpperLimitPointsAndErrorBelow0->SetLineWidth(1);
      graphUpperLimitPointsAndErrorBelow0->SetLineColor(colorcent);
      graphUpperLimitPointsAndErrorBelow0->SetMarkerStyle(1);
      graphUpperLimitPointsAndErrorBelow0->Draw("||");
   }
   if(graphArrorErrorBelow0){
      graphArrorErrorBelow0->SetLineWidth(1);
      graphArrorErrorBelow0->SetLineColor(colorcent);
      graphArrorErrorBelow0->SetMarkerStyle(1);
      graphArrorErrorBelow0->Draw(">");
   }
   if(graphArrorErrorBelow0)graphArrorErrorBelow0->Draw("Z");
   if(graphArrowPointErrorBelow0){
      graphArrowPointErrorBelow0->SetLineWidth(1);
      graphArrowPointErrorBelow0->SetLineColor(colorcent);
      graphArrowPointErrorBelow0->SetMarkerStyle(1);
      graphArrowPointErrorBelow0->Draw(">");
   }


   if(grSystFull)grSystFull->Write(Form("DirectPhotons_points_syst_%s",centralityCutNumber.Data()),TObject::kOverwrite);
   if(grSystFullA)grSystFullA->Write(Form("DirectPhotons_points_systA_%s",centralityCutNumber.Data()),TObject::kOverwrite);
   if(grSystFullB)grSystFullB->Write(Form("DirectPhotons_points_systB_%s",centralityCutNumber.Data()),TObject::kOverwrite);
   if(grSystFullC)grSystFullC->Write(Form("DirectPhotons_points_systC_%s",centralityCutNumber.Data()),TObject::kOverwrite);
   if(grStatFull)grStatFull->Write(Form("DirectPhotons_points_stat_%s",centralityCutNumber.Data()));
   if(graphUpperLimitForPointAboveButError)graphUpperLimitForPointAboveButError->Write(Form("DirectPhotons_upperLimits_withPoint_%s",centralityCutNumber.Data()),TObject::kOverwrite);
   // //if(graphArrowConfi)graphArrowConfi->Draw(">");
   if(graphUpperLimitPointsAndErrorBelow0)graphUpperLimitPointsAndErrorBelow0->Write(Form("DirectPhotons_upperLimits_noPoint_%s",centralityCutNumber.Data()),TObject::kOverwrite);
   if(graphArrorErrorBelow0)graphArrorErrorBelow0->Write(Form("DirectPhotons_Arrow_withPoint_%s",centralityCutNumber.Data()),TObject::kOverwrite);
   // //if(graphArrorErrorBelow0)graphArrorErrorBelow0->Draw("Z");
   if(graphArrowPointErrorBelow0)graphArrowPointErrorBelow0->Write(Form("DirectPhotons_Arrow_noPoint_%s",centralityCutNumber.Data()),TObject::kOverwrite);
   if(NLO)NLO->Write(Form("NLO_%s",centralityCutNumber.Data()),TObject::kOverwrite);

   if(grSystFullCombined)grSystFullCombined->Write(Form("Combined_DirectPhotons_points_syst_%s",centralityCutNumber.Data()),TObject::kOverwrite);
   if(grStatFullCombined)grStatFullCombined->Write(Form("Combined_DirectPhotons_points_stat_%s",centralityCutNumber.Data()),TObject::kOverwrite);
   if(graphUpperLimitForPointAboveButErrorCombined)graphUpperLimitForPointAboveButErrorCombined->Write(Form("Combined_DirectPhotons_upperLimits_withPoint_%s",centralityCutNumber.Data()),TObject::kOverwrite);
   // //if(graphArrowConfi)graphArrowConfi->Draw(">");
   if(graphUpperLimitPointsAndErrorBelow0Combined)graphUpperLimitPointsAndErrorBelow0Combined->Write(Form("Combined_DirectPhotons_upperLimits_noPoint_%s",centralityCutNumber.Data()),TObject::kOverwrite);
   if(graphArrorErrorBelow0Combined)graphArrorErrorBelow0Combined->Write(Form("Combined_DirectPhotons_Arrow_withPoint_%s",centralityCutNumber.Data()),TObject::kOverwrite);
   // //if(graphArrorErrorBelow0)graphArrorErrorBelow0->Draw("Z");
   if(graphArrowPointErrorBelow0Combined)graphArrowPointErrorBelow0Combined->Write(Form("Combined_DirectPhotons_Arrow_noPoint_%s",centralityCutNumber.Data()),TObject::kOverwrite);


   FileGammaSpectrum->Close();

   TFile *FileCharged;
   if(isCharged){
      FileCharged = new TFile("GammaSpectrumCharged.root","UPDATE");
      if(grSystFull)grSystFull->Write("DirectPhotons_points_syst_charged");
      if(grStatFull)grStatFull->Write("DirectPhotons_points_stat_charged");
      if(graphUpperLimitForPointAboveButError)graphUpperLimitForPointAboveButError->Write("DirectPhotons_upperLimits_withPoint_charge");
      // //if(graphArrowConfi)graphArrowConfi->Draw(">");
      if(graphUpperLimitPointsAndErrorBelow0)graphUpperLimitPointsAndErrorBelow0->Write("DirectPhotons_upperLimits_noPoint_charged");
      if(graphArrorErrorBelow0)graphArrorErrorBelow0->Write("DirectPhotons_Arrow_withPoint_charged");
      // //if(graphArrorErrorBelow0)graphArrorErrorBelow0->Draw("Z");
      if(graphArrowPointErrorBelow0)graphArrowPointErrorBelow0->Write("DirectPhotons_Arrow_noPoint_charged");

   }

   if(drawCharged){
      TGraphAsymmErrors *grSystFullcharged = (TGraphAsymmErrors*)FileCharged->Get("DirectPhotons_points_syst_charged");
      TGraphAsymmErrors *grStatFullcharged = (TGraphAsymmErrors*)FileCharged->Get("DirectPhotons_points_stat_charged");
      TGraphAsymmErrors *graphUpperLimitForPointAboveButErrorcharged = (TGraphAsymmErrors*)FileCharged->Get("DirectPhotons_upperLimits_withPoint_charged");
      TGraphAsymmErrors *graphUpperLimitPointsAndErrorBelow0charged = (TGraphAsymmErrors*)FileCharged->Get("DirectPhotons_upperLimits_noPoint_charged");
      TGraphAsymmErrors *graphArrorErrorBelow0charged = (TGraphAsymmErrors*)FileCharged->Get("DirectPhotons_Arrow_withPoint_charged");
      TGraphAsymmErrors *graphArrowPointErrorBelow0charged = (TGraphAsymmErrors*)FileCharged->Get("DirectPhotons_Arrow_noPoint_charged");

      if(graphArrorErrorBelow0charged)graphArrorErrorBelow0charged->SetLineColor(1);
      if(graphUpperLimitForPointAboveButErrorcharged)graphUpperLimitForPointAboveButErrorcharged->SetLineColor(1);
      if(graphUpperLimitPointsAndErrorBelow0charged)graphUpperLimitPointsAndErrorBelow0charged->SetLineColor(1);
      if(graphUpperLimitPointsAndErrorBelow0charged)grStatFullcharged->SetLineColor(1);
      if(graphUpperLimitPointsAndErrorBelow0charged)grStatFullcharged->SetMarkerColor(1);
      if(graphUpperLimitPointsAndErrorBelow0charged)grStatFullcharged->SetFillColor(0);
      if(graphArrorErrorBelow0charged)graphArrorErrorBelow0charged->SetLineColor(1);
      if(grStatFullcharged)grStatFullcharged->SetMarkerStyle(20);
      //if(grSystFullcharged)grSystFullcharged->Draw("Z2");
      if(grStatFullcharged)grStatFullcharged->Draw("PZe1");
      if(graphUpperLimitForPointAboveButErrorcharged)graphUpperLimitForPointAboveButErrorcharged->Draw("||");
      //if(graphArrowConfi)graphArrowConfi->Draw(">");
      if(graphUpperLimitPointsAndErrorBelow0charged)graphUpperLimitPointsAndErrorBelow0charged->Draw("||");
      //if(graphArrorErrorBelow0charged)graphArrorErrorBelow0charged->Draw(">");
      //if(graphArrorErrorBelow0)graphArrorErrorBelow0->Draw("Z");
      if(graphArrowPointErrorBelow0charged)graphArrowPointErrorBelow0charged->Draw(">");
   }

   if(graAllShifted)graAllShifted->SetFillColor(kGray);
   if(graAllShifted)graAllShifted->SetLineColor(0);

   //AliHEPDataParser *HEPData = new AliHEPDataParser(grStatFull,grSystFull);

   //HEPData->SaveHEPDataFile("DirectPhotons_2760GeV_00-40.hep",1);

   cout<<centralityCutNumber<<endl;
   if(!(!centralityCutNumber.CompareTo("524") || !centralityCutNumber.CompareTo("548") || !centralityCutNumber.CompareTo("612")))
      exponetial->Draw("same");


   //exponetialPrelim->Draw("same");
   TLegend* legendDirectPhotonWithXshift = GetAndSetLegend(0.33,0.58,8.0);
   legendDirectPhotonWithXshift->AddEntry(GraphErrorsDummy,Form("#gamma_{dir} ALICE"),"pf");
   //legendDirectPhotonWithXshift->AddEntry(grSystFullp,Form("Prelim. #gamma_{dir} ALICE"),"pf");
   //legendDirectPhotonWithXshift->AddEntry(grStatFullcharged,Form("#gamma_{dir} with #pi^{#pm} ALICE"),"pe");
   legendDirectPhotonWithXshift->AddEntry(NLO,"PDF: CTEQ6M5 FF: DSS","l");
   legendDirectPhotonWithXshift->AddEntry((TObject*)0,"JETPHOX:","");
   legendDirectPhotonWithXshift->AddEntry(pp276CT10BFG2_sum_pdferr_InvYield,"PDF: CT10, FF: BFG2","f");
   legendDirectPhotonWithXshift->AddEntry(PbPb276CTEQ61EPS09BFG2_sum_pdferr_InvYield,"nPDF: EPS09, FF: BFG2","f");
   legendDirectPhotonWithXshift->AddEntry((TObject*)0,"(all N_{coll} scaled pp)","");
   if(!(!centralityCutNumber.CompareTo("524") || !centralityCutNumber.CompareTo("548") || !centralityCutNumber.CompareTo("612"))){
      legendDirectPhotonWithXshift->AddEntry(exponetial,Form("exp(-#it{p}_{T}/T) (%.1f < #it{p}_{T} < %.1f GeV/#it{c})",fitRangeLow,fitRangeHigh),"l");
      legendDirectPhotonWithXshift->AddEntry((TObject*)0,Form("T = %.0f #pm %.0f^{stat} #pm %0.0f^{syst} MeV",exponetialStat->GetParameter(1)*1000,exponetialStat->GetParError(1)*1000,exponetial->GetParError(1)*1000),"");
   }

   legendDirectPhotonWithXshift->Draw();

   PbPb276CTEQ61EPS09BFG2_sum_pdferr_InvYield->Draw("E3");
   pp276CT10BFG2_sum_pdferr_InvYield->Draw("E3");
   NLO->Draw("E3");

   if(graphArrorErrorBelow0)graphArrorErrorBelow0->SetLineColor(colorcent);
   if(graphUpperLimitForPointAboveButError)graphUpperLimitForPointAboveButError->SetLineColor(colorcent);
   if(graphUpperLimitPointsAndErrorBelow0)graphUpperLimitPointsAndErrorBelow0->SetLineColor(colorcent);
   if(grSystFull)grSystFull->Draw("Z2");
   if(grStatFull)grStatFull->Draw("PZe1");
   if(graphUpperLimitForPointAboveButError){
      //graphUpperLimitForPointAboveButError->Draw("||");
   }
   if(graphUpperLimitPointsAndErrorBelow0){
      graphUpperLimitPointsAndErrorBelow0->Draw("||");
   }
   if(graphArrorErrorBelow0){
      graphArrorErrorBelow0->Draw(">");
   }
   if(graphArrorErrorBelow0)graphArrorErrorBelow0->Draw("Z");
   if(graphArrowPointErrorBelow0){
      graphArrowPointErrorBelow0->Draw(">");
   }


   PbPb276CTEQ61EPS09BFG2_sum_pdferr_InvYield->SetFillColor(kAzure+6);

   DrawSystem(0.19,0.15,option.Data(),(GetCentralityStringC(cutSel)).Data());


   gPad->RedrawAxis();   
   
   canvasDirectPhotonWithShift->Print(Form("%s/DirectPhotonSpectrum_%s.eps",outputDir.Data(),cent.Data()));
   canvasDirectPhotonWithShift->Print(Form("%s/DirectPhotonSpectrum_%s.C",outputDir.Data(),cent.Data()));



   TGraphAsymmErrors *grSystFullp = (TGraphAsymmErrors*)PrelimFile->Get("DirectPhotons_points_syst");
   TGraphAsymmErrors *grStatFullp = (TGraphAsymmErrors*)PrelimFile->Get("DirectPhotons_points_stat");
   TGraphAsymmErrors *graphUpperLimitForPointAboveButErrorp = (TGraphAsymmErrors*)PrelimFile->Get("DirectPhotons_upperLimits_withPoint");
   TGraphAsymmErrors *graphUpperLimitPointsAndErrorBelow0p = (TGraphAsymmErrors*)PrelimFile->Get("DirectPhotons_upperLimits_noPoint");
   TGraphAsymmErrors *graphArrorErrorBelow0p = (TGraphAsymmErrors*)PrelimFile->Get("DirectPhotons_Arrow_withPoint");
   TGraphAsymmErrors *graphArrowPointErrorBelow0p = (TGraphAsymmErrors*)PrelimFile->Get("DirectPhotons_Arrow_noPoint");
   TF1 *exponetialPrelim = new TF1("Exponential Fit","[0]*exp(-x/[1])",0.8,2.2);



   if(!centralityCutNumber.CompareTo("504")){

      //dummy->DrawCopy();

      exponetialPrelim->SetParameters(100,0.300);
      exponetialPrelim->SetLineColor(colorComb0020);
      exponetialPrelim->FixParameter(1,0.304);
      for(Int_t i = 0;i<grStatFullp->GetN();i++){
         grStatFullp->SetPointEXhigh(i,0);
         grStatFullp->SetPointEXlow(i,0);
      }
      grStatFullp->Fit(exponetialPrelim,"NRE+","",0.8,2.2);
      if(grSystFull)grSystFullp->Draw("Z2");
      PbPb276CTEQ61EPS09BFG2_sum_pdferr_InvYield->Draw("E3");
      pp276CT10BFG2_sum_pdferr_InvYield->Draw("E3");
      NLO->Draw("3l");
      if(grStatFull)grStatFullp->Draw("PZe1");
      if(graphUpperLimitForPointAboveButError)graphUpperLimitForPointAboveButErrorp->Draw("||");
      //if(graphArrowConfi)graphArrowConfi->Draw(">");
      if(graphUpperLimitPointsAndErrorBelow0)graphUpperLimitPointsAndErrorBelow0p->Draw("||");
      if(graphArrorErrorBelow0)graphArrorErrorBelow0p->Draw(">");
      //if(graphArrorErrorBelow0)graphArrorErrorBelow0->Draw("Z");
      if(graphArrowPointErrorBelow0)graphArrowPointErrorBelow0p->Draw(">");


     
      if(grSystFull)grSystFull->Draw("Z2");
      if(grStatFull){
         grStatFull->RemovePoint(16);
         grStatFull->Draw("PZe1");
      }
      if(graphUpperLimitForPointAboveButError){
         graphUpperLimitForPointAboveButError->Draw("||");
      }
      if(graphUpperLimitPointsAndErrorBelow0){
         graphUpperLimitPointsAndErrorBelow0->Draw("||");
      }
      if(graphArrorErrorBelow0){
         graphArrorErrorBelow0->Draw(">");
      }
      if(graphArrorErrorBelow0)graphArrorErrorBelow0->Draw("Z");
      if(graphArrowPointErrorBelow0){
         graphArrowPointErrorBelow0->Draw(">");
      }
      grSystFullp->SetLineColor(0);
      exponetial->Draw("same");
      legendDirectPhotonWithXshift->AddEntry(grSystFullp,Form("Prelim. #gamma_{dir} ALICE"),"pf");

      legendDirectPhotonWithXshift->Draw();

      
      TF1* directPhotonFit = FitObject("h","Name","Pi0",grStatFullp,.8,5.,NULL,"QNRME+");

      PrintParameter(directPhotonFit,"h","Name",6,"blup");
      grStatFullp->Fit(directPhotonFit);
grStatFullp->Fit(directPhotonFit);
grStatFullp->Fit(directPhotonFit);
grStatFullp->Fit(directPhotonFit);
grStatFullp->Fit(directPhotonFit);
grStatFullp->Fit(directPhotonFit);
grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);
grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);
grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);
grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);
grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);
grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);
grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);


 grStatFullp->Fit(directPhotonFit);grStatFullp->Fit(directPhotonFit);
 
      //Create a histogram to hold the confidence intervals
      TH1D *hint = new TH1D("hint", 
                            "Fitted gaussian with .95 conf.band", 1000,0.8 ,5. );
      (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint);
      cout<<TVirtualFitter::GetFitter()<<endl;
      //Now the "hint" histogram has the fitted function values as the 
      //bin contents and the confidence intervals as bin errors
      hint->SetStats(kFALSE);
      hint->SetMarkerSize(0);
      hint->SetFillColor(4);
      hint->Draw("samee3");
      //directPhotonFit->Draw("same");
      
      
      canvasDirectPhotonWithShift->Print(Form("%s/DirectPhotonSpectrumPlusPrelim_%s.eps",outputDir.Data(),cent.Data()));
      canvasDirectPhotonWithShift->Print(Form("%s/DirectPhotonSpectrumPlusPrelim_%s.C",outputDir.Data(),cent.Data()));


      TCanvas *FitToSpectrum = GetAndSetCanvas("FitToSpectrum");
      dummy->GetXaxis()->SetRangeUser(0.0,5);      
      dummy->DrawCopy();



      TGraphAsymmErrors *ratio = grStatFullp->Clone("ratio");
      
      for(Int_t i = 0;i<=grStatFullp->GetN(); i++){
         grStatFullp->SetPoint(i,grStatFullp->GetX()[i],grStatFullp->GetY()[i]/directPhotonFit->Eval(grStatFullp->GetX()[i]));
         grStatFullp->SetPointError(i,0,0,grStatFullp->GetEYlow()[i]/directPhotonFit->Eval(grStatFullp->GetX()[i]),grStatFullp->GetEYhigh()[i]/directPhotonFit->Eval(grStatFullp->GetX()[i]));
         if(grStatFullp->GetX()[i]>5)grStatFullp->RemovePoint(i);
         
      }
      grStatFullp->RemovePoint(grStatFullp->GetN()-1);
      grStatFullp->Draw("ap");
      
      
      One->Draw("same");
      
      FitToSpectrum->Print("SpectrumToFitRatio.pdf");
   }



   // *********************** Subtracted NLO Plots ***************************************

   TCanvas *canvasDirectPhotonNLOSubtracted = GetAndSetCanvas("canvasDirectPhotonNLOSubtracted",0.15,0.1, 1000 ,1250);
   canvasDirectPhotonNLOSubtracted->SetLogy();
   dummy->GetYaxis()->SetRangeUser(2e-6,2e1);
   dummy->GetXaxis()->SetRangeUser(0.80,4.);
   dummy->GetYaxis()->SetTitleOffset(1.5);
   dummy->GetXaxis()->SetLabelOffset(0);
   dummy->DrawCopy();



   if(grStatFullSubNLO && grSystFullSubNLO){
      if(!centralityCutNumber.CompareTo("601")){
         // grStatFullSubNLO->RemovePoint(5);
         // grStatFullSubNLO->RemovePoint(6);
         // grStatFullSubNLO->RemovePoint(6);
         // grStatFullSubNLO->RemovePoint(6);
         // grSystFullSubNLO->RemovePoint(0);
         // grSystFullSubNLO->RemovePoint(5);
      }
      if(!centralityCutNumber.CompareTo("612")){
         grStatFullSubNLO->RemovePoint(0);
         grSystFullSubNLO->RemovePoint(0);
         grStatFullSubNLO->RemovePoint(4);
         //grStatFullSubNLO->RemovePoint(3);
         //grSystFullSubNLO->RemovePoint(3);
         grSystFullSubNLO->RemovePoint(4);
         // grSystFullSubNLO->RemovePoint(3);


      }
      if(!centralityCutNumber.CompareTo("502")){
         grStatFullSubNLO->RemovePoint(0);
         grStatFullSubNLO->RemovePoint(7);
         //grSystFullSubNLO->RemovePoint(7);
         grStatFullSubNLO->RemovePoint(8);
         grStatFullSubNLO->RemovePoint(8);
         grStatFullSubNLO->RemovePoint(8);
         //grStatFullSubNLO->RemovePoint(9);
         // grStatFullSubNLO->RemovePoint(10);
         grSystFullSubNLO->RemovePoint(7);
      }
      if(!centralityCutNumber.CompareTo("504")){
         grStatFullSubNLO->RemovePoint(0);
         grStatFullSubNLO->RemovePoint(3);
         grStatFullSubNLO->RemovePoint(10);
         grSystFullSubNLO->RemovePoint(3);

      }
      if(!centralityCutNumber.CompareTo("501")){
         //grStatFullSubNLO->RemovePoint(0);
         grStatFullSubNLO->RemovePoint(8);
         grStatFullSubNLO->RemovePoint(8);
         grStatFullSubNLO->RemovePoint(8);
         grSystFullSubNLO->RemovePoint(8);
         grSystFullSubNLO->RemovePoint(8);
         grSystFullSubNLO->RemovePoint(8);

         // grStatFullSubNLO->RemovePoint(7);
         // grStatFullSubNLO->RemovePoint(7);
         // grStatFullSubNLO->RemovePoint(7);
         //grSystFullSubNLO->RemovePoint(0);
      }
      if(!centralityCutNumber.CompareTo("512")){
         grStatFullSubNLO->RemovePoint(9);
         grStatFullSubNLO->RemovePoint(10);
         // grStatFullSubNLO->RemovePoint(5);
         // grStatFullSubNLO->RemovePoint(5);
         // grStatFullSubNLO->RemovePoint(6);
         // grSystFullSubNLO->RemovePoint(3);
         // grSystFullSubNLO->RemovePoint(5);
      }
      if(!centralityCutNumber.CompareTo("524")){
      }
   }

   if(graphArrorErrorBelow0SubNLO)graphArrorErrorBelow0SubNLO->SetLineColor(colorcent);
   if(graphUpperLimitForPointAboveButErrorSubNLO)graphUpperLimitForPointAboveButErrorSubNLO->SetLineColor(colorcent);
   if(graphUpperLimitPointsAndErrorBelow0SubNLO)graphUpperLimitPointsAndErrorBelow0SubNLO->SetLineColor(colorcent);
   if(grSystFullSubNLO){
      DrawGammaSetMarkerTGraphAsym(grSystFullSubNLO , 20,1, colorcent, colorcent, 1, kTRUE);
      grSystFullSubNLO->Draw("Z2");
   }
   if(grStatFullSubNLO){
      grStatFullSubNLO->SetMarkerColor(colorcent);
      grStatFullSubNLO->SetLineColor(colorcent);
      grStatFullSubNLO->Draw("PZe1");
   }
   if(graphUpperLimitForPointAboveButErrorSubNLO){
      graphUpperLimitForPointAboveButErrorSubNLO->SetLineColor(colorcent);
      graphUpperLimitForPointAboveButErrorSubNLO->SetMarkerStyle(1);
      graphUpperLimitForPointAboveButErrorSubNLO->Draw("||");
   }
   if(graphUpperLimitPointsAndErrorBelow0SubNLO){
      graphUpperLimitPointsAndErrorBelow0SubNLO->SetLineColor(colorcent);
      graphUpperLimitPointsAndErrorBelow0SubNLO->SetMarkerStyle(1);
      graphUpperLimitPointsAndErrorBelow0SubNLO->Draw("||");
   }
   if(graphArrorErrorBelow0SubNLO){
      graphArrorErrorBelow0SubNLO->SetLineColor(colorcent);
      graphArrorErrorBelow0SubNLO->SetMarkerStyle(1);
      graphArrorErrorBelow0SubNLO->Draw(">");
   }
   if(graphArrorErrorBelow0SubNLO)graphArrorErrorBelow0SubNLO->Draw("Z");
   if(graphArrowPointErrorBelow0SubNLO){
      graphArrowPointErrorBelow0SubNLO->SetLineColor(colorcent);
      graphArrowPointErrorBelow0SubNLO->SetMarkerStyle(1);
      //graphArrowPointErrorBelow0SubNLO->Draw(">");
   }


   TF1 *exponetialSubNLOSyst = new TF1("Exponential Fit","[0]*exp(-x/[1])",fitRangeLow,fitRangeHigh);
   TF1 *exponetialSubNLOStat = new TF1("Exponential Fit","[0]*exp(-x/[1])",fitRangeLow,fitRangeHigh);

   exponetialSubNLOSyst->SetParameters(100,0.300);
   exponetialSubNLOStat->SetParameters(100,0.300);
   if(grSystFullSubNLO)grSystFullSubNLO->Fit(exponetialSubNLOSyst,"NRE+","",fitRangeLow,fitRangeHigh);
   if(grStatFullSubNLO)grStatFullSubNLO->Fit(exponetialSubNLOStat,"NRE+","",fitRangeLow,fitRangeHigh);

   exponetialSubNLOStat->SetLineColor(colorComb4060);
   //exponetialSubNLOSyst->Draw("same");
   exponetialSubNLOStat->Draw("same");



   if(grStatFullSubNLO)grSystFullSubNLO->SetMarkerSize(2);
   TLegend* legendDirectPhotonNLOSubtracted = GetAndSetLegend(0.17,0.13,6.5);
   legendDirectPhotonNLOSubtracted->AddEntry(grSystFullSubNLO,Form("#gamma_{dir} ALICE"),"pf");
   legendDirectPhotonNLOSubtracted->AddEntry((TObject*)0,"NLO subtracted","");
   legendDirectPhotonNLOSubtracted->AddEntry((TObject*)0,"(nPDF: EPS09, FF: BFG2)","");
   if(!(!centralityCutNumber.CompareTo("548"))){
      legendDirectPhotonNLOSubtracted->AddEntry(exponetialSubNLOStat,Form("exp(-#it{p}_{T}/T) (%.1f < #it{p}_{T} < %.1f GeV/#it{c})",fitRangeLow,fitRangeHigh),"l");
      legendDirectPhotonNLOSubtracted->AddEntry((TObject*)0,Form("T = %.0f #pm %.0f^{stat} #pm %0.0f^{syst} MeV",exponetialSubNLOStat->GetParameter(1)*1000,exponetialSubNLOStat->GetParError(1)*1000,exponetialSubNLOSyst->GetParError(1)*1000),"");
   }
   legendDirectPhotonNLOSubtracted->Draw();

   DrawSystem(0.25,0.90,option.Data(),(GetCentralityStringC(cutSel)).Data());

   canvasDirectPhotonNLOSubtracted->Print(Form("%s/DirectPhotonSpectrumSubtractedNLO_%s.eps",outputDir.Data(),cent.Data()));
   canvasDirectPhotonNLOSubtracted->Print(Form("%s/DirectPhotonSpectrumSubtractedNLO_%s.C",outputDir.Data(),cent.Data()));









   TCanvas *canvasDirectPhotonCombined = GetAndSetCanvas("canvasDirectPhotonCombined",0.15,0.1, 1100 ,1250);
   canvasDirectPhotonCombined->SetLogy();
   dummy->GetYaxis()->SetRangeUser(1e-7,4e1);
   if(!centralityCutNumber.CompareTo("548")){
      dummy->GetYaxis()->SetRangeUser(2e-10,2e1);
   }
   dummy->GetYaxis()->SetTitle("#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c}^{2})");
   dummy->GetXaxis()->SetRangeUser(0.50,11.4);
   dummy->GetYaxis()->SetTitleOffset(1.5);
   dummy->GetXaxis()->SetLabelOffset(0);
   dummy->DrawCopy();


   PbPb276CTEQ61EPS09BFG2_sum_pdferr_InvYield->Draw("E3");
   pp276CT10BFG2_sum_pdferr_InvYield->Draw("E3");
   NLO->Draw("E3");

   
   if(grSystFullCombined){
      DrawGammaSetMarkerTGraphAsym(grSystFullCombined , 20,2, colorcent, colorcent, 1, kTRUE);
      // if(!centralityCutNumber.CompareTo("548"))grSystFullCombined->RemovePoint(0);
      // if(!centralityCutNumber.CompareTo("548"))grSystFullCombined->RemovePoint(0);
      if(!centralityCutNumber.CompareTo("548"))grSystFullCombined->RemovePoint(12);
      if(!centralityCutNumber.CompareTo("504"))grSystFullCombined->RemovePoint(16);
      if(!centralityCutNumber.CompareTo("502"))grSystFullCombined->RemovePoint(0);
      if(!centralityCutNumber.CompareTo("501"))grSystFullCombined->RemovePoint(0);
      if(!centralityCutNumber.CompareTo("512"))grSystFullCombined->RemovePoint(16);
      grSystFullCombined->Draw("Z2");
   }
   if(grStatFullCombined){
      DrawGammaSetMarkerTGraphAsym(grStatFullCombined , 20,2, colorcent, colorcent, 1, kTRUE);
      if(!centralityCutNumber.CompareTo("524"))grStatFullCombined->RemovePoint(0);
      if(!centralityCutNumber.CompareTo("548"))grStatFullCombined->RemovePoint(0);
      if(!centralityCutNumber.CompareTo("548"))grStatFullCombined->RemovePoint(0);
      if(!centralityCutNumber.CompareTo("548"))grStatFullCombined->RemovePoint(0);
      if(!centralityCutNumber.CompareTo("548"))grStatFullCombined->RemovePoint(12);
      if(!centralityCutNumber.CompareTo("502"))grStatFullCombined->RemovePoint(0);
      if(!centralityCutNumber.CompareTo("504"))grStatFullCombined->RemovePoint(0);
      if(!centralityCutNumber.CompareTo("504"))grStatFullCombined->RemovePoint(16);
      if(!centralityCutNumber.CompareTo("501"))grStatFullCombined->RemovePoint(0);
      if(!centralityCutNumber.CompareTo("512"))grStatFullCombined->RemovePoint(16);
      grStatFullCombined->Draw("PZe1");
   }
   if(graphUpperLimitForPointAboveButErrorCombined){
      graphUpperLimitForPointAboveButErrorCombined->SetLineWidth(1);
      graphUpperLimitForPointAboveButErrorCombined->SetLineColor(colorcent);
      graphUpperLimitForPointAboveButErrorCombined->SetMarkerStyle(1);
      graphUpperLimitForPointAboveButErrorCombined->Draw("||");
   }
   if(graphUpperLimitPointsAndErrorBelow0Combined){
      DrawGammaSetMarkerTGraphAsym( graphUpperLimitPointsAndErrorBelow0Combined, 20,2, colorcent, colorcent, 1, kTRUE);
      graphUpperLimitPointsAndErrorBelow0Combined->Draw("||");
   }
   if(graphArrorErrorBelow0Combined){
      graphArrorErrorBelow0Combined->SetLineColor(colorcent);
      graphArrorErrorBelow0Combined->SetMarkerStyle(1);
      graphArrorErrorBelow0Combined->Draw(">");
   }
   if(graphArrorErrorBelow0Combined){
      graphArrorErrorBelow0Combined->SetLineColor(colorcent);
      graphArrorErrorBelow0Combined->SetMarkerStyle(1);
      graphArrorErrorBelow0Combined->Draw("Z");
   }
   if(graphArrowPointErrorBelow0Combined){
      DrawGammaSetMarkerTGraphAsym(graphArrowPointErrorBelow0Combined , 20,2, colorcent, colorcent, 1, kTRUE);
      if(!centralityCutNumber.CompareTo("524"))graphArrowPointErrorBelow0Combined->RemovePoint(0);
      if(!centralityCutNumber.CompareTo("548"))graphArrowPointErrorBelow0Combined->RemovePoint(0);
      if(!centralityCutNumber.CompareTo("548"))graphArrowPointErrorBelow0Combined->RemovePoint(0);
      if(!centralityCutNumber.CompareTo("548"))graphArrowPointErrorBelow0Combined->RemovePoint(0);
      graphArrowPointErrorBelow0Combined->Draw(">");
   }

   // fitRangeLow = 0.9;
   // fitRangeHigh = 3.0;

   TF1 *exponetialCombinedSyst = new TF1("Exponential Fit","[0]*exp(-x/[1])",fitRangeLow,fitRangeHigh);
   TF1 *exponetialCombinedStat = new TF1("Exponential Fit","[0]*exp(-x/[1])",fitRangeLow,fitRangeHigh);

   exponetialCombinedSyst->SetParameters(100,0.300);
   exponetialCombinedStat->SetParameters(100,0.300);
   for(Int_t i = 0;i<grSystFullCombinedA->GetN();i++){
      grSystFullCombinedA->SetPointEXhigh(i,0);
      grSystFullCombinedA->SetPointEXlow(i,0);
   }
   for(Int_t i = 0;i<grStatFullCombined->GetN();i++){
      grStatFullCombined->SetPointEXhigh(i,0);
      grStatFullCombined->SetPointEXlow(i,0);
   }

   grSystFullCombinedA->Fit(exponetialCombinedSyst,"NRE+","",fitRangeLow,fitRangeHigh);
   grStatFullCombined->Fit(exponetialCombinedStat,"NRE+","",fitRangeLow,fitRangeHigh);

   TH1D *hintB = new TH1D("hint", 
                         "Fitted gaussian with .95 conf.band", 1000,fitRangeLow ,fitRangeHigh );
   (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hintB,0.999);
   cout<<TVirtualFitter::GetFitter()<<endl;
   //Now the "hint" histogram has the fitted function values as the 
   //bin contents and the confidence intervals as bin errors
   hintB->SetStats(kFALSE);
   hintB->SetMarkerSize(0);
   hintB->SetFillColor(kYellow);
   hintB->Draw("samee3");
   
   exponetialCombinedStat->SetLineColor(colorComb4060);
   exponetialCombinedStat->Draw("same");


   TLegend* legendDirectPhotonCombined = GetAndSetLegend(0.33,0.6,8.0);
   legendDirectPhotonCombined->AddEntry(GraphErrorsDummy,Form("#gamma_{dir} ALICE"),"pf");
   //legendDirectPhotonCombined->AddEntry(grSystFullp,Form("Prelim. #gamma_{dir} ALICE"),"pf");
   //legendDirectPhotonCombined->AddEntry(grStatFullcharged,Form("#gamma_{dir} with #pi^{#pm} ALICE"),"pe");
   legendDirectPhotonCombined->AddEntry(NLO,"PDF: CTEQ6M5 FF: DSS","l");
   legendDirectPhotonCombined->AddEntry((TObject*)0,"JETPHOX:","");
   legendDirectPhotonCombined->AddEntry(pp276CT10BFG2_sum_pdferr_InvYield,"PDF: CT10, FF: BFG2","f");
   legendDirectPhotonCombined->AddEntry(PbPb276CTEQ61EPS09BFG2_sum_pdferr_InvYield,"nPDF: EPS09, FF: BFG2","f");
   legendDirectPhotonCombined->AddEntry((TObject*)0,"(all N_{coll} scaled pp)","");
   // if(!(!centralityCutNumber.CompareTo("524") || !centralityCutNumber.CompareTo("548") || !centralityCutNumber.CompareTo("612"))){
   legendDirectPhotonCombined->AddEntry(exponetial,Form("exp(-#it{p}_{T}/T) (%.1f < #it{p}_{T} < %.1f GeV/#it{c})",fitRangeLow,fitRangeHigh),"l");
   legendDirectPhotonCombined->AddEntry((TObject*)0,Form("T = %.0f #pm %.0f^{stat} #pm %0.0f^{syst} MeV",exponetialCombinedStat->GetParameter(1)*1000,exponetialCombinedStat->GetParError(1)*1000,exponetialCombinedSyst->GetParError(1)*1000),"");
   //}
   legendDirectPhotonCombined->Draw();

   if(grSystFullCombined)
      grSystFullCombined->Draw("Z2");
   if(grStatFullCombined)
      grStatFullCombined->Draw("PZe1");
   if(graphUpperLimitForPointAboveButErrorCombined)
      graphUpperLimitForPointAboveButErrorCombined->Draw("||");
   if(graphUpperLimitPointsAndErrorBelow0Combined)
      graphUpperLimitPointsAndErrorBelow0Combined->Draw("||");
   if(graphArrorErrorBelow0Combined)
      graphArrorErrorBelow0Combined->Draw(">");
   if(graphArrorErrorBelow0Combined)
      graphArrorErrorBelow0Combined->Draw("Z");
   if(graphArrowPointErrorBelow0Combined)
      graphArrowPointErrorBelow0Combined->Draw(">");

   gPad->RedrawAxis();
   canvasDirectPhotonCombined->Print(Form("%s/DirectPhotonSpectrumCombined_%s.eps",outputDir.Data(),cent.Data()));

   // Ratio NLO / Combined

   TCanvas *canvasDirectPhotonNLORaa = GetAndSetCanvas("canvasDirectPhotonNLORaa");
   //canvasDirectPhotonNLORaa->SetLogy();
   dummy->GetYaxis()->SetRangeUser(0,12);
   dummy->GetXaxis()->SetRangeUser(0.,15);
   dummy->GetYaxis()->SetTitleOffset(.8);
   dummy->GetXaxis()->SetLabelOffset(0);
   dummy->GetYaxis()->SetTitle("#it{R}_{AA,NLO}^{#gamma_{dir}}");
   dummy->DrawCopy();




   TGraphAsymmErrors *grSpectrumToNLORatioSyst = (TGraphAsymmErrors*)grSystFullCombined->Clone("grSystFullCombinedSubNLO");
   TGraphAsymmErrors *grSpectrumToNLORatioStat = (TGraphAsymmErrors*)grStatFullCombined->Clone("grSystFullCombinedSubNLO");

   if(!centralityCutNumber.CompareTo("502") || !centralityCutNumber.CompareTo("501") ||  !centralityCutNumber.CompareTo("504")){
      grSpectrumToNLORatioStat->RemovePoint(0);
      grSpectrumToNLORatioSyst->RemovePoint(0);
      grSpectrumToNLORatioStat->RemovePoint(0);
      grSpectrumToNLORatioSyst->RemovePoint(0);
      //grSpectrumToNLORatioStat->RemovePoint(0);
      //grSpectrumToNLORatioSyst->RemovePoint(0);
   }
   if(!centralityCutNumber.CompareTo("524")){
      grSpectrumToNLORatioStat->RemovePoint(0);
      grSpectrumToNLORatioSyst->RemovePoint(0);
   }
   if(!centralityCutNumber.CompareTo("512")){
      grSpectrumToNLORatioStat->RemovePoint(0);
      grSpectrumToNLORatioStat->RemovePoint(0);
      //grSpectrumToNLORatioStat->RemovePoint(0);
      grSpectrumToNLORatioSyst->RemovePoint(0);
      grSpectrumToNLORatioSyst->RemovePoint(0);
   }


   for(Int_t i = 0;i<grSpectrumToNLORatioSyst->GetN();i++){
      Double_t x1,x2,y1,y2;
      grSpectrumToNLORatioSyst->GetPoint(i,x1,y1);
      Double_t xerrlow = grSpectrumToNLORatioSyst->GetErrorXlow(i);
      Double_t xerrhigh = grSpectrumToNLORatioSyst->GetErrorXhigh(i);
      Double_t yerrlow = grSpectrumToNLORatioSyst->GetErrorYlow(i);
      Double_t yerrhigh = grSpectrumToNLORatioSyst->GetErrorYhigh(i);
      for(Int_t j = 0;j<PbPb276CTEQ61EPS09BFG2_sum_pdferr_InvYield->GetN();j++){
         PbPb276CTEQ61EPS09BFG2_sum_pdferr_InvYield->GetPoint(j,x2,y2);
         if(x2>x1-xerrlow && x2<x1+xerrhigh){
            Double_t yerrlowNLO = PbPb276CTEQ61EPS09BFG2_sum_pdferr_InvYield->GetErrorYlow(j);
            Double_t yerrhighNLO = PbPb276CTEQ61EPS09BFG2_sum_pdferr_InvYield->GetErrorYhigh(j);
            grSpectrumToNLORatioSyst->SetPoint(i,x1,y1/y2);
            grSpectrumToNLORatioSyst->SetPointError(i,xerrlow,xerrhigh,sqrt( (yerrlow*yerrlow*y2*y2+yerrlowNLO*yerrlowNLO*y1*y1)/ (y2*y2*y2*y2) ) ,sqrt( (yerrhigh*yerrhigh*y2*y2+yerrhighNLO*yerrhighNLO*y1*y1)/ (y2*y2*y2*y2) ) );
         }
      }
   }

   for(Int_t i = 0;i<grSpectrumToNLORatioStat->GetN();i++){
      Double_t x1,x2,y1,y2;
      grSpectrumToNLORatioStat->GetPoint(i,x1,y1);
      Double_t xerrlow = grSpectrumToNLORatioSyst->GetErrorXlow(i);
      Double_t xerrhigh = grSpectrumToNLORatioSyst->GetErrorXhigh(i);
      Double_t yerrlow = grSpectrumToNLORatioStat->GetErrorYlow(i);
      Double_t yerrhigh = grSpectrumToNLORatioStat->GetErrorYhigh(i);
      for(Int_t j = 0;j<PbPb276CTEQ61EPS09BFG2_sum_pdferr_InvYield->GetN();j++){
         PbPb276CTEQ61EPS09BFG2_sum_pdferr_InvYield->GetPoint(j,x2,y2);
         if(x2>x1-xerrlow && x2<x1+xerrhigh){
            cout<<x1-xerrlow<<"  "<<x1+xerrhigh<<endl;
            Double_t yerrlowNLO = 0;
            Double_t yerrhighNLO = 0;
            grSpectrumToNLORatioStat->SetPoint(i,x1,y1/y2);
            //cout<<x1<<"  "<<x2<<"  "<<y1<<"  "<<y2<<endl;
            grSpectrumToNLORatioStat->SetPointError(i,0,0,sqrt( (yerrlow*yerrlow*y2*y2+yerrlowNLO*yerrlowNLO*y1*y1)/ (y2*y2*y2*y2) ) ,sqrt( (yerrhigh*yerrhigh*y2*y2+yerrhighNLO*yerrhighNLO*y1*y1)/ (y2*y2*y2*y2) ) );
         }
      }

   }

   ONE2->Draw("same");
   grSpectrumToNLORatioSyst->Draw("Z2");
   grSpectrumToNLORatioStat->Draw("PZe1");

   DrawSystem(0.3,0.86,option.Data(),(GetCentralityStringA(cutSel)).Data());
   canvasDirectPhotonNLORaa->Print(Form("%s/DirectPhotonCombined_NLORaa_%s.eps",outputDir.Data(),cent.Data()));


   /// Subtract Combined NLO

   TCanvas *canvasDirectPhotonCombinedNLOSubtracted = GetAndSetCanvas("canvasDirectPhotonCombinedNLOSubtracted",0.15,0.1, 1000 ,1250);
   canvasDirectPhotonCombinedNLOSubtracted->SetLogy();
   dummy->GetYaxis()->SetRangeUser(2e-6,2e1);
   dummy->GetXaxis()->SetRangeUser(0.80,4.7);
   dummy->GetYaxis()->SetTitleOffset(1.5);
   dummy->GetXaxis()->SetLabelOffset(0);
   dummy->GetYaxis()->SetTitle("#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c}^{2})");
   dummy->DrawCopy();


   
   TGraphAsymmErrors *graphUpperLimitForPointAboveButErrorCombinedSubNLO;
   TGraphAsymmErrors *graphArrorErrorBelow0CombinedSubNLO;
   TGraphAsymmErrors *grSystFullCombinedSubNLO;

   if(graphUpperLimitForPointAboveButErrorCombinedSubNLO){
      graphUpperLimitForPointAboveButErrorCombinedSubNLO = (TGraphAsymmErrors*) graphUpperLimitForPointAboveButErrorCombined->Clone();
      graphArrorErrorBelow0CombinedSubNLO = (TGraphAsymmErrors*)graphArrorErrorBelow0Combined->Clone();
      grSystFullCombinedSubNLO = (TGraphAsymmErrors*)grSystFullCombined->Clone("grSystFullCombinedSubNLO");
      for(Int_t i = 0;i<grSystFullCombined->GetN();i++){
         Double_t x1,x2,y1,y2;
         grSystFullCombined->GetPoint(i,x1,y1);
         Double_t xerrlow = grSystFullCombined->GetErrorXlow(i);
         Double_t xerrhigh = grSystFullCombined->GetErrorXhigh(i);
         Double_t yerrlow = grSystFullCombined->GetErrorYlow(i);
         Double_t yerrhigh = grSystFullCombined->GetErrorYhigh(i);
         //cout<<x1-xerrlow<<"  "<<x1+xerrhigh<<endl;
         for(Int_t j = 0;j<PbPb276CTEQ61EPS09BFG2_sum_pdferr_InvYield->GetN();j++){
            PbPb276CTEQ61EPS09BFG2_sum_pdferr_InvYield->GetPoint(j,x2,y2);
            if(x2>x1-xerrlow && x2<x1+xerrhigh){
               Double_t yerrlowNLO = PbPb276CTEQ61EPS09BFG2_sum_pdferr_InvYield->GetErrorYlow(j);
               Double_t yerrhighNLO = PbPb276CTEQ61EPS09BFG2_sum_pdferr_InvYield->GetErrorYhigh(j);
               grSystFullCombinedSubNLO->SetPoint(i,x1,y1-y2);
               grSystFullCombinedSubNLO->SetPointError(i,xerrlow,xerrhigh,sqrt(yerrlow*yerrlow+yerrlowNLO*yerrlowNLO),sqrt(yerrhigh*yerrhigh+yerrhighNLO*yerrhighNLO));
               if((y1-y2-sqrt(yerrlow*yerrlow+yerrlowNLO*yerrlowNLO))<0){
                  graphUpperLimitForPointAboveButErrorCombinedSubNLO->SetPoint(graphUpperLimitForPointAboveButErrorCombinedSubNLO->GetN(),x1,y1-y2);
                  graphUpperLimitForPointAboveButErrorCombinedSubNLO->SetPointError(graphUpperLimitForPointAboveButErrorCombinedSubNLO->GetN()-1,0,0,0,1.28*abs(sqrt(yerrhigh*yerrhigh+yerrhighNLO*yerrhighNLO)));
                  graphArrorErrorBelow0CombinedSubNLO->SetPoint(graphArrorErrorBelow0CombinedSubNLO->GetN(),x1,1.28*abs(sqrt(yerrhigh*yerrhigh+yerrhighNLO*yerrhighNLO))+(y1-y2));
                  graphArrorErrorBelow0CombinedSubNLO->SetPointError(graphArrorErrorBelow0CombinedSubNLO->GetN()-1,0,0,1.28*abs(sqrt(yerrhigh*yerrhigh+yerrhighNLO*yerrhighNLO))+((y1-y2)*0.5*1.3),0);
               }
            }
         }
      }
   }

   TGraphAsymmErrors *grSystFullCombinedSubNLOA = (TGraphAsymmErrors*)grSystFullCombinedA->Clone("grSystFullCombinedSubNLO");
   for(Int_t i = 0;i<grSystFullCombinedA->GetN();i++){
      Double_t x1,x2,y1,y2;
      grSystFullCombinedA->GetPoint(i,x1,y1);
      Double_t xerrlow = grSystFullCombinedA->GetErrorXlow(i);
      Double_t xerrhigh = grSystFullCombinedA->GetErrorXhigh(i);
      Double_t yerrlow = grSystFullCombinedA->GetErrorYlow(i);
      Double_t yerrhigh = grSystFullCombinedA->GetErrorYhigh(i);
      cout<<x1-xerrlow<<"  "<<x1+xerrhigh<<endl;
      for(Int_t j = 0;j<PbPb276CTEQ61EPS09BFG2_sum_pdferr_InvYield->GetN();j++){
         PbPb276CTEQ61EPS09BFG2_sum_pdferr_InvYield->GetPoint(j,x2,y2);
         if(x2>x1-xerrlow && x2<x1+xerrhigh){
            Double_t yerrlowNLO = PbPb276CTEQ61EPS09BFG2_sum_pdferr_InvYield->GetErrorYlow(j);
            Double_t yerrhighNLO = PbPb276CTEQ61EPS09BFG2_sum_pdferr_InvYield->GetErrorYhigh(j);
            grSystFullCombinedSubNLOA->SetPoint(i,x1,y1-y2);
            grSystFullCombinedSubNLOA->SetPointError(i,xerrlow,xerrhigh,sqrt(yerrlow*yerrlow+yerrlowNLO*yerrlowNLO),sqrt(yerrhigh*yerrhigh+yerrhighNLO*yerrhighNLO));
         }
      }

   }


   TGraphAsymmErrors *grStatFullCombinedSubNLO = (TGraphAsymmErrors*)grStatFullCombined->Clone("grStatFullCombinedSubNLO");
   for(Int_t i = 0;i<grStatFullCombined->GetN();i++){
      Double_t x1,x2,y1,y2;
      grStatFullCombined->GetPoint(i,x1,y1);
      Double_t xerrlow = grSystFullCombined->GetErrorXlow(i);
      Double_t xerrhigh = grSystFullCombined->GetErrorXhigh(i);
      cout<<x1-xerrlow<<"  "<<x1+xerrhigh<<endl;
      for(Int_t j = 0;j<PbPb276CTEQ61EPS09BFG2_sum_pdferr_InvYield->GetN();j++){
         PbPb276CTEQ61EPS09BFG2_sum_pdferr_InvYield->GetPoint(j,x2,y2);
         if(x2>x1-xerrlow && x2<x1+xerrhigh){
            grStatFullCombinedSubNLO->SetPoint(i,x1,y1-y2);
         }
      }

   }

   TF1 *exponetialCombinedSystSubNLO = new TF1("Exponential Fit","[0]*exp(-x/[1])",fitRangeLow,fitRangeHigh);
   TF1 *exponetialCombinedStatSubNLO = new TF1("Exponential Fit","[0]*exp(-x/[1])",fitRangeLow,fitRangeHigh);

   exponetialCombinedSystSubNLO->SetParameters(100,0.300);
   exponetialCombinedStatSubNLO->SetParameters(100,0.300);

 for(Int_t i = 0;i<grSystFullCombinedSubNLOA->GetN();i++){
      grSystFullCombinedSubNLOA->SetPointEXhigh(i,0);
      grSystFullCombinedSubNLOA->SetPointEXlow(i,0);
   }
   for(Int_t i = 0;i<grStatFullCombinedSubNLO->GetN();i++){
      grStatFullCombinedSubNLO->SetPointEXhigh(i,0);
      grStatFullCombinedSubNLO->SetPointEXlow(i,0);
   }
   grSystFullCombinedSubNLOA->Fit(exponetialCombinedSystSubNLO,"NRE+","",fitRangeLow,fitRangeHigh);
   grStatFullCombinedSubNLO->Fit(exponetialCombinedStatSubNLO,"NRE+","",fitRangeLow,fitRangeHigh);

   exponetialCombinedStatSubNLO->SetLineColor(colorComb4060);
   exponetialCombinedStatSubNLO->Draw("same");
   
   if(grSystFullCombinedSubNLO){
      if(!centralityCutNumber.CompareTo("501"))grSystFullCombinedSubNLO->RemovePoint(11);
      if(!centralityCutNumber.CompareTo("501"))grSystFullCombinedSubNLO->RemovePoint(11);
      if(!centralityCutNumber.CompareTo("502"))grSystFullCombinedSubNLO->RemovePoint(11);
      if(!centralityCutNumber.CompareTo("502"))grSystFullCombinedSubNLO->RemovePoint(11);
      
      if(!centralityCutNumber.CompareTo("524"))grSystFullCombinedSubNLO->RemovePoint(1);
      if(!centralityCutNumber.CompareTo("524"))grSystFullCombinedSubNLO->RemovePoint(1);
      if(!centralityCutNumber.CompareTo("524"))grSystFullCombinedSubNLO->RemovePoint(9);

      grSystFullCombinedSubNLO->Draw("Z2");
   }

   grStatFullCombinedSubNLO->Draw("PE1");
   if(graphArrorErrorBelow0CombinedSubNLO)graphArrorErrorBelow0CombinedSubNLO->Draw(">same");
   if(graphUpperLimitForPointAboveButErrorCombinedSubNLO)graphUpperLimitForPointAboveButErrorCombinedSubNLO->Draw("||same");
   graphUpperLimitPointsAndErrorBelow0->Draw("||same");

   // if(graphUpperLimitForPointAboveButError0020CombinedSubNLO)graphUpperLimitForPointAboveButError0020CombinedSubNLO->Draw("||same");
   // if(graphUpperLimitPointsAndErrorBelow00020Combined)graphUpperLimitPointsAndErrorBelow00020Combined->Draw("||same");
   // if(graphArrorErrorBelow00020CombinedSubNLO)graphArrorErrorBelow00020CombinedSubNLO->Draw(">same");



   TLegend* legendDirectPhotonNLOSubtractedCombined = GetAndSetLegend(0.17,0.13,6.5);
   if(grSystFullCombinedSubNLO)legendDirectPhotonNLOSubtractedCombined->AddEntry(grSystFullCombinedSubNLO,Form("#gamma_{dir} ALICE"),"pf");
   legendDirectPhotonNLOSubtractedCombined->AddEntry((TObject*)0,"NLO subtracted","");
   legendDirectPhotonNLOSubtractedCombined->AddEntry((TObject*)0,"(nPDF: EPS09, FF: BFG2)","");
   if(!(!centralityCutNumber.CompareTo("548"))){
      legendDirectPhotonNLOSubtractedCombined->AddEntry(exponetialSubNLOStat,Form("exp(-#it{p}_{T}/T) (%.1f < #it{p}_{T} < %.1f GeV/#it{c})",fitRangeLow,fitRangeHigh),"l");
      legendDirectPhotonNLOSubtractedCombined->AddEntry((TObject*)0,Form("T = %.0f #pm %.0f^{stat} #pm %0.0f^{syst} MeV",exponetialCombinedStatSubNLO->GetParameter(1)*1000,exponetialCombinedStatSubNLO->GetParError(1)*1000,exponetialCombinedSystSubNLO->GetParError(1)*1000),"");
   }
   legendDirectPhotonNLOSubtractedCombined->Draw();

   if(grSystFullCombinedSubNLO)grSystFullCombinedSubNLO->Draw("Z2");
   grStatFullCombinedSubNLO->Draw("PZe1");
   //grSystFullCombinedSubNLOA->Draw("PZe1");

   DrawSystem(0.25,0.9,option.Data(),(GetCentralityString(cutSel)).Data());
   canvasDirectPhotonCombinedNLOSubtracted->Print(Form("%s/DirectPhotonSpectrumSubtractedNLOCombined_%s.eps",outputDir.Data(),cent.Data()));


   TFile *CombinedDirectPhotons = new TFile("CombinedDiretPhotons.root","UPDATE");
   if(grSystFullCombined)grSystFullCombined->Write(Form("DP_syst_%s",centralityCutNumber.Data()),TObject::kOverwrite);
   if(grStatFullCombined)grStatFullCombined->Write(Form("DP_stat_%s",centralityCutNumber.Data()),TObject::kOverwrite);
   if(graphUpperLimitForPointAboveButErrorCombined)graphUpperLimitForPointAboveButErrorCombined->Write(Form("DP_UL1_%s",centralityCutNumber.Data()),TObject::kOverwrite);
   if(graphUpperLimitForPointAboveButErrorCombinedSubNLO)graphUpperLimitForPointAboveButErrorCombinedSubNLO->Write(Form("DP_UL1_subNLO_%s",centralityCutNumber.Data()),TObject::kOverwrite);
   if(graphUpperLimitPointsAndErrorBelow0Combined)graphUpperLimitPointsAndErrorBelow0Combined->Write(Form("DP_UL2_%s",centralityCutNumber.Data()),TObject::kOverwrite);
   if(graphArrorErrorBelow0Combined)graphArrorErrorBelow0Combined->Write(Form("DP_ARR1_%s",centralityCutNumber.Data()),TObject::kOverwrite);
   if(graphArrorErrorBelow0CombinedSubNLO)graphArrorErrorBelow0CombinedSubNLO->Write(Form("DP_ARR1_subNLO_%s",centralityCutNumber.Data()),TObject::kOverwrite);
   if(graphArrowPointErrorBelow0Combined)graphArrowPointErrorBelow0Combined->Write(Form("DP_ARR2_%s",centralityCutNumber.Data()),TObject::kOverwrite);
   PbPb276CTEQ61EPS09BFG2_sum_pdferr_InvYield->Write(Form("NLO_PbPb_%s",centralityCutNumber.Data()),TObject::kOverwrite);
   pp276CT10BFG2_sum_pdferr_InvYield->Write(Form("NLO_pp1_%s",centralityCutNumber.Data()),TObject::kOverwrite);
   NLO->Write(Form("NLO_pp2_%s",centralityCutNumber.Data()),TObject::kOverwrite);
   exponetialCombinedSyst->Write(Form("Exp_Syst_%s",centralityCutNumber.Data()),TObject::kOverwrite);
   exponetialCombinedStat->Write(Form("Exp_Stat_%s",centralityCutNumber.Data()),TObject::kOverwrite);

   DoubleRatioCombinedSyst->Write(Form("DR_Comb_syst_%s",centralityCutNumber.Data()),TObject::kOverwrite);
   DoubleRatioCombinedStat->Write(Form("DR_Comb_stat_%s",centralityCutNumber.Data()),TObject::kOverwrite);
   DoubleRatioPi0FitRebined->Write(Form("DR_PCM_stat_%s",centralityCutNumber.Data()),TObject::kOverwrite);
   DoubleRatioPi0FitErrorsAsym->Write(Form("DR_PCM_syst_%s",centralityCutNumber.Data()),TObject::kOverwrite);
   DoubleRatioPHOSStat->Write(Form("DR_PHOS_stat_%s",centralityCutNumber.Data()),TObject::kOverwrite);
   GraphDoubleRatioPHOSSystRebinned->Write(Form("DR_PHOS_syst_%s",centralityCutNumber.Data()),TObject::kOverwrite);
   if(grSystFullCombinedSubNLO)grSystFullCombinedSubNLO->Write(Form("DP_sub_syst_%s",centralityCutNumber.Data()),TObject::kOverwrite);
   grStatFullCombinedSubNLO->Write(Form("DP_sub_stat_%s",centralityCutNumber.Data()),TObject::kOverwrite);
   exponetialCombinedSystSubNLO->Write(Form("Exp_sub_Syst_%s",centralityCutNumber.Data()),TObject::kOverwrite);
   exponetialCombinedStatSubNLO->Write(Form("Exp_sub_Stat_%s",centralityCutNumber.Data()),TObject::kOverwrite);

   CombinedDirectPhotons->Close();


   Double_t arrayBoundsXDPCombined[7];
   Double_t arrayBoundsYDPCombined[7];
   Double_t relativeMarginsDPCombinedRatioX[7];
   Double_t relativeMarginsDPCombinedRatioY[7];
   ReturnCorrectValuesForCanvasScaling(1200,1300, 1, 1, 0.12, 0.005, 0.005,0.07,arrayBoundsXDPCombined,arrayBoundsYDPCombined,relativeMarginsDPCombinedRatioX,relativeMarginsDPCombinedRatioY);

   TCanvas *canvasDirectPhotonCombinedAllCent = GetAndSetCanvas("canvasDirectPhotonCombinedAllCent",0.12,0.07, 1200 ,1300);
   canvasDirectPhotonCombinedAllCent->SetLogy();




   Int_t textSizeLabelsPixel = 40;
   Double_t margin = 0.4*1200;
   Double_t textsizeLabelsInter = 0;
   Double_t textsizeFacInter = 0;
   Double_t textsizeLabelsTop = 0;
   Double_t textsizeFacTop = 0;
   Double_t textsizeLabelsBottom = 0;
   Double_t textsizeFacBottom = 0;
   Double_t markerSize = 2;

   if (gPad->XtoPixel(gPad->GetX2()) < gPad->YtoPixel(gPad->GetY1())){
      textsizeLabelsTop = (Double_t)textSizeLabelsPixel/gPad->XtoPixel(gPad->GetX2()) ;
      textsizeFacTop = (Double_t)1./gPad->XtoPixel(gPad->GetX2()) ;
   } else {
      textsizeLabelsTop = (Double_t)textSizeLabelsPixel/gPad->YtoPixel(gPad->GetY1());
      textsizeFacTop = (Double_t)1./gPad->YtoPixel(gPad->GetY1());
   }

   SetStyleHistoTH1ForGraphs(dummy, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c}^{2})", 0.85*textsizeLabelsTop,textsizeLabelsTop, 0.85*textsizeLabelsTop,textsizeLabelsTop, 1.,0.6/(textsizeFacTop*margin), 515, 409);
   dummy->GetXaxis()->SetRangeUser(0.50,11.4);
   dummy->GetYaxis()->SetRangeUser(1e-11,4e3);
   dummy->DrawCopy();


   TFile *CombinedDirectPhotonsPlot = new TFile("CombinedDiretPhotons.root");

   TGraphAsymmErrors *grSystFull0020Combined = (TGraphAsymmErrors*)CombinedDirectPhotonsPlot->Get("DP_syst_502");
   TGraphAsymmErrors *grSystFull2040Combined = (TGraphAsymmErrors*)CombinedDirectPhotonsPlot->Get("DP_syst_524");
   grSystFull2040Combined = ScaleGraph(grSystFull2040Combined,0.1);
   TGraphAsymmErrors *grSystFull4080Combined = (TGraphAsymmErrors*)CombinedDirectPhotonsPlot->Get("DP_syst_548");
   grSystFull4080Combined= ScaleGraph(grSystFull4080Combined,0.01);
   TGraphAsymmErrors *grStatFull0020Combined = (TGraphAsymmErrors*)CombinedDirectPhotonsPlot->Get("DP_stat_502");
   TGraphAsymmErrors *grStatFull2040Combined = (TGraphAsymmErrors*)CombinedDirectPhotonsPlot->Get("DP_stat_524");
   grStatFull2040Combined= ScaleGraph(grStatFull2040Combined,0.1);
   TGraphAsymmErrors *grStatFull4080Combined = (TGraphAsymmErrors*)CombinedDirectPhotonsPlot->Get("DP_stat_548");
   grStatFull4080Combined= ScaleGraph(grStatFull4080Combined,0.01);
   TGraphAsymmErrors *graphUpperLimitForPointAboveButError0020Combined = (TGraphAsymmErrors*)CombinedDirectPhotonsPlot->Get("DP_UL1_502");
   TGraphAsymmErrors *graphUpperLimitForPointAboveButError0020CombinedSubNLO = (TGraphAsymmErrors*)CombinedDirectPhotonsPlot->Get("DP_UL1_subNLO_502");
   TGraphAsymmErrors *graphUpperLimitForPointAboveButError2040Combined = (TGraphAsymmErrors*)CombinedDirectPhotonsPlot->Get("DP_UL1_524");
   TGraphAsymmErrors *graphUpperLimitForPointAboveButError2040CombinedSubNLO = (TGraphAsymmErrors*)CombinedDirectPhotonsPlot->Get("DP_UL1_subNLO_524");
   graphUpperLimitForPointAboveButError2040Combined= ScaleGraph(graphUpperLimitForPointAboveButError2040Combined,0.1);
   TGraphAsymmErrors *graphUpperLimitForPointAboveButError4080Combined = (TGraphAsymmErrors*)CombinedDirectPhotonsPlot->Get("DP_UL1_548");
   graphUpperLimitForPointAboveButError4080Combined= ScaleGraph(graphUpperLimitForPointAboveButError4080Combined,0.01);
   TGraphAsymmErrors *graphUpperLimitPointsAndErrorBelow00020Combined = (TGraphAsymmErrors*)CombinedDirectPhotonsPlot->Get("DP_UL2_502");
   TGraphAsymmErrors *graphUpperLimitPointsAndErrorBelow02040Combined = (TGraphAsymmErrors*)CombinedDirectPhotonsPlot->Get("DP_UL2_524");
   graphUpperLimitPointsAndErrorBelow02040Combined = ScaleGraph(graphUpperLimitPointsAndErrorBelow02040Combined,0.1);
   TGraphAsymmErrors *graphUpperLimitPointsAndErrorBelow04080Combined = (TGraphAsymmErrors*)CombinedDirectPhotonsPlot->Get("DP_UL2_548");
   graphUpperLimitPointsAndErrorBelow04080Combined= ScaleGraph(graphUpperLimitPointsAndErrorBelow04080Combined,0.01);
   TGraphAsymmErrors *graphArrorErrorBelow00020Combined = (TGraphAsymmErrors*)CombinedDirectPhotonsPlot->Get("DP_ARR1_502");
   TGraphAsymmErrors *graphArrorErrorBelow00020CombinedSubNLO = (TGraphAsymmErrors*)CombinedDirectPhotonsPlot->Get("DP_ARR1_subNLO_502");
   TGraphAsymmErrors *graphArrorErrorBelow02040Combined = (TGraphAsymmErrors*)CombinedDirectPhotonsPlot->Get("DP_ARR1_524");
   TGraphAsymmErrors *graphArrorErrorBelow02040CombinedSubNLO = (TGraphAsymmErrors*)CombinedDirectPhotonsPlot->Get("DP_ARR1_subNLO_524");
   graphArrorErrorBelow02040Combined= ScaleGraph(graphArrorErrorBelow02040Combined,0.1);
   TGraphAsymmErrors *graphArrorErrorBelow04080Combined = (TGraphAsymmErrors*)CombinedDirectPhotonsPlot->Get("DP_ARR1_548");
   graphArrorErrorBelow04080Combined= ScaleGraph(graphArrorErrorBelow04080Combined,0.01);
   TGraphAsymmErrors *graphArrowPointErrorBelow00020Combined = (TGraphAsymmErrors*)CombinedDirectPhotonsPlot->Get("DP_ARR2_502");
   TGraphAsymmErrors *graphArrowPointErrorBelow02040Combined = (TGraphAsymmErrors*)CombinedDirectPhotonsPlot->Get("DP_ARR2_524");
   graphArrowPointErrorBelow02040Combined= ScaleGraph(graphArrowPointErrorBelow02040Combined,0.1);
   TGraphAsymmErrors *graphArrowPointErrorBelow04080Combined = (TGraphAsymmErrors*)CombinedDirectPhotonsPlot->Get("DP_ARR2_548");
   graphArrowPointErrorBelow04080Combined= ScaleGraph(graphArrowPointErrorBelow04080Combined,0.01);
   TGraphAsymmErrors *NLO0020Combined = (TGraphAsymmErrors*)CombinedDirectPhotonsPlot->Get("NLO_pp1_502");
   TGraphAsymmErrors *NLO2040Combined = (TGraphAsymmErrors*)CombinedDirectPhotonsPlot->Get("NLO_pp1_524");
   NLO2040Combined= ScaleGraph(NLO2040Combined,0.1);
   TGraphAsymmErrors *NLO4080Combined = (TGraphAsymmErrors*)CombinedDirectPhotonsPlot->Get("NLO_pp1_548");
   NLO4080Combined= ScaleGraph(NLO4080Combined,0.01);
   TGraphAsymmErrors *NLOpp0020Combined = (TGraphAsymmErrors*)CombinedDirectPhotonsPlot->Get("NLO_pp2_502");
   TGraphAsymmErrors *NLOpp2040Combined = (TGraphAsymmErrors*)CombinedDirectPhotonsPlot->Get("NLO_pp2_524");
   NLOpp2040Combined= ScaleGraph(NLOpp2040Combined,0.1);
   TGraphAsymmErrors *NLOpp4080Combined = (TGraphAsymmErrors*)CombinedDirectPhotonsPlot->Get("NLO_pp2_548");
   NLOpp4080Combined= ScaleGraph(NLOpp4080Combined,0.01);
   TGraphAsymmErrors *NLOPb0020Combined = (TGraphAsymmErrors*)CombinedDirectPhotonsPlot->Get("NLO_PbPb_502");
   TGraphAsymmErrors *NLOPb2040Combined = (TGraphAsymmErrors*)CombinedDirectPhotonsPlot->Get("NLO_PbPb_524");
   NLOPb2040Combined= ScaleGraph(NLOPb2040Combined,0.1);
   TGraphAsymmErrors *NLOPb4080Combined = (TGraphAsymmErrors*)CombinedDirectPhotonsPlot->Get("NLO_PbPb_548");
   NLOPb4080Combined= ScaleGraph(NLOPb4080Combined,0.01);

   TF1 *ExpStat0020 = (TF1*) CombinedDirectPhotonsPlot->Get("Exp_Stat_502");
   TF1 *ExpStat2040 = (TF1*) CombinedDirectPhotonsPlot->Get("Exp_Stat_524");
   if(ExpStat2040)ExpStat2040->SetParameter(0,ExpStat2040->GetParameter(0)*0.1);
   TF1 *ExpStat4080 = (TF1*) CombinedDirectPhotonsPlot->Get("Exp_Stat_548");
   if(ExpStat4080)ExpStat4080->SetParameter(0,ExpStat4080->GetParameter(0)*0.01);
   TF1 *ExpSyst0020 = (TF1*) CombinedDirectPhotonsPlot->Get("Exp_Syst_502");
   TF1 *ExpSyst2040 = (TF1*) CombinedDirectPhotonsPlot->Get("Exp_Syst_524");
   if(ExpSyst2040)ExpSyst2040->SetParameter(0,ExpSyst2040->GetParameter(0)*0.1);
   TF1 *ExpSyst4080 = (TF1*) CombinedDirectPhotonsPlot->Get("Exp_Syst_548");
   if(ExpSyst4080)ExpSyst4080->SetParameter(0,ExpSyst4080->GetParameter(0)*0.01);

   TLegend* legendDirectPhotonCombinedAllCent = GetAndSetLegend(0.24,0.64,5.0,2,"Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV");
   legendDirectPhotonCombinedAllCent->SetTextSize(textsizeLabelsTop);
   legendDirectPhotonCombinedAllCent->AddEntry(grSystFull0020Combined,Form("  0-20%"),"pf");
   legendDirectPhotonCombinedAllCent->AddEntry(NLO,"PDF: CTEQ6M5 FF: DSS","l");
   legendDirectPhotonCombinedAllCent->AddEntry(grSystFull2040Combined,Form("20-40%"),"pf");
   legendDirectPhotonCombinedAllCent->AddEntry((TObject*)0,"JETPHOX:","");
   legendDirectPhotonCombinedAllCent->AddEntry(grStatFull4080Combined,Form("40-80%"),"pf");
   legendDirectPhotonCombinedAllCent->AddEntry(pp276CT10BFG2_sum_pdferr_InvYield,"PDF: CT10, FF: BFG2","f");
   legendDirectPhotonCombinedAllCent->AddEntry((TObject*)0,"","");
   legendDirectPhotonCombinedAllCent->AddEntry(PbPb276CTEQ61EPS09BFG2_sum_pdferr_InvYield,"nPDF: EPS09, FF: BFG2","f");
   legendDirectPhotonCombinedAllCent->AddEntry((TObject*)0,"","");
   legendDirectPhotonCombinedAllCent->AddEntry((TObject*)0,"(all N_{coll} scaled pp)","");
   legendDirectPhotonCombinedAllCent->Draw();

   TLegend* legendDirectPhotonCombinedAllCentNLO = GetAndSetLegend(0.15,0.10,1.7,1,Form("exp(-#it{p}_{T}/T) (fit range %.1f < #it{p}_{T} < %.1f GeV/#it{c})",fitRangeLow,fitRangeHigh));
   legendDirectPhotonCombinedAllCentNLO->SetTextSize(textsizeLabelsTop);
   //legendDirectPhotonCombinedAllCentNLO->AddEntry((TObject*)0,Form("exp(-#it{p}_{T}/T) (fit range %.1f < #it{p}_{T} < %.1f GeV/#it{c})",fitRangeLow,fitRangeHigh),"");
   legendDirectPhotonCombinedAllCentNLO->AddEntry(ExpStat0020,Form("T = %.0f #pm %.0f^{stat} #pm %0.0f^{syst} MeV (0-20%)",ExpStat0020->GetParameter(1)*1000,ExpStat0020->GetParError(1)*1000,ExpSyst0020->GetParError(1)*1000),"l");
   legendDirectPhotonCombinedAllCentNLO->AddEntry(ExpStat2040,Form("T = %.0f #pm %.0f^{stat} #pm %0.0f^{syst} MeV (20-40%)",ExpStat2040->GetParameter(1)*1000,ExpStat2040->GetParError(1)*1000,ExpSyst2040->GetParError(1)*1000),"l");
   legendDirectPhotonCombinedAllCentNLO->Draw();

   //DrawSystem(0.14,0.92,option.Data(),"",textsizeLabelsTop);




   if(NLOPb0020Combined)NLOPb0020Combined->Draw("E3same");
   if(NLOPb2040Combined)NLOPb2040Combined->Draw("E3same");
   if(NLOPb4080Combined)NLOPb4080Combined->Draw("E3same");

   if(NLO0020Combined)NLO0020Combined->Draw("E3same");
   if(NLO2040Combined)NLO2040Combined->Draw("E3same");
   if(NLO4080Combined)NLO4080Combined->Draw("E3same");

   if(NLOpp0020Combined)NLOpp0020Combined->Draw("E3same");
   if(NLOpp2040Combined)NLOpp2040Combined->Draw("E3same");
   if(NLOpp4080Combined)NLOpp4080Combined->Draw("E3same");

   if(grSystFull0020Combined)grSystFull0020Combined->Draw("Z2same");
   if(grSystFull2040Combined)grSystFull2040Combined->Draw("Z2same");
   if(grSystFull4080Combined)grSystFull4080Combined->Draw("Z2same");

   if(grStatFull0020Combined)grStatFull0020Combined->Draw("PZe1same");
   if(grStatFull2040Combined)grStatFull2040Combined->Draw("PZe1same");
   if(grStatFull4080Combined)grStatFull4080Combined->Draw("PZe1same");

   if(graphUpperLimitForPointAboveButError0020Combined)graphUpperLimitForPointAboveButError0020Combined->Draw("||same");
   if(graphUpperLimitForPointAboveButError2040Combined)graphUpperLimitForPointAboveButError2040Combined->Draw("||same");
   if(graphUpperLimitForPointAboveButError4080Combined)graphUpperLimitForPointAboveButError4080Combined->Draw("||same");

   if(graphUpperLimitPointsAndErrorBelow00020Combined)graphUpperLimitPointsAndErrorBelow00020Combined->Draw("||same");
   if(graphUpperLimitPointsAndErrorBelow02040Combined)graphUpperLimitPointsAndErrorBelow02040Combined->Draw("||same");
   if(graphUpperLimitPointsAndErrorBelow04080Combined)graphUpperLimitPointsAndErrorBelow04080Combined->Draw("||same");

   if(graphArrorErrorBelow00020Combined)graphArrorErrorBelow00020Combined->Draw(">same");
   if(graphArrorErrorBelow02040Combined)graphArrorErrorBelow02040Combined->Draw(">same");
   if(graphArrorErrorBelow04080Combined)graphArrorErrorBelow04080Combined->Draw(">same");

   if(graphArrowPointErrorBelow00020Combined)graphArrowPointErrorBelow00020Combined->Draw("Zsame");
   if(graphArrowPointErrorBelow02040Combined)graphArrowPointErrorBelow02040Combined->Draw("Zsame");
   if(graphArrowPointErrorBelow04080Combined)graphArrowPointErrorBelow04080Combined->Draw("Zsame");

   if(ExpStat0020){
      ExpStat0020->SetLineColor(kCyan-8);
      ExpStat0020->Draw("same");
   }
   if(ExpStat2040){
      ExpStat2040->SetLineColor(kYellow-8);
      ExpStat2040->Draw("same");
   }


   TLatex *Scaled2040 = new TLatex(0.88,0.33,"#times 0.1");
   Scaled2040->SetNDC();
   Scaled2040->SetTextColor(colorComb2040);
   Scaled2040->SetTextFont(42);
   Scaled2040->SetTextSize(0.03);
   Scaled2040->SetLineWidth(2);
   Scaled2040->Draw();

   TLatex *Scaled4080 = new TLatex(0.88,0.22,"#times 0.01");
   Scaled4080->SetNDC();
   Scaled4080->SetTextColor(colorComb4060);
   Scaled4080->SetTextFont(42);
   Scaled4080->SetTextSize(0.03);
   Scaled4080->SetLineWidth(2);
   Scaled4080->Draw();

   gPad->RedrawAxis();
   canvasDirectPhotonCombinedAllCent->Print(Form("%s/DirectPhotonSpectrumCombined_allCent.eps",outputDir.Data()));

  
   TH1D *DR_PCM_stat_502 = (TH1D*)CombinedDirectPhotonsPlot->Get("DR_PCM_stat_502");
   TH1D *DR_PCM_stat_524 = (TH1D*)CombinedDirectPhotonsPlot->Get("DR_PCM_stat_524");
   TH1D *DR_PCM_stat_548 = (TH1D*)CombinedDirectPhotonsPlot->Get("DR_PCM_stat_548");

   TGraphAsymmErrors *DR_PCM_syst_502 = (TGraphAsymmErrors*)CombinedDirectPhotonsPlot->Get("DR_PCM_syst_502");
   TGraphAsymmErrors *DR_PCM_syst_524 = (TGraphAsymmErrors*)CombinedDirectPhotonsPlot->Get("DR_PCM_syst_524");
   TGraphAsymmErrors *DR_PCM_syst_548 = (TGraphAsymmErrors*)CombinedDirectPhotonsPlot->Get("DR_PCM_syst_548");

   TH1D *DR_PHOS_stat_502 = (TH1D*)CombinedDirectPhotonsPlot->Get("DR_PHOS_stat_502");
   TH1D *DR_PHOS_stat_524 = (TH1D*)CombinedDirectPhotonsPlot->Get("DR_PHOS_stat_524");
   TH1D *DR_PHOS_stat_548 = (TH1D*)CombinedDirectPhotonsPlot->Get("DR_PHOS_stat_548");

   TGraphAsymmErrors *DR_PHOS_syst_502 = (TGraphAsymmErrors*)CombinedDirectPhotonsPlot->Get("DR_PHOS_syst_502");
   TGraphAsymmErrors *DR_PHOS_syst_524 = (TGraphAsymmErrors*)CombinedDirectPhotonsPlot->Get("DR_PHOS_syst_524");
   TGraphAsymmErrors *DR_PHOS_syst_548 = (TGraphAsymmErrors*)CombinedDirectPhotonsPlot->Get("DR_PHOS_syst_548");

   TGraphAsymmErrors *DR_Comb_syst_502 = (TGraphAsymmErrors*)CombinedDirectPhotonsPlot->Get("DR_Comb_syst_502");
   TGraphAsymmErrors *DR_Comb_syst_524 = (TGraphAsymmErrors*)CombinedDirectPhotonsPlot->Get("DR_Comb_syst_524");
   TGraphAsymmErrors *DR_Comb_syst_548 = (TGraphAsymmErrors*)CombinedDirectPhotonsPlot->Get("DR_Comb_syst_548");

   TGraphAsymmErrors *DR_Comb_stat_502 = (TGraphAsymmErrors*)CombinedDirectPhotonsPlot->Get("DR_Comb_stat_502");
   TGraphAsymmErrors *DR_Comb_stat_524 = (TGraphAsymmErrors*)CombinedDirectPhotonsPlot->Get("DR_Comb_stat_524");
   TGraphAsymmErrors *DR_Comb_stat_548 = (TGraphAsymmErrors*)CombinedDirectPhotonsPlot->Get("DR_Comb_stat_548");



   Double_t arrayBoundsXDRCombined[7];
   Double_t arrayBoundsYDRCombined[7];
   Double_t relativeMarginsDRCombinedRatioX[7];
   Double_t relativeMarginsDRCombinedRatioY[7];
   ReturnCorrectValuesForCanvasScaling(1200,1500, 1, 3, 0.1, 0.005, 0.005,0.06,arrayBoundsXDRCombined,arrayBoundsYDRCombined,relativeMarginsDRCombinedRatioX,relativeMarginsDRCombinedRatioY);
   TCanvas *canvasDoubleRatioCombinedAllCent = GetAndSetCanvas("canvasDoubleRatioCombinedAllCent",0.15,0.1, 1200 ,1500);
   TPad* pad0020combined = new TPad("pad0020", "",arrayBoundsXDRCombined[0], arrayBoundsYDRCombined[1],arrayBoundsXDRCombined[1], arrayBoundsYDRCombined[0],-1, -1, -2);
   DrawGammaPadSettings( pad0020combined, relativeMarginsDRCombinedRatioX[0], relativeMarginsDRCombinedRatioX[2], relativeMarginsDRCombinedRatioY[0], relativeMarginsDRCombinedRatioY[1]);
   pad0020combined->Draw();
   TPad* pad2040combined = new TPad("pad2040", "",arrayBoundsXDRCombined[0], arrayBoundsYDRCombined[2],arrayBoundsXDRCombined[1], arrayBoundsYDRCombined[1],-1, -1, -2);
   DrawGammaPadSettings( pad2040combined, relativeMarginsDRCombinedRatioX[0], relativeMarginsDRCombinedRatioX[2], relativeMarginsDRCombinedRatioY[1], relativeMarginsDRCombinedRatioY[1]);
   pad2040combined->Draw();
   TPad* pad4080combined = new TPad("pad4080","",arrayBoundsXDRCombined[0], arrayBoundsYDRCombined[3],arrayBoundsXDRCombined[1], arrayBoundsYDRCombined[2],-1, -1, -2);
   DrawGammaPadSettings( pad4080combined,relativeMarginsDRCombinedRatioX[0], relativeMarginsDRCombinedRatioX[2], relativeMarginsDRCombinedRatioY[1], relativeMarginsDRCombinedRatioY[2]);
   pad4080combined->Draw();

   textSizeLabelsPixel = 40;
   margin = 0.4*1200;
   textsizeLabelsInter = 0;
   textsizeFacInter = 0;
   textsizeLabelsTop = 0;
   textsizeFacTop = 0;
   textsizeLabelsBottom = 0;
   textsizeFacBottom = 0;
   markerSize = 2;

   if (pad0020combined->XtoPixel(pad0020combined->GetX2()) < pad0020combined->YtoPixel(pad0020combined->GetY1())){
      textsizeLabelsTop = (Double_t)textSizeLabelsPixel/pad0020combined->XtoPixel(pad0020combined->GetX2()) ;
      textsizeFacTop = (Double_t)1./pad0020combined->XtoPixel(pad0020combined->GetX2()) ;
   } else {
      textsizeLabelsTop = (Double_t)textSizeLabelsPixel/pad0020combined->YtoPixel(pad0020combined->GetY1());
      textsizeFacTop = (Double_t)1./pad0020combined->YtoPixel(pad0020combined->GetY1());
   }
   if (pad2040combined->XtoPixel(pad2040combined->GetX2()) < pad2040combined->YtoPixel(pad2040combined->GetY1())){
      textsizeLabelsInter = (Double_t)textSizeLabelsPixel/pad2040combined->XtoPixel(pad2040combined->GetX2()) ;
      textsizeFacInter = (Double_t)1./pad2040combined->XtoPixel(pad2040combined->GetX2()) ;
   } else {
      textsizeLabelsInter = (Double_t)textSizeLabelsPixel/pad2040combined->YtoPixel(pad2040combined->GetY1());
      textsizeFacInter = (Double_t)1./pad2040combined->YtoPixel(pad2040combined->GetY1());
   }
   if (pad4080combined->XtoPixel(pad4080combined->GetX2()) < pad4080combined->YtoPixel(pad4080combined->GetY1())){
      textsizeLabelsBottom = (Double_t)textSizeLabelsPixel/pad4080combined->XtoPixel(pad4080combined->GetX2()) ;
      textsizeFacBottom = (Double_t)1./pad4080combined->XtoPixel(pad4080combined->GetX2()) ;
   } else {
      textsizeLabelsBottom = (Double_t)textSizeLabelsPixel/pad4080combined->YtoPixel(pad4080combined->GetY1());
      textsizeFacBottom = (Double_t)1./pad4080combined->YtoPixel(pad4080combined->GetY1());
   }


   dummy->GetYaxis()->SetRangeUser(0.75, 1.95);
   dummy->GetXaxis()->SetRangeUser(doubleRatioX[0],12);

   pad0020combined->cd();
   pad0020combined->SetLogx();
   SetStyleHistoTH1ForGraphs(dummy, "#it{p}_{T} (GeV/#it{c})","(#gamma_{inc}/#pi^{0})/(#gamma_{decay}/#pi^{0})", 0.85*textsizeLabelsTop,textsizeLabelsTop, 0.85*textsizeLabelsTop,textsizeLabelsTop, 1,0.5/(textsizeFacTop*margin), 515, 409);
   dummy->DrawCopy();
   ONE2->Draw("same");
   DR_Comb_syst_502->SetMarkerSize(markerSize);
   DR_Comb_syst_502->Draw("pZ2");
   DR_Comb_stat_502->Draw("pe1");
   DrawSystem(0.15,0.83,option.Data(),"0-20%",textsizeLabelsTop);

   pad2040combined->cd();
   pad2040combined->SetLogx();
   SetStyleHistoTH1ForGraphs(dummy, "#it{p}_{T} (GeV/#it{c})","(#gamma_{inc}/#pi^{0})/(#gamma_{decay}/#pi^{0})", 0.85*textsizeLabelsInter,textsizeLabelsInter, 0.85*textsizeLabelsInter,textsizeLabelsInter, 1,0.5/(textsizeFacInter*margin), 515, 409);
   dummy->DrawCopy();
   ONE2->Draw("same");
   DR_Comb_syst_524->SetMarkerSize(markerSize);
   DR_Comb_syst_524->Draw("pZ2");
   DR_Comb_stat_524->Draw("pe1");
   DrawSystem(0.15,0.85,option.Data(),"20-40%",textsizeLabelsInter);

   pad4080combined->cd();
   pad4080combined->SetLogx();
SetStyleHistoTH1ForGraphs(dummy, "#it{p}_{T} (GeV/#it{c})","(#gamma_{inc}/#pi^{0})/(#gamma_{decay}/#pi^{0})", 0.85*textsizeLabelsBottom,textsizeLabelsBottom, 0.85*textsizeLabelsBottom,textsizeLabelsBottom, 1,0.5/(textsizeFacBottom*margin), 515, 409);
 dummy->DrawCopy();
   ONE2->Draw("same");
   DR_Comb_syst_548->Draw("pZ2");
   DR_Comb_syst_548->SetMarkerSize(markerSize);
   DR_Comb_stat_548->Draw("pe1");
   DrawSystem(0.15,0.88,option.Data(),"40-80%",textsizeLabelsBottom);
   canvasDoubleRatioCombinedAllCent->Print(Form("%s/DoubleRatioCombined_allCent.eps",outputDir.Data()));



   Double_t arrayBoundsXDRSingle[7];
   Double_t arrayBoundsYDRSingle[7];
   Double_t relativeMarginsDRSingleRatioX[7];
   Double_t relativeMarginsDRSingleRatioY[7];
   ReturnCorrectValuesForCanvasScaling(1200,1500, 1, 3, 0.1, 0.005, 0.005,0.06,arrayBoundsXDRSingle,arrayBoundsYDRSingle,relativeMarginsDRSingleRatioX,relativeMarginsDRSingleRatioY);

   TCanvas *canvasDoubleRatioSingleAllCent = GetAndSetCanvas("canvasDoubleRatioSingleAllCent",0.1,0.06, 1200 ,1500);
   TPad* pad0020single = new TPad("pad0020", "",arrayBoundsXDRSingle[0], arrayBoundsYDRSingle[1],arrayBoundsXDRSingle[1], arrayBoundsYDRSingle[0],-1, -1, -2);
   DrawGammaPadSettings( pad0020single, relativeMarginsDRSingleRatioX[0], relativeMarginsDRSingleRatioX[2], relativeMarginsDRSingleRatioY[0], relativeMarginsDRSingleRatioY[1]);
   pad0020single->Draw();
   TPad* pad2040single = new TPad("pad2040", "",arrayBoundsXDRSingle[0], arrayBoundsYDRSingle[2],arrayBoundsXDRSingle[1], arrayBoundsYDRSingle[1],-1, -1, -2);
   DrawGammaPadSettings( pad2040single, relativeMarginsDRSingleRatioX[0], relativeMarginsDRSingleRatioX[2], relativeMarginsDRSingleRatioY[1], relativeMarginsDRSingleRatioY[1]);
   pad2040single->Draw();
   TPad* pad4080single = new TPad("pad4080","",arrayBoundsXDRSingle[0], arrayBoundsYDRSingle[3],arrayBoundsXDRSingle[1], arrayBoundsYDRSingle[2],-1, -1, -2);
   DrawGammaPadSettings( pad4080single,relativeMarginsDRSingleRatioX[0], relativeMarginsDRSingleRatioX[2], relativeMarginsDRSingleRatioY[1], relativeMarginsDRSingleRatioY[2]);
   pad4080single->Draw();

   textSizeLabelsPixel = 40;
   margin = 0.4*1200;
   textsizeLabelsInter = 0;
   textsizeFacInter = 0;
   textsizeLabelsTop = 0;
   textsizeFacTop = 0;
   textsizeLabelsBottom = 0;
   textsizeFacBottom = 0;

   if (pad0020single->XtoPixel(pad0020single->GetX2()) < pad0020single->YtoPixel(pad0020single->GetY1())){
      textsizeLabelsTop = (Double_t)textSizeLabelsPixel/pad0020single->XtoPixel(pad0020single->GetX2()) ;
      textsizeFacTop = (Double_t)1./pad0020single->XtoPixel(pad0020single->GetX2()) ;
   } else {
      textsizeLabelsTop = (Double_t)textSizeLabelsPixel/pad0020single->YtoPixel(pad0020single->GetY1());
      textsizeFacTop = (Double_t)1./pad0020single->YtoPixel(pad0020single->GetY1());
   }
   if (pad2040single->XtoPixel(pad2040single->GetX2()) < pad2040single->YtoPixel(pad2040single->GetY1())){
      textsizeLabelsInter = (Double_t)textSizeLabelsPixel/pad2040single->XtoPixel(pad2040single->GetX2()) ;
      textsizeFacInter = (Double_t)1./pad2040single->XtoPixel(pad2040single->GetX2()) ;
   } else {
      textsizeLabelsInter = (Double_t)textSizeLabelsPixel/pad2040single->YtoPixel(pad2040single->GetY1());
      textsizeFacInter = (Double_t)1./pad2040single->YtoPixel(pad2040single->GetY1());
   }
   if (pad4080single->XtoPixel(pad4080single->GetX2()) < pad4080single->YtoPixel(pad4080single->GetY1())){
      textsizeLabelsBottom = (Double_t)textSizeLabelsPixel/pad4080single->XtoPixel(pad4080single->GetX2()) ;
      textsizeFacBottom = (Double_t)1./pad4080single->XtoPixel(pad4080single->GetX2()) ;
   } else {
      textsizeLabelsBottom = (Double_t)textSizeLabelsPixel/pad4080single->YtoPixel(pad4080single->GetY1());
      textsizeFacBottom = (Double_t)1./pad4080single->YtoPixel(pad4080single->GetY1());
   }

   pad0020single->cd();
   pad0020single->SetLogx();
   SetStyleHistoTH1ForGraphs(dummy, "#it{p}_{T} (GeV/#it{c})","(#gamma_{inc}/#pi^{0})/(#gamma_{decay}/#pi^{0})", 0.85*textsizeLabelsTop,textsizeLabelsTop, 0.85*textsizeLabelsTop,textsizeLabelsTop, 1,0.5/(textsizeFacTop*margin), 515, 409);
   dummy->DrawCopy();
   ONE2->Draw("same");
   DR_PCM_syst_502->SetMarkerSize(markerSize);
   DR_PCM_syst_502->Draw("pZ2");
   DR_PCM_stat_502->Draw("e1same");
   DR_PHOS_syst_502->SetMarkerSize(markerSize);
   DR_PHOS_syst_502->Draw("pZ2");
   DR_PHOS_stat_502->Draw("e1same");

   TLegend* legendDoubleRatioSingleAllCent = GetAndSetLegend(0.15,0.62,5,1,"");
   legendDoubleRatioSingleAllCent->SetTextSize(textsizeLabelsTop);
   legendDoubleRatioSingleAllCent->AddEntry(DR_PCM_syst_502,"Photon Conversion Method" ,"pf");
   legendDoubleRatioSingleAllCent->AddEntry(DR_PHOS_syst_502,"Calorimteter PHOS","pf");
   legendDoubleRatioSingleAllCent->Draw();

   DrawSystem(0.15,0.83,option.Data(),"0-20%",textsizeLabelsTop);

   pad2040single->cd();
   pad2040single->SetLogx();
   SetStyleHistoTH1ForGraphs(dummy, "#it{p}_{T} (GeV/#it{c})","(#gamma_{inc}/#pi^{0})/(#gamma_{decay}/#pi^{0})", 0.85*textsizeLabelsInter,textsizeLabelsInter, 0.85*textsizeLabelsInter,textsizeLabelsInter, 1,0.5/(textsizeFacInter*margin), 515, 409);
   dummy->DrawCopy();
   ONE2->Draw("same");
   DR_PCM_syst_524->SetMarkerSize(markerSize);
   DR_PCM_syst_524->Draw("pZ2");
   DR_PCM_stat_524->Draw("e1same");
   DR_PHOS_syst_524->SetMarkerSize(markerSize);
   DR_PHOS_syst_524->Draw("peZ2");
   DR_PHOS_stat_524->Draw("e1same");

   DrawSystem(0.15,0.85,option.Data(),"20-40%",textsizeLabelsInter);

   pad4080single->cd();
   pad4080single->SetLogx();
   SetStyleHistoTH1ForGraphs(dummy, "#it{p}_{T} (GeV/#it{c})","(#gamma_{inc}/#pi^{0})/(#gamma_{decay}/#pi^{0})", 0.85*textsizeLabelsBottom,textsizeLabelsBottom, 0.85*textsizeLabelsBottom,textsizeLabelsBottom, 1,0.5/(textsizeFacBottom*margin), 515, 409);
   dummy->DrawCopy();
   ONE2->Draw("same");
   DR_PCM_syst_548->SetMarkerSize(markerSize);
   DR_PCM_syst_548->Draw("pZ2");
   DR_PCM_stat_548->Draw("e1same");
   DR_PHOS_syst_548->SetMarkerSize(markerSize);
   DR_PHOS_syst_548->Draw("peZ2");
   DR_PHOS_stat_548->Draw("e1same");
   DrawSystem(0.15,0.88,option.Data(),"40-80%",textsizeLabelsBottom);

   canvasDoubleRatioSingleAllCent->Print(Form("%s/DoubleRatioSingle_allCent.eps",outputDir.Data()));


   Double_t arrayBoundsXIndMeasRatio[7];
   Double_t arrayBoundsYIndMeasRatio[7];
   Double_t relativeMarginsIndMeasRatioX[7];
   Double_t relativeMarginsIndMeasRatioY[7];
   ReturnCorrectValuesForCanvasScaling(1200,2800, 1, 6,0.1, 0.005, 0.005,0.04,arrayBoundsXIndMeasRatio,arrayBoundsYIndMeasRatio,relativeMarginsIndMeasRatioX,relativeMarginsIndMeasRatioY);

   TCanvas *canvasDoubleRatioCombinedAllCentPlusIndividual = GetAndSetCanvas("canvasDoubleRatioCombinedAllCentPlusIndividual",0.15,0.1, 1200 ,2800);
   TPad* pad0020Ind = new TPad("pad0020", "",arrayBoundsXIndMeasRatio[0], arrayBoundsYIndMeasRatio[1],arrayBoundsXIndMeasRatio[1], arrayBoundsYIndMeasRatio[0],-1, -1, -2);
   DrawGammaPadSettings( pad0020Ind, relativeMarginsIndMeasRatioX[0], relativeMarginsIndMeasRatioX[2], relativeMarginsIndMeasRatioY[0], relativeMarginsIndMeasRatioY[1]);
   pad0020Ind->Draw();
   TPad* pad2040Ind = new TPad("pad2040", "",arrayBoundsXIndMeasRatio[0], arrayBoundsYIndMeasRatio[2],arrayBoundsXIndMeasRatio[1], arrayBoundsYIndMeasRatio[1],-1, -1, -2);
   DrawGammaPadSettings( pad2040Ind, relativeMarginsIndMeasRatioX[0], relativeMarginsIndMeasRatioX[2], relativeMarginsIndMeasRatioY[1], relativeMarginsIndMeasRatioY[1]);
   pad2040Ind->Draw();
   TPad* pad4080Ind = new TPad("pad4080","",arrayBoundsXIndMeasRatio[0], arrayBoundsYIndMeasRatio[3],arrayBoundsXIndMeasRatio[1], arrayBoundsYIndMeasRatio[2],-1, -1, -2);
   DrawGammaPadSettings( pad4080Ind,relativeMarginsIndMeasRatioX[0], relativeMarginsIndMeasRatioX[2], relativeMarginsIndMeasRatioY[1], relativeMarginsIndMeasRatioY[1]);
   pad4080Ind->Draw();
   TPad* pad0020Comb = new TPad("pad0020","",arrayBoundsXIndMeasRatio[0], arrayBoundsYIndMeasRatio[4],arrayBoundsXIndMeasRatio[1], arrayBoundsYIndMeasRatio[3],-1, -1, -2);
   DrawGammaPadSettings( pad0020Comb,relativeMarginsIndMeasRatioX[0], relativeMarginsIndMeasRatioX[2], relativeMarginsIndMeasRatioY[1], relativeMarginsIndMeasRatioY[1]);
   pad0020Comb->Draw();
   TPad* pad2040Comb = new TPad("pad2040","",arrayBoundsXIndMeasRatio[0], arrayBoundsYIndMeasRatio[5],arrayBoundsXIndMeasRatio[1], arrayBoundsYIndMeasRatio[4],-1, -1, -2);
   DrawGammaPadSettings( pad2040Comb,relativeMarginsIndMeasRatioX[0], relativeMarginsIndMeasRatioX[2], relativeMarginsIndMeasRatioY[1], relativeMarginsIndMeasRatioY[1]);
   pad2040Comb->Draw();
   TPad* pad4080Comb = new TPad("pad4080","",arrayBoundsXIndMeasRatio[0], arrayBoundsYIndMeasRatio[6],arrayBoundsXIndMeasRatio[1], arrayBoundsYIndMeasRatio[5],-1, -1, -2);
   DrawGammaPadSettings( pad4080Comb,relativeMarginsIndMeasRatioX[0], relativeMarginsIndMeasRatioX[2], relativeMarginsIndMeasRatioY[1], relativeMarginsIndMeasRatioY[2]);
   pad4080Comb->Draw();

   textSizeLabelsPixel = 40;
   margin = 0.4*1200;
   textsizeLabelsInter = 0;
   textsizeFacInter = 0;
   textsizeLabelsTop = 0;
   textsizeFacTop = 0;
   textsizeLabelsBottom = 0;
   textsizeFacBottom = 0;

   if (pad0020Ind->XtoPixel(pad0020Ind->GetX2()) < pad0020Ind->YtoPixel(pad0020Ind->GetY1())){
      textsizeLabelsTop = (Double_t)textSizeLabelsPixel/pad0020Ind->XtoPixel(pad0020Ind->GetX2()) ;
      textsizeFacTop = (Double_t)1./pad0020Ind->XtoPixel(pad0020Ind->GetX2()) ;
   } else {
      textsizeLabelsTop = (Double_t)textSizeLabelsPixel/pad0020Ind->YtoPixel(pad0020Ind->GetY1());
      textsizeFacTop = (Double_t)1./pad0020Ind->YtoPixel(pad0020Ind->GetY1());
   }
   if (pad2040Ind->XtoPixel(pad2040Ind->GetX2()) < pad2040Ind->YtoPixel(pad2040Ind->GetY1())){
      textsizeLabelsInter = (Double_t)textSizeLabelsPixel/pad2040Ind->XtoPixel(pad2040Ind->GetX2()) ;
      textsizeFacInter = (Double_t)1./pad2040Ind->XtoPixel(pad2040Ind->GetX2()) ;
   } else {
      textsizeLabelsInter = (Double_t)textSizeLabelsPixel/pad2040Ind->YtoPixel(pad2040Ind->GetY1());
      textsizeFacInter = (Double_t)1./pad2040Ind->YtoPixel(pad2040Ind->GetY1());
   }
   if (pad4080Comb->XtoPixel(pad4080Comb->GetX2()) < pad4080Comb->YtoPixel(pad4080Comb->GetY1())){
      textsizeLabelsBottom = (Double_t)textSizeLabelsPixel/pad4080Comb->XtoPixel(pad4080Comb->GetX2()) ;
      textsizeFacBottom = (Double_t)1./pad4080Comb->XtoPixel(pad4080Comb->GetX2()) ;
   } else {
      textsizeLabelsBottom = (Double_t)textSizeLabelsPixel/pad4080Comb->YtoPixel(pad4080Comb->GetY1());
      textsizeFacBottom = (Double_t)1./pad4080Comb->YtoPixel(pad4080Comb->GetY1());
   }



   cout << textsizeLabelsInter << endl;

   gStyle->SetEndErrorSize(4);

   dummy->GetYaxis()->SetRangeUser(0.76, 1.94);
   dummy->GetXaxis()->SetRangeUser(doubleRatioX[0],12);

   pad0020Ind->cd();
   pad0020Ind->SetLogx();
   SetStyleHistoTH1ForGraphs(dummy, "#it{p}_{T} (GeV/#it{c})","(#gamma_{inc}/#pi^{0})/(#gamma_{decay}/#pi^{0})", 0.85*textsizeLabelsTop,textsizeLabelsTop, 0.85*textsizeLabelsTop,textsizeLabelsTop, 1,0.5/(textsizeFacTop*margin), 515, 409);

   dummy->DrawCopy();
   ONE2->Draw("same");
   DR_PCM_syst_502->SetMarkerSize(markerSize);
   DR_PCM_stat_502->SetMarkerSize(markerSize);
   DR_PCM_syst_502->Draw("pZ2");
   DR_PCM_stat_502->Draw("E1same");
   DR_PHOS_syst_502->SetMarkerSize(markerSize);
   DR_PHOS_stat_502->SetMarkerSize(markerSize);
   DR_PHOS_syst_502->Draw("pZ2");
   DR_PHOS_stat_502->Draw("E1same");

   TLegend* legendDoubleRatioIndAllCent = GetAndSetLegend(0.15,0.62,5,1,"");
   legendDoubleRatioIndAllCent->SetTextSize(textsizeLabelsTop);
   legendDoubleRatioIndAllCent->AddEntry(DR_PCM_syst_502,"Photon Conversion Method" ,"pf");
   legendDoubleRatioIndAllCent->AddEntry(DR_PHOS_syst_502,"Calorimteter PHOS","pf");
   legendDoubleRatioIndAllCent->Draw();

   DrawSystem(0.15,0.83,option.Data(),"0-20%",textsizeLabelsTop);

   pad2040Ind->cd();
   pad2040Ind->SetLogx();
   SetStyleHistoTH1ForGraphs(dummy, "#it{p}_{T} (GeV/#it{c})","(#gamma_{inc}/#pi^{0})/(#gamma_{decay}/#pi^{0})", 0.85*textsizeLabelsInter,textsizeLabelsInter, 0.85*textsizeLabelsInter,textsizeLabelsInter, 1,0.5/(textsizeFacInter*margin), 515, 409);
   dummy->DrawCopy();
   ONE2->Draw("same");
   DR_PCM_syst_524->SetMarkerSize(markerSize);
   DR_PCM_stat_524->SetMarkerSize(markerSize);
   DR_PCM_syst_524->Draw("pZ2");
   DR_PCM_stat_524->Draw("E1same");
   DR_PHOS_syst_524->SetMarkerSize(markerSize);
   DR_PHOS_stat_524->SetMarkerSize(markerSize);
   DR_PHOS_syst_524->Draw("peZ2");
   DR_PHOS_stat_524->Draw("E1same");

   DrawSystem(0.15,0.85,option.Data(),"20-40%",textsizeLabelsInter);

   pad4080Ind->cd();
   pad4080Ind->SetLogx();
   dummy->DrawCopy();
   ONE2->Draw("same");
   DR_PCM_syst_548->SetMarkerSize(markerSize);
   DR_PCM_stat_548->SetMarkerSize(markerSize);
   DR_PCM_syst_548->Draw("pZ2");
   DR_PCM_stat_548->Draw("E1same");
   DR_PHOS_syst_548->SetMarkerSize(markerSize);
   DR_PHOS_stat_548->SetMarkerSize(markerSize);
   DR_PHOS_syst_548->Draw("peZ2");
   DR_PHOS_stat_548->Draw("E1same");
   DrawSystem(0.15,0.88,option.Data(),"40-80%",textsizeLabelsInter);

   pad0020Comb->cd();
   pad0020Comb->SetLogx();
   dummy->DrawCopy();
   ONE2->Draw("same");
   DR_Comb_syst_502->SetMarkerSize(markerSize);
   DR_Comb_stat_502->SetMarkerSize(markerSize);
   DR_Comb_syst_502->Draw("pZ2");
   DR_Comb_stat_502->Draw("pe1");
   DrawSystem(0.15,0.83,option.Data(),"Combined 0-20%",textsizeLabelsInter);

   pad2040Comb->cd();
   pad2040Comb->SetLogx();
   dummy->DrawCopy();
   ONE2->Draw("same");
   DR_Comb_syst_524->SetMarkerSize(markerSize);
   DR_Comb_stat_524->SetMarkerSize(markerSize);
   DR_Comb_syst_524->Draw("pZ2");
   DR_Comb_stat_524->Draw("pe1");
   DrawSystem(0.15,0.85,option.Data(),"Combined 20-40%",textsizeLabelsInter);

   pad4080Comb->cd();
   pad4080Comb->SetLogx();
   SetStyleHistoTH1ForGraphs(dummy, "#it{p}_{T} (GeV/#it{c})","(#gamma_{inc}/#pi^{0})/(#gamma_{decay}/#pi^{0})", 0.85*textsizeLabelsBottom,textsizeLabelsBottom, 0.85*textsizeLabelsBottom,textsizeLabelsBottom, 1,0.5/(textsizeFacBottom*margin), 509, 409);
   dummy->DrawCopy();
   ONE2->Draw("same");
   DR_Comb_syst_548->SetMarkerSize(markerSize);
   DR_Comb_stat_548->SetMarkerSize(markerSize);
   DR_Comb_syst_548->Draw("pZ2");
   DR_Comb_stat_548->Draw("pe1");
   DrawSystem(0.15,0.88,option.Data(),"Combined 40-80%",textsizeLabelsBottom);
   canvasDoubleRatioCombinedAllCentPlusIndividual->Print(Form("%s/DoubleRatioCombinedPlusIndividual_allCent.eps",outputDir.Data()));


   TGraphAsymmErrors *DP_sub_syst_502 = (TGraphAsymmErrors*)CombinedDirectPhotonsPlot->Get("DP_sub_syst_502");
   TGraphAsymmErrors *DP_sub_stat_502 = (TGraphAsymmErrors*)CombinedDirectPhotonsPlot->Get("DP_sub_stat_502");
   TGraphAsymmErrors *DP_sub_syst_524 = (TGraphAsymmErrors*)CombinedDirectPhotonsPlot->Get("DP_sub_syst_524");
   TGraphAsymmErrors *DP_sub_stat_524 = (TGraphAsymmErrors*)CombinedDirectPhotonsPlot->Get("DP_sub_stat_524");

   TF1 *Exp_sub_Syst_502 = (TF1*)CombinedDirectPhotonsPlot->Get("Exp_sub_Syst_502");
   TF1 *Exp_sub_Stat_502 = (TF1*)CombinedDirectPhotonsPlot->Get("Exp_sub_Stat_502");
   TF1 *Exp_sub_Syst_524 = (TF1*)CombinedDirectPhotonsPlot->Get("Exp_sub_Syst_524");
   TF1 *Exp_sub_Stat_524 = (TF1*)CombinedDirectPhotonsPlot->Get("Exp_sub_Stat_524");



   Double_t arrayBoundsXSubNLORatio[7];
   Double_t arrayBoundsYSubNLORatio[7];
   Double_t relativeMarginsSubNLORatioX[7];
   Double_t relativeMarginsSubNLORatioY[7];
   ReturnCorrectValuesForCanvasScaling(2400,1200, 2, 1,0.06, 0.005, 0.005,0.10,arrayBoundsXSubNLORatio,arrayBoundsYSubNLORatio,relativeMarginsSubNLORatioX,relativeMarginsSubNLORatioY);

   TCanvas *canvasDirectPhotonSubAllCent = GetAndSetCanvas("canvasDirectPhotonSubAllCent",0.06,0.10, 2400 ,1200);
   TPad* pad0020sub = new TPad("pad0020", "",arrayBoundsXSubNLORatio[0], arrayBoundsYSubNLORatio[0],arrayBoundsXSubNLORatio[1], arrayBoundsYSubNLORatio[1],-1, -1, -2);
   DrawGammaPadSettings( pad0020sub, relativeMarginsSubNLORatioX[0], relativeMarginsSubNLORatioX[1], relativeMarginsSubNLORatioY[0], relativeMarginsSubNLORatioY[2]);
   pad0020sub->Draw();
   TPad* pad2040sub = new TPad("pad2040", "",arrayBoundsXSubNLORatio[1], arrayBoundsYSubNLORatio[0],arrayBoundsXSubNLORatio[2], arrayBoundsYSubNLORatio[1],-1, -1, -2);
   DrawGammaPadSettings( pad2040sub, relativeMarginsSubNLORatioX[1], relativeMarginsSubNLORatioX[2], relativeMarginsSubNLORatioY[0], relativeMarginsSubNLORatioY[2]);
   pad2040sub->Draw();

   textSizeLabelsPixel = 40;
   margin = 0.17*1200;
   textsizeLabelsInter = 0;
   textsizeFacInter = 0;
   textsizeLabelsTop = 0;
   textsizeFacTop = 0;
   textsizeLabelsBottom = 0;
   textsizeFacBottom = 0;

   if (pad0020sub->XtoPixel(pad0020sub->GetX2()) < pad0020sub->YtoPixel(pad0020sub->GetY1())){
      textsizeLabelsTop = (Double_t)textSizeLabelsPixel/pad0020sub->XtoPixel(pad0020sub->GetX2()) ;
      textsizeFacTop = (Double_t)1./pad0020sub->XtoPixel(pad0020sub->GetX2()) ;
   } else {
      textsizeLabelsTop = (Double_t)textSizeLabelsPixel/pad0020sub->YtoPixel(pad0020sub->GetY1());
      textsizeFacTop = (Double_t)1./pad0020sub->YtoPixel(pad0020sub->GetY1());
   }
   if (pad2040sub->XtoPixel(pad2040sub->GetX2()) < pad2040sub->YtoPixel(pad2040sub->GetY1())){
      textsizeLabelsBottom = (Double_t)textSizeLabelsPixel/pad2040sub->XtoPixel(pad2040sub->GetX2()) ;
      textsizeFacBottom = (Double_t)1./pad2040sub->XtoPixel(pad2040sub->GetX2()) ;
   } else {
      textsizeLabelsBottom = (Double_t)textSizeLabelsPixel/pad2040sub->YtoPixel(pad2040sub->GetY1());
      textsizeFacBottom = (Double_t)1./pad2040sub->YtoPixel(pad2040sub->GetY1());
   }

   dummy->GetYaxis()->SetRangeUser(2e-5, 13);
   dummy->GetXaxis()->SetRangeUser(.6,4.2);


   pad0020sub->cd();
   pad0020sub->SetLogy();
   SetStyleHistoTH1ForGraphs(dummy, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c}^{2})", 0.85*textsizeLabelsTop,textsizeLabelsTop, 0.85*textsizeLabelsTop,textsizeLabelsTop, 1./*0.21/(textsizeFacTop*margin)*/,0.23/(textsizeFacTop*margin), 509, 409);
   dummy->GetXaxis()->SetLabelOffset(0.007);
   dummy->DrawCopy();
   // TLegend* legendDirectPhotonNLOSubtractedCombinedSingle = GetAndSetLegend(0.23,0.81,2.0,1,"");
   // legendDirectPhotonNLOSubtractedCombinedSingle->SetTextSize(textsizeLabelsTop);
   // legendDirectPhotonNLOSubtractedCombinedSingle->AddEntry(Exp_sub_Stat_502,Form("T = %.0f #pm %.0f^{stat} #pm %0.0f^{syst} MeV",Exp_sub_Stat_502->GetParameter(1)*1000,Exp_sub_Stat_502->GetParError(1)*1000,Exp_sub_Syst_502->GetParError(1)*1000),"l");
   // legendDirectPhotonNLOSubtractedCombinedSingle->AddEntry((TObject*)0,Form("%.1f < #it{p}_{T} < %.1f GeV/#it{c}",fitRangeLow,fitRangeHigh),"");
   // legendDirectPhotonNLOSubtractedCombinedSingle->Draw();

   TLegend* legendDirectPhotonNLOSubtractedCombinedSingleA = GetAndSetLegend(0.15,0.17,2.0,1,"");
   legendDirectPhotonNLOSubtractedCombinedSingleA->SetTextSize(textsizeLabelsTop);
   //legendDirectPhotonNLOSubtractedCombinedSingleC->AddEntry(DP_sub_syst_502,"  0-20%","pf");
   legendDirectPhotonNLOSubtractedCombinedSingleA->AddEntry(Exp_sub_Stat_502,Form("T = %.0f #pm %.0f^{ stat} #pm %0.0f^{ syst} MeV",Exp_sub_Stat_502->GetParameter(1)*1000,Exp_sub_Stat_502->GetParError(1)*1000,Exp_sub_Syst_502->GetParError(1)*1000),"l");
   legendDirectPhotonNLOSubtractedCombinedSingleA->AddEntry((TObject*)0,Form("%.1f < #it{p}_{T} < %.1f GeV/#it{c}",fitRangeLow,fitRangeHigh),"");
   
   legendDirectPhotonNLOSubtractedCombinedSingleA->Draw();
   DrawSystem(0.55,0.92,option.Data(),"0-20%",textsizeLabelsTop);
   DP_sub_stat_502->Draw("pe1");
   DP_sub_syst_502->Draw("Z2");


   if(graphUpperLimitForPointAboveButError0020CombinedSubNLO)graphUpperLimitForPointAboveButError0020CombinedSubNLO->Draw("||same");
   if(graphUpperLimitPointsAndErrorBelow00020Combined)graphUpperLimitPointsAndErrorBelow00020Combined->Draw("||same");
   if(graphArrorErrorBelow00020CombinedSubNLO)graphArrorErrorBelow00020CombinedSubNLO->Draw(">same");
   Exp_sub_Stat_502->Draw("same");
   
   gPad->RedrawAxis();
   pad2040sub->cd();
   pad2040sub->SetLogy();
   SetStyleHistoTH1ForGraphs(dummy, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c}^{2})", 0.85*textsizeLabelsBottom,textsizeLabelsBottom, 0.85*textsizeLabelsBottom,textsizeLabelsBottom, 1./*0.21/(textsizeFacBottom*margin)*/,0.45/(textsizeFacBottom*margin), 509, 409);
   dummy->GetXaxis()->SetLabelOffset(0.005);
   cout<<(textsizeFacBottom*margin)<<endl;
   dummy->DrawCopy();

   DP_sub_stat_524->Draw("pe1");
   DP_sub_syst_524->Draw("Z2");

   if(graphUpperLimitForPointAboveButError2040CombinedSubNLO){
      graphUpperLimitForPointAboveButError2040CombinedSubNLO->Draw("||same");
   }
   if(graphUpperLimitPointsAndErrorBelow02040Combined){
      graphUpperLimitPointsAndErrorBelow02040Combined->Draw("||same");
   }
   if(graphArrorErrorBelow02040CombinedSubNLO){
      graphArrorErrorBelow02040CombinedSubNLO->Draw(">same");
   }
   Exp_sub_Stat_524->Draw("same");
   TLegend* legendDirectPhotonNLOSubtractedCombinedSingleB = GetAndSetLegend(0.05,0.17,2.0,1,"");
   legendDirectPhotonNLOSubtractedCombinedSingleB->SetTextSize(textsizeLabelsTop);
   //legendDirectPhotonNLOSubtractedCombinedSingleB->AddEntry(DP_sub_syst_524,"20-40%","pf");
   legendDirectPhotonNLOSubtractedCombinedSingleB->AddEntry(Exp_sub_Stat_524,Form("T = %.0f #pm %.0f^{ stat} #pm %0.0f^{ syst} MeV",Exp_sub_Stat_524->GetParameter(1)*1000,Exp_sub_Stat_524->GetParError(1)*1000,Exp_sub_Syst_524->GetParError(1)*1000),"l");
   legendDirectPhotonNLOSubtractedCombinedSingleB->AddEntry((TObject*)0,Form("1.1 < #it{p}_{T} < %.1f GeV/#it{c}",fitRangeHigh),"");
   legendDirectPhotonNLOSubtractedCombinedSingleB->Draw();
   DrawSystem(0.45,0.92,option.Data(),"20-40%",textsizeLabelsBottom);

   // TLegend* legendDirectPhotonNLOSubtractedCombinedSingleB = GetAndSetLegend(0.01,0.8,3.5,1,"");
   // legendDirectPhotonNLOSubtractedCombinedSingleB->SetTextSize(textsizeLabelsBottom);
   // legendDirectPhotonNLOSubtractedCombinedSingleB->AddEntry(Exp_sub_Stat_524,Form("T = %.0f #pm %.0f^{stat} #pm %0.0f^{syst} MeV",Exp_sub_Stat_524->GetParameter(1)*1000,Exp_sub_Stat_524->GetParError(1)*1000,Exp_sub_Syst_524->GetParError(1)*1000),"l");
   // legendDirectPhotonNLOSubtractedCombinedSingleB->AddEntry((TObject*)0,Form("1.3 < #it{p}_{T} < %.1f GeV/#it{c}",fitRangeHigh),"");
   // legendDirectPhotonNLOSubtractedCombinedSingleB->Draw();
   
   gPad->RedrawAxis();
   canvasDirectPhotonSubAllCent->Print(Form("%s/AAADoubleRatioSingle_allCent.eps",outputDir.Data()));
   // grSystFullCombinedSubNLO->Write(Form("DP_sub_syst_%s",centralityCutNumber.Data()),TObject::kOverwrite);
   // grStatFullCombinedSubNLO->Write(Form("DP_sub_stat_%s",centralityCutNumber.Data()),TObject::kOverwrite);
   // exponetialCombinedSystSubNLO->Write(Form("Exp_sub_Syst_%s",centralityCutNumber.Data()),TObject::kOverwrite);
   // exponetialCombinedStatSubNLO->Write(Form("Exp_sub_Stat_%s",centralityCutNumber.Data()),TObject::kOverwrite);



   //FileGammaSpectrum->Open();

   /// ALL SPECTRA In One Plot

   // FileGammaSpectrum = new TFile("GammaSpectrum.root","UPDATE");

   // TGraphAsymmErrors *grSystFull0005 = (TGraphAsymmErrors*)FileGammaSpectrum->Get("DirectPhotons_points_syst_601");
   // TGraphAsymmErrors *grSystFull0510 = (TGraphAsymmErrors*)FileGammaSpectrum->Get("DirectPhotons_points_syst_612");
   // TGraphAsymmErrors *grSystFull0010 = (TGraphAsymmErrors*)FileGammaSpectrum->Get("DirectPhotons_points_syst_501");
   // TGraphAsymmErrors *grSystFull1020 = (TGraphAsymmErrors*)FileGammaSpectrum->Get("DirectPhotons_points_syst_512");
   // TGraphAsymmErrors *grSystFull0020 = (TGraphAsymmErrors*)FileGammaSpectrum->Get("DirectPhotons_points_syst_502");
   // TGraphAsymmErrors *grSystFull2040 = (TGraphAsymmErrors*)FileGammaSpectrum->Get("DirectPhotons_points_syst_524");

   // TGraphAsymmErrors *grStatFull0005 = (TGraphAsymmErrors*)FileGammaSpectrum->Get("DirectPhotons_points_stat_601");
   // TGraphAsymmErrors *grStatFull0510 = (TGraphAsymmErrors*)FileGammaSpectrum->Get("DirectPhotons_points_stat_612");
   // TGraphAsymmErrors *grStatFull0010 = (TGraphAsymmErrors*)FileGammaSpectrum->Get("DirectPhotons_points_stat_501");
   // TGraphAsymmErrors *grStatFull1020 = (TGraphAsymmErrors*)FileGammaSpectrum->Get("DirectPhotons_points_stat_512");
   // TGraphAsymmErrors *grStatFull0020 = (TGraphAsymmErrors*)FileGammaSpectrum->Get("DirectPhotons_points_stat_502");
   // TGraphAsymmErrors *grStatFull2040 = (TGraphAsymmErrors*)FileGammaSpectrum->Get("DirectPhotons_points_stat_524");

   // TGraphAsymmErrors *graphUpperLimitForPointAboveButError0005 = (TGraphAsymmErrors*)FileGammaSpectrum->Get("DirectPhotons_upperLimits_withPoint_601");
   // TGraphAsymmErrors *graphUpperLimitForPointAboveButError0510 = (TGraphAsymmErrors*)FileGammaSpectrum->Get("DirectPhotons_upperLimits_withPoint_612");
   // TGraphAsymmErrors *graphUpperLimitForPointAboveButError0010 = (TGraphAsymmErrors*)FileGammaSpectrum->Get("DirectPhotons_upperLimits_withPoint_501");
   // TGraphAsymmErrors *graphUpperLimitForPointAboveButError1020 = (TGraphAsymmErrors*)FileGammaSpectrum->Get("DirectPhotons_upperLimits_withPoint_512");
   // TGraphAsymmErrors *graphUpperLimitForPointAboveButError0020 = (TGraphAsymmErrors*)FileGammaSpectrum->Get("DirectPhotons_upperLimits_withPoint_502");
   // TGraphAsymmErrors *graphUpperLimitForPointAboveButError2040 = (TGraphAsymmErrors*)FileGammaSpectrum->Get("DirectPhotons_upperLimits_withPoint_524");

   // TGraphAsymmErrors *graphUpperLimitPointsAndErrorBelow00005 = (TGraphAsymmErrors*)FileGammaSpectrum->Get("DirectPhotons_upperLimits_noPoint_601");
   // TGraphAsymmErrors *graphUpperLimitPointsAndErrorBelow00510 = (TGraphAsymmErrors*)FileGammaSpectrum->Get("DirectPhotons_upperLimits_noPoint_612");
   // TGraphAsymmErrors *graphUpperLimitPointsAndErrorBelow00010 = (TGraphAsymmErrors*)FileGammaSpectrum->Get("DirectPhotons_upperLimits_noPoint_501");
   // TGraphAsymmErrors *graphUpperLimitPointsAndErrorBelow01020 = (TGraphAsymmErrors*)FileGammaSpectrum->Get("DirectPhotons_upperLimits_noPoint_512");
   // TGraphAsymmErrors *graphUpperLimitPointsAndErrorBelow00020 = (TGraphAsymmErrors*)FileGammaSpectrum->Get("DirectPhotons_upperLimits_noPoint_502");
   // TGraphAsymmErrors *graphUpperLimitPointsAndErrorBelow02040 = (TGraphAsymmErrors*)FileGammaSpectrum->Get("DirectPhotons_upperLimits_noPoint_524");

   // TGraphAsymmErrors *graphArrorErrorBelow00005 = (TGraphAsymmErrors*)FileGammaSpectrum->Get("DirectPhotons_Arrow_withPoint_601");
   // TGraphAsymmErrors *graphArrorErrorBelow00510 = (TGraphAsymmErrors*)FileGammaSpectrum->Get("DirectPhotons_Arrow_withPoint_612");
   // TGraphAsymmErrors *graphArrorErrorBelow00010 = (TGraphAsymmErrors*)FileGammaSpectrum->Get("DirectPhotons_Arrow_withPoint_501");
   // TGraphAsymmErrors *graphArrorErrorBelow01020 = (TGraphAsymmErrors*)FileGammaSpectrum->Get("DirectPhotons_Arrow_withPoint_512");
   // TGraphAsymmErrors *graphArrorErrorBelow00020 = (TGraphAsymmErrors*)FileGammaSpectrum->Get("DirectPhotons_Arrow_withPoint_502");
   // TGraphAsymmErrors *graphArrorErrorBelow02040 = (TGraphAsymmErrors*)FileGammaSpectrum->Get("DirectPhotons_Arrow_withPoint_524");

   // TGraphAsymmErrors *graphArrowPointErrorBelow00005 = (TGraphAsymmErrors*)FileGammaSpectrum->Get("DirectPhotons_Arrow_noPoint_601");
   // TGraphAsymmErrors *graphArrowPointErrorBelow00510 = (TGraphAsymmErrors*)FileGammaSpectrum->Get("DirectPhotons_Arrow_noPoint_612");
   // TGraphAsymmErrors *graphArrowPointErrorBelow00010 = (TGraphAsymmErrors*)FileGammaSpectrum->Get("DirectPhotons_Arrow_noPoint_501");
   // TGraphAsymmErrors *graphArrowPointErrorBelow01020 = (TGraphAsymmErrors*)FileGammaSpectrum->Get("DirectPhotons_Arrow_noPoint_512");
   // TGraphAsymmErrors *graphArrowPointErrorBelow00020 = (TGraphAsymmErrors*)FileGammaSpectrum->Get("DirectPhotons_Arrow_noPoint_502");
   // TGraphAsymmErrors *graphArrowPointErrorBelow02040 = (TGraphAsymmErrors*)FileGammaSpectrum->Get("DirectPhotons_Arrow_noPoint_524");

   // TGraphAsymmErrors *NLO0005 = (TGraphAsymmErrors*)FileGammaSpectrum->Get("NLO_601");
   // TGraphAsymmErrors *NLO0510 = (TGraphAsymmErrors*)FileGammaSpectrum->Get("NLO_612");
   // TGraphAsymmErrors *NLO0010 = (TGraphAsymmErrors*)FileGammaSpectrum->Get("NLO_501");
   // TGraphAsymmErrors *NLO1020 = (TGraphAsymmErrors*)FileGammaSpectrum->Get("NLO_512");
   // TGraphAsymmErrors *NLO0020 = (TGraphAsymmErrors*)FileGammaSpectrum->Get("NLO_502");
   // TGraphAsymmErrors *NLO2040 = (TGraphAsymmErrors*)FileGammaSpectrum->Get("NLO_524");


   // TCanvas *canvasDirectPhotonAllCent = GetAndSetCanvas("canvasDirectPhotonAllCent",0.15,0.1, 1000 ,1250);
   // canvasDirectPhotonAllCent->SetLogy();
   // dummy->GetYaxis()->SetRangeUser(2e-6,2e1);
   // dummy->GetXaxis()->SetRangeUser(0.50,11.2);
   // dummy->GetYaxis()->SetTitleOffset(1.5);
   // dummy->GetXaxis()->SetLabelOffset(0);
   // dummy->DrawCopy();


   // // if(grSystFull0005)grSystFull0005->Draw("Z2");
   // // if(grStatFull0005)grStatFull0005->Draw("PZe1");
   // // if(graphUpperLimitForPointAboveButError0005) graphUpperLimitForPointAboveButError0005->Draw("||");
   // // if(graphUpperLimitPointsAndErrorBelow00005)  graphUpperLimitPointsAndErrorBelow00005->Draw("||");
   // // if(graphArrorErrorBelow00005) graphArrorErrorBelow00005->Draw(">");
   // // if(graphArrorErrorBelow00005)graphArrorErrorBelow00005->Draw("Z");
   // // if(graphArrowPointErrorBelow00005) graphArrowPointErrorBelow00005->Draw(">");
   // // if(NLO0005) NLO0005->Draw("3l");

   // // if(grSystFull0510)grSystFull0510->Draw("Z2");
   // // if(grStatFull0510)grStatFull0510->Draw("PZe1");
   // // if(graphUpperLimitForPointAboveButError0510) graphUpperLimitForPointAboveButError0510->Draw("||");
   // // if(graphUpperLimitPointsAndErrorBelow00510)  graphUpperLimitPointsAndErrorBelow00510->Draw("||");
   // // if(graphArrorErrorBelow00510) graphArrorErrorBelow00510->Draw(">");
   // // if(graphArrorErrorBelow00510)graphArrorErrorBelow00510->Draw("Z");
   // // if(graphArrowPointErrorBelow00510) graphArrowPointErrorBelow00510->Draw(">");
   // // if(NLO0510) NLO0510->Draw("3l");


   // if(NLO0010) NLO0010->Draw("3l");
   // if(grSystFull0010)grSystFull0010->Draw("Z2");
   // if(grStatFull0010)grStatFull0010->Draw("PZe1");
   // if(graphUpperLimitForPointAboveButError0010) graphUpperLimitForPointAboveButError0010->Draw("||");
   // if(graphUpperLimitPointsAndErrorBelow00010)  graphUpperLimitPointsAndErrorBelow00010->Draw("||");
   // if(graphArrorErrorBelow00010) graphArrorErrorBelow00010->Draw(">");
   // if(graphArrorErrorBelow00010)graphArrorErrorBelow00010->Draw("Z");
   // if(graphArrowPointErrorBelow00010) graphArrowPointErrorBelow00010->Draw(">");


   // if(NLO1020) NLO1020->Draw("3l");
   // if(grSystFull1020)grSystFull1020->Draw("Z2");
   // if(grStatFull1020)grStatFull1020->Draw("PZe1");
   // if(graphUpperLimitForPointAboveButError1020) graphUpperLimitForPointAboveButError1020->Draw("||");
   // if(graphUpperLimitPointsAndErrorBelow01020)  graphUpperLimitPointsAndErrorBelow01020->Draw("||");
   // if(graphArrorErrorBelow01020) graphArrorErrorBelow01020->Draw(">");
   // if(graphArrorErrorBelow01020)graphArrorErrorBelow01020->Draw("Z");
   // if(graphArrowPointErrorBelow01020) graphArrowPointErrorBelow01020->Draw(">");


   // // if(grSystFull0020)grSystFull0020->Draw("Z2");
   // // if(grStatFull0020)grStatFull0020->Draw("PZe1");
   // // if(graphUpperLimitForPointAboveButError0020) graphUpperLimitForPointAboveButError0020->Draw("||");
   // // if(graphUpperLimitPointsAndErrorBelow00020)  graphUpperLimitPointsAndErrorBelow00020->Draw("||");
   // // if(graphArrorErrorBelow00020) graphArrorErrorBelow00020->Draw(">");
   // // if(graphArrorErrorBelow00020)graphArrorErrorBelow00020->Draw("Z");
   // // if(graphArrowPointErrorBelow00020) graphArrowPointErrorBelow00020->Draw(">");
   // // if(NLO0020) NLO0020->Draw("3l");

   // if(NLO2040) NLO2040->Draw("3l");
   // if(grSystFull2040)grSystFull2040->Draw("Z2");
   // if(grStatFull2040)grStatFull2040->Draw("PZe1");
   // if(graphUpperLimitForPointAboveButError2040) graphUpperLimitForPointAboveButError2040->Draw("||");
   // if(graphUpperLimitPointsAndErrorBelow02040)  graphUpperLimitPointsAndErrorBelow02040->Draw("||");
   // if(graphArrorErrorBelow02040) graphArrorErrorBelow02040->Draw(">");
   // if(graphArrorErrorBelow02040)graphArrorErrorBelow02040->Draw("Z");
   // if(graphArrowPointErrorBelow02040) graphArrowPointErrorBelow02040->Draw(">");



   // canvasDirectPhotonAllCent->Print(Form("%s/DirectPhotonSpectrumAllCent_%s.eps",outputDir.Data(),cent.Data()));
   // canvasDirectPhotonAllCent->Print(Form("%s/DirectPhotonSpectrumAllCent_%s.C",outputDir.Data(),cent.Data()));





   // return;











































   // TCanvas *canvasDirectPhotonWithShiftLogX = GetAndSetCanvas("canvasDirectPhotonWithShiftLogX",0.135);
   // canvasDirectPhotonWithShiftLogX->SetLogy();
   // //canvasDirectPhotonWithShiftLogX->SetLogx();
   // dummy->GetXaxis()->SetRangeUser(0.0,14.2);
   // dummy->GetYaxis()->SetTitleSize(0.05);
   // dummy->GetYaxis()->SetTitleOffset(1.2);
   // dummy->GetXaxis()->SetTitleSize(0.05);
   // dummy->GetXaxis()->SetTitleOffset(.85);
   // dummy->GetYaxis()->SetLabelSize(0.045);
   // dummy->GetXaxis()->SetLabelSize(0.045);
   // dummy->DrawCopy();
   // if(grSystFull)grSystFull->Draw("Z2");
   // if(grStatFull)grStatFull->Draw("PZe1");
   // if(graphUpperLimitForPointAboveButError)graphUpperLimitForPointAboveButError->Draw("||");
   // // //if(graphArrowConfi)graphArrowConfi->Draw(">");
   // if(graphUpperLimitPointsAndErrorBelow0)graphUpperLimitPointsAndErrorBelow0->Draw("||");
   // if(graphArrorErrorBelow0)graphArrorErrorBelow0->Draw(">");
   // // //if(graphArrorErrorBelow0)graphArrorErrorBelow0->Draw("Z");
   // if(graphArrowPointErrorBelow0)graphArrowPointErrorBelow0->Draw(">");


   // //NLO
   // NLO->Draw("p3l");
   // exponetial->SetLineStyle(2);
   // exponetial->SetLineColor(kRed+1);
   // exponetial->Draw("same");
   // TLegend* legendDirectPhotonWithXshiftLogx = GetAndSetLegend(0.4,0.50,6.5);
   // legendDirectPhotonWithXshiftLogx->AddEntry(GraphErrorsDummy,"Direct photons","pf");
   // //legendDirectPhotonWithXshiftLogx->AddEntry(NLO,"Direct photon NLO for #mu = 0.5,1.0,2.0 #it{p}_{T} (scaled pp)","l");
   // legendDirectPhotonWithXshiftLogx->AddEntry(NLO,"Direct photon NLO","l");
   // legendDirectPhotonWithXshiftLogx->AddEntry((TObject*)0, "for #mu = 0.5 to 2.0 #it{p}_{T} (scaled pp)", "");
   // //legendDirectPhotonWithXshiftLogx->AddEntry(exponetial,Form("Exponential fit: A #times exp(-#it{p}_{T}/T), T = 304 #pm 52 MeV"),"l");
   // if(GetCentralityStringA(cutSel).CompareTo("0-40%") == 0){
   //    legendDirectPhotonWithXshiftLogx->AddEntry(exponetial,"Exponential fit: A #times exp(-#it{p}_{T}/T),","l");
   //    legendDirectPhotonWithXshiftLogx->AddEntry((TObject*)0, Form("T = %.0f #pm 51^{stat+syst} MeV",exponetial->GetParameter(1)*1000,exponetial->GetParError(1)*1000), "");
   // }


   // DrawSystem(0.3,0.86,option.Data(),(GetCentralityStringA(cutSel)).Data());


   // canvasDirectPhotonWithShiftLogX->cd();
   // legendDirectPhotonWithXshiftLogx->Draw();


   // canvasDirectPhotonWithShiftLogX->Print(Form("%s/DirectPhotonSpectrumXShiftLogX_%s.eps",outputDir.Data(),cent.Data()));
   // canvasDirectPhotonWithShiftLogX->Print(Form("%s/DirectPhotonSpectrumXShiftLogX_%s.C",outputDir.Data(),cent.Data()));


   // TCanvas *canvasDirectPhotonWithShiftNoNLO = GetAndSetCanvas("canvasDirectPhotonWithShiftNoNLO");
   // canvasDirectPhotonWithShiftNoNLO->SetLogy();
   // //canvasDirectPhotonWithShiftNoNLO->SetLogx();
   // dummy->GetXaxis()->SetRangeUser(0.,14.2);
   // dummy->DrawCopy();

   // if(grSystFull)grSystFull->Draw("Z2");
   // if(grStatFull)grStatFull->Draw("PZe1");
   // if(graphUpperLimitForPointAboveButError)graphUpperLimitForPointAboveButError->Draw("||");
   // // //if(graphArrowConfi)graphArrowConfi->Draw(">");
   // if(graphUpperLimitPointsAndErrorBelow0)graphUpperLimitPointsAndErrorBelow0->Draw("||");
   // if(graphArrorErrorBelow0)graphArrorErrorBelow0->Draw(">");
   // // //if(graphArrorErrorBelow0)graphArrorErrorBelow0->Draw("Z");
   // if(graphArrowPointErrorBelow0)graphArrowPointErrorBelow0->Draw(">");


   // TLegend* legendDirectPhotonWithXshiftNoNLO = GetAndSetLegend(0.29,0.69,2.5);//0.16,0.22,1.5);
   // legendDirectPhotonWithXshiftNoNLO->AddEntry(GraphErrorsDummy,"Direct photons","pf");
   // legendDirectPhotonWithXshiftNoNLO->Draw();



   // DrawSystem(0.3,0.86,option.Data(),(GetCentralityStringA(cutSel)).Data());

   // canvasDirectPhotonWithShiftNoNLO->Print(Form("%s/DirectPhotonSpectrumXShiftNoNLO_%s.eps",outputDir.Data(),cent.Data()));
   // canvasDirectPhotonWithShiftNoNLO->Print(Form("%s/DirectPhotonSpectrumXShiftNoNLO_%s.C",outputDir.Data(),cent.Data()));



   // TFile *upperLimitFile = new TFile("DirectPhotons_0-40.root");


   // TCanvas *canvasDirectPhotonWithChargedPions = GetAndSetCanvas("canvasDirectPhotonWithChargedPions");
   // canvasDirectPhotonWithChargedPions->SetLogy();
   // dummy->GetXaxis()->SetRangeUser(0.,8.0);
   // //   dummy->GetXaxis()->SetRangeUser(0.5,5.1);
   // //dummy->GetYaxis()->SetRangeUser(5e-5,5e3);
   // dummy->DrawCopy();

   // TGraphAsymmErrors *grSystFull_cock = (TGraphAsymmErrors*) upperLimitFile->Get("1_chargecock");
   // grSystFull_cock->SetLineColor(4);
   // grSystFull_cock->SetMarkerColor(4);
   // TGraphAsymmErrors *grStatFull_cock = (TGraphAsymmErrors*) upperLimitFile->Get("2_chargecock");
   // grStatFull_cock->SetLineColor(4);
   // grStatFull_cock->SetMarkerColor(4);
   // grStatFull_cock->SetMarkerStyle(20);
   // TGraphAsymmErrors *graphUpperLimitForPointAboveButError_cock = (TGraphAsymmErrors*) upperLimitFile->Get("3_chargecock");
   // graphUpperLimitForPointAboveButError_cock->SetLineColor(4);
   // graphUpperLimitForPointAboveButError_cock->SetMarkerColor(4);
   // TGraphAsymmErrors *graphUpperLimitPointsAndErrorBelow0_cock = (TGraphAsymmErrors*) upperLimitFile->Get("4_chargecock");
   // graphUpperLimitPointsAndErrorBelow0_cock->SetLineColor(4);
   // graphUpperLimitPointsAndErrorBelow0_cock->SetMarkerColor(4);
   // TGraphAsymmErrors *graphArrorErrorBelow0_cock = (TGraphAsymmErrors*) upperLimitFile->Get("5_chargecock");
   // graphArrorErrorBelow0_cock->SetLineColor(4);
   // graphArrorErrorBelow0_cock->SetMarkerColor(4);
   // TGraphAsymmErrors *graphArrowPointErrorBelow0_cock = (TGraphAsymmErrors*) upperLimitFile->Get("6_chargecock");
   // graphArrowPointErrorBelow0_cock->SetLineColor(4);
   // graphArrowPointErrorBelow0_cock->SetMarkerColor(4);
   // TF1 *exponential_cock = (TF1*)  upperLimitFile->Get("exponential_cock");


   // if(grSystFull)grSystFull->Draw("Z2");
   // if(grStatFull)grStatFull->Draw("PZe1");
   // if(graphUpperLimitForPointAboveButError)graphUpperLimitForPointAboveButError->Draw("||");
   // // //if(graphArrowConfi)graphArrowConfi->Draw(">");
   // if(graphUpperLimitPointsAndErrorBelow0)graphUpperLimitPointsAndErrorBelow0->Draw("||");
   // if(graphArrorErrorBelow0)graphArrorErrorBelow0->Draw(">");
   // // //if(graphArrorErrorBelow0)graphArrorErrorBelow0->Draw("Z");
   // if(graphArrowPointErrorBelow0)graphArrowPointErrorBelow0->Draw(">");



   // //grSystFull_cock->Draw("Z2");
   // // grStatFull_cock->Draw("PZe1");
   // // graphUpperLimitForPointAboveButError_cock->Draw("||");
   // // graphUpperLimitPointsAndErrorBelow0_cock->Draw("||");
   // // graphArrorErrorBelow0_cock->Draw(">");
   // // graphArrowPointErrorBelow0_cock->Draw(">");


   // TGraphAsymmErrors *grSystFull_both = (TGraphAsymmErrors*) upperLimitFile->Get("1_chargeboth");
   // grSystFull_both->SetLineColor(2);
   // grSystFull_both->SetMarkerColor(2);
   // TGraphAsymmErrors *grStatFull_both = (TGraphAsymmErrors*) upperLimitFile->Get("2_chargeboth");
   // grStatFull_both->SetLineColor(2);
   // grStatFull_both->SetMarkerColor(2);
   // grStatFull_both->SetMarkerStyle(20);
   // TGraphAsymmErrors *graphUpperLimitForPointAboveButError_both = (TGraphAsymmErrors*) upperLimitFile->Get("3_chargeboth");
   // graphUpperLimitForPointAboveButError_both->SetLineColor(2);
   // graphUpperLimitForPointAboveButError_both->SetMarkerColor(2);
   // TGraphAsymmErrors *graphUpperLimitPointsAndErrorBelow0_both = (TGraphAsymmErrors*) upperLimitFile->Get("4_chargeboth");
   // graphUpperLimitPointsAndErrorBelow0_both->SetLineColor(2);
   // graphUpperLimitPointsAndErrorBelow0_both->SetMarkerColor(2);
   // TGraphAsymmErrors *graphArrorErrorBelow0_both = (TGraphAsymmErrors*) upperLimitFile->Get("5_chargeboth");
   // graphArrorErrorBelow0_both->SetLineColor(2);
   // graphArrorErrorBelow0_both->SetMarkerColor(2);
   // TGraphAsymmErrors *graphArrowPointErrorBelow0_both = (TGraphAsymmErrors*) upperLimitFile->Get("6_chargeboth");
   // graphArrowPointErrorBelow0_both->SetLineColor(2);
   // graphArrowPointErrorBelow0_both->SetMarkerColor(2);
   // //TF1 *exponential_both = (TF1*)  upperLimitFile->Get("exponential_both");

   // TF1 *exponential_both = new TF1("Exponential Fit both","[0]*exp(-x/[1])",1.2,2.2);
   // exponential_both->SetParameters(1.,300.);
   // grStatFull_both->Fit(exponential_both,"NRME+","",1.2,2.2);



   // //grSystFull_both->Draw("Z2");
   // grStatFull_both->Draw("PZe1");
   // graphUpperLimitForPointAboveButError_both->Draw("||");
   // graphUpperLimitPointsAndErrorBelow0_both->Draw("||");
   // graphArrorErrorBelow0_both->Draw(">");
   // graphArrowPointErrorBelow0_both->Draw(">");


   // exponetial->SetLineColor(1);
   // exponetial->Draw("same");
   // exponential_both->SetLineColor(2);
   // exponential_both->SetLineStyle(2);
   // exponential_both->Draw("same");
   // // exponential_cock->SetLineColor(4);
   // // exponential_cock->SetLineStyle(3);
   // // exponential_cock->Draw("same");


   // TLegend* legendDirectPhotonWithChargedPions = GetAndSetLegend(0.30,0.58,7.5);
   // legendDirectPhotonWithChargedPions->AddEntry(GraphErrorsDummy,"Direct photons","pf");
   // legendDirectPhotonWithChargedPions->AddEntry(exponetial,Form("Exponential fit: A #times exp(-#it{p}_{T}/T), T = %.0f #pm 51^{stat+syst} MeV",exponetial->GetParameter(1)*1000,exponetial->GetParError(1)*1000),"l");
   // //legendDirectPhotonWithChargedPions->AddEntry(grStatFull_cock,"Direct photon with cocktail from charged pions","pl");
   // // legendDirectPhotonWithChargedPions->AddEntry(exponential_cock,Form("Exponential Fit (cocktail from charged pions)"),"l");
   // //legendDirectPhotonWithChargedPions->AddEntry((TObject*)0,Form("T = %.0f #pm %.0f MeV",exponential_cock->GetParameter(1)*1000, exponential_cock->GetParError(1)*1000),"");
   // legendDirectPhotonWithChargedPions->AddEntry(grStatFull_both,"Direct photon with charged pions","pl");
   // legendDirectPhotonWithChargedPions->AddEntry(exponential_both,Form("Exponential Fit (charged pions)"),"l");
   // legendDirectPhotonWithChargedPions->AddEntry((TObject*)0,Form("T = %.0f #pm 49^{stat+syst} MeV",exponential_both->GetParameter(1)*1000, exponential_both->GetParError(1)*1000),"");
   // legendDirectPhotonWithChargedPions->Draw();







   // canvasDirectPhotonWithChargedPions->Print(Form("%s/DirectPhotonSpectrumChargedPions_%s.eps",outputDir.Data(),cent.Data()));








   // TGraphAsymmErrors * spectrum = (TGraphAsymmErrors*)graAllShifted->Clone("Spectrum");

   // gStyle->SetOptFit(0);
   // gStyle->SetOptTitle(0);

   // TF1 *fitHag = new TF1("Hagedorn",Hagedorn,3.5,14.,3);
   // fitHag->SetParameters(1.7e4,0.05,5.6) ;
   // fitHag->SetLineColor(kBlue) ;
   // fitHag->SetLineWidth(2) ;
   // fitHag->SetLineStyle(2) ;
   // fitHag->SetParName(0,"A") ;
   // fitHag->SetParName(1,"C") ;
   // fitHag->SetParName(2,"n") ;
   // spectrum->Fit(fitHag,"E EX0","",3.5,14.) ;

   // TF1 *fitExp = new TF1("Exponent",Exponent,0.8,2.6,2);
   // fitExp->SetParameters(50.,0.3) ;
   // fitExp->SetLineColor(kGreen+1) ;
   // fitExp->SetLineWidth(2) ;
   // fitExp->SetLineStyle(2) ;
   // fitExp->SetParName(0,"B") ;
   // fitExp->SetParName(1,"T") ;
   // spectrum->Fit(fitExp,"E EX0 +","",0.8,2.6) ;

   // TF1 *fit = new TF1("Hagedorn+exponent",FitSpectrum,0.5,14.,5);
   // fit->SetParameters(fitHag->GetParameter(0),fitHag->GetParameter(1),fitHag->GetParameter(2),
   //                    fitExp->GetParameter(0),fitExp->GetParameter(1)) ;
   // // fit->SetParameters(100.,0.1,6.,1.,0.3) ;
   // fit->SetLineColor(kRed) ;
   // fit->SetLineWidth(2) ;
   // fit->SetParName(0,"A") ;
   // fit->SetParName(1,"C") ;
   // fit->SetParName(2,"n") ;
   // fit->SetParName(3,"B") ;
   // fit->SetParName(4,"T") ;
   // spectrum->Fit(fit,"E EX0 +","",0.8,14.) ;

   // TH1F *box = spectrum->GetHistogram();
   // box->SetXTitle("#it{p}_{T} (GeV/#it{c})");
   // box->SetYTitle("#frac{1}{2#pi} #frac{d^{ 2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} ((GeV/#it{c})^{-2})");
   // box->SetTitleOffset(1.3,"Y");

   // TPaveText *fitText = new TPaveText(0.50,0.50,0.96,0.95,"NDC");
   // char title[512] ;
   // sprintf(title,"#frac{1}{2#pi} #frac{d^{ 2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} = #frac{A}{2#pi} #frac{1}{(1+#it{p}_{T}/nC)^{n}} + B exp(-#it{p}_{T}/T)");
   // fitText->AddText(title) ;
   // fitText->AddText(" ") ;
   // sprintf(title,"#chi^{2}/NDF\t = %4.1f/%d"     ,fit->GetChisquare(),fit->GetNDF()) ;
   // fitText->AddText(title) ;
   // sprintf(title,"A\t = %.3g  #pm %.3g (GeV/#it{c})^{-2}",fit->GetParameter(0),fit->GetParError(0)) ;
   // fitText->AddText(title) ;
   // sprintf(title,"C\t = %4.2f  #pm %3.2f GeV"       ,fit->GetParameter(1),fit->GetParError(1)) ;
   // fitText->AddText(title) ;
   // sprintf(title,"n\t = %4.1f  #pm %3.1f "          ,fit->GetParameter(2),fit->GetParError(2)) ;
   // fitText->AddText(title) ;
   // sprintf(title,"B\t = %.2g  #pm %.2g (GeV/#it{c})^{-2}",fit->GetParameter(3),fit->GetParError(3)) ;
   // fitText->AddText(title) ;
   // sprintf(title,"T\t =%5.2f  #pm %5.2f GeV"        ,fit->GetParameter(4),fit->GetParError(4)) ;
   // fitText->AddText(title) ;
   // fitText->SetFillColor(0) ;
   // fitText->SetTextAlign(12) ;
   // fitText->SetBorderSize(0);

   // TLegend *fitLegend = new TLegend(0.15,0.15,0.5,0.30);
   // fitLegend->SetFillColor(kWhite);
   // fitLegend->SetBorderSize(0);
   // fitLegend->AddEntry(fitHag,"Hagedorn","l");
   // fitLegend->AddEntry(fitExp,"exponent","l");
   // fitLegend->AddEntry(fit,"Hagedorn+exponent","l");

   // TCanvas * cM  = GetAndSetCanvas("cocktailCanvasSpecDirectPhtonCombinedFit");
   // cM->SetRightMargin(0.02);
   // cM->SetTopMargin(0.02);
   // cM->SetLeftMargin(0.12);
   // cM->SetFillColor(0) ;
   // cM->SetFillStyle(0) ;
   // cM->Range(0,0,1,1);
   // cM->SetBorderSize(0);
   // cM->SetLogy();
   // dummy->Draw();


   // if(grSystFull)grSystFull->Draw("Z2");
   // if(grStatFull)grStatFull->Draw("PZe1");
   // if(graphUpperLimitForPointAboveButError)graphUpperLimitForPointAboveButError->Draw("||");
   // // //if(graphArrowConfi)graphArrowConfi->Draw(">");
   // if(graphUpperLimitPointsAndErrorBelow0)graphUpperLimitPointsAndErrorBelow0->Draw("||");
   // if(graphArrorErrorBelow0)graphArrorErrorBelow0->Draw(">");
   // // //if(graphArrorErrorBelow0)graphArrorErrorBelow0->Draw("Z");
   // if(graphArrowPointErrorBelow0)graphArrowPointErrorBelow0->Draw(">");






   // // box->Draw();
   // //spectrum->Draw("AP") ;
   // fit->Draw("same");
   // fitHag->Draw("same");
   // fitExp->Draw("same");

   // fitText->Draw() ;
   // fitLegend->Draw();


   // cM->Print(Form("%s/DirectPhotonSpectrumXShiftCombinedFit_%s.eps",outputDir.Data(),cent.Data()));
   // cM->Print(Form("%s/DirectPhotonSpectrumXShiftCombinedFit_%s.C",outputDir.Data(),cent.Data()));




   // TCanvas *canvasDirectPhotonWithShiftCombined = GetAndSetCanvas("canvasDirectPhotonWithShiftCombined");
   // canvasDirectPhotonWithShiftCombined->SetLogy();
   // //canvasDirectPhotonWithShiftCombined->SetLogx();
   // dummy->GetXaxis()->SetRangeUser(0.0,8);
   // dummy->DrawCopy();



   // if(grSystFull)grSystFull->Draw("Z2");
   // if(grStatFull)grStatFull->Draw("PZe1");
   // if(graphUpperLimitForPointAboveButError)graphUpperLimitForPointAboveButError->Draw("||");
   // // //if(graphArrowConfi)graphArrowConfi->Draw(">");
   // if(graphUpperLimitPointsAndErrorBelow0)graphUpperLimitPointsAndErrorBelow0->Draw("||");
   // if(graphArrorErrorBelow0)graphArrorErrorBelow0->Draw(">");
   // // //if(graphArrorErrorBelow0)graphArrorErrorBelow0->Draw("Z");
   // if(graphArrowPointErrorBelow0)graphArrowPointErrorBelow0->Draw(">");




   // fit->Draw("same");
   // fitHag->Draw("same");
   // fitExp->Draw("same");



   // TLegend* legendDirectPhotonWithXshiftCombined = GetAndSetLegend(0.42,0.40,9.5);
   // legendDirectPhotonWithXshiftCombined->AddEntry(GraphErrorsDummy,"Direct photons","pf");
   // legendDirectPhotonWithXshiftCombined->AddEntry(fitHag,"Hagedorn","l");
   // legendDirectPhotonWithXshiftCombined->AddEntry(fitExp,"Exponential","l");
   // legendDirectPhotonWithXshiftCombined->AddEntry(fit,"Hagedorn+Exponential","l");
   // legendDirectPhotonWithXshiftCombined->AddEntry((TObject*)0, "#frac{1}{2#pi} #frac{d^{ 2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} = #frac{A}{2#pi} #frac{1}{(1+#it{p}_{T}/nC)^{n}} + B exp(-#it{p}_y{T}/T)", "");
   // //legendDirectPhotonWithXshiftCombined->AddEntry((TObject*)0, Form("T = %.0f #pm %.0f MeV",fit->GetParameter(4)*1000, fit->GetParError(4)*1000),"");
   // legendDirectPhotonWithXshiftCombined->AddEntry((TObject*)0, Form("T = 299 #pm 51^{stat+syst} MeV"),"");


   // legendDirectPhotonWithXshiftCombined->Draw();

   // if(grSystFull)grSystFull->Draw("Z2");
   // if(grStatFull)grStatFull->Draw("PZe1");
   // if(graphUpperLimitForPointAboveButError)graphUpperLimitForPointAboveButError->Draw("||");
   // // //if(graphArrowConfi)graphArrowConfi->Draw(">");
   // if(graphUpperLimitPointsAndErrorBelow0)graphUpperLimitPointsAndErrorBelow0->Draw("||");
   // if(graphArrorErrorBelow0)graphArrorErrorBelow0->Draw(">");
   // // //if(graphArrorErrorBelow0)graphArrorErrorBelow0->Draw("Z");
   // if(graphArrowPointErrorBelow0)graphArrowPointErrorBelow0->Draw(">");


   // fit->Draw("same");
   // fitHag->Draw("same");
   // fitExp->Draw("same");


   // DrawSystem(0.3,0.86,option.Data(),(GetCentralityStringA(cutSel)).Data());

   // canvasDirectPhotonWithShiftCombined->Print(Form("%s/DirectPhotonSpectrumXShiftCombined_%s.eps",outputDir.Data(),cent.Data()));
   // canvasDirectPhotonWithShiftCombined->Print(Form("%s/DirectPhotonSpectrumXShiftCombined_%s.C",outputDir.Data(),cent.Data()));



   // TCanvas *canvasDirectPhotonWithShiftBlastWave = GetAndSetCanvas("canvasDirectPhotonWithShiftBlastWave");
   // Int_t nParam = 3;
   // Double_t NN = 0.001;      // normalization
   // Double_t TT = 90.;     // T, in MeV
   // Double_t BB = 0.6;     // beta --> "velocity"
   // Double_t xMin = 0.8;
   // Double_t xMax = 2.4;

   // TF1 *FuncBWdndpt = new TF1("Blast Wave",BWdndpt,xMin, xMax, nParam);
   // FuncBWdndpt->SetParameters(NN, TT, BB);
   // FuncBWdndpt->SetParNames("norm","T","beta");
   // FuncBWdndpt->FixParameter(2,0.66);
   // //FuncBWdndpt->FixParameter(1,300);
   // FuncBWdndpt->SetLineColor(kSpring);
   // spectrum->Fit(FuncBWdndpt,"INRME+", "");

   // canvasDirectPhotonWithShiftBlastWave->SetLogy();
   // dummy->GetXaxis()->SetRangeUser(0.0,7.2);
   // dummy->DrawCopy();
   // if(grSystFull)grSystFull->Draw("Z2");
   // if(grStatFull)grStatFull->Draw("PZe1");
   // if(graphUpperLimitForPointAboveButError)graphUpperLimitForPointAboveButError->Draw("||");
   // // //if(graphArrowConfi)graphArrowConfi->Draw(">");
   // if(graphUpperLimitPointsAndErrorBelow0)graphUpperLimitPointsAndErrorBelow0->Draw("||");
   // if(graphArrorErrorBelow0)graphArrorErrorBelow0->Draw(">");
   // // //if(graphArrorErrorBelow0)graphArrorErrorBelow0->Draw("Z");
   // if(graphArrowPointErrorBelow0)graphArrowPointErrorBelow0->Draw(">");



   // FuncBWdndpt->Draw("same");


   // TLegend* legendDirectPhotonWithXshiftBlastWave = GetAndSetLegend(0.14,0.12,2.5);
   // legendDirectPhotonWithXshiftBlastWave->AddEntry(GraphErrorsDummy,"Direct photons","pf");
   // legendDirectPhotonWithXshiftBlastWave->AddEntry(FuncBWdndpt, Form("Blast wave fit, T = %.0f #pm %.0f^{stat+syst} MeV",FuncBWdndpt->GetParameter(1), FuncBWdndpt->GetParError(1)), "l");
   // legendDirectPhotonWithXshiftBlastWave->Draw();

   // canvasDirectPhotonWithShiftBlastWave->Print(Form("%s/DirectPhotonSpectrumXShiftBlastWave_%s.eps",outputDir.Data(),cent.Data()));


   // TFile *DirectPhotonsOutput = new TFile("DirectPhotons_00-40.root","RECREATE");

   // DoubleRatioRebinedFit->Write("DoubleRatio0-40_stat");
   // DoubleRatioWithFitErrors->Write("DoubleRatio0-40_syst");

   // if(grSystFull)grSystFull->Write("DirectPhotons_points_syst");
   // if(grStatFull)grStatFull->Write("DirectPhotons_points_stat");
   // if(graphUpperLimitForPointAboveButError)graphUpperLimitForPointAboveButError->Write("DirectPhotons_upperLimits_withPoint");
   // // //if(graphArrowConfi)graphArrowConfi->Draw(">");
   // if(graphUpperLimitPointsAndErrorBelow0)graphUpperLimitPointsAndErrorBelow0->Write("DirectPhotons_upperLimits_noPoint");
   // if(graphArrorErrorBelow0)graphArrorErrorBelow0->Write("DirectPhotons_Arrow_withPoint");
   // // //if(graphArrorErrorBelow0)graphArrorErrorBelow0->Draw("Z");
   // if(graphArrowPointErrorBelow0)graphArrowPointErrorBelow0->Write("DirectPhotons_Arrow_noPoint");

   // // gammaAssymError->Write("Gamma_0-40_stat");
   // // GammaErrors->Write("Gamma_0-40_syst");
   // // perGamma->Write("Gamma_40-80_stat");
   // // perGammaErrors->Write("Gamma_40-80_syst");


   // DirectPhotonsOutput->Close();

   // TFile *Pi0Output = new TFile("Conversion_Pi0_Yields_PbPb.root","UPDATE");



   // if(grSystFull)grSystFull->Write("DirectPhotons_points_syst");
   // if(grStatFull)grStatFull->Write("DirectPhotons_points_stat");
   // if(graphUpperLimitForPointAboveButError)graphUpperLimitForPointAboveButError->Write("DirectPhotons_upperLimits_withPoint");
   // // //if(graphArrowConfi)graphArrowConfi->Draw(">");
   // if(graphUpperLimitPointsAndErrorBelow0)graphUpperLimitPointsAndErrorBelow0->Write("DirectPhotons_upperLimits_noPoint");
   // if(graphArrorErrorBelow0)graphArrorErrorBelow0->Write("DirectPhotons_Arrow_withPoint");
   // // //if(graphArrorErrorBelow0)graphArrorErrorBelow0->Draw("Z");
   // if(graphArrowPointErrorBelow0)graphArrowPointErrorBelow0->Write("DirectPhotons_Arrow_noPoint");


   // fit->Write("ExpPlusHag_0040");
   // exponetial->Write("Exponential_0040");



   // Pi0Output->Close();

}

