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
#include "TRandom3.h"
#include "Math/WrappedTF1.h"
#include "Math/BrentRootFinder.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TGenPhaseSpace.h"
#include "../CommonHeaders/PlottingGammaConversionHistos.h"
#include "../CommonHeaders/PlottingGammaConversionAdditional.h"
#include "../CommonHeaders/FittingGammaConversion.h"
#include "../CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "../CommonHeaders/ConversionFunctions.h"
#include "TThread.h"

// Global definition of inputs accessible for all threads
Int_t fNEvtsGlobal;
TThread *hx[10];
TH1D* h1_ptdistSmeared[4];
TH1D* h1_ptdistSmearedReb[4];
TProfile* p1_fractionsSmeared[4];
TProfile* p1_fractionsSmearedReb[4];
TSpline3* spFractions[10][4];

TH2D *h2_ptMothervsDaughter;
TH2D *h2_asym_geom;
TH2D *h2_OpenAngle;
TH2D *h2_ClusterDistance;

TH1D *h1_ptdistribution;
TH1D *h1_ptdistributionReb;

TH1D *h1_phidistribution;
TH1D *h1_etadistribution;
TH1D *h1_ptDaughter[10];
TH1D *h1_ptGammaDaughters;
TH1D *h1_ptPiPlDaughters;
TH1D *h1_ptPiMiDaughters;
TH1D *h1_ptPiZeroDaughters;

TH1D* h1ResolProjection[10][5][1000];
TF1* ptDistribution[10];
Bool_t haveResol                 = kFALSE;
Bool_t haveFractions             = kFALSE;
TH2F* histoResolutionInput[5];
TH1F* histoFractionsInput[5];
TString labelResolHist[5];
Int_t nResolHist                 = 4;
Double_t massParticle            = 0;
TString outputlabel              = "";
TString plotLabel                = "";
TString daughterLabels[10]       = {"","","","","","","","","",""};
Color_t colors[10]               = { kRed+1, kGreen+2, kOrange+7, kBlue+2, kAzure,
                                    kYellow-8, kViolet, kCyan+2, kGray+1, kPink+2 };
Int_t nDaughters                 = 0;
Double_t* masses                 = NULL;
Int_t* pdgCodesDaughters         = NULL;

TH1D* h1_ptReweightedRecMCRebModFrac[3];
TH1D* h1_ptReweightedRecMCRebModFracRatioToStandard[3];


void *threadhandle(void *ptr)
{
   long nr = (long) ptr;
   TRandom3 randy1;
   TRandom3 randy2;
   TRandom3 randy3;
   TGenPhaseSpace event;
   Double_t rEMC               = 440.;
   for(Int_t n=0; n<fNEvtsGlobal; n++){ // this is the important loop (nEvents)
      if (n%500000 == 0)
      TThread::Printf("generated: %4.1f Mio events in thread: %ld", (Double_t)n/1e6, nr);
      // TThread::Lock();
      Double_t ptcurrent      = ptDistribution[nr]->GetRandom();
      // TThread::UnLock();

      Double_t phiCurrent     = randy1.Uniform(2*TMath::Pi());
      Double_t etaCurrent     = randy2.Uniform(-1,1);

      h1_ptdistribution->Fill(ptcurrent);
      h1_ptdistributionReb->Fill(ptcurrent);
      h1_phidistribution->Fill(phiCurrent);
      h1_etadistribution->Fill(etaCurrent);

      Double_t shift          = 0;
      for (Int_t j = 0; (j < 4 && j < 5); j++){
         Double_t ptcurrentDist      = ptcurrent;
         if (haveResol){
            // TThread::Lock();
            shift                   = h1ResolProjection[nr][j][histoResolutionInput[j]->GetXaxis()->FindBin(ptcurrent)-1]->GetRandom();
            // TThread::UnLock();
            ptcurrentDist           = (shift*ptcurrent)+ptcurrent;
         }
         h1_ptdistSmeared[j]->Fill(ptcurrentDist);
         h1_ptdistSmearedReb[j]->Fill(ptcurrentDist);
         Double_t fraction           = 0;
         if (haveFractions){
            fraction                = spFractions[nr][j]->Eval(ptcurrentDist);
         }
         p1_fractionsSmeared[j]->Fill(ptcurrentDist,fraction);
         p1_fractionsSmearedReb[j]->Fill(ptcurrentDist,fraction);
      }

      TLorentzVector particle(0.0, 0.0, 0, massParticle);
      particle.SetPtEtaPhiM(ptcurrent, etaCurrent, phiCurrent, massParticle);
      event.SetDecay(particle, nDaughters, masses);

      Double_t weight         = event.Generate();
      Double_t ptDaughter[nDaughters];
      TLorentzVector* pDaughter[nDaughters];
      for (Int_t i = 0; i < nDaughters; i++){
         pDaughter[i]                = event.GetDecay(i); // these are my daughters !!
         ptDaughter[i]               = pDaughter[i]->Pt();

         h1_ptDaughter[i]->Fill(ptDaughter[i]);
         if (pdgCodesDaughters[i] == 22)
            h1_ptGammaDaughters->Fill(ptDaughter[i]);
         if (pdgCodesDaughters[i] == 111)
            h1_ptPiZeroDaughters->Fill(ptDaughter[i]);
         if (pdgCodesDaughters[i] == 211)
            h1_ptPiPlDaughters->Fill(ptDaughter[i]);
         if (pdgCodesDaughters[i] == -211)
            h1_ptPiMiDaughters->Fill(ptDaughter[i]);

         if (i == 0)
            h2_ptMothervsDaughter->Fill(ptcurrent,ptDaughter[i]);
      }
      if (nDaughters > 1){
         h2_asym_geom->Fill((ptDaughter[0]-ptDaughter[1])/(ptDaughter[0]+ptDaughter[1]), ptcurrent);
         Double_t openangle      = pDaughter[0]->Angle(pDaughter[1]->Vect());
         Double_t distanceCl     = TMath::Tan(openangle/2)*rEMC*2;

         h2_OpenAngle->Fill(ptcurrent,openangle);
         h2_ClusterDistance->Fill(ptcurrent,distanceCl);
      }
   }
   return 0;
}

// exemplary call for this macro (runs 5 threads which each generates 2e8 events)
// root -b -l -q -x 'ToyModels/threads.C(2e8,5,1,"pPb_8TeV",10,200,"pdf","GridOutput/LHC18b9bc_GammaCaloMerged_3430_all.root","80010123_4117947050032200000_4117947050022700001_0163300000000000",10,"","8008d123_4117947050032200000_4117947050022700001_0163300000000000/pPb_8TeV/Pi0_MC_GammaMergedCorrection_8008d123_4117947050032200000_4117947050022700001_0163300000000000.root")'
void NeutralMesonDecay_Threaded(
   Int_t nEvts                 = 1000000,
   const Int_t numThreads      = 5,
   Int_t particle              = 1,
   TString energy              = "2.76TeV",
   Double_t minPt              = 2,
   Double_t maxPt              = 50,
   TString suffix              = "eps",
   TString fileName            = "",
   TString cutSting            = "",
   Int_t mode                  = 10,
   TString subMode             = "",
   TString fileNameFractions   = ""
){
   fNEvtsGlobal = nEvts;
   //*************************************************************************************************
   //******************************** Set Style settings globally ************************************
   //*************************************************************************************************
   gROOT->Reset();
   gROOT->SetStyle("Plain");

   StyleSettingsThesis(suffix);
   SetPlotStyle();
   TString fCollisionSystem            = ReturnFullCollisionsSystem(energy);
   TString fCollisionSystenWrite       = ReturnCollisionEnergyOutputString(energy);
   TString fOutputDir                  = Form("%s/%s/FractionsCheck/",fCollisionSystenWrite.Data(),suffix.Data());
   gSystem->Exec("mkdir -p "+fOutputDir);

   //*************************************************************************************************
   //*************************** Initialize pt spectra for particles**********************************
   //*************************************************************************************************

   if (particle == 1){
      massParticle        = TDatabasePDG::Instance()->GetParticle(111)->Mass();
      if (energy.CompareTo("2.76TeV") == 0){
         for (Int_t tt = 0; tt < 10; tt++){ // one entry per thread for safety
            ptDistribution[tt]      = FitObject("tcm","ptDistribution","Pi0",NULL,0.4,50.);
            ptDistribution[tt]->SetParameter(0,3.04120116616420654e+11);
            ptDistribution[tt]->SetParameter(1,1.33176621422733510e-01);
            ptDistribution[tt]->SetParameter(2,2.73953300141987877e+10);
            ptDistribution[tt]->SetParameter(3,5.64798456522145997e-01);
            ptDistribution[tt]->SetParameter(4,3.19909580984167752e+00);
            ptDistribution[tt]->SetRange(minPt,maxPt);
         }
      } else if (energy.CompareTo("pPb_8TeV") == 0){
         for (Int_t tt = 0; tt < 10; tt++){ // one entry per thread for safety
            ptDistribution[tt]      = FitObject("tcm",Form("ptDistribution%d",tt),"Pi0",NULL,10,200.);
            ptDistribution[tt]->SetParameter(0,59.7673972756);
            ptDistribution[tt]->SetParameter(1,0.1133121714);
            ptDistribution[tt]->SetParameter(2,4.2802668896);
            ptDistribution[tt]->SetParameter(3,0.5676437241);
            ptDistribution[tt]->SetParameter(4,2.9785817202);
            ptDistribution[tt]->SetRange(minPt,maxPt);
         }
      } else if (energy.CompareTo("8TeV") == 0){
         for (Int_t tt = 0; tt < 10; tt++){ // one entry per thread for safety
            ptDistribution[tt]      = FitObject("tcm",Form("ptDistribution%d",tt),"Pi0",NULL,10,200.);
            ptDistribution[tt]->SetParameter(0,709.7622755997);
            ptDistribution[tt]->SetParameter(1,0.0022980359);
            ptDistribution[tt]->SetParameter(2,0.0758656254);
            ptDistribution[tt]->SetParameter(3,0.8831101616);
            ptDistribution[tt]->SetParameter(4,3.0568183147);
            ptDistribution[tt]->SetRange(minPt,maxPt);
         }
      } else {
         cout << "ERROR: undefined energy for pi0" << endl;
         return;
      }
      outputlabel         = "Pi0";
      plotLabel           = "#pi^{0}";
      daughterLabels[0]   = "#gamma_{1}";
      daughterLabels[1]   = "#gamma_{2}";
      nDaughters          = 2;
      pdgCodesDaughters   = new Int_t[nDaughters];
      pdgCodesDaughters[0]= 22;
      pdgCodesDaughters[1]= 22;
      masses              = new Double_t[nDaughters];
      masses[0]           = TDatabasePDG::Instance()->GetParticle(pdgCodesDaughters[0])->Mass();;
      masses[1]           = TDatabasePDG::Instance()->GetParticle(pdgCodesDaughters[1])->Mass();;
   } else if (particle == 2){
      massParticle        = TDatabasePDG::Instance()->GetParticle(221)->Mass();
      if (energy.CompareTo("2.76TeV") == 0){
         for (Int_t tt = 0; tt < 10; tt++){ // one entry per thread for safety
            ptDistribution[tt]      = FitObject("tcm","ptDistribution","Eta",NULL,0.4,50.);
            ptDistribution[tt]->SetParameter(0,1.71597713217546654e+10);
            ptDistribution[tt]->SetParameter(1,1.51262659810865535e-01);
            ptDistribution[tt]->SetParameter(2,1.41711438002461696e+09);
            ptDistribution[tt]->SetParameter(3,8.49829150135323452e-01);
            ptDistribution[tt]->SetParameter(4,3.32752667783400913e+00);
            ptDistribution[tt]->SetRange(minPt,maxPt);
         }
      } else {
         cout << "ERROR: undefined energy for eta" << endl;
         return;
      }
      outputlabel         = "Eta";
      plotLabel           = "#eta";
      daughterLabels[0]   = "#gamma_{1}";
      daughterLabels[1]   = "#gamma_{2}";
      nDaughters          = 2;
      pdgCodesDaughters   = new Int_t[nDaughters];
      pdgCodesDaughters[0]= 22;
      pdgCodesDaughters[1]= 22;
      masses              = new Double_t[nDaughters];
      masses[0]           = TDatabasePDG::Instance()->GetParticle(pdgCodesDaughters[0])->Mass();;
      masses[1]           = TDatabasePDG::Instance()->GetParticle(pdgCodesDaughters[1])->Mass();;
   } else {
      cout << "particle not defined, aborting" << endl;
   }

   Color_t colorSmear[5]               = {kAzure, kRed+2, kGreen+2, kViolet+2, kCyan+2};
   Style_t markerSmear[5]              = {21, 20, 25, 24, 33};
   Color_t colorScaled[3]              = {kRed-6, kGreen-6, kBlue-6};
   Style_t markerScaled[3]             = {25, 24, 27};

   Double_t scaleFactors[4][5] = {
      {0.00, 0.00, 0.0, 0.0, 0.0},
      {0.00, -0.2, 0.2, 0.0, 0.0},
      {0.20, -0.2, 0.0, 0.0, 0.0},
      {-0.2, 0.20, 0.0, 0.0, 0.0}
   };
   Double_t ptBinning[75] = {
      0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
      1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,
      2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9,
      3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8,
      5.0, 5.4, 5.8, 6.2, 6.6, 7.0, 7.5, 8.0, 8.5, 9.0,
      10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0, 25.0,
      30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0,
      80., 100.,125.,150.,175.,200
   };
   Int_t nBinsX = 74;

   //*************************************************************************************************
   //************************** Loading the resolutions **********************************************
   //*************************************************************************************************

   if (fileName.CompareTo("")!= 0 && cutSting.CompareTo("") != 0){
      TFile* resolutionFile           = new TFile(fileName.Data());
      TString autoDetectedMainDir     = AutoDetectMainTList(mode , resolutionFile);
      if (autoDetectedMainDir.CompareTo("") == 0){
         cout << "ERROR: trying to read file, which is incompatible with mode selected" << endl;;
         return;
      }

      //************************** Container Loading ********************************************************************
      TList* TopDir                   = (TList*)resolutionFile->Get(autoDetectedMainDir.Data());
      if(TopDir == NULL){
         cout<<"ERROR: TopDir not Found"<<endl;
         return;
      }
      TList *HistosGammaConversion    = (TList*)TopDir->FindObject(Form("Cut Number %s",cutSting.Data()));
      TList *TrueConversionContainer  = (TList*)HistosGammaConversion->FindObject(Form("%s True histograms",cutSting.Data()));

      //****************************** Resolution dPt vs Pt **********************************************
      TString nameResolHist[5];
      nameResolHist[0]                = Form("ESD_TruePrimary%s_MCPt_ResolPt",outputlabel.Data());
      labelResolHist[0]               = Form("%s", plotLabel.Data());
      if (mode == 10){
         nameResolHist[0]            = Form("ESD_TruePrimary%sPureMerged_MCPt_ResolPt",outputlabel.Data());
         nameResolHist[1]            = Form("ESD_TruePrimary%sMergedPartConv_MCPt_ResolPt",outputlabel.Data());
         nameResolHist[2]            = Form("ESD_TruePrimary%s1Gamma_MCPt_ResolPt",outputlabel.Data());
         nameResolHist[3]            = Form("ESD_TruePrimary%s1Electron_MCPt_ResolPt",outputlabel.Data());
         labelResolHist[0]               = Form("%s fully merged", plotLabel.Data());
         labelResolHist[1]               = Form("%s merged part conv.", plotLabel.Data());
         labelResolHist[2]               = Form("%s only 1 #gamma", plotLabel.Data());
         labelResolHist[3]               = Form("%s only 1 e^{#pm}", plotLabel.Data());
      }
      for (Int_t i = 0; (i < 5) && (i < nResolHist); i++){
         histoResolutionInput[i]     = (TH2F*)TrueConversionContainer->FindObject(nameResolHist[i].Data()); //
         if (histoResolutionInput[i]){
            histoResolutionInput[i]->Sumw2();
            haveResol               = kTRUE;
         } else {
            haveResol               = kFALSE;
         }
      }
   }

   if (fileNameFractions.CompareTo("") != 0 && cutSting.CompareTo("") != 0 && haveResol){
      TFile* fractionsFile   = new TFile(fileNameFractions.Data());
      TString nameFracHist[5];
      nameFracHist[0]                 = "";
      if (mode == 10){
         nameFracHist[0]             = "RatioMergedPure";
         nameFracHist[1]             = "RatioPi0MergedPartConv";
         nameFracHist[2]             = "RatioPi0MergedOneGamma";
         nameFracHist[3]             = "RatioPi0MergedOneElectron";
      }
      for (Int_t i = 0; (i < 5) && (i < nResolHist); i++){
         cout << "trying to load: " << nameFracHist[i].Data() << endl;
         histoFractionsInput[i]      = (TH1F*)fractionsFile->Get(nameFracHist[i].Data()); //
         if (histoFractionsInput[i]){
            histoFractionsInput[i]->Sumw2();
            histoFractionsInput[i]->Scale(1./100);
            haveFractions           = kTRUE;
            cout << "succeeded." << endl;
         } else {
            haveFractions           = kFALSE;
            cout << "failed." << endl;
         }
      }
   }


   h2_ptMothervsDaughter = new TH2D("h2_ptMothervsDaughter","", 500,0,maxPt,500,0,maxPt);  
   h2_ptMothervsDaughter->Sumw2();
   h2_asym_geom          = new TH2D("h2_asym_geom","", 500,-1,1, 500,0,maxPt);
   h2_ptMothervsDaughter->Sumw2();
   h2_OpenAngle          = new TH2D("h2_OpenAngle","", 500,0,maxPt, 400,0,4);
   h2_OpenAngle->SetXTitle(Form("#it{p}_{T,%s} (GeV/#it{c})",plotLabel.Data()));
   h2_OpenAngle->SetYTitle("#theta_{open}");
   h2_OpenAngle->Sumw2();
   h2_ClusterDistance    = new TH2D("h2_ClusterDistance","", 500,0,maxPt, 5000,0,100);
   h2_ClusterDistance->SetXTitle(Form("#it{p}_{T,%s} (GeV/#it{c})",plotLabel.Data()));
   h2_ClusterDistance->SetYTitle("#it{R}_{#gamma's at r=440cm} (cm)");
   h2_ClusterDistance->Sumw2();

   h1_ptdistribution     = new TH1D("h1_ptdistribution","", 500,0,maxPt);
   h1_ptdistribution->Sumw2();
   h1_ptdistributionReb  = new TH1D("h1_ptdistributionReb","", nBinsX, ptBinning);
   h1_ptdistributionReb->Sumw2();


   TSpline3* spFractionsDummy[4];
   for (Int_t i = 0; i < nResolHist; i++){
      h1_ptdistSmeared[i]         = new TH1D(Form("h1_ptdistSmeared_%d",i),"", 500,0,maxPt);
      h1_ptdistSmeared[i]->Sumw2();
      h1_ptdistSmearedReb[i]      = new TH1D(Form("h1_ptdistSmearedReb_%d",i),"", nBinsX, ptBinning);
      h1_ptdistSmearedReb[i]->Sumw2();
      p1_fractionsSmeared[i]      = new TProfile(Form("p1_fractionsSmeared_%d",i),"", 500,0,maxPt);
      p1_fractionsSmeared[i]->Sumw2();
      p1_fractionsSmearedReb[i]   = new TProfile(Form("p1_fractionsSmearedReb_%d",i),"", nBinsX, ptBinning);
      p1_fractionsSmearedReb[i]->Sumw2();
      if (haveFractions){
         spFractionsDummy[i]              = new TSpline3(histoFractionsInput[i]);
         for (Int_t tt = 0; tt < 10; tt++){ // one entry per thread for safety
            spFractions[tt][i]              =(TSpline3*)spFractionsDummy[i]->Clone(Form("spFractionsThreaded_%d_%d",tt,i));
         }
      }
   }

   h1_phidistribution    = new TH1D("h1_phidistribution","", 100,0,2*TMath::Pi());
   h1_phidistribution->Sumw2();
   h1_etadistribution    = new TH1D("h1_etadistribution","", 200,-1,1);
   h1_etadistribution->Sumw2();
   for (Int_t i = 0; i < nDaughters; i++){
      h1_ptDaughter[i]        = new TH1D(Form("h1_ptDaughter%d",i),"", 500,0,maxPt);
   }
   h1_ptGammaDaughters   = new TH1D("h1_ptGammaDaughters","", 500,0,maxPt);
   h1_ptPiPlDaughters    = new TH1D("h1_ptPiPlDaughters","", 500,0,maxPt);
   h1_ptPiMiDaughters    = new TH1D("h1_ptPiMiDaughters","", 500,0,maxPt);
   h1_ptPiZeroDaughters  = new TH1D("h1_ptPiZeroDaughters","", 500,0,maxPt);

   if (haveResol){
      for (Int_t j=0; (j < nResolHist && j < 5); j++){
         for (Int_t i=0; (i < histoResolutionInput[j]->GetNbinsX()+1 && i < 1000); i++){
            for (Int_t tt = 0; tt < 10; tt++){ // one entry per thread for safety
               h1ResolProjection[tt][j][i] = (TH1D*)histoResolutionInput[j]->ProjectionY(Form("dummy_%d_%d_%d",j,i,tt), i+1,i+1,"e");
            }
         }
      }
   }



   //*************************************************************************************************
   //******************************** Thread creation and running ************************************
   //*************************************************************************************************

   // create and run individual threads
   for (Int_t cthr=1; cthr < numThreads+1; cthr++){
      printf("Starting Thread %d\n",cthr);
      hx[cthr] = new TThread(Form("h%d",cthr), threadhandle, (void*) cthr);
      hx[cthr]->Run();
   }

   // print info of threads and join them
   TThread::Ps();
   for (Int_t cthr=1; cthr < numThreads+1; cthr++){
      hx[cthr]->Join();
      TThread::Ps();
   }

   //*************************************************************************************************
   //******************************** Plotting with thread output ************************************
   //*************************************************************************************************




    //****************************** Plot input pT and decay distributions ***********************************************************************
    TCanvas *canvasQA = new TCanvas("canvasQA","canvasQA",1000,800);
    DrawGammaCanvasSettings( canvasQA, 0.07, 0.02, 0.02, 0.08);
    canvasQA->cd();
    canvasQA->SetLogy(1);

    TLegend* legendSpectra = GetAndSetLegend2(0.86, 0.70, 0.95, 0.93, 32,1); 
    DrawAutoGammaMesonHistos(   h1_ptdistribution, 
                                "", "#it{p}_{T} (GeV/#it{c})", "#it{N}_{X}", // (%)", 
                                kTRUE, 10, 1e-1, kFALSE,
                                kFALSE, 0., 0.7, 
                                kFALSE, 0., 10.);
    h1_ptdistribution->GetYaxis()->SetTitleOffset(0.85);
    DrawGammaSetMarker(h1_ptdistribution, 20, 1.5, kBlack, kBlack);
    h1_ptdistribution->DrawClone("pe");
    legendSpectra->AddEntry(h1_ptdistribution,plotLabel.Data(),"pe");

    for (Int_t i = 0; i < nDaughters; i++){
      DrawGammaSetMarker(h1_ptDaughter[i], 21+i, 1., colors[i], colors[i]);
      h1_ptDaughter[i]->Draw("same,pe");
      legendSpectra->AddEntry(h1_ptDaughter[i],daughterLabels[i].Data(),"pe");
    }
    if (h1_ptGammaDaughters->GetEntries() > 0){
      DrawGammaSetMarker(h1_ptGammaDaughters, 24, 1.5, kGray+2, kGray+2);
      h1_ptGammaDaughters->Draw("same,pe");
      legendSpectra->AddEntry(h1_ptGammaDaughters,"#gamma's","pe");
    }

    legendSpectra->Draw();
    canvasQA->SaveAs(Form("%s%s_PtDistribution_Input.%s",fOutputDir.Data(), outputlabel.Data(), suffix.Data()));

    //****************************** Plot smeared and input pT ***********************************************************************
    canvasQA->cd();
    TLegend* legendSpectra2 = GetAndSetLegend2(0.45, 0.93-(nResolHist+1)*0.035, 0.65, 0.93, 32,1); 
    DrawGammaSetMarker(h1_ptdistribution, 20, 1.5, kBlack, kBlack);
    h1_ptdistribution->DrawClone("pe");
    legendSpectra2->AddEntry(h1_ptdistribution,plotLabel.Data(),"pe");
    for (Int_t j = 0; (j < 5 && j < nResolHist); j++ ){
        DrawGammaSetMarker(h1_ptdistSmeared[j], markerSmear[j], 1.5, colorSmear[j], colorSmear[j]);
        h1_ptdistSmeared[j]->Draw("same,pe");
        legendSpectra2->AddEntry(h1_ptdistSmeared[j],Form("%s smeared",labelResolHist[j].Data()),"pe");
    }
    legendSpectra2->Draw();

    canvasQA->SaveAs(Form("%s%s_PtDistribution_InputVsSmeared%s.%s",fOutputDir.Data(), outputlabel.Data(), subMode.Data(), suffix.Data()));

    //**************************** Draw fractions as input ***************************************************************************
    if (haveFractions){
      canvasQA->cd();
      canvasQA->SetLogy(0);
      TH2F* dummyDrawingHist  = new TH2F("dummyDrawingHist","dummyDrawingHist",5000,0,maxPt,10000, 0., 0.8); 
      SetStyleHistoTH2ForGraphs(  dummyDrawingHist, "#it{p}_{T} (GeV/#it{c})", "L_{X} = clus. rec from X / all clus.", 0.028, 0.04, 
                                 0.028, 0.04, 0.86, 0.82, 510, 505);
      dummyDrawingHist->Draw();

      TLegend* legendFractions = GetAndSetLegend2(0.15, 0.93-nResolHist*0.035, 0.35, 0.93, 32,1); 
      for (Int_t j = 0; (j < 5 && j < nResolHist); j++ ){
         DrawGammaSetMarkerProfile(p1_fractionsSmearedReb[j], markerSmear[j], 1.5, colorSmear[j], colorSmear[j]);
         p1_fractionsSmearedReb[j]->DrawClone("same,pe");
         histoFractionsInput[j]->SetLineColor(colorSmear[j]);
         histoFractionsInput[j]->Draw("same,hist,c");

         legendFractions->AddEntry(p1_fractionsSmearedReb[j],labelResolHist[j].Data(),"pe");
      }
      legendFractions->Draw();

      DrawGammaLines(0, maxPt , 0.60, 0.60 ,1, kGray, 8);
      DrawGammaLines(0, maxPt , 0.40, 0.40 ,1, kGray, 8);
      DrawGammaLines(0, maxPt , 0.20, 0.20 ,1, kGray, 8);

      canvasQA->SaveAs(Form("%s%s_FractionsInputVsPt.%s",fOutputDir.Data(), outputlabel.Data(), suffix.Data()));
    }

    //**************************** Calculate and draw ratio of smeared and input dist ************************************************
    TH1D* histoRatioSmearedDivInput[nResolHist];

    canvasQA->cd();
    canvasQA->SetLogy(0);

    TLegend* legendRatio1 = GetAndSetLegend2(0.15, 0.93-nResolHist*0.035, 0.35, 0.93, 32,1); 
    for (Int_t j = 0; (j < 5 && j < nResolHist); j++ ){
      histoRatioSmearedDivInput[j]    = (TH1D*)h1_ptdistSmeared[j]->Clone(Form("histoRatioSmearedDivInput%d",j));
      histoRatioSmearedDivInput[j]->Divide(histoRatioSmearedDivInput[j],h1_ptdistribution,1,1,"B");
      if (j == 0){
         DrawAutoGammaMesonHistos(   histoRatioSmearedDivInput[j], 
                                       "", "#it{p}_{T} (GeV/#it{c})", "smeared/input", // (%)", 
                                       kFALSE, 10, 1e-1, kFALSE,
                                       kTRUE, 0., 5., 
                                       kFALSE, 0., 10.);
         histoRatioSmearedDivInput[j]->GetYaxis()->SetTitleOffset(0.85);
         DrawGammaSetMarker(histoRatioSmearedDivInput[j], markerSmear[j], 1.5, colorSmear[j], colorSmear[j]);
         histoRatioSmearedDivInput[j]->DrawClone("pe");

      }  else {
         DrawGammaSetMarker(histoRatioSmearedDivInput[j], markerSmear[j], 1.5, colorSmear[j], colorSmear[j]);
         histoRatioSmearedDivInput[j]->DrawClone("same,pe");
      }
      legendRatio1->AddEntry(histoRatioSmearedDivInput[j],labelResolHist[j].Data(),"pe");
    }
    legendRatio1->Draw();

    DrawGammaLines(0, maxPt , 1, 1 ,1, kGray, 7);
    DrawGammaLines(0, maxPt , 1.2, 1.2 ,1, kGray, 8);
    DrawGammaLines(0, maxPt , 0.8, 0.8 ,1, kGray, 8);

    canvasQA->SaveAs(Form("%s%s_Ratio_SmearedDivInputVsPt.%s",fOutputDir.Data(), outputlabel.Data(), suffix.Data()));

    //**************************** Calculate and draw ratio of smeared and input dist ************************************************
    TH1D* histoRatioSmearedDivInputReb[nResolHist];
    canvasQA->cd();
    canvasQA->SetLogy(0);

    for (Int_t j = 0; (j < 5 && j < nResolHist); j++ ){
      histoRatioSmearedDivInputReb[j]    = (TH1D*)h1_ptdistSmearedReb[j]->Clone(Form("histoRatioSmearedDivInputReb%d",j));
      histoRatioSmearedDivInputReb[j]->Divide(histoRatioSmearedDivInputReb[j],h1_ptdistributionReb,1,1,"B");
      if (j == 0){
         DrawAutoGammaMesonHistos(   histoRatioSmearedDivInputReb[j],
                                       "", "#it{p}_{T} (GeV/#it{c})", "smeared/input", // (%)",
                                       kFALSE, 10, 1e-1, kFALSE,
                                       kTRUE, 0., 5,
                                       kTRUE, 0., maxPt);
         histoRatioSmearedDivInputReb[j]->GetYaxis()->SetTitleOffset(0.85);
         DrawGammaSetMarker(histoRatioSmearedDivInputReb[j], markerSmear[j], 1.5, colorSmear[j], colorSmear[j]);
         histoRatioSmearedDivInputReb[j]->DrawClone("pe");
      }  else {
         DrawGammaSetMarker(histoRatioSmearedDivInputReb[j], markerSmear[j], 1.5, colorSmear[j], colorSmear[j]);
         histoRatioSmearedDivInputReb[j]->DrawClone("same,pe");
      }
    }
    legendRatio1->Draw();

    DrawGammaLines(0, maxPt , 1, 1 ,1, kGray, 7);
    DrawGammaLines(0, maxPt , 1.2, 1.2 ,1, kGray, 8);
    DrawGammaLines(0, maxPt , 0.8, 0.8 ,1, kGray, 8);

    canvasQA->SaveAs(Form("%s%s_Ratio_SmearedDivInputVsPtRebined.%s",fOutputDir.Data(), outputlabel.Data(), suffix.Data()));

    if (haveResol && haveFractions){
      //********************************************************************************************************************************
      //**************************** Calculate reconstructed spectrum ******************************************************************
      //********************************************************************************************************************************
      TH1D* h1_ptdistSmearedReweighted[4][nResolHist];
      TH1D* h1_ptdistSmearedReweightedReb[4][nResolHist];

      for (Int_t k = 0; k < 4; k++){
         for (Int_t j = 0; (j < 5 && j < nResolHist); j++ ){
            h1_ptdistSmearedReweighted[k][j] = (TH1D*)h1_ptdistSmeared[j]->Clone(Form("h1_ptdistSmeared_%d",j));
            h1_ptdistSmearedReweightedReb[k][j] = (TH1D*)h1_ptdistSmearedReb[j]->Clone(Form("h1_ptdistSmearedReb_%d",j));
            for (Int_t i = 1; i < h1_ptdistSmearedReweighted[k][j]->GetNbinsX()+1; i++){
               if (h1_ptdistSmearedReweighted[k][j]->GetBinCenter(i) > 10 && h1_ptdistSmearedReweighted[k][j]->GetBinCenter(i) < maxPt){
                  h1_ptdistSmearedReweighted[k][j]->SetBinContent(i, h1_ptdistSmearedReweighted[k][j]->GetBinContent(i)*(p1_fractionsSmeared[j]->GetBinContent(i)+scaleFactors[k][j])/
                                                                     h1_ptdistSmearedReweighted[k][j]->GetBinWidth(i));
                  h1_ptdistSmearedReweighted[k][j]->SetBinError(i, h1_ptdistSmearedReweighted[k][j]->GetBinError(i)*(p1_fractionsSmeared[j]->GetBinContent(i)+scaleFactors[k][j])/
                                                                     h1_ptdistSmearedReweighted[k][j]->GetBinWidth(i));
               } else {
                  h1_ptdistSmearedReweighted[k][j]->SetBinContent(i, 0);
                  h1_ptdistSmearedReweighted[k][j]->SetBinError(i, 0);
               }
            }
            for (Int_t i = 1; i < h1_ptdistSmearedReweightedReb[k][j]->GetNbinsX()+1; i++){
               if (h1_ptdistSmearedReweightedReb[k][j]->GetBinCenter(i) > 10 && h1_ptdistSmearedReweightedReb[k][j]->GetBinCenter(i) < maxPt){
                  h1_ptdistSmearedReweightedReb[k][j]->SetBinContent(i, h1_ptdistSmearedReweightedReb[k][j]->GetBinContent(i)*(p1_fractionsSmearedReb[j]->GetBinContent(i)+scaleFactors[k][j])/
                                                                        h1_ptdistSmearedReweightedReb[k][j]->GetBinWidth(i));
                  h1_ptdistSmearedReweightedReb[k][j]->SetBinError(i, h1_ptdistSmearedReweightedReb[k][j]->GetBinError(i)*(p1_fractionsSmearedReb[j]->GetBinContent(i)+scaleFactors[k][j])/
                                                                        h1_ptdistSmearedReweightedReb[k][j]->GetBinWidth(i));
               } else {
                  h1_ptdistSmearedReweightedReb[k][j]->SetBinContent(i, 0);
                  h1_ptdistSmearedReweightedReb[k][j]->SetBinError(i, 0);
               }
            }
         }
      }

      for (Int_t i = 1; i < h1_ptdistributionReb->GetNbinsX()+1; i++){
         if (h1_ptdistributionReb->GetBinCenter(i) > 10 && h1_ptdistributionReb->GetBinCenter(i) < maxPt){
            h1_ptdistributionReb->SetBinContent(i, h1_ptdistributionReb->GetBinContent(i)/h1_ptdistributionReb->GetBinWidth(i));
            h1_ptdistributionReb->SetBinError(i, h1_ptdistributionReb->GetBinError(i)/h1_ptdistributionReb->GetBinWidth(i));
         } else {
            h1_ptdistributionReb->SetBinContent(i, 0);
            h1_ptdistributionReb->SetBinError(i, 0);
         }
      }

      TH1D* h1_ptReweightedRecMCReb       = (TH1D*)h1_ptdistSmearedReweightedReb[0][0]->Clone("h1_ptReweightedRecMCReb");
      h1_ptReweightedRecMCReb->Sumw2();
      for (Int_t j = 1; (j < 5 && j < nResolHist); j++ ){
         h1_ptReweightedRecMCReb->Add(h1_ptdistSmearedReweightedReb[0][j]);
      }

      for (Int_t i = 0; i< 3; i++){
         h1_ptReweightedRecMCRebModFrac[i]      = new TH1D(Form("h1_ptReweightedRecMCRebModFrac_%d",i),"", nBinsX, ptBinning);
         h1_ptReweightedRecMCRebModFrac[i]->Sumw2();
         for (Int_t j = 0; (j < 5 && j < nResolHist); j++ ){
            h1_ptReweightedRecMCRebModFrac[i]->Add(h1_ptdistSmearedReweightedReb[i+1][j]);
         }
         h1_ptReweightedRecMCRebModFracRatioToStandard[i]      = (TH1D*)h1_ptReweightedRecMCRebModFrac[i]->Clone(Form("h1_ptReweightedRecMCRebModFracRatioToStandard_%d",i));
         h1_ptReweightedRecMCRebModFracRatioToStandard[i]->Divide(h1_ptReweightedRecMCRebModFracRatioToStandard[i],h1_ptReweightedRecMCReb,1,1,"B");
      }

      //****************************** Plot smeared and input pT reweighted************************************************************
      canvasQA->cd();
      canvasQA->SetLogy(1);
      DrawAutoGammaMesonHistos(   h1_ptdistributionReb,
                                 "", "#it{p}_{T} (GeV/#it{c})", "#it{N}_{X}",
                                 kTRUE, 10, 1e0, kFALSE,
                                 kFALSE, 0., 0.7,
                                 kFALSE, 0., 10.);
      h1_ptdistributionReb->GetYaxis()->SetTitleOffset(0.85);
      DrawGammaSetMarker(h1_ptdistributionReb, 20, 1.5, kBlack, kBlack);
      h1_ptdistributionReb->DrawClone("pe");
      for (Int_t j = 0; (j < 5 && j < nResolHist); j++ ){
         DrawGammaSetMarker(h1_ptdistSmearedReweightedReb[0][j], markerSmear[j], 1.5, colorSmear[j], colorSmear[j]);
         h1_ptdistSmearedReweightedReb[0][j]->Draw("same,pe");
      }
      legendSpectra2->Draw();
      canvasQA->SaveAs(Form("%s%s_PtDistribution_InputVsSmearedReweighted%s.%s",fOutputDir.Data(), outputlabel.Data(), subMode.Data(), suffix.Data()));

      //****************************** Plot reconstructed spectrum ********************************************************************
      canvasQA->cd();
      canvasQA->SetLogy(1);

      DrawGammaSetMarker(h1_ptdistributionReb, 20, 1.5, kBlack, kBlack);
      h1_ptdistributionReb->DrawClone("pe");

      DrawGammaSetMarker(h1_ptReweightedRecMCReb, 20, 1.5, kGray+2, kGray+2);
      h1_ptReweightedRecMCReb->DrawClone("same,pe");

      TLegend* legendSpectra3 = GetAndSetLegend2(0.55, 0.93-(2+3)*0.035, 0.75, 0.93, 32,1); 
      legendSpectra3->AddEntry(h1_ptdistributionReb,plotLabel.Data(),"pe");
      legendSpectra3->AddEntry(h1_ptReweightedRecMCReb,Form("rec. %s",plotLabel.Data()),"pe");

      for (Int_t i = 0; i<3; i++){
         DrawGammaSetMarker(h1_ptReweightedRecMCRebModFrac[i], markerScaled[i], 1.5, colorScaled[i], colorScaled[i]);
         h1_ptReweightedRecMCRebModFrac[i]->Draw("same,pe");
         legendSpectra3->AddEntry(h1_ptReweightedRecMCRebModFrac[i],Form("rec. %s, var %d",plotLabel.Data(),i+1),"pe");
      }
      legendSpectra3->Draw();
      canvasQA->SaveAs(Form("%s%s_PtDistribution_RecSpec%s.%s",fOutputDir.Data(), outputlabel.Data(), subMode.Data(), suffix.Data()));

      //**************************** Draw ratio of modified reweighted to MC reweighted ************************************************
      canvasQA->cd();
      canvasQA->SetLogy(0);

      TLegend* legendRatio2 = GetAndSetLegend2(0.55, 0.93-(2+3)*0.035, 0.75, 0.93, 32,1); 
      for (Int_t j = 0; j < 3 ; j++ ){
         if (j == 0){
            DrawAutoGammaMesonHistos(   h1_ptReweightedRecMCRebModFracRatioToStandard[j], 
                                       "", "#it{p}_{T} (GeV/#it{c})", "mod/standard", // (%)", 
                                       kFALSE, 10, 1e-1, kFALSE,
                                       kTRUE, 0.6, 1.4, 
                                       kTRUE, 0., maxPt);
            h1_ptReweightedRecMCRebModFracRatioToStandard[j]->GetYaxis()->SetTitleOffset(0.85);
            DrawGammaSetMarker(h1_ptReweightedRecMCRebModFracRatioToStandard[j], markerScaled[j], 2.5, colorScaled[j], colorScaled[j]);
            h1_ptReweightedRecMCRebModFracRatioToStandard[j]->DrawClone("pe");
         }  else {
            DrawGammaSetMarker(h1_ptReweightedRecMCRebModFracRatioToStandard[j], markerScaled[j], 2.5, colorScaled[j], colorScaled[j]);
            h1_ptReweightedRecMCRebModFracRatioToStandard[j]->DrawClone("same,pe");
         }
         legendRatio2->AddEntry(h1_ptReweightedRecMCRebModFracRatioToStandard[j],Form("rec. %s, var %d",plotLabel.Data(),j),"pe");
      }
      legendRatio2->Draw();

      DrawGammaLines(0, maxPt , 1, 1 ,1, kGray, 7);
      DrawGammaLines(0, maxPt , 1.1, 1.1 ,1, kGray, 8);
      DrawGammaLines(0, maxPt , 0.9, 0.9 ,1, kGray, 8);

      canvasQA->SaveAs(Form("%s%s_Ratio_ReweightedModDivStandard.%s",fOutputDir.Data(), outputlabel.Data(), suffix.Data()));
   }

   //******************************** Plot input phi phi distribution ***************************************************************
   canvasQA->SetLogy(1);
   DrawAutoGammaMesonHistos(   h1_phidistribution,
                              "", "#varphi (rad)", Form("#it{N}_{%s}",plotLabel.Data()),
                              kTRUE, 10, 1e-1, kFALSE,
                              kFALSE, 0., 0.7,
                              kFALSE, 0., 10.);
   h1_phidistribution->GetYaxis()->SetTitleOffset(0.85);
   h1_phidistribution->DrawClone();
   canvasQA->SaveAs(Form("%s_PhiDistribution_Input.%s",outputlabel.Data(),suffix.Data()));

   DrawAutoGammaMesonHistos(   h1_etadistribution,
                              "", "#eta", Form("#it{N}_{%s}",plotLabel.Data()),
                              kTRUE, 10, 1e-1, kFALSE,
                              kFALSE, 0., 0.7,
                              kFALSE, 0., 10.);
   h1_etadistribution->GetYaxis()->SetTitleOffset(0.85);
   h1_etadistribution->DrawClone();
   canvasQA->SaveAs(Form("%s%s_EtaDistribution_Input.%s",fOutputDir.Data(), outputlabel.Data(),suffix.Data()));

   TCanvas *canvasQA2D = new TCanvas("canvasQA2D","canvasQA2D",1000,800);
   DrawGammaCanvasSettings( canvasQA2D, 0.1, 0.09, 0.02, 0.11);
   canvasQA2D->cd();
   canvasQA2D->SetLogy(0);
   canvasQA2D->SetLogz(1);
   DrawAutoGammaHisto(h2_ptMothervsDaughter,
                              "",
                              Form("#it{p}_{T,%s} (GeV/#it{c})",plotLabel.Data()),
                              Form("#it{p}_{T,%s} (GeV/#it{c})",daughterLabels[0].Data()),
                              0,0,0,
                              0,0,5,
                              0,0, 50,0.85);

   h2_ptMothervsDaughter->Draw("colz");
   canvasQA2D->SaveAs(Form("%s%s_PtMothervsDaughter1.%s",fOutputDir.Data(), outputlabel.Data(),suffix.Data()));
   DrawGammaCanvasSettings( canvasQA2D, 0.1, 0.09, 0.02, 0.11);
   DrawAutoGammaHisto(h2_asym_geom,
                              "",
                              "#it{A} = (#it{p}_{T,1}-#it{p}_{T,2})/(#it{p}_{T,1}+#it{p}_{T,2})",
                              Form("#it{p}_{T,%s} (GeV/#it{c})",plotLabel.Data()),
                              0,0,0,
                              0,0,5,
                              0,0, 50,0.85);
   h2_asym_geom->Draw("colz");
   canvasQA2D->SaveAs(Form("%s%s_AsymmetryDaughtersVsMotherPt.%s",fOutputDir.Data(),outputlabel.Data(),suffix.Data()));

   DrawGammaCanvasSettings( canvasQA2D, 0.1, 0.09, 0.02, 0.11);
   DrawAutoGammaHisto(h2_OpenAngle,
                              "",
                              Form("#it{p}_{T,%s} (GeV/#it{c})",plotLabel.Data()),
                              "#theta_{open}",
                              0,0,0,
                              0,0,5,
                              0,0, 50,0.85);
   h2_OpenAngle->Draw("colz");
   canvasQA2D->SaveAs(Form("%s%s_MotherPtVsOpeningAngleDaughters.%s",fOutputDir.Data(),outputlabel.Data(),suffix.Data()));

   DrawAutoGammaHisto(h2_ClusterDistance,
                              "",
                              Form("#it{p}_{T,%s} (GeV/#it{c})",plotLabel.Data()),
                              "#it{R}_{#gamma's at r=440cm} (cm)",
                              0,0,0,
                              0,0,5,
                              0,0, 50,0.85);
   h2_ClusterDistance->Draw("colz");
   canvasQA2D->SaveAs(Form("%s%s_MotherPtVsClusterDistance.%s",fOutputDir.Data(),outputlabel.Data(),suffix.Data()));

   if (haveResol){
      canvasQA2D->SetLeftMargin(0.115);
      canvasQA2D->SetRightMargin(0.11);
      for (Int_t j = 0; (j < 5 && j < nResolHist); j++ ){
         DrawAutoGammaHisto(  histoResolutionInput[j],
                                       "",
                                       "#it{p}_{T,MC} (GeV/#it{c})", 
                                       Form ("(#it{p}^{%s}_{T,rec} -#it{p}^{%s}_{T,MC})/#it{p}^{%s}_{T,MC}", plotLabel.Data(), plotLabel.Data(), plotLabel.Data()),
                                       0,0,0,
                                       0,0,5,
                                       0,0, 50,0.85);
         histoResolutionInput[j]->GetZaxis()->SetRangeUser(1e-7, histoResolutionInput[j]->GetMaximum());
         histoResolutionInput[j]->GetYaxis()->SetTitleOffset(1.05);
         histoResolutionInput[j]->Draw("colz");
         PutProcessLabelAndEnergyOnPlot(0.15, 0.25, 28, fCollisionSystem.Data(), "mEMC", labelResolHist[j].Data(), 63, 0.03);

         canvasQA2D->SaveAs(Form("%s%s_Resolutions%d.%s",fOutputDir.Data(),outputlabel.Data(),j,suffix.Data()));
      }
   }


   //*************************************************************************************************
   //********************** Write histograms to file *************************************************
   //*************************************************************************************************
   TFile* fileOutput = new TFile(Form("%s/ToyMCOutput_%s.root",fCollisionSystenWrite.Data(), outputlabel.Data()),"RECREATE");

   h1_ptdistribution->Write();
   h1_phidistribution->Write();
   h1_etadistribution->Write();
   for (Int_t i = 0; i < nDaughters; i++){
      h1_ptDaughter[i]->Write();
   }
   for (Int_t j = 0; j < 3 ; j++ ){
      if( h1_ptReweightedRecMCRebModFracRatioToStandard[j])h1_ptReweightedRecMCRebModFracRatioToStandard[j]->Write();
   }
   h2_ptMothervsDaughter->Write();
   h2_ClusterDistance->Write();
   h2_OpenAngle->Write();
   h2_asym_geom->Write();

   fileOutput->Write();
   fileOutput->Close();
}
