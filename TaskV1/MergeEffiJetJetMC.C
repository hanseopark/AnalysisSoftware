/****************************************************************************************************************************
 ******     provided by Gamma Conversion Group, PWGGA,                                                       *****
 ******     Friederike Bock, friederike.bock@cern.ch                                                    *****
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

void MergeEffiJetJetMC(  TString fCutSelection, 
                         TString mesonType, 
                         TString fSuffix, 
                         TString fEnergyFlag, 
                         TString nameFileCorrectionFileFull, 
                         TString nameFileCorrectionFileMBMC, 
                         TString nameFileCorrectionFileJetJetMC,
                         Double_t minPt = 6.
                      ){
    
    gSystem->Exec("cp "+nameFileCorrectionFileFull+" "+nameFileCorrectionFileMBMC);
    TString outputDir                       = Form("%s/%s/%s/ExtractSignal", fCutSelection.Data(), fEnergyFlag.Data(), fSuffix.Data());
    gSystem->Exec("mkdir "+outputDir);
   
    TFile* file                             = new TFile (nameFileCorrectionFileFull);
    TH1D *histoTrueEffiPt                   = (TH1D*)file->Get("TrueMesonEffiPt");
    TH1D *histoEffiPt                       = (TH1D*)file->Get("MesonEffiPt");
    TH1D *histoEffiLeftPt                   = (TH1D*)file->Get("MesonLeftEffiPt");
    TH1D *histoTrueYield                    = (TH1D*)file->Get("histoYieldTrueMeson");
    TH1D *histoTrueEffiNarrowPt             = (TH1D*)file->Get("TrueMesonNarrowEffiPt");
    TH1D *histoEffiNarrowPt                 = (TH1D*)file->Get("MesonNarrowEffiPt");
    TH1D *histoEffiLeftNarrowPt             = (TH1D*)file->Get("MesonLeftNarrowEffiPt");
    TH1D *histoTrueEffiWidePt               = (TH1D*)file->Get("TrueMesonWideEffiPt");
    TH1D *histoEffiWidePt                   = (TH1D*)file->Get("MesonWideEffiPt");
    TH1D *histoEffiLeftWidePt               = (TH1D*)file->Get("MesonLeftWideEffiPt");
    TH1D* histoTrueMassMeson                = (TH1D*)file->Get("histoTrueMassMeson");
    TH1D* histoMassMeson                    = (TH1D*)file->Get("histoMassMeson");
    TH1D* histoTrueFWHMMeson                = (TH1D*)file->Get("histoTrueFWHMMeson");
    TH1D* histoFWHMMeson                    = (TH1D*)file->Get("histoFWHMMeson");
    TH1D* histoAcceptance                   = (TH1D*)file->Get("fMCMesonAccepPt");
    TH1D* histoAcceptanceWOWeights          = (TH1D*)file->Get("fMCMesonAccepPtWOWeights");
    TH1D* histoTrueEffiPtWOWeights          = (TH1D*)file->Get("TrueMesonEffiPtUnweighted");
    TH1D* histoTrueEffiNarrowPtWOWeights    = (TH1D*)file->Get("TrueMesonNarrowEffiPtUnweighted");
    TH1D* histoTrueEffiWidePtWOWeights      = (TH1D*)file->Get("TrueMesonWideEffiPtUnweighted");
    TH1D* histoTrueSecFrac                  = (TH1D*)file->Get("TrueSecFrac");
    TH1D* histoTrueSecFracWide              = (TH1D*)file->Get("TrueSecFracWide");
    TH1D* histoTrueSecFracNarrow            = (TH1D*)file->Get("TrueSecFracNarrow");
    TH1D* histoTrueSecFracFromK0S           = (TH1D*)file->Get("TrueSecFracFromK0S");
    TH1D* histoTrueSecFracFromLambda        = (TH1D*)file->Get("TrueSecFracFromLambda");
    TH1D* histoTrueSecFracFromK0SWide       = (TH1D*)file->Get("TrueSecFracFromK0SWide");
    TH1D* histoTrueSecFracFromK0SNarrow     = (TH1D*)file->Get("TrueSecFracFromK0SNarrow");
    
    
    TFile* fileMinBias                      = new TFile (nameFileCorrectionFileMBMC);
    TH1D *histoTrueEffiPtMinBias            = (TH1D*)fileMinBias->Get("TrueMesonEffiPt");
    TH1D *histoEffiPtMinBias                = (TH1D*)fileMinBias->Get("MesonEffiPt");
    TH1D *histoEffiLeftPtMinBias            = (TH1D*)fileMinBias->Get("MesonLeftEffiPt");
    TH1D *histoTrueYieldMinBias             = (TH1D*)fileMinBias->Get("histoYieldTrueMeson");
    TH1D *histoTrueEffiNarrowPtMinBias      = (TH1D*)fileMinBias->Get("TrueMesonNarrowEffiPt");
    TH1D *histoEffiNarrowPtMinBias          = (TH1D*)fileMinBias->Get("MesonNarrowEffiPt");
    TH1D *histoEffiLeftNarrowPtMinBias      = (TH1D*)fileMinBias->Get("MesonLeftNarrowEffiPt");
    TH1D *histoTrueEffiWidePtMinBias        = (TH1D*)fileMinBias->Get("TrueMesonWideEffiPt");
    TH1D *histoEffiWidePtMinBias            = (TH1D*)fileMinBias->Get("MesonWideEffiPt");
    TH1D *histoEffiLeftWidePtMinBias        = (TH1D*)fileMinBias->Get("MesonLeftWideEffiPt");
    TH1D* histoTrueMassMesonMinBias         = (TH1D*)fileMinBias->Get("histoTrueMassMeson");
    TH1D* histoMassMesonMinBias             = (TH1D*)fileMinBias->Get("histoMassMeson");
    TH1D* histoTrueFWHMMesonMinBias         = (TH1D*)fileMinBias->Get("histoTrueFWHMMeson");
    TH1D* histoFWHMMesonMinBias             = (TH1D*)fileMinBias->Get("histoFWHMMeson");
    TH1D* histoAcceptanceMinBias            = (TH1D*)fileMinBias->Get("fMCMesonAccepPt");
    TH1D* histoAcceptanceWOWeightsMinBias   = (TH1D*)fileMinBias->Get("fMCMesonAccepPtWOWeights");
    TH1D* histoTrueEffiPtWOWeightsMinBias   = (TH1D*)fileMinBias->Get("TrueMesonEffiPtUnweighted");
    TH1D* histoTrueEffiNarrowPtWOWeightsMinBias= (TH1D*)fileMinBias->Get("TrueMesonNarrowEffiPtUnweighted");
    TH1D* histoTrueEffiWidePtWOWeightsMinBias  = (TH1D*)fileMinBias->Get("TrueMesonWideEffiPtUnweighted");
    TH1D* histoTrueSecFracMinBias              = (TH1D*)fileMinBias->Get("TrueSecFrac");
    TH1D* histoTrueSecFracWideMinBias          = (TH1D*)fileMinBias->Get("TrueSecFracWide");
    TH1D* histoTrueSecFracNarrowMinBias        = (TH1D*)fileMinBias->Get("TrueSecFracNarrow");
    TH1D* histoTrueSecFracFromK0SMinBias       = (TH1D*)fileMinBias->Get("TrueSecFracFromK0S");
    TH1D* histoTrueSecFracFromLambdaMinBias    = (TH1D*)fileMinBias->Get("TrueSecFracFromLambda");
    TH1D* histoTrueSecFracFromK0SWideMinBias   = (TH1D*)fileMinBias->Get("TrueSecFracFromK0SWide");
    TH1D* histoTrueSecFracFromK0SNarrowMinBias = (TH1D*)fileMinBias->Get("TrueSecFracFromK0SNarrow");
        
    TFile* fileJetJetMC                     = new TFile (nameFileCorrectionFileJetJetMC);
    TH1F *histoEventQualityJetJetMC         = NULL;
    TH1D *histoMCInputJetJetMC              = NULL;
    TH1D *histoMCInputWOWeightingJetJetMC   = NULL;
    TH1D *histoMCInputWeightsJetJetMC       = NULL;
    histoEventQualityJetJetMC               = (TH1F*)fileJetJetMC->Get("NEvents");
    histoMCInputJetJetMC                    = (TH1D*)fileJetJetMC->Get("MC_Meson_genPt_oldBin");
    histoMCInputWOWeightingJetJetMC         = (TH1D*)fileJetJetMC->Get("MC_Meson_genPt_WOWeights");
    histoMCInputWeightsJetJetMC             = (TH1D*)fileJetJetMC->Get("MC_Meson_genPt_Weights");
    TH1D *histoTrueEffiPtJetJetMC           = (TH1D*)fileJetJetMC->Get("TrueMesonEffiPt");
    TH1D *histoEffiPtJetJetMC               = (TH1D*)fileJetJetMC->Get("MesonEffiPt");
    TH1D *histoEffiLeftPtJetJetMC           = (TH1D*)fileJetJetMC->Get("MesonLeftEffiPt");
    TH1D *histoTrueYieldJetJetMC            = (TH1D*)fileJetJetMC->Get("histoYieldTrueMeson");
    TH1D *histoTrueEffiNarrowPtJetJetMC     = (TH1D*)fileJetJetMC->Get("TrueMesonNarrowEffiPt");
    TH1D *histoEffiNarrowPtJetJetMC         = (TH1D*)fileJetJetMC->Get("MesonNarrowEffiPt");
    TH1D *histoEffiLeftNarrowPtJetJetMC     = (TH1D*)fileJetJetMC->Get("MesonLeftNarrowEffiPt");
    TH1D *histoTrueEffiWidePtJetJetMC       = (TH1D*)fileJetJetMC->Get("TrueMesonWideEffiPt");
    TH1D *histoEffiWidePtJetJetMC           = (TH1D*)fileJetJetMC->Get("MesonWideEffiPt");
    TH1D *histoEffiLeftWidePtJetJetMC       = (TH1D*)fileJetJetMC->Get("MesonLeftWideEffiPt");
    TH1D* histoTrueMassMesonJetJetMC        = (TH1D*)fileJetJetMC->Get("histoTrueMassMeson");
    TH1D* histoMassMesonJetJetMC            = (TH1D*)fileJetJetMC->Get("histoMassMeson");
    TH1D* histoTrueFWHMMesonJetJetMC        = (TH1D*)fileJetJetMC->Get("histoTrueFWHMMeson");
    TH1D* histoFWHMMesonJetJetMC            = (TH1D*)fileJetJetMC->Get("histoFWHMMeson");
    TH1D* histoAcceptanceJetJetMC           = (TH1D*)fileJetJetMC->Get("fMCMesonAccepPt");
    TH1D* histoAcceptanceWOWeightsJetJetMC  = (TH1D*)fileJetJetMC->Get("fMCMesonAccepPtWOWeights");
    TH1D* histoTrueEffiPtWOWeightsJetJetMC  = (TH1D*)fileJetJetMC->Get("TrueMesonEffiPtUnweighted");
    TH1D* histoTrueEffiNarrowPtWOWeightsJetJetMC= (TH1D*)fileJetJetMC->Get("TrueMesonNarrowEffiPtUnweighted");
    TH1D* histoTrueEffiWidePtWOWeightsJetJetMC  = (TH1D*)fileJetJetMC->Get("TrueMesonWideEffiPtUnweighted");
    TH1D* histoTrueSecFracJetJetMC              = (TH1D*)fileJetJetMC->Get("TrueSecFrac");
    TH1D* histoTrueSecFracWideJetJetMC          = (TH1D*)fileJetJetMC->Get("TrueSecFracWide");
    TH1D* histoTrueSecFracNarrowJetJetMC        = (TH1D*)fileJetJetMC->Get("TrueSecFracNarrow");
    TH1D* histoTrueSecFracFromK0SJetJetMC       = (TH1D*)fileJetJetMC->Get("TrueSecFracFromK0S");
    TH1D* histoTrueSecFracFromLambdaJetJetMC    = (TH1D*)fileJetJetMC->Get("TrueSecFracFromLambda");
    TH1D* histoTrueSecFracFromK0SWideJetJetMC   = (TH1D*)fileJetJetMC->Get("TrueSecFracFromK0SWide");
    TH1D* histoTrueSecFracFromK0SNarrowJetJetMC = (TH1D*)fileJetJetMC->Get("TrueSecFracFromK0SNarrow");
    
    TH1D *histoTrueEffiPtWeighted           = (TH1D*)histoTrueEffiPtJetJetMC->Clone("histoTrueEffiPtWeighted");
   
    for (Int_t i = 1;i < histoTrueEffiPtMinBias->GetNbinsX()+1 ;i++){
        if (histoTrueEffiPtMinBias->GetBinCenter(i) < minPt) continue;
        Double_t relErrMinBias      = histoTrueEffiPtMinBias->GetBinError(i)/histoTrueEffiPtMinBias->GetBinContent(i)*100;
        Double_t relErrJetJetMC     = histoTrueEffiPtJetJetMC->GetBinError(i)/histoTrueEffiPtJetJetMC->GetBinContent(i)*100;
                
        Double_t weightMinBias      = 1/TMath::Power(histoTrueEffiPtMinBias->GetBinError(i),2);
        Double_t weightJetJetMC     = 1/TMath::Power(histoTrueEffiPtJetJetMC->GetBinError(i),2);
        Double_t weightSum          = weightMinBias + weightJetJetMC;
        Double_t weightedEffi       = (weightMinBias*histoTrueEffiPtMinBias->GetBinContent(i) + weightJetJetMC * histoTrueEffiPtJetJetMC->GetBinContent(i))/weightSum;
        Double_t weightedEffiErr    = pow((weightMinBias +  weightJetJetMC),-0.5);
        
        if (isfinite(weightedEffi) && isfinite(weightedEffiErr)){
            histoTrueEffiPt->SetBinContent(i, weightedEffi);
            histoTrueEffiPt->SetBinError(i, weightedEffiErr);
        }
        cout << histoTrueEffiPtMinBias->GetBinContent(i) << "\t" << histoTrueEffiPtJetJetMC->GetBinContent(i)<<"\t" << histoTrueEffiPtWeighted->GetBinContent(i) << "\t" << histoTrueEffiPtWeighted->GetBinError(i) << endl;
        cout << "NORMAL: pt: "<< histoTrueEffiPtMinBias->GetBinCenter(i) <<"\t error min Bias: " << relErrMinBias << "\t error added sig: "  << relErrJetJetMC << endl;
        
        if (histoTrueEffiPtWOWeightsMinBias){
            cout << histoTrueEffiPtWOWeightsMinBias->GetBinError(i) << endl;
            cout << histoTrueEffiPtWOWeightsJetJetMC->GetBinError(i) << endl;
            Double_t weightMinBiasW0We      = 1/TMath::Power(histoTrueEffiPtWOWeightsMinBias->GetBinError(i),2);
            Double_t weightJetJetMCW0We     = 1/TMath::Power(histoTrueEffiPtWOWeightsJetJetMC->GetBinError(i),2);
            Double_t weightSumW0We          = weightMinBiasW0We + weightJetJetMCW0We;
            Double_t weightedEffiW0We       = (weightMinBiasW0We*histoTrueEffiPtWOWeights->GetBinContent(i) + weightJetJetMCW0We * histoTrueEffiPtWOWeightsJetJetMC->GetBinContent(i))/weightSumW0We;
            Double_t weightedEffiErrW0We    = pow((weightMinBiasW0We +  weightJetJetMCW0We),-0.5);
            
            if (isfinite(weightedEffiW0We) && isfinite(weightedEffiErrW0We)){
                histoTrueEffiPtWOWeights->SetBinContent(i, weightedEffiW0We);
                histoTrueEffiPtWOWeights->SetBinError(i, weightedEffiErrW0We);
            }
        }
        
        if (histoTrueEffiNarrowPtWOWeightsMinBias){
            Double_t weightMinBiasW0WeN      = 1/TMath::Power(histoTrueEffiNarrowPtWOWeightsMinBias->GetBinError(i),2);
            Double_t weightJetJetMCW0WeN     = 1/TMath::Power(histoTrueEffiNarrowPtWOWeightsJetJetMC->GetBinError(i),2);
            Double_t weightSumW0WeN          = weightMinBiasW0WeN + weightJetJetMCW0WeN;
            Double_t weightedEffiW0WeN       = (weightMinBiasW0WeN*histoTrueEffiNarrowPtWOWeightsMinBias->GetBinContent(i) 
                                                + weightJetJetMCW0WeN * histoTrueEffiNarrowPtWOWeightsJetJetMC->GetBinContent(i))/weightSumW0WeN;
            Double_t weightedEffiErrW0WeN    = pow((weightMinBiasW0WeN +  weightJetJetMCW0WeN),-0.5);
            
            if (isfinite(weightedEffiErrW0WeN) && isfinite(weightedEffiW0WeN)){
                histoTrueEffiNarrowPtWOWeights->SetBinContent(i, weightedEffiW0WeN);
                histoTrueEffiNarrowPtWOWeights->SetBinError(i, weightedEffiErrW0WeN);
            }
        }
        
        if (histoTrueEffiWidePtWOWeightsMinBias){
            Double_t weightMinBiasW0WeW      = 1/TMath::Power(histoTrueEffiWidePtWOWeightsMinBias->GetBinError(i),2);
            Double_t weightJetJetMCW0WeW     = 1/TMath::Power(histoTrueEffiWidePtWOWeightsJetJetMC->GetBinError(i),2);
            Double_t weightSumW0WeW          = weightMinBiasW0WeW + weightJetJetMCW0WeW;
            Double_t weightedEffiW0WeW       = (weightMinBiasW0WeW*histoTrueEffiWidePtWOWeightsMinBias->GetBinContent(i) 
                                                + weightJetJetMCW0WeW * histoTrueEffiWidePtWOWeightsJetJetMC->GetBinContent(i))/weightSumW0WeW;
            Double_t weightedEffiErrW0WeW    = pow((weightMinBiasW0WeW +  weightJetJetMCW0WeW),-0.5);
            
            if (isfinite(weightedEffiW0WeW) && isfinite(weightedEffiErrW0WeW)){
                histoTrueEffiWidePtWOWeights->SetBinContent(i, weightedEffiW0WeW);
                histoTrueEffiWidePtWOWeights->SetBinError(i, weightedEffiErrW0WeW);
            }
        }
        
        Double_t weightMinBiasAcc      = 1/TMath::Power(histoAcceptanceMinBias->GetBinError(i),2);
        Double_t weightJetJetMCAcc     = 1/TMath::Power(histoAcceptanceJetJetMC->GetBinError(i),2);
        Double_t weightSumAcc          = weightMinBiasAcc + weightJetJetMCAcc;
        Double_t weightedAcc       = (weightMinBiasAcc*histoAcceptanceMinBias->GetBinContent(i) + weightJetJetMCAcc * histoAcceptanceJetJetMC->GetBinContent(i))/weightSumAcc;
        Double_t weightedAccErr    = pow((weightMinBiasAcc +  weightJetJetMCAcc),-0.5);
        
        if (isfinite(weightedAcc) && isfinite(weightedAccErr)){
            histoAcceptance->SetBinContent(i, weightedAcc);
            histoAcceptance->SetBinError(i, weightedAccErr);
        }
        
        if (histoAcceptanceWOWeightsMinBias){
            Double_t weightMinBiasAccW0We      = 1/TMath::Power(histoAcceptanceWOWeightsMinBias->GetBinError(i),2);
            Double_t weightJetJetMCAccW0We     = 1/TMath::Power(histoAcceptanceWOWeightsJetJetMC->GetBinError(i),2);
            Double_t weightSumAccW0We          = weightMinBiasAccW0We + weightJetJetMCAccW0We;
            Double_t weightedAccW0We           = (weightMinBiasAccW0We*histoAcceptanceWOWeightsMinBias->GetBinContent(i) 
                                                + weightJetJetMCAccW0We * histoAcceptanceWOWeightsJetJetMC->GetBinContent(i))/weightSumAccW0We;
            Double_t weightedAccW0WeErr    = pow((weightMinBiasAccW0We +  weightJetJetMCAccW0We),-0.5);
            
            if (isfinite(weightedAccW0We) && isfinite(weightedAccW0WeErr)){
                histoAcceptanceWOWeights->SetBinContent(i, weightedAccW0We);
                histoAcceptanceWOWeights->SetBinError(i, weightedAccW0WeErr);
            }
        }
        
        Double_t weightMinBiasNEff      = 1/TMath::Power(histoEffiPtMinBias->GetBinError(i),2);
        Double_t weightJetJetMCNEff     = 1/TMath::Power(histoEffiPtJetJetMC->GetBinError(i),2);
        Double_t weightSumNEff          = weightMinBiasNEff + weightJetJetMCNEff;
        Double_t weightedNEff       = (weightMinBiasNEff*histoEffiPtMinBias->GetBinContent(i) + weightJetJetMCNEff * histoEffiPtJetJetMC->GetBinContent(i))/weightSumNEff;
        Double_t weightedNEffErr    = pow((weightMinBiasNEff +  weightJetJetMCNEff),-0.5);
        
        if (isfinite(weightedNEff) && isfinite(weightedNEffErr)){
            histoEffiPt->SetBinContent(i, weightedNEff);
            histoEffiPt->SetBinError(i, weightedNEffErr);
        }

        Double_t weightMinBiasNEffL      = 1/TMath::Power(histoEffiLeftPtMinBias->GetBinError(i),2);
        Double_t weightJetJetMCNEffL     = 1/TMath::Power(histoEffiLeftPtJetJetMC->GetBinError(i),2);
        Double_t weightSumNEffL          = weightMinBiasNEffL + weightJetJetMCNEffL;
        Double_t weightedNEffL       = (weightMinBiasNEffL*histoEffiLeftPtMinBias->GetBinContent(i) + weightJetJetMCNEffL * histoEffiLeftPtJetJetMC->GetBinContent(i))/weightSumNEffL;
        Double_t weightedNEffLErr    = pow((weightMinBiasNEffL +  weightJetJetMCNEffL),-0.5);
        
        if (isfinite(weightedNEffL) && isfinite(weightedNEffLErr)){
            histoEffiLeftPt->SetBinContent(i, weightedNEffL);
            histoEffiLeftPt->SetBinError(i, weightedNEffLErr);
        }
        
        
        if (relErrJetJetMC < relErrMinBias && isfinite(relErrJetJetMC) ){
//          histoTrueEffiPt->SetBinContent(i,histoTrueEffiPtJetJetMC->GetBinContent(i));
//          histoTrueEffiPt->SetBinError(i,histoTrueEffiPtJetJetMC->GetBinError(i));
          histoTrueMassMeson->SetBinContent(i,histoTrueMassMesonJetJetMC->GetBinContent(i));
          histoTrueMassMeson->SetBinError(i,histoTrueMassMesonJetJetMC->GetBinError(i));
          histoTrueFWHMMeson->SetBinContent(i,histoTrueFWHMMesonJetJetMC->GetBinContent(i));
          histoTrueFWHMMeson->SetBinError(i,histoTrueFWHMMesonJetJetMC->GetBinError(i));
          histoMassMeson->SetBinContent(i,histoMassMesonJetJetMC->GetBinContent(i));
          histoMassMeson->SetBinError(i,histoMassMesonJetJetMC->GetBinError(i));
          histoFWHMMeson->SetBinContent(i,histoFWHMMesonJetJetMC->GetBinContent(i));
          histoFWHMMeson->SetBinError(i,histoFWHMMesonJetJetMC->GetBinError(i));
          if (mesonType.Contains("Pi0")){
            histoTrueSecFrac->SetBinContent(i,histoTrueSecFracJetJetMC->GetBinContent(i));
            histoTrueSecFrac->SetBinError(i,histoTrueSecFracJetJetMC->GetBinError(i));
            histoTrueSecFracWide->SetBinContent(i,histoTrueSecFracWideJetJetMC->GetBinContent(i));
            histoTrueSecFracWide->SetBinError(i,histoTrueSecFracWideJetJetMC->GetBinError(i));
            histoTrueSecFracNarrow->SetBinContent(i,histoTrueSecFracNarrowJetJetMC->GetBinContent(i));
            histoTrueSecFracNarrow->SetBinError(i,histoTrueSecFracNarrowJetJetMC->GetBinError(i));
            histoTrueSecFracFromK0S->SetBinContent(i,histoTrueSecFracFromK0SJetJetMC->GetBinContent(i));
            histoTrueSecFracFromK0S->SetBinError(i,histoTrueSecFracFromK0SJetJetMC->GetBinError(i));
            histoTrueSecFracFromK0SWide->SetBinContent(i,histoTrueSecFracFromK0SWideJetJetMC->GetBinContent(i));
            histoTrueSecFracFromK0SWide->SetBinError(i,histoTrueSecFracFromK0SWideJetJetMC->GetBinError(i));
            histoTrueSecFracFromK0SNarrow->SetBinContent(i,histoTrueSecFracFromK0SNarrowJetJetMC->GetBinContent(i));
            histoTrueSecFracFromK0SNarrow->SetBinError(i,histoTrueSecFracFromK0SNarrowJetJetMC->GetBinError(i));
            histoTrueSecFracFromLambda->SetBinContent(i,histoTrueSecFracFromLambdaJetJetMC->GetBinContent(i));
            histoTrueSecFracFromLambda->SetBinError(i,histoTrueSecFracFromLambdaJetJetMC->GetBinError(i));
          }  
        } else {
//          histoTrueEffiPt->SetBinContent(i,histoTrueEffiPtMinBias->GetBinContent(i));
//          histoTrueEffiPt->SetBinError(i,histoTrueEffiPtMinBias->GetBinError(i));
          histoTrueMassMeson->SetBinContent(i,histoTrueMassMesonMinBias->GetBinContent(i));
          histoTrueMassMeson->SetBinError(i,histoTrueMassMesonMinBias->GetBinError(i));
          histoTrueFWHMMeson->SetBinContent(i,histoTrueFWHMMesonMinBias->GetBinContent(i));
          histoTrueFWHMMeson->SetBinError(i,histoTrueFWHMMesonMinBias->GetBinError(i));
          histoMassMeson->SetBinContent(i,histoMassMesonMinBias->GetBinContent(i));
          histoMassMeson->SetBinError(i,histoMassMesonMinBias->GetBinError(i));
          histoFWHMMeson->SetBinContent(i,histoFWHMMesonMinBias->GetBinContent(i));
          histoFWHMMeson->SetBinError(i,histoFWHMMesonMinBias->GetBinError(i));
          if (mesonType.Contains("Pi0")){
            histoTrueSecFrac->SetBinContent(i,histoTrueSecFracMinBias->GetBinContent(i));
            histoTrueSecFrac->SetBinError(i,histoTrueSecFracMinBias->GetBinError(i));
            histoTrueSecFracWide->SetBinContent(i,histoTrueSecFracWideMinBias->GetBinContent(i));
            histoTrueSecFracWide->SetBinError(i,histoTrueSecFracWideMinBias->GetBinError(i));
            histoTrueSecFracNarrow->SetBinContent(i,histoTrueSecFracNarrowMinBias->GetBinContent(i));
            histoTrueSecFracNarrow->SetBinError(i,histoTrueSecFracNarrowMinBias->GetBinError(i));
            histoTrueSecFracFromK0S->SetBinContent(i,histoTrueSecFracFromK0SMinBias->GetBinContent(i));
            histoTrueSecFracFromK0S->SetBinError(i,histoTrueSecFracFromK0SMinBias->GetBinError(i));
            histoTrueSecFracFromK0SWide->SetBinContent(i,histoTrueSecFracFromK0SWideMinBias->GetBinContent(i));
            histoTrueSecFracFromK0SWide->SetBinError(i,histoTrueSecFracFromK0SWideMinBias->GetBinError(i));
            histoTrueSecFracFromK0SNarrow->SetBinContent(i,histoTrueSecFracFromK0SNarrowMinBias->GetBinContent(i));
            histoTrueSecFracFromK0SNarrow->SetBinError(i,histoTrueSecFracFromK0SNarrowMinBias->GetBinError(i));
            histoTrueSecFracFromLambda->SetBinContent(i,histoTrueSecFracFromLambdaMinBias->GetBinContent(i));
            histoTrueSecFracFromLambda->SetBinError(i,histoTrueSecFracFromLambdaMinBias->GetBinError(i));
          }            
        }
        
        Double_t relErrMinBiasNarrow = histoTrueEffiNarrowPtMinBias->GetBinError(i)/histoTrueEffiNarrowPtMinBias->GetBinContent(i)*100;
        Double_t relErrJetJetMCNarrow = histoTrueEffiNarrowPtJetJetMC->GetBinError(i)/histoTrueEffiNarrowPtJetJetMC->GetBinContent(i)*100;
        cout << "NARROW: pt: "<< histoTrueEffiPtMinBias->GetBinCenter(i) <<"\t error min Bias: " << relErrMinBiasNarrow << "\t error added sig: "  << relErrJetJetMCNarrow << endl;
        
        Double_t weightMinBiasNarrow = 1/TMath::Power(histoTrueEffiNarrowPtMinBias->GetBinError(i),2);
        Double_t weightJetJetMCNarrow = 1/TMath::Power(histoTrueEffiNarrowPtJetJetMC->GetBinError(i),2);
        Double_t weightSumNarrow = weightMinBiasNarrow + weightJetJetMCNarrow;
        Double_t weightedEffiNarrow = (weightMinBiasNarrow*histoTrueEffiNarrowPtMinBias->GetBinContent(i) + weightJetJetMCNarrow * histoTrueEffiNarrowPtJetJetMC->GetBinContent(i))/weightSumNarrow;
        Double_t weightedEffiErrNarrow = pow((weightMinBiasNarrow +  weightJetJetMCNarrow),-0.5);

        if (isfinite(weightedEffiNarrow) && isfinite(weightedEffiErrNarrow)){
            histoTrueEffiNarrowPt->SetBinContent(i, weightedEffiNarrow);
            histoTrueEffiNarrowPt->SetBinError(i, weightedEffiErrNarrow);
        }

        Double_t weightMinBiasNEffN      = 1/TMath::Power(histoEffiNarrowPtMinBias->GetBinError(i),2);
        Double_t weightJetJetMCNEffN     = 1/TMath::Power(histoEffiNarrowPtJetJetMC->GetBinError(i),2);
        Double_t weightSumNEffN          = weightMinBiasNEffN + weightJetJetMCNEffN;
        Double_t weightedNEffN       = (weightMinBiasNEffN*histoEffiNarrowPtMinBias->GetBinContent(i) + weightJetJetMCNEffN * histoEffiNarrowPtJetJetMC->GetBinContent(i))/weightSumNEffN;
        Double_t weightedNEffNErr    = pow((weightMinBiasNEffN +  weightJetJetMCNEffN),-0.5);
        
        if (isfinite(weightedNEffN) && isfinite(weightedNEffNErr)){
            histoEffiNarrowPt->SetBinContent(i, weightedNEffN);
            histoEffiNarrowPt->SetBinError(i, weightedNEffNErr);
        }

        Double_t weightMinBiasNEffLN      = 1/TMath::Power(histoEffiLeftNarrowPtMinBias->GetBinError(i),2);
        Double_t weightJetJetMCNEffLN     = 1/TMath::Power(histoEffiLeftNarrowPtJetJetMC->GetBinError(i),2);
        Double_t weightSumNEffLN          = weightMinBiasNEffLN + weightJetJetMCNEffLN;
        Double_t weightedNEffLN       = (weightMinBiasNEffLN*histoEffiLeftNarrowPtMinBias->GetBinContent(i) + weightJetJetMCNEffLN * histoEffiLeftNarrowPtJetJetMC->GetBinContent(i))/weightSumNEffLN;
        Double_t weightedNEffLNErr    = pow((weightMinBiasNEffLN +  weightJetJetMCNEffLN),-0.5);
        
        if (isfinite(weightedNEffLN) && isfinite(weightedNEffLNErr)){
            histoEffiLeftNarrowPt->SetBinContent(i, weightedNEffLN);
            histoEffiLeftNarrowPt->SetBinError(i, weightedNEffLNErr);
        }
       
//      if ( relErrJetJetMCNarrow< relErrMinBiasNarrow && isfinite(relErrJetJetMCNarrow) ){
//          histoTrueEffiNarrowPt->SetBinContent(i,histoTrueEffiNarrowPtJetJetMC->GetBinContent(i));
//          histoTrueEffiNarrowPt->SetBinError(i,histoTrueEffiNarrowPtJetJetMC->GetBinError(i));
//      } else {
//          histoTrueEffiNarrowPt->SetBinContent(i,histoTrueEffiNarrowPtMinBias->GetBinContent(i));
//          histoTrueEffiNarrowPt->SetBinError(i,histoTrueEffiNarrowPtMinBias->GetBinError(i));
//      }

        Double_t relErrMinBiasWide = histoTrueEffiWidePtMinBias->GetBinError(i)/histoTrueEffiWidePtMinBias->GetBinContent(i)*100;
        Double_t relErrJetJetMCWide = histoTrueEffiWidePtJetJetMC->GetBinError(i)/histoTrueEffiWidePtJetJetMC->GetBinContent(i)*100;
        cout << "WIDE: : pt: "<< histoTrueEffiPtMinBias->GetBinCenter(i) <<"\t error min Bias: " <<  relErrMinBiasWide << "\t error added sig: "  << relErrJetJetMCWide << endl;
        
        Double_t weightMinBiasWide = 1/TMath::Power(histoTrueEffiWidePtMinBias->GetBinError(i),2);
        Double_t weightJetJetMCWide = 1/TMath::Power(histoTrueEffiWidePtJetJetMC->GetBinError(i),2);
        Double_t weightSumWide = weightMinBiasWide + weightJetJetMCWide;
        Double_t weightedEffiWide = (weightMinBiasWide*histoTrueEffiWidePtMinBias->GetBinContent(i) + weightJetJetMCWide * histoTrueEffiWidePtJetJetMC->GetBinContent(i))/weightSumWide;
        Double_t weightedEffiErrWide = pow((weightMinBiasWide +  weightJetJetMCWide),-0.5);

        if (isfinite(weightedEffiWide) && isfinite(weightedEffiErrWide)){
            histoTrueEffiWidePt->SetBinContent(i, weightedEffiWide);
            histoTrueEffiWidePt->SetBinError(i, weightedEffiErrWide);
        }

        Double_t weightMinBiasNEffW      = 1/TMath::Power(histoEffiWidePtMinBias->GetBinError(i),2);
        Double_t weightJetJetMCNEffW     = 1/TMath::Power(histoEffiWidePtJetJetMC->GetBinError(i),2);
        Double_t weightSumNEffW          = weightMinBiasNEffW + weightJetJetMCNEffW;
        Double_t weightedNEffW       = (weightMinBiasNEffW*histoEffiWidePtMinBias->GetBinContent(i) + weightJetJetMCNEffW * histoEffiWidePtJetJetMC->GetBinContent(i))/weightSumNEffW;
        Double_t weightedNEffWErr    = pow((weightMinBiasNEffW +  weightJetJetMCNEffW),-0.5);
        
        if (isfinite(weightedNEffW) && isfinite(weightedNEffWErr)){
            histoEffiWidePt->SetBinContent(i, weightedNEffW);
            histoEffiWidePt->SetBinError(i, weightedNEffWErr);
        }

        Double_t weightMinBiasNEffLW      = 1/TMath::Power(histoEffiLeftWidePtMinBias->GetBinError(i),2);
        Double_t weightJetJetMCNEffLW     = 1/TMath::Power(histoEffiLeftWidePtJetJetMC->GetBinError(i),2);
        Double_t weightSumNEffLW          = weightMinBiasNEffLW + weightJetJetMCNEffLW;
        Double_t weightedNEffLW       = (weightMinBiasNEffLW*histoEffiLeftWidePtMinBias->GetBinContent(i) + weightJetJetMCNEffLW * histoEffiLeftWidePtJetJetMC->GetBinContent(i))/weightSumNEffLW;
        Double_t weightedNEffLWErr    = pow((weightMinBiasNEffLW +  weightJetJetMCNEffLW),-0.5);
        
        if (isfinite(weightedNEffLW) && isfinite(weightedNEffLWErr)){
            histoEffiLeftWidePt->SetBinContent(i, weightedNEffLW);
            histoEffiLeftWidePt->SetBinError(i, weightedNEffLWErr);
        }
       
        
//      if (relErrJetJetMCWide < relErrMinBiasWide  && isfinite(relErrJetJetMCWide) ){
//          histoTrueEffiWidePt->SetBinContent(i,histoTrueEffiWidePtJetJetMC->GetBinContent(i));
//          histoTrueEffiWidePt->SetBinError(i,histoTrueEffiWidePtJetJetMC->GetBinError(i));
//      } else {
//          histoTrueEffiWidePt->SetBinContent(i,histoTrueEffiWidePtMinBias->GetBinContent(i));
//          histoTrueEffiWidePt->SetBinError(i,histoTrueEffiWidePtMinBias->GetBinError(i));
//          
//      }
    }
    
    TH1D* ratioMinBiasFinal = (TH1D*) histoTrueEffiPtMinBias->Clone("ratioMinBiasFinal");
    ratioMinBiasFinal->Divide(histoTrueEffiPtMinBias,histoTrueEffiPt , 1.,1.,"B");
    TH1D* ratioJetJetMCFinal = (TH1D*) histoTrueEffiPtJetJetMC->Clone("ratioJetJetMCFinal");
    ratioJetJetMCFinal->Divide(histoTrueEffiPtJetJetMC,histoTrueEffiPt , 1.,1.,"B");
    TH1D* ratioMinBiasJetJetMC = (TH1D*) histoTrueEffiPtMinBias->Clone("ratioJetJetMCFinal");
    ratioMinBiasJetJetMC->Divide(histoTrueEffiPtMinBias,histoTrueEffiPtJetJetMC , 1.,1.,"B");

    
    TCanvas* canvasFraction2 = new TCanvas("canvasFraction2","",1550,1200);// gives the page size
    canvasFraction2->SetTickx();
    canvasFraction2->SetTicky();
    canvasFraction2->SetGridx(0);
    canvasFraction2->SetGridy(0);
    canvasFraction2->SetLogy(0);
    canvasFraction2->SetLeftMargin(0.13);
    canvasFraction2->SetRightMargin(0.02);
    canvasFraction2->SetTopMargin(0.02);
    canvasFraction2->SetFillColor(0);


    DrawGammaSetMarker(ratioMinBiasFinal, 24, 1., kBlue+2, kBlue+2);
    DrawAutoGammaMesonHistos( ratioMinBiasFinal,
                    "", "#it{p}_{T} (GeV/#it{c})", "effi A/ effi B",
                    kFALSE, 5., 10e-10, kTRUE,
                    kTRUE, 0.5, 1.5,
                    kFALSE, 0., 7.9);
    
    DrawGammaSetMarker(ratioJetJetMCFinal, 25, 1., kRed, kRed);
    ratioJetJetMCFinal->Draw("same");
    DrawGammaSetMarker(ratioMinBiasJetJetMC, 26, 1., kGreen+2, kGreen+2);
    ratioMinBiasJetJetMC->Draw("same");
    canvasFraction2->Update();
    DrawGammaLines(0., 8.,1., 1.,0.1);

    TLegend* legendMultDataPP = new TLegend(0.15,0.8,0.4,0.95);
    legendMultDataPP->SetFillColor(0);
    legendMultDataPP->SetLineColor(0);
    legendMultDataPP->SetTextSize(0.04);
    legendMultDataPP->AddEntry(ratioJetJetMCFinal,"Jet Jet MC / final","p");
    legendMultDataPP->AddEntry(ratioMinBiasJetJetMC,"Min Bias/ Jet Jet MC","p");
    legendMultDataPP->AddEntry(ratioMinBiasFinal, "Min Bias/ final","p");
    legendMultDataPP->Draw();

    
    canvasFraction2->SaveAs(Form("%s/%s_MC_RatioComparisonEffiMerged_%s.%s",outputDir.Data(),mesonType.Data(),fCutSelection.Data(),fSuffix.Data()));

    
    //Drawing different Efficiencys
    TCanvas* canvasEffSimple = new TCanvas("canvasEffSimple","",200,10,1350,900);// gives the page size
    DrawGammaCanvasSettings( canvasEffSimple, 0.13, 0.02, 0.02, 0.09);
    canvasEffSimple->SetLogy(1);
                    
    DrawAutoGammaMesonHistos( histoTrueEffiPt, 
                                    "", "#it{p}_{T} (GeV/#it{c})", "True Efficiency", 
                                    kTRUE, 2., 3e-6, kFALSE,
                                    kFALSE, 0., 0.7, 
                                    kFALSE, 0., 10.);
            
    DrawGammaSetMarker(histoTrueEffiPt, 24, 1., kBlack, kBlack);
    histoTrueEffiPt->DrawCopy("e1");
    
    DrawGammaSetMarker(histoTrueEffiPtMinBias, 25, 1., kBlue, kBlue);
    histoTrueEffiPtMinBias->DrawCopy("e1,same");
    
    
    DrawGammaSetMarker(histoTrueEffiPtJetJetMC, 26, 1., kRed, kRed);
    histoTrueEffiPtJetJetMC->DrawCopy("e1,same");
    
    TLegend* legendYield = new TLegend(0.6,0.15,0.95,0.3);
    legendYield->SetTextSize(0.02);
    legendYield->SetFillColor(0);
    legendYield->AddEntry(histoTrueEffiPt,"merged ");
    legendYield->AddEntry(histoTrueEffiPtMinBias,"min Bias");
    legendYield->AddEntry(histoTrueEffiPtJetJetMC,"Jet Jet MC");
    legendYield->Draw();
        
        
    
    canvasEffSimple->Update();
    canvasEffSimple->SaveAs(Form("%s/%s_MC_ComparisonTrueEffiMerged_%s.%s",outputDir.Data(),mesonType.Data(),fCutSelection.Data(),fSuffix.Data()));

    DrawAutoGammaMesonHistos( histoEffiPt, 
                                    "", "#it{p}_{T} (GeV/#it{c})", "#epsilon_{rec}", 
                                    kTRUE, 2., 3e-6, kFALSE,
                                    kFALSE, 0., 0.7, 
                                    kFALSE, 0., 10.);
            
    DrawGammaSetMarker(histoEffiPt, 24, 1., kBlack, kBlack);
    histoEffiPt->DrawCopy("e1");
    
    DrawGammaSetMarker(histoEffiPtMinBias, 25, 1., kBlue, kBlue);
    histoEffiPtMinBias->DrawCopy("e1,same");
    
    
    DrawGammaSetMarker(histoEffiPtJetJetMC, 26, 1., kRed, kRed);
    histoEffiPtJetJetMC->DrawCopy("e1,same");
    
    legendYield->Draw();
            
    canvasEffSimple->Update();
    canvasEffSimple->SaveAs(Form("%s/%s_MC_ComparisonEffiMerged_%s.%s",outputDir.Data(),mesonType.Data(),fCutSelection.Data(),fSuffix.Data()));
    
        //Drawing different Efficiencys
    canvasEffSimple->SetLogy(0);
                    
    DrawAutoGammaMesonHistos( histoTrueEffiPt, 
                                    "", "#it{p}_{T} (GeV/#it{c})", "True Efficiency", 
                                    kTRUE, 0.75, 3e-6, kFALSE,
                                    kFALSE, 0., 0.7, 
                                    kFALSE, 0., 10.);
            
    histoTrueEffiPt->DrawCopy("e1");
    histoTrueEffiPtMinBias->DrawCopy("e1,same");
    histoTrueEffiPtJetJetMC->DrawCopy("e1,same");
    legendYield->Draw();
        
    canvasEffSimple->Update();
    canvasEffSimple->SaveAs(Form("%s/%s_MC_ComparisonTrueEffiMerged_linY_%s.%s",outputDir.Data(),mesonType.Data(),fCutSelection.Data(),fSuffix.Data()));

    DrawAutoGammaMesonHistos( histoEffiPt, 
                                    "", "#it{p}_{T} (GeV/#it{c})", "#epsilon_{rec}", 
                                    kTRUE, 0.75, 3e-6, kFALSE,
                                    kFALSE, 0., 0.7, 
                                    kFALSE, 0., 10.);
            
    histoEffiPt->DrawCopy("e1");
    histoEffiPtMinBias->DrawCopy("e1,same");
    histoEffiPtJetJetMC->DrawCopy("e1,same");
    legendYield->Draw();
        
    canvasEffSimple->Update();
    canvasEffSimple->SaveAs(Form("%s/%s_MC_ComparisonEffiMerged_linY_%s.%s",outputDir.Data(),mesonType.Data(),fCutSelection.Data(),fSuffix.Data()));

    delete canvasEffSimple;
    
    TFile* fileSum = new TFile (nameFileCorrectionFileFull,"UPDATE");
      histoTrueEffiPt->Write("TrueMesonEffiPt",TObject::kOverwrite);
      histoTrueEffiNarrowPt->Write("TrueMesonNarrowEffiPt",TObject::kOverwrite);
      histoTrueEffiWidePt->Write("TrueMesonWideEffiPt",TObject::kOverwrite);
      if(histoTrueEffiPtWOWeights) histoTrueEffiPtWOWeights->Write("TrueMesonEffiPtUnweighted",TObject::kOverwrite);
      if(histoTrueEffiNarrowPtWOWeights) histoTrueEffiNarrowPtWOWeights->Write("TrueMesonNarrowEffiPtUnweighted",TObject::kOverwrite);
      if(histoTrueEffiWidePtWOWeights) histoTrueEffiWidePtWOWeights->Write("TrueMesonWideEffiPtUnweighted",TObject::kOverwrite);

      histoEffiPt->Write("MesonEffiPt",TObject::kOverwrite);
      histoEffiNarrowPt->Write("MesonNarrowEffiPt",TObject::kOverwrite);
      histoEffiWidePt->Write("MesonWideEffiPt",TObject::kOverwrite);
      histoEffiLeftPt->Write("MesonLeftEffiPt",TObject::kOverwrite);
      histoEffiLeftNarrowPt->Write("MesonLeftNarrowEffiPt",TObject::kOverwrite);
      histoEffiLeftWidePt->Write("MesonLeftWideEffiPt",TObject::kOverwrite);
      histoAcceptance->Write("fMCMesonAccepPt",TObject::kOverwrite);
      if(histoAcceptanceWOWeights) histoAcceptanceWOWeights->Write("fMCMesonAccepPtWOWeights",TObject::kOverwrite);

      histoTrueMassMeson->Write("histoTrueMassMeson",TObject::kOverwrite);
      histoTrueFWHMMeson->Write("histoTrueFWHMMeson",TObject::kOverwrite);
      histoMassMeson->Write("histoMassMeson",TObject::kOverwrite);
      histoFWHMMeson->Write("histoFWHMMeson",TObject::kOverwrite);
      if (mesonType.Contains("Pi0")){
        histoTrueSecFrac->Write("TrueSecFrac",TObject::kOverwrite);
        histoTrueSecFracWide->Write("TrueSecFracWide",TObject::kOverwrite);
        histoTrueSecFracNarrow->Write("TrueSecFracNarrow",TObject::kOverwrite);
        histoTrueSecFracFromK0S->Write("TrueSecFracFromK0S",TObject::kOverwrite);
        histoTrueSecFracFromLambda->Write("TrueSecFracFromLambda",TObject::kOverwrite);
        histoTrueSecFracFromK0SWide->Write("TrueSecFracFromK0SWide",TObject::kOverwrite);
        histoTrueSecFracFromK0SNarrow->Write("TrueSecFracFromK0SNarrow",TObject::kOverwrite);
      }
      histoEventQualityJetJetMC->Write("NEvents_JetJetMC",TObject::kOverwrite);
      histoMCInputJetJetMC->Write("MC_Meson_genPt_oldBin_JetJetMC",TObject::kOverwrite);
      histoMCInputWOWeightingJetJetMC->Write("MC_Meson_genPt_WOWeights_JetJetMC",TObject::kOverwrite);
      histoMCInputWeightsJetJetMC->Write("MC_Meson_genPt_Weights_JetJetMC",TObject::kOverwrite);
      
      // adding pure MinBias histos
      histoTrueEffiPtMinBias->Write("TrueMesonEffiPt_MB",TObject::kOverwrite);
      histoTrueEffiNarrowPtMinBias->Write("TrueMesonNarrowEffiPt_MB",TObject::kOverwrite);
      histoTrueEffiWidePtMinBias->Write("TrueMesonWideEffiPt_MB",TObject::kOverwrite);
      histoEffiPtMinBias->Write("MesonEffiPt_MB",TObject::kOverwrite);
      histoEffiNarrowPtMinBias->Write("MesonNarrowEffiPt_MB",TObject::kOverwrite);
      histoEffiWidePtMinBias->Write("MesonWideEffiPt_MB",TObject::kOverwrite);
      histoEffiLeftPtMinBias->Write("MesonLeftEffiPt_MB",TObject::kOverwrite);
      histoEffiLeftNarrowPtMinBias->Write("MesonLeftNarrowEffiPt_MB",TObject::kOverwrite);
      histoEffiLeftWidePtMinBias->Write("MesonLeftWideEffiPt_MB",TObject::kOverwrite);
      histoAcceptance->Write("fMCMesonAccepPt_MB",TObject::kOverwrite);
      histoTrueMassMesonMinBias->Write("histoTrueMassMeson_MB",TObject::kOverwrite);
      histoTrueFWHMMesonMinBias->Write("histoTrueFWHMMeson_MB",TObject::kOverwrite);
      histoMassMesonMinBias->Write("histoMassMeson_MB",TObject::kOverwrite);
      histoFWHMMesonMinBias->Write("histoFWHMMeson_MB",TObject::kOverwrite);
      if (mesonType.Contains("Pi0")){
          histoTrueSecFracMinBias->Write("TrueSecFrac_MB",TObject::kOverwrite);
          histoTrueSecFracWideMinBias->Write("TrueSecFracWide_MB",TObject::kOverwrite);
          histoTrueSecFracNarrowMinBias->Write("TrueSecFracNarrow_MB",TObject::kOverwrite);
          histoTrueSecFracFromK0SMinBias->Write("TrueSecFracFromK0S_MB",TObject::kOverwrite);
          histoTrueSecFracFromLambdaMinBias->Write("TrueSecFracFromLambda_MB",TObject::kOverwrite);
          histoTrueSecFracFromK0SWideMinBias->Write("TrueSecFracFromK0SWide_MB",TObject::kOverwrite);
          histoTrueSecFracFromK0SNarrowMinBias->Write("TrueSecFracFromK0SNarrow_MB",TObject::kOverwrite);
      }
      // adding pure JetJetMC histos
      histoTrueEffiPtJetJetMC->Write("TrueMesonEffiPt_JJ",TObject::kOverwrite);
      histoTrueEffiNarrowPtJetJetMC->Write("TrueMesonNarrowEffiPt_JJ",TObject::kOverwrite);
      histoTrueEffiWidePtJetJetMC->Write("TrueMesonWideEffiPt_JJ",TObject::kOverwrite);
      histoEffiPtJetJetMC->Write("MesonEffiPt_JJ",TObject::kOverwrite);
      histoEffiNarrowPtJetJetMC->Write("MesonNarrowEffiPt_JJ",TObject::kOverwrite);
      histoEffiWidePtJetJetMC->Write("MesonWideEffiPt_JJ",TObject::kOverwrite);
      histoEffiLeftPtJetJetMC->Write("MesonLeftEffiPt_JJ",TObject::kOverwrite);
      histoEffiLeftNarrowPtJetJetMC->Write("MesonLeftNarrowEffiPt_JJ",TObject::kOverwrite);
      histoEffiLeftWidePtJetJetMC->Write("MesonLeftWideEffiPt_JJ",TObject::kOverwrite);
      histoAcceptance->Write("fMCMesonAccepPt_JJ",TObject::kOverwrite);
      histoTrueMassMesonJetJetMC->Write("histoTrueMassMeson_JJ",TObject::kOverwrite);
      histoTrueFWHMMesonJetJetMC->Write("histoTrueFWHMMeson_JJ",TObject::kOverwrite);
      histoMassMesonJetJetMC->Write("histoMassMeson_JJ",TObject::kOverwrite);
      histoFWHMMesonJetJetMC->Write("histoFWHMMeson_JJ",TObject::kOverwrite);
      if (mesonType.Contains("Pi0")){
          histoTrueSecFracJetJetMC->Write("TrueSecFrac_JJ",TObject::kOverwrite);
          histoTrueSecFracWideJetJetMC->Write("TrueSecFracWide_JJ",TObject::kOverwrite);
          histoTrueSecFracNarrowJetJetMC->Write("TrueSecFracNarrow_JJ",TObject::kOverwrite);
          histoTrueSecFracFromK0SJetJetMC->Write("TrueSecFracFromK0S_JJ",TObject::kOverwrite);
          histoTrueSecFracFromLambdaJetJetMC->Write("TrueSecFracFromLambda_JJ",TObject::kOverwrite);
          histoTrueSecFracFromK0SWideJetJetMC->Write("TrueSecFracFromK0SWide_JJ",TObject::kOverwrite);
          histoTrueSecFracFromK0SNarrowJetJetMC->Write("TrueSecFracFromK0SNarrow_JJ",TObject::kOverwrite);
      }
      
      fileSum->Write();
    fileSum->Close();
    
}
