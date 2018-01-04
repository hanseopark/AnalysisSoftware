/*****************************************************************************
******         provided by Gamma Conversion Group, PWGGA,               ******
******        Friederike Bock, friederike.bock@cern.ch                  ******
*****************************************************************************/

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

//**************************************************************************************************
//**************************** Main function *******************************************************
//**************************************************************************************************
void AdjustOutputFormat(    TString energy          = "",
                            TString centrality      = "",
                            Int_t mode              = -1,
                            Int_t optMode           = -1,
                            TString meson           = "",
                            TString inputFileName1  = "",
                            TString inputFileName2  = "",
                            TString outputFileName  = ""
                        ){

    StyleSettingsThesis();
    SetPlotStyle();

    cout << energy.Data() << "\t" << centrality << endl;
    cout << "File 1: " << inputFileName1.Data() << endl;
    cout << "File 2: " << inputFileName2.Data() << endl;
    cout << "Output file name: " << outputFileName.Data() << endl;

//     mode: 0....13 (PCM, PCM-EMC, PCM-PHOS, EMC, PHOS ...)
//     optMode:
//         0 - combine,
//         1 - copy and add hists to first folder struct
//         2 - bring to proper format

    //**********************************************************************************************************************
    TFile fileOutputFile(outputFileName.Data(),"RECREATE");

    TString outputDirName           = Form("%s%s%s", meson.Data(), energy.Data(), centrality.Data());
    fileOutputFile.mkdir(outputDirName.Data());
    TDirectoryFile* directoryWrite  = (TDirectoryFile*)fileOutputFile.Get(outputDirName.Data());

    //**********************************************************************************************************************
    if (optMode == 2){
        if (mode == 5 && energy.CompareTo("pPb_5.023TeV") == 0 && meson.CompareTo("Pi0") == 0){
            TFile* fileInput1                   = new TFile(inputFileName1.Data());
            TH1F* histoSysSources[7]            = {NULL, NULL, NULL, NULL, NULL, NULL};
            TString nameSourcesRead[7]          = {"Ext", "Acc", "Con", "Bge", "Non", "ESc", "TOF"};
            TString nameSourcesWrite[7]         = {"SignalExtraction", "Acceptance", "Material", "BGEstimate", "NonLinearity", "EnergyScale", "TimeOfFlight"};
            for (Int_t i= 0; i<7; i++ ){
                cout << "reading: " << Form("hSys%s",nameSourcesRead[i].Data()) << endl;
                histoSysSources[i]              = (TH1F*)fileInput1->Get(Form("hSys%s",nameSourcesRead[i].Data()));
                histoSysSources[i]->SetName(Form("histoSysSources_%s",nameSourcesWrite[i].Data()));
                directoryWrite->Add(histoSysSources[i]);
            }
            TH1F* histoSysTot                   = (TH1F*)fileInput1->Get("hSystTotal");
            histoSysTot->SetName("histoSysSources_Total");
            directoryWrite->Add(histoSysTot);
            TF1* fNonLinFunct                   = (TF1*)fileInput1->Get("fNon9");
            fNonLinFunct->SetName("NonLinearityCorrectionFunction");
            directoryWrite->Add(fNonLinFunct);
            TH1F* histoInvYieldStat             = (TH1F*)fileInput1->Get("hCor_stat");
            histoInvYieldStat->SetName("CorrectedYieldPi0");
            directoryWrite->Add(histoInvYieldStat);
            TH1F* histoInvYieldSys              = (TH1F*)fileInput1->Get("hCor_syst");
            histoInvYieldSys->SetName("histPi0SystError");
            directoryWrite->Add(histoInvYieldSys);
            TGraphAsymmErrors* graphInvYieldSys = new TGraphAsymmErrors(histoInvYieldSys);
            graphInvYieldSys->SetName("Pi0SystError");
            directoryWrite->Add(graphInvYieldSys);

            TH1F* histoMass                     = (TH1F*)fileInput1->Get("Gmean_Real");
            histoMass->SetName("Pi0_Mass_data");
            directoryWrite->Add(histoMass);
            TH1F* histoMassMC                   = (TH1F*)fileInput1->Get("Gmean_MC");
            histoMassMC->SetName("Pi0_Mass_MC");
            directoryWrite->Add(histoMassMC);
            TH1F* histoWidth                    = (TH1F*)fileInput1->Get("Gsigma_Real");
            histoWidth->SetName("Pi0_Width_data");
            TH1F* histoWidthMC                  = (TH1F*)fileInput1->Get("Gsigma_MC");
            directoryWrite->Add(histoWidth);
            histoWidthMC->SetName("Pi0_Width_MC");
            directoryWrite->Add(histoWidthMC);

            TFile* fileInput2                   = new TFile(inputFileName2.Data());
            TH1F* histoEffi                     = (TH1F*)fileInput2->Get("hEfficiency");
            histoEffi->SetName("EffTimesAccPi0");
            directoryWrite->Add(histoEffi);
            TF1* fitEffi                        = (TF1*)fileInput2->Get("fEfficiency_BC");
            fitEffi->SetName("fitEfficiencyPi0");
            directoryWrite->Add(fitEffi);

            fstream fileSysOutputCompilation;
            TString fFileSysErrName             = Form("SystematicErrorAveragedSinglePHOSPHOS_Pi0_pPb_5_023TeV_createdFromRootFile.dat");
            fileSysOutputCompilation.open(fFileSysErrName, ios::out);

            for (Int_t i = 1;i< histoInvYieldStat->GetNbinsX()+1; i++ ){
                if (i == 1) {
                    fileSysOutputCompilation << "Pt" << "\t";
                    for (Int_t k= 0; k<7; k++ ){
                        fileSysOutputCompilation << nameSourcesWrite[k] << "\t";
                    }
                    fileSysOutputCompilation << endl;
                }
                if (histoInvYieldStat->GetBinContent(i) > 0){
                    fileSysOutputCompilation <<  histoInvYieldStat->GetBinCenter(i) << "\t";
                    Double_t currentPt      = histoInvYieldStat->GetBinCenter(i);
                    Double_t currentTotSys  = 0;
                    for (Int_t k= 0; k<7; k++ ){
                        fileSysOutputCompilation << histoSysSources[k]->GetBinContent(histoSysSources[k]->FindBin(currentPt))*100 << "\t";
                        currentTotSys       = currentTotSys+TMath::Power(histoSysSources[k]->GetBinContent(histoSysSources[k]->FindBin(currentPt))*100,2);
                    }
                    fileSysOutputCompilation << TMath::Sqrt(currentTotSys) << endl;
                }
            }
        }
    }

    //**********************************************************************************************************************
    directoryWrite->Write();
    fileOutputFile.Close();

}
