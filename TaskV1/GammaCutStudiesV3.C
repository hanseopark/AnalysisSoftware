#include <Riostream.h>
#include <fstream>
#include "TMath.h"
#include <stdio.h>
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

struct SysErrorConversion {
   Double_t value;
   Double_t error;
};


void PlotCanvas( Int_t i, 
                 Int_t number, 
                 TCanvas *canvas, 
                 TH1D *spectrum, 
                 TLegend *legend, 
                 TCanvas *ratiocanvas, 
                 TH1D *ratio, 
                 TLegend *ratiolegend, 
                 TString cutSelection, 
                 TF1 *OneA = NULL, 
                 TF1 *OneB = NULL){
   

    if(spectrum->GetMinimum()<=0.0)
        spectrum->SetMinimum(spectrum->GetBinContent(spectrum->GetNbinsX())/10.);

    canvas->cd();
    spectrum->GetYaxis()->SetTitleSize(0.035);
    spectrum->GetXaxis()->SetTitleSize(0.035);
    spectrum->GetYaxis()->SetLabelSize(0.03);
    spectrum->GetXaxis()->SetLabelSize(0.03);
    spectrum->GetYaxis()->SetTitleOffset(1.);
    spectrum->GetXaxis()->SetTitleOffset(1.);
    if(!i){
        spectrum->DrawCopy();
        if(OneB)OneB->Draw("same");
    }
    else spectrum->DrawCopy("same");
    legend->AddEntry(spectrum,cutSelection,"p");
    if(i==number-1) legend->Draw();
    ratiocanvas->cd();
    //    ratiocanvas->SetGridx();
    //    ratiocanvas->SetGridy();
    if(!i){
        ratio->DrawCopy("e1][");
        OneA->Draw("same");
    }
    else ratio->DrawCopy("e1same][");
    ratiolegend->AddEntry(ratio,cutSelection,"p");
    if(i==number-1) ratiolegend->Draw();
   
}

void CalculateSystematicsGraphs( TH1D** histoArray, 
                                 TH1D**histoArrayBound, 
                                 TString* cutSelectionOut,
                                 Int_t numberOfCuts, 
                                 TString prefix, TString prefixCommonSysFile, 
                                 TString outputDir, TString outputFileDir,
                                 TString cutVariationName
    
){
    if (numberOfCuts<2) return;
    
    Int_t NBinsPt = histoArray[0]->GetNbinsX();
    const Int_t NBinstPtConst       = NBinsPt+1;

    const Int_t ConstNumberOfCuts   = numberOfCuts;

    Double_t  BinsXCenter[NBinstPtConst];
    Double_t  BinsXWidth[NBinstPtConst];
    Double_t BinValue[NBinstPtConst];
    BinsXCenter[0] = 0;
    BinsXWidth[0]=0.;
    BinValue[0]=0.;
    for (Int_t i = 1; i < NBinsPt +1; i++){
        BinsXCenter[i] = histoArray[0]->GetBinCenter(i);
        BinsXWidth[i]= histoArray[0]->GetBinWidth(i)/2.;
    }

    Double_t BinValueGamma[NBinstPtConst];
    BinValueGamma[0]=0.;

    SysErrorConversion SysErrCutGamma[ConstNumberOfCuts][NBinstPtConst];
    SysErrorConversion SysErrCutGammaRaw[ConstNumberOfCuts][NBinstPtConst];

    for (Int_t j = 0; j < numberOfCuts; j++){
        for (Int_t i = 1; i < NBinsPt +1; i++){
            BinValueGamma[i]= histoArray[0]->GetBinContent(i);
            SysErrCutGamma[j][i].value = histoArray[j]->GetBinContent(i);
            SysErrCutGamma[j][i].error = histoArray[j]->GetBinError(i);
            SysErrCutGammaRaw[j][i].value = histoArrayBound[j]->GetBinContent(i);
            SysErrCutGammaRaw[j][i].error = histoArrayBound[j]->GetBinError(i);
        }
    }

    Double_t DifferenceCutGamma[ConstNumberOfCuts][NBinstPtConst];
    Double_t DifferenceErrorCutGamma[ConstNumberOfCuts][NBinstPtConst];

    Double_t LargestDiffGammaNeg[NBinstPtConst];
    Double_t LargestDiffGammaPos[NBinstPtConst];
    Double_t LargestDiffGammaErrorNeg[NBinstPtConst];
    Double_t LargestDiffGammaErrorPos[NBinstPtConst];

    Double_t LargestDiffGammaRelNeg[NBinstPtConst];
    Double_t LargestDiffGammaRelPos[NBinstPtConst];
    Double_t LargestDiffGammaRelErrorNeg[NBinstPtConst];
    Double_t LargestDiffGammaRelErrorPos[NBinstPtConst];

    Double_t RelDifferenceCutGamma[ConstNumberOfCuts][NBinstPtConst];
    Double_t RelDifferenceErrorCutGamma[ConstNumberOfCuts][NBinstPtConst];
    Double_t RelDifferenceRawCut[ConstNumberOfCuts][NBinstPtConst];
            
    for (Int_t j = 1; j < numberOfCuts; j++){
        for ( Int_t i = 1; i < NBinstPtConst; i++) {
            DifferenceCutGamma[j][i]=0.;
            DifferenceErrorCutGamma[j][i]=0.;
            LargestDiffGammaNeg[i]=0.;
            LargestDiffGammaPos[i]=0.;
            LargestDiffGammaErrorNeg[i]=0.;
            LargestDiffGammaErrorPos[i]=0.;
            RelDifferenceCutGamma[j][i]=0.;
            RelDifferenceRawCut[j][i]=0.;
            RelDifferenceErrorCutGamma[j][i]=0.;
        }
    }

    for(Int_t j = 1; j < numberOfCuts; j++){
        for (Int_t i = 1; i < NBinsPt +1; i++){
            //Calculate differences
            DifferenceCutGamma[j][i] = SysErrCutGamma[j][i].value - SysErrCutGamma[0][i].value;
            DifferenceErrorCutGamma[j][i] = TMath::Sqrt(TMath::Abs(TMath::Power(SysErrCutGamma[j][i].error,2)-TMath::Power(SysErrCutGamma[0][i].error,2)));
            if(SysErrCutGamma[0][i].value != 0){
                RelDifferenceCutGamma[j][i] = DifferenceCutGamma[j][i]/SysErrCutGamma[0][i].value*100. ;
                RelDifferenceErrorCutGamma[j][i] = DifferenceErrorCutGamma[j][i]/SysErrCutGamma[0][i].value*100. ;
            } else {
                RelDifferenceCutGamma[j][i] = -10000.;
                RelDifferenceErrorCutGamma[j][i] = 100. ;
            }
            if(SysErrCutGammaRaw[0][i].value != 0){
                RelDifferenceRawCut[j][i] = (SysErrCutGammaRaw[j][i].value - SysErrCutGammaRaw[0][i].value)/SysErrCutGammaRaw[0][i].value*100. ;
            } else {
                RelDifferenceRawCut[j][i] = -10000.;
            }
                    
            if(DifferenceCutGamma[j][i] < 0){
                if (TMath::Abs(LargestDiffGammaNeg[i]) < TMath::Abs(DifferenceCutGamma[j][i]) && RelDifferenceRawCut[j][i] > -75.){
                LargestDiffGammaNeg[i] = DifferenceCutGamma[j][i];
                LargestDiffGammaErrorNeg[i] = DifferenceErrorCutGamma[j][i];
                }
            }else{
                if (TMath::Abs(LargestDiffGammaPos[i]) < TMath::Abs(DifferenceCutGamma[j][i]) && RelDifferenceRawCut[j][i] > -75.){
                LargestDiffGammaPos[i] = DifferenceCutGamma[j][i];
                LargestDiffGammaErrorPos[i] = DifferenceErrorCutGamma[j][i];
                }
            }
        }
    }

    cout << "done filling" << endl;
    const char *SysErrDatnameGamma = Form("%s/%s_SystematicErrorCutStudies.dat",prefix.Data(),outputDir.Data());
    fstream SysErrDatGamma;
    SysErrDatGamma.open(SysErrDatnameGamma, ios::out);
    SysErrDatGamma << "Calculation of the systematic error due to the yield cuts" << endl;

    for (Int_t l=0; l< numberOfCuts; l++){
        if (l == 0) {
            SysErrDatGamma << endl <<"Bin" << "\t" << cutSelectionOut[l] << "\t" <<endl;
            for(Int_t i = 1; i < (NBinsPt +1); i++){
                SysErrDatGamma << BinsXCenter[i] << "\t" << SysErrCutGamma[l][i].value << "\t" << SysErrCutGamma[l][i].error << endl;   
            }
        } else{
            SysErrDatGamma << endl <<"Bin" << "\t" << cutSelectionOut[l] << "\t" << "Error " << "\t Dif to Cut1" << endl;
            for(Int_t i = 1; i < (NBinsPt +1); i++){
                if (RelDifferenceRawCut[l][i] > -75.){
                    SysErrDatGamma  << BinsXCenter[i] << "\t" << SysErrCutGamma[l][i].value << "\t" << SysErrCutGamma[l][i].error << "\t" <<  DifferenceCutGamma[l][i] << "\t"
                                    << DifferenceErrorCutGamma[l][i] << "\t"<< RelDifferenceCutGamma[l][i] <<  "\t" << RelDifferenceErrorCutGamma[l][i] <<"\t" 
                                    << RelDifferenceRawCut[l][i]<< endl;
                } else {
                    SysErrDatGamma  << BinsXCenter[i] << "\t" << SysErrCutGamma[l][i].value << "\t" << SysErrCutGamma[l][i].error << "\t" <<  DifferenceCutGamma[l][i] << "\t"
                                    << DifferenceErrorCutGamma[l][i] << "\t"<< RelDifferenceCutGamma[l][i] <<  "\t" << RelDifferenceErrorCutGamma[l][i] <<"\t" 
                                    << RelDifferenceRawCut[l][i]  <<"\t not considered in largest dev" <<endl;
                }
            }
        }
    }


    SysErrDatGamma << endl;
    SysErrDatGamma << endl;
    SysErrDatGamma << "Bin" << "\t" << "Largest Dev Neg" << "\t" << "Largest Dev Pos"  << endl;
    for(Int_t i = 1; i < (NBinsPt +1); i++){
        SysErrDatGamma << BinsXCenter[i]  << "\t" << LargestDiffGammaNeg[i] << "\t" <<LargestDiffGammaErrorNeg[i]<< "\t" << LargestDiffGammaPos[i] << "\t" << LargestDiffGammaErrorPos[i]<<endl;
    }
    SysErrDatGamma << endl << endl <<"Bin" << "\t" << "Largest Dev Neg rel" << "\t" << "Largest Dev Pos rel"  << endl;
    for(Int_t i = 0; i < (NBinsPt +1); i++){
        if ( SysErrCutGamma[0][i].value != 0.){
            LargestDiffGammaRelNeg[i] = - LargestDiffGammaNeg[i]/SysErrCutGamma[0][i].value*100.;
            LargestDiffGammaRelPos[i] = LargestDiffGammaPos[i]/SysErrCutGamma[0][i].value*100.;
            LargestDiffGammaRelErrorNeg[i] = - LargestDiffGammaErrorNeg[i]/SysErrCutGamma[0][i].value*100.;
            LargestDiffGammaRelErrorPos[i] = LargestDiffGammaErrorPos[i]/SysErrCutGamma[0][i].value*100.;
            if (i > 0) SysErrDatGamma   << BinsXCenter[i] << "\t" << LargestDiffGammaNeg[i]/SysErrCutGamma[0][i].value*100. << "\t" 
                                        << LargestDiffGammaErrorNeg[i]/SysErrCutGamma[0][i].value*100. << "\t" << LargestDiffGammaPos[i]/SysErrCutGamma[0][i].value*100. << "\t" 
                                        << LargestDiffGammaErrorPos[i]/SysErrCutGamma[0][i].value*100.<<endl;
        } else {
            LargestDiffGammaRelNeg[i] = 0.;
            LargestDiffGammaRelPos[i] = 0.;
            LargestDiffGammaRelErrorNeg[i] = 0.;
            LargestDiffGammaRelErrorPos[i] = 0.;
        }
    }

    SysErrDatGamma.close();
                
    TGraphAsymmErrors* SystErrGraphNegGamma = new TGraphAsymmErrors(NBinsPt+1, BinsXCenter, LargestDiffGammaRelNeg, BinsXWidth, BinsXWidth, LargestDiffGammaRelErrorNeg, LargestDiffGammaRelErrorNeg);
    SystErrGraphNegGamma->SetName(Form("%s_SystErrorRelNeg_%s", prefix.Data(), cutVariationName.Data()));
    TGraphAsymmErrors* SystErrGraphPosGamma = new TGraphAsymmErrors(NBinsPt+1, BinsXCenter, LargestDiffGammaRelPos, BinsXWidth, BinsXWidth, LargestDiffGammaRelErrorPos, LargestDiffGammaRelErrorPos);
    SystErrGraphPosGamma->SetName(Form("%s_SystErrorRelPos_%s", prefix.Data(), cutVariationName.Data()));
    const char* Outputname = Form("%s/%s_SystematicErrorCuts.root", outputFileDir.Data(), prefixCommonSysFile.Data());
    TFile* SystematicErrorFile = new TFile(Outputname,"UPDATE");
    SystErrGraphPosGamma->Write(Form("%s_SystErrorRelPos_%s_%s",prefix.Data(), cutVariationName.Data(),(GetCentralityString(cutSelectionOut[0])).Data()),TObject::kOverwrite);
    SystErrGraphNegGamma->Write(Form("%s_SystErrorRelNeg_%s_%s",prefix.Data(), cutVariationName.Data(),(GetCentralityString(cutSelectionOut[0])).Data()),TObject::kOverwrite);
    SystematicErrorFile->Write();
    SystematicErrorFile->Close();
        
    delete SystErrGraphNegGamma;
    delete SystErrGraphPosGamma;
    delete SystematicErrorFile;

    return;
}    

void CalculateSystematicsGraphsWOBound( TH1D** histoArray, 
                                        TString* cutSelectionOut,
                                        Int_t numberOfCuts, 
                                        TString prefix, TString prefixCommonSysFile, 
                                        TString outputDir, TString outputFileDir,
                                        TString cutVariationName
){
    if (numberOfCuts<2) return;
    Int_t NBinsPt = histoArray[0]->GetNbinsX();
    const Int_t NBinstPtConst       = NBinsPt+1;

    const Int_t ConstNumberOfCuts   = numberOfCuts;

    Double_t  BinsXCenter[NBinstPtConst];
    Double_t  BinsXWidth[NBinstPtConst];
    Double_t BinValue[NBinstPtConst];
    BinsXCenter[0] = 0;
    BinsXWidth[0]=0.;
    BinValue[0]=0.;
    for (Int_t i = 1; i < NBinsPt +1; i++){
        BinsXCenter[i] = histoArray[0]->GetBinCenter(i);
        BinsXWidth[i]= histoArray[0]->GetBinWidth(i)/2.;
    }

    Double_t BinValueGamma[NBinstPtConst];
    BinValueGamma[0]=0.;

    SysErrorConversion SysErrCutGamma[ConstNumberOfCuts][NBinstPtConst];

    for (Int_t j = 0; j < numberOfCuts; j++){
        for (Int_t i = 1; i < NBinsPt +1; i++){
            BinValueGamma[i]= histoArray[0]->GetBinContent(i);
            SysErrCutGamma[j][i].value = histoArray[j]->GetBinContent(i);
            SysErrCutGamma[j][i].error = histoArray[j]->GetBinError(i);
        }
    }

    Double_t DifferenceCutGamma[ConstNumberOfCuts][NBinstPtConst];
    Double_t DifferenceErrorCutGamma[ConstNumberOfCuts][NBinstPtConst];

    Double_t LargestDiffGammaNeg[NBinstPtConst];
    Double_t LargestDiffGammaPos[NBinstPtConst];
    Double_t LargestDiffGammaErrorNeg[NBinstPtConst];
    Double_t LargestDiffGammaErrorPos[NBinstPtConst];

    Double_t LargestDiffGammaRelNeg[NBinstPtConst];
    Double_t LargestDiffGammaRelPos[NBinstPtConst];
    Double_t LargestDiffGammaRelErrorNeg[NBinstPtConst];
    Double_t LargestDiffGammaRelErrorPos[NBinstPtConst];

    Double_t RelDifferenceCutGamma[ConstNumberOfCuts][NBinstPtConst];
    Double_t RelDifferenceErrorCutGamma[ConstNumberOfCuts][NBinstPtConst];
            
    for (Int_t j = 1; j < numberOfCuts; j++){
        for ( Int_t i = 1; i < NBinstPtConst; i++) {
            DifferenceCutGamma[j][i]=0.;
            DifferenceErrorCutGamma[j][i]=0.;
            LargestDiffGammaNeg[i]=0.;
            LargestDiffGammaPos[i]=0.;
            LargestDiffGammaErrorNeg[i]=0.;
            LargestDiffGammaErrorPos[i]=0.;
            RelDifferenceCutGamma[j][i]=0.;
            RelDifferenceErrorCutGamma[j][i]=0.;
        }
    }

    for(Int_t j = 1; j < numberOfCuts; j++){
        for (Int_t i = 1; i < NBinsPt +1; i++){
            //Calculate differences
            DifferenceCutGamma[j][i] = SysErrCutGamma[j][i].value - SysErrCutGamma[0][i].value;
            DifferenceErrorCutGamma[j][i] = TMath::Sqrt(TMath::Abs(TMath::Power(SysErrCutGamma[j][i].error,2)-TMath::Power(SysErrCutGamma[0][i].error,2)));
            if(SysErrCutGamma[0][i].value != 0){
                RelDifferenceCutGamma[j][i] = DifferenceCutGamma[j][i]/SysErrCutGamma[0][i].value*100. ;
                RelDifferenceErrorCutGamma[j][i] = DifferenceErrorCutGamma[j][i]/SysErrCutGamma[0][i].value*100. ;
            } else {
                RelDifferenceCutGamma[j][i] = -10000.;
                RelDifferenceErrorCutGamma[j][i] = 100. ;
            }
                    
            if(DifferenceCutGamma[j][i] < 0){
                if (TMath::Abs(LargestDiffGammaNeg[i]) < TMath::Abs(DifferenceCutGamma[j][i])){
                LargestDiffGammaNeg[i] = DifferenceCutGamma[j][i];
                LargestDiffGammaErrorNeg[i] = DifferenceErrorCutGamma[j][i];
                }
            }else{
                if (TMath::Abs(LargestDiffGammaPos[i]) < TMath::Abs(DifferenceCutGamma[j][i])){
                LargestDiffGammaPos[i] = DifferenceCutGamma[j][i];
                LargestDiffGammaErrorPos[i] = DifferenceErrorCutGamma[j][i];
                }
            }
        }
    }

    cout << "done filling" << endl;
    const char *SysErrDatnameGamma = Form("%s/%s_SystematicErrorCutStudies.dat",outputDir.Data(),prefix.Data());
    fstream SysErrDatGamma;
    SysErrDatGamma.open(SysErrDatnameGamma, ios::out);
    SysErrDatGamma << "Calculation of the systematic error due to the yield cuts" << endl;

    cout << "works" << endl;
    for (Int_t l=0; l< numberOfCuts; l++){
        if (l == 0) {
            SysErrDatGamma << endl <<"Bin" << "\t" << cutSelectionOut[l] << "\t" <<endl;
            for(Int_t i = 1; i < (NBinsPt +1); i++){
                SysErrDatGamma << BinsXCenter[i] << "\t" << SysErrCutGamma[l][i].value << "\t" << SysErrCutGamma[l][i].error << endl;   
            }
        } else{
            SysErrDatGamma << endl <<"Bin" << "\t" << cutSelectionOut[l] << "\t" << "Error " << "\t Dif to Cut1" << endl;
            for(Int_t i = 1; i < (NBinsPt +1); i++){
                SysErrDatGamma  << BinsXCenter[i] << "\t" << SysErrCutGamma[l][i].value << "\t" << SysErrCutGamma[l][i].error << "\t" <<  DifferenceCutGamma[l][i] << "\t"
                                << DifferenceErrorCutGamma[l][i] << "\t"<< RelDifferenceCutGamma[l][i] <<  "\t" << RelDifferenceErrorCutGamma[l][i] << endl;
            }
        }
    }


    SysErrDatGamma << endl;
    SysErrDatGamma << endl;
    SysErrDatGamma << "Bin" << "\t" << "Largest Dev Neg" << "\t" << "Largest Dev Pos"  << endl;
    for(Int_t i = 1; i < (NBinsPt +1); i++){
        SysErrDatGamma << BinsXCenter[i]  << "\t" << LargestDiffGammaNeg[i] << "\t" <<LargestDiffGammaErrorNeg[i]<< "\t" << LargestDiffGammaPos[i] << "\t" << LargestDiffGammaErrorPos[i]<<endl;
    }
    SysErrDatGamma << endl << endl <<"Bin" << "\t" << "Largest Dev Neg rel" << "\t" << "Largest Dev Pos rel"  << endl;
    for(Int_t i = 0; i < (NBinsPt +1); i++){
        if ( SysErrCutGamma[0][i].value != 0.){
            LargestDiffGammaRelNeg[i] = - LargestDiffGammaNeg[i]/SysErrCutGamma[0][i].value*100.;
            LargestDiffGammaRelPos[i] = LargestDiffGammaPos[i]/SysErrCutGamma[0][i].value*100.;
            LargestDiffGammaRelErrorNeg[i] = - LargestDiffGammaErrorNeg[i]/SysErrCutGamma[0][i].value*100.;
            LargestDiffGammaRelErrorPos[i] = LargestDiffGammaErrorPos[i]/SysErrCutGamma[0][i].value*100.;
            if (i > 0) SysErrDatGamma   << BinsXCenter[i] << "\t" << LargestDiffGammaNeg[i]/SysErrCutGamma[0][i].value*100. << "\t" 
                                        << LargestDiffGammaErrorNeg[i]/SysErrCutGamma[0][i].value*100. << "\t" << LargestDiffGammaPos[i]/SysErrCutGamma[0][i].value*100. << "\t" 
                                        << LargestDiffGammaErrorPos[i]/SysErrCutGamma[0][i].value*100.<<endl;
        } else {
            LargestDiffGammaRelNeg[i] = 0.;
            LargestDiffGammaRelPos[i] = 0.;
            LargestDiffGammaRelErrorNeg[i] = 0.;
            LargestDiffGammaRelErrorPos[i] = 0.;
        }
    }

    SysErrDatGamma.close();
                
    TGraphAsymmErrors* SystErrGraphNegGamma = new TGraphAsymmErrors(NBinsPt+1, BinsXCenter, LargestDiffGammaRelNeg, BinsXWidth, BinsXWidth, LargestDiffGammaRelErrorNeg, LargestDiffGammaRelErrorNeg);
    SystErrGraphNegGamma->SetName(Form("%s_SystErrorRelNeg_%s", prefix.Data(), cutVariationName.Data()));
    TGraphAsymmErrors* SystErrGraphPosGamma = new TGraphAsymmErrors(NBinsPt+1, BinsXCenter, LargestDiffGammaRelPos, BinsXWidth, BinsXWidth, LargestDiffGammaRelErrorPos, LargestDiffGammaRelErrorPos);
    SystErrGraphPosGamma->SetName(Form("%s_SystErrorRelPos_%s", prefix.Data(), cutVariationName.Data()));
    const char* Outputname = Form("%s/%s_SystematicErrorCuts.root", prefixCommonSysFile.Data(),outputFileDir.Data());
    TFile* SystematicErrorFile = new TFile(Outputname,"UPDATE");
    SystErrGraphPosGamma->Write(Form("%s_SystErrorRelPos_%s_%s",prefix.Data(), cutVariationName.Data(),(GetCentralityString(cutSelectionOut[0])).Data()),TObject::kOverwrite);
    SystErrGraphNegGamma->Write(Form("%s_SystErrorRelNeg_%s_%s",prefix.Data(), cutVariationName.Data(),(GetCentralityString(cutSelectionOut[0])).Data()),TObject::kOverwrite);
    SystematicErrorFile->Write();
    SystematicErrorFile->Close();
        
    delete SystErrGraphNegGamma;
    delete SystErrGraphPosGamma;
    delete SystematicErrorFile;

    return;
}    

//***********************************************************************************************************************************************
//*********************************** CutVariation studies for photons- main routing ************************************************************
//***********************************************************************************************************************************************
void GammaCutStudiesV3(TString cutFile = "CombineCuts.dat",TString energy="",TString cutVariationName ="",TString suffix = "eps", Int_t mode = 9){
    
    //****************************************************************************************
    //*************************** Catch old versions *****************************************
    //****************************************************************************************
    if (mode == 1 || mode > 5){
        cout << "ERROR: This macro isn't designed to run for Dalitz or the old software version" << endl;
        return;
    } else if (mode == 2 || mode == 3 ){    
        cout << "WARNING: running hybrid mode, this macro is still under construction for this mode" << endl;    
    } else if ( mode == 4 || mode == 5){
        cout << "This macro can't yet deal with these modi" << endl;
        return;
    }    

    
    //*****************************************************************************************
    //*************************** Flag for processing *****************************************
    //*****************************************************************************************
    Bool_t haveOutputGammaToPi0 = 1;
    
    
    //*****************************************************************************************
    //************************** Set general style settings ***********************************
    //*****************************************************************************************
    StyleSettingsThesis();  
    SetPlotStyle();

    //*****************************************************************************************
    //************************** Setting output directories ***********************************
    //*****************************************************************************************
    TString outputDir           = Form("CutStudies/%s/%s",energy.Data(),cutVariationName.Data());
    TString outputFileDir       = Form("CutStudies/%s",energy.Data());
    gSystem->Exec("mkdir -p "+outputDir);


    //*******************************************************************************************
    //************************ More specific styling and labeling issues ************************
    //*******************************************************************************************
    TString collisionSystem     = ReturnFullCollisionsSystem(energy);   
    if (collisionSystem.CompareTo("") == 0){
        cout << "No correct collision system specification, has been given" << endl;
        return;     
    }
    TString detectionProcess    = ReturnFullTextReconstructionProcess(mode);
    TString process             = "#gamma";

    Color_t color[20]           = { kBlack, kAzure, kGreen+2, kOrange+2, kRed, 
                                    kViolet, kBlue-9, kSpring+10, kCyan+3, kCyan-10, 
                                    kCyan, kGreen+4, kGreen-9, kGreen,  kYellow+4, 
                                    kYellow+3, kMagenta+4, kMagenta-8, kGray, kGray+3
        
                                  };
    
//                                     {kBlack, kRed+1, kBlue+1, kGreen+1, kMagenta+1, 
//                                    809, kCyan+1, kGreen-2, kAzure+2, kGray};
    Int_t markerType            = 24;
    

    TF1 *One                    = new TF1("One","1",0,25);
    One->SetLineWidth(1.2);
    One->SetLineColor(1);

    
    //*******************************************************************************************
    //************************** Initialization of variables ************************************
    //*******************************************************************************************
    cout << endl << endl <<endl;
    cout << "Processing '" << cutVariationName <<"' Study for " << collisionSystem <<endl;
    cout<<endl;
    TString cutStringsName          [50];
    TString cutSelection            [50];
    
    
    //*******************************************************************************************
    //**************************** Reading and deciphering cutstrings ***************************
    //*******************************************************************************************
    ifstream in(cutFile);
    cout<<"=========================="<<endl;
    cout<<"Available Cuts:"<<endl;
    string currentCutNumber;
    Int_t number                = 0;

    while(getline(in, currentCutNumber)){
        cutSelection[number]    = currentCutNumber;
        cout<<"---> " << currentCutNumber << endl;
        number++;
    }

    cout<<"=========================="<<endl;


    //*******************************************************************************************
    //***************************** Initialization of histo arrays ******************************
    //*******************************************************************************************
    TString nameCurrentCorrectedFile            = "";
    TH1D **histoIncGamma                        = new TH1D*[number];
    TH1D **histoIncGammaRatio                   = new TH1D*[number];
    TH1D **histoPi0Spectrum                     = new TH1D*[number];
    TH1D **histoPi0SpectrumRatio                = new TH1D*[number];
    TH1D **histoPi0SpectrumFit                  = new TH1D*[number];
    TH1D **histoPi0SpectrumFitRatio             = new TH1D*[number];
    TH1D **histoIncGammaToPi0Ratio              = new TH1D*[number];
    TH1D **histoIncGammaToPi0RatioRatio         = new TH1D*[number];
    TH1D **histoIncGammaToPi0RatioFit           = new TH1D*[number];
    TH1D **histoIncGammaToPi0RatioFitRatio      = new TH1D*[number];
    TH1D **histoDR                              = new TH1D*[number];
    TH1D **histoDRRatio                         = new TH1D*[number];
    TH1D **histoDRFit                           = new TH1D*[number];
    TH1D **histoDRFitRatio                      = new TH1D*[number];

    TH1D **histoPurity                          = new TH1D*[number];
    TH1D **histoPurityRatio                     = new TH1D*[number];
    TH1D **histoGammaEff                        = new TH1D*[number];
    TH1D **histoGammaEffRatio                   = new TH1D*[number];
    TH1D **histoGammaResolCorr                  = new TH1D*[number];
    TH1D **histoGammaResolCorrRatio             = new TH1D*[number];
    TH1D **histoRawGamma                        = new TH1D*[number];
    TH1D **histoRawGammaRatio                   = new TH1D*[number];
    TH1D **histoGammaEffRecPt                   = new TH1D*[number];
    TH1D **histoGammaEffRecPtRatio              = new TH1D*[number];
    TH1D **histoGammaConvProb                   = new TH1D*[number];
    TH1D **histoGammaConvProbRatio              = new TH1D*[number];
    TH1D **histoGammaCorrFac                    = new TH1D*[number];
    TH1D **histoGammaCorrFacRatio               = new TH1D*[number];
    
    TFile **fileCurrentFinal                    = new TFile*[number];
    TFile **fileCurrentCorrection               = new TFile*[number];

    //*******************************************************************************************
    //*****************************Initialization of Canvases ***********************************
    //*******************************************************************************************    
    TCanvas *canvasGammaSpectrum                = GetAndSetCanvas("GammaSpectra");
    TCanvas *canvasPi0Spectrum                  = GetAndSetCanvas("Pi0Spectra");canvasPi0Spectrum->SetLogy();
    TCanvas *canvasPi0SpectrumFit               = GetAndSetCanvas("Pi0FitSpectra");canvasPi0SpectrumFit->SetLogy();
    TCanvas *canvasIncGammaToPi0Ratio           = GetAndSetCanvas("InclusiveRatios");
    TCanvas *canvasIncGammaToPi0RatioFit        = GetAndSetCanvas("InclusiveRatiosFit");
    TCanvas *canvasDR                           = GetAndSetCanvas("DoubleRatios");
    TCanvas *canvasDRFit                        = GetAndSetCanvas("DoubleRatiosFit");

    TCanvas *canvasGammaSpectrumRatio           = GetAndSetCanvas("GammaSpectraRatio");
    TCanvas *canvasPi0SpectrumRatio             = GetAndSetCanvas("Pi0SpectraRatio");
    TCanvas *canvasPi0SpectrumFitRatio          = GetAndSetCanvas("Pi0FitSpectraRatio");
    TCanvas *canvasIncGammaToPi0RatioRatio      = GetAndSetCanvas("InclusiveRatiosRatio");
    TCanvas *canvasIncGammaToPi0RatioFitRatio   = GetAndSetCanvas("InclusiveRatiosFitRatio");
    TCanvas *canvasDRRatio                      = GetAndSetCanvas("DoubleRatiosRatio");
    TCanvas *canvasDRFitRatio                   = GetAndSetCanvas("DoubleRatiosFitRatio");
    
    TCanvas *canvasRawGamma                     = GetAndSetCanvas("RawGamma");canvasRawGamma->SetLogy();
    TCanvas *canvasRawGammaRatio                = GetAndSetCanvas("RawGammaRatio");
    TCanvas *canvasPurity                       = GetAndSetCanvas("Purity");
    TCanvas *canvasPurityRatio                  = GetAndSetCanvas("PurityRatio");
    TCanvas *canvasGammaEff                     = GetAndSetCanvas("GammaEff");
    TCanvas *canvasGammaEffRatio                = GetAndSetCanvas("GammaEffRatio");
    TCanvas *canvasGammaEffRecPt                = GetAndSetCanvas("GammaEffRecPt");
    TCanvas *canvasGammaEffRecPtRatio           = GetAndSetCanvas("GammaEffRecPtRatio");
    TCanvas *canvasGammaConvProb                = GetAndSetCanvas("GammaConvProb");
    TCanvas *canvasGammaConvProbRatio           = GetAndSetCanvas("GammaConvProbRatio");
    TCanvas *canvasGammaResolCorr               = GetAndSetCanvas("GammaResolCorr");
    TCanvas *canvasGammaResolCorrRatio          = GetAndSetCanvas("GammaResolCorrRatio");
    TCanvas *canvasGammaCorrFac                 = GetAndSetCanvas("GammaCorrFac");
    TCanvas *canvasGammaCorrFacRatio            = GetAndSetCanvas("GammaCorrFacRatio");
    
    //*******************************************************************************************
    //*****************************Initialization of Canvases ***********************************
    //*******************************************************************************************        
    TLegend *legendGammaSpectrum                = GetAndSetLegend(0.3,0.7,number);
    TLegend *legendGammaSpectrumRatio           = GetAndSetLegend(0.15,0.75,number);
    TLegend *legendPi0Spectrum                  = GetAndSetLegend(0.3,0.7,number);
    TLegend *legendPi0SpectrumRatio             = GetAndSetLegend(0.15,0.75,number);
    TLegend *legendPi0SpectrumFit               = GetAndSetLegend(0.3,0.7,number);
    TLegend *legendPi0SpectrumFitRatio          = GetAndSetLegend(0.15,0.75,number);
    TLegend *legendPurity                       = GetAndSetLegend(0.2,0.3,number);
    TLegend *legendPurityRatio                  = GetAndSetLegend(0.15,0.15,number);
    TLegend *legendGammaEff                     = GetAndSetLegend(0.3,0.15,number);
    TLegend *legendGammaEffRatio                = GetAndSetLegend(0.15,0.75,number);
    TLegend *legendGammaEffRecPt                = GetAndSetLegend(0.3,0.15,number);
    TLegend *legendGammaEffRecPtRatio           = GetAndSetLegend(0.15,0.75,number);
    TLegend *legendGammaConvProb                = GetAndSetLegend(0.3,0.15,number);
    TLegend *legendGammaConvProbRatio           = GetAndSetLegend(0.15,0.75,number);
    TLegend *legendGammaResolCorr               = GetAndSetLegend(0.15,0.75,number);
    TLegend *legendGammaResolCorrRatio          = GetAndSetLegend(0.15,0.75,number);
    TLegend *legendGammaCorrFac                 = GetAndSetLegend(0.15,0.75,number);
    TLegend *legendGammaCorrFacRatio            = GetAndSetLegend(0.15,0.75,number);
    
    TLegend *legendIncGammaToPi0Ratio           = GetAndSetLegend(0.3,0.7,number);
    TLegend *legendIncGammaToPi0RatioRatio      = GetAndSetLegend(0.15,0.75,number);
    TLegend *legendIncGammaToPi0RatioFit        = GetAndSetLegend(0.3,0.7,number);
    TLegend *legendIncGammaToPi0RatioFitRatio   = GetAndSetLegend(0.15,0.75,number);
    TLegend *legendDR                           = GetAndSetLegend(0.3,0.7,number);
    TLegend *legendDRRatio                      = GetAndSetLegend(0.15,0.75,number);
    TLegend *legendDRFit                        = GetAndSetLegend(0.3,0.7,number);
    TLegend *legendDRFitRatio                   = GetAndSetLegend(0.15,0.75,number);
    TLegend *legendRawGamma                     = GetAndSetLegend(0.45,0.75,number);
    TLegend *legendRawGammaRatio                = GetAndSetLegend(0.15,0.15,number);

    //********************************************************************************************
    //***************************** Reading histos from file *************************************
    //********************************************************************************************
    for(Int_t i = 0; i<number; i++){    
        TString fEventCutSelection;
        TString fGammaCutSelection;
        TString fElectronCutSelection;
        TString fMesonCutSelection;    
        TString fClusterCutSelection;
        ReturnSeparatedCutNumberAdvanced(cutSelection[i].Data(),fEventCutSelection, fGammaCutSelection, fClusterCutSelection, fElectronCutSelection, fMesonCutSelection, mode);

        cout << "entered reading" << endl;
        //Putting correct labels
        if (cutVariationName.Contains("SpecialTrigg")){
            TString fTrigger = fEventCutSelection(GetEventSelectSpecialTriggerCutPosition(),1);
            cutStringsName[i] = AnalyseSpecialTriggerCut(fTrigger.Atoi());      
        } else if (cutVariationName.Contains("MultiplicityPP")){
            TString minMult                                     = fEventCutSelection(GetEventCentralityMinCutPosition(),1);
            TString maxMult                                     = fEventCutSelection(GetEventCentralityMaxCutPosition(),1);
            cutStringsName[i]                                   = AnalysePPMultiplicityCut(minMult.Atoi(),maxMult.Atoi());      
        } else if (cutVariationName.Contains("V0Reader")){
            TString fV0Reader = fGammaCutSelection(GetPhotonV0FinderCutPosition(fGammaCutSelection),1);
            cutStringsName[i] = AnalyseV0ReaderCut(fV0Reader.Atoi());      
        } else if (cutVariationName.Contains("Eta")){
            TString fEtaCut = fGammaCutSelection(GetPhotonEtaCutPosition(fGammaCutSelection),1);    
            cout << fGammaCutSelection.Data() << "\t"<<fEtaCut.Data() << endl;
            cutStringsName[i] = AnalyseEtaCut(fEtaCut.Atoi());
        } else if (cutVariationName.Contains("RCut")){
            TString fRCut = fGammaCutSelection(GetPhotonMinRCutPosition(fGammaCutSelection),1);
            cutStringsName[i] = AnalyseRCut(fRCut.Atoi());
        } else if (cutVariationName.Contains("SinglePt")){
            TString fSinglePtCut = fGammaCutSelection(GetPhotonSinglePtCutPosition(fGammaCutSelection),1);
            cutStringsName[i] = AnalyseSinglePtCut(fSinglePtCut.Atoi());
        } else if (cutVariationName.Contains("Cluster")){
            TString fClusterCut = fGammaCutSelection(GetPhotonClsTPCCutPosition(fGammaCutSelection),1);
            cutStringsName[i] = AnalyseTPCClusterCut(fClusterCut.Atoi());
            cout << i << "\t" << fClusterCut.Data() << "\t" << cutStringsName[i].Data()<< endl;
        } else if (cutVariationName.Contains("dEdxE")){  
            TString fdEdxCut = fGammaCutSelection(GetPhotonEDedxSigmaCutPosition(fGammaCutSelection),1);
            cutStringsName[i] = AnalyseTPCdEdxCutElectronLine(fdEdxCut.Atoi());
        } else if (cutVariationName.Contains("dEdxPi")){    
            TString fdEdxCut = fGammaCutSelection(GetPhotonPiDedxSigmaCutPosition(fGammaCutSelection),3);
            cutStringsName[i] = AnalyseTPCdEdxCutPionLine(fdEdxCut.Data());      
        } else if (cutVariationName.Contains("Qt")){
            TString fQtCut = fGammaCutSelection(GetPhotonQtMaxCutPosition(fGammaCutSelection),1);
            cutStringsName[i] = AnalyseQtMaxCut(fQtCut.Atoi());
        } else if (cutVariationName.Contains("Chi2")){
            TString fChi2Cut = fGammaCutSelection(GetPhotonChi2GammaCutPosition(fGammaCutSelection),1);
            TString fPsiPairCut = fGammaCutSelection(GetPhotonPsiPairCutPosition(fGammaCutSelection),1);
            cutStringsName[i] = AnalyseChi2GammaCut(fChi2Cut.Atoi(),fPsiPairCut.Atoi());
        } else if (cutVariationName.Contains("PsiPairAndR")){
            TString fPsiPairCut = fGammaCutSelection(GetPhotonPsiPairCutPosition(fGammaCutSelection),1);
            TString fRCut = fGammaCutSelection(GetPhotonMinRCutPosition(fGammaCutSelection),1);
            cutStringsName[i] = AnalysePsiPairAndR(fPsiPairCut.Atoi(),fRCut.Atoi());      
        } else if (cutVariationName.Contains("PsiPair")){
            TString fPsiPairCut = fGammaCutSelection(GetPhotonPsiPairCutPosition(fGammaCutSelection),1);
            TString fChi2Cut = fGammaCutSelection(GetPhotonChi2GammaCutPosition(fGammaCutSelection),1);
            cutStringsName[i] = AnalysePsiPair(fPsiPairCut.Atoi(),fChi2Cut.Atoi());   
        } else if (cutVariationName.Contains("DCAZPhoton")){   
            TString fDCAZCut = fGammaCutSelection(GetPhotonDcaZPrimVtxCutPosition(fGammaCutSelection),1);
            cutStringsName[i] = AnalyseDCAZPhotonCut(fDCAZCut.Atoi());
        } else if (cutVariationName.Contains("CosPoint")){
            TString fCosPoint = fGammaCutSelection(GetPhotonCosinePointingAngleCutPosition(fGammaCutSelection),1);
            cutStringsName[i] = AnalyseCosPointCut(fCosPoint.Atoi());               
        } else if (cutVariationName.Contains("PhotonQuality")){
            TString fPhotonQuality = fGammaCutSelection(GetPhotonSharedElectronCutPosition(fGammaCutSelection),1);
            cutStringsName[i] = AnalysePhotonQuality(fPhotonQuality.Atoi());              
        } else if (cutVariationName.Contains("BG")){
            TString fBGCut = fMesonCutSelection(GetMesonBGSchemeCutPosition(),3);
            cutStringsName[i] = AnalyseBackgroundScheme(fBGCut.Data());   
        } else if (cutVariationName.Contains("Rapidity")){
            TString fRapidityCut = fMesonCutSelection(GetMesonRapidityCutPosition(),1);
            cutStringsName[i] = AnalyseRapidityMesonCut(fRapidityCut.Atoi());      
        } else if (cutVariationName.Contains("Alpha")){
            TString fAlphaCut = fMesonCutSelection(GetMesonAlphaCutPosition(),1);
            cutStringsName[i] = AnalyseAlphaMesonCut(fAlphaCut.Atoi());
        } else if (cutVariationName.Contains("Cent")){
            cutStringsName[i] = GetCentralityString(fGammaCutSelection.Data());
        } else if (cutVariationName.Contains("DiffRapWindow")){
            TString fRapidityCut = fMesonCutSelection(GetMesonRapidityCutPosition(),1);
            cutStringsName[i] = AnalyseRapidityMesonCutpPb(fRapidityCut.Atoi());      
        } else if (cutVariationName.Contains("MCSmearing")){
            TString fMCSmearing = fMesonCutSelection(GetMesonUseMCPSmearingCutPosition(),1);
            cutStringsName[i] = AnalyseMCSmearingCut(fMCSmearing.Atoi());      
        } else if (cutVariationName.Contains("IntRange")){
            if (i==0) cutStringsName[i] = "standard #pi^{0} integration range";
            if (i==1) cutStringsName[i] = "wide #pi^{0} integration range";
            if (i==2) cutStringsName[i] = "narrow #pi^{0} integration range";
        } else {
            cutStringsName[i] = cutSelection[i].Data();
        }

        
        nameCurrentCorrectedFile = Form("%s/%s/Gamma_Pi0_data_GammaConvV1_InclusiveRatio_0-100.root",cutSelection[i].Data(),energy.Data());
        fileCurrentFinal[i] = new TFile(nameCurrentCorrectedFile);
        if (fileCurrentFinal[i]->IsZombie()) haveOutputGammaToPi0=0;
        
        if (haveOutputGammaToPi0){
            histoIncGamma[i] = (TH1D*) fileCurrentFinal[i]->Get("histoGammaSpecCorrPurity");
            histoIncGamma[i]->SetTitle("");
            DrawGammaSetMarker(histoIncGamma[i], markerType, 2.0, color[i], color[i]);
            histoIncGammaRatio[i] = (TH1D*) histoIncGamma[i]->Clone(Form("histoIncGammaRatio_%s/%s",cutSelection[i].Data(),cutSelection[0].Data()));
            histoIncGammaRatio[i]->Divide(histoIncGammaRatio[i],histoIncGamma[0],1,1,"b");
            SetHistogramm(histoIncGammaRatio[i],"#it{p}_{T} (GeV/c)","Ratios of #gamma Spectra",0.0,2);
            
            histoIncGammaToPi0Ratio[i] = (TH1D*) fileCurrentFinal[i]->Get("IncRatioPurity_trueEff");
            if(i == 1 && cutVariationName.Contains("IntRange")) histoIncGammaToPi0Ratio[i] = (TH1D*) fileCurrentFinal[i]->Get("IncRatioPurity_trueEffWide");
            if(i == 2 && cutVariationName.Contains("IntRange")) histoIncGammaToPi0Ratio[i] = (TH1D*) fileCurrentFinal[i]->Get("IncRatioPurity_trueEffNarrow");
            histoIncGammaToPi0Ratio[i]->SetTitle("");
            DrawGammaSetMarker(histoIncGammaToPi0Ratio[i], markerType, 2.0, color[i], color[i]);
            histoIncGammaToPi0RatioRatio[i] = (TH1D*) histoIncGammaToPi0Ratio[i]->Clone(Form("histoIncGammaToPi0RatioRatio_%s/%s",cutSelection[i].Data(),cutSelection[0].Data()));
            histoIncGammaToPi0RatioRatio[i]->Divide(histoIncGammaToPi0RatioRatio[i],histoIncGammaToPi0Ratio[0],1,1,"b");
            SetHistogramm(histoIncGammaToPi0RatioRatio[i],"#it{p}_{T} (GeV/c)","Ratios of #gamma/#pi^{0} Ratios",0,2);
            
            histoPi0Spectrum[i] = (TH1D*) fileCurrentFinal[i]->Get("CorrectedYieldTrueEff");
            if(i == 1 && cutVariationName.Contains("IntRange")) histoPi0Spectrum[i] = (TH1D*) fileCurrentFinal[i]->Get("CorrectedYieldTrueEffWide");
            if(i == 2 && cutVariationName.Contains("IntRange")) histoPi0Spectrum[i] = (TH1D*) fileCurrentFinal[i]->Get("CorrectedYieldTrueEffNarrow");
            histoPi0Spectrum[i]->SetTitle("");
            DrawGammaSetMarker(histoPi0Spectrum[i], markerType, 2.0, color[i], color[i]);
            histoPi0SpectrumRatio[i] = (TH1D*) histoPi0Spectrum[i]->Clone(Form("histoPi0SpectrumRatio_%s/%s",cutSelection[i].Data(),cutSelection[0].Data()));
            histoPi0SpectrumRatio[i]->Divide(histoPi0SpectrumRatio[i],histoPi0Spectrum[0],1,1,"b");
            SetHistogramm(histoPi0SpectrumRatio[i],"#it{p}_{T} (GeV/c)","Ratios of #pi^{0} Spectra",0.0,2.0);

            histoPi0SpectrumFit[i] = (TH1D*) fileCurrentFinal[i]->Get("CorrectedYieldTrueEffPi0Fit");
            if(i == 1 && cutVariationName.Contains("IntRange")) histoPi0SpectrumFit[i] = (TH1D*) fileCurrentFinal[i]->Get("CorrectedYieldTrueEffPi0FitWide");
            if(i == 2 && cutVariationName.Contains("IntRange")) histoPi0SpectrumFit[i] = (TH1D*) fileCurrentFinal[i]->Get("CorrectedYieldTrueEffPi0FitNarrow");
            histoPi0SpectrumFit[i]->SetTitle("");
            DrawGammaSetMarker(histoPi0SpectrumFit[i], markerType, 2.0, color[i], color[i]);
            histoPi0SpectrumFitRatio[i] = (TH1D*) histoPi0SpectrumFit[i]->Clone(Form("histoPi0SpectrumFitRatio_%s/%s",cutSelection[i].Data(),cutSelection[0].Data()));
            histoPi0SpectrumFitRatio[i]->Divide(histoPi0SpectrumFitRatio[i],histoPi0SpectrumFit[0],1,1,"b");
            SetHistogramm(histoPi0SpectrumFitRatio[i],"#it{p}_{T} (GeV/c)","Ratios of #pi^{0} Spectra",0.0,2.0);
            
            histoIncGammaToPi0Ratio[i] = (TH1D*) fileCurrentFinal[i]->Get("IncRatioPurity_trueEff");
            if(i == 1 && cutVariationName.Contains("IntRange")) histoIncGammaToPi0Ratio[i] = (TH1D*) fileCurrentFinal[i]->Get("IncRatioPurity_trueEffWide");
            if(i == 2 && cutVariationName.Contains("IntRange")) histoIncGammaToPi0Ratio[i] = (TH1D*) fileCurrentFinal[i]->Get("IncRatioPurity_trueEffNarrow");
            histoIncGammaToPi0Ratio[i]->SetTitle("");
            DrawGammaSetMarker(histoIncGammaToPi0Ratio[i], markerType, 2.0, color[i], color[i]);
            histoIncGammaToPi0RatioRatio[i] = (TH1D*) histoIncGammaToPi0Ratio[i]->Clone(Form("histoIncGammaToPi0RatioRatio_%s/%s",cutSelection[i].Data(),cutSelection[0].Data()));
            histoIncGammaToPi0RatioRatio[i]->Divide(histoIncGammaToPi0RatioRatio[i],histoIncGammaToPi0Ratio[0],1,1,"b");
            SetHistogramm(histoIncGammaToPi0RatioRatio[i],"#it{p}_{T} (GeV/c)","Ratios of #gamma/#pi^{0} Ratios",0,2);

            histoIncGammaToPi0RatioFit[i] = (TH1D*) fileCurrentFinal[i]->Get("histoIncRatioFitPurity");
            if(i == 1 && cutVariationName.Contains("IntRange")) histoIncGammaToPi0RatioFit[i] = (TH1D*) fileCurrentFinal[i]->Get("histoIncRatioFitPurityWide");
            if(i == 2 && cutVariationName.Contains("IntRange")) histoIncGammaToPi0RatioFit[i] = (TH1D*) fileCurrentFinal[i]->Get("histoIncRatioFitPurityNarrow");
            if(i == 1 && cutVariationName.Contains("Fit")) histoIncGammaToPi0RatioFit[i] = (TH1D*) fileCurrentFinal[i]->Get("histoIncRatioLowFitPurity");
            if(i == 2 && cutVariationName.Contains("Fit")) histoIncGammaToPi0RatioFit[i] = (TH1D*) fileCurrentFinal[i]->Get("histoIncRatioHighFitPurity");

            histoIncGammaToPi0RatioFit[i]->SetTitle("");
            DrawGammaSetMarker(histoIncGammaToPi0RatioFit[i], markerType, 2.0, color[i], color[i]);
            histoIncGammaToPi0RatioFitRatio[i] = (TH1D*) histoIncGammaToPi0RatioFit[i]->Clone(Form("histoIncGammaToPi0RatioFitRatio_%s/%s",cutSelection[i].Data(),cutSelection[0].Data()));
            histoIncGammaToPi0RatioFitRatio[i]->Divide(histoIncGammaToPi0RatioFitRatio[i],histoIncGammaToPi0RatioFit[0],1,1,"b");
            SetHistogramm(histoIncGammaToPi0RatioFitRatio[i],"#it{p}_{T} (GeV/c)","Ratios of #gamma/#pi^{0}_{Fit} Ratios",0,2);

            histoDR[i] = (TH1D*) fileCurrentFinal[i]->Get("DoubleRatioConversionTrueEffPurity");
            if(i == 1 && cutVariationName.Contains("IntRange")) histoDR[i] = (TH1D*) fileCurrentFinal[i]->Get("DoubleRatioConversionTrueEffPurityWide");
            if(i == 2 && cutVariationName.Contains("IntRange")) histoDR[i] = (TH1D*) fileCurrentFinal[i]->Get("DoubleRatioConversionTrueEffPurityNarrow");
            if(i == 1 && cutVariationName.CompareTo("CocktailEta") == 0) histoDR[i] = (TH1D*) fileCurrentFinal[i]->Get("DoubleRatioConversionTrueEffPurityHigh");
            if(i == 2 && cutVariationName.CompareTo("CocktailEta") == 0) histoDR[i] = (TH1D*) fileCurrentFinal[i]->Get("DoubleRatioConversionTrueEffPurityLow");
            if(i == 1 && cutVariationName.CompareTo("CocktailEtaNorm") == 0) histoDR[i] = (TH1D*) fileCurrentFinal[i]->Get("DoubleRatioConversionTrueEffPurityEtaHigh");
            if(i == 2 && cutVariationName.CompareTo("CocktailEtaNorm") == 0) histoDR[i] = (TH1D*) fileCurrentFinal[i]->Get("DoubleRatioConversionTrueEffPurityEtaLow");
            if(i == 1 && cutVariationName.Contains("CocktailParam")) histoDR[i] = (TH1D*) fileCurrentFinal[i]->Get("DoubleRatioConversionTrueEffPurityModA");
            if(i == 2 && cutVariationName.Contains("CocktailParam")) histoDR[i] = (TH1D*) fileCurrentFinal[i]->Get("DoubleRatioConversionTrueEffPurityModB");
            histoDR[i]->SetTitle("");
            DrawGammaSetMarker(histoDR[i], markerType, 2.0, color[i], color[i]);
            histoDRRatio[i] = (TH1D*) histoDR[i]->Clone(Form("histoDRRatio_%s/%s",cutSelection[i].Data(),cutSelection[0].Data()));
            histoDRRatio[i]->Divide(histoDRRatio[i],histoDR[0],1,1,"b");
            SetHistogramm(histoDRRatio[i],"#it{p}_{T} (GeV/c)","Ratios of Double Ratios",0,2);

            histoDRFit[i] = (TH1D*) fileCurrentFinal[i]->Get("DoubleRatioConversionFitPurity");
            if(i == 1 && cutVariationName.Contains("IntRange")) histoDRFit[i] = (TH1D*) fileCurrentFinal[i]->Get("DoubleRatioConversionFitPurityWide");
            if(i == 2 && cutVariationName.Contains("IntRange")) histoDRFit[i] = (TH1D*) fileCurrentFinal[i]->Get("DoubleRatioConversionFitPurityNarrow");
            if(i == 1 && cutVariationName.Contains("Fit")) histoDRFit[i] = (TH1D*) fileCurrentFinal[i]->Get("DoubleRatioConversionLowFitPurity");
            if(i == 2 && cutVariationName.Contains("Fit")) histoDRFit[i] = (TH1D*) fileCurrentFinal[i]->Get("DoubleRatioConversionHighFitPurity");
            if(i == 1 && cutVariationName.CompareTo("CocktailEta") == 0) histoDRFit[i] = (TH1D*) fileCurrentFinal[i]->Get("DoubleRatioConversionFitPurityHigh");
            if(i == 2 && cutVariationName.CompareTo("CocktailEta") == 0) histoDRFit[i] = (TH1D*) fileCurrentFinal[i]->Get("DoubleRatioConversionFitPurityLow");
            if(i == 1 && cutVariationName.CompareTo("CocktailEtaNorm") == 0) histoDRFit[i] = (TH1D*) fileCurrentFinal[i]->Get("DoubleRatioConversionFitPurityEtaHigh");
            if(i == 2 && cutVariationName.CompareTo("CocktailEtaNorm") == 0) histoDRFit[i] = (TH1D*) fileCurrentFinal[i]->Get("DoubleRatioConversionFitPurityEtaLow");
            if(i == 1 && cutVariationName.Contains("CocktailParam")) histoDRFit[i] = (TH1D*) fileCurrentFinal[i]->Get("DoubleRatioConversionFitPurityModA");
            if(i == 2 && cutVariationName.Contains("CocktailParam")) histoDRFit[i] = (TH1D*) fileCurrentFinal[i]->Get("DoubleRatioConversionFitPurityModB");
            histoDRFit[i]->SetTitle("");
            DrawGammaSetMarker(histoDRFit[i], markerType, 2.0, color[i], color[i]);
            histoDRFitRatio[i] = (TH1D*) histoDRFit[i]->Clone(Form("histoDRFitRatio_%s/%s",cutSelection[i].Data(),cutSelection[0].Data()));
            histoDRFitRatio[i]->Divide(histoDRFitRatio[i],histoDRFit[0],1,1,"b");
            SetHistogramm(histoDRFitRatio[i],"#it{p}_{T} (GeV/c)","Ratios of Double Ratios",0,2);

            PlotCanvas(i,number,canvasGammaSpectrum,histoIncGamma[i],legendGammaSpectrum,canvasGammaSpectrumRatio,histoIncGammaRatio[i],legendGammaSpectrumRatio,cutStringsName[i],One);
            PlotCanvas(i,number,canvasPi0Spectrum,histoPi0Spectrum[i],legendPi0Spectrum,canvasPi0SpectrumRatio,histoPi0SpectrumRatio[i],legendPi0SpectrumRatio,cutStringsName[i],One);
            PlotCanvas(i,number,canvasPi0SpectrumFit,histoPi0SpectrumFit[i],legendPi0SpectrumFit,canvasPi0SpectrumFitRatio,histoPi0SpectrumFitRatio[i],legendPi0SpectrumFitRatio,cutStringsName[i],One);
            PlotCanvas(i,number,canvasIncGammaToPi0Ratio,histoIncGammaToPi0Ratio[i],legendIncGammaToPi0Ratio,canvasIncGammaToPi0RatioRatio,histoIncGammaToPi0RatioRatio[i],
                       legendIncGammaToPi0RatioRatio,cutStringsName[i],One);
            PlotCanvas(i,number,canvasIncGammaToPi0RatioFit,histoIncGammaToPi0RatioFit[i],legendIncGammaToPi0RatioFit,canvasIncGammaToPi0RatioFitRatio,histoIncGammaToPi0RatioFitRatio[i],
                       legendIncGammaToPi0RatioFitRatio,cutStringsName[i],One);
            PlotCanvas(i,number,canvasDR,histoDR[i],legendDR,canvasDRRatio,histoDRRatio[i],legendDRRatio,cutStringsName[i],One,One);
            PlotCanvas(i,number,canvasDRFit,histoDRFit[i],legendDRFit,canvasDRFitRatio,histoDRFitRatio[i],legendDRFitRatio,cutStringsName[i],One,One);
        } else {
            cout << "WARNING: I had to shrink the output to basics, the output of CalculateGammaToPi0 was missing" << endl;
        }   
        
        nameCurrentCorrectedFile = Form("%s/%s/Gamma_Pi0_data_GammaConvV1Correction_%s.root",cutSelection[i].Data(),energy.Data(),cutSelection[i].Data());
        fileCurrentCorrection[i] = new TFile(nameCurrentCorrectedFile);
        if (fileCurrentCorrection[i]->IsZombie()){
            cout << "ERROR: you are missing even the basic files" << endl;
            return;
        }
        
        histoRawGamma[i] = (TH1D*) fileCurrentCorrection[i]->Get("GammaRaw_Pt");
        histoRawGamma[i]->SetTitle("");
        DrawGammaSetMarker(histoRawGamma[i], markerType, 2.0, color[i], color[i]);
        histoRawGammaRatio[i] = (TH1D*) histoRawGamma[i]->Clone(Form("histoRawGammaRatio_%s/%s",cutSelection[i].Data(),cutSelection[0].Data()));
        histoRawGammaRatio[i]->Divide(histoRawGammaRatio[i],histoRawGamma[0],1,1,"b");
        SetHistogramm(histoRawGammaRatio[i],"#it{p}_{T} (GeV/c)","Ratios of Raw #gamma",0,2);

        histoPurity[i] = (TH1D*) fileCurrentCorrection[i]->Get("GammaPurityWOSec_Pt");
        histoPurity[i]->SetTitle("");
        DrawGammaSetMarker(histoPurity[i], markerType, 2.0, color[i], color[i]);
        histoPurityRatio[i] = (TH1D*) histoPurity[i]->Clone(Form("histoPurityRatio_%s/%s",cutSelection[i].Data(),cutSelection[0].Data()));
        histoPurityRatio[i]->Divide(histoPurityRatio[i],histoPurity[0],1,1,"b");
        SetHistogramm(histoPurityRatio[i],"#it{p}_{T} (GeV/c)","Ratios of Purities",0.8,1.2);

        histoGammaConvProb[i] = (TH1D*) fileCurrentCorrection[i]->Get("GammaConvProb_Pt");
        histoGammaConvProb[i]->SetTitle("");
        DrawGammaSetMarker(histoGammaConvProb[i], markerType, 2.0, color[i], color[i]);
        histoGammaConvProbRatio[i] = (TH1D*) histoGammaConvProb[i]->Clone(Form("histoConvProbRatio_%s/%s",cutSelection[i].Data(),cutSelection[0].Data()));
        histoGammaConvProbRatio[i]->Divide(histoGammaConvProbRatio[i],histoGammaConvProb[0],1,1,"b");
        SetHistogramm(histoGammaConvProbRatio[i],"#it{p}_{T} (GeV/c)","Ratios of Purities",0.8,1.2);
        
        histoGammaEff[i] = (TH1D*) fileCurrentCorrection[i]->Get("GammaRecoEff_MCPt");
        histoGammaEff[i]->SetTitle("");
        histoGammaEff[i]->GetYaxis()->SetRangeUser(0,1);
        SetHistogramm(histoGammaEff[i],"p_{T,MC} (GeV/c)","#epsilon_{reco,#gamma}",0.,1.);
        DrawGammaSetMarker(histoGammaEff[i], markerType, 2.0, color[i], color[i]);
        histoGammaEffRatio[i] = (TH1D*) histoGammaEff[i]->Clone(Form("histoGammaEffRatio_%s/%s",cutSelection[i].Data(),cutSelection[0].Data()));
        histoGammaEffRatio[i]->Divide(histoGammaEffRatio[i],histoGammaEff[0],1,1,"b");
        SetHistogramm(histoGammaEffRatio[i],"p_{T,MC} (GeV/c)","Ratios of Efficiencies",0.5,1.8);

        histoGammaEffRecPt[i] = (TH1D*) fileCurrentCorrection[i]->Get("GammaRecoEff_WithResolCorr_Pt");
        histoGammaEffRecPt[i]->SetTitle("");
        histoGammaEffRecPt[i]->GetYaxis()->SetRangeUser(0,1);
        SetHistogramm(histoGammaEffRecPt[i],"#it{p}_{T} (GeV/c)","#epsilon_{reco,#gamma} inc resol corr",0.,1.);
        DrawGammaSetMarker(histoGammaEffRecPt[i], markerType, 2.0, color[i], color[i]);
        histoGammaEffRecPtRatio[i] = (TH1D*) histoGammaEffRecPt[i]->Clone(Form("histoGammaEffRecPtRatio_%s/%s",cutSelection[i].Data(),cutSelection[0].Data()));
        histoGammaEffRecPtRatio[i]->Divide(histoGammaEffRecPtRatio[i],histoGammaEffRecPt[0],1,1,"b");
        SetHistogramm(histoGammaEffRecPtRatio[i],"#it{p}_{T} (GeV/c)","Ratios of Efficiencies inc resol corr",0.5,1.8);

        
        histoGammaResolCorr[i] = (TH1D*) fileCurrentCorrection[i]->Get("GammaResolCorrUnfold_Pt");
        histoGammaResolCorr[i]->SetTitle("");
        histoGammaResolCorr[i]->GetXaxis()->SetRangeUser(0,histoGammaEff[0]->GetXaxis()->GetBinUpEdge(histoGammaEff[0]->GetNbinsX()));
        histoGammaResolCorr[i]->GetYaxis()->SetRangeUser(0,2);
        DrawGammaSetMarker(histoGammaResolCorr[i], markerType, 2.0, color[i], color[i]);
        histoGammaResolCorrRatio[i] = (TH1D*) histoGammaResolCorr[i]->Clone(Form("histoGammaResolCorrRatio_%s/%s",cutSelection[i].Data(),cutSelection[0].Data()));
        histoGammaResolCorrRatio[i]->Divide(histoGammaResolCorrRatio[i],histoGammaResolCorr[0],1,1,"b");
        SetHistogramm(histoGammaResolCorrRatio[i],"#it{p}_{T} (GeV/c)","Ratios of resolution corrections",0,2);

        histoGammaCorrFac[i] = (TH1D*) fileCurrentCorrection[i]->Get("GammaCorrFac_Pt");
        histoGammaCorrFac[i]->SetTitle("");
        SetHistogramm(histoGammaCorrFac[i],"#it{p}_{T} (GeV/c)","#epsilon_{reco,#gamma} (#it{p}_{T}) #it{P}_{conv} (#it{p}_{T,MC})",0.,0.1);
        histoGammaCorrFac[i]->GetXaxis()->SetRangeUser(0,histoGammaEff[0]->GetXaxis()->GetBinUpEdge(histoGammaEff[0]->GetNbinsX()));
        DrawGammaSetMarker(histoGammaCorrFac[i], markerType, 2.0, color[i], color[i]);
        histoGammaCorrFacRatio[i] = (TH1D*) histoGammaCorrFac[i]->Clone(Form("histoGammaCorrFacRatio_%s/%s",cutSelection[i].Data(),cutSelection[0].Data()));
        histoGammaCorrFacRatio[i]->Divide(histoGammaCorrFacRatio[i],histoGammaCorrFac[0],1,1,"b");
        SetHistogramm(histoGammaCorrFacRatio[i],"#it{p}_{T} (GeV/c)","Ratios of correction factors",0,2);

        
        PlotCanvas(i,number,canvasRawGamma,histoRawGamma[i],legendRawGamma,canvasRawGammaRatio,histoRawGammaRatio[i],legendRawGammaRatio,cutStringsName[i],One);
        PlotCanvas(i,number,canvasPurity,histoPurity[i],legendPurity,canvasPurityRatio,histoPurityRatio[i],legendPurityRatio,cutStringsName[i],One);
        PlotCanvas(i,number,canvasGammaConvProb,histoGammaConvProb[i],legendGammaConvProb,canvasGammaConvProbRatio,histoGammaConvProbRatio[i],legendGammaConvProbRatio,cutStringsName[i],One);
        PlotCanvas(i,number,canvasGammaEff,histoGammaEff[i],legendGammaEff,canvasGammaEffRatio,histoGammaEffRatio[i],legendGammaEffRatio,cutStringsName[i],One);
        PlotCanvas(i,number,canvasGammaEffRecPt,histoGammaEffRecPt[i],legendGammaEffRecPt,canvasGammaEffRecPtRatio,histoGammaEffRecPtRatio[i],legendGammaEffRecPtRatio,cutStringsName[i],One);
        PlotCanvas(i,number,canvasGammaResolCorr,histoGammaResolCorr[i],legendGammaResolCorr,canvasGammaResolCorrRatio,histoGammaResolCorrRatio[i],legendGammaResolCorrRatio,
                   cutStringsName[i],One);
        PlotCanvas(i,number,canvasGammaCorrFac,histoGammaCorrFac[i],legendGammaCorrFac,canvasGammaCorrFacRatio,histoGammaCorrFacRatio[i],legendGammaCorrFacRatio,
                   cutStringsName[i],One);
        
        if (!haveOutputGammaToPi0){
            histoIncGamma[i] = (TH1D*) fileCurrentCorrection[i]->Get("GammaCorrUnfold_Pt");
            histoIncGamma[i]->SetTitle("");
            DrawGammaSetMarker(histoIncGamma[i], markerType, 2.0, color[i], color[i]);
            histoIncGammaRatio[i] = (TH1D*) histoIncGamma[i]->Clone(Form("histoIncGammaRatio_%s/%s",cutSelection[i].Data(),cutSelection[0].Data()));
            histoIncGammaRatio[i]->Divide(histoIncGammaRatio[i],histoIncGamma[0],1,1,"b");
            SetHistogramm(histoIncGammaRatio[i],"#it{p}_{T} (GeV/c)","Ratios of #gamma Spectra",0.0,2);
            
            PlotCanvas(i,number,canvasGammaSpectrum,histoIncGamma[i],legendGammaSpectrum,canvasGammaSpectrumRatio,histoIncGammaRatio[i],legendGammaSpectrumRatio,cutStringsName[i],One);
        }    
        
    }

    if (haveOutputGammaToPi0){
        canvasPi0Spectrum->SaveAs(Form("%s/Pi0Spectra_%s.eps",outputDir.Data(), (GetCentralityString(cutSelection[0])).Data()));
        canvasPi0SpectrumRatio->SaveAs(Form("%s/Pi0SpectraRatios_%s.eps",outputDir.Data(), (GetCentralityString(cutSelection[0])).Data()));
        canvasPi0SpectrumFit->SaveAs(Form("%s/Pi0SpectraFit_%s.eps",outputDir.Data(), (GetCentralityString(cutSelection[0])).Data()));
        canvasPi0SpectrumFitRatio->SaveAs(Form("%s/Pi0SpectraFitRatios_%s.eps",outputDir.Data(), (GetCentralityString(cutSelection[0])).Data()));
        canvasIncGammaToPi0Ratio->SaveAs(Form("%s/InclusiveRatios_%s.eps",outputDir.Data(), (GetCentralityString(cutSelection[0])).Data()));
        canvasIncGammaToPi0RatioRatio->SaveAs(Form("%s/InclusiveRatiosRatios_%s.eps",outputDir.Data(), (GetCentralityString(cutSelection[0])).Data()));
        canvasIncGammaToPi0RatioFit->SaveAs(Form("%s/InclusiveRatiosFit_%s.eps",outputDir.Data(), (GetCentralityString(cutSelection[0])).Data()));
        canvasIncGammaToPi0RatioFitRatio->SaveAs(Form("%s/InclusiveRatiosFitRatios_%s.eps",outputDir.Data(), (GetCentralityString(cutSelection[0])).Data()));
        canvasDR->SaveAs(Form("%s/DoubleRatios_%s.eps",outputDir.Data(), (GetCentralityString(cutSelection[0])).Data()));
        canvasDRRatio->SaveAs(Form("%s/DoubleRatiosRatios_%s.eps",outputDir.Data(), (GetCentralityString(cutSelection[0])).Data()));
        canvasDRFit->SaveAs(Form("%s/DoubleRatiosFit_%s.eps",outputDir.Data(), (GetCentralityString(cutSelection[0])).Data()));
        canvasDRFitRatio->SaveAs(Form("%s/DoubleRatiosFitRatios_%s.eps",outputDir.Data(), (GetCentralityString(cutSelection[0])).Data()));
    }
    
    canvasGammaSpectrum->SaveAs(Form("%s/GammaSpectra_%s.eps",outputDir.Data(),(GetCentralityString(cutSelection[0])).Data()));
    canvasGammaSpectrumRatio->SaveAs(Form("%s/GammaSpectraRatios_%s.eps",outputDir.Data(),(GetCentralityString(cutSelection[0])).Data()));
    canvasRawGamma->SaveAs(Form("%s/GammaRaw_%s.eps",outputDir.Data(),(GetCentralityString(cutSelection[0])).Data()));
    canvasRawGammaRatio->SaveAs(Form("%s/GammaRawRatios_%s.eps",outputDir.Data(),(GetCentralityString(cutSelection[0])).Data()));    
    canvasPurity->SaveAs(Form("%s/GammaPurity_%s.eps",outputDir.Data(),(GetCentralityString(cutSelection[0])).Data()));
    canvasPurityRatio->SaveAs(Form("%s/GammaPurityRatios_%s.eps",outputDir.Data(),(GetCentralityString(cutSelection[0])).Data()));
    canvasGammaEff->SaveAs(Form("%s/GammaEff_%s.eps",outputDir.Data(),(GetCentralityString(cutSelection[0])).Data()));
    canvasGammaEffRatio->SaveAs(Form("%s/GammaEffRatios_%s.eps",outputDir.Data(),(GetCentralityString(cutSelection[0])).Data()));
    canvasGammaEffRecPt->SaveAs(Form("%s/GammaEffRecPt_%s.eps",outputDir.Data(),(GetCentralityString(cutSelection[0])).Data()));
    canvasGammaEffRecPtRatio->SaveAs(Form("%s/GammaEffRatiosRecPt_%s.eps",outputDir.Data(),(GetCentralityString(cutSelection[0])).Data()));
    canvasGammaResolCorr->SaveAs(Form("%s/GammaResolCorr_%s.eps",outputDir.Data(),(GetCentralityString(cutSelection[0])).Data()));
    canvasGammaResolCorrRatio->SaveAs(Form("%s/GammaResolCorrRatios_%s.eps",outputDir.Data(),(GetCentralityString(cutSelection[0])).Data()));
    canvasGammaConvProb->SaveAs(Form("%s/GammaConvProb_%s.eps",outputDir.Data(),(GetCentralityString(cutSelection[0])).Data()));
    canvasGammaConvProbRatio->SaveAs(Form("%s/GammaConvProbRatios_%s.eps",outputDir.Data(),(GetCentralityString(cutSelection[0])).Data()));
    canvasGammaCorrFac->SaveAs(Form("%s/GammaCorrFac_%s.eps",outputDir.Data(),(GetCentralityString(cutSelection[0])).Data()));
    canvasGammaCorrFacRatio->SaveAs(Form("%s/GammaCorrFacRatios_%s.eps",outputDir.Data(),(GetCentralityString(cutSelection[0])).Data()));
    
    CalculateSystematicsGraphs( histoIncGamma, histoRawGamma, cutSelection, number, "Gamma", "Gamma", outputDir, outputFileDir, cutVariationName);

    if (haveOutputGammaToPi0){
        //*************************************************************************************************
        //******************** Output of the systematic Error due to Signal extraction for DoubleRatio ****
        //*************************************************************************************************
        Int_t NumberOfCutsDoubleRatio = number;
        if  (   cutVariationName.CompareTo("CocktailEtaNorm") == 0 || 
                cutVariationName.CompareTo("CocktailEta") == 0 || 
                cutVariationName.Contains("CocktailParam") || 
                cutVariationName.CompareTo("Yield") == 0 || 
                cutVariationName.Contains("Fit") || 
                cutVariationName.CompareTo("Purity") == 0
            ){
            NumberOfCutsDoubleRatio = 3;
        }
        if(cutVariationName.CompareTo("Charged") == 0){
            NumberOfCutsDoubleRatio = 2;
        }
        CalculateSystematicsGraphsWOBound (histoDR, cutSelection, NumberOfCutsDoubleRatio, "DoubleRatio", "Gamma", outputDir, outputFileDir, cutVariationName);
        
        //*************************************************************************************************
        //***************** Output of the systematic Error due to Signal extraction for DoubleRatioFit ****
        //*************************************************************************************************
        CalculateSystematicsGraphsWOBound (histoDRFit, cutSelection, NumberOfCutsDoubleRatio, "DoubleRatioFit", "Gamma", outputDir, outputFileDir, cutVariationName);
        
        //*************************************************************************************************
        //******************** Output of the systematic Error due to Signal extraction for IncRatio ****
        //*************************************************************************************************
        Int_t NumberOfCutsIncRatio = number;
        if  (   cutVariationName.CompareTo("Yield") == 0 || 
                cutVariationName.Contains("Fit") || 
                cutVariationName.CompareTo("Purity") == 0
            ){
            NumberOfCutsIncRatio = 3;
        }
        CalculateSystematicsGraphsWOBound (histoIncGammaToPi0Ratio, cutSelection, NumberOfCutsIncRatio, "IncRatio", "Gamma", outputDir, outputFileDir, cutVariationName);

        //*************************************************************************************************
        //******************** Output of the systematic Error due to Signal extraction for IncRatioFit ****
        //*************************************************************************************************

        CalculateSystematicsGraphsWOBound (histoIncGammaToPi0RatioFit, cutSelection, NumberOfCutsIncRatio, "IncRatioFit", "Gamma", outputDir, outputFileDir, cutVariationName);
        
        //*************************************************************************************************
        //******************** Output of the systematic Error due to Signal extraction for Pi0 ****
        //*************************************************************************************************

        Int_t NumberOfCutsPi0 = number;
        if  (   cutVariationName.CompareTo("Yield") == 0 || 
                cutVariationName.Contains("Fit") || 
                cutVariationName.CompareTo("Purity") == 0
            ){
            NumberOfCutsPi0 = 3;
        }
        CalculateSystematicsGraphsWOBound (histoPi0Spectrum, cutSelection, NumberOfCutsPi0, "Pi0", "Gamma", outputDir, outputFileDir, cutVariationName);
        
        //*************************************************************************************************
        //******************** Output of the systematic Error due to Signal extraction for Pi0Fit ****
        //*************************************************************************************************

        CalculateSystematicsGraphsWOBound (histoPi0SpectrumFit, cutSelection, NumberOfCutsPi0, "Pi0Fit", "Gamma", outputDir, outputFileDir, cutVariationName);

    }
    return;
}

