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
#include "THnSparse.h"
#include "TParallelCoord.h"
#include "TParallelCoordVar.h"
#include "TTree.h"
#include "TLeaf.h"
#include "THnSparse.h"
#include "TAxis.h"
#include "TRandom.h"
#include "../CommonHeaders/ConversionFunctions.h"
#include "../CommonHeaders/PlottingGammaConversionHistos.h"
#include "../CommonHeaders/AdjustHistRange.h"
#include "../CommonHeaders/PlottingGammaConversionAdditional.h"
#include "../CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"

/*
- This macro reads a data root file from train running with THnSparse option switched on (needs to contain <cut> Back histograms and <cut> Mother histograms)
- It's purpose is to crosscheck the event mixing classes
- The THnSparse is a 4D object containing info about the inv mass of gamma pairs, the pT of the tracks, and the vertex Z and (photon/track) multiplicity of the event
- The macro produces invariant mass histograms in bins of pT, M and Z (In four versions: Same events, mixed events, same events normalized and mixed events normalized)
- It also cretes 2D histograms (z vs m) of the integrals of invariant mass plot (=> number of photon pairs). One 2D histo for each pT bin
- provided by Meike Danisch, meike.charlotte.danisch@cern.ch, July 2019, Gamma Conversion Group, PWGGA
 */

TString GetBinLabel(TString variable, Int_t bin, TString cent = "" ){

    TString binLabel = "";

    Double_t fBinLimitsArrayPT[6] = {0.,5.,10.,20.,50.,250.};  // here
    Double_t fBinLimitsArrayZPbPb[8] =  {-50, -5.5, -2.9, -0.65, 1.45, 3.65, 6.15, 50};
    Double_t fBinLimitsArrayMultiplicityPbPb0010[5] =   {1., 50., 65., 90., 200.};             
    Double_t fBinLimitsArrayMultiplicityPbPb1020[5] =   {1., 30., 40., 60., 200.};
    Double_t fBinLimitsArrayMultiplicityPbPb0020[5] =   {1., 30., 50., 70., 200.};
    Double_t fBinLimitsArrayMultiplicityPbPb2040[5] =   {1., 12., 20., 30., 200.};
    Double_t fBinLimitsArrayMultiplicityPbPb4060[5] =   {1., 4., 7., 13., 200.};
    Double_t fBinLimitsArrayMultiplicityPbPb6080[5] =   {1., 3., 4., 5., 200.};
    Double_t fBinLimitsArrayMultiplicityPbPb2050[5] =   {1., 8., 18., 30., 200.};

    if(variable.CompareTo("z")==0){
        binLabel = Form("%.1f < z < %.1f", fBinLimitsArrayZPbPb[bin],fBinLimitsArrayZPbPb[bin+1]);
    } else if(variable.CompareTo("m")==0){
        if(cent.CompareTo("0-10%")==0){
            binLabel = Form("%.1f < m < %.1f", fBinLimitsArrayMultiplicityPbPb0010[bin],fBinLimitsArrayMultiplicityPbPb0010[bin+1]);
        } else if(cent.CompareTo("10-20%")==0){
            binLabel = Form("%.1f < m < %.1f", fBinLimitsArrayMultiplicityPbPb1020[bin],fBinLimitsArrayMultiplicityPbPb1020[bin+1]);
        } else if(cent.CompareTo("0-20%")==0){
            binLabel = Form("%.1f < m < %.1f", fBinLimitsArrayMultiplicityPbPb0020[bin],fBinLimitsArrayMultiplicityPbPb0020[bin+1]);
        } else if(cent.CompareTo("20-40%")==0){
            binLabel = Form("%.1f < m < %.1f", fBinLimitsArrayMultiplicityPbPb2040[bin],fBinLimitsArrayMultiplicityPbPb2040[bin+1]);
        } else if(cent.CompareTo("40-60%")==0){
            binLabel = Form("%.1f < m < %.1f", fBinLimitsArrayMultiplicityPbPb4060[bin],fBinLimitsArrayMultiplicityPbPb4060[bin+1]);
        } else if(cent.CompareTo("60-80%")==0){
            binLabel = Form("%.1f < m < %.1f", fBinLimitsArrayMultiplicityPbPb6080[bin],fBinLimitsArrayMultiplicityPbPb6080[bin+1]);
        } else if(cent.CompareTo("20-50%")==0){
            binLabel = Form("%.1f < m < %.1f", fBinLimitsArrayMultiplicityPbPb2050[bin],fBinLimitsArrayMultiplicityPbPb2050[bin+1]);
        } else
            cout << "GetBinLabel: cent not defined correctly" << endl;
    } else if(variable.CompareTo("p")==0){
        binLabel = Form("%.1f < p < %.1f GeV/c", 0.1*fBinLimitsArrayPT[bin],0.1*fBinLimitsArrayPT[bin+1]);
    } else
        cout << "GetBinLabel: variable not defined correctly" << endl;
    return binLabel;

}

void DrawCompareHistoTH1(   TCanvas* canvas,Double_t leftMargin,Double_t rightMargin,Double_t topMargin,Double_t bottomMargin,Bool_t logX, Bool_t logY, Bool_t logZ,
                            std::vector<TH1D*> vecHist, TString title,TString xTitle,TString yTitle, Double_t xOffset, Double_t yOffset,
                            Color_t *color, Bool_t adjustRange, Double_t rangeMin, Double_t rangeMax, Bool_t useErrors,
                            Double_t p1,Double_t p2,Double_t p3, TString text1, TString* plotDataSets, TString text2, TString text3,
                            Int_t textAlign = 31)
{
    if((Int_t) vecHist.size()<1) {
        cout << "ERROR: DrawCompareHistoTH1, vector size <1! Returning...";
        return;
    }

    canvas->cd();
    canvas->SetLeftMargin(leftMargin);
    canvas->SetRightMargin(rightMargin);
    canvas->SetTopMargin(topMargin);
    canvas->SetBottomMargin(bottomMargin);
    canvas->SetLogx(logX);
    canvas->SetLogy(logY);
    canvas->SetLogz(logZ);

    TLegend *leg1 = new TLegend(0.2,0.92,0.99,0.99);
    leg1->SetNColumns(4);
    leg1->SetFillColor(0);
    leg1->SetLineColor(0);
    leg1->SetTextSize(0.03);
    leg1->SetTextFont(42);

    if(adjustRange) AdjustHistRange(vecHist,rangeMin,rangeMax,useErrors);

    Int_t iEnd = (Int_t) vecHist.size();
    for(Int_t i=0; i<iEnd; i++){
        TH1D* fHist = (TH1D*) vecHist.at(i);
        DrawGammaSetMarker(fHist, 1, 0.5, color[i], color[i]);
        leg1->AddEntry(fHist,plotDataSets[i].Data());
    }

    Int_t iStart = (Int_t) vecHist.size()-1;
    for(Int_t i=iStart; i>=0; i--){
        TH1D* fHist = (TH1D*) vecHist.at(i);
        //x-axis
        fHist->GetXaxis()->SetTitleOffset(xOffset);
        fHist->SetXTitle(xTitle.Data());
        fHist->GetXaxis()->SetLabelFont(42);
        fHist->GetXaxis()->SetTitleFont(62);
        fHist->GetXaxis()->SetTitleSize(0.04);
        fHist->GetXaxis()->SetLabelSize(0.035);
        if (logX){
            fHist->GetXaxis()->SetLabelOffset(-0.01);
        }
        //y-axis
        fHist->GetYaxis()->SetTitleOffset(yOffset);
        fHist->SetYTitle(yTitle.Data());
        fHist->GetYaxis()->SetLabelFont(42);
        fHist->GetYaxis()->SetTitleFont(62);
        fHist->GetYaxis()->SetTitleSize(0.04);
        fHist->GetYaxis()->SetLabelSize(0.035);
        fHist->GetYaxis()->SetDecimals();

        fHist->SetTitle(title.Data());

        fHist->SetLineWidth(3);
        
        if(i==iStart) fHist->DrawCopy("hist");
        else fHist->DrawCopy("hist,same");
    }

    leg1->Draw();

    PutProcessLabelAndEnergyOnPlot(p1,p2,p3, text1.Data(), text2.Data(), text3.Data(), 42, 0.03, "", 1, 1.25, textAlign);

    return;
}

void AnalyseBackgroundFromTHnSparse(TString fileData, Int_t mode, TString cutString, TString suffix, TString meson = "Pi0"){

    // read THnSparse objects from fileData

    TFile* f = new TFile(fileData.Data());
    TString autoDetectedMainDir = "";
    autoDetectedMainDir = AutoDetectMainTList(mode,f,"","");
    if (autoDetectedMainDir.CompareTo("") == 0){
        cout << "ERROR: trying to read file, which is incompatible with mode selected (mode " << mode << ")" << endl;;
        return;
    }
    TList *TopDir                 = (TList*)f->Get(autoDetectedMainDir.Data());
    if(TopDir == NULL){
        cout<<"ERROR: TopDir not Found"<<endl;
        return;
    }
    cout << "Reading from TopDir \"" << autoDetectedMainDir << "\"" << endl;

    TList *HistosGammaConversion = (TList*)TopDir->FindObject(Form("Cut Number %s",cutString.Data()));
    TList *BackgroundContainer   = (TList*) HistosGammaConversion->FindObject(Form("%s Back histograms",cutString.Data()));
    TList *MotherContainer       = (TList*) HistosGammaConversion->FindObject(Form("%s Mother histograms",cutString.Data()));
    THnSparseF* fSparseMotherZM  = NULL;
    THnSparseF* fSparseBckZM     = NULL;
    fSparseMotherZM              = (THnSparseF*)MotherContainer->FindObject("Back_Mother_InvMass_Pt_z_m");
    fSparseBckZM                 = (THnSparseF*)BackgroundContainer->FindObject("Back_Back_InvMass_Pt_z_m");
    if(fSparseMotherZM == NULL) {
        cout << "ERROR: could not get Back_Mother_InvMass_Pt_z_m from Mother histograms" << endl;
        return;
    }
    if(fSparseBckZM == NULL){
        cout << "ERROR: could not get Back_Back_InvMass_Pt_z_m from Back histograms" << endl;
        return;
    }

    // create histograms integrated over all bins

    TH1D* fHistoInvMassMotherFull;
    TH1D* fHistoInvMassBckFull;

    fHistoInvMassMotherFull = (TH1D*)fSparseMotherZM->Projection(0);
    fHistoInvMassMotherFull->SetName("fHistoInvMassMother_Full");
    fHistoInvMassBckFull    = (TH1D*)fSparseBckZM->Projection(0);
    fHistoInvMassBckFull->SetName("fHistoInvMassBckFull_Full");

    // plotting settings and globally used objects

    StyleSettingsThesis(suffix);
    SetPlotStyle();
    TString outputDir   = Form("AnalyseBck_%s_%s",meson.Data(), cutString.Data());
    gSystem->Exec("mkdir -p "+outputDir);
    cout<<"Pictures are saved as "<< suffix.Data()<< endl;

    TCanvas *canvasInvMass = new TCanvas("canvasInvMass","",1550,1200);    
    DrawGammaCanvasSettings( canvasInvMass, 0.1, 0.05, 0.05, 0.10);
    TLegend* legendInvMass = GetAndSetLegend2(0.4,0.15,0.8,0.25, 0.04*1200,1);

    TCanvas *canvas2D = new TCanvas("canvas2D","",1700,1500);    
    DrawGammaCanvasSettings( canvas2D, 0.16, 0.11, 0.06, 0.1);  // left right top bottom

    std::vector<TH1D*> vecInvMass;

    // get centrality
    TString centCut = GetCentralityString(cutString);

    // define x range around meson mass
    Double_t xMin = 0.05;
    Double_t xMax = 0.2;
    if (meson.CompareTo("Eta")==0){
        xMin = 0.4;
        xMax = 0.7;
    }

    // plot integrated over all bins, compare background to mother histo
    canvasInvMass->cd();
    canvasInvMass->SetLogy();
    DrawGammaSetMarker(fHistoInvMassMotherFull, 20, 1.0, kRed, kRed); 
    DrawGammaSetMarker(fHistoInvMassBckFull, 24, 1.5, kBlack, kBlack);
    fHistoInvMassMotherFull->SetTitle("");
    fHistoInvMassBckFull->SetTitle("");

    // set x range
    Int_t xMinBin = fHistoInvMassMotherFull->FindBin(xMin);
    Int_t xMaxBin = fHistoInvMassMotherFull->FindBin(xMax);
    fHistoInvMassMotherFull->GetXaxis()->SetRange(xMinBin,xMaxBin);
    fHistoInvMassBckFull->GetXaxis()->SetRange(xMinBin,xMaxBin);

    fHistoInvMassBckFull->Draw("p,same");
    fHistoInvMassMotherFull->Draw("p,same");

    // set y range
    vecInvMass.push_back(fHistoInvMassBckFull);
    vecInvMass.push_back(fHistoInvMassMotherFull);
    AdjustHistRange(vecInvMass, 5, 5, kFALSE, 0, 0, kTRUE, xMinBin, xMaxBin);

    legendInvMass->AddEntry(fHistoInvMassMotherFull, "same events","p");
    legendInvMass->AddEntry(fHistoInvMassBckFull, "mixed events","p");
    legendInvMass->Draw();
    canvasInvMass->Update();
    canvasInvMass->SaveAs(Form("%s/InvMassProjection.%s",outputDir.Data(),suffix.Data()));

    // plot integrated over all bins, compare background to mother histo, normalized
    canvasInvMass->SetLogy(0);

    // normalize
    fHistoInvMassBckFull->Scale(1./fHistoInvMassBckFull->Integral());
    fHistoInvMassMotherFull->Scale(1./fHistoInvMassMotherFull->Integral());

    // set y range (should be the same or both histos due to normalization)
    AdjustHistRange(fHistoInvMassMotherFull, 1.1, 1.1, kFALSE, 0, 0, xMinBin, xMaxBin);
    AdjustHistRange(fHistoInvMassBckFull, 1.1, 1.1, kFALSE, 0, 0, xMinBin, xMaxBin);

    fHistoInvMassBckFull->Draw("p,same");
    fHistoInvMassMotherFull->Draw("p,same");
    legendInvMass->Draw();
    canvasInvMass->Update();
    canvasInvMass->SaveAs(Form("%s/InvMassProjectionNormalized.%s",outputDir.Data(),suffix.Data()));

    // create histograms in bins

    const Int_t nBinsPt     = 250;
    const Int_t nBinsPtProj = 5;
    const Int_t nBinsZ      = 7;
    const Int_t nBinsM      = 4;

    Double_t fBinLimitsArrayPT[6] = {0.,5.,10.,20.,50.,250.};  // here        bin numbers correspond to 0. - 0.5 - 1 - 2 - 5 - 25 GeV

    if(fSparseMotherZM->GetAxis(1)->GetNbins() != nBinsPt || fSparseMotherZM->GetAxis(2)->GetNbins() != nBinsZ || fSparseMotherZM->GetAxis(3)->GetNbins() != nBinsM){
        cout << "wrong number of bins" << endl;
        return;
    }

    // bin value labels

    TString* plotLabelZ = new TString[nBinsZ];
    for(Int_t z=0;z<nBinsZ;z++){
        plotLabelZ[z] =  GetBinLabel("z", z);  
    }
    TString* plotLabelM = new TString[nBinsM];
    for(Int_t m=0;m<nBinsM;m++){
        plotLabelM[m] =  GetBinLabel("m", m, centCut);  
    }
    TString* plotLabelP = new TString[nBinsPtProj]; 
    for(Int_t p=0;p<nBinsPtProj;p++){
        plotLabelP[p] =  GetBinLabel("p", p);  
    }

    TH1D* fHistoInvMassMother[nBinsPt][nBinsZ][nBinsM];
    TH1D* fHistoInvMassBck[nBinsPt][nBinsZ][nBinsM];

    TH2D* fHistoMotherIntegralMapZM[nBinsPt];
    TH2D* fHistoBckIntegralMapZM[nBinsPt];
    for(Int_t p=0; p<nBinsPtProj; p++) {
        fHistoMotherIntegralMapZM[p] = new TH2D("fHistoMotherIntegralMapZM",Form("# of #gamma pairs from same events for %s",plotLabelP[p].Data()),nBinsZ,0,nBinsZ,nBinsM,0,nBinsM); // bin 1 from 0 to 1 etc.
        fHistoBckIntegralMapZM[p] = new TH2D("fHistoBckIntegralMapZM",Form("# of #gamma pairs from mixed events for %s",plotLabelP[p].Data()),nBinsZ,0,nBinsZ,nBinsM,0,nBinsM);
    }

    cout << "Projection for ..." << endl;

    for(Int_t m=0; m<nBinsM; m++) {

        for(Int_t z=0;z<nBinsZ;z++) {

            for(Int_t p=0; p<nBinsPtProj; p++) {

                cout <<  "m = " << m << ", z = " << z << ", p = " << fBinLimitsArrayPT[p] << " - " << fBinLimitsArrayPT[p+1] << endl;   // here

                fSparseMotherZM->GetAxis(1)->SetRange(fBinLimitsArrayPT[p],fBinLimitsArrayPT[p+1]);  // here
                fSparseMotherZM->GetAxis(2)->SetRange(z+1,z+1);    // bins are counted starting from 1
                fSparseMotherZM->GetAxis(3)->SetRange(m+1,m+1);
                fHistoInvMassMother[p][z][m] = (TH1D*)fSparseMotherZM->Projection(0);
                fHistoInvMassMother[p][z][m]->SetName(Form("fHistoInvMassMother_Pt%d_Z%d_M%d",p,z,m));
                fHistoInvMassMother[p][z][m]->GetXaxis()->SetRange(xMinBin,xMaxBin);

                fSparseBckZM->GetAxis(1)->SetRange(fBinLimitsArrayPT[p],fBinLimitsArrayPT[p+1]);  // here
                fSparseBckZM->GetAxis(2)->SetRange(z+1,z+1);
                fSparseBckZM->GetAxis(3)->SetRange(m+1,m+1);
                fHistoInvMassBck[p][z][m] = (TH1D*)fSparseBckZM->Projection(0);
                fHistoInvMassBck[p][z][m]->SetName(Form("fHistoInvMassBck_Pt%d_Z%d_M%d",p,z,m));
                fHistoInvMassBck[p][z][m]->GetXaxis()->SetRange(xMinBin,xMaxBin);

                // fill 2D histo                
                fHistoMotherIntegralMapZM[p]->Fill(z,m,fHistoInvMassMother[p][z][m]->Integral());  // z = 0 => bin 1 [0-1] is filled, z = 3 => bin 4 [3-4] is filled
                fHistoBckIntegralMapZM[p]->Fill(z,m,fHistoInvMassBck[p][z][m]->Integral());      
            }
        }
    }

    //=========================================================================================  
    //========== plotting Z bins comparisons ==================================================
    //=========================================================================================

    cout << "Plotting z bin comparisons..." << endl;

    Color_t colorsZ[nBinsZ] = {kBlack, kBlue, kRed, kGreen-2, kMagenta+1, kOrange+7, kCyan-2};

    for(Int_t p=0; p<nBinsPtProj; p++) {

        for(Int_t m=0; m<nBinsM; m++) {

            // compare background between z bins
            canvasInvMass->Clear();
            vecInvMass.clear();
            for(Int_t z=0;z<nBinsZ;z++) {
                vecInvMass.push_back(fHistoInvMassBck[p][z][m]);
            }
            DrawCompareHistoTH1(canvasInvMass,0.14, 0.02, 0.12, 0.1,kFALSE,kFALSE,kFALSE,  // left right top bottom margin
                                vecInvMass,"","M_{#gamma#gamma}","#frac{dN}{dM_{#gamma#gamma}}",1,1.4,
                                colorsZ, kTRUE, 1.1, 1.1, kFALSE, 
                                0.17,0.85,0.03,plotLabelP[p],plotLabelZ,plotLabelM[m],"mixed events",11);  // x,y,height
            canvasInvMass->Update();
            canvasInvMass->SaveAs(Form("%s/InvMassBGProjection_p%d_z_m%d.%s",outputDir.Data(),p,m,suffix.Data()));

            // compare background between z bins, normalized
            canvasInvMass->Clear();
            for(Int_t iVec=0; iVec<(Int_t)vecInvMass.size(); iVec++){
                TH1D* temp = vecInvMass.at(iVec);
                temp->Sumw2();
                temp->Scale(1./temp->Integral());
            }
            DrawCompareHistoTH1(canvasInvMass,0.14, 0.02, 0.12, 0.1,kFALSE,kFALSE,kFALSE,
                                vecInvMass,"","M_{#gamma#gamma}","#frac{1}{N} #frac{dN}{dM_{#gamma#gamma}}",1,1.4,
                                colorsZ, kFALSE, 1, 1, kFALSE,
                                0.17,0.85,0.03,plotLabelP[p],plotLabelZ,plotLabelM[m],"mixed events",11);
            canvasInvMass->Update();
            canvasInvMass->SaveAs(Form("%s/InvMassBGProjectionNormalized_p%d_z_m%d.%s",outputDir.Data(),p,m,suffix.Data()));

            // compare mother histo between z bins
            canvasInvMass->Clear();
            vecInvMass.clear();
            for(Int_t z=0;z<nBinsZ;z++) {
                vecInvMass.push_back(fHistoInvMassMother[p][z][m]);
            }
            DrawCompareHistoTH1(canvasInvMass,0.14, 0.02, 0.12, 0.1,kFALSE,kFALSE,kFALSE,
                                vecInvMass,"","M_{#gamma#gamma}","#frac{dN}{dM_{#gamma#gamma}}",1,1.4,
                                colorsZ, kTRUE, 1.1, 1.1, kFALSE, 
                                0.17,0.85,0.03,plotLabelP[p],plotLabelZ,plotLabelM[m],"same events",11);
            canvasInvMass->Update();
            canvasInvMass->SaveAs(Form("%s/InvMassMotherProjection_p%d_z_m%d.%s",outputDir.Data(),p,m,suffix.Data()));

            // compare mother histo between z bins, normalized
            canvasInvMass->Clear();
            for(Int_t iVec=0; iVec<(Int_t)vecInvMass.size(); iVec++){
                TH1D* temp = vecInvMass.at(iVec);
                temp->Sumw2();
                temp->Scale(1./temp->Integral());
            }

            DrawCompareHistoTH1(canvasInvMass,0.14, 0.02, 0.12, 0.1,kFALSE,kFALSE,kFALSE,
                                vecInvMass,"","M_{#gamma#gamma}","#frac{1}{N} #frac{dN}{dM_{#gamma#gamma}}",1,1.4,
                                colorsZ, kFALSE, 1, 1, kFALSE,
                                0.17,0.85,0.03,plotLabelP[p],plotLabelZ,plotLabelM[m],"same events",11);
            canvasInvMass->Update();
            canvasInvMass->SaveAs(Form("%s/InvMassMotherProjectionNormalized_p%d_z_m%d.%s",outputDir.Data(),p,m,suffix.Data()));
        }
    }

    //=========================================================================================  
    //========== plotting M bins comparisons ==================================================
    //=========================================================================================

    cout << "Plotting m bin comparisons..." << endl;

    Color_t colorsM[nBinsM] = {kBlack, kRed, kGreen+2, kMagenta+1};

    for(Int_t p=0; p<nBinsPtProj; p++) {

        for(Int_t z=0; z<nBinsM; z++) {

            // compare background between m bins
            canvasInvMass->Clear();
            vecInvMass.clear();
            for(Int_t m=0;m<nBinsM;m++) {
                vecInvMass.push_back(fHistoInvMassBck[p][z][m]);
            }
            DrawCompareHistoTH1(canvasInvMass,0.14, 0.02, 0.12, 0.1,kFALSE,kFALSE,kFALSE,  // left right top bottom margin
                                vecInvMass,"","M_{#gamma#gamma}","#frac{dN}{dM_{#gamma#gamma}}",1,1.4,
                                colorsM, kTRUE, 1.1, 1.1, kFALSE, 
                                0.17,0.85,0.03,plotLabelP[p],plotLabelM,plotLabelZ[z],"mixed events",11);  // x,y,height
            canvasInvMass->Update();
            canvasInvMass->SaveAs(Form("%s/InvMassBGProjection_p%d_z%d_m.%s",outputDir.Data(),p,z,suffix.Data()));

            // compare background between m bins, normalized
            canvasInvMass->Clear();
            for(Int_t iVec=0; iVec<(Int_t)vecInvMass.size(); iVec++){
                TH1D* temp = vecInvMass.at(iVec);
                temp->Sumw2();
                temp->Scale(1./temp->Integral());
            }
            DrawCompareHistoTH1(canvasInvMass,0.14, 0.02, 0.12, 0.1,kFALSE,kFALSE,kFALSE,
                                vecInvMass,"","M_{#gamma#gamma}","#frac{1}{N} #frac{dN}{dM_{#gamma#gamma}}",1,1.4,
                                colorsM, kFALSE, 1, 1, kFALSE,
                                0.17,0.85,0.03,plotLabelP[p],plotLabelM,plotLabelZ[z],"mixed events",11);
            canvasInvMass->Update();
            canvasInvMass->SaveAs(Form("%s/InvMassBGProjectionNormalized_p%d_z%d_m.%s",outputDir.Data(),p,z,suffix.Data()));

            // compare mother histo between m bins
            canvasInvMass->Clear();
            vecInvMass.clear();
            for(Int_t m=0;m<nBinsM;m++) {
                vecInvMass.push_back(fHistoInvMassMother[p][z][m]);
            }
            DrawCompareHistoTH1(canvasInvMass,0.14, 0.02, 0.12, 0.1,kFALSE,kFALSE,kFALSE,
                                vecInvMass,"","M_{#gamma#gamma}","#frac{dN}{dM_{#gamma#gamma}}",1,1.4,
                                colorsM, kTRUE, 1.1, 1.1, kFALSE, 
                                0.17,0.85,0.03,plotLabelP[p],plotLabelM,plotLabelZ[z],"same events",11);
            canvasInvMass->Update();
            canvasInvMass->SaveAs(Form("%s/InvMassMotherProjection_p%d_z%d_m.%s",outputDir.Data(),p,z,suffix.Data()));

            // compare mother histo between m bins, normalized
            canvasInvMass->Clear();
            for(Int_t iVec=0; iVec<(Int_t)vecInvMass.size(); iVec++){
                TH1D* temp = vecInvMass.at(iVec);
                temp->Sumw2();
                temp->Scale(1./temp->Integral());
            }

            DrawCompareHistoTH1(canvasInvMass,0.14, 0.02, 0.12, 0.1,kFALSE,kFALSE,kFALSE,
                                vecInvMass,"","M_{#gamma#gamma}","#frac{1}{N} #frac{dN}{dM_{#gamma#gamma}}",1,1.4,
                                colorsM, kFALSE, 1, 1, kFALSE,
                                0.17,0.85,0.03,plotLabelP[p],plotLabelM,plotLabelZ[z],"same events",11);
            canvasInvMass->Update();
            canvasInvMass->SaveAs(Form("%s/InvMassMotherProjectionNormalized_p%d_z%d_m.%s",outputDir.Data(),p,z,suffix.Data()));
        }
    }

    //=========================================================================================  
    //========== plotting 2D histograms =======================================================
    //=========================================================================================

    canvas2D->cd();

    for(Int_t p=0; p<nBinsPtProj; p++) {

        canvas2D->Clear();

        for(Int_t z=0;z<nBinsZ;z++) {
             fHistoMotherIntegralMapZM[p]->GetXaxis()->SetBinLabel(z+1,plotLabelZ[z]);
             fHistoBckIntegralMapZM[p]->GetXaxis()->SetBinLabel(z+1,plotLabelZ[z]);
         }

        for(Int_t m=0;m<nBinsM;m++) {
             fHistoMotherIntegralMapZM[p]->GetYaxis()->SetBinLabel(m+1,plotLabelM[m]);
             fHistoBckIntegralMapZM[p]->GetYaxis()->SetBinLabel(m+1,plotLabelM[m]);
         }

        DrawAutoGammaHisto2D(fHistoMotherIntegralMapZM[p], "", "", "", "", kFALSE, 0, 0, kFALSE, 0, 0, 1, 1.2);
        canvas2D->Update();
        canvas2D->SaveAs(Form("%s/InvMassMotherProjectionIntegral2D_p%d.%s",outputDir.Data(),p,suffix.Data()));

        canvas2D->Clear();
        DrawAutoGammaHisto2D(fHistoBckIntegralMapZM[p], "", "", "", "", kFALSE, 0, 0, kFALSE, 0, 0, 1, 1.2);
        canvas2D->Update();
        canvas2D->SaveAs(Form("%s/InvMassBckProjectionIntegral2D_p%d.%s",outputDir.Data(),p,suffix.Data()));
   
    }

    //=========================================================================================  
    //========== saving and cleanup ===========================================================
    //=========================================================================================

    TString outputFileName = Form("%s/AnalyseBck_%s.root",outputDir.Data(),cutString.Data());
    TFile * outputFile = new TFile(outputFileName, "RECREATE");

    cout << "Write to file..." << endl;
    
    fHistoInvMassMotherFull->Write(); 
    fHistoInvMassBckFull->Write();
    
    for (Int_t p=0;p<0;p++) {
        for (Int_t z=0;z<nBinsZ;z++) {
            for (Int_t m=0;m<nBinsM;m++) {
                fHistoInvMassMother[p][z][m]->Write();
                fHistoInvMassBck[p][z][m]->Write();
            }
        }
    }

    outputFile->Close();
    delete outputFile;
    delete canvasInvMass;
    delete canvas2D;
    f->Close();
    delete f;

}
