#include <stdlib.h>
#include <iostream>
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
#include "TPad.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TDatabasePDG.h"
#include "TMinuit.h"
#include "TLatex.h"
#include "TASImage.h"
#include "TGraphErrors.h"
#include "TArrow.h"
#include "TGraphAsymmErrors.h" 
#include "TGaxis.h"
#include "TMath.h"
#include "TMarker.h"
#include "../CommonHeaders/PlottingMeson.h"
#include "../CommonHeaders/PlottingGammaConversionHistos.h"
#include "../CommonHeaders/PlottingGammaConversionAdditional.h"
#include "../CommonHeaders/FittingGammaConversion.h"
#include "../CommonHeaders/ConversionFunctions.h"
#include "../CommonHeaders/ExtractSignalBinning.h"

//****************************************************************************
//******************** Projection out of 2D in X *****************************
//****************************************************************************
TH1D* FillProjectionX (TH2F* fDummy2D, TString name, Double_t minY, Double_t maxY, Int_t rebin){
    TH1D* dummy1D           = new TH1D(name.Data(), name.Data(), fDummy2D->GetNbinsX(), 0., fDummy2D->GetXaxis()->GetBinUpEdge(fDummy2D->GetNbinsX()));
    dummy1D->Sumw2();
    Int_t startBin          = fDummy2D->GetYaxis()->FindBin(minY+0.001);
    Int_t endBin            = fDummy2D->GetYaxis()->FindBin(maxY-0.001);
    fDummy2D->ProjectionX(name.Data(),startBin,endBin,"e");
    dummy1D                 = (TH1D*)gDirectory->Get(name.Data());
    if(rebin>1){
        dummy1D->Rebin(rebin);
    }
    return dummy1D; 
}


//****************************************************************************
//******************** Correct Bin Contents for Chi2 Test ********************
//****************************************************************************
TH1D* CorrectBinContentForChi2Test ( TH1D* histo ){
    TH1D* histoDummy    = (TH1D*)histo->Clone(Form("%s_Chi2Test",histo->GetName())); 
    for (Int_t j = 1; j < histo->GetNbinsX()+1; j++){
        if (histoDummy->GetBinContent(j) == 0){
            histoDummy->SetBinError(j,1.14);   
        }    
    }    
    return histoDummy;
}    

//__________________________________________ Plotting _______________________________________________
void PlotMergedCompDataMCInPtBins(  TH1D**      fHistoMergedRecDataPtBinPlot, 
                                    TH1D**      fHistoMergedRecMCPtBinPlot, 
                                    TH1D**      fHistoMergedValMCPtBinPlot, 
                                    Double_t*   fChi2Value,
                                    Double_t*   fMeanData,
                                    Double_t*   fMeanMC,
                                    TString     namePlot, 
                                    TString     nameCanvas, 
                                    TString     namePad, 
                                    Double_t*   fPlottingRangeMeson, 
                                    TString     dateDummy, 
                                    TString     fMesonType, 
                                    Int_t       fRowPlot, 
                                    Int_t       fColumnPlot, 
                                    Int_t       fStartBinPtRange, 
                                    Int_t       fNumberPtBins, 
                                    Double_t*   fRangeBinsPt, 
                                    TString     fDecayChannel, 
                                    TString     decayChannel,  
                                    TString     fDetectionChannel, 
                                    TString     fEnergy,
                                    Int_t       plottingMode
                                        ){
    TGaxis::SetMaxDigits(3);
    TCanvas * canvasDataSpectra         = new TCanvas(nameCanvas.Data(),"",1400,900);  // gives the page size
    DrawGammaCanvasSettings(canvasDataSpectra, 0.0, 0., 0., 0.);

    TPad * padDataSpectra               = new TPad(namePad.Data(),"",0.0,0.0,1.,1.,0);   // gives the size of the histo areas
    padDataSpectra->SetFillColor(0);
    padDataSpectra->GetFrame()->SetFillColor(0);
    padDataSpectra->SetBorderMode(0);
    padDataSpectra->SetLogy(0);
    padDataSpectra->Divide(fColumnPlot,fRowPlot,0.0,0.0);
    padDataSpectra->Draw();

    cout<<"fColumnPlot: "<<fColumnPlot<<" fRowPlot: "<<fRowPlot<<endl;
    cout << fMesonType.Data() << endl;
    Int_t place                         = 0;
    for(Int_t iPt=fStartBinPtRange;iPt<fNumberPtBins;iPt++){
        cout<<"Pt: "<<iPt<<" of "<<fNumberPtBins<<endl;
        Double_t startPt                = fRangeBinsPt[iPt];
        Double_t endPt                  = fRangeBinsPt[iPt+1];

        place                           = place + 1;  //give the right place in the page
        if (place == fColumnPlot){
            iPt--;
            padDataSpectra->cd(place);

            cout << "entered ALICE plotting" << endl;
            TString textAlice           = "ALICE performance";
            Double_t nPixels            = 13;
            Double_t textHeight         = 0.08;
            if (padDataSpectra->cd(place)->XtoPixel(padDataSpectra->cd(place)->GetX2()) < padDataSpectra->cd(place)->YtoPixel(padDataSpectra->cd(place)->GetY1())){
                textHeight              = (Double_t)nPixels/padDataSpectra->cd(place)->XtoPixel(padDataSpectra->cd(place)->GetX2()) ;
            } else {
                textHeight              = (Double_t)nPixels/padDataSpectra->cd(place)->YtoPixel(padDataSpectra->cd(place)->GetY1());
            }   
            
            Double_t startTextX         = 0.10;
            Double_t startTextY         = 0.9;
            Double_t differenceText     = textHeight*1.25;

            TLatex *alice               = new TLatex(startTextX, startTextY, Form("%s",textAlice.Data()));
            TLatex *latexDate           = new TLatex(startTextX, (startTextY-1.25*differenceText), dateDummy.Data());
            TLatex *energy              = new TLatex(startTextX, (startTextY-2.25*differenceText), fEnergy);
            TLatex *process             = new TLatex(startTextX, (startTextY-3.25*differenceText), fDecayChannel);
            TLatex *detprocess          = new TLatex(startTextX, (startTextY-4.25*differenceText), fDetectionChannel);

            alice->SetNDC();
            alice->SetTextColor(1);
            alice->SetTextSize(textHeight*1.3);
            alice->Draw();

            latexDate->SetNDC();
            latexDate->SetTextColor(1);
            latexDate->SetTextSize(textHeight);
            latexDate->Draw();

            energy->SetNDC();
            energy->SetTextColor(1);
            energy->SetTextSize(textHeight);
            energy->Draw();

            process->SetNDC(); 
            process->SetTextColor(1);
            process->SetTextSize(textHeight);
            process->Draw();

            detprocess->SetNDC(); 
            detprocess->SetTextColor(1);
            detprocess->SetTextSize(textHeight);
            detprocess->Draw();

            Int_t nLegend               = 0;
            if (fHistoMergedRecDataPtBinPlot) nLegend++;
            if (fHistoMergedRecMCPtBinPlot) nLegend++;
            if (fHistoMergedValMCPtBinPlot) nLegend++;
            
            TLegend* legend             = new TLegend(startTextX,startTextY-5.75*differenceText,1,startTextY-(5.75+nLegend)*differenceText);
            legend->SetTextSize(textHeight);
            legend->SetTextFont(62);
            legend->SetFillColor(0);
            legend->SetFillStyle(0);
            legend->SetLineWidth(0);
            legend->SetLineColor(0);
            legend->SetMargin(0.15);
            Size_t markersize           = 0;
            Size_t linesize             = 0;
            if (fHistoMergedRecDataPtBinPlot){
                if (fHistoMergedRecDataPtBinPlot[iPt]){
                    markersize          = fHistoMergedRecDataPtBinPlot[iPt]->GetMarkerSize();
                    fHistoMergedRecDataPtBinPlot[iPt]->SetMarkerSize(1.5*markersize);
                    legend->AddEntry(fHistoMergedRecDataPtBinPlot[iPt],"merged clus. Data","ep");
                }
            }
            if (fHistoMergedRecMCPtBinPlot){
                if (fHistoMergedRecMCPtBinPlot[iPt]){
                    markersize          = fHistoMergedRecMCPtBinPlot[iPt]->GetMarkerSize();
                    fHistoMergedRecMCPtBinPlot[iPt]->SetMarkerSize(3*markersize);
                    legend->AddEntry(fHistoMergedRecMCPtBinPlot[iPt],"merged clus. MC","l");
                }
            }
            if (fHistoMergedValMCPtBinPlot){
                if (fHistoMergedValMCPtBinPlot[iPt]){
                    markersize          = fHistoMergedValMCPtBinPlot[iPt]->GetMarkerSize();
                    fHistoMergedValMCPtBinPlot[iPt]->SetMarkerSize(3*markersize);
                    legend->AddEntry(fHistoMergedValMCPtBinPlot[iPt],"merged clus. val. MC","l");
                }
            }
            legend->Draw(); 
            
        } else {            
            padDataSpectra->cd(place);
            padDataSpectra->cd(place)->SetTopMargin(0.12);
            padDataSpectra->cd(place)->SetBottomMargin(0.15);
            padDataSpectra->cd(place)->SetRightMargin(0.02);

            Double_t nPixels            = 13;
            Double_t textHeight         = 0.08;
            if (padDataSpectra->cd(place)->XtoPixel(padDataSpectra->cd(place)->GetX2()) < padDataSpectra->cd(place)->YtoPixel(padDataSpectra->cd(place)->GetY1())){
                textHeight              = (Double_t)nPixels/padDataSpectra->cd(place)->XtoPixel(padDataSpectra->cd(place)->GetX2()) ;
            } else {
                textHeight              = (Double_t)nPixels/padDataSpectra->cd(place)->YtoPixel(padDataSpectra->cd(place)->GetY1());
            }   

            //             cout << "place: " << place << "\t "<< startPt << " < #it{p}_{T} < " << endPt  << endl;
            int remaining               = (place-1)%fColumnPlot;
            if (remaining > 0) padDataSpectra->cd(place)->SetLeftMargin(0.15);
            else padDataSpectra->cd(place)->SetLeftMargin(0.25);
//             cout << "here" << endl;
            TString labelAbove  = Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}",startPt,endPt);
            TString xaxisLabel  = "";
            TString yaxisLabel  = "";
            if (plottingMode == 0){
                xaxisLabel  = Form("#it{M}_{%s} (GeV/#it{c}^{2})",decayChannel.Data());
                yaxisLabel  = Form("dN_{%s}/d#it{M}_{%s} norm int",decayChannel.Data(), decayChannel.Data());
            } else if  (plottingMode == 1){
                xaxisLabel  = "#lambda_{0}^{2}";
                yaxisLabel  = Form("dN_{%s}/d#lambda_{0}^{2} norm int", decayChannel.Data());
            } else if (plottingMode == 2){
                labelAbove  = Form("%3.2f GeV/#it{c} < #it{E} < %3.2f GeV/#it{c}",startPt,endPt);
            }
            
            if (fHistoMergedRecDataPtBinPlot){
                if (fHistoMergedRecDataPtBinPlot[iPt]){
                    if (plottingMode != 2 ){
                        DrawGammaHistoColored( fHistoMergedRecDataPtBinPlot[iPt], labelAbove, xaxisLabel, yaxisLabel, fPlottingRangeMeson[0],fPlottingRangeMeson[1],1,kBlack,20,0.8);
                    } else {
                        DrawGammaHistoColored( fHistoMergedRecDataPtBinPlot[iPt], labelAbove, xaxisLabel, yaxisLabel, fPlottingRangeMeson[0],fPlottingRangeMeson[1],0,kBlack,20,0.8);
                        fHistoMergedRecDataPtBinPlot[iPt]->Draw("pe");
                    }    
                } else {
                    continue;
                }
            } else {
                cout << "ABORT: plotting no default histos" << endl;
                return;
            }
            Double_t startLabels        = 0.7;
            if (remaining > 0)          
                startLabels             = 0.66;
            if (fColumnPlot > 4){
                startLabels             = 0.59;
                if (remaining > 0)          
                    startLabels         = 0.55;
            }
            if (fChi2Value){
               if (fChi2Value[iPt] != 0 && fChi2Value[iPt] != 10000){
                    TLatex *chi2Value           = new TLatex(startLabels, 0.8, Form("#chi^{2}/ndf = %2.2f",fChi2Value[iPt]));
                    chi2Value->SetNDC();
                    chi2Value->SetTextColor(1);
                    chi2Value->SetTextSize(textHeight*1.3);
                    chi2Value->Draw();   
                }    
            }    
            if (fMeanData){
               if (fMeanData[iPt] != 10000){
                    TLatex *meanData           = new TLatex(startLabels, 0.8-textHeight*1.3, Form("#mu_{Data} = %2.2f",fMeanData[iPt]));
                    meanData->SetNDC();
                    meanData->SetTextColor(1);
                    meanData->SetTextSize(textHeight*1.3);
                    meanData->Draw();   
                }    
            }    
            if (fMeanMC){
               if (fMeanMC[iPt] != 10000){
                    TLatex *meanMC           = new TLatex(startLabels, 0.8-2*textHeight*1.3, Form("#mu_{MC} = %2.2f",fMeanMC[iPt]));
                    meanMC->SetNDC();
                    meanMC->SetTextColor(1);
                    meanMC->SetTextSize(textHeight*1.3);
                    meanMC->Draw();   
                }    
            }    
            
            if (fHistoMergedRecMCPtBinPlot){
                if (fHistoMergedRecMCPtBinPlot[iPt]){
                    DrawGammaHistoColored(  fHistoMergedRecMCPtBinPlot[iPt], labelAbove, xaxisLabel, yaxisLabel, fPlottingRangeMeson[0],fPlottingRangeMeson[1],0,kBlue+1);
                }
            }
            if (fHistoMergedValMCPtBinPlot){
                if (fHistoMergedValMCPtBinPlot[iPt]){
                    DrawGammaHistoColored(  fHistoMergedValMCPtBinPlot[iPt], labelAbove, xaxisLabel, yaxisLabel, fPlottingRangeMeson[0],fPlottingRangeMeson[1],0,kGreen+2);
                }
            }
            cout << iPt << endl;
        }
//         cout << "here " << endl;
    }
    canvasDataSpectra->Print(namePlot.Data());
//     cout << "here " << endl;
    delete padDataSpectra;
    delete canvasDataSpectra;
//     cout << "here " << endl;
}


//__________________________________________ Plotting _______________________________________________
void PlotMergedCompDataMCSingleBin( TH1D**      fHistoMergedRecDataPtBinPlot, 
                                    TH1D**      fHistoMergedRecMCPtBinPlot, 
                                    TH1D**      fHistoMergedValMCPtBinPlot, 
                                    Double_t*   fChi2Value,
                                    Double_t*   fMeanData,
                                    Double_t*   fMeanMC,
                                    TString     namePlot, 
                                    TString     nameCanvas, 
                                    TString     namePad, 
                                    Double_t*   fPlottingRangeMeson, 
                                    TString     dateDummy, 
                                    TString     fMesonType, 
                                    Int_t       iPtBin, 
                                    Double_t*   fRangeBinsPt, 
                                    Int_t       detMode, 
                                    TString     fEnergy,
                                    Int_t       trigger,
                                    Int_t       plotttingMode   = 0
                                    ){
    
    TString methodStr               = ReturnTextReconstructionProcess(detMode);
    TString triggerStr2             = ReturnTriggerName(trigger,fEnergy);
    TString triggerStr              = Form("%s triggered", triggerStr2.Data());

    if (!fHistoMergedRecDataPtBinPlot[iPtBin]){
        cout << "I am missing the data hist" << endl;
        return;  
    } 
    
    Double_t marginUp       = 0.015;
    if (plotttingMode == 2)
        marginUp       = 0.035;
    TCanvas* canvasMergedSamplePlot     = new TCanvas("canvasMergedSamplePlotNew","",0,0,1500,1500);  // gives the page size
    DrawGammaCanvasSettings( canvasMergedSamplePlot,  0.09, 0.015, marginUp, 0.08);
    
    Double_t startPt                    = fRangeBinsPt[iPtBin];
    Double_t endPt                      = fRangeBinsPt[iPtBin+1];
    Double_t startTextY                 = 0.945-marginUp;
    
    Double_t textSizeLabelsPixel        = 100*3/5;
    Double_t marginSample               = 0.1*1500;
    Double_t textsizeLabelsSample       = 0;
    Double_t textsizeFacSample          = 0;
    if (canvasMergedSamplePlot->XtoPixel(canvasMergedSamplePlot->GetX2()) < canvasMergedSamplePlot->YtoPixel(canvasMergedSamplePlot->GetY1())){
        textsizeLabelsSample       = (Double_t)textSizeLabelsPixel/canvasMergedSamplePlot->XtoPixel(canvasMergedSamplePlot->GetX2()) ;
        textsizeFacSample          = (Double_t)1./canvasMergedSamplePlot->XtoPixel(canvasMergedSamplePlot->GetX2()) ;
    } else {
        textsizeLabelsSample       = (Double_t)textSizeLabelsPixel/canvasMergedSamplePlot->YtoPixel(canvasMergedSamplePlot->GetY1());
        textsizeFacSample          = (Double_t)1./canvasMergedSamplePlot->YtoPixel(canvasMergedSamplePlot->GetY1());
    }

    TString xaxisLabel  = fHistoMergedRecDataPtBinPlot[iPtBin]->GetXaxis()->GetTitle();
    TString yaxisLabel  = fHistoMergedRecDataPtBinPlot[iPtBin]->GetYaxis()->GetTitle();
    
    TH1F * histo1DSampleDummy;
    histo1DSampleDummy             = new TH1F("histo1DSample2","histo1DSample2",11000,fPlottingRangeMeson[0],fPlottingRangeMeson[1]);
    SetStyleHistoTH1ForGraphs(histo1DSampleDummy, xaxisLabel, yaxisLabel, 0.85*textsizeLabelsSample, textsizeLabelsSample,
                            0.85*textsizeLabelsSample, textsizeLabelsSample,0.88, 0.115/(textsizeFacSample*marginSample));
    histo1DSampleDummy->GetYaxis()->SetLabelOffset(0.008);
    histo1DSampleDummy->GetXaxis()->SetLabelOffset(0.005);

    TString ptLabel         =  "#it{p}_{T} ";
    TString unit            =  "GeV/#it{c}";
    if (plotttingMode == 2){
        unit                =  "GeV";
        ptLabel             =  "#it{E} ";
    }    
    // Set range for fits and labels
    TLatex *labelSamplePtRange;

    if(fMesonType.CompareTo("Pi0") == 0 || fMesonType.CompareTo("Pi0EtaBinning") == 0){
        labelSamplePtRange = new TLatex(0.95,startTextY, Form("#pi^{0}: %3.1f %s < %s< %3.1f %s",startPt, unit.Data(), ptLabel.Data(),endPt, unit.Data() ));
    } else {
        labelSamplePtRange = new TLatex(0.95,startTextY, Form("#eta:  %3.1f %s < %s< %3.1f %s",startPt, unit.Data(), ptLabel.Data(),endPt, unit.Data() ));
    }

    Double_t minimum = fHistoMergedRecDataPtBinPlot[iPtBin]->GetMinimum();
    if (minimum < 0) minimum = 1.1*minimum;
    else if (minimum == 0) minimum = -0.001;
    else minimum = 0.9*minimum;

    Double_t scale = 1;
    if (plotttingMode == 1)
        scale = 1.3;
    canvasMergedSamplePlot->cd();
    histo1DSampleDummy->GetYaxis()->SetRangeUser(minimum,scale*fHistoMergedRecDataPtBinPlot[iPtBin]->GetMaximum());
    histo1DSampleDummy->Draw("AXIS");
  
    DrawGammaSetMarker(fHistoMergedRecDataPtBinPlot[iPtBin], 20, 2, kBlack, kBlack);
    fHistoMergedRecDataPtBinPlot[iPtBin]->Draw("same");
    if (fHistoMergedRecMCPtBinPlot){
        if (fHistoMergedRecMCPtBinPlot[iPtBin]){
            DrawGammaSetMarker(fHistoMergedRecMCPtBinPlot[iPtBin], 24, 2, kRed+2, kRed+2);
            fHistoMergedRecMCPtBinPlot[iPtBin]->Draw("same");
        }
    }

    // labels left side
    TLatex *labelALICE      = new TLatex(0.135,startTextY,"ALICE this work");
    SetStyleTLatex( labelALICE, 0.85*textSizeLabelsPixel,4, 1, 43);
    labelALICE->Draw();    
    TLatex *labelSampleEnergy      = new TLatex(0.135,startTextY-0.9*textsizeLabelsSample,fEnergy.Data());
    SetStyleTLatex( labelSampleEnergy, 0.85*textSizeLabelsPixel,4, 1, 43);
    labelSampleEnergy->Draw();
    TLatex *labelTrigger  = new TLatex(0.135,startTextY-2*0.9*textsizeLabelsSample,triggerStr.Data());
    SetStyleTLatex( labelTrigger, 0.85*textSizeLabelsPixel,4, 1, 43);
    labelTrigger->Draw();
    TLatex *labelSampleReco  = new TLatex(0.135,startTextY-3*0.9*textsizeLabelsSample, methodStr);
    SetStyleTLatex( labelSampleReco, 0.85*textSizeLabelsPixel,4, 1, 43);
    labelSampleReco->Draw();

    // labels right side
    SetStyleTLatex( labelSamplePtRange, 0.85*textSizeLabelsPixel,4, 1, 43, kTRUE, 31);
    labelSamplePtRange->Draw();
    TLatex *chi2Value  = new TLatex(0.95, startTextY-4*0.9*textsizeLabelsSample, Form("#chi^{2}/ndf = %2.2f",fChi2Value[iPtBin]));
    SetStyleTLatex( chi2Value, 0.85*textSizeLabelsPixel,4, 1, 43, kTRUE, 31);
    chi2Value->Draw();
    TLatex *meanData  = new TLatex(0.95, startTextY-5*0.9*textsizeLabelsSample, Form("#mu_{Data} = %2.2f",fMeanData[iPtBin]));
    SetStyleTLatex( meanData, 0.85*textSizeLabelsPixel,4, 1, 43, kTRUE, 31);
    meanData->Draw();
    TLatex *meanMC  = new TLatex(0.95, startTextY-6*0.9*textsizeLabelsSample, Form("#mu_{MC} = %2.2f",fMeanMC[iPtBin]));
    SetStyleTLatex( meanMC, 0.85*textSizeLabelsPixel,4, 1, 43, kTRUE, 31);
    meanMC->Draw();
    
    Int_t nLegendLines               = 0;
    if (fHistoMergedRecDataPtBinPlot) nLegendLines++;
    if (fHistoMergedRecMCPtBinPlot) nLegendLines++;

    TLegend* legendSample  = GetAndSetLegend2(0.82, startTextY-0.01-nLegendLines*0.9*textsizeLabelsSample, 0.95, startTextY-0.01, 0.85*textSizeLabelsPixel);
    legendSample->SetMargin(0.4);
    if( fHistoMergedRecDataPtBinPlot[iPtBin]) legendSample->AddEntry(fHistoMergedRecDataPtBinPlot[iPtBin],"Data","p");
    if( fHistoMergedRecMCPtBinPlot[iPtBin]) legendSample->AddEntry(fHistoMergedRecMCPtBinPlot[iPtBin],"MC","p");
    legendSample->Draw();
    histo1DSampleDummy->Draw("AXIS,same");

    cout << "done plotting single" << endl;
    canvasMergedSamplePlot->Print(namePlot.Data());
    delete canvasMergedSamplePlot;
}



//************************************************************************************************************
//******************************* Main routine for comparing MC & data ***************************************
//************************************************************************************************************
void CompareShapeMergedClusterQuantities(   TString dataFileName    = "rawSignalData", 
                                            TString mcFileName      = "rawSignalMC", 
                                            TString fCutSelection   = "", 
                                            TString mesonType       = "Pi0", 
                                            TString fSuffix         = "", 
                                            TString energyFlag      = "" ,
                                            Int_t numberOfBins      = 25,
                                            Int_t mode              = 10
                                        )
{
    gROOT->Reset();
    // mode:    10 // merged cluster EMC
    //          11 // merged cluster PHOS
    
    if (!(mode == 10 || mode == 11)){
        cout << "ERROR: You are running this macro in the wrong mode, aborting ..." << endl;
        return;
    }    

    StyleSettingsThesis(fSuffix);   
    SetPlotStyle();
    TFile fileRawSignalData(dataFileName.Data());
    TFile fileRawSignalMC(mcFileName.Data());

    TString fdate                       = ReturnDateString();
    cout << dataFileName.Data() << endl;
    cout << mcFileName.Data() << endl;
    cout << fCutSelection.Data() << endl;
    cout << mesonType.Data() << endl;
    cout << fSuffix.Data()<< endl;
    cout << energyFlag.Data() << endl;
    cout << numberOfBins << endl;
    
    Bool_t fAdvanceClusterQA            = kFALSE;
    
    Double_t *fMesonMassRange           = NULL;
    Double_t *fMesonM02Range            = NULL;
    TString outputDir                   = Form("%s/%s/%s/ExtractSignalMergedMeson",fCutSelection.Data(),energyFlag.Data(),fSuffix.Data());

    TString fEventCutSelection          = "";
    TString fClusterCutSelection        = "";
    TString fClusterMergedCutSelection  = "";
    TString fMesonCutSelection          = "";
    TString dummyString                 = "";
    ReturnSeparatedCutNumberAdvanced( fCutSelection, fEventCutSelection, fClusterCutSelection, fClusterMergedCutSelection, dummyString, fMesonCutSelection, mode);
    Int_t fTrigger                      = ((TString)fEventCutSelection(3,2)).Atoi();
    
    //****************************** Specification of collision system ************************************************
    TString textProcess                 = ReturnMesonString (mesonType);
    if(textProcess.CompareTo("") == 0 ){
        cout << "Meson unknown" << endl;
        return ;
    }
    
    TString fTextMeasurement            = Form("%s #rightarrow #gamma#gamma", textProcess.Data());
    TString fCollisionSystem            = ReturnFullCollisionsSystem(energyFlag);
    TString fDecayChannel               = "#gamma#gamma";
    if (fCollisionSystem.CompareTo("") == 0){
        cout << "No correct collision system specification, has been given" << endl;
        return;
    }
    TString fDetectionProcess           = ReturnFullTextReconstructionProcess(mode,0,textProcess.Data());
    TString fDetectionProcessPtBins     = ReturnFullTextReconstructionProcess(mode,0,textProcess.Data(), fClusterMergedCutSelection);

    TString fNLMString                  = "";
    Int_t fNLMmin                       = ReturnClusterNLM(fClusterMergedCutSelection);
    if (ReturnClusterNLM(fClusterMergedCutSelection) == 1)
        fNLMString                      = Form("%i local maximum", ReturnClusterNLM(fClusterMergedCutSelection));        
    else 
        fNLMString                      = Form("%i local maxima", ReturnClusterNLM(fClusterMergedCutSelection));
    
    InitializeBinning(mesonType, numberOfBins, energyFlag, "", mode, fEventCutSelection, fClusterMergedCutSelection);
    
    if (mesonType.CompareTo("Pi0") == 0 || mesonType.CompareTo("Pi0EtaBinning") == 0){
        fMesonMassRange         = new Double_t[2]; 
        fMesonMassRange[0]      = 0.; 
        fMesonMassRange[1]      = 0.3;
        fMesonM02Range          = new Double_t[2]; 
        fMesonM02Range[0]       = 0.; 
        fMesonM02Range[1]       = 4.8;
        if (fNLMmin == 1)
            fMesonM02Range[1]   = 2.;
    } else if (mesonType.CompareTo("Eta") == 0){
        fMesonMassRange         = new Double_t[2]; 
        fMesonMassRange[0]      = 0.35; 
        fMesonMassRange[1]      = 0.79;
        fMesonM02Range          = new Double_t[2]; 
        fMesonM02Range[0]       = 0.; 
        fMesonM02Range[1]       = 4.8;
    } 

    cout << fStartPtBin << endl;
    
    //******************************* Reading histograms **************************************************************

    TString histonameSignal;
    TString histonameMCTruth;
    Double_t scaleFacData;
    Double_t scaleFacMC;

    Double_t chi2NDFInvMass                 [70];
    Double_t chi2NDFM02                     [70];
    for (Int_t i = 0; i < 50; i++){
        chi2NDFInvMass[i]                   = 10000;
        chi2NDFM02[i]                       = 10000;
    }
    TH1D* histoDataInvMassPtBin             [70];
    TH1D* histoMCInvMassPtBin               [70];
    TH1D* histoValidatedInvMassMCPtBin      [70];
    TH1D* histoDataM02PtBin                 [70];
    TH1D* histoMCM02PtBin                   [70];
    TH1D* histoValidatedM02MCPtBin          [70];
        
    for(Int_t iPt = fStartPtBin; iPt < fNBinsPt; iPt++){
        scaleFacData                        = 1;
        scaleFacMC                          = 1;
        cout << "reading bin: " << iPt << endl;
        histonameSignal                     = Form("InvMass_PtBin%02d", iPt);
        histonameMCTruth                    = Form("ValidatedMerged_InvMass_PtBin%02d", iPt);
        TH1D* histoDataInvMassChi2Test      = NULL;
        TH1D* histoMCInvMassChi2Test        = NULL;
        
        histoDataInvMassPtBin[iPt]          = (TH1D*)fileRawSignalData.Get(histonameSignal);
        if (histoDataInvMassPtBin[iPt]){
           histoDataInvMassChi2Test         = CorrectBinContentForChi2Test(histoDataInvMassPtBin[iPt]);
        }    
        scaleFacData                        = histoDataInvMassPtBin[iPt]->Integral();
        if (scaleFacData > 0 && histoDataInvMassPtBin[iPt]){
            histoDataInvMassPtBin[iPt]->Scale(1./scaleFacData);
            if (histoDataInvMassChi2Test) histoDataInvMassChi2Test->Scale(1./scaleFacData);
        }
        
        
        histoMCInvMassPtBin[iPt]            = (TH1D*)fileRawSignalMC.Get(histonameSignal);
        if (histoMCInvMassPtBin[iPt]){
           histoMCInvMassChi2Test           = CorrectBinContentForChi2Test(histoMCInvMassPtBin[iPt]);
        }    

        if (histoMCInvMassPtBin[iPt]){
            scaleFacMC                      = histoMCInvMassPtBin[iPt]->Integral();
        }    
        if (scaleFacMC > 0 && histoMCInvMassPtBin[iPt]){
            histoMCInvMassPtBin[iPt]->Scale(1./scaleFacMC);    
            if (histoMCInvMassChi2Test) histoMCInvMassChi2Test->Scale(1./scaleFacMC);
        }    
        
        histoValidatedInvMassMCPtBin[iPt]   = (TH1D*)fileRawSignalMC.Get(histonameMCTruth);
        if (scaleFacMC > 0 && histoValidatedInvMassMCPtBin[iPt]){
            histoValidatedInvMassMCPtBin[iPt]->Scale(1./scaleFacMC);
        }
        if (histoDataInvMassPtBin[iPt] && histoMCInvMassPtBin[iPt]){
            chi2NDFInvMass[iPt]             = histoDataInvMassChi2Test->Chi2Test(histoMCInvMassChi2Test,"NORM,CHI2/NDF");
            cout << "Inv Mass: "<< chi2NDFInvMass[iPt] << endl;
        }
        
        histonameSignal                     = Form("M02_PtBin%02d", iPt);
        histonameMCTruth                    = Form("ValidatedMerged_M02_PtBin%02d", iPt);
        TH1D* histoDataM02Chi2Test          = NULL;
        TH1D* histoMCM02Chi2Test            = NULL;
        
        histoDataM02PtBin[iPt]              = (TH1D*)fileRawSignalData.Get(histonameSignal);
        if (histoDataM02PtBin[iPt]){
           histoDataM02Chi2Test             = CorrectBinContentForChi2Test(histoDataM02PtBin[iPt]);
        }    
        scaleFacData                        = histoDataM02PtBin[iPt]->Integral();
        if (scaleFacData > 0 && histoDataM02PtBin[iPt]){
            histoDataM02PtBin[iPt]->Scale(1./scaleFacData);
            if (histoDataM02Chi2Test) histoDataM02Chi2Test->Scale(1./scaleFacData);
        }
        
        
        histoMCM02PtBin[iPt]                = (TH1D*)fileRawSignalMC.Get(histonameSignal);
        if (histoMCM02PtBin[iPt]){
           histoMCM02Chi2Test               = CorrectBinContentForChi2Test(histoMCM02PtBin[iPt]);
        }    

        if (histoMCM02PtBin[iPt]){
            scaleFacMC                      = histoMCM02PtBin[iPt]->Integral();
        }    
        if (scaleFacMC > 0 && histoMCM02PtBin[iPt]){
            histoMCM02PtBin[iPt]->Scale(1./scaleFacMC);    
            if (histoMCM02Chi2Test) histoMCM02Chi2Test->Scale(1./scaleFacMC);
        }    
        
        histoValidatedM02MCPtBin[iPt]       = (TH1D*)fileRawSignalMC.Get(histonameMCTruth);
        if (scaleFacMC > 0 && histoValidatedM02MCPtBin[iPt]){
            histoValidatedM02MCPtBin[iPt]->Scale(1./scaleFacMC);
        }
        if (histoDataM02PtBin[iPt] && histoMCM02PtBin[iPt]){
            chi2NDFM02[iPt]                 = histoDataM02Chi2Test->Chi2Test(histoMCM02Chi2Test,"NORM,CHI2/NDF");
            cout << "M02: "<< chi2NDFM02[iPt] << endl;
        }
        
    }   
    TString nameCompPlotInvMass             = Form("%s/%s_MesonInvMassCompared_%s.%s",outputDir.Data(),mesonType.Data(),fCutSelection.Data(),fSuffix.Data());
    PlotMergedCompDataMCInPtBins(   histoDataInvMassPtBin, histoMCInvMassPtBin, histoValidatedInvMassMCPtBin,  chi2NDFInvMass, NULL, NULL,
                                    nameCompPlotInvMass, "DOESNTMATTER_1", "IDONTCARE_1", fMesonMassRange, 
                                    fdate, mesonType, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fDecayChannel, fDetectionProcessPtBins, fCollisionSystem, 0);

    TString nameCompPlotM02                 = Form("%s/%s_MesonM02Compared_%s.%s",outputDir.Data(),mesonType.Data(),fCutSelection.Data(),fSuffix.Data()) ; 
    PlotMergedCompDataMCInPtBins(   histoDataM02PtBin, histoMCM02PtBin, histoValidatedM02MCPtBin, chi2NDFM02, NULL, NULL,
                                    nameCompPlotM02,  "bla", "bla", fMesonM02Range, 
                                    fdate, mesonType, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fDecayChannel, fDetectionProcessPtBins, fCollisionSystem, 1);
    
    TH2F* histoDataMergedNCellsPt           = NULL;
    TH2F* histoMCMergedNCellsPt             = NULL;
    histoDataMergedNCellsPt                 = (TH2F*)fileRawSignalData.Get("ClusterMergedNCellsPt");
    histoMCMergedNCellsPt                   = (TH2F*)fileRawSignalMC.Get("ClusterMergedNCellsPt");
    TH2F* histoDataMergedNCellsAroundPt     = NULL;
    TH2F* histoMCMergedNCellsAroundPt       = NULL;
    histoDataMergedNCellsAroundPt           = (TH2F*)fileRawSignalData.Get("ClusterMergedNCellsAroundClusterPt");
    histoMCMergedNCellsAroundPt             = (TH2F*)fileRawSignalMC.Get("ClusterMergedNCellsAroundClusterPt");
    TH2F* histoDataMergedNCellsAroundAndInPt= NULL;
    TH2F* histoMCMergedNCellsAroundAndInPt  = NULL;
    histoDataMergedNCellsAroundAndInPt      = (TH2F*)fileRawSignalData.Get("ClusterMergedNCellsAroundAndInClusterPt");
    histoMCMergedNCellsAroundAndInPt        = (TH2F*)fileRawSignalMC.Get("ClusterMergedNCellsAroundAndInClusterPt");
    TH2F* histoDataMergedEAroundE           = NULL;
    TH2F* histoMCMergedEAroundE             = NULL;
    histoDataMergedEAroundE                 = (TH2F*)fileRawSignalData.Get("ClusterMergedEAroundClusterE");
    histoMCMergedEAroundE                   = (TH2F*)fileRawSignalMC.Get("ClusterMergedEAroundClusterE");

    if (histoDataMergedNCellsPt && histoMCMergedNCellsPt) 
        fAdvanceClusterQA                   = kTRUE;
    
    if (fAdvanceClusterQA){
        TH1D* histoDataMergedNCellsPtBin            [70];
        TH1D* histoMCMergedNCellsPtBin              [70];
        TH1D* histoDataMergedNCellsAroundPtBin      [70];
        TH1D* histoMCMergedNCellsAroundPtBin        [70];
        TH1D* histoDataMergedNCellsAroundAndInPtBin [70];
        TH1D* histoMCMergedNCellsAroundAndInPtBin   [70];
        TH1D* histoDataMergedEAroundEBin            [70];
        TH1D* histoMCMergedEAroundEBin              [70];
        Double_t chi2NDFNCells                      [70];
        Double_t chi2NDFNCellsAround                [70];
        Double_t chi2NDFNCellsAroundAndIn           [70];
        Double_t chi2NDFEAround                     [70];
        Double_t meanDataNCells                     [70];
        Double_t meanDataNCellsAround               [70];
        Double_t meanDataNCellsAroundAndIn          [70];
        Double_t meanDataEAround                    [70];
        Double_t meanMCNCells                       [70];
        Double_t meanMCNCellsAround                 [70];
        Double_t meanMCNCellsAroundAndIn            [70];
        Double_t meanMCEAround                      [70];
        for (Int_t i = 0; i < 50; i++){
            chi2NDFNCells[i]                        = 10000;
            chi2NDFNCellsAround[i]                  = 10000;
            chi2NDFNCellsAroundAndIn[i]             = 10000;
            chi2NDFEAround[i]                       = 10000;
            meanDataNCells[i]                       = 10000;
            meanDataNCellsAround[i]                 = 10000;
            meanDataNCellsAroundAndIn[i]            = 10000;
            meanDataEAround[i]                      = 10000;
            meanMCNCells[i]                         = 10000;
            meanMCNCellsAround[i]                   = 10000;
            meanMCNCellsAroundAndIn[i]              = 10000;
            meanMCEAround[i]                        = 10000;
        }

        
        Double_t fMergedNCellsRange[2]              = {-0.5,55.5};
        Double_t fMergedNCellsAroundRange[2]        = {1,55.5};
        Double_t fMergedERange[2]                   = {10.5,50};
        Double_t fMergedNCellsAroundAndInRange[2]   = {1,75.5};
        
        for (Int_t iPt = fStartPtBin; iPt < fNBinsPt; iPt++){
            // histos for NCells in merged cluster
            TH1D* histoDataNCellsChi2Test       = NULL;
            TH1D* histoMCNCellsChi2Test         = NULL;

            TString fNameDataHistoMergedNCells  = Form("NCellsMerged_in_Pt_Bin%02d", iPt);
            if (histoDataMergedNCellsPt){
                histoDataMergedNCellsPtBin[iPt] = FillProjectionX(histoDataMergedNCellsPt, fNameDataHistoMergedNCells, fBinsPt[iPt], fBinsPt[iPt+1], 1);
                if (histoDataMergedNCellsPtBin[iPt]){
                    meanDataNCells[iPt]         = histoDataMergedNCellsPtBin[iPt]->GetMean();
                    histoDataNCellsChi2Test     = CorrectBinContentForChi2Test(histoDataMergedNCellsPtBin[iPt]);
                }    
                scaleFacData                    = histoDataMergedNCellsPtBin[iPt]->Integral();
                if (scaleFacData > 0){
                    histoDataMergedNCellsPtBin[iPt]->Scale(1./scaleFacData);
                    if (histoDataNCellsChi2Test) histoDataNCellsChi2Test->Scale(1./scaleFacData);
                }    
            }
            
            TString fNameMCHistoMergedNCells    = Form("NCellsMerged_in_Pt_Bin%02d", iPt);
            if (histoMCMergedNCellsPt){
                histoMCMergedNCellsPtBin[iPt]   = FillProjectionX(histoMCMergedNCellsPt, fNameMCHistoMergedNCells, fBinsPt[iPt], fBinsPt[iPt+1], 1);
                if (histoMCMergedNCellsPtBin[iPt]){
                    meanMCNCells[iPt]           = histoMCMergedNCellsPtBin[iPt]->GetMean();
                    histoMCNCellsChi2Test       = CorrectBinContentForChi2Test(histoMCMergedNCellsPtBin[iPt]);
                }    
                scaleFacMC                      = histoMCMergedNCellsPtBin[iPt]->Integral();
                if (scaleFacMC > 0){
                    histoMCMergedNCellsPtBin[iPt]->Scale(1./scaleFacMC);    
                    if (histoMCNCellsChi2Test) histoMCNCellsChi2Test->Scale(1./scaleFacMC);
                }    
            }
            if (histoDataNCellsChi2Test && histoMCNCellsChi2Test){
                 chi2NDFNCells[iPt]             = histoDataNCellsChi2Test->Chi2Test(histoMCNCellsChi2Test,"NORM,CHI2/NDF");   
            }    

            // histos for NCells around merged cluster delta R < 0.15
            TH1D* histoDataNCellsAroundChi2Test         = NULL;
            TH1D* histoMCNCellsAroundChi2Test           = NULL;

            TString fNameDataHistoMergedNCellsAround    = Form("NCellsAroundMerged_in_Pt_Bin%02d", iPt);
            if (histoDataMergedNCellsAroundPt){
                histoDataMergedNCellsAroundPtBin[iPt]   = FillProjectionX(histoDataMergedNCellsAroundPt, fNameDataHistoMergedNCellsAround, fBinsPt[iPt], fBinsPt[iPt+1], 1);
                if (histoDataMergedNCellsAroundPtBin[iPt]){
                    meanDataNCellsAround[iPt]           = histoDataMergedNCellsAroundPtBin[iPt]->GetMean();
                    histoDataNCellsAroundChi2Test       = CorrectBinContentForChi2Test(histoDataMergedNCellsAroundPtBin[iPt]);
                }
                scaleFacData                            = histoDataMergedNCellsAroundPtBin[iPt]->GetMaximum();
                if (scaleFacData > 0){
                    histoDataMergedNCellsAroundPtBin[iPt]->Scale(1./scaleFacData);
                    if (histoDataNCellsAroundChi2Test) histoDataNCellsAroundChi2Test->Scale(1./scaleFacData);
                }    
            }
            
            TString fNameMCHistoMergedNCellsAround      = Form("NCellsAroundMerged_in_Pt_Bin%02d", iPt);
            if (histoMCMergedNCellsAroundPt){
                histoMCMergedNCellsAroundPtBin[iPt]     = FillProjectionX(histoMCMergedNCellsAroundPt, fNameMCHistoMergedNCellsAround, fBinsPt[iPt], fBinsPt[iPt+1], 1);
                if (histoMCMergedNCellsAroundPtBin[iPt]){
                    meanMCNCellsAround[iPt]             = histoMCMergedNCellsAroundPtBin[iPt]->GetMean();
                    histoMCNCellsAroundChi2Test         = CorrectBinContentForChi2Test(histoMCMergedNCellsAroundPtBin[iPt]);
                }    
                scaleFacMC                              = histoMCMergedNCellsAroundPtBin[iPt]->GetMaximum();
                if (scaleFacMC > 0){
                    histoMCMergedNCellsAroundPtBin[iPt]->Scale(1./scaleFacMC);    
                    if (histoMCNCellsAroundChi2Test) histoMCNCellsAroundChi2Test->Scale(1./scaleFacMC);
                }                    
            }
            if (histoDataNCellsAroundChi2Test && histoMCNCellsAroundChi2Test){
                 chi2NDFNCellsAround[iPt]               = histoDataNCellsAroundChi2Test->Chi2Test(histoMCNCellsAroundChi2Test,"NORM,CHI2/NDF");   
            }    

            // histos for NCells around and in merged cluster delta R < 0.15
            TH1D* histoDataNCellsAroundAndInChi2Test        = NULL;
            TH1D* histoMCNCellsAroundAndInChi2Test          = NULL;
            
            TString fNameDataHistoMergedNCellsAroundAndIn   = Form("NCellsAroundAndInMerged_in_Pt_Bin%02d", iPt);
            if (histoDataMergedNCellsAroundAndInPt){
                histoDataMergedNCellsAroundAndInPtBin[iPt]  = FillProjectionX(histoDataMergedNCellsAroundAndInPt, fNameDataHistoMergedNCellsAroundAndIn, fBinsPt[iPt], fBinsPt[iPt+1], 1);
                if (histoDataMergedNCellsAroundAndInPtBin[iPt]){
                    meanDataNCellsAroundAndIn[iPt]          = histoDataMergedNCellsAroundAndInPtBin[iPt]->GetMean();
                    histoDataNCellsAroundAndInChi2Test      = CorrectBinContentForChi2Test(histoDataMergedNCellsAroundAndInPtBin[iPt]);
                }
                scaleFacData                                = histoDataMergedNCellsAroundAndInPtBin[iPt]->Integral();
                if (scaleFacData > 0){
                    histoDataMergedNCellsAroundAndInPtBin[iPt]->Scale(1./scaleFacData);
                    if (histoDataNCellsAroundAndInChi2Test) histoDataNCellsAroundAndInChi2Test->Scale(1./scaleFacData);
                }    

            }
            TString fNameMCHistoMergedNCellsAroundAndIn     = Form("NCellsAroundAndInMerged_in_Pt_Bin%02d", iPt);
            if (histoMCMergedNCellsAroundAndInPt){
                histoMCMergedNCellsAroundAndInPtBin[iPt]    = FillProjectionX(histoMCMergedNCellsAroundAndInPt, fNameMCHistoMergedNCellsAroundAndIn, fBinsPt[iPt], fBinsPt[iPt+1], 1);
                if (histoMCMergedNCellsAroundAndInPtBin[iPt]){
                    meanMCNCellsAroundAndIn[iPt]            = histoMCMergedNCellsAroundAndInPtBin[iPt]->GetMean();
                    histoMCNCellsAroundAndInChi2Test        = CorrectBinContentForChi2Test(histoMCMergedNCellsAroundAndInPtBin[iPt]);
                }    
                scaleFacMC                                  = histoMCMergedNCellsAroundAndInPtBin[iPt]->Integral();
                if (scaleFacMC > 0){
                    histoMCMergedNCellsAroundAndInPtBin[iPt]->Scale(1./scaleFacMC);    
                    if (histoMCNCellsAroundAndInChi2Test) histoMCNCellsAroundAndInChi2Test->Scale(1./scaleFacMC);
                }                    
            }
            if (histoDataNCellsAroundAndInChi2Test && histoMCNCellsAroundAndInChi2Test){
                 chi2NDFNCellsAroundAndIn[iPt]              = histoDataNCellsAroundAndInChi2Test->Chi2Test(histoMCNCellsAroundAndInChi2Test,"NORM,CHI2/NDF");   
            }    

            
            // histos for E around merged cluster delta R < 0.15
            TH1D* histoDataEAroundChi2Test          = NULL;
            TH1D* histoMCEAroundChi2Test            = NULL;

            TString fNameDataHistoMergedEAroundE    = Form("EAroundEMerged_in_Pt_Bin%02d", iPt);
            if (histoDataMergedEAroundE){
                histoDataMergedEAroundEBin[iPt]     = FillProjectionX(histoDataMergedEAroundE, fNameDataHistoMergedEAroundE, fBinsPt[iPt], fBinsPt[iPt+1], 1);
                histoDataMergedEAroundEBin[iPt]->Rebin(5);
                if (histoDataMergedEAroundEBin[iPt]){
                    meanDataEAround[iPt]            = histoDataMergedEAroundEBin[iPt]->GetMean();
                    histoDataEAroundChi2Test        = CorrectBinContentForChi2Test(histoDataMergedEAroundEBin[iPt]);
                }
                scaleFacData                        = histoDataMergedEAroundEBin[iPt]->GetMaximum();
                if (scaleFacData > 0){
                    histoDataMergedEAroundEBin[iPt]->Scale(1./scaleFacData);
                    if (histoDataEAroundChi2Test) histoDataEAroundChi2Test->Scale(1./scaleFacData);
                }    
            }
            
            TString fNameMCHistoMergedEAroundE      = Form("EAroundMerged_in_Pt_Bin%02d", iPt);
            if (histoMCMergedEAroundE){
                histoMCMergedEAroundEBin[iPt]       = FillProjectionX(histoMCMergedEAroundE, fNameMCHistoMergedEAroundE, fBinsPt[iPt], fBinsPt[iPt+1], 1);
                histoMCMergedEAroundEBin[iPt]->Rebin(5);
                if (histoMCMergedEAroundEBin[iPt]){
                    meanMCEAround[iPt]              = histoMCMergedEAroundEBin[iPt]->GetMean();
                    histoMCEAroundChi2Test          = CorrectBinContentForChi2Test(histoMCMergedEAroundEBin[iPt]);
                }
                scaleFacMC                          = histoMCMergedEAroundEBin[iPt]->GetMaximum();
                if (scaleFacMC > 0){
                    histoMCMergedEAroundEBin[iPt]->Scale(1./scaleFacMC);    
                    if (histoMCEAroundChi2Test) histoMCEAroundChi2Test->Scale(1./scaleFacMC);
                }    
            }
            if (histoDataEAroundChi2Test && histoMCEAroundChi2Test){
                 chi2NDFEAround[iPt]                = histoDataEAroundChi2Test->Chi2Test(histoMCEAroundChi2Test,"NORM,CHI2/NDF");   
            }    

        }
        
        if (histoDataMergedNCellsPt && histoMCMergedNCellsPt) {
            TString nameCompPlotNCellsPlot          = Form("%s/%s_MergedClusterNCellsCompared_%s.%s",outputDir.Data(),mesonType.Data(),fCutSelection.Data(),fSuffix.Data()) ;
            PlotMergedCompDataMCInPtBins(   histoDataMergedNCellsPtBin, histoMCMergedNCellsPtBin, NULL, chi2NDFNCells, meanDataNCells, meanMCNCells,
                                            nameCompPlotNCellsPlot,  "bla", "bla", fMergedNCellsRange, 
                                            fdate, mesonType, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fDecayChannel, fDetectionProcessPtBins, fCollisionSystem, 3);
            nameCompPlotNCellsPlot          = Form("%s/%s_MergedClusterNCellsComparedSingleBin1.%s",outputDir.Data(),mesonType.Data(),fSuffix.Data()) ;
            PlotMergedCompDataMCSingleBin(  histoDataMergedNCellsPtBin, histoMCMergedNCellsPtBin, NULL, chi2NDFNCells, meanDataNCells, meanMCNCells,
                                            nameCompPlotNCellsPlot,  "bla", "bla", fMergedNCellsAroundRange, 
                                            fdate, mesonType, fExampleBin, fBinsPt, mode, fCollisionSystem, fTrigger,1); 
            nameCompPlotNCellsPlot          = Form("%s/%s_MergedClusterNCellsComparedSingleBin2.%s",outputDir.Data(),mesonType.Data(),fSuffix.Data()) ;
            PlotMergedCompDataMCSingleBin(  histoDataMergedNCellsPtBin, histoMCMergedNCellsPtBin, NULL, chi2NDFNCells, meanDataNCells, meanMCNCells,
                                            nameCompPlotNCellsPlot,  "bla", "bla", fMergedNCellsAroundRange, 
                                            fdate, mesonType, fExampleBin, fBinsPt+2, mode, fCollisionSystem, fTrigger,1); 
        }
        if (histoDataMergedNCellsAroundPt && histoMCMergedNCellsAroundPt) {
            TString nameCompPlotNCellsAroundPlot    = Form("%s/%s_MergedClusterNCellsAroundClusterCompared_%s.%s",outputDir.Data(),mesonType.Data(),fCutSelection.Data(),fSuffix.Data()) ;
            PlotMergedCompDataMCInPtBins(   histoDataMergedNCellsAroundPtBin, histoMCMergedNCellsAroundPtBin, NULL, chi2NDFNCellsAround, meanDataNCellsAround, meanMCNCellsAround,
                                            nameCompPlotNCellsAroundPlot,  "bla", "bla", fMergedNCellsAroundRange, 
                                            fdate, mesonType, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fDecayChannel, fDetectionProcessPtBins, fCollisionSystem, 3);
            nameCompPlotNCellsAroundPlot    = Form("%s/%s_MergedClusterNCellsAroundClusterComparedSingleBin.%s",outputDir.Data(),mesonType.Data(),fSuffix.Data()) ;
            PlotMergedCompDataMCSingleBin(  histoDataMergedNCellsAroundPtBin, histoMCMergedNCellsAroundPtBin, NULL, chi2NDFNCellsAround, meanDataNCellsAround, meanMCNCellsAround,
                                            nameCompPlotNCellsAroundPlot,  "bla", "bla", fMergedNCellsAroundRange, 
                                            fdate, mesonType, fExampleBin, fBinsPt, mode, fCollisionSystem, fTrigger); 
        }
        if (histoDataMergedNCellsAroundAndInPt && histoMCMergedNCellsAroundAndInPt) {
            TString nameCompPlotNCellsAroundAndInPlot   = Form("%s/%s_MergedClusterNCellsAroundAndInClusterCompared_%s.%s",outputDir.Data(),mesonType.Data(),fCutSelection.Data(),fSuffix.Data()) ;
            PlotMergedCompDataMCInPtBins(   histoDataMergedNCellsAroundAndInPtBin, histoMCMergedNCellsAroundAndInPtBin, NULL, chi2NDFNCellsAroundAndIn, meanDataNCellsAroundAndIn, meanMCNCellsAroundAndIn,
                                            nameCompPlotNCellsAroundAndInPlot,  "bla", "bla", fMergedNCellsAroundAndInRange, 
                                            fdate, mesonType, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fDecayChannel, fDetectionProcessPtBins, fCollisionSystem, 3);
            nameCompPlotNCellsAroundAndInPlot   = Form("%s/%s_MergedClusterNCellsAroundAndInClusterComparedSingleBin.%s",outputDir.Data(),mesonType.Data(),fSuffix.Data()) ;
            PlotMergedCompDataMCSingleBin( histoDataMergedNCellsAroundAndInPtBin, histoMCMergedNCellsAroundAndInPtBin, NULL, chi2NDFNCellsAroundAndIn, meanDataNCellsAroundAndIn, meanMCNCellsAroundAndIn,
                                            nameCompPlotNCellsAroundAndInPlot,  "bla", "bla", fMergedNCellsAroundAndInRange, 
                                            fdate, mesonType, fExampleBin, fBinsPt, mode, fCollisionSystem, fTrigger); 
        }                            

        
        if (histoDataMergedEAroundE && histoMCMergedEAroundE) {
            TString nameCompPlotEAroundPlot = Form("%s/%s_MergedClusterEAroundClusterCompared_%s.%s",outputDir.Data(),mesonType.Data(),fCutSelection.Data(),fSuffix.Data()) ;
            for (Int_t i= fStartPtBin; i< fNBinsPt; i++){
                cout << histoDataMergedEAroundEBin[i]->GetBinContent(histoDataMergedEAroundEBin[i]->FindBin(12.0)) << endl;
                histoDataMergedEAroundEBin[i]->GetYaxis()->SetRangeUser(-0.001,histoMCMergedEAroundEBin[i]->GetBinContent(histoMCMergedEAroundEBin[i]->FindBin(10.5)*10));
            }    
            PlotMergedCompDataMCInPtBins(   histoDataMergedEAroundEBin, histoMCMergedEAroundEBin, NULL, chi2NDFEAround, meanDataEAround, meanMCEAround,
                                            nameCompPlotEAroundPlot,  "bla", "bla", fMergedERange, 
                                            fdate, mesonType, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fDecayChannel, fDetectionProcessPtBins, fCollisionSystem, 2);
            nameCompPlotEAroundPlot         = Form("%s/%s_MergedClusterEAroundClusterComparedSingleBin.%s",outputDir.Data(),mesonType.Data(),fSuffix.Data()) ;
            PlotMergedCompDataMCSingleBin(  histoDataMergedEAroundEBin, histoMCMergedEAroundEBin, NULL, chi2NDFEAround, meanDataEAround, meanMCEAround,
                                            nameCompPlotEAroundPlot,  "bla", "bla", fMergedERange, 
                                            fdate, mesonType, fExampleBin, fBinsPt, mode, fCollisionSystem, fTrigger, 2); 
        }
    }    
}



