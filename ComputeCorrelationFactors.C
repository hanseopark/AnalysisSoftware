/****************************************************************************************************************************
 ******     provided by Gamma Conversion Group, PWGGA,                                                                  *****
 ******        Daniel Mühlheim, d.muehlheim@cern.ch                                                                     *****
 *****************************************************************************************************************************/
#include <TString.h>
#include <TObject.h>
#include <TH1.h>
#include <TH2.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TGaxis.h>
#include <TLegend.h>
#include <TFrame.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TF1.h>
#include <TASImage.h>
#include <TMarker.h>
#include <TDatime.h>
#include <TArrow.h>
#include <TROOT.h>
#include <TObjString.h>
#include <TFile.h>
#include "TMath.h"


#include <Riostream.h>
#include <string>
#include <sstream>
#include <fstream>
#include <stdlib.h>
#include <fstream>
#include <math.h>

#include "CommonHeaders/PlottingGammaConversionHistos.h"
#include "CommonHeaders/PlottingGammaConversionAdditional.h"
#include "CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"

using namespace std;

Int_t GetPosInVec( std::vector<TString> vec, TString lookup){
    for(Int_t i=0; i<(Int_t)vec.size(); i++){
        if(vec.at(i).CompareTo(lookup) == 0) return i;
    }
    cout << "ERROR in GetPosInVec, did not find '" << lookup.Data() << "'" << endl;
    return -1;
}

void ComputeCorrelationFactors(
                                    TString fileInput   = "input.txt",
                                    TString combMode    = "triggers",
                                    TString meson       = "",
                                    TString energy      = "",
                                    Int_t mode          = -1,
                                    TString suffix      = "eps",
                                    Bool_t isStatCorr   = kFALSE,
                                    TString centrality  = "",
                                    TString eventCut    = ""
                                   ){

    //---------------------------------------------------------------------------------------------------------------
    //---------------------------- General setting & global variables -----------------------------------------------
    //---------------------------------------------------------------------------------------------------------------
    gROOT->Reset();
    gROOT->SetStyle("Plain");
    TH1::AddDirectory(kFALSE);

    StyleSettingsThesis();
    SetPlotStyle();

    TString mesonPlot       = "";
    if(meson.CompareTo("Pi0") == 0 || meson.CompareTo("Pi0RpPb") == 0)
        mesonPlot           = "#pi^{0}";
    else if(meson.CompareTo("Eta") == 0 || meson.CompareTo("EtaRpPb") == 0)
        mesonPlot           = "#eta";
    else if(meson.CompareTo("Pi0EtaBinning") == 0 || meson.CompareTo("EtaToPi0") == 0)
        mesonPlot           = "#eta/#pi^{0}";
    else if(meson.CompareTo("Omega") == 0)
        mesonPlot           = "#omega";
    else if(meson.CompareTo("GammaInc") == 0 )
        mesonPlot           = "#gamma_{inc}";
    else if(meson.CompareTo("IncGammaToPi0") == 0 )
        mesonPlot           = "#gamma_{inc}/#pi^{0}";
    else if(meson.CompareTo("RGamma") == 0 )
        mesonPlot           = "R_{#gamma}";

    TString modeOutput      = Form("%d",mode);
    if(combMode.CompareTo("systems")==0)
        modeOutput                = "Systems";

    if (isStatCorr) mesonPlot       = "stat. corr "+mesonPlot;
    else mesonPlot       = "sys. corr "+mesonPlot;

    TString dateForOutput           = ReturnDateStringForOutput();
    TString fCollisionSystemWrite   = ReturnCollisionEnergyOutputString(energy);
    TString fCollisionSystemAndCent = fCollisionSystemWrite;
    TString centralityForOutput     = GetCentralityStringOutput(eventCut);
    if (energy.CompareTo("pPb_5.023TeVCent") == 0)
        centralityForOutput         = centralityForOutput+"Cent";
    TString collisionSystem         = ReturnFullCollisionsSystem(energy);
    if (centrality.CompareTo("") != 0){
        collisionSystem             = centrality+" "+collisionSystem;
        fCollisionSystemAndCent     = centralityForOutput+"_"+fCollisionSystemWrite;
    }
    TString outputDir               = Form("%s/%s/ComputeCorrelationFactors_%s",suffix.Data(),dateForOutput.Data(), fCollisionSystemWrite.Data());
    TString detectionProcess        = ReturnFullTextReconstructionProcess(mode);

    fstream fLog;
    fLog.open(Form("%s/%s_%s_%s_mode%s.log",outputDir.Data(),meson.Data(),combMode.Data(),fCollisionSystemAndCent.Data(),modeOutput.Data()), ios::out);
    fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    fLog << "Energy: " << fCollisionSystemAndCent.Data() << endl;
    fLog <<  collisionSystem.Data() << endl;
    fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

    cout << "creating outpur dir: " << outputDir.Data() << endl;
    gSystem->Exec("mkdir -p "+outputDir);
    gSystem->Exec("mkdir -p "+outputDir+"/sys");
    gSystem->Exec("mkdir -p "+outputDir+"/corrFactors");
    gSystem->Exec("mkdir -p "+outputDir+"/corrPlotting");
    gSystem->Exec(Form("cp %s %s/%s_%s_%s_mode%s.config", fileInput.Data(), outputDir.Data(), meson.Data(), combMode.Data(), fCollisionSystemAndCent.Data(), modeOutput.Data()));

    cout << "running combination mode: " << combMode.Data() << endl;
    fLog << "running combination mode: " << combMode.Data() << endl;
    if(combMode.CompareTo("triggers") && combMode.CompareTo("systems")){ cout << "ERROR in line '" << __LINE__ << "' combMode can only be 'triggers' or 'systems', returning..."<< endl; return;}

    Int_t nMeasTot     = 0;            // number of measurements
    vector<TString>  vecComb;       // vector with names of measurements
    vector<TString>  vecSysFiles;   // vector with names of sytematics files for measurements
    vector<Int_t>    vecNBinsPt;      // vector with number of pT bins for measurements
    vector<Double_t> vecUseMinPt;   // vector with minimum pT for measurements
    vector<Double_t> vecUseMaxPt;   // vector with maximum pT for measurements

    //---------------------------------------------------------------------------------------------------------------
    //--------------------------------- read config file
    //---------------------------------------------------------------------------------------------------------------
    cout << "reading input file" << endl;

    // reading how many measurements should be considered
    ifstream in(fileInput.Data());
    in >> nMeasTot;
    cout << Form("number of %s available: ",combMode.Data()) << nMeasTot << "\n" << endl;
    fLog << Form("number of %s available: ",combMode.Data()) << nMeasTot << "\n" << endl;
    if(nMeasTot<1){ cout << "ERROR in line '" << __LINE__ << "', returning..."<< endl; return;}

    // reading measurement considered for correlations
    for(Int_t iMeas=0; iMeas<nMeasTot; iMeas++) {
        TString temp            = "";
        TString tempSysFile     = "";
        Int_t nBinsTemp         = 0;
        Double_t minTemp        = 0;
        Double_t maxTemp        = 0;
        // format in file should be:
        // name meas \t number of bins \t minimum pt \t maximum pt \t sytematics file \n
        in >> temp >> nBinsTemp >> minTemp >> maxTemp >> tempSysFile;
        if(!temp.IsNull() && !tempSysFile.IsNull() && nBinsTemp > 0 && minTemp > 0 && maxTemp > 0){
            vecComb.push_back(temp);
            vecSysFiles.push_back(tempSysFile);
            vecNBinsPt.push_back(nBinsTemp);
            vecUseMinPt.push_back(minTemp);
            vecUseMaxPt.push_back(maxTemp);
            TObjArray *rToken     = tempSysFile.Tokenize("/");
            TObjString *rString   = (TObjString*)rToken->At(rToken->GetLast());
            gSystem->Exec(Form("cp %s %s/sys/%s", tempSysFile.Data(), outputDir.Data(), ((TString)rString->GetString()).Data()));
            fLog << "read: " << temp.Data() << ", " << minTemp << ", " << maxTemp << ", " << tempSysFile.Data() << endl;
        }else{ cout << "ERROR in line '" << __LINE__ << "', returning..."<< endl; return;}
    }
    fLog << endl;


    vector<TString>** corrName            = new vector<TString>*[nMeasTot];      // nMeasTot x nMeasTot matrix with names for the measurements
    vector<TString>** corr                = new vector<TString>*[nMeasTot];      // nMeasTot x nMeasTot matrix with correlated fraction for the measurements as strings
    vector<Double_t> **corrFactorsBins    = new vector<Double_t>*[nMeasTot];     // nMeasTot x nMeasTot matrix with pt values
    vector<Double_t> **corrFactors        = new vector<Double_t>*[nMeasTot];     // nMeasTot x nMeasTot matrix with rho_AB's
    for(Int_t iMeas=0; iMeas<nMeasTot; iMeas++){
        corrName[iMeas]        = new vector<TString>[nMeasTot];
        corr[iMeas]            = new vector<TString>[nMeasTot];
        corrFactorsBins[iMeas] = new vector<Double_t>[nMeasTot];
        corrFactors[iMeas]     = new vector<Double_t>[nMeasTot];
    }

    // reading assumed correlation settings
    while(!in.eof()){
        TString tempA = ""; TString tempB = ""; TString tempC = ""; TString tempD = "";
        // format in the file should be:
        // name meas A \t name meas B \t name systematics component X \t uncorrelated error X in meas A with respect to B \n
        // if last value given with % at the end:
        //      => will calculate relative error for each pt bin
        //      => ranges from 0% - 100% , 0% meaning they are treated as fully correlated, 100% meaning they are treated as fully uncorrelated,
        // if last value given without % at the end:
        //      => absolute error according to that number is taken out of total error for meas A with respect to B
        in >> tempA >> tempB >> tempC >> tempD;
        if(!tempA.IsNull() && !tempB.IsNull() && !tempC.IsNull() && !tempD.IsNull()){
            corrName[GetPosInVec(vecComb,tempA)][GetPosInVec(vecComb,tempB)].push_back(tempC);
            corr    [GetPosInVec(vecComb,tempA)][GetPosInVec(vecComb,tempB)].push_back(tempD);
        }
    }
    in.close();

    //---------------------------------------------------------------------------------------------------------------
    //--------------------------------- read detailed systematics for each input
    //---------------------------------------------------------------------------------------------------------------
    // nMeasTot x nPtbins x nError sources matrix of names of systematics
    // ptSys[iMeas][iPt].push_back(iESource)
    // if (iPt > 0 && iESource == 0) ptSys[iMeas][iPt].push_back(0)         = pTValue for iMeas's measurement
    // if (iPt == 0 && iESource > 0) ptSys[iMeas][0].push_back(iESource)     = name of iESource's error for iMeas's measurement
    // if (iPt > 0 && iESource > 0)  ptSys[iMeas][iPt].push_back(iESource)   = error for iMeas's measurement for iESource's  at ptSys[iMeas][iPt].push_back(0)
    vector<TString>** ptSys   = new vector<TString>*[nMeasTot];

    for(Int_t iMeas=0; iMeas<nMeasTot; iMeas++)
        ptSys[iMeas]             = new vector<TString>[vecNBinsPt.at(iMeas)+1];

    // loop over all measurements with respective systematics files
    for(Int_t iMeas=0; iMeas<nMeasTot; iMeas++){
        ifstream fileSysErr;
        fileSysErr.open(vecSysFiles.at(iMeas).Data(),ios_base::in);
        cout << "opening: " << vecSysFiles.at(iMeas).Data() << endl;
        fLog << "opening: " << vecSysFiles.at(iMeas).Data() << endl;
        Int_t iPt = 0;
        string line;
        while (getline(fileSysErr, line) && iPt < 100) {
            istringstream ss(line);
            TString temp        = "";
            TString tempBin     = "";
            Int_t iESource = 0;
            fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
            fLog << "\t*** reading line " << iPt << " ***" << endl;
            while(ss && iESource < 100){
                ss >> temp;
                if(iPt > 0 && iESource == 0){
                    tempBin = temp;
                    if( !tempBin.IsNull() && (tempBin.Atof() < vecUseMinPt.at(iMeas) || tempBin.Atof() > vecUseMaxPt.at(iMeas)) ){
                        cout << "INFO: skipping pT " << tempBin.Atof() << ", outside of specified interval of " << vecUseMinPt.at(iMeas) << "-" << vecUseMaxPt.at(iMeas) << endl;
                        break;
                    }
                }
                if( !(iPt==0 && temp.CompareTo("bin")==0) && !temp.IsNull() && iPt < (vecNBinsPt.at(iMeas)+1)){
                    ptSys[iMeas][iPt].push_back(temp);
                    fLog << temp.Data() << ", ";
                    iESource++;
                }
            }

            if(iPt == 0){
                ptSys[iMeas][iPt++].push_back("TotalError");
                fLog << "TotalError";
                iESource++;
            } else if(!tempBin.IsNull() && (tempBin.Atof() >= vecUseMinPt.at(iMeas) && tempBin.Atof() <= vecUseMaxPt.at(iMeas)))
                iPt++;
            fLog << "\n\t***" << iESource << " columns read ***" << endl;
            fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n" << endl;
        }
        vecNBinsPt.at(iMeas) = --iPt;
        cout << "\t---" << vecNBinsPt.at(iMeas)+1 << " pT bins read ---" << endl;
        fLog << "\t---" << vecNBinsPt.at(iMeas)+1 << " pT bins read ---" << endl;
        fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        if(iPt < 1){
            cout << "ERROR: read " << vecNBinsPt.at(iMeas)+1<< " pT-bins, returning!" << endl;
            fLog << "ERROR: read " << vecNBinsPt.at(iMeas)+1 << " pT-bins, returning!" << endl;
            return;
        }
        fileSysErr.close();
    }

    //---------------------------------------------------------------------------------------------------------------
    //-------------------------------- Define variable for plotting -------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------
    Int_t textSizeLabelsPixel   = 900*0.04;

    // canvas
    TCanvas* canvasWeights      = new TCanvas("canvasWeights","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasWeights, 0.07, 0.005, 0.05, 0.09);
    canvasWeights->SetLogx();

    Int_t plotErr               = 0;
    TH2F * histo2DPi0Weights    = NULL;
    Double_t minCorrYaxis       = 0.45;
    Double_t maxCorrYaxis       = 1.05;
    if(combMode.CompareTo("systems") == 0 && energy.Contains("pPb_5.023TeV")  && !meson.Contains("Gamma")){
        minCorrYaxis            = 0.0;
    } else if(combMode.CompareTo("systems") == 0 && (energy.Contains("pPb_5.023TeV")  || energy.CompareTo("2.76TeV") == 0 || energy.CompareTo("8TeV") == 0) &&  meson.Contains("Gamma") && !isStatCorr){
        minCorrYaxis            = 0.005;
        maxCorrYaxis            = 1.5;
        canvasWeights->SetLogy();
    } else if(  combMode.CompareTo("systems") == 0 && (energy.Contains("pPb_5.023TeV")  || energy.CompareTo("2.76TeV") == 0 || energy.CompareTo("8TeV") == 0) &&
                (( meson.Contains("RGamma") ||  meson.Contains("IncGammaToPi0")) && isStatCorr ) ){
        minCorrYaxis            = 0.005;
        maxCorrYaxis            = 1.05;
    }else if(combMode.CompareTo("systems") == 0 ){
        minCorrYaxis            = 0.01;
        if(mode >= 40 && mode <= 45){
          minCorrYaxis = 0.2;
          maxCorrYaxis = 0.82;
        }
    } else if (combMode.CompareTo("triggers") == 0 && (energy.CompareTo("2.76TeV") == 0  || energy.Contains("5TeV2017") )&& ( mode == 2 || mode == 4 ) ){
        minCorrYaxis            = 0.15;
    }
    Double_t maxPt              = 0;
    for (Int_t iMeas = 0; iMeas < nMeasTot; iMeas++){
        if (maxPt < vecUseMaxPt.at(iMeas))
            maxPt       = vecUseMaxPt.at(iMeas);
    }
    maxPt                       = maxPt*1.2;

    Double_t minPt              = 10000;
    for (Int_t iMeas = 0; iMeas < nMeasTot; iMeas++){
        if (minPt > vecUseMinPt.at(iMeas))
            minPt               = vecUseMinPt.at(iMeas);
    }

    // dummy hist
    histo2DPi0Weights           = new TH2F("histo2DPi0Weights","histo2DPi0Weights",11000,minPt,maxPt,1000,minCorrYaxis,maxCorrYaxis);
    SetStyleHistoTH2ForGraphs(histo2DPi0Weights, "#it{p}_{T} (GeV/#it{c})","#rho_{ij}",0.04,0.045, 0.04,0.045, 0.85,0.73);
    histo2DPi0Weights->GetXaxis()->SetMoreLogLabels();
    histo2DPi0Weights->GetXaxis()->SetLabelOffset(-0.01);
    histo2DPi0Weights->Draw("copy");

    Int_t nRowsLabels           = 2;
    if (modeOutput.CompareTo("Systems") != 0)
        nRowsLabels++;

    // labels
    TLatex *labelWeightsEnergy  = new TLatex(0.95,0.14+((nRowsLabels-1)*0.05),collisionSystem.Data());
    SetStyleTLatex( labelWeightsEnergy, textSizeLabelsPixel,4);
    labelWeightsEnergy->SetTextFont(43);
    labelWeightsEnergy->SetTextAlign(31);
    labelWeightsEnergy->Draw();

    TLatex *labelWeightsPi0     = 0x0;
    if(meson.CompareTo("Pi0EtaBinning") == 0 || meson.CompareTo("EtaToPi0") == 0 || meson.Contains("Gamma") )
        labelWeightsPi0         = new TLatex(0.95,0.14+((nRowsLabels-2)*0.05),Form("%s",mesonPlot.Data()));
    else if(mode >= 40 && mode <= 45)
        labelWeightsPi0         = new TLatex(0.95,0.14+((nRowsLabels-2)*0.05),Form("%s #rightarrow #pi^{+}#pi^{-}#pi^{0}",mesonPlot.Data()));
    else
        labelWeightsPi0         = new TLatex(0.95,0.14+((nRowsLabels-2)*0.05),Form("%s #rightarrow #gamma#gamma",mesonPlot.Data()));
    SetStyleTLatex( labelWeightsPi0, textSizeLabelsPixel,4);
    labelWeightsPi0->SetTextFont(43);
    labelWeightsPi0->SetTextAlign(31);
    labelWeightsPi0->Draw();

    TLatex *labelWeightsDetectionProcess  = new TLatex(0.95,0.14,detectionProcess.Data());
    SetStyleTLatex( labelWeightsDetectionProcess, textSizeLabelsPixel,4);
    labelWeightsDetectionProcess->SetTextFont(43);
    labelWeightsDetectionProcess->SetTextAlign(31);
    if (modeOutput.CompareTo("Systems") != 0)
        labelWeightsDetectionProcess->Draw();


    // colors, markers
    Color_t colorTrigg[30]      = { kBlack, kGray+1, kRed+2, kBlue+2, kGreen+3, kCyan+2, kViolet+1, kMagenta+2,  kRed-2, kBlue-2,
                                    807, kAzure+2, kGreen-2, kMagenta-2, kYellow+2, kAzure-5, kPink-2, kCyan-5,  kRed-2, kBlue-2,
                                    kBlack, kGray+1, kRed+2, kBlue+2, kGreen+3, kCyan+2, kViolet, kMagenta+2,  kRed-2, kBlue-2 };
    Marker_t markerTrigg[30]    = { 20, 24, 21, 34, 29, 33, 21, 27, 28, 30,
                                    24, 25, 20, 21, 34, 29, 33, 21, 27, 28,
                                    30, 24, 25, 20, 21, 34, 29, 33, 21, 27  };
    Int_t iTrigg                = 0;
    Int_t iColumnsLegend        = 1;

    // legend
    Int_t maxNLegendColumn      = 4;
    TLegend* legendWeights      = GetAndSetLegend2(0.11, 0.11, 0.42, 0.1+(0.04*maxNLegendColumn), textSizeLabelsPixel);
    legendWeights->SetMargin(0.15);
    TLegend* legendWeightsTop   = GetAndSetLegend2(0.055, 0.95, 0.995, 0.995, textSizeLabelsPixel, nMeasTot+1, "", 43, 0.02);
    legendWeightsTop->AddEntry((TObject*)0, "i,j:","");
    //---------------------------------------------------------------------------------------------------------------
    //--------------------------------- calculate correlation factors
    //---------------------------------------------------------------------------------------------------------------
    fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n" << endl;

    fstream fCorr;
    if (!isStatCorr)
        fCorr.open(Form("%s/corrFactors/%s_%s_%s_mode%s.log",outputDir.Data(),meson.Data(),combMode.Data(),fCollisionSystemAndCent.Data(),modeOutput.Data()), ios::out);
    else
        fCorr.open(Form("%s/corrFactors/%s_%s_Stat_%s_mode%s.log",outputDir.Data(),meson.Data(),combMode.Data(),fCollisionSystemAndCent.Data(),modeOutput.Data()), ios::out);
    TFile *fOutput              = new TFile(Form("%s/%s.root",outputDir.Data(),fCollisionSystemWrite.Data()),"UPDATE");
    if (centrality.CompareTo("") != 0){
        TDirectoryFile* directoryCentOut  = NULL;
        directoryCentOut                  = (TDirectoryFile*)fOutput->Get(centralityForOutput.Data());

        if (!directoryCentOut){
            fOutput->mkdir(centralityForOutput.Data());
            directoryCentOut              = (TDirectoryFile*)fOutput->Get(centralityForOutput.Data());
        }
        fOutput->cd(centralityForOutput.Data());
    }
    // Calculate rho_AB for A = iMeasA and B = iMeasB
    for(Int_t iMeasA=0; iMeasA<nMeasTot; iMeasA++){
        for(Int_t iMeasB=0; iMeasB<nMeasTot; iMeasB++){
            // don't calculate it for the correlation with itself
            if (iMeasA == 0){
                legendWeightsTop->AddEntry((TObject*)0, Form("%d-%s",iMeasB,vecComb.at(iMeasB).Data()),"");
            }
            if(iMeasA==iMeasB) continue;
            Int_t commonBins    = 0;
            Int_t iPtA          = 1;
            Int_t iPtB          = 1;
            Double_t pTA         = -1;
            Double_t pTB        = -1;

            Double_t factor     = 0;
            TString tempCorr    = "";
            TString tempPlot    = "";
            if(iMeasA<iMeasB){
                tempCorr = Form("%s_%s-%s",vecComb.at(iMeasA).Data(),vecComb.at(iMeasA).Data(),vecComb.at(iMeasB).Data());
                tempPlot = Form("i=%d,j=%d",iMeasA,iMeasB);
                cout << "\n\n\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
                cout << "correlation factors " << vecComb.at(iMeasA) << "_" << vecComb.at(iMeasA) << "-" << vecComb.at(iMeasB) << endl;
                cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
                fLog << "\n\n\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
                fLog << "correlation factors " << vecComb.at(iMeasA) << "_" << vecComb.at(iMeasA) << "-" << vecComb.at(iMeasB) << endl;
                fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
            }else{
                tempCorr = Form("%s_%s-%s",vecComb.at(iMeasA).Data(),vecComb.at(iMeasB).Data(),vecComb.at(iMeasA).Data());
                tempPlot = Form("i=%d,j=%d",iMeasA,iMeasB);
                cout << "\n\n\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
                cout << "correlation factors " << vecComb.at(iMeasA) << "_" << vecComb.at(iMeasB) << "-" << vecComb.at(iMeasA) << endl;
                cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
                fLog << "\n\n\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
                fLog << "correlation factors " << vecComb.at(iMeasA) << "_" << vecComb.at(iMeasB) << "-" << vecComb.at(iMeasA) << endl;
                fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
            }
            fCorr << endl;
            fCorr << tempCorr.Data() << endl;

            // create correlations histogram
            TH1F* histoCorr         = new TH1F(Form("%s_%s_%s",modeOutput.Data(),meson.Data(),tempCorr.Data()),Form("%s_%s_%s",modeOutput.Data(),meson.Data(),tempCorr.Data()),2000,-0.025,100-0.025);
            // set all bin 0 for correlations histogram
            for(Int_t iBin=1; iBin<histoCorr->GetNbinsX()+1; iBin++)
                histoCorr->SetBinContent(iBin,0.);

            for(;; iPtA++, iPtB++){
                // check if bin count goes above defined value in config for both measurements
                if( iPtA == vecNBinsPt.at(iMeasA)+1 || iPtB == vecNBinsPt.at(iMeasB)+1 ) break;

                // get pt values at beginning
                pTA     = ((TString)ptSys[iMeasA][iPtA].at(0)).Atof();
                pTB     = ((TString)ptSys[iMeasB][iPtB].at(0)).Atof();

                // find overlapping bin
                if( iPtA < vecNBinsPt.at(iMeasA) && iPtB < vecNBinsPt.at(iMeasB) ){
                    do {
                        // increase bin for measurement A if pT meas A < pT meas B
                        if(pTA < pTB){
                            if( iPtA == vecNBinsPt.at(iMeasA) ) break;
                            pTA = ((TString)ptSys[iMeasA][++iPtA].at(0)).Atof();
                        // increase bin for measurement B if pT meas B < pT meas C
                        }else if(pTB < pTA){
                            if( iPtB == vecNBinsPt.at(iMeasB) ) break;
                            pTB = ((TString)ptSys[iMeasB][++iPtB].at(0)).Atof();
                        }
                    // continue doing this until no bins are left
                    } while (pTA!=pTB && iPtA <= vecNBinsPt.at(iMeasA) && iPtB <= vecNBinsPt.at(iMeasB));
                }

                if( pTA!=pTB ) break;
                // increase number of common bins
                commonBins++;

                cout << "\n\t---------- pT: " << pTA << " ----------" << endl;
                fLog << "\n\t---------- pT: " << pTA << " ----------" << endl;

                if( (Int_t) corrName[iMeasA][iMeasB].size() > 0){
                    // get total error for measurement A at ptBin = iPtA
                    Double_t totalErr   = ((TString)ptSys[iMeasA][iPtA].at(GetPosInVec(ptSys[iMeasA][0],"TotalError"))).Atof();
                    cout << vecComb.at(iMeasA) << " totErr: " << totalErr << endl;
                    fLog << vecComb.at(iMeasA) << " totErr: " << totalErr << endl;
                    Double_t unCorrErr = 0;         // uncorrelated error with respect to measurement B

                    // run through the array given for the respective measurement combination
                    for(Int_t iSourceUnCorr=0; iSourceUnCorr<(Int_t) corrName[iMeasA][iMeasB].size(); iSourceUnCorr++){
                        TString tempName        = corrName[iMeasA][iMeasB].at(iSourceUnCorr);
                        TString tempD           = corr[iMeasA][iMeasB].at(iSourceUnCorr);
                        // get error for respective source of measurement A
                        Double_t tempUnCorrErr  = ((TString)ptSys[iMeasA][iPtA].at(GetPosInVec(ptSys[iMeasA][0],tempName))).Atof();
                        Double_t tempErr        = 0;
                        // if tempD contains % at end calculate relative error in pt based on fraction given
                        if(tempD.Contains("%")){
                            tempD.Resize(tempD.Length()-1);
                            tempErr     = tempUnCorrErr*tempD.Atof()/100.;
                        // if tempD does not contain %, take out given value from total error
                        } else
                            tempErr     = tempD.Atof();
                        unCorrErr += tempErr*tempErr;
                        cout << corrName[iMeasA][iMeasB].at(iSourceUnCorr) << "(" << tempUnCorrErr << "), unCorrErr: " << tempErr << endl;
                        fLog << corrName[iMeasA][iMeasB].at(iSourceUnCorr) << "(" << tempUnCorrErr << "), unCorrErr: " << tempErr << endl;
                    }
                    // calculate rho_AB at respective pTA
                    factor              = sqrt(totalErr*totalErr-unCorrErr)/totalErr;
                    if ( TMath::IsNaN(factor)) factor = 0;
                    if ( !TMath::Finite(factor)) factor = 0;
                }
                // set pt values and correlation factors rho_AB
                corrFactorsBins[iMeasA][iMeasB].push_back(pTA);
                corrFactors[iMeasA][iMeasB].push_back(factor);

                cout << "\n\t|-----------" << endl;
                cout << "\t|corrFactor:\t" << factor << "|" << endl;
                cout << "\t|------------" << endl;
                fLog << "\n\t|------------" << endl;
                fLog << "\n\tcorrFactor:\t" << factor << "|" << endl;
                fLog << "\t|------------" << endl;
                fCorr << pTA << " \t-\t " << factor << endl;
                histoCorr->SetBinContent(histoCorr->FindBin(pTA),factor);
                if(factor>0 && factor<=plotErr){
                    cout << "\n\n\tWARNING: point out of range (" << factor <<  ") for generated plots, please adjust y-range of histograms!!!\n\n" << endl;
                    fLog << "\n\n\tWARNING: point out of range (" << factor <<  ") for generated plots, please adjust y-range of histograms!!!\n\n" << endl;
                    plotErr++;
                }
            }
            cout << commonBins << endl;
            // plot if there are overlapping bins among the triggers and they are not completely uncorrelated
            if (commonBins > 0 && histoCorr->GetMaximum() > 0 ){
                histoCorr->SetMarkerSize(2.);
                histoCorr->SetMarkerColor(colorTrigg[iTrigg]);
                histoCorr->SetMarkerStyle(markerTrigg[iTrigg]);
                if (iTrigg > maxNLegendColumn)     iColumnsLegend  = 2;
                if (iTrigg > 2*maxNLegendColumn)   iColumnsLegend  = 3;
                if (iTrigg > 3*maxNLegendColumn)   iColumnsLegend  = 4;
                if (iTrigg > 4*maxNLegendColumn)   iColumnsLegend  = 5;
                iTrigg++;
                histoCorr->Draw("p,same");
                legendWeights->SetNColumns(iColumnsLegend);
                legendWeights->AddEntry(histoCorr,tempPlot.Data(),"p");
            } else {
               cout << "not plotting " << tempCorr.Data() << ", no correlation factors different from 0 found!" << endl;
            }
            // write to output
            if (!isStatCorr)
                histoCorr->Write(Form("%s_%s_%s",modeOutput.Data(),meson.Data(),tempCorr.Data()),TObject::kOverwrite);
            else
                histoCorr->Write(Form("%s_Stat_%s_%s",modeOutput.Data(),meson.Data(),tempCorr.Data()),TObject::kOverwrite);
        }
    }

    legendWeights->SetX2(0.11+0.10*iColumnsLegend);
    if(iTrigg+1 < maxNLegendColumn)
        legendWeights->SetY2(0.1+(0.035*iTrigg));
    legendWeightsTop->Draw();

    legendWeights->Draw();
    canvasWeights->RedrawAxis();

    labelWeightsEnergy->Draw();
    labelWeightsPi0->Draw();

    DrawGammaLines(minPt,maxPt, 1., 1.,0.1, kGray, 1);
    if (!meson.Contains("Gamma")){
        DrawGammaLines(minPt,maxPt, 0.98, 0.98,0.1, kGray, 3);
        DrawGammaLines(minPt,maxPt, 0.96, 0.96,0.1, kGray, 7);
        DrawGammaLines(minPt,maxPt, 0.94, 0.94,0.1, kGray, 3);
        DrawGammaLines(minPt,maxPt, 0.92, 0.92,0.1, kGray, 7);
    } else if (!isStatCorr) {
        DrawGammaLines(minPt,maxPt, 0.1, 0.1,0.1, kGray, 3);
    }
    DrawGammaLines(minPt,maxPt, 0.9, 0.9,0.1, kGray, 1);
    DrawGammaLines(minPt,maxPt, 0.8, 0.8,0.1, kGray, 3);
    DrawGammaLines(minPt,maxPt, 0.7, 0.7,0.1, kGray, 7);
    if(mode >= 40 && mode <= 45){
      DrawGammaLines(minPt,maxPt, 0.6, 0.6,0.1, kGray, 1);
      DrawGammaLines(minPt,maxPt, 0.5, 0.5,0.1, kGray, 3);
    }

    if (!isStatCorr)
        canvasWeights->SaveAs(Form("%s/corrPlotting/%s_%s_%s_corrFactors.%s",outputDir.Data(),fCollisionSystemAndCent.Data(),modeOutput.Data(),meson.Data(),suffix.Data()));
    else
        canvasWeights->SaveAs(Form("%s/corrPlotting/%s_%s_Stat_%s_corrFactors.%s",outputDir.Data(),fCollisionSystemAndCent.Data(),modeOutput.Data(),meson.Data(),suffix.Data()));

    if(plotErr>0){
        cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        cout << "\n\t|-----------" << endl;
        cout << "\t|total plotting errors:\t" << plotErr << "| -> adjust y-range in ComputeCorrelationFactors!" << endl;
        cout << "\t|------------" << endl;
        fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        fLog << "\n\t|------------" << endl;
        fLog << "\n\ttotal plotting errors:\t" << plotErr << "| -> adjust y-range in ComputeCorrelationFactors!" << endl;
        fLog << "\t|------------" << endl;
    }

    fOutput->Write();
    fOutput->Close();
    fCorr.close();
    fLog.close();

    return;
}
