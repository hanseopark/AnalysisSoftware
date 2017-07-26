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
#include "TGraphAsymmErrors.h" 
#include "TGaxis.h"
#include "TMarker.h"
#include "CommonHeaders/PlottingGammaConversionHistos.h"
#include "CommonHeaders/PlottingGammaConversionAdditional.h"
#include "CommonHeaders/FittingGammaConversion.h"
#include "CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "CommonHeaders/ConversionFunctions.h"

void FinaliseSystematicErrorsConvPbPb_LHC11h(   const char* nameDataFileErrors  = "", 
                                                TString energy                  = "", 
                                                TString meson                   = "", 
                                                Int_t numberOfPtBins            = 1 ,
                                                Int_t numberCutStudies          = 1, 
                                                Int_t offSetBeginning           = 0, 
                                                TString centPercent             = "0-5%",
                                                TString cent                    = "0005",
                                                TString suffix                  = "pdf", 
                                                Int_t smooth                    = 2, 
                                                Int_t mode                      = 0
                                            ){
	
	StyleSettingsThesis();	
	SetPlotStyle();
	
	TString date = ReturnDateString();
	TString dateForOutput = ReturnDateStringForOutput();
	TString collisionSystem= ReturnFullCollisionsSystem(energy);
	TString centralityRead = cent;

    Int_t 	numberOfEntriesPos = 0;
	Int_t 	numberOfEntriesNeg = 0;
	
	TFile* fileErrorInput= new TFile(nameDataFileErrors);
	const Int_t nPtBins = numberOfPtBins;
	const Int_t nCuts = numberCutStudies;
	Double_t* ptBins = NULL;
	Double_t* ptBinsErr = NULL;
	Double_t* newPoint = NULL;
	
	TString nameCutVariation[11] = {"Yield extract.", "dE/dx e-line", "dE/dx #pi-line", 
									"TPC cluster", "Single e^{#pm} #it{p}_{T}", "2D #chi^{2} #gamma, #psi_{pair} #gamma",
									"2D q_{T}",	"#alpha meson", "#varphi_{conv}",
									"Pile-up", "Yield extract. #pi^{0}"
// 									"#eta",	"Max #pi momentum", "TOF", "cosine point. angle"
                                    };  
	if (meson.CompareTo("EtaToPi0") == 0) nameCutVariation[0]          = "Yield extract. #eta";
	TString nameCutVariationSC[11] = {"YieldExtraction", "dEdxE", "dEdxPi", 
										"TPCCluster", "SinglePt", "Chi2", 
										"Qt", "Alpha", "ConvPhi",
										"Pileup", "YieldExtraction"
// 										"Eta", "pdEdxPi", "TOF", "CosPoint"
                                    };
    // Set colors and markers
    Color_t color[20];
    Color_t markerStyle[20];	
    for (Int_t k =0; k<nCuts; k++ ){
        color[k]                                = GetColorSystematics( nameCutVariationSC[k] ); 
        if(k==10) color[k] = GetColorSystematics("Rapidity"); 
        markerStyle[k]                          = GetMarkerStyleSystematics( nameCutVariationSC[k] );     
//         nameCutVariation[k]                     = GetSystematicsName(nameCutVariationSC[k]);
    }
		
    Bool_t benable[11]              = { 0, 0, 0, 0, 0, 
                                        0, 0, 0, 0, 0,
                                        0};
    Bool_t benableMeson[11]         = { 1, 1, 1, 1, 1, 
                                        1, 1, 1, 1, 0,
                                        0};
    Bool_t benableRatio[11]         = { 1, 1, 1, 1, 1, 
                                        1, 1, 1, 1, 0,
                                        1};
    for (Int_t i = 0; i < numberCutStudies; i++){
        if(meson.CompareTo("EtaToPi0")==0){
            benable[i] = benableRatio[i];    
        } else {
            benable[i] = benableMeson[i];
        }
    }
    if(centPercent.CompareTo("20-40%") || centPercent.CompareTo("20-50%")) benable[9] = 0;
    
	TString outputdir = Form("SystematicErrorsCalculatedGammaPbPb_LHC11h_%s/%s",dateForOutput.Data(),meson.Data());
	gSystem->Exec("mkdir -p "+outputdir);
	TString outputsmooth = Form("%s/Smoothing",outputdir.Data());
	gSystem->Exec("mkdir -p "+outputsmooth);
	
	Double_t* errorsNeg[nCuts];
	Double_t errorsNegCorr[nCuts][nPtBins];
	Double_t errorsNegSummed[nPtBins];
	Double_t errorsNegCorrSummed[nPtBins];
	Double_t errorsNegCorrMatSummed[nPtBins];
	
	Double_t* errorsNegErr[nCuts];
	Double_t errorsNegErrCorr[nCuts][nPtBins];
	Double_t errorsNegErrSummed[nPtBins];
	Double_t errorsNegErrCorrSummed[nPtBins];
	
	Double_t* errorsPos[nCuts];
	Double_t errorsPosCorr[nCuts][nPtBins];
	Double_t errorsPosSummed[nPtBins];
	Double_t errorsPosCorrSummed[nPtBins];
	Double_t errorsPosCorrMatSummed[nPtBins];
	
	Double_t* errorsPosErr[nCuts];
	Double_t errorsPosErrSummed[nPtBins];
	Double_t errorsPosErrCorr[nCuts][nPtBins];
	Double_t errorsPosErrCorrSummed[nPtBins];
	
	Double_t errorsMean[nCuts][nPtBins];
	Double_t errorsMeanCorr[nCuts][nPtBins];
	Double_t errorsMeanSummed[nPtBins];
	Double_t errorsMeanCorrSummed[nPtBins];
	Double_t errorsMeanCorrMatSummed[nPtBins];

	Double_t* errorsMeanCorrSmoothed[nCuts];
	Double_t errorsMeanCorrSummedSmoothed[nPtBins];
	Double_t errorsMeanCorrMatSummedSmoothed[nPtBins];
	Double_t errorsMeanErrCorrSummedSmoothed[nPtBins];
	Double_t errorsMeanErrCorrMatSummedSmoothed[nPtBins];
	
	Double_t errorsMeanErr[nCuts][nPtBins];
	Double_t errorsMeanErrCorr[nCuts][nPtBins];
	Double_t errorsMeanErrSummed[nPtBins];
	Double_t errorsMeanErrCorrSummed[nPtBins];
	Double_t errorsMeanErrCorrMatSummed[nPtBins];

    Double_t errorsMeanNotSmooth[nCuts][nPtBins];
    Double_t errorsMeanErrNotSmooth[nCuts][nPtBins];
    Double_t errorsMeanCorrNotSmooth[nCuts][nPtBins];
    Double_t errorsMeanErrCorrNotSmooth[nCuts][nPtBins];
    
	TGraphErrors* negativeErrors[nCuts];
	TGraphErrors* negativeErrorsSummed;
	TGraphErrors* positiveErrors[nCuts];
	TGraphErrors* positiveErrorsSummed;
	TGraphErrors* negativeErrorsCorr[nCuts];
	TGraphErrors* negativeErrorsCorrSummed;
	TGraphErrors* positiveErrorsCorr[nCuts];
	TGraphErrors* positiveErrorsCorrSummed;
	TGraphErrors* meanErrors[nCuts];
	TGraphErrors* meanErrorsSummed;
	TGraphErrors* meanErrorsCorr[nCuts];
	TGraphErrors* meanErrorsCorrSummed;
	TGraphErrors* meanErrorsCorrSummedIncMat;
	TGraphErrors* meanErrorsCorrSummedSmoothed;
	TGraphErrors* meanErrorsCorrSummedIncMatSmoothed = NULL;
    TGraphErrors* graphMaterialError;
    TGraphErrors* graphMassResolError;
	TGraphErrors* meanErrorsCorrSmoothed[nCuts];
    TGraphErrors* meanErrorsCorrNotSmooth[nCuts];
	
	for (Int_t l = 0; l < nPtBins; l++){
        errorsPosSummed[l] = 0.;
        errorsNegSummed[l] = 0.;
        errorsMeanSummed[l] = 0.;
        errorsPosCorrSummed[l] = 0.;
        errorsNegCorrSummed[l] = 0.;
        errorsMeanCorrSummed[l] = 0.;
        errorsMeanCorrSummedSmoothed[l] = 0.;
    } 
	
	for (Int_t i = 0; i < nCuts; i++){
        TGraphAsymmErrors* graphPosErrors       = NULL;
        TGraphAsymmErrors* graphNegErrors       = NULL;
        TString nameGraphPos    = "";
        TString nameGraphNeg    = "";
        if(benable[i]){
            if (meson.CompareTo("EtaToPi0") == 0){
                if(i==0){// special for treatment eta to pi0 ratio syst errors
                    nameGraphPos        = Form("Eta_SystErrorRelPos_YieldExtraction_%s",centPercent.Data() );
                    nameGraphNeg        = Form("Eta_SystErrorRelNeg_YieldExtraction_%s",centPercent.Data() );         
                } else if(i==10){
                    nameGraphPos    = Form("Pi0EtaBinning_SystErrorRelPos_YieldExtraction_%s",centPercent.Data() );
                    nameGraphNeg    = Form("Pi0EtaBinning_SystErrorRelNeg_YieldExtraction_%s",centPercent.Data() );
                } else if(i==9){ //pile-up
                    nameGraphPos        = Form("Eta_SystErrorRel_BGEstimate_%s",centPercent.Data() );
                    nameGraphNeg        = Form("Eta_SystErrorRel_BGEstimate_%s",centPercent.Data() );   
                } else {
                    nameGraphPos = Form("%s_SystErrorRelPos_%s%s",meson.Data(),nameCutVariationSC[i].Data(),centPercent.Data()  );
                    nameGraphNeg = Form("%s_SystErrorRelNeg_%s%s",meson.Data(),nameCutVariationSC[i].Data(),centPercent.Data()  );
                }
            } else {
                if (i==0){
                    nameGraphPos        = Form("%s_SystErrorRelPos_YieldExtraction_%s",meson.Data(),centPercent.Data() );
                    nameGraphNeg        = Form("%s_SystErrorRelNeg_YieldExtraction_%s",meson.Data(),centPercent.Data() );
                } else if (i==9) {
                    nameGraphPos        = Form("%s_SystErrorRel_BGEstimate_%s",meson.Data(),centPercent.Data() );
                    nameGraphNeg        = Form("%s_SystErrorRel_BGEstimate_%s",meson.Data(),centPercent.Data() );
                } else {
                    nameGraphPos = Form("%s_SystErrorRelPos_%s%s",meson.Data(),nameCutVariationSC[i].Data(),centPercent.Data()  );
                    nameGraphNeg = Form("%s_SystErrorRelNeg_%s%s",meson.Data(),nameCutVariationSC[i].Data(),centPercent.Data()  );
                }
            }
            cout << "Cutstudies " << i<< "\t" <<nameGraphPos.Data() << "\t" << nameGraphNeg.Data()<<  endl;
            graphPosErrors = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphPos.Data());
            graphNegErrors = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphNeg.Data());
            
            for (Int_t j = 0; j < offSetBeginning; j++){
                graphPosErrors->RemovePoint(0);
                graphNegErrors->RemovePoint(0);
            }        
            if (i == 0) {
                ptBins = graphNegErrors->GetX();
                ptBinsErr = graphNegErrors->GetEXhigh();
            }
            errorsNeg[i] = graphNegErrors->GetY();
            errorsNegErr[i] = graphNegErrors->GetEYhigh();
            errorsPos[i] = graphPosErrors->GetY();
            errorsPosErr[i] = graphPosErrors->GetEYhigh();
                
            // Calculate systematic error from input spectrum
            CalculateMeanSysErr(errorsMean[i], errorsMeanErr[i], errorsPos[i], errorsNeg[i], nPtBins);	
            CorrectSystematicErrorsWithMean(errorsPos[i],errorsPosErr[i], errorsPosCorr[i], errorsPosErrCorr[i], nPtBins);
            CorrectSystematicErrorsWithMean(errorsNeg[i],errorsNegErr[i], errorsNegCorr[i], errorsNegErrCorr[i], nPtBins);
            CorrectSystematicErrorsWithMean(errorsMean[i], errorsMeanErr[i], errorsMeanCorr[i], errorsMeanErrCorr[i], nPtBins);

            // ***************************************************************************************************
            // ************************ Adjust errors if requested to fixed values *******************************
            // ***************************************************************************************************
            Double_t errorFixed                 = -1;
            Bool_t adjustPtDependent            = kFALSE;

            for (Int_t k = 0; k < nPtBins; k++){
                errorsMeanNotSmooth[i][k]                = 0;
                errorsMeanErrNotSmooth[i][k]             = 0.0;
                errorsMeanCorrNotSmooth[i][k]            = 0;
                errorsMeanErrCorrNotSmooth[i][k]         = 0.0;            
                
                errorsMeanNotSmooth[i][k]        = errorsMean[i][k];
                errorsMeanErrNotSmooth[i][k]     = errorsMeanErr[i][k];
                errorsMeanCorrNotSmooth[i][k]    = errorsMeanCorr[i][k];
                errorsMeanErrCorrNotSmooth[i][k] = errorsMeanErrCorr[i][k];
            }

            if(smooth==2){
                if(nameCutVariationSC[i].Contains("Alpha")){
                    for (Int_t k = 0; k < nPtBins; k++){
                        if(meson.CompareTo("Pi0")==0){
                            if(!centPercent.CompareTo("0-5%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 1.+pow(ptBins[k],2)*0.023;
                            } 
                            else if (!centPercent.CompareTo("0-10%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 0.98+pow(ptBins[k],2)*0.02;
                            }
                            else if(!centPercent.CompareTo("5-10%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 1.1+pow(ptBins[k],2)*0.023;
                            } 
                            else if(!centPercent.CompareTo("20-40%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 0.9+pow(ptBins[k],2)*0.023;
                            }
                            else if(!centPercent.CompareTo("20-50%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 0.88+pow(ptBins[k],2)*0.023;
                            }
                            
                            if (errorFixed != -1){
                                errorsMean[i][k]        = errorFixed;
                                errorsMeanErr[i][k]     = errorFixed*0.01;
                                errorsMeanCorr[i][k]    = errorFixed;
                                errorsMeanErrCorr[i][k] = errorFixed*0.01;
                            }
                        } 
                        else {
                            if(!centPercent.CompareTo("0-5%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 1+2.*(ptBins[k]-2)+20/pow(ptBins[k]+0.2,2);
                            } 
                            else if(!centPercent.CompareTo("5-10%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 1+2.*(ptBins[k]-2)+20/pow(ptBins[k]+0.2,2);
                            } 
                            else if (!centPercent.CompareTo("0-10%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 1+2.*(ptBins[k]-2)+20/pow(ptBins[k]+0.2,2);
                            }
                            else if(!centPercent.CompareTo("20-40%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 1+2.*(ptBins[k]-2)+20/pow(ptBins[k]+0.2,2);
                            }
                            else if(!centPercent.CompareTo("20-50%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 1+2.*(ptBins[k]-2)+20/pow(ptBins[k]+0.2,2);
                            }
                            
                            if (errorFixed != -1){
                                errorsMean[i][k]        = errorFixed;
                                errorsMeanErr[i][k]     = errorFixed*0.01;
                                errorsMeanCorr[i][k]    = errorFixed;
                                errorsMeanErrCorr[i][k] = errorFixed*0.01;
                            }
                        }
                    }    
                
                }
                if(nameCutVariationSC[i].Contains("Chi2")){
                    if(meson.CompareTo("EtaToPi0")==0){
                        for (Int_t k = 0; k < nPtBins; k++){
                            
                            if(!centPercent.CompareTo("0-5%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 9+1.3*(ptBins[k]-6)+15/pow(ptBins[k],1.5);
                            } 
                            else if(!centPercent.CompareTo("5-10%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 8+1.3*(ptBins[k]-6)+15/pow(ptBins[k],1.5);
                            } 
                            else if (!centPercent.CompareTo("0-10%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 9+1.3*(ptBins[k]-6)+20/pow(ptBins[k],1.5);
                            }
                            else if(!centPercent.CompareTo("20-40%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 8+1.3*(ptBins[k]-6)+10/pow(ptBins[k],1.5);
                            }
                            else if(!centPercent.CompareTo("20-50%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 8+1.3*(ptBins[k]-6)+10/pow(ptBins[k],1.5);
                            }
                                                    
                            if (errorFixed != -1){
                                errorsMean[i][k]        = errorFixed;
                                errorsMeanErr[i][k]     = errorFixed*0.01;
                                errorsMeanCorr[i][k]    = errorFixed;
                                errorsMeanErrCorr[i][k] = errorFixed*0.01;
                            }
                        }

                    }
                    else if(meson.CompareTo("Pi0")==0){
                        for (Int_t k = 0; k < nPtBins; k++){
                            if(!centPercent.CompareTo("0-5%")){
                                adjustPtDependent               = kTRUE;
                                if(ptBins[k] <= 7) 
                                    errorFixed = -0.6+0.3*ptBins[k]+20/(pow(ptBins[k]+1.5,1.5));
                                else errorFixed = 0.96+pow(ptBins[k],2)*0.03;
                            } 
                            else if(!centPercent.CompareTo("5-10%")){
                                adjustPtDependent               = kTRUE;
                                if(ptBins[k] <= 6) 
                                    errorFixed = -0.6+0.3*ptBins[k]+20/(pow(ptBins[k]+1.5,1.5));
                                else errorFixed = 0.98+pow(ptBins[k],2)*0.03;
                            } 
                            else if (!centPercent.CompareTo("0-10%")){
                                adjustPtDependent               = kTRUE;
                                if(ptBins[k] <= 6) 
                                    errorFixed = -0.6+0.21*ptBins[k]+20/(pow(ptBins[k]+1.5,1.5));
                                else errorFixed = 0.98+pow(ptBins[k],2)*0.02;
                            }
                            else if(!centPercent.CompareTo("20-40%")){
                                adjustPtDependent               = kTRUE;
                                if(ptBins[k] <= 6) 
                                    errorFixed = -0.6+0.23*ptBins[k]+16/(pow(ptBins[k]+1.8,1.5));
                                else errorFixed = 0.98+pow(ptBins[k],2)*0.02;
                            }
                            else if(!centPercent.CompareTo("20-50%")){
                                adjustPtDependent               = kTRUE;
                                if(ptBins[k] <= 6) 
                                    errorFixed = -0.6+0.23*ptBins[k]+16/(pow(ptBins[k]+1.8,1.5));
                                else errorFixed = 0.98+pow(ptBins[k],2)*0.02;
                            }
                            
                            if (errorFixed != -1){
                                errorsMean[i][k]        = errorFixed;
                                errorsMeanErr[i][k]     = errorFixed*0.01;
                                errorsMeanCorr[i][k]    = errorFixed;
                                errorsMeanErrCorr[i][k] = errorFixed*0.01;
                            }
                        }
                    } 
                    else {
                        for (Int_t k = 0; k < nPtBins; k++){
                            
                            if(!centPercent.CompareTo("0-5%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 10+1.3*(ptBins[k]-6)+20/pow(ptBins[k],1.5);
                            } 
                            else if(!centPercent.CompareTo("5-10%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 8+1.3*(ptBins[k]-6)+20/pow(ptBins[k],1.5);
                            } 
                            else if (!centPercent.CompareTo("0-10%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 9+1.3*(ptBins[k]-6)+20/pow(ptBins[k],1.5);
                            }
                            else if(!centPercent.CompareTo("20-40%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 10+1.3*(ptBins[k]-6)+10/pow(ptBins[k],1.5);
                            }
                            else if(!centPercent.CompareTo("20-50%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 8+1.3*(ptBins[k]-6)+20/pow(ptBins[k],1.5);
                            }
                            
                            if (errorFixed != -1){
                                errorsMean[i][k]        = errorFixed;
                                errorsMeanErr[i][k]     = errorFixed*0.01;
                                errorsMeanCorr[i][k]    = errorFixed;
                                errorsMeanErrCorr[i][k] = errorFixed*0.01;
                            }
                        }
                    }
                } 
                if(nameCutVariationSC[i].Contains("ConvPhi")){
                    if(meson.CompareTo("Pi0")==0){
                        for (Int_t k = 0; k < nPtBins; k++){

                            if(!centPercent.CompareTo("0-5%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 1.1+pow(ptBins[k],1.5)*0.015;
                            } 
                            else if(!centPercent.CompareTo("5-10%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 1.1+pow(ptBins[k],1.8)*0.023;
                            } 
                            else if (!centPercent.CompareTo("0-10%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 0.98+pow(ptBins[k],1.5)*0.02;
                            }
                            else if(!centPercent.CompareTo("20-40%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 0.9+pow(ptBins[k],1.5)*0.02;
                            }
                            else if(!centPercent.CompareTo("20-50%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 0.88+pow(ptBins[k],1.5)*0.018;
                            }
                            
                            if (errorFixed != -1){
                                errorsMean[i][k]        = errorFixed;
                                errorsMeanErr[i][k]     = errorFixed*0.01;
                                errorsMeanCorr[i][k]    = errorFixed;
                                errorsMeanErrCorr[i][k] = errorFixed*0.01;
                            }
                        }
                    } else {
                        for (Int_t k = 0; k < nPtBins; k++){
                            
                            if(!centPercent.CompareTo("0-5%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 4.+pow(ptBins[k],1.5)*0.015;
                            } 
                            else if(!centPercent.CompareTo("5-10%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 4.+pow(ptBins[k],1.8)*0.023;
                            } 
                            else if (!centPercent.CompareTo("0-10%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 4.+pow(ptBins[k],1.5)*0.02;
                            }
                            else if(!centPercent.CompareTo("20-40%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 5.+pow(ptBins[k],1.5)*0.02;
                            }
                            else if(!centPercent.CompareTo("20-50%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 4+pow(ptBins[k],1.5)*0.018;
                            }
                            
                            if (errorFixed != -1){
                                errorsMean[i][k]        = errorFixed;
                                errorsMeanErr[i][k]     = errorFixed*0.01;
                                errorsMeanCorr[i][k]    = errorFixed;
                                errorsMeanErrCorr[i][k] = errorFixed*0.01;
                            }
                        }
                    }
                }
                if(nameCutVariationSC[i].Contains("TPCCluster")){
                    if(meson.CompareTo("EtaToPi0")==0){
                        for (Int_t k = 0; k < nPtBins; k++){
                            if(!centPercent.CompareTo("0-5%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 4.8+pow(ptBins[k],1.5)*0.023;
                            } 
                            else if(!centPercent.CompareTo("5-10%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 4.8+pow(ptBins[k],1.5)*0.023;
                            } 
                            else if (!centPercent.CompareTo("0-10%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 3+pow(ptBins[k],1.5)*0.02;
                            }
                            else if(!centPercent.CompareTo("20-40%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 2.5+pow(ptBins[k],1.5)*0.023;
                            }
                            else if(!centPercent.CompareTo("20-50%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 2.5+pow(ptBins[k],1.5)*0.023;
                            }
                            
                            if (errorFixed != -1){
                                errorsMean[i][k]        = errorFixed;
                                errorsMeanErr[i][k]     = errorFixed*0.01;
                                errorsMeanCorr[i][k]    = errorFixed;
                                errorsMeanErrCorr[i][k] = errorFixed*0.01;
                            }
                        }
                    } 
                    else if(meson.CompareTo("Pi0")==0){
                        for (Int_t k = 0; k < nPtBins; k++){
                            
                            if(!centPercent.CompareTo("0-5%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 1.+pow(ptBins[k],1.5)*0.023;
                            } 
                            else if (!centPercent.CompareTo("0-10%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 0.98+pow(ptBins[k],1.5)*0.02;
                            }
                            else if(!centPercent.CompareTo("5-10%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 1.1+pow(ptBins[k],1.5)*0.023;
                            } 
                            else if(!centPercent.CompareTo("20-40%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 0.9+pow(ptBins[k],1.5)*0.023;
                            }
                            else if(!centPercent.CompareTo("20-50%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 0.88+pow(ptBins[k],1.5)*0.023;
                            }
                            
                            if (errorFixed != -1){
                                errorsMean[i][k]        = errorFixed;
                                errorsMeanErr[i][k]     = errorFixed*0.01;
                                errorsMeanCorr[i][k]    = errorFixed;
                                errorsMeanErrCorr[i][k] = errorFixed*0.01;
                            }
                        }
                    } 
                    else {
                        for (Int_t k = 0; k < nPtBins; k++){
                            if(!centPercent.CompareTo("0-5%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 4.8+pow(ptBins[k],1.5)*0.023;
                            } 
                            else if(!centPercent.CompareTo("5-10%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 4.8+pow(ptBins[k],1.5)*0.023;
                            } 
                            else if (!centPercent.CompareTo("0-10%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 2.7+pow(ptBins[k],1.5)*0.02;
                            }
                            else if(!centPercent.CompareTo("20-40%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 2.+pow(ptBins[k],1.5)*0.023;
                            }
                            else if(!centPercent.CompareTo("20-50%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 2.5+pow(ptBins[k],1.5)*0.023;
                            }
                            
                            if (errorFixed != -1){
                                errorsMean[i][k]        = errorFixed;
                                errorsMeanErr[i][k]     = errorFixed*0.01;
                                errorsMeanCorr[i][k]    = errorFixed;
                                errorsMeanErrCorr[i][k] = errorFixed*0.01;
                            }
                        }
                    }
                }
                if(nameCutVariationSC[i].Contains("dEdxE")){
                    if(meson.CompareTo("EtaToPi0")==0){
                        for (Int_t k = 0; k < nPtBins; k++){
                            if(!centPercent.CompareTo("0-5%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed = 3.+0.3*ptBins[k]+20/(pow(ptBins[k]+0.3,2));
                            } 
                            else if(!centPercent.CompareTo("5-10%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed = 5+0.3*ptBins[k]+10/(pow(ptBins[k]+0.3,2));
                            } 
                            else if (!centPercent.CompareTo("0-10%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed = 3+0.3*ptBins[k]+20/(pow(ptBins[k]+0.3,2));
                            }
                            else if(!centPercent.CompareTo("20-40%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed = 2+0.3*ptBins[k]+20/(pow(ptBins[k]+0.3,2));
                            }
                            else if(!centPercent.CompareTo("20-50%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed = 2+0.3*ptBins[k]+20/(pow(ptBins[k]+0.3,2));
                            }
                            
                            if (errorFixed != -1){
                                errorsMean[i][k]        = errorFixed;
                                errorsMeanErr[i][k]     = errorFixed*0.01;
                                errorsMeanCorr[i][k]    = errorFixed;
                                errorsMeanErrCorr[i][k] = errorFixed*0.01;
                            }
                        }
                    } 
                    else if(meson.CompareTo("Pi0")==0){
                        for (Int_t k = 0; k < nPtBins; k++){
                            if(!centPercent.CompareTo("0-5%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed = -1.5+0.3*ptBins[k]+20/(pow(ptBins[k]+1.7,1.5));
                            } 
                            else if(!centPercent.CompareTo("5-10%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed = -1.5+0.3*ptBins[k]+20/(pow(ptBins[k]+1.7,1.5));
                            } 
                            else if (!centPercent.CompareTo("0-10%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed = -1.7+0.28*ptBins[k]+20/(pow(ptBins[k]+1.8,1.5));
                            }
                            else if(!centPercent.CompareTo("20-40%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed = -1.7+0.3*ptBins[k]+20/(pow(ptBins[k]+2.1,1.5));
                            }
                            else if(!centPercent.CompareTo("20-50%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed = -1.7+0.3*ptBins[k]+20/(pow(ptBins[k]+2.2,1.5));
                            }
                            
                            if (errorFixed != -1){
                                errorsMean[i][k]        = errorFixed;
                                errorsMeanErr[i][k]     = errorFixed*0.01;
                                errorsMeanCorr[i][k]    = errorFixed;
                                errorsMeanErrCorr[i][k] = errorFixed*0.01;
                            }
                        }
                    } 
                    else {
                        for (Int_t k = 0; k < nPtBins; k++){
                            if(!centPercent.CompareTo("0-5%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed = 4+0.3*ptBins[k]+20/(pow(ptBins[k]+0.3,2));
                            } 
                            else if(!centPercent.CompareTo("5-10%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed = 4+0.3*ptBins[k]+20/(pow(ptBins[k]+0.3,2));
                            } 
                            else if (!centPercent.CompareTo("0-10%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed = 3.7+0.3*ptBins[k]+20/(pow(ptBins[k]+0.3,2));
                            }
                            else if(!centPercent.CompareTo("20-40%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed = 2.8+0.3*ptBins[k]+20/(pow(ptBins[k]+0.3,2));
                            }
                            else if(!centPercent.CompareTo("20-50%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed = 4+0.3*ptBins[k]+20/(pow(ptBins[k]+0.3,2));
                            }

                            if (errorFixed != -1){
                                errorsMean[i][k]        = errorFixed;
                                errorsMeanErr[i][k]     = errorFixed*0.01;
                                errorsMeanCorr[i][k]    = errorFixed;
                                errorsMeanErrCorr[i][k] = errorFixed*0.01;
                            }
                        }
                    }
                }
                if(nameCutVariationSC[i].Contains("dEdxPi")){
                    if(meson.CompareTo("EtaToPi0")==0){
                        for (Int_t k = 0; k < nPtBins; k++){
                            if(!centPercent.CompareTo("0-5%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 7+(6+ptBins[k])*0.1+1./pow(ptBins[k]-0.7,2);
                            } 
                            else if(!centPercent.CompareTo("5-10%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 7+(6+ptBins[k])*0.1+1./pow(ptBins[k]-0.7,2);
                            } 
                            else if (!centPercent.CompareTo("0-10%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 5+(6+ptBins[k])*0.1+1./pow(ptBins[k]-0.7,2);
                            }
                            else if(!centPercent.CompareTo("20-40%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 5+(6+ptBins[k])*0.1+1./pow(ptBins[k]-0.7,2);
                            }
                            else if(!centPercent.CompareTo("20-50%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 0.3+(2+ptBins[k])*0.5+1./pow(ptBins[k]-0.7,2);
                            }

                            if (errorFixed != -1){
                                errorsMean[i][k]        = errorFixed;
                                errorsMeanErr[i][k]     = errorFixed*0.01;
                                errorsMeanCorr[i][k]    = errorFixed;
                                errorsMeanErrCorr[i][k] = errorFixed*0.01;
                            }
                        }
                    } 
                    else if(meson.CompareTo("Pi0")==0){
                        for (Int_t k = 0; k < nPtBins; k++){
                            if(!centPercent.CompareTo("0-5%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 1.3+pow(ptBins[k],2)*0.023;
                            } 
                            else if(!centPercent.CompareTo("5-10%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 1.1+pow(ptBins[k],2)*0.023;
                            } 
                            else if (!centPercent.CompareTo("0-10%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 0.98+pow(ptBins[k],2)*0.017;
                            }
                            else if(!centPercent.CompareTo("20-40%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 0.95+pow(ptBins[k],2)*0.03;
                            }
                            else if(!centPercent.CompareTo("20-50%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 0.88+pow(ptBins[k],2)*0.03;
                            }
                            
                            if (errorFixed != -1){
                                errorsMean[i][k]        = errorFixed;
                                errorsMeanErr[i][k]     = errorFixed*0.01;
                                errorsMeanCorr[i][k]    = errorFixed;
                                errorsMeanErrCorr[i][k] = errorFixed*0.01;
                            }
                        }
                    } 
                    else {
                        for (Int_t k = 0; k < nPtBins; k++){
                            if(!centPercent.CompareTo("0-5%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 7+(6+ptBins[k])*0.1+1./pow(ptBins[k]-0.7,2);
                            } 
                            else if(!centPercent.CompareTo("5-10%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 7+(6+ptBins[k])*0.1+1./pow(ptBins[k]-0.7,2);
                            } 
                            else if (!centPercent.CompareTo("0-10%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 7+(6+ptBins[k])*0.1+1./pow(ptBins[k]-0.7,2);
                            }
                            else if(!centPercent.CompareTo("20-40%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 6+(6+ptBins[k])*0.1+1./pow(ptBins[k]-0.7,2);
                            }
                            else if(!centPercent.CompareTo("20-50%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 6+(6+ptBins[k])*0.1+1./pow(ptBins[k]-0.7,2);
                            }
                            
                            if (errorFixed != -1){
                                errorsMean[i][k]        = errorFixed;
                                errorsMeanErr[i][k]     = errorFixed*0.01;
                                errorsMeanCorr[i][k]    = errorFixed;
                                errorsMeanErrCorr[i][k] = errorFixed*0.01;
                            }
                        }
                    }
                }
                if(nameCutVariationSC[i].Contains("Qt")){
                    if(meson.CompareTo("EtaToPi0")==0){
                        for (Int_t k = 0; k < nPtBins; k++){
                            if(!centPercent.CompareTo("0-5%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 7+pow(ptBins[k],2)*0.02;
                            } 
                            else if(!centPercent.CompareTo("5-10%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 5+pow(ptBins[k],2)*0.07;
                            } 
                            else if (!centPercent.CompareTo("0-10%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 3+pow(ptBins[k],2)*0.07;
                            }
                            else if(!centPercent.CompareTo("20-40%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 4+pow(ptBins[k],2)*0.07;
                            }
                            else if(!centPercent.CompareTo("20-50%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 2+pow(ptBins[k],2)*0.07;
                            }
                            
                            if (errorFixed != -1){
                                errorsMean[i][k]        = errorFixed;
                                errorsMeanErr[i][k]     = errorFixed*0.01;
                                errorsMeanCorr[i][k]    = errorFixed;
                                errorsMeanErrCorr[i][k] = errorFixed*0.01;
                            }
                        }
                    } 
                    else if(meson.CompareTo("Pi0")==0){
                        for (Int_t k = 0; k < nPtBins; k++){
                            
                            if(!centPercent.CompareTo("0-5%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 1.45+pow(ptBins[k],2)*0.02;
                            } 
                            else if(!centPercent.CompareTo("5-10%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 1.1+pow(ptBins[k],2)*0.023;
                            } 
                            else if (!centPercent.CompareTo("0-10%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 0.98+pow(ptBins[k],2)*0.02;
                            }
                            else if(!centPercent.CompareTo("20-40%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 1.1+pow(ptBins[k],2)*0.02;
                            }
                            else if(!centPercent.CompareTo("20-50%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 1.+pow(ptBins[k],2)*0.02;
                            }
                            
                            if (errorFixed != -1){
                                errorsMean[i][k]        = errorFixed;
                                errorsMeanErr[i][k]     = errorFixed*0.01;
                                errorsMeanCorr[i][k]    = errorFixed;
                                errorsMeanErrCorr[i][k] = errorFixed*0.01;
                            }
                        }
                    } 
                    else {
                        for (Int_t k = 0; k < nPtBins; k++){
                            if(!centPercent.CompareTo("0-5%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 5+pow(ptBins[k],2)*0.1;
                            } 
                            else if(!centPercent.CompareTo("5-10%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 5+pow(ptBins[k],2)*0.1;
                            } 
                            else if (!centPercent.CompareTo("0-10%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 3+pow(ptBins[k],2)*0.08;
                            }
                            else if(!centPercent.CompareTo("20-40%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 4+pow(ptBins[k],2)*0.1;
                            }
                            else if(!centPercent.CompareTo("20-50%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 2+pow(ptBins[k],2)*0.1;
                            }
                            
                            if (errorFixed != -1){
                                errorsMean[i][k]        = errorFixed;
                                errorsMeanErr[i][k]     = errorFixed*0.01;
                                errorsMeanCorr[i][k]    = errorFixed;
                                errorsMeanErrCorr[i][k] = errorFixed*0.01;
                            }
                        }
                    }
                }
                if(nameCutVariationSC[i].Contains("SinglePt")){
                    if(meson.CompareTo("EtaToPi0")==0){
                        for (Int_t k = 0; k < nPtBins; k++){

                            if(!centPercent.CompareTo("0-5%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed = 4.8+0.15*ptBins[k]+18/(pow(ptBins[k]+1.5,1.5));
                            } 
                            else if(!centPercent.CompareTo("5-10%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed = 2+0.15*ptBins[k]+20/(pow(ptBins[k]+1.5,1.5));
                            } 
                            else if (!centPercent.CompareTo("0-10%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed = 2+0.15*ptBins[k]+17/(pow(ptBins[k]+1.5,1.5));
                            }
                            else if(!centPercent.CompareTo("20-40%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed = 1.+0.15*ptBins[k]+20/(pow(ptBins[k]+1.5,1.5));
                            }
                            else if(!centPercent.CompareTo("20-50%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed = 1.7+0.15*ptBins[k]+17/(pow(ptBins[k]+1.5,1.5));
                            }
                            
                            if (errorFixed != -1){
                                errorsMean[i][k]        = errorFixed;
                                errorsMeanErr[i][k]     = errorFixed*0.01;
                                errorsMeanCorr[i][k]    = errorFixed;
                                errorsMeanErrCorr[i][k] = errorFixed*0.01;
                            }
                        }
                    } 
                    else if(meson.CompareTo("Pi0")==0){
                        for (Int_t k = 0; k < nPtBins; k++){
                            
                            if(!centPercent.CompareTo("0-5%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed = -1+0.15*ptBins[k]+18/(pow(ptBins[k]+1.5,1.5));
                            } 
                            else if(!centPercent.CompareTo("5-10%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed = -1+0.15*ptBins[k]+20/(pow(ptBins[k]+1.5,1.5));
                            } 
                            else if (!centPercent.CompareTo("0-10%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed = -1+0.15*ptBins[k]+17/(pow(ptBins[k]+1.5,1.5));
                            }
                            else if(!centPercent.CompareTo("20-40%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed = -1.1+0.15*ptBins[k]+20/(pow(ptBins[k]+1.5,1.5));
                            }
                            else if(!centPercent.CompareTo("20-50%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed = -1.+0.15*ptBins[k]+17/(pow(ptBins[k]+1.5,1.5));
                            }
                            
                            if (errorFixed != -1){
                                errorsMean[i][k]        = errorFixed;
                                errorsMeanErr[i][k]     = errorFixed*0.01;
                                errorsMeanCorr[i][k]    = errorFixed;
                                errorsMeanErrCorr[i][k] = errorFixed*0.01;
                            }
                        }
                    } 
                    else {
                        for (Int_t k = 0; k < nPtBins; k++){
                            
                            if(!centPercent.CompareTo("0-5%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed = 4.8+0.15*ptBins[k]+18/(pow(ptBins[k]+1.5,1.5));
                            } 
                            else if(!centPercent.CompareTo("5-10%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed = 2+0.15*ptBins[k]+20/(pow(ptBins[k]+1.5,1.5));
                            } 
                            else if (!centPercent.CompareTo("0-10%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed = 2+0.15*ptBins[k]+17/(pow(ptBins[k]+1.5,1.5));
                            }
                            else if(!centPercent.CompareTo("20-40%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed = 1.+0.15*ptBins[k]+20/(pow(ptBins[k]+1.5,1.5));
                            }
                            else if(!centPercent.CompareTo("20-50%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed = 1.5+0.15*ptBins[k]+17/(pow(ptBins[k]+1.5,1.5));
                            }
                            
                            if (errorFixed != -1){
                                errorsMean[i][k]        = errorFixed;
                                errorsMeanErr[i][k]     = errorFixed*0.01;
                                errorsMeanCorr[i][k]    = errorFixed;
                                errorsMeanErrCorr[i][k] = errorFixed*0.01;
                            }
                        }
                    }
                }            
                // put fixed values for pt independent errors, which were adjusted
                if (!adjustPtDependent && errorFixed != -1){
                    for (Int_t k = 0; k < nPtBins; k++){
                        errorsMean[i][k]        = errorFixed;
                        errorsMeanErr[i][k]     = errorFixed*0.01;
                        errorsMeanCorr[i][k]    = errorFixed;
                        errorsMeanErrCorr[i][k] = errorFixed*0.01;
                    }
                }    

                
            } // end smoothing if
            
            
            cout << "errors added quadratically for variation " << nameCutVariationSC[i].Data() << endl;
            for (Int_t l = 0; l < nPtBins; l++){
                errorsPosSummed[l]      = errorsPosSummed[l]+pow(errorsPos[i][l],2);
                errorsNegSummed[l]      = errorsNegSummed[l]+pow(errorsNeg[i][l],2);
                errorsMeanSummed[l]     = errorsMeanSummed[l]+pow(errorsMean[i][l],2);
                errorsPosCorrSummed[l]  = errorsPosCorrSummed[l]+pow(errorsPosCorr[i][l],2);
                errorsNegCorrSummed[l]  = errorsNegCorrSummed[l]+pow(errorsNegCorr[i][l],2);
                errorsMeanCorrSummed[l] = errorsMeanCorrSummed[l]+pow(errorsMeanCorr[i][l],2);
            }

            negativeErrors[i]             = new TGraphErrors(nPtBins,ptBins ,errorsNeg[i] ,ptBinsErr ,errorsNegErr[i] );
            meanErrors[i]                 = new TGraphErrors(nPtBins,ptBins ,errorsMean[i] ,ptBinsErr ,errorsMeanErr[i] );
            positiveErrors[i]             = new TGraphErrors(nPtBins,ptBins ,errorsPos[i] ,ptBinsErr ,errorsPosErr[i] );
            negativeErrorsCorr[i]         = new TGraphErrors(nPtBins,ptBins ,errorsNegCorr[i] ,ptBinsErr ,errorsNegErrCorr[i] );
            meanErrorsCorr[i]             = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorr[i] ,ptBinsErr ,errorsMeanErrCorr[i] );
            positiveErrorsCorr[i]         = new TGraphErrors(nPtBins,ptBins ,errorsPosCorr[i] ,ptBinsErr ,errorsPosErrCorr[i] );
            meanErrorsCorrNotSmooth[i]    = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrNotSmooth[i] ,ptBinsErr ,errorsMeanErrCorrNotSmooth[i] );
        }
	}
	
	//smoothing with the fit function:
	if(smooth==1){
		
		TCanvas* canvasCheckSmooth = new TCanvas("canvasCheckSmooth","",200,10,1350,900);  // gives the page size
		DrawGammaCanvasSettings( canvasCheckSmooth, 0.08, 0.01, 0.015, 0.09);
	
		TH2D *histo2DCheckSmooth;
		if (meson.Contains("Pi0")){
			histo2DCheckSmooth = new TH2D("histo2DCheckSmooth", "histo2DCheckSmooth", 20,0.,ptBins[nPtBins-1]+1,1000.,-0.5,30.);
		} else { 
			histo2DCheckSmooth = new TH2D("histo2DCheckSmooth", "histo2DCheckSmooth", 20,0.,ptBins[nPtBins-1]+1,1000.,-0.5,65.);
		}
		histo2DCheckSmooth->SetYTitle("mean smoothed systematic Err %");
		histo2DCheckSmooth->SetXTitle("#it{p}_{T} (GeV/#it{c})");
		histo2DCheckSmooth->GetYaxis()->SetNdivisions(510,kTRUE);
		histo2DCheckSmooth->GetYaxis()->SetDecimals();
		histo2DCheckSmooth->GetYaxis()->SetTitleSize(0.04);
		histo2DCheckSmooth->GetXaxis()->SetTitleSize(0.04);
		histo2DCheckSmooth->GetYaxis()->SetTitleOffset(0.9);
		histo2DCheckSmooth->GetXaxis()->SetTitleOffset(1.);
		histo2DCheckSmooth->GetYaxis()->SetLabelSize(0.03);
		histo2DCheckSmooth->GetXaxis()->SetLabelSize(0.03);
		histo2DCheckSmooth->SetTitle("");

		Double_t minPt = 1;
		Double_t maxPt = 14;
		TF1* pol2 = new TF1("pol2","[0]+[1]*x",0.4,14.);
		TF1* pol3 = new TF1("pol3","[0]+[1]*x*x+[2]*x*x*x",0.4,14.); 
		TF1* pol4 = new TF1("pol4","[0]+[1]*x+[2]*x*x+[3]*x*x*x*x",0.4,14.); 
		pol4->SetLineColor(kRed+2);
		pol3->SetLineColor(kBlue+1);
		pol2->SetLineColor(kGreen+2);

		Double_t newPoint;
		
		for(Int_t cut = 0; cut< numberCutStudies ; cut++){
		
			if(meson.Contains("Pi0")){
				minPt = 0.8;
				maxPt = ptBins[nPtBins-1]+1;
			} else if(meson.Contains("Eta")){
				minPt = 1.0;
				maxPt = ptBins[nPtBins-1]+1;
			}

			meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",minPt,maxPt);
			meanErrorsCorr[cut]->Fit(pol3,"QNRMEX0+","",minPt,maxPt);
			meanErrorsCorr[cut]->Fit(pol4,"QNRMEX0+","",minPt,maxPt);

			meanErrorsCorrSmoothed[cut] = (TGraphErrors*)meanErrorsCorr[cut]->Clone();

			for (Int_t l=0; l< nPtBins; l++){ //
				if (!(meson.CompareTo("EtaToPi0") == 0)){
//////////////////////////////////////////////////////////////////////////////
					if(nameCutVariationSC[cut].Contains("Alpha")){
						if(meson.Contains("Pi0")){
							if(cent.Contains("0010") || cent.Contains("0005") || cent.Contains("0510")){
								meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",0.8,10.);
                                newPoint = pol2->Eval(ptBins[l]);
                                meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
							} else if(cent.Contains("2040")){
								meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",0.8,12.);
								if(ptBins[l] < 10.){
                                    newPoint = pol2->Eval(ptBins[l]);
                                    meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                                } else if(ptBins[l] == 11.){
                                    cout << "taking graph point" << endl;
                                } else {
                                    newPoint = pol3->Eval(ptBins[l]);
                                    meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                                }
                            } else if(cent.Contains("2050")){
								meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",0.8,12.);
								if(ptBins[l] < 4.){
                                    newPoint = pol3->Eval(ptBins[l]);
                                    meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                                } else {
                                    newPoint = pol2->Eval(ptBins[l]);
                                    meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                                }                            
                            }
						} else if(meson.Contains("Eta")){
							if(cent.Contains("0010") ){
                                meanErrorsCorr[cut]->Fit(pol3,"QNRMEX0+","",1.,8.);
                                newPoint = pol3->Eval(ptBins[l]);
                                meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint*0.7);
                            } else if(cent.Contains("0005") ){
                                meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",1.,8.);
                                meanErrorsCorr[cut]->Fit(pol3,"QNRMEX0+","",1.,6.);
                                if(ptBins[l] < 6.){
                                    newPoint = pol3->Eval(ptBins[l]);
                                    meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint*0.8);
                                } else {
                                    newPoint = pol2->Eval(ptBins[l]);
                                    meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint*0.8);
                                }
                            } else if( cent.Contains("0510") ){
                                meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",1.5,4.);
                                meanErrorsCorr[cut]->Fit(pol3,"QNRMEX0+","",1.5,6.);
                                if(ptBins[l] < 8.){
                                    newPoint = pol3->Eval(ptBins[l]);
                                    meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                                } else {
                                    newPoint = pol2->Eval(ptBins[l]);
                                    meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint*2.);
                                }
							} else if(cent.Contains("2040")){
                                if(ptBins[l] < 8.){
                                    newPoint = pol3->Eval(ptBins[l]);
                                    meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint*0.9);
                                } else {
                                    newPoint = pol3->Eval(7.5);
                                    meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                                }
                            } else if(cent.Contains("2050")){
                                if(ptBins[l] < 6.){
                                    newPoint = pol3->Eval(ptBins[l]);
                                    meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                                } else if(ptBins[l] >= 6.){
                                    meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",1.,8.);
                                    newPoint = pol2->Eval(ptBins[l]);
                                    meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                                } else {
                                    cout << "taking graph point" << endl;
                                }
							}
						}
//////////////////////////////////////////////////////////////////////////////
					} else if(nameCutVariationSC[cut].Contains("Chi2")){
							if(meson.Contains("Pi0")){
								if(cent.Contains("0010")){
									meanErrorsCorr[cut]->Fit(pol4,"QNRMEX0+","",1.,maxPt);
 									meanErrorsCorr[cut]->Fit(pol3,"QNRMEX0+","",1.,6.);
                                    if(ptBins[l] < 5.){
                                      newPoint = pol3->Eval(ptBins[l]);
                                      meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                                    } else if(ptBins[l] > 5. && ptBins[l] < 10){
                                      newPoint = pol4->Eval(ptBins[l]);
                                      meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                                    } else {
                                      newPoint = pol2->Eval(ptBins[l]);
                                      meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                                    }
                                } else if(cent.Contains("0005") ){
                                    meanErrorsCorr[cut]->Fit(pol4,"QNRMEX0+","",1.,12.);
   									meanErrorsCorr[cut]->Fit(pol3,"QNRMEX0+","",1.,6.);
                                    if(ptBins[l] < 3.){
                                      newPoint = pol3->Eval(ptBins[l]);
                                      meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                                    } else if(ptBins[l] > 3. && ptBins[l] < 5){
                                      newPoint = pol4->Eval(ptBins[l]);
                                      meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                                    } else if(ptBins[l] > 12.){
                                      newPoint = pol2->Eval(ptBins[l]);
                                      meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                                    } else {
                                        cout << "taking graph point" << endl;
                                    }
                                } else if( cent.Contains("0510") ){
                                    meanErrorsCorr[cut]->Fit(pol4,"QNRMEX0+","",1.,12.);
									meanErrorsCorr[cut]->Fit(pol3,"QNRMEX0+","",1.,6.);
                                    if(ptBins[l] < 5.){
                                      newPoint = pol3->Eval(ptBins[l]);
                                      meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                                    } else if(ptBins[l] > 12.){
                                      newPoint = pol2->Eval(ptBins[l]);
                                      meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                                    } else {
                                        cout << "taking graph point" << endl;
                                    }
								} else if(cent.Contains("2040")){
									meanErrorsCorr[cut]->Fit(pol4,"QNRMEX0+","",1.,14.);
									meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",2.,12.);
									if(ptBins[l] < 10.){
										newPoint = pol4->Eval(ptBins[l]);
										meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                                    } else if(ptBins[l] == 11.){
                                        cout << "taking graph point" << endl;
									} else {
										newPoint = pol2->Eval(ptBins[l]);
										meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);										
									}
								} else if(cent.Contains("2050")){
									meanErrorsCorr[cut]->Fit(pol4,"QNRMEX0+","",1.,maxPt);
									meanErrorsCorr[cut]->Fit(pol3,"QNRMEX0+","",1.,8.);
                                    if(ptBins[l] < 6.){
                                        newPoint = pol4->Eval(ptBins[l]);
                                        meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                                    } else if(ptBins[l] == 9.){
                                        newPoint = pol3->Eval(ptBins[l]);
                                        meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);                                        
                                    } else {
                                        cout << "taking graph point" << endl;
                                    }
								}
							} else if(meson.Contains("Eta")){
								if(cent.Contains("0010")){
                                    if(ptBins[l] < 1.5){
                                        meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",1.,9.);
										newPoint = pol2->Eval(ptBins[l]);
										meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                                    } else if( ptBins[l] > 6.){
										newPoint = pol4->Eval(6.);
										meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                                    } else {
                                        cout << "taking graph point" << endl;
                                    }
                                } else if( cent.Contains("0005") ){
                                    meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",2,maxPt);
                                    newPoint = pol2->Eval(ptBins[l]);
                                    meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                                } else if(cent.Contains("0510") ){
                                    meanErrorsCorr[cut]->Fit(pol3,"QNRMEX0+","",1.5,maxPt);
                                    if(ptBins[l] < 8){    
                                        newPoint = pol3->Eval(ptBins[l]);
                                        meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                                    } else {
                                        cout << "taking graph point" << endl;
                                    }
                                } else if(cent.Contains("2040")){
                                    meanErrorsCorr[cut]->Fit(pol3,"QNRMEX0+","",1,maxPt);
                                    meanErrorsCorr[cut]->Fit(pol4,"QNRMEX0+","",1,maxPt);
                                    meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",2.,8.);
                                    if(ptBins[l] < 6){    
                                        newPoint = pol4->Eval(ptBins[l]);
                                        meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint*0.75);  
                                    } else {
                                        newPoint = pol2->Eval(ptBins[l]);
                                        meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint*0.8);
                                    }
                                } else if(cent.Contains("2050")){
                                    newPoint = pol3->Eval(ptBins[l]);
                                    meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
								}
							}
//////////////////////////////////////////////////////////////////////////////
					} else if( nameCutVariationSC[cut].Contains("ConvPhi")){
						if(meson.Contains("Pi0")){
							if(cent.Contains("0010") || cent.Contains("0005")){
                                meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",1.,10);
                                newPoint = pol2->Eval(ptBins[l]);
                                meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                            } else if(cent.Contains("0510")){
                                meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",0.8,8.);
                                newPoint = pol2->Eval(ptBins[l]);
                                meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);                                
							} else if(cent.Contains("2050") || cent.Contains("2040")){
                                meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",1.,maxPt);
                                newPoint = pol2->Eval(ptBins[l]);
                                meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
							}
						} else if(meson.Contains("Eta")){
							if(cent.Contains("0010") ){
                                if(ptBins[l] < 2 || ptBins[l] > 8.){
                                    meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",1.5,9.);
									newPoint = pol2->Eval(ptBins[l]);
									meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                                } else {
                                    cout << "taking graph point" << endl;
                                }
                            } else if(cent.Contains("0005") ){
                                meanErrorsCorr[cut]->Fit(pol4,"QNRMEX0+","",3,maxPt);
                                meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",3,maxPt);
                                if( ptBins[l] > 1.5 && ptBins[l] < 3.){
                                    newPoint = pol4->Eval(ptBins[l]);
                                    meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                                } else if( ptBins[l] > 6.){
                                    newPoint = pol4->Eval(5);
                                    meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                                } else {
                                    cout << "taking graph point" << endl;
                                }
                            } else if(cent.Contains("0510") ){
                                meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",minPt,maxPt);
                                newPoint = pol2->Eval(ptBins[l]);
                                meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                            } else if(cent.Contains("2040")){
                                meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",3,maxPt);
                                meanErrorsCorr[cut]->Fit(pol3,"QNRMEX0+","",minPt,maxPt);
                                if( ptBins[l] < 4.){
                                    newPoint = pol3->Eval(ptBins[l]);
                                    meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint*0.6);
                                } else {
                                    cout << "taking graph point" << endl;
                                }
							} else if(cent.Contains("2050")){
                                meanErrorsCorr[cut]->Fit(pol3,"QNRMEX0+","",2.,8.);
								if(ptBins[l] < 1.5 || ptBins[l] == 3.5){
									newPoint = pol3->Eval(ptBins[l]);
									meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
								} else if(ptBins[l] > 8){
                                    newPoint = pol3->Eval(6.);
                                    meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                                } else {
                                    cout << "taking graph point" << endl;
                                }
                            }
						}
//////////////////////////////////////////////////////////////////////////////
					} else if(nameCutVariationSC[cut].Contains("dEdxE")){
						if(meson.Contains("Pi0")){
							if(cent.Contains("0010")){
                                meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",2.,8.);
                                meanErrorsCorr[cut]->Fit(pol3,"QNRMEX0+","",1.,maxPt);
                                if(ptBins[l] < 2.5){
                                    newPoint = pol3->Eval(ptBins[l]);
                                    meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);                                
                                } else {
                                    cout << "taking graph point" << endl;
                                }
                            } else if(cent.Contains("0005") ){
                                meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",2.,8.);
                                meanErrorsCorr[cut]->Fit(pol3,"QNRMEX0+","",1.,maxPt);
                                if(ptBins[l] < 3){
                                    newPoint = pol3->Eval(ptBins[l]);
                                    meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);                                
                                } else {
                                    newPoint = pol2->Eval(ptBins[l]);
                                    meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                                }
                            } else if(cent.Contains("0510")){                                
                                meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",2.,8.);
                                meanErrorsCorr[cut]->Fit(pol3,"QNRMEX0+","",1.,maxPt);
                                if(ptBins[l] < 3.5){
                                    newPoint = pol3->Eval(ptBins[l]);
                                    meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);                                
                                } else if(ptBins[l] > 10){   
                                    newPoint = pol2->Eval(ptBins[l]);
                                    meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);                                
                                } else {
                                    cout << "taking graph point" << endl;
                                }
							} else if(cent.Contains("2050")){
								if(ptBins[l] < 5.){
                                    meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",0.8,8.);
									newPoint = pol2->Eval(ptBins[l]);
									meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
								} else if(ptBins[l] > 8.){
                                    meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",0.8,10.);
									newPoint = pol2->Eval(ptBins[l]);
									meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);	
                                } else {
                                    cout << "taking graph point" << endl;
                                }
							}  else if(cent.Contains("2040")){
                                if(ptBins[l] < 5){
                                    meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",0.8,8.);
                                    newPoint = pol2->Eval(ptBins[l]);
                                    meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                                } else if(ptBins[l] > 8.){
                                    meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",0.8,12.);
                                    newPoint = pol2->Eval(ptBins[l]);
                                    meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                                } else {
                                    cout << "taking graph point" << endl;
                                }
							}	
						} else if(meson.Contains("Eta")){
							if(cent.Contains("0010") ){
								meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",1.5,10.);
                                if(ptBins[l] < 2.){
									newPoint = pol2->Eval(ptBins[l]);
									meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
								} else {
									cout << "taking graph point" << endl;
								}					
                            } else if(cent.Contains("0005")  ){
                                newPoint = pol3->Eval(ptBins[l]);
                                meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint*0.8);
                            } else if(cent.Contains("0510") ){
                                meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",1.5,10.);
                                newPoint = pol2->Eval(ptBins[l]);
                                meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
							}  else if(cent.Contains("2040")){
                                if(ptBins[l] < 2.){
                                  meanErrorsCorr[cut]->Fit(pol3,"QNRMEX0+","",2.,6.);
                                  newPoint = pol3->Eval(ptBins[l]);
                                  meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                                } else if(ptBins[l] == 5.){
                                  newPoint = pol4->Eval(ptBins[l]);
                                  meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                                } else {
                                    cout << "taking graph point" << endl;
                                }
							}  else if(cent.Contains("2050")){
                                meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",1.5,8.);
                                newPoint = pol2->Eval(ptBins[l]);
                                meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                            }
						}
//////////////////////////////////////////////////////////////////////////////
					} else if(nameCutVariationSC[cut].Contains("dEdxPi")){
						if(meson.Contains("Pi0")){
							if(cent.Contains("0010") || cent.Contains("0005") || cent.Contains("0510") ){
                                meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",1.,8);
                                if(ptBins[l] < 1.){
//                                     meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",1.,3.);
                                    newPoint = pol2->Eval(ptBins[l]);
                                    meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                                } else if(ptBins[l] >= 1.){
//                                     meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",1.,8.);
									newPoint = pol2->Eval(ptBins[l]);
									meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                                } else {
                                    cout << "taking graph point" << endl;
                                }
							} else if(cent.Contains("2050")){
                                meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",1.,4.);
                                meanErrorsCorr[cut]->Fit(pol3,"QNRMEX0+","",minPt,10.);
                                meanErrorsCorr[cut]->Fit(pol4,"QNRMEX0+","",minPt,10.);
								if(ptBins[l] < 4.){
									newPoint = pol2->Eval(ptBins[l]);
									meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                                } else if(ptBins[l] == 7.){
                                    newPoint = pol3->Eval(ptBins[l]);
                                    meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                                } else if(ptBins[l] == 13){
                                    newPoint = pol3->Eval(11.);
                                    meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                                } else {
                                    cout << "taking graph point" << endl;
								}
                            } else if(cent.Contains("2040")){
                                meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",1.,4.);
                                newPoint = pol2->Eval(ptBins[l]);
                                meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
							}
						} else if(meson.Contains("Eta")){
							if(cent.Contains("0010") ){
                                if(ptBins[l] < 2.){
                                    newPoint = pol4->Eval(ptBins[l]);
                                    meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint*0.5);
                                } else if(ptBins[l] == 3.5){
                                    newPoint = pol4->Eval(ptBins[l]);
                                    meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);                                    
                                } else if(ptBins[l] == 7.){
                                    newPoint = pol3->Eval(ptBins[l]);
                                    meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);                                    
                                } else {
                                    cout << "taking graph point" << endl;
								}
                            } else if( cent.Contains("0005") ){
                                if(ptBins[l] == 1.75){
                                    newPoint = pol4->Eval(ptBins[l]);
                                    meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                                } else if(ptBins[l] == 5.){
                                    newPoint = pol4->Eval(4.);
                                    meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                                } else if(ptBins[l] > 6){
                                    newPoint = pol2->Eval(ptBins[l]);
                                    meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                                } else {
                                    cout << "taking graph point" << endl;
								}
                            } else if( cent.Contains("0510") ){
                                newPoint = pol3->Eval(ptBins[l]);
                                meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint*0.8);
                            }  else if( cent.Contains("2040")){
                                if(ptBins[l] < 3.){
                                    meanErrorsCorr[cut]->Fit(pol3,"QNRMEX0+","",2.,10.);
                                    newPoint = pol3->Eval(ptBins[l]);
                                    meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                                } else if(ptBins[l] >= 3.){
                                    meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",2,8.);
                                    newPoint = pol2->Eval(ptBins[l]);
                                    meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                                }
                            }  else if(cent.Contains("2050")){
//                                 if(ptBins[l] < 2.){
//                                   meanErrorsCorr[cut]->Fit(pol3,"QNRMEX0+","",1.5,8.);
//                                   newPoint = pol3->Eval(ptBins[l]);
//                                   meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
//                                 } else if(ptBins[l] >= 2.){
//                                   meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",1.5,8.);
//                                   newPoint = pol2->Eval(ptBins[l]);
//                                   meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
//                                 }
                                cout << "taking graph point" << endl;
                            }
						}
//////////////////////////////////////////////////////////////////////////////
					} else if(nameCutVariationSC[cut].Contains("Qt")){
						if(meson.Contains("Pi0")){
							if(cent.Contains("0010") || cent.Contains("0510") ){
                                if(ptBins[l] <= 6.){
                                    meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",1.,8.);
                                    newPoint = pol2->Eval(ptBins[l]);
                                    meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                                } else if(ptBins[l] > 6.){
                                    meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",1.,12.);
                                    newPoint = pol2->Eval(ptBins[l]);
                                    meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                                }
                            } else if(cent.Contains("0005")){
                                    meanErrorsCorr[cut]->Fit(pol3,"QNRMEX0+","",1.,maxPt);
                                    meanErrorsCorr[cut]->Fit(pol4,"QNRMEX0+","",1.,maxPt);
                                if(ptBins[l] <= 2.){
                                    meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",0.8,8.);
                                    newPoint = pol2->Eval(ptBins[l]);
                                    meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                                } else if(ptBins[l] > 2. && ptBins[l] < 6){
                                    newPoint = pol4->Eval(ptBins[l]);
                                    meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                                } else {
                                    meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",1.,14.);
                                    newPoint = pol2->Eval(ptBins[l]);
                                    meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                                }
							} else if(cent.Contains("2050")){
								meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",1.,12.);
								newPoint = pol2->Eval(ptBins[l]);
								meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
							} else if(cent.Contains("2040")){
								meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",1.,8.);
								newPoint = pol2->Eval(ptBins[l]);
								meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
							}
						} else if(meson.Contains("Eta")){
							if(cent.Contains("0010") ){
								if(ptBins[l] == 1.75 || ptBins[l] == 7){
                                    newPoint = pol4->Eval(1.75);
                                    meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                                } else {
                                    cout << "taking graph point" << endl;
                                }
                            } else if(cent.Contains("0005")){
                                meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",1.5,10.);
                                newPoint = pol2->Eval(ptBins[l]);
                                meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                            } else if(cent.Contains("0510") ){
								if(ptBins[l] == 1.75 || ptBins[l] > 6.){
                                    newPoint = pol4->Eval(5.);
                                    meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                                } else {
                                    cout << "taking graph point" << endl;
                                }
							} else if(cent.Contains("2040")){
								if(ptBins[l] < 6.){
									newPoint = pol3->Eval(ptBins[l]);
									meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                                } else if(ptBins[l] > 6.){
                                    meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",1.,8.);
                                    newPoint = pol2->Eval(4.);
                                    meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                                } else {
                                    cout << "taking graph point" << endl;
                                }
                            } else if(cent.Contains("2050")){
                                meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",4,8.);
                                if(ptBins[l] < 2|| ptBins[l] > 8.){
                                    newPoint = pol4->Eval(3.5);
                                    meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                                } else {
                                    cout << "taking graph point" << endl;
                                }
							}
						}
//////////////////////////////////////////////////////////////////////////////
					} else if(nameCutVariationSC[cut].Contains("SinglePt")){
						if(meson.Contains("Pi0")){
							if(cent.Contains("0010") ){
								meanErrorsCorr[cut]->Fit(pol4,"QNRMEX0+","",.8,maxPt);
                                meanErrorsCorr[cut]->Fit(pol3,"QNRMEX0+","",.8,maxPt);
                                meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",1.,maxPt);
                                if(ptBins[l] <= 1.5 ){
                                    newPoint = pol3->Eval(ptBins[l]);
                                    meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                                } else if(ptBins[l] > 12 ){
                                    newPoint = pol4->Eval(ptBins[l]);
                                    meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                                } else {
                                    cout << "taking graph point" << endl;
                                }
                            } else if(cent.Contains("0005")){
								meanErrorsCorr[cut]->Fit(pol4,"QNRMEX0+","",.8,maxPt);
                                meanErrorsCorr[cut]->Fit(pol3,"QNRMEX0+","",.8,maxPt);
                                meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",1.,maxPt);
                                if(ptBins[l] <= 1.5 ){
                                    newPoint = pol3->Eval(ptBins[l]);
                                    meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                                } else if(ptBins[l] > 12 ){
                                    newPoint = pol2->Eval(ptBins[l]);
                                    meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                                } else {
                                    cout << "taking graph point" << endl;
                                }
                            } else if(cent.Contains("0510")){
								meanErrorsCorr[cut]->Fit(pol4,"QNRMEX0+","",.8,maxPt);
                                meanErrorsCorr[cut]->Fit(pol3,"QNRMEX0+","",.8,maxPt);
                                meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",1.,maxPt);
                                if(ptBins[l] < 2 ){
                                    newPoint = pol3->Eval(ptBins[l]);
                                    meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                                } else if(ptBins[l] > 2 && ptBins[l] < 8){
                                    newPoint = pol4->Eval(ptBins[l]);
                                    meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);                                    
                                } else if(ptBins[l] > 8){
                                    newPoint = pol2->Eval(ptBins[l]);
                                    meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);                                    
                                }
							} else if(cent.Contains("2050")){
								meanErrorsCorr[cut]->Fit(pol4,"QNRMEX0+","",1.5,10);
                                if(ptBins[l] < 1.5 ){
                                    newPoint = pol4->Eval(ptBins[l]);
                                    meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                                } else if(ptBins[l] > 8. ){
                                    meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",1.5,12.);
                                    newPoint = pol2->Eval(ptBins[l]);
                                    meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                                } else {
                                    cout << "taking graph point" << endl;
                                }
							} else if(cent.Contains("2040")){
								meanErrorsCorr[cut]->Fit(pol4,"QNRMEX0+","",1,maxPt);
                                meanErrorsCorr[cut]->Fit(pol3,"QNRMEX0+","",1,12);
                                meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",2,12);
                                if(ptBins[l] <= 2.5 ){
									newPoint = pol4->Eval(ptBins[l]);
									meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);									
                                } else if(ptBins[l] > 6. &&  ptBins[l] < 8){
                                    newPoint = pol3->Eval(ptBins[l]);
                                    meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                                } else if(ptBins[l] > 10.){
                                    newPoint = pol2->Eval(ptBins[l]);
                                    meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                                } else {
                                    cout << "taking graph point" << endl;
                                }
                            }
						} else if(meson.Contains("Eta")){
							if(cent.Contains("0010") ){
								if(ptBins[l] < 4.){
                                    meanErrorsCorr[cut]->Fit(pol4,"QNRMEX0+","",1.5,maxPt);
									newPoint = pol4->Eval(ptBins[l]);
									meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                                } else {
                                    cout << "taking graph point" << endl;
                                }
                            } else if( cent.Contains("0005")){
                                newPoint = pol3->Eval(ptBins[l]);
                                meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                            } else if( cent.Contains("0510") ){
                                  newPoint = pol3->Eval(ptBins[l]);
                                  meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                            } else if(cent.Contains("2040")){
								if(ptBins[l] < 2.){
									newPoint = pol4->Eval(ptBins[l]);
									meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                                } else {
                                    cout << "taking graph point" << endl;
                                }
							} else if(cent.Contains("2050")){
                              if(ptBins[l] <= 8.){
									newPoint = pol3->Eval(ptBins[l]);
									meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                                } else if(ptBins[l] > 8.){
                                    meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",1.2,maxPt);
                                    newPoint = pol2->Eval(ptBins[l]);
                                    meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                                } else {
                                    cout << "taking graph point" << endl;
                                }
                            }
						}
//////////////////////////////////////////////////////////////////////////////
					} else if(nameCutVariationSC[cut].Contains("TPCCluster")){
						if(meson.Contains("Pi0")){
                            if(cent.Contains("0010")){
                                if(ptBins[l] < 2 || ptBins[l] >12 ){
                                    meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",2,5);
                                    newPoint = pol2->Eval(ptBins[l]);
                                    meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);                            
                                } else {
                                    cout << "taking graph point" << endl;
                                }
                            } else if(cent.Contains("0005")){
                                meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",1.5,12);  
                                if(ptBins[l] < 3 || ptBins[l] > 10 ){
                                    newPoint = pol2->Eval(ptBins[l]);
                                    meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);                                
                                } else {
                                    cout << "taking graph point" << endl;
                                }
                            } else if(cent.Contains("0510") ){
                                meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",1.,maxPt);
                                newPoint = pol2->Eval(ptBins[l]);
                                meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
							} else if(cent.Contains("2050")){
                                meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",1.,12.);
                                newPoint = pol2->Eval(ptBins[l]);
                                meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
							} else if(cent.Contains("2040")){
                                meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",1.,10.);
                                newPoint = pol2->Eval(ptBins[l]);
                                meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
							}
						} else if(meson.Contains("Eta")){
							if(cent.Contains("0010")){
                                meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",1.5,10);
                                newPoint = pol2->Eval(ptBins[l]);
                                meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                            } else if (cent.Contains("0510") ){
                                newPoint = pol2->Eval(ptBins[l]);
                                meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                            } else if(cent.Contains("0005") ){
                                if(ptBins[l] < 1.5 || ptBins[l] == 2.5 ){
                                    meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",3.,maxPt);
                                    newPoint = pol2->Eval(ptBins[l]);
                                    meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);                                    
                                } else {
                                    cout << "taking graph point" << endl;
                                }
							} else if(cent.Contains("2040")){
								meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",2.,8.);
								newPoint = pol2->Eval(ptBins[l]);
								meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);									
							} else if(cent.Contains("2050")){
                                if(ptBins[l] < 2 ||  ptBins[l] == 5. || ptBins[l] > 8 ){
                                    newPoint = pol3->Eval(5.);
                                    meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);                                    
                                } else {
                                    cout << "taking graph point" << endl;
                                }
							}
						}
					}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

				} else 	if (meson.CompareTo("EtaToPi0") == 0 && nCuts > 9){
				
					if(nameCutVariationSC[cut].Contains("Alpha")){
						if(cent.Contains("0010")){
                            newPoint = pol2->Eval(ptBins[l]);
                            meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                        } else if(cent.Contains("0005") ){
                            newPoint = pol3->Eval(4.);
                            meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                        } else if( cent.Contains("0510") ){
							meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",2.,6.);
                            newPoint = pol2->Eval(ptBins[l]);
                            meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                        } else if( cent.Contains("2040")){
                            if( ptBins[l] < 6.){
                                newPoint = pol3->Eval(ptBins[l]);
                                meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint*0.6);
                            } else {
                                newPoint = pol2->Eval(ptBins[l]);
                                meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint*0.5);
                            }
						} else if(cent.Contains("2050") ){
							if( ptBins[l] < 6.){
                                newPoint = pol3->Eval(ptBins[l]);
                                meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
							} else {
                                newPoint = pol3->Eval(6.);
                                meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
							}
						}
//////////////////////////////////////////////////////////////////////////////
					} else if(nameCutVariationSC[cut].Contains("Chi2")){
						if(cent.Contains("0010") ){
                            newPoint = pol3->Eval(5);
                            meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                        } else if(cent.Contains("0005") || cent.Contains("0510") ){
                            newPoint = pol4->Eval(5.);
                            meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
						} else if(cent.Contains("2050") || cent.Contains("2040")){
                            newPoint = pol3->Eval(ptBins[l]);
                            meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                        }
//////////////////////////////////////////////////////////////////////////////
					} else if( nameCutVariationSC[cut].Contains("ConvPhi")){
						if(cent.Contains("0010") ){
                            if( ptBins[l] < 2.){
                                newPoint = pol3->Eval(ptBins[l]);
                                meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint*0.5);
							} else {
								cout << "taking graph point" << endl;					
							}                            
                        } else if( cent.Contains("0005")  ){
                            meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",3,6.);
                            newPoint = pol2->Eval(ptBins[l]);
                            meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                        } else if(cent.Contains("0510") ){
                            newPoint = pol4->Eval(4.);
                            meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
						} else if(cent.Contains("2050")){
							meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",2.,6.);
							newPoint = pol2->Eval(ptBins[l]);
							meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                        } else if( cent.Contains("2040")){
                            newPoint = pol3->Eval(4.5);
                            meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                        }
//////////////////////////////////////////////////////////////////////////////
					} else if(nameCutVariationSC[cut].Contains("dEdxE")){
						if(cent.Contains("0010")){
							meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",2.5,maxPt);
                            newPoint = pol2->Eval(ptBins[l]);
                            meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                        } else if(cent.Contains("0510") ){
                            newPoint = pol2->Eval(ptBins[l]);
                            meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                        } else if(cent.Contains("0005") ){
                            newPoint = pol2->Eval(ptBins[l]);
                            meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                        }  else if(  cent.Contains("2040")){
                            meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",2.,maxPt);
                            newPoint = pol2->Eval(ptBins[l]);
                            meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
						}  else if(cent.Contains("2050")){
                            meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",1.5,8.);
                            newPoint = pol2->Eval(ptBins[l]);
                            meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                        }
//////////////////////////////////////////////////////////////////////////////
					} else if(nameCutVariationSC[cut].Contains("dEdxPi")){
						if(cent.Contains("0010")  ){
                            if(ptBins[l] < 8.){
                                newPoint = pol3->Eval(ptBins[l]);
                                meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint*0.8);
							} else {
                                newPoint = pol2->Eval(ptBins[l]);
                                meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint*0.8);
							}
                        } else if( cent.Contains("0005") ){
                            newPoint = pol4->Eval(4.);
                            meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                        } else if( cent.Contains("0510") ){
                            newPoint = pol3->Eval(2.);
                            meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
						}  else if(cent.Contains("2050") ){
                            if(ptBins[l] < 6.){
                                meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",1.5,8.);
                                newPoint = pol2->Eval(ptBins[l]);
                                meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                            } else if(ptBins[l] > 6.){
                                meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",2.,8.);
                                newPoint = pol2->Eval(ptBins[l]);
                                meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                            }
                        }  else if( cent.Contains("2040")){
                            if(ptBins[l] > 2. && ptBins[l] < 6.){
                                newPoint = pol4->Eval(ptBins[l]);
                                meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                            } else if(ptBins[l] == 7){
                                newPoint = pol3->Eval(ptBins[l]);
                                meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                            } else {
                                newPoint = pol2->Eval(6.);
                                meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);                                
                            }
                        }
//////////////////////////////////////////////////////////////////////////////
					} else if(nameCutVariationSC[cut].Contains("Qt")){
						if(cent.Contains("0010")  ){
                            if(ptBins[l] == 1.75 || ptBins[l] == 7.){
                                newPoint = pol4->Eval(2.);
                                meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);                                
							} else {
								cout << "taking graph point" << endl;					
							}
                        } else if( cent.Contains("0005") ){
                            newPoint = pol4->Eval(4.5);
                            meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                        } else if(cent.Contains("0510") ){
                            if(ptBins[l] == 1.75 || ptBins[l] > 6.){
                                newPoint = pol4->Eval(4.);
                                meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);                                
							} else {
								cout << "taking graph point" << endl;					
							}
                        } else if(cent.Contains("2040")){
                            meanErrorsCorr[cut]->Fit(pol3,"QNRMEX0+","",1.5,8.);
                            newPoint = pol3->Eval(3.5);
                            meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
						} else if(cent.Contains("2050")){
                            if(ptBins[l] < 2. || ptBins[l] > 4.){
								newPoint = pol3->Eval(3.);
								meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
							} else {
								cout << "taking graph point" << endl;					
							}
						}
//////////////////////////////////////////////////////////////////////////////
					} else if(nameCutVariationSC[cut].Contains("SinglePt")){
						if(cent.Contains("0010") ){
                            if(ptBins[l] < 3.){
                                newPoint = pol3->Eval(ptBins[l]);
                                meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
							} else {
								cout << "taking graph point" << endl;					
							}
                        } else if(cent.Contains("0005") ){
                            newPoint = pol3->Eval(ptBins[l]);
                            meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint*0.9);
                        } else if( cent.Contains("0510") ){
                            if(ptBins[l] < 8.){
                                newPoint = pol3->Eval(ptBins[l]);
                                meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                            } else {
                                newPoint = pol2->Eval(ptBins[l]);
                                meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);                                
                            }
                        } else if(cent.Contains("2040")){
                            if(ptBins[l] < 2.){
                                meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",2.,6.);
                                newPoint = pol2->Eval(ptBins[l]);
                                meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                            } else if(ptBins[l] > 8.){
                                newPoint = pol2->Eval(7.);
                                meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
							} else {
								cout << "taking graph point" << endl;					
							}
						} else if(cent.Contains("2050")){
                            if(ptBins[l] < 6.){
                                newPoint = pol3->Eval(ptBins[l]);
                                meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                            } else if(ptBins[l] > 8.){
                                newPoint = pol3->Eval(6.5);
                                meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
							} else {
								cout << "taking graph point" << endl;					
							}
						}
//////////////////////////////////////////////////////////////////////////////
					} else if(nameCutVariationSC[cut].Contains("TPCCluster")){
						if(cent.Contains("0010")){                            
                            if(ptBins[l] < 4.){
                                meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",4.,maxPt);
                                newPoint = pol2->Eval(ptBins[l]);
                                meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
							} else {
								cout << "taking graph point" << endl;					
							}
                        } else if(cent.Contains("0005")  ){
                            if(ptBins[l] < 3. || ptBins[l] > 8.){
                                meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",3.,8.);
                                newPoint = pol2->Eval(ptBins[l]);
                                meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
							} else {
								cout << "taking graph point" << endl;					
							}
                        } else if( cent.Contains("0510") ){
                            newPoint = pol2->Eval(ptBins[l]);
                            meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
                        } else if( cent.Contains("2040")){
                            if(ptBins[l] < 2. || ptBins[l] > 8.){
                                meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",2.,8.);
                                newPoint = pol2->Eval(ptBins[l]);
                                meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
							} else {
								cout << "taking graph point" << endl;					
							}
						} else if(cent.Contains("2050") ){
                            if(ptBins[l] < 1.5 || ptBins[l] == 5. || ptBins[l] > 8.){
                                newPoint = pol4->Eval(5.);
                                meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBins[l],newPoint);
							} else {
								cout << "taking graph point" << endl;					
							}
						}
					}
				}// else cout <<"bla"<<endl;
				
			} //end of pt bins loop	
			
			histo2DCheckSmooth->DrawCopy();

			TLegend* legendSingleSyst;
			if (meson.Contains("Pi0")){legendSingleSyst= new TLegend(0.13,0.7,0.65,0.95);} 
			else {legendSingleSyst= new TLegend(0.20,0.7,0.7,0.95);}
			legendSingleSyst->SetTextSize(0.035);
			legendSingleSyst->SetFillColor(0);
			legendSingleSyst->SetFillStyle(0);	
			legendSingleSyst->SetBorderSize(0);
			legendSingleSyst->SetMargin(0.1);
			legendSingleSyst->AddEntry(meanErrorsCorr[cut], Form("%s",nameCutVariationSC[cut].Data()),"p");
			legendSingleSyst->AddEntry(pol3,"pol3, [0]+[1]*x^{2}+[2]*x^{3}","l");
			legendSingleSyst->AddEntry(pol4,"pol4, [0]+[1]*x+[2]*x^{2}+[3]*x^{4}","l");
						
			DrawGammaSetMarkerTGraphErr(meanErrorsCorr[cut], 20+cut, 1.,color[cut],color[cut]);
			meanErrorsCorr[cut]->Draw("p,csame");
			pol2->Draw("same");
			pol3->Draw("same");
			pol4->Draw("same");

			DrawGammaSetMarkerTGraphErr(meanErrorsCorrSmoothed[cut], 20, 1.5,kGray+2,kGray+2);
			meanErrorsCorrSmoothed[cut]->Draw("p,csame");
			legendSingleSyst->AddEntry(meanErrorsCorrSmoothed[cut], Form("%s smoothed",nameCutVariationSC[cut].Data()),"p");
					
			legendSingleSyst->Draw("same");
			
			canvasCheckSmooth->SaveAs(Form("%s/SysMeanNewWithMeanSingle%s_%s_%s%s_%s.%s",outputdir.Data(),nameCutVariationSC[cut].Data(),meson.Data(), energy.Data(),cent.Data(),dateForOutput.Data(),suffix.Data()));
		} //cuts loop
	}//end smoothing if

	
	Double_t errorMaterial;
	Double_t errorMassResolution;
	if (meson.CompareTo("EtaToPi0")==0){
		errorMaterial = 0.;
		errorMassResolution = 0.;
	} else {
		errorMaterial = 4.50;
		errorMassResolution = 1.4; // was 1.2;
	}
	
	if(smooth==1){	
		for(Int_t cut = 0; cut< numberCutStudies; cut++){
			errorsMeanCorrSmoothed[cut] = meanErrorsCorrSmoothed[cut]->GetY();
			for (Int_t l = 0; l < nPtBins; l++){
				errorsMeanCorrSummedSmoothed[l] = errorsMeanCorrSummedSmoothed[l]+ pow(errorsMeanCorrSmoothed[cut][l],2);
	// 			cout << l << " errorsMeanCorrSummedSmoothed[l] " << errorsMeanCorrSummedSmoothed[l] << endl;
			}
		}
	}
	
	for (Int_t l = 0; l < nPtBins; l++){
		errorsPosSummed[l]                    = pow(errorsPosSummed[l],0.5);
		errorsMeanSummed[l]                   = pow(errorsMeanSummed[l],0.5);
        errorsPosErrSummed[l]                 = errorsPosSummed[l]*0.001;
		errorsMeanErrSummed[l]                = errorsMeanSummed[l]*0.001;
		errorsNegSummed[l]                    = -pow(errorsNegSummed[l],0.5);
		errorsNegErrSummed[l]                 = errorsNegSummed[l]*0.001;
		errorsPosCorrSummed[l]                = errorsPosCorrSummed[l]+ pow(2*errorMassResolution ,2.);
		errorsMeanCorrSummed[l]               = errorsMeanCorrSummed[l]+ pow(2*errorMassResolution ,2.);
		errorsNegCorrSummed[l]                = errorsNegCorrSummed[l]+ pow(2*errorMassResolution ,2.);
		errorsPosCorrMatSummed[l]             = pow(errorsPosCorrSummed[l]+ pow(2*errorMaterial ,2.),0.5);
		errorsMeanCorrMatSummed[l]            = pow(errorsMeanCorrSummed[l]+ pow(2*errorMaterial ,2.),0.5);
		errorsNegCorrMatSummed[l]             = -pow(errorsNegCorrSummed[l]+ pow(2*errorMaterial ,2.),0.5);

		if(smooth==1){
			errorsMeanCorrMatSummedSmoothed[l] = pow(errorsMeanCorrSummedSmoothed[l]+ pow(2*errorMaterial ,2.),0.5);
			errorsMeanCorrSummedSmoothed[l]    = pow(errorsMeanCorrSummedSmoothed[l],0.5);
			errorsMeanErrCorrSummedSmoothed[l] = errorsMeanCorrSummedSmoothed[l]*0.001;
		}

		errorsPosCorrSummed[l]                = pow(errorsPosCorrSummed[l],0.5);
		errorsMeanCorrSummed[l]               = pow(errorsMeanCorrSummed[l],0.5);
		errorsPosErrCorrSummed[l]             = errorsPosCorrSummed[l]*0.001;
		errorsMeanErrCorrSummed[l]            = errorsMeanCorrSummed[l]*0.001;
		errorsMeanErrCorrMatSummed[l]         = errorsMeanCorrMatSummed[l]*0.001;
		errorsMeanErrCorrMatSummedSmoothed[l] = errorsMeanCorrMatSummedSmoothed[l]*0.001;
		errorsNegCorrSummed[l]                = -pow(errorsNegCorrSummed[l],0.5);
		errorsNegErrCorrSummed[l]             = errorsNegCorrSummed[l]*0.001;
	}
	
	Double_t errorsMat[nPtBins];
	Double_t errorMassRes[nPtBins];
	for (Int_t l = 0; l < nPtBins; l++){
		errorsMat[l] = 2*errorMaterial;
		errorMassRes[l] = errorMassResolution;
		
	}
	graphMaterialError                     = new TGraphErrors(nPtBins,ptBins ,errorsMat ,ptBinsErr ,errorsMeanErrSummed );
	graphMassResolError                    = new TGraphErrors(nPtBins,ptBins ,errorMassRes ,ptBinsErr ,errorsMeanErrSummed );
	negativeErrorsSummed 		           = new TGraphErrors(nPtBins,ptBins ,errorsNegSummed ,ptBinsErr ,errorsNegErrSummed );
	negativeErrorsCorrSummed 	           = new TGraphErrors(nPtBins,ptBins ,errorsNegCorrSummed ,ptBinsErr ,errorsNegErrCorrSummed );
	positiveErrorsSummed                   = new TGraphErrors(nPtBins,ptBins ,errorsPosSummed ,ptBinsErr ,errorsPosErrSummed );
	positiveErrorsCorrSummed 	           = new TGraphErrors(nPtBins,ptBins ,errorsPosCorrSummed ,ptBinsErr ,errorsPosErrCorrSummed );
	meanErrorsSummed 			           = new TGraphErrors(nPtBins,ptBins ,errorsMeanSummed ,ptBinsErr ,errorsMeanErrSummed );
	meanErrorsCorrSummed		           = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrSummed ,ptBinsErr ,errorsMeanErrCorrSummed );
	meanErrorsCorrSummedIncMat 	           = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrMatSummed ,ptBinsErr ,errorsMeanErrCorrMatSummed );
	if(smooth==1){
		meanErrorsCorrSummedSmoothed 	   = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrSummedSmoothed ,ptBinsErr ,errorsMeanErrCorrSummed );
		meanErrorsCorrSummedIncMatSmoothed = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrMatSummedSmoothed,ptBinsErr ,errorsMeanErrCorrMatSummed );
	}
	
	
	
	/////////////////////////////////////////////////////////////////////////
	///////////////////////////////// Plotting //////////////////////////////
	/////////////////////////////////////////////////////////////////////////
	cout << "Plotting the systematic: " << endl;
    Double_t textSizeLabelsPixelmean        = 55;
    Double_t textSizeLabelsRelmean          = 55./1300;
	
	TCanvas* canvasSystErrors = new TCanvas("canvasSystErrors","",200,10,1300,1000);  // gives the page size
	DrawGammaCanvasSettings( canvasSystErrors, 0.075, 0.03, 0.015, 0.095);
        TH2D *histo2DSystErrors;
        if (meson.CompareTo("Pi0")==0)histo2DSystErrors = new TH2D("histo2DSystErrors", "histo2DSystErrors", 20,0.,15,1000.,0.,30.);
        else histo2DSystErrors = new TH2D("histo2DSystErrors", "histo2DSystErrors", 20,0.,11,1000.,0.,65.);
        SetStyleHistoTH2ForGraphs( histo2DSystErrors, "#it{p}_{T} (GeV/#it{c})", "mean systematic Err %",
                    0.85*textSizeLabelsRelmean, textSizeLabelsRelmean, 0.85*textSizeLabelsRelmean, textSizeLabelsRelmean,0.9,0.75);
        
        histo2DSystErrors->DrawCopy();
	
        TLegend* legendMean; 
        if (meson.Contains("Pi0")) legendMean = new TLegend(0.15,0.6,0.57,0.95);
        else legendMean= new TLegend(0.20,0.6,0.62,0.95);
        legendMean->SetTextSize(0.035);
        legendMean->SetFillColor(0);
        legendMean->SetBorderSize(0);
        for(Int_t i = 0; i< numberCutStudies ; i++){
            if(benable[i]){
                DrawGammaSetMarkerTGraphErr(meanErrors[i], markerStyle[i], 1.,color[i],color[i]);
                meanErrors[i]->Draw("pE0,csame");
                legendMean->AddEntry(meanErrors[i],nameCutVariation[i].Data(),"p");
            }
        }
        DrawGammaSetMarkerTGraphErr(graphMaterialError, 20, 1.,color[8]+1,color[8]+1);
        DrawGammaSetMarkerTGraphErr(graphMassResolError, 21, 1.,color[8]+2,color[8]+2);
        if(!meson.CompareTo("EtaToPi0")==0){
            graphMaterialError->Draw("pX0,csame");
            legendMean->AddEntry(graphMaterialError,"Material","p");
            graphMassResolError->Draw("pX0,csame");
            legendMean->AddEntry(graphMassResolError,"Mass resolution","p");
        }
    // 	DrawGammaSetMarkerTGraphErr(meanErrorsSummed, 20, 1.,1,1);
    // 	meanErrorsSummed->Draw("pE0,csame");
    // 	legendMean->AddEntry(meanErrorsSummed,"quad. sum.","p");
        legendMean->Draw();
        
    canvasSystErrors->Update();
//     canvasSystErrors->SaveAs(Form("%s/SysMean_%s_%s%s_%s.%s",outputdir.Data(),meson.Data(), energy.Data(),cent.Data(),dateForOutput.Data(),suffix.Data()));
	
    canvasSystErrors->cd();
        histo2DSystErrors->DrawCopy();
    
        TLegend* legendMeanNew;
        if (meson.Contains("Pi0")) legendMeanNew= new TLegend(0.11,0.65,0.5,0.95);
        else legendMeanNew= new TLegend(0.11,0.65,0.55,0.95);
        legendMeanNew->SetTextSize(0.035);
        legendMeanNew->SetMargin(0.1);
        legendMeanNew->SetFillColor(0);
        legendMeanNew->SetFillStyle(0);
        legendMeanNew->SetBorderSize(0);
        legendMeanNew->SetNColumns(2);
        for(Int_t i = 0; i< numberCutStudies ; i++){
            if(benable[i]){
                DrawGammaSetMarkerTGraphErr(meanErrorsCorr[i], markerStyle[i], 1.,color[i],color[i]);
                meanErrorsCorr[i]->Draw("pX0,csame");
                legendMeanNew->AddEntry(meanErrorsCorr[i],nameCutVariation[i].Data(),"p");	
            }
        }
        if(!meson.CompareTo("EtaToPi0")==0){
            graphMaterialError->Draw("pX0,csame");
            legendMeanNew->AddEntry(graphMaterialError,"Material","p");
            graphMassResolError->Draw("pX0,csame");
            legendMeanNew->AddEntry(graphMassResolError,"Mass resolution","p");
        }
        
        DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummedIncMat, 20, 1.,kBlack,kBlack);
        meanErrorsCorrSummedIncMat->Draw("p,csame");
        legendMeanNew->AddEntry(meanErrorsCorrSummedIncMat,"quad. sum.","p");
        legendMeanNew->Draw();
	
        TLatex *labelThesis = new TLatex(0.6,0.83,"This thesis");
        SetStyleTLatex( labelThesis, 0.035,4);
        labelThesis->Draw();        
        TLatex *labelMeson;
        if (meson.Contains("Pi0")) labelMeson= new TLatex(0.6,0.89,Form("#pi^{0} #rightarrow #gamma_{conv}#gamma_{conv}"));
        else labelMeson= new TLatex(0.6,0.89,Form("#eta #rightarrow #gamma_{conv}#gamma_{conv}"));
        SetStyleTLatex( labelMeson, 0.035,4);
        labelMeson->Draw();
        TLatex *labelCentrality = new TLatex(0.6,0.93,Form("%s %s",centPercent.Data(),collisionSystem.Data()	));
        SetStyleTLatex( labelCentrality, 0.035,4);
        labelCentrality->Draw();

	canvasSystErrors->Update();
	canvasSystErrors->SaveAs(Form("%s/SysMeanNewWithMean_%s_%s%s_%s.%s",outputdir.Data(),meson.Data(), energy.Data(),cent.Data(),dateForOutput.Data(),suffix.Data()));
	
    if(smooth==2){
        for(Int_t i = 0; i< numberCutStudies ; i++){
            if(benable[i]){
                canvasSystErrors->cd();
                TH2D *histo2DCheckSmooth;
                if (meson.Contains("Pi0"))histo2DCheckSmooth = new TH2D("", "", 20,0.,ptBins[nPtBins-1]+1,1000.,0.,30.);
                else histo2DCheckSmooth = new TH2D("", "", 20,0.,ptBins[nPtBins-1]+1,1000.,0.,65.);
                SetStyleHistoTH2ForGraphs( histo2DCheckSmooth, "#it{p}_{T} (GeV/#it{c})", "mean smoothed systematic Err %",
                                    0.85*55./1200, 55./1200, 0.85*55./1200, 55./1200,0.9,0.75);
//                 histo2DCheckSmooth->GetYaxis()->SetLabelOffset(0.001);
//                 histo2DCheckSmooth->GetXaxis()->SetLabelOffset(-0.01);
        //      histo2DCheckSmooth->GetXaxis()->SetMoreLogLabels(kTRUE);
                histo2DCheckSmooth->DrawCopy();
                    
                DrawGammaSetMarkerTGraphErr(meanErrorsCorr[i], markerStyle[i], 1.,color[i],color[i]);
                meanErrorsCorr[i]->Draw("pX0,csame");

                DrawGammaSetMarkerTGraphErr(meanErrorsCorrNotSmooth[i], markerStyle[i], 1.,kGray+2,kGray+2);
                meanErrorsCorrNotSmooth[i]->Draw("pX0,csame");
                    
                canvasSystErrors->SaveAs(Form("%s/%sSysBeforeSmoothing_%s_%s_%s_%s.%s",outputsmooth.Data(), nameCutVariationSC[i].Data(), meson.Data(), energy.Data(),cent.Data(),dateForOutput.Data(),suffix.Data()));
            }
        }
    }
    
	if(smooth==1){
		canvasSystErrors->cd();
		histo2DSystErrors->Draw();
			
		for(Int_t i = 0; i< numberCutStudies ; i++){
            if(benable[i]){
                DrawGammaSetMarkerTGraphErr(meanErrorsCorrSmoothed[i], markerStyle[i], 1.,color[i],color[i]);
                meanErrorsCorrSmoothed[i]->Draw("pX0,csame");
        // 		legendMeanNew->AddEntry(meanErrorsCorr[i],nameCutVariation[i].Data(),"p");	
            }
		}
        if(meson.CompareTo("EtaToPi0")){
            graphMassResolError->Draw("pX0,csame");
            graphMaterialError->Draw("pX0,csame");
        }
		DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummedIncMatSmoothed, 20, 1.,kBlack,kBlack);
		meanErrorsCorrSummedIncMatSmoothed->Draw("p,csame");

		labelMeson->Draw();
		labelCentrality->Draw();
        legendMeanNew->Draw();
        
        TLatex *labelThesis = new TLatex(0.7,0.83,"This thesis");
        SetStyleTLatex( labelThesis, 0.038,4);
        labelThesis->Draw();

		canvasSystErrors->Update();
		canvasSystErrors->SaveAs(Form("%s/SysMeanNewWithMeanSmoothed_%s_%s%s_%s.%s",outputsmooth.Data(),meson.Data(), energy.Data(),cent.Data(),dateForOutput.Data(),suffix.Data()));
	}	
// 	
	Double_t errorsMeanCorrPID[nPtBins];	
	Double_t errorsMeanCorrSignalExtraction[nPtBins];	
	Double_t errorsMeanCorrTrackReco[nPtBins];	
	Double_t errorsMeanCorrPhotonReco[nPtBins];	
	Double_t errorsMeanCorrAcceptance[nPtBins];			
	Double_t errorsMeanCorrOther[nPtBins];	
    Double_t errorsMeanCorrPileUp[nPtBins];
	
    //"YieldExtraction"-0,"dEdxE"-1,"dEdxPi"-2, "TPCCluster"-3, "SinglePt"-4, "Chi2"-5, "Qt"-6, "Alpha"-7, "ConvPhi"-8, "YieldsfotEtaToPi0" -9
	for (Int_t l=0; l< nPtBins; l++){
        if(smooth==1){
            if (meson.CompareTo("EtaToPi0")==0){
                // Signal extraction: Yield extraction 0, Alpha 7, and secondyields extraction for eta to pi0
                errorsMeanCorrSignalExtraction[l] = TMath::Sqrt(errorsMeanCorr[0][l]*errorsMeanCorr[0][l]+errorsMeanCorrSmoothed[7][l]*errorsMeanCorrSmoothed[7][l]+errorsMeanCorr[10][l]*errorsMeanCorr[10][l]+errorMassRes[l]*errorMassRes[l]);	
            } else {
                // Signal extraction: Yield extraction 0, Alpha 7, mass resolution 
                errorsMeanCorrSignalExtraction[l] = TMath::Sqrt(errorsMeanCorr[0][l]*errorsMeanCorr[0][l]+errorsMeanCorrSmoothed[7][l]*errorsMeanCorrSmoothed[7][l]+errorMassRes[l]*errorMassRes[l]);	
            }
            // PID: dEdxE 1, dEdxPi 2
            errorsMeanCorrPID[l] =TMath::Sqrt(errorsMeanCorrSmoothed[1][l]*errorsMeanCorrSmoothed[1][l]+ errorsMeanCorrSmoothed[2][l]*errorsMeanCorrSmoothed[2][l]);
            // photon reco: Chi2 5, Qt 6
            errorsMeanCorrPhotonReco[l] =TMath::Sqrt( errorsMeanCorrSmoothed[5][l]* errorsMeanCorrSmoothed[5][l]+errorsMeanCorrSmoothed[6][l]*errorsMeanCorrSmoothed[6][l]);
            // track reconstruction: TPCCluster 3, Single pt 4, ConvPhi 8
            errorsMeanCorrTrackReco[l] = TMath::Sqrt(errorsMeanCorrSmoothed[3][l]*errorsMeanCorrSmoothed[3][l]+errorsMeanCorrSmoothed[4][l]*errorsMeanCorrSmoothed[4][l]+errorsMeanCorrSmoothed[8][l]*errorsMeanCorrSmoothed[8][l]);
            // Pile-up
            if(benable[9]) errorsMeanCorrPileUp[l] = TMath::Sqrt(errorsMeanCorrSmoothed[9][l]*errorsMeanCorrSmoothed[9][l]);            
        } 
        if(smooth!=1){
            // Signal extraction: Yield extraction 0, Alpha 7, and secondyields extraction for eta to pi0
            if (meson.CompareTo("EtaToPi0")==0){
                errorsMeanCorrSignalExtraction[l] = TMath::Sqrt(errorsMeanCorr[0][l]*errorsMeanCorr[0][l]+errorsMeanCorr[7][l]*errorsMeanCorr[7][l]+errorsMeanCorr[10][l]*errorsMeanCorr[10][l]+errorMassRes[l]*errorMassRes[l]);	
            } else {
                errorsMeanCorrSignalExtraction[l] = TMath::Sqrt(errorsMeanCorr[0][l]*errorsMeanCorr[0][l]+errorsMeanCorr[7][l]*errorsMeanCorr[7][l]+errorMassRes[l]*errorMassRes[l]);	
            }
            // PID: dEdxE 1, dEdxPi 2
            errorsMeanCorrPID[l] =TMath::Sqrt(errorsMeanCorr[1][l]*errorsMeanCorr[1][l]+ errorsMeanCorr[2][l]*errorsMeanCorr[2][l]);
            // photon reco: Chi2 5, Qt 6
            errorsMeanCorrPhotonReco[l] =TMath::Sqrt( errorsMeanCorr[5][l]* errorsMeanCorr[5][l]+errorsMeanCorr[6][l]*errorsMeanCorr[6][l]);
            // track reconstruction: TPCCluster 3, Single pt 4, ConvPhi 8
            errorsMeanCorrTrackReco[l] = TMath::Sqrt(errorsMeanCorr[3][l]*errorsMeanCorr[3][l]+errorsMeanCorr[4][l]*errorsMeanCorr[4][l]+errorsMeanCorr[8][l]*errorsMeanCorr[8][l]);
            // Pile-up
            if(benable[9]) errorsMeanCorrPileUp[l] = TMath::Sqrt(errorsMeanCorr[9][l]*errorsMeanCorr[9][l]);
        }
	}
		
		
	TGraphErrors* meanErrorsPID = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrPID ,ptBinsErr ,errorsMeanErrCorrSummed );
	TGraphErrors* meanErrorsPhotonReco = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrPhotonReco ,ptBinsErr ,errorsMeanErrCorrSummed );
	TGraphErrors* meanErrorsSignalExtraction = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrSignalExtraction ,ptBinsErr ,errorsMeanErrCorrSummed );
	TGraphErrors* meanErrorsTrackReco = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrTrackReco ,ptBinsErr ,errorsMeanErrCorrSummed );
	TGraphErrors* meanErrorsPileUp = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrPileUp ,ptBinsErr ,errorsMeanErrCorrSummed );
// 	TGraphErrors* meanErrorsAcceptance = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrAcceptance ,ptBinsErr ,errorsMeanErrCorrSummed );
// 	TGraphErrors* meanErrorsOther = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrOther ,ptBinsErr ,errorsMeanErrCorrSummed );

    canvasSystErrors->cd();
        histo2DSystErrors->DrawCopy();
	
        TLegend* legendSummedMeanNew;
        if (meson.Contains("Pi0")) legendSummedMeanNew= new TLegend(0.13,0.7,0.65,0.95);
        else legendSummedMeanNew= new TLegend(0.20,0.7,0.7,0.95);
        legendSummedMeanNew->SetTextSize(0.035);
        legendSummedMeanNew->SetFillColor(0);
        legendSummedMeanNew->SetFillStyle(0);	
        legendSummedMeanNew->SetBorderSize(0);
        legendSummedMeanNew->SetNColumns(2);
        legendSummedMeanNew->SetMargin(0.1);
	
        // Signal extraction error
        DrawGammaSetMarkerTGraphErr(meanErrorsSignalExtraction, 20, 1.,color[0],color[0]);
        meanErrorsSignalExtraction->Draw("pX0,csame");
        // electron PID error
        DrawGammaSetMarkerTGraphErr(meanErrorsPID, 21, 1.,color[1],color[1]);
        meanErrorsPID->Draw("pX0,csame");
        // track reconstruction error
        DrawGammaSetMarkerTGraphErr(meanErrorsTrackReco, 22, 1.,color[2],color[2]);
        meanErrorsTrackReco->Draw("pX0,csame");
        // Photon reco error
        DrawGammaSetMarkerTGraphErr(meanErrorsPhotonReco, 23, 1.,color[3],color[3]);
        meanErrorsPhotonReco->Draw("pX0,csame");
        // Eta error
    // 	DrawGammaSetMarkerTGraphErr(meanErrorsAcceptance, 28, 1.,color[5],color[5]);
    // 	meanErrorsAcceptance->Draw("pX0,csame");
        // Other error
    // 	DrawGammaSetMarkerTGraphErr(meanErrorsOther, 26, 1.,color[6],color[6]);
    // 	meanErrorsOther->Draw("pX0,csame");
        //pile-up
        DrawGammaSetMarkerTGraphErr(meanErrorsPileUp, 26, 1.,color[6],color[6]);
        if(benable[9]) meanErrorsPileUp->Draw("pX0,csame");
        // PCM Material budget 
        if(!meson.CompareTo("EtaToPi0")==0){
            DrawGammaSetMarkerTGraphErr(graphMaterialError, 24, 1.,color[4],color[4]);
            graphMaterialError->Draw("pX0,csame");
        }
        // PCM Mass Resolution
    // 	DrawGammaSetMarkerTGraphErr(graphMassResolError, 21, 1.,color[11],color[11]);
    // 	graphMassResolError->Draw("pX0,csame");
        
        legendSummedMeanNew->AddEntry(meanErrorsSignalExtraction,"Signal Extraction","p");
        legendSummedMeanNew->AddEntry(meanErrorsPID,"Electron PID","p");
        legendSummedMeanNew->AddEntry(meanErrorsTrackReco,"Track Reconstruction","p");
        legendSummedMeanNew->AddEntry(meanErrorsPhotonReco,"Photon Reconstruction","p");
    // 	legendSummedMeanNew->AddEntry(meanErrorsAcceptance,"Acceptance (#eta)","p");
        if(!meson.CompareTo("EtaToPi0")==0) legendSummedMeanNew->AddEntry(graphMaterialError,"Material","p");
    // 	legendSummedMeanNew->AddEntry(graphMassResolError,"Mass resolution","p");	
    // 	legendSummedMeanNew->AddEntry(meanErrorsOther,"TOF, cos(P.A.)","p");
        if(benable[9]) legendSummedMeanNew->AddEntry(meanErrorsPileUp,"Pile-up","p");
        
        DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummedIncMat, 20, 1.,kBlack,kBlack);
        meanErrorsCorrSummedIncMat->Draw("p,csame");
        legendSummedMeanNew->AddEntry(meanErrorsCorrSummedIncMat,"quad. sum.","p");
        legendSummedMeanNew->Draw();
        
        labelMeson->Draw();
        labelCentrality->Draw();
        
    canvasSystErrors->Update();
    canvasSystErrors->SaveAs(Form("%s/SysErrorSummedVisu_%s_%s%s_%s.%s",outputdir.Data(),meson.Data(), energy.Data(),cent.Data(),dateForOutput.Data(),suffix.Data()));	
	delete canvasSystErrors;

    cout << " \n meanErrorsTrackReco:" << endl;
    meanErrorsTrackReco->Print();
    cout << " \n meanErrorsSignalExtraction:" << endl;
    meanErrorsSignalExtraction->Print();
    cout << " \n meanErrorsPID:" << endl;
    meanErrorsPID->Print();
    cout << " \n meanErrorsPhotonReco:" << endl;
    meanErrorsPhotonReco->Print();
    if(benable[9]){
        cout << " \n meanErrorsPileUp:" << endl;
        meanErrorsPileUp->Print();
    }
    if(smooth==1){
        cout << " \n meanErrorsCorrSummedIncMatSmoothed:" << endl;
        meanErrorsCorrSummedIncMatSmoothed->Print();        
    } else {
        cout << " \n meanErrorsCorrSummedIncMat:" << endl;
        meanErrorsCorrSummedIncMat->Print();        
    }


	const char *SysErrDatname = Form("SystematicErrorsCalculatedGammaPbPb_LHC11h_%s/SystematicError_%s_%s%s_%s.dat",dateForOutput.Data(),meson.Data(),energy.Data(),cent.Data(),dateForOutput.Data());
	fstream SysErrDat;
	cout << SysErrDatname << endl;
	SysErrDat.open(SysErrDatname, ios::out);
	for (Int_t l=0; l< nPtBins; l++){
		SysErrDat <<errorsNegCorrMatSummed[l] << "\t" <<errorsPosCorrMatSummed[l] << "\t"  <<errorsNegCorrSummed[l] << "\t" <<errorsPosCorrSummed[l]  << endl;
	}
	SysErrDat.close();

	const char *SysErrDatnameMean = Form("SystematicErrorsCalculatedGammaPbPb_LHC11h_%s/SystematicErrorAveraged_%s_%s%s_%s.dat",dateForOutput.Data(),meson.Data(),energy.Data(),cent.Data(),dateForOutput.Data());
	fstream SysErrDatAver;
	cout << SysErrDatnameMean << endl;
	SysErrDatAver.open(SysErrDatnameMean, ios::out);
	for (Int_t l=0; l< nPtBins; l++){
		if(smooth){
			SysErrDatAver << "-"<< errorsMeanCorrMatSummedSmoothed[l] << "\t" <<errorsMeanCorrMatSummedSmoothed[l] << "\t"  << "-"<< errorsMeanCorrSummedSmoothed[l] << "\t" <<errorsMeanCorrSummedSmoothed[l]  << endl;
		} else {
			SysErrDatAver << "-"<< errorsMeanCorrMatSummed[l] << "\t" <<errorsMeanCorrMatSummed[l] << "\t"  << "-"<< errorsMeanCorrSummed[l] << "\t" <<errorsMeanCorrSummed[l]  << endl;
		}
	}
	SysErrDatAver.close();

}