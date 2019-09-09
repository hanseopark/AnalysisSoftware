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
#include "CommonHeaders/ExtractSignalBinning.h"

void FinaliseSystematicErrorsConvPbPb(   const char* nameDataFileErrors  = "", 
                                         TString energy                  = "", 
                                         TString meson                   = "", 
                                         Int_t numberOfPtBins            = 1 ,   // for size of arrays
                                         Int_t offSetBeginning           = 0,    // remove number of points at the beginning of SystError graphs
                                         TString centPercent             = "0-5%",
                                         TString cent                    = "0005",
                                         TString suffix                  = "pdf", 
                                         Int_t smooth                    = 2, 
                                         Int_t mode                      = 0,
                                         Int_t verbose                   = 0
    ){
	
	StyleSettingsThesis();	
	SetPlotStyle();
	
	TString date = ReturnDateString();
	TString dateForOutput = ReturnDateStringForOutput();
	TString collisionSystem= ReturnFullCollisionsSystem(energy);
	TString centralityRead = cent;

	Int_t 	numberOfEntriesPos = 0;
	Int_t 	numberOfEntriesNeg = 0;
	
	TFile* fileErrorInput = new TFile(nameDataFileErrors);
	const Int_t nPtBins   = numberOfPtBins;
	const Int_t nCuts     = 11;  // for size of arrays   
	Double_t* ptBins      = NULL;
	Double_t* ptBinsErr   = NULL;
    
    // get desired pT binning
	Double_t ptBinsNew[50];
    Double_t ptBinsCentersNew[50];
    Double_t ptBinsNewErr[50];
    Int_t ptBinsMax = 0;
    TString mesonForGetBinning;
    Int_t startBin = 1;
    if (meson.Contains("Eta")) mesonForGetBinning = "Eta";
    else mesonForGetBinning = "Pi0";
    GetBinning( ptBinsNew, ptBinsMax, mesonForGetBinning, energy, mode, -1, kFALSE, centPercent, kFALSE);
    if(numberOfPtBins>ptBinsMax) {
      cout << "ERROR: too many pT bins requested" << endl;
      return;
    }
    const Int_t nPtBinsNew = ptBinsMax-startBin;
    cout << "will evaluate smoothing in the following binning: " << endl;
    for(Int_t i=0; i<(ptBinsMax+1); i++){
        cout << ptBinsNew[i] << endl;
    }
    cout << "with bin centers and widhts:" << endl;
    for(Int_t i=0; i<nPtBinsNew; i++){
        ptBinsNewErr[i]     = (ptBinsNew[i+1+startBin]-ptBinsNew[i+startBin])/2;
        ptBinsCentersNew[i] = (ptBinsNew[i+1+startBin]+ptBinsNew[i+startBin])/2;
        cout << ptBinsCentersNew[i] << " " << ptBinsNewErr[i] << endl;
    }

    // cut variations
	TString nameCutVariation[nCuts] = {"Yield extract.", "dE/dx e-line",              "dE/dx #pi-line", 
                                       "TPC cluster",    "Single e^{#pm} #it{p}_{T}", "2D #chi^{2} #gamma, #psi_{pair} #gamma",
                                       "2D q_{T}",	  "#alpha meson",              "#varphi_{conv}",
                                       "Pile-up",        "Yield extract. #pi^{0}"   
                                       //"#eta",	  "Max #pi momentum",           "TOF", "cosine point. angle"
	};
	// in case of EtaToPi0: [0]: Eta, [10]: Pi0EtaBinning yield extraction
	
	if (meson.CompareTo("EtaToPi0") == 0) nameCutVariation[0]          = "Yield extract. #eta";

	// names for default colors
	TString nameCutVariationSC[nCuts] = {"YieldExtraction", "dEdxE", "dEdxPi", 
                                         "TPCCluster", "SinglePt", "Chi2", 
                                         "Qt", "Alpha", "ConvPhi",
                                         "Pileup", "YieldExtraction"
                                         //"Eta", "pdEdxPi", "TOF", "CosPoint"
	};
	
	// Set colors and markers
	Color_t color[nCuts];
	Color_t markerStyle[nCuts];	
	for (Int_t k =0; k<nCuts; k++ ){
        color[k]                                = GetColorSystematics( nameCutVariationSC[k] ); 
        if(k==10) color[k]                      = GetColorSystematics("Rapidity"); 
        markerStyle[k]                          = GetMarkerStyleSystematics( nameCutVariationSC[k] );     
        // nameCutVariation[k]                  = GetSystematicsName(nameCutVariationSC[k]);
	}
		
	Bool_t benable[nCuts]            = { 0, 0, 0, 0, 0, 
                                         0, 0, 0, 0, 0,
                                         0};

	Bool_t benablePi0[nCuts]         = { 1, 1, 1, 1, 1,   
                                         1, 1, 0, 0, 0,
                                         0};

    Bool_t benableEta[nCuts]         = { 1, 1, 1, 1, 1,   
                                         1, 1, 0, 0, 0,
                                         0};
	
	Bool_t benableRatio[nCuts]       = { 1, 1, 1, 1, 0,  // disabled SinglePt for now    
                                         1, 1, 0, 0, 0,
                                         1};


	// END OF SETTINGS
	
	
    if(meson.CompareTo("EtaToPi0")==0){
        if( (centPercent.CompareTo("0-10%") == 0) || (centPercent.CompareTo("60-80%") == 0) ){
            for (Int_t i = 0; i < nCuts; i++){
                benable[i] = benableRatio[i];
            }
        } else {
            benable[0] = 1;                       // for all cent classes except 0-10% and 60-80% enable only yield extraction for Eta & Pi0EtaBinning
            benable[nCuts-1] = 1;
            for (Int_t i = 1; i < (nCuts-1); i++){
                benable[i] = 0;
            }
        }
    } else if(meson.CompareTo("Pi0")==0){
        if( (centPercent.CompareTo("0-10%") == 0) || (centPercent.CompareTo("60-80%") == 0) ){
            for (Int_t i = 0; i < nCuts; i++){
                benable[i] = benablePi0[i];
            }
        } else {
            benable[0] = 1;                       // for all cent classes except 0-10% and 60-80% enable only yield extraction
            for (Int_t i = 1; i < nCuts; i++){
                benable[i] = 0;
            }
        }
    } else if(meson.CompareTo("Eta")==0){
        if( (centPercent.CompareTo("0-10%") == 0) || (centPercent.CompareTo("60-80%") == 0) ){
            for (Int_t i = 0; i < nCuts; i++){
                benable[i] = benableEta[i];
            }
        } else {
            benable[0] = 1;                       
            for (Int_t i = 1; i < nCuts; i++){
                benable[i] = 0;
            }
        }
    }


	// PLOTTING SETTINGS
	Double_t maxPtPlot  = 15.0;   
	Double_t maxErrPlot = 30.0;
    Double_t maxErrPlotIncMat = 30.0;
    Double_t maxErrPlotEta = 65.;
    Double_t maxErrPlotEtaIncMat = 80.;
	if(!centPercent.CompareTo("60-80%")){
        maxErrPlot = 10.0;
        maxErrPlotIncMat = 20;
    }
    
	// FITTING SETTINGS
	Double_t minPtFit = 0.;   
	Double_t maxPtFit = 0.;
	if(meson.Contains("Pi0")){
        minPtFit = 0.4;
        maxPtFit = 15.;
	} else if(meson.Contains("Eta")){
        minPtFit = 1.0;
        maxPtFit = 12.;
	}
	
	TString outputdirName = "MesonSystematicErrorsCalculatedConv"; // SystematicErrorsCalculatedGammaPbPb
	TString outputdir = Form("%s_%s/%s",outputdirName.Data(), dateForOutput.Data(),meson.Data()); 
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
	Double_t errorsMeanCorrSummedSmoothed[nPtBinsNew];  // Meike
	Double_t errorsMeanCorrMatSummedSmoothed[nPtBinsNew];
	Double_t errorsMeanErrCorrSummedSmoothed[nPtBinsNew];
	Double_t errorsMeanErrCorrMatSummedSmoothed[nPtBinsNew];
	
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

    }
    for (Int_t l = 0; l < nPtBinsNew; l++){
        errorsMeanCorrSummedSmoothed[l] = 0.;
    }
    
	// go throught cut variations 
	for (Int_t i = 0; i < nCuts; i++){
        TGraphAsymmErrors* graphPosErrors       = NULL;
        TGraphAsymmErrors* graphNegErrors       = NULL;
        TString nameGraphPos    = "";
        TString nameGraphNeg    = "";

        // names of graphs to be read from input root file
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
                    nameGraphPos = Form("%s_SystErrorRelPos_%s_%s",meson.Data(),nameCutVariationSC[i].Data(),centPercent.Data()  );
                    nameGraphNeg = Form("%s_SystErrorRelNeg_%s_%s",meson.Data(),nameCutVariationSC[i].Data(),centPercent.Data()  );
                }
            } else {
                if (i==0){
                    nameGraphPos        = Form("%s_SystErrorRelPos_YieldExtraction_%s",meson.Data(),centPercent.Data() );
                    nameGraphNeg        = Form("%s_SystErrorRelNeg_YieldExtraction_%s",meson.Data(),centPercent.Data() );
                } else if (i==9) {
                    nameGraphPos        = Form("%s_SystErrorRel_BGEstimate_%s",meson.Data(),centPercent.Data() );
                    nameGraphNeg        = Form("%s_SystErrorRel_BGEstimate_%s",meson.Data(),centPercent.Data() );
                } else {
                    nameGraphPos = Form("%s_SystErrorRelPos_%s_%s",meson.Data(),nameCutVariationSC[i].Data(),centPercent.Data()  );
                    nameGraphNeg = Form("%s_SystErrorRelNeg_%s_%s",meson.Data(),nameCutVariationSC[i].Data(),centPercent.Data()  );
                }
            }
            cout << "Cutstudies " << i<< "\t" <<nameGraphPos.Data() << "\t" << nameGraphNeg.Data()<<  endl;

            // read error graph from root file
            graphPosErrors = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphPos.Data());
            graphNegErrors = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphNeg.Data());
            
            for (Int_t j = 0; j < offSetBeginning; j++){
                Double_t x, y;
                graphPosErrors->GetPoint(0,x,y);
                cout << "Remove point " << " " << x << "\t " << y <<  endl;
                graphPosErrors->RemovePoint(0);
                graphNegErrors->RemovePoint(0);
            } cout << "done" << endl;
            // fill arrays with points from graphs
            if (i == 0) {
                ptBins = graphNegErrors->GetX();
                ptBinsErr = graphNegErrors->GetEXhigh();
            }
	    
            errorsNeg[i] = graphNegErrors->GetY();
            errorsNegErr[i] = graphNegErrors->GetEYhigh();
            errorsPos[i] = graphPosErrors->GetY();
            errorsPosErr[i] = graphPosErrors->GetEYhigh();
                
            // Calculate mean systematic error from pos and neg
            CalculateMeanSysErr(errorsMean[i], errorsMeanErr[i], errorsPos[i], errorsNeg[i], nPtBins);
            // if 0, use mean of previous and next pT bin
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
                    if((meson.CompareTo("EtaToPi0")==0)){
                        for (Int_t k = 0; k < nPtBins; k++){
                            if(!centPercent.CompareTo("0-5%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 1+1.5*(ptBins[k]-2)+20/pow(ptBins[k]+0.2,2);
                            } 
                            else if(!centPercent.CompareTo("5-10%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 1+1.5*(ptBins[k]-2)+20/pow(ptBins[k]+0.2,2);
                            } 
                            else if (!centPercent.CompareTo("0-10%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 1+1.5*(ptBins[k]-2)+20/pow(ptBins[k]+0.2,2);
                            }
                            else if(!centPercent.CompareTo("20-40%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 1+1.5*(ptBins[k]-2)+20/pow(ptBins[k]+0.2,2);
                            }
                            else if(!centPercent.CompareTo("20-50%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 1+1.5*(ptBins[k]-2)+20/pow(ptBins[k]+0.2,2);
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
                    }
                    else {
                        for (Int_t k = 0; k < nPtBins; k++){
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
                                errorFixed      = 9+1.3*(ptBins[k]-6)+10/pow(ptBins[k],1.5);
                            } 
                            else if(!centPercent.CompareTo("5-10%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 8+1.3*(ptBins[k]-6)+15/pow(ptBins[k],1.5);
                            } 
                            else if (!centPercent.CompareTo("0-10%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 9+1.3*(ptBins[k]-6)+10/pow(ptBins[k],1.5);
                            }
                            else if(!centPercent.CompareTo("20-40%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 8+1.3*(ptBins[k]-6)+10/pow(ptBins[k],1.5);
                            }
                            else if(!centPercent.CompareTo("20-50%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 9+1.*(ptBins[k]-6)+10/pow(ptBins[k],1.5);
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
                                errorFixed      = 10+1.3*(ptBins[k]-6)+10/pow(ptBins[k],1.5);
                            } 
                            else if(!centPercent.CompareTo("5-10%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 8+1.2*(ptBins[k]-6)+20/pow(ptBins[k],1.5);
                            } 
                            else if (!centPercent.CompareTo("0-10%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 9+1.3*(ptBins[k]-6)+10/pow(ptBins[k],1.5);
                            }
                            else if(!centPercent.CompareTo("20-40%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 7+1.1*(ptBins[k]-6)+20/pow(ptBins[k],1.5);
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
                                errorFixed = 4.+0.3*ptBins[k]+15/(pow(ptBins[k]+0.3,2));
                            } 
                            else if(!centPercent.CompareTo("5-10%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed = 5.5+0.3*ptBins[k]+8/(pow(ptBins[k]+0.3,2));
                            } 
                            else if (!centPercent.CompareTo("0-10%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed = 2+0.3*ptBins[k]+20/(pow(ptBins[k]+0.3,2));
                            }
                            else if(!centPercent.CompareTo("20-40%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed = 2+0.3*ptBins[k]+18/(pow(ptBins[k]+0.3,2));
                            }
                            else if(!centPercent.CompareTo("20-50%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed = 1.8+0.3*ptBins[k]+20/(pow(ptBins[k]+0.3,2));
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
                                errorFixed      = 4+pow(ptBins[k],2)*0.07;
                            } 
                            else if (!centPercent.CompareTo("0-10%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 3+pow(ptBins[k],2)*0.07;
                            }
                            else if(!centPercent.CompareTo("20-40%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 3+pow(ptBins[k],2)*0.07;
                            }
                            else if(!centPercent.CompareTo("20-50%")){
                                adjustPtDependent               = kTRUE;
                                errorFixed      = 2+pow(ptBins[k],2)*0.05;
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
                                errorFixed = 1.+0.2*ptBins[k]+23/(pow(ptBins[k]+1.5,1.5));
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
                } else if(i==0 && meson.CompareTo("Pi0")==0 && !centPercent.CompareTo("0-10%")){
                    for (Int_t k = 0; k < nPtBins; k++){
                        if(ptBins[k]==9){
                            adjustPtDependent               = kFALSE;
                            errorFixed = 6.;
                        } else if(ptBins[k]==11){
                            adjustPtDependent               = kFALSE;
                            errorFixed = 7.;
                        } else errorFixed = -1;
                            
                        if (errorFixed != -1){
                            errorsMean[i][k]        = errorFixed;
                            errorsMeanErr[i][k]     = errorFixed*0.01;
                            errorsMeanCorr[i][k]    = errorFixed;
                            errorsMeanErrCorr[i][k] = errorFixed*0.01;
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

            // create graphs from arrays
            negativeErrors[i]             = new TGraphErrors(nPtBins,ptBins ,errorsNeg[i] ,ptBinsErr ,errorsNegErr[i] );
            meanErrors[i]                 = new TGraphErrors(nPtBins,ptBins ,errorsMean[i] ,ptBinsErr ,errorsMeanErr[i] );
            positiveErrors[i]             = new TGraphErrors(nPtBins,ptBins ,errorsPos[i] ,ptBinsErr ,errorsPosErr[i] );
            negativeErrorsCorr[i]         = new TGraphErrors(nPtBins,ptBins ,errorsNegCorr[i] ,ptBinsErr ,errorsNegErrCorr[i] );
            meanErrorsCorr[i]             = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorr[i] ,ptBinsErr ,errorsMeanErrCorr[i] );
            positiveErrorsCorr[i]         = new TGraphErrors(nPtBins,ptBins ,errorsPosCorr[i] ,ptBinsErr ,errorsPosErrCorr[i] );
            meanErrorsCorrNotSmooth[i]    = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrNotSmooth[i] ,ptBinsErr ,errorsMeanErrCorrNotSmooth[i] );

            meanErrorsCorrSmoothed[i]     = new TGraphErrors(nPtBinsNew, ptBinsCentersNew, errorsMeanCorr[i], ptBinsNewErr, errorsMeanErrCorr[i]); // values just for initialization
	    
        }
	}
	
	//smoothing with the fit function:
	if(smooth==1){

        // plot to check smoothing procedure
        TCanvas* canvasCheckSmooth = new TCanvas("canvasCheckSmooth","",200,10,1350,900);  // gives the page size
        DrawGammaCanvasSettings( canvasCheckSmooth, 0.08, 0.01, 0.015, 0.09);
	  
        TH2D *histo2DCheckSmooth;
        if (meson.Contains("Pi0")){
            histo2DCheckSmooth = new TH2D("histo2DCheckSmooth", "histo2DCheckSmooth", maxPtPlot,0.,maxPtPlot,1000.,-0.5,maxErrPlot);
        } else { 
            histo2DCheckSmooth = new TH2D("histo2DCheckSmooth", "histo2DCheckSmooth", maxPtPlot,0.,maxPtPlot,1000.,-0.5,maxErrPlotEta);
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


        cout << "Fit between " << minPtFit<< " and " << maxPtFit << endl;
        TF1* pol1  = new TF1("pol1", "[0]+[1]*x",minPtFit,maxPtFit);
        TF1* pol2  = new TF1("pol2", "[0]+[1]*x+[2]*x*x",minPtFit,maxPtFit);
        TF1* pol3  = new TF1("pol3", "[0]+[1]*x+[2]*x*x+[3]*x*x*x",minPtFit,maxPtFit); 
        TF1* pol4  = new TF1("pol4", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x",minPtFit,maxPtFit);
        TF1* pol5  = new TF1("pol5", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x+[5]*x*x*x*x*x",minPtFit,maxPtFit);
        TF1* pol6  = new TF1("pol6", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x+[5]*x*x*x*x*x+[6]*x*x*x*x*x*x",minPtFit,maxPtFit);
        TF1* pol7  = new TF1("pol7", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x+[5]*x*x*x*x*x+[6]*x*x*x*x*x*x+[7]*x*x*x*x*x*x*x",minPtFit,maxPtFit);
        TF1* pol8  = new TF1("pol8", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x+[5]*x*x*x*x*x+[6]*x*x*x*x*x*x+[7]*x*x*x*x*x*x*x+[8]*x*x*x*x*x*x*x*x",minPtFit,maxPtFit);
        TF1* pol2N = new TF1("pol2N","[0]+[1]*x+[2]*x*x+[3]/x",minPtFit,maxPtFit);
        TF1* pol3N = new TF1("pol3N","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]/x",minPtFit,maxPtFit);
        TF1* pol4N = new TF1("pol4N","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x+[5]/x",minPtFit,maxPtFit);

        TString pol1Label  = "c+a#upointx";
        TString pol2Label  = "c+a#upointx+b#upointx^{2}";
        TString pol3Label  = "c+a#upointx+b#upointx^{2}+c#upointx^{3}";
        TString pol4Label  = "c+a#upointx+b#upointx^{2}+c#upointx^{3}+d#upointx^{4}";
        TString pol5Label  = "c+a#upointx+b#upointx^{2}+c#upointx^{3}+d#upointx^{4}+e#upointx^{5}";
        TString pol6Label  = "c+a#upointx+b#upointx^{2}+c#upointx^{3}+d#upointx^{4}+e#upointx^{5}+f#upointx^{6}";
        TString pol7Label  = "c+a#upointx+b#upointx^{2}+c#upointx^{3}+d#upointx^{4}+e#upointx^{5}+f#upointx^{6}+g#upointx^{7}";
        TString pol8Label  = "c+a#upointx+b#upointx^{2}+c#upointx^{3}+d#upointx^{4}+e#upointx^{5}+f#upointx^{6}+g#upointx^{7}+h#upointx^{8}";
        TString pol2NLabel = "c+a#upointx+b#upointx^{2}+c#upointx^{-2}";
        TString pol3NLabel = "c+a#upointx+b#upointx^{2}+c#upointx^{3}+d#upointx^{-1}";
        TString pol4NLabel = "c+a#upointx+b#upointx^{2}+c#upointx^{3}+d#upointx^{4}+e#upointx^{-1}";

        pol1->SetLineColor(kMagenta+2);
        pol2->SetLineColor(kGreen+2);
        pol3->SetLineColor(kBlue+1);
        pol4->SetLineColor(kRed+2);
        pol5->SetLineColor(kGreen+2);
        pol6->SetLineColor(kOrange+7);
        pol7->SetLineColor(kGray+1);
        pol8->SetLineColor(kCyan+2);
        pol2N->SetLineColor(kGreen+2);
        pol3N->SetLineColor(kMagenta+2);
        pol4N->SetLineColor(kGray+3);

        Bool_t usePol1[nCuts]  = {kFALSE};
        Bool_t usePol2[nCuts]  = {kFALSE};
        Bool_t usePol3[nCuts]  = {kFALSE};
        Bool_t usePol4[nCuts]  = {kFALSE};
        Bool_t usePol5[nCuts]  = {kFALSE};
        Bool_t usePol6[nCuts]  = {kFALSE};
        Bool_t usePol7[nCuts]  = {kFALSE};
        Bool_t usePol8[nCuts]  = {kFALSE};
        Bool_t usePol2N[nCuts] = {kFALSE};
        Bool_t usePol3N[nCuts] = {kFALSE};
        Bool_t usePol4N[nCuts] = {kFALSE};
	  
        Double_t newPoint = 0;
	  
        for(Int_t cut = 0; cut< nCuts ; cut++){
            if(benable[cut]){
	      
                cout << "Smoothing with fit function for " << nameCutVariationSC[cut].Data() << endl;

                meanErrorsCorr[cut]->Fit(pol1, "QNRMEX0+","",minPtFit,maxPtFit);
                meanErrorsCorr[cut]->Fit(pol2, "QNRMEX0+","",minPtFit,maxPtFit);  
                meanErrorsCorr[cut]->Fit(pol3, "QNRMEX0+","",minPtFit,maxPtFit);
                meanErrorsCorr[cut]->Fit(pol4, "QNRMEX0+","",minPtFit,maxPtFit);
                meanErrorsCorr[cut]->Fit(pol5, "QNRMEX0+","",minPtFit,maxPtFit);
                meanErrorsCorr[cut]->Fit(pol6, "QNRMEX0+","",minPtFit,maxPtFit);
                meanErrorsCorr[cut]->Fit(pol7, "QNRMEX0+","",minPtFit,maxPtFit);
                meanErrorsCorr[cut]->Fit(pol8, "QNRMEX0+","",minPtFit,maxPtFit);
                meanErrorsCorr[cut]->Fit(pol2N,"QNRMEX0+","",minPtFit,maxPtFit);
                meanErrorsCorr[cut]->Fit(pol3N,"QNRMEX0+","",minPtFit,maxPtFit);
                meanErrorsCorr[cut]->Fit(pol4N,"QNRMEX0+","",minPtFit,maxPtFit);
	      
                //meanErrorsCorrSmoothed[cut] = (TGraphErrors*)meanErrorsCorr[cut]->Clone();  // Meike
                  
                for (Int_t l=0; l<nPtBinsNew; l++){     // loop until nPtBinsNew
                    
                        if(nameCutVariationSC[cut].Contains("Alpha")){                    // Alpha

                            cout << "Smoothing not yet implemented..." << endl;
                            
                            //////////////////////////////////////////////////////////////////////////////
                        } else if(nameCutVariationSC[cut].Contains("YieldExtraction")){   // Yield Extraction
                            if(meson.CompareTo("Pi0")==0){                                       // Pi0
                                if(cent.Contains("6080")){                                    // 60-80%
                                    usePol3N[cut] = kTRUE;
                                    newPoint = pol3N->Eval(ptBinsCentersNew[l]);
                                    
                                } else if(cent.Contains("0010")){                             // 0-10%
                                    usePol3N[cut] = kTRUE;
                                    newPoint = pol3N->Eval(ptBinsCentersNew[l]);
                                    
                                } else if(cent.Contains("1020")){                             // 10-20%
                                    usePol3N[cut] = kTRUE;
                                    newPoint = pol3N->Eval(ptBinsCentersNew[l]);
                                    
                                } else if(cent.Contains("2040")){                             // 20-40%
                                    usePol3N[cut] = kTRUE;
                                    newPoint = pol3N->Eval(ptBinsCentersNew[l]);
                                    
                                } else if(cent.Contains("4060")){                             // 40-60%
                                    usePol2N[cut] = kTRUE;
                                    newPoint = pol2N->Eval(ptBinsCentersNew[l]);
                                    
                                }

                            } else if(meson.CompareTo("EtaToPi0")==0){
                                if(cent.Contains("6080")){                                    // 60-80%
                                    usePol3N[cut] = kTRUE;
                                    newPoint = pol3N->Eval(ptBinsCentersNew[l]);
                                    
                                } else if(cent.Contains("0010")){                             // 0-10%
                                    usePol2N[cut] = kTRUE;
                                    newPoint = pol2N->Eval(ptBinsCentersNew[l]);
                                    
                                } else if(cent.Contains("1020")){                             // 10-20%
                                    usePol3N[cut] = kTRUE;
                                    newPoint = pol3N->Eval(ptBinsCentersNew[l]);
                                    if(newPoint<0)  newPoint = pol3N->Eval(ptBinsCentersNew[5]);
                                    
                                } else if(cent.Contains("2040")){                             // 20-40%
                                    usePol2N[cut] = kTRUE;
                                    newPoint = pol2N->Eval(ptBinsCentersNew[l]);
                                    
                                } else if(cent.Contains("4060")){                             // 40-60%
                                    usePol2N[cut] = kTRUE;
                                    newPoint = pol2N->Eval(ptBinsCentersNew[l]);
                                    
                                }
                                                                
                            } else if(meson.CompareTo("Eta")==0){                               // Eta

                                if(cent.Contains("0010")){                                   // 0-10%
                                    usePol3N[cut] = kTRUE;
                                    newPoint = pol3N->Eval(ptBinsCentersNew[l]);
                                    if(newPoint<0)  newPoint = pol3N->Eval(ptBinsCentersNew[5]);
                                    
                                } else if(cent.Contains("1020")){                             // 10-20%
                                    usePol3[cut] = kTRUE;
                                    newPoint = pol3->Eval(ptBinsCentersNew[l]);
                                    if(newPoint<0)  newPoint = pol3N->Eval(ptBinsCentersNew[5]);
                                    
                                } else if(cent.Contains("2040")){                             // 20-40%
                                    usePol2N[cut] = kTRUE;
                                    newPoint = pol2N->Eval(ptBinsCentersNew[l]);
                                    if(newPoint<0)  newPoint = pol3N->Eval(ptBinsCentersNew[5]);
                                    
                                } else if(cent.Contains("4060")){                             // 40-60%
                                    usePol3N[cut] = kTRUE;
                                    newPoint = pol3N->Eval(ptBinsCentersNew[l]);
                                    if(newPoint<0)  newPoint = pol3N->Eval(ptBinsCentersNew[6]);
                                    
                                } else if(cent.Contains("6080")){                             // 60-80%
                                    usePol2[cut] = kTRUE;
                                    newPoint = pol2->Eval(ptBinsCentersNew[l]);
                                }
                            }
                            meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBinsCentersNew[l],newPoint);
                            
                            //////////////////////////////////////////////////////////////////////////////
                        } else if(nameCutVariationSC[cut].Contains("Chi2")){              // Chi2
                            if(meson.Contains("Pi0")){                                       // Pi0
			  
                                if(cent.Contains("6080")){                                    // 60-80 %
                                    usePol3N[cut] = kTRUE;
                                    newPoint = pol3N->Eval(ptBinsCentersNew[l]);
                                    
                                } else if(cent.Contains("0010")){                            // 0-10%
                                    usePol3N[cut] = kTRUE;			    
                                    newPoint = pol3N->Eval(ptBinsCentersNew[l]);
                                    
								}
                                
							} else if(meson.Contains("Eta")){                                 // Eta
								if(cent.Contains("0010")){                                      // 0-10%
                                    
                                     usePol2N[cut] = kTRUE;
                                     newPoint = pol2N->Eval(ptBinsCentersNew[l]);
                                     
                                } else if(cent.Contains("6080")){                                // 60-80%
                                     usePol3N[cut] = kTRUE;
                                     newPoint = pol3N->Eval(ptBinsCentersNew[l]);
                                     
                                }
							}
                            
                            meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBinsCentersNew[l],newPoint);
                            
                            //////////////////////////////////////////////////////////////////////////////
                        } else if(nameCutVariationSC[cut].Contains("dEdxE")){                   // dEdxE
                            if(meson.Contains("Pi0")){

                                usePol3N[cut] = kTRUE;
                                newPoint = pol3N->Eval(ptBinsCentersNew[l]);
                                
                            } else if(meson.Contains("Eta")){
                                if(cent.Contains("0010")){ 
                                    usePol2N[cut] = kTRUE;
                                    newPoint = pol2N->Eval(ptBinsCentersNew[l]);
                                    
                                } else if(cent.Contains("6080")){
                                    usePol1[cut] = kTRUE;
                                    newPoint = pol1->Eval(ptBinsCentersNew[l]);
                                }
                            }


                            meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBinsCentersNew[l],newPoint);
                            
                            //////////////////////////////////////////////////////////////////////////////
                        } else if(nameCutVariationSC[cut].Contains("dEdxPi")){
                            if(meson.Contains("Pi0")){

                                usePol1[cut] = kTRUE;
                                newPoint = pol1->Eval(ptBinsCentersNew[l]);
                                
                            } else if(meson.Contains("Eta")){
                                if(cent.Contains("0010")){
                                    usePol3N[cut] = kTRUE;
                                    newPoint = pol3N->Eval(ptBinsCentersNew[l]);
                                } else if(cent.Contains("6080")){
                                    usePol1[cut] = kTRUE;
                                    newPoint = pol1->Eval(ptBinsCentersNew[l]);
                                }
                            }
                            
                            meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBinsCentersNew[l],newPoint);
                            
                            //////////////////////////////////////////////////////////////////////////////
                        } else if(nameCutVariationSC[cut].Contains("Qt")){
                            if(meson.Contains("Pi0")){

                                if(cent.Contains("0010")){
                                    usePol1[cut] = kTRUE;
                                    newPoint = pol1->Eval(ptBinsCentersNew[l]);
                                    
                                } else if(cent.Contains("6080")){
                                    usePol3[cut] = kTRUE;
                                    newPoint = pol3->Eval(ptBinsCentersNew[l]);
                                
                                }
                                
                            } else if(meson.Contains("Eta")){

                                usePol4[cut] = kTRUE;
                                newPoint = pol4->Eval(ptBinsCentersNew[l]);

                            }

                            meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBinsCentersNew[l],newPoint);
                            
                            ////////////////////////////////////////////////////////////////////////////////////
                        } else if(nameCutVariationSC[cut].Contains("SinglePt")){                     // SinglePt
                            if(meson.Contains("Pi0")){                                               // Pi0
                                
                                if(cent.Contains("0010") ){                                          // 0-10%
                                    usePol3N[cut] = kTRUE;                                    
                                    newPoint = pol3N->Eval(ptBinsCentersNew[l]);
                                    
                                } else if(cent.Contains("6080")){                                    // 60-80%
                                    usePol4N[cut] = kTRUE;
                                    newPoint = pol4N->Eval(ptBinsCentersNew[l]);
                                                                        
                                } 
                            } else if(meson.Contains("Eta")){                                             // Eta

                                if(cent.Contains("0010") ){


                                    usePol4[cut] = kTRUE;   
                                    newPoint = pol4->Eval(ptBinsCentersNew[l]);
                                    
                                                                    
                                } else if(cent.Contains("6080")){

                                    usePol1[cut] = kTRUE;
                                    newPoint = pol1->Eval(ptBinsCentersNew[l]);

                                }   
                            }

                            meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBinsCentersNew[l],newPoint);
                            
                            //////////////////////////////////////////////////////////////////////////////
                        } else if(nameCutVariationSC[cut].Contains("TPCCluster")){                    // TPCCluster
                            if(meson.Contains("Pi0")){                                                // Pi0
                                
                                if(cent.Contains("0010")){                                            // 0-10%
                                    if(meson.CompareTo("Pi0")==0){
                                        usePol3N[cut] = kTRUE;
                                        if(ptBinsCentersNew[l] < 5){
                                            newPoint = pol3N->Eval(ptBinsCentersNew[l]);
                                        } else
                                            newPoint = 0;
                                    } else { // EtaToPi0
                                        usePol3N[cut] = kTRUE;
                                        if(ptBinsCentersNew[l] < 3){
                                            newPoint = pol3N->Eval(ptBinsCentersNew[l]);
                                        } else
                                            newPoint = 0;

                                    }
                                } else if(cent.Contains("6080")){                                     // 60-80%
                                    newPoint = 0;
                                }
                            } else if(meson.Contains("Eta")){                                         // Eta

                                if(cent.Contains("0010")){
                                    usePol3N[cut] = kTRUE;
                                    if(ptBinsCentersNew[l] < 3){
                                        newPoint = pol3N->Eval(ptBinsCentersNew[l]);
                                    } else
                                        newPoint = 0;

                                } else if(cent.Contains("6080")){                                     // 60-80%
                                    newPoint = 0;
                                }
                            }

                            meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBinsCentersNew[l],newPoint);
                        }

                        // in case fit function goes to negative values
                        if(newPoint<0) newPoint=0;
                        meanErrorsCorrSmoothed[cut]->SetPoint(l,ptBinsCentersNew[l],newPoint);
                        
                } //end of pt bins loop
	      
                cout << "done smoothing" << endl;
		  
                // *********** plot smoothing check ******************
		  
                histo2DCheckSmooth->DrawCopy();

                TLegend* legendSingleSyst;
                if (meson.Contains("Pi0")){legendSingleSyst= new TLegend(0.13,0.6,0.65,0.95);} 
                else {legendSingleSyst= new TLegend(0.20,0.7,0.7,0.95);}
                legendSingleSyst->SetTextSize(0.035);
                legendSingleSyst->SetFillColor(0);
                legendSingleSyst->SetFillStyle(0);	
                legendSingleSyst->SetBorderSize(0);
                legendSingleSyst->SetMargin(0.1);
                legendSingleSyst->AddEntry(meanErrorsCorr[cut], Form("%s",nameCutVariationSC[cut].Data()),"p");
                if (usePol1[cut])  legendSingleSyst->AddEntry(pol1, pol1Label,"l");
                if (usePol2[cut])  legendSingleSyst->AddEntry(pol2, pol2Label,"l");
                if (usePol3[cut])  legendSingleSyst->AddEntry(pol3, pol3Label,"l");
                if (usePol4[cut])  legendSingleSyst->AddEntry(pol4, pol4Label,"l");
                if (usePol5[cut])  legendSingleSyst->AddEntry(pol5, pol5Label,"l");
                if (usePol2N[cut]) legendSingleSyst->AddEntry(pol2N,pol2NLabel,"l");
                if (usePol3N[cut]) legendSingleSyst->AddEntry(pol3N,pol3NLabel,"l");
                if (usePol4N[cut]) legendSingleSyst->AddEntry(pol4N,pol4NLabel,"l");
                if (usePol6[cut]) legendSingleSyst->AddEntry(pol6, pol6Label,"l");
                if (usePol7[cut]) legendSingleSyst->AddEntry(pol7, pol7Label,"l");
                if (usePol8[cut]) legendSingleSyst->AddEntry(pol8, pol8Label,"l");
		  
                DrawGammaSetMarkerTGraphErr(meanErrorsCorr[cut], 20+cut, 1.,color[cut],color[cut]);
                meanErrorsCorr[cut]->Draw("p,csame");
                if (usePol1[cut])  pol1->Draw("same");
                if (usePol2[cut])  pol2->Draw("same");
                if (usePol3[cut])  pol3->Draw("same");
                if (usePol4[cut])  pol4->Draw("same");
                if (usePol5[cut])  pol5->Draw("same");
                if (usePol2N[cut]) pol2N->Draw("same");
                if (usePol3N[cut]) pol3N->Draw("same");
                if (usePol4N[cut]) pol4N->Draw("same");
                if (usePol6[cut]) pol6->Draw("same");
                if (usePol7[cut]) pol7->Draw("same");
                if (usePol8[cut]) pol8->Draw("same");

                DrawGammaSetMarkerTGraphErr(meanErrorsCorrSmoothed[cut], 24, 1.5,kGray+2,kGray+2);
                meanErrorsCorrSmoothed[cut]->Draw("p,csame");
                legendSingleSyst->AddEntry(meanErrorsCorrSmoothed[cut], Form("%s smoothed",nameCutVariationSC[cut].Data()),"p");
		  
                legendSingleSyst->Draw("same");
                canvasCheckSmooth->SaveAs(Form("%s/SysMeanNewWithMeanSingle_%s_%s_%s%s.%s",outputdir.Data(),nameCutVariationSC[cut].Data(),meson.Data(), energy.Data(),cent.Data(),suffix.Data()));
		  
            } 
        } //cuts loop
	}//end smoothing if

	//****************************************************
	// ***** values for material and mass resulution error
	//****************************************************
	Double_t errorMaterial;
	Double_t errorMassResolution;
	if (meson.CompareTo("EtaToPi0")==0){
		errorMaterial = 0.;
		errorMassResolution = 0.;
	} else {
		errorMaterial = 4.50;
		errorMassResolution = 1.4; // was 1.2;
	}
	//****************************************************

	// quadratically sum errors from all variations for each pT bin
	if(smooth==1){
        for(Int_t cut = 0; cut< nCuts; cut++){
            if(benable[cut]){
                errorsMeanCorrSmoothed[cut] = meanErrorsCorrSmoothed[cut]->GetY();
                for (Int_t l = 0; l < nPtBinsNew; l++){  
                    errorsMeanCorrSummedSmoothed[l] = errorsMeanCorrSummedSmoothed[l]+ pow(errorsMeanCorrSmoothed[cut][l],2);
                }
            } else {
                errorsMeanCorrSmoothed[cut] = NULL;      
            }
        }
	}

	// quadratically add material budget and mass resolution errors to the summed errors from all cut variations
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
		errorsPosCorrSummed[l]                = pow(errorsPosCorrSummed[l],0.5);
		errorsMeanCorrSummed[l]               = pow(errorsMeanCorrSummed[l],0.5);
		errorsPosErrCorrSummed[l]             = errorsPosCorrSummed[l]*0.001;
		errorsMeanErrCorrSummed[l]            = errorsMeanCorrSummed[l]*0.001;
		errorsMeanErrCorrMatSummed[l]         = errorsMeanCorrMatSummed[l]*0.001;
		errorsNegCorrSummed[l]                = -pow(errorsNegCorrSummed[l],0.5);
		errorsNegErrCorrSummed[l]             = errorsNegCorrSummed[l]*0.001;
	}

    if(smooth==1){  // Meike
        for (Int_t l = 0; l < nPtBinsNew; l++){
			errorsMeanCorrMatSummedSmoothed[l]    = pow(errorsMeanCorrSummedSmoothed[l]+ pow(2*errorMaterial ,2.),0.5);
			errorsMeanCorrSummedSmoothed[l]       = pow(errorsMeanCorrSummedSmoothed[l],0.5);
			errorsMeanErrCorrSummedSmoothed[l]    = errorsMeanCorrSummedSmoothed[l]*0.001;
            errorsMeanErrCorrMatSummedSmoothed[l] = errorsMeanCorrMatSummedSmoothed[l]*0.001;
        }
    }
    
	Double_t errorsMat[nPtBins];
	Double_t errorMassRes[nPtBins];
	for (Int_t l = 0; l < nPtBins; l++){
		errorsMat[l] = 2*errorMaterial;
		errorMassRes[l] = errorMassResolution;
		
	}

	// make graphs from arrays
	graphMaterialError                     = new TGraphErrors(nPtBins,ptBins ,errorsMat ,ptBinsErr ,errorsMeanErrSummed );
	graphMassResolError                    = new TGraphErrors(nPtBins,ptBins ,errorMassRes ,ptBinsErr ,errorsMeanErrSummed );
	negativeErrorsSummed 		           = new TGraphErrors(nPtBins,ptBins ,errorsNegSummed ,ptBinsErr ,errorsNegErrSummed );
	negativeErrorsCorrSummed 	           = new TGraphErrors(nPtBins,ptBins ,errorsNegCorrSummed ,ptBinsErr ,errorsNegErrCorrSummed );
	positiveErrorsSummed                   = new TGraphErrors(nPtBins,ptBins ,errorsPosSummed ,ptBinsErr ,errorsPosErrSummed );
	positiveErrorsCorrSummed 	           = new TGraphErrors(nPtBins,ptBins ,errorsPosCorrSummed ,ptBinsErr ,errorsPosErrCorrSummed );
	meanErrorsSummed 		               = new TGraphErrors(nPtBins,ptBins ,errorsMeanSummed ,ptBinsErr ,errorsMeanErrSummed );
	meanErrorsCorrSummed		           = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrSummed ,ptBinsErr ,errorsMeanErrCorrSummed );
	meanErrorsCorrSummedIncMat 	           = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrMatSummed ,ptBinsErr ,errorsMeanErrCorrMatSummed );
    
	if(smooth==1){
		meanErrorsCorrSummedSmoothed 	   = new TGraphErrors(nPtBinsNew,ptBinsCentersNew ,errorsMeanCorrSummedSmoothed ,ptBinsNewErr ,errorsMeanErrCorrSummedSmoothed );  // Meike
		meanErrorsCorrSummedIncMatSmoothed = new TGraphErrors(nPtBinsNew,ptBinsCentersNew ,errorsMeanCorrMatSummedSmoothed,ptBinsNewErr ,errorsMeanErrCorrMatSummedSmoothed );
	}
	
	
	/////////////////////////////////////////////////////////////////////////
	///////////////////////////////// Plotting //////////////////////////////
	/////////////////////////////////////////////////////////////////////////
	
	cout << "Plotting the systematic errors... " << endl;
	Double_t textSizeLabelsPixelmean        = 55;
	Double_t textSizeLabelsRelmean          = 55./1300;
	
	// canvas
	TCanvas* canvasSystErrors = new TCanvas("canvasSystErrors","",200,10,1300,1000);  // gives the page size
	DrawGammaCanvasSettings( canvasSystErrors, 0.075, 0.03, 0.015, 0.095);
	
	// dummy histo
    TH2D *histo2DSystErrors;
    if (meson.CompareTo("Pi0")==0)histo2DSystErrors = new TH2D("histo2DSystErrors", "histo2DSystErrors", 20,-0.1,maxPtPlot,1000.,0.,maxErrPlotIncMat);
    else histo2DSystErrors = new TH2D("histo2DSystErrors", "histo2DSystErrors", 20,0.,11,1000.,0.,maxErrPlotEtaIncMat);
    SetStyleHistoTH2ForGraphs( histo2DSystErrors, "#it{p}_{T} (GeV/#it{c})", "mean systematic Err %",
                               0.85*textSizeLabelsRelmean, textSizeLabelsRelmean, 0.85*textSizeLabelsRelmean, textSizeLabelsRelmean,0.9,0.75);
    histo2DSystErrors->DrawCopy();
    
	// *********** SysMean plot *******************************************
	
	// legend
    TLegend* legendMean; 
    if (meson.Contains("Pi0")) legendMean = new TLegend(0.11,0.65,0.5,0.95);
    else legendMean= new TLegend(0.11,0.65,0.55,0.95);
    legendMean->SetTextSize(0.035);
    legendMean->SetMargin(0.1);
    legendMean->SetFillColor(0);
    legendMean->SetFillStyle(0);
    legendMean->SetBorderSize(0);
    legendMean->SetNColumns(2);

	
	// cut variations systematic errors
    for(Int_t i = 0; i< nCuts ; i++){
        if(benable[i]){
            DrawGammaSetMarkerTGraphErr(meanErrors[i], markerStyle[i], 1.,color[i],color[i]);
            meanErrors[i]->Draw("pX0,csame");
            legendMean->AddEntry(meanErrors[i],nameCutVariation[i].Data(),"p");
        }
    }
	
	// summed error without material budget error
	DrawGammaSetMarkerTGraphErr(meanErrorsSummed, 20, 1.,1,1);
	meanErrorsSummed->Draw("pE0,csame");
	legendMean->AddEntry(meanErrorsSummed,"quad. sum.","p");

	
	// labels
    TLatex *labelThesis = new TLatex(0.6,0.83,"This thesis");
    SetStyleTLatex( labelThesis, 0.035,4);
    //labelThesis->Draw();        
    TLatex *labelMeson;
    if (meson.Contains("Pi0")) labelMeson= new TLatex(0.6,0.89,Form("#pi^{0} #rightarrow #gamma#gamma"));
    else labelMeson= new TLatex(0.6,0.89,Form("#eta #rightarrow #gamma_{conv}#gamma_{conv}"));
    SetStyleTLatex( labelMeson, 0.035,4);
    labelMeson->Draw();
    TLatex *labelCentrality = new TLatex(0.6,0.93,Form("%s %s",centPercent.Data(),collisionSystem.Data()));
    SetStyleTLatex( labelCentrality, 0.035,4);
    labelCentrality->Draw();
	
	// save canvas
    legendMean->Draw();
	canvasSystErrors->SaveAs(Form("%s/SysMean_%s_%s%s.%s",outputdir.Data(),meson.Data(), energy.Data(),cent.Data(),suffix.Data()));

	// *********** SysMeanNew plot *******************************************
	
	canvasSystErrors->cd();
    histo2DSystErrors->DrawCopy();
	
	// legend
    TLegend* legendMeanNew;
    if (meson.Contains("Pi0")) legendMeanNew= new TLegend(0.11,0.65,0.5,0.95);
    else legendMeanNew= new TLegend(0.11,0.65,0.55,0.95);
    legendMeanNew->SetTextSize(0.035);
    legendMeanNew->SetMargin(0.1);
    legendMeanNew->SetFillColor(0);
    legendMeanNew->SetFillStyle(0);
    legendMeanNew->SetBorderSize(0);
    legendMeanNew->SetNColumns(2);
	
	// cut variations systematic errors
    for(Int_t i = 0; i< nCuts ; i++){
        if(benable[i]){
            DrawGammaSetMarkerTGraphErr(meanErrorsCorr[i], markerStyle[i], 1.,color[i],color[i]);
            meanErrorsCorr[i]->Draw("pX0,csame");
            legendMeanNew->AddEntry(meanErrorsCorr[i],nameCutVariation[i].Data(),"p");	
        }
    }

	// material budget and mass resolution error only if NOT EtaToPi0
    if(!(meson.CompareTo("EtaToPi0")==0)){
        // material
        DrawGammaSetMarkerTGraphErr(graphMaterialError, 20, 1., kMagenta+2,kMagenta+2);
        graphMaterialError->Draw("pX0,csame");
        legendMeanNew->AddEntry(graphMaterialError,"Material","p");
        // resolution
        DrawGammaSetMarkerTGraphErr(graphMassResolError, 24, 1.,kGray+1,kGray+1);
        graphMassResolError->Draw("pX0,csame");
        legendMeanNew->AddEntry(graphMassResolError,"Mass resolution","p");
    }

	// summed systematic error
    DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummedIncMat, 20, 1.,kBlack,kBlack);  
    meanErrorsCorrSummedIncMat->Draw("p,csame");
    legendMeanNew->AddEntry(meanErrorsCorrSummedIncMat,"quad. sum.","p");

	// labels
    //labelThesis->Draw();        
    labelMeson->Draw();
    labelCentrality->Draw();

	// save canvas
	canvasSystErrors->Update();
	legendMeanNew->Draw();
	canvasSystErrors->SaveAs(Form("%s/SysMeanNew_%s_%s%s.%s",outputdir.Data(),meson.Data(), energy.Data(),cent.Data(),suffix.Data()));

	// smoothing
	if(smooth==2){
        for(Int_t i = 0; i< nCuts ; i++){
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
              
                canvasSystErrors->SaveAs(Form("%s/%sSysBeforeSmoothing_%s_%s_%s.%s",outputsmooth.Data(), nameCutVariationSC[i].Data(), meson.Data(), energy.Data(),cent.Data(),suffix.Data()));
            }
        }
	}
	
	if(smooth==1){
	  
        canvasSystErrors->cd();

        histo2DSystErrors->Draw();
              
        for(Int_t i = 0; i< nCuts ; i++){
            if(benable[i]){
                DrawGammaSetMarkerTGraphErr(meanErrorsCorrSmoothed[i], markerStyle[i], 1.,color[i],color[i]);
                meanErrorsCorrSmoothed[i]->Draw("pX0,csame");
                //legendMeanNew->AddEntry(meanErrorsCorr[i],nameCutVariation[i].Data(),"p");	
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
        //labelThesis->Draw();
	  
        canvasSystErrors->Update();
        canvasSystErrors->SaveAs(Form("%s/SysMeanNewWithMeanSmoothed_%s_%s%s.%s",outputsmooth.Data(),meson.Data(), energy.Data(),cent.Data(),suffix.Data()));
	  
	}	
    //

	/*
      Double_t errorsMeanCorrPID[nPtBins];	
      Double_t errorsMeanCorrSignalExtraction[nPtBins];	
      Double_t errorsMeanCorrTrackReco[nPtBins];	
      Double_t errorsMeanCorrPhotonReco[nPtBins];	
      Double_t errorsMeanCorrAcceptance[nPtBins];			
      Double_t errorsMeanCorrOther[nPtBins];	
      Double_t errorsMeanCorrPileUp[nPtBins];

      // sum up different categories of error sources
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
		
      // make graphs from arrays
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
      canvasSystErrors->SaveAs(Form("%s/SysErrorSummedVisu_%s_%s%s.%s",outputdir.Data(),meson.Data(), energy.Data(),cent.Data(),suffix.Data()));	

      // print to stdout
      if(verbose>0){
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
      }
	*/

	// cleanup
	delete canvasSystErrors;
	
	// save results in .dat files
	cout << "write results" << endl;

	TString energyForOutput                 = energy;
	energyForOutput.ReplaceAll(".","_");

    
	const char *SysErrDatnameMean = Form("%s_%s/SystematicErrorAveragedPCM_%s_%s%s_%s.dat",outputdirName.Data(), dateForOutput.Data(),meson.Data(),energyForOutput.Data(),cent.Data(), dateForOutput.Data());
	fstream SysErrDatAver;
	cout << SysErrDatnameMean << endl;
	SysErrDatAver.open(SysErrDatnameMean, ios::out);
	
    if(smooth==1){
        for (Int_t l=0; l<nPtBinsNew; l++){
            SysErrDatAver << ptBinsCentersNew[l] <<"\t" << "-"<< std::setprecision(4) <<errorsMeanCorrMatSummedSmoothed[l] << "\t" << errorsMeanCorrMatSummedSmoothed[l] << "\t"  << "-" << errorsMeanCorrSummedSmoothed[l] << "\t"  << errorsMeanCorrSummedSmoothed[l]  << endl;
        }
    } else {
        for (Int_t l=0; l< nPtBins; l++){
            SysErrDatAver  << ptBins[l] <<"\t"<< "-"<< std::setprecision(4) <<errorsMeanCorrMatSummed[l] << "\t"  << errorsMeanCorrMatSummed[l] << "\t"  << "-" << errorsMeanCorrSummed[l] << "\t"  << errorsMeanCorrSummed[l]  << endl;
        }
    }
   
	SysErrDatAver.close();


	const char *SysErrDatnameMeanSingleErr = Form("%s_%s/SystematicErrorAveragedSinglePCM_%s_%s%s_%s.dat", outputdirName.Data(), dateForOutput.Data(), meson.Data(), energyForOutput.Data(), cent.Data(), dateForOutput.Data());
	fstream SysErrDatAverSingle;
	cout << SysErrDatnameMeanSingleErr << endl;
	SysErrDatAverSingle.open(SysErrDatnameMeanSingleErr, ios::out);
	SysErrDatAverSingle << "Pt bin\t" ;
	for (Int_t i= 0; i< nCuts; i++){
        if(i==10) SysErrDatAverSingle << "YEPi0\t";
        else SysErrDatAverSingle << nameCutVariationSC[i] << "\t";
	}
	SysErrDatAverSingle << "Material\t";
    SysErrDatAverSingle << "Resolution\t";
    SysErrDatAverSingle << "Summed" << endl;
    
    if(smooth==1){
        for (Int_t l=0;l< nPtBinsNew;l++){ 
            SysErrDatAverSingle << ptBinsCentersNew[l] << "\t";
            for (Int_t i= 0; i< nCuts; i++){
                if(benable[i]){
                    SysErrDatAverSingle <<  std::setprecision(4) << errorsMeanCorrSmoothed[i][l] << "\t";      
                } else{
                    SysErrDatAverSingle <<  0 << "\t";
                }
            }
            SysErrDatAverSingle << 2*errorMaterial << "\t";
            SysErrDatAverSingle << errorMassResolution << "\t" ;
            SysErrDatAverSingle << errorsMeanCorrMatSummedSmoothed[l] << endl; 
        }
    } else{
        for (Int_t l=0;l< nPtBins;l++){ 
            SysErrDatAverSingle << ptBins[l] << "\t";
            for (Int_t i= 0; i< nCuts; i++){
                if(benable[i]){
                    SysErrDatAverSingle <<  std::setprecision(4) << errorsMeanCorr[i][l] << "\t";
                } else{
                    SysErrDatAverSingle <<  0 << "\t";
                }
            }
            SysErrDatAverSingle << 2*errorMaterial << "\t";
            SysErrDatAverSingle << errorMassResolution << "\t";
            SysErrDatAverSingle << errorsMeanCorrMatSummed[l] << endl;
        }
    }

	SysErrDatAverSingle.close();   
}


