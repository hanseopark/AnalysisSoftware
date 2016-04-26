/*****************************************************************************
******         provided by Gamma Conversion Group, PWG4,                ******
******        Ana Marin, marin@physi.uni-heidelberg.de                  ******
******           Kathrin Koch, kkoch@physi.uni-heidelberg.de            ******
******        Friederike Bock, friederike.bock@cern.ch                  ******
******        Lucia Leardini, lucia.leardini@cern.ch                    ******
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

extern TRandom*    gRandom;
extern TBenchmark*    gBenchmark;
extern TSystem*    gSystem;
extern TMinuit*      gMinuit;

Double_t    xSection900GeV =     47.78*1e-3;
Double_t    xSection2760GeV =    55.416*1e-3;
Double_t    xSection7000GeV =    62.22*1e-3;
Double_t    xSection8000GeV =    55.74*1e-3;
Double_t    recalcBarn =      1e12; //NLO in pbarn!!!!


TGraph* ScaleGraph (TGraph* graph, Double_t scaleFac){
	TGraph* dummyGraph = (TGraph*)graph->Clone(Form("%s_Scaled",graph->GetName()));
	Double_t * xValue = dummyGraph->GetX();
	Double_t * yValue = dummyGraph->GetY();

	Int_t nPoints = dummyGraph->GetN();
	for (Int_t i = 0; i < nPoints; i++){
		yValue[i] = yValue[i]*scaleFac;
	}
	TGraph* returnGraph = new TGraph(nPoints,xValue,yValue);
	return returnGraph;
}

TGraphAsymmErrors* ScaleGraphAsym (TGraphAsymmErrors* graph, Double_t scaleFac){
	TGraphAsymmErrors* dummyGraph = (TGraphAsymmErrors*)graph->Clone(Form("%s_Scaled",graph->GetName()));

	Double_t * xValue = dummyGraph->GetX(); 
	Double_t * yValue = dummyGraph->GetY();
	Double_t* xErrorLow = dummyGraph->GetEXlow();
	Double_t* xErrorHigh = dummyGraph->GetEXhigh();
	Double_t* yErrorLow = dummyGraph->GetEYlow();
	Double_t* yErrorHigh = dummyGraph->GetEYhigh();
	Int_t nPoints = dummyGraph->GetN();
	for (Int_t i = 0; i < nPoints; i++){
		yValue[i] = yValue[i]*scaleFac;
		yErrorLow[i] = yErrorLow[i]*scaleFac;
		yErrorHigh[i] = yErrorHigh[i]*scaleFac;
	}
	TGraphAsymmErrors* returnGraph =  new TGraphAsymmErrors(nPoints,xValue,yValue,xErrorLow,xErrorHigh,yErrorLow,yErrorHigh); 
	return returnGraph;
}

void ProduceTheoryGraphsDirectPhotons(){    
    
	StyleSettingsThesis();    
	SetPlotStyle();


    //******************************************************************************************************************
    //***************************** Direct Photon calculations *********************************************************
    //******************************************************************************************************************
    
    //******************************************************************************************************************
    //*********************************Old MC Gill *********************************************************************
    //******************************************************************************************************************
//     Double_t theoryMcGillPbPb0020_x[15]     = {    2.000000000000000111e-01, 4.000000000000000222e-01, 5.999999999999999778e-01, 8.000000000000000444e-01, 1.000000000000000000e+00, 
//                                                 1.199999999999999956e+00, 1.399999999999999911e+00, 1.600000000000000089e+00, 1.800000000000000044e+00, 2.000000000000000000e+00, 
//                                                 2.200000000000000178e+00, 2.399999999999999911e+00, 2.600000000000000089e+00, 2.799999999999999822e+00, 3.000000000000000000e+00};
//     Double_t theoryMcGillPbPb0020_xerr[15]     = {    0, 0, 0, 0, 0,
//                                                 0, 0, 0, 0, 0,
//                                                 0, 0, 0, 0, 0};
//     Double_t theoryMcGillPbPb0020_y[15]     = {    1.540859385329232794e+02, 2.044997387338719363e+01, 5.878499923344167044e+00, 2.250958012412696885e+00, 9.991863147603015083e-01,
//                                                 4.845094433818496471e-01, 2.497293734693046829e-01, 1.349175639294052931e-01, 7.576424621951352578e-02, 4.395437070531096890e-02,
//                                                 2.638121048976970612e-02, 1.632128463217997344e-02, 1.037788835831024437e-02, 6.782615199120981507e-03, 4.547945132616882172e-03 }; 
//     Double_t theoryMcGillPbPb0020_yerr[15]     = {    2.416453821976370708e+00, 3.301167073906846050e-01, 9.581685975141712719e-02, 3.707312930049495858e-02, 1.673623054763240595e-02, 
//                                                 8.318531597347338796e-03, 4.417758588599256936e-03, 2.462353352140851659e-03, 1.423226057051350064e-03, 8.454263715662204250e-04, 
//                                                 5.166073785903992467e-04, 3.233128079544251221e-04, 2.066024310156839574e-04, 1.349119265242332950e-04, 8.990109325466884435e-05}; 
//     TGraphErrors* graphTheoryMcGill0020     = new TGraphErrors(15,theoryMcGillPbPb0020_x,theoryMcGillPbPb0020_y,theoryMcGillPbPb0020_xerr, theoryMcGillPbPb0020_yerr);  
//     TGraphErrors* graphTheoryMcGill0020Plot    = ScaleGraph(graphTheoryMcGill0020,100);
//     graphTheoryMcGill0020Plot->RemovePoint(0);                                           
//     graphTheoryMcGill0020Plot->RemovePoint(0);
//     graphTheoryMcGill0020Plot->RemovePoint(0);
//     Double_t theoryMcGillPbPb2040_x[15]     = {    2.000000000000000111e-01, 4.000000000000000222e-01, 5.999999999999999778e-01, 8.000000000000000444e-01, 1.000000000000000000e+00, 
//                                                 1.199999999999999956e+00, 1.399999999999999911e+00, 1.600000000000000089e+00, 1.800000000000000044e+00, 2.000000000000000000e+00, 
//                                                 2.200000000000000178e+00, 2.399999999999999911e+00, 2.600000000000000089e+00, 2.799999999999999822e+00, 3.000000000000000000e+00}; 
//     Double_t theoryMcGillPbPb2040_xerr[15]     = {    0, 0, 0, 0, 0,
//                                                 0, 0, 0, 0, 0,
//                                                 0, 0, 0, 0, 0};
//     Double_t theoryMcGillPbPb2040_y[15]     = {    5.996162411828655081e+01, 7.682871584213885718e+00, 2.170336260707449672e+00, 8.182166469176795909e-01, 3.583236917664220922e-01,
//                                                 1.718390778093310256e-01, 8.781549983209124832e-02, 4.716549515653831182e-02, 2.639893678394553828e-02, 1.530106534174231064e-02,
//                                                 9.200643808420352898e-03, 5.715716866426879920e-03, 3.655916602971762304e-03, 2.407585539368926539e-03, 1.628794054128191544e-03 }; 
//     Double_t theoryMcGillPbPb2040_yerr[15]     = {    1.293765778931855293e+00, 1.700623864527457396e-01, 4.903841383176194696e-02, 1.908670016811140485e-02, 8.698233566001163999e-03,
//                                                 4.353428213951922483e-03, 2.319044963732616853e-03, 1.292520172768387397e-03, 7.456386913794923587e-04, 4.415830500157004066e-04, 
//                                                 2.693008616739346649e-04, 1.683700730118350920e-04, 1.075702856825843344e-04, 7.032156976682658054e-05, 4.696434116417641080e-05 }; 
//     TGraphErrors* graphTheoryMcGill2040     = new TGraphErrors(15,theoryMcGillPbPb2040_x,theoryMcGillPbPb2040_y,theoryMcGillPbPb2040_xerr, theoryMcGillPbPb2040_yerr);  
//     TGraphErrors* graphTheoryMcGill2040Plot    = ScaleGraph(graphTheoryMcGill2040,10);
//     graphTheoryMcGill2040Plot->RemovePoint(0);                                           
//     graphTheoryMcGill2040Plot->RemovePoint(0);
//     graphTheoryMcGill2040Plot->RemovePoint(0);
    
    
    //******************************************************************************************************************
    //*********************************MC Gill (Paquett, IPGlasma, 14.08.2015 - physi-mail) ****************************
    //******************************************************************************************************************    
    Double_t ptMCGill0020        [100];
    Double_t yieldMCGill0020    [100];
    Double_t errYieldMCGill0020    [100];
    Double_t errYieldXMCGill0020    [100];
    Int_t nlinesMCGill0020                     = 0;
    
    TString fileNameMCGill0020 = "ExternalInputPbPb/Theory/McGill/direct_photons_LHC2760_cent0020_Paquet_McGill.dat";
    ifstream  fileMCGill0020;
    fileMCGill0020.open(fileNameMCGill0020,ios_base::in);
    cout << fileNameMCGill0020 << endl;
    
    while(!fileMCGill0020.eof() && nlinesMCGill0020< 100){
        fileMCGill0020 >> ptMCGill0020[nlinesMCGill0020] >> yieldMCGill0020[nlinesMCGill0020] >> errYieldMCGill0020[nlinesMCGill0020] ; 
        cout << nlinesMCGill0020 << "\t"  << ptMCGill0020[nlinesMCGill0020] << "\t"  << yieldMCGill0020[nlinesMCGill0020] << "\t"  << errYieldMCGill0020[nlinesMCGill0020] << endl;;
        errYieldMCGill0020[nlinesMCGill0020]     = 0;
        errYieldXMCGill0020[nlinesMCGill0020]     = 0;
        nlinesMCGill0020++;
    }
    fileMCGill0020.close();
    TGraphErrors* graphDirectPhotonMCGill0020 = new TGraphErrors(nlinesMCGill0020-1,ptMCGill0020,yieldMCGill0020, errYieldXMCGill0020, errYieldMCGill0020); 

    Double_t ptMCGillPrompt0020        [100];
    Double_t yieldMCGillPrompt0020    [100];
    Double_t errYieldMCGillPrompt0020    [100];
    Double_t errYieldXMCGillPrompt0020    [100];
    Int_t nlinesMCGillPrompt0020                     = 0;
        
    TString fileNameMCGillPrompt0020 = "ExternalInputPbPb/Theory/McGill/prompt_photons_LHC2760_cent0020_Paquet_McGill.dat";
    ifstream  fileMCGillPrompt0020;
    fileMCGillPrompt0020.open(fileNameMCGillPrompt0020,ios_base::in);
    cout << fileNameMCGillPrompt0020 << endl;
    
    while(!fileMCGillPrompt0020.eof() && nlinesMCGillPrompt0020< 100){
        fileMCGillPrompt0020 >> ptMCGillPrompt0020[nlinesMCGillPrompt0020] >> yieldMCGillPrompt0020[nlinesMCGillPrompt0020] >> errYieldMCGillPrompt0020[nlinesMCGillPrompt0020] ; 
        cout << nlinesMCGillPrompt0020 << "\t"  << ptMCGillPrompt0020[nlinesMCGillPrompt0020] << "\t"  << yieldMCGillPrompt0020[nlinesMCGillPrompt0020] << "\t"  <<
        errYieldMCGillPrompt0020[nlinesMCGillPrompt0020] << endl;;
        errYieldMCGillPrompt0020[nlinesMCGillPrompt0020]    = 0;
        errYieldXMCGillPrompt0020[nlinesMCGillPrompt0020]     = 0;
        nlinesMCGillPrompt0020++;
    }
    fileMCGillPrompt0020.close();
    TGraphErrors* graphPromptPhotonMCGill0020 = new TGraphErrors(nlinesMCGillPrompt0020-1,ptMCGillPrompt0020,yieldMCGillPrompt0020, errYieldXMCGillPrompt0020, errYieldMCGillPrompt0020); 

    
    Double_t ptMCGill2040        [100];
    Double_t yieldMCGill2040    [100];
    Double_t errYieldMCGill2040    [100];
    Double_t errYieldXMCGill2040    [100];
    Int_t nlinesMCGill2040                     = 0;
    
    TString fileNameMCGill2040 = "ExternalInputPbPb/Theory/McGill/direct_photons_LHC2760_cent2040_Paquet_McGill.dat";
    ifstream  fileMCGill2040;
    fileMCGill2040.open(fileNameMCGill2040,ios_base::in);
    cout << fileNameMCGill2040 << endl;
    
    while(!fileMCGill2040.eof() && nlinesMCGill2040< 100){
        fileMCGill2040 >> ptMCGill2040[nlinesMCGill2040] >> yieldMCGill2040[nlinesMCGill2040] >> errYieldMCGill2040[nlinesMCGill2040] ; 
        cout << nlinesMCGill2040 << "\t"  << ptMCGill2040[nlinesMCGill2040] << "\t"  << yieldMCGill2040[nlinesMCGill2040] << "\t"  << errYieldMCGill2040[nlinesMCGill2040] << endl;;
        errYieldMCGill2040[nlinesMCGill2040]     = 0;
        errYieldXMCGill2040[nlinesMCGill2040]     = 0;
        nlinesMCGill2040++;
    }
    fileMCGill2040.close();
    TGraphErrors* graphDirectPhotonMCGill2040 = new TGraphErrors(nlinesMCGill2040-1,ptMCGill2040,yieldMCGill2040, errYieldXMCGill2040, errYieldMCGill2040); 

    Double_t ptMCGillPrompt2040        [100];
    Double_t yieldMCGillPrompt2040    [100];
    Double_t errYieldMCGillPrompt2040    [100];
    Double_t errYieldXMCGillPrompt2040    [100];
    Int_t nlinesMCGillPrompt2040                     = 0;
        
    TString fileNameMCGillPrompt2040 = "ExternalInputPbPb/Theory/McGill/prompt_photons_LHC2760_cent2040_Paquet_McGill.dat";
    ifstream  fileMCGillPrompt2040;
    fileMCGillPrompt2040.open(fileNameMCGillPrompt2040,ios_base::in);
    cout << fileNameMCGillPrompt2040 << endl;
    
    while(!fileMCGillPrompt2040.eof() && nlinesMCGillPrompt2040< 100){
        fileMCGillPrompt2040 >> ptMCGillPrompt2040[nlinesMCGillPrompt2040] >> yieldMCGillPrompt2040[nlinesMCGillPrompt2040] >> errYieldMCGillPrompt2040[nlinesMCGillPrompt2040] ; 
        cout << nlinesMCGillPrompt2040 << "\t"  << ptMCGillPrompt2040[nlinesMCGillPrompt2040] << "\t"  << yieldMCGillPrompt2040[nlinesMCGillPrompt2040] << "\t"  <<
        errYieldMCGillPrompt2040[nlinesMCGillPrompt2040] << endl;;
        errYieldMCGillPrompt2040[nlinesMCGillPrompt2040]    = 0;
        errYieldXMCGillPrompt2040[nlinesMCGillPrompt2040]     = 0;
        nlinesMCGillPrompt2040++;
    }
    fileMCGillPrompt2040.close();
    TGraphErrors* graphPromptPhotonMCGill2040 = new TGraphErrors(nlinesMCGillPrompt2040-1,ptMCGillPrompt2040,yieldMCGillPrompt2040, errYieldXMCGillPrompt2040, errYieldMCGillPrompt2040); 

    
    //******************************************************************************************************************
    //*********************************PHSD arXiv:1504.05699 [nucl-th], O. Linnyk **************************************
    //******************************************************************************************************************
    Double_t ptPHSD            [100];
    Double_t yieldPHSD0020    [100];
    Double_t yieldPHSD0040    [100];
    Double_t yieldPHSD2040    [100];
    Double_t yieldPHSD4080    [100];
    Double_t errPHSD0020    [100];
    Double_t errPHSD0040    [100];
    Double_t errPHSD2040    [100];
    Double_t errPHSD4080    [100];
    Double_t errXPHSD        [100];
    
    Int_t nlinesPHSD                     = 0;
    
    TString fileNamePHSD     = "ExternalInputPbPb/Theory/PHSD/Direct_Spcetra_Centralities_prepForRead.dat";
    ifstream  filePHSD;
    filePHSD.open(fileNamePHSD,ios_base::in);
    cout << fileNamePHSD << endl;
    
    while(!filePHSD.eof() && nlinesPHSD< 100){
        filePHSD >> ptPHSD[nlinesPHSD] >> yieldPHSD0020[nlinesPHSD] >> yieldPHSD2040[nlinesPHSD] >> yieldPHSD0040[nlinesPHSD]>> yieldPHSD4080[nlinesPHSD]; 
        cout << nlinesPHSD << "\t"  << ptPHSD[nlinesPHSD] << "\t"  << yieldPHSD0020[nlinesPHSD]<< "\t"  << yieldPHSD2040[nlinesPHSD] << "\t"  << yieldPHSD4080[nlinesPHSD]<< endl;;
        errPHSD0020[nlinesPHSD] = 0;//yieldPHSD0020[nlinesPHSD]*0.3;
        errPHSD0040[nlinesPHSD] = 0;//yieldPHSD0040[nlinesPHSD]*0.3;
        errPHSD2040[nlinesPHSD] = 0;//yieldPHSD2040[nlinesPHSD]*0.3;
        errPHSD4080[nlinesPHSD] = 0;//yieldPHSD4080[nlinesPHSD]*0.3;
        errXPHSD[nlinesPHSD]    = 0;
        nlinesPHSD++;
    }
    filePHSD.close();
    TGraphErrors* graphDirectPhotonPHSD0020 = new TGraphErrors(nlinesPHSD-1,ptPHSD,yieldPHSD0020, errXPHSD, errPHSD0020); 
    TGraphErrors* graphDirectPhotonPHSD2040 = new TGraphErrors(nlinesPHSD-1,ptPHSD,yieldPHSD2040, errXPHSD, errPHSD2040); 
    TGraphErrors* graphDirectPhotonPHSD0040 = new TGraphErrors(nlinesPHSD-1,ptPHSD,yieldPHSD0040, errXPHSD, errPHSD0040); 
    TGraphErrors* graphDirectPhotonPHSD4080 = new TGraphErrors(nlinesPHSD-1,ptPHSD,yieldPHSD4080, errXPHSD, errPHSD4080); 

    //******************************************************************************************************************
    //*********************************Fireball He et al. PRC85(2012)044911 ********************************************
    //******************************************************************************************************************
    Double_t ptHe            [100];
    Double_t yieldHe0020    [100];
//     Double_t yieldHe0040    [100];
    Double_t yieldHe2040    [100];
    Double_t yieldHe4080    [100];
    Double_t errHe0020    [100];
//     Double_t errHe0040    [100];
    Double_t errHe2040    [100];
    Double_t errHe4080    [100];
    Double_t errXHe        [100];
    
    Int_t nlinesHe                     = 0;
    
    TString fileNameHe     = "ExternalInputPbPb/Theory/Rapp/YieldFireball.txt";
    ifstream  fileHe;
    fileHe.open(fileNameHe,ios_base::in);
    cout << fileNameHe << endl;
    
    while(!fileHe.eof() && nlinesHe< 100){
        fileHe >> ptHe[nlinesHe] >> yieldHe0020[nlinesHe] >> yieldHe2040[nlinesHe] >> yieldHe4080[nlinesHe]; 
        cout << nlinesHe << "\t"  << ptHe[nlinesHe] << "\t"  << yieldHe0020[nlinesHe]<< "\t"  << yieldHe2040[nlinesHe] << "\t"  << yieldHe4080[nlinesHe]<< endl;;
        errHe0020[nlinesHe] = 0;//yieldHe0020[nlinesHe]*0.3;
//         errHe0040[nlinesHe] = 0;//yieldHe0040[nlinesHe]*0.3;
        errHe2040[nlinesHe] = 0;//yieldHe2040[nlinesHe]*0.3;
        errHe4080[nlinesHe] = 0;//yieldHe4080[nlinesHe]*0.3;
        errXHe[nlinesHe]    = 0;
        nlinesHe++;
    }
    fileHe.close();
    TGraphErrors* graphDirectPhotonHe0020 = new TGraphErrors(nlinesHe-1,ptHe,yieldHe0020, errXHe, errHe0020); 
    while (graphDirectPhotonHe0020->GetX()[0] < 0.7) graphDirectPhotonHe0020->RemovePoint(0);

    TGraphErrors* graphDirectPhotonHe2040 = new TGraphErrors(nlinesHe-1,ptHe,yieldHe2040, errXHe, errHe2040); 
    while (graphDirectPhotonHe2040->GetX()[0] < 0.7) graphDirectPhotonHe2040->RemovePoint(0);

//     TGraphErrors* graphDirectPhotonHe0040 = new TGraphErrors(nlinesHe-1,ptHe,yieldHe0040, errXHe, errHe0040); 
    TGraphErrors* graphDirectPhotonHe4080 = new TGraphErrors(nlinesHe-1,ptHe,yieldHe4080, errXHe, errHe4080); 
    while (graphDirectPhotonHe4080->GetX()[0] < 0.7) graphDirectPhotonHe4080->RemovePoint(0);
    
    //******************************************************************************************************************
    //*********************************Chatterjee http://journals.aps.org/prc/pdf/10.1103/PhysRevC.85.064910 ***********
    //******************************************************************************************************************
    Double_t ptChatterjee0020        [100];
    Double_t yieldThermalChatterjee0020    [100];
    Double_t yieldPromptChatterjee0020    [100];
    Double_t yieldChatterjee0020    [100];
    Double_t errYieldChatterjee0020    [100];
    Double_t errYieldXChatterjee0020    [100];
    Int_t nlinesChatterjee0020                     = 0;
    
    TString fileNameChatterjee0020 = "ExternalInputPbPb/Theory/Chatterjee/direct_0_20_forReading.dat";
    ifstream  fileChatterjee0020;
    fileChatterjee0020.open(fileNameChatterjee0020,ios_base::in);
    cout << fileNameChatterjee0020 << endl;
    
    while(!fileChatterjee0020.eof() && nlinesChatterjee0020< 100){
        fileChatterjee0020 >> ptChatterjee0020[nlinesChatterjee0020] >> yieldThermalChatterjee0020[nlinesChatterjee0020] >> yieldPromptChatterjee0020[nlinesChatterjee0020] 
                           >> yieldChatterjee0020[nlinesChatterjee0020] ; 
        cout << nlinesChatterjee0020 << "\t"  << ptChatterjee0020[nlinesChatterjee0020] << "\t"  << yieldThermalChatterjee0020[nlinesChatterjee0020] << "\t"  
             << yieldPromptChatterjee0020[nlinesChatterjee0020] << "\t"  << yieldChatterjee0020[nlinesChatterjee0020] << endl;;
        errYieldChatterjee0020[nlinesChatterjee0020]     = 0;
        errYieldXChatterjee0020[nlinesChatterjee0020]     = 0;
        nlinesChatterjee0020++;
    }
    fileChatterjee0020.close();
    TGraphErrors* graphDirectPhotonChatterjee0020 = new TGraphErrors(nlinesChatterjee0020-1,ptChatterjee0020,yieldChatterjee0020, errYieldXChatterjee0020, errYieldChatterjee0020); 
    
    TGraphErrors* graphDirectPhotonThermalChatterjee0020 = new TGraphErrors(nlinesChatterjee0020-1,ptChatterjee0020,yieldThermalChatterjee0020, errYieldXChatterjee0020, errYieldChatterjee0020); 
    while (graphDirectPhotonThermalChatterjee0020->GetY()[graphDirectPhotonThermalChatterjee0020->GetN()-1] == 0) 
        graphDirectPhotonThermalChatterjee0020->RemovePoint(graphDirectPhotonThermalChatterjee0020->GetN()-1); 
    
    TGraphErrors* graphDirectPhotonPromptChatterjee0020 = new TGraphErrors(nlinesChatterjee0020-1,ptChatterjee0020,yieldPromptChatterjee0020, errYieldXChatterjee0020, errYieldChatterjee0020); 
    while (graphDirectPhotonPromptChatterjee0020->GetY()[0] == 0) 
        graphDirectPhotonPromptChatterjee0020->RemovePoint(0);

    Double_t ptChatterjee2040        [100];
    Double_t yieldThermalChatterjee2040    [100];
    Double_t yieldPromptChatterjee2040    [100];
    Double_t yieldChatterjee2040    [100];
    Double_t errYieldChatterjee2040    [100];
    Double_t errYieldXChatterjee2040    [100];
    Int_t nlinesChatterjee2040                     = 0;
    
    TString fileNameChatterjee2040 = "ExternalInputPbPb/Theory/Chatterjee/direct_20_40_forReading.dat";
    ifstream  fileChatterjee2040;
    fileChatterjee2040.open(fileNameChatterjee2040,ios_base::in);
    cout << fileNameChatterjee2040 << endl;
    
    while(!fileChatterjee2040.eof() && nlinesChatterjee2040< 100){
        fileChatterjee2040 >> ptChatterjee2040[nlinesChatterjee2040] >> yieldThermalChatterjee2040[nlinesChatterjee2040] >> yieldPromptChatterjee2040[nlinesChatterjee2040] 
                           >> yieldChatterjee2040[nlinesChatterjee2040] ; 
        cout << nlinesChatterjee2040 << "\t"  << ptChatterjee2040[nlinesChatterjee2040] << "\t"  << yieldThermalChatterjee2040[nlinesChatterjee2040] << "\t"  
             << yieldPromptChatterjee2040[nlinesChatterjee2040] << "\t"  << yieldChatterjee2040[nlinesChatterjee2040] << endl;;
        errYieldChatterjee2040[nlinesChatterjee2040]     = 0;
        errYieldXChatterjee2040[nlinesChatterjee2040]     = 0;
        nlinesChatterjee2040++;
    }
    fileChatterjee2040.close();
    TGraphErrors* graphDirectPhotonChatterjee2040 = new TGraphErrors(nlinesChatterjee2040-1,ptChatterjee2040,yieldChatterjee2040, errYieldXChatterjee2040, errYieldChatterjee2040); 
    
    TGraphErrors* graphDirectPhotonThermalChatterjee2040 = new TGraphErrors(nlinesChatterjee2040-1,ptChatterjee2040,yieldThermalChatterjee2040, errYieldXChatterjee2040, errYieldChatterjee2040); 
    while (graphDirectPhotonThermalChatterjee2040->GetY()[graphDirectPhotonThermalChatterjee2040->GetN()-1] == 0) 
        graphDirectPhotonThermalChatterjee2040->RemovePoint(graphDirectPhotonThermalChatterjee2040->GetN()-1); 
    
    TGraphErrors* graphDirectPhotonPromptChatterjee2040 = new TGraphErrors(nlinesChatterjee2040-1,ptChatterjee2040,yieldPromptChatterjee2040, errYieldXChatterjee2040, errYieldChatterjee2040); 
    while (graphDirectPhotonPromptChatterjee2040->GetY()[0] == 0) 
        graphDirectPhotonPromptChatterjee2040->RemovePoint(0);

    Double_t ptChatterjee4060        [100];
    Double_t yieldThermalChatterjee4060    [100];
    Double_t yieldPromptChatterjee4060    [100];
    Double_t yieldChatterjee4060    [100];
    Double_t errYieldChatterjee4060    [100];
    Double_t errYieldXChatterjee4060    [100];
    Int_t nlinesChatterjee4060                     = 0;
    
    TString fileNameChatterjee4060 = "ExternalInputPbPb/Theory/Chatterjee/direct_40_60_forReading.dat";
    ifstream  fileChatterjee4060;
    fileChatterjee4060.open(fileNameChatterjee4060,ios_base::in);
    cout << fileNameChatterjee4060 << endl;
    
    while(!fileChatterjee4060.eof() && nlinesChatterjee4060< 100){
        fileChatterjee4060 >> ptChatterjee4060[nlinesChatterjee4060] >> yieldThermalChatterjee4060[nlinesChatterjee4060] >> yieldPromptChatterjee4060[nlinesChatterjee4060] 
                           >> yieldChatterjee4060[nlinesChatterjee4060] ; 
        cout << nlinesChatterjee4060 << "\t"  << ptChatterjee4060[nlinesChatterjee4060] << "\t"  << yieldThermalChatterjee4060[nlinesChatterjee4060] << "\t"  
             << yieldPromptChatterjee4060[nlinesChatterjee4060] << "\t"  << yieldChatterjee4060[nlinesChatterjee4060] << endl;;
        errYieldChatterjee4060[nlinesChatterjee4060]     = 0;
        errYieldXChatterjee4060[nlinesChatterjee4060]     = 0;
        nlinesChatterjee4060++;
    }
    fileChatterjee4060.close();
    TGraphErrors* graphDirectPhotonChatterjee4060 = new TGraphErrors(nlinesChatterjee4060-1,ptChatterjee4060,yieldChatterjee4060, errYieldXChatterjee4060, errYieldChatterjee4060); 
    
    TGraphErrors* graphDirectPhotonThermalChatterjee4060 = new TGraphErrors(nlinesChatterjee4060-1,ptChatterjee4060,yieldThermalChatterjee4060, errYieldXChatterjee4060, errYieldChatterjee4060); 
    while (graphDirectPhotonThermalChatterjee4060->GetY()[graphDirectPhotonThermalChatterjee4060->GetN()-1] == 0) 
        graphDirectPhotonThermalChatterjee4060->RemovePoint(graphDirectPhotonThermalChatterjee4060->GetN()-1); 
    
    TGraphErrors* graphDirectPhotonPromptChatterjee4060 = new TGraphErrors(nlinesChatterjee4060-1,ptChatterjee4060,yieldPromptChatterjee4060, errYieldXChatterjee4060, errYieldChatterjee4060); 
    while (graphDirectPhotonPromptChatterjee4060->GetY()[0] == 0) 
        graphDirectPhotonPromptChatterjee4060->RemovePoint(0);

    //******************************************************************************************************************
    //*********************************Chatterjee thermal from Rupa, pQCD from JHEP 1305 (2013) 030 -0-20% *************
    //******************************************************************************************************************
    Double_t ptChatterjee0020_2[6]                = { 1.3, 2, 3, 4, 5,
                                                        6};
    Double_t yieldThermalChatterjee0020_2[6]      = { 0.231531, 0.0343742, 0.00357161, 0.000518726, 9.30806e-05,
                                                        1.91335e-05};
    Double_t yieldPromptChatterjee0020_2[6]       = { 0.0831538, 0.0160724, 0.00253054, 0.00063178, 0.000216407,
                                                        8.72282e-05};
    Double_t yieldChatterjee0020_2[6]             = { 0.314685, 0.0504466, 0.00610215, 0.00115051, 0.000309488,
                                                        0.000106362};
    Double_t errYieldXChatterjee0020_2[6]         = { 0, 0, 0, 0, 0,
                                                        0  };
    Double_t errYieldYChatterjee0020_2[6]         = { 0, 0, 0, 0, 0,
                                                        0  };
    Int_t nlinesChatterjee0020_2                  = 6;

    TGraphErrors* graphDirectPhotonChatterjee0020_2         = new TGraphErrors(nlinesChatterjee0020_2,ptChatterjee0020_2,yieldChatterjee0020_2, errYieldXChatterjee0020_2, errYieldYChatterjee0020_2); 
    TGraphErrors* graphDirectPhotonThermalChatterjee0020_2  = new TGraphErrors(nlinesChatterjee0020_2,ptChatterjee0020_2,yieldThermalChatterjee0020_2, errYieldXChatterjee0020_2, errYieldYChatterjee0020_2);   
    TGraphErrors* graphDirectPhotonPromptChatterjee0020_2   = new TGraphErrors(nlinesChatterjee0020_2,ptChatterjee0020_2,yieldPromptChatterjee0020_2, errYieldXChatterjee0020_2, errYieldYChatterjee0020_2); 

    //******************************************************************************************************************
    //*********************************Chatterjee thermal from Rupa, pQCD from JHEP 1305 (2013) 030 -20-40% ************
    //******************************************************************************************************************
    Double_t ptChatterjee2040_2[6]                = { 1.3, 2, 3, 4, 5,
                                                        6};
    Double_t yieldThermalChatterjee2040_2[6]      = { 0.0920386, 0.013742, 0.00140182, 0.000194282, 3.29015e-05,
                                                        6.38523e-06};
    Double_t yieldPromptChatterjee2040_2[6]       = { 0.0302371, 0.00574401, 0.000894192, 0.000222269, 7.58425e-05,
                                                        3.05298e-05};
    Double_t yieldChatterjee2040_2[6]             = { 0.122276, 0.019486, 0.00229601, 0.000416551, 0.000108744,
                                                        3.6915e-05};
    Double_t errYieldXChatterjee2040_2[6]         = { 0, 0, 0, 0, 0,
                                                        0  };
    Double_t errYieldYChatterjee2040_2[6]         = { 0, 0, 0, 0, 0,
                                                        0  };
    Int_t nlinesChatterjee2040_2                  = 6;                                               

    TGraphErrors* graphDirectPhotonChatterjee2040_2         = new TGraphErrors(nlinesChatterjee2040_2, ptChatterjee2040_2, yieldChatterjee2040_2, 
                                                                               errYieldXChatterjee2040_2, errYieldYChatterjee2040_2); 
    TGraphErrors* graphDirectPhotonThermalChatterjee2040_2  = new TGraphErrors(nlinesChatterjee2040_2, ptChatterjee2040_2, yieldThermalChatterjee2040_2, 
                                                                               errYieldXChatterjee2040_2, errYieldYChatterjee2040_2);   
    TGraphErrors* graphDirectPhotonPromptChatterjee2040_2   = new TGraphErrors(nlinesChatterjee2040_2, ptChatterjee2040_2, yieldPromptChatterjee2040_2, 
                                                                               errYieldXChatterjee2040_2, errYieldYChatterjee2040_2); 
    
    //******************************************************************************************************************
    //*********************************v.Hees, Rapp NPA933(2015)256 ****************************************************
    //******************************************************************************************************************

    Double_t ptHees0020        [100];
    Double_t yieldRhoSF0020    [100];
    Double_t yieldQGPHees0020    [100];
    Double_t yieldOmegaHees0020    [100];
    Double_t yieldMesonGasHees0020    [100];
    Double_t yieldPrimordialHees0020    [100];
    Double_t yieldHees0020    [100];
    Double_t errYieldHees0020    [100];
    Double_t errYieldXHees0020    [100];
    Int_t nlinesHees0020                     = 0;
    
    TString fileNameHees0020 = "ExternalInputPbPb/Theory/Rapp/Yield0020.txt";
    ifstream  fileHees0020;
    fileHees0020.open(fileNameHees0020,ios_base::in);
    cout << fileNameHees0020 << endl;
    
    while(!fileHees0020.eof() && nlinesHees0020< 100){
        fileHees0020 >> ptHees0020[nlinesHees0020] >> yieldRhoSF0020[nlinesHees0020] >>yieldQGPHees0020[nlinesHees0020] >>yieldOmegaHees0020[nlinesHees0020] >> yieldMesonGasHees0020[nlinesHees0020]
                     >> yieldPrimordialHees0020[nlinesHees0020] >> yieldHees0020[nlinesHees0020] ; 
        cout << nlinesHees0020 << "\t"  << ptHees0020[nlinesHees0020] << "\t"  << yieldQGPHees0020[nlinesHees0020] << "\t"  
             << yieldPrimordialHees0020[nlinesHees0020] << "\t"  << yieldHees0020[nlinesHees0020] << endl;;
        errYieldHees0020[nlinesHees0020]     = 0;
        errYieldXHees0020[nlinesHees0020]     = 0;
        nlinesHees0020++;
    }
    fileHees0020.close();
    TGraphErrors* graphDirectPhotonHees0020 = new TGraphErrors(nlinesHees0020-1,ptHees0020,yieldHees0020, errYieldXHees0020, errYieldHees0020); 
    while (graphDirectPhotonHees0020->GetX()[0] < 0.7) graphDirectPhotonHees0020->RemovePoint(0);
    TGraphErrors* graphDirectPhotonRhoSFHees0020 = new TGraphErrors(nlinesHees0020-1,ptHees0020,yieldRhoSF0020, errYieldXHees0020, errYieldHees0020);
    TGraphErrors* graphDirectPhotonQGPHees0020 = new TGraphErrors(nlinesHees0020-1,ptHees0020,yieldQGPHees0020, errYieldXHees0020, errYieldHees0020);
    TGraphErrors* graphDirectPhotonOmegaHees0020 = new TGraphErrors(nlinesHees0020-1,ptHees0020,yieldOmegaHees0020, errYieldXHees0020, errYieldHees0020);
    TGraphErrors* graphDirectPhotonMesonGasHees0020 = new TGraphErrors(nlinesHees0020-1,ptHees0020,yieldMesonGasHees0020, errYieldXHees0020, errYieldHees0020);
    TGraphErrors* graphDirectPhotonPrimordialHees0020 = new TGraphErrors(nlinesHees0020-1,ptHees0020,yieldPrimordialHees0020, errYieldXHees0020, errYieldHees0020); 
    
    Double_t ptHees2040        [100];
    Double_t yieldRhoSF2040    [100];
    Double_t yieldQGPHees2040    [100];
    Double_t yieldOmegaHees2040    [100];
    Double_t yieldMesonGasHees2040    [100];
    Double_t yieldPrimordialHees2040    [100];
    Double_t yieldHees2040    [100];
    Double_t errYieldHees2040    [100];
    Double_t errYieldXHees2040    [100];
    Int_t nlinesHees2040                     = 0;
    
    TString fileNameHees2040 = "ExternalInputPbPb/Theory/Rapp/Yield2040.txt";
    ifstream  fileHees2040;
    fileHees2040.open(fileNameHees2040,ios_base::in);
    cout << fileNameHees2040 << endl;
    
    while(!fileHees2040.eof() && nlinesHees2040< 100){
        fileHees2040 >> ptHees2040[nlinesHees2040] >> yieldRhoSF2040[nlinesHees2040] >>yieldQGPHees2040[nlinesHees2040] >>yieldOmegaHees2040[nlinesHees2040] >> yieldMesonGasHees2040[nlinesHees2040]
                     >> yieldPrimordialHees2040[nlinesHees2040] >> yieldHees2040[nlinesHees2040] ; 
        cout << nlinesHees2040 << "\t"  << ptHees2040[nlinesHees2040] << "\t"  << yieldQGPHees2040[nlinesHees2040] << "\t"  
             << yieldPrimordialHees2040[nlinesHees2040] << "\t"  << yieldHees2040[nlinesHees2040] << endl;;
        errYieldHees2040[nlinesHees2040]     = 0;
        errYieldXHees2040[nlinesHees2040]     = 0;
        nlinesHees2040++;
    }
    fileHees2040.close();
    TGraphErrors* graphDirectPhotonHees2040 = new TGraphErrors(nlinesHees2040-1,ptHees2040,yieldHees2040, errYieldXHees2040, errYieldHees2040); 
    while (graphDirectPhotonHees2040->GetX()[0] < 0.7) graphDirectPhotonHees2040->RemovePoint(0);
    TGraphErrors* graphDirectPhotonRhoSFHees2040 = new TGraphErrors(nlinesHees2040-1,ptHees2040,yieldRhoSF2040, errYieldXHees2040, errYieldHees2040);
    TGraphErrors* graphDirectPhotonQGPHees2040 = new TGraphErrors(nlinesHees2040-1,ptHees2040,yieldQGPHees2040, errYieldXHees2040, errYieldHees2040);
    TGraphErrors* graphDirectPhotonOmegaHees2040 = new TGraphErrors(nlinesHees2040-1,ptHees2040,yieldOmegaHees2040, errYieldXHees2040, errYieldHees2040);
    TGraphErrors* graphDirectPhotonMesonGasHees2040 = new TGraphErrors(nlinesHees2040-1,ptHees2040,yieldMesonGasHees2040, errYieldXHees2040, errYieldHees2040);
    TGraphErrors* graphDirectPhotonPrimordialHees2040 = new TGraphErrors(nlinesHees2040-1,ptHees2040,yieldPrimordialHees2040, errYieldXHees2040, errYieldHees2040); 

    Double_t ptHees4080        [100];
    Double_t yieldRhoSF4080    [100];
    Double_t yieldQGPHees4080    [100];
    Double_t yieldOmegaHees4080    [100];
    Double_t yieldMesonGasHees4080    [100];
    Double_t yieldPrimordialHees4080    [100];
    Double_t yieldHees4080    [100];
    Double_t errYieldHees4080    [100];
    Double_t errYieldXHees4080    [100];
    Int_t nlinesHees4080                     = 0;
    
    TString fileNameHees4080 = "ExternalInputPbPb/Theory/Rapp/Yield4080.txt";
    ifstream  fileHees4080;
    fileHees4080.open(fileNameHees4080,ios_base::in);
    cout << fileNameHees4080 << endl;
    
    while(!fileHees4080.eof() && nlinesHees4080< 100){
        fileHees4080 >> ptHees4080[nlinesHees4080] >> yieldRhoSF4080[nlinesHees4080] >>yieldQGPHees4080[nlinesHees4080] >>yieldOmegaHees4080[nlinesHees4080] >> yieldMesonGasHees4080[nlinesHees4080]
                     >> yieldPrimordialHees4080[nlinesHees4080] >> yieldHees4080[nlinesHees4080] ; 
        cout << nlinesHees4080 << "\t"  << ptHees4080[nlinesHees4080] << "\t"  << yieldQGPHees4080[nlinesHees4080] << "\t"  
             << yieldPrimordialHees4080[nlinesHees4080] << "\t"  << yieldHees4080[nlinesHees4080] << endl;;
        errYieldHees4080[nlinesHees4080]     = 0;
        errYieldXHees4080[nlinesHees4080]     = 0;
        nlinesHees4080++;
    }
    fileHees4080.close();
    TGraphErrors* graphDirectPhotonHees4080 = new TGraphErrors(nlinesHees4080-1,ptHees4080,yieldHees4080, errYieldXHees4080, errYieldHees4080); 
    while (graphDirectPhotonHees4080->GetX()[0] < 0.7) graphDirectPhotonHees4080->RemovePoint(0);
    TGraphErrors* graphDirectPhotonRhoSFHees4080 = new TGraphErrors(nlinesHees4080-1,ptHees4080,yieldRhoSF4080, errYieldXHees4080, errYieldHees4080);
    TGraphErrors* graphDirectPhotonQGPHees4080 = new TGraphErrors(nlinesHees4080-1,ptHees4080,yieldQGPHees4080, errYieldXHees4080, errYieldHees4080);
    TGraphErrors* graphDirectPhotonOmegaHees4080 = new TGraphErrors(nlinesHees4080-1,ptHees4080,yieldOmegaHees4080, errYieldXHees4080, errYieldHees4080);
    TGraphErrors* graphDirectPhotonMesonGasHees4080 = new TGraphErrors(nlinesHees4080-1,ptHees4080,yieldMesonGasHees4080, errYieldXHees4080, errYieldHees4080);
    TGraphErrors* graphDirectPhotonPrimordialHees4080 = new TGraphErrors(nlinesHees4080-1,ptHees4080,yieldPrimordialHees4080, errYieldXHees4080, errYieldHees4080); 
    
    //******************************************************************************************************************
    //*********************************Holopainen **********************************************************************
    //******************************************************************************************************************
    
    Double_t ptHolopainen0020        [100];
    Double_t v2Holopainen0020    [100];
    Double_t yieldHolopainen0020    [100];
    Double_t errYieldHolopainen0020    [100];
    Double_t errYieldXHolopainen0020    [100];
    Int_t nlinesHolopainen0020                     = 0;
    
    TString fileNameHolopainen0020 = "ExternalInputPbPb/Theory/Holopainen/yieldAndv2_0020_toRead.txt";
    ifstream  fileHolopainen0020;
    fileHolopainen0020.open(fileNameHolopainen0020,ios_base::in);
    cout << fileNameHolopainen0020 << endl;
    
    while(!fileHolopainen0020.eof() && nlinesHolopainen0020< 100){
        fileHolopainen0020 >> ptHolopainen0020[nlinesHolopainen0020] >> yieldHolopainen0020[nlinesHolopainen0020] >> v2Holopainen0020[nlinesHolopainen0020] ; 
        cout << nlinesHolopainen0020 << "\t"  << ptHolopainen0020[nlinesHolopainen0020] << "\t"  << v2Holopainen0020[nlinesHolopainen0020] << "\t"  << yieldHolopainen0020[nlinesHolopainen0020] << endl;    
        errYieldHolopainen0020[nlinesHolopainen0020]     = 0;
        errYieldXHolopainen0020[nlinesHolopainen0020]     = 0;
        nlinesHolopainen0020++;
    }
    fileHolopainen0020.close();
    TGraphErrors* graphDirectPhotonHolopainen0020 = new TGraphErrors(nlinesHolopainen0020-1,ptHolopainen0020,yieldHolopainen0020, errYieldXHolopainen0020, errYieldHolopainen0020); 
    while (graphDirectPhotonHolopainen0020->GetX()[0] < 0.7) graphDirectPhotonHolopainen0020->RemovePoint(0);
    TGraphErrors* graphDirectPhotonV2Holopainen0020 = new TGraphErrors(nlinesHolopainen0020-1,ptHolopainen0020,v2Holopainen0020, errYieldXHolopainen0020, errYieldHolopainen0020); 
    while (graphDirectPhotonV2Holopainen0020->GetX()[0] < 0.7) graphDirectPhotonV2Holopainen0020->RemovePoint(0);
    
    Double_t ptHolopainen2040        [100];
    Double_t v2Holopainen2040    [100];
    Double_t yieldHolopainen2040    [100];
    Double_t errYieldHolopainen2040    [100];
    Double_t errYieldXHolopainen2040    [100];
    Int_t nlinesHolopainen2040                     = 0;
    
    TString fileNameHolopainen2040 = "ExternalInputPbPb/Theory/Holopainen/yieldAndv2_2040_toRead.txt";
    ifstream  fileHolopainen2040;
    fileHolopainen2040.open(fileNameHolopainen2040,ios_base::in);
    cout << fileNameHolopainen2040 << endl;
    
    while(!fileHolopainen2040.eof() && nlinesHolopainen2040< 100){
        fileHolopainen2040 >> ptHolopainen2040[nlinesHolopainen2040] >> yieldHolopainen2040[nlinesHolopainen2040] >> v2Holopainen2040[nlinesHolopainen2040] ; 
        cout << nlinesHolopainen2040 << "\t"  << ptHolopainen2040[nlinesHolopainen2040] << "\t"  << v2Holopainen2040[nlinesHolopainen2040] << "\t"  << yieldHolopainen2040[nlinesHolopainen2040] << endl;    
        errYieldHolopainen2040[nlinesHolopainen2040]     = 0;
        errYieldXHolopainen2040[nlinesHolopainen2040]     = 0;
        nlinesHolopainen2040++;
    }
    fileHolopainen2040.close();
    TGraphErrors* graphDirectPhotonHolopainen2040 = new TGraphErrors(nlinesHolopainen2040-1,ptHolopainen2040,yieldHolopainen2040, errYieldXHolopainen2040, errYieldHolopainen2040); 
    while (graphDirectPhotonHolopainen2040->GetX()[0] < 0.7) graphDirectPhotonHolopainen2040->RemovePoint(0);
    TGraphErrors* graphDirectPhotonV2Holopainen2040 = new TGraphErrors(nlinesHolopainen2040-1,ptHolopainen2040,v2Holopainen2040, errYieldXHolopainen2040, errYieldHolopainen2040); 
    while (graphDirectPhotonV2Holopainen2040->GetX()[0] < 0.7) graphDirectPhotonV2Holopainen2040->RemovePoint(0);
    
    
    
	TFile *fileTheoryGraphsPbPb = new TFile("ExternalInputPbPb/Theory/TheoryCompilationPbPb.root","RECREATE");

	
        fileTheoryGraphsPbPb->mkdir("DirectPhoton");
        TDirectoryFile* directoryGamma = (TDirectoryFile*)fileTheoryGraphsPbPb->Get("DirectPhoton"); 
        fileTheoryGraphsPbPb->cd("DirectPhoton");

        graphDirectPhotonMCGill0020->Write("graphDirectPhotonYield_McGill_0020");
        graphDirectPhotonMCGill2040->Write("graphDirectPhotonYield_McGill_2040");
        graphPromptPhotonMCGill0020->Write("graphPromptPhotonYield_McGill_0020");
        graphPromptPhotonMCGill2040->Write("graphPromptPhotonYield_McGill_2040");
        
        graphDirectPhotonPHSD0020->Write("graphDirectPhotonYield_PHSD_0020");
        graphDirectPhotonPHSD0040->Write("graphDirectPhotonYield_PHSD_0040");
        graphDirectPhotonPHSD2040->Write("graphDirectPhotonYield_PHSD_2040");
        graphDirectPhotonPHSD4080->Write("graphDirectPhotonYield_PHSD_4080");

        graphDirectPhotonHe0020->Write("graphDirectPhotonYield_He_0020");
//         graphDirectPhotonHe0040->Write("graphDirectPhotonYield_He_0040");
        graphDirectPhotonHe2040->Write("graphDirectPhotonYield_He_2040");
        graphDirectPhotonHe4080->Write("graphDirectPhotonYield_He_4080");
        
        graphDirectPhotonChatterjee0020->Write("graphDirectPhotonYield_Chatterjee_0020");
        graphDirectPhotonThermalChatterjee0020->Write("graphDirectPhotonThermalYield_Chatterjee_0020");
        graphDirectPhotonPromptChatterjee0020->Write("graphDirectPhotonPromptYield_Chatterjee_0020");
        graphDirectPhotonChatterjee2040->Write("graphDirectPhotonYield_Chatterjee_2040");
        graphDirectPhotonThermalChatterjee2040->Write("graphDirectPhotonThermalYield_Chatterjee_2040");
        graphDirectPhotonPromptChatterjee2040->Write("graphDirectPhotonPromptYield_Chatterjee_2040");
        graphDirectPhotonChatterjee4060->Write("graphDirectPhotonYield_Chatterjee_4060");
        graphDirectPhotonThermalChatterjee4060->Write("graphDirectPhotonThermalYield_Chatterjee_4060");
        graphDirectPhotonPromptChatterjee4060->Write("graphDirectPhotonPromptYield_Chatterjee_4060");

        graphDirectPhotonChatterjee0020_2->Write("graphDirectPhotonYield_Chatterjee_0020_2");
        graphDirectPhotonThermalChatterjee0020_2->Write("graphDirectPhotonThermalYield_Chatterjee_0020_2");
        graphDirectPhotonPromptChatterjee0020_2->Write("graphDirectPhotonPromptYield_Chatterjee_0020_2");
        graphDirectPhotonChatterjee2040_2->Write("graphDirectPhotonYield_Chatterjee_2040_2");
        graphDirectPhotonThermalChatterjee2040_2->Write("graphDirectPhotonThermalYield_Chatterjee_2040_2");
        graphDirectPhotonPromptChatterjee2040_2->Write("graphDirectPhotonPromptYield_Chatterjee_2040_2");
        
        graphDirectPhotonHees0020->Write("graphDirectPhotonYield_VanHees_0020");
        graphDirectPhotonQGPHees0020->Write("graphDirectPhotonQGPYield_VanHees_0020");
        graphDirectPhotonPrimordialHees0020->Write("graphDirectPhotonPrimordialYield_VanHees_0020");
        graphDirectPhotonOmegaHees0020->Write("graphDirectPhotonOmegaYield_VanHees_0020");
        graphDirectPhotonMesonGasHees0020->Write("graphDirectPhotonMesonGasYield_VanHees_0020");
        graphDirectPhotonRhoSFHees0020->Write("graphDirectPhotonRhoSFYield_VanHees_0020");

        graphDirectPhotonHees2040->Write("graphDirectPhotonYield_VanHees_2040");
        graphDirectPhotonQGPHees2040->Write("graphDirectPhotonQGPYield_VanHees_2040");
        graphDirectPhotonPrimordialHees2040->Write("graphDirectPhotonPrimordialYield_VanHees_2040");
        graphDirectPhotonOmegaHees2040->Write("graphDirectPhotonOmegaYield_VanHees_2040");
        graphDirectPhotonMesonGasHees2040->Write("graphDirectPhotonMesonGasYield_VanHees_2040");
        graphDirectPhotonRhoSFHees2040->Write("graphDirectPhotonRhoSFYield_VanHees_2040");

        graphDirectPhotonHees4080->Write("graphDirectPhotonYield_VanHees_4080");
        graphDirectPhotonQGPHees4080->Write("graphDirectPhotonQGPYield_VanHees_4080");
        graphDirectPhotonPrimordialHees4080->Write("graphDirectPhotonPrimordialYield_VanHees_4080");
        graphDirectPhotonOmegaHees4080->Write("graphDirectPhotonOmegaYield_VanHees_4080");
        graphDirectPhotonMesonGasHees4080->Write("graphDirectPhotonMesonGasYield_VanHees_4080");
        graphDirectPhotonRhoSFHees4080->Write("graphDirectPhotonRhoSFYield_VanHees_4080");
        
        graphDirectPhotonHolopainen0020->Write("graphDirectPhotonYield_Holopainen_0020");
        graphDirectPhotonV2Holopainen0020->Write("graphDirectPhotonV2_Holopainen_0020");
        graphDirectPhotonHolopainen2040->Write("graphDirectPhotonYield_Holopainen_2040");
        graphDirectPhotonV2Holopainen2040->Write("graphDirectPhotonV2_Holopainen_2040");
        
    fileTheoryGraphsPbPb->Close();

	
}