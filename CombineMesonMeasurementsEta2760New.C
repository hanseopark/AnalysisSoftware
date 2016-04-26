/****************************************************************************************************************************
******          provided by Gamma Conversion Group, PWG4,                                                                                       
                *****
******          Ana Marin, marin@physi.uni-heidelberg.de                                                                                        
                *****
******          Kathrin Koch, kkoch@physi.uni-heidelberg.de                                                                                     
                *****
******          Friederike Bock, friederike.bock@cern.ch                                                                                        
                *****
*****************************************************************************************************************************/
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
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "CommonHeaders/PlottingGammaConversionHistos.h"
#include "CommonHeaders/PlottingGammaConversionAdditional.h"
#include "TaskV1/ProduceFinalResults.h"
#include "CommonHeaders/FittingGammaConversion.h"
#include "CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "CommonHeaders/ConversionFunctions.h"
#include "CommonHeaders/CombinationFunctions.h"
#include "CombineMesonMeasurementsEta2760.h"

extern TRandom* gRandom;
extern TBenchmark*      gBenchmark;
extern TSystem* gSystem;
extern TMinuit*         gMinuit;        

struct SysErrorConversion {
        Double_t value;
        Double_t error;
        //      TString name;
};

void CombineMesonMeasurementsEta2760New(TString fileNamePCM = "", 
					//TString fileNameEMCAL = "", 
					TString suffix = "eps", 
					TString isMC= "", 
					TString thesisPlots = "", 
					TString bWCorrection="X"){   

        date = ReturnDateString();
        
        gROOT->Reset(); 
        gROOT->SetStyle("Plain");
        
        StyleSettingsThesis();  
        SetPlotStyle();
        Double_t xSection2760GeVppINEL 	= 62.8*1e9;
        TString dateForOutput                   = ReturnDateStringForOutput();
        cout << dateForOutput.Data() << endl;
        //___________________________________ Declaration of files _____________________________________________

        collisionSystem2760GeV                  = "pp, #sqrt{#it{s}} = 2.76 TeV";               
        collisionSystemCombinedReallyAll = "pp #sqrt{#it{s}} = 0.9, 2.76, 7 TeV";
        fileNamePublishedPi0Allpp="CombinedResultsPaperX_18_Feb_2014.root";

	TString outputDir                               = Form("%s/%s/CombineMesonMeasurements%s",suffix.Data(),dateForOutput.Data(),bWCorrection.Data());
	gSystem->Exec("mkdir -p "+outputDir);
        nameFinalResDat                                 = Form("%s/CombinedResults%s_FitResults.dat",dateForOutput.Data(),bWCorrection.Data());
        cout << outputDir.Data() << endl;
        cout << fileNamePCM.Data() << endl;
        cout<< fileNamePublishedPi0Allpp.Data()<<endl;

        Bool_t thesis                                                   = kFALSE;
        if(thesisPlots.CompareTo("thesis") == 0){
                thesis                                                          = kTRUE;
        }


       if (isMC.CompareTo("kTRUE")==0){ 
                prefix2 = "MC";
                pictDrawingOptions[1] = kTRUE;
        } else {        
                prefix2 = "data";
                pictDrawingOptions[1] = kFALSE;
        }
        
        mesonMassExpectPi0 = TDatabasePDG::Instance()->GetParticle(111)->Mass();
        mesonMassExpectEta = TDatabasePDG::Instance()->GetParticle(221)->Mass();

 
    //************************** Read Theory*******************************************
	cout<< " Reading theory "<<endl;

	TFile* fileTheoryCompilation = new TFile("ExternalInput/TheoryCompilationPP.root");
	TH1F* histoPythia8InvXSection = (TH1F*) fileTheoryCompilation->Get("histoInvSecPythia8Spec2760GeV");
	TH1F* histoPythia8InvXSection_VarBinning = (TH1F*) fileTheoryCompilation->Get("histoInvSecPythia8Spec2760GeVVarBinning");
	graphNLOCalcMuHalf2760GeV=    (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecPi0MuHalf2760GeV");
	graphNLOCalcMuOne2760GeV=     (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecPi0MuOne2760GeV");
	graphNLOCalcMuTwo2760GeV=     (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecPi0MuTwo2760GeV");
	graphNLOCalcEtaMuHalf2760GeV=   (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecEtaMuHalf2760GeV");
	graphNLOCalcEtaMuOne2760GeV=    (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecEtaMuOne2760GeV");
	graphNLOCalcEtaMuTwo2760GeV=    (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecEtaMuTwo2760GeV");
	graphEtaToPi0NLOMuHalf2760GeV = (TGraph*)fileTheoryCompilation->Get("graphNLOCalcEtaOverPi0MuHalf2760GeV");
	graphEtaToPi0NLOMuOne2760GeV =  (TGraph*)fileTheoryCompilation->Get("graphNLOCalcEtaOverPi0MuOne2760GeV");
	graphEtaToPi0NLOMuTwo2760GeV =  (TGraph*)fileTheoryCompilation->Get("graphNLOCalcEtaOverPi0MuTwo2760GeV");
	graphNLODSSCalcMuTwo2760GeV=    (TGraph*)fileTheoryCompilation->Get("graphNLOCalcDSSInvSecPi0MuTwo2760GeV");





	//************************** Ratio + World data ******************************************** 
	cout << "Before reading WorldDataFile " << endl;
	fileWorldDataPi0Eta = new TFile("ExternalInput/WorldDataPi0Eta.root");
	graphDonaldson100GeV = (TGraphErrors*)fileWorldDataPi0Eta->Get("donaldson100GeV");
	DrawGammaSetMarkerTGraphErr(graphDonaldson100GeV, 24, 0.8, kGray+3, kGray+3);
	graphDonaldson200GeV = (TGraphErrors*)fileWorldDataPi0Eta->Get("donaldson200GeV");
	DrawGammaSetMarkerTGraphErr(graphDonaldson200GeV, 25, 0.8, kGray+3, kGray+3);
	graphAntille87pp = (TGraphErrors*)fileWorldDataPi0Eta->Get("Antille24.3GeVpp");
	DrawGammaSetMarkerTGraphErr(graphAntille87pp, 26, 0.8, kGray+3, kGray+3);
	graphAguilar400GeV = (TGraphErrors*)fileWorldDataPi0Eta->Get("Aguilar400GeV");
	DrawGammaSetMarkerTGraphErr(graphAguilar400GeV, 27, 0.8, kGray+3, kGray+3);
	graphKourkou79pp = (TGraphErrors*)fileWorldDataPi0Eta->Get("Kourkou30.6GeVpp");
	DrawGammaSetMarkerTGraphErr(graphKourkou79pp, 28, 0.8, kGreen+2, kGreen+2);
	graphApana530GeV = (TGraphErrors*)fileWorldDataPi0Eta->Get("Apana530GeV");
	DrawGammaSetMarkerTGraphErr(graphApana530GeV, 29, 0.8, kGray+3, kGray+3);
	graphKourkou79pp52 = (TGraphErrors*)fileWorldDataPi0Eta->Get("Kourkou52.7GeVpp");
	DrawGammaSetMarkerTGraphErr(graphKourkou79pp52, 30, 0.8, kGreen+2, kGreen+2);
	graphAkesson53GeVpp = (TGraphErrors*)fileWorldDataPi0Eta->Get("Akesson53GeVpp");
	DrawGammaSetMarkerTGraphErr(graphAkesson53GeVpp, 3, 0.8, kGreen+2, kGreen+2);
	graphKourkou79pp62 = (TGraphErrors*)fileWorldDataPi0Eta->Get("Kourkou79pp62");
	DrawGammaSetMarkerTGraphErr(graphKourkou79pp62, 22, 0.6, kGreen+2, kGreen+2);
	graphPhenix200GeV = (TGraphErrors*)fileWorldDataPi0Eta->Get("Phenix200GeV");
	DrawGammaSetMarkerTGraphErr(graphPhenix200GeV, 5, 0.8, kGreen+2, kGreen+2);
	cout << " done reading WorldDataFile" << endl;
	
	//---------------------------------------------------------------------------------



	//************************** Read data for PCM **************************************************
	//
	// Published pi0 pp 2.76TeV
	// --------------------------------------------
        filePublishedPi0Allpp = new TFile(fileNamePublishedPi0Allpp.Data());
	TGraphAsymmErrors * graphInvCrossSectionPi0PCM2760GeVStatErr = (TGraphAsymmErrors *)filePublishedPi0Allpp->Get("graphInvCrossSectionPi0PCM2760GeVStatErr");
	TGraphAsymmErrors * graphInvCrossSectionPi0PCM2760GeVSysErr  = (TGraphAsymmErrors *)filePublishedPi0Allpp->Get("graphInvCrossSectionPi0PCM2760GeVSysErr");
	TGraphAsymmErrors * graphInvCrossSectionPi0PCMStat2760GeVYShifted = (TGraphAsymmErrors *)filePublishedPi0Allpp->Get("graphInvCrossSectionPi0PCMStat2760GeV_YShifted");
        TGraphAsymmErrors * graphInvCrossSectionPi0PCMSys2760GeVYShifted  = (TGraphAsymmErrors *)filePublishedPi0Allpp->Get("graphInvCrossSectionPi0PCMSys2760GeV_YShifted");
	TGraphAsymmErrors * graphEtaToPi0Comb7TeVPublished  = (TGraphAsymmErrors *)filePublishedPi0Allpp->Get("graphEtaToPi0Comb7TeV");
	TGraphAsymmErrors * graphEtaToPi0Comb7TeVPublishedStat  = (TGraphAsymmErrors *)filePublishedPi0Allpp->Get("graphEtaToPi0Comb7TeVStat");
	TGraphAsymmErrors * graphEtaToPi0Comb7TeVPublishedSys  = (TGraphAsymmErrors *)filePublishedPi0Allpp->Get("graphEtaToPi0Comb7TeVSys");

	TGraphAsymmErrors * graphInvCrossSectionPi0PCMStat2760GeVYShiftedForFit = (TGraphAsymmErrors *)graphInvCrossSectionPi0PCMStat2760GeVYShifted->Clone();
	TGraphAsymmErrors * graphInvCrossSectionPi0PCMSys2760GeVYShiftedForFit  = (TGraphAsymmErrors *)graphInvCrossSectionPi0PCMSys2760GeVYShifted->Clone();

	TGraphAsymmErrors * graphInvCrossSectionEtaComb2760GeV_PrelimQM2011 = (TGraphAsymmErrors *)filePublishedPi0Allpp->Get("graphInvCrossSectionEtaComb2760GeV");
	TGraphAsymmErrors * graphInvCrossSectionEtaComb2760GeVStatErr_PrelimQM2011 = (TGraphAsymmErrors *)filePublishedPi0Allpp->Get("graphInvCrossSectionEtaComb2760GeVStatErr_PrelimQM2011");	
	TGraphAsymmErrors * graphInvCrossSectionEtaComb2760GeVSysErr_PrelimQM2011 = (TGraphAsymmErrors *)filePublishedPi0Allpp->Get("graphInvCrossSectionEtaComb2760GeVStatErr_PrelimQM2011");


	// cout<< "Published pi0 stat/sys"<<endl;
	// graphInvCrossSectionPi0PCMStat2760GeVYShiftedForFit->Print();
	// graphInvCrossSectionPi0PCMSys2760GeVYShiftedForFit->Print();
	// cout<<endl;
	// cout<<endl;
	Double_t paramPi0pp2760GeV[3]={1.1e11,5.8,0.13};
	fitInvCrossSectionPi0PCM2760GeVStat = FitObject("l","fitLevyInvCrossSectionPi0PCM2760GeVStat","Pi0",graphInvCrossSectionPi0PCMStat2760GeVYShiftedForFit,0.4,8., paramPi0pp2760GeV,"QNRMEX0+");
        DrawGammaSetMarkerTF1(fitInvCrossSectionPi0PCM2760GeVStat, 1, 1.5, kRed);


	fitInvCrossSectionPi0PCM2760GeVSys = FitObject("l","fitLevyInvCrossSectionPi0PCM2760GeVSys","Pi0",graphInvCrossSectionPi0PCMSys2760GeVYShiftedForFit,0.4,8.,paramPi0pp2760GeV,"QNRMEX0+");
	DrawGammaSetMarkerTF1(fitInvCrossSectionPi0PCM2760GeVSys, 1, 1.5, kBlue);


	//---------------------------------------------------------
	//Eta pp 2.76 new analysis with pile-up correction

        filePCM =               new TFile(fileNamePCM.Data());
        cout<< "Start reading histos from files-eta"<<endl;

        //Eta from new analysis
        histoPCMNumberOfEventsEta2760GeV=                    (TH1D*)filePCM->Get("histoNumberOfEvents2.76TeV");
        directoryPCMEta2760GeV =                             (TDirectory*)filePCM->Get("Eta2.76TeV"); 
        histoPCMCorrectedYieldEta2760GeV=                    (TH1D*)directoryPCMEta2760GeV->Get("CorrectedYieldEta");
        histoPCMUnCorrectedYieldEta2760GeV=                  (TH1D*)directoryPCMEta2760GeV->Get("RAWYieldPerEventsEta");
        histoPCMMassMesonEta2760GeV=                         (TH1D*)directoryPCMEta2760GeV->Get("MassEta");
        histoPCMFWHMMesonEtaMeV2760GeV =                     (TH1D*)directoryPCMEta2760GeV->Get("FWHMEtaMeV");
        histoPCMTrueMassMesonEta2760GeV =                    (TH1D*)directoryPCMEta2760GeV->Get("TrueMassEta");
        histoPCMTrueFWHMMesonEtaMeV2760GeV =                 (TH1D*)directoryPCMEta2760GeV->Get("TrueFWHMEtaMeV");
        histoPCMEtaCorrYieldBinShifted2760GeV=               (TH1D*)directoryPCMEta2760GeV->Get("CorrectedYieldEtaBinShifted");
        histoPCMEtaCorrYieldMtBinShifted2760GeV=             (TH1D*)directoryPCMEta2760GeV->Get("CorrectedYieldEtaMtBinShifted");
        histoPCMInvCrossSectionEta2760GeV=                   (TH1D*)directoryPCMEta2760GeV->Get("InvCrossSectionEta");     // non binshifted
        histoPCMRatioEtaPi02760GeV=                          (TH1D*)directoryPCMEta2760GeV->Get("EtatoPi0RatioConversionBinShifted");          
        graphPCMInvCrossSectionSysAEta2760GeV=               (TGraphAsymmErrors*)directoryPCMEta2760GeV->Get("InvCrossSectionEtaSysA");
        graphPCMInvCrossSectionSysEta2760GeV=                (TGraphAsymmErrors*)directoryPCMEta2760GeV->Get("InvCrossSectionEtaSys");
        graphPCMCorrectedYieldEtaSysErr2760GeV=              (TGraphAsymmErrors*)directoryPCMEta2760GeV->Get("EtaSystError");
        graphPCMCorrectedYieldEtaSysErrBinShifted2760GeV=    (TGraphAsymmErrors*)directoryPCMEta2760GeV->Get("EtaSystErrorBinShifted");
	//        graphPCMCorrectedYieldEtaSysErrMtBinShifted2760GeV=  (TGraphAsymmErrors*)directoryPCMEta2760GeV->Get("EtaSystErrorBinShiftedMt");
        graphPCMSystErrRatio2760GeV=                         (TGraphAsymmErrors*)directoryPCMEta2760GeV->Get("EtatoPi0RatioConversionBinShiftedSys");

	maxPtPCMEta2760GeV = histoPCMInvCrossSectionEta2760GeV->GetXaxis()->GetBinUpEdge(histoPCMInvCrossSectionEta2760GeV->GetNbinsX());

	graphPCMInvCrossSectionSysEta2760GeV->Print();
        cout<< "Done reading histos from eta-PCM file"<<endl;


     //************************** Read data for EMCAL **************************************************
	//	fileEMCAL =               new TFile("ExternalInput/EMCAL/2.76TeV/Eta2760_ppReferenceSpectraEMCal7August2015-TH1.root");
	//        fileEMCAL =               new TFile("ExternalInput/EMCAL/2.76TeV/pp-Ref-RAA-EMCal21August2015-TH1.root");
	//fileEMCALforComp = new TFile("ExternalInput/EMCAL/2.76TeV/ppReferenceSpectraEMCal24August2015LowandHighPt.root");

	fileEMCAL =               new TFile("ExternalInput/EMCAL/2.76TeV/ppReferenceSpectraEMCal27August2015.root");
	fileEMCALforComp = new TFile("ExternalInput/EMCAL/2.76TeV/ppReferenceSpectraEMCal27LowandHighPt.root");

    //************************************************************************************************
	histoEMCALCorrectedYieldEta2760GeV = (TH1D*)fileEMCAL->Get("RefStatEta");
	histoEMCALCorrectedYieldSysEta2760GeV = (TH1D*)fileEMCAL->Get("RefSysEta");

	graphEMCALCorrectedYieldStatEta2760GeV = new TGraphAsymmErrors(histoEMCALCorrectedYieldEta2760GeV); 
	graphEMCALCorrectedYieldSysEta2760GeV = new TGraphAsymmErrors(histoEMCALCorrectedYieldSysEta2760GeV);

	maxPtEta2760GeV = histoEMCALCorrectedYieldEta2760GeV->GetXaxis()->GetBinUpEdge(histoEMCALCorrectedYieldEta2760GeV->GetNbinsX()-1);


	histoEMCALCorrectedYieldEta2760GeVforComp = (TH1D*)fileEMCALforComp->Get("RefStatEta");
	histoEMCALCorrectedYieldSysEta2760GeVforComp = (TH1D*)fileEMCALforComp->Get("RefSysEta");
	graphEMCALCorrectedYieldStatEta2760GeVforComp = new TGraphAsymmErrors(histoEMCALCorrectedYieldEta2760GeVforComp); 
	graphEMCALCorrectedYieldSysEta2760GeVforComp = new TGraphAsymmErrors(histoEMCALCorrectedYieldSysEta2760GeVforComp);

	cout<< "Done reading histos from eta-EMCAL file"<<endl;

	graphEMCALStatErrInvCrossSectionEta2760GeV = (TGraphAsymmErrors*) graphEMCALCorrectedYieldStatEta2760GeV->Clone("graphEMCALStatErrInvCrossSectionEta2760GeV");
	graphEMCALStatErrInvCrossSectionEta2760GeV = ScaleGraph(graphEMCALCorrectedYieldStatEta2760GeV,xSection2760GeVppINEL);

	graphEMCALSysErrInvCrossSectionEta2760GeV = (TGraphAsymmErrors*) graphEMCALCorrectedYieldSysEta2760GeV->Clone("graphEMCALSysErrInvCrossSectionEta2760GeV");
	graphEMCALSysErrInvCrossSectionEta2760GeV = ScaleGraph(graphEMCALCorrectedYieldSysEta2760GeV,xSection2760GeVppINEL);
	graphEMCALInvCrossSectionSysAEta2760GeV = (TGraphAsymmErrors*) graphEMCALSysErrInvCrossSectionEta2760GeV ->Clone("graphEMCALInvCrossSectionSysAEta2760GeV");




	histoEMCALInvCrossSectionStatEta2760GeV=(TH1D*)histoEMCALCorrectedYieldEta2760GeV->Clone();
	histoEMCALInvCrossSectionStatEta2760GeV->Scale(xSection2760GeVppINEL);

	histoEMCALInvCrossSectionSysEta2760GeV=(TH1D*)histoEMCALCorrectedYieldSysEta2760GeV->Clone();
	histoEMCALInvCrossSectionSysEta2760GeV->Scale(xSection2760GeVppINEL);


	graphEMCALStatErrInvCrossSectionEta2760GeVforComp = (TGraphAsymmErrors*) graphEMCALCorrectedYieldStatEta2760GeVforComp->Clone("graphEMCALStatErrInvCrossSectionEta2760GeVforComp");
	graphEMCALStatErrInvCrossSectionEta2760GeVforComp = ScaleGraph(graphEMCALCorrectedYieldStatEta2760GeVforComp,xSection2760GeVppINEL);

	graphEMCALSysErrInvCrossSectionEta2760GeVforComp = (TGraphAsymmErrors*) graphEMCALCorrectedYieldSysEta2760GeVforComp->Clone("graphEMCALSysErrInvCrossSectionEta2760GeVforComp");
	graphEMCALSysErrInvCrossSectionEta2760GeVforComp = ScaleGraph(graphEMCALCorrectedYieldSysEta2760GeVforComp,xSection2760GeVppINEL);



        nEvtPCMEta2760GeV =  histoPCMNumberOfEventsEta2760GeV->GetBinContent(1);
        cout<<"number of events eta-file::"<<nEvtPCMEta2760GeV<<endl; 

	//**************** CombinePoints 2760 GeV***************************************************//
	//**************** Eta cross section  sqrt(s)=2760GeV  should be combined with EMCal********//
        //**************** At present taken as PCM only************************* *******************//
        //**************** Graph created from Cross section Histogram that had not binshift********//

	cout<< "before shifts in x for eta"<<endl;
	graphInvCrossSectionEtaComb2760GeV = CombinePtPointsSpectra( histoPCMInvCrossSectionEta2760GeV, graphPCMInvCrossSectionSysEta2760GeV,
								     histoEMCALInvCrossSectionStatEta2760GeV, graphEMCALSysErrInvCrossSectionEta2760GeV,
     	 				                             graphInvCrossSectionEtaComb2760GeVStatErr, graphInvCrossSectionEtaComb2760GeVSysErr,
								     xPtLimits2760GeVNew, 13, 1, 0,6,kFALSE,1,kTRUE);    
	                                                             //xPtLimits2760GeVNew, 13, 1, 0,6,kTRUE,1,kTRUE);    // AM-To be used with 2013 eta syst. errors 
								     //xPtLimits2760GeVNew, 13, 1, 0,5,kTRUE,1,kFALSE); // AM-To be used in the 4-6 GeV bin is joined PCM-EMCAL    


	// graphInvCrossSectionEtaComb2760GeV->Print();
	// graphInvCrossSectionEtaComb2760GeVStatErr->Print();
	// graphInvCrossSectionEtaComb2760GeVSysErr->Print();

	//AM august to be checked
	// graphInvCrossSectionEtaComb2760GeVStatErr = new TGraphAsymmErrors(histoPCMInvCrossSectionEta2760GeV);	
	// graphInvCrossSectionEtaComb2760GeVStatErr->Print();
	// graphInvCrossSectionEtaComb2760GeVStatErr->RemovePoint(0);

	// cout<< "Reading Eta syst. errors"<< endl;

        // relSystErrorEta2760GeVDown = ExtractRelErrDownAsymmGraph(graphPCMCorrectedYieldEtaSysErr2760GeV);
        // relSystErrorEta2760GeVUp = ExtractRelErrUpAsymmGraph(graphPCMCorrectedYieldEtaSysErr2760GeV);
	// graphPCMCorrectedYieldEtaSysErr2760GeV->RemovePoint(graphPCMCorrectedYieldEtaSysErr2760GeV->GetN()-1);

        // nPointsPCMEta2760GeV = graphPCMCorrectedYieldEtaSysErr2760GeV->GetN();


	// graphInvCrossSectionEtaComb2760GeV = CalculateSysErrFromRelSysHistoComplete( histoPCMInvCrossSectionEta2760GeV, 
	// 									     "graphInvCrossSectionEtaComb2760GeV",
	// 									     relSystErrorEta2760GeVDown, 
	// 									     relSystErrorEta2760GeVUp, 2, nPointsPCMEta2760GeV);
	// graphInvCrossSectionEtaComb2760GeVSysErr =  (TGraphAsymmErrors*)graphPCMInvCrossSectionSysEta2760GeV->Clone("graphInvCrossSectionEtaComb2760GeVSysErr");
	// graphInvCrossSectionEtaComb2760GeVSysErr->RemovePoint(graphInvCrossSectionEtaComb2760GeVSysErr->GetN()-1);
	// graphInvCrossSectionEtaComb2760GeVSysErr->Print();



  	TGraphAsymmErrors* graphInvCrossSectionEtaComb2760GeVUnShifted = (TGraphAsymmErrors*)graphInvCrossSectionEtaComb2760GeV->Clone("Unshifted"); 
	TGraphAsymmErrors* graphInvCrossSectionEtaComb2760GeVStatErrUnShifted = (TGraphAsymmErrors*)graphInvCrossSectionEtaComb2760GeVStatErr->Clone("UnshiftedStat"); 
	TGraphAsymmErrors* graphInvCrossSectionEtaComb2760GeVSysErrUnShifted = (TGraphAsymmErrors*)graphInvCrossSectionEtaComb2760GeVSysErr->Clone("UnshiftedSys"); 


	cout<< "New graph PCM "<<endl;
 	graphPCMXSectionEta2760GeV = new TGraphAsymmErrors(histoPCMInvCrossSectionEta2760GeV);
 	graphPCMXSectionEta2760GeV->RemovePoint(0);
	// 	graphPCMXSectionEta2760GeV->RemovePoint(graphPCMXSectionEta2760GeV->GetN()-1);

 	TGraphAsymmErrors* graphPCMXSectionEta2760GeVUnShifted = (TGraphAsymmErrors*)graphPCMXSectionEta2760GeV->Clone("UnshiftedPCM"); 
 	//graphPCMInvCrossSectionSysEta2760GeV->RemovePoint(graphPCMInvCrossSectionSysEta2760GeV->GetN()-1);
 	TGraphAsymmErrors* graphPCMXSectionSysEta2760GeVUnShifted = (TGraphAsymmErrors*)graphPCMInvCrossSectionSysEta2760GeV->Clone("UnshiftedSysPCM"); 
	graphPCMXSectionEta2760GeV->Print();

	TGraphAsymmErrors* graphPCMXSectionStatEta2760GeVXShifted = (TGraphAsymmErrors*)graphPCMXSectionEta2760GeV->Clone("graphPCMXSectionEta2760GeVXShifted"); 
 	TGraphAsymmErrors* graphPCMXSectionSysEta2760GeVXShifted = (TGraphAsymmErrors*)graphPCMInvCrossSectionSysEta2760GeV->Clone("graphPCMXSectionSysEta2760GeVXShifted"); 


        cout<< endl;

	//*************EMCAL should come here******************
	cout<< "New graph EMCAL "<<endl;

 	graphEMCALXSectionEta2760GeV = new TGraphAsymmErrors(histoEMCALInvCrossSectionStatEta2760GeV); 
 	TGraphAsymmErrors* graphEMCALXSectionEta2760GeVUnShifted = (TGraphAsymmErrors*)graphEMCALXSectionEta2760GeV->Clone("UnshiftedEMCAL"); 
 	TGraphAsymmErrors* graphEMCALXSectionSysEta2760GeVUnShifted = (TGraphAsymmErrors*)graphEMCALSysErrInvCrossSectionEta2760GeV ->Clone("UnshiftedSysEMCAL"); 

	TGraphAsymmErrors* graphEMCALXSectionStatEta2760GeVXShifted = (TGraphAsymmErrors*)graphEMCALXSectionEta2760GeV->Clone("graphEMCALXSectionStatEta2760GeVXShifted"); 
 	TGraphAsymmErrors* graphEMCALXSectionSysEta2760GeVXShifted = (TGraphAsymmErrors*)graphEMCALSysErrInvCrossSectionEta2760GeV ->Clone("graphEMCALXSectionSysEta2760GeVXShifted"); 


	cout<< "after"<< graphInvCrossSectionEtaComb2760GeVStatErr->GetN()<< " " << graphInvCrossSectionEtaComb2760GeVSysErr->GetN()<< endl;     
	cout<<"before x shift-stat errors/total"<<endl;
	graphInvCrossSectionEtaComb2760GeVStatErr->Print();
	graphInvCrossSectionEtaComb2760GeV->Print();

	Double_t paramGraphEta[3] = {1.0e10,7.,0.13};
	Double_t ParametersEta2760GeV[5] = { 2.4e+10, 0.3, 1.,0.3,3.8};
	//Double_t paramGraphEta[3] = {2.4e10,4.7,0.11};
	if(bWCorrection.CompareTo("X")==0 ){
	  //	  cout<< "doing xshift::"<<graphInvCrossSectionEtaComb2760GeV->GetName()<< endl;
	  TF1* fitTsallisEta2760GeVPt = FitObject("tmpt","tmptTsallisWithPtEta2760GeV","Eta");
	  fitTsallisEta2760GeVPt->SetParameters(paramGraphEta[0],paramGraphEta[1], paramGraphEta[2]) ; // standard parameter optimize if necessary
	  
	  graphInvCrossSectionEtaComb2760GeV   = ApplyXshift(graphInvCrossSectionEtaComb2760GeV ,fitTsallisEta2760GeVPt,"Eta");

	  cout<< endl;
	  graphInvCrossSectionEtaComb2760GeV->Print();
	  cout<< endl;

	  graphInvCrossSectionEtaComb2760GeVStatErr = ApplyXshiftIndividualSpectra (graphInvCrossSectionEtaComb2760GeV, 
	  									    graphInvCrossSectionEtaComb2760GeVStatErr, 
	  									    fitTsallisEta2760GeVPt, 0, graphInvCrossSectionEtaComb2760GeVStatErr->GetN(),"Eta");
	  cout<<endl;
	  graphInvCrossSectionEtaComb2760GeVStatErr->Print(),
	    cout<< endl;

	  cout<<"before shift"<<endl;
	  graphInvCrossSectionEtaComb2760GeVSysErr->Print();
	  graphInvCrossSectionEtaComb2760GeVSysErr = ApplyXshiftIndividualSpectra (graphInvCrossSectionEtaComb2760GeV, 
	  									   graphInvCrossSectionEtaComb2760GeVSysErr, 
	  									   fitTsallisEta2760GeVPt, 0, graphInvCrossSectionEtaComb2760GeVSysErr->GetN(),"Eta");
	  
	  cout<<"after shift"<<endl;
	  graphInvCrossSectionEtaComb2760GeVSysErr->Print();
	  cout<< endl;

	  // Apply X-shift to PCM

 	  
	  graphPCMXSectionStatEta2760GeVXShifted = ApplyXshiftIndividualSpectra (graphInvCrossSectionEtaComb2760GeV, 
	  									 graphPCMXSectionStatEta2760GeVXShifted, 
	  									 fitTsallisEta2760GeVPt, 0, graphPCMXSectionStatEta2760GeVXShifted->GetN(),"Eta");

	  graphPCMXSectionSysEta2760GeVXShifted = ApplyXshiftIndividualSpectra (graphInvCrossSectionEtaComb2760GeV, 
	  									 graphPCMXSectionSysEta2760GeVXShifted, 
	  									 fitTsallisEta2760GeVPt, 0, graphPCMXSectionSysEta2760GeVXShifted->GetN(),"Eta");

	  // Apply X-shift to EMCAL

	  cout<< " going to apply x shift to EMCAL spectrum"<< endl;
	  graphEMCALXSectionStatEta2760GeVXShifted->Print();
	  graphEMCALXSectionStatEta2760GeVXShifted = ApplyXshiftIndividualSpectra (graphInvCrossSectionEtaComb2760GeV, 
	  									 graphEMCALXSectionStatEta2760GeVXShifted, 
										   fitTsallisEta2760GeVPt, 5, 12,"Eta");
	  graphEMCALXSectionStatEta2760GeVXShifted->Print();
	  cout<<endl;

	  graphEMCALXSectionSysEta2760GeVXShifted = ApplyXshiftIndividualSpectra (graphInvCrossSectionEtaComb2760GeV, 
	  									 graphEMCALXSectionSysEta2760GeVXShifted, 
										  fitTsallisEta2760GeVPt, 5, 12,"Eta");
	  graphEMCALXSectionSysEta2760GeVXShifted->Print();
	  cout<<endl;
	  /*
	   TF1* fitBylinkinEta2760GeVPt = FitObject("tcm","BylinkinWithPtEta2760GeVtcm","Eta");
	   fitBylinkinEta2760GeVPt->SetParameters(ParametersEta2760GeV[0],ParametersEta2760GeV[1], ParametersEta2760GeV[2],ParametersEta2760GeV[3],ParametersEta2760GeV[4]) ; // standard parameter optimize if necessary
	   graphInvCrossSectionEtaComb2760GeV   = ApplyXshift(graphInvCrossSectionEtaComb2760GeV ,fitBylinkinEta2760GeVPt,"Eta");
	  */
	}

        cout<< " after applying X-shift"<<endl;


	graphInvCrossSectionEtaComb2760GeV->SetMarkerStyle(markerStyleConv);
	graphInvCrossSectionEtaComb2760GeV->SetMarkerColor(1);
	graphInvCrossSectionEtaComb2760GeV->SetLineColor(1);
	graphInvCrossSectionEtaComb2760GeV->SetLineWidth(5);
	graphInvCrossSectionEtaComb2760GeV->SetFillColor(0);
	graphInvCrossSectionEtaComb2760GeV->SetFillStyle(0);
	graphInvCrossSectionEtaComb2760GeVUnscaled = (TGraphAsymmErrors*)graphInvCrossSectionEtaComb2760GeV->Clone("graphInvCrossSectionEtaComb2760GeVUnscaled");
	graphInvCrossSectionEtaComb2760GeVStatErrUnscaled = (TGraphAsymmErrors*)graphInvCrossSectionEtaComb2760GeVStatErr->Clone("graphInvCrossSectionEtaComb2760GeVStatErrUnscaled");
	graphInvCrossSectionEtaComb2760GeVSysErrUnscaled = (TGraphAsymmErrors*)graphInvCrossSectionEtaComb2760GeVSysErr->Clone("graphInvCrossSectionEtaComb2760GeVSysErrUnscaled");


	TGraphAsymmErrors * graphInvCrossSectionEtaComb2760GeVForFit=(TGraphAsymmErrors*) graphInvCrossSectionEtaComb2760GeV->Clone("graphInvCrossSectionEtaComb2760GeVForFit");

	//	fitInvCrossSectionEta2760GeV = FitObject("l","fitInvCrossSectionEta2760GeV","Eta",histoPCMInvCrossSectionEta2760GeV,0.6,20.,paramGraphEta,"NQRME+","fixn");

	fitInvCrossSectionEta2760GeV = FitObject("l","fitInvCrossSectionEta2760GeV","Eta",graphInvCrossSectionEtaComb2760GeVForFit,0.5,20.,paramGraphEta,"QNRMEX0+");
	//	graphInvCrossSectionEtaComb2760GeV->Fit(fitInvCrossSectionEta2760GeV,"QNRMEX0+","",0.6,20.);
	paramGraphEta[1] =fitInvCrossSectionEta2760GeV->GetParameter(1);	

	cout << WriteParameterToFile(fitInvCrossSectionEta2760GeV)<< endl;
	forOutput= WriteParameterToFile(fitInvCrossSectionEta2760GeV);
	fileFinalResults<< forOutput.Data()<< endl;

	//        cout << "Hola 1 after xshift "<< endl;
	//---------------------Yshifts for eta@2760-------------------------//
	
       	TF1* fitTsallisEta2760GeVPtYShift = ApplyYShift(graphInvCrossSectionEtaComb2760GeVUnShifted,
							&graphInvCrossSectionEtaComb2760GeVYShifted, "l","InvXsectionY",0.3, paramGraphEta,0.00001,kFALSE,"Eta");
	// cout << "Unshifted/yshifted" << endl;
	// graphInvCrossSectionEtaComb2760GeVUnShifted->Print();
        // cout << endl;  
	// graphInvCrossSectionEtaComb2760GeVYShifted->Print();
	// cout << endl;  
	// cout << endl;  

	graphInvCrossSectionEtaComb2760GeVYShiftedSys = (TGraphAsymmErrors*)graphInvCrossSectionEtaComb2760GeVSysErrUnShifted->Clone("YShiftedCombSys");
	graphInvCrossSectionEtaComb2760GeVYShiftedSys= ApplyYshiftIndividualSpectra( graphInvCrossSectionEtaComb2760GeVYShiftedSys, fitTsallisEta2760GeVPtYShift);
        
	graphInvCrossSectionEtaComb2760GeVYShiftedStat = (TGraphAsymmErrors*)graphInvCrossSectionEtaComb2760GeVStatErrUnShifted->Clone("YShiftedCombStat");
	graphInvCrossSectionEtaComb2760GeVYShiftedStat= ApplyYshiftIndividualSpectra( graphInvCrossSectionEtaComb2760GeVYShiftedStat, fitTsallisEta2760GeVPtYShift);
        
	//apply y-shifts to PCM
	graphPCMXSectionEta2760GeVYShifted = (TGraphAsymmErrors*)graphPCMXSectionEta2760GeVUnShifted->Clone("YShiftedPCM");
	graphPCMXSectionEta2760GeVYShifted= ApplyYshiftIndividualSpectra( graphPCMXSectionEta2760GeVYShifted, fitTsallisEta2760GeVPtYShift);
        
	graphPCMXSectionEta2760GeVYShiftedSys = (TGraphAsymmErrors*)graphPCMXSectionSysEta2760GeVUnShifted->Clone("YShiftedPCMSys");
	graphPCMXSectionEta2760GeVYShiftedSys = ApplyYshiftIndividualSpectra( graphPCMXSectionEta2760GeVYShiftedSys, fitTsallisEta2760GeVPtYShift);
        
	graphPCMXSectionEta2760GeVYShiftedSysRAA = (TGraphAsymmErrors*)graphPCMInvCrossSectionSysAEta2760GeV->Clone("YShiftedPCMSysRAA");
	graphPCMXSectionEta2760GeVYShiftedSysRAA = ApplyYshiftIndividualSpectra( graphPCMXSectionEta2760GeVYShiftedSysRAA, fitTsallisEta2760GeVPtYShift);

	// cout<< " PCM stat/sys y shifted"<<endl;
	// graphPCMXSectionEta2760GeVYShifted->Print();
	// graphPCMXSectionEta2760GeVYShiftedSys->Print();
	// cout<< " after PCM sys/stat"<<endl;


	//apply y-shifts to EMCAL
                                                                   
	graphEMCALXSectionEta2760GeVYShifted = (TGraphAsymmErrors*)graphEMCALXSectionEta2760GeVUnShifted->Clone("YShiftedEMCAL");
	graphEMCALXSectionEta2760GeVYShifted= ApplyYshiftIndividualSpectra( graphEMCALXSectionEta2760GeVYShifted, fitTsallisEta2760GeVPtYShift);
        
	graphEMCALXSectionEta2760GeVYShiftedSys = (TGraphAsymmErrors*)graphEMCALXSectionSysEta2760GeVUnShifted->Clone("YShiftedEMCALSys");
	graphEMCALXSectionEta2760GeVYShiftedSys = ApplyYshiftIndividualSpectra( graphEMCALXSectionEta2760GeVYShiftedSys, fitTsallisEta2760GeVPtYShift);
        
	graphEMCALXSectionEta2760GeVYShiftedSysRAA = (TGraphAsymmErrors*)graphEMCALInvCrossSectionSysAEta2760GeV->Clone("YShiftedEMCALSysRAA");
	graphEMCALXSectionEta2760GeVYShiftedSysRAA = ApplyYshiftIndividualSpectra( graphEMCALXSectionEta2760GeVYShiftedSysRAA, fitTsallisEta2760GeVPtYShift);


	graphEMCALXSectionEta2760GeVYShifted->SetMarkerStyle(20);
	graphPCMXSectionEta2760GeVYShifted->SetMarkerStyle(21);
	graphEMCALXSectionEta2760GeVYShifted->SetMarkerColor(kBlue);
	graphEMCALXSectionEta2760GeVYShifted->SetLineColor(kBlue);
	graphPCMXSectionEta2760GeVYShifted->SetMarkerColor(kGreen);
	graphPCMXSectionEta2760GeVYShifted->SetLineColor(kGreen);
	graphPCMXSectionEta2760GeVYShifted->SetLineWidth(3);

	graphEMCALXSectionEta2760GeVUnShifted->SetMarkerStyle(20);
	graphPCMXSectionEta2760GeVUnShifted->SetMarkerStyle(21);
	graphEMCALXSectionEta2760GeVUnShifted->SetMarkerColor(kBlue);
	graphEMCALXSectionEta2760GeVUnShifted->SetLineColor(kBlue);
	graphPCMXSectionEta2760GeVUnShifted->SetMarkerColor(kGreen);
	graphPCMXSectionEta2760GeVUnShifted->SetLineColor(kGreen);
	graphPCMXSectionEta2760GeVUnShifted->SetLineWidth(3);

	graphEMCALSysErrInvCrossSectionEta2760GeVforComp->SetMarkerStyle(20);
	graphEMCALSysErrInvCrossSectionEta2760GeVforComp->SetMarkerColor(kBlue-7);
	graphEMCALSysErrInvCrossSectionEta2760GeVforComp->SetLineColor(kBlue-7);

	//--------------------EtaPi0------------------------------------//

	TGraphAsymmErrors* graphRatioEtaPCMtoPi0Published2760GeVSys= (TGraphAsymmErrors*)graphPCMXSectionEta2760GeVYShiftedSys->Clone();
	graphRatioEtaPCMtoPi0Published2760GeVSys = CalculateGraphErrRatioToFit (graphPCMXSectionEta2760GeVYShiftedSys, fitInvCrossSectionPi0PCM2760GeVStat); 

        TGraphAsymmErrors* graphRatioEtaPCMtoPi0Published2760GeVStat= (TGraphAsymmErrors*)graphPCMXSectionEta2760GeVYShifted->Clone();
	graphRatioEtaPCMtoPi0Published2760GeVStat = CalculateGraphErrRatioToFit (graphPCMXSectionEta2760GeVYShifted, fitInvCrossSectionPi0PCM2760GeVStat); 

	graphRatioEtaPCMtoPi0Published2760GeVSys->RemovePoint(6);
	graphRatioEtaPCMtoPi0Published2760GeVStat->RemovePoint(6);
	// cout<<"stat errors"<<endl;
	// graphRatioEtaPCMtoPi0Published2760GeVStat->Print();	
	// cout<<"sys errors"<<endl;
	// graphRatioEtaPCMtoPi0Published2760GeVSys->Print();
	
	graphPCMSystErrRatio2760GeV->RemovePoint(graphPCMSystErrRatio2760GeV->GetN()-1);
        TGraphAsymmErrors* dummyGraphSystErr = (TGraphAsymmErrors*)graphPCMSystErrRatio2760GeV->Clone("dummyGraphSystErr");
	//	dummyGraphSystErr->Print(); 
	relSystErrorEtaPi02760GeVDown = ExtractRelErrDownAsymmGraph(dummyGraphSystErr);
        relSystErrorEtaPi02760GeVUp = ExtractRelErrUpAsymmGraph(dummyGraphSystErr);
        nPointsEta2760GeV = dummyGraphSystErr->GetN();
        relSystErrorEtaPi02760GeVUp[nPointsEta2760GeV-1]=relSystErrorEtaPi02760GeVUp[nPointsEta2760GeV-2];
        relSystErrorEtaPi02760GeVDown[nPointsEta2760GeV-1]=relSystErrorEtaPi02760GeVDown[nPointsEta2760GeV-2];
	
	//	graphRatioEtaPi0StatErr2760GeV->RemovePoint(nPointsEta2760GeV-1);
        cout << "systematic errors Eta Conv 2.76TeV" << endl;
        for (Int_t i = 0; i < nPointsEta2760GeV; i++){
                cout << relSystErrorEtaPi02760GeVDown[i] << "\t" << relSystErrorEtaPi02760GeVUp[i] << endl;
        }
	
	
        graphRatioEtaPi0ComplErr2760GeV = CalculateSysErrFromRelSysHistoComplete( histoPCMRatioEtaPi02760GeV , "EtaPi0ComplError2760GeV",relSystErrorEtaPi02760GeVDown , relSystErrorEtaPi02760GeVDown, 2, nPointsEta2760GeV);
        graphRatioEtaPi0StatErr2760GeV = new TGraphAsymmErrors(histoPCMRatioEtaPi02760GeV);
        graphRatioEtaPi0StatErr2760GeV->RemovePoint(0);
	//        graphRatioEtaPi0StatErr2760GeV->RemovePoint(graphRatioEtaPi0StatErr2760GeV->GetN()-1);
        // cout << "stat Eta/Pi0 2.76 TeV" << endl;
        // graphRatioEtaPi0StatErr2760GeV->Print();

        graphRatioEtaPi0SysErr2760GeV = (TGraphAsymmErrors*)graphPCMSystErrRatio2760GeV->Clone("graphRatioEtaPi0SysErr2760GeV");
	

	
       //----------------------- Eta 2760GeV NLO mu = pt/2 ------------------------------
        graphNLOMuHalfEta2760GeV = (TGraph*)graphNLOCalcEtaMuHalf2760GeV->Clone();
        Double_t* dummyArrayXValuesNLO = graphNLOMuHalfEta2760GeV->GetX();
        Int_t  dummyArrayN = graphNLOMuHalfEta2760GeV->GetN();
        Int_t numberOfBins = 0;
        Int_t n = 0;
        while (n!=dummyArrayN ){
	  //                if (dummyArrayXValuesNLO[n] <= 7.5){
                if (dummyArrayXValuesNLO[n] <= 30){
                        n++;
                        numberOfBins++;
                        //cout << n << endl;
                } else {
                        n = dummyArrayN;
                } 
        }
	//  cout << "Number of bins Eta 2.76TeV "<<numberOfBins << endl;
	
	
       //----------------------- Eta 2760GeV NLO mu = pt/2 ------------------------------
	//        graphNLOMuHalfEta2760GeV = ScaleGraph(graphNLOMuHalfEta2760GeV,xSection2760GeV*recalcBarn);
        DrawGammaSetMarkerTGraph(graphNLOMuHalfEta2760GeV, styleMarkerNLOMuHalf, sizeMarkerNLO, colorNLOEta2760GeVMuHalf, colorNLOEta2760GeVMuHalf );

        for (Int_t i = numberOfBins; i < dummyArrayN; i++){
                graphNLOMuHalfEta2760GeV->RemovePoint(numberOfBins);
        }

        graphNLOMuHalfEta2760GeV->RemovePoint(0);
        graphNLOMuHalfEta2760GeV->RemovePoint(0);
        //------------------------- Eta 2760GeV NLO mu = pt ----------------------------
        graphNLOMuOneEta2760GeV = (TGraph*)graphNLOCalcEtaMuOne2760GeV->Clone();
        DrawGammaSetMarkerTGraph(graphNLOMuOneEta2760GeV,  styleMarkerNLOMuOne, sizeMarkerNLO, colorNLOEta2760GeVMuOne, colorNLOEta2760GeVMuOne);
        for (Int_t i = numberOfBins; i < dummyArrayN; i++){
                graphNLOMuOneEta2760GeV->RemovePoint(numberOfBins);
        }
	
	
        //------------------------- Eta 2760GeV NLO mu = 2pt -----------------------------
        graphNLOMuTwoEta2760GeV = (TGraph*)graphNLOCalcEtaMuTwo2760GeV->Clone();
        DrawGammaSetMarkerTGraph(graphNLOMuTwoEta2760GeV, styleMarkerNLOMuTwo, sizeMarkerNLO, colorNLOEta2760GeVMuTwo, colorNLOEta2760GeVMuTwo);
        for (Int_t i = numberOfBins; i < dummyArrayN; i++){
                graphNLOMuTwoEta2760GeV->RemovePoint(numberOfBins);
        }
        
	// ******************** Pi0 pp at 2760 published  PCM part ****************************
	TCanvas* canvasPi0ppPCM2760Published = new TCanvas("canvasPi0ppPCM2760Published","",200,10,1350,1350*1.15);  // gives the page size
	DrawGammaCanvasSettings( canvasPi0ppPCM2760Published,0.14, 0.02, 0.02, 0.09);
	canvasPi0ppPCM2760Published->SetLogy(1);
	canvasPi0ppPCM2760Published->SetLogx(1);
        TH2D *histo2DPi0ppPCM2760Published = new TH2D("histo2DPi0ppPCM2760Published", "histo2DPi0ppPCM2760Published", 20,0.3,10.,1000.,4e4,1e11);
	SetStyleHistoTH2ForGraphs(histo2DPi0ppPCM2760Published , "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )", 0.02,0.02, 0.02,0.02, 0.8,0.65, 510, 510);
 
	//         histo2DPi0ppPCM2760Published ->GetYaxis()->SetRangeUser(0.,1.02);
	histo2DPi0ppPCM2760Published ->Draw("copy");
	//	graphInvCrossSectionPi0PCMSys2760GeVYShifted->Draw("p,E2same");
	DrawGammaSetMarkerTGraphAsym(graphInvCrossSectionPi0PCMSys2760GeVYShifted, markerStyleCommmonSpectrum2760GeV,1.6*markerSizeCommonSpectrum, colorCommonSpectrumEta2760GeV, colorCommonSpectrumEta2760GeV, widthCommonSpectrumBoxes, kTRUE);

	graphInvCrossSectionPi0PCMStat2760GeVYShifted->Draw("p,same");
	graphInvCrossSectionPi0PCMSys2760GeVYShifted->Draw("E2,same");
	//	fitInvCrossSectionPi0PCM2760GeVStat->Draw("same");
	fitInvCrossSectionPi0PCM2760GeVSys->Draw("same");
	canvasPi0ppPCM2760Published->Update();
	canvasPi0ppPCM2760Published->SaveAs(Form("%s/%s_Pi0ppPCM2760Published.%s",outputDir.Data(), prefix2.Data(), suffix.Data()));


       //*********************** Eta/pi0 at 2.76 TeV ***************************************************
	
	TCanvas* canvasRatioEtaPi0DiffEnergies = new TCanvas("canvasRatioEtaPi0DiffEnergies","",200,10,1350,1350*1.15);  // gives the page size
        DrawGammaCanvasSettings( canvasRatioEtaPi0DiffEnergies,0.14, 0.02, 0.02, 0.09);
   
        TH2D *histo2DRatioEtaPi0DiffEnergies = new TH2D("histo2DRatioEtaPi0DiffEnergies", "histo2DRatioEtaPi0DiffEnergies", 20,0.,maxPtPCMEta2760GeV+10.,1000.,-0.4,1.2);
	SetStyleHistoTH2ForGraphs(histo2DRatioEtaPi0DiffEnergies, "#it{p}_{T} (GeV/#it{c})","#eta/#pi^{0}", 0.04,0.04, 0.04,0.04, 0.8,0.65, 510, 510);

	histo2DRatioEtaPi0DiffEnergies->GetXaxis()->SetRangeUser(0.,maxPtPCMEta2760GeV+1);
	histo2DRatioEtaPi0DiffEnergies->GetYaxis()->SetRangeUser(0.,1.02);
	histo2DRatioEtaPi0DiffEnergies->Draw("copy");
	DrawGammaSetMarkerTGraphAsym(graphRatioEtaPCMtoPi0Published2760GeVSys, markerStyleCommmonSpectrum2760GeV,1.6*markerSizeCommonSpectrum, colorCommonSpectrumEta2760GeV, colorCommonSpectrumEta2760GeV, widthCommonSpectrumBoxes, kTRUE);
	DrawGammaSetMarkerTGraphAsym(graphRatioEtaPCMtoPi0Published2760GeVStat, markerStyleCommmonSpectrum2760GeV,1.6*markerSizeCommonSpectrum, colorCommonSpectrumEta2760GeV, colorCommonSpectrumEta2760GeV, widthCommonSpectrumBoxes, kFALSE);
	DrawGammaSetMarkerTGraphAsym(graphEtaToPi0Comb7TeVPublished, markerStyleCommmonSpectrum7TeV,1.6*markerSizeCommonSpectrum, colorCommonSpectrumEta7TeV, colorCommonSpectrumEta7TeV, widthCommonSpectrumBoxes, kTRUE);
	DrawGammaSetMarkerTGraphAsym(graphEtaToPi0Comb7TeVPublishedStat, markerStyleCommmonSpectrum7TeV,1.6*markerSizeCommonSpectrum, colorCommonSpectrumEta7TeV, colorCommonSpectrumEta7TeV, widthCommonSpectrumBoxes, kTRUE);
	DrawGammaSetMarkerTGraphAsym(graphEtaToPi0Comb7TeVPublishedSys, markerStyleCommmonSpectrum7TeV,1.6*markerSizeCommonSpectrum, colorCommonSpectrumEta7TeV, colorCommonSpectrumEta7TeV, widthCommonSpectrumBoxes, kTRUE);
	graphRatioEtaPCMtoPi0Published2760GeVSys->Draw("p,E2same");
	graphRatioEtaPCMtoPi0Published2760GeVStat->Draw("p,same");

	 //graphEtaToPi0Comb7TeVPublished->Print();
	 //	 graphRatioEtaPCMtoPi0Published2760GeVSys->Draw("p,same");

	DrawGammaNLOTGraph(graphEtaToPi0NLOMuHalf2760GeV, widthLineNLO, styleLineNLOMuHalf, colorNLOPi02760GeVMuHalf);
	DrawGammaNLOTGraph(graphEtaToPi0NLOMuOne2760GeV, widthLineNLO, styleLineNLOMuOne, colorNLOPi02760GeVMuOne);
	DrawGammaNLOTGraph(graphEtaToPi0NLOMuTwo2760GeV, widthLineNLO, styleLineNLOMuTwo, colorNLOPi02760GeVMuTwo);

	TLegend* legendRatioEtaToPi0Theory2760GeV = new TLegend(0.5,0.16,0.977,0.34);
	legendRatioEtaToPi0Theory2760GeV->SetTextSize(0.03);			
	legendRatioEtaToPi0Theory2760GeV->SetFillColor(0);
	legendRatioEtaToPi0Theory2760GeV->SetBorderSize(0);
	legendRatioEtaToPi0Theory2760GeV->AddEntry(graphRatioEtaPCMtoPi0Published2760GeVSys,Form("PCM, #sqrt{#it{s}} = 2.76 TeV "),"fp");

	legendRatioEtaToPi0Theory2760GeV->AddEntry(graphEtaToPi0NLOMuHalf2760GeV,"NLO #mu = 0.5 #it{p}_{T}","l");
	legendRatioEtaToPi0Theory2760GeV->AddEntry(graphEtaToPi0NLOMuOne2760GeV,"NLO #mu = #it{p}_{T}","l");
	legendRatioEtaToPi0Theory2760GeV->AddEntry(graphEtaToPi0NLOMuTwo2760GeV,"NLO #mu = 2 #it{p}_{T}","l");
	legendRatioEtaToPi0Theory2760GeV->Draw();


	graphEtaToPi0NLOMuHalf2760GeV->Draw("same,c");
	graphEtaToPi0NLOMuOne2760GeV->Draw("same,c");
	graphEtaToPi0NLOMuTwo2760GeV->Draw("same,c");
	canvasRatioEtaPi0DiffEnergies->Update();
	canvasRatioEtaPi0DiffEnergies->SaveAs(Form("%s/%s_Pi02760PublishedEtaNewRatio.%s",outputDir.Data(), prefix2.Data(), suffix.Data()));


	histo2DRatioEtaPi0DiffEnergies->GetXaxis()->SetRangeUser(0.,maxPtPCMEta2760GeV+10.);
	histo2DRatioEtaPi0DiffEnergies->Draw("copy");
	graphRatioEtaPCMtoPi0Published2760GeVSys->Draw("p,E2same");
	graphRatioEtaPCMtoPi0Published2760GeVStat->Draw("p,same");
	graphEtaToPi0NLOMuHalf2760GeV->Draw("same,c");
	graphEtaToPi0NLOMuOne2760GeV->Draw("same,c");
	graphEtaToPi0NLOMuTwo2760GeV->Draw("same,c");


	legendRatioEtaToPi0Theory2760GeV->AddEntry(graphEtaToPi0Comb7TeVPublishedSys,Form("Phys. Lett. B 717 (2012) 162 , #sqrt{#it{s}} = 7 TeV "),"fp");
	legendRatioEtaToPi0Theory2760GeV->Draw();

	graphEtaToPi0Comb7TeVPublishedStat->Draw("p,same");
	graphEtaToPi0Comb7TeVPublishedSys->Draw("p,E2same");
	canvasRatioEtaPi0DiffEnergies->Update();
	canvasRatioEtaPi0DiffEnergies->SaveAs(Form("%s/%s_Pi02760PublishedEtaNewRatio_Compared7TeV.%s",outputDir.Data(), prefix2.Data(), suffix.Data()));

	//-------------------Comparison to Theory-2760GeV------------------------------/
	/*	
	TCanvas* canvasRatioEtaPi0NLO2760GeV = new TCanvas("canvasRatioEtaPi0NLO2760GeV","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings( canvasRatioEtaPi0NLO2760GeV,0.09, 0.01, 0.015, 0.115);
	
	TH2D *histo2DRatioEtaPi0Theory2760GeV = new TH2D("histo2DRatioEtaPi0Theory2760GeV", "histo2DRatioEtaPi0Theory2760GeV", 20,0.,maxPtEta2760GeV,1000.,-0.4,1.2);
	SetStyleHistoTH2ForGraphs(histo2DRatioEtaPi0Theory2760GeV, "#it{p}_{T} (GeV/#it{c})","#eta/#pi^{0}", 0.046,0.058, 0.046,0.058, 0.8,0.65, 510, 510);
	histo2DRatioEtaPi0Theory2760GeV->GetYaxis()->SetRangeUser(0.,1.02);
	histo2DRatioEtaPi0Theory2760GeV->Draw();


	DrawGammaNLOTGraph(graphEtaToPi0NLOMuHalf2760GeV, widthLineNLO, styleLineNLOMuHalf, colorNLOPi02760GeVMuHalf);
	DrawGammaNLOTGraph(graphEtaToPi0NLOMuOne2760GeV, widthLineNLO, styleLineNLOMuOne, colorNLOPi02760GeVMuOne);
	DrawGammaNLOTGraph(graphEtaToPi0NLOMuTwo2760GeV, widthLineNLO, styleLineNLOMuTwo, colorNLOPi02760GeVMuTwo);

	graphEtaToPi0NLOMuHalf2760GeV->Draw("same,c");
	graphEtaToPi0NLOMuOne2760GeV->Draw("same,c");
	graphEtaToPi0NLOMuTwo2760GeV->Draw("same,c");

 	histo2DRatioEtaPi0Theory2760GeV->Draw("AXIS,same");
	
	TLegend* legendRatioEtaToPi0Theory2760GeV = new TLegend(0.5,0.16,0.977,0.34);
	legendRatioEtaToPi0Theory2760GeV->SetTextSize(0.036);			
	legendRatioEtaToPi0Theory2760GeV->SetFillColor(0);
	legendRatioEtaToPi0Theory2760GeV->SetBorderSize(0);
	legendRatioEtaToPi0Theory2760GeV->AddEntry(graphRatioEtaPi0ComplErr2760GeVNewXError,Form("stat + syst, PCM, #sqrt{#it{s}} = 2.76 TeV "),"fp");
	legendRatioEtaToPi0Theory2760GeV->AddEntry(graphEtaToPi0NLOMuHalf2760GeV,"NLO #mu = 0.5 #it{p}_{T}","l");
	legendRatioEtaToPi0Theory2760GeV->AddEntry(graphEtaToPi0NLOMuOne2760GeV,"NLO #mu = #it{p}_{T}","l");
	legendRatioEtaToPi0Theory2760GeV->AddEntry(graphEtaToPi0NLOMuTwo2760GeV,"NLO #mu = 2 #it{p}_{T}","l");
	legendRatioEtaToPi0Theory2760GeV->Draw();
	
	canvasRatioEtaPi0NLO2760GeV->Update();
	canvasRatioEtaPi0NLO2760GeV->SaveAs(Form("%s/%s_Pi0EtaRatioTheory2760GeV_Paper.%s",outputDir.Data(), prefix2.Data(), suffix.Data()));
	*/
	
	/*
	TCanvas* canvasRatioEtaPi03 = new TCanvas("canvasRatioEtaPi03","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings( canvasRatioEtaPi03,0.09, 0.01, 0.015, 0.115);
	canvasRatioEtaPi03->SetLogy(1);
	
	TH2D *histo2DRatioEtaPi0 = new TH2D("histo2DRatioEtaPi0", "histo2DRatioEtaPi0", 20,0.,maxPtEta2760GeV,1000.,-0.4,6.0);
	SetStyleHistoTH2ForGraphs(histo2DRatioEtaPi0, "#it{p}_{T} (GeV/#it{c})","#eta/#pi^{0}", 0.046,0.058, 0.046,0.058, 0.8,0.65, 510, 510);
	histo2DRatioEtaPi0->GetYaxis()->SetRangeUser(0.03,2.8);
	histo2DRatioEtaPi0->Draw();

	cout<<"printout pp2760 eta data"<<endl;
	graphRatioEtaPi0ComplErr2760GeVNewXError->Print();
	graphRatioEtaPi0ComplErr2760GeVNewXError->Draw("p,E2same");  

 
	graphDonaldson100GeV->Draw("p,same");
	graphDonaldson200GeV->Draw("p,same");
	graphAntille87pp->Draw("p,same");
	graphAguilar400GeV->Draw("p,same");
	graphKourkou79pp->Draw("p,same");
	graphApana530GeV->Draw("p,same");
	graphKourkou79pp52->Draw("p,same");
	graphAkesson53GeVpp->Draw("p,same");
	graphKourkou79pp62->Draw("p,same");
	graphPhenix200GeV->Draw("p,same");


	graphRatioEtaPi0ComplErr2760GeVNoError->Draw("p,same");  

		
	TLegend* legendRatio3 = new TLegend(0.22,0.14,0.96,0.28);
	legendRatio3->SetTextSize(0.027);			
	legendRatio3->SetFillColor(0);
	legendRatio3->SetBorderSize(0);
	legendRatio3->SetNColumns(2);
//legendRatio3->AddEntry(graphCombinedEtaToPi0NewXError,Form("p+p ALICE (#sqrt{#it{s}} = 7 TeV, PLB 717 (2012) 162-172)"),"fp");
	legendRatio3->AddEntry(graphRatioEtaPi0ComplErr2760GeVNewXError,Form("p+p ALICE (#sqrt{#it{s}} = 2.76 TeV)"),"fp");
//	legendRatio3->AddEntry(graphRatioEtaPi0ComplErr900GeVNewXError,Form("p+p ALICE (#sqrt{#it{s}} = 0.9 TeV)"),"fp");
	legendRatio3->AddEntry(graphDonaldson100GeV,graphDonaldson100GeV->GetTitle(),"p");
	legendRatio3->AddEntry(graphDonaldson200GeV,graphDonaldson200GeV->GetTitle(),"p");
	legendRatio3->AddEntry(graphAntille87pp,graphAntille87pp->GetTitle(),"p");
	legendRatio3->AddEntry(graphAguilar400GeV,graphAguilar400GeV->GetTitle(),"p");
	legendRatio3->AddEntry(graphKourkou79pp,graphKourkou79pp->GetTitle(),"p");
	legendRatio3->AddEntry(graphApana530GeV,graphApana530GeV->GetTitle(),"p");
	legendRatio3->AddEntry(graphKourkou79pp52,graphKourkou79pp52->GetTitle(),"p");
	legendRatio3->AddEntry(graphAkesson53GeVpp,graphAkesson53GeVpp->GetTitle(),"p");
	legendRatio3->AddEntry(graphKourkou79pp62,graphKourkou79pp62->GetTitle(),"p");
	legendRatio3->AddEntry(graphPhenix200GeV,graphPhenix200GeV->GetTitle(),"p");	
	legendRatio3->Draw();
	
	canvasRatioEtaPi03->Update();
	canvasRatioEtaPi03->SaveAs(Form("%s/%s_Pi0EtaRatioWorld_Paper.%s",outputDir.Data(), prefix2.Data(), suffix.Data()));
	
	
	//	if(!thesis)DrawAliceLogoCombinedPreliminary(pictDrawingCoordinatesPi0Eta[0], pictDrawingCoordinatesPi0Eta[1], pictDrawingCoordinatesPi0Eta[2], pictDrawingCoordinatesPi0Eta[3], pictDrawingCoordinatesPi0Eta[4], pictDrawingCoordinatesPi0Eta[5], pictDrawingCoordinatesPi0Eta[6], pictDrawingCoordinatesPi0Eta[7], pictDrawingCoordinates[8],collisionSystemCombinedReallyAll, pictDrawingOptions[1], pictDrawingOptions[2], pictDrawingOptions[3],1350,900);
	
	canvasRatioEtaPi03->Update();
	canvasRatioEtaPi03->SaveAs(Form("%s/%s_Pi0EtaRatioWorld.%s",outputDir.Data(), prefix2.Data(), suffix.Data()));
	
	delete canvasRatioEtaPi03;
	delete legendRatio3;
	*/

	//-----------------------------Ratios to NLO
	
	graphRatioCombNLOEta2760GeVMuHalf= (TGraph*)graphNLOMuHalfEta2760GeV->Clone();
	graphRatioCombNLOEta2760GeVMuOne= (TGraph*)graphNLOMuOneEta2760GeV->Clone();
	graphRatioCombNLOEta2760GeVMuTwo= (TGraph*)graphNLOMuTwoEta2760GeV->Clone();
	// graphRatioCombNLOEta2760GeVMuHalf->Print();
	// graphRatioCombNLOEta2760GeVMuOne->Print();
	// graphRatioCombNLOEta2760GeVMuTwo->Print();



	histoRatioCombConvEta2760GeV = (TH1D*) histoPCMInvCrossSectionEta2760GeV->Clone();		
	histoRatioCombConvEta2760GeV = CalculateHistoRatioToFit (histoRatioCombConvEta2760GeV, fitInvCrossSectionEta2760GeV); 
	graphRatioCombNLOEta2760GeVMuHalf = CalculateGraphRatioToFit (graphRatioCombNLOEta2760GeVMuHalf, fitInvCrossSectionEta2760GeV); 
	graphRatioCombNLOEta2760GeVMuOne = CalculateGraphRatioToFit (graphRatioCombNLOEta2760GeVMuOne, fitInvCrossSectionEta2760GeV); 
	graphRatioCombNLOEta2760GeVMuTwo = CalculateGraphRatioToFit (graphRatioCombNLOEta2760GeVMuTwo, fitInvCrossSectionEta2760GeV); 
	
	graphRatioCombCombFitEta2760GeV = (TGraphAsymmErrors*)graphInvCrossSectionEtaComb2760GeV->Clone();
	graphRatioCombCombFitEta2760GeV = CalculateGraphErrRatioToFit(graphRatioCombCombFitEta2760GeV, fitInvCrossSectionEta2760GeV); 
	graphRatioCombCombFitEta2760GeVStat = (TGraphAsymmErrors*)graphInvCrossSectionEtaComb2760GeVStatErr->Clone();
	graphRatioCombCombFitEta2760GeVStat = CalculateGraphErrRatioToFit(graphRatioCombCombFitEta2760GeVStat, fitInvCrossSectionEta2760GeV); 
	graphRatioCombCombFitEta2760GeVSys = (TGraphAsymmErrors*)graphInvCrossSectionEtaComb2760GeVSysErr->Clone();
	graphRatioCombCombFitEta2760GeVSys = CalculateGraphErrRatioToFit(graphRatioCombCombFitEta2760GeVSys, fitInvCrossSectionEta2760GeV); 
	histoFitInvCrossSectionEta2760GeV = (TH1D*)fitInvCrossSectionEta2760GeV->GetHistogram();

	TGraphAsymmErrors *graphRatioPrelimQM2011CombFitEta2760GeV = (TGraphAsymmErrors*)graphInvCrossSectionEtaComb2760GeV_PrelimQM2011->Clone();
	graphRatioPrelimQM2011CombFitEta2760GeV = CalculateGraphErrRatioToFit(graphRatioPrelimQM2011CombFitEta2760GeV, fitInvCrossSectionEta2760GeV); 

	TGraphAsymmErrors *graphRatioPrelimQM2011CombFitEta2760GeVStatErr = (TGraphAsymmErrors*)graphInvCrossSectionEtaComb2760GeVStatErr_PrelimQM2011->Clone();
	graphRatioPrelimQM2011CombFitEta2760GeVStatErr = CalculateGraphErrRatioToFit(graphRatioPrelimQM2011CombFitEta2760GeVStatErr, fitInvCrossSectionEta2760GeV); 

	TGraphAsymmErrors *graphRatioPrelimQM2011CombFitEta2760GeVSysErr = (TGraphAsymmErrors*)graphInvCrossSectionEtaComb2760GeVSysErr_PrelimQM2011->Clone();
	graphRatioPrelimQM2011CombFitEta2760GeVSysErr = CalculateGraphErrRatioToFit(graphRatioPrelimQM2011CombFitEta2760GeVSysErr, fitInvCrossSectionEta2760GeV); 


	// for(Int_t i=0;i< 200;i++){
	//   cout<< i<< " "<< histoFitInvCrossSectionEta2760GeV->GetBinCenter(i)<< " " <<histoFitInvCrossSectionEta2760GeV->GetBinContent(i)<< endl;
  	// }
	// histoFitInvCrossSectionEta2760GeV->SetAxisRange(0.6,19.5);
	if (suffix.CompareTo("eps")==0){
		widthLinesBoxes				= 1.4;
		widthCommonFit					= 2.;
		widthStatErrBars				= 1.5;
		widthCommonErrors				= 1.1;
		widthCommonSpectrumBoxes			= 0.99;
	} else {
		widthLinesBoxes				= 2.3;
		widthCommonFit					= 2.6;
		widthStatErrBars				= 2.6;
		widthCommonErrors				= 2.;
		widthCommonSpectrumBoxes			= 2.3;
	}


	//***************************************************************************************************************
	//************************ Eta Combined + NLO *******************************************************************
	//***************************************************************************************************************
	
	TCanvas* canvasInvXSectionEtaALLEnergies = new TCanvas("canvasInvXSectionEtaALLEnergies","",200,10,1200,2000);  // gives the page size
	DrawGammaCanvasSettings( canvasInvXSectionEtaALLEnergies,  0.15, 0.02, 0.03, 0.06);
	
	TPad* padComparisonXSectionEtaALLEnergies = new TPad("padComparisonXSectionEtaALLEnergies", "", 0., 0.42, 1., 1.,-1, -1, -2);
	DrawGammaPadSettings( padComparisonXSectionEtaALLEnergies, 0.15, 0.02, 0.03, 0.);
	padComparisonXSectionEtaALLEnergies->Draw();

	TPad* padXSectionEtaALLEnergiesRatioEta2760GeV = new TPad("padXSectionEtaALLEnergiesRatioEta2760GeV", "", 0., 0.18, 1., 0.30,-1, -1, -2);
	DrawGammaPadSettings( padXSectionEtaALLEnergiesRatioEta2760GeV, 0.15, 0.02, 0., 0.);
	padXSectionEtaALLEnergiesRatioEta2760GeV->Draw();
	padComparisonXSectionEtaALLEnergies->cd();
	padComparisonXSectionEtaALLEnergies->SetLogy();		
	padComparisonXSectionEtaALLEnergies->SetLogx();		
	
	//-------------- Plotting ------------------------------------------------------
	TH2F * histo2DInvXSectionEtaALLEnergies = new TH2F("histo2DInvXSectionEtaALLEnergies","histo2DInvXSectionEtaALLEnergies",1000,0.4,25.,1000,20,1e10);
	SetStyleHistoTH2ForGraphs(histo2DInvXSectionEtaALLEnergies, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )", 0.032,0.04, 0.04,0.035, 1,1.55);
	histo2DInvXSectionEtaALLEnergies->DrawCopy(); 
	
	//	graphInvCrossSectionEtaComb2760GeV = ScaleGraph(graphInvCrossSectionEtaComb2760GeV,1e-1);
        DrawGammaSetMarkerTGraphAsym(graphInvCrossSectionEtaComb2760GeV, markerStyleCommmonSpectrum2760GeV,1.6*markerSizeCommonSpectrum, colorCommonSpectrumEta2760GeV, colorCommonSpectrumEta2760GeV, widthCommonSpectrumBoxes, kTRUE);
        graphInvCrossSectionEtaComb2760GeV->SetLineWidth(widthCommonErrors);
	// //graphInvCrossSectionEtaComb2760GeV->Draw("same,p");
        graphInvCrossSectionEtaComb2760GeV->Draw("p,E2same");

	//histoFitInvCrossSectionEta2760GeV->Scale(1e-1);
	SetStyleHisto(histoFitInvCrossSectionEta2760GeV, widthCommonFit, styleFitCommonSpectrum, colorCommonSpectrumPi02760GeV);
	//histoFitInvCrossSectionEta2760GeV->Draw("same,c");

	//	graphNLOMuHalfEta2760GeV= ScaleGraph(graphNLOMuHalfEta2760GeV,1e-1);
	DrawGammaNLOTGraph( graphNLOMuHalfEta2760GeV, widthCommonFit, styleLineNLOMuHalf, colorNLOPi02760GeVMuHalf);
	graphNLOMuHalfEta2760GeV->Draw("same,c");
	//graphNLOMuOneEta2760GeV= ScaleGraph(graphNLOMuOneEta2760GeV,1e-1);
	DrawGammaNLOTGraph( graphNLOMuOneEta2760GeV, widthCommonFit, styleLineNLOMuOne, colorNLOPi02760GeVMuOne);
	graphNLOMuOneEta2760GeV->Draw("same,c");
	//graphNLOMuTwoEta2760GeV = ScaleGraph(graphNLOMuTwoEta2760GeV,1e-1);
	DrawGammaNLOTGraph( graphNLOMuTwoEta2760GeV, widthCommonFit, styleLineNLOMuTwo, colorNLOPi02760GeVMuTwo);
	graphNLOMuTwoEta2760GeV->Draw("same,c");

	//graphEMCALXSectionEta2760GeVYShifted->Draw("same,c");
	//graphPCMXSectionEta2760GeVYShifted->Draw("same,c");
	graphEMCALXSectionEta2760GeVYShifted->Draw("p,same");
	graphPCMXSectionEta2760GeVYShifted->Draw("p,same");

	graphEMCALXSectionEta2760GeVUnShifted->Draw("same,c");
	graphEMCALSysErrInvCrossSectionEta2760GeVforComp->Draw("same,c");
	graphPCMXSectionEta2760GeVUnShifted->Draw("same,c");


	graphEMCALSysErrInvCrossSectionEta2760GeVforComp->Draw("same,p");
	graphEMCALXSectionEta2760GeVUnShifted->Draw("p,same");
	graphPCMXSectionEta2760GeVUnShifted->Draw("p,same");

	graphInvCrossSectionEtaComb2760GeV->Draw("p,E2same");
 	TLatex *labelScalingEta2760GeVAllEnergiesSp = new TLatex(0.46,3E8,"x 10^{-1}");
 	SetStyleTLatex( labelScalingEta2760GeVAllEnergiesSp, 0.028,4,colorCommonSpectrumPi02760GeV,kFALSE);
	// 	labelScalingEta2760GeVAllEnergiesSp->Draw();

	TLatex *labelScalingEta2760GeVALLEnergies = new TLatex(0.46,3E8,"x 10^{-1}");
	SetStyleTLatex( labelScalingEta2760GeVALLEnergies, 0.025,4,histoFitInvCrossSectionEta2760GeV->GetLineColor(),62,kFALSE);
	//labelScalingEta2760GeVALLEnergies->Draw();
	// DrawNormalizationErrorText(normalizationInvX1MesonALLEn[0],normalizationInvX1MesonALLEn[1],normalizationInvX1MesonALLEn[2],
	// 					  normalizationInvX1MesonALLEn[3],normalizationInvX1MesonALLEn[4],"all"); 
			
	/*
	TPad* padXSectionEtaALLEnergiesLegend = new TPad("padXSectionEtaALLEnergiesLegend", "", 0.17, 0.005, 0.95, 0.21,-1, -1, -2); 
	DrawGammaPadSettings( padXSectionEtaALLEnergiesLegend, 0., 0., 0., 0.);
	padXSectionEtaALLEnergiesLegend->Draw();
	padXSectionEtaALLEnergiesLegend->cd();
	*/
	//*************** third Column **********************************************************
	// TLatex *textEta2760GeVNLOEtaALLEnergies = new TLatex(columnsNLOLegendPrelAndFinal[2],rowsNLOLegendPrelAndFinal[0],"#eta, #sqrt{#it{s}} = 2.76 TeV");
	// SetStyleTLatex( textEta2760GeVNLOEtaALLEnergies, textSizeTopLablesPrelAndFinal,4);
        // textEta2760GeVNLOEtaALLEnergies->Draw();

	//************************************************* End Legend ***************************************************
	padXSectionEtaALLEnergiesRatioEta2760GeV->cd();
	padXSectionEtaALLEnergiesRatioEta2760GeV->SetLogx();
	TH2F * ratio2DInvXSectionEtaALLEnergiesEta2760GeV;
	ratio2DInvXSectionEtaALLEnergiesEta2760GeV = new TH2F("ratio2DInvXSectionEtaALLEnergiesEta2760GeV","ratio2DInvXSectionEtaALLEnergiesEta2760GeV",1000,0.23,30.,1000,0.4,3.55);
	SetStyleHistoTH2ForGraphs(ratio2DInvXSectionEtaALLEnergiesEta2760GeV, "#it{p}_{T} (GeV/#it{c})","#frac{NLO}{fit}", 0.13,0.13, 0.175,0.18, 1,0.32, 512, 505);
	ratio2DInvXSectionEtaALLEnergiesEta2760GeV->DrawCopy(); 



	DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitEta2760GeV, markerStyleCommmonSpectrum2760GeV,1.6*markerSizeCommonSpectrum, colorCommonSpectrumPi02760GeV, colorCommonSpectrumPi02760GeV, widthCommonSpectrumBoxes, kTRUE);
	graphRatioCombCombFitEta2760GeV->SetLineWidth(widthCommonErrors);
	graphRatioCombCombFitEta2760GeV->Draw("same,p");
       	graphRatioCombCombFitEta2760GeV->Draw("p,E2same");
	// cout<< "Ratio NLO to fit; helf , one two"<<endl;

	// graphRatioCombNLOEta2760GeVMuHalf->Print();
        // cout<< endl;
 	// graphRatioCombNLOEta2760GeVMuOne->Print();
        // cout<< endl;
	// graphRatioCombNLOEta2760GeVMuTwo->Print();
        // cout<< endl;
        // cout<< endl;
	DrawGammaNLOTGraph( graphRatioCombNLOEta2760GeVMuHalf, widthCommonFit, styleLineNLOMuHalf, colorNLOPi02760GeVMuHalf);
	graphRatioCombNLOEta2760GeVMuHalf->Draw("same,c");
	DrawGammaNLOTGraph( graphRatioCombNLOEta2760GeVMuOne, widthCommonFit, styleLineNLOMuOne, colorNLOPi02760GeVMuOne);
	graphRatioCombNLOEta2760GeVMuOne->Draw("same,c");
	DrawGammaNLOTGraph( graphRatioCombNLOEta2760GeVMuTwo, widthCommonFit, styleLineNLOMuTwo, colorNLOPi02760GeVMuTwo);
	graphRatioCombNLOEta2760GeVMuTwo->Draw("same,c");
	
	// boxErrorSigmaPi02760GeVRatio->Draw();
	
	TLatex *labelRatioNLOEta2760GeV = new TLatex(0.18,0.75,"#eta, #sqrt{#it{s}} = 2.76 TeV");
	SetStyleTLatex( labelRatioNLOEta2760GeV, 0.17,4);
	labelRatioNLOEta2760GeV->Draw();
	
	DrawGammaLines(0., 20.,1., 1.,0.1);

	canvasInvXSectionEtaALLEnergies->Update();
	canvasInvXSectionEtaALLEnergies->Print(Form("%s/Eta_InvXSectionALLEnergies_Paper.%s",outputDir.Data(),suffix.Data()));

	
	padComparisonXSectionEtaALLEnergies->cd();
	//	if(!thesis)DrawAliceLogoPi0WithPHOSPrelim(pictDrawingCoordinatesCombineMeasCross[0], pictDrawingCoordinatesCombineMeasCross[1], pictDrawingCoordinatesCombineMeasCross[2], pictDrawingCoordinatesCombineMeasCross[3], pictDrawingCoordinatesCombineMeasCross[4], pictDrawingCoordinatesCombineMeasCross[5], pictDrawingCoordinatesCombineMeasCross[6], pictDrawingCoordinatesCombineMeasCross[7], pictDrawingCoordinates[8],collisionSystemCombinedReallyAll, pictDrawingOptions[1], pictDrawingOptions[2], kFALSE,1200,1220);

	canvasInvXSectionEtaALLEnergies->Update();
	canvasInvXSectionEtaALLEnergies->Print(Form("%s/Eta_InvXSectionALLEnergies.%s",outputDir.Data(),suffix.Data()));



	//*****************************************************************************************************************
	// Plotting only ratio----	
	//*****************************************************************************************************************
	TCanvas* canvasInvXSectionOnlyRatioEta = new TCanvas("canvasInvXSectionOnlyRatioEta","",200,10,1200,600);  // gives the page size
	DrawGammaCanvasSettings( canvasInvXSectionOnlyRatioEta,  0.075, 0.05, 0.05, 0.12);
	canvasInvXSectionOnlyRatioEta->SetLogx(1);
	TH2F * ratio2DInvXSectionEtaOnlyEta2760GeV;
	ratio2DInvXSectionEtaOnlyEta2760GeV = new TH2F("ratio2DInvXSectionEtaOnlyEta2760GeV","ratio2DInvXSectionEtaOnlyEta2760GeV",1000,0.4,30.,1000,0.4,3.55);
	SetStyleHistoTH2ForGraphs(ratio2DInvXSectionEtaOnlyEta2760GeV, "#it{p}_{T} (GeV/#it{c})","#frac{Data,NLO}{fit}", 0.04,0.04, 0.03,0.03, 1,0.32, 512, 505);
	ratio2DInvXSectionEtaOnlyEta2760GeV->DrawCopy(); 
	DrawGammaLines(0., 30.,1., 1.,0.1);
	DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitEta2760GeV, markerStyleCommmonSpectrum2760GeV,1.6*markerSizeCommonSpectrum, colorCommonSpectrumPi02760GeV, colorCommonSpectrumPi02760GeV, widthCommonSpectrumBoxes, kTRUE);
	graphRatioCombCombFitEta2760GeV->SetLineWidth(widthCommonErrors);
	graphRatioCombCombFitEta2760GeV->Draw("same,p");
       	graphRatioCombCombFitEta2760GeV->Draw("p,E2same");
	
	DrawGammaNLOTGraph( graphRatioCombNLOEta2760GeVMuHalf, widthCommonFit, styleLineNLOMuHalf, colorNLOPi02760GeVMuHalf);
	graphRatioCombNLOEta2760GeVMuHalf->Draw("same,c");
	DrawGammaNLOTGraph( graphRatioCombNLOEta2760GeVMuOne, widthCommonFit, styleLineNLOMuOne, colorNLOPi02760GeVMuOne);
	graphRatioCombNLOEta2760GeVMuOne->Draw("same,c");
	DrawGammaNLOTGraph( graphRatioCombNLOEta2760GeVMuTwo, widthCommonFit, styleLineNLOMuTwo, colorNLOPi02760GeVMuTwo);
	graphRatioCombNLOEta2760GeVMuTwo->Draw("same,c");

	TLegend* legendOnlyRatioEta2760GeV = new TLegend(0.65,0.72,0.9,0.92);
	legendOnlyRatioEta2760GeV->SetTextSize(0.03);			
	legendOnlyRatioEta2760GeV->SetFillColor(0);
	legendOnlyRatioEta2760GeV->SetBorderSize(0);
	legendOnlyRatioEta2760GeV->AddEntry(graphRatioCombCombFitEta2760GeV,Form(" #sqrt{#it{s}} = 2.76 TeV "),"fp");

	legendOnlyRatioEta2760GeV->AddEntry(graphRatioCombNLOEta2760GeVMuHalf,"NLO #mu = 0.5 #it{p}_{T}","l");
	legendOnlyRatioEta2760GeV->AddEntry(graphRatioCombNLOEta2760GeVMuOne,"NLO #mu = #it{p}_{T}","l");
	legendOnlyRatioEta2760GeV->AddEntry(graphRatioCombNLOEta2760GeVMuTwo,"NLO #mu = 2 #it{p}_{T}","l");
	legendOnlyRatioEta2760GeV->Draw();

	canvasInvXSectionOnlyRatioEta->Update();
	canvasInvXSectionOnlyRatioEta->Print(Form("%s/Eta_InvXSectionOnlyRatioEta.%s",outputDir.Data(),suffix.Data()));


	//*****************************************************************************************************************
	//************************ Eta Meson Combined + compared to Prelim QM2011******************************************************
	//**************************************************************************

	TCanvas* canvasInvXSectionOnlyRatioEtaPrelimQM2011 = new TCanvas("canvasInvXSectionOnlyRatioEtaPrelimQM2011","",200,10,1200,600);  // gives the page size
	DrawGammaCanvasSettings( canvasInvXSectionOnlyRatioEtaPrelimQM2011 ,  0.075, 0.05, 0.05, 0.12);
	canvasInvXSectionOnlyRatioEtaPrelimQM2011 ->SetLogx(1);

	ratio2DInvXSectionEtaOnlyEta2760GeV->DrawCopy(); 
	DrawGammaLines(0., 30.,1., 1.,0.1);

	graphRatioCombCombFitEta2760GeV->SetLineWidth(widthCommonErrors);
	graphRatioCombCombFitEta2760GeV->Draw("same,p");
       	graphRatioCombCombFitEta2760GeV->Draw("p,E2same");

	DrawGammaSetMarkerTGraphAsym(graphRatioPrelimQM2011CombFitEta2760GeVStatErr, markerStyleCommmonSpectrum2760GeV,1.6*markerSizeCommonSpectrum, kBlue, kBlue, widthCommonSpectrumBoxes, kFALSE);
	graphRatioPrelimQM2011CombFitEta2760GeVStatErr->Draw("p,same");
	DrawGammaSetMarkerTGraphAsym(graphRatioPrelimQM2011CombFitEta2760GeVSysErr, markerStyleCommmonSpectrum2760GeV,1.6*markerSizeCommonSpectrum, kBlue, kBlue, widthCommonSpectrumBoxes, kTRUE);
	graphRatioPrelimQM2011CombFitEta2760GeVSysErr->Draw("p,E2same");

	TLegend* legendOnlyRatioEta2760GeVPrelimQM2011 = new TLegend(0.65,0.72,0.9,0.92);
	legendOnlyRatioEta2760GeVPrelimQM2011->SetTextSize(0.03);			
	legendOnlyRatioEta2760GeVPrelimQM2011->SetFillColor(0);
	legendOnlyRatioEta2760GeVPrelimQM2011->SetBorderSize(0);
	legendOnlyRatioEta2760GeVPrelimQM2011->AddEntry(graphRatioCombCombFitEta2760GeV,Form(" #sqrt{#it{s}} = 2.76 TeV "),"fp");
	legendOnlyRatioEta2760GeVPrelimQM2011->AddEntry(graphRatioPrelimQM2011CombFitEta2760GeVSysErr,Form("PCM- Prelim QM2011, #sqrt{#it{s}} = 2.76 TeV "),"fp");
	legendOnlyRatioEta2760GeVPrelimQM2011->Draw();
	canvasInvXSectionOnlyRatioEtaPrelimQM2011->Update();
	canvasInvXSectionOnlyRatioEtaPrelimQM2011->Print(Form("%s/Eta_InvXSectionOnlyRatioEtaPrelimQM2011.%s",outputDir.Data(),suffix.Data()));




	//*****************************************************************************************************************
	//************************ Eta Meson Combined + NLO just spectrum******************************************************
	//*****************************************************************************************************************
	
	TCanvas* canvasInvXSectionNLOOnlySpectraEta = new TCanvas("canvasInvXSectionNLOOnlySpectraEta","",200,10,1200,1200);  // gives the page size
	DrawGammaCanvasSettings( canvasInvXSectionNLOOnlySpectraEta,  0.15, 0.02, 0.03, 0.06);
	
	TPad* padComparisonXSectionNLOOnlySpectraEta = new TPad("padComparisonXSectionNLOOnlySpectraEta", "", 0., 0., 1., 0.8,-1, -1, -2);
	DrawGammaPadSettings( padComparisonXSectionNLOOnlySpectraEta, 0.15, 0.02, 0.02, 0.09);
	padComparisonXSectionNLOOnlySpectraEta->Draw();
	
	TPad* padXSectionNLOLegendOnlySpectraEta = new TPad("padXSectionNLOLegendOnlySpectraEta", "", 0.15, 0.8, 0.98, 1.,-1, -1, -2);
	DrawGammaPadSettings( padXSectionNLOLegendOnlySpectraEta, 0.15, 0.02, 0.03, 0.);
	padXSectionNLOLegendOnlySpectraEta->Draw();
	
	padComparisonXSectionNLOOnlySpectraEta->cd();
	padComparisonXSectionNLOOnlySpectraEta->SetLogy();		
	padComparisonXSectionNLOOnlySpectraEta->SetLogx();		
	
	//-------------- Plotting ------------------------------------------------------
	TH2F * histo2DInvXSectionNLOOnlySpectraEta = new TH2F("histo2DInvXSectionNLOOnlySpectraEta","histo2DInvXSectionNLOOnlySpectraEta",1000,0.4,25.,1000,20.,1e10);
	SetStyleHistoTH2ForGraphs(histo2DInvXSectionNLOOnlySpectraEta, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )", 0.032,0.04, 0.04,0.025, 1,1.55);
	histo2DInvXSectionNLOOnlySpectraEta->DrawCopy(); 
	
	// graphInvCrossSectionEtaComb900GeV->Draw("p,E2same");
	// graphInvCrossSectionEtaComb7TeV->Draw("p,E2same");

	// fitInvCrossSectionEta7TeV->Draw("same");

       	graphInvCrossSectionEtaComb2760GeV->Draw("p,E2same");
	//histoFitInvCrossSectionEta2760GeV ->Draw("same,c");

	
	graphEMCALXSectionEta2760GeVUnShifted->Draw("same,c");
	graphPCMXSectionEta2760GeVUnShifted->Draw("same,c");
	// graphInvCrossSectionEtaComb2760GeV_PrelimQM2011 ->Draw("same,p");
	// graphInvCrossSectionEtaComb2760GeVStatErr_PrelimQM2011 ->Draw("same,p");
	// graphInvCrossSectionEtaComb2760GeVSysErr_PrelimQM2011 ->Draw("p,E2same");

	graphEMCALXSectionEta2760GeVUnShifted->Draw("same,p");
	graphPCMXSectionEta2760GeVUnShifted->Draw("same,p");

	graphInvCrossSectionEtaComb2760GeV->Draw("p,E2same");
	
	graphNLOMuHalfEta2760GeV->Draw("same,c");
	graphNLOMuOneEta2760GeV->Draw("same,c");
	graphNLOMuTwoEta2760GeV->Draw("same,c");

//labelScalingEta7TeVALLEnergies->Draw();
//	labelScalingEta900GeVALLEnergies->Draw();
//	labelScalingEta2760GeVALLEnergies->Draw();

	DrawNormalizationErrorText(normalizationInvXOnlySpec[0]+0.1,normalizationInvXOnlySpec[1],normalizationInvXOnlySpec[2],
						  normalizationInvXOnlySpec[3],normalizationInvXOnlySpec[4],"all"); 
	
	//************************************************* Begin Legend ***************************************************

	// padXSectionNLOLegendOnlySpectraEta->cd();
	// DrawGammaPadSettings( padXSectionNLOLegendOnlySpectraEta, 0., 0., 0., 0.);
	// padXSectionNLOLegendOnlySpectraEta->SetBorderMode(-1);
	// padXSectionNLOLegendOnlySpectraEta->SetBorderSize(3);
	// padXSectionNLOLegendOnlySpectraEta->Draw();
	// padXSectionNLOLegendOnlySpectraEta->cd();

	//*************** first Column **********************************************************
	// textSpectrumALLEnergies->Draw();
	// textFitCombALLEnergies->Draw();
	// textNLOMuHalfALLEnergies->Draw();
	// textNLOMuOneALLEnergies->Draw();
	// textNLOMuTwoALLEnergies->Draw();
	
	//*************** second Column **********************************************************
	// textEta7TeVNLOEtaALLEnergies->Draw();
	// textPi07TeVNLOsysALLEnergies->Draw();
	// boxCombinedPi07TeVALLEnergies->Draw("l");
	// markerCombinedPi07TeVALLEnergies->DrawMarker(columnsNLOLegendPrelAndFinal[1]+offsetNLOLegendPrelAndFinalMarkerX,rowsNLOLegendPrelAndFinal[2]+offsetNLOLegendPrelAndFinalMarkerY);
	// lineFit7TeVNLOALLEnergies->Draw("same");
	// lineNLOPi07TeVMuHalfALLEnergies->Draw("same");
	// lineNLOPi07TeVMuOneALLEnergies->Draw("same");
	// lineNLOPi07TeVMuTwoALLEnergies->Draw("same");
	
	//*************** third Column **********************************************************
	// textEta2760GeVNLOEtaALLEnergies->Draw();
	// textPi02760GeVNLOsysALLEnergies->Draw();
	// boxCombinedPi02760GeVALLEnergies->Draw("l");
	// markerCombinedPi02760GeVALLEnergies->DrawMarker(columnsNLOLegendPrelAndFinal[2]+offsetNLOLegendPrelAndFinalMarkerX,rowsNLOLegendPrelAndFinal[2]+offsetNLOLegendPrelAndFinalMarkerY);
	// lineFitPi02760GeVNLOALLEnergies->Draw("same");
	// lineNLOPi02760GeVMuHalfALLEnergies->Draw("same");
	// lineNLOPi02760GeVMuOneALLEnergies->Draw("same");
	// lineNLOPi02760GeVMuTwoALLEnergies->Draw("same");
		
	//**************** forth Column **********************************************************
	// textEta900GeVNLOEtaALLEnergies->Draw();
	// textPi0900GeVNLOsysALLEnergies->Draw();
	// boxCombinedPi0900GeVALLEnergies->Draw("l");
	// markerCombinedPi0900GeVALLEnergies->DrawMarker(columnsNLOLegendPrelAndFinal[3]+offsetNLOLegendPrelAndFinalMarkerX,rowsNLOLegendPrelAndFinal[2]+offsetNLOLegendPrelAndFinalMarkerY);
	// lineFitPi0900GeVNLOALLEnergies->Draw("same");
	// lineNLOPi0900GeVMuHalfALLEnergies->Draw("same");
	// lineNLOPi0900GeVMuOneALLEnergies->Draw("same");
	// lineNLOPi0900GeVMuTwoALLEnergies->Draw("same");
	// textArxivALLEnergies->Draw();
	
	
	canvasInvXSectionNLOOnlySpectraEta->Update();
	canvasInvXSectionNLOOnlySpectraEta->Print(Form("%s/Eta_InvXSectionNLO_OnlySpectrum_Paper.%s",outputDir.Data(),suffix.Data()));
	
	padComparisonXSectionNLOOnlySpectraEta->cd();

	//	if(!thesis)DrawAliceLogoPi0WithPHOSPrelim(pictDrawingCoordinatesOnlySpectrum[0]-0.05, pictDrawingCoordinatesOnlySpectrum[1], pictDrawingCoordinatesOnlySpectrum[2], pictDrawingCoordinatesOnlySpectrum[3], pictDrawingCoordinatesOnlySpectrum[4], pictDrawingCoordinatesOnlySpectrum[5], pictDrawingCoordinatesOnlySpectrum[6], pictDrawingCoordinatesOnlySpectrum[7], pictDrawingCoordinates[8],collisionSystemCombinedReallyAll, pictDrawingOptions[1], pictDrawingOptions[2], kFALSE,1200,960);
	
	canvasInvXSectionNLOOnlySpectraEta->Update();
	canvasInvXSectionNLOOnlySpectraEta->Print(Form("%s/Eta_InvXSectionNLO_OnlySpectrum.%s",outputDir.Data(),suffix.Data()));

	
	//-------------- Plotting ------------------------------------------------------
	TCanvas* canvasInvXSectionEta = new TCanvas("canvasInvXSectionEta","",200,10,1200,1200);  // gives the page size
	DrawGammaCanvasSettings( canvasInvXSectionEta,  0.15, 0.02, 0.03, 0.06);

	canvasInvXSectionEta->SetLogy(1);
	canvasInvXSectionEta->SetLogx(1);
	TH2F * histo2DInvXSectionEtaALLEnergiesSep = new TH2F("histo2DInvXSectionEtaALLEnergiesSep","histo2DInvXSectionEtaALLEnergiesSep",1000,0.4,25.,1000,20,1e10);
	SetStyleHistoTH2ForGraphs(histo2DInvXSectionEtaALLEnergiesSep, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )", 0.032,0.03, 0.04,0.025, 1,1.55);
	histo2DInvXSectionEtaALLEnergiesSep->DrawCopy(); 
	
	//	graphInvCrossSectionEtaComb2760GeVStatErr = ScaleGraph(graphInvCrossSectionEtaComb2760GeVStatErr,1e-1);
	DrawGammaSetMarkerTGraphAsym(graphInvCrossSectionEtaComb2760GeVStatErr, markerStyleCommmonSpectrum2760GeV,1.6*markerSizeCommonSpectrum, colorCommonSpectrumPi02760GeV, colorCommonSpectrumPi02760GeV, widthCommonSpectrumBoxes, kFALSE);
	graphInvCrossSectionEtaComb2760GeVStatErr->SetLineWidth(widthCommonErrors);
	//	graphInvCrossSectionEtaComb2760GeVSysErr = ScaleGraph(graphInvCrossSectionEtaComb2760GeVSysErr,1e-1);
	DrawGammaSetMarkerTGraphAsym(graphInvCrossSectionEtaComb2760GeVSysErr, markerStyleCommmonSpectrum2760GeV,1.6*markerSizeCommonSpectrum, colorCommonSpectrumPi02760GeV, colorCommonSpectrumPi02760GeV, widthCommonSpectrumBoxes, kTRUE,colorCommonSpectrumPi02760GeVBox );
	graphInvCrossSectionEtaComb2760GeVSysErr->SetLineWidth(widthCommonErrors);
	graphInvCrossSectionEtaComb2760GeVSysErr->Draw("E2same");
	graphInvCrossSectionEtaComb2760GeVStatErr->Draw("p,same");
	//histoFitInvCrossSectionEta2760GeV->Draw("same,c");
	graphNLOMuHalfEta2760GeV->Draw("same,c");
	graphNLOMuOneEta2760GeV->Draw("same,c");
	graphNLOMuTwoEta2760GeV->Draw("same,c");
	fitInvCrossSectionEta2760GeV->SetLineColor(colorCommonSpectrumPi02760GeV);
	fitInvCrossSectionEta2760GeV->SetLineWidth(2);
	fitInvCrossSectionEta2760GeV->Draw("c,same");
	//	labelScalingEta2760GeVALLEnergies->Draw();
	
	TLegend* legendEtaComb2760GeV = new TLegend(0.25,0.16,0.477,0.34);
	legendEtaComb2760GeV->SetTextSize(0.03);			
	legendEtaComb2760GeV->SetFillColor(0);
	legendEtaComb2760GeV->SetBorderSize(0);
	legendEtaComb2760GeV->AddEntry(graphInvCrossSectionEtaComb2760GeVSysErr,Form("pp #sqrt{#it{s}} = 2.76 TeV "),"fp");
	legendEtaComb2760GeV->AddEntry(fitInvCrossSectionEta2760GeV,"Fit Tsallis");
	legendEtaComb2760GeV->AddEntry(graphNLOMuHalfEta2760GeV,"NLO #mu = 0.5 #it{p}_{T}","l");
	legendEtaComb2760GeV->AddEntry(graphNLOMuOneEta2760GeV,"NLO #mu = #it{p}_{T}","l");
	legendEtaComb2760GeV->AddEntry(graphNLOMuTwoEta2760GeV,"NLO #mu = 2 #it{p}_{T}","l");

	legendEtaComb2760GeV->Draw();


	canvasInvXSectionEta->Update();
	canvasInvXSectionEta->Print(Form("%s/Eta_InvXSection.%s",outputDir.Data(),suffix.Data()));
	
	//******************************Going to save results*************************//
        TFile *fCombResults= new TFile(Form("CombinedResultsEta2760%s_%s.root",bWCorrection.Data(),dateForOutput.Data()),"RECREATE");
	//	graphRatioEtaPi0ComplErr2760GeV->Write("graphEtaToPi0PCM2760GeV");
	//	graphRatioEtaPi0SysErr2760GeV->Write("graphRatioEtaPi0SysErr2760GeV");   //  PCM binshifted

	graphRatioEtaPCMtoPi0Published2760GeVSys->Write("graphRatioEtaPCMtoPi0Published2760GeVSys");
	graphRatioEtaPCMtoPi0Published2760GeVStat->Write("graphRatioEtaPCMtoPi0Published2760GeVStat");

	  //graphInvCrossSectionEtaComb2760GeVYShiftedSys->Write("");
	  //graphInvCrossSectionEtaComb2760GeVYShiftedSys->Write("");
	  //graphInvCrossSectionEtaComb2760GeVYShiftedStat->Write("");
	  //graphInvCrossSectionEtaComb2760GeVYShiftedStat->Write("");
	  //graphInvCrossSectionEtaComb2760GeV->Write();

	//Saving Combined Eta 2760 GeV

        fitTsallisEta2760GeVPtYShift->Write("fitInvCrossSectionEtaComb2760GeV_YShift");

	graphInvCrossSectionEtaComb2760GeV->Write("graphInvCrossSectionEtaComb2760GeV");
	graphInvCrossSectionEtaComb2760GeVSysErr->Write("graphInvCrossSectionEtaComb2760GeVSysErr");
	graphInvCrossSectionEtaComb2760GeVStatErr->Write("graphInvCrossSectionEtaComb2760GeVStatErr");

	graphInvCrossSectionEtaComb2760GeVYShifted->Write("graphInvCrossSectionEtaComb2760GeV_YShifted");
	graphInvCrossSectionEtaComb2760GeVYShiftedSys->Write("graphInvCrossSectionEtaComb2760GeVSysErr_YShifted");
	graphInvCrossSectionEtaComb2760GeVYShiftedStat->Write("graphInvCrossSectionEtaComb2760GeVStatErr_YShifted");

	graphInvCrossSectionEtaComb2760GeVUnShifted->Write("graphInvCrossSectionEtaComb2760GeVUnShifted");
	graphInvCrossSectionEtaComb2760GeVSysErrUnShifted->Write("graphInvCrossSectionEtaComb2760GeVSysErrUnShifted");
	graphInvCrossSectionEtaComb2760GeVStatErrUnShifted->Write("graphInvCrossSectionEtaComb2760GeVStatErrUnShifted");


	// saving PCM contribution
	graphPCMXSectionStatEta2760GeVXShifted ->Write("graphPCMXSectionStatEta2760GeV_XShifted");
	graphPCMXSectionSysEta2760GeVXShifted ->Write("graphPCMXSectionSysEta2760GeV_XShifted");

	graphPCMXSectionEta2760GeVUnShifted->Write("graphInvCrossSectionEtaPCMStat2760GeVUnshifted");
 	graphPCMXSectionSysEta2760GeVUnShifted ->Write("graphInvCrossSectionEtaPCMSys2760GeVUnshifted");


	graphPCMXSectionEta2760GeVYShifted->Write("graphInvCrossSectionEtaPCMStat2760GeV_YShifted");
	graphPCMXSectionEta2760GeVYShiftedSys->Write("graphInvCrossSectionEtaPCMSys2760GeV_YShifted");
	graphPCMXSectionEta2760GeVYShiftedSysRAA->Write("graphInvCrossSectionEtaPCMSysForRAA2760GeV_YShifted");


	//saving EMCal contribution 

	graphEMCALXSectionStatEta2760GeVXShifted ->Write("graphEMCALXSectionStatEta2760GeV_XShifted");
	graphEMCALXSectionSysEta2760GeVXShifted ->Write("graphEMCALXSectionSysEta2760GeV_XShifted");



	graphEMCALXSectionEta2760GeVYShifted->Write("graphInvCrossSectionEtaEMCALStat2760GeV_YShifted");
	graphEMCALXSectionEta2760GeVYShiftedSys->Write("graphInvCrossSectionEtaEMCALSys2760GeV_YShifted");

	graphEMCALXSectionEta2760GeVUnShifted->Write("graphEMCALXSectionEta2760GeVUnShifted");
	graphEMCALXSectionSysEta2760GeVUnShifted->Write("graphEMCALXSectionSysEta2760GeVUnShifted");


        fCombResults->Close();
}


 
