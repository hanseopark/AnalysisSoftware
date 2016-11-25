/*******************************************************************************
 ******  provided by Gamma Conversion Group, PWGGA,                        *****
 ******     Daniel Muehlheim, d.muehlheim@cern.ch                          ***** 
 ******     Friederike Bock, fbock@cern.ch                                 ***** 
 *******************************************************************************/

#include "TaskQA/QA.h"

void EventQA_CompareCentralities_Runwise( TString suffix  = "eps", 
			                   Int_t mode      = 0
                      ){
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "EventQA_CompareCentralities_Runwise" << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

    gROOT->Reset();

    StyleSettingsThesis();
    SetPlotStyle();                                  

    //**************************************************************************************************************
    // LHC15o
      const Int_t nSets           = 2;
      Size_t constMarkerSize      = 1;
      TString fEnergyFlag         = "PbPb_5.02TeV";
      TString period              = "LHC15o";                // used for marker style and color
      TString dataSet             = "LHC15o";                // used for filename and histo names inside 
      TString centralities[nSets] = {"0-10%","50-90%"};      // used for legend and marker style and color
      TString cutsDataSets[nSets] = {
	"10110013_00200009247602008250404000_0652501500000000",
	//"11210013_00200009247602008250404000_0652501500000000",
	//"12510013_00200009247602008250404000_0652501500000000",
	"15910013_00200009247602008250404000_0652501500000000"  
      };
                      
//**************************************************************************************************************
    Style_t hMarkerStyle[nSets];
    Size_t hMarkerSize[nSets];
    Color_t hMarkerColor[nSets];
    Color_t hLineColor[nSets];

    TString pathDataSets[nSets];

    TString fCutSelection[nSets];
    TString fEventCutSelection[nSets], fGammaCutSelection[nSets], fMesonCutSelection[nSets];
    TString dummyCutSelection[nSets];

    for(Int_t i=0; i<nSets; i++)
    {
        pathDataSets[i] = Form("%s/%s/EventQA", cutsDataSets[i].Data(),fEnergyFlag.Data());
	hMarkerStyle[i] = GetDefaultMarkerStyle(fEnergyFlag, period.Data() , centralities[i].Data());
	hMarkerSize[i]  = constMarkerSize;
	hMarkerColor[i] = GetColorDefaultColor(fEnergyFlag, period.Data() , centralities[i].Data());
	hLineColor[i]   = GetColorDefaultColor(fEnergyFlag, period.Data() , centralities[i].Data());

	ReturnSeparatedCutNumberAdvanced(cutsDataSets[i], fEventCutSelection[i], fGammaCutSelection[i], dummyCutSelection[i],dummyCutSelection[i], fMesonCutSelection[i], mode);
	cout << "set " << i << ": " << fEventCutSelection[i].Data() << ", " << fGammaCutSelection[i].Data() << ", " << fMesonCutSelection[i].Data() << endl;

    }

    TString fCollisionSystem = ReturnFullCollisionsSystem(fEnergyFlag);
    if (fCollisionSystem.CompareTo("") == 0){
        cout << "No correct collision system specification, has been given" << endl;
        return;
    }

    TString fDetectionProcess = ReturnFullTextReconstructionProcess(mode);

    std::vector<TH1D*> vecHistosRunwise[nSets];
    std::vector<TString> vecHistosRunwiseName[nSets];

    TH1D* hNEvents[nSets];
    TH1D* hNEventsMinBias[nSets];
    TH1D* hNEventsAll[nSets];
    TH1D* hNEventsFracGood[nSets];
    TH1D* hNEventsFracNorm[nSets];
    TH1D* hNEventsFracMinBias[nSets];

    TH1D* hTracksMeanGood[nSets];
    TH1D* hTracksRMSGood[nSets];
    TH1D* hVertexZMean[nSets];
    TH1D* hVertexZRMS[nSets];
    TH1D* hFracWtxOut[nSets];
    TH1D* hFracSPDClusTrack[nSets];
    TH1D* hFracPileUp[nSets];
    TH1D* hFracWOVtx[nSets];    
    TH1D* hConvNCandidates[nSets];
    TH1D* hConvNCandidatesQA[nSets];
    TH1D* hCentrality[nSets];
    TH1D* hEventPlaneAngle[nSets];

    Int_t nRangeRunwise = 0;

    for(Int_t i=0; i<nSets; i++){
        
        vecHistosRunwise[i].push_back(hNEvents[i]);
        vecHistosRunwiseName[i].push_back("hNEvents");
        if(i==0)nRangeRunwise++;

        vecHistosRunwise[i].push_back(hNEventsMinBias[i]);
        vecHistosRunwiseName[i].push_back("hNEventsMinBias");
        if(i==0)nRangeRunwise++;

        vecHistosRunwise[i].push_back(hNEventsAll[i]);
        vecHistosRunwiseName[i].push_back("hNEventsAll");
        if(i==0)nRangeRunwise++;

        vecHistosRunwise[i].push_back(hNEventsFracGood[i]);
        vecHistosRunwiseName[i].push_back("hNEventsFracGoodEvents");
        if(i==0)nRangeRunwise++;
                                          
        vecHistosRunwise[i].push_back(hNEventsFracGood[i]);
        vecHistosRunwiseName[i].push_back("hNEventsFracGoodEvents");
        if(i==0)nRangeRunwise++;

        vecHistosRunwise[i].push_back(hNEventsFracNorm[i]);
        vecHistosRunwiseName[i].push_back("hNEventsFracNormAll");
        if(i==0)nRangeRunwise++;
                                          
        vecHistosRunwise[i].push_back(hNEventsFracMinBias[i]);
        vecHistosRunwiseName[i].push_back("hNEventsFracMinBias");
        if(i==0)nRangeRunwise++;
                                          
        vecHistosRunwise[i].push_back(hTracksMeanGood[i]);
        vecHistosRunwiseName[i].push_back("hTracksGood-Mean");
        if(i==0)nRangeRunwise++;

        vecHistosRunwise[i].push_back(hTracksRMSGood[i]);
        vecHistosRunwiseName[i].push_back("hTracksGood-RMS");
        if(i==0)nRangeRunwise++;

        vecHistosRunwise[i].push_back(hVertexZMean[i]);
        vecHistosRunwiseName[i].push_back("hVertexZ-Mean");
        if(i==0)nRangeRunwise++;

        vecHistosRunwise[i].push_back(hVertexZRMS[i]);
        vecHistosRunwiseName[i].push_back("hVertexZ-RMS");
        if(i==0)nRangeRunwise++;

        vecHistosRunwise[i].push_back(hFracWtxOut[i]);
        vecHistosRunwiseName[i].push_back("hFracWVtxOutside10cm");
        if(i==0)nRangeRunwise++;

        vecHistosRunwise[i].push_back(hFracSPDClusTrack[i]);
        vecHistosRunwiseName[i].push_back("hFracSPDClusTrack");
        if(i==0)nRangeRunwise++;
  
        vecHistosRunwise[i].push_back(hFracPileUp[i]);
        vecHistosRunwiseName[i].push_back("hFracPileUp");
        if(i==0)nRangeRunwise++;

        vecHistosRunwise[i].push_back(hFracWOVtx[i]);
        vecHistosRunwiseName[i].push_back("hFracWOVtx");
        if(i==0)nRangeRunwise++;
        
	vecHistosRunwise[i].push_back(hConvNCandidates[i]);
	vecHistosRunwiseName[i].push_back("hConvNCandidates");
	if(i==0)nRangeRunwise++;

	vecHistosRunwise[i].push_back(hConvNCandidatesQA[i]);
	vecHistosRunwiseName[i].push_back("hConvNCandidatesQA");
	if(i==0)nRangeRunwise++;

	vecHistosRunwise[i].push_back(hCentrality[i]);
	vecHistosRunwiseName[i].push_back("hCentrality-Mean");
	if(i==0)nRangeRunwise++;

	vecHistosRunwise[i].push_back(hEventPlaneAngle[i]);
	vecHistosRunwiseName[i].push_back("hEventPlaneAngle-Mean");
	if(i==0)nRangeRunwise++;
        
    }

    for(Int_t i=0; i<nSets; i++)    {
        TString fFileRunwise = Form("%s/%s_EventQARunwise.root", pathDataSets[i].Data(), dataSet.Data());
        TFile* FileRunwise = new TFile(fFileRunwise.Data(),"READ");
        if(FileRunwise->IsZombie()) {cout << "Warning: ROOT file '" << fFileRunwise.Data() << "' could not be openend!" << endl;}

        cout << "\t\t----------------------------------------------------------------------------" << endl;
        cout << "\t\t\tNumber of Runwise Histograms being processed: " << nRangeRunwise;
        cout << "\t\t----------------------------------------------------------------------------" << endl;

        for(Int_t j=0; j<nRangeRunwise; j++)
        {
	    cout << Form("%s_%s", vecHistosRunwiseName[i].at(j).Data(), dataSet.Data());
            TH1D* temp = (TH1D*) FileRunwise->Get(Form("%s_%s", vecHistosRunwiseName[i].at(j).Data(), dataSet.Data()));
            cout << ", ";
            OnlyEditTH1(temp, hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
            temp->SetTitle("");
            vecHistosRunwise[i].at(j) = new TH1D(*temp);
        }
    }

    gSystem->Exec("mkdir -p EventQA_CompareCentralities/"+dataSet);
    gSystem->Exec("mkdir -p EventQA_CompareCentralities/"+dataSet+"/Runwise");

    TString nameOutput = Form("EventQA_CompareCentralities/%s/Output_Runwise.root",dataSet.Data());
    TFile* fOutput = new TFile(nameOutput,"RECREATE");
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "Output file: " << nameOutput << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

    Int_t nTopRows = (nSets-1)/(3) + 1;
    TCanvas* canvas = new TCanvas("canvas","",200,10,1350,(846.+nTopRows*54.));  // gives the page size
    Double_t leftMar = 0.09; Double_t rightMar = 0.02; Double_t topMargin = (nTopRows*54.)/(846.+nTopRows*54.); Double_t bottomMargin = 0.09;
    DrawGammaCanvasSettings(canvas, leftMar, rightMar, topMargin, bottomMargin);

//****************************** Combined Trending Histograms ************************************************

    TLegend *legend = new TLegend(0.12,0.99 - (nTopRows*0.04),0.95,0.98);
    legend->SetNColumns(3);
    legend->SetFillColor(0);
    legend->SetLineColor(0);
    legend->SetTextSize(0.03);
    legend->SetTextFont(42);
//---------

        for(Int_t h=0; h<(Int_t)vecHistosRunwise[0].size(); h++){
            cout << h << " " << vecHistosRunwiseName[0].at(h).Data() <<  ", " ;
            Double_t scaleFactor    = 1.2;
            if ( h==0 || vecHistosRunwiseName[0].at(h).Contains("hCaloNClusters") || vecHistosRunwiseName[0].at(h).Contains("hCaloMergedNClusters") ||  
                vecHistosRunwiseName[0].at(h).Contains("hConvNCandidates") ) 
                scaleFactor    = 12.;
            if (vecHistosRunwiseName[0].at(h).Contains("hTracksGood-Mean"))
                scaleFactor    = 2.5;
            
            AdjustHistRange(vecHistosRunwise, scaleFactor, scaleFactor, h, nSets, kTRUE);
            for(Int_t i=0; i<nSets; i++){
                TString draw;
                if(h==0) draw = (i==0)?"p":"p, same";
                else draw = (i==0)?"px0e1":"px0e1, same";
                ((TH1D*) vecHistosRunwise[i].at(h))->Draw(draw.Data());
                legend->AddEntry(((TH1D*) vecHistosRunwise[i].at(h)),centralities[i].Data(),"p");
            }
            legend->Draw();           
	    
	    PutProcessLabelAndEnergyOnPlot(0.75, 0.97-nTopRows*0.06, 0.03, fCollisionSystem.Data(), "", "");
            if ( h==0 || vecHistosRunwiseName[0].at(h).Contains("hCaloNClusters") || vecHistosRunwiseName[0].at(h).Contains("hCaloMergedNClusters") ||  
                vecHistosRunwiseName[0].at(h).Contains("hConvNCandidates") || vecHistosRunwiseName[0].at(h).Contains("hTracksGood-Mean")  
               ) 
               SaveWriteCanvas(canvas, Form("EventQA_CompareCentralities/%s/Runwise/%s.%s", dataSet.Data(), vecHistosRunwiseName[0].at(h).Data(),suffix.Data()), kFALSE, kTRUE);
            else SaveWriteCanvas(canvas, Form("EventQA_CompareCentralities/%s/Runwise/%s.%s", dataSet.Data(), vecHistosRunwiseName[0].at(h).Data(),suffix.Data()));
            legend->Clear();
        }
        cout << endl;

    fOutput->Write();
    fOutput->Close();

    cout << "Done with EventQA_CompareCentralities_Runwise" << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;


    return;
}
