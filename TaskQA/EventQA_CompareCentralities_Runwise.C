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
      const Int_t nSets           = 4;
      Size_t constMarkerSize      = 1;
      TString fEnergyFlag         = "PbPb_5.02TeV";
      TString period              = "LHC15o";                // used for marker style and color
      TString dataSet             = "LHC15o";                // used for filename and histo names inside 
      TString centralities[nSets] = {"0-10%","10-20%","20-50%","50-90%"};      // used for legend and marker style and color
      TString cutsDataSets[nSets] = {
	"10110013_00200009247602008250404000_0652501500000000",
	"11210013_00200009247602008250404000_0652501500000000",
	"12510013_00200009247602008250404000_0652501500000000",
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

    Bool_t drawVerticalLines = kFALSE;
    Int_t nLines;        // number of vertical lines
    Int_t runRanges[10]; // array of bin numbers where to draw vertical lines
    TLine* verticalLines[10];
    if (fEnergyFlag.CompareTo("PbPb_5.02TeV") == 0){
      drawVerticalLines = kTRUE;
      nLines = 6;
      runRanges[0] = 10; runRanges[1] = 16; runRanges[2] = 21;
      runRanges[3] = 24; runRanges[4] = 60; runRanges[5] = 63;
    }
    if (nLines > 10) cout << "ERROR: nLines cannot be larger than 10. Increase size of runRanges[10] and verticalLines[10]" << endl;

    TString fDetectionProcess = ReturnFullTextReconstructionProcess(mode);

    std::vector<TH1D*> vecHistosRunwise[nSets];
    std::vector<TString> vecHistosRunwiseName[nSets];
    std::vector<TString> vecHistosRunwiseNameSave[nSets];

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
        vecHistosRunwiseNameSave[i].push_back("NEvents");
        if(i==0)nRangeRunwise++;

        vecHistosRunwise[i].push_back(hNEventsMinBias[i]);
        vecHistosRunwiseName[i].push_back("hNEventsMinBias");
        vecHistosRunwiseNameSave[i].push_back("NEventsMinBias");
        if(i==0)nRangeRunwise++;

        vecHistosRunwise[i].push_back(hNEventsAll[i]);
        vecHistosRunwiseName[i].push_back("hNEventsAll");
        vecHistosRunwiseNameSave[i].push_back("NEventsAll");
        if(i==0)nRangeRunwise++;

        vecHistosRunwise[i].push_back(hNEventsFracGood[i]);
        vecHistosRunwiseName[i].push_back("hNEventsFracGoodEvents");
        vecHistosRunwiseNameSave[i].push_back("NEventsFracGoodEvents");
        if(i==0)nRangeRunwise++;
                                          
        vecHistosRunwise[i].push_back(hNEventsFracGood[i]);
        vecHistosRunwiseName[i].push_back("hNEventsFracGoodEvents");
        vecHistosRunwiseNameSave[i].push_back("NEventsFracGoodEvents");
        if(i==0)nRangeRunwise++;

        vecHistosRunwise[i].push_back(hNEventsFracNorm[i]);
        vecHistosRunwiseName[i].push_back("hNEventsFracNormAll");
        vecHistosRunwiseNameSave[i].push_back("NEventsFracNormAll");
        if(i==0)nRangeRunwise++;
                                          
        vecHistosRunwise[i].push_back(hNEventsFracMinBias[i]);
        vecHistosRunwiseName[i].push_back("hNEventsFracMinBias");
        vecHistosRunwiseNameSave[i].push_back("NEventsFracMinBias");
        if(i==0)nRangeRunwise++;
                                          
        vecHistosRunwise[i].push_back(hTracksMeanGood[i]);
        vecHistosRunwiseName[i].push_back("hTracksGood-Mean");
        vecHistosRunwiseNameSave[i].push_back("TracksGood-Mean");
        if(i==0)nRangeRunwise++;

        vecHistosRunwise[i].push_back(hTracksRMSGood[i]);
        vecHistosRunwiseName[i].push_back("hTracksGood-RMS");
        vecHistosRunwiseNameSave[i].push_back("TracksGood-RMS");
        if(i==0)nRangeRunwise++;

        vecHistosRunwise[i].push_back(hVertexZMean[i]);
        vecHistosRunwiseName[i].push_back("hVertexZ-Mean");
        vecHistosRunwiseNameSave[i].push_back("VertexZ-Mean");
        if(i==0)nRangeRunwise++;

        vecHistosRunwise[i].push_back(hVertexZRMS[i]);
        vecHistosRunwiseName[i].push_back("hVertexZ-RMS");
        vecHistosRunwiseNameSave[i].push_back("VertexZ-RMS");
        if(i==0)nRangeRunwise++;

        vecHistosRunwise[i].push_back(hFracWtxOut[i]);
        vecHistosRunwiseName[i].push_back("hFracWVtxOutside10cm");
        vecHistosRunwiseNameSave[i].push_back("FracWVtxOutside10cm");
        if(i==0)nRangeRunwise++;

        vecHistosRunwise[i].push_back(hFracSPDClusTrack[i]);
        vecHistosRunwiseName[i].push_back("hFracSPDClusTrack");
        vecHistosRunwiseNameSave[i].push_back("FracSPDClusTrack");
        if(i==0)nRangeRunwise++;
  
        vecHistosRunwise[i].push_back(hFracPileUp[i]);
        vecHistosRunwiseName[i].push_back("hFracPileUp");
        vecHistosRunwiseNameSave[i].push_back("FracPileUp");
        if(i==0)nRangeRunwise++;

        vecHistosRunwise[i].push_back(hFracWOVtx[i]);
        vecHistosRunwiseName[i].push_back("hFracWOVtx");
        vecHistosRunwiseNameSave[i].push_back("FracWOVtx");
        if(i==0)nRangeRunwise++;
        
	vecHistosRunwise[i].push_back(hConvNCandidates[i]);
	vecHistosRunwiseName[i].push_back("hConvNCandidates");
	vecHistosRunwiseNameSave[i].push_back("ConvNCandidates");
	if(i==0)nRangeRunwise++;

	vecHistosRunwise[i].push_back(hConvNCandidatesQA[i]);
	vecHistosRunwiseName[i].push_back("hConvNCandidatesQA");
	vecHistosRunwiseNameSave[i].push_back("ConvNCandidatesQA");
	if(i==0)nRangeRunwise++;

	vecHistosRunwise[i].push_back(hCentrality[i]);
	vecHistosRunwiseName[i].push_back("hCentrality-Mean");
	vecHistosRunwiseNameSave[i].push_back("Centrality-Mean");
	if(i==0)nRangeRunwise++;

	vecHistosRunwise[i].push_back(hEventPlaneAngle[i]);
	vecHistosRunwiseName[i].push_back("hEventPlaneAngle-Mean");
	vecHistosRunwiseNameSave[i].push_back("EventPlaneAngle-Mean");
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

    Int_t nTopRows = 1;
    TCanvas* canvas = new TCanvas("canvas","",200,10,1350,(846.+nTopRows*54.));  // gives the page size
    Double_t leftMar = 0.09; Double_t rightMar = 0.02; Double_t topMargin = (nTopRows*54.)/(846.+nTopRows*54.); Double_t bottomMargin = 0.09;
    DrawGammaCanvasSettings(canvas, leftMar, rightMar, topMargin, bottomMargin);
    Double_t yMax;  // min of y range for vertical lines
    Double_t yMin;  // max of y range for vertical lines
    Bool_t adjustedRange = kFALSE;
//****************************** Combined Trending Histograms ************************************************

    TLegend *legend = new TLegend(0.12,0.99 - (nTopRows*0.04),0.95,0.98);
    legend->SetNColumns(4);
    legend->SetFillColor(0);
    legend->SetLineColor(0);
    legend->SetTextSize(0.03);
    legend->SetTextFont(42);

//---------

        for(Int_t h=0; h<(Int_t)vecHistosRunwise[0].size(); h++){
            cout << h << " " << vecHistosRunwiseName[0].at(h).Data() <<  ", " ;
            Double_t scaleFactor    = 1.2;
            if ( h==0 || vecHistosRunwiseName[0].at(h).Contains("hConvNCandidates") )
                scaleFactor    = 12.;
            if (vecHistosRunwiseName[0].at(h).Contains("hTracksGood-Mean"))
                scaleFactor    = 2.5;
            
            adjustedRange = AdjustHistRange(vecHistosRunwise, scaleFactor, scaleFactor, h, nSets, kTRUE, &yMin, &yMax);
            for(Int_t i=0; i<nSets; i++){
                TString draw;
                if(h==0) draw = (i==0)?"p":"p, same";
                else draw = (i==0)?"px0e1":"px0e1, same";
                ((TH1D*) vecHistosRunwise[i].at(h))->Draw(draw.Data());
                legend->AddEntry(((TH1D*) vecHistosRunwise[i].at(h)),centralities[i].Data(),"p");
            }
            legend->Draw();           
	    if (drawVerticalLines){
	      if(!adjustedRange){
		canvas->Update();
		yMax = gPad->GetUymax();
		yMin = gPad->GetUymin();
	      }
	      for(Int_t lineBin=0; lineBin<nLines; lineBin++){
		verticalLines[lineBin] = new TLine(runRanges[lineBin],yMin,runRanges[lineBin],yMax);
		verticalLines[lineBin]->SetLineWidth(1);
		verticalLines[lineBin]->SetLineStyle(2);
		verticalLines[lineBin]->Draw("same");
	      }
	    }
	    
	    PutProcessLabelAndEnergyOnPlot(0.75, 0.97-nTopRows*0.06, 0.03, fCollisionSystem.Data(), "", "");
            if ( h==0 ||  vecHistosRunwiseName[0].at(h).Contains("hConvNCandidates") || vecHistosRunwiseName[0].at(h).Contains("hTracksGood-Mean")
               ) 
               SaveWriteCanvas(canvas, Form("EventQA_CompareCentralities/%s/Runwise/%s.%s", dataSet.Data(), vecHistosRunwiseNameSave[0].at(h).Data(),suffix.Data()), kFALSE, kTRUE);
            else SaveWriteCanvas(canvas, Form("EventQA_CompareCentralities/%s/Runwise/%s.%s", dataSet.Data(), vecHistosRunwiseNameSave[0].at(h).Data(),suffix.Data()));
            legend->Clear();
	    if(drawVerticalLines){
	      for(Int_t lineBin=0; lineBin<nLines; lineBin++){
		delete verticalLines[lineBin];
	      }
	    }
        }
        cout << endl;

    fOutput->Write();
    fOutput->Close();

    cout << "Done with EventQA_CompareCentralities_Runwise" << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;


    return;
}
