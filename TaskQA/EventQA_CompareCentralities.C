/*******************************************************************************
 ******  provided by Gamma Conversion Group, PWGGA,                        *****
 ******     Daniel Muehlheim, d.muehlheim@cern.ch                          ***** 
 ******     Friederike Bock, fbock@cern.ch                                 ***** 
 *******************************************************************************/

#include "QA.h"

void EventQA_CompareCentralities( TString suffix  = "eps", 
                                  Int_t mode      = 0
                      ){
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "EventQA_CompareCentralities" << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

    gROOT->Reset();

    StyleSettingsThesis();
    SetPlotStyle();

    //choose which data sets to process
    //**************************************************************************************************************
    // LHC15o
    const Int_t nSets           = 4;
    Size_t constMarkerSize      = 1;
    TString fEnergyFlag         = "PbPb_5.02TeV";
    TString dataSet             = "LHC15o";                                 // used for filename, GetDefaultMarkerStyle 
    TString MCSet               = "LHC16g1";
    TString centralities[nSets] = {"0-10%","10-20%","20-50%","50-90%"};     // used for plot and marker style and color
    TString cutsDataSets[nSets] = {
      "10110013_00200009247602008250404000_0652501500000000",
      "11210013_00200009247602008250404000_0652501500000000",
      "12510013_00200009247602008250404000_0652501500000000",
      "15910013_00200009247602008250404000_0652501500000000"  
    };
                     
//**************************************************************************************************************
    Style_t hMarkerStyle[nSets];
    Style_t hMarkerStyleMC[nSets];
    Size_t hMarkerSize[nSets];
    Size_t hMarkerSizeMC[nSets];
    Color_t hMarkerColor[nSets];
    Color_t hMarkerColorMC[nSets];
    Color_t hLineColor[nSets];
    Color_t hLineColorMC[nSets];

    TString pathDataSets[nSets];

    TString fCutSelection[nSets];
    TString fEventCutSelection[nSets], fGammaCutSelection[nSets], fMesonCutSelection[nSets];
    TString dummyCutSelection[nSets];

    for(Int_t i=0; i<nSets; i++)
    {
        pathDataSets[i] = Form("%s/%s/EventQA", cutsDataSets[i].Data(),fEnergyFlag.Data());
	hMarkerStyle[i] = GetDefaultMarkerStyle(fEnergyFlag, dataSet.Data() , centralities[i].Data());
	hMarkerStyleMC[i] = GetDefaultMarkerStyle(fEnergyFlag, MCSet.Data() , centralities[i].Data());
	hMarkerSize[i]  = constMarkerSize;
	hMarkerColor[i] = GetColorDefaultColor(fEnergyFlag, dataSet.Data() , centralities[i].Data());
	hMarkerColorMC[i] = GetColorDefaultColor(fEnergyFlag, MCSet.Data() , centralities[i].Data());
	hLineColor[i]   = GetColorDefaultColor(fEnergyFlag, dataSet.Data() , centralities[i].Data());
	hLineColorMC[i]   = GetColorDefaultColor(fEnergyFlag, MCSet.Data() , centralities[i].Data());

	ReturnSeparatedCutNumberAdvanced(cutsDataSets[i], fEventCutSelection[i], fGammaCutSelection[i], dummyCutSelection[i],dummyCutSelection[i], fMesonCutSelection[i], mode);
	cout << "set " << i << ": " << fEventCutSelection[i].Data() << ", " << fGammaCutSelection[i].Data() << ", " << fMesonCutSelection[i].Data() << endl;
  }


    TString fCollisionSystem = ReturnFullCollisionsSystem(fEnergyFlag);
    if (fCollisionSystem.CompareTo("") == 0){
        cout << "No correct collision system specification, has been given" << endl;
        return;
    }

    TString fDetectionProcess = ReturnFullTextReconstructionProcess(mode);

    std::vector<TH1D*> vecHistos[nSets];                       // data histos
    std::vector<TH1D*> vecHistosMC[nSets];                     // MC histos
    std::vector<TString> vecHistosName[nSets];                 // data and MC histos name for reading in
    std::vector<TString> vecHistosNameForSaving[nSets];        // names for saving data histos and data-MC-comparison histos

    TH1D* hcVertexZ[nSets];
    TH1D* hcGoodESDTracks[nSets];
    TH1D* hcV0Mult[nSets];
    TH1D* hcGammaCandidates[nSets];
    TH1D* hcCentrality[nSets];
    TH1D* hcEventPlaneAngle[nSets];
    TH1D* hcNEvents[nSets];

    TH1D* hcVertexZMC[nSets];
    TH1D* hcGoodESDTracksMC[nSets];
    TH1D* hcV0MultMC[nSets];
    TH1D* hcGammaCandidatesMC[nSets];
    TH1D* hcCentralityMC[nSets];
    TH1D* hcEventPlaneAngleMC[nSets];
    TH1D* hcNEventsMC[nSets];

    Int_t nRange = 0;

    for(Int_t i=0; i<nSets; i++){
        
        vecHistos[i].push_back(hcVertexZ[i]);
	vecHistosMC[i].push_back(hcVertexZMC[i]);
        vecHistosName[i].push_back("VertexZ");
        vecHistosNameForSaving[i].push_back("hVertex_Z");
        if(i==0)nRange++;

        vecHistos[i].push_back(hcGoodESDTracks[i]);
	vecHistosMC[i].push_back(hcGoodESDTracksMC[i]);
        vecHistosName[i].push_back("GoodESDTracks");
        vecHistosNameForSaving[i].push_back("hNGoodTracks");
        if(i==0)nRange++;

        
        vecHistos[i].push_back(hcV0Mult[i]);
	vecHistosMC[i].push_back(hcV0MultMC[i]);
        vecHistosName[i].push_back("V0 Multiplicity");
        vecHistosNameForSaving[i].push_back("hV0Mult");
        if(i==0)nRange++;
    
        vecHistos[i].push_back(hcGammaCandidates[i]);
        vecHistosMC[i].push_back(hcGammaCandidatesMC[i]);
        vecHistosName[i].push_back("GammaCandidates");
        vecHistosNameForSaving[i].push_back("hGammaCandidates");
        if(i==0)nRange++;

        vecHistos[i].push_back(hcEventPlaneAngle[i]);
        vecHistosMC[i].push_back(hcEventPlaneAngleMC[i]);
        vecHistosName[i].push_back(Form("EventPlaneAngle %s",fEventCutSelection[i].Data()));
        vecHistosNameForSaving[i].push_back("hEventPlaneAngle");
        if(i==0)nRange++;
       
        vecHistos[i].push_back(hcNEvents[i]);
	vecHistosMC[i].push_back(hcNEventsMC[i]);
        vecHistosName[i].push_back("NEvents");
        vecHistosNameForSaving[i].push_back("hNoEvents");
        if(i==0)nRange++;

    }

    for(Int_t i=0; i<nSets; i++)    {
        TString fFile = Form("%s/EventQA_%s.root", pathDataSets[i].Data(), dataSet.Data());
        TString fFileMC = Form("%s/EventQA_%s.root", pathDataSets[i].Data(), MCSet.Data());
        TFile* File = new TFile(fFile.Data(),"READ");
        TFile* FileMC = new TFile(fFileMC.Data(),"READ");

        if(File->IsZombie()) {cout << "Warning: ROOT file '" << fFile.Data() << "' could not be openend!" << endl; return;}
        if(FileMC->IsZombie()) {cout << "Warning: ROOT file '" << fFile.Data() << "' could not be openend!" << endl; return;}

        cout << endl;
        cout << "\t\t----------------------------------------------------------------------------" << endl;
        cout << "\t\tDataSet: " << dataSet.Data() << ", " << centralities[i].Data() << ":" << endl;
        cout << "\t\tProcessing file: " << fFile.Data() << endl;
        cout << "\t\tMC:      " << MCSet.Data() << ", " << centralities[i].Data() << ":" << endl;
        cout << "\t\tProcessing file: " << fFileMC.Data() << endl;
        cout << "\t\t----------------------------------------------------------------------------" << endl;
        cout << endl;

        cout << "\n\t\t----------------------------------------------------------------------------" << endl;
        cout << "\t\t\tNumber of Histograms being processed: " << nRange << endl;
        cout << "\t\t----------------------------------------------------------------------------" << endl;

        for(Int_t j=0; j<nRange; j++) // go through histograms, scale and add to vector
        {
            cout << ((TString)vecHistosName[i].at(j)).Data();
            TH1D* temp = (TH1D*) File->Get(((TString)vecHistosName[i].at(j)).Data());
            TH1D* tempMC = (TH1D*) FileMC->Get(((TString)vecHistosName[i].at(j)).Data());
            cout << ", ";
            temp->Sumw2();
            tempMC->Sumw2();
	    if(j!=5){ // not for NoEvents histo
	      // scale histograms with number of entries:
	      Double_t tempEntries = temp->GetEntries();
	      Double_t tempEntriesMC = tempMC->GetEntries();
	      temp->Scale(1./tempEntries);
	      tempMC->Scale(1./tempEntriesMC);
	      temp->GetYaxis()->SetTitle(Form("#frac{1}{N} %s",temp->GetYaxis()->GetTitle()));
	      tempMC->GetYaxis()->SetTitle(Form("#frac{1}{N} %s",tempMC->GetYaxis()->GetTitle()));
	    }
	    OnlyEditTH1(temp, hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
	    OnlyEditTH1(tempMC, hMarkerStyleMC[i], hMarkerSize[i], hMarkerColorMC[i], hLineColorMC[i]);
            temp->GetXaxis()->SetTitleOffset(1);
            tempMC->GetXaxis()->SetTitleOffset(1);
            temp->GetYaxis()->SetTitleOffset(1.1);
            tempMC->GetYaxis()->SetTitleOffset(1.1);
            temp->SetTitle("");
            tempMC->SetTitle("");
            vecHistos[i].at(j) = new TH1D(*temp);
            vecHistosMC[i].at(j) = new TH1D(*tempMC);
        }

        cout << "done!" << endl;
    }

    gSystem->Exec("mkdir -p EventQA_CompareCentralities/"+dataSet);
    gSystem->Exec("mkdir -p EventQA_CompareCentralities/"+dataSet+"/Periodwise");
    gSystem->Exec("mkdir -p EventQA_CompareCentralities/Comparison/"+dataSet+"_"+MCSet+"/Periodwise");

    TString nameOutput = Form("EventQA_CompareCentralities/%s/Output_Periodwise.root",dataSet.Data());
    TFile* fOutput = new TFile(nameOutput,"RECREATE");
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "Output file: " << nameOutput << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

    Int_t nTopRows = 1;
    TCanvas* canvas = new TCanvas("canvas","",200,10,1350,(846.+nTopRows*54.));  // gives the page size
    canvas->SetLogy(1);
    Double_t leftMar = 0.09; Double_t rightMar = 0.02; Double_t topMargin = (nTopRows*54.)/(846.+nTopRows*54.); Double_t bottomMargin = 0.09;
    DrawGammaCanvasSettings(canvas, leftMar, rightMar, topMargin, bottomMargin);
    Double_t scaleFactor = 5.0;

//************************************* Data histograms ************************************************
    TLegend *legend = new TLegend(0.12,0.99 - (nTopRows*0.04),0.95,0.98);
    legend->SetNColumns(4);
    legend->SetFillColor(0);
    legend->SetLineColor(0);
    legend->SetTextSize(0.03);
    legend->SetTextFont(42);

    for(Int_t h=0; h<(Int_t)vecHistos[0].size(); h++){
        cout << h << ", ";
	if ( vecHistosNameForSaving[0].at(h).Contains("VertexZ") || vecHistosNameForSaving[0].at(h).Contains("EventPlaneAngle")) scaleFactor = 1.2;
        AdjustHistRange(vecHistos,scaleFactor,scaleFactor,h,nSets,kFALSE);
	if ( vecHistosNameForSaving[0].at(h).Contains("V0Mult") ) ((TH1D*) vecHistos[0].at(h))->GetXaxis()->SetRangeUser(0,(((TH1D*) vecHistos[0].at(h))->GetXaxis()->GetXmax()));
        for(Int_t i=0; i<nSets; i++) {
            TString draw;
            if(h==0) draw = (i==0)?"p":"p, same";
            if(h==5) draw = (i==0)?"histo":"histo, same";
            else draw = (i==0)?"px0e1":"px0e1, same";
            ((TH1D*) vecHistos[i].at(h))->Draw(draw.Data());
	    if(h==5)legend->AddEntry(((TH1D*) vecHistos[i].at(h)),centralities[i],"l");
            else legend->AddEntry(((TH1D*) vecHistos[i].at(h)),centralities[i],"p");
        }
        legend->Draw();
        if ( vecHistosNameForSaving[0].at(h).Contains("Candidates") ) 
            PutProcessLabelAndEnergyOnPlot(0.75, 0.97-nTopRows*0.06, 0.03, fCollisionSystem.Data(),"","");
        else 
            PutProcessLabelAndEnergyOnPlot(0.75, 0.97-nTopRows*0.06, 0.03, fCollisionSystem.Data(), "", "");
	if ( vecHistosNameForSaving[0].at(h).Contains("VertexZ") || vecHistosNameForSaving[0].at(h).Contains("EventPlaneAngle")) SaveWriteCanvas(canvas, Form("EventQA_CompareCentralities/%s/Periodwise/%s.%s", dataSet.Data(), vecHistosNameForSaving[0].at(h).Data(),suffix.Data()));
        else SaveWriteCanvas(canvas, Form("EventQA_CompareCentralities/%s/Periodwise/%s.%s", dataSet.Data(), vecHistosNameForSaving[0].at(h).Data(),suffix.Data()), kFALSE, kTRUE);
	legend->Clear();
        // if(h==9) TGaxis::SetExponentOffset(0, 0, "x");
    }
    cout << endl;

//******************************************* Data & MC histograms ***************************************************

    nTopRows = 2;
    TCanvas* canvasComparison = new TCanvas("canvasComparison","",200,10,1350,(846.+nTopRows*54.));  // gives the page size
    canvasComparison->SetLogy(1);
    topMargin = (nTopRows*54.)/(846.+nTopRows*54.);
    DrawGammaCanvasSettings(canvasComparison, leftMar, rightMar, topMargin, bottomMargin);

    TLegend *legendComparison = new TLegend(0.12,0.99 - (nTopRows*0.04),0.95,0.98);
    legendComparison->SetNColumns(4);
    legendComparison->SetFillColor(0);
    legendComparison->SetLineColor(0);
    legendComparison->SetTextSize(0.03);
    legendComparison->SetTextFont(42);

    for(Int_t h=0; h<(Int_t)vecHistos[0].size(); h++){
        cout << h << ", ";
	if ( vecHistosNameForSaving[0].at(h).Contains("VertexZ") || vecHistosNameForSaving[0].at(h).Contains("EventPlaneAngle")) scaleFactor = 1.2;
        AdjustHistRange(vecHistos,scaleFactor,scaleFactor,h,nSets,kFALSE);
	if ( vecHistosNameForSaving[0].at(h).Contains("V0Mult") ) ((TH1D*) vecHistos[0].at(h))->GetXaxis()->SetRangeUser(0,(((TH1D*) vecHistos[0].at(h))->GetXaxis()->GetXmax()));
        for(Int_t i=0; i<nSets; i++) {  // Draw Data histograms
            TString draw;
            if(h==0) draw = (i==0)?"p":"p, same";
            else draw = (i==0)?"px0e1":"px0e1, same";
            ((TH1D*) vecHistos[i].at(h))->Draw(draw.Data());
            legendComparison->AddEntry(((TH1D*) vecHistos[i].at(h)),Form("Data %s", centralities[i].Data()),"p");
        }
       for(Int_t i=0; i<nSets; i++) {  // Draw MC histograms
            TString draw;
            if(h==0) draw = "p, same";
            else draw = "px0e1, same";
            ((TH1D*) vecHistosMC[i].at(h))->Draw(draw.Data());
            legendComparison->AddEntry(((TH1D*) vecHistosMC[i].at(h)),Form("MC %s", centralities[i].Data()),"p");
        }
        legendComparison->Draw();
        if ( vecHistosNameForSaving[0].at(h).Contains("Candidates") )
            PutProcessLabelAndEnergyOnPlot(0.75, 0.97-nTopRows*0.06, 0.03, fCollisionSystem.Data(),"","");
        else
            PutProcessLabelAndEnergyOnPlot(0.75, 0.97-nTopRows*0.06, 0.03, fCollisionSystem.Data(), "", "");
	if ( vecHistosNameForSaving[0].at(h).Contains("VertexZ") || vecHistosNameForSaving[0].at(h).Contains("EventPlaneAngle")) SaveWriteCanvas(canvasComparison, Form("EventQA_CompareCentralities/Comparison/%s_%s/Periodwise/%s.%s", dataSet.Data(), MCSet.Data(), vecHistosNameForSaving[0].at(h).Data(),suffix.Data()));
        else SaveWriteCanvas(canvasComparison, Form("EventQA_CompareCentralities/Comparison/%s_%s/Periodwise/%s.%s", dataSet.Data(), MCSet.Data(), vecHistosNameForSaving[0].at(h).Data(),suffix.Data()), kFALSE, kTRUE);
	legendComparison->Clear();
        // if(h==9) TGaxis::SetExponentOffset(0, 0, "x");
    }
    cout << endl;

//*************************************************************************************************************

    fOutput->Write();
    fOutput->Close();

    cout << "Done with EventQA_CompareCentralities" << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;


    return;
}
