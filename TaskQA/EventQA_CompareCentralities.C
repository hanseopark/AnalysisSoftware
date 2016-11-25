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
    const Int_t nSets           = 2;
    Size_t constMarkerSize      = 1;
    TString fEnergyFlag         = "PbPb_5.02TeV";
    TString period              = "LHC15o";                                 // for GetDefaultMarkerStyle and legend
    TString dataSet             = "LHC15o_HighIR_pass1_AOD";                // used for filename 
    //TString centralities[nSets] = {"0-10%","10-20%","20-50%","50-90%"};     // used for plot and marker style and color
    TString centralities[nSets] = {"0-10%","50-90%"};  
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

    std::vector<TH1D*> vecHistos[nSets];
    std::vector<TString> vecHistosName[nSets];
    std::vector<TString> vecHistosNameForSaving[nSets];

    TH1D* hcVertexZ[nSets];
    TH1D* hcGoodESDTracks[nSets];
    TH1D* hcV0Mult[nSets];
    TH1D* hcGammaCandidates[nSets];
    TH1D* hcCentrality[nSets];
    TH1D* hcEventPlaneAngle[nSets];

    Int_t nRange = 0;

    for(Int_t i=0; i<nSets; i++){
        
        vecHistos[i].push_back(hcVertexZ[i]);
        vecHistosName[i].push_back("VertexZ");
        vecHistosNameForSaving[i].push_back("VertexZ");
        if(i==0)nRange++;

        vecHistos[i].push_back(hcGoodESDTracks[i]);
        vecHistosName[i].push_back("GoodESDTracks");
        vecHistosNameForSaving[i].push_back("GoodESDTracks");
        if(i==0)nRange++;

        
        vecHistos[i].push_back(hcV0Mult[i]);
        vecHistosName[i].push_back("V0 Multiplicity");
        vecHistosNameForSaving[i].push_back("V0Mult");
        if(i==0)nRange++;
    
        vecHistos[i].push_back(hcGammaCandidates[i]);
        vecHistosName[i].push_back("GammaCandidates");
        vecHistosNameForSaving[i].push_back("GammaCandidates");
        if(i==0)nRange++;

        vecHistos[i].push_back(hcEventPlaneAngle[i]);
        vecHistosName[i].push_back(Form("EventPlaneAngle %s",fEventCutSelection[i].Data()));
        vecHistosNameForSaving[i].push_back("EventPlaneAngle");
        if(i==0)nRange++;
       
    }

    Double_t nEvents[nSets];

    for(Int_t i=0; i<nSets; i++)    {
        TString fFile = Form("%s/EventQA_%s.root", pathDataSets[i].Data(), dataSet.Data());
        TFile* File = new TFile(fFile.Data(),"READ");

        if(File->IsZombie()) {cout << "Warning: ROOT file '" << fFile.Data() << "' could not be openend!" << endl; return;}

        nEvents[i] = GetNEvents((TH1*) File->Get("NEvents"),kFALSE);

        cout << endl;
        cout << "\t\t----------------------------------------------------------------------------" << endl;
        cout << "\t\tDataSet: " << dataSet.Data() << ", " << centralities[i].Data() << ":" << endl;
        cout << "\t\tTotal NEvents: " << nEvents[i] << endl;
        cout << "\t\tProcessing file: " << fFile.Data() << endl;
        cout << "\t\t----------------------------------------------------------------------------" << endl;
        cout << endl;

        cout << "\n\t\t----------------------------------------------------------------------------" << endl;
        cout << "\t\t\tNumber of Histograms being processed: " << nRange << endl;
        cout << "\t\t----------------------------------------------------------------------------" << endl;

        for(Int_t j=0; j<nRange; j++)
        {
            cout << ((TString)vecHistosName[i].at(j)).Data();
            TH1D* temp = (TH1D*) File->Get(((TString)vecHistosName[i].at(j)).Data());
            cout << ", ";
            temp->Sumw2();
	    // scale histograms with number of entries:
	    Double_t tempEntries = temp->GetEntries();
	    temp->Scale(1./tempEntries);
	    temp->GetYaxis()->SetTitle(Form("#frac{1}{N} %s",temp->GetYaxis()->GetTitle()));
	    // with number of events: 
	    // temp->Scale(1./nEvents[i]);
	    // temp->GetYaxis()->SetTitle(Form("#frac{1}{N_{Events}} %s",temp->GetYaxis()->GetTitle()));
	    OnlyEditTH1(temp, hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
            temp->GetXaxis()->SetTitleOffset(1);
            temp->GetYaxis()->SetTitleOffset(1.1);
            temp->SetTitle("");
            vecHistos[i].at(j) = new TH1D(*temp);
        }

        cout << "done!" << endl;
    }

    gSystem->Exec("mkdir -p EventQA_CompareCentralities/"+dataSet);
    gSystem->Exec("mkdir -p EventQA_CompareCentralities/"+dataSet+"/Periodwise");

    TString nameOutput = Form("EventQA_CompareCentralities/%s/Output_Periodwise.root",dataSet.Data());
    TFile* fOutput = new TFile(nameOutput,"RECREATE");
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "Output file: " << nameOutput << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

    Int_t nTopRows = (nSets-1)/(3) + 1;
    TCanvas* canvas = new TCanvas("canvas","",200,10,1350,(846.+nTopRows*54.));  // gives the page size
    canvas->SetLogy(1);
    Double_t leftMar = 0.09; Double_t rightMar = 0.02; Double_t topMargin = (nTopRows*54.)/(846.+nTopRows*54.); Double_t bottomMargin = 0.09;
    DrawGammaCanvasSettings(canvas, leftMar, rightMar, topMargin, bottomMargin);


//****************************** Combined Comparison Histograms ************************************************
    TLegend *legend = new TLegend(0.12,0.99 - (nTopRows*0.04),0.95,0.98);
    legend->SetNColumns(3);
    legend->SetFillColor(0);
    legend->SetLineColor(0);
    legend->SetTextSize(0.03);
    legend->SetTextFont(42);

    canvas->SetLogy(1);
    for(Int_t h=0; h<(Int_t)vecHistos[0].size(); h++){
        cout << h << ", ";    
        AdjustHistRange(vecHistos,5,5,h,nSets,kFALSE);
        for(Int_t i=0; i<nSets; i++) {
            TString draw;
            if(h==0) draw = (i==0)?"p":"p, same";
            else draw = (i==0)?"px0e1":"px0e1, same";
	    canvas->SetLogy(1);
            ((TH1D*) vecHistos[i].at(h))->Draw(draw.Data());
            legend->AddEntry(((TH1D*) vecHistos[i].at(h)),centralities[i],"p");
        }
	canvas->SetLogy(1);
        legend->Draw();
        if ( vecHistosNameForSaving[0].at(h).Contains("Candidates") ) 
            PutProcessLabelAndEnergyOnPlot(0.75, 0.97-nTopRows*0.06, 0.03, fCollisionSystem.Data(),"","");
        else 
            PutProcessLabelAndEnergyOnPlot(0.75, 0.97-nTopRows*0.06, 0.03, fCollisionSystem.Data(), "", "");
	canvas->SetLogy(1);
        SaveWriteCanvas(canvas, Form("EventQA_CompareCentralities/%s/Periodwise/%s.%s", dataSet.Data(), vecHistosNameForSaving[0].at(h).Data(),suffix.Data()), kFALSE, kTRUE);
        legend->Clear();
        // if(h==9) TGaxis::SetExponentOffset(0, 0, "x");
    }
    cout << endl;

    fOutput->Write();
    fOutput->Close();

    cout << "Done with EventQA_CompareCentraliries" << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;


    return;
}
