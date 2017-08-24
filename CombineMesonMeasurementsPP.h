

void drawLatexAdd(TString latextext, Double_t textcolumn, Double_t textrow, Double_t textSizePixel,Bool_t setFont = kFALSE, Bool_t setFont2 = kFALSE, Bool_t alignRight = kFALSE, Color_t textcolor = kBlack){
    TLatex *latexDummy                              = new TLatex(textcolumn ,textrow,latextext);
    SetStyleTLatex( latexDummy, textSizePixel,4);
    if(setFont)
        latexDummy->SetTextFont(62);
    if(setFont2)
        latexDummy->SetTextFont(43);
    if(alignRight)
        latexDummy->SetTextAlign(31);
    latexDummy->SetTextColor(textcolor);
    latexDummy->Draw();
}
void drawMarker(TH1D* histoDummy, Double_t column, Double_t row, Double_t markerScale){
    TMarker* markerdummy;
    markerdummy                                     = CreateMarkerFromHisto(histoDummy,column ,row ,markerScale);
    markerdummy->DrawMarker(column ,row);
}

//______________________________________________ Definition of files and directories
    TFile* inputFile[6];
    TDirectory* directoryPi0[6];
    TDirectory* directoryEta[6];

//______________________________________________ Definition of strings
    TString energyLatex[6]                      = {"#sqrt{s} = 900 GeV", "#sqrt{s} = 2.76 TeV", "#sqrt{s} = 5 TeV", "#sqrt{s} = 7 TeV", "#sqrt{s} = 8 TeV", "#sqrt{s} = 13 TeV"};
    TString nameMeasGlobal[11]                  = {"PCM", "PHOS", "EMCal", "PCM-PHOS", "PCM-EMCal", "PCM-Dalitz", "PHOS-Dalitz", "EMCal-Dalitz", "EMCal high pT", "EMCal merged", "Comb"};
    TString nameMeasGlobalshort[11]             = {"PCM", "PHOS", "EMCAL", "PCMPHOS", "PCMEMCAL", "PCM-Dalitz", "PHOS-Dalitz", "EMCal-Dalitz", "EMCal high pT", "EMCal merged", "Comb"};
    TString graphNameModifier[11]               = {"", "", "", "", "", "", "", "", "", "", "A"};
    TString nameEnergyGlobal[6]                 = {"900GeV", "2.76TeV",   "5TeV",   "7TeV",   "8TeV",  "13TeV"};
    TString nameEnergyGlobal2[6]                = {"900GeV", "2760GeV",   "5TeV",   "7TeV",   "8TeV",  "13TeV"};
    TString nameGeneratorGlobal[6]              = {"Pythia",  "Pythia", "Pythia", "Pythia", "Pythia", "Pythia"};
    TString nameGeneratorGlobal2[6]             = {"Phojet",  "Phojet", "Phojet", "Phojet", "Phojet", "Phojet"};

//______________________________________________ Definition of colors and markers
    Color_t colorMeas[11];
    Color_t colorMeasMC[11];
    Style_t markerStyleMeas[11];
    Style_t markerStyleMeasMC[11];
    Size_t  markerSizeMeas[11];
    Size_t  markerSizeMeasMC[11];

    Color_t colorEnergy[6];
    Color_t colorEnergyMC[6];
    Color_t colorEnergyMC2[6];
    Color_t colorEnergyBox[6];
    Style_t markerStyleEnergy[6];
    Style_t markerStyleEnergyMC[6];
    Size_t  markerSizeEnergy[6];
    Size_t  markerSizeEnergyMC[6];

    Color_t pythia8color = kRed+2;
    Color_t colorNLO     = kAzure-4;

    Width_t widthLinesBoxes                     = 1.4;
    Width_t widthCommonFit                      = 1.7;

void LoadColorsMarkersAndSizes(){
    for (Int_t i = 0; i < 11; i++){
        colorMeas[i]                             = GetDefaultColorDiffDetectors(        nameMeasGlobal[i].Data(), kFALSE, kFALSE, kTRUE );
        colorMeasMC[i]                           = GetDefaultColorDiffDetectors(        nameMeasGlobal[i].Data(), kTRUE, kFALSE, kTRUE  );
        markerStyleMeas[i]                       = GetDefaultMarkerStyleDiffDetectors(  nameMeasGlobal[i].Data(), kFALSE                );
        markerStyleMeasMC[i]                     = GetDefaultMarkerStyleDiffDetectors(  nameMeasGlobal[i].Data(), kTRUE                 );
        markerSizeMeas[i]                        = GetDefaultMarkerSizeDiffDetectors(   nameMeasGlobal[i].Data(), kFALSE                );
        markerSizeMeasMC[i]                      = GetDefaultMarkerSizeDiffDetectors(   nameMeasGlobal[i].Data(), kTRUE                 );
    }

    for (Int_t i = 0; i < 6; i++){
        colorEnergy[i]                          = GetColorDefaultColor( nameEnergyGlobal[i].Data(),                              "", "", kFALSE  );
        colorEnergyMC[i]                        = GetColorDefaultColor( nameEnergyGlobal[i].Data(),   nameGeneratorGlobal[i].Data(), "", kFALSE  );
        colorEnergyMC2[i]                       = GetColorDefaultColor( nameEnergyGlobal[i].Data(),  nameGeneratorGlobal2[i].Data(), "", kFALSE  );
        colorEnergyBox[i]                       = GetColorDefaultColor( nameEnergyGlobal[i].Data(),                              "", "", kTRUE   );
        markerStyleEnergy[i]                    = GetDefaultMarkerStyle(nameEnergyGlobal[i].Data(),   "","");
        markerStyleEnergyMC[i]                  = GetDefaultMarkerStyle(nameEnergyGlobal[i].Data(), "MC","");
        markerSizeEnergy[i]                     = GetDefaultMarkerSize( nameEnergyGlobal[i].Data(),   "","");
        markerSizeEnergyMC[i]                   = GetDefaultMarkerSize( nameEnergyGlobal[i].Data(), "MC","");
    }
}

 void getPadTextSizes(TPad* padDummy,Double_t textSizeLabelsPixel, Double_t textsizeLabelDummy, Double_t textsizeFracDummy)
{
    if (padDummy->XtoPixel(padDummy->GetX2()) < padDummy->YtoPixel(padDummy->GetY1())){
        textsizeLabelDummy                  = (Double_t)textSizeLabelsPixel/padDummy->XtoPixel(padDummy->GetX2()) ;
        textsizeFracDummy                   = (Double_t)1./padDummy->XtoPixel(padDummy->GetX2()) ;
        cout << "here " << textsizeFracDummy << "\t" << padDummy->GetX2() << "\t" << padDummy->XtoPixel(padDummy->GetX2()) << endl;
    } else {
        textsizeLabelDummy                  = (Double_t)textSizeLabelsPixel/padDummy->YtoPixel(padDummy->GetY1());
        textsizeFracDummy                   = (Double_t)1./padDummy->YtoPixel(padDummy->GetY1());
        cout << "now here " << textsizeFracDummy << "\t" << padDummy->GetY1() << "\t" << padDummy->YtoPixel(padDummy->GetY1()) << endl;
    }
}

