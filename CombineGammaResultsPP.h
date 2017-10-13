// define number of sets
    Int_t nSets                                 = 13;
    Int_t nSetsMC                               = 4;
TString outputDir;
TString suffix = "pdf";
    TDirectory* directoryPi0[3];
    TDirectory* directoryPi02[3];
    TDirectory* directoryEta[3];
    
    Width_t widthLinesBoxes                     = 1.4;
    Width_t widthCommonFit                      = 2;

    
    TString combinatorics[17]                               = { "Elec+Elec","Elec+Pion","Elec+Kaon","Elec+Proton","Elec+Muon","Pion+Pion","Pion+Kaon","Pion+Proton",
                                                                "Pion+Muon","Kaon+Kaon","Kaon+Proton","Kaon+Muon","Proton+Proton","Proton+Muon","Muon+Muon","Rest","All"};
    TString combinatoricsCalo[11]                           = { "Electron","Pion","Proton","Kaon","Neutron","K0s","Lambda","Muon","K0l","Rest","All"};
    
    Color_t colorsCombinatorics[17]                         = { kAzure-1, kRed+1, kOrange+7, kMagenta+2, kRed-9,
                                                                kBlue-6, kAzure+5, kPink, kCyan-3, kGreen-3,
                                                                kSpring+9, kGreen+2, kBlue+2, kMagenta-6, kSpring+4,
                                                                kCyan+2, 809};
    Color_t colorsCombinatoricsCalo[11]                     = { kAzure-1, kRed+1, kOrange+7, kMagenta+2, kRed-9, kBlue,
                                                                kBlue-6, kAzure+5, kPink, kCyan-3, kGreen-3};
    
    Style_t markersCombinatorics[17]                        = { 20, 21, 24, 25, 27,
                                                                28, 29, 30, 33, 34, 
                                                                20, 21, 24, 25, 27,
                                                                28, 29};
    Style_t markersCombinatoricsCalo[11]                    = { 20, 21, 24, 25, 27,
                                                                28, 29, 30, 33, 34, 20};
    
    TString decays[9]                                       = { "Pi0","Eta","Etap","Omega","Rho",
                                                                "Phi","Sigma", "All decays","direct #gamma"
                                                              };
    TString decaysLabels[9]                                  = { "#gamma from #pi^{0}","#gamma from #eta","#gamma from #eta'","#gamma from #omega","#gamma from #rho^{0}",
                                                                "#gamma from #phi", "#gamma from #Sigma^{0}", "All decays","direct #gamma"
                                                              };

    Color_t colorsDecay[9]                                  = { kRed+1, kAzure-1, kOrange+7, kGreen+2, kCyan-3, 809, 
                                                                kBlue+2, kBlack, kRed-1
                                                              };
    TH1D*  histoCombinatorialSpecies_Pt[13][17];
    TH1D*  histoSignalToCombBackgroundRatio[13][17];
    Color_t colorComb                               = kBlack;
    Style_t markerStyleComb                         = 20;
    Size_t markerSizeComb                           = 1.2;
    // DATA: Definition of colors, styles and markers sizes
    Float_t nEvtMC[13];
    TH1F*   histoEventQualityMC[13];
    TH1D*   BinShiftRatio[13];
    TH1D* TruePrimaryConvGamma_Pt[13];
    TH1D* histoGammaPurity_Pt[13];
    TH1D* GammaRecoEff_MCPt[13];
    TH1D* GammaRecoEff_WithResolCorr_Pt[13];
    TH1D* GammaConvProb_Pt[13];
    TH1D* PileUpCorrectionFactor[13];
    
    TH1D* ESDGammaDCAzAllMC[13];
    TH1D* ESDGammaDCAzAll[13];
    TH1D* ESDGammaDCAzBack[13];
    TH1D* ESDGammaDCAzAllSub[13];
    TH1D* ESDGammaDCAzAll3MC[13];
    TH1D* ESDGammaDCAzAll3[13];
    TH1D* ESDGammaDCAzBack3[13];
    TH1D* ESDGammaDCAzAllSub3[13];

    
    TH1F* histoPythia8InvXSectionPi08TeV;
    TH1F* histoPythia8InvXSectionPi0276TeV;
        TH1F* histoPythia8InvXSectionPi07TeV;
        TH1F* histoPythia8InvXSectionPi0900GeV;      
    
    TH1D* histoDirectPhotonSpectrum[13];
    TGraphAsymmErrors* graphDirectPhotonSpectrum[13];
    TGraphAsymmErrors* graphNLODirGamma[13];
    
    TH1D* DoubleRatioConversionTrueEffPurity[13];
    TGraphAsymmErrors* graphDoubleRatioConversionTrueEffPurity[13];
    TGraphAsymmErrors* DoubleRatioConversionTrueEffPurity_SystErr[13];
    
    TH1D* histoGammaSpecCorrPurity[13];
    TGraphAsymmErrors* graphGammaSpecCorrPurity[13];
    TGraphAsymmErrors* ESD_TruePrimaryConvGammaESD_PtMCPt_SystErr[13];
    
    TH1D* IncRatioPurity_trueEff[13];
    TGraphAsymmErrors* graphIncRatioPurity_trueEff[13];
    TGraphAsymmErrors* IncRatioPurity_trueEff_SystErr[13];
    
    TH1D* DoubleRatioStatError;
    TGraphAsymmErrors* graphDoubleRatioStatError;
    TGraphAsymmErrors* DoubleRatioSystError;
    
    TH1D* EtaToPi0RatioConversionBinShifted;
    TGraphAsymmErrors* graphEtaToPi0RatioConversionBinShifted;
    TGraphAsymmErrors* EtaToPi0RatioConversionBinShiftedSys;
    
    TH1D* SecondaryGammaFromXFromK0sRecoEff_Pt[13];
    TH1D* SecondaryGammaFromXFromK0lRecoEff_Pt[13];
    TH1D* SecondaryGammaFromXFromLambdaRecoEff_Pt[13];
    TH1D* SecondaryGammaFromXFromRestRecoEff_Pt[13];
    
    TH1D* SecondaryGammaFromXFromK0sConvProb_MCPt[13];
    TH1D* SecondaryGammaFromXFromK0lConvProb_MCPt[13];
    TH1D* SecondaryGammaFromXFromLambdaConvProb_MCPt[13];
    TH1D* SecondaryGammaFromXFromRestConvProb_MCPt[13];
    
    TH1D* GammaRaw_Pt[13];
    TH1D* histoGammaTrueSecConvGammaFromXFromK0s_Cocktail_Raw_Pt[13];
    TH1D* histoGammaTrueSecConvGammaFromXFromK0l_Cocktail_Raw_Pt[13];
    TH1D* histoGammaTrueSecConvGammaFromXFromLambda_Cocktail_Raw_Pt[13];
    TH1D* histoGammaTrueSecCocktailGammaRest_Pt[13];

    TH1D* histoFracAllGammaToSecFromXFromK0s_Cocktail_Pt[13];
    TH1D* histoFracAllGammaToSecFromXFromK0l_Cocktail_Pt[13];
    TH1D* histoFracAllGammaToSecFromXFromLambda_Cocktail_Pt[13];
    TH1D* histoFracAllGammaToSecRest_Cocktail_Pt[13];
    
    TH1D* hNEventsFracGoodEvents[13];
    TH1D* hPi0Mass[13];
    TH1D* hEtaMass[13];
    TH1D* hTracksGoodMean[13];
    TH1D* hConvNCandidatesQA[13];
    TH1D* hFracPileup[13];
    TH1D* hFracWVtxOutside10cm[13];
    TH1D* hPi0Frac[13];
    TH1D* hVertexZMean[13];

    TH1D* hNEventsFracGoodEventsMC[13];
    TH1D* hPi0MassMC[13];
    TH1D* hEtaMassMC[13];
    TH1D* hTracksGoodMeanMC[13];
    TH1D* hConvNCandidatesQAMC[13];
    TH1D* hFracPileupMC[13];
    TH1D* hFracWVtxOutside10cmMC[13];
    TH1D* hPi0FracMC[13];
    TH1D* hVertexZMeanMC[13];
    
    TH1D* TrueMesonEffiPt[10];
    TH1D* TrueSecFromK0SEffiPt[10];
    TH1D* TrueSecFromK0LEffiPt[10];
    TH1D* TrueSecFromLambdaEffiPt[10];
    TH1D* TrueSecFromRestEffiPt[10];
    
    TH1D* fMCMesonAcceptPt[10];
    TH1D* fMCSecPi0FromK0SAccepPt[10];
    TH1D* fMCSecPi0FromK0LAccepPt[10];
    TH1D* fMCSecPi0FromLambdaAccepPt[10];
    TH1D* fMCSecPi0FromRestAccepPt[10];
    
    TH1D* histoYieldMesonPerEvent[10];
    TH1D* SecYieldFromK0SMesonFromCocktail[10];
    TH1D* SecYieldFromK0LMesonFromCocktail[10];
    TH1D* SecYieldFromLambdaMesonFromCocktail[10];
    TH1D* SecYieldFromRestMeson[10];
    
    TH1D* RatioToRaw[10];
    TH1D* RatioSecYieldFromK0SMesonFromCocktailToRaw[10];
    TH1D* RatioSecYieldFromK0LMesonFromCocktailToRaw[10];
    TH1D* RatioSecYieldFromLambdaMesonFromCocktailToRaw[10];
    TH1D* RatioSecYieldFromRestMesonToRaw[10];
    
    TH1D* histoNumberOfEvents[10];
    TH1D* histoPi0RawYields[10];
    TH1D* histoEtaRawYields[10];
    TH1D* histoPi0Mass[10];
    TH1D* histoPi0FWHMMeV[10];
    TH1D* histoPi0TrueMass[10];
    TH1D* histoPi0TrueFWHMMeV[10];
    TH1D* histoEtaMass[10];
    TH1D* histoEtaFWHMMeV[10];
    TH1D* histoEtaTrueMass[10];
    TH1D* histoEtaTrueFWHMMeV[10];
    TH1D* histoPi0AccEff[10];
    TH1D* histoPi0Acc[10];
    TH1D* histoPi0TrueEffPt[10];
    TH1D* histoEtaAcc[10];
    TH1D* histoEtaAccEff[10];
    TH1D* histoEtaTrueEffPt[10];
    TH1D* histoPi0AccTimesEff[10];
    TH1D* histoPi0InvCrossSection[10];
    TH1D* BGEstimateFromPileup[10];
    TH1D* PileupContamination[10];
    TH1D* fHistFracCat_1_vsPt[10];
    TH1D* fHistFracCat_2_vsPt[10];
    TH1D* fHistFracCat_3_vsPt[10];
    TH1D* fHistFracCat_4_vsPt[10];
    TH1D* fHistFracCat_5_vsPt[10];
    TH1D* fHistFracCat_6_vsPt[10];
    TH1D* fHistFracIntHistBGvsPt_Cat_1_Variant_1[10];
    TH1D* fHistFracIntHistBGvsPt_Cat_2_Variant_1[10];
    TH1D* fHistFracIntHistBGvsPt_Cat_3_Variant_1[10];
    TH1D* fHistFracIntHistBGvsPt_Cat_4_Variant_1[10];
    TH1D* fHistFracIntHistBGvsPt_Cat_5_Variant_1[10];
    TH1D* fHistFracIntHistBGvsPt_Cat_6_Variant_1[10];
    TH1D* HistDCAZUnderMesonCat_1_MesonPt_0304MCTRUE[10];
    TH1D* HistDCAZUnderMesonCat_1_MesonPt_0304MC[10];
    TH1D* HistDCAZUnderMesonCat_1_MesonPt_0304[10];
    TH1D* fHistDCAZUnderMesonBGEstimateCat_1_MesonPt_0304[10];
    TH1D* histoChargedPionStat[10];
    TH1D* histoChargedPionSys[10];
    TH1D* fDCAEventQualityMC[10];
    Int_t fNEventsDCAMC[10];
    TH1D* fDCAEventQuality[10];
    Int_t fNEventsDCA[10];
    TGraphAsymmErrors* graphChargedPionStat[10];
    TGraphAsymmErrors* graphChargedPionSys[10];
    TGraphAsymmErrors* histCorrectedYieldPi0[10];
    TGraphAsymmErrors* graphCorrectedYieldPi0[10];
    TGraphAsymmErrors* Pi0SystError[10];
    TGraphAsymmErrors* graphPi0InvCrossSectionSys[10];
    TGraphAsymmErrors* graphPi0InvCrossSectionStat[10];
    TGraphAsymmErrors* graphEtaInvCrossSectionSys[10];
    TGraphAsymmErrors* graphEtaInvCrossSectionStat[10];
    TString strTrigName[3]                       = {"INT1", "INT1", "MB"};

    TH1D* histoEtaToPi0Stat[10];
    TGraphAsymmErrors* graphEtaToPi0Stat[10]        = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TGraphAsymmErrors* graphEtaToPi0Sys[10]         = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TGraphAsymmErrors* graphEtaToPi0ShiftStat[10]        = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TGraphAsymmErrors* graphEtaToPi0ShiftSys[10]         = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};

    Double_t rapidityMeas[3]                       = {1.6, 1.6,1.6};


Style_t markerstyles[4]                         ={34,20,33,29};
    Style_t markerstylesMC[4]                   ={28,24,27,30};
    Size_t markersize[4]                        ={2.3,2.3,3.4,3.4};
    Size_t markersizeMC[4]                      ={2.3,2.3,3.4,2.2};

    Color_t colorData[4]                        ={kRed+2,kBlue+2,kGreen+2,kMagenta+2};
    Color_t colorMC[4]                          ={kRed-8,kBlue-8,kGreen-8,kMagenta-8};
    
    TString nameMeasGlobal[11]                  = {"#sqrt{s} = 900 GeV", "#sqrt{s} = 7 TeV", "#sqrt{s} = 8 TeV", "#sqrt{s} = 2.76 TeV"};

    TString nameDataset[13]                     = { "LHC10c_900GeV",
                                                    "LHC10b", "LHC10c", "LHC10d", "LHC10e", "LHC10f",
                                                    "LHC12a", "LHC12b", "LHC12c", "LHC12d", "LHC12f", "LHC12h", "LHC12i"};
    TString nameDatasetAdd[13]                  = { "",
                                                    "_pass4", "_pass4", "_pass4", "_pass4", "_pass4",
                                                    "", "", "", "", "", "", ""};
    TString fEnergyFlag[13]                     = { "900GeV",
                                                    "7TeV", "7TeV", "7TeV", "7TeV", "7TeV",
                                                    "8TeV", "8TeV", "8TeV", "8TeV", "8TeV", "8TeV", "8TeV"};
    
    // MC: Definition of colors, styles and markers sizes
    TString nameDatasetMCLoad[4]                = { "LHC14j4c_900GeV",
                                                    "LHC14j4",
                                                    "LHC15h1", "LHC15h2"};
    TString nameDatasetMC[4]                    = { "Pythia6_900GeV",
                                                    "Pythia6",
                                                    "Pythia8", "Phojet"};
    TString fEnergyFlagMC[4]                    = { "900GeV",
                                                    "7TeV",
                                                    "8TeV", "8TeV"};