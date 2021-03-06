// provided by Gamma Conversion Group, $ALICE_ROOT/PWG4/GammaConv ;https://twiki.cern.ch/twiki/bin/view/ALICE/PWG4GammaConversion

Int_t       fMode                                                   = 0;
Int_t       fIsMC                                                   = 0;
Double_t    fNEvnt                                                  = 0;
Bool_t      fEnablePCM                                              = 0;
Bool_t      fEnableCalo                                             = 0;
Bool_t      fUseAnalysisBinning                                     = 0;
Bool_t      fUseFitConstraints                                      = kTRUE;

TString     fDate                                                   = "";
TString     fDirectPhoton                                           = "";

TString     fOutputDir                                              = "";
TString     fSuffix                                                 = "";
TString     fMeson                                                  = "";

TString     fEventCutNumber                                         = "";
TString     fCutSelection                                           = "";
TString     fGammaCutNumber                                         = "";
TString     fClusterCutNumber                                       = "";
TString     fElectronCutNumber                                      = "";
TString     fMesonCutNumber                                         = "";

TString     fCutSelectionRead                                       = "";
TString     fEventCutSelection                                      = "";
TString     fGammaCutSelection                                      = "";
TString     fMesonCutSelection                                      = "";

TString     fEnergyFlag                                             = "";
TString     fPrefix                                                 = "";
TString     fPrefix2                                                = "";
TString     fPeriodFlag                                             = "";
TString     fTextMeasurement                                        = "";
TString     fCollisionSystem                                        = "";
TString     fDetectionProcess                                       = "";


Double_t NEntriesData;
Int_t status;
Double_t signalPurityErr;
Double_t signalPurityNumErr;
Double_t signalPurityDenomErr;

Double_t nSigmaLow  = -20.0;
Double_t nSigmaHigh = 20.0;
Double_t constrainLow  = 0.;
Double_t constrainHigh = 1.;

TString SigmaStarForm = "#Kappa = #frac{|#kappa^{+}|+|#kappa^{-}|}{2}+2(#kappa^{+}+#kappa^{-})";

TH1D*       fDeltaPt                                              = NULL;
// 2D-plots
TH2D*       hKappaTPCPt                                           = NULL;
TH2D*       hKappaTPCPtElEl                                       = NULL;
TH2D*       hKappaTPCPtElPi                                       = NULL;
TH2D*       hKappaTPCPtPiPi                                       = NULL;
TH2D*       hKappaTPCPtRest                                       = NULL;
// projections
TH1D*       hKappaTPC                                             = NULL;
TH1D**      hKappaTPCAfterCut                                     = NULL;
TH1D**      hKappaTPCElEl                                         = NULL;
TH1D**      hKappaTPCElPi                                         = NULL;
TH1D**      hKappaTPCPiPi                                         = NULL;
TH1D**      hKappaTPCRest                                         = NULL;
TH1D**      hKappaTPCElK                                          = NULL;
TH1D**      hKappaTPCElP                                          = NULL;
TH1D**      hKappaTPCPiK                                          = NULL;
TH1D**      hKappaTPCPiP                                          = NULL;
TH1D**      hKappaTPCKK                                           = NULL;
TH1D**      hKappaTPCHad                                          = NULL;
TH1D**      hKappaTPCLeftoverRest                                 = NULL;
TH1D**      hKappaTPCPartialRest                                  = NULL;
TH1D**      hKappaTPCElPK                                         = NULL;
TH1D**      hKappaTPCPiPK                                         = NULL;

TH1D**      hKappaTPCSum                                          = NULL;
TH1D**      hKappaTPCSumBkg                                       = NULL;
TH1D**      hKappaTPCTotalSum                                     = NULL;

TH1D*       hValuesElEl                                           = NULL;
TH1D*       hFractionElEl                                         = NULL;
TH1D*       hValuesElPi                                           = NULL;
TH1D*       hFractionElPi                                         = NULL;
TH1D*       hValuesPiPi                                           = NULL;
TH1D*       hFractionPiPi                                         = NULL;
TH1D*       hValuesRest                                           = NULL;
TH1D*       hFractionRest                                         = NULL;
TH1D*       hValuesElPK                                           = NULL;
TH1D*       hFractionElPK                                         = NULL;
TH1D*       hValuesPiPK                                           = NULL;
TH1D*       hFractionPiPK                                         = NULL;
TH1D*       hValuesPartialRest                                    = NULL;
TH1D*       hFractionPartialRest                                  = NULL;

TH1D*       hSignalPurity1                                        = NULL;
TH1D*       hSignalPurity2                                        = NULL;
TH1D*       hSignalPurity3                                        = NULL;
TH1D*       hBackgroundPiPi1                                      = NULL;
TH1D*       hBackgroundPiPi2                                      = NULL;
TH1D*       hBackgroundPiPi3                                      = NULL;
TH1D*       hBackgroundElPi1                                      = NULL;
TH1D*       hBackgroundElPi2                                      = NULL;
TH1D*       hBackgroundElPi3                                      = NULL;
TH1D*       hBackgroundRest1                                      = NULL;
TH1D*       hBackgroundRest2                                      = NULL;
TH1D*       hBackgroundRest3                                      = NULL;
TH1D*       hBackgroundElPK1                                      = NULL;
TH1D*       hBackgroundElPK2                                      = NULL;
TH1D*       hBackgroundElPK3                                      = NULL;
TH1D*       hBackgroundPiPK1                                      = NULL;
TH1D*       hBackgroundPiPK2                                      = NULL;
TH1D*       hBackgroundPiPK3                                      = NULL;
TH1D*       hBackgroundPartialRest1                               = NULL;
TH1D*       hBackgroundPartialRest2                               = NULL;
TH1D*       hBackgroundPartialRest3                               = NULL;

TH1D**      hTemplateElEl                                         = NULL;
TH1D**      hTemplateElPi                                         = NULL;
TH1D**      hTemplatePiPi                                         = NULL;
TH1D**      hTemplateRest                                         = NULL;
TH1D**      hTemplateSum                                          = NULL;
TH1D**      hTemplateElPK                                         = NULL;
TH1D**      hTemplatePiPK                                         = NULL;
TH1D**      hTemplatePartialRest                                  = NULL;

Double_t signalPurity1;
Double_t signalPurityNum1;
Double_t signalPurityDenom1;
Double_t signalPurity2;
Double_t signalPurityNum2;
Double_t signalPurityDenom2;
Double_t signalPurity3;
Double_t signalPurityNum3;
Double_t signalPurityDenom3;
Double_t backgroundPiPi1;
Double_t backgroundPiPiNum1;
Double_t backgroundPiPi2;
Double_t backgroundPiPiNum2;
Double_t backgroundPiPi3;
Double_t backgroundPiPiNum3;
Double_t backgroundElPi1;
Double_t backgroundElPiNum1;
Double_t backgroundElPi2;
Double_t backgroundElPiNum2;
Double_t backgroundElPi3;
Double_t backgroundElPiNum3;
Double_t backgroundRest1;
Double_t backgroundRestNum1;
Double_t backgroundRest2;
Double_t backgroundRestNum2;
Double_t backgroundRest3;
Double_t backgroundRestNum3;
Double_t backgroundElPK1;
Double_t backgroundElPKNum1;
Double_t backgroundElPK2;
Double_t backgroundElPKNum2;
Double_t backgroundElPK3;
Double_t backgroundElPKNum3;
Double_t backgroundPiPK1;
Double_t backgroundPiPKNum1;
Double_t backgroundPiPK2;
Double_t backgroundPiPKNum2;
Double_t backgroundPiPK3;
Double_t backgroundPiPKNum3;


//******************** Definition of functions ****************************
void SaveHistos(TString cutString, Bool_t moreBG);

void HistoPlotSettings(TH1* histo1,
                       TString XTitle,
                       TString YTitle,
                       Double_t YMin,
                       Double_t YMax,
                       Double_t XMin,
                       Double_t XMax,
                       Color_t lineColor,
                       Color_t fillColor,
                       Style_t fillStyle);

void HistoMarkerPlotSettings(TH1* histo1,
                             TString XTitle,
                             TString YTitle,
                             Double_t YMin,
                             Double_t YMax,
                             Double_t XMin,
                             Double_t XMax,
                             Color_t lineColor,
                             Color_t markerColor,
                             Style_t markerStyle);

void SetStylePave(TPaveText *pave, Style_t textFontStyle,
                  Style_t textAlign,
                  Size_t textSize);

void Initialize(TString setPi0,
                TString Energy,
                Int_t   numberOfBins,
                Int_t   mode,
                Bool_t  addSig);
