// provided by Gamma Conversion Group, $ALICE_ROOT/PWG4/GammaConv ;https://twiki.cern.ch/twiki/bin/view/ALICE/PWG4GammaConversion

// Double, Int, etc.
TDatime now;
fstream fFileErrLog;
fstream fFileDataLog;

Int_t fCrysFitting = 0;
Int_t fIsMC = 0;
Int_t fMesonId = 0;
Float_t fBackgroundMultNumber;

Double_t fYMaxMeson = 0.9;
Double_t fMesonMassExpect = 0;   // expected meson mass
Double_t fNEvents = 0;

Double_t *fPeakRange;
Double_t *fFitRange;
Double_t *fBGFitRange    = nullptr;
Double_t *fMesonIntRange = nullptr;

// TSring 
TString fTextCent;
TString fEnergyFlag;
TString fPrefix;
TString fPrefix2;
TString fPeriodFlag;
TString ftextDayth;
TString fdate;
TString fCollisionSystem;
TString fTextMeasurement;
TString fCutSelection;
TString fBackgroundMultCutNumber;

TString fcMonth[12]={"Jan","Feb","Mar","Apr","May","Jun",
                    "Jul","Aug","Sep","Oct","Nov","Dec"};

TString fDirectPhoton;
TString fThesis ="";

// Output Files
TFile* fOutput                                          = nullptr;
TFile* fOutput1                                         = nullptr;
TFile* fOutput2                                         = nullptr;

//TH1D
TH1D*  fBckNorm                                         = nullptr;
TH1D*  fSignal                                          = nullptr;
TH1D*  fRatioSB                                         = nullptr;
TH1D*  fFittingHistMidPtSignal                          = nullptr;
TH1D*  fFittingHistMidPtBackground                      = nullptr;
TH1D*  fFittingHistMidPtSignalSub                       = nullptr;
TH1D*  fMesonFullPtSignal                               = nullptr;
TH1D*  fMesonFullPtBackground                           = nullptr;
TH1D*  fMesonFullPtBackNorm                             = nullptr;
TH1D*  fHistoMotherZMProj;
TH1D*  fHistoBckZMProj;
TH1D*  fNumberOfGoodESDTracks;
TH1D*  fEventQuality;
TH1D*  fHistoMappingBackNormInvMass                     = nullptr;
TH1D*  fHistoMappingSignalInvMass                       = nullptr;
TH1D*  fHistoYieldMeson                                 = nullptr;
TH1D*  fHistoYieldMesonBackFit                          = nullptr;
TH1D*  fHistoYieldMesonPerEvent                         = nullptr;
TH1D*  fHistoYieldMesonPerEventBackFit                  = nullptr;
TH1D*  fHistoSignMeson                                  = nullptr;
TH1D*  fHistoSBMeson                                    = nullptr;
TH1D*  fHistoTrueSignMeson                              = nullptr;
TH1D*  fHistoTrueSBMeson                                = nullptr;
TH1D*  fHistoYieldMesonNarrow                           = nullptr;
TH1D*  fHistoYieldMesonPerEventNarrow                   = nullptr;
TH1D*  fHistoSignMesonNarrow                            = nullptr;
TH1D*  fHistoSBMesonNarrow                              = nullptr;
TH1D*  fHistoYieldMesonWide                             = nullptr;
TH1D*  fHistoYieldMesonPerEventWide                     = nullptr;
TH1D*  fHistoSignMesonWide                              = nullptr;
TH1D*  fHistoSBMesonWide                                = nullptr;
TH1D*  fHistoMassPosition                               = nullptr;
TH1D*  fHistoMassMeson                                  = nullptr;
TH1D*  fHistoWidthMeson                                 = nullptr;
TH1D*  fHistoFWHMMeson                                  = nullptr;
TH1D*  fHistoTrueMassMeson                              = nullptr;
TH1D*  fHistoTrueMassMesonReweighted                    = nullptr;
TH1D*  fHistoTrueFWHMMeson                              = nullptr;
TH1D*  fHistoTrueFWHMMesonReweighted                    = nullptr;
TH1D*  fHistoFWHMMesonAlpha01                           = nullptr;
TH1D*  fDeltaPt                                         = nullptr;
TH1D*  fHistoYieldMesonLeft                             = nullptr;
TH1D*  fHistoYieldMesonLeftPerEvent                     = nullptr;
TH1D*  fHistoSignMesonLeft                              = nullptr;
TH1D*  fHistoSBMesonLeft                                = nullptr;
TH1D*  fHistoYieldMesonLeftNarrow                       = nullptr;
TH1D*  fHistoYieldMesonLeftPerEventNarrow               = nullptr;
TH1D*  fHistoSignMesonLeftNarrow                        = nullptr;
TH1D*  fHistoSBMesonLeftNarrow                          = nullptr;
TH1D*  fHistoYieldMesonLeftWide                         = nullptr;
TH1D*  fHistoYieldMesonLeftPerEventWide                 = nullptr;
TH1D*  fHistoSignMesonLeftWide                          = nullptr;
TH1D*  fHistoSBMesonLeftWide;
TH1D*  fHistoMassMesonLeft                              = nullptr;
TH1D*  fHistoWidthMesonLeft                             = nullptr;
TH1D*  fHistoFWHMMesonLeft                              = nullptr;
TH1D*  fHistoMCMesonPtWithinAcceptance                  = nullptr;
TH1D*  fHistoMCMesonPt                                  = nullptr;
TH1D*  fHistoMCMesonPtWOWeights                         = nullptr;
TH1D*  fHistoMCMesonPtWeights                           = nullptr;
TH1D*  fHistoMCMesonWithinAccepPt                       = nullptr; // Proper bins in Pt
TH1D*  fHistoMCMesonPt1                                 = nullptr; // Proper bins in Pt
TH1D*  fHistoYieldTrueMeson                             = nullptr;
TH1D*  fHistoYieldTrueMesonWide                         = nullptr;
TH1D*  fHistoYieldTrueMesonNarrow                       = nullptr;
TH1D*  fHistoYieldTrueMesonReweighted                   = nullptr;
TH1D*  fHistoYieldTrueMesonReweightedWide               = nullptr;
TH1D*  fHistoYieldTrueMesonReweightedNarrow             = nullptr;

TH1D*  fHistoYieldTrueGGMeson                           = nullptr;
TH1D*  fHistoYieldTrueGGCont                            = nullptr;
TH1D*  fHistoYieldTrueDalitzMeson                       = nullptr;
TH1D*  fHistoYieldTrueDalitzCont                        = nullptr;
TH1D*  fHistoYieldTrueGGMesonWide                       = nullptr;
TH1D*  fHistoYieldTrueGGContWide                        = nullptr;
TH1D*  fHistoYieldTrueDalitzMesonWide                   = nullptr;
TH1D*  fHistoYieldTrueDalitzContWide                    = nullptr;
TH1D*  fHistoYieldTrueGGMesonNarrow                     = nullptr;
TH1D*  fHistoYieldTrueGGContNarrow                      = nullptr;
TH1D*  fHistoYieldTrueDalitzMesonNarrow                 = nullptr;
TH1D*  fHistoYieldTrueDalitzContNarrow                  = nullptr;
TH1D*  fHistoMCMesonAcceptPt                            = nullptr;
TH1D*  fHistoMCMesonEffiPt                              = nullptr;
TH1D*  fHistoMCMesonEffiFitPt                           = nullptr;
TH1D*  fHistoTrueMesonEffiPt                            = nullptr;
TH1D*  fHistoMonteMesonEffiPt                           = nullptr;
TH1D*  fHistoMonteMesonEffiBackFitPt                    = nullptr;
TH1D*  fHistoMonteMesonNarrowEffiPt                     = nullptr;
TH1D*  fHistoMonteMesonWideEffiPt                       = nullptr;
TH1D*  fHistoMonteMesonLeftEffiPt                       = nullptr;
TH1D*  fHistoMonteMesonLeftNarrowEffiPt                 = nullptr;
TH1D*  fHistoMonteMesonLeftWideEffiPt                   = nullptr;
TH1D*  fHistoMonteMesonEffiMCAll                        = nullptr;
TH1D * fHistoMCTrueMesonEffiPt                          = nullptr;
TH1D * fHistoMCTrueMesonEffiPtReweighted                = nullptr;
TH1D * fHistoMCTrueMesonNarrowEffiPt                    = nullptr;
TH1D * fHistoMCTrueMesonNarrowEffiPtReweighted          = nullptr;
TH1D * fHistoMCTrueMesonWideEffiPt                      = nullptr;
TH1D * fHistoMCTrueMesonWideEffiPtReweighted            = nullptr;
TH1D*  fHistoTrueMesonEffiFitPt                         = nullptr;
TH1D*  fHistoMonteMesonEffiFitPt                        = nullptr;
TH1D*  fHistoMonteMesonNarrowEffiFitPt                  = nullptr;
TH1D*  fHistoMonteMesonWideEffiFitPt                    = nullptr;
TH1D*  fHistoMonteMesonLeftEffiFitPt                    = nullptr;
TH1D*  fHistoMonteMesonLeftNarrowEffiFitPt              = nullptr;
TH1D*  fHistoMonteMesonLeftWideEffiFitPt                = nullptr;
TH1D*  fHistoMonteMesonEffiFitMCAll                     = nullptr;
TH1D * fHistoMCTrueMesonEffiFitPt                       = nullptr;
TH1D * fHistoMCTrueMesonNarrowEffiFitPt                 = nullptr;
TH1D * fHistoMCTrueMesonWideEffiFitPt                   = nullptr;

TH1D** fHistoMappingTrueMesonInvMassPtBins              = nullptr;    
TH1D** fHistoMappingTrueMesonInvMassPtReweightedBins    = nullptr;    
TH1D** fHistoMappingTrueGGMesonInvMassPtBins            = nullptr;    
TH1D** fHistoMappingTrueDalitzInvMassPtBins             = nullptr;     
TH1D** fHistoMappingPiPlPiMiGammaInvMassPtBin           = nullptr;    
TH1D** fHistoMappingGGInvMassBackFitPtBin               = nullptr;    
TH1D** fHistoMappingGGInvMassBackFitWithoutSignalPtBin  = nullptr;    
TH1D** fHistoMappingBackInvMassPtBin                    = nullptr;
TH1D** fHistoMappingBackNormInvMassPtBin                = nullptr;
TH1D** fHistoMappingRatioSBInvMassPtBin                 = nullptr;
TH1D** fHistoMappingSignalInvMassPtBin                  = nullptr;
TH1D** fHistoMappingPeakPosInvMassPtBin                 = nullptr;
// Histograms for normalization on the left of the peak
TH1D**  fHistoMappingBackNormInvMassLeftPtBin           = nullptr;
TH1D**  fHistoMappingSignalInvMassLeftPtBin             = nullptr;


TF1**  fBackgroundFitPol                                = nullptr;
TF1 *  fFitReco                                         = nullptr;
TF1 *  fFitGausExp                                      = nullptr;
TF1 *  fFitLinearBck                                    = nullptr;
TF1 *  fFitLinearBckOut                                 = nullptr;
TF1 *  fFitSignalInvMassMidPt                           = nullptr;

Double_t  fYields;
Double_t  fYieldsError;

Double_t  fFWHMFunc;
Double_t  fFWHMFuncError;
Double_t  fYieldsFunc;
Double_t  fYieldsFuncError;
Double_t  fIntLinearBck;
Double_t  fIntLinearBckError;
Double_t  fIntLinearBckOut;
Double_t  fIntLinearBckErrorOut;

Float_t  pictDrawingCoordinatesFWHM[9] =   {0.6, 0.8, 0.30, 0.04, 0.15,0.7, 0.1, 0.035,0};
Int_t  fNRebinGlobal =      2;

// common meson analysis variables (So it is easy to switch between the pi0 and the eta without changing the code)
Double_t  *fBGFitRangeLeft                              = nullptr;
Double_t  *fMesonPlotRange                              = nullptr;

Double_t  *fMesonIntRangeWide                           = nullptr;
Double_t  *fMesonIntRangeNarrow                         = nullptr;
Double_t  *fMesonCurIntRange                            = nullptr;
Double_t  *fMesonCurIntRangeBackFit                     = nullptr;
Double_t  *fMesonCurIntRangeWide                        = nullptr;
Double_t  *fMesonCurIntRangeNarrow                      = nullptr;
Double_t  *fMesonCurLeftIntRange                        = nullptr;
Double_t  *fMesonCurLeftIntRangeWide                    = nullptr;
Double_t  *fMesonCurLeftIntRangeNarrow                  = nullptr;
Double_t  *fMesonTrueIntRange                           = nullptr;
Double_t  *fMesonTrueIntRangeWide                       = nullptr;
Double_t  *fMesonTrueIntRangeNarrow                     = nullptr;
Double_t    *fMesonTrueIntReweightedRange               = nullptr;
Double_t    *fMesonTrueIntReweightedRangeWide           = nullptr;
Double_t    *fMesonTrueIntReweightedRangeNarrow         = nullptr;

Double_t  *fMesonMassRange                              = nullptr;
Double_t  *fMesonFitRange                               = nullptr;
Double_t  *fMesonFitRangeWithoutPeak                    = nullptr;

Double_t  fMesonWidthExpect=     0;
Double_t  fMesonLambdaTail=     0;
Double_t  *fMesonWidthRange                             = nullptr;
Double_t  *fMesonLambdaTailRange                        = nullptr;
Double_t  *fMidPt                                       = nullptr;
Double_t  *fFullPt                                      = nullptr;
// end common meson analysis variables

//Background histograms in different M and Z bins

TH2D* fHistoMotherZM;
TH2D* fHistoBckZM;

Double_t  fScalingFactorBck[7][4];
Double_t  fCBAlpha;
Double_t  fCBn;

///////////////////////////
TString  fNameHistoGG;
TString  fNameHistoBack;
TString  fNameHistoPP;
TString  fNameHistoBackNormOut;
TString  fNameHistoBackNormLeftOut;
TString  fNameHistoTrue;
TString  fNameHistoTrueGG;
TString  fNameHistoTrueDalitz;
TString  fNameHistoEffi;
TString  fNameHistoFrac;
TString  fNameHistoMotherZM;
TString  fNameHistoBckZM;
TString  fNameFitSignalPos;

Double_t* fPiPlPiMiGammaYields                          = nullptr;
Double_t* fBckYields                                    = nullptr;
Double_t* fTotalBckYields                               = nullptr;
Double_t* fMesonYields                                  = nullptr;
Double_t* fMesonYieldsBackFit                           = nullptr;
Double_t* fMesonTrueYields                              = nullptr;
Double_t* fMesonTrueYieldsReweighted                    = nullptr;
Double_t* fMesonTrueGGContYields                        = nullptr;
Double_t* fMesonTrueDalitzContYields                    = nullptr;
Double_t* fMesonTrueGGContYieldsWide                    = nullptr;
Double_t* fMesonTrueDalitzContYieldsWide                = nullptr;
Double_t* fMesonTrueGGContYieldsNarrow                  = nullptr;
Double_t* fMesonTrueDalitzContYieldsNarrow              = nullptr;
Double_t* fMesonYieldsFunc                              = nullptr;
Double_t* fMesonYieldsResidualBckFunc                   = nullptr;
Double_t* fMesonYieldsResidualBckFuncBackFit            = nullptr;
Double_t* fMesonYieldsCorResidualBckFunc                = nullptr;
Double_t* fMesonYieldsCorResidualBckFuncBackFit         = nullptr;
Double_t* fMesonYieldsPerEvent                          = nullptr;
Double_t* fMesonYieldsPerEventBackFit                   = nullptr;
Double_t* fMesonMass                                    = nullptr;
Double_t* fMesonMassBackFit                             = nullptr;
Double_t* fMesonWidth                                   = nullptr;
Double_t* fMesonWidthBackFit                            = nullptr;
Double_t* fMesonSB                                      = nullptr;
Double_t* fMesonTrueSB                                  = nullptr;
Double_t* fMesonTrueSign                                = nullptr;
Double_t* fMesonSign                                    = nullptr;
Double_t* fMesonFWHM                                    = nullptr;
Double_t* fMesonTrueMass                                = nullptr;
Double_t* fMesonTrueMassReweighted                      = nullptr;
Double_t* fMesonTrueFWHM                                = nullptr;
Double_t* fMesonTrueFWHMReweighted                      = nullptr;
Double_t* fMesonFWHMAlpha01                             = nullptr;

// Normalization at the left of the peak
Double_t* fPiPlPiMiGammaYieldsLeft                      = nullptr;
Double_t* fBckYieldsLeft                                = nullptr;
Double_t* fTotalBckYieldsLeft                           = nullptr;
Double_t* fMesonYieldsLeft                              = nullptr;
Double_t* fMesonYieldsFuncLeft                          = nullptr;
Double_t* fMesonYieldsResidualBckFuncLeft               = nullptr;
Double_t* fMesonYieldsCorResidualBckFuncLeft            = nullptr;
Double_t* fMesonYieldsLeftPerEvent                      = nullptr;
Double_t* fMesonMassLeft                                = nullptr;
Double_t* fMesonWidthLeft                               = nullptr;
Double_t* fMesonSBLeft                                  = nullptr;
Double_t* fMesonSignLeft                                = nullptr;
Double_t* fMesonFWHMLeft                                = nullptr;

// Narrow Integration Window
Double_t* fPiPlPiMiGammaYieldsNarrow                    = nullptr;
Double_t* fBckYieldsNarrow                              = nullptr;
Double_t* fTotalBckYieldsNarrow                         = nullptr;
Double_t* fMesonYieldsNarrow                            = nullptr;
Double_t* fMesonTrueYieldsNarrow                        = nullptr;
Double_t* fMesonTrueYieldsReweightedNarrow              = nullptr;
Double_t* fMesonYieldsFuncNarrow                        = nullptr;
Double_t* fMesonYieldsResidualBckFuncNarrow             = nullptr;
Double_t* fMesonYieldsCorResidualBckFuncNarrow          = nullptr;
Double_t* fMesonYieldsPerEventNarrow                    = nullptr;
Double_t* fMesonSBNarrow                                = nullptr;
Double_t* fMesonSignNarrow                              = nullptr;

Double_t* fPiPlPiMiGammaYieldsLeftNarrow                = nullptr;
Double_t* fBckYieldsLeftNarrow                          = nullptr;
Double_t* fTotalBckYieldsLeftNarrow                     = nullptr;
Double_t* fMesonYieldsLeftNarrow                        = nullptr;
Double_t* fMesonYieldsFuncLeftNarrow                    = nullptr;
Double_t* fMesonYieldsResidualBckFuncLeftNarrow         = nullptr;
Double_t* fMesonYieldsCorResidualBckFuncLeftNarrow      = nullptr;
Double_t* fMesonYieldsLeftPerEventNarrow                = nullptr;
Double_t* fMesonSBLeftNarrow                            = nullptr;
Double_t* fMesonSignLeftNarrow                          = nullptr;

// Wide Integration Window
Double_t* fPiPlPiMiGammaYieldsWide                      = nullptr;
Double_t* fBckYieldsWide                                = nullptr;
Double_t* fTotalBckYieldsWide                           = nullptr;
Double_t* fMesonYieldsWide                              = nullptr;
Double_t* fMesonTrueYieldsWide                          = nullptr;
Double_t* fMesonTrueYieldsReweightedWide                = nullptr;
Double_t* fMesonYieldsFuncWide                          = nullptr;
Double_t* fMesonYieldsResidualBckFuncWide               = nullptr;
Double_t* fMesonYieldsCorResidualBckFuncWide            = nullptr;
Double_t* fMesonYieldsPerEventWide                      = nullptr;
Double_t* fMesonSBWide                                  = nullptr;
Double_t* fMesonSignWide                                = nullptr;

Double_t* fPiPlPiMiGammaYieldsLeftWide                  = nullptr;
Double_t* fBckYieldsLeftWide                            = nullptr;
Double_t* fTotalBckYieldsLeftWide                       = nullptr;
Double_t* fMesonYieldsLeftWide                          = nullptr;
Double_t* fMesonYieldsFuncLeftWide                      = nullptr;
Double_t* fMesonYieldsResidualBckFuncLeftWide           = nullptr;
Double_t* fMesonYieldsCorResidualBckFuncLeftWide        = nullptr;
Double_t* fMesonYieldsLeftPerEventWide                  = nullptr;
Double_t* fMesonSBLeftWide                              = nullptr;
Double_t* fMesonSignLeftWide                            = nullptr;

Double_t* fPiPlPiMiGammaYieldsError                     = nullptr;
Double_t* fBckYieldsError                               = nullptr;
Double_t* fTotalBckYieldsError                          = nullptr;
Double_t* fMesonYieldsError                             = nullptr;
Double_t* fMesonYieldsBackFitError                      = nullptr;
Double_t* fMesonTrueYieldsError                         = nullptr;
Double_t* fMesonTrueYieldsReweightedError               = nullptr;
Double_t* fMesonTrueGGContYieldsError                   = nullptr;
Double_t* fMesonTrueDalitzContYieldsError               = nullptr;
Double_t* fMesonTrueGGContYieldsWideError               = nullptr;
Double_t* fMesonTrueDalitzContYieldsWideError           = nullptr;
Double_t* fMesonTrueGGContYieldsNarrowError             = nullptr;
Double_t* fMesonTrueDalitzContYieldsNarrowError         = nullptr;
Double_t* fMesonYieldsFuncError                         = nullptr;
Double_t* fMesonYieldsResidualBckFuncError              = nullptr;
Double_t* fMesonYieldsResidualBckFuncBackFitError       = nullptr;
Double_t* fMesonYieldsCorResidualBckFuncError           = nullptr;
Double_t* fMesonYieldsCorResidualBckFuncBackFitError    = nullptr;
Double_t* fMesonYieldsPerEventError                     = nullptr;
Double_t* fMesonYieldsPerEventBackFitError              = nullptr;
Double_t* fMesonMassError                               = nullptr;
Double_t* fMesonMassBackFitError                        = nullptr;
Double_t* fMesonWidthError                              = nullptr;
Double_t* fMesonWidthBackFitError                       = nullptr;
Double_t* fMesonTrueMassError                           = nullptr;
Double_t* fMesonTrueMassReweightedError                 = nullptr;
Double_t* fMesonTrueFWHMReweightedError                 = nullptr;
Double_t* fMesonTrueFWHMError                           = nullptr;

Double_t* fMesonSBError                                 = nullptr;
Double_t* fMesonSignError                               = nullptr;
Double_t* fMesonTrueSBError                             = nullptr;
Double_t* fMesonTrueSignError                           = nullptr;
Double_t* fMesonFWHMError                               = nullptr;

Double_t* fPiPlPiMiGammaYieldsLeftError                 = nullptr;
Double_t* fBckYieldsLeftError                           = nullptr;
Double_t* fTotalBckYieldsLeftError                      = nullptr;
Double_t* fMesonYieldsLeftError                         = nullptr;
Double_t* fMesonYieldsFuncLeftError                     = nullptr;
Double_t* fMesonYieldsResidualBckFuncLeftError          = nullptr;
Double_t* fMesonYieldsCorResidualBckFuncLeftError       = nullptr;
Double_t* fMesonYieldsLeftPerEventError                 = nullptr;
Double_t* fMesonMassLeftError                           = nullptr;
Double_t* fMesonWidthLeftError                          = nullptr;
Double_t* fMesonSBLeftError                             = nullptr;
Double_t* fMesonSignLeftError                           = nullptr;
Double_t* fMesonFWHMLeftError                           = nullptr;
Double_t* fMesonFWHMAlpha01Error                        = nullptr;

// Narrow integration Window
Double_t* fPiPlPiMiGammaYieldsNarrowError               = nullptr;
Double_t* fBckYieldsNarrowError                         = nullptr;
Double_t* fTotalBckYieldsNarrowError                    = nullptr;
Double_t* fMesonYieldsNarrowError                       = nullptr;
Double_t* fMesonTrueYieldsNarrowError                   = nullptr;
Double_t* fMesonTrueYieldsReweightedNarrowError         = nullptr;
Double_t* fMesonYieldsFuncNarrowError                   = nullptr;
Double_t* fMesonYieldsResidualBckFuncNarrowError        = nullptr;
Double_t* fMesonYieldsCorResidualBckFuncNarrowError     = nullptr;
Double_t* fMesonYieldsPerEventNarrowError               = nullptr;
Double_t* fMesonSBNarrowError                           = nullptr;
Double_t* fMesonSignNarrowError                         = nullptr;

Double_t* fPiPlPiMiGammaYieldsLeftNarrowError           = nullptr;
Double_t* fBckYieldsLeftNarrowError                     = nullptr;
Double_t* fTotalBckYieldsLeftNarrowError                = nullptr;
Double_t* fMesonYieldsLeftNarrowError                   = nullptr;
Double_t* fMesonYieldsFuncLeftNarrowError               = nullptr;
Double_t* fMesonYieldsResidualBckFuncLeftNarrowError    = nullptr;
Double_t* fMesonYieldsCorResidualBckFuncLeftNarrowError = nullptr;
Double_t* fMesonYieldsLeftPerEventNarrowError           = nullptr;
Double_t* fMesonSBLeftNarrowError                       = nullptr;
Double_t* fMesonSignLeftNarrowError                     = nullptr;

// Wide integration Window
Double_t* fPiPlPiMiGammaYieldsWideError                 = nullptr;
Double_t* fBckYieldsWideError                           = nullptr;
Double_t* fTotalBckYieldsWideError                      = nullptr;
Double_t* fMesonYieldsWideError                         = nullptr;
Double_t* fMesonTrueYieldsWideError                     = nullptr;
Double_t* fMesonTrueYieldsReweightedWideError           = nullptr;
Double_t* fMesonYieldsFuncWideError                     = nullptr;
Double_t* fMesonYieldsResidualBckFuncWideError          = nullptr;
Double_t* fMesonYieldsCorResidualBckFuncWideError       = nullptr;
Double_t* fMesonYieldsPerEventWideError                 = nullptr;
Double_t* fMesonSBWideError                             = nullptr;
Double_t* fMesonSignWideError                           = nullptr;

Double_t* fPiPlPiMiGammaYieldsLeftWideError             = nullptr;
Double_t* fBckYieldsLeftWideError                       = nullptr;
Double_t* fTotalBckYieldsLeftWideError                  = nullptr;
Double_t* fMesonYieldsLeftWideError                     = nullptr;
Double_t* fMesonYieldsFuncLeftWideError                 = nullptr;
Double_t* fMesonYieldsResidualBckFuncLeftWideError      = nullptr;
Double_t* fMesonYieldsCorResidualBckFuncLeftWideError   = nullptr;
Double_t* fMesonYieldsLeftPerEventWideError             = nullptr;
Double_t* fMesonSBLeftWideError                         = nullptr;
Double_t* fMesonSignLeftWideError                       = nullptr;

TF1**  fFitSignalInvMassPtBin                           = nullptr;
TF1**  fFitSignalInvMassBackFitPtBin                    = nullptr;
TF1**  fFitTrueSignalInvMassPtBin                       = nullptr;
TF1**    fFitTrueSignalInvMassPtReweightedBin           = nullptr;
TF1**  fFitSignalPeakPosInvMassPtBin                    = nullptr;
TF1**  fFitSignalPeakPosInvMassBackFitPtBin             = nullptr;
TF1**  fFitBckInvMassPtBin                              = nullptr;
TF1**  fFitBckInvMassBackFitPtBin                       = nullptr;
TF1**  fFitRatioInvMassPtBin                            = nullptr;
TF1** fFitWithPol2ForBG                                 = nullptr;
TF1*     fFitDCAZPhotonTrainBGH                         = nullptr;
TF1*     fFitDCAZPhotonTrainBG                          = nullptr;
TH1F*    fHistDCAZPhotonTrainBG                         = nullptr;
TF1*     fFitDCAZPhotonInterTrainBG                     = nullptr;
TF1*     fFitDCAZPhotonInterTrainBGH                    = nullptr;
TF1*     fFitDCAZPhotonBGH                              = nullptr;
TF1*     fFitDCAZPhotonBG                               = nullptr;
TF1**  fFitPeakPosPtBin                                 = nullptr;
Double_t* fMesonMassPeakPos                             = nullptr;
Double_t* fMesonMassPeakPosError                        = nullptr;

TF1**  fFitInvMassLeftPtBin                             = nullptr;
TF1**  fFitSignalPeakPosInvMassLeftPtBin                = nullptr;
TF1**  fFitBckInvMassLeftPtBin                          = nullptr;

TH2D*  fHistoTrueMesonInvMassVSPt                       = nullptr;
TH2D*    fHistoTrueMesonInvMassVSPtWOWeights            = nullptr;
TH2D*    fHistoTrueMesonInvMassVSPtReweighted           = nullptr;
TProfile2D*    fProfileTrueMesonInvMassVSPtWeights      = nullptr;
TH2D*  fHistoTrueMesonInvMassVSPtSec                    = nullptr;
TH2D*  fHistoTrueGGMesonInvMassVSPt                     = nullptr;
TH2D*  fHistoTrueDalitzMesonInvMassVSPt                 = nullptr;
TH2D*   fGammaGammaInvMassVSPt                          = nullptr;
TH2D*   fBckInvMassVSPt                                 = nullptr;

TH2F**    fHistoWeightsBGZbinVsMbin = 0x0;
TH2F**    fHistoFillPerEventBGZbinVsMbin = 0x0;


TString ObjectNameTrue;
TString ObjectNameTrueWOWeights;
TString ObjectNameProfileWeights;
TString ObjectNameContaminationGG;
TString ObjectNameContaminationDalitz;
TString ObjectNameMCEtaAcc;
TString ObjectNameMCEta;
TString ObjectNameMCEtaWOWeights;
   
Bool_t fAdvancedMesonQA = kFALSE;
Bool_t fEstimateTrainPileUp = kFALSE;

void ProcessEM(TH1D*,TH1D*,Double_t *);
void ProcessBckFitSubtraction(TH1D *fGammaGamma, Int_t i, Double_t * fPeakRange, Double_t *fFitRange);
void ProcessRatioSignalBackground(TH1D* , TH1D* );
void FillMassHistosArray(TH2D*);
void FillMassMCTrueMesonHistosArray(TH2D*);
void FillMassMCTrueGGMesonHistosArray(TH2D*);
void FillMassMCTrueDalitzMesonHistosArray(TH2D*);
TH1D* CalulateContaminationFraction(TH1D* histoRawYield, TH1D* histoRawYieldSec, TString nameHistoFrac);
void CreatePtHistos();
void FillPtHistos();
void FitSubtractedInvMassInPtBins(TH1D * ,Double_t *, Int_t, Bool_t );  // Fits the Invariant Mass histos with a given function
void FitTrueInvMassInPtBins(TH1D * ,Double_t *, Int_t);  // Fits the Invariant Mass histos with a given function
void FitWithPol2ForBG(TH1D*,Double_t*,Int_t,Bool_t);
void ProduceBckProperWeighting(TList*, TList* );
void ProduceBckWithoutWeighting(TH2D *);
void IntegrateHistoInvMassStream(TH1D * , Double_t *);
void IntegrateHistoInvMass(TH1D * , Double_t *);
void IntegrateFitFunc(TF1 * , TH1D *, Double_t *);
void FillHistosArrayMC(TH1D* , TH1D*, TH1D* ,TString);
void CalculateMesonAcceptance();
void CalculateMesonEfficiency(TH1D*, TH1D*,TString,TString);
void SaveHistos(Int_t, TString, TString);
void SaveCorrectionHistos(TString , TString);
void Initialize(TString setPi0, Int_t);
void CalculateFWHM(TF1 *);
void Delete();
void SetCorrectMCHistogrammNames(TString, TString);

TString centralityString = "";