// provided by Gamma Conversion Group, $ALICE_PHYSICS/PWGGA/GammaConv ;https://twiki.cern.ch/twiki/bin/view/ALICE/PWG4GammaConversion

//************************** Some general definitions *************************************
TString     fAnalyzedMeson                              = "";
TString     fEnergyFlag                                 = "";
TString     fPeriodFlag                                 = "";
TString     fdirectphoton                               = "";
TString     fSuffix                                     = "";
TString     fCutSelection                               = "";
TString     fEventCutSelection                          = "";
TString     fGammaCutSelection                          = "";
TString     fClusterCutSelection                        = "";
TString     fElectronCutSelection                       = "";
TString     fMesonCutSelection                          = "";
Double_t    ptMin                                       = 0;
Double_t    ptMax                                       = 20;
Double_t    ptPlotMin                                   = 0;
Double_t    ptPlotMax                                   = 20;
Double_t    fRapidity                                   = 0;
Int_t       fMode                                       = 0;
Float_t     nEvents                                     = 0;

const Int_t nMotherParticles                            = 16;
TString     motherParticles[nMotherParticles]           = {"Eta","K0s","K0l","Lambda","rho0","EtaPrim","omega","rho+","rho-","phi","J/psi","Delta-","Delta0","Delta+","Delta++","Sigma0"};
TString     motherParticlesPDG[nMotherParticles]        = {"221","310","130","3122","113","331","223","213","-213","333","443","1114","2114","2214","2224","3212"};
TString     motherParticlesLatex[nMotherParticles]      = { "#eta","K^{0}_{S}","K^{0}_{L}","#Lambda","#rho^{0}","#eta'","#omega","#rho^{+}","#rho^{-}","#phi",
                                                            "J/#psi","#Delta^{-}","#Delta^{0}","#Delta^{+}","#Delta^{++}","#Sigma^{0}"};
//ctau given in (cm) below! - zero means: ctau < 0.01cm
Double_t    motherParticles_ctau[nMotherParticles]      = {0,2.68,1533.74,7.89,0,0,0,0,0,0,0,0,0,0,0,0};
Double_t    motherFactorDecayLength[nMotherParticles]   = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
Int_t       motherParticleDec[nMotherParticles]         = {2,256,65536,131072,4,16,8,8192,16384,32,64,2048,4096,1024,512,128};
Color_t     cocktailColor[nMotherParticles]             = {kRed+2,kBlue+1,kOrange+1,kYellow+2,kAzure-2,kGreen+2,kRed-2,kViolet,kMagenta,kViolet+2,kBlue-3,kTeal+9,kCyan+2,kMagenta+2,kCyan+4,kViolet+4};
Style_t     cocktailMarker[nMotherParticles]            = {20,21,24,25,20,21,24,25,20,21,24,25,20,21,24,25};
TString     decayChannelsLatex[nMotherParticles][18];
Double_t    decayChannelsBR[nMotherParticles][18];

// additional scaling factor (param per inel. event vs. MB event)
Double_t    eventNormScalingFactor                      = 1.;

//************************** cocktail settings ********************************************
Double_t ptGenMin                                       = 0;
Double_t ptGenMax                                       = 20;
Int_t    nParticles                                     = 1000;
Int_t    selectedMothers                                = 62591;
Bool_t   hasMother[nMotherParticles]                    = {0};
Double_t mtScaleFactor[nMotherParticles]                = {1.};

//************************** Declaration of histograms ************************************
TH1F*  fDeltaPt                                         = NULL;
TH1F*  histoNEvents                                     = NULL;
TH1F*  histMtScalingFactors                             = NULL;
TH1F** histoRatioPi0FromXToPi0Param                     = NULL;
TH1F** histoDecayChannels                               = NULL;
TH1F** histoDecayChannelsBR                             = NULL;
TH2F** histoMesonDaughterPtY                            = NULL;
TH2F** histoMesonDaughterPtYCorr                        = NULL;
TH2F** histoMesonDaughterPtPhi                          = NULL;
TH2F** histoMesonMotherPtY                              = NULL;
TH2F** histoMesonMotherPtPhi                            = NULL;
TH1F** histoMesonDaughterPtOrBin                        = NULL;
TH1F** histoMesonDaughterYOrBin                         = NULL;
TH1F** histoMesonDaughterPhiOrBin                       = NULL;
TH1F** histoMesonMotherPtOrBin                          = NULL;
TH1F** histoMesonMotherYOrBin                           = NULL;
TH1F** histoMesonMotherPhiOrBin                         = NULL;
TH2F** histoGammaFromXFromMotherPtY                     = NULL;
TH2F** histoGammaFromXFromMotherPtYCorr                 = NULL;
TH2F** histoGammaFromXFromMotherPtPhi                   = NULL;
TH1F** histoGammaFromXFromMotherPtOrBin                 = NULL;
TH1F** histoGammaFromXFromMotherYOrBin                  = NULL;
TH1F** histoGammaFromXFromMotherPhiOrBin                = NULL;

//************************** Declaration of lists *****************************************
TList** listYSlicesMesonDaughter                        = NULL;
TList** listYSlicesGammaFromXFromMother                 = NULL;

//************************** Cocktail input ***********************************************
TFile* cocktailInputFile                                = NULL;
TList* cocktailInputList                                = NULL;
TF1**  cocktailInputParametrizations                    = NULL;
TF1**  cocktailInputParametrizationsMtScaled            = NULL;
TF1*   cocktailInputParametrizationPi0                  = NULL;
TF1*   paramScaleBase                                   = NULL;

//************************** Methods ******************************************************
void        Initialize                              (   TString meson,
                                                        TString energy,
                                                        Int_t   numberOfBins    );
void        RebinSpectrum                           (   TH1F*   Spectrum,
                                                        TString NewName         );
void        SaveMesonHistos                         (                           );
void        SavePhotonHistos                        (                           );
TF1*        MtScaledParam                           (   TF1*    param,
                                                        Int_t   particleNumber  );
Double_t    GetMass                                 (   TString particle        );
void        SetHistogramTitles                      (   TH1F*   input,
                                                        TString title,
                                                        TString xTitle,
                                                        TString yTitle          );
TH1D*       CalculateRatioToTF1                     (   TH1D*   hist,
                                                        TF1*    func            );
void        CorrectForNonFlatRapidity               (   TH2F*   histCorr,
                                                        TH2F*   histOr,
                                                        TList*  list            );
void        CreateBRTableLatex                      (                           );
Int_t       GetMinimumBinAboveThreshold             (   TH1F*       hist,
                                                        Double_t    thres       );

