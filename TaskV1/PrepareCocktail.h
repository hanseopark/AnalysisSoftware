// provided by Gamma Conversion Group, $ALICE_ROOT/PWG4/GammaConv ;https://twiki.cern.ch/twiki/bin/view/ALICE/PWG4GammaConversion

//************************** Some general definitions *************************************
//TString fDate                                         = "";
TString fEnergyFlag                                     = "";
TString fPeriodFlag                                     = "";
TString fSuffix                                         = "";
TString fCutSelection                                   = "";
TString fEventCutSelection                              = "";
TString fGammaCutSelection                              = "";
TString fClusterCutSelection                            = "";
TString fElectronCutSelection                           = "";
TString fMesonCutSelection                              = "";
TString centralityString                                = "";
TString fTextCent                                       = "";
TString fCollisionSystem                                = "";
Double_t ptMin                                          = 0;
Double_t ptMax                                          = 20;
Double_t fRapidity                                      = 0;
Int_t fMode                                             = 0;
Float_t nEvents                                         = 0;

const Int_t nMotherParticles                            = 14;
TString motherParticles[nMotherParticles]               = {"Pi0", "Eta", "EtaPrim", "omega", "rho0", "rho+", "rho-", "phi", "J/psi", "Delta-", "Delta0", "Delta+", "Delta++", "Sigma0"};
TString motherParticlesPDG[nMotherParticles]            = {"111", "221", "331", "223", "113", "213", "-213", "333", "443", "1114", "2114", "2214", "2224", "3212"};
TString motherParticlesLatex[nMotherParticles]          = {"#pi^{0}", "#eta", "#eta'", "#omega", "#rho^{0}", "#rho^{+}", "#rho^{-}", "#phi", "J/#psi", "#Delta^{-}", "#Delta^{0}", "#Delta^{+}", "#Delta^{++}",
                                                            "#Sigma^{0}"};
Int_t motherParticleDec[nMotherParticles]               = {1, 2, 16, 8, 4, 8192, 16384, 32, 64, 2048, 4096, 1024, 512, 128};
TString decayChannelsLatex[nMotherParticles]            = {""};
Color_t cocktailColor[nMotherParticles]                 = {kBlue+3, kBlue, kBlue-3, kRed-3, kRed, kRed+3, kAzure+3, kAzure, kAzure-3, kGreen+3, kGreen, kCyan+3, kCyan+1, kCyan-3};
Style_t cocktailMarker[nMotherParticles]                = {20, 21, 24, 25, 20, 21, 24, 25, 20, 21, 24, 25, 20, 21};

//************************** cocktail settings ********************************************
Double_t ptGenMin                                       = 0;
Double_t ptGenMax                                       = 20;
Int_t nParticles                                        = 1000;
Int_t selectedMothers                                   = 62591;
Bool_t hasMother[nMotherParticles]                      = {0};

//************************** mT scale factor (AliGenEMlib) ********************************
const Double_t mtScaleFactor[3][14] = {
    // {1.0, 0.5, 1.0, 0.9, 0.4, 0.23, 0.054},  // factor for pp from arXiv:1110.3929
    // {1.0, 0.55, 1.0, 0.9, 0.4, 0.25, 0.004}    // factor for PbPb from arXiv:1110.3929
    //{1., 0.48, 1.0, 0.9, 0.25, 0.4}, (old values)
    //{1., 0.48, 1.0, 0.9, 0.4, 0.25}, (nlo values)
    //{1., 0.48, 1.0, 0.8, 0.4, 0.2, 0.06} (combination of nlo and LHC measurements)
    //https://aliceinfo.cern.ch/Figure/node/2634
    //https://aliceinfo.cern.ch/Figure/node/2788
    //https://aliceinfo.cern.ch/Figure/node/4403
    //https://aliceinfo.cern.ch/Figure/node/5842
    //https://aliceinfo.cern.ch/Notes/node/87
    /*best guess:
     - pp values for eta/pi0 [arXiv:1205.5724], omega/pi0 [arXiv:1210.5749], phi/(pi+/-) [arXiv:1208.5717] from measured 7 Tev data
     */
    {1.0, 0.476, 0.4, 0.85, 1.0, 1.0, 1.0, 0.13, 1.0, 1.0, 1.0, 1.0, 1.0, 0.49}, //pp
    {1.0, 0.476, 0.4, 0.85, 1.0, 1.0, 1.0, 0.25, 1.0, 1.0, 1.0, 1.0, 1.0, 0.49}, //pPb
    {1.0, 0.476, 0.4, 0.85, 1.0, 1.0, 1.0, 0.25, 1.0, 1.0, 1.0, 1.0, 1.0, 0.49}, //PbPb
};


//************************** Declaration of histograms ************************************
TH1F* fDeltaPt                                          = NULL;
TH1F* histoNEvents                                      = NULL;
TH1F* histMtScalingFactors                              = NULL;
TH1D* histoPi0YieldData                                 = NULL;
TH1F** histoDecayChannels                               = NULL;
TH1F** histoPythiaBRs                                   = NULL;
TH2F* histoGammaSumPtY                                  = NULL;
TH2F* histoGammaSumPtPhi                                = NULL;
TH1F* histoGammaSumPtOrBin                              = NULL;
TH1F* histoGammaSumPt                                   = NULL;
TH1F* histoGammaSumYOrBin                               = NULL;
TH1F* histoGammaSumPhiOrBin                             = NULL;
TH2F** histoGammaPtY                                    = NULL;
TH2F** histoGammaPtPhi                                  = NULL;
TH1F** histoGammaPtOrBin                                = NULL;
TH1F** histoGammaPt                                     = NULL;
TH1F** histoGammaYOrBin                                 = NULL;
TH1F** histoGammaPhiOrBin                               = NULL;
TH2F** histoGammaMotherPtY                              = NULL;
TH1F** histoGammaMotherPtOrBin                          = NULL;
TH1F** histoGammaMotherPt                               = NULL;
TH1F** histoGammaMotherYOrBin                           = NULL;
TH2F** histoGammaMotherPtPhi                            = NULL;
TH1F** histoGammaMotherPhiOrBin                         = NULL;


//************************** Cocktail input ***********************************************
TFile* cocktailInputFile                                = NULL;
TList* cocktailInputList                                = NULL;
TF1** cocktailInputParametrizations                     = NULL;
TF1** cocktailInputParametrizationsMtScaled             = NULL;
TF1* paramScaleBase                                     = NULL;


//************************** Methods ******************************************************
void Initialize                                     (   TString energy,
                                                        Int_t numberOfBins      );
void RebinSpectrum                                  (   TH1F *Spectrum,
                                                        TString NewName         );
TH1F* ConvertYieldHisto                             (   TH1F* input             );
void SaveHistos                                     (                           );
TF1* MtScaledParam                                  (   TF1* param,
                                                        Int_t particleNumber    );
Double_t GetMass                                    (   TString particle        );
void SetHistogramTitles                             (   TH1F* input,
                                                        TString title,
                                                        TString xTitle,
                                                        TString yTitle          );
void DeleteObjects                                  (                           );
