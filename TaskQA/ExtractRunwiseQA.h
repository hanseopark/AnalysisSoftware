TFile *File;
TDirectoryFile *folderInFile;

TH1F *VertexZ;
Double_t Events;

TH1D *elNSigmaElectronCut= NULL;
TH1D *posNSigmaElectronCut= NULL;
Double_t meanelNSigmaElectronCut;
Double_t meanErrelNSigmaElectronCut;
Double_t meanposNSigmaElectronCut;
Double_t meanErrposNSigmaElectronCut;

TH1D *elNSigmaProtonCut= NULL;
TH1D *posNSigmaProtonCut= NULL;
Double_t meanelNSigmaProtonCut;
Double_t meanErrelNSigmaProtonCut;
Double_t meanposNSigmaProtonCut;
Double_t meanErrposNSigmaProtonCut;

TH3F *NSigmadEdxElectron= NULL;
TH3F *NSigmadEdxPositron= NULL;

//nsigma
TH2D *projxNSigmadEdxElectron= NULL;
TH2D *projxNSigmadEdxPositron= NULL;
Double_t meanNSigmadEdxElectron;
Double_t meanErrNSigmadEdxElectron;
Double_t meanNSigmadEdxPositron;
Double_t meanErrNSigmadEdxPositron;
//eta
TH1D *projyNSigmadEdxElectron= NULL;
TH1D *projyNSigmadEdxPositron= NULL;
Double_t meanEtaElectron;
Double_t meanErrEtaElectron;
Double_t meanEtaPositron;
Double_t meanErrEtaPositron;
Double_t meanEtaNegElectron;
Double_t meanErrEtaNegElectron;
Double_t meanEtaNegPositron;
Double_t meanErrEtaNegPositron;
Double_t meanEtaPosElectron;
Double_t meanErrEtaPosElectron;
Double_t meanEtaPosPositron;
Double_t meanErrEtaPosPositron;

//gamma
TH2F *GammaEtaPt= NULL;
TH1D *GammaEta= NULL;
Double_t meanGammaEta;
Double_t meanErrGammaEta;
TH1D *GammaPt= NULL;
Double_t meanGammaPt;
Double_t meanErrGammaPt;
Double_t meanGammaEtaNeg;
Double_t meanErrGammaEtaNeg;
Double_t meanGammaEtaPos;
Double_t meanErrGammaEtaPos;

Double_t sigmaelNSigmaElectronCut;
Double_t sigmaErrelNSigmaElectronCut;
Double_t sigmaposNSigmaElectronCut;
Double_t sigmaErrposNSigmaElectronCut;

Int_t bin1el;
Double_t bin1elErr;
Int_t bin2el;
Double_t bin2elErr;
Double_t widthelNSigmaElectronCut;
Double_t widthErrelNSigmaElectronCut;
Int_t bin1pos;
Double_t bin1posErr;
Int_t bin2pos;
Double_t bin2posErr;
Double_t widthposNSigmaElectronCut;
Double_t widthErrposNSigmaElectronCut;

Double_t sigmaelNSigmaProtonCut;
Double_t sigmaErrelNSigmaProtonCut;
Double_t sigmaposNSigmaProtonCut;
Double_t sigmaErrposNSigmaProtonCut;
Int_t bin3el;
Double_t bin3elErr;
Int_t bin4el;
Double_t bin4elErr;
Double_t widthelNSigmaProtonCut;
Double_t widthErrelNSigmaProtonCut;
Int_t bin3pos;
Double_t bin3posErr;
Int_t bin4pos;
Double_t bin4posErr;
Double_t widthposNSigmaProtonCut; 
Double_t widthErrposNSigmaProtonCut; 

Double_t ErrNumGammaAllSectorEtaNeg; 
Double_t ErrNumGammaAllSectorEtaPos; 
Double_t ErrNumGammaSector0EtaNeg; 
Double_t ErrNumGammaSector0EtaPos;
Double_t ErrNumGammaSector1EtaNeg; 
Double_t ErrNumGammaSector1EtaPos; 
Double_t ErrNumGammaSector2EtaNeg;
Double_t ErrNumGammaSector2EtaPos; 
Double_t ErrNumGammaSector3EtaNeg; 
Double_t ErrNumGammaSector3EtaPos; 
Double_t ErrNumGammaSector4EtaNeg; 
Double_t ErrNumGammaSector4EtaPos; 
Double_t ErrNumGammaSector5EtaNeg; 
Double_t ErrNumGammaSector5EtaPos; 
Double_t ErrNumGammaSector6EtaNeg; 
Double_t ErrNumGammaSector6EtaPos; 
Double_t ErrNumGammaSector7EtaNeg; 
Double_t ErrNumGammaSector7EtaPos; 
Double_t ErrNumGammaSector8EtaNeg; 
Double_t ErrNumGammaSector8EtaPos; 
Double_t ErrNumGammaSector9EtaNeg; 
Double_t ErrNumGammaSector9EtaPos; 
Double_t ErrNumGammaSector10EtaNeg; 
Double_t ErrNumGammaSector10EtaPos; 
Double_t ErrNumGammaSector11EtaNeg; 
Double_t ErrNumGammaSector11EtaPos; 
Double_t ErrNumGammaSector12EtaNeg; 
Double_t ErrNumGammaSector12EtaPos; 
Double_t ErrNumGammaSector13EtaNeg; 
Double_t ErrNumGammaSector13EtaPos; 
Double_t ErrNumGammaSector14EtaNeg; 
Double_t ErrNumGammaSector14EtaPos; 
Double_t ErrNumGammaSector15EtaNeg; 
Double_t ErrNumGammaSector15EtaPos; 
Double_t ErrNumGammaSector16EtaNeg; 
Double_t ErrNumGammaSector16EtaPos; 
Double_t ErrNumGammaSector17EtaNeg; 
Double_t ErrNumGammaSector17EtaPos;
















