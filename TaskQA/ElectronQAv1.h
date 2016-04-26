const Double_t  fPi0GGBRDPG     = 0.98823;  //DPG
const Double_t  fPi0DalitzBRDPG = 0.01174;  //DPG 0.01174
const Int_t kNbins = 42;
Double_t xbins[kNbins]={0.0,  0.1,  0.2,  0.3,  0.4,  0.5,  0.6,  0.7,  0.8,  0.9,
				1.0,  1.1,  1.2,  1.3,  1.4,  1.5,  1.6,  1.7,  1.8,  1.9,
				2.0,  2.2,  2.4,  2.6,  2.8,  3.0,  3.4,  3.8,  4.2,  4.6,
				5.0,  6.0,  7.0,  8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0,
				22.0, 25.0};
const Int_t kNRbins = 109;
Double_t xRbins[kNRbins] = {  0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 
			      10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 
			      20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 
			      30, 31, 32, 33, 34, 35, 37, 39, 41, 43, 
			      45, 47, 49, 51, 53, 55, 57, 59, 61, 63, 
			      65, 67, 69, 71, 73, 75, 77, 79, 81, 83, 
			      85, 87, 89, 91, 93, 95, 97, 99, 101, 103, 
			      105, 107, 109, 111, 113, 115, 117, 119, 121, 123, 
			      125, 127, 129, 131, 133, 135, 137, 139, 141, 143, 
			      145, 147, 149, 151, 153, 155, 157, 159, 161, 163, 
			      165, 167, 169, 171, 173, 175, 177, 179, 180};
			      
Double_t nEventsData 	= 0.0;
Double_t nEventsMC   	= 0.0;
Double_t mean_Data   	= 0.0;
Double_t mean_MC     	= 0.0;
Double_t NormFactorData = 1.0;
Double_t NormFactorMC   = 1.0;
TString fEventCutSelection    = "";
TString fGammaCutSelection    = "";
TString fClusterCutSelection  = "";
TString fElectronCutSelection = "";
TString fMesonCutSelection    = "";
TString fCutSelectionMC;
TString fEventCutSelectionMC    = "";
TString fGammaCutSelectionMC    = "";
TString fClusterCutSelectionMC  = "";
TString fElectronCutSelectionMC = "";
TString fMesonCutSelectionMC    = "";

TH1F* fEventQualityData = 0;
      
TH1F* fESD_NumberOfGoodESDTracksData = 0;

//TString collisionSystem = "";
TString textYAxisPhiP      = 	"#frac{dN_{e^{+}}}{N_{e^{+}} d#phi} (rad)^{-1}";
TString textYAxisPhiE      = 	"#frac{dN_{e^{-}}}{N_{e^{-}} d#phi} (rad)^{-1}";
TString  textYAxisPhiEP    =    "#frac{dN_{e}}{N_{e} d#phi} (rad)^{-1}";
TString textYAxisPtDistP   = 	"#frac{dN_{e^{+}}}{N_{e^{+}} dp_{T}} (GeV/c)^{-1}";
TString textYAxisPtDistE   = 	"#frac{dN_{e^{-}}}{N_{e^{-}} dp_{T}} (GeV/c)^{-1}";
TString textYAxisPtDistEP  =    "#frac{dN_{e^{#pm}}}{N_{e^{#pm}} dp_{T}} (GeV/c)^{-1}";
TString textYAxisEtaDistP  = 	"#frac{dN_{e^{+}}}{N_{e^{+}} d#eta}";
TString textYAxisEtaDistE  = 	"#frac{dN_{e^{-}}}{N_{e^{-}} d#eta}";
TString textYAxisEtaDistEP =    "#frac{dN_{e^{#pm}}}{N_{e^{#pm}} d#eta}";
TString textYAxisEtaDistGamma = "#frac{dN_{#gamma}}{N_{#gamma} d#gamma} (GeV/c)^{-1}";
TString textYAxisDCAxy = 	"#frac{dN_{e^{#pm}}}{N_{e^{#pm}} d DCA_{xy}} (cm^{-1})";
TString textYAxisDCAz =  	"#frac{dN_{e^{#pm}}}{N_{e^{#pm}} d DCA_{z}} (cm^{-1})";

TString textYAxisITSclsE      = "#frac{dN_{e^{-}}}{N_{e^{-}} dN_{ITS_{cls}}}";
TString textYAxisITSclsP      = "#frac{dN_{e^{+}}}{N_{e^{+}} dN_{ITS_{cls}}}";
TString textYAxisITSclsEP     = "#frac{dN_{e^{#pm}}}{N_{e^{#pm}} dN_{ITS_{cls}}}";


TString textYAxisTPCclsP      = "#frac{dN_{e^{+}}}{N_{e^{+}} dN_{TPCcls}}";
TString textYAxisTPCclsE      = "#frac{dN_{e^{-}}}{N_{e^{-}} dN_{TPCcls}}";
TString textYAxisTPCclsEP     = "#frac{dN_{e}}{N_{e} dN_{TPCcls}}";

TString textYAxisTPCcrossedRowsP     = "#frac{dN_{e^{+}}}{N_{e^{+}} dN_{TPCcrossedRows}}";
TString textYAxisTPCcrossedRowsE     = "#frac{dN_{e^{-}}}{N_{e^{-}} dN_{TPCcrossedRows}}";
TString textYAxisTPCcrossedRowsEP    = "#frac{dN_{e}}{N_{e} dN_{TPCcrossedRows}}";
TString pTLabel 		     = "#it{p}_{T} (GeV/#it{c})";
TString pLabel = "p (GeV/#it{c})";



//TString 

TH1F* hESDConvGammaEta_data               		= 0;
TH1F* hESDDalitzElectronAfterPt_data      		= 0;
TH1F* hESDDalitzPositronAfterPt_data      		= 0;
TH1F* hESDDalitzElectronAfterEta_data     		= 0;
TH1F* hESDDalitzPositronAfterEta_data     		= 0;
TH1F* hESDDalitzElectronAfterEtaPCut_data 		= 0;
TH1F* hESDDalitzPositronAfterEtaPCut_data 		= 0;
TH1F* hESDDalitzElectronAfterPhi_data 	  		= 0;
TH1F* hESDDalitzPositronAfterPhi_data     		= 0; 
TH1F* hESDDalitzElectronAfterNClsITS_data 		= 0;
TH1F* hESDDalitzPositronAfterNClsITS_data 		= 0;
TH2F* hESDDalitzElectronAfterNFindClsTPC_data		= 0;
TH2F* hESDDalitzPositronAfterNFindClsTPC_data		= 0;
TH1F* hESDDalitzElectronAfterNFindClsTPCPCut_data	= 0;
TH1F* hESDDalitzPositronAfterNFindClsTPCPCut_data	= 0;
TH2F* hESDDalitzElectronAfterNClsTPC_data		= 0;
TH2F* hESDDalitzPositronAfterNClsTPC_data		= 0;
TH1F* hESDDalitzElectronAfterNClsTPCPCut_data		= 0;
TH1F* hESDDalitzPositronAfterNClsTPCPCut_data		= 0;
TH2F* hESDDalitzElectronAfterNCrossedRowsTPC_data	= 0;
TH2F* hESDDalitzPositronAfterNCrossedRowsTPC_data	= 0;
TH1F* hESDDalitzElectronAfterNCrossedRowsTPCPCut_data   = 0;
TH1F* hESDDalitzPositronAfterNCrossedRowsTPCPCut_data   = 0;
TH2F* hESDDalitzPosEleAfterDCAxy_data			= 0;
TH2F* hESDDalitzPosEleAfterDCAz_data  			= 0;
TH2F* hESDDalitzElectronAfterTPCdEdxVsP_data		= 0;
TH2F* hESDDalitzPositronAfterTPCdEdxVsP_data		= 0;
TH2F* hESDDalitzElectronAfterTPCdEdxSignalVsP_data	= 0;
TH2F* hESDDalitzPositronAfterTPCdEdxSignalVsP_data	= 0;
TH2F* hESDDalitzElectronAfterTPCdEdxVsEta_data		= 0;
TH2F* hESDDalitzPositronAfterTPCdEdxVsEta_data		= 0;
TH2F* hESDDalitzElectronAfterTPCdEdxVsPhi_data		= 0;
TH2F* hESDDalitzPositronAfterTPCdEdxVsPhi_data		= 0;
TH1F* hESDMotherPhi_data				= 0;
TH2F* hESDEposEnegPsiPairDPhi_data			= 0;
TH2F* hESDEposEnegInvMassPt_data			= 0;
TH2F* hESDEposEnegLikeSignBackInvMassPt_data		= 0;
TH2F* hTPCdEdxSigafter_data 				= 0;
TH2F* hTPCdEdxafter_data    				= 0;
TH2F* hTPCdEdxSigbefore_data                            = 0;
TH2F* hTPCdEdxbefore_data				= 0;
TH2F* hTPCEledEdxSigafter_data				= 0;
TH2F* hTPCEledEdxafter_data				= 0;



TH1F* fEventQualityMC 					= 0;
TH1F* fESD_NumberOfGoodESDTracksMC  			= 0;
TH1F * hESDConvGammaEta_mc                              = 0;
TH2F * hMCPi0EposEnegInvMassPt_mc			= 0;
TH1F* hMCConvGammaR	                                = 0;
TH1F* hMCConvGammaPt                                    = 0;
TH1F* hMCAllGammaPt					= 0;
TH2F* hMCConvGammaPtR					= 0;
TH1F* hMCAllPositronsPt					= 0;
TH1F* hMCAllElectronsPt					= 0;
TH1F* hMCDecayPositronPi0Pt				= 0;
TH1F* hMCDecayElectronPi0Pt				= 0;
TH1F* hMCPi0DalitzGammaPt 				= 0;
TH1F* hMCConvGammaPi0Pt					= 0;
TH1F* hMCAllGammaPi0Pt					= 0;
TH1F* hESDDalitzElectronAfterPt_mc 			= 0;
TH1F* hESDDalitzPositronAfterPt_mc  			= 0;
TH1F* hESDDalitzElectronAfterEta_mc 			= 0;
TH1F* hESDDalitzPositronAfterEta_mc 			= 0;
TH1F* hESDDalitzElectronAfterEtaPCut_mc               	= 0;
TH1F* hESDDalitzPositronAfterEtaPCut_mc               	= 0;
TH1F* hESDDalitzElectronAfterPhi_mc 			= 0;
TH1F* hESDDalitzPositronAfterPhi_mc 			= 0;
TH1F* hESDDalitzElectronAfterNClsITS_mc			= 0;
TH1F* hESDDalitzPositronAfterNClsITS_mc			= 0;
TH2F* hESDDalitzElectronAfterNFindClsTPC_mc		= 0;
TH2F* hESDDalitzPositronAfterNFindClsTPC_mc		= 0;
TH1F* hESDDalitzElectronAfterNFindClsTPCPCut_mc		= 0;
TH1F* hESDDalitzPositronAfterNFindClsTPCPCut_mc		= 0;
TH2F* hESDDalitzElectronAfterNClsTPC_mc			= 0;
TH2F* hESDDalitzPositronAfterNClsTPC_mc			= 0;
TH1F* hESDDalitzElectronAfterNClsTPCPCut_mc		= 0;
TH1F* hESDDalitzPositronAfterNClsTPCPCut_mc		= 0;
TH2F* hESDDalitzElectronAfterNCrossedRowsTPC_mc		= 0;
TH2F* hESDDalitzPositronAfterNCrossedRowsTPC_mc		= 0;
TH1F* hESDDalitzElectronAfterNCrossedRowsTPCPCut_mc	= 0;
TH1F* hESDDalitzPositronAfterNCrossedRowsTPCPCut_mc     = 0;
TH2F* hESDDalitzPosEleAfterDCAxy_mc			= 0;
TH2F* hESDDalitzPosEleAfterDCAz_mc  			= 0;
TH2F* hESDDalitzElectronAfterTPCdEdxVsP_mc		= 0;
TH2F* hESDDalitzPositronAfterTPCdEdxVsP_mc		= 0;
TH2F* hESDDalitzElectronAfterTPCdEdxSignalVsP_mc	= 0;
TH2F* hESDDalitzPositronAfterTPCdEdxSignalVsP_mc	= 0;
TH2F* hESDDalitzElectronAfterTPCdEdxVsEta_mc		= 0;
TH2F* hESDDalitzPositronAfterTPCdEdxVsEta_mc		= 0;
TH2F* hESDDalitzElectronAfterTPCdEdxVsPhi_mc		= 0;
TH2F* hESDDalitzPositronAfterTPCdEdxVsPhi_mc		= 0;
TH1F* hESDMotherPhi_mc					= 0;
TH2F* hESDEposEnegPsiPairDPhi_mc			= 0;
TH2F* hESDEposEnegInvMassPt_mc				= 0;
TH2F* hESDEposEnegLikeSignBackInvMassPt_mc		= 0;
TH2F* hTPCdEdxSigafter_mc 				= 0;
TH2F* hTPCdEdxafter_mc    				= 0;
TH2F* hESDEposEnegTruePi0DalitzInvMassPt   		= 0;
TH2F* hESDEposEnegTrueEtaDalitzInvMassPt                = 0;
TH2F* hESDEposEnegTruePhotonInvMassPt			= 0;
TH2F* hESDEposEnegInvMassPt				= 0;
TH2F* hESDEposEnegTrueInvMassPt				= 0;

TH1F* hESDEposEnegTruePi0DalitzPt 			= 0;
TH1F* hESDEposEnegTrueEtaDalitzPt 			= 0; 
TH1F* hESDEposEnegTruePhotonPt    			= 0; 
TH1F* hESDEposEnegPt              			= 0;
TH1F* hESDTrueEposEnegPt				= 0;

TH1F* hESDEposEnegTruePi0DalitzInvMass 			= 0;
TH1F* hESDEposEnegTrueEtaDalitzInvMass 			= 0; 
TH1F* hESDEposEnegTruePhotonInvMass    			= 0; 
TH1F* hESDEposEnegTrueInvMass				= 0;
TH1F* hESDEposEnegInvMass 				= 0;
TH1F* hESDTrueEposEnegInvMass              		= 0;
TH2F* hTPCdEdxSigbefore_mc      	                = 0;
TH2F* hTPCdEdxbefore_mc					= 0;

	



TH1F* hESDTrueConvGammaRMC 		  		= 0;
TH1F* hESDTrueConvGammaR 				= 0;
TH1F* hESDTrueConvGammaRvsEffi  			= 0;
TH1F* hESDTrueConvGammaPt				= 0;
TH1F* hESDTrueConvGammaPtvsEffi     			= 0;
TH1F* hMCGammaPtvsConvProb 	   			= 0;
TH2F* hESDEposEnegInvMassPi0Pt_mc			= 0;
TH1F* hESDEposEnegInvMassPi0Below			= 0;
TH1F* hESDEposEnegInvMassPi0Above                       = 0;
TH1F* hESDDalitzElectronPositronAfterClsITS_mc          = 0;
TH1F* hESDDalitzElectronPositronAfterClsITS_data        = 0;
TH2F* hTPCEledEdxSigafter_mc				= 0;
TH2F* hTPCEledEdxafter_mc 				= 0;
    

const Int_t nBins = 7;
TH1F*  hMCConvGammaEffiVsPtProb[nBins];
TString gammaCutEffiLegend = "|#eta|<0.9, p_{T,min}=0.05 GeV/c, R_{#gamma conv,min}=5.0cm}";
TString electronCutEffiLegend = "|#eta|<0.9";
		
TH1F* hESDTruePositronPt				= 0;
TH1F* hESDTrueElectronPt				= 0;
       
TH1F* hESDTruePositronPtvsEffi  	= 0;
TH1F* hESDTrueElectronPtvsEffi  	= 0;
TH1F* hMCAllPositronsPtScaled   	= 0;
TH1F* hMCAllElectronsPtScaled   	= 0;
TH1F* hESDTruePositronPtScaled  	= 0;
TH1F* hESDTrueElectronPtScaled  	= 0;
       
TH1F* hMCAllPositronsPtScaledBR   	= 0;
TH1F* hMCAllElectronsPtScaledBR   	= 0;
TH1F* hESDTruePositronPtScaledBR  	= 0;
TH1F* hESDTrueElectronPtScaledBR  	= 0;
     
TH1D* fHistoMCMesonPt	       		= 0;   
TH1D* fHistoMCMesonDalitzPt     	= 0;
TH1D* fHistoMCMesonGGPt         	= 0;
      
TH2F* hESDEposEnegTruePi0DalitzPsiPairDPhi	= 0;
TH2F* hESDEposEnegTrueEtaDalitzPsiPairDPhi	= 0; 
TH2F* hESDEposEnegTruePhotonPsiPairDPhi		= 0; 
       
TH1F*	hESDTruePi0DalitzElectronPt		= 0;
TH1F*	hESDTruePi0DalitzPositronPt		= 0;
TH1F* 	hESDTruePi0DalitzPositronPtvsEffi 	= 0;
TH1F* 	hESDTruePi0DalitzElectronPtvsEffi 	= 0;  
TH1F*   hESDTruePi0DalitzConvGammaPt            = 0;
TH1F*   hHistoTruePi0DalitzClusGammaPt          = 0;
TH1F*   hESDTruePi0DalitzConvGammaPtvsEffi      = 0;
TH1F*   hHistoTruePi0DalitzClusGammaPtvsEffi    = 0;
TH1F*   hESDTruePi0DalitzGammaPtvsEffi          = 0;
TH1F*   hMCPi0EposEnegInvMass                   = 0;
TH2D*   hESDConvGammaZR_data			= 0;
TH2D*   hESDConvGammaXY_data                    = 0;
THnSparseF* sESDConvGammaZR_mc			= 0;
THnSparseF* sESDConvGammaXY_mc			= 0;
THnSparseF* sESDConvGammaZR_data		= 0;
THnSparseF* sESDConvGammaXY_data		= 0;

TString 	arrayNamesRBins[13]=	{"Beam Pipe", 
									"SPD 1", 
									"SPD 2",
									"Thermal shield/Support between SPD/SDD", 
									"SDD 1", 
									"SDD 2", 
									"Thermal shield/Support between SDD/SSD", 
									"SSD 1", 
									"SSD 2", 
									"Air + TPC in. cont. vessel + CO_{2}", 
									"CO_{2} + TPC in. field cage vessel", 
									"TPC rods + Ne: CO_{2}: N_{2}", 
									"Ne: CO_{2}: N_{2}"};


		

void DrawAliceLogoPerformancePlot(Float_t startX, Float_t startY, Float_t widthLogo, Float_t textHeight, Float_t decrease, TString date, TString collisionSystem, TString textGenerator, TString textPeriod, Double_t xLengthCanvas, Double_t yLengthCanvas);
//TGraph2D* ConvertTH2FtoTGraph(TH2F* h2);

/*TGraph2D* ConvertTH2FtoTGraph(TH2F* h2){
  
   const Int_t xNBin = h2->GetNbinsX();
   const Int_t yNBin = h2->GetNbinsY();
   
   Double_t xBinVal[xNBin];
   Double_t yBinVal[yNBin];
  
  for(Int_t iBin=0; iBin < xNBin; iBin++){
     
		 xBinVal[iBin] = h2->GetXaxis()->GetBinCenter(iBin+1);
		 yBinVal[iBin] = h2->GetBinContent(iBin+1);
		 
     
  }
   
   //TGraph* graph  = new TGraph(xNBin,xBinVal,yBinVal);
   TGraph2D* graph = new TGraph2D(h2);
   
   return graph;
  
}*/


       