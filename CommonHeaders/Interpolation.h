extern TRandom*         gRandom;	
extern TBenchmark*      gBenchmark;
extern TSystem*         gSystem;
extern TMinuit*         gMinuit;











void CalcRpPb(TGraphAsymmErrors* PPSpectrumSystErr, TGraphAsymmErrors*  PPSpectrumStatErr, TGraphAsymmErrors* pPbSpectrumSystErr, TGraphAsymmErrors* pPbSpectrumStatErr,
	      TGraphAsymmErrors** graphRpPbSystErr, TGraphAsymmErrors** graphRpPbStatErr);

TGraphAsymmErrors* CalculateSystErrors(TGraphErrors* spectrum,TGraphAsymmErrors* g1,TGraphAsymmErrors* g2);

TF1* RebinWithFitToTGraph(TGraphAsymmErrors *spectrum, TGraphAsymmErrors** newSpectrum, TGraphAsymmErrors *newBins, TString FitType,Double_t minPt, Double_t maxPt, Double_t* parameters,Double_t probability);
TGraphErrors *GetInterpolSpectrum2D(TGraphErrors *g1, TGraphErrors *g2,Double_t d1, Double_t d2,Double_t dSqrts,TString method);
TGraphAsymmErrors* GetChargeParticlesRpPb(TString typeErr);
TGraphErrors *ConvertTGraphAsymmErrorstoTGraphErrors(TGraphAsymmErrors* g1);
TGraphAsymmErrors* ConvertTGraphErrorstoTGraphAsymmErrors(TGraphErrors* g1);
TH1F* GetPPReferenceFromPythia(TString fileName);
void SavingFiles(TString outputDir);



Double_t	xSection2760GeVppINEL    =  62.8*1e-3;
Double_t        xSectionpPb5023GeVINEL   =  70*1e-3;
Double_t	xSection7TeVINEL         =  73.2*1e-3;
//Double_t        recalcBarn               =  1e12; 
Double_t 	fNcoll = 6.9;
Double_t 	fTpPb     =   0.0983e3*(1/recalcBarn);
Double_t 	fTpPbErr  =   0.0035e3*(1/recalcBarn);
	 

TGraphAsymmErrors* graphErrosInterPolation5023GeVBinDalitzStatErr;
TGraphAsymmErrors* graphErrosInterPolation5023GeVBinDalitzSystErr;
TGraphAsymmErrors* graphErrosInterPolation5023GeVBinDalitzSystWOMatErr;
TGraphAsymmErrors* graphInvYieldPi0DalitzpPb5023GeVXShiftedSystWOMatErr;
TGraphAsymmErrors* graphInvYieldPi0DalitzpPb5023GeVYShiftedSystWOMatErr;

TGraphAsymmErrors* graphErrosInterPolation5023GeVBinPCMStatErr;
TGraphAsymmErrors* graphErrosInterPolation5023GeVBinPCMSystErr;	
TGraphAsymmErrors* graphErrosInterPolation5023GeVBinPCMSystWOMatErr = NULL;
TGraphAsymmErrors* graphInvYieldPi0PCMpPb5023GeVXShiftedSystWOMatErr = NULL;
TGraphAsymmErrors* graphInvYieldPi0PCMpPb5023GeVYShiftedSystWOMatErr = NULL;
TGraphAsymmErrors* graphRpPbDalitzSystErr;
TGraphAsymmErrors* graphRpPbDalitzStatErr;
TGraphAsymmErrors* graphRpPbPCMStatErr;
TGraphAsymmErrors* graphRpPbPCMSystErr;
TGraphAsymmErrors* graphRpPbPHOSStatErr;
TGraphAsymmErrors* graphRpPbPHOSSystErr;
TGraphErrors*      graphAlphaDalitz;
TGraphErrors*      graphAlphaPCM;
TGraphErrors**     graphPtvsSqrtsDalitz;
TGraphErrors**     graphPtvsSqrtsPCM; 
TGraphErrors**     gPtvsEnergiesPCM;
TGraphErrors**     gPtvsEnergiesDalitz;
TF1**              fPowerlawDalitz;
TF1**		   fPowerlawPCM;

TGraphAsymmErrors* graphAsymmChargedParticlesRpPbSystErr;
TGraphAsymmErrors* graphAsymmChargedParticlesRpPbStatErr;
TGraphAsymmErrors*  graphChargedPionRpPbSystErr;
TGraphAsymmErrors*  graphChargedPionRpPbStatErr;
TGraph*      graphPi0ESP09sPi0KKP;
TGraph*      graphPi0ESP09sPi0AKK;
TGraph*      graphPi0DSS5000;
TGraphAsymmErrors*      graphAsymmErrorsPi0DSS5000;
TGraph*      graphPi0CGC;
//TGraphErrors* 		t[];

/*Labels*/

TString invYieldLabel = "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}N_{#pi^{0}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c}^{2})";
TString pTLabel = "#it{p}_{T} (GeV/#it{c})";
TString RpPbLabel= "#it{R}_{pPb}";
TString energyLabel = "p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV";
TString DalitzLabel = "#pi^{0} #rightarrow e^{+}e^{-}#gamma";
TString PCMLabel    = "#pi^{0} #rightarrow #gamma#gamma #rightarrow e^{+}e^{-}e^{+}e^{-}";
TString ChargedPionsLabel = "#pi^{+} + #pi^{-}";
TString thesisPlotLabel = "";

enum{kDalitz,kPCM,kPHOS,kPion,kKaon,kProton,kNAna};
Int_t colorsArray[kNAna];
      










TGraphAsymmErrors* GetChargeParticlesRpPb(TString typeErr){

  // Plot: p8424_d3x1y1
  double p8424_d3x1y1_xval[] = { 0.525, 0.575, 0.625, 0.675, 0.725, 0.775, 0.825, 0.875, 0.925, 
    0.975, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 
    1.95, 2.1, 2.3, 2.5, 2.7, 2.9, 3.1, 3.3, 3.5, 3.7, 
    3.9, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 
    10.5, 11.5, 12.5, 13.5, 15.0, 18.0 };
  double p8424_d3x1y1_xerrminus[] = { 0.025000000000000022, 0.02499999999999991, 0.025000000000000022, 0.025000000000000022, 0.025000000000000022, 0.025000000000000022, 0.02499999999999991, 0.025000000000000022, 0.025000000000000022, 
    0.025000000000000022, 0.050000000000000044, 0.04999999999999982, 0.050000000000000044, 0.050000000000000044, 0.050000000000000044, 0.050000000000000044, 0.04999999999999982, 0.050000000000000044, 0.050000000000000044, 
    0.050000000000000044, 0.10000000000000009, 0.09999999999999964, 0.10000000000000009, 0.10000000000000009, 0.10000000000000009, 0.10000000000000009, 0.09999999999999964, 0.10000000000000009, 0.10000000000000009, 
    0.10000000000000009, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.5, 0.5, 0.5, 
    0.5, 0.5, 0.5, 0.5, 1.0, 2.0 };
  double p8424_d3x1y1_xerrplus[] = { 0.025000000000000022, 0.025000000000000022, 0.025000000000000022, 0.02499999999999991, 0.025000000000000022, 0.025000000000000022, 0.025000000000000022, 0.025000000000000022, 0.02499999999999991, 
    0.025000000000000022, 0.050000000000000044, 0.050000000000000044, 0.050000000000000044, 0.04999999999999982, 0.050000000000000044, 0.050000000000000044, 0.050000000000000044, 0.050000000000000044, 0.04999999999999982, 
    0.050000000000000044, 0.10000000000000009, 0.10000000000000009, 0.10000000000000009, 0.09999999999999964, 0.10000000000000009, 0.10000000000000009, 0.10000000000000009, 0.10000000000000009, 0.09999999999999964, 
    0.10000000000000009, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.5, 0.5, 0.5, 
    0.5, 0.5, 0.5, 0.5, 1.0, 2.0 };

    double p8424_d3x1y1_yval[] = { 0.5824, 0.5987, 0.6139, 0.6337, 0.6514, 0.6716, 0.6878, 0.6989, 0.7118, 
    0.722, 0.7485, 0.7822, 0.8046, 0.8306, 0.8539, 0.8811, 0.8947, 0.9145, 0.9431, 
    0.9584, 0.9844, 1.01, 1.034, 1.051, 1.051, 1.08, 1.083, 1.104, 1.102, 
    1.111, 1.1, 1.091, 1.111, 1.058, 1.035, 1.087, 1.064, 1.064, 1.11, 
    0.9183, 1.144, 1.061, 1.17, 0.9558, 1.121 };
  double p8424_d3x1y1_yerrminus[] = { 0.05600723167591842, 0.05670881765651617, 0.056910631695668255, 0.058712264476853564, 0.060314011639087645, 0.06221575363201831, 0.06371765846294103, 0.06472232999514155, 0.0659245781177248, 
    0.06692697512961422, 0.06931623186527092, 0.07251992829560713, 0.07452684080249208, 0.07693438762997987, 0.07914271918502674, 0.08165512843661445, 0.08296969326205805, 0.08478519918004557, 0.08740583504549339, 
    0.08892429364352579, 0.09128335006998813, 0.09413288479590966, 0.09618731725128839, 0.09725224933131367, 0.09732933781753578, 0.10031948963187563, 0.1004987562112089, 0.10259142264341595, 0.10270345661174213, 
    0.1039471019317037, 0.10270345661174213, 0.10196568050084304, 0.10932977636490435, 0.10413932974625868, 0.10235721762533408, 0.1078563859954523, 0.10412012293500235, 0.1077032961426901, 0.11764777940955791, 
    0.10736125930707034, 0.13917614738165443, 0.14509651959988565, 0.1776091213873882, 0.14606854555310667, 0.17462244987400674 };
  double p8424_d3x1y1_yerrplus[] = { 0.05600723167591842, 0.05670881765651617, 0.056910631695668255, 0.058712264476853564, 0.060314011639087645, 0.06221575363201831, 0.06371765846294103, 0.06472232999514155, 0.0659245781177248, 
    0.06692697512961422, 0.06931623186527092, 0.07251992829560713, 0.07452684080249208, 0.07693438762997987, 0.07914271918502674, 0.08165512843661445, 0.08296969326205805, 0.08478519918004557, 0.08740583504549339, 
    0.08892429364352579, 0.09128335006998813, 0.09413288479590966, 0.09618731725128839, 0.09725224933131367, 0.09732933781753578, 0.10031948963187563, 0.1004987562112089, 0.10259142264341595, 0.10270345661174213, 
    0.1039471019317037, 0.10270345661174213, 0.10196568050084304, 0.10932977636490435, 0.10413932974625868, 0.10235721762533408, 0.1078563859954523, 0.10412012293500235, 0.1077032961426901, 0.11764777940955791, 
    0.10736125930707034, 0.13917614738165443, 0.14509651959988565, 0.1776091213873882, 0.14606854555310667, 0.17462244987400674 };
  double p8424_d3x1y1_ystatminus[] = { 9.0E-4, 0.001, 0.0011, 0.0012, 0.0013, 0.0014, 0.0015, 0.0017, 0.0018, 
    0.0019, 0.0015, 0.0017, 0.002, 0.0023, 0.0026, 0.003, 0.0034, 0.0038, 0.0043, 
    0.0047, 0.0039, 0.005, 0.006, 0.007, 0.008, 0.008, 0.01, 0.011, 0.012, 
    0.014, 0.012, 0.014, 0.017, 0.021, 0.026, 0.032, 0.029, 0.04, 0.055, 
    0.064, 0.089, 0.107, 0.141, 0.1159, 0.138 };
  double p8424_d3x1y1_ystatplus[] = { 9.0E-4, 0.001, 0.0011, 0.0012, 0.0013, 0.0014, 0.0015, 0.0017, 0.0018, 
    0.0019, 0.0015, 0.0017, 0.002, 0.0023, 0.0026, 0.003, 0.0034, 0.0038, 0.0043, 
    0.0047, 0.0039, 0.005, 0.006, 0.007, 0.008, 0.008, 0.01, 0.011, 0.012, 
    0.014, 0.012, 0.014, 0.017, 0.021, 0.026, 0.032, 0.029, 0.04, 0.055, 
    0.064, 0.089, 0.107, 0.141, 0.1159, 0.138 };
  int p8424_d3x1y1_numpoints = 45;
  
  TGraphAsymmErrors* p8424_d3x1y1;
  
  if( typeErr.CompareTo("Syst") == 0 ) {
	p8424_d3x1y1 = new TGraphAsymmErrors(p8424_d3x1y1_numpoints, p8424_d3x1y1_xval, p8424_d3x1y1_yval, p8424_d3x1y1_xerrminus, p8424_d3x1y1_xerrplus, p8424_d3x1y1_yerrminus, p8424_d3x1y1_yerrplus);
  } else if ( typeErr.CompareTo("Stat") == 0 ) {
	  
        p8424_d3x1y1 = new TGraphAsymmErrors(p8424_d3x1y1_numpoints, p8424_d3x1y1_xval, p8424_d3x1y1_yval, p8424_d3x1y1_xerrminus, p8424_d3x1y1_xerrplus, p8424_d3x1y1_ystatminus, p8424_d3x1y1_ystatplus);
      
  }
  p8424_d3x1y1->SetName("/HepData/8424/d3x1y1");
  p8424_d3x1y1->SetTitle("/HepData/8424/d3x1y1");
  
  return p8424_d3x1y1;
  
  
}


	
TGraphErrors* RemovePointsFromGraph(TGraphErrors *graph, Int_t 
NbPoints){

   TGraphErrors* dummyGraph = (TGraphErrors*)	graph->Clone();
	Double_t * xValue = dummyGraph->GetX();
	Double_t * yValue = dummyGraph->GetY();
	Double_t* xError = dummyGraph->GetEX();
	Double_t* yError = dummyGraph->GetEY();

	Int_t nPoints = dummyGraph->GetN();
	Int_t nPoints_new = nPoints  - NbPoints;
	for (Int_t i = NbPoints; i < nPoints; i++){
		yValue[i-NbPoints] = yValue[i];
		xValue[i-NbPoints] = xValue[i];
		xError[i-NbPoints] = xError[i];
		yError[i-NbPoints] = yError[i];
	}
	TGraphErrors* returnGraph = new TGraphErrors(nPoints_new,xValue,yValue, xError,yError);
	return returnGraph;
}

	
	
	

TF1* RebinWithFitToTGraph(TGraphAsymmErrors *spectrum, TGraphAsymmErrors** newSpectrum, TGraphAsymmErrors *newBins, TString FitType,Double_t minPt, Double_t maxPt, Double_t* parameters){
  
  

        TF1 * CurrentFit = new TF1();
	
	
	Int_t     nOldPoints       = spectrum->GetN();
	Double_t *xOldValue 	   = spectrum->GetX();
	Double_t *yOldValueErrlow  = spectrum->GetEYlow();
	Double_t *yOldValueErrhigh = spectrum->GetEYhigh();
	Double_t *yOldValue        = spectrum->GetY();
	
	///////////////////////////////////////////////////////////
	
	Int_t nNewPoints            = newBins->GetN();
	Double_t *xNewValue         = newBins->GetX();
	Double_t *xNewValueErrlow   = newBins->GetEXlow();
	Double_t *xNewValueErrhigh  = newBins->GetEXhigh();
	
	
	//////////////////////////////////////////////////////////////    
	
	
	Double_t  yNewValueRelErrlow[nNewPoints];
	Double_t  yNewValueRelErrhigh[nNewPoints];
	
			
	
	(*newSpectrum)        = new TGraphAsymmErrors (nNewPoints);
	
	
	
	
	
	
	if( FitType.BeginsWith("l") || FitType.BeginsWith("L") ) {
	
	       
	       CurrentFit = FitObject("l","fitInvCrossSectionPi0","Pi0");
	       
	       

	       CurrentFit->SetRange(minPt,maxPt);
 	       CurrentFit->SetParameters(parameters[0],parameters[1],parameters[2]); // standard
	
	
	         //TFitResultPtr resultfit = spectrum->Fit(CurrentFit,"SQNRME+","",minPt,maxPt);
	         //spectrum->Fit(CurrentFit,"QNRMEX0+","",minPt,maxPt);  //One time
	         spectrum->Fit(CurrentFit,"SQNRME+","",minPt,maxPt);  //One time
	       
	     
	} else if (FitType.BeginsWith("tcm") || FitType.BeginsWith("TCM") ){
	  
	  
	      
		
	       CurrentFit = FitObject("tcm","fitInvCrossSectionPi0","Pi0");
	       
	       

	       CurrentFit->SetRange(minPt,maxPt);
 	       CurrentFit->SetParameters(parameters[0],parameters[1],parameters[2],parameters[3],parameters[4]); // standard
	
	
	         //TFitResultPtr resultfit = spectrum->Fit(CurrentFit,"SQNRME+","",minPt,maxPt);
	         //spectrum->Fit(CurrentFit,"QNRMEX0+","",minPt,maxPt);  //One time
	       spectrum->Fit(CurrentFit,"SQNRME+","",minPt,maxPt);  //One time
	         
	    
	}
	
	cout<<"iPoint"<<"\t"<<"x"<<"\t"<<"y"<<"\t"<<"xElow"<<"\t"<<"xEhigh"<<"\t"<<"yElow"<<"\t"<<"yEhigh"<<endl;
    
        Double_t decisionBoundary = 0.0000001;
        
	for( Int_t iPoint = 0; iPoint < nNewPoints; iPoint++){
	  
	  
		    Double_t yEval;
		    
		    Int_t index = 0;
		    
		    while(  ( xNewValue[iPoint] - xOldValue[index] ) > decisionBoundary  && index < nOldPoints-1 ) {index++;}
		    
		    cout<<"xNewValue "<<xNewValue[iPoint]<<" xOldValue "<<xOldValue[index]<<"  "<<(xNewValue[iPoint] - xOldValue[index])<<" nOldPoints "<<nOldPoints-1<<" index "<<index<<endl;
		    
		    
		     if( index == 0 ){
		
			yNewValueRelErrlow[iPoint]  = yOldValueErrlow[index]  / yOldValue[index];
			yNewValueRelErrhigh[iPoint] = yOldValueErrhigh[index] / yOldValue[index];
			
			
		      } else if (  TMath::Abs(  xNewValue[iPoint] - xOldValue[index-1] ) < TMath::Abs(  xNewValue[iPoint] - xOldValue[index] ) ) {
		
		
	      		yNewValueRelErrlow[iPoint]  = yOldValueErrlow[index-1] / yOldValue[index-1];
			yNewValueRelErrhigh[iPoint] = yOldValueErrhigh[index-1] / yOldValue[index-1];
			
			
	      
		     } else {
		       
			yNewValueRelErrlow[iPoint]  = yOldValueErrlow[index]  / yOldValue[index]; 
			yNewValueRelErrhigh[iPoint] = yOldValueErrhigh[index] / yOldValue[index];
			
		     }
		     
	      
		     yEval = CurrentFit->Eval(xNewValue[iPoint]);
		     
		     Double_t yNewValueErrlow     = yEval * yNewValueRelErrlow[iPoint];
		     Double_t yNewValueErrhigh    = yEval * yNewValueRelErrhigh[iPoint];
		    
		    
		    (*newSpectrum)->SetPoint(iPoint,xNewValue[iPoint],yEval);
		    (*newSpectrum)->SetPointError(iPoint,xNewValueErrlow[iPoint],xNewValueErrhigh[iPoint],yNewValueErrlow,yNewValueErrhigh);
		    
		    cout<<iPoint<<"\t"<<xNewValue[iPoint]<<"\t"<<yEval<<"\t"<<xNewValueErrlow[iPoint]<<"\t"<<xNewValueErrhigh[iPoint]<<"\t"<<yNewValueErrlow<<"\t"<<yNewValueErrhigh<<endl;
    
		   
		    
		      
	}
	
	
	  
	
	return CurrentFit;
	
	
	
	
}

TGraphAsymmErrors* CalculateSystErrors(TGraphErrors* spectrum,TGraphAsymmErrors* g1,TGraphAsymmErrors* g2){
  
  
	Int_t nPoints = spectrum->GetN();
	Double_t* xBins = spectrum->GetX();
	Double_t* yBins = spectrum->GetY();
	
	
	Double_t *g1yErrlow  = g1->GetEYlow();
	Double_t *g1xErrlow  = g1->GetEXlow();
	Double_t *g1xErrhigh = g1->GetEXhigh();
	Double_t *g1yVal     = g1->GetY();
	
	
	Double_t *g2yErrlow  = g2->GetEYlow();
	Double_t *g2yVal  = g2->GetY();
	
	
	

	
	TGraphAsymmErrors* graphErr = new TGraphAsymmErrors(nPoints);
    
	
	cout<<"iPoin\tx\ty\tXerrLow\tXerrHigh\tYerrLow\tYerrHigh"<<endl;;
	
	for(Int_t iPoint = 0; iPoint < nPoints; iPoint++){
	  
		Double_t g1yrelErr = g1yErrlow[iPoint]/g1yVal[iPoint];
		Double_t g2yrelErr = g2yErrlow[iPoint]/g2yVal[iPoint];
		Double_t ySysErr = 0.0;
		
		
		if(  g1yrelErr >= g2yrelErr ){
		  ySysErr = yBins[iPoint]*g1yrelErr;
		} else {
		  ySysErr = yBins[iPoint]*g2yrelErr;
		}
	  
		graphErr->SetPoint(iPoint,xBins[iPoint],yBins[iPoint]);
		graphErr->SetPointError(iPoint,g1xErrlow[iPoint],g1xErrhigh[iPoint],ySysErr,ySysErr);
		
		cout<<iPoint<<"\t"<<xBins[iPoint]<<"\t"<<yBins[iPoint]<<"\t"<<g1xErrlow[iPoint]<<"\t"<<g1xErrhigh[iPoint]<<"\t"<<ySysErr<<"\t"<<ySysErr<<endl;
		
	  
	}
      
  return graphErr;
}




void CalcRpPb(TGraphAsymmErrors* PPSpectrumSystErr, TGraphAsymmErrors*  PPSpectrumStatErr, TGraphAsymmErrors* pPbSpectrumSystErr, TGraphAsymmErrors* pPbSpectrumStatErr,
	      TGraphAsymmErrors** graphRpPbSystErr, TGraphAsymmErrors** graphRpPbStatErr){
	      


	
	 //Computing RpPb using the interpolated pp@5.023 TeV and PCM pp@7 TeV spectra
	 
	//Double_t        xSectionpPb5023GeVINEL   =  70*1e-3;  //Dangerous to should be put in a library
	
	cout<<"material Error"<<endl;
	
	Int_t nPoints = pPbSpectrumStatErr->GetN();

	
	(*graphRpPbStatErr) = new TGraphAsymmErrors( nPoints);
	(*graphRpPbSystErr)  = new TGraphAsymmErrors( nPoints);
	
	
	Double_t *xBins      =  pPbSpectrumStatErr->GetX();
	Double_t *xErrlow    =  pPbSpectrumStatErr->GetEXlow();
	Double_t *xErrhigh   =  pPbSpectrumStatErr->GetEXhigh();
	
	Double_t *ypPbBins   =  pPbSpectrumStatErr->GetY();
	Double_t *yPPBins    =  PPSpectrumStatErr->GetY();
	
	
	Double_t *ypPbStatErrlow      = pPbSpectrumStatErr->GetEYlow();
	//Double_t *ypPbStatErrhigh     = pPbSpectrumStatErr->GetEYlow();
	
	
	Double_t *ypPbSystErrlow    = pPbSpectrumSystErr->GetEYlow();
	Double_t *ypPbSystErrhigh   = pPbSpectrumSystErr->GetEYlow();
	
	
	Double_t *yPPStatErrlow  = PPSpectrumStatErr->GetEYlow();
	//Double_t *yPPStatErrhigh = PPSpectrumStatErr->GetEYlow();
	
	
	Double_t *yPPSystErrlow  = PPSpectrumSystErr->GetEYlow();
	Double_t *yPPSystErrhigh = PPSpectrumSystErr->GetEYlow();
	
	Double_t *RpPb   = new Double_t[nPoints];
	
	
	
	
	
	 Double_t fNcoll = 6.9;
	
	 Double_t fTpPb     =   0.0983e3*(1/recalcBarn);
	 Double_t fTpPbErr  =   0.0035e3*(1/recalcBarn);
	 
	
	
	
	for(Int_t iPoint = 0; iPoint < nPoints; iPoint++){
	  
	  RpPb[iPoint] =  ypPbBins[iPoint] / (fTpPb* yPPBins[iPoint]);
	  
	 (*graphRpPbSystErr)->SetPoint( iPoint, xBins[iPoint], RpPb[iPoint]);
	 (*graphRpPbStatErr)->SetPoint( iPoint, xBins[iPoint], RpPb[iPoint] );
	
	 
	 Double_t errYStat = 	 pow( pow( ypPbStatErrlow[iPoint]/ypPbBins[iPoint], 2. ) + pow( yPPStatErrlow[iPoint]/yPPBins[iPoint],  2.), 0.5)* RpPb[iPoint];
	 Double_t errYSystlow =  pow( pow( ypPbSystErrlow[iPoint]/ypPbBins[iPoint], 2. ) + pow( yPPSystErrlow[iPoint]/yPPBins[iPoint]  ,2. ) + pow ((fTpPbErr/fTpPb) , 2),   0.5) * RpPb[iPoint];
	 Double_t errYSysthigh = pow( pow( ypPbSystErrhigh[iPoint]/ypPbBins[iPoint],2. ) + pow( yPPSystErrhigh[iPoint]/yPPBins[iPoint], 2. ) + pow ((fTpPbErr/fTpPb) , 2),   0.5) * RpPb[iPoint];
	 
	 
	 
	 
	 (*graphRpPbStatErr)->SetPointError(iPoint, xErrlow[iPoint],xErrhigh[iPoint],errYStat,errYStat);
	 (*graphRpPbSystErr)->SetPointError(iPoint, xErrlow[iPoint],xErrhigh[iPoint],errYSystlow,errYSysthigh);
		 
	 
	}
	
	
	
}

TGraphAsymmErrors* ConvertTGraphErrorstoTGraphAsymmErrors(TGraphErrors* g1){
  
  Int_t nPoints = g1->GetN();
  Double_t *xValue = g1->GetX();
  Double_t *yValue = g1->GetY();
  
  TGraphAsymmErrors* gNew = new TGraphAsymmErrors(nPoints);
  
  
    for(Int_t iPoint = 0; iPoint < nPoints; iPoint++){
      
	    Double_t xError,yError;
	    
	    
	      xError = g1->GetErrorX(iPoint);
	      yError = g1->GetErrorY(iPoint);
	      
	      gNew->SetPoint(iPoint,xValue[iPoint],yValue[iPoint]);
	      gNew->SetPointError(iPoint,xError,xError,yError,yError);
	      
      
    }
  
  return gNew;
  
}


TGraphErrors *ConvertTGraphAsymmErrorstoTGraphErrors(TGraphAsymmErrors* g1){
  
  Int_t nPoints = g1->GetN();
  Double_t *xValue = g1->GetX();
  Double_t *yValue = g1->GetY();
  
  TGraphErrors* gNew = new TGraphErrors(nPoints);
  
  
    for(Int_t iPoint = 0; iPoint < nPoints; iPoint++){
      
	    Double_t xError,yError;
	    
	    
	      xError = g1->GetErrorX(iPoint);
	      yError = g1->GetErrorY(iPoint);
	      
	      gNew->SetPoint(iPoint,xValue[iPoint],yValue[iPoint]);
	      gNew->SetPointError(iPoint,xError,yError);
	      
      
    }
  
  return gNew;
  
}
TGraphErrors *GetInterpolSpectrum2D(TGraphErrors *g1, TGraphErrors *g2, Double_t d1, Double_t d2,Double_t dSqrts,TString method)
{
         
         if(!g1) return 0x0;
         if(!g2) return 0x0;

         TGraphErrors  *gInterpol     = new TGraphErrors(g1->GetN());
	 TGraphErrors  *gAlpha        = new TGraphErrors(g1->GetN());
	 TGraphErrors** gPtvsSqrts    = new TGraphErrors*[g1->GetN()];
	 TGraphErrors** gPtvsEnergies = new TGraphErrors*[g1->GetN()];
	 TF1**          fPowerlawFits = new TF1*[g1->GetN()];
	 
	 

         for(Int_t i = 0; i < g1->GetN(); i++)
         {

                 TGraphErrors *grint = new TGraphErrors(1);
                 grint->SetPoint(0, dSqrts, 0);
                 TGraphErrors *gToFit = new TGraphErrors(2);
		 
		 
                 gToFit->SetPoint(0, d1, g1->GetY()[i]);
                 gToFit->SetPointError(0, 0, g1->GetEY()[i]);
		 
		
                 gToFit->SetPoint(1, d2, g2->GetY()[i]);
                 gToFit->SetPointError(1, 0, g2->GetEY()[i]);
		 
		 
		 
		 
                 TF1 *fPowerlaw = new TF1("fPowerlaw","[0]*x^([1])", 0,10000);
                 fPowerlaw->SetParameters(0, 0.1);
                 fPowerlaw->SetParameters(1, 2.0);

	         for(Int_t l = 0; l < 10; l++) gToFit->Fit(fPowerlaw,"Q");

		
                 (TVirtualFitter::GetFitter())->GetConfidenceIntervals(grint, 0.68);
		 
		 Double_t alpha  = fPowerlaw->GetParameter(1);
		 Double_t alphaE = fPowerlaw->GetParError(1);
		 
                 gInterpol->SetPoint(i, g1->GetX()[i],fPowerlaw->Eval(dSqrts));
		 gInterpol->SetPointError(i, 0, grint->GetEY()[0]);
		 
		 gAlpha->SetPoint(i, g1->GetX()[i],alpha);
		 //gAlpha->SetPointError(i, 0,alphaE);
		 
		 gPtvsSqrts[i]= new TGraphErrors(1);
		 gPtvsSqrts[i]->SetPoint(0,dSqrts,gInterpol->GetY()[i]);
		 //gPtvsSqrts[i]->Sort();
		 
		 fPowerlawFits[i] = fPowerlaw;
		 
		 gPtvsEnergies[i] = gToFit;
		 
		 
		 
                 delete grint;
                 //delete fPowerlaw;
		 
		 
         }

         if ( method.CompareTo("PCM") == 0 ){
	   
		  graphAlphaPCM = gAlpha;
		  graphPtvsSqrtsPCM = new TGraphErrors*[g1->GetN()];
		  fPowerlawPCM	    = new TF1*[g1->GetN()];
		  gPtvsEnergiesPCM  = new TGraphErrors*[g1->GetN()];
		  
		  for ( Int_t i = 0; i < g1->GetN(); i++ ){
		    
			graphPtvsSqrtsPCM[i] = gPtvsSqrts[i];
			fPowerlawPCM[i] = fPowerlawFits[i];
			gPtvsEnergiesPCM[i] = gPtvsEnergies[i];
			
		  }
	   
	 } else {
	   
		  graphAlphaDalitz = gAlpha;
		  graphPtvsSqrtsDalitz = new TGraphErrors*[g1->GetN()];
		  fPowerlawDalitz      = new TF1*[g1->GetN()];
		  gPtvsEnergiesDalitz  = new TGraphErrors*[g1->GetN()];
		  
		  for ( Int_t i = 0; i < g1->GetN(); i++ ){
		    
			graphPtvsSqrtsDalitz[i] = gPtvsSqrts[i];
			fPowerlawDalitz[i] = fPowerlawFits[i];
			gPtvsEnergiesDalitz[i] = gPtvsEnergies[i];
			
		  }
	 }
         
         
         
         return gInterpol;
}

TGraphErrors *GetInterpolSpectrum3D(TGraphErrors *g1, TGraphErrors *g2,TGraphErrors *g3, Double_t d1, Double_t d2, Double_t d3, Double_t dSqrts= 10000)
{
         if(!g1) return 0x0;
         if(!g2) return 0x0;
         if(!g3) return 0x0;

         TGraphErrors *gInterpol = new TGraphErrors(g1->GetN());

         for(Int_t i = 0; i < g1->GetN(); i++)
         {
                 TGraphErrors *grint = new TGraphErrors(1);
                 grint->SetPoint(0, dSqrts, 0);
                 TGraphErrors *gToFit = new TGraphErrors(3);
                 gToFit->SetPoint(0, d1, g1->GetY()[i]);
                 gToFit->SetPointError(0, 0, g1->GetEY()[i]);
                 gToFit->SetPoint(1, d2, g2->GetY()[i]);
                 gToFit->SetPointError(1, 0, g2->GetEY()[i]);
                 gToFit->SetPoint(2, d3, g3->GetY()[i]);
                 gToFit->SetPointError(2, 0, g3->GetEY()[i]);

                 TF1 *fPowerlaw = new TF1("fPowerlaw","[0]*x^([1])", 0,10000);
                 fPowerlaw->SetParameters(0, 0.1);
                 fPowerlaw->SetParameters(1, 2.0);

                 for(Int_t l = 0; l < 10; l++) gToFit->Fit(fPowerlaw,"0Q");

                    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(grint, 0.68);
                    gInterpol->SetPoint(i, g1->GetX()[i],
                    fPowerlaw->Eval(dSqrts));
                    gInterpol->SetPointError(i, 0, grint->GetEY()[0]);

                 delete grint;
                 delete fPowerlaw;
         }

         return gInterpol;
}


TGraphAsymmErrors* CancelOutMaterialError(TGraphAsymmErrors* graphYieldPi0, TString method){
  
  
     Double_t* valueX     = graphYieldPi0->GetX();
     Double_t* valueY     = graphYieldPi0->GetY();
     Double_t* errorXlow  = graphYieldPi0->GetEXlow();
     Double_t* errorXhigh = graphYieldPi0->GetEXhigh();
     Double_t* errorYlow  = graphYieldPi0->GetEYlow();
     Double_t* errorYhigh = graphYieldPi0->GetEYhigh();

     Int_t nPoints = graphYieldPi0->GetN();
  
     
          
      if( method.CompareTo("PCMPCM") == 0 ){  //For RpPb PCM/PCM cancelout the error of the material
	
	
	
	    Double_t  errorMaterial = 4.5;
       
	 
	    for(Int_t i = 0; i < nPoints; i++){
	  

                Double_t materialTwoE   =  ( 2*errorMaterial * valueY[i] ) / 100.0;
		
		//cout<<"Before "<<i<<" "<<valueY[i]<<" "<<errorYlow[i]<<" errorYLow/valueY "<<errorYlow[i]/valueY[i]<<endl;
                errorYlow[i]  = TMath::Sqrt(   ( errorYlow[i] * errorYlow[i]  )  - ( materialTwoE * materialTwoE ) );
		//cout<<"After  "<<i<<" "<<valueY[i]<<" "<<errorYlow[i]<<" errorYLow/valueY "<<errorYlow[i]/valueY[i]<<endl;
              
		//cout<<"Before "<<i<<" "<<valueY[i]<<" "<<errorYhigh[i]<<" errorYhigh/valueY "<<errorYhigh[i]/valueY[i]<<endl;
               	errorYhigh[i] = TMath::Sqrt(   ( errorYhigh[i]* errorYhigh[i] )  - ( materialTwoE * materialTwoE ) );
		//cout<<"After  "<<i<<" "<<valueY[i]<<" "<<errorYhigh[i]<<" errorYhigh/valueY "<<errorYhigh[i]/valueY[i]<<endl;
               
		

	    }
        
	
      } else if( method.CompareTo("PCMDalitz") == 0) {
	
	    

	    Double_t  errorMaterial = 4.5;
       
	 
	    for(Int_t i = 0; i < nPoints; i++){
	  

                Double_t materialE 	=  (  errorMaterial   * valueY[i] ) / 100.0;
		Double_t materialTwoE   =  (  2*errorMaterial * valueY[i] ) / 100.0;

                errorYlow[i]  = TMath::Sqrt(   ( errorYlow[i] * errorYlow[i]  )  - ( materialTwoE * materialTwoE ) + ( materialE*materialE ) );
                errorYhigh[i] = TMath::Sqrt(   ( errorYhigh[i]* errorYhigh[i] )  - ( materialTwoE * materialTwoE ) + ( materialE*materialE ) );

	    }
        
	   
	
	
      } else if ( method.CompareTo("Dalitz") == 0 ) {
	
	cout<<"Entro Dalitz"<<endl;
	
	
	    Double_t  errorMaterial = 4.5;
       


	    for(Int_t i = 0; i < nPoints; i++){
	  

                Double_t materialE 	=  (  errorMaterial   * valueY[i] ) / 100.0;
		
		//cout<<valueX[i]<<" "<<"errorLow:  "<<errorYlow[i] * errorYlow[i]<<"  material:    "<<materialE*materialE<<endl;
	        //cout<<valueX[i]<<" "<<"errorHigh: "<<errorYhigh[i] * errorYhigh[i]<<"  material:  "<<materialE*materialE<<endl;

                errorYlow[i]  = TMath::Sqrt(   ( errorYlow[i] * errorYlow[i]  )  - ( materialE*materialE ) );
                errorYhigh[i] = TMath::Sqrt(   ( errorYhigh[i]* errorYhigh[i] )  - ( materialE*materialE ) );

	    }
	
	
      }
      
      //cout<<"LLego al fin"<<endl;
      
      TGraphAsymmErrors* graphInvYieldPi0CancelOutMaterial = new TGraphAsymmErrors(nPoints,valueX,valueY,errorXlow,errorXhigh,errorYlow,errorYhigh);
      
      
      //cout<<"Paso la inicializacion"<<endl;
	
      
      return graphInvYieldPi0CancelOutMaterial;
}
