/****************************************************************************************************************************
******  provided by Gamma Conversion Group, PWGGA,                                                                      *****
******    Lucia Leardini, leardini@cern.ch                                                                              *****
*****************************************************************************************************************************/

Double_t xSection2760GeVpp =        55.416*1e-3;
Double_t xSection2760GeVErrpp =     3.9;
Double_t xSection2760GeVppINEL = 62.8*1e9;


TGraph* ScaleGraph (TGraph* graph, Double_t scaleFac){
    TGraph* dummyGraph      = (TGraph*)graph->Clone(Form("%s_Scaled",graph->GetName()));
    Double_t * xValue       = dummyGraph->GetX();
    Double_t * yValue       = dummyGraph->GetY();
    Int_t nPoints           = dummyGraph->GetN();
    for (Int_t i = 0; i < nPoints; i++){
        yValue[i]           = yValue[i]*scaleFac;
    }
    TGraph* returnGraph     = new TGraph(nPoints,xValue,yValue);
    return returnGraph;
}

TGraphAsymmErrors* ScaleGraph (TGraphAsymmErrors* graph, Double_t scaleFac){
   TGraphAsymmErrors* dummyGraph    = (TGraphAsymmErrors*)graph->Clone(Form("%s_Scaled",graph->GetName()));

    Double_t* xValue                = dummyGraph->GetX();
    Double_t* yValue                = dummyGraph->GetY();
    Double_t* xErrorLow             = dummyGraph->GetEXlow();
    Double_t* xErrorHigh            = dummyGraph->GetEXhigh();
    Double_t* yErrorLow             = dummyGraph->GetEYlow();
    Double_t* yErrorHigh            = dummyGraph->GetEYhigh();
    Int_t nPoints                   = dummyGraph->GetN();
    for (Int_t i = 0; i < nPoints; i++){
        yValue[i]                   = yValue[i]*scaleFac;
        yErrorLow[i]                = yErrorLow[i]*scaleFac;
        yErrorHigh[i]               = yErrorHigh[i]*scaleFac;
    }
    TGraphAsymmErrors* returnGraph  = new TGraphAsymmErrors(nPoints,xValue,yValue,xErrorLow,xErrorHigh,yErrorLow,yErrorHigh);
    return returnGraph;
}

TGraphErrors* ScaleGraph (TGraphErrors* graph, Double_t scaleFac){
    TGraphErrors* dummyGraph    = (TGraphErrors*)graph->Clone(Form("%s_Scaled",graph->GetName()));
    Double_t* xValue            = dummyGraph->GetX();
    Double_t* yValue            = dummyGraph->GetY();
    Double_t* xError            = dummyGraph->GetEX();
    Double_t* yError            = dummyGraph->GetEY();
    Int_t nPoints               = dummyGraph->GetN();
    for (Int_t i = 0; i < nPoints; i++){
        yValue[i]               = yValue[i]*scaleFac;
        yError[i]               = yError[i]*scaleFac;
    }
    TGraphErrors* returnGraph   = new TGraphErrors(nPoints,xValue,yValue,xError,yError);

    return returnGraph;
}



TH1D *GraphAsymErrorsToHist(TGraphAsymmErrors *graph,Int_t maxPt = 50, TString name = ""){
    Double_t* xValue        =  graph->GetX();
    Double_t* yValue        = graph->GetY();
    Double_t* Exhigh        = graph->GetEXhigh();
    Double_t* Exlow         = graph->GetEXlow();
    Int_t nPoints           = graph->GetN();
    Int_t maxPoints         = 0;

    for(Int_t i = 0; i<nPoints; i++){
        if(xValue[i]<=maxPt) maxPoints++;
    }

    Double_t *newBinningX   = new Double_t[maxPoints];
    for(Int_t i = 0;i<maxPoints;i++)
        newBinningX[i]      = xValue[i]-Exlow[i];

    TH1D *hist              = new TH1D(name,"",maxPoints-1,newBinningX);

    for(Int_t i = 1;i<maxPoints;i++) hist->SetBinContent(i,yValue[i-1]);

    return hist;
}


TF1* DivideTF1(TF1* f1, TF1* f2, TString name) {

        if (!f1 || !f2) return NULL;

        Double_t xmin, xmax;
        f1->GetRange(xmin, xmax);
        Int_t nPar1                         = f1->GetNpar();
        Int_t nPar2                         = f2->GetNpar();
        TString formula1                    = f1->GetExpFormula();
        TString formula2                    = f2->GetExpFormula();

        for (Int_t i = 0; i< nPar2; i++){
            formula2.ReplaceAll(Form("[%d]",i), Form("[%d]",i+nPar1));
        }

        TF1* result = new TF1(name.Data(),Form("(%s)/(%s)",formula1.Data(), formula2.Data()), xmin, xmax);
        for (Int_t i = 0; i < nPar1; i++ ){
            result->SetParameter(i, f1->GetParameter(i));
        }
        for (Int_t j = 0; j < nPar2; j++ ){
            result->SetParameter(nPar1+j, f2->GetParameter(j));
        }

        return result;
}


TF1* MtScaledParam(TF1* param, Int_t particlePDG, Int_t particleBasePDG, Double_t scaleFactor, Bool_t isInvYield = kTRUE, Bool_t doAdditionalScaling = kFALSE) {

        if (!param || particlePDG==0 || particleBasePDG==0 || !scaleFactor || scaleFactor<0) return NULL;

        Double_t mass                   = TDatabasePDG::Instance()->GetParticle(particlePDG)->Mass();
        Double_t massBase               = TDatabasePDG::Instance()->GetParticle(particleBasePDG)->Mass();

        if (!mass || !massBase)
            return NULL;

        Double_t xMin, xMax;
        param->GetRange(xMin, xMax);
        TString paramPi0Formula         = param->GetExpFormula();
        //cout << "input parametrization : " << paramPi0Formula.Data() << endl;

        // check for cut off when m(particleBasePDG) > m(particlePDG)
        if ( (xMin*xMin + mass*mass - massBase*massBase) < 0 ) xMin = TMath::Sqrt(xMin*xMin + massBase*massBase - mass*mass);

        TString mT                      = Form("TMath::Sqrt(x*x + %f * %f - %f * %f)",mass,mass,massBase,massBase);
        TString pTovermT                = Form("x/TMath::Sqrt(x*x + %f * %f - %f * %f)",mass,mass,massBase,massBase);
        TString mTScaledFormula         = paramPi0Formula.ReplaceAll("exp", "placeholder");
        TString dummyFormula            = mTScaledFormula.ReplaceAll("x",mT.Data() );
        mTScaledFormula                 = dummyFormula.ReplaceAll("placeholder","exp");
        //cout << "output parametrization in mT: " << mTScaledFormula.Data() << endl;

        Double_t paramEvaluated         = param->Eval(5.)/param->Eval(TMath::Sqrt(25. + mass*mass - massBase*massBase));
        if (doAdditionalScaling)
            scaleFactor                 = scaleFactor * paramEvaluated;

        TString         outputFormula   = "";
        if (isInvYield) outputFormula   = Form("%f * (%s)",scaleFactor,mTScaledFormula.Data());
        else            outputFormula   = Form("%f * (%s) * (%s)",scaleFactor,pTovermT.Data(),mTScaledFormula.Data());
        //cout << "output formula : " << outputFormula.Data() << endl;

        TF1* scaledParam                = new TF1("scaledParam",outputFormula.Data(),xMin, xMax);
        Int_t nPar = param->GetNpar();
        for (Int_t i = 0; i< nPar; i++){
            scaledParam->SetParameter(i,param->GetParameter(i));
        }

        return scaledParam;
}

TF1* MtScaledParam(TF1* param, Int_t particlePDG, Double_t scaleFactor, Bool_t isInvYield = kTRUE, Bool_t doAdditionalScaling = kFALSE) {

        // wrapper for direct use with pi0 as a basis for the scaling (implemented to prevent break due to use of old function status before 22.03.2017)

        return MtScaledParam(param, particlePDG, 111, scaleFactor, isInvYield, doAdditionalScaling);
}



//WITHOUT material budget
TGraphAsymmErrors* graphCombInvYieldStatPbPb2760GeVA_0010   = NULL;
TGraphAsymmErrors* graphCombInvYieldSysPbPb2760GeVA_0010    = NULL;
TGraphAsymmErrors* graphCombInvYieldTotPbPb2760GeVA_0010    = NULL;
TGraphAsymmErrors* graphCombInvYieldStatPbPb2760GeVA_2050   = NULL;
TGraphAsymmErrors* graphCombInvYieldSysPbPb2760GeVA_2050    = NULL;
TGraphAsymmErrors* graphCombInvYieldTotPbPb2760GeVA_2050    = NULL;
//WITH material budget
TGraphAsymmErrors* graphCombInvYieldStatPbPb2760GeV_0010    = NULL;
TGraphAsymmErrors* graphCombInvYieldSysPbPb2760GeV_0010     = NULL;
TGraphAsymmErrors* graphCombInvYieldTotPbPb2760GeV_0010     = NULL;
TGraphAsymmErrors* graphCombInvYieldStatPbPb2760GeV_2050    = NULL;
TGraphAsymmErrors* graphCombInvYieldSysPbPb2760GeV_2050     = NULL;
TGraphAsymmErrors* graphCombInvYieldTotPbPb2760GeV_2050     = NULL;

//unshifted graphs
TGraphAsymmErrors* graphPCMInvYieldStatPbPb2760GeVUnShifted_0010;
TGraphAsymmErrors* graphPCMInvYieldSysPbPb2760GeVUnShifted_0010;
TGraphAsymmErrors* graphPCMInvYieldSysPbPb2760GeVforRAAUnShifted_0010;
TGraphAsymmErrors* graphPHOSInvYieldStatPbPb2760GeVUnshifted_0010;
TGraphAsymmErrors* graphPHOSInvYieldSysPbPb2760GeVUnshifted_0010;
TGraphAsymmErrors* graphPHOSInvYieldSysPbPb2760GeVforRAAUnshifted_0010;
TGraphAsymmErrors* graphEMCalInvYieldStatPbPb2760GeVUnshifted_0010;
TGraphAsymmErrors* graphEMCalInvYieldSysPbPb2760GeVUnshifted_0010;
TGraphAsymmErrors* graphEMCalInvYieldSysPbPb2760GeVforRAAUnshifted_0010;

TGraphAsymmErrors* graphPCMInvYieldStatPbPb2760GeVUnShifted_2050;
TGraphAsymmErrors* graphPCMInvYieldSysPbPb2760GeVUnShifted_2050;
TGraphAsymmErrors* graphPCMInvYieldSysPbPb2760GeVforRAAUnShifted_2050;
TGraphAsymmErrors* graphEMCalInvYieldStatPbPb2760GeVUnshifted_2050;
TGraphAsymmErrors* graphEMCalInvYieldSysPbPb2760GeVUnshifted_2050;
TGraphAsymmErrors* graphEMCalInvYieldSysPbPb2760GeVforRAAUnshifted_2050;

TGraphAsymmErrors* graphRatioPCMCombFitStat2760GeV_0010;
TGraphAsymmErrors* graphRatioPCMCombFitSys2760GeV_0010;
TGraphAsymmErrors* graphRatioPHOSCombFitStat2760GeV_0010;
TGraphAsymmErrors* graphRatioPHOSCombFitSys2760GeV_0010;
TGraphAsymmErrors* graphRatioEMCalCombFitStat2760GeV_0010;
TGraphAsymmErrors* graphRatioEMCalCombFitSys2760GeV_0010;
TGraphAsymmErrors* graphRatioPCMCombFitStat2760GeV_2050;
TGraphAsymmErrors* graphRatioPCMCombFitSys2760GeV_2050;
TGraphAsymmErrors* graphRatioEMCalCombFitStat2760GeV_2050;
TGraphAsymmErrors* graphRatioEMCalCombFitSys2760GeV_2050;

// TGraphAsymmErrors*  graphCombInvYieldTotPbPb2760GeVYShifted_0010;
TGraphAsymmErrors*  graphCombInvYieldSysPbPb2760GeVYShifted_0010;
TGraphAsymmErrors*  graphCombInvYieldStatPbPb2760GeVYShifted_0010;
TGraphAsymmErrors*  graphPCMInvYieldStatPbPb2760GeVYShifted_0010;
TGraphAsymmErrors*  graphPCMInvYieldSysPbPb2760GeVYShifted_0010;
TGraphAsymmErrors*  graphPHOSInvYieldStatPbPb2760GeVYShifted_0010;
TGraphAsymmErrors*  graphPHOSInvYieldSysPbPb2760GeVYShifted_0010;
TGraphAsymmErrors*  graphEMCalInvYieldStatPbPb2760GeVYShifted_0010;
TGraphAsymmErrors*  graphEMCalInvYieldSysPbPb2760GeVYShifted_0010;

// TGraphAsymmErrors*  graphCombInvYieldTotPbPb2760GeVYShifted_2050;
TGraphAsymmErrors*  graphCombInvYieldSysPbPb2760GeVYShifted_2050;
TGraphAsymmErrors*  graphCombInvYieldStatPbPb2760GeVYShifted_2050;
TGraphAsymmErrors*  graphPCMInvYieldStatPbPb2760GeVYShifted_2050;
TGraphAsymmErrors*  graphPCMInvYieldSysPbPb2760GeVYShifted_2050;
TGraphAsymmErrors*  graphEMCalInvYieldStatPbPb2760GeVYShifted_2050;
TGraphAsymmErrors*  graphEMCalInvYieldSysPbPb2760GeVYShifted_2050;

TF1* fitBylinkinPbPb2760GeVPtLHC11h_0010;
TF1* fitBylinkinPbPb2760GeVPtLHC11h_2050;
TF1* fitQCDPbPb2760GeVPtLHC11h_0010;
TF1* fitQCDPbPb2760GeVPtLHC11h_2050;
TF1* fitTsallisPbPb2760GeVPtLHC11h_0010;
TF1* fitTsallisPbPb2760GeVPtLHC11h_2050;

TF1 *fitNormBylinkinPbPb2760GeVPtLHC11h_0010;
TF1 *fitNormBylinkinPbPb2760GeVPtLHC11h_2050;

// TF1* fitBylinkinPbPb2760GeVPtLHC11hYshift_0010;
// TF1* fitBylinkinPbPb2760GeVPtLHC11hYshift_2050;

TF1* fitQCDInvYield2760GeVLHC11h_0010;
TF1* fitQCDInvYield2760GeVLHC11h_2050;
TF1* fitQCDInvYieldPbPb2760GeV_0010;
TF1* fitQCDInvYieldPbPb2760GeV_2050;

TF1 *fitLowPtBylinkin_0010;
TF1 *fitHighPtBylinkin_0010;
TF1 *fitLowPtBylinkin_2050;
TF1 *fitHighPtBylinkin_2050;

// RAA directly from measurement (non calculated)
TGraphAsymmErrors*  graphPCMPi0RAAStatPbPb2760GeVYShifted_0010;
TGraphAsymmErrors*  graphPCMPi0RAASysPbPb2760GeVYShifted_0010;
TGraphAsymmErrors*  graphEMCalPi0RAAStatPbPb2760GeVYShifted_0010;
TGraphAsymmErrors*  graphEMCalPi0RAASysPbPb2760GeVYShifted_0010;
TGraphAsymmErrors*  graphPCMPi0RAAStatPbPb2760GeVYShifted_2050;
TGraphAsymmErrors*  graphPCMPi0RAASysPbPb2760GeVYShifted_2050;
TGraphAsymmErrors*  graphEMCalPi0RAAStatPbPb2760GeVYShifted_2050;
TGraphAsymmErrors*  graphEMCalPi0RAASysPbPb2760GeVYShifted_2050;
TGraphAsymmErrors*  graphPCMEtaRAAStatPbPb2760GeVYShifted_0010;
TGraphAsymmErrors*  graphPCMEtaRAASysPbPb2760GeVYShifted_0010;
TGraphAsymmErrors*  graphEMCalEtaRAAStatPbPb2760GeVYShifted_0010;
TGraphAsymmErrors*  graphEMCalEtaRAASysPbPb2760GeVYShifted_0010;
TGraphAsymmErrors*  graphPCMEtaRAAStatPbPb2760GeVYShifted_2050;
TGraphAsymmErrors*  graphPCMEtaRAASysPbPb2760GeVYShifted_2050;
TGraphAsymmErrors*  graphEMCalEtaRAAStatPbPb2760GeVYShifted_2050;
TGraphAsymmErrors*  graphEMCalEtaRAASysPbPb2760GeVYShifted_2050;

TGraphAsymmErrors *graphPHOSPi0InvYieldStatPbPb2760GeV_0010;
TGraphAsymmErrors *graphPHOSPi0InvYieldSysPbPb2760GeV_0010;

//Raa calculated for each meas with comb fit
TGraphAsymmErrors* graphRAAPCM0010;
TGraphAsymmErrors* graphRAASysPCM0010;
TGraphAsymmErrors* graphRAAPHOS0010;
TGraphAsymmErrors* graphRAASysPHOS0010;
TGraphAsymmErrors* graphRAAEMCal0010;
TGraphAsymmErrors* graphRAASysEMCal0010;
TGraphAsymmErrors* graphRAAPCM2050;
TGraphAsymmErrors* graphRAASysPCM2050;
TGraphAsymmErrors* graphRAAEMCal2050;
TGraphAsymmErrors* graphRAASysEMCal2050;
//Raa stat error histos for meas Raa combination
TH1D* histoRAAStatPCM0010;
TH1D* histoRAAStatPCM2050;
TH1D* histoRAAStatEMCal0010;
TH1D* histoRAAStatEMCal2050;
TH1D* histoRAAStatPHOS0010;
// Raa graphs from calcraa with comb fit (see above)
TGraphAsymmErrors* graphCombRAATotPbPb2760GeV_0010 = NULL;
TGraphAsymmErrors* graphCombRAAStatPbPb2760GeV_0010 = NULL;
TGraphAsymmErrors* graphCombRAASysPbPb2760GeV_0010 = NULL;
TGraphAsymmErrors* graphCombRAATotPbPb2760GeV_2050 = NULL;
TGraphAsymmErrors* graphCombRAAStatPbPb2760GeV_2050 = NULL;
TGraphAsymmErrors* graphCombRAASysPbPb2760GeV_2050 = NULL;


//loading file for both meson plots
TGraphAsymmErrors *graphRAAPi0StatBothMeson_0010;
TGraphAsymmErrors *graphRAAPi0SysBothMeson_0010;
TGraphAsymmErrors *graphRAAPi0StatBothMeson_2050;
TGraphAsymmErrors *graphRAAPi0SysBothMeson_2050;
TGraphAsymmErrors *graphRAAEtaStatBothMeson_0010;
TGraphAsymmErrors *graphRAAEtaSysBothMeson_0010;
TGraphAsymmErrors *graphRAAEtaStatBothMeson_2050;
TGraphAsymmErrors *graphRAAEtaSysBothMeson_2050;


TGraphAsymmErrors* graphCombEtatoPi0StatPbPb2760GeV_0010 = NULL;
TGraphAsymmErrors* graphCombEtatoPi0SysPbPb2760GeV_0010 = NULL;
TGraphAsymmErrors* graphCombEtatoPi0TotPbPb2760GeV_0010 = NULL;
TGraphAsymmErrors* graphCombEtatoPi0StatPbPb2760GeV_2050 = NULL;
TGraphAsymmErrors* graphCombEtatoPi0SysPbPb2760GeV_2050 = NULL;
TGraphAsymmErrors* graphCombEtatoPi0TotPbPb2760GeV_2050 = NULL;

TFile*  filePHOSPbPb;
TDirectory* directoryPHOSPi0PbPb0010;
TH1D*   histoPi0PHOSPbPb0010;
TH1D*   histoPi0PHOSSysPbPb0010;
TGraphAsymmErrors* graphPHOSYieldPi0SysErrPbPb0010;
TH1D*   histoPi0PHOSSysRAAPbPb0010;
TGraphAsymmErrors* graphSysErrRAAYieldPi0PHOSPbPb0010;
TGraphAsymmErrors* graphPHOSYieldPi0SysErrPbPb0010Red;
TGraphAsymmErrors* graphSysErrRAAYieldPi0PHOSPbPb0010Red;
TH1D*   histoPi0PHOSPbPb0010Red;
TH1D*   histoPHOSMassData0010;
TH1D*   histoPHOSMassMC0010;
TH1D*   histoPHOSWidthData0010;
TH1D*   histoPHOSWidthMC0010;



TFile*  fileFinalResultsPP;
// PP combined yields (x-shifted)
TGraphAsymmErrors* graphInvSectionCombStatPi02760GeV;
TGraphAsymmErrors* graphInvSectionCombStatPi02760GeVPlot;
TGraphAsymmErrors* graphInvSectionCombSysPi02760GeV;
TGraphAsymmErrors* graphInvSectionCombSysPi02760GeVPlot;
TGraphAsymmErrors* graphInvSectionPCMStatPi02760GeV;
TGraphAsymmErrors* graphInvSectionPCMSysPi02760GeV;
TGraphAsymmErrors* graphInvSectionPHOSStatPi02760GeV;
TGraphAsymmErrors* graphInvSectionPHOSSysPi02760GeV;
TGraphAsymmErrors* graphInvSectionEMCalStatPi02760GeV;
TGraphAsymmErrors* graphInvSectionEMCalSysPi02760GeV;
TGraphAsymmErrors* graphInvSectionCombStatEta2760GeV;
TGraphAsymmErrors* graphInvSectionCombStatEta2760GeVPlot;
TGraphAsymmErrors* graphInvSectionCombSysEta2760GeV;
TGraphAsymmErrors* graphInvSectionCombSysEta2760GeVPlot;
TGraphAsymmErrors* graphInvSectionPCMStatEta2760GeV;
TGraphAsymmErrors* graphInvSectionPCMSysEta2760GeV;
TGraphAsymmErrors* graphInvSectionEMCalStatEta2760GeV;
TGraphAsymmErrors* graphInvSectionEMCalSysEta2760GeV;

TGraphAsymmErrors* graphRatioEtaToPi0Comb2760GeVStatErr;
TGraphAsymmErrors* graphRatioEtaToPi0Comb2760GeVSysErr;

// PP combined yields (y-shifted)
TGraphAsymmErrors*  graphInvSectionCombStatPi02760GeVforRAA;
TGraphAsymmErrors*  graphInvSectionCombSysPi02760GeVforRAA;
TGraphAsymmErrors*  graphInvSectionPCMStatPi02760GeVforRAA;
TGraphAsymmErrors*  graphInvSectionPCMSysPi02760GeV_yShifted;
TGraphAsymmErrors*  graphInvSectionPCMSysPi02760GeVforRAA;
TGraphAsymmErrors*  graphInvSectionPHOSStatPi02760GeVforRAA;
TGraphAsymmErrors*  graphInvSectionPHOSSysPi02760GeV_yShifted;
TGraphAsymmErrors*  graphInvSectionPHOSSysPi02760GeVforRAA;
TGraphAsymmErrors*  graphInvSectionEMCalStatPi02760GeVforRAA;
TGraphAsymmErrors*  graphInvSectionEMCalSysPi02760GeV_yShifted;
TGraphAsymmErrors*  graphInvSectionEMCalSysPi02760GeVforRAA;
TGraphAsymmErrors*  graphInvSectionCombStatEta2760GeVforRAA;
TGraphAsymmErrors*  graphInvSectionCombSysEta2760GeVforRAA;
TGraphAsymmErrors*  graphInvSectionPCMStatEta2760GeVforRAA;
TGraphAsymmErrors*  graphInvSectionPCMSysEta2760GeV_yShifted;
TGraphAsymmErrors*  graphInvSectionPCMSysEta2760GeVforRAA;
TGraphAsymmErrors*  graphInvSectionEMCalStatEta2760GeVforRAA;
TGraphAsymmErrors*  graphInvSectionEMCalSysEta2760GeV_yShifted;
TGraphAsymmErrors*  graphInvSectionEMCalSysEta2760GeVforRAA;



TGraph *graphPi0JetQuenching18;
TGraph *graphPi0JetQuenching22;
TGraph *graphPi0JetQuenching26;

TGraph *graphEtaJetQuenching18;
TGraph *graphEtaJetQuenching22;
TGraph *graphEtaJetQuenching26;

TGraph *graphEtatoPi0RatioJetQuenching18;
TGraph *graphEtatoPi0RatioJetQuenching22;
TGraph *graphEtatoPi0RatioJetQuenching26;

TGraphErrors* Xiao_Raa_0020;
TGraphErrors* Xiao_Raa_0005;
TGraphErrors* Xiao_Raa_0510;
TGraphErrors* Xiao_Raa_1020;
TGraphErrors*   Xiao_Raa_2040;
TGraphErrors* Xiao_Raa_4060;
TGraphErrors* Xiao_Raa_6080;
TGraphAsymmErrors*  gWHDG_Raa_0020;
TGraphAsymmErrors*  gWHDG_Raa_0005;
TGraphAsymmErrors*  gWHDG_Raa_0510;
TGraphAsymmErrors*  gWHDG_Raa_0010;
TGraphAsymmErrors*  gWHDG_Raa_1020;
TGraphAsymmErrors*  gWHDG_Raa_2040;
TGraphAsymmErrors*  gWHDG_Raa_2050;
TGraphAsymmErrors*  gWHDG_Raa_4060;
TGraphAsymmErrors*  gWHDG_Raa_6080;
TGraphAsymmErrors*  gWHDG_Eta_Raa_0010;
TGraphAsymmErrors*  gWHDG_Eta_Raa_2050;

TGraph* gEPOS_Spec_0005;
TGraph* gEPOS_Spec_0510;
TGraph* gEPOS_Spec_1020;
TGraph* gEPOS_Spec_0020;
TGraph* gEPOS_Spec_2040;
TGraph* gEPOS_Spec_4060;
TGraph* gEPOS_Spec_6080;
TGraph* gKopeliovichELoss_Spec_0005;
TGraph* gKopeliovichELoss_Spec_0510;
TGraph* gKopeliovichELoss_Spec_1020;
TGraph* gKopeliovichELoss_Spec_0020;
TGraph* gKopeliovichELoss_Spec_2040;
TGraph* gKopeliovichELoss_Spec_4060;
TGraph* gKopeliovichELoss_Spec_6080;
TGraph* gKopeliovichTotal_Spec_0005;
TGraph* gKopeliovichTotal_Spec_0510;
TGraph* gKopeliovichTotal_Spec_1020;
TGraph* gKopeliovichTotal_Spec_0020;
TGraph* gKopeliovichTotal_Spec_2040;
TGraph* gKopeliovichTotal_Spec_4060;
TGraph* gKopeliovichTotal_Spec_6080;

TGraph* gKopeliovichHydro_Spec_0005;
TGraph* gKopeliovichHydro_Spec_0510;
TGraph* gKopeliovichHydro_Spec_1020;
TGraph* gKopeliovichHydro_Spec_0020;
TGraph* gKopeliovichHydro_Spec_2040;
TGraph* gKopeliovichHydro_Spec_4060;
TGraph* gKopeliovichHydro_Spec_6080;

TGraph* gEPOS_RAA_0005;
TGraph* gEPOS_RAA_0510;
TGraph* gEPOS_RAA_1020;
TGraph* gEPOS_RAA_0020;
TGraph* gEPOS_RAA_2040;
TGraph* gEPOS_RAA_4060;
TGraph* gEPOS_RAA_6080;
TGraph* gRatioEPOSToFit0005;
TGraph* gRatioEPOSToFit0510;
TGraph* gRatioEPOSToFit1020;
TGraph* gRatioEPOSToFit2040;
TGraph* gRatioEPOSToFit4060;
TGraph* gRatioEPOSToFit6080;
TGraph* gRatioKopeliovichELossToFit0005;
TGraph* gRatioKopeliovichELossToFit0510;
TGraph* gRatioKopeliovichELossToFit1020;
TGraph* gRatioKopeliovichELossToFit2040;
TGraph* gRatioKopeliovichELossToFit4060;
TGraph* gRatioKopeliovichELossToFit6080;
TGraph* gRatioKopeliovichTotalToFit0005;
TGraph* gRatioKopeliovichTotalToFit0510;
TGraph* gRatioKopeliovichTotalToFit1020;
TGraph* gRatioKopeliovichTotalToFit2040;
TGraph* gRatioKopeliovichTotalToFit4060;
TGraph* gRatioKopeliovichTotalToFit6080;
TGraph* gRatioKopeliovichHydroToFit0005;
TGraph* gRatioKopeliovichHydroToFit0510;
TGraph* gRatioKopeliovichHydroToFit1020;
TGraph* gRatioKopeliovichHydroToFit2040;
TGraph* gRatioKopeliovichHydroToFit4060;
TGraph* gRatioKopeliovichHydroToFit6080;

TGraphErrors*   Vitev_Bas_Raa_0020;
TGraphErrors*   Vitev_Bas_Raa_0005;
TGraphErrors*   Vitev_Bas_Raa_0510;
TGraphErrors*   Vitev_Bas_Raa_1020;
TGraphErrors*   Vitev_Bas_Raa_2040;
TGraphErrors*   Vitev_Bas_Raa_4060;
TGraphErrors*   Vitev_Bas_Raa_6080;
TGraphErrors*   Vitev_ShlSel_Raa_0020;
TGraphErrors*   Vitev_ShlSel_Raa_0005;
TGraphErrors*   Vitev_ShlSel_Raa_0510;
TGraphErrors*   Vitev_ShlSel_Raa_1020;

TGraphErrors*   graphPHENIX200GeVRAA_0010;
TGraphErrors*   graphPHENIX200GeVRAA_1020;
TGraphErrors*   graphPHENIX200GeVRAA_0020;
TGraphErrors*   graphPHENIX200GeVRAA_2040;
TGraphErrors*   graphPHENIX200GeVRAA_4060;
TGraphErrors*   graphPHENIX200GeVRAA_6080;
TGraphErrors*   graphPHENIX39GeVRAA_0010;
TGraphErrors*   graphPHENIX39GeVRAA_1020;
TGraphErrors*   graphPHENIX39GeVRAA_0020;
TGraphErrors*   graphPHENIX39GeVRAA_2040;
TGraphErrors*   graphPHENIX39GeVRAA_4060;
TGraphErrors*   graphPHENIX39GeVRAA_6080;

TGraphErrors*   graphPHENIX62GeVRAA_0010;
TGraphErrors*   graphPHENIX62GeVRAA_1020;
TGraphErrors*   graphPHENIX62GeVRAA_0020;
TGraphErrors*   graphPHENIX62GeVRAA_2040;
TGraphErrors*   graphPHENIX62GeVRAA_4060;
TGraphErrors*   graphPHENIX62GeVRAA_6080;

TGraphErrors*   graphPHENIX200GeVEtaRAA_0010;
TGraphErrors*   graphPHENIX200GeVEtaRAA_2040;
TGraphErrors*   graphPHENIX200GeVEtaRAA_2060;

TGraphErrors*   graphWA98_17_3GeVRAA_0013;

TH2F * histo2DInvYieldSectionPi0LHC11h;
TLatex *labelRatioToFitPi0;
TLatex *labelRelSysErrPi0;
TLatex *labelWeightsPi0;
TGraphErrors *graphEPOSpred_0010;
TGraphErrors *graphEPOSpred_2050;
TH1D *histoEPOSpred_0010;
TH1D *histoEPOSpred_2050;
TLatex *labelFactorLower;
TLatex *labelFactorLowerOnlyPbPb;
TLatex *labelFactorUpper;
TLatex *labelSystOnlyPbPb;
TLatex *labelSyst;
TLatex *labelSystRaa;
TLatex *labelDetSysInvYieldSectionPi0LHC11h;
TLatex *labelDetSysInvYieldSectionPi0LHC11hnoPrelim;
TLatex *labelDetSysInvYieldSectionPi0LHC11hwithPP;
TLatex *labelRelStatErrPi0;
TLegend* legendXsectionPaperOnlyRatios;
TLatex *labelDetSysRatioToModelsLHC11h;
TH2F * histo2DInvYieldSectionLHC11hwithPP;

TGraphAsymmErrors* sysErrorRel2010and2011_0010[2];
TGraphAsymmErrors* sysErrorRel2010and2011_2040[2];
TGraphAsymmErrors* sysErrorRelRAA2010and2011_0010[2];
TGraphAsymmErrors* sysErrorRelRAA2010and2011_2040[2];
TGraphAsymmErrors* statErrorRel2010and2011_0010[2];
TGraphAsymmErrors* statErrorRel2010and2011_2040[2];
TGraphAsymmErrors* statErrorRelRAA2010and2011_0010[2];
TGraphAsymmErrors* statErrorRelRAA2010and2011_2040[2];
