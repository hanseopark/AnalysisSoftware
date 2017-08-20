#ifndef GAMMACONV_ConversionFunctions
#define GAMMACONV_ConversionFunctions

    #include "Math/WrappedTF1.h"
    #include "Math/BrentRootFinder.h"
    #include "TFitResultPtr.h"
    #include "TFitResult.h"
    #include "TObjString.h"
    #include "TRandom.h"
    #include "TSpline.h"
    #include "TKey.h"

    // ****************************************************************************************************************
    // *********************** declaration of functions defined in this header ****************************************
    // ****************************************************************************************************************
    Float_t             CalculateMeanPt(const TF1* );
    void                CalculateMeanPtWithError(const TF1* , Float_t& , Float_t& );
    // TH1D*            CalculateHistoRatioToFit (TH1D*, TF1*);
    // TH1F*            CalculateHistoRatioToFit (TH1F*, TF1*);
    TF1*                ScaleTF1(TF1* func, Double_t constant, TString name);
    TF1*                MultiplyTF1(TF1* f1, TF1* f2, TString name);
    TF1*                DivideTF1(TF1* f1, TF1* f2, TString name);
    TH1D*               DivideTF1IntoHisto(TF1* f1, TF1* f2, TString name, TH1D *dummy);
    TH1D*               CorrectHistoToBinCenter (TH1D*);
    TGraphErrors* 	    CalculateGraphRatioToGraph(TGraphErrors*, TGraphErrors*);
    TGraphAsymmErrors*  CalculateAsymGraphRatioToGraph(TGraphAsymmErrors* graphA, TGraphAsymmErrors* graphB);
    TGraphErrors*       CalculateGraphErrRatioToFit (TGraphErrors* , TF1* );
    TGraphAsymmErrors*  CalculateGraphErrRatioToFit (TGraphAsymmErrors* , TF1* );
    TGraph*             CalculateGraphRatioToFit (TGraph* , TF1* );
    TH1D*               CalculateHistoRatioToFitNLO (TH1D* , TF1* , Double_t );
    TGraphAsymmErrors*  CalculateSysErrFromRelSysHistoWithPtBins( TH1D* , TString , Double_t* , Double_t* , Double_t* , const Int_t  );
    TGraphAsymmErrors*  CalculateSysErrFromRelSysHisto( TH1D* , TString , Double_t* , Double_t* , Int_t , const Int_t  );
    TGraphAsymmErrors*  CalculateSysErrAFromRelSysHisto( TH1D* , TString , Double_t* , Double_t* , Int_t , const Int_t  );
    TGraphAsymmErrors*  CalculateSysErrFromRelSysHistoComplete( TH1D* , TString , Double_t* , Double_t* , Int_t , const Int_t  );
    TGraphAsymmErrors*  CalculateCombinedSysAndStatError( TGraphAsymmErrors* , TGraphAsymmErrors* );
    TGraphErrors*  	    CalculateCombinedSysAndStatError( TGraphErrors*, 	   TGraphErrors* );

    TH1D*               CalculateShiftedSpectrumInXtOrMt(TH1D* , TH1D* , TString );
    TH1D*               CalculateMassMinusExpectedMass(TH1D* , Double_t );
    void                RemoveScalingWithPt(TH1D* );
    void                RemoveScalingWithPtGraph(TGraphAsymmErrors* );
    void                ScaleWithPtGraph(TGraphAsymmErrors* );
    TGraph*             ScaleGraph (TGraph* , Double_t );
    TGraphAsymmErrors*  ScaleGraph (TGraphAsymmErrors* , Double_t );
    TGraphErrors*       ScaleGraph (TGraphErrors* , Double_t );
    Double_t*           ExtractRelErrDownAsymmGraph(TGraphAsymmErrors* );
    Double_t*           ExtractRelErrUpAsymmGraph(TGraphAsymmErrors* );
    TGraphAsymmErrors*  RebinCombPi0Graph(TGraphAsymmErrors* , TH1D* );
    TGraphErrors*       RebinNLOGraph(TGraphErrors* );
    TH1D*               CalculateFitToFitRatio(TF1 *, TF1 *);
    TH1D*               CalculateWeightedAveragePCM( TH1D* , TH1D* );
    void                ReadOutFromSysErrVector(Double_t* , Double_t* , Double_t* , Int_t , Int_t ,Int_t , Int_t ,Double_t );
    void                CalculateMeanSysErr(Double_t* , Double_t* , Double_t* , Double_t* , Int_t );
    void                CorrectSystematicErrorsWithMean(Double_t* ,Double_t* , Double_t* , Double_t* , Int_t );
    void                ProduceGraphAsymmWithoutXErrors(TGraphAsymmErrors* );
    void                ProduceGraphPartialXErrors(TGraphErrors* , Double_t );
    void                ProduceGraphAsymmPartialXErrors(TGraphAsymmErrors* , Double_t );
    void                ProduceGraphFixedXErrors(TGraphErrors* , Double_t );
    void                ProduceGraphAsymmFixedXErrors(TGraphAsymmErrors* , Double_t );
    void                ProduceGraphAsymmWithoutXYErrors(TGraphAsymmErrors* );
    void                ProduceGraphErrWithoutXErrors(TGraphErrors* );
    void                ProduceGraphErrDisplacedX(TGraphErrors* ,Double_t );
    TGraphAsymmErrors*  Add2TGraphAsymmErrorsSameBinning(TGraphAsymmErrors* ,TGraphAsymmErrors* );
    TGraphAsymmErrors*  CorrectTGraphAsymmErrorsToBinCenter(TGraphAsymmErrors* );
    TH1D*               ShortChargedHadronHisto(TH1D* );
    TH1D*               ConvertChargedHadronHisto(TH1D* , TH1D* );
    TGraphAsymmErrors*  ApplyYshiftIndividualSpectra(TGraphAsymmErrors * , TF1 *);
    Double_t            bin_shift_x(TF1 *, Double_t , Double_t , Double_t );
    Int_t               GetBinning(TObject *, Double_t* );
    TString             AutoDetectMainTList(Int_t);

    // ****************************************************************************************************************
    // ********************** definition of functions defined in this header ******************************************
    // ****************************************************************************************************************
    Float_t CalculateMeanPt(const TF1* fit){

        // fit is the invariant cross section (or proportinal to it):
        // fit ~ 1/pT dN/dpT

        TF1* fs = (TF1*) fit->Clone();
        fs->SetName("fs");

        TF1* f1 = new TF1("f1","x*fs",0.,16.);
        TF1* f2 = new TF1("f2","x^2*fs",0.,16.);

        Float_t i1 = f1->Integral(0.,16.);
        Float_t i2 = f2->Integral(0.,16.);


        delete fs;
        delete f1;
        delete f2;

        return i2/i1;
    }

    //**********************************************************************************************************
    // Calculates the mean momentum of the function given
    //**********************************************************************************************************
    void CalculateMeanPtWithError(const TF1* fSigma, Float_t& meanpt, Float_t& meanpt_error){

        meanpt = CalculateMeanPt(fSigma);

        // quadratic summ of errors
        Float_t sum_err2    = 0.;
        // loop over the parameters of the function fSigma
        for (Int_t i=0; i<fSigma->GetNpar(); i++) {
            Float_t par                 = fSigma->GetParameter(i);
            Float_t sigma               = fSigma->GetParError(i);

            // mean pt for a function with par(i) = par(i) + sigma(par(i))
            TF1* fmod_up                = (TF1*) fSigma->Clone();
            fmod_up->SetName("fmod_up");
            fmod_up->SetParameter(i, par+sigma);
            Float_t delta_meanpt_up     = fabs(meanpt - CalculateMeanPt(fmod_up));

            // mean pt for a function with par(i) = par(i) - sigma(par(i))
            TF1* fmod_down              = (TF1*) fSigma->Clone();
            fmod_up->SetName("fmod_down");
            fmod_down->SetParameter(i, par-sigma);
            Float_t delta_meanpt_down   = fabs(meanpt - CalculateMeanPt(fmod_down));

            // take the larger one of the two deltas
            Float_t delta_meanpt    = delta_meanpt_up;
            if (delta_meanpt_up < delta_meanpt_down) delta_meanpt = delta_meanpt_down;

            sum_err2                += delta_meanpt*delta_meanpt;

            delete fmod_up;
            delete fmod_down;

        }

        meanpt_error        = sqrt(sum_err2);
    }

    //**********************************************************************************************************
    // Calculates the inverse of the histo handed to the function and returns it
    //**********************************************************************************************************
    TH1D* InvertHisto (TH1D* histo){
        TH1D* histo2        = (TH1D*)histo->Clone("Dummy");
        for( Int_t I = 1; I < histo->GetNbinsX() +1 ;I++){
            Double_t yValue = 1./histo2->GetBinContent(I);
            Double_t yError = 1./histo2->GetBinContent(I)*histo2->GetBinError(I)/histo2->GetBinContent(I);
            histo2->SetBinContent(I,yValue);
            histo2->SetBinError(I,yError);
        }
        return histo2;
    }

    //**********************************************************************************************************
    // Calculates the ratio of a histogram and a spline
    //**********************************************************************************************************
    TH1D* CalculateHistoRatioToSpline (TH1D* histo, TSpline* fit){
        TH1D* histo2                = (TH1D*)histo->Clone("Dummy");
        for( Int_t I = 1; I < histo->GetNbinsX() +1 ;I++){
            Double_t xValue         = histo2->GetBinCenter(I);
            Double_t yValue         = fit->Eval(xValue);
            Double_t formerYValue   = histo2->GetBinContent(I);
            if (yValue != 0){
                histo2->SetBinContent(I,formerYValue/yValue);
                histo2->SetBinError(I,histo2->GetBinError(I)/yValue);
            }
        }
        return histo2;
    }

    //**********************************************************************************************************
    // Calculates the ratio of a histogram and a spline
    //**********************************************************************************************************
    TH1F* CalculateHistoRatioToSpline (TH1F* histo, TSpline* fit){
        TH1F* histo2                = (TH1F*)histo->Clone("Dummy");
        for( Int_t I = 1; I < histo->GetNbinsX() +1 ;I++){
            Double_t xValue         = histo2->GetBinCenter(I);
            Double_t yValue         = fit->Eval(xValue);
            Double_t formerYValue   = histo2->GetBinContent(I);
            if (yValue != 0){
                histo2->SetBinContent(I,formerYValue/yValue);
                histo2->SetBinError(I,histo2->GetBinError(I)/yValue);
            }
        }
        return histo2;
    }

    //================================================================================================================
    //Scale TF1 with constant
    //================================================================================================================
    TF1* ScaleTF1(TF1* func, Double_t constant, TString name) {

        if (!func) return NULL;

        Double_t    xMin, xMax;
        TString     formula         = func->GetExpFormula();
        func->GetRange(xMin, xMax);

        for (Int_t i=0; i<func->GetNpar(); i++) {
            formula.ReplaceAll(Form("[%d]", i), Form("[placeholder%d]",i+1));
        }
        for (Int_t i=1; i<func->GetNpar()+1; i++) {
            formula.ReplaceAll(Form("[placeholder%d]", i), Form("[%d]",i));
        }

        TF1* result                 = new TF1(name.Data(), Form("[0] * (%s)", formula.Data()), xMin, xMax);
        for (Int_t i=0; i<func->GetNpar()+1; i++) {
            if (i==0)   result->SetParameter(i, constant);
            else        result->SetParameter(i, func->GetParameter(i-1));
        }

        return result;
    }

    //================================================================================================================
    //Multiply two TF1s
    //================================================================================================================
    TF1* MultiplyTF1(TF1* f1, TF1* f2, TString name) {

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

        TF1* result = new TF1(name.Data(),Form("(%s)*(%s)",formula1.Data(), formula2.Data()), xmin, xmax);
        for (Int_t i = 0; i < nPar1; i++ ){
            result->SetParameter(i, f1->GetParameter(i));
        }
        for (Int_t j = 0; j < nPar2; j++ ){
            result->SetParameter(nPar1+j, f2->GetParameter(j));
        }

        return result;
    }

    //================================================================================================================
    //Divide two TF1s
    //================================================================================================================
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

    //================================================================================================================
    //Divide two TF1s
    //================================================================================================================
    TH1D* DivideTF1IntoHisto(TF1* f1, TF1* f2, TString name, TH1D *dummy) {

        if (!f1 || !f2) return NULL;

        TH1D *result = (TH1D*)dummy->Clone(name.Data());
        for(Int_t i=1; i<dummy->GetNbinsX()+1; i++){

        Double_t x = dummy->GetBinCenter(i);

        Double_t ratio = f1->Integral(dummy->GetXaxis()->GetBinLowEdge(i), dummy->GetXaxis()->GetBinUpEdge(i))/f2->Integral(dummy->GetXaxis()->GetBinLowEdge(i), dummy->GetXaxis()->GetBinUpEdge(i));

        result->SetBinContent(i, ratio);

        }


        return result;
    }



    //************************** Routine to calculate mt scaled params **************************************************
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

    //**********************************************************************************************************
    // Calculates spectrum in mt from pt (y-shifted) spectrum
    //**********************************************************************************************************
    TGraphAsymmErrors* CalculateSpectrumInMt(TH1D* ptSpectrum, Int_t particlePDG, TString newSpectrumName){

        if(!ptSpectrum) return NULL;
        TGraphAsymmErrors *dummy = NULL;
        TString type = ptSpectrum->ClassName();
        if(type.BeginsWith("TH1"))
            dummy = new TGraphAsymmErrors(ptSpectrum);
        TGraphAsymmErrors *mtspectrum = (TGraphAsymmErrors*)dummy->Clone();

        Int_t numberPoints     = dummy->GetN();
        Double_t mass          = TDatabasePDG::Instance()->GetParticle(particlePDG)->Mass();

        Double_t xMtSpectrum[numberPoints];
        Double_t yMtSpectrum[numberPoints];
        Double_t errorMtXlow[numberPoints];
        Double_t errorMtXhigh[numberPoints];
        Double_t errorMtYlow[numberPoints];
        Double_t errorMtYhigh[numberPoints];

        Double_t *xPtSpectrum  = dummy->GetX();
        Double_t *yPtSpectrum  = dummy->GetY();
        Double_t *errorPtXlow  = dummy->GetEXlow();
        Double_t *errorPtXhigh = dummy->GetEXhigh();
        Double_t *errorPtYlow  = dummy->GetEYlow();
        Double_t *errorPtYhigh = dummy->GetEYhigh();
        for (Int_t i = 0; i < numberPoints; i++){
            if(yPtSpectrum[i]!=0){
                xMtSpectrum[i] = TMath::Power(xPtSpectrum[i]*xPtSpectrum[i] + mass*mass,0.5);
                errorMtXlow[i] = TMath::Power(xPtSpectrum[i]*xPtSpectrum[i] + mass*mass,-0.5)*xPtSpectrum[i]*errorPtXlow[i];

                yMtSpectrum[i] = yPtSpectrum[i];

                mtspectrum->SetPoint(i,xMtSpectrum[i],yMtSpectrum[i]);
                mtspectrum->SetPointEXlow (i,errorMtXlow[i]);
                mtspectrum->SetPointEXhigh(i,errorMtXlow[i]);
                mtspectrum->SetPointEYlow (i,errorPtYlow[i]);
                mtspectrum->SetPointEYhigh(i,errorPtYhigh[i]);
            }
        }

        return mtspectrum;
    }

    TGraphAsymmErrors* CalculateSpectrumInMt(TGraphAsymmErrors* ptSpectrum, Int_t particlePDG, TString newSpectrumName){

        if(!ptSpectrum) return NULL;
        TGraphAsymmErrors *dummy      = (TGraphAsymmErrors*)ptSpectrum->Clone();
        TGraphAsymmErrors *mtspectrum = (TGraphAsymmErrors*)dummy->Clone();

        Int_t numberPoints     = dummy->GetN();
        Double_t mass          = TDatabasePDG::Instance()->GetParticle(particlePDG)->Mass();

        Double_t xMtSpectrum[numberPoints];
        Double_t yMtSpectrum[numberPoints];
        Double_t errorMtXlow[numberPoints];
        Double_t errorMtXhigh[numberPoints];
        Double_t errorMtYlow[numberPoints];
        Double_t errorMtYhigh[numberPoints];

        Double_t *xPtSpectrum  = dummy->GetX();
        Double_t *yPtSpectrum  = dummy->GetY();
        Double_t *errorPtXlow  = dummy->GetEXlow();
        Double_t *errorPtXhigh = dummy->GetEXhigh();
        Double_t *errorPtYlow  = dummy->GetEYlow();
        Double_t *errorPtYhigh = dummy->GetEYhigh();
        for (Int_t i = 0; i < numberPoints; i++){
            if(yPtSpectrum[i]!=0){
                xMtSpectrum[i] = TMath::Power(xPtSpectrum[i]*xPtSpectrum[i] + mass*mass,0.5);
                errorMtXlow[i] = TMath::Power(xPtSpectrum[i]*xPtSpectrum[i] + mass*mass,-0.5)*xPtSpectrum[i]*errorPtXlow[i];

                yMtSpectrum[i] = yPtSpectrum[i];

                mtspectrum->SetPoint(i,xMtSpectrum[i],yMtSpectrum[i]);
                mtspectrum->SetPointEXlow (i,errorMtXlow[i]);
                mtspectrum->SetPointEXhigh(i,errorMtXlow[i]);
                mtspectrum->SetPointEYlow (i,errorPtYlow[i]);
                mtspectrum->SetPointEYhigh(i,errorPtYhigh[i]);
            }
        }

        return mtspectrum;
    }

    //**********************************************************************************************************
    // Calculates the ratio of a histogram and a fit, with the posibility to integrate the function in the same
    // bin width as the data
    //**********************************************************************************************************
    TH1D* CalculateHistoRatioToFit (TH1D* histo, TF1* fit, Bool_t integrateFunction=kFALSE){
        TH1D* histo2                = (TH1D*)histo->Clone("Dummy");
        for( Int_t I = 1; I < histo->GetNbinsX() +1 ;I++){
            Double_t xValue         = histo2->GetBinCenter(I);
            Double_t yValue         = fit->Eval(xValue);
            if (integrateFunction){
                Double_t xMin       = histo2->GetXaxis()->GetBinLowEdge(I);
                Double_t xMax       = histo2->GetXaxis()->GetBinUpEdge(I);
                yValue              = fit->Integral(xMin,xMax)/(xMax-xMin);
            }
            Double_t formerYValue   = histo2->GetBinContent(I);
            if (yValue != 0){
                histo2->SetBinContent(I,formerYValue/yValue);
                histo2->SetBinError(I,histo2->GetBinError(I)/yValue);
            }
        }
        return histo2;
    }

    //**********************************************************************************************************
    // Calculates the ratio of a histogram and a fit, with the posibility to integrate the function in the same
    // bin width as the data
    //**********************************************************************************************************
    TH1F* CalculateHistoRatioToFit (TH1F* histo, TF1* fit, Bool_t integrateFunction=kFALSE){
        TH1F* histo2                = (TH1F*)histo->Clone("Dummy");
        for( Int_t I = 1; I < histo->GetNbinsX() +1 ;I++){
            Double_t xValue         = histo2->GetBinCenter(I);
            Double_t yValue         = fit->Eval(xValue);
            Double_t formerYValue   = histo2->GetBinContent(I);
            if (integrateFunction){
                Double_t xMin       = histo2->GetXaxis()->GetBinLowEdge(I);
                Double_t xMax       = histo2->GetXaxis()->GetBinUpEdge(I);
                yValue = fit->Integral(xMin,xMax)/(xMax-xMin);
            }

            if (yValue != 0){
                histo2->SetBinContent(I,formerYValue/yValue);
                histo2->SetBinError(I,histo2->GetBinError(I)/yValue);
            }
        }
        return histo2;
    }

    //**********************************************************************************************************
    // Devides the content of a histogram by the respective bin center values
    //**********************************************************************************************************
    TH1D* CorrectHistoToBinCenter (TH1D* histo){
        histo->Sumw2();
        for( Int_t I = 1; I < histo->GetNbinsX() +1 ;I++){
            Double_t xValue         = histo->GetBinCenter(I);
            Double_t yValue         = histo->GetBinContent(I);
            Double_t yValueError    = histo->GetBinError(I);
            histo->SetBinContent(I,yValue/xValue);
            histo->SetBinError(I,yValueError/xValue);
        }
        return histo;
    }

    //**********************************************************************************************************
    // Calculates the ratio of a graph and a fit
    //**********************************************************************************************************
    TGraphAsymmErrors* MultiplyGraphAsymmErrs (TGraphAsymmErrors* graph_1, TGraphAsymmErrors* graph_2){

        // check for same "binning" in x
        Int_t nPoints1          = graph_1->GetN();
        Int_t nPoints2          = graph_2->GetN();
        if (nPoints1 != nPoints2) return NULL;
        Double_t* xValue1       = graph_1->GetX();
        Double_t* xValue2       = graph_2->GetX();
        for (Int_t i=0; i<nPoints1; i++) {
            if (xValue1[i] != xValue2[i]) return NULL;
        }
        Double_t* xErrorLow1    = graph_1->GetEXlow();
        Double_t* xErrorHigh1   = graph_1->GetEXhigh();
        Double_t* xErrorLow2    = graph_2->GetEXlow();
        Double_t* xErrorHigh2   = graph_2->GetEXhigh();
        for (Int_t i=0; i<nPoints1; i++) {
            if (xErrorLow1[i] != xErrorLow2[i])     return NULL;
            if (xErrorHigh1[i] != xErrorHigh2[i])   return NULL;
        }

        // get y values and errors
        Double_t* yValue1       = graph_1->GetY();
        Double_t* yValue2       = graph_2->GetY();
        Double_t* yErrorLow1    = graph_1->GetEYlow();
        Double_t* yErrorHigh1   = graph_1->GetEYhigh();
        Double_t* yErrorLow2    = graph_2->GetEYlow();
        Double_t* yErrorHigh2   = graph_2->GetEYhigh();

        // calculate new yValues and errors
        Double_t* yValueNew     = new Double_t[nPoints1];
        Double_t* yErrorLowNew  = new Double_t[nPoints1];
        Double_t* yErrorHighNew = new Double_t[nPoints1];
        for (Int_t i=0; i<nPoints1; i++) {
            yValueNew[i]        = yValue1[i] * yValue2[i];
            yErrorLowNew[i]     = TMath::Sqrt( TMath::Power(yErrorLow1[i]*yValue2[i],2) + TMath::Power(yErrorLow2[i]*yValue1[i],2) );
            yErrorHighNew[i]    = TMath::Sqrt( TMath::Power(yErrorHigh1[i]*yValue2[i],2) + TMath::Power(yErrorHigh2[i]*yValue1[i],2) );
        }

        TGraphAsymmErrors* graph = new TGraphAsymmErrors(nPoints1, xValue1, yValueNew, xErrorLow1, xErrorHigh1, yErrorLowNew, yErrorHighNew);
        graph->SetName(Form("%sTimes%s", graph_1->GetName(), graph_2->GetName()));
        return graph;
    }

    //**********************************************************************************************************
    // Calculates the ratio of two error graphs
    //**********************************************************************************************************
    TGraphErrors* CalculateGraphRatioToGraph(TGraphErrors* graphA, TGraphErrors* graphB){

        TGraphErrors* graphACopy        = (TGraphErrors*)graphA->Clone("GraphCopy");
        Double_t* xValue                = graphACopy->GetX();
        Double_t* yValue                = graphACopy->GetY();
        Double_t* xError                = graphACopy->GetEX();
        Double_t* yError                = graphACopy->GetEY();
        Int_t nPoints                   = graphACopy->GetN();

        for (Int_t i = 0; i < nPoints; i++){
            yValue[i]                   = yValue[i]/graphB->GetY()[i];
            cout<<i<<" A EY "<< graphA->GetEY()[i]<<" A Y "<<graphA->GetY()[i]<<" B EY "<<graphB->GetEY()[i]<<" B Y"<<graphB->GetY()[i]<<endl;
            Double_t yErrorRatio        = yValue[i]*TMath::Sqrt( TMath::Power(graphA->GetEY()[i]/graphA->GetY()[i],2) + TMath::Power(graphB->GetEY()[i]/graphB->GetY()[i],2));
            yError[i] = TMath::Abs(yErrorRatio);
        }
        TGraphErrors* returnGraph       = new TGraphErrors(nPoints,xValue,yValue,xError,yError);
        return returnGraph;
    }

    //**********************************************************************************************************
    // Calculates the ratio of two asymmerror graphs
    //**********************************************************************************************************
    TGraphAsymmErrors* CalculateAsymGraphRatioToGraph(TGraphAsymmErrors* graphA, TGraphAsymmErrors* graphB){

        TGraphErrors* graphACopy        = (TGraphErrors*)graphA->Clone("GraphCopy");
        Double_t* xValue                = graphACopy->GetX();
        Double_t* yValue                = graphACopy->GetY();
        Double_t* xErrorHigh            = graphACopy->GetEXhigh();
        Double_t* xErrorLow             = graphACopy->GetEXlow();
        Double_t* yErrorHigh            = graphACopy->GetEYhigh();
        Double_t* yErrorLow             = graphACopy->GetEYlow();

        Int_t nPoints                   = graphACopy->GetN();
        for (Int_t i = 0; i < nPoints; i++){
            if (abs(xValue[i]-graphB->GetX()[i]) < 0.0001){
                if (graphB->GetY()[i] != 0){
                    yValue[i]                   = yValue[i]/graphB->GetY()[i];
                    Double_t yErrorRatioHigh    = yValue[i]*TMath::Sqrt( TMath::Power(graphA->GetEYhigh()[i]/graphA->GetY()[i],2) + TMath::Power(graphB->GetEYhigh()[i]/graphB->GetY()[i],2));
                    yErrorHigh[i]               = TMath::Abs(yErrorRatioHigh);
                    Double_t yErrorRatioLow     = yValue[i]*TMath::Sqrt( TMath::Power(graphA->GetEYlow()[i]/graphA->GetY()[i],2) + TMath::Power(graphB->GetEYlow()[i]/graphB->GetY()[i],2));
                    yErrorLow[i]                = TMath::Abs(yErrorRatioLow);
                } else {
                    yValue[i]                   = 0;
                    yErrorHigh[i]               = 0;
                    yErrorLow[i]                = 0;
                }
            } else {
                cout << "ERROR: graphs don't have same binning in bin " << i << " with x1 = " << xValue[i] << " and x2 = " << graphB->GetX()[i] << endl;
                return NULL;
            }
        }
        TGraphAsymmErrors* returnGraph      = new TGraphAsymmErrors(nPoints,xValue,yValue,xErrorLow,xErrorHigh,yErrorLow,yErrorHigh);
        return returnGraph;
    }

    //**********************************************************************************************************
    // Calculates the difference of two asymmerror graphs
    //**********************************************************************************************************
    TGraphAsymmErrors* CalculateAsymGraphDifferenceToGraph(TGraphAsymmErrors* graphA, TGraphAsymmErrors* graphB){

        TGraphErrors* graphACopy        = (TGraphErrors*)graphA->Clone("GraphCopy");
        Double_t* xValue                = graphACopy->GetX();
        Double_t* yValue                = graphACopy->GetY();
        Double_t* xErrorHigh            = graphACopy->GetEXhigh();
        Double_t* xErrorLow             = graphACopy->GetEXlow();
        Double_t* yErrorHigh            = graphACopy->GetEYhigh();
        Double_t* yErrorLow             = graphACopy->GetEYlow();

        Int_t nPoints                   = graphACopy->GetN();
        for (Int_t i = 0; i < nPoints; i++){
            if (abs(xValue[i]-graphB->GetX()[i]) < 0.0001){
                yValue[i]                   = yValue[i]-graphB->GetY()[i];
                Double_t yErrorRatioHigh    = TMath::Sqrt( TMath::Power(graphA->GetEYhigh()[i],2) + TMath::Power(graphB->GetEYhigh()[i],2));
                yErrorHigh[i]               = yErrorRatioHigh;
                Double_t yErrorRatioLow     = TMath::Sqrt( TMath::Power(graphA->GetEYlow()[i],2) + TMath::Power(graphB->GetEYlow()[i],2));
                yErrorLow[i]                = yErrorRatioLow;
            } else {
                cout << "ERROR: graphs don't have same binning " << endl;
                return NULL;
            }
        }
        TGraphAsymmErrors* returnGraph      = new TGraphAsymmErrors(nPoints,xValue,yValue,xErrorLow,xErrorHigh,yErrorLow,yErrorHigh);
        return returnGraph;
    }

    //**********************************************************************************************************
    // Calculates the ratio of a graph and a fit
    //**********************************************************************************************************
    TGraphErrors* CalculateGraphErrRatioToFit (TGraphErrors* graph_Org, TF1* fit){
        TGraphErrors* graph         = (TGraphErrors*)graph_Org->Clone("Dummy");
        Double_t * xValue           = graph->GetX();
        Double_t * yValue           = graph->GetY();
        Double_t* xError            = graph->GetEX();
        Double_t* yError            = graph->GetEY();
        Int_t nPoints               = graph->GetN();
        for (Int_t i = 0; i < nPoints; i++){
            yValue[i]               = yValue[i]/fit->Eval(xValue[i]);
            yError[i]               = yError[i]/fit->Eval(xValue[i]);
        }
        TGraphErrors* returnGraph   = new TGraphErrors(nPoints,xValue,yValue,xError,yError);
        return returnGraph;
    }

    //**********************************************************************************************************
    // Calculates the ratio of a graph and a spline
    //**********************************************************************************************************
    TGraphErrors* CalculateGraphErrRatioToSpline (TGraphErrors* graph_Org, TSpline* fit){
        TGraphErrors* graph         = (TGraphErrors*)graph_Org->Clone("Dummy");
        Double_t * xValue           = graph->GetX();
        Double_t * yValue           = graph->GetY();
        Double_t* xError            = graph->GetEX();
        Double_t* yError            = graph->GetEY();
        Int_t nPoints               = graph->GetN();
        for (Int_t i = 0; i < nPoints; i++){
            yValue[i]               = yValue[i]/fit->Eval(xValue[i]);
            yError[i]               = yError[i]/fit->Eval(xValue[i]);
        }
        TGraphErrors* returnGraph   = new TGraphErrors(nPoints,xValue,yValue,xError,yError);
        return returnGraph;
    }


    //**********************************************************************************************************
    // Calculates the ratio of a graph and a fit
    //**********************************************************************************************************
    TGraphAsymmErrors* CalculateGraphErrRatioToFit (TGraphAsymmErrors* graph_Org, TF1* fit){
        TGraphAsymmErrors* graph        = (TGraphAsymmErrors*)graph_Org->Clone("Dummy");
        Double_t * xValue               = graph->GetX();
        Double_t * yValue               = graph->GetY();
        Double_t* xErrorLow             = graph->GetEXlow();
        Double_t* xErrorHigh            = graph->GetEXhigh();
        Double_t* yErrorLow             = graph->GetEYlow();
        Double_t* yErrorHigh            = graph->GetEYhigh();
        Int_t nPoints                   = graph->GetN();
        for (Int_t i = 0; i < nPoints; i++){
            yValue[i]                   = yValue[i]/fit->Eval(xValue[i]);
            yErrorLow[i]                = yErrorLow[i]/fit->Eval(xValue[i]);
            yErrorHigh[i]               = yErrorHigh[i]/fit->Eval(xValue[i]);
        }
        TGraphAsymmErrors* returnGraph  = new TGraphAsymmErrors(nPoints,xValue,yValue,xErrorLow,xErrorHigh,yErrorLow,yErrorHigh);
        return returnGraph;
    }

    //**********************************************************************************************************
    // Calculates the ratio of a graph and a fit
    //**********************************************************************************************************
    TGraph* CalculateGraphRatioToFit (TGraph* graph, TF1* fit){
        Double_t * xValue       = graph->GetX();
        Double_t * yValue       = graph->GetY();
        Int_t nPoints           = graph->GetN();
        for (Int_t i = 0; i < nPoints; i++){
            yValue[i]           = yValue[i]/fit->Eval(xValue[i]);
        }
        TGraph* returnGraph     = new TGraph(nPoints,xValue,yValue);
        return returnGraph;
    }

    //**********************************************************************************************************
    // Calculates the ratio of a graph and a NLO fit
    //**********************************************************************************************************
    TH1D* CalculateHistoRatioToFitNLO (TH1D* histo, TF1* fit, Double_t startX){
        Double_t beginBin           = histo->GetXaxis()->FindBin(startX);
        for (Int_t L = 0; L < beginBin; L++){
            histo->SetBinContent(L+1,-1);
        }
        for (Int_t L = beginBin; L < histo->GetNbinsX()+1; L++){
            Double_t xValue         = histo->GetBinCenter(L);
            Double_t yValue         = fit->Eval(xValue);
            Double_t formerYValue   = histo->GetBinContent(L);
            if (formerYValue > 0){
                histo->SetBinContent(L,formerYValue/yValue);
                histo->SetBinError(L,0);
            }
        }
        return histo;
    }

    //**********************************************************************************************************
    // Calculates the ratio of a graph and a fit
    //**********************************************************************************************************
    TGraphAsymmErrors* CalculateGraphErrMultiplicationOfFit (TGraphAsymmErrors* graph_Org, TF1* fit){
        TGraphAsymmErrors* graph        = (TGraphAsymmErrors*)graph_Org->Clone("Dummy");
        Double_t * xValue               = graph->GetX();
        Double_t * yValue               = graph->GetY();
        Double_t* xErrorLow             = graph->GetEXlow();
        Double_t* xErrorHigh            = graph->GetEXhigh();
        Double_t* yErrorLow             = graph->GetEYlow();
        Double_t* yErrorHigh            = graph->GetEYhigh();
        Int_t nPoints                   = graph->GetN();
        for (Int_t i = 0; i < nPoints; i++){
            yValue[i]                   = yValue[i]*fit->Eval(xValue[i]);
            yErrorLow[i]                = yErrorLow[i]*fit->Eval(xValue[i]);
            yErrorHigh[i]               = yErrorHigh[i]*fit->Eval(xValue[i]);
        }
        TGraphAsymmErrors* returnGraph  = new TGraphAsymmErrors(nPoints,xValue,yValue,xErrorLow,xErrorHigh,yErrorLow,yErrorHigh);
        return returnGraph;
    }

    //**********************************************************************************************************
    // Calculates a graph with systematic errors based on an input histogram and two arrays of doubles
    // containing the relative systematic errors
    //**********************************************************************************************************
    TGraphAsymmErrors* CalculateSysErrFromRelSysHisto( TH1D* histo, TString nameGraph, Double_t* relSystErrorDown, Double_t* relSystErrorUp, Int_t offset, const Int_t nPoints ){
        Double_t xValueCorr[nPoints];
        Double_t yValueCorr[nPoints];
        //    Double_t yErrorCorr[nPoints];
        Double_t xErrorCorr[nPoints];
        Double_t systErrorDown[nPoints];
        Double_t systErrorUp[nPoints];
        for (Int_t i = 0; i < nPoints; i++){
            xValueCorr[i]           = histo->GetBinCenter(i+offset);
            yValueCorr[i]           = histo->GetBinContent(i+offset);
            xErrorCorr[i]           = histo->GetBinWidth(i+offset)/2.;
            systErrorDown[i]        = relSystErrorDown[i]*yValueCorr[i]/100.*(-1);
            systErrorUp[i]          = relSystErrorUp[i]*yValueCorr[i]/100.;
            cout << i << "\t" << i+offset <<"\t" <<xValueCorr[i] << "\t" << yValueCorr[i] << "\t" << relSystErrorDown[i] << endl;
        }
        cout << "create graph" << endl;
        TGraphAsymmErrors* graphSys = new TGraphAsymmErrors(nPoints,xValueCorr,yValueCorr,xErrorCorr,xErrorCorr,systErrorDown,systErrorUp);
        graphSys->SetName(nameGraph);
        cout << "graph could be generated" << endl;
        return graphSys;
    }

    //**********************************************************************************************************
    // Calculates a graph with systematic errors based on an input histogram and two arrays of doubles
    // containing the relative systematic errors
    //**********************************************************************************************************
    TGraphAsymmErrors* CalculateSysErrFromRelSysHistoWithPtBins( TH1D* histo, TString nameGraph, Double_t* relSystErrorDown, Double_t* relSystErrorUp, Double_t* ptCenter, Int_t nMaxPtSys ){
        TGraphAsymmErrors* graphSysDummy  = new TGraphAsymmErrors(histo);

        Bool_t kGraphLower          = kFALSE;
        while (graphSysDummy->GetX()[0] < ptCenter[0] && graphSysDummy->GetN() > 1){
            graphSysDummy->RemovePoint(0);
            kGraphLower             = kTRUE;
        }
        if (graphSysDummy->GetN() == 0) {
            cout << "ERROR: could not find a common bin" << endl;
            return NULL;
        }
        Int_t startBinPt            = 0;
        Int_t counter               = nMaxPtSys-1;
        if (!kGraphLower && (TMath::Abs(graphSysDummy->GetX()[0] - ptCenter[startBinPt])>0.001)){
            while( TMath::Abs(graphSysDummy->GetX()[0] - ptCenter[startBinPt])>0.001 && counter > 0){
                startBinPt++;
                counter--;
            }
        }
        if (counter == 0){
            cout << "ERROR: could not find a common bin" << endl;
            return NULL;
        }
        TGraphAsymmErrors* graphSys      = (TGraphAsymmErrors*)graphSysDummy->Clone(nameGraph.Data());

        for (Int_t i = 0; i < graphSys->GetN() && i < (nMaxPtSys-startBinPt); i++){
            Double_t xErrorCorrUp   = graphSysDummy->GetEXhigh()[i];
            Double_t xErrorCorrDown = graphSysDummy->GetEXlow()[i];
            Double_t yvalue         = graphSysDummy->GetY()[i];
            Double_t yErrorUp       = 0;
            Double_t yErrorDown     = 0;
            if ( TMath::Abs(graphSysDummy->GetX()[i] - ptCenter[startBinPt+i])<0.001) {
                yErrorUp            = relSystErrorUp[startBinPt+i]*yvalue/100.;
                yErrorDown          = relSystErrorDown[startBinPt+i]*yvalue/100.*(-1);
            } else {
            cout << "couldn't find matching bin: "    << graphSysDummy->GetX()[i] << "\t" << ptCenter[startBinPt+i] <<endl;
            }
            graphSys->SetPointError(i,xErrorCorrDown, xErrorCorrUp, yErrorDown, yErrorUp);
        }
        cout << "return graph" << endl;
        return graphSys;
    }

    //**********************************************************************************************************
    // Calculates a graph with systematic errors based on an input histogram and two arrays of doubles
    // containing the relative systematic errors
    //**********************************************************************************************************
    TGraphAsymmErrors* CalculateSysErrAFromRelSysHisto( TH1D* histo, TString nameGraph, Double_t* relSystErrorDown, Double_t* relSystErrorUp, Int_t offset, const Int_t nPoints ){
        Double_t xValueCorr[nPoints];
        Double_t yValueCorr[nPoints];
        //    Double_t yErrorCorr[nPoints];
        Double_t xErrorCorr[nPoints];
        Double_t systErrorADown[nPoints];
        Double_t systErrorAUp[nPoints];

        Double_t relSystErrorADown[nPoints];
        Double_t relSystErrorAUp[nPoints];
        Double_t relSysMatDown          = 2.*6.2;
        Double_t relSysMatUp            = 2.*3.4;
        for (Int_t i = 0; i < nPoints; i++){
            if( TMath::Abs(relSystErrorDown[i])>=relSysMatDown &&  TMath::Abs(relSystErrorUp[i])>=relSysMatUp){
                relSystErrorADown[i]    = -TMath::Sqrt(relSystErrorDown[i]*relSystErrorDown[i]-relSysMatDown*relSysMatDown);
                relSystErrorAUp[i]      = TMath::Sqrt(relSystErrorUp[i]*relSystErrorUp[i]-relSysMatUp*relSysMatUp);
            }else{
                relSystErrorADown[i]    = -TMath::Sqrt(relSystErrorDown[i]*relSystErrorDown[i]);
                relSystErrorAUp[i]      = TMath::Sqrt(relSystErrorUp[i]*relSystErrorUp[i]);
            }
            xValueCorr[i]               = histo->GetBinCenter(i+offset);
            yValueCorr[i]               = histo->GetBinContent(i+offset);
            xErrorCorr[i]               = histo->GetBinWidth(i+offset)/2.;
            systErrorADown[i]           = relSystErrorADown[i]*yValueCorr[i]/100.*(-1);
            systErrorAUp[i]             = relSystErrorAUp[i]*yValueCorr[i]/100.;
            cout << xValueCorr[i] << "\t" << yValueCorr[i] << endl;

        }

        TGraphAsymmErrors* graphSys     = new TGraphAsymmErrors(nPoints,xValueCorr,yValueCorr,xErrorCorr,xErrorCorr,systErrorADown,systErrorAUp);
        graphSys->SetName(nameGraph);
        return graphSys;
    }

    //**********************************************************************************************************
    // Calculates a graph with systematic errors based on an input histogram and two arrays of doubles
    // containing the relative systematic errors
    //**********************************************************************************************************
    TGraphAsymmErrors* CalculateSysErrFromRelSysHistoComplete( TH1D* histo, TString nameGraph, Double_t* relSystErrorDown, Double_t* relSystErrorUp, Int_t offset, const Int_t nPoints ){
        Double_t xValueCorr[nPoints];
        Double_t yValueCorr[nPoints];
        Double_t yErrorCorr[nPoints];
        Double_t xErrorCorr[nPoints];
        Double_t systErrorDown[nPoints];
        Double_t systErrorUp[nPoints];
        for (Int_t i = 0; i < nPoints; i++){
            xValueCorr[i]               = histo->GetBinCenter(i+offset);
            xErrorCorr[i]               = histo->GetBinWidth(i+offset)/2.;
            yValueCorr[i]               = histo->GetBinContent(i+offset);
            yErrorCorr[i]               = histo->GetBinError(i+offset);
            systErrorDown[i]            = TMath::Sqrt((TMath::Power(yErrorCorr[i],2) + TMath::Power(relSystErrorDown[i]*yValueCorr[i]/100.,2)));
            systErrorUp[i]              = TMath::Sqrt(TMath::Power(yErrorCorr[i],2) + TMath::Power(relSystErrorUp[i]*yValueCorr[i]/100.,2));
            cout << xValueCorr[i] << "\t" << yValueCorr[i] << endl;
        }

        TGraphAsymmErrors* graphSys     = new TGraphAsymmErrors(nPoints,xValueCorr,yValueCorr,xErrorCorr,xErrorCorr,systErrorDown,systErrorUp);

        graphSys->SetName(nameGraph);
        return graphSys;
    }

    //**********************************************************************************************************
    // Calculates a graph with total errors based on an 2 input graphs containing systematic and statistical
    // errors
    //**********************************************************************************************************


    TGraphErrors* CalculateCombinedSysAndStatError( TGraphErrors* graphStat, TGraphErrors* graphSys){
    TGraphErrors* graphStatCopy          = (TGraphErrors*)graphStat->Clone("graphStatCopy");
    TGraphErrors* graphSysCopy           = (TGraphErrors*)graphSys->Clone("graphSysCopy");
    Double_t* xValue                     = graphStatCopy->GetX();
    Double_t* xError                     = graphStatCopy->GetEX();
    Double_t* yValueStat                 = graphStatCopy->GetY();
    Double_t* yErrorStat                 = graphStatCopy->GetEY();
    Double_t* yErrorSys                  = graphSysCopy->GetEY();
    Int_t nPoints                        = graphStatCopy->GetN();

    Double_t yErrorComb[nPoints];
    // Double_t yErrorComb[nPoints];
    for (Int_t i = 0; i < nPoints; i++){
        yErrorComb[i]                 = TMath::Sqrt((TMath::Power(yErrorStat[i],2) + TMath::Power(yErrorSys[i],2)));
    }

    TGraphErrors* graphSysAndComb   = new TGraphErrors(nPoints,xValue,yValueStat,xError,yErrorComb);

    graphSysAndComb->SetName(graphSys->GetName());
    return graphSysAndComb;
    }

    TGraphAsymmErrors* CalculateCombinedSysAndStatError( TGraphAsymmErrors* graphStat, TGraphAsymmErrors* graphSys){
    TGraphAsymmErrors* graphStatCopy     = (TGraphAsymmErrors*)graphStat->Clone("graphStatCopy");
    TGraphAsymmErrors* graphSysCopy      = (TGraphAsymmErrors*)graphSys->Clone("graphSysCopy");
    Double_t* xValue                     = graphStatCopy->GetX();
    Double_t* xErrorLow                  = graphStatCopy->GetEXlow();
    Double_t* xErrorHigh                 = graphStatCopy->GetEXhigh();
    Double_t* yValueStat                 = graphStatCopy->GetY();
    Double_t* yErrorLowStat              = graphStatCopy->GetEYlow();
    Double_t* yErrorHighStat             = graphStatCopy->GetEYhigh();
    Double_t* yErrorLowSys               = graphSysCopy->GetEYlow();
    Double_t* yErrorHighSys              = graphSysCopy->GetEYhigh();
    Int_t nPoints                        = graphStatCopy->GetN();
    Double_t yErrorLowComb[nPoints];
    Double_t yErrorHighComb[nPoints];
    for (Int_t i = 0; i < nPoints; i++){
        yErrorLowComb[i]                  = TMath::Sqrt((TMath::Power(yErrorLowStat[i],2) + TMath::Power(yErrorLowSys[i],2)));
        yErrorHighComb[i]                 = TMath::Sqrt((TMath::Power(yErrorHighStat[i],2) + TMath::Power(yErrorHighSys[i],2)));

    }

    TGraphAsymmErrors* graphSysAndComb   = new TGraphAsymmErrors(nPoints,xValue,yValueStat,xErrorLow,xErrorHigh,yErrorLowComb,yErrorHighComb);

    graphSysAndComb->SetName(graphSys->GetName());
    return graphSysAndComb;
    }

    //**********************************************************************************************************
    // Calculates the ratio of two graphs in the same binning
    //**********************************************************************************************************
    TGraphAsymmErrors* CalculateGraphErrRatioToOtherTGraphErr(TGraphAsymmErrors* graph_Org,
                                                            TGraph* graphToDivBy,
                                                            Bool_t calCombinedErr = kFALSE
                                                            ){
        TGraphAsymmErrors* graph            = (TGraphAsymmErrors*)graph_Org->Clone("Dummy");
        TGraphAsymmErrors* graphToDivBy2    = (TGraphAsymmErrors*)graphToDivBy->Clone("graphToDivBy2");
        Int_t j                             = 0;
        for (Int_t i = 0; i < graph->GetN(); i++){
    //         cout << i << "\t" << j<< endl;
    //         cout << "graph A: " << graph->GetX()[i] << endl;
    //         cout << "graph B: " << graphToDivBy2->GetX()[j] << endl;
            while (abs(graph->GetX()[i] - graphToDivBy2->GetX()[j]) > 0.0001 && j < graphToDivBy2->GetN()){
    //             cout << "graph B: " << graphToDivBy2->GetX()[j] << endl;
                if (graph->GetX()[i] > graphToDivBy2->GetX()[j]){
                    graphToDivBy2->RemovePoint(j);
                }else {
                    cout << "I failed" << endl;
                    return 0x0;
                }
    //             cout << "next graph B: " << graphToDivBy2->GetX()[j] << endl;
            }
    //         cout << "correct" << endl;
            j++;
        }
        Double_t* xValue                    = graph->GetX();
        Double_t* yValue                    = graph->GetY();
        Double_t* xErrorLow                 = graph->GetEXlow();
        Double_t* xErrorHigh                = graph->GetEXhigh();
        Double_t* yErrorLow                 = graph->GetEYlow();
        Double_t* yErrorHigh                = graph->GetEYhigh();
        Int_t nPoints                       = graph->GetN();
        for (Int_t i = 0; i < nPoints; i++){
            yValue[i]                       = yValue[i]/graphToDivBy2->GetY()[i];
            if (!calCombinedErr){
                yErrorLow[i]                = yErrorLow[i]/graphToDivBy2->GetY()[i];
                yErrorHigh[i]               = yErrorHigh[i]/graphToDivBy2->GetY()[i];
            } else {
                cout << i << endl;
                cout<< yValue[i] << "\t" << graph_Org->GetEYlow()[i]/graph_Org->GetY()[i] << "\t" << graphToDivBy2->GetEYlow()[i]/graphToDivBy2->GetY()[i] << endl;
                Double_t intermediateLow    = yValue[i]* TMath::Sqrt( TMath::Power(graph_Org->GetEYlow()[i]/graph_Org->GetY()[i],2)+ TMath::Power(graphToDivBy2->GetEYlow()[i]/graphToDivBy2->GetY()[i],2));
                cout << intermediateLow << endl;
                Double_t intermediateHigh   = yValue[i]* TMath::Sqrt( TMath::Power(graph_Org->GetEYhigh()[i]/graph_Org->GetY()[i],2)+ TMath::Power(graphToDivBy2->GetEYhigh()[i]/graphToDivBy2->GetY()[i],2));
                yErrorLow[i]                = intermediateLow;
                yErrorHigh[i]               = intermediateHigh;
            }
        }
        TGraphAsymmErrors* returnGraph      = new TGraphAsymmErrors(nPoints,xValue,yValue,xErrorLow,xErrorHigh,yErrorLow,yErrorHigh);
        return returnGraph;
    }

    //**********************************************************************************************************
    // Calculates shifted spectrum in xt or mt
    //**********************************************************************************************************
    TH1D* CalculateShiftedSpectrumInXtOrMt(TH1D* ptSpectrum, TH1D* spectrum, TString nameNewSpectrum){
        Double_t dummyValue     = 0;
        Double_t dummyError     = 0;
        TH1D* shiftedSpectrum   = (TH1D*)spectrum->Clone(nameNewSpectrum);
        for (Int_t i = 0; i < shiftedSpectrum->GetNbinsX()+1; i++){
            dummyValue          = ptSpectrum->GetBinContent(i);
            dummyError          = ptSpectrum->GetBinError(i);
            shiftedSpectrum->SetBinContent(i,dummyValue);
            shiftedSpectrum->SetBinError(i,dummyError);
        }
        return shiftedSpectrum;
    }

    //**********************************************************************************************************
    // Calculates mean mass from fit - expected mass
    //**********************************************************************************************************
    TH1D* CalculateMassMinusExpectedMass(TH1D* histo, Double_t mass){
        TH1D* histoNew                  = (TH1D*) histo->Clone();
        for ( Int_t l=0; l < histoNew->GetNbinsX()+1; l++){
            Double_t intermediateValue  = histoNew->GetBinContent(l);
            Double_t intermediateError  = histoNew->GetBinError(l);
            if (intermediateValue != 0) {
                histoNew->SetBinContent(l,(intermediateValue - mass)*1000);
                histoNew->SetBinError(l,intermediateError*1000);
            }
        }
        return histoNew;
    }

    //**********************************************************************************************************
    // Takes out scaling with pT bin by bin for histograms
    //**********************************************************************************************************
    void RemoveScalingWithPt(TH1D* histo){
        for (Int_t i = 0; i < histo->GetNbinsX()+1;i++){
            histo->SetBinContent(i, histo->GetBinContent(i)*histo->GetBinCenter(i));
            histo->SetBinError(i, histo->GetBinError(i)*histo->GetBinCenter(i));
        }
        return;
    }

    //**********************************************************************************************************
    // Takes out scaling with pT bin by bin for TGraphAsymmErrors
    //**********************************************************************************************************
    void RemoveScalingWithPtGraph(TGraphAsymmErrors* graph){
        Double_t * xValue       = graph->GetX();
        Double_t * yValue       = graph->GetY();
        Double_t* yErrorHigh    = graph->GetEYhigh();
        Double_t* yErrorLow     = graph->GetEYlow();
        Int_t nPoints           = graph->GetN();
        for (Int_t i = 0; i < nPoints; i++){
            yValue[i]           = yValue[i]*xValue[i];
            yErrorHigh[i]       = yErrorHigh[i]*xValue[i];
            yErrorLow[i]        = yErrorLow[i]*xValue[i];
        }
        return;
    }

    //**********************************************************************************************************
    // Scales every bin with for TGraphAsymmErrors
    //**********************************************************************************************************
    void ScaleWithPtGraph(TGraphAsymmErrors* graph){
        Double_t * xValue       = graph->GetX();
        Double_t * yValue       = graph->GetY();
        Double_t* yErrorHigh    = graph->GetEYhigh();
        Double_t* yErrorLow     = graph->GetEYlow();
        Int_t nPoints           = graph->GetN();
        for (Int_t i = 0; i < nPoints; i++){
            yValue[i]           = yValue[i]/xValue[i];
            yErrorHigh[i]       = yErrorHigh[i]/xValue[i];
            yErrorLow[i]        = yErrorLow[i]/xValue[i];
        }
        return;
    }

    // ****************************************************************************************************************
    // ****************************************************************************************************************
    // ****************************************************************************************************************
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

    // ****************************************************************************************************************
    // ****************************************************************************************************************
    // ****************************************************************************************************************
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

    // ****************************************************************************************************************
    // ****************************************************************************************************************
    // ****************************************************************************************************************
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

    // ****************************************************************************************************************
    // ****************************************************************************************************************
    // ****************************************************************************************************************
    Double_t* ExtractRelErrDownAsymmGraph(TGraphAsymmErrors* graph){
        Double_t* yValue        = graph->GetY();
        Double_t* yErrorLow     = graph->GetEYlow();
        Double_t* yErrorLow2    = yErrorLow;
        Int_t nPoints           = graph->GetN();
        for (Int_t i = 0; i < nPoints; i++){
            yErrorLow2[i]       = -yErrorLow2[i]/yValue[i]*100;
        }
        return yErrorLow2;
    }

    // ****************************************************************************************************************
    // ****************************************************************************************************************
    // ****************************************************************************************************************
    Double_t* ExtractRelErrUpAsymmGraph(TGraphAsymmErrors* graph){
        Double_t* yValue        = graph->GetY();
        Double_t* yErrorHigh    = graph->GetEYhigh();
        Double_t* yErrorHigh2   = yErrorHigh;
        Int_t nPoints           = graph->GetN();

        for (Int_t i = 0; i < nPoints; i++){
            yErrorHigh2[i]      = yErrorHigh2[i]/yValue[i]*100;
        }
        return yErrorHigh2;
    }

    // ****************************************************************************************************************
    // ****************************************************************************************************************
    // ****************************************************************************************************************
    TGraphAsymmErrors* CalculateRelErrUpAsymmGraph( TGraphAsymmErrors* graph,
                                                    TString nameNewGraph            = "relativeError"
                                                ){
        TGraphAsymmErrors* returnGraph  = (TGraphAsymmErrors*)graph->Clone(nameNewGraph.Data());
        Double_t* yValue                = returnGraph->GetY();
        Double_t* yErrorHigh            = returnGraph->GetEYhigh();
        Double_t* yErrorLow             = returnGraph->GetEYlow();
        Int_t nPoints                   = returnGraph->GetN();

        for (Int_t i = 0; i < nPoints; i++){
            yValue[i]                   = yErrorHigh[i]/yValue[i]*100;
            yErrorHigh[i]               = 0.;
            yErrorLow[i]                = 0.;
        }
        return returnGraph;
    }

    // ****************************************************************************************************************
    // ****************************************************************************************************************
    // ****************************************************************************************************************
    TGraphAsymmErrors* CalculateRelErrAsymmGraphAround1(    TGraphAsymmErrors* graph,
                                                            TString nameNewGraph        = "relativeError"
                                                    ){
        TGraphAsymmErrors* returnGraph  = (TGraphAsymmErrors*)graph->Clone(nameNewGraph.Data());
        Double_t* yValue                = returnGraph->GetY();
        Double_t* yErrorHigh            = returnGraph->GetEYhigh();
        Double_t* yErrorLow             = returnGraph->GetEYlow();
        Int_t nPoints                   = returnGraph->GetN();

        for (Int_t i = 0; i < nPoints; i++){
            yErrorHigh[i]               = yErrorHigh[i]/yValue[i];
            yErrorLow[i]                = yErrorLow[i]/yValue[i];
            yValue[i]                   = 1;
        }
        return returnGraph;
    }

    // ****************************************************************************************************************
    // ****************************************************************************************************************
    // ****************************************************************************************************************
    TGraphErrors* CalculateRelErrGraphAround1(  TGraphErrors* graph,
                                                TString nameNewGraph    = "relativeError"
                                            ){
        TGraphErrors* returnGraph       = (TGraphErrors*)graph->Clone(nameNewGraph.Data());
        Double_t* yValue                = returnGraph->GetY();
        Double_t* yError                = returnGraph->GetEY();
        Int_t nPoints                   = returnGraph->GetN();

        for (Int_t i = 0; i < nPoints; i++){
            yError[i]                   = yError[i]/yValue[i];
            yValue[i]                   = 1;
        }
        return returnGraph;
    }

    // ****************************************************************************************************************
    // ****************************************************************************************************************
    // ****************************************************************************************************************
    TH1D* CalculateRelErrUpTH1D(    TH1D* histo,
                                    TString nameNewGraph    = "relativeError"
                            ){
        TH1D* returnHisto               = (TH1D*)histo->Clone(nameNewGraph.Data());
        cout << nameNewGraph.Data() << endl;
        for (Int_t i = 1; i < returnHisto->GetNbinsX()+1; i++){
            if (abs(returnHisto->GetBinContent(i)) > 0){
                cout <<i << "\t" << returnHisto->GetBinContent(i) << "\t"<< returnHisto->GetBinError(i) << "\t" << returnHisto->GetBinError(i)/returnHisto->GetBinContent(i)*100 << endl;
                returnHisto->SetBinContent(i,returnHisto->GetBinError(i)/returnHisto->GetBinContent(i)*100);
                returnHisto->SetBinError(i,0);
            } else {
                returnHisto->SetBinContent(i,-1);
                returnHisto->SetBinError(i,0);
            }
        }
        return returnHisto;
    }

    // ****************************************************************************************************************
    // ****************************************************************************************************************
    // ****************************************************************************************************************
    TGraphAsymmErrors* CalculateGraphAsymErrRatioToGraphErr (   TGraphAsymmErrors* graphAsymA_v1,
                                                                TGraphAsymmErrors* graphAsymB_v1,
                                                                Bool_t binomial                     = kFALSE,
                                                                Double_t commonPercentageADown      = 0.,
                                                                Double_t commonPercentageAUp        = 0.,
                                                                Double_t commonPercentageBDown      = 0.,
                                                                Double_t commonPercentageBUp        = 0.
                                                            ){
        TGraphAsymmErrors *graphAsymA       = (TGraphAsymmErrors*) graphAsymA_v1->Clone("graphAsymA");

        Double_t* xValueAsymA               = graphAsymA->GetX();
        Double_t* yValueAsymA               = graphAsymA->GetY();
        Double_t* xErrorLowAsymA            = graphAsymA->GetEXlow();
        Double_t* xErrorHighAsymA           = graphAsymA->GetEXhigh();
        Double_t* yErrorLowAsymA            = graphAsymA->GetEYlow();
        Double_t* yErrorHighAsymA           = graphAsymA->GetEYhigh();

        TGraphAsymmErrors *graphAsymB       = (TGraphAsymmErrors*) graphAsymB_v1->Clone("graphAsymB");
        //    Double_t * xValueAsymB            = graphAsymB->GetX();
        Double_t * yValueAsymB              = graphAsymB->GetY();
        Double_t* yErrorLowAsymB            = graphAsymB->GetEYlow();

        //  Double_t* xErrorHighAsymB           = graphAsymB->GetEXhigh();
        //  Double_t* yErrorLowAsymB            = graphAsymB->GetEYlow();
        Double_t* yErrorHighAsymB           = graphAsymB->GetEYhigh();

        Int_t nPointsB                      = graphAsymB->GetN();
        Int_t nPoints                       = graphAsymA->GetN();
        if (nPoints > nPointsB) nPoints     = nPointsB;
        Double_t yValueRatio[nPoints];
        for (Int_t i = 0; i < nPoints; i++){
            yValueRatio[i]                  = yValueAsymA[i]/yValueAsymB[i];
        //       cout << "before low: "<<yErrorLowAsymA[i]/yValueAsymA[i] << "\t" << commonPercentageADown << "\t"  <<  yErrorLowAsymB[i]/yValueAsymB[i] << "\t"  << commonPercentageBDown << endl;
        //       cout << "before up: "<<yErrorHighAsymA[i]/yValueAsymA[i] << "\t" << commonPercentageAUp << "\t"  <<  yErrorHighAsymB[i]/yValueAsymB[i] << "\t"  <<  commonPercentageBUp << endl;
            yErrorLowAsymA[i]               = TMath::Sqrt(pow(yErrorLowAsymA[i],2) - pow(commonPercentageADown*yValueAsymA[i]/100,2));
            yErrorHighAsymA[i]              = TMath::Sqrt(pow(yErrorHighAsymA[i],2) - pow(commonPercentageAUp*yValueAsymA[i]/100,2));
            yErrorLowAsymB[i]               = TMath::Sqrt(pow(yErrorLowAsymB[i],2) - pow(commonPercentageBDown*yValueAsymB[i]/100,2));
            yErrorHighAsymB[i]              = TMath::Sqrt(pow(yErrorHighAsymB[i],2) - pow(commonPercentageBUp*yValueAsymB[i]/100,2));
        //       cout << "after low: "<< yErrorLowAsymA[i]/yValueAsymA[i] << "\t"  <<  yErrorLowAsymB[i]/yValueAsymB[i] << endl;
        //       cout << "after up: "<< yErrorHighAsymA[i]/yValueAsymA[i] << "\t"  <<  yErrorHighAsymB[i]/yValueAsymB[i] << endl;
            if (binomial){
                Double_t w                  = yValueAsymA[i]/yValueAsymB[i];
                yErrorLowAsymA[i]           = TMath::Abs( ( (1.-2.*w)*yErrorLowAsymA[i]*yErrorLowAsymA[i] + w*w*yErrorLowAsymB[i]*yErrorLowAsymB[i] )/(yValueAsymB[i]*yValueAsymB[i]) );
                yErrorHighAsymA[i]          = TMath::Abs( ( (1.-2.*w)*yErrorHighAsymA[i]*yErrorHighAsymA[i] + w*w*yErrorHighAsymB[i]*yErrorHighAsymB[i] )/(yValueAsymB[i]*yValueAsymB[i]) );
            } else {

                yErrorLowAsymA[i]           = TMath::Sqrt( pow(yErrorLowAsymA[i]/yValueAsymB[i],2)  + pow( yErrorLowAsymB[i]*yValueAsymA[i]/pow(yValueAsymB[i],2),2) );
                yErrorHighAsymA[i]          = TMath::Sqrt( pow(yErrorHighAsymA[i]/yValueAsymB[i],2) + pow( yErrorHighAsymB[i]*yValueAsymA[i]/pow(yValueAsymB[i],2),2) );
            }
        }
        TGraphAsymmErrors* returnGraph      = new TGraphAsymmErrors(nPoints,xValueAsymA,yValueRatio,xErrorLowAsymA,xErrorHighAsymA,yErrorLowAsymA,yErrorHighAsymA);
        return returnGraph;
    }

    // ****************************************************************************************************************
    // ****************************************************************************************************************
    // ****************************************************************************************************************
    TGraphAsymmErrors* RebinCombPi0Graph(   TGraphAsymmErrors* graphAsymV1,
                                            TH1D* histoGraph
                                        ){
        TGraphAsymmErrors *graphAsym    = (TGraphAsymmErrors*) graphAsymV1->Clone("graphAsym");

        Double_t * xValueAsym           = graphAsym->GetX();
        Double_t * yValueAsym           = graphAsym->GetY();
        //Double_t* xErrorLowAsym       = graphAsym->GetEXlow();
        //Double_t* xErrorHighAsym      = graphAsym->GetEXhigh();
        Double_t* yErrorLowAsym         = graphAsym->GetEYlow();
        Double_t* yErrorHighAsym        = graphAsym->GetEYhigh();
        Int_t nPoints                   = graphAsym->GetN();

        Int_t nBins                     = histoGraph->GetNbinsX();
        Double_t *newBinningY           = new Double_t[nBins];
        Double_t *newBinningX           = new Double_t[nBins];
        Double_t *newBinningErrorLowY   = new Double_t[nBins];
        Double_t *newBinningErrorLowX   = new Double_t[nBins];
        Double_t *newBinningErrorHighY  = new Double_t[nBins];
        Double_t *newBinningErrorHighX  = new Double_t[nBins];

        Int_t notusedBins               = 0;

        for(Int_t bin = 1; bin <nBins+1; bin++){
            Double_t yValue     = 0;
            Double_t yErrorLow  = 0;
            Double_t yErrorHigh = 0;
            Int_t nMergedBins   = 0;
            Double_t binCenter  = histoGraph->GetBinCenter(bin);
            Double_t binWidth   = histoGraph->GetBinWidth(bin);

            for(Int_t graphBin  = 0; graphBin<nPoints; graphBin++){
                if(histoGraph->FindBin(xValueAsym[graphBin]) == bin){
                    yValue      = yValue+yValueAsym[graphBin];
                    yErrorLow   = yErrorLow+yErrorLowAsym[graphBin];
                    yErrorHigh  = yErrorHigh+yErrorHighAsym[graphBin];
                    nMergedBins++;
                }
            }
            if(nMergedBins == 0){
                notusedBins++;
                continue;
            }

            newBinningY[bin-1-notusedBins]          = yValue/nMergedBins;
            newBinningX[bin-1-notusedBins]          = binCenter;
            newBinningErrorLowX[bin-1-notusedBins]  = binWidth/2;
            newBinningErrorHighX[bin-1-notusedBins] = binWidth/2;
            newBinningErrorLowY[bin-1-notusedBins]  = yErrorLow/nMergedBins;
            newBinningErrorHighY[bin-1-notusedBins] = yErrorHigh/nMergedBins;

            histoGraph->SetBinContent(bin,newBinningY[bin-1-notusedBins]);
            if(newBinningErrorLowY[bin-1-notusedBins]>=newBinningErrorHighY[bin-1-notusedBins]) histoGraph->SetBinError(bin,newBinningErrorLowY[bin-1-notusedBins]);
            else histoGraph->SetBinError(bin,newBinningErrorHighY[bin-1-notusedBins]);
        }

        Double_t *newnewBinningY                = new Double_t[nBins-notusedBins];
        Double_t *newnewBinningX                = new Double_t[nBins-notusedBins];
        Double_t *newnewBinningErrorLowY        = new Double_t[nBins-notusedBins];
        Double_t *newnewBinningErrorLowX        = new Double_t[nBins-notusedBins];
        Double_t *newnewBinningErrorHighY       = new Double_t[nBins-notusedBins];
        Double_t *newnewBinningErrorHighX       = new Double_t[nBins-notusedBins];


        for(Int_t graphBin = 0; graphBin<nBins-notusedBins; graphBin++){
            newnewBinningY[graphBin]            = newBinningY[graphBin];
            newnewBinningX[graphBin]            = newBinningX[graphBin];
            newnewBinningErrorLowY[graphBin]    = newBinningErrorLowY[graphBin];
            newnewBinningErrorLowX[graphBin]    = newBinningErrorLowX[graphBin];
            newnewBinningErrorHighY[graphBin]   = newBinningErrorHighY[graphBin];
            newnewBinningErrorHighX[graphBin]   = newBinningErrorHighX[graphBin];
        }

        TGraphAsymmErrors* returnGraph          = new TGraphAsymmErrors(nBins-notusedBins, newnewBinningX, newnewBinningY, newnewBinningErrorLowX, newnewBinningErrorHighX,
                                                                        newnewBinningErrorLowY, newnewBinningErrorHighY );

        delete graphAsym;
        delete[] newBinningY;
        delete[] newBinningX;
        delete[] newBinningErrorLowY;
        delete[] newBinningErrorLowX;
        delete[] newBinningErrorHighY;
        delete[] newBinningErrorHighX;
        delete[] newnewBinningY;
        delete[] newnewBinningX;
        delete[] newnewBinningErrorLowY;
        delete[] newnewBinningErrorLowX;
        delete[] newnewBinningErrorHighY;
        delete[] newnewBinningErrorHighX;
        return returnGraph;
    }

    // ****************************************************************************************************************
    // ****************************************************************************************************************
    // ****************************************************************************************************************
    TGraphErrors* RebinNLOGraph(TGraphErrors* graphV1){
        TGraphErrors *graph         = (TGraphErrors*) graphV1->Clone("graph");

        Double_t* xValue            = graph->GetX();
        Double_t* yValue            = graph->GetY();
        Int_t nPoints               = graph->GetN();

        Int_t newBinning            = nPoints*0.5;
        Double_t *newBinningY       = new Double_t[newBinning];
        Double_t *newBinningX       = new Double_t[newBinning];
        Double_t *newBinningErrorY  = new Double_t[newBinning];
        Double_t *newBinningErrorX  = new Double_t[newBinning];
        Int_t newBin                = 0;

        for(Int_t i = 2; i<nPoints; i=i+2){
            Double_t yValueOne          = 0;
            Double_t yValueTwo          = 0;
            Double_t yValueThree        = 0;
            Double_t yValueNew          = 0;
            Double_t xValueOne          = 0;
            Double_t xValueTwo          = 0;
            Double_t xValueThree        = 0;
            Double_t xValueNew          = 0;

            yValueOne                   = yValue[i-1];
            yValueTwo                   = yValue[i];
            yValueThree                 = yValue[i+1];

            xValueOne                   = xValue[i-1];
            xValueTwo                   = xValue[i];
            xValueThree                 = xValue[i+1];

            yValueNew                   = (yValueOne + yValueTwo+yValueThree)/3;
            xValueNew                   = (xValueOne + xValueTwo+xValueThree)/3;

            //   cout<<i<<"    "<<xValueOne<<"  "<<xValueTwo<<"   "<<xValueThree<<"   "<<  xValueNew<<endl;
            newBinningY[newBin]         = yValueNew;
            newBinningX[newBin]         = xValueNew;
            newBinningErrorY[newBin]    = 1.;
            newBinningErrorX[newBin]    = xValueTwo-xValueOne;
            //cout<<-(xValueOne-xValueTwo)*0.5<<endl;

            newBin++;
        }

        TGraphErrors* returnGraph       = new TGraphErrors(newBinning,newBinningX,newBinningY,newBinningErrorX,newBinningErrorY);

        delete graph;
        return returnGraph;
    }

    // ****************************************************************************************************************
    // ********************* Converts TGraphErrors to TH1D* ***********************************************************
    // ****************************************************************************************************************
    TH1D *GraphToHist(  TGraphErrors *graph,
                        Int_t maxPt         = 50,
                        TString name        = ""
                    ){

        //     graph->Print();
        Double_t* xValue        = graph->GetX();
        Double_t* yValue        = graph->GetY();
        Double_t* Ex            = graph->GetEX();
        Int_t  nPoints          = graph->GetN();
        Int_t maxPoints         = 0;

        for(Int_t i = 0; i<nPoints; i++){
            if(xValue[i]<=maxPt) maxPoints++;
        }

        Double_t *newBinningX   = new Double_t[maxPoints];
        for(Int_t i = 0;i<maxPoints;i++)
            newBinningX[i]      = xValue[i]-Ex[i];

        TH1D *hist              = new TH1D(name,"",maxPoints-1,newBinningX);

        for(Int_t i = 1;i<maxPoints;i++){
    //         cout << i << "\t pt: " << xValue[i] << "\t y: "<<yValue[i] << endl;
            hist->SetBinContent(i,yValue[i-1]);
        }
        return hist;
    }

    // ****************************************************************************************************************
    // ****************************************************************************************************************
    // ****************************************************************************************************************
    TH1D *GraphAsymErrorsToHist(    TGraphAsymmErrors *graph,
                                    Int_t maxPt                 = 50,
                                    TString name                = ""
                            ){
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

    // ****************************************************************************************************************
    // ****************************************************************************************************************
    // ****************************************************************************************************************
    TH1D *GraphAsymErrorsToHist_withErrors( TGraphAsymmErrors *graph,
                                            TString name                = ""
                                        ){
        Double_t* xValue        =  graph->GetX();
        Double_t* yValue        = graph->GetY();
        Double_t* Exhigh        = graph->GetEXhigh();
        Double_t* Exlow         = graph->GetEXlow();
        Double_t* Eyhigh        = graph->GetEYhigh();
        Double_t* Eylow         = graph->GetEYlow();
        Int_t nPoints           = graph->GetN();

        Double_t *newBinningX   = new Double_t[nPoints+1];
        for(Int_t i = 0;i<nPoints;i++)
            newBinningX[i]      = xValue[i]-Exlow[i];
            newBinningX[nPoints] = xValue[nPoints-1]+Exhigh[nPoints-1];
        TH1D *hist              = new TH1D(name,"",nPoints,newBinningX);

        for(Int_t i = 1;i<=nPoints;i++){
        hist->SetBinContent(i,yValue[i-1]);
        if (Eyhigh[i-1]<Eylow[i-1]) Eyhigh[i-1]=Eylow[i-1];
        hist->SetBinError(i,Eyhigh[i-1]);
        }

        return hist;
    }

    // ****************************************************************************************************************
    // ****************************************************************************************************************
    // ****************************************************************************************************************
    TH1D* RebinTH1D(    TH1D* histoV2,
                        TH1D* histoBinningV2,
                        Bool_t deltaPt          = kFALSE
                ){

        Double_t binArray[60];

        TH1D *histo             = (TH1D*) histoV2->Clone(histoV2->GetName());
        TH1D *histoBinning      = (TH1D*) histoBinningV2->Clone(histoBinningV2->GetName());

        Int_t nBinsPi0          = histoBinning->GetNbinsX();
        Bool_t sameStartBin     = kFALSE;
        Int_t missingStartBins  = 0;

        for(Int_t bin = 0; bin<nBinsPi0; bin++){
            if(!sameStartBin && (histoBinning->GetBinLowEdge(bin+1) < histo->GetBinLowEdge(1))){
                missingStartBins++;
                continue;
            }
            else if(!sameStartBin && histoBinning->GetBinLowEdge(bin+1) >= histo->GetBinLowEdge(1)){
                sameStartBin    = kTRUE;
            }

            binArray[bin-missingStartBins]          = histoBinning->GetBinLowEdge(bin+1);
            if(bin == (nBinsPi0-1))
                binArray[bin+1-missingStartBins]    = histoBinning->GetXaxis()->GetBinUpEdge(bin+1);
        }

        nBinsPi0                = nBinsPi0 - missingStartBins;
        TH1D *Rebin             = new TH1D(Form("Rebin_%s",histoV2->GetName()),"",nBinsPi0,binArray);
        Rebin->Sumw2();
        TH1D *DeltaPt           = new TH1D(Form("DeltaPt_%s",histoV2->GetName()),"",nBinsPi0,binArray);
        DeltaPt->Sumw2();

        for(Int_t iPt = 0; iPt < nBinsPi0; iPt++){
            Int_t startBin      = histo->GetXaxis()->FindBin(binArray[iPt]+0.001);
            Int_t endBin        = histo->GetXaxis()->FindBin(binArray[iPt+1]-0.001);
            Int_t nMergedBins   = 1+endBin - startBin;
            Rebin->SetBinContent(iPt+1,nMergedBins);
            Rebin->SetBinError(iPt+1,0);
            Double_t diffPt     = binArray[iPt+1]-binArray[iPt];
            DeltaPt->SetBinContent(iPt+1,diffPt);
            DeltaPt->SetBinError(iPt+1,0);
        }

        histo                   = (TH1D*) histo->Rebin(nBinsPi0,histo->GetName(),binArray);
        if(!deltaPt) histo->Divide(Rebin);
        if(deltaPt) histo->Divide(DeltaPt);

        delete Rebin;
        delete DeltaPt;
        delete histoBinning;

        return histo;
    }

    // ****************************************************************************************************************
    // ****************************************************************************************************************
    // ****************************************************************************************************************
    TH1D* CalculateFitToFitRatio(   TF1 *funcA,
                                    TF1 *funcB
                                ){

        Double_t maxA           = funcA->GetXmax();
        Double_t maxB           = funcB->GetXmax();
        Double_t newMax         = 0;
        if(maxA>maxB) newMax    = maxB;
        else newMax             = maxA;

        Int_t nBinsX            = newMax/0.1;
        TH1D *returnHist        = new TH1D(Form("%s/%s",funcA->GetName(),funcB->GetName()),"",nBinsX,0,newMax);

        for(Int_t i = 1; i<nBinsX+1; i++){
            Double_t binCenter      = returnHist->GetBinCenter(i);
            Double_t newBinContent  = funcA->Eval(binCenter)/funcB->Eval(binCenter);
            returnHist->SetBinContent(i,newBinContent);
            returnHist->SetBinError(i,0);
        }
        return returnHist;
    }

    // ****************************************************************************************************************
    // ****************************************************************************************************************
    // ****************************************************************************************************************
    TH1D* CalculateWeightedAveragePCM( TH1D* histoConvOnfly,
                                    TH1D* histoConvOffline
                                    ){

        TH1D* returnHisto                   = (TH1D*)histoConvOnfly->Clone("returnHisto");
        TGraphErrors* graphStatErrOffline   = new TGraphErrors(histoConvOffline);
        TGraphErrors* graphStatErrOnfly     = new TGraphErrors(histoConvOnfly);

        Int_t nOffline                      = graphStatErrOffline->GetN();
        Double_t* xOffline                  = graphStatErrOffline->GetX();
        Double_t* yOffline                  = graphStatErrOffline->GetY();
        Double_t* eyStaOffline              = graphStatErrOffline->GetEY();

        for(Int_t i=0;i<nOffline;i++){
            cout<< "Offline:: "<< i<< "\t" <<xOffline[i]<< " "<<yOffline[i]<< " " <<  eyStaOffline[i]<< endl;
        }

        Int_t nOnfly = graphStatErrOnfly->GetN();
        Double_t* xOnfly                    = graphStatErrOnfly->GetX();
        Double_t* yOnfly                    = graphStatErrOnfly->GetY();
        Double_t* eyStaOnfly                = graphStatErrOnfly->GetEY();

        for(Int_t i=0;i<nOnfly;i++){
            cout<< "Onfly::"<< i<< "\t" << xOnfly[i]<< " "<<yOnfly[i]<< " " <<  eyStaOnfly[i]<< endl;
        }

        cout<<endl;

        Bool_t okOffline,okOnfly;
        for (Int_t i=0;i<nOnfly;i++){
            okOffline                       = kFALSE;
            okOnfly                         = kFALSE;

            if (  yOffline[i]!= 0.){
                okOffline                   = kTRUE;
            }

            if (  yOnfly[i]!= 0.){
                okOnfly                     = kTRUE;
            }

            if ( okOffline && okOnfly ){
                if (eyStaOffline[i]== 0. && eyStaOnfly[i]>0. ){
                    returnHisto->SetBinContent(i+1, yOnfly[i]);
                    returnHisto->SetBinError(i+1, eyStaOnfly[i]);
                    cout<< " Combined::"<< i << "\t empty Offline " << yOnfly[i] << "\t" << eyStaOnfly[i]<<  endl;
                } else if( eyStaOnfly[i]!=0. &&  eyStaOffline[i]!=0. && eyStaOnfly[i]>0. && eyStaOffline[i]>0.){
                    Double_t wOffline       = 1./eyStaOffline[i];
                    Double_t wOnfly         = 1./eyStaOnfly[i];
                    Double_t wSum           = wOnfly+wOffline;
                    Double_t xSectionComb   = (wOnfly*yOnfly[i] +  wOffline*yOffline[i])/ wSum;
                    Double_t xSectionCombErr= pow(1./2.* wOnfly/wSum* eyStaOnfly[i]*eyStaOnfly[i] + 1./2.* wOffline/wSum* eyStaOffline[i]*eyStaOffline[i],0.5);
                    cout<< " Combined::"<< i << "\t" <<xSectionComb<< " " << xSectionCombErr << " " << yOnfly[i]<< " "<< eyStaOnfly[i] << " "
                    << yOffline[i]<<" "<< eyStaOffline[i]<< endl;
                    returnHisto->SetBinContent(i+1, xSectionComb);
                    returnHisto->SetBinError(i+1, xSectionCombErr);
                } else {
                    returnHisto->SetBinContent(i+1, 0.);
                    returnHisto->SetBinError(i+1, 0.);
                    cout<< " Combined::"<< i << "\t empty" << endl;
                    cout << "here" << endl;
                }
            } else {
                returnHisto->SetBinContent(i+1, 0.);
                returnHisto->SetBinError(i+1, 0.);
                cout<< " Combined::"<< i << "\t empty" << endl;
            }
        }
        return returnHisto;
    }

    // ****************************************************************************************************************
    // ****************************************************************************************************************
    // ****************************************************************************************************************
    void CalculateFitResults(   TF1* fitstat,
                                TF1* fitsys,
                                Double_t* fitresults,
                                TString fitName         = "",
                                Double_t sigma          = 1
                            ){
        Int_t nPar                  = fitstat->GetNpar();
        cout << nPar << endl;
        for (Int_t i =0; i < nPar; i++){
            if (fitName.CompareTo("Levy") == 0 && i == 0){
                fitresults[i*3]     = fitsys->GetParameter(i)/sigma;
                fitresults[i*3+1]   = fitstat->GetParError(i)/sigma;
                fitresults[i*3+2]   = TMath::Sqrt(TMath::Abs(pow((fitsys->GetParError(i)/sigma),2)-pow((fitstat->GetParError(i)/sigma),2)));
            } else {
                fitresults[i*3]     = fitsys->GetParameter(i);
                fitresults[i*3+1]   = fitstat->GetParError(i);
                fitresults[i*3+2]   = TMath::Sqrt(TMath::Abs(pow((fitsys->GetParError(i)),2)-pow((fitstat->GetParError(i)),2)));
            }
            cout << fitresults[i*3] << "\t" << fitresults[i*3+1] << "\t" << fitresults[i*3+2] << endl;
        }
        return;
    }

    // ****************************************************************************************************************
    // ****************************************************************************************************************
    // ****************************************************************************************************************
    void ReadOutFromSysErrVector(   Double_t* readoutVector,
                                    Double_t* fillInVector,
                                    Double_t* fillInVectorError,
                                    Int_t cutNr,
                                    Int_t nrCutStudies,
                                    Int_t offsetAtEnd,
                                    Int_t numberOfLines,
                                    Double_t decisionBoundary
                                ){
        Double_t dummyValue;
        Double_t dummyError;
        Int_t l                         = cutNr*3-2;
        Int_t totalLineEntries          = 1+ nrCutStudies*3+offsetAtEnd;
        for (Int_t i = 0; i < numberOfLines; i++){
            dummyValue                  = readoutVector[l];
            dummyError                  = readoutVector[l+1];
            cout << "readout \t " << l << "\t" << dummyValue << "\t" << dummyError << endl;
            if (dummyError != 0 &&TMath::Abs(dummyValue)/dummyError > decisionBoundary){
                fillInVector[i]         = dummyValue/TMath::Sqrt(2);
                fillInVectorError[i]    = fillInVector[i]*0.001;
            } else {
                fillInVector[i]         = 0;
                fillInVectorError[i]    = 0;
            }
            l                           = l + totalLineEntries;
        }
        return;
    }

    // ****************************************************************************************************************
    // ****************************************************************************************************************
    // ****************************************************************************************************************
    void CalculateMeanSysErr(   Double_t* fillInVector,
                                Double_t* fillInVectorError,
                                Double_t* posError,
                                Double_t* negError,
                                Int_t numberOfLines
                            ){
        Double_t dummyValue;
        Double_t dummyError;
        for (Int_t i = 0; i < numberOfLines; i++){
            dummyValue              = (posError[i] + TMath::Abs(negError[i]))/2;
            cout << posError[i] << "\t" << negError[i] << "\t" << dummyValue << endl;
            dummyError              = dummyValue*0.001;
            fillInVector[i]         = dummyValue;
            fillInVectorError[i]    = dummyError;
        }
        return;
    }

    // ****************************************************************************************************************
    // ****************************************************************************************************************
    // ****************************************************************************************************************
    void CorrectSystematicErrorsWithMean(   Double_t* oldErrorVector,
                                            Double_t* oldErrorVectorsErrors,
                                            Double_t* newErrorVector,
                                            Double_t* newErrorVectorError,
                                            Int_t numberOfPtBins
                                        ){
        if (numberOfPtBins > 2) {
            for (Int_t i = 0; i < numberOfPtBins; i++){
                if (i == 0){ // bin 1
                    if (oldErrorVector[i] == 0 && oldErrorVector[i+1] == 0){
                        newErrorVector[i]           = oldErrorVector[i+2];
                        newErrorVectorError[i]      = oldErrorVectorsErrors[i+2];
                    } else if (oldErrorVector[i] == 0){
                        newErrorVector[i]           = oldErrorVector[i+1];
                        newErrorVectorError[i]      = oldErrorVectorsErrors[i+1];
                    } else {
                        newErrorVector[i]           = oldErrorVector[i];
                        newErrorVectorError[i]      = newErrorVector[i]*0.001;
                    }
                } else if (i == 1){ // bin 2
                    if (oldErrorVector[i] == 0 && oldErrorVector[i+1] != 0){
                        newErrorVector[i]           = oldErrorVector[i-1]*0.2 + oldErrorVector[i+1]*0.8;
                        newErrorVectorError[i]      = newErrorVector[i]*0.001;
                        if (oldErrorVector[i-1] == 0 && newErrorVector[i-1] != 0){
                            newErrorVector[i]       = newErrorVector[i-1]*0.2 + oldErrorVector[i+1]*0.8 ;
                            newErrorVectorError[i]  = newErrorVector[i]*0.001;
                        }
                    } else if (oldErrorVector[i] == 0 && oldErrorVector[i+1] == 0){
                        newErrorVector[i]           = oldErrorVector[i-1]*0.2 + oldErrorVector[i+2]*0.8;
                        newErrorVectorError[i]      = newErrorVector[i]*0.001;
                        if (oldErrorVector[i-1] == 0 && newErrorVector[i-1] != 0){
                            newErrorVector[i]       = (newErrorVector[i-1] + oldErrorVector[i+1])/2. ;
                            newErrorVectorError[i]  = newErrorVector[i]*0.001;
                        }
                    } else {
                        newErrorVector[i]           = oldErrorVector[i];
                        newErrorVectorError[i]      = newErrorVector[i]*0.001;
                    }
                }else if (i == numberOfPtBins - 2 ){ // previous to the last bin
                    if (oldErrorVector[i] == 0 && oldErrorVector[i+1] != 0){
                        newErrorVector[i]           = oldErrorVector[i-1]*0.8 + oldErrorVector[i+1]*0.2;
                        newErrorVectorError[i]      = newErrorVector[i]*0.001;
                        if (oldErrorVector[i-1] == 0 && newErrorVector[i-1] != 0){
                            newErrorVector[i]       = newErrorVector[i-1]*0.8 + oldErrorVector[i+1]*0.2 ;
                            newErrorVectorError[i]  = newErrorVector[i]*0.001;
                        }
                    } else if (oldErrorVector[i] == 0 && oldErrorVector[i+1] == 0){
                        newErrorVector[i]           = (oldErrorVector[i-1] + oldErrorVector[i+1])/2. ;
                        newErrorVectorError[i]      = newErrorVector[i]*0.001;
                        if (oldErrorVector[i-1] == 0 && newErrorVector[i-1] != 0){
                            newErrorVector[i]       = (newErrorVector[i-1] + oldErrorVector[i+1])/2. ;
                            newErrorVectorError[i]  = newErrorVector[i]*0.001;
                        }
                    } else {
                        newErrorVector[i]           = oldErrorVector[i];
                        newErrorVectorError[i]      = newErrorVector[i]*0.001;
                    }
                } else if (i == numberOfPtBins - 3){ //3rd form end
                    if (oldErrorVector[i] == 0 && oldErrorVector[i+1] != 0){
                        if (TMath::Abs(oldErrorVector[i-1]) < TMath::Abs(oldErrorVector[i+1])){
                            newErrorVector[i]       = oldErrorVector[i-1];
                            newErrorVectorError[i]  = newErrorVector[i]*0.001;
                        } else {
                            newErrorVector[i]       = oldErrorVector[i+1];
                            newErrorVectorError[i]  = newErrorVector[i]*0.001;
                        }
                    } else if (oldErrorVector[i] == 0 && oldErrorVector[i+1] == 0 && oldErrorVector[i+2] != 0){
                        if (TMath::Abs(oldErrorVector[i-1]) < TMath::Abs(oldErrorVector[i+2])){
                            newErrorVector[i]       = oldErrorVector[i-1];
                            newErrorVectorError[i]  = newErrorVector[i]*0.001;
                        } else {
                            newErrorVector[i]       = oldErrorVector[i+2];
                            newErrorVectorError[i]  = newErrorVector[i]*0.001;
                        }
                    } else {
                        newErrorVector[i]           = oldErrorVector[i];
                        newErrorVectorError[i]      = oldErrorVectorsErrors[i];
                    }
                } else if (i == numberOfPtBins - 1 ){ //last
                    if (oldErrorVector[i] == 0 && oldErrorVector[i-1] != 0){
                        newErrorVector[i]           = oldErrorVector[i-1];
                        newErrorVectorError[i]      = oldErrorVectorsErrors[i-1];
                    } else if (oldErrorVector[i] == 0 && oldErrorVector[i-1] == 0 && newErrorVector[i-1] != 0){
                        newErrorVector[i]           = newErrorVector[i-1];
                        newErrorVectorError[i]      = newErrorVectorError[i-1];
                    } else {
                        newErrorVector[i]           = oldErrorVector[i];
                        newErrorVectorError[i]      = newErrorVector[i]*0.001;
                    }
                } else {     //every where else
                    if (oldErrorVector[i] == 0 && oldErrorVector[i+1] == 0 && i+2 != numberOfPtBins){
                        newErrorVector[i]           = (oldErrorVector[i-1] + oldErrorVector[i+2])/2. ;
                        newErrorVectorError[i]      = newErrorVector[i]*0.001;
                        if (oldErrorVector[i-1] == 0 && newErrorVector[i-1] != 0){
                            newErrorVector[i]       = (newErrorVector[i-1] + oldErrorVector[i+2])/2. ;
                            newErrorVectorError[i]  = newErrorVector[i]*0.001;
                        }
                    } else if (oldErrorVector[i] == 0){
                        newErrorVector[i]           = (oldErrorVector[i-1] + oldErrorVector[i+1])/2. ;
                        newErrorVectorError[i]      = newErrorVector[i]*0.001;
                        if (oldErrorVector[i-1] == 0 && newErrorVector[i-1] != 0){
                            newErrorVector[i]       = (newErrorVector[i-1] + oldErrorVector[i+1])/2. ;
                            newErrorVectorError[i]  = newErrorVector[i]*0.001;
                        }
                    } else {
                        newErrorVector[i]           = oldErrorVector[i];
                        newErrorVectorError[i]      = newErrorVector[i]*0.001;
                    }
                }
            }
        } else {
            newErrorVector[0]                       = (oldErrorVector[0] + oldErrorVector[1])/2.;
            newErrorVector[1]                       = newErrorVector[0];
            newErrorVectorError[0]                  = newErrorVector[0]*0.001;
            newErrorVectorError[1]                  = newErrorVectorError[0];
        }
    }

    // ****************************************************************************************************************
    // ****************************************************************************************************************
    // ****************************************************************************************************************
    void ProduceGraphAsymmWithoutXErrors(TGraphAsymmErrors* inputgraph){
        Int_t n                 = inputgraph->GetN();
        Double_t* xValue        = inputgraph->GetX();
        Double_t* xErrorHigh    = inputgraph->GetEXhigh();
        Double_t* xErrorLow     = inputgraph->GetEXlow();
        Double_t* yValue        = inputgraph->GetY();
        Double_t* yErrorLow     = inputgraph->GetEYlow();
        Double_t* yErrorHigh    = inputgraph->GetEYhigh();
        for (Int_t i= 0; i < n; i++){
            xErrorHigh[i]       = 0.;
            xErrorLow[i]        = 0.;
        }
        inputgraph              = new TGraphAsymmErrors(n,xValue,yValue,xErrorLow,xErrorHigh,yErrorLow,yErrorHigh);
    //     inputgraph->Print();
    }

    // ****************************************************************************************************************
    // ****************************************************************************************************************
    // ****************************************************************************************************************
    void ProduceGraphPartialXErrors(TGraphErrors* inputgraph, Double_t part){  // kk
        Int_t n             = inputgraph->GetN();
        Double_t* xValue    = inputgraph->GetX();
        Double_t* xError    = inputgraph->GetEX();
        Double_t* yValue    = inputgraph->GetY();
        Double_t* yError    = inputgraph->GetEY();
        for (Int_t i= 0; i < n; i++){
            xError[i]       = part*xError[i];
            xError[i]       = part*xError[i];
        }
        inputgraph          = new TGraphErrors(n,xValue,yValue,xError,yError);
    //     inputgraph->Print();
    }

    // ****************************************************************************************************************
    // ****************************************************************************************************************
    // ****************************************************************************************************************
    void ProduceGraphAsymmPartialXErrors(TGraphAsymmErrors* inputgraph, Double_t part){  // kk
        Int_t n                 = inputgraph->GetN();
        Double_t* xValue        = inputgraph->GetX();
        Double_t* xErrorHigh    = inputgraph->GetEXhigh();
        Double_t* xErrorLow     = inputgraph->GetEXlow();
        Double_t* yValue        = inputgraph->GetY();
        Double_t* yErrorLow     = inputgraph->GetEYlow();
        Double_t* yErrorHigh    = inputgraph->GetEYhigh();
        for (Int_t i= 0; i < n; i++){
            xErrorHigh[i]       = part*xErrorHigh[i];
            xErrorLow[i]        = part*xErrorLow[i];
        }
        inputgraph              = new TGraphAsymmErrors(n,xValue,yValue,xErrorLow,xErrorHigh,yErrorLow,yErrorHigh);
    //     inputgraph->Print();
    }

    // ****************************************************************************************************************
    // ****************************************************************************************************************
    // ****************************************************************************************************************
    void ProduceGraphFixedXErrors(TGraphErrors* inputgraph, Double_t widthXError){
        Int_t n             = inputgraph->GetN();
        Double_t* xValue    = inputgraph->GetX();
        Double_t* xError    = inputgraph->GetEX();
        Double_t* yValue    = inputgraph->GetY();
        Double_t* yError    = inputgraph->GetEY();
        for (Int_t i= 0; i < n; i++){
            xError[i]       = widthXError/2;
            xError[i]       = widthXError/2;
        }
        inputgraph          = new TGraphErrors(n,xValue,yValue,xError,yError);
    //     inputgraph->Print();
    }

    // ****************************************************************************************************************
    // ****************************************************************************************************************
    // ****************************************************************************************************************
    void ProduceGraphAsymmFixedXErrors(TGraphAsymmErrors* inputgraph, Double_t widthXError){
        Int_t n                 = inputgraph->GetN();
        Double_t* xValue        = inputgraph->GetX();
        Double_t* xErrorHigh    = inputgraph->GetEXhigh();
        Double_t* xErrorLow     = inputgraph->GetEXlow();
        Double_t* yValue        = inputgraph->GetY();
        Double_t* yErrorLow     = inputgraph->GetEYlow();
        Double_t* yErrorHigh    = inputgraph->GetEYhigh();
        for (Int_t i= 0; i < n; i++){
            xErrorHigh[i]       = widthXError/2;
            xErrorLow[i]        = widthXError/2;
        }
        inputgraph              = new TGraphAsymmErrors(n,xValue,yValue,xErrorLow,xErrorHigh,yErrorLow,yErrorHigh);
    //     inputgraph->Print();
    }

    // ****************************************************************************************************************
    // ****************************************************************************************************************
    // ****************************************************************************************************************
    void ProduceGraphAsymmWithoutXYErrors(TGraphAsymmErrors* inputgraph){  // kk
        Int_t n                 = inputgraph->GetN();
        Double_t* xValue        = inputgraph->GetX();
        Double_t* xErrorHigh    = inputgraph->GetEXhigh();
        Double_t* xErrorLow     = inputgraph->GetEXlow();
        Double_t* yValue        = inputgraph->GetY();
        Double_t* yErrorLow     = inputgraph->GetEYlow();
        Double_t* yErrorHigh    = inputgraph->GetEYhigh();
        for (Int_t i= 0; i < n; i++){
            xErrorHigh[i]       = 0.;
            xErrorLow[i]        = 0.;
            yErrorHigh[i]       = 0.;
            yErrorLow[i]        = 0.;
        }
        inputgraph              = new TGraphAsymmErrors(n,xValue,yValue,xErrorLow,xErrorHigh,yErrorLow,yErrorHigh);
    //     inputgraph->Print();
    }

    // ****************************************************************************************************************
    // ****************************************************************************************************************
    // ****************************************************************************************************************
    void ProduceGraphErrWithoutXErrors(TGraphErrors* inputgraph){
        Int_t n                 = inputgraph->GetN();
        Double_t* xValue        = inputgraph->GetX();
        Double_t* xErrorHigh    = inputgraph->GetEX();
        Double_t* yValue        = inputgraph->GetY();
        Double_t* yErrorHigh    = inputgraph->GetEY();
        for (Int_t i= 0; i < n; i++){
            xErrorHigh[i]       = 0.;
        }
        inputgraph              = new TGraphErrors(n,xValue,yValue,xErrorHigh,yErrorHigh);
    //     inputgraph->Print();
    }

    // ****************************************************************************************************************
    // ****************************************************************************************************************
    // ****************************************************************************************************************
    void ProduceGraphErrDisplacedX(TGraphErrors* inputgraph,Double_t displaceX){
        Int_t n                 = inputgraph->GetN();
        Double_t* xValue        = inputgraph->GetX();
        Double_t* xErrorHigh    = inputgraph->GetEX();
        Double_t* yValue        = inputgraph->GetY();
        Double_t* yErrorHigh    = inputgraph->GetEY();
        for (Int_t i= 0; i < n; i++){
            xValue[i]           = xValue[i]+displaceX;
            xErrorHigh[i]       = 0.;
        }
        inputgraph              = new TGraphErrors(n,xValue,yValue,xErrorHigh,yErrorHigh);
    //     inputgraph->Print();
    }

    // ****************************************************************************************************************
    // ****************************************************************************************************************
    // ****************************************************************************************************************
    TGraphAsymmErrors* Add2TGraphAsymmErrorsSameBinning(TGraphAsymmErrors* inputgraph1,TGraphAsymmErrors* inputgraph2){
        Double_t* xValue1       = inputgraph1->GetX();
        Double_t* yValue1       = inputgraph1->GetY();
        Double_t* xErrorLow1    = inputgraph1->GetEXlow();
        Double_t* xErrorHigh1   = inputgraph1->GetEXhigh();
        Double_t* yErrorLow1    = inputgraph1->GetEYlow();
        Double_t* yErrorHigh1   = inputgraph1->GetEYhigh();

        Double_t* yValue2       = inputgraph2->GetY();
        Double_t* yErrorLow2    = inputgraph2->GetEYlow();
        Double_t* yErrorHigh2   = inputgraph2->GetEYhigh();

        Double_t* yValue        = inputgraph2->GetY();
        Double_t* yErrorLow     = inputgraph2->GetEYlow();
        Double_t* yErrorHigh    = inputgraph2->GetEYhigh();

        Int_t nPoints           = inputgraph1->GetN();
        for (Int_t i = 0; i < nPoints; i++){
            yValue[i]           = yValue1[i]+ yValue2[i];
            yErrorLow[i]        = TMath::Sqrt(TMath::Power(yErrorLow1[i],2) + TMath::Power(yErrorLow2[i],2));
            yErrorHigh[i]       = TMath::Sqrt(TMath::Power(yErrorHigh1[i],2) + TMath::Power(yErrorHigh2[i],2));
        }
        TGraphAsymmErrors* returnGraph = new TGraphAsymmErrors(nPoints,xValue1,yValue,xErrorLow1,xErrorHigh1,yErrorLow,yErrorHigh);
        return returnGraph;
    }

    // ****************************************************************************************************************
    // ****************************************************************************************************************
    // ****************************************************************************************************************
    TGraphAsymmErrors* CorrectTGraphAsymmErrorsToBinCenter(TGraphAsymmErrors* inputgraph1){
        Double_t* xValue1           = inputgraph1->GetX();
        Double_t* yValue1           = inputgraph1->GetY();
        Double_t* xErrorLow1        = inputgraph1->GetEXlow();
        Double_t* xErrorHigh1       = inputgraph1->GetEXhigh();
        Double_t* yErrorLow1        = inputgraph1->GetEYlow();
        Double_t* yErrorHigh1       = inputgraph1->GetEYhigh();

        Double_t* yValueNew         = inputgraph1->GetY();
        Double_t* yErrorLowNew      = inputgraph1->GetEYlow();
        Double_t* yErrorHighNew     = inputgraph1->GetEYhigh();

        Int_t nPoints               = inputgraph1->GetN();
        for (Int_t i = 0; i < nPoints; i++){
            yValueNew[i]            = yValue1[i]/xValue1[i];
            yErrorLowNew[i]         = yErrorLow1[i]/xValue1[i];
            yErrorHighNew[i]        = yErrorHigh1[i]/xValue1[i];
        }
        TGraphAsymmErrors* returnGraph = new TGraphAsymmErrors(nPoints,xValue1,yValueNew,xErrorLow1,xErrorHigh1,yErrorLowNew,yErrorHighNew);
        return returnGraph;
    }

    // ****************************************************************************************************************
    // ****************************************************************************************************************
    // ****************************************************************************************************************
    TH1D* ShortChargedHadronHisto(TH1D* histIn) {
        const Int_t nPtBins         = 57;
        Double_t xBins[nPtBins+1]   = { 0. ,  0.05, 0.1,  0.15, 0.2,  0.25, 0.3,  0.35, 0.4,  0.45,
                                        0.5,  0.55, 0.6,  0.65, 0.7,  0.75, 0.8,  0.85, 0.9,  0.95,
                                        1.0,  1.1 , 1.2,  1.3 , 1.4,  1.5 , 1.6,  1.7 , 1.8,  1.9 ,
                                        2.0,  2.2 , 2.4,  2.6 , 2.8,  3.0 , 3.2,  3.4 , 3.6,  3.8 ,
                                        4.0,  4.5 , 5.0,  5.5 , 6.0,  6.5 , 7.0,  8.0 , 9.0,  10.0,
                                        11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0};

        TH1D* hist                  = new TH1D("hist", "", nPtBins, xBins);
        hist->SetTitle(histIn->GetTitle());
        hist->GetXaxis()->SetTitle(histIn->GetXaxis()->GetTitle());
        hist->GetYaxis()->SetTitle(histIn->GetYaxis()->GetTitle());

        const Double_t deltapt      = 0.0001;

        for(Int_t bin = 1; bin <= nPtBins; bin++) {

            // check bin size
            if(TMath::Abs(hist->GetXaxis()->GetBinLowEdge(bin) -
                histIn->GetXaxis()->GetBinLowEdge(bin)) > deltapt) {
                cout << "pt edge low does not agree!!!!!!!!!!!!!!!" << endl;
                return 0;
            }
            if(TMath::Abs(hist->GetXaxis()->GetBinUpEdge(bin) -
                histIn->GetXaxis()->GetBinUpEdge(bin)) > deltapt) {
                cout << "pt edge high does not agree!!!!!!!!!!!!!!!" << endl;
                return 0;
            }
            if(TMath::Abs(hist->GetXaxis()->GetBinCenter(bin) -
                histIn->GetXaxis()->GetBinCenter(bin)) > deltapt) {
                cout << "pt center does not agree!!!!!!!!!!!!!!!" << endl;
                return 0;
            }

            hist->SetBinContent(bin, histIn->GetBinContent(bin));
            cout << "\n converting n\n" << endl;
            cout << "pt in   " << histIn->GetXaxis()->GetBinCenter(bin) << "    bin content in  " << histIn->GetBinContent(bin) << endl;
            cout << " pt out   " << hist->GetXaxis()->GetBinCenter(bin) << "    bin content out   " << hist->GetBinContent(bin) << endl;
            hist->SetBinError(bin, histIn->GetBinError(bin));
        }

        return hist;
    }

    // ****************************************************************************************************************
    // ****************************************************************************************************************
    // ****************************************************************************************************************
    TH1D* ConvertChargedHadronHisto(    TH1D* histIn,
                                        TH1D* fraction
                                ) {
        const Int_t nPtBins         = 57;
        Double_t xBins[nPtBins+1]   = { 0. ,  0.05, 0.1,  0.15, 0.2,  0.25, 0.3,  0.35, 0.4,  0.45,
                                        0.5,  0.55, 0.6,  0.65, 0.7,  0.75, 0.8,  0.85, 0.9,  0.95,
                                        1.0,  1.1 , 1.2,  1.3 , 1.4,  1.5 , 1.6,  1.7 , 1.8,  1.9 ,
                                        2.0,  2.2 , 2.4,  2.6 , 2.8,  3.0 , 3.2,  3.4 , 3.6,  3.8 ,
                                        4.0,  4.5 , 5.0,  5.5 , 6.0,  6.5 , 7.0,  8.0 , 9.0,  10.0,
                                        11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0};

        TH1D* hist                  = new TH1D("hist", "", nPtBins, xBins);
        hist->SetTitle(histIn->GetTitle());
        hist->GetXaxis()->SetTitle(histIn->GetXaxis()->GetTitle());
        hist->GetYaxis()->SetTitle(histIn->GetYaxis()->GetTitle());

        const Double_t deltapt      = 0.0001;

        for(Int_t bin = 1; bin <= nPtBins; bin++) {

            // check bin size
            if(TMath::Abs(hist->GetXaxis()->GetBinLowEdge(bin) -
                histIn->GetXaxis()->GetBinLowEdge(bin)) > deltapt) {
                cout << "pt edge low does not agree!!!!!!!!!!!!!!!" << endl;
                return 0;
            }
            if(TMath::Abs(hist->GetXaxis()->GetBinUpEdge(bin) -
                histIn->GetXaxis()->GetBinUpEdge(bin)) > deltapt) {
                cout << "pt edge high does not agree!!!!!!!!!!!!!!!" << endl;
                return 0;
            }
            if(TMath::Abs(hist->GetXaxis()->GetBinCenter(bin) -
                histIn->GetXaxis()->GetBinCenter(bin)) > deltapt) {
                cout << "pt center does not agree!!!!!!!!!!!!!!!" << endl;
                return 0;
            }

            if(TMath::Abs(hist->GetXaxis()->GetBinLowEdge(bin) -
                fraction->GetXaxis()->GetBinLowEdge(bin)) > deltapt) {
                cout << "pt edge fraction low does not agree!!!!!!!!!!!!!!!" << endl;
                return 0;
            }
            if(TMath::Abs(hist->GetXaxis()->GetBinUpEdge(bin) -
                fraction->GetXaxis()->GetBinUpEdge(bin)) > deltapt) {
                cout << "pt edge fraction high does not agree!!!!!!!!!!!!!!!" << endl;
                return 0;
            }
            if(TMath::Abs(hist->GetXaxis()->GetBinCenter(bin) -
                fraction->GetXaxis()->GetBinCenter(bin)) > deltapt) {
                cout << "pt center fraction does not agree!!!!!!!!!!!!!!!" << endl;
                return 0;
            }

            Double_t content        = histIn->GetBinContent(bin)*fraction->GetBinContent(bin);
            Double_t error          = TMath::Sqrt(pow(histIn->GetBinError(bin)*fraction->GetBinContent(bin),2)+ pow(histIn->GetBinContent(bin)*fraction->GetBinError(bin),2) );

            hist->SetBinContent(bin, content);
            cout << "\n converting n\n" << endl;
            cout << "pt in   " << histIn->GetXaxis()->GetBinCenter(bin) << "    bin content in  " << histIn->GetBinContent(bin) << endl;
            cout << "pt frac   " << fraction->GetXaxis()->GetBinCenter(bin) << "    bin content fraction  " << fraction->GetBinContent(bin) << endl;
            cout << " pt out   " << hist->GetXaxis()->GetBinCenter(bin) << "    bin content out   " << hist->GetBinContent(bin) << endl;
            hist->SetBinError(bin, error);
        }

        return hist;
    }

    // ****************************************************************************************************************
    // ****************************************************************************************************************
    // ****************************************************************************************************************
    TF1 *BinShiftTH1D(  TH1D *UnshiftedYield,
                        TH1D **ShiftedYield,
                        TString mesonType           = "Pi0",
                        TString BinShiftType        = "L",
                        TString NameBinShiftHist    = "DummyFit",
                        Double_t minPtForFitsDummy  = 0.6,
                        Double_t* parameters        = 0x00
                    ){
        TF1 *CurrentFit             = new TF1();
        TH1D *CurrentHist           = (TH1D*)UnshiftedYield->Clone("");
        (*ShiftedYield)             = (TH1D*)CurrentHist->Clone(NameBinShiftHist);

        Int_t binNumber             = CurrentHist->GetXaxis()->GetNbins();
        Double_t globalRatioBS      = 0;
        Double_t testGlobalRatioBS  = 1;
        Int_t colorBS               = 1;
        Double_t maxPtPi0BS         = CurrentHist->GetXaxis()->GetBinUpEdge(CurrentHist->GetNbinsX());

        if(BinShiftType.BeginsWith("h") || BinShiftType.BeginsWith("H")){
            CurrentFit              = FitObject("h","fitBinShifting",mesonType.Data(),CurrentHist,minPtForFitsDummy,maxPtPi0BS,parameters,"QNRME+");
        } else if(BinShiftType.BeginsWith("l") || BinShiftType.BeginsWith("L")){
            CurrentFit              = FitObject("l","fitBinShifting",mesonType.Data(),CurrentHist,minPtForFitsDummy,maxPtPi0BS,parameters,"QNRME+");
        } else if(BinShiftType.BeginsWith("p") || BinShiftType.BeginsWith("P")){
            CurrentFit              = FitObject("p","fitBinShifting",mesonType.Data(),CurrentHist,minPtForFitsDummy,maxPtPi0BS,parameters,"QNRME+");
        } else if(BinShiftType.BeginsWith("rad") || BinShiftType.BeginsWith("RAD")){
            CurrentFit              = FitObject("rad","fitBinShifting",mesonType.Data(),CurrentHist,minPtForFitsDummy,maxPtPi0BS,parameters,"QNRME+");
        } else if(BinShiftType.BeginsWith("qcd") || BinShiftType.BeginsWith("QCD")){
            CurrentFit              = FitObject("qcd","fitBinShifting",mesonType.Data(),CurrentHist,minPtForFitsDummy,maxPtPi0BS,parameters,"QNRME+");
        } else if(BinShiftType.BeginsWith("xqcd") || BinShiftType.BeginsWith("XQCD")){
            CurrentFit              = FitObject("xqcd","fitBinShifting",mesonType.Data(),CurrentHist,minPtForFitsDummy,maxPtPi0BS,parameters,"QNRME+");
        } else if(BinShiftType.BeginsWith("tcm") || BinShiftType.BeginsWith("TCM")){
            CurrentFit              = FitObject("tcm","fitBinShifting",mesonType.Data(),CurrentHist,minPtForFitsDummy,maxPtPi0BS,parameters,"QNRME+");
        }
        CurrentFit->SetRange(minPtForFitsDummy,maxPtPi0BS);

        while(globalRatioBS != testGlobalRatioBS){
            if(colorBS == 200) break;
            testGlobalRatioBS       = globalRatioBS;
            globalRatioBS           = 0.;
            colorBS++;
            // fit acc+eff corrected yield with Hagedorn function
            (*ShiftedYield)->Fit(CurrentFit,"QRM0");

            // apply bin shift correction
            for (Int_t ib=2; ib<=binNumber; ib++) {
                Double_t ptMin          = CurrentHist->GetBinLowEdge(ib);
                Double_t ptMax          = ptMin + CurrentHist->GetBinWidth(ib);
                // the bin shift affected value of the fit function in current bin
                Double_t shiftedValue   = CurrentFit->Integral(ptMin,ptMax) / (ptMax - ptMin);
                // the correct value at the bin center
                Double_t trueValue      = CurrentFit->Eval((ptMax + ptMin)/2.);


                // the bin shift correction factor
                Double_t ratio          = shiftedValue / trueValue;

                (*ShiftedYield)->SetBinContent(ib,  CurrentHist->GetBinContent(ib) / ratio);
                (*ShiftedYield)->SetBinError(ib,  CurrentHist->GetBinError(ib) / ratio);
                globalRatioBS           = globalRatioBS + ratio;
            }

            globalRatioBS           = globalRatioBS/(binNumber-1);

            (*ShiftedYield)->SetMarkerStyle(24);
            (*ShiftedYield)->SetMarkerSize(0.9);
            (*ShiftedYield)->SetMarkerColor(kBlack);

            DrawGammaSetMarkerTF1( CurrentFit, 1, 0.4, kBlue-4);

            Double_t parameter[10];
            CurrentFit->GetParameters(parameter);
            CurrentFit->SetParameter(0,parameter[0]);
            CurrentFit->SetParameter(1,parameter[1]);
            CurrentFit->SetParameter(2,parameter[2]);
            CurrentFit->SetParameter(3,parameter[3]);
            CurrentFit->SetParameter(4,parameter[4]);
            CurrentFit->SetParameter(5,parameter[5]);
            CurrentFit->SetParameter(6,parameter[6]);
            CurrentFit->SetParameter(7,parameter[7]);
            CurrentFit->SetParameter(8,parameter[8]);
            CurrentFit->SetParameter(9,parameter[9]);

            cout<<colorBS<<" ";
        }
    //     cout << WriteParameterToFile(CurrentFit)<< endl;
        return CurrentFit;

    }

    // ****************************************************************************************************************
    // ****************************************************************************************************************
    // ****************************************************************************************************************
    TF1 *ApplyYShift(   TGraphAsymmErrors *UnshiftedYield,
                        TGraphAsymmErrors **ShiftedYield,
                        TString BinShiftType,
                        TString NameBinShiftHist,
                        Double_t minPtForFitsDummy          = 0.6,
                        Double_t* parameters                = 0x00,
                        Double_t accuracy                   = 0.0001,
                        Bool_t incParlimits                 = kFALSE,
                        TString mesonType                   =  "Pi0"
                    ){
    //     cout << "entered Y shift" << endl;
        TF1 *CurrentFit                     = new TF1();
        TGraphAsymmErrors *CurrentGraph     = (TGraphAsymmErrors*)UnshiftedYield->Clone("");
        (*ShiftedYield)                     = (TGraphAsymmErrors*)CurrentGraph->Clone(NameBinShiftHist);

        Int_t binNumber                     = CurrentGraph->GetN();
        Double_t* xvalue                    = UnshiftedYield->GetX();
        Double_t* xvalueErr                 = UnshiftedYield->GetEXlow();
        Double_t* yValue                    = UnshiftedYield->GetY();
        Double_t* yValueErrlow              = UnshiftedYield->GetEYlow();
    //     Double_t *yValueErrhigh          = UnshiftedYield->GetEYhigh();
    //     Double_t xbins[binNumber+1];

        Double_t globalRatioBS              = 0;
        Double_t testGlobalRatioBS          = 1;
        Int_t colorBS                       = 1;
        Double_t maxPtPi0BS                 = xvalue[binNumber] + xvalueErr[binNumber];

        if(BinShiftType.BeginsWith("h") || BinShiftType.BeginsWith("H")){
            CurrentFit                      = FitObject("h","fitBinShiftingPi0",mesonType.Data(),CurrentGraph,minPtForFitsDummy,maxPtPi0BS,parameters);
        } else if(BinShiftType.BeginsWith("l") || BinShiftType.BeginsWith("L")){
            CurrentFit                      = FitObject("l","fitBinShiftingPi0",mesonType.Data(),CurrentGraph,minPtForFitsDummy,maxPtPi0BS,parameters);
        } else if(BinShiftType.BeginsWith("p") || BinShiftType.BeginsWith("P")){
            CurrentFit                      = FitObject("p","fitBinShiftingPi0",mesonType.Data(),CurrentGraph,minPtForFitsDummy,maxPtPi0BS,parameters);
        } else if(BinShiftType.BeginsWith("rad") || BinShiftType.BeginsWith("RAD")){
            if (!incParlimits){
                CurrentFit                  = FitObject("rad","fitBinShiftingPi0",mesonType.Data(),CurrentGraph,minPtForFitsDummy,maxPtPi0BS,parameters);
            } else {
                CurrentFit                  = FitObject("rad","fitBinShiftingPi0",mesonType.Data(),CurrentGraph,minPtForFitsDummy,maxPtPi0BS);
                SetParametersLimitsForFit(CurrentFit, 5, parameters);
            }
        } else if(BinShiftType.BeginsWith("qcd") || BinShiftType.BeginsWith("QCD")){
            CurrentFit                      = FitObject("qcd","fitBinShiftingPi0",mesonType.Data(),CurrentGraph,minPtForFitsDummy,maxPtPi0BS,parameters);
        } else if(BinShiftType.BeginsWith("xqcd") || BinShiftType.BeginsWith("XQCD")){
            CurrentFit                      = FitObject("xqcd","fitBinShiftingPi0",mesonType.Data(),CurrentGraph,minPtForFitsDummy,maxPtPi0BS,parameters);
        } else if(BinShiftType.BeginsWith("tcm") || BinShiftType.BeginsWith("TCM")){
            CurrentFit                      = FitObject("tcm","fitBinShiftingPi0",mesonType.Data(),CurrentGraph,minPtForFitsDummy,maxPtPi0BS,parameters);
        }

        CurrentFit->SetRange(minPtForFitsDummy,maxPtPi0BS);


        Double_t ratio[binNumber+1];
        while(TMath::Abs(globalRatioBS - testGlobalRatioBS) > accuracy){
            if(colorBS == 200) break;
            testGlobalRatioBS               = globalRatioBS;
            globalRatioBS                   = 0.;
            colorBS++;
            (*ShiftedYield)->Fit(CurrentFit,"NRMEX0+");

    //         cout << WriteParameterToFile(CurrentFit) << endl;
            Double_t *xPoint                = (*ShiftedYield)->GetX();
            Double_t *yPoint                = (*ShiftedYield)->GetY();

            Double_t *errorXlow             = (*ShiftedYield)->GetEXlow();
    //         Double_t *errorXhigh         = (*ShiftedYield)->GetEXhigh();
            Double_t *errorYlow             = (*ShiftedYield)->GetEYlow();
    //         Double_t *errorYhigh         = (*ShiftedYield)->GetEYhigh();

            // apply bin shift correction
            for (Int_t ib=0; ib<binNumber; ib++) {
                Double_t ptMin              = xPoint[ib]-errorXlow[ib];
                Double_t ptMax              = xPoint[ib]+errorXlow[ib];

                // the bin shift affected value of the fit function in current bin
                Double_t shiftedValue       = CurrentFit->Integral(ptMin,ptMax) / (ptMax - ptMin);
                // the correct value at the bin center
                Double_t trueValue          = CurrentFit->Eval((ptMax + ptMin)/2.);
                // the bin shift correction factor
                ratio[ib]                   = shiftedValue / trueValue;
                errorYlow[ib]               = yValueErrlow[ib]/ ratio[ib];
                errorYlow[ib]               = yValueErrlow[ib]/ ratio[ib];
                yPoint[ib]                  = yValue[ib]/ratio[ib];
                (*ShiftedYield)->SetMarkerStyle(24);
                (*ShiftedYield)->SetMarkerSize(0.9);
                (*ShiftedYield)->SetMarkerColor(kBlack);
                globalRatioBS               = globalRatioBS + ratio[ib];
            }

            globalRatioBS                   = globalRatioBS/(binNumber);

            DrawGammaSetMarkerTF1( CurrentFit, 1, 0.4, kBlue-4);

            Double_t parameter[10];
            CurrentFit->GetParameters(parameter);

            CurrentFit->SetParameter(0,parameter[0]);
            CurrentFit->SetParameter(1,parameter[1]);
            CurrentFit->SetParameter(2,parameter[2]);
            CurrentFit->SetParameter(3,parameter[3]);
            CurrentFit->SetParameter(4,parameter[4]);
            CurrentFit->SetParameter(5,parameter[5]);
            CurrentFit->SetParameter(6,parameter[6]);
            CurrentFit->SetParameter(7,parameter[7]);
            CurrentFit->SetParameter(8,parameter[8]);
            CurrentFit->SetParameter(9,parameter[9]);
    //         cout<<colorBS<<" " << globalRatioBS << "  " << testGlobalRatioBS << endl;

        }
    //     for (Int_t ib=0; ib<=binNumber; ib++) {
    //         cout << ratio[ib] << "\t" << endl;
    //     }
        return CurrentFit;

    }

    // ****************************************************************************************************************
    // ****************************************************************************************************************
    // ****************************************************************************************************************
    TGraphAsymmErrors* ApplyYshiftIndividualSpectra(TGraphAsymmErrors *IndividualSpectum , TF1 *commonFit) {
    //     cout << "entered bin shifting individual spectra y  "<<IndividualSpectum->GetName()  << endl;
        TGraphAsymmErrors* indSpecShifted           = (TGraphAsymmErrors*)IndividualSpectum->Clone();
        TGraphAsymmErrors* dummySpec                = (TGraphAsymmErrors*)IndividualSpectum->Clone();

        Double_t* xvalueI                           = indSpecShifted->GetX();
        Double_t* xvalueIErrlow                     = indSpecShifted->GetEXlow();
        Double_t* xvalueIErrhigh                    = indSpecShifted->GetEXhigh();
        Double_t* yvalueI                           = indSpecShifted->GetY();
        Double_t* yvalueIErrlow                     = dummySpec->GetEYlow();
        Double_t* yvalueIErrhigh                    = dummySpec->GetEYhigh();
        Int_t numberPointsI                         = indSpecShifted->GetN();

        for (Int_t ip = 0; ip < numberPointsI; ip++){
    //       cout << "before:  " << yvalueIErrlow[ip]/yvalueI[ip] << "\t" << yvalueIErrhigh[ip]/yvalueI[ip] << endl;
            Double_t ptMin                          = xvalueI[ip]-xvalueIErrlow[ip];
            Double_t ptMax                          = xvalueI[ip]+xvalueIErrhigh[ip];
            yvalueIErrlow[ip]                       = yvalueIErrlow[ip]/yvalueI[ip];
            yvalueIErrhigh[ip]                      = yvalueIErrhigh[ip]/yvalueI[ip];
            // the bin shift affected value of the fit function in current bin
            Double_t shiftedValue                   = commonFit->Integral(ptMin,ptMax) / (ptMax - ptMin);
            // the correct value at the bin center
            Double_t trueValue                      = commonFit->Eval((ptMax + ptMin)/2.);
            // the bin shift correction factor
            Double_t ratio                          = shiftedValue / trueValue;
            Double_t output                         = yvalueI[ip]/ratio;
            indSpecShifted->SetPoint(ip,xvalueI[ip],yvalueI[ip]/ratio);
            indSpecShifted->SetPointEYlow(ip,yvalueIErrlow[ip]*output);
            indSpecShifted->SetPointEYhigh(ip,yvalueIErrhigh[ip]*output);
        }

    //     indSpecShifted->Print();
        return indSpecShifted;
    }

    // ****************************************************************************************************************
    // ****************************************************************************************************************
    // ****************************************************************************************************************
    TH1D* ApplyYshiftIndividualSpectra(TH1D *IndividualSpectum , TF1 *commonFit) {
    //     cout << "entered bin shifting individual spectra y  "<<IndividualSpectum->GetName()  << endl;
        TH1D* indSpecShifted            = (TH1D*)IndividualSpectum->Clone(IndividualSpectum->GetName());
        TH1D* dummySpec                 = (TH1D*)IndividualSpectum->Clone("b");

        for (Int_t ip = 1; ip < indSpecShifted->GetNbinsX()+1; ip++){
            cout << ip <<":\t" << dummySpec->GetBinCenter(ip) << "\t+-" << dummySpec->GetBinWidth(ip)/2 << "\t "<< dummySpec->GetBinContent(ip) << endl;
            if (!(dummySpec->GetBinContent(ip) > 0)) continue;
            Double_t ptMin                          = dummySpec->GetBinCenter(ip)-dummySpec->GetBinWidth(ip)/2;
            Double_t ptMax                          = dummySpec->GetBinCenter(ip)+dummySpec->GetBinWidth(ip)/2;
            Double_t yvalueIErr                     = dummySpec->GetBinError(ip)/dummySpec->GetBinContent(ip);
            // the bin shift affected value of the fit function in current bin
            Double_t shiftedValue                   = commonFit->Integral(ptMin,ptMax) / dummySpec->GetBinWidth(ip);
            // the correct value at the bin center
            Double_t trueValue                      = commonFit->Eval(dummySpec->GetBinCenter(ip));
            // the bin shift correction factor
            Double_t ratio                          = shiftedValue / trueValue;
            Double_t output                         = dummySpec->GetBinContent(ip)/ratio;

            cout << ratio << "\t" << output<< endl;

            indSpecShifted->SetBinContent(ip,output);
            indSpecShifted->SetBinError(ip,yvalueIErr*output);
        }

        return indSpecShifted;
    }


    // ****************************************************************************************************************
    // ****************************************************************************************************************
    // ****************************************************************************************************************
    Double_t bin_shift_x(TF1 *fYield, Double_t ptMin, Double_t ptMax, Double_t yValue){

        //
        // sample code for bin shift in x direction
        // Klaus Reygers
        //

        // define function describing yield dN/dpT
        // parameter [3] needed for root finding (f(pTlw) - meas = 0)

        // define the bin limits
        Float_t pT1         = ptMin; // lower bin limit in GeV/c
        Float_t pT2         = ptMax; // upper bin limit in GeV/c

        // integrate yield of that bin
        Float_t meas        = 1./(pT2 - pT1) * fYield->Integral(pT1,pT2);

        // Set offset to integrated yield (variable "meas")
        // in order to solve the eq. f(pTlw) - meas = 0
        TString funcName    = fYield->GetName();
        Int_t nPar          = fYield->GetNpar();
    //     cout << "Number of parameter for function " << funcName.Data() << " is " << nPar << endl;
        fYield->SetParameter(nPar-1, meas);

        // prepare root finding with root's BrentRootFinder
        ROOT::Math::WrappedTF1 wfYield(*fYield); // create wrapper function
        ROOT::Math::BrentRootFinder brf; // create root finder
        brf.SetFunction(wfYield, pT1, pT2);   // set range for function wfYield

        // solve the equation f(pTlw) - meas = 0
        brf.Solve();
        Double_t pt0        = brf.Root();

        // correct position of the point according to
        // Lafferty, Wyatt., Nucl. Instr. and Meth. A 355, 541, 1995
        // printf("%4.1f < pt < %4.1f GeV/c, pt0 = %6.3f GeV/c, y=%g, int=%g\n",ptMin,ptMax,pt0,yValue,meas);
        return pt0;

    if(yValue == 0) {}
    }

    // ****************************************************************************************************************
    // *********************************** ApplyXshift ****************************************************************
    // ****************************************************************************************************************
    TGraphAsymmErrors *ApplyXshift(TGraphAsymmErrors *spectrum, TF1 *tsallis, TString mesonType = "Pi0", Bool_t useRangesFromTF1 = kFALSE){
        //----------------------------------------------------------------------
        // This function takes a spectrum, fits it by a function tsallis
        // and calculates pt-shift iteratively.
        // Return value: a new spectrum with pt-shifts for each point.
        // Yuri Kharlov. 19.12.2011
        // modified by Friederike Bock : 19.12.2011
        //----------------------------------------------------------------------

    //     cout << "entered bin shifting x" << endl;
        TGraphAsymmErrors *spectrumShift    = (TGraphAsymmErrors*)spectrum->Clone();
        spectrumShift->SetName(Form("%s_xshift",spectrum->GetName()));
        Int_t nIter                         = 10;
    //     cout << "loaded spectrum" << endl;
        RemoveScalingWithPtGraph(spectrumShift);

        // Define a function Int(f)-f0

    //  Double_t* xvalue                    = spectrum->GetX();
    //  Double_t* xvalueErr                 = spectrum->GetEXlow();
        Int_t numberPoints                  = spectrum->GetN();
    //     cout << "number of bin in x" << numberPoints<< endl;
    //     cout << "red x values from spectrum" << endl;

        TString formulaName                 = Form("%s",tsallis->GetName());
        Int_t nPar                          = tsallis->GetNpar();
        TString formula                     = Form("%s - [%d]",(tsallis->GetExpFormula()).Data(),nPar);
    //     cout << "Formula used for bin shift: " << formula << endl;
        TF1 * fYield;
        if(useRangesFromTF1) fYield        = new TF1(formulaName.Data(), formula, tsallis->GetMinimumX(),tsallis->GetMaximumX());
        else fYield                        = new TF1(formulaName.Data(), formula, 0.2,25.0);

        for (Int_t iter=0; iter<nIter; iter++) {
            Double_t *xPoint     = spectrumShift->GetX();
            Double_t *yPoint     = spectrumShift->GetY();

            Double_t *errorXlow  = spectrumShift->GetEXlow();
            Double_t *errorXhigh = spectrumShift->GetEXhigh();
            Double_t *errorYlow  = spectrumShift->GetEYlow();
            Double_t *errorYhigh = spectrumShift->GetEYhigh();

            spectrumShift->Fit(tsallis, "NRMEXQ0","", xPoint[0]-errorXlow[0],xPoint[numberPoints-1]+errorXhigh[numberPoints-1]);
            for(Int_t p=0; p<nPar; p++){
                fYield->SetParameter(p,tsallis->GetParameter(p));
            }

    //         cout << "iteration " <<iter << endl;
            for (Int_t i=0; i<numberPoints; i++) {
                Double_t ptMin  = xPoint[i]-errorXlow[i];
                Double_t ptMax  = xPoint[i]+errorXhigh[i];
                fYield->SetParameter(nPar,0.);
                Double_t pt0    = bin_shift_x(fYield,ptMin,ptMax,yPoint[i]);
                Double_t dX     = xPoint[i] - pt0;
                spectrumShift->SetPoint(i,pt0,yPoint[i]);
                spectrumShift->SetPointEXlow (i,errorXlow[i] -dX);
                spectrumShift->SetPointEXhigh(i,errorXhigh[i]+dX);
                if (iter == nIter-1) {
                    cout << Form("%4.1f<pt<%4.1f GeV/c, pt0=%6.3f GeV/c, Ed3s/dp3 = %.3g + %.3g - %.3g pb/GeV2",
                    ptMin,ptMax,pt0,yPoint[i],errorYhigh[i],errorYlow[i]) << endl;
                }
            }
        }
    //     cout<<"%Integral measured::"<<100.*tsallis->Integral(xvalue[0]-xvalueErr[0],xvalue[numberPoints-1]+xvalueErr[numberPoints-1])/tsallis->Integral(0.,100.)<<endl;
    //     cout<<"%Integral NOT measured::"<<100.*(1.-tsallis->Integral(xvalue[0]-xvalueErr[0],xvalue[numberPoints-1]+xvalueErr[numberPoints-1])/tsallis->Integral(0.,100.)) <<endl;

        ScaleWithPtGraph(spectrumShift);
    //     spectrumShift->Print();
        return spectrumShift;
    }

    // ****************************************************************************************************************
    // ****************************************************************************************************************
    // ****************************************************************************************************************
    TGraphAsymmErrors* ApplyXshiftIndividualSpectra( TGraphAsymmErrors *CombinedSpectrum,
                                                    TGraphAsymmErrors *IndividualSpectum,
                                                    TF1 *shiftFunction,
                                                    Int_t startBin =0, Int_t numberOfCommonBins = 0, TString mesonType = "Pi0"){
        //----------------------------------------------------------------------
        // function bin shift the individual spectra based the commond fit.
        // startBin must be the bin in the common spectra at which the individual spectra starts
        // numberOfCommonBins must be the numbers of bins where the individual and common spectra have the same binning
        //         -> where they are not the same, a manual shift is done for the individual spectra
        //----------------------------------------------------------------------

    //     CombinedSpectrum->Print();
    //     cout << "entered bin shifting individual spectra x" << endl;
        TGraphAsymmErrors *spectrumShifted  = (TGraphAsymmErrors*)CombinedSpectrum->Clone();
        Double_t* xvalue                    = spectrumShifted->GetX();
        Double_t* yvalue                    = spectrumShifted->GetY();
        Double_t* xvalueErrUp               = spectrumShifted->GetEXhigh();
        Double_t* xvalueErrLow              = spectrumShifted->GetEXlow();

        TGraphAsymmErrors* indSpecShifted   = (TGraphAsymmErrors*)IndividualSpectum->Clone();
    //     IndividualSpectum->Print();
        Double_t* xvalueI                   = indSpecShifted->GetX();

        Double_t* xvalueIErrlow             = indSpecShifted->GetEXlow();
        Double_t* xvalueIErrhigh            = indSpecShifted->GetEXhigh();
        Double_t* yvalueI                   = indSpecShifted->GetY();
        Double_t* yvalueIErrlow             = indSpecShifted->GetEYlow();
        Double_t* yvalueIErrhigh            = indSpecShifted->GetEYhigh();
        Int_t numberPointsI                 = indSpecShifted->GetN();

    //     indSpecShifted->Print();
        cout  << "startBin: " << startBin << "\t largest common bin: " << numberOfCommonBins << "\t number of real common bins: " <<  numberOfCommonBins-startBin  << endl;
        Int_t iBinsInd                      = 0;
        for (Int_t iBins = startBin ; iBins < numberOfCommonBins; iBins++){
            cout << "Common bins: " << iBinsInd << "\t"<< xvalueI[iBinsInd] << "\t"<<xvalue[iBins] << endl;
            //cout << "Before Shift" << yvalue[iBinsInd] << "\t" << yvalueI[iBinsInd] << endl;
            yvalueI[iBinsInd]               = yvalueI[iBinsInd]*xvalueI[iBinsInd]/xvalue[iBins];
                //cout << "after Shift" << yvalue[iBinsInd] << "\t" << yvalueI[iBinsInd] << endl;
            yvalueIErrhigh[iBinsInd]        = yvalueIErrhigh[iBinsInd]*xvalueI[iBinsInd]/xvalue[iBins];
            yvalueIErrlow[iBinsInd]         = yvalueIErrlow[iBinsInd]*xvalueI[iBinsInd]/xvalue[iBins];
            xvalueI[iBinsInd]               = xvalue[iBins];
            xvalueIErrhigh[iBinsInd]        = xvalueErrUp[iBins];
            xvalueIErrlow[iBinsInd]         = xvalueErrLow[iBins];
            iBinsInd++;
        }
        //cout << iBinsInd << "\t" << numberPointsI << endl;
        TString formulaName                 = Form("%s",shiftFunction->GetName());
        Int_t nPar                          = shiftFunction->GetNpar();
        TString formula                     = Form("%s - [%d]",(shiftFunction->GetExpFormula()).Data(),nPar);
        cout << "Formula used for bin shift: " << formula << endl;
        TF1 * fYield                        = new TF1(formulaName.Data(),formula, 0.2,25.);
        for(Int_t p=0; p<nPar; p++){
            fYield->SetParameter(p,shiftFunction->GetParameter(p));
        }

        for (Int_t ip = (iBinsInd); ip < numberPointsI; ip++){
            cout << "Individual shift: " << ip << "\t" <<xvalueI[ip] ;
            Double_t ptMin                  = xvalueI[ip]-xvalueIErrlow[ip];
            Double_t ptMax                  = xvalueI[ip]+xvalueIErrhigh[ip];
            fYield->SetParameter(nPar,0.);
            Double_t pt0                    = bin_shift_x(fYield,ptMin,ptMax,yvalueI[ip]*xvalueI[ip]);
            Double_t dX                     = xvalueI[ip] - pt0;
            indSpecShifted->SetPoint(ip,pt0,yvalueI[ip]*xvalueI[ip]/pt0);
            indSpecShifted->SetPointEXlow(ip,xvalueIErrlow[ip] -dX);
            indSpecShifted->SetPointEXhigh(ip,xvalueIErrhigh[ip]+dX);
            cout << "\t" << pt0 << endl;
        }

        return indSpecShifted;
    }

    // ****************************************************************************************************************
    // ***************** Calculation of RAA with data *****************************************************************
    // ****************************************************************************************************************
    void CalcRaa(   TGraphAsymmErrors* graphPPSpectrum,
                    TGraphAsymmErrors* graphPPSpectrumSystNoMat,
                    TGraphAsymmErrors* graphPPCombinedSpectrum,
                    TF1* fitPP,
                    TGraphAsymmErrors* graphPbPbSpectrum,
                    TGraphAsymmErrors* graphPbPbSpectrumSysNoMat,  //PbPb Yields
                    TGraphAsymmErrors** graphRAA,
                    TGraphAsymmErrors** graphRAASys,
                    Double_t fNcoll,
                    Double_t fNcollError,
                    TString meson                                   = "Pi0",
                    Double_t maxPtPP                                = 0.,
                    Int_t nSubEnd                                   = 0,
                    TString functionInterpolation                   = "powPure",
                    Bool_t quiet                                    = kTRUE
                ){
        //--------------------------------------- extrapolating the pp CONV yield in pT --------------------------------------------------

        TString nameSystem = graphPbPbSpectrum->GetName();
        cout<< nameSystem.Data() << endl;
        TString labelSystem;
        if(nameSystem.Contains("PCM"))
            labelSystem     = "PCM";
        else if(nameSystem.Contains("PHOS"))
            labelSystem     = "PHOS";
        else if(nameSystem.Contains("EMCal"))
            labelSystem     = "EMCal";

        if(!quiet){
            cout << Form("PbPb %s spectrum: ",labelSystem.Data()) << endl;
            graphPbPbSpectrum->Print();
            cout << "\n";
            graphPbPbSpectrumSysNoMat->Print();

            cout << Form("\n PP %s spectrum: ",labelSystem.Data()) << endl;
            graphPPSpectrum->Print();
            cout << "\n";
            graphPPSpectrumSystNoMat->Print();

            cout << "\n combined PP spectrum (for extrapolation): " << endl;
            graphPPCombinedSpectrum->Print();
            cout << "\n";
        }
        TGraphAsymmErrors* dummyPPSpectrum  = (TGraphAsymmErrors*)graphPPCombinedSpectrum->Clone("dummyPPSpectrum");
        Double_t* xBinsPPFit                = dummyPPSpectrum->GetX();
        Double_t* xBinsPPErrFitLow          = dummyPPSpectrum->GetEXlow();
        Double_t* xBinsPPErrFitHigh         = dummyPPSpectrum->GetEXhigh();
        Int_t nBinsPPFit                    = dummyPPSpectrum->GetN();
        if(!quiet){
            cout << "PP pt bins and errors: " << endl;
            for (Int_t i = 0; i < nBinsPPFit; i++){
                // xBinsPPErrFitLow[i] = 0.5*2*(xBinsPPErrFitLow[i]);
                // xBinsPPErrFitHigh[i] = 0.5*2*(xBinsPPErrFitHigh[i]);
                // xBinsPPErrFitLow[i] = 0.001*(xBinsPPFit[i]);
                // xBinsPPErrFitHigh[i] = 0.001*(xBinsPPFit[i]);
                cout << xBinsPPFit[i] << "\t" << xBinsPPErrFitLow[i] << "\t" <<  xBinsPPErrFitHigh[i] << endl;
            }
        }

        TString fittingOptions              = "QNSMEX0+";
        Double_t* xBinsPP                   = graphPPSpectrum->GetX();
        Double_t* xBinsErrPP                = graphPPSpectrum->GetEXlow();
        Int_t nBinsPP                       = graphPPSpectrum->GetN();
        Double_t* yPP                       = graphPPSpectrum->GetY();
        Double_t* yErrLowPP                 = graphPPSpectrum->GetEYlow();
        Double_t* yErrHighPP                = graphPPSpectrum->GetEYhigh();
        Double_t* yErrLowPPSys              = graphPPSpectrumSystNoMat->GetEYlow();
        Double_t* yErrHighPPSys             = graphPPSpectrumSystNoMat->GetEYhigh();

        Double_t maxPtPPForSpec;
        if (maxPtPP !=0.){
            maxPtPPForSpec                  = maxPtPP;
        } else {
            maxPtPPForSpec                  = xBinsPP[nBinsPP-1];
        }

        Double_t* xBinsPbPb                 = graphPbPbSpectrum->GetX();
        Double_t* xBinsErrPbPb              = graphPbPbSpectrum->GetEXlow();
        Double_t* yPbPb                     = graphPbPbSpectrum->GetY();

    //     Double_t* xBinsPPComb            =  graphPPCombinedSpectrum->GetX();
        Int_t nBinsPPComb                   = graphPPSpectrum->GetN();
        Int_t firstBinPbPb                  = 0;
        Double_t decisionBoundary           = 0.0000001;
        while (TMath::Abs(xBinsPP[firstBinPbPb] - xBinsPbPb[0]) >decisionBoundary && firstBinPbPb < nBinsPPComb){
            cout << xBinsPP[firstBinPbPb] << "\t" << xBinsPbPb[0] << endl;
            firstBinPbPb++;
        }

        TF1* fitPPPowerlaw;
        TF1* fitPPPowerlaw2;
        TFitResultPtr resultPP;
        TFitResultPtr resultPPPowerlaw;
        TFitResultPtr resultPPPowerlaw2;

        if(functionInterpolation.Contains("h")){
            if(meson.Contains("Eta")){
                fitPPPowerlaw               = FitObject("h","fitRAARefLevy",meson.Data(),dummyPPSpectrum,0.8, 20.);
                fitPPPowerlaw2              = FitObject("h","fitRAARefLevy2",meson.Data(),dummyPPSpectrum, 0.8, 20.);

                //fit fixed as Pi0:
                //     fitRAARef
                // C_{H}  2360.8849828565      1320.2816544575
                // n      6.3976643035      0.1485156261
                // p_{0} (GeV/c)      0.2725372973      0.0406270874
                // Chi2      3.6183191741      ndf      6.0000000000
                // Chi2/ndf      0.6030531957

                // fitPP->FixParameter(0,2360.8849828565);
                // fitPP->FixParameter(1,6.3976643035);
                // fitPP->FixParameter(2,0.2725372973);
                resultPP                    = dummyPPSpectrum->Fit(fitPP,fittingOptions.Data() , "", 1., maxPtPPForSpec);
                resultPPPowerlaw            = dummyPPSpectrum->Fit(fitPPPowerlaw,fittingOptions.Data() , "", 1., maxPtPPForSpec);
                resultPPPowerlaw2           = dummyPPSpectrum->Fit(fitPPPowerlaw2,fittingOptions.Data() , "", 0.5, maxPtPPForSpec);
            } else {
                fitPPPowerlaw               = FitObject("h","fitRAARefLevy",meson.Data(),dummyPPSpectrum,0.4, 40.);
                fitPPPowerlaw2              = FitObject("h","fitRAARefLevy2",meson.Data(),dummyPPSpectrum, 0.4, 40.);
                resultPP                    = dummyPPSpectrum->Fit(fitPP,fittingOptions.Data() , "", 2., maxPtPPForSpec);
                resultPPPowerlaw            = dummyPPSpectrum->Fit(fitPPPowerlaw,fittingOptions.Data() , "", 1., maxPtPPForSpec);
                resultPPPowerlaw2           = dummyPPSpectrum->Fit(fitPPPowerlaw2,fittingOptions.Data() , "", 2.5, maxPtPPForSpec);
            }
        } else {
            fitPPPowerlaw               = FitObject(functionInterpolation.Data(),"fitRAARefPowerlaw2",meson.Data(),dummyPPSpectrum);
            fitPPPowerlaw2              = FitObject(functionInterpolation.Data(),"fitRAARefPowerlaw4",meson.Data(),dummyPPSpectrum);
            resultPP                    = dummyPPSpectrum->Fit(fitPP,fittingOptions.Data() , "", 3., maxPtPPForSpec);
            resultPPPowerlaw            = dummyPPSpectrum->Fit(fitPPPowerlaw,fittingOptions.Data() , "", 2., maxPtPPForSpec);
            resultPPPowerlaw2           = dummyPPSpectrum->Fit(fitPPPowerlaw2,fittingOptions.Data() , "", 4., maxPtPPForSpec);
        }

        cout << WriteParameterToFile(fitPP)<< endl;
        cout << WriteParameterToFile(fitPPPowerlaw)<< endl;
        cout << WriteParameterToFile(fitPPPowerlaw2)<< endl;

        TGraphAsymmErrors* graphPPSpectrumExtended      = (TGraphAsymmErrors*) graphPbPbSpectrum->Clone("graphPPSpectrumExtended");
        TGraphAsymmErrors* graphPPSpectrumExtendedSys   = (TGraphAsymmErrors*) graphPbPbSpectrumSysNoMat->Clone("graphPPSpectrumExtendedSys");
        TGraphAsymmErrors *relSysErrorFuncPP            = new TGraphAsymmErrors(graphPbPbSpectrum->GetN()-nSubEnd);
        TGraphAsymmErrors *relSysErrorPowerLaw1         = new TGraphAsymmErrors(graphPbPbSpectrum->GetN()-nSubEnd);
        TGraphAsymmErrors *relSysErrorPowerLaw2         = new TGraphAsymmErrors(graphPbPbSpectrum->GetN()-nSubEnd);

        (*graphRAA)                         = new TGraphAsymmErrors(graphPbPbSpectrum->GetN()-nSubEnd);
        (*graphRAASys)                      = new TGraphAsymmErrors(graphPbPbSpectrum->GetN()-nSubEnd);

        Double_t lastSystematicRel          = 0;
        for (Int_t i=0;i<graphPbPbSpectrum->GetN()-nSubEnd ;i++){
    //         Double_t xCenter = 0.5*(xPtLimits[i+1]+xPtLimits[i]);
    //         cout <<  xBinsPP[i+firstBinPbPb] << "\t" << xBinsPbPb[i] << "\t" << xBinsErrPP[i+firstBinPbPb] << "\t" << xBinsErrPbPb[i] << endl;
            if (TMath::Abs(xBinsPP[i+firstBinPbPb] - xBinsPbPb[i]) < decisionBoundary && TMath::Abs(xBinsErrPP[i+firstBinPbPb] - xBinsErrPbPb[i]) < decisionBoundary && xBinsPP[i] < maxPtPPForSpec ){
                cout << "loop with pp points" << endl;
    //          cout<< "xPP: "<< xBinsPP[i+firstBinPbPb]<< " xPPErr: " << xBinsErrPP[i+firstBinPbPb] << " yPP: " << yPP[i+firstBinPbPb]<< endl;
    //          cout << "xPbPb: "<<xBinsPbPb[i]<< " xPbPbErr: " <<xBinsErrPbPb[i] <<"  yPbPb/Ncoll: " << yPbPb[i]/fNcoll<< " Raa: " << yPbPb[i] /(fNcoll*yPP[i+firstBinPbPb])<< endl;
                graphPPSpectrumExtended->SetPoint(i,xBinsPbPb[i],yPP[i+firstBinPbPb]);
                graphPPSpectrumExtendedSys->SetPoint(i,xBinsPbPb[i],yPP[i+firstBinPbPb]);
                graphPPSpectrumExtended->SetPointError(i, xBinsErrPbPb[i], xBinsErrPbPb[i], yErrLowPP[i+firstBinPbPb], yErrHighPP[i+firstBinPbPb]);
                graphPPSpectrumExtendedSys->SetPointError(i, xBinsErrPbPb[i], xBinsErrPbPb[i], yErrLowPPSys[i+firstBinPbPb], yErrHighPPSys[i+firstBinPbPb]);
                if (yPP[i+firstBinPbPb] != 0){
                    lastSystematicRel       = yErrLowPPSys[i+firstBinPbPb]/yPP[i+firstBinPbPb]*100;
                    if(!quiet) cout << "rel.syst.err. PPext : (" << yErrLowPPSys[i+firstBinPbPb] << " / " << yPP[i+firstBinPbPb] << ")*100 = " << lastSystematicRel << "\n" << endl;
                }

                relSysErrorFuncPP->SetPoint(i,xBinsPbPb[i],lastSystematicRel);

            } else {
                cout << "loop with pp interpolation" << endl;
                Double_t ptStart                = xBinsPbPb[i] - xBinsErrPbPb[i];
                Double_t ptEnd                  = xBinsPbPb[i] + xBinsErrPbPb[i];
                Double_t binWidth               = ptEnd-ptStart;

                Double_t yieldPP                = fitPP->Integral(ptStart, ptEnd, resultPP->GetParams())/binWidth;
                Double_t errorYieldPP           = fitPP->IntegralError(ptStart, ptEnd, resultPP->GetParams(), resultPP->GetCovarianceMatrix().GetMatrixArray())/binWidth;
                Double_t relErrPP               = abs(errorYieldPP)/yieldPP *100;

                cout << "yieldPP: " << yieldPP << " relErrPP: " << relErrPP << endl;
                Double_t yieldPPPowerlaw        = fitPPPowerlaw->Integral(ptStart, ptEnd, resultPPPowerlaw->GetParams())/binWidth;
                Double_t errorYieldPPPowerlaw   = fitPPPowerlaw->IntegralError(ptStart, ptEnd, resultPPPowerlaw->GetParams(), resultPPPowerlaw->GetCovarianceMatrix().GetMatrixArray())/binWidth;
                if(!quiet) cout << "pt: " << ptStart << " " << ptEnd << " , yieldPPvar1 " << yieldPPPowerlaw << "+-" << errorYieldPPPowerlaw << endl;
                Double_t relErrPow1             = abs(yieldPPPowerlaw-yieldPP)/yieldPP*100;
                if(!quiet) cout << "relErrPP var1: " << relErrPow1 << endl;

                Double_t yieldPPPowerlaw2       = fitPPPowerlaw2->Integral(ptStart, ptEnd, resultPPPowerlaw2->GetParams())/binWidth;
                Double_t errorYieldPPPowerlaw2  = fitPPPowerlaw2->IntegralError(ptStart, ptEnd, resultPPPowerlaw2->GetParams(), resultPPPowerlaw2->GetCovarianceMatrix().GetMatrixArray())/binWidth;
                if(!quiet) cout<< "pt: "<< ptStart << " " << ptEnd << " , yieldPPvar2 " << yieldPPPowerlaw2 << "+-" << errorYieldPPPowerlaw2 << endl;
                Double_t relErrPow2             = abs(yieldPPPowerlaw2-yieldPP)/yieldPP *100;
                if(!quiet) cout << "relErrPP var2: " << relErrPow2 << endl;

                graphPPSpectrumExtended->SetPoint(i,xBinsPbPb[i],yieldPP);
                graphPPSpectrumExtended->SetPointError(i, xBinsErrPbPb[i], xBinsErrPbPb[i], errorYieldPP, errorYieldPP);
                graphPPSpectrumExtendedSys->SetPoint(i,xBinsPbPb[i],yieldPP);

                relSysErrorFuncPP->SetPoint(i,xBinsPbPb[i],relErrPP);
                relSysErrorPowerLaw1->SetPoint(i,xBinsPbPb[i],relErrPow1);
                relSysErrorPowerLaw2->SetPoint(i,xBinsPbPb[i],relErrPow2);

                Double_t syst                   = 0;
                if(!quiet) cout << "lastSystematicRel: " << lastSystematicRel << endl;
                if (abs(relErrPow1) > abs(relErrPow2)){
                    syst                    = TMath::Sqrt(relErrPow1*relErrPow1 + lastSystematicRel*lastSystematicRel)*yieldPP/100;
                    if(!quiet) cout << "total syst for PP: " << TMath::Sqrt(relErrPow1*relErrPow1 + lastSystematicRel*lastSystematicRel) << endl;
                } else {
                    syst                    = TMath::Sqrt(relErrPow2*relErrPow2 + lastSystematicRel*lastSystematicRel)*yieldPP/100;
                    if(!quiet) cout << "total syst for PP: " << TMath::Sqrt(relErrPow2*relErrPow2 + lastSystematicRel*lastSystematicRel) << endl;
                }
                graphPPSpectrumExtendedSys->SetPointError(i, xBinsErrPbPb[i], xBinsErrPbPb[i],syst, syst); //-syst, +syst);

            }
        }

        TCanvas* canvasDummy6 = new TCanvas("canvasDummy2","",200,10,1200,1100);  // gives the page size
        DrawGammaCanvasSettings( canvasDummy6,  0.1, 0.01, 0.015, 0.08);
        canvasDummy6->SetLogy();
        canvasDummy6->SetLogx();
        TH2F * histo2DDummy6 = new TH2F("histo2DDummy6","histo2DDummy5",1000,0.3,40.,1000,1e-10,10);
        SetStyleHistoTH2ForGraphs(histo2DDummy6, "#it{p}_{T} (GeV/#it{c})","", 0.032,0.04, 0.04,0.04, 1,1.55);
        histo2DDummy6->DrawCopy();

        DrawGammaSetMarkerTGraphAsym(graphPPSpectrum, 24,2.5, kBlue, kBlue);
        graphPPSpectrum->Draw("pEsame");
        DrawGammaSetMarkerTGraphAsym(graphPPCombinedSpectrum, 24,2.5, kRed, kRed);
        graphPPCombinedSpectrum->Draw("pEsame");
        DrawGammaSetMarkerTGraphAsym(graphPPSpectrumExtendedSys, 20,2, kGray+1, kGray+1,1,kTRUE,kGray+1);
        graphPPSpectrumExtendedSys->Draw("2same");
        DrawGammaSetMarkerTGraphAsym(graphPPSpectrumExtended, 20,2, kBlack, kBlack);
        graphPPSpectrumExtended->Draw("pEsame");

        fitPPPowerlaw->SetLineColor(kOrange+2);
        fitPPPowerlaw2->SetLineColor(kGreen+2);
        fitPP->SetLineColor(kBlue+2);
        fitPPPowerlaw->Draw("same");
        fitPPPowerlaw2->Draw("same");
        fitPP->Draw("same");

        TLegend* legend = new TLegend(0.15,0.15,0.5,0.25);
        legend->SetFillColor(0);
        legend->SetLineColor(0);
        legend->SetTextSize(0.03);
        legend->AddEntry(graphPPSpectrum,labelSystem.Data(),"p");
        legend->AddEntry(graphPPCombinedSpectrum,"comb","p");
        legend->AddEntry(graphPPSpectrumExtended,"extrapolated","p");
        legend->Draw();

        canvasDummy6->Update();
        canvasDummy6->Print(Form("debugWithInterpolation_%s%s.eps",labelSystem.Data(), meson.Data()));

        if(!quiet) {
        graphPPSpectrumExtendedSys->Print();
        graphPPSpectrumExtended->Print();
        }

        Double_t* yPPExtended                   =  graphPPSpectrumExtended->GetY();
        Double_t errYlow;
        Double_t errYhigh;

        for(Int_t i=0; i<(*graphRAA)->GetN(); i++){
            (*graphRAA)->SetPoint(i,xBinsPbPb[i],yPbPb[i]/(fNcoll* yPPExtended[i]));
            Double_t errYStat                   = pow( pow(graphPbPbSpectrum->GetErrorYlow(i)/yPbPb[i],2.) +
                                                        pow(graphPPSpectrumExtended->GetErrorYlow(i)/yPPExtended[i],2.),
                                                        0.5)*yPbPb[i]/(fNcoll*yPPExtended[i]);
            (*graphRAA)->SetPointError(i, xBinsErrPbPb[i], xBinsErrPbPb[i], errYStat, errYStat);

            (*graphRAASys)->SetPoint(i,xBinsPbPb[i],yPbPb[i] /(fNcoll* yPPExtended[i]));
            if (TMath::Abs(xBinsPP[i+firstBinPbPb]-xBinsPbPb[i]) < decisionBoundary && TMath::Abs(xBinsErrPP[i+firstBinPbPb]-xBinsErrPbPb[i]) < decisionBoundary && xBinsPP[i]<maxPtPPForSpec ){
            if(!quiet) {
                cout << "loop raa with pp points" << endl;
                cout << "rel syst err PbPb: " << graphPbPbSpectrumSysNoMat->GetErrorYlow(i)/yPbPb[i] << endl;
                cout << "rel syst err pp: "   << graphPPSpectrumExtendedSys->GetErrorYlow(i)/yPPExtended[i] << endl;
            }
                errYlow                         = pow(  pow(graphPbPbSpectrumSysNoMat->GetErrorYlow(i)/yPbPb[i],2.) +
                                                        pow(graphPPSpectrumExtendedSys->GetErrorYlow(i)/yPPExtended[i],2.)
                                                        , 0.5)*yPbPb[i] /(fNcoll* yPPExtended[i]); //+ pow(fNcollError/fNcoll,2.) fNcollError not taking into account
                errYhigh                        = pow(  pow(graphPbPbSpectrumSysNoMat->GetErrorYhigh(i)/yPbPb[i],2.) +
                                                        pow(graphPPSpectrumExtendedSys->GetErrorYhigh(i)/yPPExtended[i],2.)
                                                        , 0.5)*yPbPb[i] /(fNcoll* yPPExtended[i]); //+ pow(fNcollError/fNcoll,2.) fNcollError not taking into account
                if(!quiet) cout << "rel syst err Raa: " << errYlow/(yPbPb[i]/(fNcoll* yPPExtended[i])) << endl;
                if(!quiet) cout << "---------------" << endl;

            } else {
            if(!quiet) {
                cout << "loop raa with extr points" << endl;
                cout << "rel syst err PbPb: " << graphPbPbSpectrumSysNoMat->GetErrorYlow(i)/yPbPb[i] << endl;
                cout << "rel syst err pp: "   << graphPPSpectrumExtendedSys->GetErrorYlow(i)/yPPExtended[i] << endl;
            }
                errYlow                         = pow(  pow(graphPbPbSpectrumSysNoMat->GetErrorYlow(i)/yPbPb[i],2.) +
                                                        pow(graphPPSpectrumExtendedSys->GetErrorYlow(i)/yPPExtended[i],2.)
                                                        , 0.5)*yPbPb[i] /(fNcoll* yPPExtended[i]); //+ pow(fNcollError/fNcoll,2.) fNcollError not taking into account
                errYhigh                        = pow(  pow(graphPbPbSpectrumSysNoMat->GetErrorYhigh(i)/yPbPb[i],2.) +
                                                        pow(graphPPSpectrumExtendedSys->GetErrorYhigh(i)/yPPExtended[i],2.)
                                                        , 0.5)*yPbPb[i] /(fNcoll* yPPExtended[i]); //+ pow(fNcollError/fNcoll,2.) fNcollError not taking into account
                if(!quiet) cout << "rel syst err Raa: " << errYlow/(yPbPb[i]/(fNcoll* yPPExtended[i])) << endl;
                if(!quiet) cout << "---------------" << endl;
            }
            (*graphRAASys)->SetPointError(i, xBinsErrPbPb[i], xBinsErrPbPb[i], errYlow, errYhigh);
        }

        Int_t b = 0;
        while (yPbPb[b] == 0){
            (*graphRAA)->RemovePoint(0);
            (*graphRAASys)->RemovePoint(0);
            b++;
        }
        delete dummyPPSpectrum;
        if(fNcollError){}
    }

    // ****************************************************************************************************************
    // ********************** Calculation of RAA to a theory graph as replacement for pp reference ********************
    // ****************************************************************************************************************
    void CalcRaaWithTheoryGraph(    TGraphErrors* graphTheory,                 // Theory graph
                                    TGraphAsymmErrors* graphPbPbSpectrum,   // PbPb Yields
                                    TGraphAsymmErrors** graphRAA             // RAA return graphs
                            ){

        cout << "PbPb spectrum: " << endl;
        graphPbPbSpectrum->Print();

        TH1D* dummyPPSpectrum       = GraphToHist(graphTheory,graphTheory->GetN(),"dummyPPSpectrum");
        Double_t startTheory        = graphTheory->GetX()[0];
        for (Int_t i=1; i< graphTheory->GetN(); i++){
            cout << dummyPPSpectrum->GetBinCenter(i) << "\t" << dummyPPSpectrum->GetBinContent(i) << endl;
        }
        cout << "theory starts at: " << startTheory << endl;

        Double_t* xBinsPbPb         = graphPbPbSpectrum->GetX();
        Double_t* xBinsErrPbPb      = graphPbPbSpectrum->GetEXlow();
        Double_t* yPbPb             = graphPbPbSpectrum->GetY();

        Int_t firstBinPbPb          = 0;
        Double_t decisionBoundary   = 0.0000001;
        while ( (xBinsPbPb[firstBinPbPb]- startTheory) < 0 && firstBinPbPb < graphPbPbSpectrum->GetN()){
            cout << startTheory << "\t" << xBinsPbPb[firstBinPbPb] << endl;
            firstBinPbPb++;
        }


        (*graphRAA)                 = new TGraphAsymmErrors(graphPbPbSpectrum->GetN()-firstBinPbPb);

        for(Int_t i=firstBinPbPb; i<(*graphRAA)->GetN(); i++){
            Double_t pp             = dummyPPSpectrum->Interpolate(xBinsPbPb[i]);
            Double_t rAA            = yPbPb[i]/pp;

            (*graphRAA)->SetPoint(i,xBinsPbPb[i],rAA);
            Double_t errYStat       = graphPbPbSpectrum->GetErrorYlow(i)/yPbPb[i]*rAA ;

            (*graphRAA)->SetPointError(i, xBinsErrPbPb[i], xBinsErrPbPb[i], errYStat, errYStat);
        }

        Int_t b                     = 0;
        while (yPbPb[b] == 0){
            (*graphRAA)->RemovePoint(0);
            b++;
        }

        delete dummyPPSpectrum;
    }

    // ****************************************************************************************************************
    // ********************** Calculation of RAA to a theory function as replacement for pp reference *****************
    // ****************************************************************************************************************
    void CalcRaaWithTheoryFit(  TF1* fitTheory,                                                 // Theory graph
                                TGraphAsymmErrors* graphPbPbSpectrum,                           // PbPb Yields
                                TGraphAsymmErrors** graphRAA,                                   // RAA return graphs
                                Double_t* binningX                      = NULL,                 // binning
                                Int_t nBins                             = 0                     // number of bins
                            ){

        cout << "PbPb spectrum: " << endl;
        graphPbPbSpectrum->Print();

        Double_t* xBinsPbPb         = graphPbPbSpectrum->GetX();
        Double_t* xBinsErrPbPb      = graphPbPbSpectrum->GetEXlow();
        Double_t* yPbPb             = graphPbPbSpectrum->GetY();
        Int_t firstBinPbPb          = 0;
        (*graphRAA)                 = new TGraphAsymmErrors(graphPbPbSpectrum->GetN()-firstBinPbPb);

        for(Int_t i=firstBinPbPb; i<(*graphRAA)->GetN(); i++){
            Double_t pp = 0;
            if (xBinsErrPbPb[i] != 0.) {
                pp                  = fitTheory->Integral(xBinsPbPb[i]-xBinsErrPbPb[i], xBinsPbPb[i]+xBinsErrPbPb[i])/(2*xBinsErrPbPb[i]);
            } else {
                if (binningX){
                    Int_t currentBin    = 0;
                    while (binningX[currentBin] < xBinsPbPb[i] && currentBin< nBins )  currentBin++;
                    cout << binningX[currentBin-1] << "\t" << binningX[currentBin] << "\t" << xBinsPbPb[i] << endl;
                    pp                  = fitTheory->Integral(binningX[currentBin-1], binningX[currentBin])/(binningX[currentBin]-binningX[currentBin-1]);
                } else {
                    cout << "calculation failed, no binning given!" << endl;
                    return;
                }
            }
            Double_t rAA            = yPbPb[i]/pp;
            (*graphRAA)->SetPoint(i,xBinsPbPb[i],rAA);
            Double_t errYStat       = graphPbPbSpectrum->GetErrorYlow(i)/yPbPb[i]*rAA ;
            (*graphRAA)->SetPointError(i, xBinsErrPbPb[i], xBinsErrPbPb[i], errYStat, errYStat);

        }

        Int_t b                     = 0;
        while (yPbPb[b] == 0){
            (*graphRAA)->RemovePoint(0);
            b++;
        }
    }


    // ****************************************************************************************************************
    // ********************** Calculation of RpA with Graphs in same binning ******************************************
    // ****************************************************************************************************************
    TGraphAsymmErrors* CalcRpPbV2(  TGraphAsymmErrors* PPSpectrumStatErr,
                                    TGraphAsymmErrors* PPSpectrumSystErr,
                                    TGraphAsymmErrors* pPbSpectrumStatErr,
                                    TGraphAsymmErrors* pPbSpectrumSystErr,
                                    TGraphAsymmErrors** graphRpPbStatErr,
                                    TGraphAsymmErrors** graphRpPbSystErr,
                                    Double_t tAA,
                                    Double_t tAAErr,
                                    vector<TString>nameSysTakeOut,
                                    TString fileNameSysPPb                  = "",
                                    TString fileNameSysPP                   = "",
                                    TString fileNameSysOut                  = ""
                                ){

        Int_t nPoints                           = pPbSpectrumStatErr->GetN();

        (*graphRpPbStatErr)                     = new TGraphAsymmErrors( nPoints);
        (*graphRpPbSystErr)                     = new TGraphAsymmErrors( nPoints);
        TGraphAsymmErrors* graphRpPbCombErr     = new TGraphAsymmErrors( nPoints);

        Double_t *xBins                         = pPbSpectrumStatErr->GetX();
        Double_t *xErrlow                       = pPbSpectrumStatErr->GetEXlow();
        Double_t *xErrhigh                      = pPbSpectrumStatErr->GetEXhigh();
        Double_t *ypPbBins                      = pPbSpectrumStatErr->GetY();
        Double_t *yPPBins                       = PPSpectrumStatErr->GetY();
        Double_t *ypPbStatErrlow                = pPbSpectrumStatErr->GetEYlow();
        //Double_t *ypPbStatErrhigh             = pPbSpectrumStatErr->GetEYlow();
        Double_t *ypPbSystErrlow                = pPbSpectrumSystErr->GetEYlow();
        Double_t *ypPbSystErrhigh               = pPbSpectrumSystErr->GetEYlow();
        Double_t *yPPStatErrlow                 = PPSpectrumStatErr->GetEYlow();
        //Double_t *yPPStatErrhigh = PPSpectrumStatErr->GetEYlow();
        Double_t *yPPSystErrlow                 = PPSpectrumSystErr->GetEYlow();
        Double_t *yPPSystErrhigh                = PPSpectrumSystErr->GetEYlow();
        Double_t *RpPb                          = new Double_t[nPoints];

        Int_t nBinsSysPtPPb                     = 0;
        Int_t nSysAvailSinglePPb                = 0;
        Int_t nBinsSysPtPP                      = 0;
        Int_t nSysAvailSinglePP                 = 0;
        Bool_t haveDetailedSyspPb               = kFALSE;
        Bool_t haveDetailedSysPP                = kFALSE;
        vector<Double_t>unCorrSys;
        vector<Double_t>ptSysExternal;

        if (fileNameSysPPb.CompareTo("") != 0 && fileNameSysPP.CompareTo("") != 0){
            cout << "fileNames correctly set" << endl;
            cout << "pPb: " << fileNameSysPPb.Data() << endl;
            cout << "pp: " << fileNameSysPP.Data() << endl;
            ifstream fileSysErrDetailedPPb;
            fileSysErrDetailedPPb.open(fileNameSysPPb,ios_base::in);

            // check if the file exists
            if(fileSysErrDetailedPPb.is_open()) {
                haveDetailedSyspPb      = kTRUE;
            }
            vector<TString>nameSysPPb;
            // possibly 100 pt bins
            vector<Double_t>* ptSysSplitPPb     = new vector<Double_t>[100];
            vector<Bool_t>enablePPbSys;
            // read detailed file pPb
            if (haveDetailedSyspPb){
                Int_t iPtBin                = 0;
                Bool_t isFirstLine          = kTRUE;
                string line;
                Int_t nDiffErrContribPPb    = 0;

                while (getline(fileSysErrDetailedPPb, line) && iPtBin < 100) {
                    istringstream ss(line);
                    TString temp        ="";
                    if (isFirstLine){
                        while(ss && nDiffErrContribPPb < 100){
                            ss >> temp;
                            if( !(iPtBin==0 && temp.CompareTo("bin")==0) && !temp.IsNull()){
                                nameSysPPb.push_back(temp);
                                nDiffErrContribPPb++;
                            }
                        }
                        nameSysPPb.push_back("TotalError");
                        nDiffErrContribPPb++;
                        isFirstLine             = kFALSE;
                    } else {
                        Int_t nRunning          = 0;
                        while(ss && nRunning < nDiffErrContribPPb){
                            ss >> temp;
                            ptSysSplitPPb[iPtBin].push_back(temp.Atof());
                            nRunning++;
                        }
                        iPtBin++;
                    }
                }
                fileSysErrDetailedPPb.close();

                nBinsSysPtPPb             = iPtBin;
                nSysAvailSinglePPb        = (Int_t)nameSysPPb.size()-1;
                cout <<  nSysAvailSinglePPb << " individual errors:"<< endl;
                for (Int_t k = 0; k < nSysAvailSinglePPb+1; k++ ){
                    cout << ((TString)nameSysPPb.at(k)).Data() << "\t";
                    Bool_t enabled      = kTRUE;
                    for (Int_t m = 0; m < (Int_t)nameSysTakeOut.size() ; m++){
                        if (((TString)nameSysPPb.at(k)).CompareTo((TString)nameSysTakeOut.at(m)) == 0)
                            enabled     = kFALSE;
                    }
                    if (((TString)nameSysPPb.at(k)).CompareTo("TotalError") == 0 || ((TString)nameSysPPb.at(k)).CompareTo("Pt") == 0 || ((TString)nameSysPPb.at(k)).CompareTo("pt") == 0 || ((TString)nameSysPPb.at(k)).CompareTo("bin") == 0 )
                        enabled     = kFALSE;

                    enablePPbSys.push_back(enabled);
                }
                cout << endl;
                for (Int_t k = 0; k < nSysAvailSinglePPb+1; k++ ){
                    cout << enablePPbSys.at(k) << "\t" ;
                }
                cout << endl;
            }

            // check if the file exists
            ifstream fileSysErrDetailedPP;
            fileSysErrDetailedPP.open(fileNameSysPP,ios_base::in);

            if(fileSysErrDetailedPP.is_open()) {
                haveDetailedSysPP       = kTRUE;
            }
            vector<TString>nameSysPP;
            // possibly 100 pt bins
            vector<Double_t>* ptSysSplitPP     = new vector<Double_t>[100];
            vector<Bool_t>enablePPSys;
            // read detailed file pp
            if (haveDetailedSysPP){
                cout << fileNameSysPP.Data() << endl;
                Int_t iPtBin                = 0;
                Bool_t isFirstLine          = kTRUE;
                string line;
                Int_t nDiffErrContribPP     = 0;

                while (getline(fileSysErrDetailedPP, line) && iPtBin < 100) {
                    istringstream ss(line);
                    TString temp        ="";
                    if (isFirstLine){
                        while(ss && nDiffErrContribPP < 100){
                            ss >> temp;
                            if( !(iPtBin==0 && temp.CompareTo("bin")==0) && !temp.IsNull()){
                                nameSysPP.push_back(temp);
                                nDiffErrContribPP++;
                            }
                        }
                        isFirstLine             = kFALSE;
                    } else {
                        Int_t nRunning          = 0;
                        while(ss && nRunning < nDiffErrContribPP){
                            ss >> temp;
                            ptSysSplitPP[iPtBin].push_back(temp.Atof());
                            nRunning++;
                        }
                        iPtBin++;
                    }
                }
                fileSysErrDetailedPP.close();

                nBinsSysPtPP             = iPtBin;
                nSysAvailSinglePP        = (Int_t)nameSysPP.size()-1;
                cout <<  nSysAvailSinglePP << " individual errors pp:"<< endl;
                for (Int_t k = 0; k < nSysAvailSinglePP+1; k++ ){
                    cout << ((TString)nameSysPP.at(k)).Data() << "\t";
                    Bool_t enabled      = kTRUE;
                    for (Int_t m = 0; m < (Int_t)nameSysTakeOut.size() ; m++){
                        if (((TString)nameSysPP.at(k)).CompareTo((TString)nameSysTakeOut.at(m)) == 0)
                            enabled     = kFALSE;
                    }
                    if (((TString)nameSysPP.at(k)).CompareTo("TotalErrorUncorrPP") == 0 || ((TString)nameSysPP.at(k)).CompareTo("Pt") == 0 || ((TString)nameSysPP.at(k)).CompareTo("pt") == 0 || ((TString)nameSysPP.at(k)).CompareTo("bin") == 0 || ((TString)nameSysPP.at(k)).CompareTo("TotalErrorUncorrPP") == 0 )
                        enabled     = kFALSE;

                    enablePPSys.push_back(enabled);
                }
                cout << endl;
                for (Int_t k = 0; k < nSysAvailSinglePP+1; k++ ){
                    cout << enablePPSys.at(k) << "\t" ;
                }
                cout << endl;
            }

            fstream SysErrDatAverOut;
            cout << fileNameSysOut << endl;
            SysErrDatAverOut.open(fileNameSysOut, ios::out);

            if (haveDetailedSyspPb && haveDetailedSysPP){
                SysErrDatAverOut << nameSysPPb.at(0) << "\t" ;
                for (Int_t k = 1; k < nSysAvailSinglePPb; k++){
                    if (enablePPbSys.at(k)) {
                        SysErrDatAverOut << nameSysPPb.at(k) << "\t" ;
                    }
                }
                for (Int_t k = 1; k < nSysAvailSinglePP; k++){
                    if (enablePPSys.at(k)) {
                        SysErrDatAverOut << nameSysPP.at(k) << "\t" ;
                    }
                }
                SysErrDatAverOut << endl;

                for (Int_t iPt = 0; iPt < nBinsSysPtPPb; iPt++){
                    Double_t uncorr     = 0;
                    if ((Double_t)ptSysSplitPPb[iPt].at(0) != (Double_t)ptSysSplitPP[iPt].at(0))
                        cout << "BINNING external sources out of sync" << "\t" << (Double_t)ptSysSplitPPb[iPt].at(0) << "\t"<< (Double_t)ptSysSplitPP[iPt].at(0) << endl;

                    ptSysExternal.push_back((Double_t)ptSysSplitPPb[iPt].at(0));
                    SysErrDatAverOut << (Double_t)ptSysSplitPPb[iPt].at(0) << "\t";
                    for (Int_t k = 1; k < nSysAvailSinglePPb; k++){
                        if (enablePPbSys.at(k)) {
                            uncorr          = uncorr + (Double_t)ptSysSplitPPb[iPt].at(k)*(Double_t)ptSysSplitPPb[iPt].at(k);
                            SysErrDatAverOut << ptSysSplitPPb[iPt].at(k) << "\t" ;
                        }
                    }
                    for (Int_t k = 1; k < nSysAvailSinglePP; k++){
                        if (enablePPSys.at(k)) {
                            uncorr          = uncorr + (Double_t)ptSysSplitPP[iPt].at(k)*(Double_t)ptSysSplitPP[iPt].at(k);
                            SysErrDatAverOut << ptSysSplitPP[iPt].at(k) << "\t" ;
                        }
                    }
                    uncorr              = TMath::Sqrt(uncorr);
                    unCorrSys.push_back(uncorr);
                    SysErrDatAverOut << uncorr << endl;
                }
                cout << "made it here: " << __LINE__ << endl;
                SysErrDatAverOut.close();
            } else {
                cout << "****************************************************************" << endl;
                cout << "****************************************************************" << endl;
                cout << "WARNING: at least one of the single error files wasn't available" << endl;
                cout << "****************************************************************" << endl;
                cout << "****************************************************************" << endl;
            }
        }

        for(Int_t iPoint = 0; iPoint < nPoints; iPoint++){
            cout << "calculating point: "<< iPoint << "\t" << xBins[iPoint] << endl;
            RpPb[iPoint]            = ypPbBins[iPoint] / (tAA* yPPBins[iPoint]);

            (*graphRpPbSystErr)->SetPoint( iPoint, xBins[iPoint], RpPb[iPoint]);
            (*graphRpPbStatErr)->SetPoint( iPoint, xBins[iPoint], RpPb[iPoint] );
            graphRpPbCombErr->SetPoint( iPoint, xBins[iPoint], RpPb[iPoint] );

            Double_t errYStat       = pow( pow( ypPbStatErrlow[iPoint]/ypPbBins[iPoint], 2. ) + pow( yPPStatErrlow[iPoint]/yPPBins[iPoint],  2.), 0.5)* RpPb[iPoint];
            Double_t errYSystlow    = pow( pow( ypPbSystErrlow[iPoint]/ypPbBins[iPoint], 2. ) + pow( yPPSystErrlow[iPoint]/yPPBins[iPoint]  ,2. ) ,   0.5) * RpPb[iPoint];
            Double_t errYSysthigh   = pow( pow( ypPbSystErrhigh[iPoint]/ypPbBins[iPoint],2. ) + pow( yPPSystErrhigh[iPoint]/yPPBins[iPoint], 2. ) ,   0.5) * RpPb[iPoint];
            Double_t errYComblow    = pow( pow( ypPbStatErrlow[iPoint]/ypPbBins[iPoint], 2. ) + pow( ypPbSystErrlow[iPoint]/ypPbBins[iPoint],2. )
                                        + pow( yPPSystErrlow[iPoint]/yPPBins[iPoint], 2. ) ,   0.5) * RpPb[iPoint];
            Double_t errYCombhigh   = pow( pow( ypPbStatErrlow[iPoint]/ypPbBins[iPoint], 2. ) + pow( ypPbSystErrhigh[iPoint]/ypPbBins[iPoint],2. )
                                            + pow( yPPSystErrhigh[iPoint]/yPPBins[iPoint], 2. ) ,   0.5) * RpPb[iPoint];

            if (haveDetailedSyspPb && haveDetailedSyspPb){
                if (TMath::Abs((Double_t)ptSysExternal.at(iPoint) - xBins[iPoint]) < 0.0001){
                    errYSystlow     = (Double_t)unCorrSys.at(iPoint)/100.*RpPb[iPoint];
                    errYSysthigh    = (Double_t)unCorrSys.at(iPoint)/100.*RpPb[iPoint];
                    errYComblow     = pow( pow( ypPbStatErrlow[iPoint]/ypPbBins[iPoint], 2.) + pow( (Double_t)unCorrSys.at(iPoint)/100.,2.),   0.5) * RpPb[iPoint];
                    errYCombhigh    = pow( pow( ypPbStatErrlow[iPoint]/ypPbBins[iPoint], 2.) + pow( (Double_t)unCorrSys.at(iPoint)/100.,2.) ,   0.5) * RpPb[iPoint];

                } else {
                    cout << "BINNING of detailed syst incorrect!!" << endl<< endl<< endl;
                }
            }
            (*graphRpPbStatErr)->SetPointError(iPoint, xErrlow[iPoint],xErrhigh[iPoint],errYStat,errYStat);
            (*graphRpPbSystErr)->SetPointError(iPoint, xErrlow[iPoint],xErrhigh[iPoint],errYSystlow,errYSysthigh);
            graphRpPbCombErr->SetPointError(iPoint, xErrlow[iPoint],xErrhigh[iPoint],errYComblow,errYCombhigh);
        }
        cout << "finished RpA calc" << endl;

        return graphRpPbCombErr;

    }


    // ******************************************************************************************************************
    //****************************************** Integrate spectrum above a certain pt **********************************
    // ******************************************************************************************************************
    TGraphAsymmErrors* SubtractPromptPhotonsViaFit(     TF1* fitTheory,                                             // Theory graph
                                                        TGraphAsymmErrors* graphPbPbSpectrum,                       // PbPb Yields
                                                        Double_t* binningX                      = NULL,             // binning
                                                        Int_t nBins                             = 0                 // number of bins
                                                    ){

        cout << "PbPb spectrum: " << endl;
        graphPbPbSpectrum->Print();
        TGraphAsymmErrors* graphReturn      = (TGraphAsymmErrors*)graphPbPbSpectrum->Clone("dummy");
        Double_t* xBinsPbPb                 = graphPbPbSpectrum->GetX();
        Double_t* xBinsErrPbPb              = graphPbPbSpectrum->GetEXlow();
        Double_t* yPbPb                     = graphPbPbSpectrum->GetY();
        Int_t firstBinPbPb                  = 0;

        for(Int_t i=firstBinPbPb; i<graphReturn->GetN(); i++){
            Double_t pp                     = 0;
            if (xBinsErrPbPb[i] != 0.) {
                pp                          = fitTheory->Integral(xBinsPbPb[i]-xBinsErrPbPb[i], xBinsPbPb[i]+xBinsErrPbPb[i])/(2*xBinsErrPbPb[i]);
            } else {
                if (binningX){
                    Int_t currentBin        = 0;
                    while (binningX[currentBin] < xBinsPbPb[i] && currentBin< nBins )  currentBin++;
                    cout << binningX[currentBin-1] << "\t" << binningX[currentBin] << "\t" << xBinsPbPb[i] << endl;
                    pp                      = fitTheory->Integral(binningX[currentBin-1], binningX[currentBin])/(binningX[currentBin]-binningX[currentBin-1]);
                } else {
                    cout << "calculation failed, no binning given!" << endl;
                    return NULL;
                }
            }
            Double_t y                      = yPbPb[i]-pp;

            cout << y << "\t" << yPbPb[i] << "\t" << pp << endl;
            graphReturn->SetPoint(i,xBinsPbPb[i],y);
            Double_t errYStat               = graphPbPbSpectrum->GetErrorYlow(i)/yPbPb[i]*y ;
            graphReturn->SetPointError(i, xBinsErrPbPb[i], xBinsErrPbPb[i], errYStat, errYStat);
        }

        Int_t b                             = 0;
        while (yPbPb[b] == 0 ){
            graphReturn->RemovePoint(0);
            b++;
        }
        return graphReturn;
    }



    // ****************************************************************************************************************
    // ****************************** Calculate RCP between different centralities ************************************
    // ****************************************************************************************************************
    void CalcRcp(   TGraphAsymmErrors* graphCentralSpectrum, TGraphAsymmErrors* graphPeripheralSpectrum, TGraphAsymmErrors* graphPbPbCentralSpectrumSysNoMat,
                    TGraphAsymmErrors* graphPbPbPeripheralSpectrumSysNoMat, TGraphAsymmErrors** graphRCP, TGraphAsymmErrors** graphRCPSys,
                    Double_t fNcollCentral, Double_t fNcollPeripheral,  TString meson = "Pi0"){//, Double_t fNcollCentralError, Double_t fNcollPeripheralError){
        //Ncoll error is not considered
        cout << "Print central spectrum: " << endl;
        graphCentralSpectrum->Print();
        graphPeripheralSpectrum->Print();

        Double_t* xBinsCentralPbPb          = graphCentralSpectrum->GetX();
        Double_t* xBinsCentralErrPbPb       = graphCentralSpectrum->GetEXlow();
        Double_t* yCentralPbPb              = graphCentralSpectrum->GetY();

        Double_t* xBinsPeripheralPbPb       = graphPeripheralSpectrum->GetX();
        Double_t* xBinsPeripheralErrPbPb    = graphPeripheralSpectrum->GetEXlow();
        Double_t* yPeripheralPbPb           = graphPeripheralSpectrum->GetY();

        (*graphRCP)                         = new TGraphAsymmErrors(graphCentralSpectrum->GetN());
        (*graphRCPSys)                      = new TGraphAsymmErrors(graphCentralSpectrum->GetN());

        for(Int_t i=0; i<(*graphRCP)->GetN(); i++){
            if(xBinsCentralPbPb[i] == xBinsPeripheralPbPb[i]){

                (*graphRCP)->SetPoint(i,xBinsCentralPbPb[i],(yCentralPbPb[i]*fNcollPeripheral)/(yPeripheralPbPb[i]*fNcollCentral));

                Double_t errYStat           = pow(  pow( (graphCentralSpectrum->GetErrorYlow(i)*fNcollPeripheral)/(yPeripheralPbPb[i]*fNcollCentral),2.) +
                                                    pow( (graphPeripheralSpectrum->GetErrorYlow(i)*yCentralPbPb[i]*fNcollPeripheral)/(yPeripheralPbPb[i] *
                                                    yPeripheralPbPb[i] *fNcollCentral *fNcollCentral),2.), 0.5);
                (*graphRCP)->SetPointError(i, xBinsCentralErrPbPb[i], xBinsCentralErrPbPb[i], errYStat, errYStat);


                (*graphRCPSys)->SetPoint(i,xBinsCentralPbPb[i],(yCentralPbPb[i]*fNcollPeripheral)/(yPeripheralPbPb[i]*fNcollCentral));

                Double_t errYlow            = pow(  pow( (graphPbPbCentralSpectrumSysNoMat->GetErrorYlow(i)*fNcollPeripheral)/(yPeripheralPbPb[i]*fNcollCentral),2.) +
                                                    pow( (graphPbPbPeripheralSpectrumSysNoMat->GetErrorYlow(i)*yCentralPbPb[i]*fNcollPeripheral)/(yPeripheralPbPb[i] *
                                                    yPeripheralPbPb[i] *fNcollCentral *fNcollCentral),2.), 0.5); //+ something like pow(fNcollError/fNcoll,2.) fNcollError not taking into account

                Double_t errYhigh           = pow(  pow( (graphPbPbCentralSpectrumSysNoMat->GetErrorYhigh(i)*fNcollPeripheral)/(yPeripheralPbPb[i]*fNcollCentral),2.) +
                                                    pow( (graphPbPbPeripheralSpectrumSysNoMat->GetErrorYhigh(i)*yCentralPbPb[i]*fNcollPeripheral)/(yPeripheralPbPb[i] *
                                                    yPeripheralPbPb[i] *fNcollCentral *fNcollCentral),2.), 0.5); //+ something like pow(fNcollError/fNcoll,2.) fNcollError not taking into account

                cout << i << "\t num: " << yCentralPbPb[i]<< "\t denum: " << yPeripheralPbPb[i] << "\t RCP: " << (yCentralPbPb[i]*fNcollPeripheral)/(yPeripheralPbPb[i]*fNcollCentral) <<"\t" << errYlow << "\t" << errYhigh << endl;

                (*graphRCPSys)->SetPointError(i, xBinsCentralErrPbPb[i], xBinsCentralErrPbPb[i], errYlow, errYhigh);
            }
        }
        Int_t b                             = 0;
        while (yCentralPbPb[b] == 0){
            (*graphRCP)->RemovePoint(0);
            (*graphRCPSys)->RemovePoint(0);
            b++;
        }
    }

    // ****************************************************************************************************************
    // ******************************* Blastwave fit for pi0 **********************************************************
    // ****************************************************************************************************************
    Double_t BWdndptPi0(Double_t *x, Double_t *par) {

        Double_t pt         = x[0]*1000.; //GeV->MeV
        Double_t t          = par[1];
        Double_t beta       = abs(par[2]);
        Double_t yt         = 0.5*TMath::Log((1+beta)/(1-beta)); // atanh(beta)
        Double_t m0         = 135.; //pi0
        Double_t mt         = TMath::Sqrt(m0*m0+pt*pt);
        Double_t f          = par[0]*pt*mt*TMath::BesselI(0,pt*sinh(yt)/t)*TMath::BesselK(1,TMath::Abs(mt*cosh(yt)/t));

        return f;
    }

    // ****************************************************************************************************************
    // ******************************* Blastwave fit for eta **********************************************************
    // ****************************************************************************************************************
    Double_t BWdndptEta(Double_t *x, Double_t *par) {

        Double_t pt         = x[0]*1000.; //GeV->MeV
        Double_t t          = par[1];
        Double_t beta       = abs(par[2]);
        Double_t yt         = 0.5*TMath::Log((1+beta)/(1-beta)); // atanh(beta)
        Double_t m0         = 548.; //eta
        Double_t mt         = TMath::Sqrt(m0*m0+pt*pt);
        Double_t f          = par[0]*pt*mt*TMath::BesselI(0,pt*sinh(yt)/t)*TMath::BesselK(1,TMath::Abs(mt*cosh(yt)/t));

        return f;
    }

    // ****************************************************************************************************************
    // ****************************************************************************************************************
    // ****************************************************************************************************************
    void CorrectGammaDataUnfold(TH1D* histoGammaSpecCorr,TH1D* histoConvProb, TH1D* histoRecoEff, Double_t deltaEtaDummy, Double_t scalingDummy, Double_t nEvtDummy){
        histoGammaSpecCorr->Divide(histoGammaSpecCorr,histoConvProb,1.,1.,"");
        histoGammaSpecCorr->Divide(histoGammaSpecCorr,histoRecoEff,1.,1.,"");
        histoGammaSpecCorr->Scale(1./deltaEtaDummy);
        histoGammaSpecCorr->Scale(scalingDummy);
        histoGammaSpecCorr->Scale(1./nEvtDummy);
        for (Int_t i = 1; i < histoGammaSpecCorr->GetNbinsX()+1 ; i++){
            Double_t newBinContent  = histoGammaSpecCorr->GetBinContent(i)/histoGammaSpecCorr->GetBinCenter(i);
            Double_t newBinError    = histoGammaSpecCorr->GetBinError(i)/histoGammaSpecCorr->GetBinCenter(i);
            histoGammaSpecCorr->SetBinContent(i,newBinContent);
            histoGammaSpecCorr->SetBinError(i,newBinError);
        }
    }

    // ****************************************************************************************************************
    // ****************************************************************************************************************
    // ****************************************************************************************************************
    TH1D *FixEfficiency(TH1D *FixedEff, TH1D* Eff, TString option, TString centString){
        if(option.CompareTo("2.76TeV") == 0){
            FixedEff->SetBinContent(8,(Eff->GetBinContent(7)+Eff->GetBinContent(9))/2);
            FixedEff->SetBinContent(13,(Eff->GetBinContent(12)+Eff->GetBinContent(14))/2);
        }
        if(centString.CompareTo("0-5%") == 0){
            FixedEff->SetBinContent(7,(Eff->GetBinContent(6)+Eff->GetBinContent(8))/2);
            FixedEff->SetBinContent(10,(Eff->GetBinContent(9)+Eff->GetBinContent(11))/2);
            FixedEff->SetBinContent(12,(Eff->GetBinContent(11)+Eff->GetBinContent(13))/2);
            FixedEff->SetBinContent(15,(Eff->GetBinContent(14)+Eff->GetBinContent(16))/2);
        } else if(centString.CompareTo("5-10%") == 0){
            FixedEff->SetBinContent(15,(Eff->GetBinContent(14)+Eff->GetBinContent(16))/2);
            FixedEff->SetBinContent(18,(Eff->GetBinContent(17)+Eff->GetBinContent(19))/2);
        } else if(centString.CompareTo("0-10%") == 0){
            FixedEff->SetBinContent(10,(Eff->GetBinContent(9)+Eff->GetBinContent(11))/2);
            FixedEff->SetBinContent(13,(Eff->GetBinContent(12)+Eff->GetBinContent(14))/2);
            FixedEff->SetBinContent(16,(Eff->GetBinContent(15)+Eff->GetBinContent(17))/2);
        } else if(centString.CompareTo("0-20%") == 0){
            FixedEff->SetBinContent(10,(Eff->GetBinContent(9)+Eff->GetBinContent(11))/2);
            FixedEff->SetBinContent(15,(Eff->GetBinContent(14)+Eff->GetBinContent(16))/2);
        } else if(centString.CompareTo("0-80%") == 0){
            FixedEff->SetBinContent(11,(Eff->GetBinContent(10)+Eff->GetBinContent(12))/2);
            FixedEff->SetBinContent(14,(Eff->GetBinContent(13)+Eff->GetBinContent(15))/2);
        } else if(centString.CompareTo("10-20%") == 0){
            FixedEff->SetBinContent(10,(Eff->GetBinContent(9)+Eff->GetBinContent(11))/2);
            FixedEff->SetBinContent(13,(Eff->GetBinContent(12)+Eff->GetBinContent(16))/2);
            FixedEff->SetBinContent(14,(FixedEff->GetBinContent(13)+Eff->GetBinContent(16))/2);
            FixedEff->SetBinContent(15,(FixedEff->GetBinContent(14)+Eff->GetBinContent(16))/2);
            FixedEff->SetBinContent(17,(Eff->GetBinContent(16)+Eff->GetBinContent(18))/2);
        } else if(centString.CompareTo("20-30%") == 0){
            FixedEff->SetBinContent(14,(Eff->GetBinContent(13)+Eff->GetBinContent(15))/2);
            FixedEff->SetBinContent(16,(Eff->GetBinContent(15)+Eff->GetBinContent(17))/2);
        } else if(centString.CompareTo("20-40%") == 0){
            FixedEff->SetBinContent(12,(Eff->GetBinContent(11)+Eff->GetBinContent(13))/2);
            FixedEff->SetBinContent(14,(Eff->GetBinContent(13)+Eff->GetBinContent(15))/2);
            FixedEff->SetBinContent(16,(Eff->GetBinContent(15)+Eff->GetBinContent(17))/2);
        } else if(centString.CompareTo("20-50%") == 0){
            FixedEff->SetBinContent(12,(Eff->GetBinContent(11)+Eff->GetBinContent(13))/2);
            FixedEff->SetBinContent(14,(Eff->GetBinContent(13)+Eff->GetBinContent(15))/2);
            FixedEff->SetBinContent(16,(Eff->GetBinContent(15)+Eff->GetBinContent(17))/2);
        } else if(centString.CompareTo("30-50%") == 0){
            FixedEff->SetBinContent(11,(Eff->GetBinContent(10)+Eff->GetBinContent(12))/2);
            FixedEff->SetBinContent(16,(Eff->GetBinContent(15)+Eff->GetBinContent(17))/2);
        } else if(centString.CompareTo("40-60%") == 0){
            FixedEff->SetBinContent(8,(Eff->GetBinContent(7)+Eff->GetBinContent(9))/2);
            FixedEff->SetBinContent(12,(FixedEff->GetBinContent(13)+Eff->GetBinContent(16))/2);
            FixedEff->SetBinContent(13,(FixedEff->GetBinContent(12)+Eff->GetBinContent(16))/2);
            FixedEff->SetBinContent(14,(FixedEff->GetBinContent(13)+Eff->GetBinContent(16))/2);
            FixedEff->SetBinContent(15,(FixedEff->GetBinContent(14)+Eff->GetBinContent(16))/2);
            FixedEff->SetBinContent(17,(Eff->GetBinContent(16)+Eff->GetBinContent(18))/2);
        } else if(centString.CompareTo("50-80%") == 0){
            FixedEff->SetBinContent(11,(Eff->GetBinContent(10)+Eff->GetBinContent(12))/2);
            FixedEff->SetBinContent(10,(Eff->GetBinContent(9)+FixedEff->GetBinContent(11))/2);
            FixedEff->SetBinContent(11,(FixedEff->GetBinContent(10)+Eff->GetBinContent(12))/2);
            FixedEff->SetBinContent(15,(Eff->GetBinContent(14)+Eff->GetBinContent(16))/2);
        } else if(centString.CompareTo("60-80%") == 0){
            FixedEff->SetBinContent(11,(Eff->GetBinContent(10)+Eff->GetBinContent(12))/2);
            FixedEff->SetBinContent(10,(Eff->GetBinContent(9)+FixedEff->GetBinContent(11))/2);
            FixedEff->SetBinContent(11,(FixedEff->GetBinContent(10)+Eff->GetBinContent(12))/2);
            FixedEff->SetBinContent(10,(Eff->GetBinContent(9)+FixedEff->GetBinContent(11))/2);
        } else if(centString.CompareTo("40-80%") == 0){
            FixedEff->SetBinContent(11,(Eff->GetBinContent(10)+Eff->GetBinContent(12))/2);
            FixedEff->SetBinContent(14,(Eff->GetBinContent(13)+Eff->GetBinContent(15))/2);
        }
        return FixedEff;
    }

    // ****************************************************************************************************************
    // ****************************************************************************************************************
    // ****************************************************************************************************************
    Int_t CalcDeltaPtOverPt(    const TGraphAsymmErrors* gScaledPPStatErr, const TGraphAsymmErrors* gScaledPPSysErr,
                                const TGraphAsymmErrors* gPbPbStatErr, const TGraphAsymmErrors* gPbPbSysErr,
                                TGraphAsymmErrors* gResultStatErr, TGraphAsymmErrors* gResultSysErr,
                                bool bResultAtProtonProtonPt = true) {

        //
        // Klaus Reygers, 15-Sep-2013
        //

        //
        // Calculate relative pT shift between the Ncoll scaled pp invariant neutral pion yield and the Pb+Pb inv. yield.
        // For each Pb+Pb data points the pT,pp is calculated at which the saled pp inv. yield has the same value
        // as the Pb+Pb yield. The result of this routine is the normalized difference (pT - pT,pp)/pT,pp.
        // which cann be seen a measure of the fractional energy loss of a parton.
        // By default this result is given at pT,pp. In case of bResultAtProtonProtonPt == false it is given at the A+A pT.
        // See arXiv:1208.2254 for details.
        //
        // Klaus Reygers, 15-Sep-2013
        //
        Int_t iStatus                           = 0;

        // define Tsallis function to parameterize scaled pp invariant yield
        Double_t mass                           = 0.135; // pi0 mass
        TF1* fTsallis                           = new TF1("Tsallis", Form("[0] * pow(1.+(sqrt(x*x+%.10f*%.10f)-%.10f)/([1]*[2]), -[1])",
                                                                            mass, mass, mass),0.,30.);
        fTsallis->SetParNames("A","n","T_{Tsallis} (GeV/c)");
        fTsallis->SetParameters(1.,5.,0.18);

        //
        // fit scaled pp invariant yield
        //

        // copy scaled pp spectrum as fit would modify it (not allows due to 'const')
        TGraphAsymmErrors* gScaledPPStatErrCopy = dynamic_cast<TGraphAsymmErrors*>(gScaledPPStatErr->Clone());

        gScaledPPStatErrCopy->Fit(fTsallis, "0qN");
        Double_t prob                           = fTsallis->GetProb();
        Double_t p_A                            = fTsallis->GetParameter(0);
        Double_t p_n                            = fTsallis->GetParameter(1);
        Double_t p_T                            = fTsallis->GetParameter(2);
        delete gScaledPPStatErrCopy;

        // check fit status
        if (prob < 0.05) {
            cout << "CalcDeltaPtOverPt: Warning: bad fit of pp spectrum (p-Value = " << prob << ")" << endl;
            iStatus                             = 1;
        }

        // function used to determine value of pTppStar at which the pp yield (pT * dN/dpT = pT^2 * inv. yield)
        // takes a given value yieldPbPb
        TF1* fdNdPtScaledPP                     = new TF1(  "fdNdPtScaledPP", Form("[0] * x * x * pow(1.+(sqrt(x*x+%.10f*%.10f)-%.10f)/([1]*[2]), -[1]) - [3]",
                                                                mass, mass, mass),0.,30.);
        fdNdPtScaledPP->SetParNames("A","n","T_{Tsallis} (GeV/c)","yieldPbPb");

        // create root finder
        ROOT::Math::BrentRootFinder brf;

        //
        // loop over Pb+Pb data point and calculate pTppStar - pT for each Pb+Pb point with a given pT
        //
        for (Int_t i=0; i<gPbPbStatErr->GetN(); i++) {

            Double_t pTPbPb                     = gPbPbStatErr->GetX()[i];
            Double_t sfAA                       = pTPbPb * pTPbPb;
            Double_t yPbPb                      = sfAA * gPbPbStatErr->GetY()[i];  // yield pT * dN/dpT = pT^2 * inv. yield
            Double_t yPbPbStatErr               = sfAA * gPbPbStatErr->GetEYhigh()[i];
            Double_t yPbPbSysErr                = sfAA * gPbPbSysErr->GetEYhigh()[i];

            // determine relative errors of pp yields at the same pT
            Double_t pTpp                       = gScaledPPStatErr->GetX()[i];
            Double_t sfpp                       = pTpp * pTpp;
            Double_t ypp                        = sfpp * gScaledPPStatErr->GetY()[i]; // yield pT * dN/dpT = pT^2 * inv. yield
            Double_t yppStatErr                 = sfpp * gScaledPPStatErr->GetEYhigh()[i];
            Double_t yppRelStatErr              = yppStatErr/ypp;
            Double_t yppSysErr                  = sfpp * gScaledPPSysErr->GetEYhigh()[i];
            Double_t yppRelSysErr               = yppSysErr/ypp;

            // combine the error of the Pb+Pb yield and the pp yield at pT,pp assuming that one can neglect the
            // the change in the relative error of the pp yield between pT and pT,pp
            Double_t yCombinedStatErr           = sqrt(TMath::Power(yPbPbStatErr, 2) + TMath::Power(yPbPb * yppRelStatErr, 2));
            Double_t yCombinedSysErr            = sqrt(TMath::Power(yPbPbSysErr, 2) + TMath::Power(yPbPb * yppRelSysErr, 2));

            //
            // calculate the pTppStar at which the scaled pp spectrum has the value of the PbPb yield
            //

            // create wrapper function
            fdNdPtScaledPP->SetParameters(p_A, p_n, p_T, yPbPb);
            ROOT::Math::WrappedTF1 wf(*fdNdPtScaledPP);

            brf.SetFunction(wf, pTPbPb, 20.);
            brf.SetNpx(1000);
            brf.Solve();
            Double_t pTppStar                   = brf.Root(); // pT at which scaled pp yield takes value of Pb+Pb yield
            Double_t DeltaPt                    = pTppStar - pTPbPb;
            Double_t inverseSlope               = -1./fdNdPtScaledPP->Derivative(pTppStar);
            Double_t pTppStatErr                = inverseSlope * yCombinedStatErr;
            Double_t pTppSysErr                 = inverseSlope * yCombinedSysErr;

            // sanity check
            if (fabs(fdNdPtScaledPP->Eval(pTppStar)) > 1e-5) {
                cout << "CalcDeltaPtOverPt: Warning: root calculation incorrect: " << fabs(fdNdPtScaledPP->Eval(pTppStar)) << endl;
                iStatus                         = 1;
            }

            // determine at which pT the result is shown
            Double_t pTPlot                     = pTppStar;  // show result at p+p pT
            if (!bResultAtProtonProtonPt)
                pTPlot                          = pTPbPb;  // show result at A+A pT

            // result with statistical error
            gResultStatErr->GetX()[i]           = pTPlot;
            gResultStatErr->GetY()[i]           = DeltaPt/pTppStar;
            gResultStatErr->GetEXhigh()[i]      = 0;
            gResultStatErr->GetEXlow()[i]       = 0;
            gResultStatErr->GetEYhigh()[i]      = pTPbPb/(pTppStar*pTppStar) * pTppStatErr;  // Gaussian error propagation
            gResultStatErr->GetEYlow()[i]       = pTPbPb/(pTppStar*pTppStar) * pTppStatErr;   // Gaussian error propagation

            // result with systematic error
            gResultSysErr->GetX()[i]            = pTPlot;
            gResultSysErr->GetY()[i]            = DeltaPt/pTppStar;
            gResultSysErr->GetEXhigh()[i]       = 0.2;  // to make error boxes visible
            gResultSysErr->GetEXlow()[i]        = 0.2;  // to make error boxes visible
            gResultSysErr->GetEYhigh()[i]       = pTPbPb/(pTppStar*pTppStar) * pTppSysErr;  // Gaussian error propagation
            gResultSysErr->GetEYlow()[i]        = pTPbPb/(pTppStar*pTppStar) * pTppSysErr;   // Gaussian error propagation

        }

        delete fTsallis;
        delete fdNdPtScaledPP;

        return iStatus;

    }

    // *********************************************************************************************************
    // Create HEP data file for spectra/ RAA
    // *********************************************************************************************************
    void ProduceHEPDataFile(TGraphAsymmErrors* statErrors, TGraphAsymmErrors* sysErrors, TString header){
        cout << "*dataset: "<< endl;
        cout << header.Data() << endl;
        Double_t* currentX                  = statErrors->GetX();
        Double_t* currentXErrUp             = statErrors->GetEXhigh();
        Double_t* currentXErrDown           = statErrors->GetEXlow();
        Double_t* currentY                  = statErrors->GetY();
        Double_t* currentYErrStatDown       = statErrors->GetEYlow();
        Double_t* currentYErrSystDown       = sysErrors->GetEYlow();

        Int_t binStat                       = 0;
        Int_t binSyst                       = 0;
        Int_t binUpper                      = 0;

        for (Int_t i=0; i < statErrors->GetN(); i++){

            TString line                    = "";
            Bool_t statPresent              = kFALSE;
            Bool_t sysPresent               = kFALSE;
            //         cout << "entered loop" << endl;
            Int_t precision                 = 0;
            Int_t precisionSys              = 0;
            Int_t orderStat                 = 0;
            Int_t orderY                    = 0;
            if (currentYErrStatDown[binStat] > 1 )
                orderStat           = (Int_t)TMath::Log10(currentYErrStatDown[binStat]);
            else
                orderStat           = ((Int_t)TMath::Log10(currentYErrStatDown[binStat]))-1;
            if (currentY[binStat] > 1 )
                orderY              = (Int_t)TMath::Log10(currentY[binStat]);
            else
                orderY              = ((Int_t)TMath::Log10(currentY[binStat]))-1;
            precision               = orderY - orderStat;
            //                 cout << currentY[binStat] <<"\t" << (Int_t)TMath::Log10(currentY[binStat]) << "\t "<< currentYErrStatDown[binStat] << "\t"<< (Int_t)TMath::Log10(currentYErrStatDown[binStat])
            //                      << "\t"<< orderY-orderStat<< endl;
            statPresent             = kTRUE;
    //         cout << "stat finished" << endl;

            Int_t orderSys          = 0;
            if (currentYErrSystDown[binStat] > 1 )
                orderSys            = (Int_t)TMath::Log10(currentYErrSystDown[binStat]);
            else
                orderSys            = ((Int_t)TMath::Log10(currentYErrSystDown[binStat]))-1;

            precisionSys            = orderY - orderSys;
            sysPresent              = kTRUE;

            line                            = Form("%1.4f", currentX[i]);
            line                            = Form("%s (BIN= %1.1f TO %1.1f); \t",line.Data(), currentX[i]-currentXErrDown[i], currentX[i]+currentXErrUp[i] );
            if (statPresent){
                if (precision == 0){
                    line                    = Form("%s %2.2e +- %2.2e \t", line.Data(), currentY[binStat], currentYErrStatDown[binStat] );
                } else if (precision == 1 ){
                    line                    = Form("%s %2.3e +- %2.2e \t", line.Data(), currentY[binStat], currentYErrStatDown[binStat] );
                }  else if (precision == 2 ){
                    line                    = Form("%s %2.4e +- %2.2e \t", line.Data(), currentY[binStat], currentYErrStatDown[binStat] );
                }  else if (precision == 3 ){
                    line                    = Form("%s %2.5e +- %2.2e \t", line.Data(), currentY[binStat], currentYErrStatDown[binStat] );
                }  else if (precision == 4 ){
                    line                    = Form("%s %2.6e +- %2.2e \t", line.Data(), currentY[binStat], currentYErrStatDown[binStat] );
                } else {
                    line                    = Form("%s %2.2e +- %2.1e \t", line.Data(), currentY[binStat], currentYErrStatDown[binStat] );
                }
                binStat++;
                if (sysPresent){
                    if (precision-precisionSys == 0){
                        line                = Form("%s(DSYS= %2.2e);", line.Data(), currentYErrSystDown[binSyst] );
                    } else if (precision-precisionSys == 1){
                        line                = Form("%s(DSYS= %2.3e);", line.Data(), currentYErrSystDown[binSyst] );
                    } else if (precision-precisionSys == 2){
                        line                = Form("%s(DSYS= %2.4e);", line.Data(), currentYErrSystDown[binSyst] );
                    } else {
                        line                = Form("%s(DSYS= %2.2e);", line.Data(), currentYErrSystDown[binSyst] );
                    }
                    binSyst++;
                }
            }

            cout << line.Data() << endl;

//         for (Int_t i=0; i<statErrors->GetN(); i++){
//             cout << currentX[i] << " (BIN="<<currentX[i]-currentXErrDown[i] << " TO " << currentX[i]+currentXErrUp[i] << "); \t"
//                 << currentY[i] << " +- " <<currentYErrStatDown[i] << "\t (DSYS=" <<currentYErrSystDown[i] << ");" << endl;
//         }
        }
        cout << "*dataend:" << endl;

    }

    // *********************************************************************************************************
    // Create HEP data file for spectra with upper limits
    // *********************************************************************************************************
    void ProduceHEPDataFileWithUpperLimits(TGraphAsymmErrors* statErrors     = NULL,
                                        TGraphAsymmErrors* sysErrors     = NULL,
                                        TGraphAsymmErrors* upperLimits    = NULL,
                                        Double_t* pTBinning                = NULL,
                                        Int_t nPointsTot                    = 0,
                                        TString header                    = ""
                                        ){
        cout << "*dataset: "<< endl;
        cout << header.Data() << endl;
        Double_t* currentX                  = NULL;
        Double_t* currentY                  = NULL;
        Double_t* currentYErrStatDown       = NULL;
        Double_t* currentYErrSystDown       = NULL;
        Double_t* currentXSysErrors         = NULL;
        Double_t* currentXUpperLimits       = NULL;
        Double_t* currentYUpperLimits       = NULL;

        if (statErrors != NULL){
            currentX                        = statErrors->GetX();
            currentY                        = statErrors->GetY();
            currentYErrStatDown             = statErrors->GetEYlow();
        }
        if (sysErrors != NULL){
            currentYErrSystDown             = sysErrors->GetEYlow();
            currentXSysErrors               = sysErrors->GetX();
        }
        if (upperLimits != NULL){
            currentXUpperLimits             = upperLimits->GetX();
            currentYUpperLimits             = upperLimits->GetY();
        }
        TString clearanceLevel = "(95% CL)";
        Int_t binStat                       = 0;
        Int_t binSyst                       = 0;
        Int_t binUpper                      = 0;
        for (Int_t i=0; i < nPointsTot; i++){
            TString line                    = "";
            Bool_t statPresent              = kFALSE;
            Bool_t sysPresent               = kFALSE;
            Bool_t upperPresent             = kFALSE;
    //         cout << "entered loop" << endl;
            Int_t precision                 = 0;
            Int_t precisionSys              = 0;
            Int_t orderStat                 = 0;
            Int_t orderY                    = 0;
            if (statErrors != NULL){
                if ( currentX[binStat] > pTBinning[i] && currentX[binStat] < pTBinning[i+1] ) {
                    if (currentYErrStatDown[binStat] > 1 )
                        orderStat           = (Int_t)TMath::Log10(currentYErrStatDown[binStat]);
                    else
                        orderStat           = ((Int_t)TMath::Log10(currentYErrStatDown[binStat]))-1;
                    if (currentY[binStat] > 1 )
                        orderY              = (Int_t)TMath::Log10(currentY[binStat]);
                    else
                        orderY              = ((Int_t)TMath::Log10(currentY[binStat]))-1;
                    precision               = orderY - orderStat;
    //                 cout << currentY[binStat] <<"\t" << (Int_t)TMath::Log10(currentY[binStat]) << "\t "<< currentYErrStatDown[binStat] << "\t"<< (Int_t)TMath::Log10(currentYErrStatDown[binStat])
    //                      << "\t"<< orderY-orderStat<< endl;
                    statPresent             = kTRUE;
                }
            }
    //         cout << "stat finished" << endl;

            if (sysErrors != NULL){
                if ( currentXSysErrors[binSyst] > pTBinning[i] && currentXSysErrors[binSyst] < pTBinning[i+1] ) {
                    Int_t orderSys          = 0;
                    if (currentYErrSystDown[binStat] > 1 )
                        orderSys            = (Int_t)TMath::Log10(currentYErrSystDown[binStat]);
                    else
                        orderSys            = ((Int_t)TMath::Log10(currentYErrSystDown[binStat]))-1;

                    precisionSys            = orderY - orderSys;
                    sysPresent              = kTRUE;
                }
            }
    //         cout << "sys finished" << endl;
            if (upperLimits != NULL){
                if ( currentXUpperLimits[binUpper] > pTBinning[i] && currentXUpperLimits[binUpper] < pTBinning[i+1] )
                    upperPresent            = kTRUE;
            }
    //         cout << "here" << endl;

            line                            = Form("%s (BIN= %1.1f TO %1.1f); \t",line.Data(), pTBinning[i], pTBinning[i+1] );
            if (statPresent){
                if (precision == 0){
                    line                    = Form("%s %2.2e +- %2.2e \t", line.Data(), currentY[binStat], currentYErrStatDown[binStat] );
                } else if (precision == 1 ){
                    line                    = Form("%s %2.3e +- %2.2e \t", line.Data(), currentY[binStat], currentYErrStatDown[binStat] );
                }  else if (precision == 2 ){
                    line                    = Form("%s %2.4e +- %2.2e \t", line.Data(), currentY[binStat], currentYErrStatDown[binStat] );
                }  else if (precision == 3 ){
                    line                    = Form("%s %2.5e +- %2.2e \t", line.Data(), currentY[binStat], currentYErrStatDown[binStat] );
                }  else if (precision == 4 ){
                    line                    = Form("%s %2.6e +- %2.2e \t", line.Data(), currentY[binStat], currentYErrStatDown[binStat] );
                } else {
                    line                    = Form("%s %2.2e +- %2.1e \t", line.Data(), currentY[binStat], currentYErrStatDown[binStat] );
                }
                binStat++;
                if (sysPresent){
                    if (precision-precisionSys == 0){
                        line                = Form("%s(DSYS= %2.2e);", line.Data(), currentYErrSystDown[binSyst] );
                    } else if (precision-precisionSys == 1){
                        line                = Form("%s(DSYS= %2.3e);", line.Data(), currentYErrSystDown[binSyst] );
                    } else if (precision-precisionSys == 2){
                        line                = Form("%s(DSYS= %2.4e);", line.Data(), currentYErrSystDown[binSyst] );
                    } else {
                        line                = Form("%s(DSYS= %2.2e);", line.Data(), currentYErrSystDown[binSyst] );
                    }
                    binSyst++;
                } else if (upperPresent){
                    line                    = Form("%s(SYS = %2.2e %s);",line.Data(), currentYUpperLimits[binUpper], clearanceLevel.Data());
                    binUpper++;
                }
            } else if (upperPresent) {
                line                        = Form("%s %2.2e %s;", line.Data(), currentYUpperLimits[binUpper], clearanceLevel.Data());
                binUpper++;
            }

            cout << line.Data() << endl;
        }
        cout << "*dataend:" << endl;
    }

    // *********************************************************************************************************
    // Create HEP data file for spectra with upper limits errors split
    // *********************************************************************************************************
    void ProduceHEPDataFileWithUpperLimitsErrorsSplit(  TGraphAsymmErrors* statErrors     = NULL,
                                                        TGraphAsymmErrors* sysErrorsA     = NULL,
                                                        TGraphAsymmErrors* sysErrorsB     = NULL,
                                                        TGraphAsymmErrors* sysErrorsC     = NULL,
                                                        TGraphAsymmErrors* upperLimits    = NULL,
                                                        Double_t* pTBinning               = NULL,
                                                        Int_t nPointsTot                  = 0,
                                                        TString header                    = ""
                                                    ){
        cout << "*dataset: "<< endl;
        cout << header.Data() << endl;
        Double_t* currentX                  = NULL;
        Double_t* currentY                  = NULL;
        Double_t* currentYErrStatDown       = NULL;
        Double_t* currentYErrSystADown      = NULL;
        Double_t* currentYErrSystBDown      = NULL;
        Double_t* currentYErrSystCDown      = NULL;
        Double_t* currentXSysErrors         = NULL;
        Double_t* currentXUpperLimits       = NULL;
        Double_t* currentYUpperLimits       = NULL;

        if (statErrors != NULL){
            currentX                        = statErrors->GetX();
            currentY                        = statErrors->GetY();
            currentYErrStatDown             = statErrors->GetEYlow();
        }
        if (sysErrorsA != NULL){
            currentYErrSystADown            = sysErrorsA->GetEYlow();
            currentXSysErrors               = sysErrorsA->GetX();
        }
        if (sysErrorsB != NULL){
            currentYErrSystBDown            = sysErrorsB->GetEYlow();
        }
        if (sysErrorsC != NULL){
            currentYErrSystCDown            = sysErrorsC->GetEYlow();
        }
        if (upperLimits != NULL){
            currentXUpperLimits             = upperLimits->GetX();
            currentYUpperLimits             = upperLimits->GetY();
        }
        TString clearanceLevel = "(95% CL)";
        Int_t binStat                       = 0;
        Int_t binSyst                       = 0;
        Int_t binUpper                      = 0;
        for (Int_t i=0; i < nPointsTot; i++){
            TString line                    = "";
            Bool_t statPresent              = kFALSE;
            Bool_t sysPresent               = kFALSE;
            Bool_t upperPresent             = kFALSE;
    //         cout << "entered loop" << endl;

            Int_t precision                 = 0;
            Int_t precisionSysA             = 0;
            Int_t precisionSysB             = 0;
            Int_t precisionSysC             = 0;
            Int_t orderStat                 = 0;
            Int_t orderY                    = 0;
            if (statErrors != NULL){
                if ( currentX[binStat] > pTBinning[i] && currentX[binStat] < pTBinning[i+1] ) {
                    if (currentYErrStatDown[binStat] > 1 )
                        orderStat           = (Int_t)TMath::Log10(currentYErrStatDown[binStat]);
                    else
                        orderStat           = ((Int_t)TMath::Log10(currentYErrStatDown[binStat]))-1;
                    if (currentY[binStat] > 1 )
                        orderY              = (Int_t)TMath::Log10(currentY[binStat]);
                    else
                        orderY              = ((Int_t)TMath::Log10(currentY[binStat]))-1;
                    precision               = orderY - orderStat;
    //                 cout << currentY[binStat] <<"\t" << orderY << "\t "<< TMath::Log10(currentY[binStat]) << "\t "<< currentYErrStatDown[binStat] << "\t"<< orderStat
    //                      << "\t"<< orderY-orderStat<< endl;
                    statPresent             = kTRUE;
                }
            }
    //         cout << "stat finished" << endl;

            if (sysErrorsA != NULL && sysErrorsB != NULL && sysErrorsC != NULL){
                if ( currentXSysErrors[binSyst] > pTBinning[i] && currentXSysErrors[binSyst] < pTBinning[i+1] ) {
                    Int_t orderSysA         = 0;
                    if (currentYErrSystADown[binStat] > 1 )
                        orderSysA           = (Int_t)TMath::Log10(currentYErrSystADown[binSyst]);
                    else
                        orderSysA           = ((Int_t)TMath::Log10(currentYErrSystADown[binSyst]))-1;
                    precisionSysA           = orderY - orderSysA;

                    Int_t orderSysB         = 0;
                    if (currentYErrSystBDown[binStat] > 1 )
                        orderSysB           = (Int_t)TMath::Log10(currentYErrSystBDown[binSyst]);
                    else
                        orderSysB           = ((Int_t)TMath::Log10(currentYErrSystBDown[binSyst]))-1;
                    precisionSysB           = orderY - orderSysB;

                    Int_t orderSysC         = 0;
                    if (currentYErrSystCDown[binStat] > 1 )
                        orderSysC           = (Int_t)TMath::Log10(currentYErrSystCDown[binSyst]);
                    else
                        orderSysC           = ((Int_t)TMath::Log10(currentYErrSystCDown[binSyst]))-1;

                    precisionSysC           = orderY - orderSysC;
                    sysPresent              = kTRUE;
                }
            }
    //         cout << "sys finished" << endl;
            if (upperLimits != NULL){
                if ( currentXUpperLimits[binUpper] > pTBinning[i] && currentXUpperLimits[binUpper] < pTBinning[i+1] )
                    upperPresent            = kTRUE;
            }
    //         cout << "here" << endl;

            line                            = Form("%s (BIN= %1.1f TO %1.1f); \t",line.Data(), pTBinning[i], pTBinning[i+1] );
            if (statPresent){
                if (precision == 0){
                    line                    = Form("%s %2.2e +- %2.2e \t", line.Data(), currentY[binStat], currentYErrStatDown[binStat] );
                } else if (precision == 1 ){
                    line                    = Form("%s %2.3e +- %2.2e \t", line.Data(), currentY[binStat], currentYErrStatDown[binStat] );
                }  else if (precision == 2 ){
                    line                    = Form("%s %2.4e +- %2.2e \t", line.Data(), currentY[binStat], currentYErrStatDown[binStat] );
                }  else if (precision == 3 ){
                    line                    = Form("%s %2.5e +- %2.2e \t", line.Data(), currentY[binStat], currentYErrStatDown[binStat] );
                }  else if (precision == 4 ){
                    line                    = Form("%s %2.6e +- %2.2e \t", line.Data(), currentY[binStat], currentYErrStatDown[binStat] );
                } else {
                    line                    = Form("%s %2.2e +- %2.1e \t", line.Data(), currentY[binStat], currentYErrStatDown[binStat] );
                }
                binStat++;
                if (sysPresent){
                    if (precision-precisionSysA == 0){
                        line                = Form("%s(DSYS= %2.2e:A, ", line.Data(), currentYErrSystADown[binSyst] );
                    } else if (precision-precisionSysA == 1){
                        line                = Form("%s(DSYS= %2.3e:A, ", line.Data(), currentYErrSystADown[binSyst] );
                    } else if (precision-precisionSysA == 2){
                        line                = Form("%s(DSYS= %2.4e:A, ", line.Data(), currentYErrSystADown[binSyst] );
                    } else {
                        line                = Form("%s(DSYS= %2.2e:A, ", line.Data(), currentYErrSystADown[binSyst] );
                    }
                    if (precision-precisionSysB == 0){
                        line                = Form("%sDSYS= %2.2e:B, ", line.Data(), currentYErrSystBDown[binSyst] );
                    } else if (precision-precisionSysB == 1){
                        line                = Form("%sDSYS= %2.3e:B, ", line.Data(), currentYErrSystBDown[binSyst] );
                    } else if (precision-precisionSysB == 2){
                        line                = Form("%sDSYS= %2.4e:B, ", line.Data(), currentYErrSystBDown[binSyst] );
                    } else {
                        line                = Form("%sDSYS= %2.2e:B, ", line.Data(), currentYErrSystBDown[binSyst] );
                    }
                    if (precision-precisionSysC == 0){
                        line                = Form("%sDSYS= %2.2e:C);", line.Data(), currentYErrSystCDown[binSyst] );
                    } else if (precision-precisionSysC == 1){
                        line                = Form("%sDSYS= %2.3e:C);", line.Data(), currentYErrSystCDown[binSyst] );
                    } else if (precision-precisionSysC == 2){
                        line                = Form("%sDSYS= %2.4e:C);", line.Data(), currentYErrSystCDown[binSyst] );
                    } else {
                        line                = Form("%sDSYS= %2.2e:C);", line.Data(), currentYErrSystCDown[binSyst] );
                    }
                    binSyst++;
                } else if (upperPresent){
                    line                    = Form("%s(SYS = %2.2e %s);",line.Data(), currentYUpperLimits[binUpper], clearanceLevel.Data());
                    binUpper++;
                }
            } else if (upperPresent) {
                line                        = Form("%s %2.2e %s;", line.Data(), currentYUpperLimits[binUpper], clearanceLevel.Data());
                binUpper++;
            }
            cout << line.Data() << endl;
        }
        cout << "*dataend:" << endl;
    }

    // *********************************************************************************************************
    // Create HEP data file for NPart graphs
    // *********************************************************************************************************
    void ProduceHEPDataFileNPart(TGraphErrors* statErrors, TGraphErrors* sysErrors, TString header){
        cout << "*dataset: "<< endl;
        cout << header.Data() << endl;
        Double_t* currentX                  = statErrors->GetX();
        Double_t* currentXErr               = sysErrors->GetEX();
        Double_t* currentY                  = statErrors->GetY();
        Double_t* currentYErrStat           = statErrors->GetEY();
        Double_t* currentYErrSyst           = sysErrors->GetEY();

        for (Int_t i=0; i<statErrors->GetN(); i++){
            cout << currentX[i] << " (BIN="<<currentX[i]-currentXErr[i] << " TO " << currentX[i]+currentXErr[i] << "); \t" << currentY[i] << "+-" <<currentYErrStat[i]
                << "\t (DSYS=" <<currentYErrSyst[i] << ");" << endl;
        }
        cout << "*dataend:" << endl;
    }

    // *********************************************************************************************************
    // this function calculates the direct photon points and decides, whether
    // the point is a point or just an upper limit
    // *********************************************************************************************************
    TGraphAsymmErrors *CalculateDirectPhotonPointsAndUpperLimits(TH1* error, TH1* spectrum, float iFlag, float ScaleArrowlength, Int_t offset = 0){

        const int kMaxPoints                = 50;
        int nPoints                         = 0;
        float x[kMaxPoints],xel[kMaxPoints],xeh[kMaxPoints];
        float y[kMaxPoints],yel[kMaxPoints],yeh[kMaxPoints];

        ScaleArrowlength                    = ScaleArrowlength*1.3;
        //  const float  ScaleArrowlength = 0.5;
        //  const float  ScaleArrowlength = 0.5;

        for(int ib = 1+offset;ib <= error->GetNbinsX();ib++){
            float val                       = spectrum->GetBinContent(ib);
            float pt                        = spectrum->GetBinCenter(ib);
            float err                       = val - error->GetBinContent(ib);
            float errt                      = error->GetBinContent(ib);
            float errX                      = spectrum->GetBinWidth(ib);

            float confi                     = 1.645; // 95%
            //float confi                   = 1.28; // 90%
            //float lowbound                = val - 1.28*abs(err); // 90% CL
            float lowbound                  = val - confi*abs(err); // 95% CL
            //cout << "pt " << pt << " val " << val << " errt " << errt << " err " << err << " lowbound " << lowbound << endl;
            if(iFlag==0){ // points and errors > 0 syst
                if(val>0&&errt>0){
                    x[nPoints]              = pt;
                    xel[nPoints]            = errX/2;
                    xeh[nPoints]            = errX/2;
                    y[nPoints]              = val;
                    yel[nPoints]            = abs(err);
                    yeh[nPoints]            = abs(err);
                    nPoints++;
                }
            } else if(iFlag==1){ // points and errors > 0 stat
                if(val>0&&errt>0){
                    x[nPoints]              = pt;
                    xel[nPoints]            = 0;
                    xeh[nPoints]            = 0;
                    y[nPoints]              = val;
                    yel[nPoints]            = abs(err);
                    yeh[nPoints]            = abs(err);
                    nPoints++;
                }
            } else if(iFlag==2){ // points and errors and confi > 0 syst
                if(val>0&&errt>0&&lowbound>=0){
                    x[nPoints]              = pt;
                    xel[nPoints]            = errX/2;
                    xeh[nPoints]            = errX/2;
                    y[nPoints]              = val;
                    yel[nPoints]            = abs(err);
                    yeh[nPoints]            = abs(err);
                    nPoints++;
                }
            } else if(iFlag==3){ // points and errors and confi > 0 stat
                if(val>0&&errt>0&&lowbound>=0){
                    x[nPoints]              = pt;
                    xel[nPoints]            = 0;
                    xeh[nPoints]            = 0;
                    y[nPoints]              = val;
                    yel[nPoints]            = abs(err);
                    yeh[nPoints]            = abs(err);
                    nPoints++;
                }
            } else if(iFlag==4){ // ALL Points Processed
                // upperlimit for points but errors consistent with zero
                if(errt<=0){
                    x[nPoints]              = pt;
                    xel[nPoints]            = 0;
                    xeh[nPoints]            = 0;
                    y[nPoints]              = val;
                    yel[nPoints]            = 0;
                    yeh[nPoints]            = confi*abs(err);
                    nPoints++;
                }
            } else if(iFlag==5){ // arrow for points with zero
                if(errt<=0){
                    x[nPoints]              = pt;
                    xel[nPoints]            = 0;
                    xeh[nPoints]            = 0;
                    y[nPoints]              = confi*abs(err)+val;
                    //yel[nPoints]          = (1.575 * abs(err) + ScaleArrowlength * val);
                    //yel[nPoints]          = confi*abs(err)+(val*ScaleArrowlength);
                    yel[nPoints]            = abs(y[nPoints]-(0.08*y[nPoints]));
                    cout<<pt<<"  "<<yel[nPoints]<<"  "<<y[nPoints]<<endl;
                    yeh[nPoints]            = 0;
                    nPoints++;
                }
            } else if(iFlag==6){ // upperlimit for conf
                if(lowbound<=0){
                    x[nPoints]              = pt;
                    xel[nPoints]            = 0;
                    xeh[nPoints]            = 0;
                    y[nPoints]              = val;
                    yel[nPoints]            = 0;
                    yeh[nPoints]            = confi*abs(err);
                    nPoints++;
                }
            } else if(iFlag==7){ // arrow for confi
                if(lowbound<=0){
                    x[nPoints]              = pt;
                    xel[nPoints]            = 0;
                    xeh[nPoints]            = 0;
                    y[nPoints]              = confi*abs(err)+val;
                    //yel[nPoints]          = (1.575 * abs(err) + ScaleArrowlength * val);
                    //yel[nPoints]          = confi*abs(err)+(val*ScaleArrowlength);
                    yel[nPoints]            = abs(y[nPoints]-(0.08*y[nPoints]));
                    yeh[nPoints]            = 0;
                    nPoints++;
                }
            }
        }

        if(nPoints > 0){
            TGraphAsymmErrors *gr           = new TGraphAsymmErrors(nPoints,x,y,xel,xeh,yel,yeh);
            return gr;
        }

        return 0;
    }

    // ****************************************************************************************************************
    // ************************** Blast wave fit generalized **********************************************************
    // ****************************************************************************************************************
    Double_t BWdndpt(Double_t *x, Double_t *par) {

        Double_t pt     = x[0]*1000.; //GeV->MeV
        Double_t t      = par[1];
        Double_t beta   = abs(par[2]);
        Double_t yt     = 0.5*TMath::Log((1+beta)/(1-beta)); // atanh(beta)
        Double_t m0     = 0.; //pi0
        Double_t mt     = TMath::Sqrt(m0*m0+pt*pt);
        Double_t f      = par[0]*pt*mt*TMath::BesselI(0,pt*sinh(yt)/t)*TMath::BesselK(1,TMath::Abs(mt*cosh(yt)/t));

        return f;
    }

    // ****************************************************************************************************************
    // ****************************************************************************************************************
    // ****************************************************************************************************************
    Double_t FitSpectrum(Double_t *x, Double_t *p)
    {
        // Double_t hagd = p[0]*(p[2]-1)*(p[2]-2)/p[1]/p[2]/TMath::TwoPi()/TMath::Power((1.+x[0]/p[1]/p[2]),p[2]);
        Double_t hagd   = p[0]/TMath::TwoPi()/TMath::Power((1.+x[0]/p[1]/p[2]),p[2]);
        Double_t expo   = p[3]*TMath::Exp(-x[0]/p[4]);
        return hagd+expo;
    }

    // ****************************************************************************************************************
    // ****************************************************************************************************************
    // ****************************************************************************************************************
    Double_t Hagedorn(Double_t *x, Double_t *p)
    {
        Double_t hagd   = p[0]/TMath::TwoPi()/TMath::Power((1.+x[0]/p[1]/p[2]),p[2]);
        return hagd;
    }

    // ****************************************************************************************************************
    // ****************************************************************************************************************
    // ****************************************************************************************************************
    Double_t Exponent(Double_t *x, Double_t *p)
    {
        Double_t expo   = p[0]*TMath::Exp(-x[0]/p[1]);
        return expo;
    }

    // *********************************************************************************************************
    // Calculates the summed error of two graphs (errors added quadratically)
    // *********************************************************************************************************
    TGraphAsymmErrors* AddErrorsOfGraphsQuadratically (TGraphAsymmErrors* graphA, TGraphAsymmErrors* graphB){
        TGraphAsymmErrors* graphReturn      = (TGraphAsymmErrors*) graphA->Clone("summedErrors");
        for (Int_t i = 0; i < graphReturn->GetN(); i++){
            Double_t err1H                  = graphA->GetErrorYhigh(i)/graphA->GetY()[i];
            Double_t err2H                  = graphB->GetErrorYhigh(i)/graphB->GetY()[i];
            Double_t err1L                  = graphA->GetErrorYlow(i)/graphA->GetY()[i];
            Double_t err2L                  = graphB->GetErrorYlow(i)/graphB->GetY()[i];
            Double_t newErrorH              = TMath::Sqrt(err1H*err1H+err2H*err2H)*graphReturn->GetY()[i];
            Double_t newErrorL              = TMath::Sqrt(err1L*err1L+err2L*err2L)*graphReturn->GetY()[i];
            graphReturn->SetPointEYhigh(i, newErrorH);
            graphReturn->SetPointEYlow(i, newErrorL);
        }
        return graphReturn;
    }

    // *********************************************************************************************************
    // Calculates the summed error of 3 graphs (errors added quadratically)
    // *********************************************************************************************************
    TGraphAsymmErrors* Add3ErrorsOfGraphsQuadratically (TGraphAsymmErrors* graphA, TGraphAsymmErrors* graphB, TGraphAsymmErrors* graphC){
        TGraphAsymmErrors* graphReturn      = (TGraphAsymmErrors*) graphA->Clone("summedErrors");
        for (Int_t i = 0; i < graphReturn->GetN(); i++){
            Double_t err1H                  = graphA->GetErrorYhigh(i)/graphA->GetY()[i];
            Double_t err2H                  = graphB->GetErrorYhigh(i)/graphB->GetY()[i];
            Double_t err3H                  = graphC->GetErrorYhigh(i)/graphC->GetY()[i];
            Double_t err1L                  = graphA->GetErrorYlow(i)/graphA->GetY()[i];
            Double_t err2L                  = graphB->GetErrorYlow(i)/graphB->GetY()[i];
            Double_t err3L                  = graphC->GetErrorYlow(i)/graphC->GetY()[i];
            Double_t newErrorH              = TMath::Sqrt(err1H*err1H+err2H*err2H+err3H*err3H)*graphReturn->GetY()[i];
            Double_t newErrorL              = TMath::Sqrt(err1L*err1L+err2L*err2L+err3L*err3L)*graphReturn->GetY()[i];
            graphReturn->SetPointEYhigh(i, newErrorH);
            graphReturn->SetPointEYlow(i, newErrorL);
        }
        return graphReturn;
    }

    //*************************************************************************************************************
    //*******************************Convert p Value to nSigma ****************************************************
    //*************************************************************************************************************
    Double_t PValueToNSigma(Double_t pvalue) {

        TF1 f("f","1-TMath::Erf(x/sqrt(2.)) - [0]", -100, 100);
        f.SetParameter(0, pvalue);

        // create wrapper function
        ROOT::Math::WrappedTF1 wfrf(f);

        // create root finder
        ROOT::Math::BrentRootFinder brf;

        // set parameters of the method
        brf.SetFunction(wfrf, 0., 50.);
        brf.Solve();
        Double_t nsigma                     = brf.Root();

        return nsigma;
    }

    //*************************************************************************************************************
    // scale data points in a given pT interval according to type B error
    // in an anti-correlated manner
    //*************************************************************************************************************
    Double_t scaleValuesWithBError( Double_t& eps, Double_t& pt, Double_t pt_min = 1., Double_t pt_max = 2. ) {

        Double_t sf                         = 0;

        if (pt < pt_min) sf                 = -1;
        else if (pt > pt_max) sf            = 1;
        else sf                             = 2./(pt_max - pt_min) * (pt - (pt_max + pt_min)/2.);

        return eps * sf;
    }

    //*************************************************************************************************************
    // calculate chi2 of a set points w.r.t. the null hypothesis R=1
    //*************************************************************************************************************
    Double_t Chi2ForNullHypo(TGraphErrors* g) {

        Double_t R_true                     = 1.0;
        Double_t sumsq                      = 0;

        for (Int_t i=0; i<g->GetN(); i++) {
            Double_t y                      = g->GetY()[i];
            Double_t y_sig                  = g->GetEY()[i];
            Double_t y_sig_rel              = y_sig/y;
            Double_t nsig                   = (y-R_true)/(R_true * y_sig_rel);
            sumsq                           += nsig*nsig;
        }
        return sumsq;
    }

    //*************************************************************************************************************
    // Fill chi2 histogram for null hypothesis
    //*************************************************************************************************************
    void FillChi2HistForNullHypo    (  Int_t    n_pseudo_exp,
                                    TH1F*    &histo,
                                    TGraph*  &g_rel_stat_plus_type_a_error,
                                    TGraph*  &g_rel_type_b_error,
                                    Double_t &rel_type_c_error,
                                    Bool_t   anti_corr_type_b) {

        // null hypothesis: double ratio R = 1
        Double_t R_true                     = 1.0;
        // random numbers
        TRandom rndm;

        for (Int_t i_pseudo_exp=0; i_pseudo_exp<n_pseudo_exp; i_pseudo_exp++) {
            Double_t sumsq                  = 0;
            Double_t eps_b                  = rndm.Gaus(0,1);
            Double_t eps_c                  = rndm.Gaus(0,1);

            // create pseudo data set
            for (Int_t ip=0; ip<g_rel_type_b_error->GetN(); ip++) {
                Double_t rel_type_b_error   = g_rel_type_b_error->GetY()[ip];
                Double_t R_mod;
                if(anti_corr_type_b) {
                    Double_t pt             = g_rel_stat_plus_type_a_error->GetX()[ip];
                    R_mod                   = R_true * (1. + scaleValuesWithBError(eps_b,pt) * rel_type_b_error) * (1. + eps_c * rel_type_c_error);
                } else {
                    R_mod                   = R_true * (1. + eps_b * rel_type_b_error) * (1. + eps_c * rel_type_c_error);
                }

                Double_t rel_stat_plus_type_a_err           = g_rel_stat_plus_type_a_error->GetY()[ip];
                Double_t abs_stat_plus_type_a_err_scaled    = R_mod * rel_stat_plus_type_a_err;

                Double_t y                  = rndm.Gaus(R_mod, abs_stat_plus_type_a_err_scaled);
                Double_t nsig               = (y - R_true)/(R_true * rel_stat_plus_type_a_err);
                sumsq                       += nsig*nsig;

            }
            histo->Fill(sumsq);
        }
    }

    //*************************************************************************************************************
    //************************ Extract systematic errors from total and statistical errors ************************
    //*************************************************************************************************************
    void ExtractSystematicFromTotal(Int_t nPoints, Double_t* totals, Double_t* stat, Double_t* sysErr ){
        for (Int_t i = 0; i < nPoints; i++){
            sysErr[i] = TMath::Sqrt(totals[i]*totals[i]-stat[i]*stat[i]);
        }
        return ;
    }

    // ****************************************************************************************************************
    // ****************************************************************************************************************
    // ****************************************************************************************************************
    TF1* CalculateRatioOfTwoFunctions (TF1* fit1, TF1* fit2, TString name){
        Int_t nParFunc1  = fit1->GetNpar();
        Int_t nParFunc2  = fit2->GetNpar();
        TString formula1    = fit1->GetExpFormula();
        TString formula2    = fit2->GetExpFormula();

        for (Int_t i = 0; i< nParFunc2; i++){
            formula2.ReplaceAll(Form("[%d]",i), Form("[%d]",i+nParFunc1));
        }
        TF1* newFunction = new TF1(name.Data(),Form("(%s)/(%s)",formula1.Data(), formula2.Data()));
        for (Int_t i = 0; i < nParFunc1; i++ ){
            newFunction->SetParameter(i, fit1->GetParameter(i));
        }
        for (Int_t j = 0; j < nParFunc2; j++ ){
            newFunction->SetParameter(nParFunc1+j, fit2->GetParameter(j));
        }
        return newFunction;
    }

    // ****************************************************************************************************************
    // ****************************************************************************************************************
    // ****************************************************************************************************************
    TString AutoDetectMainTList(Int_t mode , TFile* fFile){
        TKey *key;
        TIter next(fFile->GetListOfKeys());
        TString mainDir = "";
        while ((key=(TKey*)next())){
            cout << Form("-> found TopDir: %s",key->GetName());
            mainDir = key->GetName();
        }
        cout << endl;
        TString nominalMainDir     = "";
        if (mode == 9 || mode == 0)
            nominalMainDir         = "GammaConvV1";
        else if( mode == 1 )
            nominalMainDir         = "GammaConvDalitzV1";
        else if (mode == 2 || mode == 3 || mode == 13)
            nominalMainDir         = "GammaConvCalo";
        else if (mode == 4 || mode == 12 || mode == 5)
            nominalMainDir         = "GammaCalo";
        else if( mode == 6 || mode == 7 )
            nominalMainDir         = "GammaConvDalitzCalo";
        else if (mode == 10 || mode == 11 )
            nominalMainDir         = "GammaCaloMerged";
        else if (mode == 30 )
            nominalMainDir         = "GammaConvV1";
        else if (mode == 40 || mode == 41 || mode == 42 || mode == 43|| mode == 44 || mode == 45 ||
                 mode == 46 || mode == 47 || mode == 48 || mode == 49|| mode == 50)
            nominalMainDir         = "GammaConvNeutralMesonPiPlPiMiPiZero";
        TObjArray *arr;
        arr             = mainDir.Tokenize("_");
        TString start   = ((TObjString*)arr->At(0))->GetString();;
        if (start.CompareTo(nominalMainDir.Data()) == 0)
            return mainDir;
        else
            return "";
    }

    // ****************************************************************************************************************
    // ****************************************************************************************************************
    // ****************************************************************************************************************
    TString GetDefaultMainTListName(Int_t mode){
        TString nominalMainDir     = "";
        if (mode == 9 || mode == 0)
            nominalMainDir         = "GammaConvV1";
        else if( mode == 1 )
            nominalMainDir         = "GammaConvDalitzV1";
        else if (mode == 2 || mode == 3 || mode == 13)
            nominalMainDir         = "GammaConvCalo";
        else if (mode == 4 || mode == 12 || mode == 5)
            nominalMainDir         = "GammaCalo";
        else if( mode == 6 || mode == 7 )
            nominalMainDir         = "GammaConvDalitzCalo";
        else if (mode == 10 || mode == 11 )
            nominalMainDir         = "GammaCaloMerged";
        else if (mode == 30 )
            nominalMainDir         = "GammaConvV1";
        else if (mode == 40 || mode == 41 || mode == 42 || mode == 43|| mode == 44 || mode == 45 ||
                 mode == 46 || mode == 47 || mode == 48 || mode == 49|| mode == 50)
            nominalMainDir         = "GammaConvNeutralMesonPiPlPiMiPiZero";
        return nominalMainDir;
    }

    // ****************************************************************************************************************
    // ******************** Add errors in quadrature for TGraphsAsymmErrors *******************************************
    // ****************************************************************************************************************
    TGraphAsymmErrors* AddErrorsQuadraticallyTGraph(TGraphAsymmErrors* graphA, TGraphAsymmErrors* graphB){
        Double_t* xValue            = graphA->GetX();
        Double_t* yValue            = graphA->GetY();
        Double_t* xErrorHigh        = graphA->GetEXhigh();
        Double_t* yErrorHigh        = graphA->GetEYhigh();
        Double_t* xErrorLow         = graphA->GetEXlow();
        Double_t* yErrorLow         = graphA->GetEYlow();
        Int_t nPoints               = graphA->GetN();

        TGraphAsymmErrors* graphR   = new TGraphAsymmErrors(nPoints,xValue,yValue,xErrorLow,xErrorHigh,yErrorLow, yErrorHigh);

        Double_t* yErrorRHigh       = graphR->GetEYhigh();
        Double_t* yErrorRLow        = graphR->GetEYlow();
        Double_t* yErrorSysHigh     = graphB->GetEYhigh();
        Double_t* yErrorSysLow      = graphB->GetEYlow();

        for(Int_t i=0; i<nPoints; i++){
            yErrorRHigh[i]          = TMath::Sqrt(TMath::Power(yErrorHigh[i],2)+TMath::Power(yErrorSysHigh[i],2));
            yErrorRLow[i]           = TMath::Sqrt(TMath::Power(yErrorLow[i],2)+TMath::Power(yErrorSysLow[i],2));
        }
        return graphR;
    }

    // ****************************************************************************************************************
    // ******************** Add errors in quadrature for TGraphErrors *************************************************
    // ****************************************************************************************************************
    TGraphErrors* AddErrorsQuadraticallyTGraph(TGraphErrors* graphA, TGraphErrors* graphB){
        Double_t* xValue            = graphA->GetX();
        Double_t* yValue            = graphA->GetY();
        Double_t* xError            = graphA->GetEX();
        Double_t* yError            = graphA->GetEY();
        Int_t nPoints               = graphA->GetN();

        TGraphErrors* graphR        = new TGraphErrors(nPoints,xValue,yValue,xError,yError);

        Double_t* yErrorR           = graphR->GetEY();
        Double_t* yErrorB           = graphB->GetEY();

        for(Int_t i=0; i<nPoints; i++){
            yErrorR[i]              = TMath::Sqrt(TMath::Power(yError[i],2)+TMath::Power(yErrorB[i],2));
        }
        return graphR;
    }

#endif
