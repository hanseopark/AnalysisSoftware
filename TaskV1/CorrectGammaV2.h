//*********************************************************************************************************
//*** Conversion method                                                                          **********
//*** Correct inclusive raw gamma candite spectrum from conversions with bin-by-bin efficiencies **********
// - scale with 1/nEvt first
// - subtract secodary gammas (MC approach)
// - take out impurities
// - reco efficiency (histoRecoEff)
// - conversion probability (histoConvProb)
// - scale with 1/deltaEta, scaling, 1/nEvt
// - convert to inv yield (scale by bin 1/bin center)
//*********************************************************************************************************
void CorrectGammaEffiResol(TH1D*       histoGammaCorr,
                           TH1D*       histoAllSec,
                           TH1D*       histoK0sSec,
                           TH1D*       histoPurity,
                           TH1D*       histoConvProb,
                           TH1D*       histoRecoEff,
                           Double_t    deltaEta,
                           Double_t    scaling,
                           Double_t    nEvt
                           ){
    // scale with 1/nEvt as secondaries are normalized per event
    histoGammaCorr->Scale(1./nEvt);
    // subtract secondary gammas
    histoGammaCorr->Add(histoAllSec,-1);
    histoGammaCorr->Add(histoK0sSec,-1);
    // multiply with purity for primary particles
    histoGammaCorr->Multiply(histoGammaCorr,histoPurity,1.,1.,"");
    // divide by reconstruction efficiency
    histoGammaCorr->Divide(histoGammaCorr,histoRecoEff,1.,1.,"");
    // divide by P_{conv}
    histoGammaCorr->Divide(histoGammaCorr,histoConvProb,1.,1.,"");
    // scale to with 1/deltaEta and 1/deltaPhi
    histoGammaCorr->Scale(1./deltaEta);
    histoGammaCorr->Scale(scaling);
    // convert to inv yield by dividing by pT at bin center
    for (Int_t i = 1; i < histoGammaCorr->GetNbinsX()+1 ; i++){
        Double_t newBinContent  = histoGammaCorr->GetBinContent(i)  / histoGammaCorr->GetBinCenter(i);
        Double_t newBinError    = histoGammaCorr->GetBinError(i)    / histoGammaCorr->GetBinCenter(i);
        histoGammaCorr->SetBinContent(i,newBinContent);
        histoGammaCorr->SetBinError(i,  newBinError);
    }
}

//*********************************************************************************************************
//*** Conversion method                                                                          **********
//*** Correct inclusive raw gamma candite spectrum from conversions with bin-by-bin efficiencies **********
// - scale with 1/nEvt first
// - subtract secodary gammas (cocktail approach)
// - take out impurities
// - reco efficiency (histoRecoEff)
// - conversion probability (histoConvProb)
// - scale with 1/deltaEta, scaling, 1/nEvt
// - convert to inv yield (scale by bin 1/bin center)
//*********************************************************************************************************
void CorrectGammaEffiResolCocktail(TH1D*       histoGammaCorr,
                                   TH1D**      histoHadronSec,
                                   TH1D*       histoRestSec,
                                   TH1D*       histoPurity,
                                   TH1D*       histoConvProb,
                                   TH1D*       histoRecoEff,
                                   Double_t    deltaEta,
                                   Double_t    scaling,
                                   Double_t    nEvt
                                   ){
    // scale with 1/nEvt as secondaries are normalized per event
    histoGammaCorr->Scale(1./nEvt);
    // subtract secondary gammas from decays
    for (Int_t k = 0; k < 3; k++){
        if (histoHadronSec[k]) histoGammaCorr->Add(histoHadronSec[k],-1);
    }
    // subtract secondary gammas from material interactions
    histoGammaCorr->Add(histoRestSec,-1);
    // multiply with purity for primary particles
    histoGammaCorr->Multiply(histoGammaCorr,histoPurity,1.,1.,"");
    // divide by reconstruction efficiency
    histoGammaCorr->Divide(histoGammaCorr,histoRecoEff,1.,1.,"");
    // divide by P_{conv}
    histoGammaCorr->Divide(histoGammaCorr,histoConvProb,1.,1.,"");
    // scale to with 1/deltaEta and 1/deltaPhi
    histoGammaCorr->Scale(1./deltaEta);
    histoGammaCorr->Scale(scaling);
    // convert to inv yield by dividing by pT at bin center
    for (Int_t i = 1; i < histoGammaCorr->GetNbinsX()+1 ; i++){
        Double_t newBinContent  = histoGammaCorr->GetBinContent(i)  / histoGammaCorr->GetBinCenter(i);
        Double_t newBinError    = histoGammaCorr->GetBinError(i)    / histoGammaCorr->GetBinCenter(i);
        histoGammaCorr->SetBinContent(i,newBinContent);
        histoGammaCorr->SetBinError(i,  newBinError);
    }
}

//*********************************************************************************************************
//*** Calorimeter method                                                                         **********
//*** Correct inclusive raw gamma candite spectrum with bin-by-bin efficiencies                  **********
// - scale with 1/nEvt first
// - subtract secodary gammas (MC approach)
// - take out impurities
// - reco efficiency (histoRecoEff)
// - conversion probability (histoConvProb)
// - scale with 1/deltaEta, scaling, 1/nEvt
// - convert to inv yield (scale by bin 1/bin center)
//*********************************************************************************************************
void CorrectGammaEffiResol(TH1D*       histoGammaCorr,
                           TH1D*       histoAllSec,
                           TH1D*       histoK0sSec,
                           TH1D*       histoPurity,
                           TH1D*       histoRecoEff,
                           Double_t    deltaEta,
                           Double_t    scaling,
                           Double_t    nEvt
                           ){
    // scale with 1/nEvt as secondaries are normalized per event
    histoGammaCorr->Scale(1./nEvt);
    // subtract secondary gammas
    histoGammaCorr->Add(histoAllSec,-1);
    histoGammaCorr->Add(histoK0sSec,-1);
    // multiply with purity for primary particles
    histoGammaCorr->Multiply(histoGammaCorr,histoPurity,1.,1.,"");
    // divide by reconstruction efficiency
    histoGammaCorr->Divide(histoGammaCorr,histoRecoEff,1.,1.,"");
    // scale to with 1/deltaEta and 1/deltaPhi
    histoGammaCorr->Scale(1./deltaEta);
    histoGammaCorr->Scale(scaling);
    // convert to inv yield by dividing by pT at bin center
    for (Int_t i = 1; i < histoGammaCorr->GetNbinsX()+1 ; i++){
        Double_t newBinContent  = histoGammaCorr->GetBinContent(i)  / histoGammaCorr->GetBinCenter(i);
        Double_t newBinError    = histoGammaCorr->GetBinError(i)    / histoGammaCorr->GetBinCenter(i);
        histoGammaCorr->SetBinContent(i,newBinContent);
        histoGammaCorr->SetBinError(i,  newBinError);
    }
}

//*********************************************************************************************************
//*** Calorimeter method                                                                         **********
//*** Correct inclusive raw gamma candite spectrum with bin-by-bin efficiencies                  **********
// - scale with 1/nEvt first
// - subtract secodary gammas (cocktail approach)
// - take out impurities
// - reco efficiency (histoRecoEff)
// - conversion probability (histoConvProb)
// - scale with 1/deltaEta, scaling, 1/nEvt
// - convert to inv yield (scale by bin 1/bin center)
//*********************************************************************************************************
void CorrectGammaEffiResolCocktail(TH1D*       histoGammaCorr,
                                   TH1D**      histoHadronSec,
                                   TH1D*       histoRestSec,
                                   TH1D*       histoPurity,
                                   TH1D*       histoRecoEff,
                                   Double_t    deltaEta,
                                   Double_t    scaling,
                                   Double_t    nEvt
                                   ){
    // scale with 1/nEvt as secondaries are normalized per event
    histoGammaCorr->Scale(1./nEvt);
    // subtract secondary gammas from decays
    for (Int_t k = 0; k < 3; k++){
        if (histoHadronSec[k]) histoGammaCorr->Add(histoHadronSec[k],-1);
    }
    // subtract secondary gammas from material interactions
    histoGammaCorr->Add(histoRestSec,-1);
    // multiply with purity for primary particles
    histoGammaCorr->Multiply(histoGammaCorr,histoPurity,1.,1.,"");
    // divide by reconstruction efficiency
    histoGammaCorr->Divide(histoGammaCorr,histoRecoEff,1.,1.,"");
    // scale to with 1/deltaEta and 1/deltaPhi
    histoGammaCorr->Scale(1./deltaEta);
    histoGammaCorr->Scale(scaling);
    // convert to inv yield by dividing by pT at bin center
    for (Int_t i = 1; i < histoGammaCorr->GetNbinsX()+1 ; i++){
        Double_t newBinContent  = histoGammaCorr->GetBinContent(i)  / histoGammaCorr->GetBinCenter(i);
        Double_t newBinError    = histoGammaCorr->GetBinError(i)    / histoGammaCorr->GetBinCenter(i);
        histoGammaCorr->SetBinContent(i,newBinContent);
        histoGammaCorr->SetBinError(i,  newBinError);
    }
}

//*********************************************************************************************************
//*** Correct inclusive raw gamma candite spectrum for secondaries and purity                    **********
// - subtract secodary gammas (MC approach)
// - take out impurities
//*********************************************************************************************************
void CorrectGammaSecAndPurity(TH1D* histoGammaCorr,
                              TH1D* histoSecGamma,
                              TH1D* histoSecGammaAddK0s,
                              TH1D* histoPurity
                              ){
    histoGammaCorr->Sumw2();
    histoGammaCorr->Add(histoSecGamma,-1);
    histoGammaCorr->Add(histoSecGammaAddK0s,-1);
    histoGammaCorr->Multiply(histoGammaCorr,histoPurity,1.,1.,"");
}

//*********************************************************************************************************
//*** Correct inclusive raw gamma candite spectrum for secondaries and purity                    **********
// - subtract secodary gammas (cocktail approach)
// - take out impurities
//*********************************************************************************************************
void CorrectGammaSecAndPurityCocktail(TH1D*   histoGammaCorr,
                                      TH1D**  histoSecGammaHadrons,
                                      TH1D*   histoSecGammaRest,
                                      TH1D*   histoPurity
                                      ){
    histoGammaCorr->Sumw2();
    for (Int_t k = 0; k < 3; k++){
        if (histoSecGammaHadrons[k]) histoGammaCorr->Add(histoSecGammaHadrons[k],-1);
    }
    histoGammaCorr->Add(histoSecGammaRest,-1);
    histoGammaCorr->Multiply(histoGammaCorr,histoPurity,1.,1.,"");
}

//*********************************************************************************************************
//*********** Correct inclusive raw gamma spectrum from conversions using unfolded as inputs **************
// - conversion probability (histoConvProb)
// - reco efficiency (histoRecoEff)
// - scale with 1/deltaEta, scaling, 1/nEvt
// - convert to inv yield (scale by bin 1/bin center)
//*********************************************************************************************************
void CorrectGammaUnfoldResol(TH1D*       histoGammaCorr,
                             TH1D*       histoConvProb,
                             TH1D*       histoRecoEff,
                             Double_t    deltaEta,
                             Double_t    scaling,
                             Double_t    nEvt
                             ){
    // scale with 1/reco eff
    histoGammaCorr->Divide(histoGammaCorr,histoRecoEff,1.,1.,"");
    // scale with 1/ P_{conv}
    histoGammaCorr->Divide(histoGammaCorr,histoConvProb,1.,1.,"");
    // normalize to pseudorapitdyrapidity and azimuthal range
    histoGammaCorr->Scale(1./deltaEta);
    histoGammaCorr->Scale(scaling);
    // normalize per event
    histoGammaCorr->Scale(1./nEvt);
    // convert to inv yield (divide each yield by 1/pT)
    for (Int_t i = 1; i < histoGammaCorr->GetNbinsX()+1 ; i++){
        Double_t newBinContent  = histoGammaCorr->GetBinContent(i)  / histoGammaCorr->GetBinCenter(i);
        Double_t newBinError    = histoGammaCorr->GetBinError(i)    / histoGammaCorr->GetBinCenter(i);
        histoGammaCorr->SetBinContent(i,newBinContent);
        histoGammaCorr->SetBinError(i,  newBinError);
    }
}

//*********************************************************************************************************
//*********** Correct inclusive raw gamma spectrum from calorimeter using unfolded as inputs **************
// - reco efficiency (histoRecoEff)
// - scale with 1/deltaEta, scaling, 1/nEvt
// - convert to inv yield (scale by bin 1/bin center)
//*********************************************************************************************************
void CorrectGammaUnfoldResol(TH1D*       histoGammaCorr,
                             TH1D*       histoRecoEff,
                             Double_t    deltaEta,
                             Double_t    scaling,
                             Double_t    nEvt
                             ){
    // scale with 1/reco eff
    histoGammaCorr->Divide(histoGammaCorr,histoRecoEff,1.,1.,"");
    // normalize to pseudorapitdyrapidity and azimuthal range
    histoGammaCorr->Scale(1./deltaEta);
    histoGammaCorr->Scale(scaling);
    // normalize per event
    histoGammaCorr->Scale(1./nEvt);
    // convert to inv yield (divide each yield by 1/pT)
    for (Int_t i = 1; i < histoGammaCorr->GetNbinsX()+1 ; i++){
        Double_t newBinContent  = histoGammaCorr->GetBinContent(i)  / histoGammaCorr->GetBinCenter(i);
        Double_t newBinError    = histoGammaCorr->GetBinError(i)    / histoGammaCorr->GetBinCenter(i);
        histoGammaCorr->SetBinContent(i,newBinContent);
        histoGammaCorr->SetBinError(i,  newBinError);
    }
}

//*********************************************************************************************************
//*********** Convert secondary spectrum from cocktail to raw (conversion method)            **************
// - multiply with secondary conv. prob.
// - use unfolding for resolution correction if response matrix given
// - multiply with secondary reco. eff.
// - scale with nEvts for same event norm as data spectrum
//*********************************************************************************************************
Bool_t ConvertCocktailSecondaryToRaw(TH1D*       histoGammaSec,
                                     TH1D*       histoConvProb,
                                     TH1D*       histoRecoEff,
                                     TH2D*       responseMatrix,
                                     Double_t    nEvt,
                                     Bool_t      useResponseMatrix,
                                     Int_t       nIterationsUnfolding = 5
                                     ){
    // multiply with conv. prob.
    if (!histoConvProb) return kFALSE;
    histoGammaSec->Multiply(histoConvProb);
    // correct for resolution, if response matrix given (otherwise included in reco. eff.)
    if (useResponseMatrix) {
        if (!responseMatrix) return kFALSE;
        cout << "will use unfolding for " << histoGammaSec->GetName() << endl;
        RooUnfoldResponse response(0,0,responseMatrix);
        RooUnfoldBayes unfold_SpectrumCocktail (&response, histoGammaSec, nIterationsUnfolding);
        histoGammaSec = (TH1D*)unfold_SpectrumCocktail.Hreco();
    }
    // multiply with reco. eff.
    if (!histoRecoEff) return kFALSE;
    histoGammaSec->Multiply(histoRecoEff);
    // scale with nEvt from data
    histoGammaSec->Scale(nEvt);
    
    return kTRUE;
}

//*********************************************************************************************************
//*********** Convert secondary spectrum from cocktail to raw (calorimeter method)           **************
// - use unfolding for resolution correction if response matrix given
// - multiply with secondary reco. eff.
// - scale with nEvts for same event norm as data spectrum
//*********************************************************************************************************
Bool_t ConvertCocktailSecondaryToRaw(TH1D*       histoGammaSec,
                                     TH1D*       histoRecoEff,
                                     TH2D*       responseMatrix,
                                     Double_t    nEvt,
                                     Bool_t      useResponseMatrix,
                                     Int_t       nIterationsUnfolding = 5
                                     ){
    // correct for resolution, if response matrix given (otherwise included in reco. eff.)
    if (useResponseMatrix) {
        if (!responseMatrix) return kFALSE;
        cout << "will use unfolding for " << histoGammaSec->GetName() << endl;
        RooUnfoldResponse response(0,0,responseMatrix);
        RooUnfoldBayes unfold_SpectrumCocktail (&response, histoGammaSec, nIterationsUnfolding);
        histoGammaSec = (TH1D*)unfold_SpectrumCocktail.Hreco();
    }
    // multiply with reco. eff.
    if (!histoRecoEff) return kFALSE;
    histoGammaSec->Multiply(histoRecoEff);
    // scale with nEvt from data
    histoGammaSec->Scale(nEvt);
    
    return kTRUE;
}

//*********************************************************************************************************
//*********** Convert MC inclusive gamma spectrum to same norm as data                       **************
// - normalize to pseudorapidity and azimuthal range
// - normalize to nEvts in MC
// - convert to inv yield (scale by bin 1/bin center)
//*********************************************************************************************************
void CorrectGammaMC(TH1D*       histoMCGammaSpec_MCPtCorr,
                    Double_t    deltaEta,
                    Double_t    scaling,
                    Double_t    nEvtMC
                    ){
    // normalize to pseudorapidity and azimuthal range
    histoMCGammaSpec_MCPtCorr->Scale(1./deltaEta);
    histoMCGammaSpec_MCPtCorr->Scale(scaling);
    // normalize to nEvts in MC
    histoMCGammaSpec_MCPtCorr->Scale(1./nEvtMC);
    // convert to inv yield (scale by bin 1/bin center)
    for (Int_t i = 1; i < histoMCGammaSpec_MCPtCorr->GetNbinsX()+1 ; i++){
        Double_t newBinContent = histoMCGammaSpec_MCPtCorr->GetBinContent(i)/histoMCGammaSpec_MCPtCorr->GetBinCenter(i);
        Double_t newBinError = histoMCGammaSpec_MCPtCorr->GetBinError(i)/histoMCGammaSpec_MCPtCorr->GetBinCenter(i);
        histoMCGammaSpec_MCPtCorr->SetBinContent(i,newBinContent);
        histoMCGammaSpec_MCPtCorr->SetBinError(i,newBinError);
    }
}

