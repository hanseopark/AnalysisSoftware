//
// This code routines addresses the following tasks
// - Calculation of v2,dir from v2,inv, v2,dec, and Rgamma
// - Averaging PCM and PHOS results
// - Quantifying the level of agreement with certain hypotheses about v2,dir (to be done)
//
//
// Klaus Reygers, November 2016
//

#include "TMath.h"
#include "TVectorD.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraphAsymmErrors.h"
#include "TLatex.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TROOT.h"

#include "RooArgList.h"
#include "RooRealVar.h"
#include "RooMultiVarGaussian.h"
#include "RooDataSet.h"


void weighted_average(const Int_t nMeas, const TVectorD** mu, const TMatrixDSym** cov, TVectorD* muAve, TMatrixDSym* covAve) {

  // weighted average of nMeas independent sets of data points (each with the same number of data points)
  //
  // nMeas : number of data sets to average
  // mu[i]: pointer to TVectorD of measurement i
  // cov[i]: pointer to TMatrixDSym of measurement i
  //
  // result:
  // muAve: pointer to mean (TVectorD) of the weighted average
  // covAve: pointer to covariance matrix (TMatrixDSym) of the weighted average
  //
  // Klaus Reygers, September 2016
  // reygers@physi.uni-heidelberg.de

  // Get dimension of the vectors and matrizes
  const Int_t nDim = mu[0]->GetNrows();

  // create empty vector and empty nDim x nDim matrix (all elements are 0)
  TVectorD sumMu(nDim);
  TMatrixDSym sumCov(nDim);

  // loop over the different measurements
  for (Int_t iMeas=0; iMeas<nMeas; iMeas++) {
    TMatrixDSym Vinv = *cov[iMeas];
    Vinv.Invert();
    sumMu += Vinv * (*mu[iMeas]);
    sumCov += Vinv;
  };

  // set covariance of the weighted average
  *covAve = sumCov.Invert();

  // set weighted average
  *muAve = (*covAve) * sumMu;

}


void weighted_average(const TVectorD& mu1, const TVectorD& mu2, const TMatrixDSym& cov1, const TMatrixDSym& cov2, TVectorD& muAve, TMatrixDSym& covAve) {

  // weighted average of two independent sets of measurements
  // mu1: mean values data set 1
  // mu2: mean values data set 2
  // cov1: covariance matrix data set 1
  // cov2: covariance matrix data set 2
  //
  // result:
  // muAve: weighted average (TVectorD)
  // covAve: covariance Matrix of the weighted average

  const TVectorD* mu[2];
  mu[0] = &mu1; mu[1] = &mu2;

  const TMatrixDSym* cov[2];
  cov[0] = &cov1; cov[1] = &cov2;

  weighted_average(2, mu, cov, &muAve, &covAve);

}


void likelihood(const Int_t nPtBins, const Int_t nSamples,
		const TVectorD& vnDirTrueVec,
		const TVectorD& RgamMeasVec, const TMatrixDSym& covRgam,
		const TVectorD& vnIncMeasVec, const TMatrixDSym& covVnInc,
		const TVectorD& vnDecMeasVec, const TMatrixDSym& covVnDec, Double_t& Lfinal) {

  // determine the likelihood of a hypothesis about v2,dir
  //
  // Approach:
  // Repeat these steps many times:
  // - get RgamTrueVec from MC sampling
  // - get vnDecTrueVec from MC sampling
  // - calculate vnIncTrueVec from vnDirTrueVec, vnDecTrueVec, and RgamTrueVec
  // - calculate likelihood L(measured Rgam, vnInc, VnDec | true Rgam, vnInc, VnDec)
  //
  // We then take the average of the so determined likelihoods

  //
  // RooFit variables for vnInc
  //
  RooArgList vnInc;
  RooArgList vnIncMeas;
  RooRealVar* vnIncVar[nPtBins];
  RooRealVar* vnIncMeasVar;

  // RooFit variables for inclusive photon vn
  for (Int_t i=0; i<nPtBins; i++) {
    char* nameMeas = Form("vnIncMeas%d", i);
    vnIncMeasVar = new RooRealVar(nameMeas, nameMeas, vnIncMeasVec[i], -0.5, 0.5);
    vnIncMeas.add(*vnIncMeasVar);

    char* nameTrue = Form("vnInc%d", i);
    vnIncVar[i] = new RooRealVar(nameTrue, nameTrue, 0., -0.5, 0.5);
    vnInc.add(*vnIncVar[i]);
  }

  // now make the multivariate Gaussian
  RooMultiVarGaussian mvgVnInc("mvgVnInc", "mvgVnInc", vnInc, vnIncMeas, covVnInc);

  //
  // Rgamma
  //
  RooArgList Rgam;
  RooArgList RgamMeas;
  RooRealVar* RgamVar;
  RooRealVar* RgamMeasVar;

  // RooFit variables Rgamma
  for (Int_t i=0; i<nPtBins; i++) {
    char* nameMeas = Form("RgamMeas%d", i);
    RgamMeasVar = new RooRealVar(nameMeas, nameMeas, RgamMeasVec[i], 0., 10.);
    RgamMeas.add(*RgamMeasVar);

    char* nameTrue = Form("Rgam%d", i);
    RgamVar = new RooRealVar(nameTrue, nameTrue, 1., 1., 10.); // here the prior knowledge Rgamma >= 1 enters !
    Rgam.add(*RgamVar);
  }

  // now make the multivariate Gaussian
  RooMultiVarGaussian mvgRgam("mvgRgam", "mvgRgam", Rgam, RgamMeas, covRgam);

  //
  // decay photon vn
  //
  RooArgList vnDec;
  RooArgList vnDecMeas;
  RooRealVar* vnDecVar;
  RooRealVar* vnDecMeasVar;

  // RooFit variables for decay photon vn
  for (Int_t i=0; i<nPtBins; i++) {
    char* nameMeas = Form("vnDecMeas%d", i);
    vnDecMeasVar = new RooRealVar(nameMeas, nameMeas, vnDecMeasVec[i], -0.5, 0.5);
    vnDecMeas.add(*vnDecMeasVar);

    char* nameTrue = Form("vnDec%d", i);
    vnDecVar = new RooRealVar(nameTrue, nameTrue, 0., -0.5, 0.5);
    vnDec.add(*vnDecVar);
  }

  // now make the multivariate Gaussian
  RooMultiVarGaussian mvgVnDec("mvgVnDec", "mvgVnDec", vnDec, vnDecMeas, covVnDec);

  //
  // sample the distributions
  //
  RooDataSet* dataVnDec = mvgVnDec.generate(vnDec, nSamples);
  RooDataSet* dataRgam = mvgRgam.generate(Rgam, nSamples);

  //
  // create Likelihood variable
  //
  RooArgList L;
  RooRealVar* LVar = new RooRealVar("L", "L", 0., 1.);
  L.add(*LVar);

  RooDataSet dataL("dataL","dataL", L);

  // mvgRgam.Print();
  for (Int_t k=0; k<dataVnDec->numEntries(); k++) {

      const RooArgSet* rowVnDec = dataVnDec->get(k);
      const RooArgSet* rowRgam = dataRgam->get(k);

      vnDec = *rowVnDec;
      Rgam = *rowRgam;

      // cout << mvgRgam.getVal() << endl;
      for (Int_t i=0; i<nPtBins; i++) {

	RooRealVar* vnDecRow = (RooRealVar*) rowVnDec->find(vnDec.at(i)->GetName());
	RooRealVar* RgamRow = (RooRealVar*) rowRgam->find(Rgam.at(i)->GetName());

	Double_t RgamVal =  RgamRow->getVal();
	Double_t vnDecVal = vnDecRow->getVal();
      	Double_t vnIncVal = (vnDirTrueVec(i)*(RgamVal - 1.) + vnDecVal)/RgamVal;

	// cout << "vnIncVal = " << vnIncVal << ", RgamVal = " << RgamVal << ", vnDecVal = " << vnDecVal << endl;

	// ((RooRealVar*) vnInc.at(i))->setVal(vnIncVal);
	*vnIncVar[i] =  vnIncVal;

    }

    Double_t LRgam = mvgRgam.getVal();
    Double_t LvnDec = mvgVnDec.getVal();
    Double_t LvnInc = mvgVnInc.getVal();

    // cout << "LRgam = " << LRgam << ", LvnDec = " << LvnDec << ", LvnInc = " << LvnInc << endl;

    *LVar = LRgam * LvnDec * LvnInc;

    dataL.add(*LVar);

  }

  // At the meoment this number is proportional to the likelihood,
  // because the it seems RooFit doesn't return the normalized values.
  // This doesn't matter however, as we are only interested in likelihood ratios
  Lfinal = dataL.mean(*LVar);
  cout << "likelihood: " << Lfinal << endl;

}


void v2dir(const Int_t nPtBins, const Int_t nSamples,
	   const TVectorD& RgamMeasVec, const TMatrixDSym& covRgam,
	   const TVectorD& vnIncMeasVec, const TMatrixDSym& covVnInc,
	   const TVectorD& vnDecMeasVec, const TMatrixDSym& covVnDec,
	   TH1** hVnDir) {

  //
  // fill histograms with v2,dir true for all pT bins by sampling multi-variate Gaussians
  //

  //
  // Rgamma
  //
  RooArgList Rgam;
  RooArgList RgamMeas;
  RooRealVar* RgamVar;
  RooRealVar* RgamMeasVar;

  // RooFit variables Rgamma
  for (Int_t i=0; i<nPtBins; i++) {
    char* nameMeas = Form("RgamMeas%d", i);
    RgamMeasVar = new RooRealVar(nameMeas, nameMeas, RgamMeasVec[i], 0., 10.);
    RgamMeas.add(*RgamMeasVar);

    char* nameTrue = Form("Rgam%d", i);
    RgamVar = new RooRealVar(nameTrue, nameTrue, 1., 1., 10.);
    Rgam.add(*RgamVar);
  }

  // now make the multivariate Gaussian
  RooMultiVarGaussian mvgRgam("mvgRgam", "mvgRgam", Rgam, RgamMeas, covRgam);

  //
  // inclusive photon vn
  //
  RooArgList vnInc;
  RooArgList vnIncMeas;
  RooRealVar* vnIncVar;
  RooRealVar* vnIncMeasVar;

  // RooFit variables for inclusive photon vn
  for (Int_t i=0; i<nPtBins; i++) {
    char* nameMeas = Form("vnIncMeas%d", i);
    vnIncMeasVar = new RooRealVar(nameMeas, nameMeas, vnIncMeasVec[i], -0.5, 0.5);
    vnIncMeas.add(*vnIncMeasVar);

    char* nameTrue = Form("vnInc%d", i);
    vnIncVar = new RooRealVar(nameTrue, nameTrue, 0., -0.5, 0.5);
    vnInc.add(*vnIncVar);
  }

  // now make the multivariate Gaussian
  RooMultiVarGaussian mvgVnInc("mvgVnInc", "mvgVnInc", vnInc, vnIncMeas, covVnInc);

  //
  // decay photon vn
  //
  RooArgList vnDec;
  RooArgList vnDecMeas;
  RooRealVar* vnDecVar;
  RooRealVar* vnDecMeasVar;

  // RooFit variables for decay photon vn
  for (Int_t i=0; i<nPtBins; i++) {
    char* nameMeas = Form("vnDecMeas%d", i);
    vnDecMeasVar = new RooRealVar(nameMeas, nameMeas, vnDecMeasVec[i], -0.5, 0.5);
    vnDecMeas.add(*vnDecMeasVar);

    char* nameTrue = Form("vnDec%d", i);
    vnDecVar = new RooRealVar(nameTrue, nameTrue, 0., -0.5, 0.5);
    vnDec.add(*vnDecVar);
  }

  // now make the multivariate Gaussian
  RooMultiVarGaussian mvgVnDec("mvgVnDec", "mvgVnDec", vnDec, vnDecMeas, covVnDec);

  //
  // sample the distributions
  //
  RooDataSet* dataVnInc = mvgVnInc.generate(vnInc, nSamples);
  RooDataSet* dataVnDec = mvgVnDec.generate(vnDec, nSamples);
  RooDataSet* dataRgam = mvgRgam.generate(Rgam, nSamples);

  //
  // create direct photon vn data set based on the sampled values of VnInc, VnDec, Rgam
  //
  RooArgList vnDir;
  RooRealVar* vnDirVar;

  // RooFit variables for direct photon vn
  for (Int_t i=0; i<nPtBins; i++) {
    char* name = Form("vnDir%d", i);
    vnDirVar = new RooRealVar(name, name, 0., -0.5, 0.5);
    vnDir.add(*vnDirVar);
  }

  RooDataSet dataVnDir("dataVnDir","dataVnDir", vnDir);

  for (Int_t k=0; k<dataVnInc->numEntries(); k++) {
    for (Int_t i=0; i<nPtBins; i++) {

      const RooArgSet* rowVnInc = dataVnInc->get(k);
      const RooArgSet* rowVnDec = dataVnDec->get(k);
      const RooArgSet* rowRgam = dataRgam->get(k);

      RooRealVar* vnIncRow = (RooRealVar*) rowVnInc->find(vnInc.at(i)->GetName());
      RooRealVar* vnDecRow = (RooRealVar*) rowVnDec->find(vnDec.at(i)->GetName());
      RooRealVar* RgamRow = (RooRealVar*) rowRgam->find(Rgam.at(i)->GetName());

      Double_t vnIncVal = vnIncRow->getVal();
      Double_t vnDecVal = vnDecRow->getVal();
      Double_t RgamVal = RgamRow->getVal();

      // cout << vnIncVal << " " << vnDecVal << " " << RgamVal << endl;

      Double_t vnDirVal = (RgamVal*vnIncVal - vnDecVal)/(RgamVal - 1.);

      ((RooRealVar*) vnDir.at(i))->setVal(vnDirVal);

      dataVnDir.add(vnDir);

    }

  }

  //
  // some output
  //
  // dataVnDir.covarianceMatrix()->Print();

  // RooPlot* vnDirframe = ((RooRealVar*) vnDir.at(9))->frame(1000);
  // dataVnDir.plotOn(vnDirframe);
  // vnDirframe->Draw("h");

  // TH1* hTest = dataVnDir.createHistogram("test", (*(RooRealVar*) vnDir.at(9)));
  //hTest->Draw();
  for (Int_t i=0; i<nPtBins; i++) {
    TString s = Form("PtBin_%d", i);
    hVnDir[i] = dataVnDir.createHistogram(s.Data(), (*(RooRealVar*) vnDir.at(i)));
  }


};


// ------------------------------------------------------------------
// The entry point of this macro
// ------------------------------------------------------------------

void v2dir_pcm_phos_comb() {

  //
  // calculate v2,dir(pT) with uncertainties from R_gamma(pT), v2,inc(pT),  v2,dec(pT) (with toy data)
  //
  // input:
  // - measured  R_gamma(pT) + covariance matrix
  // - measured  v2,inc(pT) + covariance matrix
  // - measured  v2,dec(pT) + covariance matrix
  //
  //
  // Klaus Reygers, November 2016
  //

  // number of pT bins
  const Int_t nPtBins = 16;

  // define the centrality class
  // TString centr = "00-20";
  TString centr = "20-40";
  // TString centr = "40-80";

  // define how many v2dir values are drawn randomely
  // const Int_t nSamples = 100000;
  const Int_t nSamples = 1000;

  //
  // PCM specific input
  //
  TFile f_pcm("data_pcm/data_pcm.root");
  TDirectory* dir_pcm = (TDirectory*) f_pcm.Get(centr.Data());

  TVectorD v2_inc_meas_values_pcm = *(TVectorD*) dir_pcm->Get("v2_inc_meas_values");
  TMatrixDSym cov_v2_inc_toterr_pcm = *(TMatrixDSym*) dir_pcm->Get("cov_v2_inc_toterr");
  TMatrixDSym cov_v2_inc_staterr_pcm = *(TMatrixDSym*) dir_pcm->Get("cov_v2_inc_staterr");

  // only needed for cross checks
  TVectorD Rgamma_meas_vec_pcm = *(TVectorD*) dir_pcm->Get("Rgamma_meas_vec");
  TMatrixDSym cov_Rgamma_syserr_pcm = *(TMatrixDSym*) dir_pcm->Get("cov_Rgamma_syserr");

  //
  // PHOS specific input
  //
  TFile f_phos("data_phos/data_phos.root");
  TDirectory* dir_phos = (TDirectory*) f_phos.Get(centr.Data());

  TVectorD v2_inc_meas_values_phos = *(TVectorD*) dir_phos->Get("v2_inc_meas_values");
  TMatrixDSym cov_v2_inc_toterr_phos = *(TMatrixDSym*) dir_phos->Get("cov_v2_inc_toterr");
  TMatrixDSym cov_v2_inc_staterr_phos = *(TMatrixDSym*) dir_phos->Get("cov_v2_inc_staterr");

  // only needed for cross checks
  TVectorD Rgamma_meas_vec_phos = *(TVectorD*) dir_phos->Get("Rgamma_meas_vec");
  TMatrixDSym cov_Rgamma_phos = *(TMatrixDSym*) dir_phos->Get("cov_Rgamma");

  //
  // PCM/PHOS combined
  //
  TVectorD Rgamma_meas_vec_comb = *(TVectorD*) dir_pcm->Get("Rgamma_meas_vec_comb");
  TMatrixDSym cov_Rgamma_comb_toterr = *(TMatrixDSym*) dir_pcm->Get("cov_Rgamma_comb_toterr");
  TMatrixDSym cov_Rgamma_comb_staterr = *(TMatrixDSym*) dir_pcm->Get("cov_Rgamma_comb_staterr");

  // decay photon v2 from PCM file
  // TVectorD v2_dec_meas_values = *(TVectorD*) dir_pcm->Get("v2_dec_meas_values");
  // TMatrixDSym cov_v2_dec = *(TMatrixDSym*) dir_pcm->Get("cov_v2_dec");

  // decay photon v2 from PHOS file
  TVectorD v2_dec_meas_values = *(TVectorD*) dir_phos->Get("v2_dec_meas_values_syserr");
  TMatrixDSym cov_v2_dec_syserr = *(TMatrixDSym*) dir_phos->Get("cov_v2_dec_syserr");
  TMatrixDSym cov_v2_dec_staterr = *(TMatrixDSym*) dir_phos->Get("cov_v2_dec_staterr");

  TVectorD pt = *(TVectorD*) dir_pcm->Get("pt");

  // create array of pointers to output histograms (vndir distributions for each pT bin)
  TH1* h_vn_dir_pcm_staterr[nPtBins];
  TH1* h_vn_dir_phos_staterr[nPtBins];
  TH1* h_vn_dir_comb_staterr[nPtBins];
  TH1* h_vn_dir_pcm_toterr[nPtBins];
  TH1* h_vn_dir_phos_toterr[nPtBins];
  TH1* h_vn_dir_comb_toterr[nPtBins];

  //
  // calculate vnDir PCM (total errors, using combined Rgamma and common v2dec)
  //

  // stat. err.
  v2dir(nPtBins, nSamples, Rgamma_meas_vec_comb, cov_Rgamma_comb_staterr, v2_inc_meas_values_pcm, cov_v2_inc_toterr_pcm,
	v2_dec_meas_values, cov_v2_dec_staterr, h_vn_dir_pcm_staterr);

  // tot. err.
  v2dir(nPtBins, nSamples, Rgamma_meas_vec_comb, cov_Rgamma_comb_toterr, v2_inc_meas_values_pcm, cov_v2_inc_toterr_pcm,
	v2_dec_meas_values, cov_v2_dec_syserr, h_vn_dir_pcm_toterr);

  //
  // calculate vnDir PHOS (using combined Rgamma and common v2dec)
  //

  // stat. err.
  v2dir(nPtBins, nSamples, Rgamma_meas_vec_comb, cov_Rgamma_comb_staterr, v2_inc_meas_values_phos, cov_v2_inc_staterr_phos,
	v2_dec_meas_values, cov_v2_dec_staterr, h_vn_dir_phos_staterr);

  // tot. err.
  v2dir(nPtBins, nSamples, Rgamma_meas_vec_comb, cov_Rgamma_comb_toterr, v2_inc_meas_values_phos, cov_v2_inc_toterr_phos,
	v2_dec_meas_values, cov_v2_dec_syserr, h_vn_dir_phos_toterr);


  //
  // calculate weighted average of v2inc
  //

  // weighted average of Rgamma
  // TVectorD Rgamma_meas_vec_comb(nPtBins);
  // TMatrixDSym cov_Rgamma_comb(nPtBins);
  // weighted_average(Rgamma_meas_vec_syserr_pcm, Rgamma_meas_vec_phos, cov_Rgamma_syserr_pcm, cov_Rgamma_phos, Rgamma_meas_vec_comb, cov_Rgamma_comb);

  //
  // weighted average of v2,inc
  //

  // using total error as weight
  TVectorD v2_inc_meas_values_comb(nPtBins);
  TMatrixDSym cov_v2_inc_toterr_comb(nPtBins);
  weighted_average(v2_inc_meas_values_pcm, v2_inc_meas_values_phos, cov_v2_inc_toterr_pcm, cov_v2_inc_toterr_phos,
		   v2_inc_meas_values_comb, cov_v2_inc_toterr_comb);

  // using statistical error as weight
  TVectorD v2_inc_meas_values_comb_statweights(nPtBins);
  TMatrixDSym cov_v2_inc_staterr_comb(nPtBins);
  weighted_average(v2_inc_meas_values_pcm, v2_inc_meas_values_phos, cov_v2_inc_staterr_pcm, cov_v2_inc_staterr_phos,
		   v2_inc_meas_values_comb_statweights, cov_v2_inc_staterr_comb);

  //
  // calculate vnDir combined
  //

  // stat. err.
  v2dir(nPtBins, nSamples, Rgamma_meas_vec_comb, cov_Rgamma_comb_staterr, v2_inc_meas_values_comb, cov_v2_inc_staterr_comb,
	v2_dec_meas_values, cov_v2_dec_staterr, h_vn_dir_comb_staterr);

  // tot. err.
  v2dir(nPtBins, nSamples, Rgamma_meas_vec_comb, cov_Rgamma_comb_toterr, v2_inc_meas_values_comb, cov_v2_inc_toterr_comb,
	v2_dec_meas_values, cov_v2_dec_syserr, h_vn_dir_comb_toterr);

  // as a check, calculate weighted average of Rgamma here and determine v2,dir
  // v2dir(nPtBins, nSamples, Rgamma_meas_vec_comb, cov_Rgamma_comb, v2_inc_meas_values_comb, cov_v2_inc_comb,
  //	v2_dec_meas_values, cov_v2_dec, h_vn_dir_comb);


  // graphs to store v2dir
  TGraphAsymmErrors g_vn_dir_pcm_toterr(nPtBins);
  TGraphAsymmErrors g_vn_dir_phos_toterr(nPtBins);
  TGraphAsymmErrors g_vn_dir_comb_toterr(nPtBins);
  TGraphAsymmErrors g_vn_dir_comb_staterr(nPtBins); // central values same as g_vn_dir_comb_toterr

  // graphs to store v2inc
  TGraphAsymmErrors g_vn_inc_pcm_toterr(nPtBins);
  TGraphAsymmErrors g_vn_inc_phos_toterr(nPtBins);
  TGraphAsymmErrors g_vn_inc_comb_toterr(nPtBins);

  // graph to store decay photon v2
  TGraphAsymmErrors g_vn_dec_comb(nPtBins);

  //
  // plot output
  //
  cout << "#entries = " << h_vn_dir_pcm_toterr[0]->GetEntries() << endl;

  // style settings
  TStyle* myStyle = new TStyle("myStyle","My root style");
  myStyle->SetOptStat(kFALSE);
  myStyle->SetLabelOffset(0.005,"x");   // 0.005 = root default
  myStyle->SetLabelOffset(0.005,"y");   // 0.005 = root default
  myStyle->SetTitleXSize(0.04);         // 0.04  = root default
  myStyle->SetTitleYSize(0.04);         // 0.04  = root default
  myStyle->SetTitleXOffset(1.2);
  myStyle->SetTitleYOffset(1.3);
  myStyle->SetPadLeftMargin(0.15);
  myStyle->SetPadRightMargin(0.12);      // 0.1 = root default
  myStyle->SetPadTopMargin(0.1);
  myStyle->SetPadBottomMargin(0.12);
  myStyle->SetCanvasColor(0);
  myStyle->SetPadColor(0);
  myStyle->SetCanvasBorderMode(0);
  myStyle->SetPadBorderMode(0);
  myStyle->SetOptTitle(0);
  gROOT->SetStyle("myStyle");

  //
  // Plot v2dir distributions for each pT bin (PCM, PHOS, combined)
  //
  TCanvas* c1 = new TCanvas("c1","c1", 700, 900);
  c1->Divide(4,4);

  // gStyle->SetOptStat(kFALSE);
  // gStyle->SetTitleFontSize(0.08);
  // gROOT->ForceStyle();

  TLatex lpt;

  for (Int_t i=0; i<nPtBins; i++) {
    c1->cd(i+1);
    c1->GetPad(i+1)->SetBottomMargin(0.16);
    c1->GetPad(i+1)->SetLeftMargin(0.16);

    h_vn_dir_comb_toterr[i]->SetXTitle("v2,dir");
    h_vn_dir_comb_toterr[i]->SetTitle(Form("pT bin %d",i));

    h_vn_dir_comb_toterr[i]->SetLineColor(kBlack); h_vn_dir_comb_toterr[i]->DrawCopy("HL");
    h_vn_dir_pcm_toterr[i]->SetLineColor(kRed); h_vn_dir_pcm_toterr[i]->DrawCopy("HL,same");
    h_vn_dir_phos_toterr[i]->SetLineColor(kBlue); h_vn_dir_phos_toterr[i]->DrawCopy("HL,same");

    lpt.DrawLatexNDC(0.2, 0.75, Form("p_{T} = %4.1f GeV/#it{c}",pt(i)));

    // determine median and errors
    const Int_t n_quant = 3;
    Double_t prob_integrals[n_quant] = {0.15865, 0.5, 0.84135};
    Double_t quantiles[n_quant];

    h_vn_dir_pcm_toterr[i]->GetQuantiles(n_quant, quantiles, prob_integrals);
    Double_t central_value = quantiles[1];
    Double_t err_high = quantiles[2] - central_value;
    Double_t err_low = central_value - quantiles[0];

    g_vn_dir_pcm_toterr.SetPoint(i, pt(i)+0.03, central_value);
    g_vn_dir_pcm_toterr.SetPointEYhigh(i, err_high);
    g_vn_dir_pcm_toterr.SetPointEYlow(i, err_low);

    h_vn_dir_phos_toterr[i]->GetQuantiles(n_quant, quantiles, prob_integrals);
    central_value = quantiles[1];
    err_high = quantiles[2] - central_value;
    err_low = central_value - quantiles[0];

    g_vn_dir_phos_toterr.SetPoint(i, pt(i)+0.06, central_value);
    g_vn_dir_phos_toterr.SetPointEYhigh(i, err_high);
    g_vn_dir_phos_toterr.SetPointEYlow(i, err_low);

    h_vn_dir_comb_toterr[i]->GetQuantiles(n_quant, quantiles, prob_integrals);
    central_value = quantiles[1];
    err_high = quantiles[2] - central_value;
    err_low = central_value - quantiles[0];

    g_vn_dir_comb_toterr.SetPoint(i, pt(i), central_value);
    g_vn_dir_comb_toterr.SetPointEYhigh(i, err_high);
    g_vn_dir_comb_toterr.SetPointEYlow(i, err_low);
    g_vn_dir_comb_toterr.SetPointEXhigh(i, 0.08); // for visual representation of total errors as boxes
    g_vn_dir_comb_toterr.SetPointEXlow(i, 0.08);

    // combined v2dir with statistical error bars
    h_vn_dir_comb_staterr[i]->GetQuantiles(n_quant, quantiles, prob_integrals);
    err_high = quantiles[2] - central_value;
    err_low = central_value - quantiles[0];
    g_vn_dir_comb_staterr.SetPoint(i, pt(i), central_value);  // central value same as for g_vn_dir_comb_toterr
    g_vn_dir_comb_staterr.SetPointEYhigh(i, err_high);
    g_vn_dir_comb_staterr.SetPointEYlow(i, err_low);

    // store inclusive v2 values
    g_vn_inc_pcm_toterr.SetPoint(i, pt(i), v2_inc_meas_values_pcm(i));
    g_vn_inc_pcm_toterr.SetPointEYhigh(i, TMath::Sqrt(cov_v2_inc_toterr_pcm(i,i)));
    g_vn_inc_pcm_toterr.SetPointEYlow(i, TMath::Sqrt(cov_v2_inc_toterr_pcm(i,i)));

    g_vn_inc_phos_toterr.SetPoint(i, pt(i), v2_inc_meas_values_phos(i));
    g_vn_inc_phos_toterr.SetPointEYhigh(i, TMath::Sqrt(cov_v2_inc_toterr_phos(i,i)));
    g_vn_inc_phos_toterr.SetPointEYlow(i, TMath::Sqrt(cov_v2_inc_toterr_phos(i,i)));

    g_vn_inc_comb_toterr.SetPoint(i, pt(i), v2_inc_meas_values_comb(i));
    g_vn_inc_comb_toterr.SetPointEYhigh(i, TMath::Sqrt(cov_v2_inc_toterr_comb(i,i)));
    g_vn_inc_comb_toterr.SetPointEYlow(i, TMath::Sqrt(cov_v2_inc_toterr_comb(i,i)));

    // store decay photon v2
    g_vn_dec_comb.SetPoint(i, pt(i), v2_dec_meas_values(i));
    g_vn_dec_comb.SetPointEYhigh(i, TMath::Sqrt(cov_v2_dec_syserr(i,i)));
    g_vn_dec_comb.SetPointEYlow(i, TMath::Sqrt(cov_v2_dec_syserr(i,i)));

  }
  TString filename_c1 = "v2dir_distr_" + centr + ".pdf";
  c1->SaveAs(filename_c1);

  //
  // Plot PCM, PHOS, and combined v2dir vs pT
  //
  TCanvas* c2 = new TCanvas("c2","c2", 700, 10, 700, 500);
  TH2F frame("frame","direct photon v_{2}", 1, 0., 6, 1, -0.15, 0.35);
  frame.SetXTitle("p_{T} (GeV/#it{c})");
  frame.SetYTitle("direct photon v_{2}");

  frame.DrawCopy();

  TLine ly0(0., 0., 6., 0.);
  ly0.DrawClone();

  g_vn_dir_comb_toterr.SetFillColor(kBlue-9);

  g_vn_dir_pcm_toterr.SetMarkerStyle(kFullCircle);
  g_vn_dir_phos_toterr.SetMarkerStyle(kFullCircle);
  g_vn_dir_comb_toterr.SetMarkerStyle(kFullCircle);

  g_vn_dir_pcm_toterr.SetLineColor(kRed);
  g_vn_dir_phos_toterr.SetLineColor(kBlue);
  g_vn_dir_comb_toterr.SetLineColor(kBlack);

  g_vn_dir_pcm_toterr.SetMarkerColor(kRed);
  g_vn_dir_phos_toterr.SetMarkerColor(kBlue);
  g_vn_dir_comb_toterr.SetMarkerColor(kBlack);

  g_vn_dir_pcm_toterr.DrawClone("p");
  g_vn_dir_phos_toterr.DrawClone("p");
  g_vn_dir_comb_toterr.DrawClone("p");

  TLatex l2;
  l2.DrawLatexNDC(0.2, 0.75, Form("%s%%",centr.Data()));

  TString filename_c2 = "v2dir_pcm_phos_comb_" + centr + ".pdf";
  c2->SaveAs(filename_c2);

  //
  // Plot only combined v2dir vs pT
  //
  TCanvas* c3 = new TCanvas("c3","c3", 700, 500, 700, 500);
  frame.DrawCopy();

  ly0.DrawClone();

  g_vn_dir_comb_toterr.SetMarkerStyle(kFullCircle);
  g_vn_dir_comb_toterr.SetLineColor(kBlack);
  g_vn_dir_comb_toterr.SetMarkerColor(kBlack);
  g_vn_dir_comb_toterr.DrawClone("pE5");

  g_vn_dir_comb_staterr.SetMarkerStyle(kFullCircle);
  g_vn_dir_comb_staterr.SetLineColor(kBlack);
  g_vn_dir_comb_staterr.SetMarkerColor(kBlack);
  g_vn_dir_comb_staterr.DrawClone("p");

  TLatex l3;
  l3.DrawLatexNDC(0.2, 0.75, Form("%s%%",centr.Data()));

  TString filename_c3 = "v2dir_comb_" + centr + ".pdf";
  c3->SaveAs(filename_c3);

  //
  // Plot inclusive and decay photon v2
  //
  TCanvas* c4 = new TCanvas("c4","c4", 900, 10, 700, 500);
  frame.SetXTitle("p_{T} (GeV/#it{c})");
  frame.SetYTitle("inclusive photon v_{2}, decay photon v_{2}");

  frame.DrawCopy();

  g_vn_dec_comb.SetMarkerStyle(kOpenCircle);
  g_vn_dec_comb.SetLineColor(kGray);
  g_vn_dec_comb.SetMarkerColor(kGray);
  g_vn_dec_comb.SetFillColor(kGray);
  g_vn_dec_comb.DrawClone("3l");

  g_vn_inc_pcm_toterr.SetMarkerStyle(kFullCircle);
  g_vn_inc_phos_toterr.SetMarkerStyle(kFullCircle);
  g_vn_inc_comb_toterr.SetMarkerStyle(kFullCircle);

  g_vn_inc_pcm_toterr.SetLineColor(kRed);
  g_vn_inc_phos_toterr.SetLineColor(kBlue);
  g_vn_inc_comb_toterr.SetLineColor(kBlack);

  g_vn_inc_pcm_toterr.SetMarkerColor(kRed);
  g_vn_inc_phos_toterr.SetMarkerColor(kBlue);
  g_vn_inc_comb_toterr.SetMarkerColor(kBlack);

  g_vn_inc_pcm_toterr.DrawClone("p");
  g_vn_inc_phos_toterr.DrawClone("p");
  g_vn_inc_comb_toterr.DrawClone("p");

  g_vn_dec_comb.DrawClone("p");

  TLatex l4;
  l4.DrawLatexNDC(0.2, 0.75, Form("%s%%",centr.Data()));

  TString filename_c4 = "v2inc_pcm_phos_comb_" + centr + ".pdf";
  c4->SaveAs(filename_c4);


}
