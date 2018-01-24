//
// This code routines addresses the following tasks
// - Calculation of v2,dir from v2,inv, v2,dec, and Rgamma
// - Averaging PCM and PHOS results
// - Quantifying the level of agreement with certain hypotheses about v2,dir (to be done)
//
//
// Klaus Reygers, November 2016
//

#include "Math/BrentRootFinder.h"
#include "Math/WrappedTF1.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TVectorD.h"
// #include "Math/Functions.h"

#include "RooArgList.h"
#include "RooDataSet.h"
#include "RooMultiVarGaussian.h"
#include "RooRealVar.h"

//
// used functions
//
void v2_dir_central_value_and_error_from_v2dir_distribution(TH1 *h_v2dir_distr, Double_t &v2dir_central_value,
                                                            Double_t &v2dir_err_up, Double_t &v2dir_err_low);

void fill_v2dir_graph(TH1 **hp, const TVectorD &pt, const Double_t &pt_offset, const Double_t &err_x,
                      TGraphAsymmErrors &g);

void weighted_average(const Int_t nMeas, const TVectorD **mu, const TMatrixDSym **cov, TVectorD *muAve,
                      TMatrixDSym *covAve);

void weighted_average(const TVectorD &mu1, const TVectorD &mu2, const TMatrixDSym &cov1, const TMatrixDSym &cov2,
                      TVectorD &muAve, TMatrixDSym &covAve);

void likelihood(const Int_t nPtBins, const Int_t nSamples, const TVectorD &vnDirTrueVec, const TVectorD &RgamMeasVec,
                const TMatrixDSym &covRgam, const TVectorD &vnIncMeasVec, const TMatrixDSym &covVnInc,
                const TVectorD &vnDecMeasVec, const TMatrixDSym &covVnDec, Double_t &Lfinal);

void v2dir(const Int_t nPtBins, const Int_t nSamples, const TVectorD &RgamMeasVec, const TMatrixDSym &covRgam,
           const TVectorD &vnIncMeasVec, const TMatrixDSym &covVnInc, const TVectorD &vnDecMeasVec,
           const TMatrixDSym &covVnDec, TH1 **hVnDir, TString prefix, TMatrixDSym &covarianceMatrix,
           TMatrixDSym &correlationMatrix);

void v2dir_classic(const Int_t nPtBins, TVectorD &pt, const TVectorD &RgamMeasVec, const TMatrixDSym &covRgam,
                   const TVectorD &vnIncMeasVec, const TMatrixDSym &covVnInc, const TVectorD &vnDecMeasVec,
                   const TMatrixDSym &covVnDec, TGraphErrors *g_v2dir_classic);

Double_t DotProduct(const TVectorD &a, const TVectorD &b);

Double_t chi2(const TVectorD &x, const TVectorD &mu, const TMatrixDSym &V);

Double_t p_value_to_n_sigma(const Double_t &pvalue);

// number of pT bins
const Int_t nPtBins = 16;

//
// the main function and entry point
//
void v2dir_pcm_phos_comb(TString centr, Int_t nSamples = 100000, TString output_dir = ".") {

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

    // print debug info?
    bool verbose = false;

    if (centr != "00-20" && centr != "20-40" && centr != "40-80") {
        cout << "Allowed centrality classes are \"00-20\", \"20-40\", and \"40-80\"" << endl;
        return;
    }

    output_dir += "/";

    // output file name
    TString fn_out(output_dir + "v2dir_pcm_phos_comb_" + centr + ".root");

    //
    // PCM specific input
    //
    TString fn_pcm = "data_pcm/data_pcm.root";
    TFile f_pcm(fn_pcm.Data());
    TDirectory *dir_pcm = (TDirectory *)f_pcm.Get(centr.Data());

    TVectorD v2_inc_meas_values_pcm = *(TVectorD *)dir_pcm->Get("v2_inc_meas_values");
    TMatrixDSym cov_v2_inc_toterr_pcm = *(TMatrixDSym *)dir_pcm->Get("cov_v2_inc_toterr");
    TMatrixDSym cov_v2_inc_staterr_pcm = *(TMatrixDSym *)dir_pcm->Get("cov_v2_inc_staterr");
    TMatrixDSym cov_v2_inc_syserr_pcm = *(TMatrixDSym *)dir_pcm->Get("cov_v2_inc_syserr");

    // only needed for cross checks
    TVectorD Rgamma_meas_vec_pcm = *(TVectorD *)dir_pcm->Get("Rgamma_meas_vec");
    TMatrixDSym cov_Rgamma_syserr_pcm = *(TMatrixDSym *)dir_pcm->Get("cov_Rgamma_syserr");
    TMatrixDSym cov_Rgamma_staterr_pcm = *(TMatrixDSym *)dir_pcm->Get("cov_Rgamma_staterr");
    TMatrixDSym cov_Rgamma_toterr_pcm = *(TMatrixDSym *)dir_pcm->Get("cov_Rgamma_toterr");

    //
    // PHOS specific input
    //
    TString fn_phos = "data_phos/data_phos.root";
    TFile f_phos(fn_phos.Data());
    TDirectory *dir_phos = (TDirectory *)f_phos.Get(centr.Data());

    TVectorD v2_inc_meas_values_phos = *(TVectorD *)dir_phos->Get("v2_inc_meas_values");
    TMatrixDSym cov_v2_inc_toterr_phos = *(TMatrixDSym *)dir_phos->Get("cov_v2_inc_toterr");
    TMatrixDSym cov_v2_inc_syserr_phos = *(TMatrixDSym *)dir_phos->Get("cov_v2_inc_syserr");
    TMatrixDSym cov_v2_inc_staterr_phos = *(TMatrixDSym *)dir_phos->Get("cov_v2_inc_staterr");

    // only needed for cross checks
    TVectorD Rgamma_meas_vec_phos = *(TVectorD *)dir_phos->Get("Rgamma_meas_vec");
    TMatrixDSym cov_Rgamma_staterr_phos = *(TMatrixDSym *)dir_phos->Get("cov_Rgamma_staterr");
    TMatrixDSym cov_Rgamma_syserr_phos = *(TMatrixDSym *)dir_phos->Get("cov_Rgamma_syserr");
    TMatrixDSym cov_Rgamma_toterr_phos = *(TMatrixDSym *)dir_phos->Get("cov_Rgamma_toterr");

    //
    // PCM/PHOS combined
    //
    TVectorD Rgamma_meas_vec_comb = *(TVectorD *)dir_pcm->Get("Rgamma_meas_vec_comb");
    TMatrixDSym cov_Rgamma_comb_toterr = *(TMatrixDSym *)dir_pcm->Get("cov_Rgamma_comb_toterr");
    TMatrixDSym cov_Rgamma_comb_staterr = *(TMatrixDSym *)dir_pcm->Get("cov_Rgamma_comb_staterr");

    // decay photon v2
    TString fn_cocktail = "cocktail/cocktail.root";
    TFile f_cocktail(fn_cocktail.Data());
    TDirectory *dir_cocktail = (TDirectory *)f_cocktail.Get(centr.Data());

    TVectorD v2_dec_meas_values = *(TVectorD *)dir_cocktail->Get("v2_dec_meas_values");
    TMatrixDSym cov_v2_dec_staterr = *(TMatrixDSym *)dir_cocktail->Get("cov_v2_dec_staterr");
    TMatrixDSym cov_v2_dec_syserr = *(TMatrixDSym *)dir_cocktail->Get("cov_v2_dec_syserr");
    TMatrixDSym cov_v2_dec_toterr = *(TMatrixDSym *)dir_cocktail->Get("cov_v2_dec_toterr");

    // decay photon v2 from PHOS file (default method)
    // TVectorD v2_dec_meas_values = *(TVectorD*) dir_phos->Get("v2_dec_meas_values");
    // TMatrixDSym cov_v2_dec_staterr = *(TMatrixDSym*) dir_phos->Get("cov_v2_dec_staterr");
    // TMatrixDSym cov_v2_dec_syserr = *(TMatrixDSym*) dir_phos->Get("cov_v2_dec_syserr");
    // TMatrixDSym cov_v2_dec_toterr = *(TMatrixDSym*) dir_phos->Get("cov_v2_dec_toterr");

    TVectorD pt = *(TVectorD *)dir_pcm->Get("pt");
    if (verbose) {
        cout << "pt central values: " << endl;
        pt.Print();
    }

    //
    // create array of pointers to output histograms (vndir distributions for each pT bin)
    //
    TH1 *h_v2_dir_pcm_Rpcm_toterr[nPtBins];
    TH1 *h_v2_dir_pcm_toterr[nPtBins];
    TH1 *h_v2_dir_phos_Rphos_toterr[nPtBins];
    TH1 *h_v2_dir_phos_toterr[nPtBins];
    TH1 *h_v2_dir_comb_staterr[nPtBins];
    TH1 *h_v2_dir_comb_toterr[nPtBins];
    TH1 *h_v2_dir_comb_toterr_uncorr[nPtBins];
    TH1 *h_v2_dir_comb_staterr_uncorr[nPtBins];
    TH1 *h_v2_dir_comb_toterr_xcheck1[nPtBins];
    TH1 *h_v2_dir_comb_toterr_xcheck2[nPtBins];

    //
    // covariance and correlation matrices
    //
    TMatrixDSym cov_v2_dir_pcm_Rpcm_toterr(nPtBins);
    TMatrixDSym cov_v2_dir_pcm_toterr(nPtBins);
    TMatrixDSym cov_v2_dir_phos_Rphos_toterr(nPtBins);
    TMatrixDSym cov_v2_dir_phos_toterr(nPtBins);
    TMatrixDSym cov_v2_dir_comb_staterr(nPtBins);
    TMatrixDSym cov_v2_dir_comb_toterr(nPtBins);
    TMatrixDSym cov_v2_dir_comb_toterr_uncorr(nPtBins);
    TMatrixDSym cov_v2_dir_comb_staterr_uncorr(nPtBins);
    TMatrixDSym cov_v2_dir_comb_toterr_xcheck1(nPtBins);
    TMatrixDSym cov_v2_dir_comb_toterr_xcheck2(nPtBins);

    TMatrixDSym corr_v2_dir_pcm_Rpcm_toterr(nPtBins);
    TMatrixDSym corr_v2_dir_pcm_toterr(nPtBins);
    TMatrixDSym corr_v2_dir_phos_Rphos_toterr(nPtBins);
    TMatrixDSym corr_v2_dir_phos_toterr(nPtBins);
    TMatrixDSym corr_v2_dir_comb_staterr(nPtBins);
    TMatrixDSym corr_v2_dir_comb_toterr(nPtBins);
    TMatrixDSym corr_v2_dir_comb_toterr_uncorr(nPtBins);
    TMatrixDSym corr_v2_dir_comb_staterr_uncorr(nPtBins);
    TMatrixDSym corr_v2_dir_comb_toterr_xcheck1(nPtBins);
    TMatrixDSym corr_v2_dir_comb_toterr_xcheck2(nPtBins);

    //
    // calculate vnDir PCM (total errors, using PCM Rgamma and common v2dec)
    //
    v2dir(nPtBins, nSamples, Rgamma_meas_vec_pcm, cov_Rgamma_toterr_pcm, v2_inc_meas_values_pcm, cov_v2_inc_toterr_pcm,
          v2_dec_meas_values, cov_v2_dec_toterr, h_v2_dir_pcm_Rpcm_toterr, "h_v2Dir_PCM_TotErr_RgamPCM_",
          cov_v2_dir_pcm_Rpcm_toterr, corr_v2_dir_pcm_Rpcm_toterr);

    //
    // calculate vnDir PCM (total errors, using combined Rgamma and common v2dec)
    //
    v2dir(nPtBins, nSamples, Rgamma_meas_vec_comb, cov_Rgamma_comb_toterr, v2_inc_meas_values_pcm,
          cov_v2_inc_toterr_pcm, v2_dec_meas_values, cov_v2_dec_toterr, h_v2_dir_pcm_toterr,
          "h_v2Dir_PCM_TotErr_RgamComb_", cov_v2_dir_pcm_toterr, corr_v2_dir_pcm_toterr);

    //
    // calculate vnDir PHOS (total errors, using PHOS Rgamma and common v2dec)
    //
    v2dir(nPtBins, nSamples, Rgamma_meas_vec_phos, cov_Rgamma_toterr_phos, v2_inc_meas_values_phos,
          cov_v2_inc_toterr_phos, v2_dec_meas_values, cov_v2_dec_toterr, h_v2_dir_phos_Rphos_toterr,
          "h_v2Dir_PHOS_TotErr_RgamPhos_", cov_v2_dir_phos_Rphos_toterr, corr_v2_dir_phos_Rphos_toterr);

    //
    // calculate vnDir PHOS (total errors, using combined Rgamma and common v2dec)
    //
    v2dir(nPtBins, nSamples, Rgamma_meas_vec_comb, cov_Rgamma_comb_toterr, v2_inc_meas_values_phos,
          cov_v2_inc_toterr_phos, v2_dec_meas_values, cov_v2_dec_toterr, h_v2_dir_phos_toterr,
          "h_v2Dir_PHOS_TotErr_RgamComb_", cov_v2_dir_phos_toterr, corr_v2_dir_phos_toterr);

    //
    // calculate weighted average of v2inc
    //

    // weighted average of Rgamma
    // TVectorD Rgamma_meas_vec_comb(nPtBins);
    // TMatrixDSym cov_Rgamma_comb(nPtBins);
    // weighted_average(Rgamma_meas_vec_syserr_pcm, Rgamma_meas_vec_phos, cov_Rgamma_syserr_pcm, cov_Rgamma_phos,
    // Rgamma_meas_vec_comb, cov_Rgamma_comb);

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
    // As a cross check calculate results assuming no correlations of systematic errors
    //

    // covariance matrix for uncorrelated total errors
    TMatrixDSym cov_v2_inc_toterr_pcm_uncorr(nPtBins);
    TMatrixDSym cov_v2_inc_toterr_phos_uncorr(nPtBins);
    TMatrixDSym cov_Rgamma_comb_toterr_uncorr(nPtBins);
    TMatrixDSym cov_v2_dec_toterr_uncorr(nPtBins);
    for (Int_t i = 0; i < nPtBins; i++) {
        cov_v2_inc_toterr_pcm_uncorr(i, i) = cov_v2_inc_toterr_pcm(i, i);
        cov_v2_inc_toterr_phos_uncorr(i, i) = cov_v2_inc_toterr_phos(i, i);
        cov_Rgamma_comb_toterr_uncorr(i, i) = cov_Rgamma_comb_toterr(i, i);
        cov_v2_dec_toterr_uncorr(i, i) = cov_v2_dec_toterr(i, i);
    }

    // v2,in weighted average using total error as weight, assuming no correlation of errors in pt
    TVectorD v2_inc_meas_values_comb_uncorr(nPtBins);
    TMatrixDSym cov_v2_inc_toterr_comb_uncorr(nPtBins);
    weighted_average(v2_inc_meas_values_pcm, v2_inc_meas_values_phos, cov_v2_inc_toterr_pcm_uncorr,
                     cov_v2_inc_toterr_phos_uncorr, v2_inc_meas_values_comb_uncorr, cov_v2_inc_toterr_comb_uncorr);

    // for debugging purposes
    if (verbose) {
        for (Int_t i = 0; i < nPtBins; i++) {
            Double_t v2_inc_staterr = TMath::Sqrt(cov_v2_inc_staterr_comb(i, i));
            Double_t v2_inc_toterr = TMath::Sqrt(cov_v2_inc_toterr_comb_uncorr(i, i));
            cout << "v2_inc_staterr: " << v2_inc_staterr;
            cout << ", v2_inc_toterr: " << v2_inc_toterr;
            cout << ", v2_inc_staterr/v2_inc_toterr: " << v2_inc_staterr / v2_inc_toterr << endl;
        }
    };

    //
    // calculate vnDir combined
    //

    // total error: the main result
    v2dir(nPtBins, nSamples, Rgamma_meas_vec_comb, cov_Rgamma_comb_toterr, v2_inc_meas_values_comb,
          cov_v2_inc_toterr_comb, v2_dec_meas_values, cov_v2_dec_toterr, h_v2_dir_comb_toterr,
          "h_v2Dir_Comb_totErr_RgamComb_", cov_v2_dir_comb_toterr, corr_v2_dir_comb_toterr);

    if (verbose) {
        cout << "--> v2,dir covariance matrix:" << endl;
        cov_v2_dir_comb_toterr.Print();

        cout << "--> v2,dir correlation matrix:" << endl;
        corr_v2_dir_comb_toterr.Print();
    }

    // stat. err.
    // note that we use  cov_Rgamma_comb_toterr rather than cov_Rgamma_comb_staterr
    // which actually makes the statistical error of v2,dir smaller(!)
    // Otherwise we would get a statistical error that is larger than the systematic error
    // which would not be usedul
    v2dir(nPtBins, nSamples, Rgamma_meas_vec_comb, cov_Rgamma_comb_toterr, v2_inc_meas_values_comb,
          cov_v2_inc_staterr_comb, v2_dec_meas_values, cov_v2_dec_staterr, h_v2_dir_comb_staterr,
          "h_v2Dir_Comb_StatErr_RgamComb_", cov_v2_dir_comb_staterr, corr_v2_dir_comb_staterr);

    //
    // As a cross check, do the calcualtion assuming uncorrelated systematic errors
    //

    // total error, uncorrelated errors (cross check)
    v2dir(nPtBins, nSamples, Rgamma_meas_vec_comb, cov_Rgamma_comb_toterr_uncorr, v2_inc_meas_values_comb_uncorr,
          cov_v2_inc_toterr_comb_uncorr, v2_dec_meas_values, cov_v2_dec_toterr_uncorr, h_v2_dir_comb_toterr_uncorr,
          "h_v2Dir_Comb_totErr_RgamComb_uncorr_", cov_v2_dir_comb_toterr_uncorr, corr_v2_dir_comb_toterr_uncorr);

    // stat. err. (cross check)
    v2dir(nPtBins, nSamples, Rgamma_meas_vec_comb, cov_Rgamma_comb_toterr, v2_inc_meas_values_comb_uncorr,
          cov_v2_inc_staterr_comb, v2_dec_meas_values, cov_v2_dec_staterr, h_v2_dir_comb_staterr_uncorr,
          "h_v2Dir_Comb_StatErr_RgamComb_uncorr", cov_v2_dir_comb_staterr_uncorr, corr_v2_dir_comb_staterr_uncorr);

    //
    // Another cross check: calculate v2,dir for a fixed given Rgamma
    //
    Double_t Rgamma_xcheck1_val = 1.03;
    TVectorD Rgamma_xcheck1(nPtBins);
    for (Int_t i = 0; i < nPtBins; i++)
        Rgamma_xcheck1(i) = Rgamma_xcheck1_val;
    v2dir(nPtBins, nSamples, Rgamma_xcheck1, cov_Rgamma_comb_toterr, v2_inc_meas_values_comb, cov_v2_inc_staterr_comb,
          v2_dec_meas_values, cov_v2_dec_staterr, h_v2_dir_comb_toterr_xcheck1, "h_v2Dir_Comb_totErr_RgamComb_xcheck1_",
          cov_v2_dir_comb_toterr_xcheck1, corr_v2_dir_comb_toterr_xcheck1);

    //
    // A further cross check:
    // Scale Rgamma by the ratio of the data-driven and the MC-based purity (Fig. 10 in Mike's note)
    //

    // purity ratio: data-driven / MC-based
    Double_t purity_ratio_020[nPtBins] = {0.944581, 0.930561, 0.940445, 0.948759, 0.956215, 0.96239,
                                          0.967351, 0.971467, 0.975393, 0.979183, 0.9837,   0.986301,
                                          0.988866, 0.991039, 0.993376, 0.994887};
    Double_t purity_ratio_2040[nPtBins] = {0.97793,  0.978613, 0.982284, 0.98585,  0.988225, 0.990048,
                                           0.991743, 0.992658, 0.993512, 0.994843, 0.995635, 0.996243,
                                           0.996793, 0.997282, 0.997844, 0.998199};

    TVectorD Rgamma_xcheck2(nPtBins);
    for (Int_t i = 0; i < nPtBins; i++) {
        Double_t sf = 1;
        if (centr == "00-20")
            sf = purity_ratio_020[i];
        else if (centr == "20-40")
            sf = purity_ratio_2040[i];
        Rgamma_xcheck2(i) = sf * Rgamma_meas_vec_comb(i);
    }

    v2dir(nPtBins, nSamples, Rgamma_xcheck2, cov_Rgamma_comb_toterr, v2_inc_meas_values_comb, cov_v2_inc_toterr_comb,
          v2_dec_meas_values, cov_v2_dec_toterr, h_v2_dir_comb_toterr_xcheck2, "h_v2Dir_Comb_totErr_RgamComb_xcheck2_",
          cov_v2_dir_comb_toterr_xcheck2, corr_v2_dir_comb_toterr_xcheck2);

    // Yet another cross check: calculate v2dir central values directly from measured values
    TGraphErrors g_v2_dir_comb_classic_xcheck(nPtBins);
    v2dir_classic(nPtBins, pt, Rgamma_meas_vec_comb, cov_Rgamma_comb_toterr, v2_inc_meas_values_comb,
                  cov_v2_inc_toterr_comb, v2_dec_meas_values, cov_v2_dec_toterr, &g_v2_dir_comb_classic_xcheck);

    // graphs to store v2dir
    TGraphAsymmErrors g_v2_dir_pcm_Rpcm_toterr(nPtBins);
    TGraphAsymmErrors g_v2_dir_pcm_toterr(nPtBins);
    TGraphAsymmErrors g_v2_dir_phos_Rphos_toterr(nPtBins);
    TGraphAsymmErrors g_v2_dir_phos_toterr(nPtBins);
    TGraphAsymmErrors g_v2_dir_comb_toterr(nPtBins);
    TGraphAsymmErrors g_v2_dir_comb_staterr(nPtBins); // central values same as g_v2_dir_comb_toterr
    TGraphAsymmErrors g_v2_dir_comb_toterr_uncorr(nPtBins);
    TGraphAsymmErrors g_v2_dir_comb_staterr_uncorr(nPtBins); // central values same as g_v2_dir_comb_toterr_uncorr
    TGraphAsymmErrors g_v2_dir_comb_toterr_xcheck1(nPtBins);
    TGraphAsymmErrors g_v2_dir_comb_toterr_xcheck2(nPtBins);

    // graphs to store v2inc
    TGraphAsymmErrors g_v2_inc_pcm_toterr(nPtBins);
    TGraphAsymmErrors g_v2_inc_pcm_syserr(nPtBins);
    TGraphAsymmErrors g_v2_inc_pcm_staterr(nPtBins);

    TGraphAsymmErrors g_v2_inc_phos_toterr(nPtBins);
    TGraphAsymmErrors g_v2_inc_phos_syserr(nPtBins);
    TGraphAsymmErrors g_v2_inc_phos_staterr(nPtBins);

    TGraphAsymmErrors g_v2_inc_comb_toterr(nPtBins);
    TGraphAsymmErrors g_v2_inc_comb_syserr(nPtBins);
    TGraphAsymmErrors g_v2_inc_comb_staterr(nPtBins);

    // graph to store decay photon v2
    TGraphAsymmErrors g_v2_dec_comb(nPtBins);

    //
    // fill v2dir graphs
    //
    // Double_t pt_offset_pcm = 0.03;
    // Double_t pt_offset_phos = 0.06;
    Double_t pt_offset_pcm = 0.;
    Double_t pt_offset_phos = 0.;
    Double_t err_pt = 0.08; // for visual representation of errors as boxes
    fill_v2dir_graph(h_v2_dir_pcm_Rpcm_toterr, pt, pt_offset_pcm, 0, g_v2_dir_pcm_Rpcm_toterr);
    fill_v2dir_graph(h_v2_dir_phos_Rphos_toterr, pt, pt_offset_phos, 0, g_v2_dir_phos_Rphos_toterr);
    fill_v2dir_graph(h_v2_dir_pcm_toterr, pt, pt_offset_pcm, 0, g_v2_dir_pcm_toterr);
    fill_v2dir_graph(h_v2_dir_phos_toterr, pt, pt_offset_phos, 0, g_v2_dir_phos_toterr);
    fill_v2dir_graph(h_v2_dir_comb_toterr, pt, 0, err_pt, g_v2_dir_comb_toterr);
    fill_v2dir_graph(h_v2_dir_comb_staterr, pt, 0, 0, g_v2_dir_comb_staterr);
    fill_v2dir_graph(h_v2_dir_comb_toterr_uncorr, pt, 0, err_pt, g_v2_dir_comb_toterr_uncorr);
    fill_v2dir_graph(h_v2_dir_comb_staterr_uncorr, pt, 0, err_pt, g_v2_dir_comb_staterr_uncorr);
    fill_v2dir_graph(h_v2_dir_comb_toterr_xcheck1, pt, 0, err_pt, g_v2_dir_comb_toterr_xcheck1);
    fill_v2dir_graph(h_v2_dir_comb_toterr_xcheck2, pt, 0, err_pt, g_v2_dir_comb_toterr_xcheck2);

    // set central values of graphs with statistical error to central values of graphs with total error
    for (Int_t i = 0; i < nPtBins; i++) {
        Double_t pt_val = g_v2_dir_comb_toterr.GetX()[i];
        Double_t v2dir_val = g_v2_dir_comb_toterr.GetY()[i];
        g_v2_dir_comb_staterr.SetPoint(i, pt_val, v2dir_val);

        pt_val = g_v2_dir_comb_toterr_uncorr.GetX()[i];
        v2dir_val = g_v2_dir_comb_toterr_uncorr.GetY()[i];
        g_v2_dir_comb_staterr_uncorr.SetPoint(i, pt_val, v2dir_val);
    }

    //
    // Yet another check:
    // calculate weighted average of v2,dir(PCM, RgammaPCM) and v2,dir(PHOS, RgammaPHOS).
    // To be compared with the combined v2mdir of the default method
    //
    TVectorD v2_dir_pcm_Rpcm_toterr(nPtBins);
    TVectorD v2_dir_phos_Rphos_toterr(nPtBins);
    for (Int_t i = 0; i < nPtBins; i++) {
        v2_dir_pcm_Rpcm_toterr(i) = g_v2_dir_pcm_Rpcm_toterr.GetY()[i];
        v2_dir_phos_Rphos_toterr(i) = g_v2_dir_phos_Rphos_toterr.GetY()[i];
    }

    // the result of the simple averaging
    TVectorD v2_dir_comb_ave_xcheck(nPtBins);
    TMatrixDSym cov_v2_dir_comb_ave_xcheck(nPtBins);

    // weighted average of v2,dir(PCM, RgammaPCM) and v2,dir(PHOS, RgammaPHOS)
    weighted_average(v2_dir_pcm_Rpcm_toterr, v2_dir_phos_Rphos_toterr, cov_v2_dir_pcm_Rpcm_toterr,
                     cov_v2_dir_phos_Rphos_toterr, v2_dir_comb_ave_xcheck, cov_v2_dir_comb_ave_xcheck);

    // store result in graph
    TGraphAsymmErrors g_v2_dir_comb_ave_xcheck(nPtBins);
    for (Int_t i = 0; i < nPtBins; i++) {
        g_v2_dir_comb_ave_xcheck.SetPoint(i, pt(i), v2_dir_comb_ave_xcheck(i));
        Double_t toterr = TMath::Sqrt(cov_v2_dir_comb_ave_xcheck(i, i));
        g_v2_dir_comb_ave_xcheck.SetPointError(i, err_pt, err_pt, toterr, toterr);
    }

    //
    // plot output
    //
    cout << "#entries = " << h_v2_dir_pcm_toterr[0]->GetEntries() << endl;

    // style settings
    TStyle *myStyle = new TStyle("myStyle", "My root style");
    myStyle->SetOptStat(kFALSE);
    myStyle->SetLabelOffset(0.005, "x"); // 0.005 = root default
    myStyle->SetLabelOffset(0.005, "y"); // 0.005 = root default
    myStyle->SetTitleXSize(0.04);        // 0.04  = root default
    myStyle->SetTitleYSize(0.04);        // 0.04  = root default
    myStyle->SetTitleXOffset(1.2);
    myStyle->SetTitleYOffset(1.3);
    myStyle->SetPadLeftMargin(0.15);
    myStyle->SetPadRightMargin(0.12); // 0.1 = root default
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
    TCanvas *c1 = new TCanvas("c1", "c1", 700, 900);
    c1->Divide(4, 4);

    TLatex lpt;

    for (Int_t i = 0; i < nPtBins; i++) {
        c1->cd(i + 1);
        c1->GetPad(i + 1)->SetBottomMargin(0.16);
        c1->GetPad(i + 1)->SetLeftMargin(0.16);

        h_v2_dir_comb_toterr[i]->SetXTitle("v2,dir");

        h_v2_dir_comb_toterr[i]->SetTitle(Form("v2dir, combined, pT bin %d, pT = %4.1f GeV/c", i, pt(i)));
        h_v2_dir_comb_toterr[i]->SetName(Form("h_v2_dir_comb_toterr_ptbin_%d", i));

        h_v2_dir_pcm_toterr[i]->SetTitle(Form("v2dir, PCM, combined Rgamma, pT bin %d, pT = %4.1f GeV/c", i, pt(i)));
        h_v2_dir_pcm_toterr[i]->SetName(Form("h_v2_dir_pcm_toterr_ptbin_%d", i));

        h_v2_dir_phos_toterr[i]->SetTitle(Form("v2dir, PHOS, combined Rgamma, pT bin %d, pT = %4.1f GeV/c", i, pt(i)));
        h_v2_dir_phos_toterr[i]->SetName(Form("h_v2_dir_phos_toterr_ptbin_%d", i));

        h_v2_dir_pcm_Rpcm_toterr[i]->SetTitle(
            Form("v2dir, PCM, Rgamma from PCM, pT bin %d, pT = %4.1f GeV/c", i, pt(i)));
        h_v2_dir_pcm_Rpcm_toterr[i]->SetName(Form("h_v2_dir_pcm_Rpcm_toterr_ptbin_%d", i));

        h_v2_dir_phos_Rphos_toterr[i]->SetTitle(
            Form("v2dir, PHOS, Rgamma from PHOS, pT bin %d, pT = %4.1f GeV/c", i, pt(i)));
        h_v2_dir_phos_Rphos_toterr[i]->SetName(Form("h_v2_dir_phos_Rphos_toterr_ptbin_%d", i));

        h_v2_dir_comb_toterr[i]->SetLineColor(kBlack);
        h_v2_dir_comb_toterr[i]->DrawCopy("HL");
        h_v2_dir_pcm_Rpcm_toterr[i]->SetLineColor(kRed);
        h_v2_dir_pcm_Rpcm_toterr[i]->DrawCopy("HL,same");
        h_v2_dir_phos_Rphos_toterr[i]->SetLineColor(kBlue);
        h_v2_dir_phos_Rphos_toterr[i]->DrawCopy("HL,same");

        lpt.DrawLatexNDC(0.2, 0.75, Form("p_{T} = %4.1f GeV/#it{c}", pt(i)));

        //
        // store inclusive v2 values
        //

        // v2,inc, PCM, total error
        g_v2_inc_pcm_toterr.SetPoint(i, pt(i), v2_inc_meas_values_pcm(i));
        g_v2_inc_pcm_toterr.SetPointEYhigh(i, TMath::Sqrt(cov_v2_inc_toterr_pcm(i, i)));
        g_v2_inc_pcm_toterr.SetPointEYlow(i, TMath::Sqrt(cov_v2_inc_toterr_pcm(i, i)));

        // v2,inc, PCM, systematic error
        g_v2_inc_pcm_syserr.SetPoint(i, pt(i), v2_inc_meas_values_pcm(i));
        g_v2_inc_pcm_syserr.SetPointEYhigh(i, TMath::Sqrt(cov_v2_inc_syserr_pcm(i, i)));
        g_v2_inc_pcm_syserr.SetPointEYlow(i, TMath::Sqrt(cov_v2_inc_syserr_pcm(i, i)));

        // v2,inc, PCM, statistical error
        g_v2_inc_pcm_staterr.SetPoint(i, pt(i), v2_inc_meas_values_pcm(i));
        g_v2_inc_pcm_staterr.SetPointEYhigh(i, TMath::Sqrt(cov_v2_inc_staterr_pcm(i, i)));
        g_v2_inc_pcm_staterr.SetPointEYlow(i, TMath::Sqrt(cov_v2_inc_staterr_pcm(i, i)));

        // v2,inc, PHOS, total error
        g_v2_inc_phos_toterr.SetPoint(i, pt(i), v2_inc_meas_values_phos(i));
        g_v2_inc_phos_toterr.SetPointEYhigh(i, TMath::Sqrt(cov_v2_inc_toterr_phos(i, i)));
        g_v2_inc_phos_toterr.SetPointEYlow(i, TMath::Sqrt(cov_v2_inc_toterr_phos(i, i)));

        // v2,inc, PHOS, systematic error
        g_v2_inc_phos_syserr.SetPoint(i, pt(i), v2_inc_meas_values_phos(i));
        g_v2_inc_phos_syserr.SetPointEYhigh(i, TMath::Sqrt(cov_v2_inc_syserr_phos(i, i)));
        g_v2_inc_phos_syserr.SetPointEYlow(i, TMath::Sqrt(cov_v2_inc_syserr_phos(i, i)));

        // v2,inc, PHOS, statistical error
        g_v2_inc_phos_staterr.SetPoint(i, pt(i), v2_inc_meas_values_phos(i));
        g_v2_inc_phos_staterr.SetPointEYhigh(i, TMath::Sqrt(cov_v2_inc_staterr_phos(i, i)));
        g_v2_inc_phos_staterr.SetPointEYlow(i, TMath::Sqrt(cov_v2_inc_staterr_phos(i, i)));

        // v2,inc, combined, total error
        g_v2_inc_comb_toterr.SetPoint(i, pt(i), v2_inc_meas_values_comb(i));
        g_v2_inc_comb_toterr.SetPointEYhigh(i, TMath::Sqrt(cov_v2_inc_toterr_comb(i, i)));
        g_v2_inc_comb_toterr.SetPointEYlow(i, TMath::Sqrt(cov_v2_inc_toterr_comb(i, i)));

        // v2,inc, combined, statistical error
        g_v2_inc_comb_staterr.SetPoint(i, pt(i), v2_inc_meas_values_comb(i));
        g_v2_inc_comb_staterr.SetPointEYhigh(i, TMath::Sqrt(cov_v2_inc_staterr_comb(i, i)));
        g_v2_inc_comb_staterr.SetPointEYlow(i, TMath::Sqrt(cov_v2_inc_staterr_comb(i, i)));

        // v2,inc, combined, systematic error
        Double_t vn_inc_syserr = TMath::Sqrt(cov_v2_inc_toterr_comb(i, i) - cov_v2_inc_staterr_comb(i, i));
        g_v2_inc_comb_syserr.SetPoint(i, pt(i), v2_inc_meas_values_comb(i));
        g_v2_inc_comb_syserr.SetPointEYhigh(i, vn_inc_syserr);
        g_v2_inc_comb_syserr.SetPointEYlow(i, vn_inc_syserr);

        // store decay photon v2
        g_v2_dec_comb.SetPoint(i, pt(i), v2_dec_meas_values(i));
        g_v2_dec_comb.SetPointEYhigh(i, TMath::Sqrt(cov_v2_dec_syserr(i, i)));
        g_v2_dec_comb.SetPointEYlow(i, TMath::Sqrt(cov_v2_dec_syserr(i, i)));
    }
    TString filename_c1 = "v2dir_distr_" + centr + ".pdf";
    c1->SaveAs(output_dir + filename_c1);

    //
    // write root output file
    //
    TFile f_out(fn_out, "recreate");

    g_v2_dir_pcm_Rpcm_toterr.SetName("g_v2_dir_pcm_Rpcm_toterr");
    g_v2_dir_pcm_Rpcm_toterr.SetTitle("PCM v2dir using PCM Rgamma, total uncertainties");
    g_v2_dir_pcm_Rpcm_toterr.Write();

    g_v2_dir_pcm_toterr.SetName("g_v2_dir_pcm_Rcomb_toterr");
    g_v2_dir_pcm_toterr.SetTitle("PCM v2dir using combined Rgamma, total uncertainties");
    g_v2_dir_pcm_toterr.Write();

    g_v2_dir_phos_Rphos_toterr.SetName("g_v2_dir_phos_Rphos_toterr");
    g_v2_dir_phos_Rphos_toterr.SetTitle("PHOS v2dir using PHOS Rgamma, total uncertainties");
    g_v2_dir_phos_Rphos_toterr.Write();

    g_v2_dir_phos_toterr.SetName("g_v2_dir_phos_Rcomb_toterr");
    g_v2_dir_phos_toterr.SetTitle("PHOS v2dir using combined Rgamma, total uncertainties");
    g_v2_dir_phos_toterr.Write();

    g_v2_dir_comb_toterr.SetName("g_v2_dir_comb_toterr");
    g_v2_dir_comb_toterr.SetTitle("PCM/PHOS combined v2dir, total uncertainties");
    g_v2_dir_comb_toterr.Write();

    g_v2_dir_comb_toterr_uncorr.SetName("g_v2_dir_comb_toterr_uncorr");
    g_v2_dir_comb_toterr_uncorr.SetTitle("PCM/PHOS combined v2dir, total uncertainties, uncorrelated syst. errors");
    g_v2_dir_comb_toterr_uncorr.Write();

    g_v2_dir_comb_staterr.SetName("g_v2_dir_comb_staterr");
    g_v2_dir_comb_staterr.SetTitle("PCM/PHOS combined v2dir, stat. uncertainties");
    g_v2_dir_comb_staterr.Write();

    // write v2,dir covariance matrices
    cov_v2_dir_pcm_Rpcm_toterr.Write("cov_v2_dir_pcm_Rpcm_toterr");
    cov_v2_dir_pcm_toterr.Write("cov_v2_dir_pcm_Rcomb_toterr");
    cov_v2_dir_phos_Rphos_toterr.Write("cov_v2_dir_phos_Rphos_toterr");
    cov_v2_dir_phos_toterr.Write("cov_v2_dir_phos_Rcomb_toterr");
    cov_v2_dir_comb_toterr.Write("cov_v2_dir_comb_toterr");
    cov_v2_dir_comb_toterr_uncorr.Write("cov_v2_dir_comb_toterr_uncorr");
    cov_v2_dir_comb_staterr.Write("cov_v2_dir_comb_staterr");

    // write v2,dir correlation matrices
    corr_v2_dir_pcm_Rpcm_toterr.Write("corr_v2_dir_pcm_Rpcm_toterr");
    corr_v2_dir_pcm_toterr.Write("corr_v2_dir_pcm_Rcomb_toterr");
    corr_v2_dir_phos_Rphos_toterr.Write("corr_v2_dir_phos_Rphos_toterr");
    corr_v2_dir_phos_toterr.Write("corr_v2_dir_phos_Rcomb_toterr");
    corr_v2_dir_comb_toterr.Write("corr_v2_dir_comb_toterr");
    corr_v2_dir_comb_toterr_uncorr.Write("corr_v2_dir_comb_toterr_uncorr");
    corr_v2_dir_comb_staterr.Write("corr_v2_dir_comb_staterr");

    // inclusive photon v2
    g_v2_inc_pcm_toterr.SetName("g_v2_inc_pcm_toterr");
    g_v2_inc_pcm_toterr.SetTitle("inclusive photon v2 PCM, total uncertainties");
    g_v2_inc_pcm_toterr.Write();

    g_v2_inc_pcm_syserr.SetName("g_v2_inc_pcm_syserr");
    g_v2_inc_pcm_syserr.SetTitle("inclusive photon v2 PCM, systematic uncertainties");
    g_v2_inc_pcm_syserr.Write();

    g_v2_inc_pcm_staterr.SetName("g_v2_inc_pcm_staterr");
    g_v2_inc_pcm_staterr.SetTitle("inclusive photon v2 PCM, statistical uncertainties");
    g_v2_inc_pcm_staterr.Write();

    g_v2_inc_phos_toterr.SetName("g_v2_inc_phos_toterr");
    g_v2_inc_phos_toterr.SetTitle("inclusive photon v2 PHOS, total uncertainties");
    g_v2_inc_phos_toterr.Write();

    g_v2_inc_phos_syserr.SetName("g_v2_inc_phos_syserr");
    g_v2_inc_phos_syserr.SetTitle("inclusive photon v2 PHOS, systematic uncertainties");
    g_v2_inc_phos_syserr.Write();

    g_v2_inc_phos_staterr.SetName("g_v2_inc_phos_staterr");
    g_v2_inc_phos_staterr.SetTitle("inclusive photon v2 PHOS, statistical uncertainties");
    g_v2_inc_phos_staterr.Write();

    g_v2_inc_comb_toterr.SetName("g_v2_inc_comb_toterr");
    g_v2_inc_comb_toterr.SetTitle("inclusive photon v2 PCM/PHOS combined, total uncertainties");
    g_v2_inc_comb_toterr.Write();

    g_v2_inc_comb_syserr.SetName("g_v2_inc_comb_syserr");
    g_v2_inc_comb_syserr.SetTitle("inclusive photon v2 PCM/PHOS combined, systematic uncertainties");
    g_v2_inc_comb_syserr.Write();

    g_v2_inc_comb_staterr.SetName("g_v2_inc_comb_staterr");
    g_v2_inc_comb_staterr.SetTitle("inclusive photon v2 PCM/PHOS combined, statistical uncertainties");
    g_v2_inc_comb_staterr.Write();

    // decay photon v2
    g_v2_dec_comb.SetName("g_v2_dec_comb");
    g_v2_dec_comb.SetTitle("decay photon v2");
    g_v2_dec_comb.Write();

    //
    // histogram with posteriori v2dir distributions (PCM, PHOS, combined)
    //
    for (Int_t i = 0; i < nPtBins; i++) {
        h_v2_dir_pcm_toterr[i]->Write();
        h_v2_dir_phos_toterr[i]->Write();
        h_v2_dir_pcm_Rpcm_toterr[i]->Write();
        h_v2_dir_phos_Rphos_toterr[i]->Write();
        h_v2_dir_comb_toterr[i]->Write();
    }

    // write also Rgamma
    TGraphErrors g_Rgamma_toterr(nPtBins);
    for (Int_t i=0; i<nPtBins; ++i) {
        g_Rgamma_toterr.SetPoint(i, pt(i), Rgamma_meas_vec_comb(i));
        g_Rgamma_toterr.SetPointError(i, 0, TMath::Sqrt(cov_Rgamma_comb_toterr(i,i)));
    }

    g_Rgamma_toterr.SetName("g_Rgamma_toterr");
    g_Rgamma_toterr.SetTitle("g_Rgamma_toterr");
    g_Rgamma_toterr.Write();

    f_out.Close();

    //
    // Plot v2dir(v2incPCM,RgammaPCM), v2dir(v2incPHOS,RgammaPHOS), and combined v2dir vs pT
    //
    TCanvas *c2 = new TCanvas("c2", "c2", 700, 10, 700, 500);
    // TH2F frame("frame","direct photon v_{2}", 1, 0., 6, 1, -0.15, 0.35);
    TH2F frame("frame", "direct photon v_{2}", 1, 0., 6, 1, -0.05, 0.30);
    frame.SetXTitle("p_{T} (GeV/#it{c})");
    frame.SetYTitle("direct photon v_{2}");

    frame.DrawCopy();

    TLine ly0(0., 0., 6., 0.);
    ly0.DrawClone();

    g_v2_dir_comb_toterr.SetFillColorAlpha(kBlue - 9, 0.2);

    g_v2_dir_pcm_Rpcm_toterr.SetMarkerStyle(kFullCircle);
    g_v2_dir_phos_Rphos_toterr.SetMarkerStyle(kFullCircle);
    g_v2_dir_comb_toterr.SetMarkerStyle(kFullCircle);

    g_v2_dir_pcm_Rpcm_toterr.SetLineColor(kRed);
    g_v2_dir_phos_Rphos_toterr.SetLineColor(kBlue);
    g_v2_dir_comb_toterr.SetLineColor(kBlack);

    g_v2_dir_pcm_Rpcm_toterr.SetMarkerColor(kRed);
    g_v2_dir_phos_Rphos_toterr.SetMarkerColor(kBlue);
    g_v2_dir_comb_toterr.SetMarkerColor(kBlack);

    g_v2_dir_pcm_Rpcm_toterr.DrawClone("p");
    g_v2_dir_phos_Rphos_toterr.DrawClone("p");
    g_v2_dir_comb_toterr.DrawClone("pE5");

    TLatex l2;
    l2.DrawLatexNDC(0.2, 0.75, Form("%s%%", centr.Data()));

    TLegend leg2(0.55, 0.72, 0.87, 0.88, NULL, "brNDC");
    leg2.AddEntry(&g_v2_dir_pcm_Rpcm_toterr,
                  "#it{v}_{2,dir}^{PCM}(#it{v}_{2,inc}^{PCM},#it{v}_{2,dec},#it{R}_{#gamma}^{PCM})", "p");
    leg2.AddEntry(&g_v2_dir_phos_Rphos_toterr,
                  "#it{v}_{2,dir}^{PHOS}(#it{v}_{2,inc}^{PHOS},#it{v}_{2,dec},#it{R}_{#gamma}^{PHOS})", "p");
    leg2.AddEntry(&g_v2_dir_comb_toterr, "#it{v}_{2,dir}^{comb}", "p");
    leg2.SetBorderSize(1);
    leg2.SetTextSize(0.032);
    leg2.DrawClone();

    TString filename_c2 = "v2dir_pcm_Rpcm_phos_Rphos_comb_" + centr + ".pdf";
    c2->SaveAs(output_dir + filename_c2);

    //
    // Plot v2dir(v2incPCM,RgammaComb), v2dir(v2incPHOS,Rgammacomb), and combined v2dir vs pT
    //
    TCanvas *c2b = new TCanvas("c2b", "c2b", 700, 10, 700, 500);
    frame.DrawCopy();

    ly0.DrawClone();

    g_v2_dir_comb_toterr.SetFillColorAlpha(kBlue - 9, 0.2);

    g_v2_dir_pcm_toterr.SetMarkerStyle(kFullCircle);
    g_v2_dir_phos_toterr.SetMarkerStyle(kFullCircle);
    g_v2_dir_comb_toterr.SetMarkerStyle(kFullCircle);

    g_v2_dir_pcm_toterr.SetLineColor(kRed);
    g_v2_dir_phos_toterr.SetLineColor(kBlue);
    g_v2_dir_comb_toterr.SetLineColor(kBlack);

    g_v2_dir_pcm_toterr.SetMarkerColor(kRed);
    g_v2_dir_phos_toterr.SetMarkerColor(kBlue);
    g_v2_dir_comb_toterr.SetMarkerColor(kBlack);

    g_v2_dir_pcm_toterr.DrawClone("p");
    g_v2_dir_phos_toterr.DrawClone("p");
    g_v2_dir_comb_toterr.DrawClone("pE5");

    TLatex l2b;
    l2b.DrawLatexNDC(0.2, 0.75, Form("%s%%", centr.Data()));

    TLegend leg2b(0.55, 0.72, 0.87, 0.88, NULL, "brNDC");
    leg2b.AddEntry(&g_v2_dir_pcm_toterr,
                   "#it{v}_{2,dir}^{PCM}(#it{v}_{2,inc}^{PCM},#it{v}_{2,dec},#it{R}_{#gamma}^{comb})", "p");
    leg2b.AddEntry(&g_v2_dir_phos_toterr,
                   "#it{v}_{2,dir}^{PHOS}(#it{v}_{2,inc}^{PHOS},#it{v}_{2,dec},#it{R}_{#gamma}^{comb})", "p");
    leg2b.AddEntry(&g_v2_dir_comb_toterr, "#it{v}_{2,dir}^{comb}", "p");
    leg2b.SetBorderSize(1);
    leg2b.SetTextSize(0.032);
    leg2b.DrawClone();

    TString filename_c2b = "v2dir_pcm_phos_comb_" + centr + ".pdf";
    c2b->SaveAs(output_dir + filename_c2b);

    //
    // Plot only combined v2dir vs pT
    //
    TCanvas *c3 = new TCanvas("c3", "c3", 700, 500, 700, 500);
    frame.DrawCopy();

    ly0.DrawClone();

    g_v2_dir_comb_toterr.SetMarkerStyle(kFullCircle);
    g_v2_dir_comb_toterr.SetLineColor(kBlack);
    g_v2_dir_comb_toterr.SetMarkerColor(kBlack);
    g_v2_dir_comb_toterr.DrawClone("pE5");

    g_v2_dir_comb_staterr.SetMarkerStyle(kFullCircle);
    g_v2_dir_comb_staterr.SetLineColor(kBlack);
    g_v2_dir_comb_staterr.SetMarkerColor(kBlack);
    g_v2_dir_comb_staterr.DrawClone("p");

    TLatex l3;
    l3.DrawLatexNDC(0.2, 0.75, Form("%s%%", centr.Data()));

    TString filename_c3 = "v2dir_comb_" + centr + ".pdf";
    c3->SaveAs(output_dir + filename_c3);

    //
    // Plot combined v2dir vs pT along with Daniel's old PCM result
    //
    TCanvas *c3b = new TCanvas("c3b", "c3b", 700, 500, 700, 500);
    frame.DrawCopy();

    ly0.DrawClone();

    g_v2_dir_comb_toterr.SetMarkerStyle(kFullCircle);
    g_v2_dir_comb_toterr.SetLineColor(kBlack);
    g_v2_dir_comb_toterr.SetMarkerColor(kBlack);
    g_v2_dir_comb_toterr.DrawClone("pE5");

    g_v2_dir_comb_staterr.SetMarkerStyle(kFullCircle);
    g_v2_dir_comb_staterr.SetLineColor(kBlack);
    g_v2_dir_comb_staterr.SetMarkerColor(kBlack);
    g_v2_dir_comb_staterr.DrawClone("p");

    TFile f_pcm_daniel("results_daniel/PCMPi0v2.root");

    Int_t centr_id_daniel = 0;
    if (centr == "00-20")
        centr_id_daniel = 6;
    else if (centr == "20-40")
        centr_id_daniel = 7;

    TGraphAsymmErrors *v2dir_daniel_staterr =
        (TGraphAsymmErrors *)f_pcm_daniel.Get(Form("%dDirectPhotonv2", centr_id_daniel));
    TGraphAsymmErrors *v2dir_daniel_syserr =
        (TGraphAsymmErrors *)f_pcm_daniel.Get(Form("%dDirectPhotonv2Sys", centr_id_daniel));

    Int_t col = kRed;
    v2dir_daniel_staterr->SetMarkerStyle(kFullSquare);
    v2dir_daniel_staterr->SetMarkerColor(col);
    v2dir_daniel_staterr->SetLineColor(col);
    v2dir_daniel_syserr->SetMarkerStyle(kOpenCircle);
    v2dir_daniel_syserr->SetMarkerColor(col);
    v2dir_daniel_syserr->SetLineColor(col);
    v2dir_daniel_syserr->SetFillColor(col);
    v2dir_daniel_syserr->SetFillStyle(3002);

    v2dir_daniel_syserr->DrawClone("pE5");
    v2dir_daniel_staterr->DrawClone("p");

    f_pcm_daniel.Close();

    TLatex l3b;
    l3b.DrawLatexNDC(0.2, 0.75, Form("%s%%", centr.Data()));

    TLegend leg3b(0.55, 0.72, 0.87, 0.88, NULL, "brNDC");
    leg3b.AddEntry(&g_v2_dir_comb_toterr, "#it{v}_{2,dir}^{comb}", "p");
    leg3b.AddEntry(v2dir_daniel_staterr, "#it{v}_{2,dir}^{Daniel}", "p");
    leg3b.SetBorderSize(1);
    leg3b.SetTextSize(0.032);
    leg3b.DrawClone();

    TString filename_c3b = "v2dir_comb_cmp_daniel_" + centr + ".pdf";
    c3b->SaveAs(output_dir + filename_c3b);

    //
    // Plot only combined v2dir vs pT assuming pT-uncorrelated systematic errors (cross check)
    //
    TCanvas *c3c = new TCanvas("c3c", "c3c", 700, 500, 700, 500);
    frame.DrawCopy();

    ly0.DrawClone();

    g_v2_dir_comb_toterr_uncorr.SetFillColor(kOrange - 3);

    g_v2_dir_comb_toterr_uncorr.SetMarkerStyle(kFullCircle);
    g_v2_dir_comb_toterr_uncorr.SetLineColor(kBlack);
    g_v2_dir_comb_toterr_uncorr.SetMarkerColor(kBlack);
    g_v2_dir_comb_toterr_uncorr.DrawClone("pE5");

    g_v2_dir_comb_staterr_uncorr.SetMarkerStyle(kFullCircle);
    g_v2_dir_comb_staterr_uncorr.SetLineColor(kBlack);
    g_v2_dir_comb_staterr_uncorr.SetMarkerColor(kBlack);
    g_v2_dir_comb_staterr_uncorr.DrawClone("p");

    TLatex l3c;
    l3c.DrawLatexNDC(0.2, 0.75, Form("%s%%", centr.Data()));
    l3c.DrawLatexNDC(0.2, 0.15, "uncorrelated errors (cross check)");

    TString filename_c3c = "v2dir_comb_uncorr_errors_" + centr + ".pdf";
    c3c->SaveAs(output_dir + filename_c3c);

    //
    // Plot only combined v2dir vs pT + v2dir for given constant Rgamma (cross check)
    //
    TCanvas *c3d = new TCanvas("c3d", "c3d", 700, 500, 700, 500);
    frame.DrawCopy();

    ly0.DrawClone();

    g_v2_dir_comb_toterr.SetMarkerStyle(kFullCircle);
    g_v2_dir_comb_toterr.SetLineColor(kBlack);
    g_v2_dir_comb_toterr.SetMarkerColor(kBlack);
    g_v2_dir_comb_toterr.DrawClone("pE5");

    g_v2_dir_comb_staterr.SetMarkerStyle(kFullCircle);
    g_v2_dir_comb_staterr.SetLineColor(kBlack);
    g_v2_dir_comb_staterr.SetMarkerColor(kBlack);
    g_v2_dir_comb_staterr.DrawClone("p");

    g_v2_dir_comb_toterr_xcheck1.SetFillColorAlpha(kOrange - 3, 0.2);

    g_v2_dir_comb_toterr_xcheck1.SetMarkerStyle(kOpenCircle);
    g_v2_dir_comb_toterr_xcheck1.SetLineColor(kBlack);
    g_v2_dir_comb_toterr_xcheck1.SetMarkerColor(kBlack);
    g_v2_dir_comb_toterr_xcheck1.DrawClone("pE5");

    TLatex l3d;
    l3d.DrawLatexNDC(0.2, 0.75, Form("%s%%", centr.Data()));
    l3d.DrawLatexNDC(0.2, 0.15, Form("constant R_{gamma} = %4.2f (cross check)", Rgamma_xcheck1_val));

    TString filename_c3d = "v2dir_comb_xcheck1_" + centr + ".pdf";
    c3d->SaveAs(output_dir + filename_c3d);

    //
    // Plot only combined v2dir vs pT + v2dir for Rgamma scaled by purity ratio data-driven/MC-based (cross check)
    //
    TCanvas *c3e = new TCanvas("c3e", "c3e", 700, 500, 700, 500);
    frame.DrawCopy();

    ly0.DrawClone();

    g_v2_dir_comb_toterr.SetMarkerStyle(kFullCircle);
    g_v2_dir_comb_toterr.SetLineColor(kBlack);
    g_v2_dir_comb_toterr.SetMarkerColor(kBlack);
    g_v2_dir_comb_toterr.DrawClone("pE5");

    g_v2_dir_comb_staterr.SetMarkerStyle(kFullCircle);
    g_v2_dir_comb_staterr.SetLineColor(kBlack);
    g_v2_dir_comb_staterr.SetMarkerColor(kBlack);
    g_v2_dir_comb_staterr.DrawClone("p");

    g_v2_dir_comb_toterr_xcheck2.SetFillColorAlpha(kOrange - 3, 0.2);

    g_v2_dir_comb_toterr_xcheck2.SetMarkerStyle(kOpenCircle);
    g_v2_dir_comb_toterr_xcheck2.SetLineColor(kBlack);
    g_v2_dir_comb_toterr_xcheck2.SetMarkerColor(kBlack);
    g_v2_dir_comb_toterr_xcheck2.DrawClone("pE5");

    TLatex l3e;
    l3e.DrawLatexNDC(0.2, 0.75, Form("%s%%", centr.Data()));
    l3e.DrawLatexNDC(0.2, 0.17, "R_{#gamma} #times p_{data-driven}/p_{MC} [p = purity] (cross check)");

    TString filename_c3e = "v2dir_comb_xcheck2_" + centr + ".pdf";
    c3e->SaveAs(output_dir + filename_c3e);

    //
    // Plot default v2,dir along with weighted average of
    // v2,dir,PCM(RgammaPCM) and v2,dir,PHOS(RgammaPHOS)
    //
    TCanvas *c3f = new TCanvas("c3f", "c3f", 700, 500, 700, 500);
    frame.DrawCopy();

    ly0.DrawClone();

    g_v2_dir_comb_toterr.SetMarkerStyle(kFullCircle);
    g_v2_dir_comb_toterr.SetLineColor(kBlack);
    g_v2_dir_comb_toterr.SetMarkerColor(kBlack);
    g_v2_dir_comb_toterr.DrawClone("pE5");

    g_v2_dir_comb_staterr.SetMarkerStyle(kFullCircle);
    g_v2_dir_comb_staterr.SetLineColor(kBlack);
    g_v2_dir_comb_staterr.SetMarkerColor(kBlack);
    // g_v2_dir_comb_staterr.DrawClone("p");

    g_v2_dir_comb_ave_xcheck.SetFillColorAlpha(kOrange - 3, 0.2);

    g_v2_dir_comb_ave_xcheck.SetMarkerStyle(kOpenCircle);
    g_v2_dir_comb_ave_xcheck.SetLineColor(kBlack);
    g_v2_dir_comb_ave_xcheck.SetMarkerColor(kBlack);
    g_v2_dir_comb_ave_xcheck.DrawClone("pE5");

    TLatex l3f;
    l3f.DrawLatexNDC(0.2, 0.75, Form("%s%%", centr.Data()));
    l3f.SetTextSize(0.04);
    l3f.DrawLatexNDC(0.2, 0.17, "Weighted average of v_{2,dir}^{PCM} and v_{2,dir}^{PHOS} (cross check)");

    TString filename_c3f = "v2dir_comb_ave_xcheck_" + centr + ".pdf";
    c3f->SaveAs(output_dir + filename_c3f);

    //
    // Plot default v2,dir along with central value directly calculated from the
    // measured values, i.e., v2dir = (Rgamma,meas * v2,inc,meas - v2,dec) / (Rgamma,meas - 1)
    //
    TCanvas *c3g = new TCanvas("c3g", "c3g", 700, 500, 700, 500);
    frame.DrawCopy();

    ly0.DrawClone();

    g_v2_dir_comb_toterr.SetMarkerStyle(kFullCircle);
    g_v2_dir_comb_toterr.SetLineColor(kBlack);
    g_v2_dir_comb_toterr.SetMarkerColor(kBlack);
    g_v2_dir_comb_toterr.DrawClone("pE5");

    g_v2_dir_comb_staterr.SetMarkerStyle(kFullCircle);
    g_v2_dir_comb_staterr.SetLineColor(kBlack);
    g_v2_dir_comb_staterr.SetMarkerColor(kBlack);
    // g_v2_dir_comb_staterr.DrawClone("p");

    g_v2_dir_comb_classic_xcheck.SetFillColorAlpha(kOrange - 3, 0.2);

    g_v2_dir_comb_classic_xcheck.SetMarkerStyle(kOpenCircle);
    g_v2_dir_comb_classic_xcheck.SetLineColor(kOrange - 3);
    g_v2_dir_comb_classic_xcheck.SetMarkerColor(kOrange - 3);
    g_v2_dir_comb_classic_xcheck.DrawClone("pE5");

    TLatex l3g;
    l3g.DrawLatexNDC(0.2, 0.75, Form("%s%%", centr.Data()));
    l3g.SetTextSize(0.03);
    l3g.DrawLatexNDC(0.2, 0.17, "Directly calculated central values + Gaussian error prop. (cross check)");

    TString filename_c3g = "v2dir_comb_classic_xcheck_" + centr + ".pdf";
    c3g->SaveAs(output_dir + filename_c3g);

    //
    // Plot only combined v2dir vs pT assuming pT-uncorrelated systematic errors (cross check)
    // along with central value directly calculated from the
    // measured values, i.e., v2dir = (Rgamma,meas * v2,inc,meas - v2,dec) / (Rgamma,meas - 1)
    //
    TCanvas *c3h = new TCanvas("c3c", "c3c", 700, 500, 700, 500);
    frame.DrawCopy();

    ly0.DrawClone();

    g_v2_dir_comb_toterr_uncorr.SetFillColor(kOrange - 3);

    g_v2_dir_comb_toterr_uncorr.SetMarkerStyle(kFullCircle);
    g_v2_dir_comb_toterr_uncorr.SetLineColor(kBlack);
    g_v2_dir_comb_toterr_uncorr.SetMarkerColor(kBlack);
    g_v2_dir_comb_toterr_uncorr.DrawClone("pE5");

    g_v2_dir_comb_staterr_uncorr.SetMarkerStyle(kFullCircle);
    g_v2_dir_comb_staterr_uncorr.SetLineColor(kBlack);
    g_v2_dir_comb_staterr_uncorr.SetMarkerColor(kBlack);
    g_v2_dir_comb_staterr_uncorr.DrawClone("p");

    g_v2_dir_comb_classic_xcheck.SetFillColorAlpha(kBlue - 3, 0.2);
    g_v2_dir_comb_classic_xcheck.SetMarkerStyle(kOpenCircle);
    g_v2_dir_comb_classic_xcheck.SetLineColor(kBlue - 3);
    g_v2_dir_comb_classic_xcheck.SetMarkerColor(kBlue - 3);
    g_v2_dir_comb_classic_xcheck.DrawClone("pE5");

    TLatex l3h;
    l3h.SetTextSize(0.03);
    l3h.DrawLatexNDC(0.2, 0.75, Form("%s%%", centr.Data()));
    l3h.DrawLatexNDC(0.2, 0.15, "uncorrelated errors (cross check) and classic v2dir (cross check)");

    TString filename_c3h = "v2dir_comb_uncorr_errors_classic_xcheck_" + centr + ".pdf";
    c3h->SaveAs(output_dir + filename_c3h);

    //
    // Plot inclusive and decay photon v2
    //
    TCanvas *c4 = new TCanvas("c4", "c4", 900, 10, 700, 500);
    frame.SetXTitle("p_{T} (GeV/#it{c})");
    frame.SetYTitle("inclusive photon v_{2}, decay photon v_{2}");

    frame.DrawCopy();

    g_v2_dec_comb.SetMarkerStyle(kOpenCircle);
    g_v2_dec_comb.SetLineColor(kGray);
    g_v2_dec_comb.SetMarkerColor(kGray);
    g_v2_dec_comb.SetFillColor(kGray);
    g_v2_dec_comb.DrawClone("3l");

    g_v2_inc_pcm_toterr.SetMarkerStyle(kFullCircle);
    g_v2_inc_phos_toterr.SetMarkerStyle(kFullCircle);
    g_v2_inc_comb_toterr.SetMarkerStyle(kFullCircle);

    g_v2_inc_pcm_toterr.SetLineColor(kRed);
    g_v2_inc_phos_toterr.SetLineColor(kBlue);
    g_v2_inc_comb_toterr.SetLineColor(kBlack);

    g_v2_inc_pcm_toterr.SetMarkerColor(kRed);
    g_v2_inc_phos_toterr.SetMarkerColor(kBlue);
    g_v2_inc_comb_toterr.SetMarkerColor(kBlack);

    g_v2_inc_pcm_toterr.DrawClone("p");
    g_v2_inc_phos_toterr.DrawClone("p");
    g_v2_inc_comb_toterr.DrawClone("p");

    g_v2_dec_comb.DrawClone("p");

    TLatex l4;
    l4.DrawLatexNDC(0.2, 0.75, Form("%s%%", centr.Data()));

    TLegend leg4(0.65, 0.68, 0.87, 0.88, NULL, "brNDC");
    leg4.AddEntry(&g_v2_inc_pcm_toterr, "#it{v}_{2,inc}^{PCM}", "p");
    leg4.AddEntry(&g_v2_inc_phos_toterr, "#it{v}_{2,inc}^{PHOS}", "p");
    leg4.AddEntry(&g_v2_inc_comb_toterr, "#it{v}_{2,inc}^{combined}", "p");
    leg4.AddEntry(&g_v2_dec_comb, "#it{v}_{2,dec}", "l");

    leg4.SetBorderSize(1);
    leg4.SetTextSize(0.032);
    leg4.DrawClone();

    TString filename_c4 = "v2inc_pcm_phos_comb_" + centr + ".pdf";
    c4->SaveAs(output_dir + filename_c4);

    //
    // Plot current v2,inc,PCM along with Daniel's result
    //
    TCanvas *c5 = new TCanvas("c5", "c5", 900, 10, 700, 500);
    frame.SetXTitle("p_{T} (GeV/#it{c})");
    frame.SetYTitle("inclusive photon v_{2}");
    frame.DrawCopy();

    // Daniel's v2,inc (0-20%)
    TGraphErrors *gre0020 = new TGraphErrors(23);
    gre0020->SetName("grDanielv2inc0020");
    gre0020->SetTitle("");
    gre0020->SetFillColor(1);
    gre0020->SetFillStyle(0);
    gre0020->SetLineWidth(2);
    gre0020->SetMarkerStyle(24);
    gre0020->SetMarkerSize(1.2);
    gre0020->SetPoint(0, 0.25, 0.02653868);
    gre0020->SetPointError(0, 0.05, 0.0004103671);
    gre0020->SetPoint(1, 0.35, 0.03666553);
    gre0020->SetPointError(1, 0.05, 0.0005929124);
    gre0020->SetPoint(2, 0.45, 0.04550414);
    gre0020->SetPointError(2, 0.05, 0.0004489667);
    gre0020->SetPoint(3, 0.55, 0.05248478);
    gre0020->SetPointError(3, 0.05, 0.0004132783);
    gre0020->SetPoint(4, 0.65, 0.05761384);
    gre0020->SetPointError(4, 0.05, 0.0002697986);
    gre0020->SetPoint(5, 0.75, 0.06367664);
    gre0020->SetPointError(5, 0.05, 0.0002163565);
    gre0020->SetPoint(6, 0.85, 0.06925035);
    gre0020->SetPointError(6, 0.05, 0.0003906399);
    gre0020->SetPoint(7, 1, 0.07546688);
    gre0020->SetPointError(7, 0.1, 0.0004597703);
    gre0020->SetPoint(8, 1.2, 0.08402089);
    gre0020->SetPointError(8, 0.1, 0.0004792006);
    gre0020->SetPoint(9, 1.4, 0.09208977);
    gre0020->SetPointError(9, 0.1, 0.0005386815);
    gre0020->SetPoint(10, 1.6, 0.09751575);
    gre0020->SetPointError(10, 0.1, 0.0006099597);
    gre0020->SetPoint(11, 1.8, 0.1044671);
    gre0020->SetPointError(11, 0.1, 0.000524053);
    gre0020->SetPoint(12, 2, 0.1053657);
    gre0020->SetPointError(12, 0.1, 0.001143335);
    gre0020->SetPoint(13, 2.2, 0.1147294);
    gre0020->SetPointError(13, 0.1, 0.000629642);
    gre0020->SetPoint(14, 2.4, 0.1080451);
    gre0020->SetPointError(14, 0.1, 0.001381279);
    gre0020->SetPoint(15, 2.6, 0.09923275);
    gre0020->SetPointError(15, 0.1, 0.001402856);
    gre0020->SetPoint(16, 2.85, 0.1067197);
    gre0020->SetPointError(16, 0.15, 0.001420079);
    gre0020->SetPoint(17, 3.15, 0.104392);
    gre0020->SetPointError(17, 0.15, 0.001612167);
    gre0020->SetPoint(18, 3.5, 0.07886321);
    gre0020->SetPointError(18, 0.2, 0.003232181);
    gre0020->SetPoint(19, 3.9, 0.09279333);
    gre0020->SetPointError(19, 0.2, 0.002346974);
    gre0020->SetPoint(20, 4.35, 0.06173879);
    gre0020->SetPointError(20, 0.25, 0.003790896);
    gre0020->SetPoint(21, 5, 0.06296306);
    gre0020->SetPointError(21, 0.4, 0.005187371);
    gre0020->SetPoint(22, 5.8, 0.03971473);
    gre0020->SetPointError(22, 0.4, 0.006876481);

    // Daniel's v2,inc (20-40%)
    TGraphErrors *gre2040 = new TGraphErrors(23);
    gre2040->SetName("grDanielv2inc2040");
    gre2040->SetTitle("");
    gre2040->SetFillColor(1);
    gre2040->SetFillStyle(0);
    gre2040->SetLineWidth(2);
    gre2040->SetMarkerStyle(24);
    gre2040->SetMarkerSize(1.2);
    gre2040->SetPoint(0, 0.25, 0.05754349);
    gre2040->SetPointError(0, 0.05, 0.0004665395);
    gre2040->SetPoint(1, 0.35, 0.07447413);
    gre2040->SetPointError(1, 0.05, 0.0007736117);
    gre2040->SetPoint(2, 0.45, 0.09129447);
    gre2040->SetPointError(2, 0.05, 0.0004530184);
    gre2040->SetPoint(3, 0.55, 0.104128);
    gre2040->SetPointError(3, 0.05, 0.0005179349);
    gre2040->SetPoint(4, 0.65, 0.1163502);
    gre2040->SetPointError(4, 0.05, 0.0004933628);
    gre2040->SetPoint(5, 0.75, 0.1269654);
    gre2040->SetPointError(5, 0.05, 0.0004371996);
    gre2040->SetPoint(6, 0.85, 0.1370555);
    gre2040->SetPointError(6, 0.05, 0.0006342705);
    gre2040->SetPoint(7, 1, 0.149815);
    gre2040->SetPointError(7, 0.1, 0.0007434604);
    gre2040->SetPoint(8, 1.2, 0.165099);
    gre2040->SetPointError(8, 0.1, 0.0009312401);
    gre2040->SetPoint(9, 1.4, 0.1774414);
    gre2040->SetPointError(9, 0.1, 0.001089122);
    gre2040->SetPoint(10, 1.6, 0.1881529);
    gre2040->SetPointError(10, 0.1, 0.001280767);
    gre2040->SetPoint(11, 1.8, 0.1952062);
    gre2040->SetPointError(11, 0.1, 0.001302265);
    gre2040->SetPoint(12, 2, 0.1956598);
    gre2040->SetPointError(12, 0.1, 0.001264462);
    gre2040->SetPoint(13, 2.2, 0.1948985);
    gre2040->SetPointError(13, 0.1, 0.001483711);
    gre2040->SetPoint(14, 2.4, 0.1946173);
    gre2040->SetPointError(14, 0.1, 0.001965538);
    gre2040->SetPoint(15, 2.6, 0.1950112);
    gre2040->SetPointError(15, 0.1, 0.001783452);
    gre2040->SetPoint(16, 2.85, 0.1866001);
    gre2040->SetPointError(16, 0.15, 0.001826611);
    gre2040->SetPoint(17, 3.15, 0.1728413);
    gre2040->SetPointError(17, 0.15, 0.002727349);
    gre2040->SetPoint(18, 3.5, 0.1501855);
    gre2040->SetPointError(18, 0.2, 0.001926224);
    gre2040->SetPoint(19, 3.9, 0.1445375);
    gre2040->SetPointError(19, 0.2, 0.004738899);
    gre2040->SetPoint(20, 4.35, 0.1158196);
    gre2040->SetPointError(20, 0.25, 0.002635487);
    gre2040->SetPoint(21, 5, 0.1144715);
    gre2040->SetPointError(21, 0.4, 0.004965387);
    gre2040->SetPoint(22, 5.8, 0.09643851);
    gre2040->SetPointError(22, 0.4, 0.005143298);

    g_v2_inc_pcm_toterr.DrawClone("pE5");

    TGraphErrors *greDaniel = 0;
    if (centr == "00-20")
        greDaniel = gre0020;
    else if (centr == "20-40")
        greDaniel = gre2040;
    greDaniel->Draw("p");

    TLatex l5;
    l5.DrawLatexNDC(0.2, 0.75, Form("%s%%", centr.Data()));

    TLegend leg5(0.65, 0.68, 0.87, 0.88, NULL, "brNDC");
    leg5.AddEntry(&g_v2_inc_pcm_toterr, "#it{v}_{2,inc}^{PCM} (Mike)", "p");
    leg5.AddEntry(greDaniel, "#it{v}_{2,inc}^{PCM} (Daniel)", "p");

    leg5.SetBorderSize(1);
    leg5.SetTextSize(0.032);
    leg5.DrawClone();

    TString filename_c5 = "v2inc_pcm_new_old_" + centr + ".pdf";
    c5->SaveAs(output_dir + filename_c5);

    //
    // Plot current v2,inc,PCM along with Daniel's result
    //
    TCanvas *c5b = new TCanvas("c5b", "c5b", 900, 10, 700, 500);
    TH2F frame2("frame2", "frame2", 1, 0., 6., 1, 0.8, 1.2);
    frame2.SetXTitle("p_{T} (GeV/#it{c})");
    frame2.SetYTitle("v_{2}^{inc}(Daniel) / v_{2}^{inc}(Mike)");
    frame2.DrawCopy();

    TGraphErrors *g_ratio_mike_daniel = new TGraphErrors(nPtBins);
    for (Int_t i = 0; i < nPtBins; i++) {
        Double_t v2_Daniel = greDaniel->GetY()[i + 7];
        Double_t v2_Daniel_err = greDaniel->GetEY()[i + 7];
        Double_t v2_Daniel_relerr = v2_Daniel_err / v2_Daniel;
        Double_t v2_Mike = v2_inc_meas_values_pcm(i);
        Double_t v2_Mike_err = TMath::Sqrt(cov_v2_inc_toterr_pcm(i, i));
        Double_t v2_Mike_relerr = v2_Mike_err / v2_Mike;
        Double_t ratio = v2_Daniel / v2_Mike;
        Double_t ratio_err = ratio * TMath::Sqrt(v2_Mike_relerr * v2_Mike_relerr + v2_Daniel_relerr * v2_Daniel_relerr);
        g_ratio_mike_daniel->SetPoint(i, pt(i), ratio);
        g_ratio_mike_daniel->SetPointError(i, 0, ratio_err);
    }

    g_ratio_mike_daniel->SetMarkerStyle(20);
    g_ratio_mike_daniel->Draw("p");

    TLatex l5b;
    l5b.DrawLatexNDC(0.2, 0.75, Form("%s%%", centr.Data()));

    TString filename_c5b = "v2inc_ratio_daniel_" + centr + ".pdf";
    c5b->SaveAs(output_dir + filename_c5b);


    //
    // p-values w.r.t. to null hypotheses
    //

    // Null hyposthis H0: v2dir = 0
    TVectorD H0_v2dir_zero(nPtBins); // null hypothesis, contains only zeros
    TVectorD v2_dir_comb_toterr(nPtBins);
    for (Int_t i = 0; i < nPtBins; ++i)
        v2_dir_comb_toterr(i) = g_v2_dir_comb_toterr.GetY()[i];
    Double_t chi2val = chi2(v2_dir_comb_toterr, H0_v2dir_zero, cov_v2_dir_comb_toterr);
    Double_t pvalue = TMath::Prob(chi2val, nPtBins);
    cout << endl << "p-values:" << endl;
    cout << centr << "%: ";
    cout << "chi2(v2dir=0) = " << chi2val;
    cout << ", p-value = " << TMath::Prob(chi2val, nPtBins);
    cout << ", n_sigma = " << p_value_to_n_sigma(pvalue) << endl;

    // Null hypothesis H0; v2dir = v2decay
    TVectorD H0_v2dir_equal_to_v2dec(nPtBins); // null hypothesis, v2dir = v2dec
    for (Int_t i = 0; i < nPtBins; ++i)
        H0_v2dir_equal_to_v2dec(i) = v2_dec_meas_values(i);
    Double_t chi2val_2 = chi2(v2_dir_comb_toterr, H0_v2dir_equal_to_v2dec, cov_v2_dir_comb_toterr);
    Double_t pvalue_2 = TMath::Prob(chi2val_2, nPtBins);
    cout << centr << "%: ";
    cout << "chi2(v2dir=v2dec) = " << chi2val_2;
    cout << ", p-value = " << TMath::Prob(chi2val_2, nPtBins);
    cout << ", n_sigma = " << p_value_to_n_sigma(pvalue_2) << endl << endl;

    //
    // Finally calculate likelihoods for certain v2dir hypothesis:
    //
    // hypothesis for vndir,true (all values zero in this case)
    //
    // CAUTION: this part is still under development
    //

    TVectorD vnDirTrueHyp(nPtBins);
    for (Int_t i = 0; i < nPtBins; i++)
        vnDirTrueHyp(i) = 0.0;

    // likelihood 1: hypothesis v2,dir = v2,decay
    Double_t L1 = 0;
    // Double_t nlh = 1000000;
    Double_t nlh = 100000;

    likelihood(nPtBins, nlh, v2_dec_meas_values, Rgamma_meas_vec_comb, cov_Rgamma_comb_toterr, v2_inc_meas_values_comb,
               cov_v2_inc_toterr_comb, v2_dec_meas_values, cov_v2_dec_toterr, L1);

    // likelihood 2: v2,dir = 0
    Double_t L2 = 0;
    likelihood(nPtBins, nlh, vnDirTrueHyp, Rgamma_meas_vec_comb, cov_Rgamma_comb_toterr, v2_inc_meas_values_comb,
               cov_v2_inc_toterr_comb, v2_dec_meas_values, cov_v2_dec_toterr, L2);

    cout << "likelihood 1 (v2,dir = v2,decay) = " << L1 << endl;
    cout << "likelihood 2 (v2,dir = 0) = " << L2 << endl;
    cout << "Bayes factor L1 / L2 = " << L1 / L2 << endl;

    //
    // Question: How consistent is v2,dir from PCM with the hypothesis that v2,dir,true = v2,cocktail
    //
    // TVectorD vnDirTrueHypPCM(nPtBins);
    // for (Int_t i=0; i < nPtBins; i++) vnDirTrueHypPCM(i) = g_v2_dir_pcm_Rpcm_toterr.GetY()[i];
    // Double_t nlh_PCM = 10000;
    // Double_t L_PCM_1 = 0;
    // likelihood(nPtBins, nlh_PCM, vnDirTrueHypPCM, Rgamma_meas_vec_pcm, cov_Rgamma_toterr_pcm,
    // v2_inc_meas_values_pcm, cov_v2_inc_toterr_pcm, v2_dec_meas_values, cov_v2_dec_toterr, L_PCM_1);
    // likelihood 2
    // Double_t L_PCM_2 = 0;
    // likelihood(nPtBins, nlh_PCM, v2_dec_meas_values, Rgamma_meas_vec_pcm, cov_Rgamma_toterr_pcm,
    //   v2_inc_meas_values_pcm, cov_v2_inc_toterr_pcm, v2_dec_meas_values, cov_v2_dec_toterr, L_PCM_2);
    // cout << "L_PCM_1 = " << L_PCM_1 << endl;
    // cout << "L_PCM_2 = " << L_PCM_2 << endl;
    // cout << "Bayes factor L_PCM_1 / L_PCM_2 = " << L_PCM_1 / L_PCM_2 << endl;
}

//
// --- functions used in main part ---
//

void v2_dir_central_value_and_error_from_v2dir_distribution(TH1 *h_v2dir_distr, Double_t &v2dir_central_value,
                                                            Double_t &v2dir_err_up, Double_t &v2dir_err_low) {

    // determine median and errors
    const Int_t n_quant = 3;
    Double_t prob_integrals[n_quant] = {0.15865, 0.5, 0.84135};
    Double_t quantiles[n_quant];

    h_v2dir_distr->GetQuantiles(n_quant, quantiles, prob_integrals);
    v2dir_central_value = quantiles[1];
    v2dir_err_up = quantiles[2] - v2dir_central_value;
    v2dir_err_low = v2dir_central_value - quantiles[0];
}

void fill_v2dir_graph(TH1 **hp, const TVectorD &pt, const Double_t &pt_offset, const Double_t &err_x,
                      TGraphAsymmErrors &g) {
    for (Int_t i = 0; i < nPtBins; i++) {

        Double_t central_value = 0;
        Double_t err_high = 0;
        Double_t err_low = 0;
        v2_dir_central_value_and_error_from_v2dir_distribution(hp[i], central_value, err_high, err_low);

        g.SetPoint(i, pt(i) + pt_offset, central_value); // central value same as for g_v2_dir_comb_toterr
        g.SetPointEYhigh(i, err_high);
        g.SetPointEYlow(i, err_low);

        // for visual representation of errors as boxes
        g.SetPointEXhigh(i, err_x);
        g.SetPointEXlow(i, err_x);
    }
}

void weighted_average(const Int_t nMeas, const TVectorD **mu, const TMatrixDSym **cov, TVectorD *muAve,
                      TMatrixDSym *covAve) {

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
    for (Int_t iMeas = 0; iMeas < nMeas; iMeas++) {
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

void weighted_average(const TVectorD &mu1, const TVectorD &mu2, const TMatrixDSym &cov1, const TMatrixDSym &cov2,
                      TVectorD &muAve, TMatrixDSym &covAve) {

    // weighted average of two independent sets of measurements
    // mu1: mean values data set 1
    // mu2: mean values data set 2
    // cov1: covariance matrix data set 1
    // cov2: covariance matrix data set 2
    //
    // result:
    // muAve: weighted average (TVectorD)
    // covAve: covariance Matrix of the weighted average

    const TVectorD *mu[2];
    mu[0] = &mu1;
    mu[1] = &mu2;

    const TMatrixDSym *cov[2];
    cov[0] = &cov1;
    cov[1] = &cov2;

    weighted_average(2, mu, cov, &muAve, &covAve);
}

void likelihood(const Int_t nPtBins, const Int_t nSamples, const TVectorD &vnDirTrueVec, const TVectorD &RgamMeasVec,
                const TMatrixDSym &covRgam, const TVectorD &vnIncMeasVec, const TMatrixDSym &covVnInc,
                const TVectorD &vnDecMeasVec, const TMatrixDSym &covVnDec, Double_t &Lfinal) {

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
    RooRealVar *vnIncVar[nPtBins];
    RooRealVar *vnIncMeasVar;

    // RooFit variables for inclusive photon vn
    for (Int_t i = 0; i < nPtBins; i++) {
        char *nameMeas = Form("vnIncMeas%d", i);
        vnIncMeasVar = new RooRealVar(nameMeas, nameMeas, vnIncMeasVec[i], -0.5, 0.5);
        vnIncMeas.add(*vnIncMeasVar);

        char *nameTrue = Form("vnInc%d", i);
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
    RooRealVar *RgamVar;
    RooRealVar *RgamMeasVar;

    // RooFit variables Rgamma
    for (Int_t i = 0; i < nPtBins; i++) {
        char *nameMeas = Form("RgamMeas%d", i);
        RgamMeasVar = new RooRealVar(nameMeas, nameMeas, RgamMeasVec[i], 0., 10.);
        RgamMeas.add(*RgamMeasVar);

        char *nameTrue = Form("Rgam%d", i);
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
    RooRealVar *vnDecVar;
    RooRealVar *vnDecMeasVar;

    // RooFit variables for decay photon vn
    for (Int_t i = 0; i < nPtBins; i++) {
        char *nameMeas = Form("vnDecMeas%d", i);
        vnDecMeasVar = new RooRealVar(nameMeas, nameMeas, vnDecMeasVec[i], -0.5, 0.5);
        vnDecMeas.add(*vnDecMeasVar);

        char *nameTrue = Form("vnDec%d", i);
        vnDecVar = new RooRealVar(nameTrue, nameTrue, 0., -0.5, 0.5);
        vnDec.add(*vnDecVar);
    }

    // now make the multivariate Gaussian
    RooMultiVarGaussian mvgVnDec("mvgVnDec", "mvgVnDec", vnDec, vnDecMeas, covVnDec);

    //
    // sample the distributions
    //
    RooDataSet *dataVnDec = mvgVnDec.generate(vnDec, nSamples);
    RooDataSet *dataRgam = mvgRgam.generate(Rgam, nSamples);

    //
    // create Likelihood variable
    //
    RooArgList L;
    RooRealVar *LVar = new RooRealVar("L", "L", 0., 1.);
    L.add(*LVar);

    RooDataSet dataL("dataL", "dataL", L);

    // mvgRgam.Print();
    for (Int_t k = 0; k < dataVnDec->numEntries(); k++) {

        const RooArgSet *rowVnDec = dataVnDec->get(k);
        const RooArgSet *rowRgam = dataRgam->get(k);

        vnDec = *rowVnDec;
        Rgam = *rowRgam;

        // cout << mvgRgam.getVal() << endl;
        for (Int_t i = 0; i < nPtBins; i++) {

            RooRealVar *vnDecRow = (RooRealVar *)rowVnDec->find(vnDec.at(i)->GetName());
            RooRealVar *RgamRow = (RooRealVar *)rowRgam->find(Rgam.at(i)->GetName());

            Double_t RgamVal = RgamRow->getVal();
            Double_t vnDecVal = vnDecRow->getVal();
            Double_t vnIncVal = (vnDirTrueVec(i) * (RgamVal - 1.) + vnDecVal) / RgamVal;

            // cout << "vnIncVal = " << vnIncVal << ", RgamVal = " << RgamVal << ", vnDecVal = " << vnDecVal << endl;

            // ((RooRealVar*) vnInc.at(i))->setVal(vnIncVal);
            *vnIncVar[i] = vnIncVal;
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

void v2dir(const Int_t nPtBins, const Int_t nSamples, const TVectorD &RgamMeasVec, const TMatrixDSym &covRgam,
           const TVectorD &vnIncMeasVec, const TMatrixDSym &covVnInc, const TVectorD &vnDecMeasVec,
           const TMatrixDSym &covVnDec, TH1 **hVnDir, TString prefix, TMatrixDSym &covarianceMatrix,
           TMatrixDSym &correlationMatrix) {

    //
    // fill histograms with v2,dir true for all pT bins by sampling multi-variate Gaussians
    //

    //
    // Rgamma
    //
    RooArgList Rgam;
    RooArgList RgamMeas;
    RooRealVar *RgamVar;
    RooRealVar *RgamMeasVar;

    // RooFit variables Rgamma
    for (Int_t i = 0; i < nPtBins; i++) {
        char *nameMeas = Form("RgamMeas%d", i);
        RgamMeasVar = new RooRealVar(nameMeas, nameMeas, RgamMeasVec[i], 0., 10.);
        RgamMeas.add(*RgamMeasVar);

        char *nameTrue = Form("Rgam%d", i);
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
    RooRealVar *vnIncVar;
    RooRealVar *vnIncMeasVar;

    // RooFit variables for inclusive photon vn
    for (Int_t i = 0; i < nPtBins; i++) {
        char *nameMeas = Form("vnIncMeas%d", i);
        vnIncMeasVar = new RooRealVar(nameMeas, nameMeas, vnIncMeasVec[i], -0.5, 0.5);
        vnIncMeas.add(*vnIncMeasVar);

        char *nameTrue = Form("vnInc%d", i);
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
    RooRealVar *vnDecVar;
    RooRealVar *vnDecMeasVar;

    // RooFit variables for decay photon vn
    for (Int_t i = 0; i < nPtBins; i++) {
        char *nameMeas = Form("vnDecMeas%d", i);
        vnDecMeasVar = new RooRealVar(nameMeas, nameMeas, vnDecMeasVec[i], -0.5, 0.5);
        vnDecMeas.add(*vnDecMeasVar);

        char *nameTrue = Form("vnDec%d", i);
        vnDecVar = new RooRealVar(nameTrue, nameTrue, 0., -0.5, 0.5);
        vnDec.add(*vnDecVar);
    }

    // now make the multivariate Gaussian
    RooMultiVarGaussian mvgVnDec("mvgVnDec", "mvgVnDec", vnDec, vnDecMeas, covVnDec);

    //
    // sample the distributions
    //
    RooDataSet *dataVnInc = mvgVnInc.generate(vnInc, nSamples);
    RooDataSet *dataVnDec = mvgVnDec.generate(vnDec, nSamples);
    RooDataSet *dataRgam = mvgRgam.generate(Rgam, nSamples);

    //
    // create direct photon vn data set based on the sampled values of VnInc, VnDec, Rgam
    //
    RooArgList vnDir;
    RooRealVar *vnDirVar;

    // RooFit variables for direct photon vn
    for (Int_t i = 0; i < nPtBins; i++) {
        char *name = Form("vnDir%d", i);
        vnDirVar = new RooRealVar(name, name, 0., -0.5, 0.5);
        vnDir.add(*vnDirVar);
    }

    RooDataSet dataVnDir("dataVnDir", "dataVnDir", vnDir);

    for (Int_t k = 0; k < dataVnInc->numEntries(); k++) {

        TVectorD vnDirTmp(nPtBins);
        bool b_all_v2dir_vals_in_allowed_range = true;

        for (Int_t i = 0; i < nPtBins; i++) {

            const RooArgSet *rowVnInc = dataVnInc->get(k);
            const RooArgSet *rowVnDec = dataVnDec->get(k);
            const RooArgSet *rowRgam = dataRgam->get(k);

            RooRealVar *vnIncRow = (RooRealVar *)rowVnInc->find(vnInc.at(i)->GetName());
            RooRealVar *vnDecRow = (RooRealVar *)rowVnDec->find(vnDec.at(i)->GetName());
            RooRealVar *RgamRow = (RooRealVar *)rowRgam->find(Rgam.at(i)->GetName());

            Double_t vnIncVal = vnIncRow->getVal();
            Double_t vnDecVal = vnDecRow->getVal();
            Double_t RgamVal = RgamRow->getVal();

            // cout << vnIncVal << " " << vnDecVal << " " << RgamVal << endl;

            Double_t vnDirVal = (RgamVal * vnIncVal - vnDecVal) / (RgamVal - 1.);

            if (fabs(vnDirVal) > 0.5) {
                b_all_v2dir_vals_in_allowed_range = false;
                break;
            }

            vnDirTmp(i) = vnDirVal;

            // ((RooRealVar*) vnDir.at(i))->setVal(vnDirVal);
            // dataVnDir.add(vnDir);
        }

        if (b_all_v2dir_vals_in_allowed_range) {
            for (Int_t i = 0; i < nPtBins; i++) {
                Double_t vnDirVal = vnDirTmp(i);
                ((RooRealVar *)vnDir.at(i))->setVal(vnDirVal);
                dataVnDir.add(vnDir);
            }
        }
    }

    // set covariance and correlation matrix
    covarianceMatrix = *dataVnDir.covarianceMatrix();
    correlationMatrix = *dataVnDir.correlationMatrix();

    //
    // some output
    //
    // dataVnDir.covarianceMatrix()->Print();

    // RooPlot* vnDirframe = ((RooRealVar*) vnDir.at(9))->frame(1000);
    // dataVnDir.plotOn(vnDirframe);
    // vnDirframe->Draw("h");

    // TH1* hTest = dataVnDir.createHistogram("test", (*(RooRealVar*) vnDir.at(9)));
    // hTest->Draw();
    for (Int_t i = 0; i < nPtBins; i++) {
        TString s = Form("%s_PtBin_%d", prefix.Data(), i);
        hVnDir[i] = dataVnDir.createHistogram(s.Data(), (*(RooRealVar *)vnDir.at(i)));
    }
};

// v2dir directly calculated from measured central values of
// R_gam, v2_inc, and v2_dec (i.e., without Bayes)
void v2dir_classic(const Int_t nPtBins, TVectorD &pt, const TVectorD &RgamMeasVec, const TMatrixDSym &covRgam,
                   const TVectorD &vnIncMeasVec, const TMatrixDSym &covVnInc, const TVectorD &vnDecMeasVec,
                   const TMatrixDSym &covVnDec, TGraphErrors *g_v2dir_classic) {

    for (Int_t i = 0; i < nPtBins; ++i) {
        if (RgamMeasVec(i) > 1) {

            Double_t R_gam = RgamMeasVec(i);
            Double_t R_gam_err = TMath::Sqrt(covRgam(i, i));
            Double_t v2_inc = vnIncMeasVec(i);
            Double_t v2_inc_err = TMath::Sqrt(covVnInc(i, i));
            Double_t v2_dec = vnDecMeasVec(i);
            Double_t v2_dec_err = TMath::Sqrt(covVnDec(i, i));

            Double_t v2dir = (R_gam * v2_inc - v2_dec) / (R_gam - 1.);
            g_v2dir_classic->SetPoint(i, pt(i), v2dir);

            // errors from Gaussian error propagation
            Double_t v2_dir_err_squared =
                (TMath::Power(R_gam - 1., 2) *
                     (TMath::Power(v2_dec_err, 2) + TMath::Power(R_gam, 2) * TMath::Power(v2_inc_err, 2)) +
                 TMath::Power(R_gam_err, 2) * TMath::Power(v2_dec - v2_inc, 2)) /
                TMath::Power(R_gam - 1., 4);
            Double_t v2_dir_err = TMath::Sqrt(v2_dir_err_squared);
            g_v2dir_classic->SetPointError(i, 0, v2_dir_err);

        } else {
            g_v2dir_classic->SetPoint(i, pt(i), -1);
            g_v2dir_classic->SetPointError(i, 0, 0);
        }
    }
}

Double_t DotProduct(const TVectorD &a, const TVectorD &b) {
    Int_t n = a.GetNrows();
    Double_t dotprod = 0;
    for (Int_t i = 0; i < n; i++)
        dotprod += a(i) * b(i);
    return dotprod;
}

Double_t chi2(const TVectorD &x, const TVectorD &mu, const TMatrixDSym &V) {
    TVectorD d = x - mu;
    TMatrixDSym V_copy = V;
    TMatrixDSym Vinv = V_copy.Invert();
    return DotProduct(d, Vinv * d);
}

Double_t p_value_to_n_sigma(const Double_t &pvalue) {

    // two-sided p-value calculation

    TF1 f("f", "1-TMath::Erf(x/sqrt(2.)) - [0]", -100, 100);
    f.SetParameter(0, pvalue);

    // create wrapper function
    ROOT::Math::WrappedTF1 wfrf(f);

    // create root finder
    ROOT::Math::BrentRootFinder brf;

    // set parameters of the method
    brf.SetFunction(wfrf, 0., 50.);
    brf.Solve();
    Double_t nsigma = brf.Root();

    return nsigma;
}
