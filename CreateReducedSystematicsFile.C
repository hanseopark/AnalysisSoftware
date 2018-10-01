
/****************************************************************************************************************************
 ******     provided by Gamma Conversion Group, PWGGA,                                                                  *****
 ******        Daniel Mühlheim, d.muehlheim@cern.ch                                                                     *****
 *****************************************************************************************************************************/
#include <TString.h>
#include <TObject.h>
#include <TH1.h>
#include <TH2.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TGaxis.h>
#include <TLegend.h>
#include <TFrame.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TF1.h>
#include <TASImage.h>
#include <TMarker.h>
#include <TDatime.h>
#include <TArrow.h>
#include <TROOT.h>
#include <TObjString.h>
#include <TFile.h>
#include "TMath.h"


#include <Riostream.h>
#include <string>
#include <sstream>
#include <fstream>
#include <stdlib.h>
#include <fstream>
#include <math.h>

#include "CommonHeaders/PlottingGammaConversionHistos.h"
#include "CommonHeaders/PlottingGammaConversionAdditional.h"
#include "CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"

using namespace std;

Int_t GetPosInVec( std::vector<TString> vec, TString lookup){
    for(Int_t i=0; i<(Int_t)vec.size(); i++){
        if(vec.at(i).CompareTo(lookup) == 0) return i;
    }
    cout << "ERROR in GetPosInVec, did not find '" << lookup.Data() << "'" << endl;
    return -1;
}

void CreateReducedSystematicsFile(  TString inputfileName               = "",
                                    TString fileNameWithErrorsToBeUsed  = "",
                                    Int_t nBinsPt                       = 0,
                                    TString appendToSysName             = "",
                                    TString addUncName                  = "",
                                    Double_t addPtIndUnc                = -1
                                  ){



    vector<TString> sysTaken;
    ifstream fileSysNames(fileNameWithErrorsToBeUsed.Data());
    string ll;
    while (getline(fileSysNames,ll)){
        istringstream lll(ll);
        TString temp        = "";
        lll >> temp;
        sysTaken.push_back(temp);
        if( temp.IsNull()  )
            break;
        cout << temp.Data() << endl;
    }

    vector<TString> ptSys[nBinsPt+1];
    vector<Bool_t> isTaken;
    ifstream fileInput(inputfileName.Data());
    Int_t iPt = 0;
    string line;
    while (getline(fileInput, line) && iPt < 100) {
        istringstream ss(line);
        TString temp        = "";
        TString tempBin     = "";
        Int_t iESource = 0;
        cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        cout << "\t*** reading line " << iPt << " ***" << endl;
        while(ss && iESource < 100){
            ss >> temp;
            if(iPt > 0 && iESource == 0){
                tempBin = temp;
                if( tempBin.IsNull()  ){
                    cout << "INFO: skipping pT " << tempBin.Atof() << endl;
                    break;
                }
            }
            if( !(iPt==0 && temp.CompareTo("bin")==0) && !temp.IsNull() && iPt < (nBinsPt+1)){
                ptSys[iPt].push_back(temp);
                cout << temp.Data() << ", ";

                Bool_t valid    = kFALSE;
                if (iESource == 0){
                    valid   = kTRUE;
                } else {
                    for (Int_t so = 0; so < (Int_t)sysTaken.size(); so++){
                        if (temp.CompareTo(sysTaken.at(so)) == 0)
                            valid   = kTRUE;
                    }
                }
                isTaken.push_back(valid);
                iESource++;
            }
        }

        if(iPt == 0){
            ptSys[iPt++].push_back("TotalError");
            cout << "TotalError";
            isTaken.push_back(kFALSE);
            iESource++;
        } else if(!tempBin.IsNull() )
            iPt++;
        cout << "\n\t***" << iESource << " columns read ***" << endl;
        cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n" << endl;
    }

    TObjArray *tempArr  = inputfileName.Tokenize(".");

    TString outputfileName  = (TString)((TObjString*)tempArr->At(0))->GetString()+"_red."+(TString)((TObjString*)tempArr->At(1))->GetString();
    TString outputfileName2 = (TString)((TObjString*)tempArr->At(0))->GetString()+"_redSimple."+(TString)((TObjString*)tempArr->At(1))->GetString();

    cout << "writing reduced file to " << outputfileName.Data() << " and "<< outputfileName2.Data() << endl;
    fstream fCorr;
    fstream fCorr2;
    fCorr.open(outputfileName.Data(), ios::out);
    fCorr2.open(outputfileName2.Data(), ios::out);
    for (Int_t cPt = 0; cPt < nBinsPt+1; cPt++){
        fCorr << ptSys[cPt].at(0) ;
        if (cPt == 0) fCorr << " bin";
        fCorr << "\t" ;
        if (cPt != 0) fCorr2 << ptSys[cPt].at(0) << "\t";
        Double_t sysErrTot  = 0;
        for (Int_t source = 1; source< (Int_t)ptSys[cPt].size(); source++){
            if (isTaken.at(source)){
                if (cPt == 0)
                    fCorr << ptSys[cPt].at(source)<< appendToSysName.Data() << "\t";
                else
                    fCorr << ptSys[cPt].at(source) << "\t";
                sysErrTot   = sysErrTot+ ((TString)ptSys[cPt].at(source)).Atof()*((TString)ptSys[cPt].at(source)).Atof();
            }
        }
        if (addUncName.CompareTo("")){
            if (cPt == 0){
                fCorr << addUncName.Data()<< appendToSysName.Data() << "\t";
            } else {
                fCorr << addPtIndUnc << "\t";
                sysErrTot   = sysErrTot+addPtIndUnc*addPtIndUnc;
            }
        }
        if (cPt == 0){
            fCorr << "TotalErrorUncorr";
            fCorr2 << "TotalErrorUncorr";
        } else {
            fCorr << TMath::Sqrt(sysErrTot) ;
            fCorr2 << TMath::Sqrt(sysErrTot) ;
        }
        fCorr << endl;
        fCorr2 << endl;
    }
    fCorr.close();
    fCorr2.close();
}