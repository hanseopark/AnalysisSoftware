/*******************************************************************************
 ******  provided by Gamma Conversion Group, PWGGA,                        *****
 ******     Daniel Muehlheim, d.muehlheim@cern.ch                          ***** 
 *******************************************************************************/

#include <iostream>
#include <vector>
#include <fstream>

#include <TFile.h>
#include <TGrid.h>
#include <TSystem.h>
#include <TString.h>
#include <TObject.h>
#include <TIterator.h>
#include <TH1.h>
#include <TKey.h>

Bool_t readin(TString fileTxt, std::vector<TString> &vec);

void Grid_CheckFiles_Missing(TString folder = "/home/daniel/data/work/pcgGit/AnalysisSoftware")
{
    cout<<"Connecting to Alien..."<<endl;
    TGrid::Connect("alien://");
    cout<<"==============================="<<endl;
    cout<<"Successfully connected to Alien"<<endl;
    cout<<"==============================="<<endl;

    const Int_t nFiles = 1;
    TString Tag = "20160523";
    TString DataSetsFile[nFiles] = {"GammaConvCalo_101.root"};

    const Int_t nSets = 4;
    TString DataSets[nSets]={
      "LHC15h1c", "LHC15h1h", "LHC15h2b", "LHC15h2d"
    };
    TString DataSetsFolder[nSets]={
      "LHC15h1", "LHC15h1",
      "LHC15h2", "LHC15h2"
    };

    std::vector<TString> vecRuns;
    std::vector<TString> vecErrors[nSets];
    TString fDataSet;
    TString fPathLocal;
    TString fileTxt;

    Int_t nErr[nSets];

    for(Int_t i=0; i<nSets; i++)
    {
        nErr[i]=0;
        vecErrors[i].clear();
        vecRuns.clear();
		fDataSet = DataSets[i];
        if(fDataSet.CompareTo("")==0) continue;
        fileTxt = Form("%s/DownloadAndDataPrep/runlists/runNumbers%s.txt", folder.Data(), fDataSet.Data());
        cout << "\n------------------------------------------------------" << endl;
        if(!readin(fileTxt, vecRuns)) cout << "\n\n\n**********************ERROR, no Run Numbers could be found!**********************\n\n\n" << endl;
        cout << "------------------------------------------------------" << endl;

		Bool_t doNormalFolder = kFALSE;
		if(DataSetsFolder[i].IsNull()) doNormalFolder = kTRUE;

        for(Int_t j=0; j<(Int_t)vecRuns.size(); j++)
        {
          for(Int_t k=0; k<nFiles; k++)
          {
            if(doNormalFolder) fPathLocal = Form("%s/DataQA/%s/%s/%s", folder.Data(), Tag.Data(), fDataSet.Data(), vecRuns.at(j).Data());
            else fPathLocal = Form("%s/DataQA/%s/%s/%s", folder.Data(), Tag.Data(), DataSetsFolder[i].Data(), vecRuns.at(j).Data());
            fPathLocal+="/"; fPathLocal+=DataSetsFile[k];

            cout << endl;
            cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
            cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
            cout << "Processing file: " << fPathLocal.Data() << endl;
            cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
            cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

            TFile *fileCheck = new TFile(fPathLocal.Data());
            if(fileCheck->IsZombie()) {cout << "\n\t\t>>>>>>>>>>>>>>>>>>Info: ROOT-File |" << fPathLocal.Data() << "| does not exist! Continue...<<<<<<<<<<<<<<\n" << endl; continue;}

            TKey *key;
            TString nameMainDir;
            TIter next(fileCheck->GetListOfKeys());
            while ((key=(TKey*)next())){
              cout << Form(" - found TopDir: %s",key->GetName());
              nameMainDir = key->GetName();
            }
            cout << endl;
            if(nameMainDir.IsNull() || !nameMainDir.BeginsWith("Gamma")){cout << "ERROR, Unable to obtain valid name of MainDir:|" << nameMainDir.Data() << "|" << endl; return;}

            TList* TopDir = (TList*) fileCheck->Get(nameMainDir.Data());
              if(TopDir == NULL) {cout << "ERROR: TopDir not Found"<<endl; return;}
              else TopDir->SetOwner(kTRUE);

            Int_t a = 0;
            TString temp="";
            while (a<TopDir->GetEntries()){
               //cout << TopDir->At(a)->GetName() << endl;
               temp = TopDir->At(a)->GetName();
               if( temp.BeginsWith("ConvCuts_") ) break;
               a++;
            }

            TList* ConvCuts = (TList*) TopDir->FindObject(temp.Data());
              if(ConvCuts == NULL) {cout << "ERROR: " << temp.Data() << " not found in File" << endl; return;}
              else ConvCuts->SetOwner(kTRUE);

            temp.Replace(0,9,"");

            TH1D* NMissing = (TH1D*)ConvCuts->FindObject(Form("IsPhotonSelected %s",temp.Data()));
            Int_t missing = NMissing->GetBinContent(3) + NMissing->GetBinContent(4);
            if( missing > 0 ){
              cout << "\n\n\t**********************************************************************************************" << endl;
              cout << "\t**********************************************************************************************" << endl;
              cout << "\t*******************Err: Missing tracks / V0 in file contained!*******************" << endl;
              cout << "\t**********************************************************************************************" << endl;
              cout << "\t**********************************************************************************************\n" << endl;
              nErr[i]++;
              vecErrors[i].push_back(vecRuns.at(j));
            }

            delete TopDir;
            fileCheck->Close();
            delete fileCheck;
          }
        }
    }

    for(Int_t i=0; i<nSets; i++){
      cout << "DataSet: " << DataSets[i].Data() << ", number of errors: " << nErr[i] << endl;
      cout << "\t\tRuns: ";
      for(Int_t iRuns=0; iRuns < (Int_t)vecErrors[i].size(); iRuns++) cout << vecErrors[i].at(iRuns).Data() << ", ";
      cout << endl;
    }

    return;
}

Bool_t readin(TString fileTxt, std::vector<TString> &vec){
    cout << Form("Reading from %s...", fileTxt.Data()) << endl;
    fstream file;
    TString fVar;
    Int_t totalN=0;
    file.open(fileTxt.Data(), ios::in);
    if(file.good())
    {
        file.seekg(0L, ios::beg);
        if(fileTxt.Contains("binsJetJet")) cout << "Processing Bins: \"";
        else cout << "Processing Runs: \"";
        while(!file.eof())
        {
            file >> fVar;
            if(fVar.Sizeof()>1)
            {
                cout << fVar.Data() << ", ";
                vec.push_back(fVar);
                totalN++;
            }
        }
        cout << "\"" << endl;
    }
    file.close();
    if(fileTxt.Contains("binsJetJet"))  cout << "...done!\n\nIn total " << totalN << " Bins will be processed!" << endl;
    else cout << "...done!\n\nIn total " << totalN << " Runs will be processed!" << endl;
    if(totalN > 0) return kTRUE;
    else return kFALSE;
}
