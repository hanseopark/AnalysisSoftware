#include <iostream>
#include <vector>
#include <fstream>

#include <TFile.h>
#include <TGrid.h>
#include <TSystem.h>
#include <TString.h>
#include <TObject.h>
#include <TH1.h>
#include <TKey.h>

Bool_t readin(TString fileTxt, std::vector<TString> &vec);

void Grid_ReadJetJet(TString folder = "/home/daniel/data/work/pcgGit/AnalysisSoftware")
{
   //JetJet
    const Int_t nFiles = 2;
    TString Tag = "20160513";
    TString DataSetsFile[nFiles] = {"GammaCalo_111.root","GammaConvCalo_31.root"};
    TString DataSetsCut[nFiles] = {"00010113_1111100063032220000_0163103100000050","00010113_00200009327000008250400000_1111100063032230000_0163103100000010"};

    const Int_t nSets = 1;
    TString DataSets[nSets]={"LHC16c2"};
    TString DataSetsFolder[nSets]={"LHC16c2"};

//  const Int_t nFiles = 2;
//  TString Tag = "20160411";
//  TString DataSetsFile[nFiles] = {"GammaConvCalo_1.root","GammaConvCalo_2.root"};
//  TString DataSetsCut[nFiles] = {"80000013_00200009327000008250400000_1111100053022230000_0163103100000010","80085013_00200009327000008250400000_1111100053022230000_0163103100000010"};

//  const Int_t nSets = 3;
//  TString DataSets[nSets]={"LHC16c3a","LHC16c3b","LHC16c3c"};
//  TString DataSetsFolder[nSets]={"LHC16c3a","LHC16c3b","LHC16c3c"};


    std::vector<TString> vecRuns;
    std::vector<TString> vecBins;
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
        vecBins.clear();
		fDataSet = DataSets[i];
        if(fDataSet.CompareTo("")==0) continue;
        fileTxt = Form("%s/DownloadAndDataPrep/runlists/runNumbers%s.txt", folder.Data(), fDataSet.Data());
        cout << "\n------------------------------------------------------" << endl;
        if(!readin(fileTxt, vecRuns)) cout << "\n\n\n**********************ERROR, no Run Numbers could be found!**********************\n\n\n" << endl;
        cout << "------------------------------------------------------" << endl;

        fileTxt = Form("%s/DownloadAndDataPrep/binsJetJet%s.txt", folder.Data(), fDataSet.Data());
        cout << "\n------------------------------------------------------" << endl;
        if(!readin(fileTxt, vecBins)) cout << "\n\n\n**********************ERROR, no Jet Jet Bins could be found!**********************\n\n\n" << endl;
        cout << "------------------------------------------------------" << endl;

		Bool_t doNormalFolder = kFALSE;
		if(DataSetsFolder[i].IsNull()) doNormalFolder = kTRUE;

//        for(Int_t j=0; j<(Int_t)vecRuns.size(); j++)
//        {
            Int_t nBins = (Int_t)vecBins.size();
            for(Int_t b=0; b<nBins; b++)
            {
              cout << "\n\n\t**********************************************************************************************" << endl;
              cout << "\t**********************************************************************************************" << endl;
              cout << "\t\tBin:\t\t" << vecBins.at(b) << endl;
              for(Int_t k=0; k<nFiles; k++)
              {
                  if(doNormalFolder) fPathLocal = Form("%s/DataQA/%s/%s/%s", folder.Data(), Tag.Data(), fDataSet.Data(), vecBins.at(b).Data()/*, vecRuns.at(j).Data()*/);
                  else fPathLocal = Form("%s/DataQA/%s/%s/%s", folder.Data(), Tag.Data(), DataSetsFolder[i].Data(), vecBins.at(b).Data()/*, vecRuns.at(j).Data()*/);

                  fPathLocal+="/"; fPathLocal+=DataSetsFile[k];

                  TFile *fileCheck = new TFile(fPathLocal.Data());
                  if(fileCheck->IsZombie()) {cout << "\n\t\t>>>>>>>>>>>>>>>>>>Info: ROOT-File |" << fPathLocal.Data() << "| does not exist! Continue...<<<<<<<<<<<<<<\n" << endl; continue;}

                  cout << "Processing file: " << fPathLocal.Data();
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
                  TList* TopContainer= (TList*) TopDir->FindObject(Form("Cut Number %s",DataSetsCut[k].Data()));
                  if(TopContainer == NULL) {cout << "ERROR: " << Form("Cut Number %s",DataSetsCut[k].Data()) << " not found in File" << endl; return;}
                  else TopContainer->SetOwner(kTRUE);
                  TList* ESDContainer = (TList*) TopContainer->FindObject(Form("%s ESD histograms",DataSetsCut[k].Data()));
                  if(ESDContainer == NULL) {cout << "ERROR: " << Form("%s ESD histograms",DataSetsCut[k].Data()) << " not found in File" << endl; return;}
                  else ESDContainer->SetOwner(kTRUE);

                  TH1D* NTrials = (TH1D*)ESDContainer->FindObject("NTrials");
                  TH1D* XSection = (TH1D*)ESDContainer->FindObject("XSection");

                  cout << "\t\t" << DataSetsFile[k].Data() << ":\t\t NTrials: |" << NTrials->GetBinContent(1) << "|," << " XSection: |" << XSection->GetBinContent(1) << "|," << " N_{genEvents}: |" << NTrials->GetEntries() << "|" << endl;
                  cout << XSection->GetBinContent(1)/(NTrials->GetBinContent(1)/NTrials->GetEntries()) << endl;

                  fileCheck->Close();
                  delete fileCheck;

                  delete TopDir;
              }
              cout << "\t**********************************************************************************************" << endl;
              cout << "\t**********************************************************************************************\n" << endl;
            }
//          }
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
