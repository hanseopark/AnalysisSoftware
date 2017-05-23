/*******************************************************************************
 ******  provided by Gamma Conversion Group, PWGGA,                        *****
 ******     Daniel Muehlheim, d.muehlheim@cern.ch                          ***** 
 *******************************************************************************/

#include <algorithm>
#include <vector>
#include <iterator>
#include <iostream>
#include <fstream>
#include <TObject.h>
#include <TString.h>
#include <TFile.h>
#include <TGrid.h>
#include <TSystem.h>
#include <TKey.h>
#include <TH1D.h>
#include <iostream>
#include <fstream>

void ChangeStrucToStd(TString nameInputFile, TString namefileOutput, TString nameInputList);
Bool_t copyAlien2Local(TString loc, TString rem);
Bool_t readin(TString fileRuns, std::vector<TString> &vec){
    cout << Form("Reading from %s...", fileRuns.Data()) << endl;
    std::fstream file;
    TString fVar;
    Int_t totalN=0;
    file.open(fileRuns.Data(), ios::in);
    if(file.good())
    {
        file.seekg(0L, ios::beg);
        while(!file.eof())
        {
            file >> fVar;
            if(fVar.Sizeof()>1)
            {
                vec.push_back(fVar);
                totalN++;
            }
        }
    }
    file.close();
    if(totalN > 0) return kTRUE;
    else return kFALSE;
}

void Grid_CopyFilesJetJet(TString system = "pp", TString type = "ESD", TString folder = "/home/daniel/data/work/Grid", TString folderSoftware = "/home/daniel/data/work/pcgGit/AnalysisSoftware")
{
    cout<<"Connecting to Alien..."<<endl;
    TGrid::Connect("alien://");
    cout<<"==============================="<<endl;
    cout<<"Successfully connected to Alien"<<endl;
    cout<<"==============================="<<endl;

    //*********************************************************************************************************************************
    //*********************************************************************************************************************************
    //*********************************************************************************************************************************

//    const Int_t nSets = 2;
//    const Int_t nData = 0;
//    TString DataSets[nSets]={
//      "LHC16c2",
//      "LHC16c2_plus"
//    };
//    TString PrefixDataSets[nSets]={
//      "/alice/sim/2016/LHC16c2/",
//      "/alice/sim/2016/LHC16c2_plus/"
//    };

//    TString AODfiltering[nSets]={
//      "/AOD185",
//      ""
//    };

//    TString train = "Legotrain-vAN-20170518-8TeV-std_EMCal";
//    Int_t trainRuns[nSets] = {413,414};
//    TString runlist[nSets] = {"merge","merge"};

//    const Int_t nFiles = 3;
//    TString Files[nFiles] = {
//      "GammaCalo_101","GammaCalo_126","GammaCalo_146"
//                            };

//    const Int_t nMerge = 1;
//    TString strMerge[nMerge]={"LHC16c2_merge"};
//    std::vector<Int_t> mergeVec[nMerge];
//    std::vector<Int_t>::iterator it;
//    for(Int_t i=0; i<nSets; i++){
//      if(0<=i && i<=1) mergeVec[0].push_back(i);
//    }

    //*********************************************************************************************************************************
    //*********************************************************************************************************************************
    //*********************************************************************************************************************************

    const Int_t nSets = 2;
    const Int_t nData = 0;
    TString DataSets[nSets]={
      "LHC16c2",
      "LHC16c2_plus"
    };
    TString PrefixDataSets[nSets]={
      "/alice/sim/2016/LHC16c2/",
      "/alice/sim/2016/LHC16c2_plus/"
    };

    TString AODfiltering[nSets]={
      "/AOD185",
      ""
    };

    TString train = "Legotrain-vAN-20170519-8TeV-EMCal_openAngleStudies";
    Int_t trainRuns[nSets] = {419,420};
    TString runlist[nSets] = {"merge","merge"};

    const Int_t nFiles = 4;
    TString Files[nFiles] = {
      "GammaCalo_169","GammaCalo_170","GammaCalo_171","GammaCalo_172"
                            };

    const Int_t nMerge = 1;
    TString strMerge[nMerge]={"LHC16c2_merge"};
    std::vector<Int_t> mergeVec[nMerge];
    std::vector<Int_t>::iterator it;
    for(Int_t i=0; i<nSets; i++){
      if(0<=i && i<=1) mergeVec[0].push_back(i);
    }

    //---------------------------------------------------------------------------------------------------

    TString alienFolder;
    TString alienFolder_MC;
    TString sAOD;

    if(type.CompareTo("ESD")==0){
      sAOD = "";
      alienFolder = Form("/alice/cern.ch/user/a/alitrain/PWGGA/GA_%s/",system.Data());
      alienFolder_MC = Form("/alice/cern.ch/user/a/alitrain/PWGGA/GA_%s_MC/",system.Data());
    }else{
      sAOD = "_AOD";
      alienFolder = Form("/alice/cern.ch/user/a/alitrain/PWGGA/GA_%s_AOD/",system.Data());
      alienFolder_MC = Form("/alice/cern.ch/user/a/alitrain/PWGGA/GA_%s_MC_AOD/",system.Data());
    }
    gSystem->Exec(Form("alien_ls %s > tempData.log",alienFolder.Data()));
    gSystem->Exec(Form("alien_ls %s > tempMC.log", alienFolder_MC.Data()));

    TString strTrain[nSets];
    std::vector<TString> vecStrTrain;
    std::vector<TString> vecStrTrainMC;
    if(!readin("tempData.log", vecStrTrain)) cout << "\n\n\n**********************ERROR!**********************\n\n\n" << endl;
    if(!readin("tempMC.log", vecStrTrainMC)) cout << "\n\n\n**********************ERROR!**********************\n\n\n" << endl;

    for(Int_t i=0; i<nSets; i++){
      TString temp;
      TString tempRuns;
      if(i<nData){
        for(Int_t j=0; j<(Int_t)vecStrTrain.size(); j++){
          tempRuns = Form("%i",trainRuns[i]);
          temp = vecStrTrain.at(j);
          if(temp.BeginsWith(tempRuns)){
            strTrain[i] = temp;
            break;
          }
        }
      }else{
        for(Int_t j=0; j<(Int_t)vecStrTrainMC.size(); j++){
          tempRuns = Form("%i",trainRuns[i]);
          temp = vecStrTrainMC.at(j);
          if(temp.BeginsWith(tempRuns)){
            strTrain[i] = temp;
            break;
          }
        }
      }
    }

    Int_t nErr[nSets];
    std::vector<TString> vecRuns;
    std::vector<TString> vecBins;
    TString fDataSet;
    TString fileTxt;
    TString fPathGrid;
    TString fPathLocal;

    fPathLocal = Form("%s/%s", folder.Data(), train.Data());
    gSystem->Exec(Form("mkdir -p %s",fPathLocal.Data()));

    fstream fLog;
    fLog.open(Form("%s/%s/DL.log",folder.Data(),train.Data()), ios::out);

    TString mergeSets[nMerge][nFiles];
    for(Int_t iFiles=0; iFiles<nFiles; iFiles++){
      for(Int_t iM=0; iM<nMerge; iM++){
        TString mergeS = Form("%s/%s_%s.root", fPathLocal.Data(), strMerge[iM].Data(), Files[iFiles].Data());
        mergeSets[iM][iFiles] = Form("hadd -f -k %s",mergeS.Data());
        //cout << mergeSets[iM][iFiles] << endl;
      }
    }

    for(Int_t i=0; i<nSets; i++)
    {
        nErr[i]=0;
        //---------------------------------------------------------------------------------------------------
        if(runlist[i].CompareTo("runwise") == 0){
          fDataSet = DataSets[i];
          if(fDataSet.CompareTo("")==0) continue;

          TString mergePeriod[nFiles];
          for(Int_t iFiles=0; iFiles<nFiles; iFiles++){
            TString mergeP = Form("%s/%s_%s.root", fPathLocal.Data(), DataSets[i].Data(), Files[iFiles].Data());
            mergePeriod[iFiles] = Form("hadd -f -k %s",mergeP.Data());
          }

          fileTxt = Form("%s/DownloadAndDataPrep/runlists/runNumbers%s.txt", folderSoftware.Data(), fDataSet.Data());
          cout << "\n------------------------------------------------------" << endl;
          if(!readin(fileTxt, vecRuns)) cout << "\n\n\n**********************ERROR, no Run Numbers could be found!**********************\n\n\n" << endl;
          cout << "------------------------------------------------------" << endl;

          fileTxt = Form("%s/DownloadAndDataPrep/binsJetJet%s.txt", folderSoftware.Data(), fDataSet.Data());
          cout << "\n------------------------------------------------------" << endl;
          if(!readin(fileTxt, vecBins)) cout << "\n\n\n**********************ERROR, no Jet Jet Bins could be found!**********************\n\n\n" << endl;
          cout << "------------------------------------------------------" << endl;

          for(Int_t j=0; j<(Int_t)vecRuns.size(); j++)
          {
              Int_t nBins = (Int_t)vecBins.size();
              TString mergeBins[nFiles];
              for(Int_t iFiles=0; iFiles<nFiles; iFiles++){
                TString mergeB = Form("%s/%s/%s_%s_%s.root", fPathLocal.Data(), strTrain[i].Data(), vecRuns.at(j).Data(), DataSets[i].Data(), Files[iFiles].Data());
                mergeBins[iFiles] = Form("hadd -f -k %s",mergeB.Data());
              }
              for(Int_t b=0; b<nBins; b++)
              {
                for(Int_t k=0; k<nFiles; k++)
                {
                    fPathGrid = Form("%s%s/%s%s/PWGGA/GA_pp_MC%s/%s/%s.root", PrefixDataSets[i].Data(), vecBins.at(b).Data(), vecRuns.at(j).Data(), AODfiltering[i].Data(), sAOD.Data(), strTrain[i].Data(), Files[k].Data());
                    TString fPathTemp = Form("%s/%s", fPathLocal.Data(), strTrain[i].Data());
                    gSystem->Exec(Form("mkdir -p %s",fPathTemp.Data()));

                    fPathTemp+=Form("/%s_%s_%s_%s.root",vecRuns.at(j).Data(),vecBins.at(b).Data(),DataSets[i].Data(),Files[k].Data());
                    cout << endl;
                    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
                    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
                    cout << "Copying from (grid): " << fPathGrid.Data() << endl;
                    cout << "Copying to (local): " << fPathTemp.Data() << endl;
                    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
                    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

                    TFile fileCheck(fPathTemp.Data());
                    if(!fileCheck.IsZombie()) {
                      cout << "\n\t\t>>>>>>>>>>>>>>>>>>Info: ROOT-File |" << fPathTemp.Data() << "| does already exist! Continue...<<<<<<<<<<<<<<\n" << endl;
                      TFile* fFile        = new TFile(fPathTemp.Data(),"READ");
                      TKey *key;
                      TIter next(fFile->GetListOfKeys());
                      TString nameMainDir = "";
                      while ((key=(TKey*)next())){
                          cout << Form(" - found TopDir: %s",key->GetName());
                          nameMainDir  = key->GetName();
                      }
                      cout << endl;
                      if(nameMainDir.IsNull() || !nameMainDir.BeginsWith("Gamma")){
                          cout << "ERROR, Unable to obtain valid name of MainDir:|" << nameMainDir.Data() << "|" << endl; return;
                      }

                      TList *listInput    = (TList*)fFile->Get(nameMainDir.Data());
                      listInput->SetOwner(kTRUE);
                      if(!listInput) {cout << "ERROR: Could not find main dir: " << nameMainDir.Data() << " in file! Returning..." << endl; return;}
                      TList *listCuts     = (TList*)listInput->At(0);
                      TString nameCuts    = listCuts->GetName();
                      TList* TopContainer             = (TList*) listInput->FindObject(nameCuts.Data());
                      if(TopContainer == NULL) {cout << "ERROR: " << nameCuts.Data() << " not found in File" << endl; return;}
                          else TopContainer->SetOwner(kTRUE);
                      TList *listCuts2     = (TList*)TopContainer->At(0);
                      TString nameCuts2    = listCuts2->GetName();
                      TList* ESDContainer             = (TList*) TopContainer->FindObject(nameCuts2.Data());
                      if(ESDContainer == NULL) {cout << "ERROR: " << nameCuts2.Data() << " not found in File" << endl; return;}
                          else ESDContainer->SetOwner(kTRUE);
                      TH1D* fHistNEvents              = (TH1D*)ESDContainer->FindObject("NEvents");
                      fLog << "Run: '" << vecRuns.at(j).Data() << "', Bin: '" << vecBins.at(b).Data() << "', File: '" << Files[k].Data() << "' - NEntries of NEvents: " << fHistNEvents->GetEntries() << ", Mean: '" << fHistNEvents->GetMean() << "'" << endl;
                      delete listInput;
                      fFile->Close();
                      delete fFile;
                      mergeBins[k] += Form(" %s",fPathTemp.Data());
                      continue;
                    }

                    if(copyAlien2Local(fPathGrid,fPathTemp)){
                      TFile* fFile        = new TFile(fPathTemp.Data(),"READ");
                      TKey *key;
                      TIter next(fFile->GetListOfKeys());
                      TString nameMainDir = "";
                      while ((key=(TKey*)next())){
                          cout << Form(" - found TopDir: %s",key->GetName());
                          nameMainDir  = key->GetName();
                      }
                      cout << endl;
                      if(nameMainDir.IsNull() || !nameMainDir.BeginsWith("Gamma")){
                          cout << "ERROR, Unable to obtain valid name of MainDir:|" << nameMainDir.Data() << "|" << endl; return;
                      }

                      TList *listInput    = (TList*)fFile->Get(nameMainDir.Data());
                      listInput->SetOwner(kTRUE);
                      if(!listInput) {cout << "ERROR: Could not find main dir: " << nameMainDir.Data() << " in file! Returning..." << endl; return;}
                      TList *listCuts     = (TList*)listInput->At(0);
                      TString nameCuts    = listCuts->GetName();
                      TList* TopContainer             = (TList*) listInput->FindObject(nameCuts.Data());
                      if(TopContainer == NULL) {cout << "ERROR: " << nameCuts.Data() << " not found in File" << endl; return;}
                          else TopContainer->SetOwner(kTRUE);
                      TList *listCuts2     = (TList*)TopContainer->At(0);
                      TString nameCuts2    = listCuts2->GetName();
                      TList* ESDContainer             = (TList*) TopContainer->FindObject(nameCuts2.Data());
                      if(ESDContainer == NULL) {cout << "ERROR: " << nameCuts2.Data() << " not found in File" << endl; return;}
                          else ESDContainer->SetOwner(kTRUE);
                      TH1D* fHistNEvents              = (TH1D*)ESDContainer->FindObject("NEvents");
                      fLog << "Run: '" << vecRuns.at(j).Data() << "', Bin: '" << vecBins.at(b).Data() << "', File: '" << Files[k].Data() << "' - NEntries of NEvents: " << fHistNEvents->GetEntries() << ", Mean: '" << fHistNEvents->GetMean() << "'" << endl;
                      delete listInput;
                      fFile->Close();
                      delete fFile;

                      ChangeStrucToStd(fPathTemp.Data(),fPathTemp.Data(),Files[k].Data());
                      mergeBins[k] += Form(" %s",fPathTemp.Data());
                    }
                    else{
                      cout << "\n\n\t**********************************************************************************************" << endl;
                      cout << "\t**********************************************************************************************" << endl;
                      cout << "\t*******************Err: copyAlien2Local()!****************************************************" << endl;
                      cout << "\t**********************************************************************************************" << endl;
                      cout << "\t**********************************************************************************************\n" << endl;
                      nErr[i]++;
                   }
                }
             }
             cout << "\n------------------------------------------------------" << endl;
             cout << "Merging Bins: " << endl;
             cout << "------------------------------------------------------\n" << endl;

             for(Int_t iFiles=0; iFiles<nFiles; iFiles++){
               cout << "Merging " << mergeBins[iFiles].Data() << " ..." << endl;
               gSystem->Exec(mergeBins[iFiles].Data());
               cout << "done!" << endl;

               mergePeriod[iFiles] += Form(" %s/%s/%s_%s_%s.root", fPathLocal.Data(), strTrain[i].Data(), vecRuns.at(j).Data(), DataSets[i].Data(), Files[iFiles].Data());
             }
          }

          cout << "\n------------------------------------------------------" << endl;
          cout << "Merging Runs: " << endl;
          cout << "------------------------------------------------------\n" << endl;

          for(Int_t iFiles=0; iFiles<nFiles; iFiles++){
            cout << "Merging " << mergePeriod[iFiles].Data() << " ..." << endl;
            gSystem->Exec(mergePeriod[iFiles].Data());
            cout << "done!" << endl;

            for(Int_t iM=0; iM<nMerge; iM++){
              it = find( mergeVec[iM].begin(), mergeVec[iM].end(), i);
              if( it!=mergeVec[iM].end() ) mergeSets[iM][iFiles] += Form("%s/%s_%s.root", fPathLocal.Data(), DataSets[i].Data(), Files[iFiles].Data());
            }
          }

        //---------------------------------------------------------------------------------------------------
        }else{
          for(Int_t k=0; k<nFiles; k++)
          {
            if(i<nData) fPathGrid = Form("%s%s/%s/%s.root", alienFolder.Data(), strTrain[i].Data(), runlist[i].Data(), Files[k].Data());
            else fPathGrid = Form("%s%s/%s/%s.root", alienFolder_MC.Data(), strTrain[i].Data(), runlist[i].Data(), Files[k].Data());

            TString fPathTemp = Form("%s/%s", fPathLocal.Data(), strTrain[i].Data());
            gSystem->Exec(Form("mkdir -p %s",fPathTemp.Data()));

            fPathTemp+=Form("/%s_",DataSets[i].Data());
            fPathTemp+=Files[k].Data();
            fPathTemp+=".root";

            cout << endl;
            cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
            cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
            cout << "Copying from (grid): " << fPathGrid.Data() << endl;
            cout << "Copying to (local): " << fPathTemp.Data() << endl;
            cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
            cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

            TFile fileCheck(fPathTemp.Data());
            if(!fileCheck.IsZombie()) {
              cout << "\n\t\t>>>>>>>>>>>>>>>>>>Info: ROOT-File |" << fPathTemp.Data() << "| does already exist! Continue...<<<<<<<<<<<<<<\n" << endl;
              for(Int_t iM=0; iM<nMerge; iM++){
                it = find( mergeVec[iM].begin(), mergeVec[iM].end(), i);
                if( it!=mergeVec[iM].end() ) mergeSets[iM][k] += Form(" %s/%s_%s.root", fPathLocal.Data(), DataSets[i].Data(), Files[k].Data());
              }
              gSystem->Exec(Form("ln -s %s %s/%s_%s.root",fPathTemp.Data(), fPathLocal.Data(), DataSets[i].Data(), Files[k].Data()));
              continue;
            }

            if(copyAlien2Local(fPathGrid,fPathTemp)){
              ChangeStrucToStd(fPathTemp.Data(),fPathTemp.Data(),Files[k].Data());
              for(Int_t iM=0; iM<nMerge; iM++){
                it = find( mergeVec[iM].begin(), mergeVec[iM].end(), i);
                if( it!=mergeVec[iM].end() ) mergeSets[iM][k] += Form(" %s/%s_%s.root", fPathLocal.Data(), DataSets[i].Data(), Files[k].Data());
              }
              gSystem->Exec(Form("ln -s %s %s/%s_%s.root",fPathTemp.Data(), fPathLocal.Data(), DataSets[i].Data(), Files[k].Data()));
              continue;
            }
            else{
              cout << "\n\n\t**********************************************************************************************" << endl;
              cout << "\t**********************************************************************************************" << endl;
              cout << "\t*******************Err: copyAlien2Local()!****************************************************" << endl;
              cout << "\t**********************************************************************************************" << endl;
              cout << "\t**********************************************************************************************\n" << endl;
              nErr[i]++;
            }
          }
        }
        //---------------------------------------------------------------------------------------------------
    }

    cout << "\n------------------------------------------------------" << endl;
    cout << "Merging: " << endl;
    cout << "------------------------------------------------------\n" << endl;

    for(Int_t iFiles=0; iFiles<nFiles; iFiles++){
      for(Int_t iM=0; iM<nMerge; iM++){
        cout << "Merging " << strMerge[iM].Data() << " - " << Files[iFiles] << " ..." << endl;
        gSystem->Exec(mergeSets[iM][iFiles].Data());
        cout << "done!" << endl;
      }
    }

    for(Int_t i=0; i<nSets; i++){
      cout << "DataSet: " << DataSets[i].Data() << ", number of errors: " << nErr[i] << endl;
    }

    gSystem->Exec("rm tempData.log");
    gSystem->Exec("rm tempMC.log");

    return;
}



void ChangeStrucToStd(TString nameInputFile, TString namefileOutput, TString nameInputList){

   TFile *fileInput = new TFile(nameInputFile.Data());
   cout << fileInput << endl;

   TList *listInput =(TList*)fileInput->Get(nameInputList.Data());
   if (listInput == NULL){
      return;
   }else listInput->SetOwner();

   TObjArray *rArr = nameInputList.Tokenize("_");
   TObjString* rString = (TObjString*)rArr->At(0);
   TString string = rString->GetString();

   TFile *fileOutput = new TFile(namefileOutput,"RECREATE");
   TList *listOutput =(TList*)fileOutput->Get(string.Data());
   Bool_t kNewList = kFALSE;
   if (!listOutput){
      kNewList = kTRUE;
      listOutput = new TList();
      listOutput->SetName(string.Data());
   }

   for(Int_t i = 0; i<listInput->GetSize(); i++){
      TList *listToSave = (TList*)listInput->At(i);
      TString dirname = listToSave->GetName();
      cout<<dirname<<endl;
      if(listToSave){
         cout<<"found"<<endl;
         listOutput->Add(listToSave);
      }
   }

   listOutput->Write("",TObject::kSingleKey);
   delete listOutput;
   fileOutput->Close();
   delete fileOutput;

   delete listInput;
   fileInput->Close();
   delete fileInput;

   return;
}

Bool_t copyAlien2Local(TString loc, TString rem){
   TString sl(Form("alien://%s", loc.Data()));
   TString sr(Form("file://%s", rem.Data()));
   Bool_t ret = TFile::Cp(sl,sr);
   if (!ret) cout << Form("ERROR: Failed to copy %s to %s", sl.Data(), sr.Data()) << endl;
   return ret;
}
