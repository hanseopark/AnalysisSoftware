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
//      "LHC16c2"
//    };
//    TString PrefixDataSets[nSets]={
//      "/alice/sim/2016/LHC16c2/",
//      "/alice/sim/2016/LHC16c2/"
//    };

////    TString train = "Legotrain-vAN-20160608-8TeV-JetJet-ECalib";
////    Int_t trainRuns[nSets] = {2161,2162};
////    TString runlist[nSets] = {"runwise","merge"};

//    TString train = "Legotrain-vAN-20160611-8TeV-JetJet-ECalib";
//    Int_t trainRuns[nSets] = {2170};
//    TString runlist[nSets] = {"runwise"};

//    const Int_t nFiles = 4;
//    TString Files[nFiles] = {"GammaCalo_129","GammaCalo_149","GammaConvCalo_140","GammaConvCalo_170"};

    //*********************************************************************************************************************************
    //*********************************************************************************************************************************
    //*********************************************************************************************************************************


//    const Int_t nSets = 1;
//    const Int_t nData = 0;
//    TString DataSets[nSets]={
//      "LHC16c2"
//    };
//    TString PrefixDataSets[nSets]={
//      "/alice/sim/2016/LHC16c2/"
//    };

//    TString train = "Legotrain-vAN-20160615-8TeV-JetJet-ECalib_std";
//    Int_t trainRuns[nSets] = {2171};
//    TString runlist[nSets] = {"merge"};

//    const Int_t nFiles = 6;
//    TString Files[nFiles] = {"GammaCalo_101","GammaCalo_129","GammaCalo_149","GammaConvCalo_101","GammaConvCalo_140","GammaConvCalo_170"};

    //*********************************************************************************************************************************
    //*********************************************************************************************************************************
    //*********************************************************************************************************************************

//        const Int_t nSets = 1;
//        const Int_t nData = 0;
//        TString DataSets[nSets]={
//          "LHC16c2"
//        };
//        TString PrefixDataSets[nSets]={
//          "/alice/sim/2016/LHC16c2/"
//        };

//        TString train = "Legotrain-vAN-20160617-8TeV-Systematics_Calo_EMC7";
//        //Int_t trainRuns[nSets] = {2197};
//        Int_t trainRuns[nSets] = {2198};
//        TString runlist[nSets] = {"merge"};

//        //const Int_t nFiles = 5;
//        //TString Files[nFiles] = {"GammaCalo_121","GammaCalo_122","GammaCalo_123","GammaCalo_124","GammaCalo_126"};
//        const Int_t nFiles = 2;
//        TString Files[nFiles] = {"GammaCalo_128","GammaCalo_129"};

    //*********************************************************************************************************************************
    //*********************************************************************************************************************************
    //*********************************************************************************************************************************

//        const Int_t nSets = 1;
//        const Int_t nData = 0;
//        TString DataSets[nSets]={
//          "LHC16c2"
//        };
//        TString PrefixDataSets[nSets]={
//          "/alice/sim/2016/LHC16c2/"
//        };

//        TString train = "Legotrain-vAN-20160617-8TeV-Systematics_Calo_EGA";
//        //Int_t trainRuns[nSets] = {2198};
//        //Int_t trainRuns[nSets] = {2199};
//        Int_t trainRuns[nSets] = {2200};
//        TString runlist[nSets] = {"merge"};

//        //const Int_t nFiles = 2;
//        //TString Files[nFiles] = {"GammaCalo_142","GammaCalo_143"};
//        //const Int_t nFiles = 3;
//        //TString Files[nFiles] = {"GammaCalo_141","GammaCalo_144","GammaCalo_146"};
//        const Int_t nFiles = 3;
//        TString Files[nFiles] = {"GammaCalo_147","GammaCalo_148","GammaCalo_149"};

    //*********************************************************************************************************************************
    //*********************************************************************************************************************************
    //*********************************************************************************************************************************

//    const Int_t nSets = 1;
//    const Int_t nData = 0;
//    TString DataSets[nSets]={
//      "LHC16c2"
//    };
//    TString PrefixDataSets[nSets]={
//      "/alice/sim/2016/LHC16c2/"
//    };

//    TString train = "Legotrain-vAN-20160617-8TeV-Systematics_ConvCalo_EMC7";
//    Int_t trainRuns[nSets] = {2205};
//    //Int_t trainRuns[nSets] = {2274};
//    //Int_t trainRuns[nSets] = {2273};
//    //Int_t trainRuns[nSets] = {2272};
//    TString runlist[nSets] = {"merge"};

//    const Int_t nFiles = 3;
//    TString Files[nFiles] = {"GammaConvCalo_148","GammaConvCalo_149","GammaConvCalo_150"};
//    //const Int_t nFiles = 4;
//    //TString Files[nFiles] = {"GammaConvCalo_144","GammaConvCalo_145","GammaConvCalo_146","GammaConvCalo_147"};
//    //const Int_t nFiles = 4;
//    //TString Files[nFiles] = {"GammaConvCalo_136","GammaConvCalo_139","GammaConvCalo_140","GammaConvCalo_143"};
//    //const Int_t nFiles = 4;
//    //TString Files[nFiles] = {"GammaConvCalo_132","GammaConvCalo_133","GammaConvCalo_134","GammaConvCalo_135"};


    //*********************************************************************************************************************************
    //*********************************************************************************************************************************
    //*********************************************************************************************************************************

//    const Int_t nSets = 1;
//    const Int_t nData = 0;
//    TString DataSets[nSets]={
//      "LHC16c2"
//    };
//    TString PrefixDataSets[nSets]={
//      "/alice/sim/2016/LHC16c2/"
//    };

//    TString train = "Legotrain-vAN-20160617-8TeV-Systematics_ConvCalo_EGA";
//    //Int_t trainRuns[nSets] = {2277};
//    //Int_t trainRuns[nSets] = {2276};
//    Int_t trainRuns[nSets] = {2275};
//    //Int_t trainRuns[nSets] = {2206};
//    TString runlist[nSets] = {"merge"};

//    //const Int_t nFiles = 3;
//    //TString Files[nFiles] = {"GammaConvCalo_178","GammaConvCalo_179","GammaConvCalo_180"};
//    //const Int_t nFiles = 4;
//    //TString Files[nFiles] = {"GammaConvCalo_174","GammaConvCalo_175","GammaConvCalo_176","GammaConvCalo_177"};
//    const Int_t nFiles = 4;
//    TString Files[nFiles] = {"GammaConvCalo_166","GammaConvCalo_169","GammaConvCalo_170","GammaConvCalo_173"};
//    //const Int_t nFiles = 4;
//    //TString Files[nFiles] = {"GammaConvCalo_162","GammaConvCalo_163","GammaConvCalo_164","GammaConvCalo_165"};


    //*********************************************************************************************************************************
    //*********************************************************************************************************************************
    //*********************************************************************************************************************************

//    const Int_t nSets = 1;
//    const Int_t nData = 0;
//    TString DataSets[nSets]={
//      "LHC16c2"
//    };
//    TString PrefixDataSets[nSets]={
//      "/alice/sim/2016/LHC16c2/"
//    };

//    TString train = "Legotrain-vAN-20160615-8TeV-JetJet-AOD";
//    Int_t trainRuns[nSets] = {122};
//    TString runlist[nSets] = {"merge"};

//    const Int_t nFiles = 4;
//    TString Files[nFiles] = {"GammaCalo_101","GammaCalo_129","GammaConvCalo_101","GammaConvCalo_140"};

    //*********************************************************************************************************************************
    //*********************************************************************************************************************************
    //*********************************************************************************************************************************

//    const Int_t nSets = 1;
//    const Int_t nData = 0;
//    TString DataSets[nSets]={
//      "LHC16c2"
//    };
//    TString PrefixDataSets[nSets]={
//      "/alice/sim/2016/LHC16c2/"
//    };

//    TString train = "Legotrain-vAN-20161127-8TeV-std_ConvCalo";
//    Int_t trainRuns[nSets] = {2679};
//    TString runlist[nSets] = {"merge"};

//    const Int_t nFiles = 2;
//    TString Files[nFiles] = {"GammaConvCalo_100","GammaConvCalo_101"};

    //*********************************************************************************************************************************
    //*********************************************************************************************************************************
    //*********************************************************************************************************************************

    const Int_t nSets = 1;
    const Int_t nData = 0;
    TString DataSets[nSets]={
      "LHC16c2"
    };
    TString PrefixDataSets[nSets]={
      "/alice/sim/2016/LHC16c2/"
    };

    TString train = "Legotrain-vAN-20161125-8TeV-mergedAndVariations";
    Int_t trainRuns[nSets] = {289};
    TString runlist[nSets] = {"merge"};

    const Int_t nFiles = 3;
    TString Files[nFiles] = {"GammaCalo_123","GammaCaloMerged_107","GammaCaloMerged_171"};

//    const Int_t nSets = 2;
//    const Int_t nData = 0;
//    TString DataSets[nSets]={
//      "LHC16c2",
//      "LHC16c2"
//    };
//    TString PrefixDataSets[nSets]={
//      "/alice/sim/2016/LHC16c2/",
//      "/alice/sim/2016/LHC16c2/"
//    };

//    TString train = "Legotrain-vAN-20161108-8TeV-merged_new_AOD";
//    Int_t trainRuns[nSets] = {260,261};
//    TString runlist[nSets] = {"merge","merge"};

//    const Int_t nFiles = 4;
//    TString Files[nFiles] = {"GammaCaloMerged_109","GammaCaloMerged_110","GammaCaloMerged_116","GammaCaloMerged_117"};

    //---------------------------------------------------------------------------------------------------

    TString alienFolder;
    TString alienFolder_MC;

    if(type.CompareTo("ESD")==0){
      alienFolder = Form("/alice/cern.ch/user/a/alitrain/PWGGA/GA_%s/",system.Data());
      alienFolder_MC = Form("/alice/cern.ch/user/a/alitrain/PWGGA/GA_%s_MC/",system.Data());
    }else{
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
                    fPathGrid = Form("%s%s/%s/PWGGA/GA_pp_MC/%s/%s.root", PrefixDataSets[i].Data(), vecBins.at(b).Data(), vecRuns.at(j).Data(), strTrain[i].Data(), Files[k].Data());
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
              gSystem->Exec(Form("ln -s %s %s/%s_%s.root",fPathTemp.Data(), fPathLocal.Data(), DataSets[i].Data(), Files[k].Data()));
              continue;
            }

            if(copyAlien2Local(fPathGrid,fPathTemp)){
              ChangeStrucToStd(fPathTemp.Data(),fPathTemp.Data(),Files[k].Data());
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
