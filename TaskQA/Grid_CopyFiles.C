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

void Grid_CopyFiles(TString system = "pp", TString type = "ESD", TString folder = "/home/daniel/data/work/Grid")
{
    cout<<"Connecting to Alien..."<<endl;
    TGrid::Connect("alien://");
    cout<<"==============================="<<endl;
    cout<<"Successfully connected to Alien"<<endl;
    cout<<"==============================="<<endl;


    //*********************************************************************************************************************************
    //*********************************************************************************************************************************
    //*********************************************************************************************************************************

//    const Int_t nSets = 10;
//    const Int_t nData = 5;
//    TString DataSets[nSets]={
//      "LHC10b", "LHC10c", "LHC10d", "LHC10e", "LHC10f",
//      "LHC14j4b", "LHC14j4c", "LHC14j4d", "LHC14j4e", "LHC14j4f"
//    };

//    TString train = "Legotrain-vAN-20180828-7TeV-rerun";

//    TString runlist[nSets] = {
//      "_child_1/merge_runlist_2","_child_2/merge_runlist_2","_child_3/merge_runlist_2","_child_4/merge_runlist_2","_child_5/merge_runlist_2",
//      "_child_1/merge_runlist_2","_child_2/merge_runlist_2","_child_3/merge_runlist_2","_child_4/merge_runlist_2","_child_5/merge_runlist_2"
//    };

//    Int_t trainRuns[nSets] = {
//      2480,2480,2480,2480,2480,
//      3494,3494,3494,3494,3494
//    };
//    const Int_t nFiles = 1;
//    TString Files[nFiles] = {
//      "GammaCalo_201"
//    };

//    const Int_t nMerge = 2;
//    TString strMerge[nMerge]={"LHC10","LHC14j4"};
//    std::vector<Int_t> mergeVec[nMerge];
//    std::vector<Int_t>::iterator it;
//    for(Int_t i=0; i<nSets; i++){
//      if(0<=i && i<=4) mergeVec[0].push_back(i);
//      if(5<=i && i<=9) mergeVec[1].push_back(i);
//    }

    //*********************************************************************************************************************************
    //*********************************************************************************************************************************
    //*********************************************************************************************************************************

//    const Int_t nSets = 2;
//    const Int_t nData = 1;
//    TString DataSets[nSets]={
//      "LHC10c_900GeV",
//      "LHC14j4c_900GeV"
//    };

//    TString train = "Legotrain-vAN-20180828-900GeV-dirGamma_rerun";

//    TString runlist[nSets] = {
//      "merge",
//      "merge"
//    };

//    Int_t trainRuns[nSets] = {
//      2479,
//      3493
//    };
//    const Int_t nFiles = 6;
//    TString Files[nFiles] = {
//      "GammaCalo_221","GammaCalo_222","GammaCalo_281","GammaConvCalo_200","GammaConvCalo_201","GammaConvCalo_281"
//    };

//    const Int_t nMerge = 2;
//    TString strMerge[nMerge]={"LHC10_900GeV","LHC14j4_900GeV"};
//    std::vector<Int_t> mergeVec[nMerge];
//    std::vector<Int_t>::iterator it;
//    for(Int_t i=0; i<nSets; i++){
//      if(i==0) mergeVec[0].push_back(i);
//      if(i==1) mergeVec[1].push_back(i);
//    }

//*********************************************************************************************************************************
//*********************************************************************************************************************************
//*********************************************************************************************************************************

//    const Int_t nSets = 21;
//    const Int_t nData = 7;
//    TString DataSets[nSets]={
//      "LHC12a", "LHC12b", "LHC12c", "LHC12d", "LHC12f", "LHC12h", "LHC12i",
//      "LHC15h1a1", "LHC15h1b", "LHC15h1c", "LHC15h1d", "LHC15h1f", "LHC15h1h", "LHC15h1i",
//      "LHC15h2a", "LHC15h2b", "LHC15h2c", "LHC15h2d", "LHC15h2f", "LHC15h2h", "LHC15h2i"
//    };

//    TString train = "Legotrain-vAN-20180828-8TeV-TM_valid";
//    Int_t trainRuns[nSets] = {
//                              2490,2490,2490,2490,2490,2490,2490,
//                              3520,3520,3520,3520,3520,3520,3520,
//                              3521,3521,3521,3521,3521,3521,3521
//                             };

//    TString runlist[nSets] = {
//      "_child_1/merge_runlist_2","_child_2/merge_runlist_2","_child_3/merge_runlist_2","_child_4/merge_runlist_2","_child_5/merge_runlist_2","_child_6/merge_runlist_2","_child_7/merge_runlist_2",
//      "_child_1/merge_runlist_2","_child_2/merge_runlist_2","_child_3/merge_runlist_2","_child_4/merge_runlist_2","_child_5/merge_runlist_2","_child_6/merge_runlist_2","_child_7/merge_runlist_2",
//      "_child_1/merge_runlist_2","_child_2/merge_runlist_2","_child_3/merge_runlist_2","_child_4/merge_runlist_2","_child_5/merge_runlist_2","_child_6/merge_runlist_2","_child_7/merge_runlist_2"
//    };

//    const Int_t nFiles = 4;
//    TString Files[nFiles] = {
//      "GammaCalo_181", "GammaCalo_182", "GammaCalo_183", "GammaCalo_184"
//                            };

//    const Int_t nMerge = 4;
//    TString strMerge[nMerge]={"LHC12","LHC15h1","LHC15h2","LHC15h"};
//    std::vector<Int_t> mergeVec[nMerge];
//    std::vector<Int_t>::iterator it;
//    for(Int_t i=0; i<nSets; i++){
//      if(0<=i && i<=6) mergeVec[0].push_back(i);
//      if(7<=i && i<=13) mergeVec[1].push_back(i);
//      if(14<=i && i<=20) mergeVec[2].push_back(i);
//      if(7<=i && i<=20) mergeVec[3].push_back(i);
//    }

    //*********************************************************************************************************************************
    //*********************************************************************************************************************************
    //*********************************************************************************************************************************

    const Int_t nSets = 14;
    const Int_t nData = 0;
    TString DataSets[nSets]={
      "LHC15h1a1", "LHC15h1b", "LHC15h1c", "LHC15h1d", "LHC15h1f", "LHC15h1h", "LHC15h1i",
      "LHC15h2a", "LHC15h2b", "LHC15h2c", "LHC15h2d", "LHC15h2f", "LHC15h2h", "LHC15h2i"
    };

    TString train = "Legotrain-vAN-20180828-8TeV-dirGamma_rerun_20180925";
    Int_t trainRuns[nSets] = {
      3532,3532,3532,3532,3532,3532,3532,
      3533,3533,3533,3533,3533,3533,3533
    };


    TString runlist[nSets] = {
      "_child_1/merge_runlist_2","_child_2/merge_runlist_2","_child_3/merge_runlist_2","_child_4/merge_runlist_2","_child_5/merge_runlist_2","_child_6/merge_runlist_2","_child_7/merge_runlist_2",
      "_child_1/merge_runlist_2","_child_2/merge_runlist_2","_child_3/merge_runlist_2","_child_4/merge_runlist_2","_child_5/merge_runlist_2","_child_6/merge_runlist_2","_child_7/merge_runlist_2"
    };

    const Int_t nFiles = 3;
    TString Files[nFiles] = {
      "GammaCalo_119", "GammaCalo_120", "GammaConvCalo_131"
    };

    const Int_t nMerge = 3;
    TString strMerge[nMerge]={"LHC15h1","LHC15h2","LHC15h"};
    std::vector<Int_t> mergeVec[nMerge];
    std::vector<Int_t>::iterator it;
    for(Int_t i=0; i<nSets; i++){
      if(0<=i && i<=6) mergeVec[0].push_back(i);
      if(7<=i && i<=13) mergeVec[1].push_back(i);
      if(0<=i && i<=13) mergeVec[2].push_back(i);
    }


//*********************************************************************************************************************************
//*********************************************************************************************************************************
//*********************************************************************************************************************************

//    const Int_t nSets = 18;
//    const Int_t nData = 18;
//    TString DataSets[nSets]={
//      "LHC12a", "LHC12b", "LHC12c", "LHC12d", "LHC12f", "LHC12h", "LHC12i",
//      "LHC12b-kEMC7", "LHC12c-kEMC7", "LHC12d-kEMC7", "LHC12f-kEMC7", "LHC12h-kEMC7", "LHC12i-kEMC7",
//      "LHC12c-kEMCEGA", "LHC12d-kEMCEGA", "LHC12f-kEMCEGA", "LHC12h-kEMCEGA", "LHC12i-kEMCEGA"
//    };

//    TString train = "Legotrain-vAN-20171016-8TeV-dirGamma_EMC";
//    Int_t trainRuns[nSets] = {
//                              2258,2258,2258,2258,2258,2258,2258,
//                              2258,2258,2258,2258,2258,2258,
//                              2258,2258,2258,2258,2258
//                             };

//    TString runlist[nSets] = {
//                              "_child_1/merge_runlist_2","_child_2/merge_runlist_2","_child_3/merge_runlist_2","_child_4/merge_runlist_2","_child_5/merge_runlist_2","_child_6/merge_runlist_2","_child_7/merge_runlist_2",
//                              "_child_2/merge_runlist_3","_child_3/merge_runlist_3","_child_4/merge_runlist_3","_child_5/merge_runlist_3","_child_6/merge_runlist_3","_child_7/merge_runlist_3",
//                              "_child_3/merge_runlist_4","_child_4/merge_runlist_4","_child_5/merge_runlist_4","_child_6/merge_runlist_4","_child_7/merge_runlist_4"
//                             };

//    const Int_t nFiles = 3;
//    TString Files[nFiles] = {
//      "GammaCalo_103", "GammaCalo_161", "GammaConvCalo_101"
//                            };

//    const Int_t nMerge = 3;
//    TString strMerge[nMerge]={"LHC12", "LHC12-kEMC7", "LHC12-kEMCEGA"};
//    std::vector<Int_t> mergeVec[nMerge];
//    std::vector<Int_t>::iterator it;
//    for(Int_t i=0; i<nSets; i++){
//      if(0<=i && i<=6) mergeVec[0].push_back(i);
//      if(7<=i && i<=12) mergeVec[1].push_back(i);
//      if(13<=i && i<=17) mergeVec[2].push_back(i);
//    }

//*********************************************************************************************************************************
//*********************************************************************************************************************************
//*********************************************************************************************************************************

//systematics DirGamma 8 TeV
//    const Int_t nSets = 7;
//    const Int_t nData = 7;
//    TString DataSets[nSets]={
//      "LHC12a", "LHC12b", "LHC12c", "LHC12d", "LHC12f", "LHC12h", "LHC12i"
//    };

//    TString train = "Legotrain-vAN-20171002-8TeV-dirGamma_Systematics";
//    Int_t trainRuns[nSets] = {
//                              2246,2246,2246,2246,2246,2246,2246
//                             };

//    TString runlist[nSets] = {
//                              "_child_1/merge_runlist_2","_child_2/merge_runlist_2","_child_3/merge_runlist_2","_child_4/merge_runlist_2","_child_5/merge_runlist_2","_child_6/merge_runlist_2","_child_7/merge_runlist_2"
//                             };

//    const Int_t nFiles = 20;
//    TString Files[nFiles] = {
//      "GammaCalo_102", "GammaCalo_103", "GammaCalo_104", "GammaCalo_105",
//      "GammaCalo_107", "GammaCalo_109", "GammaCalo_110",
//      "GammaConvCalo_102", "GammaConvCalo_103", "GammaConvCalo_104", "GammaConvCalo_105",
//      "GammaConvCalo_106", "GammaConvCalo_109", "GammaConvCalo_110", "GammaConvCalo_113",
//      "GammaConvCalo_115", "GammaConvCalo_116", "GammaConvCalo_117", "GammaConvCalo_118", "GammaConvCalo_119"
//                            };

//    const Int_t nMerge = 1;
//    TString strMerge[nMerge]={"LHC12"};
//    std::vector<Int_t> mergeVec[nMerge];
//    std::vector<Int_t>::iterator it;
//    for(Int_t i=0; i<nSets; i++){
//      if(0<=i && i<=6) mergeVec[0].push_back(i);
//    }

//    const Int_t nSets = 14;
//    const Int_t nData = 0;
//    TString DataSets[nSets]={
//      "LHC15h1a1", "LHC15h1b", "LHC15h1c", "LHC15h1d", "LHC15h1f", "LHC15h1h", "LHC15h1i",
//      "LHC15h2a", "LHC15h2b", "LHC15h2c", "LHC15h2d", "LHC15h2f", "LHC15h2h", "LHC15h2i"
//    };

//    TString train = "Legotrain-vAN-20171002-8TeV-dirGamma_Systematics";
//    Int_t trainRuns[nSets] = {
////      3147,3147,3147,3147,3147,3147,3147,
////      3148,3148,3148,3148,3148,3148,3148

////      3145,3145,3145,3145,3145,3145,3145,
////      3146,3146,3146,3146,3146,3146,3146

////      3143,3143,3143,3143,3143,3143,3143,
////      3144,3144,3144,3144,3144,3144,3144

////      3141,3141,3141,3141,3141,3141,3141,
////      3142,3142,3142,3142,3142,3142,3142

//      3139,3139,3139,3139,3139,3139,3139,
//      3140,3140,3140,3140,3140,3140,3140
//    };

//    TString runlist[nSets] = {
//      "_child_1/merge_runlist_2","_child_2/merge_runlist_2","_child_3/merge_runlist_2","_child_4/merge_runlist_2","_child_5/merge_runlist_2","_child_6/merge_runlist_2","_child_7/merge_runlist_2",
//      "_child_1/merge_runlist_2","_child_2/merge_runlist_2","_child_3/merge_runlist_2","_child_4/merge_runlist_2","_child_5/merge_runlist_2","_child_6/merge_runlist_2","_child_7/merge_runlist_2"
//    };

////    const Int_t nFiles = 4;

////    const Int_t nFiles = 3;

////    const Int_t nFiles = 4;

////    const Int_t nFiles = 5;

//    const Int_t nFiles = 4;

//    TString Files[nFiles] = {
////      "GammaCalo_105", "GammaCalo_107", "GammaCalo_109", "GammaCalo_110"

////      "GammaCalo_102", "GammaCalo_103", "GammaCalo_104"

////      "GammaConvCalo_116", "GammaConvCalo_117", "GammaConvCalo_118", "GammaConvCalo_119"

////      "GammaConvCalo_106", "GammaConvCalo_109", "GammaConvCalo_110", "GammaConvCalo_113", "GammaConvCalo_115"

//        "GammaConvCalo_102", "GammaConvCalo_103", "GammaConvCalo_104", "GammaConvCalo_105"
//    };

//    const Int_t nMerge = 0;
//    TString strMerge[nMerge]={ /*"LHC15h1","LHC15h2"*/};
//    std::vector<Int_t> mergeVec[nMerge];
//    std::vector<Int_t>::iterator it;
//    for(Int_t i=0; i<nSets; i++){
//      //if(0<=i && i<=13) mergeVec[0].push_back(i);
////      if(0<=i && i<=6) mergeVec[0].push_back(i);
////      if(7<=i && i<=13) mergeVec[1].push_back(i);
//    }

//*********************************************************************************************************************************
//*********************************************************************************************************************************
//*********************************************************************************************************************************

//        const Int_t nSets = 18;
//        const Int_t nData = 18;
//        TString DataSets[nSets]={
//          "LHC12a", "LHC12b", "LHC12c", "LHC12d", "LHC12f", "LHC12h", "LHC12i",
//          "LHC12b-kEMC7", "LHC12c-kEMC7", "LHC12d-kEMC7", "LHC12f-kEMC7", "LHC12h-kEMC7", "LHC12i-kEMC7",
//          "LHC12c-kEMCEGA", "LHC12d-kEMCEGA", "LHC12f-kEMCEGA", "LHC12h-kEMCEGA", "LHC12i-kEMCEGA"
//        };

//        TString train = "Legotrain-vAN-20180123-8TeV-dirGamma-JetJet";
//        Int_t trainRuns[nSets] = {
//                                  2300,2300,2300,2300,2300,2300,2300,
//                                  2300,2300,2300,2300,2300,2300,
//                                  2300,2300,2300,2300,2300
//                                 };

//        TString runlist[nSets] = {
//                                  "_child_1/merge_runlist_2","_child_2/merge_runlist_2","_child_3/merge_runlist_2","_child_4/merge_runlist_2","_child_5/merge_runlist_2","_child_6/merge_runlist_2","_child_7/merge_runlist_2",
//                                  "_child_2/merge_runlist_3","_child_3/merge_runlist_3","_child_4/merge_runlist_3","_child_5/merge_runlist_3","_child_6/merge_runlist_3","_child_7/merge_runlist_3",
//                                  "_child_3/merge_runlist_4","_child_4/merge_runlist_4","_child_5/merge_runlist_4","_child_6/merge_runlist_4","_child_7/merge_runlist_4"
//                                 };

//        const Int_t nFiles = 6;
//        TString Files[nFiles] = {
//          "GammaCalo_120", "GammaCalo_140", "GammaCalo_160",
//          "GammaConvCalo_130", "GammaConvCalo_159", "GammaConvCalo_181"
//                                };

//        const Int_t nMerge = 3;
//        TString strMerge[nMerge]={"LHC12", "LHC12-kEMC7", "LHC12-kEMCEGA"};
//        std::vector<Int_t> mergeVec[nMerge];
//        std::vector<Int_t>::iterator it;
//        for(Int_t i=0; i<nSets; i++){
//          if(0<=i && i<=6) mergeVec[0].push_back(i);
//          if(7<=i && i<=12) mergeVec[1].push_back(i);
//          if(13<=i && i<=17) mergeVec[2].push_back(i);
//        }

//*********************************************************************************************************************************
//*********************************************************************************************************************************
//*********************************************************************************************************************************

//    const Int_t nSets = 14;
//    const Int_t nData = 0;
//    TString DataSets[nSets]={
//      "LHC15h1a1", "LHC15h1b", "LHC15h1c", "LHC15h1d", "LHC15h1f", "LHC15h1h", "LHC15h1i",
//      "LHC15h2a", "LHC15h2b", "LHC15h2c", "LHC15h2d", "LHC15h2f", "LHC15h2h", "LHC15h2i"
//    };

//    TString train = "Legotrain-vAN-20171213-8TeV-dirGamma_EMC_SysVar";
//    Int_t trainRuns[nSets] = {
//      3205,3205,3205,3205,3205,3205,3205,
//      3206,3206,3206,3206,3206,3206,3206
//    };

//    TString runlist[nSets] = {
//      "_child_1/merge_runlist_2","_child_2/merge_runlist_2","_child_3/merge_runlist_2","_child_4/merge_runlist_2","_child_5/merge_runlist_2","_child_6/merge_runlist_2","_child_7/merge_runlist_2",
//      "_child_1/merge_runlist_2","_child_2/merge_runlist_2","_child_3/merge_runlist_2","_child_4/merge_runlist_2","_child_5/merge_runlist_2","_child_6/merge_runlist_2","_child_7/merge_runlist_2"
//    };

//    const Int_t nFiles = 3;
//    TString Files[nFiles] = {
//      "GammaCalo_165","GammaCalo_166","GammaCalo_167"
//    };

//    const Int_t nMerge = 3;
//    TString strMerge[nMerge]={ "LHC15h","LHC15h1","LHC15h2"};
//    std::vector<Int_t> mergeVec[nMerge];
//    std::vector<Int_t>::iterator it;
//    for(Int_t i=0; i<nSets; i++){
//      if(0<=i && i<=13) mergeVec[0].push_back(i);
//      if(0<=i && i<=6) mergeVec[1].push_back(i);
//      if(7<=i && i<=13) mergeVec[2].push_back(i);
//    }

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
            if(temp.Contains("_child_")) temp.Remove(temp.Length()-8,8);
            strTrain[i] = temp;
            break;
          }
        }
      }else{
        for(Int_t j=0; j<(Int_t)vecStrTrainMC.size(); j++){
          tempRuns = Form("%i",trainRuns[i]);
          temp = vecStrTrainMC.at(j);
          if(temp.BeginsWith(tempRuns)){
            if(temp.Contains("_child_")) temp.Remove(temp.Length()-8,8);
            strTrain[i] = temp;
            break;
          }
        }
      }
    }

    Int_t nErr[nSets];
    TString fPathGrid;
    TString fPathLocal;

    fPathLocal = Form("%s/%s", folder.Data(), train.Data());
    gSystem->Exec(Form("mkdir -p %s",fPathLocal.Data()));

    TString mergePeriod[nMerge][nFiles];
    for(Int_t iFiles=0; iFiles<nFiles; iFiles++){
      for(Int_t iM=0; iM<nMerge; iM++){
        TString mergeP = Form("%s/%s_%s.root", fPathLocal.Data(), strMerge[iM].Data(), Files[iFiles].Data());
        mergePeriod[iM][iFiles] = Form("hadd -f -k %s",mergeP.Data());
      }
    }

    for(Int_t i=0; i<nSets; i++)
    {
        nErr[i]=0;
        for(Int_t k=0; k<nFiles; k++)
        {
          if(Files[k].BeginsWith("OmegaToPiZeroGamma_")){
            if(i<nData) Files[k].Append("_0");
            else Files[k].Append("_1");
          }
          if(runlist[i].Contains("/")){
            if(i<nData) fPathGrid = Form("%s%s%s/%s.root", alienFolder.Data(), strTrain[i].Data(), runlist[i].Data(), Files[k].Data());
            else fPathGrid = Form("%s%s%s/%s.root", alienFolder_MC.Data(), strTrain[i].Data(), runlist[i].Data(), Files[k].Data());
          }else{
            if(i<nData) fPathGrid = Form("%s%s/%s/%s.root", alienFolder.Data(), strTrain[i].Data(), runlist[i].Data(), Files[k].Data());
            else fPathGrid = Form("%s%s/%s/%s.root", alienFolder_MC.Data(), strTrain[i].Data(), runlist[i].Data(), Files[k].Data());
          }

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
              if( it!=mergeVec[iM].end() ) mergePeriod[iM][k] += Form(" %s",fPathTemp.Data());
            }
            gSystem->Exec(Form("ln -s %s %s/%s_%s.root",fPathTemp.Data(), fPathLocal.Data(), DataSets[i].Data(), Files[k].Data()));
            if(Files[k].BeginsWith("OmegaToPiZeroGamma_")){Files[k].Resize(Files[k].Length()-2);}
            continue;
          }

          if(copyAlien2Local(fPathGrid,fPathTemp)){
            ChangeStrucToStd(fPathTemp.Data(),fPathTemp.Data(),Files[k].Data());
            for(Int_t iM=0; iM<nMerge; iM++){
              it = find( mergeVec[iM].begin(), mergeVec[iM].end(), i);
              if( it!=mergeVec[iM].end() ) mergePeriod[iM][k] += Form(" %s",fPathTemp.Data());
            }
            gSystem->Exec(Form("ln -s %s %s/%s_%s.root",fPathTemp.Data(), fPathLocal.Data(), DataSets[i].Data(), Files[k].Data()));
            if(Files[k].BeginsWith("OmegaToPiZeroGamma_")){Files[k].Resize(Files[k].Length()-2);}
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

    cout << "\n------------------------------------------------------" << endl;
    cout << "Merging: " << endl;
    cout << "------------------------------------------------------\n" << endl;

    for(Int_t iFiles=0; iFiles<nFiles; iFiles++){
      for(Int_t iM=0; iM<nMerge; iM++){
        cout << "Merging " << strMerge[iM].Data() << " ..." << endl;
        gSystem->Exec(mergePeriod[iM][iFiles].Data());
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
