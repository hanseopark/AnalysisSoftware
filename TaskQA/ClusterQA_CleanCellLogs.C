#include "QA.h"
#include <vector>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <fstream>
#include <TObject.h>
#include <TString.h>
#include <TFile.h>
#include <TSystem.h>
#include <iostream>
#include <fstream>
#include "Riostream.h"


/*
 * This macro will go through log files in ClusterQA_HotCellCompare and ClusterQA_DeadCellCompare and check if some of the cells in there were already
 * flagged as good by the user. The cleaned logs will be written to $ORIGINALFILE-Cleaned.log
 * The information for which folders were flagged as good(+maybe) by user are taken from folders containing the Cell*_EnergyComparison files.
 *
 * Usage (add to config used by QAV2.C)
 *
 * cellCleaningUseMaybe 0 : run the cleaning of log files but ignore a maybe folder
 * cellCleaningUseMaybe 1 : consider all cells in maybe folder as good aswell
 * cellCleaningUseMaybe 2 : consider all cells in maybe folder as bad aswell
 *
 * Optional:
 * userGoodCellDirName $PATH
 * userBadCellDirName $PATH
 * userMaybeCellDirName $PATH
 *
 * if the user doesn't give them explicitly, macro will search in
 *
 * $CUTNUMBER/$ENERGY/ClusterQA/$SUFFIX/Cells/Detailed/ok
 * $CUTNUMBER/$ENERGY/ClusterQA/$SUFFIX/Cells/Detailed/maybe
 * $CUTNUMBER/$ENERGY/ClusterQA/$SUFFIX/Cells/Detailed/bad
 */

void read_directory(const TString& name, vector<TString>& v){

    TString list        = gSystem->GetFromPipe(Form("ls -w 1 %s", name.Data()));
    TObjArray *tempArr  = list.Tokenize("\n");

    for (Int_t i = 0; i < tempArr->GetEntries(); i++){
        TString current = (TString)((TObjString*)tempArr->At(i))->GetString();
        if (current.CompareTo("") != 0 && current.CompareTo(" ") != 0 )
            v.push_back(current);
    }
    delete tempArr;
}

void ClusterQA_CleanCellLogs(    TString configFileName  = "configFile.txt",
                                 TString suffix          = "eps" ){
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "Running ClusterQA_CleanCellLogs" << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

    //**************************************************************************************************************
    //******************* global settings **************************************************************************
    //**************************************************************************************************************
    Int_t useMaybe     = 0;
    const Int_t maxSets = 1 ;

    TString inputDirHotCellCompare = "ClusterQA_HotCellCompare";
    TString inputDirDeadCellCompare = "ClusterQA_DeadCellCompare";
    TString inputDirUserGood = "";
    TString inputDirUserMaybe = "";
    TString inputDirUserBad = "";


    TString dataSetsCold[maxSets];
    TString dataSetsHot[maxSets];
    TString addInputDirNameHotCells         = "";
    TString addInputDirNameDeadCells        = "";
    TString hotCellDataCuts                 = "";
    TString deadCellDataCuts                = "";
    TString energy                          = "";

    vector<Int_t> CellsThatNeedCheck;

    // initialize arrays
    for (Int_t i = 0; i< maxSets; i++){
        dataSetsCold[i]                 = "";
        dataSetsHot[i]                 = "";
    }
    //**************************************************************************************************************
    //******************************* Read config file for detailed settings ***************************************
    //**************************************************************************************************************
    // ATTENTION: The data set has to be separated with either tabs or spaces a mixture of
    //            both will most likely lead to misconfigurations
    //**************************************************************************************************************

    cout << "INFO: You have chosen the given the following config file: " << configFileName.Data() << endl;
    ifstream fileConfigQA;
    fileConfigQA.open(configFileName,ios_base::in);
    if (!fileConfigQA) {
        cout << "ERROR: settings " << configFileName.Data() << " not found!" << endl;
        return;
    }

    // read settings from file
    for( TString tempLine; tempLine.ReadLine(fileConfigQA, kTRUE); ){
        // check if line should be considered
        if (tempLine.BeginsWith("%") || tempLine.BeginsWith("#")){
            continue;
        }
//         cout << tempLine.Data() << endl;

        // Separate the string according to tabulators
        TObjArray *tempArr  = tempLine.Tokenize("\t");
        if(tempArr->GetEntries()<1){
            cout << "nothing to be done" << endl;
            delete tempArr;
            continue;
        } else if (tempArr->GetEntries() == 1 ){
            // Separate the string according to space
            tempArr       = tempLine.Tokenize(" ");
            if(tempArr->GetEntries()<1){
                cout << "nothing to be done" << endl;
                delete tempArr;
                continue;
            } else if (tempArr->GetEntries() == 1 ) {
                cout << ((TString)((TObjString*)tempArr->At(0))->GetString()).Data() << " has not be reset, no value given!" << endl;
                delete tempArr;
                continue;
            }
        }

        // Put them to the correct variables
        TString tempValue   = (TString)((TObjString*)tempArr->At(0))->GetString();
        if (tempValue.BeginsWith("hotCellAdditionalOutputDirName",TString::kIgnoreCase)){
            addInputDirNameHotCells= (TString)((TObjString*)tempArr->At(1))->GetString();
        } else if (tempValue.BeginsWith("deadCellAdditionalOutputDirName",TString::kIgnoreCase)){
           addInputDirNameDeadCells= (TString)((TObjString*)tempArr->At(1))->GetString();
        } else if (tempValue.BeginsWith("hotCellDataCuts",TString::kIgnoreCase)){
           hotCellDataCuts= (TString)((TObjString*)tempArr->At(1))->GetString();
        } else if (tempValue.BeginsWith("deadCellDataCuts",TString::kIgnoreCase)){
           deadCellDataCuts= (TString)((TObjString*)tempArr->At(1))->GetString();
        } else if (tempValue.BeginsWith("energy",TString::kIgnoreCase)){
           energy= (TString)((TObjString*)tempArr->At(1))->GetString();
        } else if (tempValue.BeginsWith("userGoodCellDirName",TString::kIgnoreCase)){
           inputDirUserGood= (TString)((TObjString*)tempArr->At(1))->GetString();
        } else if (tempValue.BeginsWith("userMaybeCellDirName",TString::kIgnoreCase)){
           inputDirUserMaybe= (TString)((TObjString*)tempArr->At(1))->GetString();
        } else if (tempValue.BeginsWith("userBadCellDirName",TString::kIgnoreCase)){
           inputDirUserBad= (TString)((TObjString*)tempArr->At(1))->GetString();
        } else if (tempValue.BeginsWith("cellCleaningUseMaybe",TString::kIgnoreCase)){
            useMaybe = ((TString)((TObjString*)tempArr->At(1))->GetString()).Atoi();
        } else if (tempValue.BeginsWith("hotCellDataSetNames",TString::kIgnoreCase)){
            for(Int_t i = 1; i<tempArr->GetEntries() && i <= maxSets ; i++){
                if (((TString)((TObjString*)tempArr->At(i))->GetString()).CompareTo("stop",TString::kIgnoreCase))
                    dataSetsHot[i-1]       = ((TString)((TObjString*)tempArr->At(i))->GetString());
                else
                    i                   = tempArr->GetEntries();
            }
        } else if (tempValue.BeginsWith("deadCellDataSetNames",TString::kIgnoreCase)){
            for(Int_t i = 1; i<tempArr->GetEntries() && i <= maxSets ; i++){
                if (((TString)((TObjString*)tempArr->At(i))->GetString()).CompareTo("stop",TString::kIgnoreCase))
                    dataSetsCold[i-1]       = ((TString)((TObjString*)tempArr->At(i))->GetString());
                else
                    i                   = tempArr->GetEntries();
            }
        }
        delete tempArr;
    }

    inputDirDeadCellCompare=inputDirDeadCellCompare+"/"+addInputDirNameDeadCells;
    inputDirHotCellCompare=inputDirHotCellCompare+"/"+addInputDirNameHotCells;

    // If no input dirs for good, maybe bad folder were set by user, try building your own.
    if(!inputDirUserBad.CompareTo("")){
        cout << "INFO: inputDirUserBad was not found. Will try to use default location ..." << endl;
        inputDirUserBad = Form("%s/%s/ClusterQA/%s/Cells/Detailed/bad",hotCellDataCuts.Data(),energy.Data(),suffix.Data());
    }
    if(!inputDirUserGood.CompareTo("")){
        cout << "INFO: inputDirUserGood was not found. Will try to use default location ..." << endl;
        inputDirUserGood = Form("%s/%s/ClusterQA/%s/Cells/Detailed/ok",hotCellDataCuts.Data(),energy.Data(),suffix.Data());
    }
    if((!inputDirUserMaybe.CompareTo(""))&&(useMaybe>0)){
        cout << "INFO: inputDirUserMaybe was not found. Will try to use default location ..." << endl;
        inputDirUserMaybe = Form("%s/%s/ClusterQA/%s/Cells/Detailed/maybe",hotCellDataCuts.Data(),energy.Data(),suffix.Data());
    }

    //**************************************************************************************************************
    //******************************* Check wether settings were valid *********************************************
    //**************************************************************************************************************

    cout << "**************************************************************************" << endl;
    cout << "**************** Settings found in config file ***************************" << endl;
    cout << "**************************************************************************" << endl;
    cout << "userGoodCellDirName:           " << inputDirUserGood.Data()         << endl;
    if(useMaybe>0){
        cout << "userMaybeCellDirName:           " << inputDirUserMaybe.Data()         << endl;
        if(useMaybe==1) cout<<"Cells in maybe folder will be considered as good" << endl;
        if(useMaybe==2) cout<<"Cells in maybe folder will be considered as bad" << endl;
    } else{
        cout << "using maybe folder disabled" << endl;
    }
    cout << "userBadCellDirName:           " << inputDirUserBad.Data()         << endl;
    cout << "**************************************************************************" << endl;
    cout << "**************** Log Files that will be used******************************" << endl;
    cout << "**************************************************************************" << endl;
    cout << "Input : " << Form("%s/%s.log",inputDirHotCellCompare.Data(),dataSetsHot[0].Data()) << endl;
    cout << "Input : " << Form("%s/%s-Detailed.log",inputDirHotCellCompare.Data(),dataSetsHot[0].Data()) << endl;
    cout << "Input : " << Form("%s/%s-SortedByRun.log",inputDirHotCellCompare.Data(),dataSetsHot[0].Data()) << endl;
    cout << "Input : " << Form("%s/%s-Runwise.log",inputDirHotCellCompare.Data(),dataSetsHot[0].Data()) << endl;
    cout << "Input : " << Form("%s/%s.log",inputDirDeadCellCompare.Data(),dataSetsCold[0].Data()) << endl;
    cout << "Input : " << Form("%s/%s-Detailed.log",inputDirDeadCellCompare.Data(),dataSetsCold[0].Data()) << endl;
    cout << "Input : " << Form("%s/%s-Final.log",inputDirDeadCellCompare.Data(),dataSetsCold[0].Data()) << endl;
    cout << "Input : " << Form("%s/%s-Runwise.log",inputDirDeadCellCompare.Data(),dataSetsCold[0].Data()) << endl;
    cout << "Output: " << Form("%s/%s-Cleaned.log",inputDirHotCellCompare.Data(),dataSetsHot[0].Data()) << endl;
    cout << "Output: " << Form("%s/%s-Detailed-Cleaned.log",inputDirHotCellCompare.Data(),dataSetsHot[0].Data()) << endl;
    cout << "Output: " << Form("%s/%s-SortedByRun-Cleaned.log",inputDirHotCellCompare.Data(),dataSetsHot[0].Data()) << endl;
    cout << "Output: " << Form("%s/%s-Runwise-Cleaned.log",inputDirHotCellCompare.Data(),dataSetsHot[0].Data()) << endl;
    cout << "Output: " << Form("%s/%s-Cleaned.log",inputDirDeadCellCompare.Data(),dataSetsCold[0].Data()) << endl;
    cout << "Output: " << Form("%s/%s-Detailed-Cleaned.log",inputDirDeadCellCompare.Data(),dataSetsCold[0].Data()) << endl;
    cout << "Output: " << Form("%s/%s-Final-Cleaned.log",inputDirDeadCellCompare.Data(),dataSetsCold[0].Data()) << endl;
    cout << "Output: " << Form("%s/%s-Runwise-Cleaned.log",inputDirDeadCellCompare.Data(),dataSetsCold[0].Data()) << endl;

    //**************************************************************************************************************
    //******************************* Get cells flagged by user ****************************************************
    //**************************************************************************************************************

    // Get Cells flagged as maybe by user (if given)
    vector<TString> filenamesMaybe;
    vector<Int_t> CellNumbersOwnMaybe;
    if(useMaybe>0){
        read_directory(inputDirUserMaybe,filenamesMaybe);
        for(UInt_t i =0;i<filenamesMaybe.size();i++){
            if((filenamesMaybe.at(i).CompareTo(""))&&(filenamesMaybe.at(i).CompareTo("."))&&(filenamesMaybe.at(i).CompareTo(".."))){
                TObjArray* fileTokenMaybe = filenamesMaybe.at(i).Tokenize('_');
                TString cellNmbMaybe  = ((TObjString*)fileTokenMaybe->At(0))->GetString();
                cellNmbMaybe.Remove(0,4);
                CellNumbersOwnMaybe.push_back(cellNmbMaybe.Atoi());
            }
        }
        sort(CellNumbersOwnMaybe.begin(),CellNumbersOwnMaybe.end());
        cout << "-----------------------------------------------------------------------------------" <<endl;
        cout << "Cells flagged as maybe by user:" <<endl;
        cout << "-----------------------------------------------------------------------------------" <<endl;
        for(UInt_t i =0;i<CellNumbersOwnMaybe.size();i++){
            cout <<CellNumbersOwnMaybe.at(i) <<", ";
        }
        cout << endl;

    }

    // Get Cells flagged as bad by user
    vector<TString> filenamesBad;
    vector<Int_t> CellNumbersOwnBad;
    read_directory(inputDirUserBad,filenamesBad);

    for(UInt_t i =0;i<filenamesBad.size();i++){
        if((filenamesBad.at(i).CompareTo(""))&&(filenamesBad.at(i).CompareTo("."))&&(filenamesBad.at(i).CompareTo(".."))){
            TObjArray* fileTokenBad = filenamesBad.at(i).Tokenize('_');
            TString cellNmbBad  = ((TObjString*)fileTokenBad->At(0))->GetString();
            cellNmbBad.Remove(0,4);
            CellNumbersOwnBad.push_back(cellNmbBad.Atoi());
        }
    }
    if(useMaybe==2){
        cout << "-----------------------------------------------------------------------------------" <<endl;
        cout << "ATTENTION: Because the maybe option 1 was given, runs in maybe will be flagged as bad from now ...." <<endl;
        cout << "-----------------------------------------------------------------------------------" <<endl;
        for(UInt_t i=0; i<CellNumbersOwnMaybe.size();i++){
            CellNumbersOwnBad.push_back(CellNumbersOwnMaybe[i]);
        }
    }
    sort(CellNumbersOwnBad.begin(),CellNumbersOwnBad.end());
    if(useMaybe==2){
        cout << "-----------------------------------------------------------------------------------" <<endl;
        cout << "Cells flagged as bad or maybe by user:" <<endl;
        cout << "-----------------------------------------------------------------------------------" <<endl;
    } else{
        cout << "-----------------------------------------------------------------------------------" <<endl;
        cout << "Cells flagged as bad by user:" <<endl;
        cout << "-----------------------------------------------------------------------------------" <<endl;
    }
    for(UInt_t i =0;i<CellNumbersOwnBad.size();i++){
        cout <<CellNumbersOwnBad.at(i) << ", ";
    }
    cout << endl;

    // Get Cells flagged as good by user
    vector<TString> filenamesGood;
    vector<Int_t> CellNumbersOwnGood;
    read_directory(inputDirUserGood,filenamesGood);

    for(UInt_t i =0;i<filenamesGood.size();i++){
        if((filenamesGood.at(i).CompareTo(""))&&(filenamesGood.at(i).CompareTo("."))&&(filenamesGood.at(i).CompareTo(".."))){
            TObjArray* fileTokenGood = filenamesGood.at(i).Tokenize('_');
            TString cellNmbGood  = ((TObjString*)fileTokenGood->At(0))->GetString();
            cellNmbGood.Remove(0,4);
            CellNumbersOwnGood.push_back(cellNmbGood.Atoi());
        }
    }
    // Push all the maybe flagged ones to good, in case user gave maybes
    if(useMaybe==1){
        cout << "-----------------------------------------------------------------------------------" <<endl;
        cout << "ATTENTION: Because the maybe option 1 was given, runs in maybe will be flagged as good from now ...." <<endl;
        cout << "-----------------------------------------------------------------------------------" <<endl;
        for(UInt_t i=0; i<CellNumbersOwnMaybe.size();i++){
            CellNumbersOwnGood.push_back(CellNumbersOwnMaybe[i]);
        }
    }
    sort(CellNumbersOwnGood.begin(),CellNumbersOwnGood.end());
    if(useMaybe){
        cout << "-----------------------------------------------------------------------------------" <<endl;
        cout << "Cells flagged as good or maybe by user:" <<endl;
        cout << "-----------------------------------------------------------------------------------" <<endl;
    }else{
        cout << "-----------------------------------------------------------------------------------" <<endl;
        cout << "Cells flagged as good by user:" <<endl;
        cout << "-----------------------------------------------------------------------------------" <<endl;
    }
    for(UInt_t i =0;i<CellNumbersOwnGood.size();i++){
        cout <<CellNumbersOwnGood.at(i) << ", ";
    }
    cout << endl;

    //**************************************************************************************************************
    //******************* Reading in Hot Cell Log Files }***********************************************************
    //**************************************************************************************************************
    string fullInputString;

    // Normal log
    cout << "-----------------------------------------------------------------------------------" <<endl;
    cout << Form("Cleaning %s/%s.log",inputDirHotCellCompare.Data(),dataSetsHot[0].Data()) << endl;
    cout << " Removed the following cells, because they were flagged as good by user!" << endl;
    cout << "-----------------------------------------------------------------------------------" <<endl;

    ifstream fileHotCells(Form("%s/%s.log",inputDirHotCellCompare.Data(),dataSetsHot[0].Data()));
    ofstream fileHotCellsCleaned(Form("%s/%s-Cleaned.log",inputDirHotCellCompare.Data(),dataSetsHot[0].Data()));

    if(fileHotCells.is_open()){
        while (std::getline(fileHotCells, fullInputString)){
            TString fullInputTSting(fullInputString);
            Int_t cellNmb = ((TObjString*)fullInputTSting.Tokenize(" ")->At(0))->GetString().Atoi();
            if(find(CellNumbersOwnGood.begin(),CellNumbersOwnGood.end(),cellNmb) == CellNumbersOwnGood.end()){
                // Cell not present in cells flagged by user -> write away
                fileHotCellsCleaned << fullInputString << endl;
            }else{
                cout << cellNmb << ", ";
            }
        }
        cout << endl;
    } else{
        cout << Form("ERROR: Could not open %s/%s.log !",inputDirHotCellCompare.Data(),dataSetsHot[0].Data()) << endl;
        return;
    }
    fileHotCells.close();
    fileHotCellsCleaned.close();

    // Detailed log
    cout << "-----------------------------------------------------------------------------------" <<endl;
    cout << Form("Cleaning %s/%s-Detailed.log",inputDirHotCellCompare.Data(),dataSetsHot[0].Data()) << endl;
    cout << " Removed the following cells, because they were flagged as good by user!" << endl;
    cout << "-----------------------------------------------------------------------------------" <<endl;

    ifstream fileHotCellsDetailed(Form("%s/%s-Detailed.log",inputDirHotCellCompare.Data(),dataSetsHot[0].Data()));
    ofstream fileHotCellsDetailedCleaned(Form("%s/%s-Detailed-Cleaned.log",inputDirHotCellCompare.Data(),dataSetsHot[0].Data()));

    // We only have a vector for the detailed log, because we will use it
    // to check if there are any cells flagged as bad by user that are not in the detailed log file.
    vector<Int_t> DeadAndHotCellsFromDetailed;

    if(fileHotCellsDetailed.is_open()){
        while (std::getline(fileHotCellsDetailed, fullInputString)){
            TString fullInputTSting(fullInputString);
            TString tmpstring = ((TObjString*)fullInputTSting.Tokenize(" ")->At(1))->GetString();
            tmpstring.Remove(tmpstring.Length()-1); // remove last digit (the ',')
            Int_t cellNmb = tmpstring.Atoi();
            if(find(CellNumbersOwnGood.begin(),CellNumbersOwnGood.end(),cellNmb) == CellNumbersOwnGood.end()){
                // Cell not present in cells flagged by user -> write away
                fileHotCellsDetailedCleaned << fullInputString << endl;
                DeadAndHotCellsFromDetailed.push_back(cellNmb);
            }else{
                cout <<  cellNmb << ", ";
            }
        }
        cout<<endl;
    } else{
        cout << Form("ERROR: Could not open %s/%s-Detailed.log !",inputDirHotCellCompare.Data(),dataSetsHot[0].Data()) << endl;
        return;
    }

    fileHotCellsDetailed.close();
    fileHotCellsDetailedCleaned.close();

    // Sorted by run log
    cout << "-----------------------------------------------------------------------------------" <<endl;
    cout << Form("Cleaning %s/%s-SortedByRun.log",inputDirHotCellCompare.Data(),dataSetsHot[0].Data()) << endl;
    cout << " Removed the following cells, because they were flagged as good by user!" << endl;
    cout << "-----------------------------------------------------------------------------------" <<endl;

    ifstream fileHotCellsSortedByRun(Form("%s/%s-SortedByRun.log",inputDirHotCellCompare.Data(),dataSetsHot[0].Data()));
    ofstream fileHotCellsSortedByRunCleaned(Form("%s/%s-SortedByRun-Cleaned.log",inputDirHotCellCompare.Data(),dataSetsHot[0].Data()));

    if(fileHotCellsSortedByRun.is_open()){
        while (std::getline(fileHotCellsSortedByRun, fullInputString)){
            TString fullInputTSting(fullInputString);
            Int_t cellNmb = ((TObjString*)fullInputTSting.Tokenize(" ")->At(0))->GetString().Atoi();
            if(find(CellNumbersOwnGood.begin(),CellNumbersOwnGood.end(),cellNmb) == CellNumbersOwnGood.end()){
                // Cell not present in cells flagged by user -> write away
                fileHotCellsSortedByRunCleaned << fullInputString << endl;
            }else{
                cout << cellNmb << ", ";
            }
        }
        cout << endl;
    } else{
        cout << Form("ERROR: Could not open %s/%s-SortedByRun.log !",inputDirHotCellCompare.Data(),dataSetsHot[0].Data()) << endl;
        return;
    }
    fileHotCellsSortedByRun.close();
    fileHotCellsSortedByRunCleaned.close();

    // Clean log file Runwise
    cout << "-----------------------------------------------------------------------------------" <<endl;
    cout << Form("Cleaning %s/%s-Runwise.log",inputDirHotCellCompare.Data(),dataSetsHot[0].Data()) << endl;
    cout << " Removed the following cells, because they were flagged as good by user!" << endl;
    cout << "-----------------------------------------------------------------------------------" <<endl;

    ifstream fileHotCellsRunwise(Form("%s/%s-Runwise.log",inputDirHotCellCompare.Data(),dataSetsHot[0].Data()));
    ofstream fileHotCellsRunwiseCleaned(Form("%s/%s-Runwise-Cleaned.log",inputDirHotCellCompare.Data(),dataSetsHot[0].Data()));

    if(fileHotCellsRunwise.is_open()){
        while (std::getline(fileHotCellsRunwise, fullInputString)){
            TString fullInputTSting(fullInputString);
            Int_t cellNmb = ((TObjString*)fullInputTSting.Tokenize("-")->At(0))->GetString().Atoi();
            if(find(CellNumbersOwnGood.begin(),CellNumbersOwnGood.end(),cellNmb) == CellNumbersOwnGood.end()){
                // Cell not present in cells flagged by user -> write away
                fileHotCellsRunwiseCleaned << fullInputString << endl;
            }else{
                cout << cellNmb << ", " ;
            }
        }
        cout << endl;
    } else{
        cout << Form("ERROR: Could not open %s/%s-Runwise.log !",inputDirHotCellCompare.Data(),dataSetsHot[0].Data()) << endl;
        return;
    }
    fileHotCellsRunwise.close();
    fileHotCellsRunwiseCleaned.close();

    //**************************************************************************************************************
    //******************* Reading in Dead Cell Log Files }***********************************************************
    //**************************************************************************************************************

    // Normal log
    cout << "-----------------------------------------------------------------------------------" <<endl;
    cout << Form("Cleaning %s/%s.log",inputDirDeadCellCompare.Data(),dataSetsCold[0].Data()) << endl;
    cout << " Removed the following cells, because they were flagged as good by user!" << endl;
    cout << "-----------------------------------------------------------------------------------" <<endl;

    ifstream fileDeadCells(Form("%s/%s.log",inputDirDeadCellCompare.Data(),dataSetsCold[0].Data()));
    ofstream fileDeadCellsCleaned(Form("%s/%s-Cleaned.log",inputDirDeadCellCompare.Data(),dataSetsCold[0].Data()));

    if(fileDeadCells.is_open()){
        while (std::getline(fileDeadCells, fullInputString)){
            TString fullInputTSting(fullInputString);
            Int_t cellNmb = ((TObjString*)fullInputTSting.Tokenize(" ")->At(0))->GetString().Atoi();
            if(find(CellNumbersOwnGood.begin(),CellNumbersOwnGood.end(),cellNmb) == CellNumbersOwnGood.end()){
                // Cell not present in cells flagged by user -> write away
                fileDeadCellsCleaned << fullInputString << endl;
            }else{
                cout <<  cellNmb << ", ";
            }
        }
        cout << endl;
    } else{
        cout << Form("ERROR: Could not open %s/%s.log !",inputDirDeadCellCompare.Data(),dataSetsCold[0].Data()) << endl;
        return;
    }
    fileDeadCells.close();
    fileDeadCellsCleaned.close();

    // Detailed log
    cout << "-----------------------------------------------------------------------------------" <<endl;
    cout << Form("Cleaning %s/%s-Detailed.log",inputDirDeadCellCompare.Data(),dataSetsCold[0].Data()) << endl;
    cout << " Removed the following cells, because they were flagged as good by user!" << endl;
    cout << "-----------------------------------------------------------------------------------" <<endl;

    ifstream fileDeadCellsDetailed(Form("%s/%s-Detailed.log",inputDirDeadCellCompare.Data(),dataSetsCold[0].Data()));
    ofstream fileDeadCellsDetailedCleaned(Form("%s/%s-Detailed-Cleaned.log",inputDirDeadCellCompare.Data(),dataSetsCold[0].Data()));

    if(fileDeadCellsDetailed.is_open()){
        while (std::getline(fileDeadCellsDetailed, fullInputString)){
            TString fullInputTSting(fullInputString);
            TString tmpstring = ((TObjString*)fullInputTSting.Tokenize(",")->At(0))->GetString();
            tmpstring.Remove(0,7);
            Int_t cellNmb = tmpstring.Atoi();
            if(find(CellNumbersOwnGood.begin(),CellNumbersOwnGood.end(),cellNmb) == CellNumbersOwnGood.end()){
                // Cell not present in cells flagged by user -> write away
                fileDeadCellsDetailedCleaned << fullInputString << endl;
                DeadAndHotCellsFromDetailed.push_back(cellNmb);
                if(find(CellNumbersOwnGood.begin(),CellNumbersOwnGood.end(),cellNmb) == CellNumbersOwnGood.end()){
                    // The Cell wasn't flagged as good, b
                }
            }else{
                cout <<cellNmb << ", ";
            }
        }
        cout << endl;
    } else{
        cout << Form("ERROR: Could not open %s/%s-Detailed.log !",inputDirDeadCellCompare.Data(),dataSetsCold[0].Data()) << endl;
        return;
    }

    sort(DeadAndHotCellsFromDetailed.begin(),DeadAndHotCellsFromDetailed.end());
    for(UInt_t i = 0; i<CellNumbersOwnBad.size();i++){
        if(find(DeadAndHotCellsFromDetailed.begin(),DeadAndHotCellsFromDetailed.end(),CellNumbersOwnBad.at(i)) == DeadAndHotCellsFromDetailed.end()){
            // bad cell flagged by user was not found in log file
            CellsThatNeedCheck.push_back(CellNumbersOwnBad.at(i));
        }
    }
    fileDeadCellsDetailed.close();
    fileDeadCellsDetailedCleaned.close();

    // Final log
    cout << "-----------------------------------------------------------------------------------" <<endl;
    cout << Form("Cleaning %s/%s-Final.log",inputDirDeadCellCompare.Data(),dataSetsCold[0].Data()) << endl;
    cout << " Removed the following cells, because they were flagged as good by user!" << endl;
    cout << "-----------------------------------------------------------------------------------" <<endl;

    ifstream fileDeadCellsFinal(Form("%s/%s-Final.log",inputDirDeadCellCompare.Data(),dataSetsCold[0].Data()));
    ofstream fileDeadCellsFinalCleaned(Form("%s/%s-Final-Cleaned.log",inputDirDeadCellCompare.Data(),dataSetsCold[0].Data()));

    if(fileDeadCellsFinal.is_open()){
        while (std::getline(fileDeadCellsFinal, fullInputString)){
            TString fullInputTSting(fullInputString);
            Int_t cellNmb = ((TObjString*)fullInputTSting.Tokenize(" ")->At(0))->GetString().Atoi();
            if(find(CellNumbersOwnGood.begin(),CellNumbersOwnGood.end(),cellNmb) == CellNumbersOwnGood.end()){
                // Cell not present in cells flagged by user -> write away
                fileDeadCellsFinalCleaned << fullInputString << endl;
            }else{
                cout <<  cellNmb << ", ";
            }
        }
        cout << endl;
    } else{
        cout << Form("ERROR: Could not open %s/%s-Final.log !",inputDirDeadCellCompare.Data(),dataSetsCold[0].Data()) << endl;
        return;
    }
    fileDeadCellsFinal.close();
    fileDeadCellsFinalCleaned.close();

    // Runwise log
    cout << "-----------------------------------------------------------------------------------" <<endl;
    cout << Form("Cleaning %s/%s-Runwise.log",inputDirDeadCellCompare.Data(),dataSetsCold[0].Data()) << endl;
    cout << " Removed the following cells, because they were flagged as good by user!" << endl;
    cout << "-----------------------------------------------------------------------------------" <<endl;

    ifstream fileDeadCellsRunwise(Form("%s/%s-Runwise.log",inputDirDeadCellCompare.Data(),dataSetsCold[0].Data()));
    ofstream fileDeadCellsRunwiseCleaned(Form("%s/%s-Runwise-Cleaned.log",inputDirDeadCellCompare.Data(),dataSetsCold[0].Data()));

    if(fileDeadCellsRunwise.is_open()){
        while (std::getline(fileDeadCellsRunwise, fullInputString)){
            TString fullInputTSting(fullInputString);
            Int_t cellNmb = ((TObjString*)fullInputTSting.Tokenize("-")->At(0))->GetString().Atoi();
            if(find(CellNumbersOwnGood.begin(),CellNumbersOwnGood.end(),cellNmb) == CellNumbersOwnGood.end()){
                // Cell not present in cells flagged by user -> write away
                fileDeadCellsRunwiseCleaned << fullInputString << endl;
            }else{
                cout <<  cellNmb << ", ";
            }
        }
        cout << endl;
    } else{
        cout << Form("ERROR: Could not open %s/%s-Runwise.log !",inputDirDeadCellCompare.Data(),dataSetsCold[0].Data()) << endl;
        return;
    }
    fileDeadCellsRunwise.close();
    fileDeadCellsRunwiseCleaned.close();

    if(CellsThatNeedCheck.size()==0){
    cout << "-----------------------------------------------------------------------------------" <<endl;
    cout << "----------------------- Finished cleaning log files -------------------------------" << endl;
    cout << "-----------------------------------------------------------------------------------" <<endl;
    } else{
        cout << "-----------------------------------------------------------------------------------" <<endl;
        cout << "| ATTENTION ATTENTION ATTENTION ATTENTION ATTENTION ATTENTION ATTENTION ATTENTION |"<< endl;
        cout << "-----------------------------------------------------------------------------------" <<endl;
        cout << " Finished cleaning the log files, however there were " << CellsThatNeedCheck.size()<< " cells found that were" << endl;
        cout << " flagged as bad by you, but not by HotCellCompare or DeadCellCompare. Please check" << endl;
        cout << " these cells again and add them to the log file by hand if they are really bad!" << endl;
        cout << endl;
        gSystem->Exec(Form("mkdir -p %s/AfterRunwiseComparison/maybe", inputDirUserBad.Data()));
        gSystem->Exec(Form("mkdir -p %s/AfterRunwiseComparison/ok", inputDirUserBad.Data()));
        gSystem->Exec(Form("mkdir -p %s/AfterRunwiseComparison/bad", inputDirUserBad.Data()));
        TString TStrCurrentCellName;
        TString BadCellLocation;
        TString MaybeCellLocation;
        TString TargetCellLocation;
        TString TStrCpBad;
        TString TStrCpMaybe;
        TString TStrCopyBadOrMaybeString;
        for(UInt_t i = 0; i < CellsThatNeedCheck.size();i++){
            cout << CellsThatNeedCheck.at(i) << ", ";
            TStrCurrentCellName=Form("Cell%d_EnergyAndTimeComparison.%s", CellsThatNeedCheck.at(i), suffix.Data() );
            MaybeCellLocation=Form("%s/%s", inputDirUserMaybe.Data(), TStrCurrentCellName.Data() );
            BadCellLocation=Form("%s/%s", inputDirUserBad.Data(), TStrCurrentCellName.Data() );
            TargetCellLocation=Form("%s/AfterRunwiseComparison/%s", inputDirUserBad.Data(), TStrCurrentCellName.Data());
            TStrCpBad=Form("cp %s %s", BadCellLocation.Data(), TargetCellLocation.Data());
            TStrCpMaybe=Form("cp %s %s", MaybeCellLocation.Data(), TargetCellLocation.Data());
            TStrCopyBadOrMaybeString=Form("if [ -f %s ]; then %s; elif  [ -f %s ]; then %s; fi", BadCellLocation.Data(), TStrCpBad.Data(), MaybeCellLocation.Data(), TStrCpMaybe.Data() );
            gSystem->Exec(Form("%s", TStrCopyBadOrMaybeString.Data()));
        }
        cout << endl;
    }
}
