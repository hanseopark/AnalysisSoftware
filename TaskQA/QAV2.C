/*******************************************************************************
 ******  provided by Gamma Conversion Group, PWGGA,                        *****
 ******     Daniel Muehlheim, d.muehlheim@cern.ch                          ***** 
 *******************************************************************************/

#include "QA.h"
#include "EventQA.C"
#include "PhotonQA.C"
#include "ClusterQA.C"

void QAV2(      TString configFileName  = "config.txt",         // set selected
                Bool_t doEventQA        = kFALSE,           // switch on EventQA
                Bool_t doPhotonQA       = kFALSE,           // switch on PCM-PhotonQA
                Bool_t doClusterQA      = kFALSE,           // switch on ClusterQA
                Bool_t doMergedQA       = kFALSE,           // switch on merged ClusterQA
                Int_t doExtQA           = 2,                // 0: switched off, 1: normal extQA, 2: with Cell level plots
                TString suffix          = "eps"             // output format of plots
         ){
    
    //**************************************************************************************************************
    //******************* global settings **************************************************************************
    //**************************************************************************************************************
    const Int_t maxSets             = 40;
    Int_t cutNr                     = -1;               // if -1 & not overwritten in configFile: you have to choose number at runtime
    Int_t mode                      = -1;               // will abort if not set in file
    
    //nSets == 0 is always data!
    Int_t nSets                     = 0;
    TString fEnergyFlag             = "";
    TString labelData               = "";
    Bool_t addSubfolder             = kFALSE;
    TString select                  = "";
    
    TString DataSets        [maxSets];
    TString plotDataSets    [maxSets];
    TString pathDataSets    [maxSets];
    TString pathPhotonQA    [maxSets];
    Bool_t diffPhotonQAPath         = kFALSE;
    
    // initialize arrays
    for (Int_t i = 0; i< maxSets; i++){
        DataSets[i]                 = "";
        plotDataSets[i]             = "";
        pathDataSets[i]             = "";
        pathPhotonQA[i]             = "";
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
    for( TString tempLine; tempLine.ReadLine(fileConfigQA, kTRUE); ) {        
        // check if line should be considered
        if (tempLine.BeginsWith("%") || tempLine.BeginsWith("#")){
            continue;
        } 
        cout << tempLine.Data() << endl;
     
        // Separate the string according to tabulators
        TObjArray *tempArr  = tempLine.Tokenize("\t");
        if(tempArr->GetEntries()<1){
            cout << "nothing to be done" << endl;
            delete tempArr;
            continue;
        } else if (tempArr->GetEntries() == 1 && !((TString)((TObjString*)tempArr->At(0))->GetString()).BeginsWith("enableSubfolder")){
            // Separate the string according to space
            tempArr       = tempLine.Tokenize(" ");
            if(tempArr->GetEntries()<1){
                cout << "nothing to be done" << endl;
                delete tempArr;
                continue;
            } else if (tempArr->GetEntries() == 1 && !((TString)((TObjString*)tempArr->At(0))->GetString()).BeginsWith("enableSubfolder") ) {
                cout << ((TString)((TObjString*)tempArr->At(0))->GetString()).Data() << " has not be reset, no value given!" << endl;    
                delete tempArr;
                continue;
            }
        }

        // Put them to the correct variables
        TString tempValue   = (TString)((TObjString*)tempArr->At(0))->GetString();
        if (tempValue.BeginsWith("select",TString::kIgnoreCase)){
            select          = (TString)((TObjString*)tempArr->At(1))->GetString();
        } else if (tempValue.BeginsWith("nSets",TString::kIgnoreCase)){
            nSets           = ((TString)((TObjString*)tempArr->At(1))->GetString()).Atoi();
        } else if (tempValue.BeginsWith("energy",TString::kIgnoreCase)){
            fEnergyFlag     = (TString)((TObjString*)tempArr->At(1))->GetString();
        } else if (tempValue.BeginsWith("labelData",TString::kIgnoreCase)){
            labelData       = (TString)((TObjString*)tempArr->At(1))->GetString();
        } else if (tempValue.BeginsWith("cutNr",TString::kIgnoreCase)){
            cutNr           = ((TString)((TObjString*)tempArr->At(1))->GetString()).Atoi();
        } else if (tempValue.BeginsWith("mode",TString::kIgnoreCase)){    
            mode            = ((TString)((TObjString*)tempArr->At(1))->GetString()).Atoi();
        } else if (tempValue.BeginsWith("enableSubfolder",TString::kIgnoreCase)){       
            addSubfolder    = kTRUE;
        } else if (tempValue.BeginsWith("pathDataSets",TString::kIgnoreCase)){    
            cout << "setting paths" << endl;
            for(Int_t i = 1; i<tempArr->GetEntries() && i < maxSets ; i++){
                cout << i << "\t" <<((TString)((TObjString*)tempArr->At(i))->GetString()).Data() << endl;
                if (((TString)((TObjString*)tempArr->At(i))->GetString()).CompareTo("stop",TString::kIgnoreCase))
                    pathDataSets[i-1]   = ((TString)((TObjString*)tempArr->At(i))->GetString());
                else 
                    i                   = tempArr->GetEntries();
            }    
        } else if (tempValue.BeginsWith("pathPhotonQA",TString::kIgnoreCase)){        
            cout << "setting paths" << endl;
            diffPhotonQAPath            = kTRUE;
            for(Int_t i = 1; i<tempArr->GetEntries() && i < maxSets ; i++){
                cout << i << "\t" <<((TString)((TObjString*)tempArr->At(i))->GetString()).Data() << endl;
                if (!((TString)((TObjString*)tempArr->At(i))->GetString()).CompareTo("stop",TString::kIgnoreCase))
                    pathPhotonQA[i-1]   = ((TString)((TObjString*)tempArr->At(i))->GetString());
                else 
                    i                   = tempArr->GetEntries();
            }    
        } else if (tempValue.BeginsWith("DataSetNamesPlot",TString::kIgnoreCase)){    
            for(Int_t i = 1; i<tempArr->GetEntries() && i < maxSets ; i++){
                if (((TString)((TObjString*)tempArr->At(i))->GetString()).CompareTo("stop",TString::kIgnoreCase))
                    plotDataSets[i-1]   = ((TString)((TObjString*)tempArr->At(i))->GetString());
                else 
                    i                   = tempArr->GetEntries();
            }    
        } else if (tempValue.BeginsWith("DataSetNames",TString::kIgnoreCase)){    
            for(Int_t i = 1; i<tempArr->GetEntries() && i < maxSets ; i++){
                if (((TString)((TObjString*)tempArr->At(i))->GetString()).CompareTo("stop",TString::kIgnoreCase))
                    DataSets[i-1]       = ((TString)((TObjString*)tempArr->At(i))->GetString());
                else 
                    i                   = tempArr->GetEntries();
            }    
        }    
        
        delete tempArr;
    }
    
    //**************************************************************************************************************
    //******************************* Check wether settings were valid *********************************************
    //**************************************************************************************************************
    
    cout << "**************************************************************************" << endl;
    cout << "**************** Settings found in config file ***************************" << endl;
    cout << "**************************************************************************" << endl;
    cout << "select:\t"<< select.Data() << endl;
    cout << "nSets:\t"<< nSets << endl;
    cout << "fEnergyFlag:\t"<< fEnergyFlag.Data() << endl;
    cout << "cutNr:\t"<< cutNr << endl;
    cout << "mode:\t"<< mode << endl;
    cout << "addSubfolder:\t"<< addSubfolder << endl;
    cout << "labelData:\t" << labelData.Data() << endl;
    
    if (nSets == 0 || !fEnergyFlag.CompareTo("") || mode == -1 ){
        cout << "ABORTING: You are missing the nSets, energy or mode setting, can't continue like that" << endl;
        return;
    }
    if (!diffPhotonQAPath){
        cout << "WARNING: photonQA disabled, no path set."<< endl;
        doPhotonQA          = kFALSE;
    }
    cout << "**************************************************************************" << endl;
    cout << "Data set setup: " << endl;
    for (Int_t i = 0; i < nSets; i++){
        if (plotDataSets[i].CompareTo(""))
            plotDataSets[i]             = DataSets[i];
        cout << i << "\t" << DataSets[i].Data() << "\t" << plotDataSets[i].Data() << "\t" << pathDataSets[i].Data();
        if (doPhotonQA && diffPhotonQAPath) cout << "\t" << pathPhotonQA[i].Data();
        cout << endl;
    }
    cout << "**************************************************************************" << endl;
    cout << "**************************************************************************" << endl;
    
    //**************************************************************************************************************
    //******************************  Starting individual QA macros ***********************************************
    //**************************************************************************************************************
    if ( doEventQA )    EventQA     (nSets, fEnergyFlag, DataSets, plotDataSets, pathDataSets, mode, cutNr, doExtQA, suffix, labelData, addSubfolder);
    if ( doPhotonQA )   PhotonQA    (nSets, fEnergyFlag, DataSets, plotDataSets, pathPhotonQA, mode, cutNr, doExtQA, suffix, labelData, addSubfolder);
    if ( doClusterQA )  ClusterQA   (nSets, fEnergyFlag, DataSets, plotDataSets, pathDataSets, mode, cutNr, doExtQA, suffix, labelData, addSubfolder);
    if ( doMergedQA )   ClusterQA   (nSets, fEnergyFlag, DataSets, plotDataSets, pathDataSets, mode, cutNr, doExtQA, suffix, labelData, addSubfolder, kTRUE);
    return;

}
