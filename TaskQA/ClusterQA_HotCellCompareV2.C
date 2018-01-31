/*******************************************************************************
 ******  provided by Gamma Conversion Group, PWGGA,                        *****
 ******     Daniel Muehlheim, d.muehlheim@cern.ch                          *****
 *******************************************************************************/

#include "QA.h"

void ClusterQA_HotCellCompareV2(    TString configFileName  = "configFile.txt",
                                    TString suffix          = "eps" ){

    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "ClusterQA_HotCellCompare" << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

    //**************************************************************************************************************
    //******************* set nice plotting settings ***************************************************************
    //**************************************************************************************************************
    gROOT->Reset();

    StyleSettingsThesis();
    SetPlotStyle();

    TString outputDir               = "ClusterQA_HotCellCompare";

    //**************************************************************************************************************
    //******************* global settings **************************************************************************
    //**************************************************************************************************************
    const Int_t maxSets             = 40;
    const Int_t maxTriggers         = 10;
    Int_t nFired                    = 1;
    Int_t threshNFired              = 1;
    Int_t threshNTotalFired         = 50;
    Int_t mode                      = -1;               // will abort if not set in file
    TString addOutputDirName        = "";

    Int_t nSets                     = 0;
    Int_t nTrigger                  = 0;

    TString fEnergyFlag             = "";
    TString addLabelRunlist         = "";
    TString dataSets        [maxSets];
    TString triggers        [maxTriggers];
    TString cuts            [maxTriggers];

    // initialize arrays
    for (Int_t i = 0; i< maxSets; i++){
        dataSets[i]                 = "";
    }
    // initialize arrays
    for (Int_t i = 0; i< maxTriggers; i++){
        triggers[i]                 = "";
        cuts[i]                     = "";
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
        if (tempValue.BeginsWith("hotCellNSets",TString::kIgnoreCase)){
            nSets           = ((TString)((TObjString*)tempArr->At(1))->GetString()).Atoi();
        } else if (tempValue.BeginsWith("hotCellNTrigger",TString::kIgnoreCase)){
            nTrigger        = ((TString)((TObjString*)tempArr->At(1))->GetString()).Atoi();
        } else if (tempValue.BeginsWith("energy",TString::kIgnoreCase)){
            fEnergyFlag     = (TString)((TObjString*)tempArr->At(1))->GetString();
        } else if (tempValue.BeginsWith("mode",TString::kIgnoreCase)){
            mode            = ((TString)((TObjString*)tempArr->At(1))->GetString()).Atoi();
        } else if (tempValue.BeginsWith("hotCellAdditionalOutputDirName",TString::kIgnoreCase)){
            addOutputDirName= (TString)((TObjString*)tempArr->At(1))->GetString();
        } else if (tempValue.BeginsWith("addLabelRunlist",TString::kIgnoreCase)){
            addLabelRunlist = (TString)((TObjString*)tempArr->At(1))->GetString();
        } else if (tempValue.BeginsWith("hotCellThreshNFired",TString::kIgnoreCase)){
            threshNFired   = ((TString)((TObjString*)tempArr->At(1))->GetString()).Atof();
        } else if (tempValue.BeginsWith("hotCellThreshNTotalFired",TString::kIgnoreCase)){
            threshNTotalFired   = ((TString)((TObjString*)tempArr->At(1))->GetString()).Atof();
        } else if (tempValue.BeginsWith("hotCellDataSetNames",TString::kIgnoreCase)){
            for(Int_t i = 1; i<tempArr->GetEntries() && i < maxSets ; i++){
                if (((TString)((TObjString*)tempArr->At(i))->GetString()).CompareTo("stop",TString::kIgnoreCase))
                    dataSets[i-1]       = ((TString)((TObjString*)tempArr->At(i))->GetString());
                else
                    i                   = tempArr->GetEntries();
            }
        } else if (tempValue.BeginsWith("hotCellTriggerNames",TString::kIgnoreCase)){
            for(Int_t i = 1; i<tempArr->GetEntries() && i < maxSets ; i++){
                if (((TString)((TObjString*)tempArr->At(i))->GetString()).CompareTo("stop",TString::kIgnoreCase))
                    triggers[i-1]       = ((TString)((TObjString*)tempArr->At(i))->GetString());
                else
                    i                   = tempArr->GetEntries();
            }
        } else if (tempValue.BeginsWith("hotCellDataCuts",TString::kIgnoreCase)){
            for(Int_t i = 1; i<tempArr->GetEntries() && i < maxSets ; i++){
                if (((TString)((TObjString*)tempArr->At(i))->GetString()).CompareTo("stop",TString::kIgnoreCase))
                    cuts[i-1]           = ((TString)((TObjString*)tempArr->At(i))->GetString());
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
    cout << "nSets:                 " << nSets                      << endl;
    cout << "nTrigger:              " << nTrigger                   << endl;
    cout << "fEnergyFlag:           " << fEnergyFlag.Data()         << endl;
    cout << "mode:                  " << mode                       << endl;
    cout << "addLabelRunlist:       " << addLabelRunlist.Data()     << endl;
    if (nSets == 0 || !fEnergyFlag.CompareTo("") || mode == -1 ){
        cout << "ABORTING: You are missing the nSets, energy or mode setting, can't continue like that" << endl;
        return;
    }
    cout << "**************************************************************************" << endl;
    if (addOutputDirName.CompareTo("") != 0){
        cout << "INFO: adding complementary output name: " << addOutputDirName.Data() << "  as subfolder" << endl;
        outputDir=outputDir+"/"+addOutputDirName;
    }

    // need to introduce cell offset for PHOS
    Int_t offSetBadChannel  = 0;
    if (mode == 3 || mode == 5)
        offSetBadChannel    = 0;

    // create output directory
    gSystem->Exec("mkdir -p "+outputDir);

    fstream fLogRunwiseBadCells;
    fstream fLogOutput;
    fstream fLogOutputDetailed;
    fstream fLogOutputRunwise;
    std::map <Int_t,Int_t> ma;
    std::map <Int_t,Float_t> maNFired;
    std::map <Int_t,Float_t> maNTimes;
    std::vector<TString> uniqueCellIDRun;
    std::vector<TString>::iterator it;
    std::vector<Int_t> vec;
    TString currentRun                  = "";

    for(Int_t j=0; j<nSets; j++) {
        cout << dataSets[j].Data() << endl;
        fLogOutput.open(Form("%s/%s.log",outputDir.Data(),dataSets[j].Data()),ios::out);
        fLogOutputDetailed.open(Form("%s/%s-Detailed.log",outputDir.Data(),dataSets[j].Data()),ios::out);
        fLogOutputRunwise.open(Form("%s/%s-Runwise.log",outputDir.Data(),dataSets[j].Data()),ios::out);
        for(Int_t i=0; i<nTrigger; i++) {
            TString logFileName                     = Form("%s/%s/ClusterQA/%s/Runwise/HotCellsRunwise-%s%s.log", cuts[i].Data(), fEnergyFlag.Data(), suffix.Data(), dataSets[j].Data(),
                                                           triggers[i].Data());
            fLogRunwiseBadCells.open(logFileName.Data(), ios::in);
            if(fLogRunwiseBadCells.good()){
                fLogRunwiseBadCells.seekg(0L, ios::beg);
                TString fVar;
                while(!fLogRunwiseBadCells.eof()) {
                    fLogRunwiseBadCells >> fVar;
                    if(fVar.BeginsWith("Run-")) {
                        TString curr                = fVar;
                        currentRun                  = curr.Remove(0,4);
                    }
                    else if(fVar.BeginsWith("NoNoisy")||fVar.BeginsWith("NotEnough")){
                        continue;
                    }
                    else if(fVar.Sizeof()>1) {
                        TObjArray *rNumber          = fVar.Tokenize("-");
                        TObjString* rString         = (TObjString*)rNumber->At(0);
                        TObjString* rStringNFired   = (TObjString*)rNumber->At(1);
                        TObjString* rStringNTimes   = (TObjString*)rNumber->At(2);
                        TString vecString           = rString->GetString();
                        TString vecStringNFired     = rStringNFired->GetString();
                        TString vecStringNTimes     = rStringNTimes->GetString();
                        vecStringNFired.Remove(0,14);
                        vecStringNTimes.Remove(0,12);
                        TString currentRunVecString = vecString;
                        currentRunVecString         +="-";
                        currentRunVecString         +=currentRun;
                        it                          = find (uniqueCellIDRun.begin(), uniqueCellIDRun.end(), currentRunVecString);
                        if( it == uniqueCellIDRun.end() ){
                            vec.push_back(vecString.Atoi());
                            uniqueCellIDRun.push_back(currentRunVecString);
                            if( ma.find(vecString.Atoi()) != ma.end() ){
                                if(currentRun.Atoi()<ma[vecString.Atoi()]) ma[vecString.Atoi()]=currentRun.Atoi();
                            }
                            else{
                                ma[vecString.Atoi()] = currentRun.Atoi();
                            }
                        }
                        if( maNFired.find(vecString.Atoi()) != maNFired.end() ){
                            maNFired[vecString.Atoi()] += vecStringNFired.Atof();
                            maNTimes[vecString.Atoi()] += vecStringNTimes.Atof();
                        }
                        else{
                            maNFired[vecString.Atoi()] = vecStringNFired.Atof();
                            maNTimes[vecString.Atoi()] = vecStringNTimes.Atof();
                        }
                        // Add the runwise output. For each occurence should be CellID-Run-Sigma
                        fLogOutputRunwise << currentRunVecString << "-" << (Float_t) (vecStringNFired.Atof()/vecStringNTimes.Atof()) << endl; 
                    }
                }
            }else{
                cout << "\tCould not open HotCellsRunwise file for: " << dataSets[j].Data() << triggers[i].Data() << endl;
                cout << "tried to open: " << logFileName.Data() << endl;
            }
            fLogRunwiseBadCells.close();
        }

      
        selection_sort(vec.begin(), vec.end());

        Int_t beforeCellID=vec.at(0);
        Int_t currentCellID=0;
        for(Int_t iVec=1; iVec<(Int_t)vec.size(); iVec++) {
            //cout << "iVec: " << iVec << endl;
            currentCellID               = vec.at(iVec);

            if(beforeCellID!=currentCellID){
                //cout << beforeCellID << " " << nFired << " " << ma[beforeCellID] << endl;
                if( (nFired>threshNFired && maNFired[beforeCellID]>threshNTotalFired) || (maNFired[beforeCellID]>threshNTotalFired && (maNFired[beforeCellID]/maNTimes[beforeCellID]>5)) ) {
                    fLogOutput << beforeCellID+offSetBadChannel << " " << ma[beforeCellID] << endl;
                    fLogOutputDetailed  << "CellID: " << beforeCellID+offSetBadChannel << ", occurred for first time in run: " << ma[beforeCellID] << ", in total noisy in - " << nFired
                                        << " - differentRuns, where it occurred in clusters " << maNFired[beforeCellID] << " times! That is "
                                        << ((Float_t)maNFired[beforeCellID])/maNTimes[beforeCellID] << " times than the average of neighboring cells!" << endl;
                }else{
                    cout << "discard " << beforeCellID+offSetBadChannel << ", occurred " << nFired << " times and total EFrac: " << maNFired[beforeCellID] << " which is "
                         << (Float_t)maNFired[beforeCellID]/maNTimes[beforeCellID] << " times than the average of neighboring cells!" << endl;
                }
                beforeCellID=currentCellID;
                nFired=1;
            }else nFired++;
        }
        //cout << beforeCellID << " " << nFired << " " << ma[beforeCellID] << endl;
        if(nFired>threshNFired && maNFired[beforeCellID]>threshNTotalFired) {
            fLogOutput << beforeCellID+offSetBadChannel << " " << ma[beforeCellID] << endl;
            fLogOutputDetailed  << "CellID: " << beforeCellID+offSetBadChannel << ", occurred for first time in run: " << ma[beforeCellID] << ", in total noisy in - " << nFired
                                << " - differentRuns, where it occurred in clusters " << maNFired[beforeCellID] << " times! That is "
                                << ((Float_t)maNFired[beforeCellID])/maNTimes[beforeCellID] << " times than the average of neighboring cells!" << endl;
        }

        gSystem->Exec(Form("touch %s/%s-SortedByRun.log",outputDir.Data(),dataSets[j].Data()));
        gSystem->Exec(Form("sort -k 2 %s/%s.log > %s/%s-SortedByRun.log",outputDir.Data(),dataSets[j].Data(),outputDir.Data(),dataSets[j].Data()));

        vec.clear();
        ma.clear();
        uniqueCellIDRun.clear();
        fLogOutput.close();
        fLogOutputDetailed.close();
        fLogOutputRunwise.close(); 
    }

    cout << "Done with ClusterQA_HotCellCompare" << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    return;
}
