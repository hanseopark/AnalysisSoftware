/*******************************************************************************
 ******  provided by Gamma Conversion Group, PWGGA,                        *****
 ******     Daniel Muehlheim, d.muehlheim@cern.ch                          *****
 *******************************************************************************/

#include "QA.h"


// This macro is meant to compare the log files from the
void ClusterQA_DeadCellCompareV2(   TString configFileName  = "",
                                    TString suffix          = "eps"
                                )
{
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "ClusterQA_DeadCellCompare" << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;


    //**************************************************************************************************************
    //******************* set nice plotting style ******************************************************************
    //**************************************************************************************************************
    gROOT->Reset();
    StyleSettingsThesis();
    SetPlotStyle();
    TString outputDir               = "ClusterQA_DeadCellCompare";

    //**************************************************************************************************************
    //******************* global settings **************************************************************************
    //**************************************************************************************************************
    const Int_t maxSets             = 40;
    const Int_t maxTriggers         = 10;
    const Int_t maxCaloNCells       = 20000;
    Double_t fractionThesh          = 0.5;
    Int_t mode                      = -1;               // will abort if not set in file
    TString addOutputDirName        = "";

    Int_t nCaloCells                = 0;
    Int_t nSets                     = 0;
    Int_t nTrigger                  = 0;
    Int_t nMCSets                   = 0;

    TString fEnergyFlag             = "";
    TString addLabelRunlist         = "";
    TString pathRunLists            = "DownloadAndDataPrep/runlists/";
    TString dataSets        [maxSets];
    TString mcSets          [maxSets];
    TString triggers        [maxTriggers];
    TString mcCuts          [maxSets];
    TString cuts            [maxTriggers];


    // initialize arrays
    for (Int_t i = 0; i< maxSets; i++){
        dataSets[i]                 = "";
        mcSets[i]                   = "";
        mcCuts[i]                   = "";
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
        if (tempValue.BeginsWith("deadCellNSets",TString::kIgnoreCase)){
            nSets           = ((TString)((TObjString*)tempArr->At(1))->GetString()).Atoi();
        } else if (tempValue.BeginsWith("deadCellNMCSets",TString::kIgnoreCase)){
            nMCSets         = ((TString)((TObjString*)tempArr->At(1))->GetString()).Atoi();
        } else if (tempValue.BeginsWith("deadCellNTrigger",TString::kIgnoreCase)){
            nTrigger        = ((TString)((TObjString*)tempArr->At(1))->GetString()).Atoi();
        } else if (tempValue.BeginsWith("energy",TString::kIgnoreCase)){
            fEnergyFlag     = (TString)((TObjString*)tempArr->At(1))->GetString();
        } else if (tempValue.BeginsWith("mode",TString::kIgnoreCase)){
            mode            = ((TString)((TObjString*)tempArr->At(1))->GetString()).Atoi();
        } else if (tempValue.BeginsWith("deadCellAdditionalOutputDirName",TString::kIgnoreCase)){
            addOutputDirName= (TString)((TObjString*)tempArr->At(1))->GetString();
        } else if (tempValue.BeginsWith("deadCellPathRunLists",TString::kIgnoreCase)){
            pathRunLists    = (TString)((TObjString*)tempArr->At(1))->GetString();
        } else if (tempValue.BeginsWith("nCaloCells",TString::kIgnoreCase)){
            nCaloCells      = ((TString)((TObjString*)tempArr->At(1))->GetString()).Atoi();
        } else if (tempValue.BeginsWith("addLabelRunlist",TString::kIgnoreCase)){
            addLabelRunlist = (TString)((TObjString*)tempArr->At(1))->GetString();
        } else if (tempValue.BeginsWith("deadCellFractionThesh",TString::kIgnoreCase)){
            fractionThesh   = ((TString)((TObjString*)tempArr->At(1))->GetString()).Atof();
        } else if (tempValue.BeginsWith("deadCellDataSetNames",TString::kIgnoreCase)){
            for(Int_t i = 1; i<tempArr->GetEntries() && i < maxSets ; i++){
                if (((TString)((TObjString*)tempArr->At(i))->GetString()).CompareTo("stop",TString::kIgnoreCase))
                    dataSets[i-1]       = ((TString)((TObjString*)tempArr->At(i))->GetString());
                else
                    i                   = tempArr->GetEntries();
            }
        } else if (tempValue.BeginsWith("deadCellMCSetNames",TString::kIgnoreCase)){
            for(Int_t i = 1; i<tempArr->GetEntries() && i < maxSets ; i++){
                if (((TString)((TObjString*)tempArr->At(i))->GetString()).CompareTo("stop",TString::kIgnoreCase))
                    mcSets[i-1]       = ((TString)((TObjString*)tempArr->At(i))->GetString());
                else
                    i                   = tempArr->GetEntries();
            }
        } else if (tempValue.BeginsWith("deadCellMCCuts",TString::kIgnoreCase)){
            for(Int_t i = 1; i<tempArr->GetEntries() && i < maxSets ; i++){
                if (((TString)((TObjString*)tempArr->At(i))->GetString()).CompareTo("stop",TString::kIgnoreCase))
                    mcCuts[i-1]       = ((TString)((TObjString*)tempArr->At(i))->GetString());
                else
                    i                   = tempArr->GetEntries();
            }
        } else if (tempValue.BeginsWith("deadCellTriggerNames",TString::kIgnoreCase)){
            for(Int_t i = 1; i<tempArr->GetEntries() && i < maxSets ; i++){
                if (((TString)((TObjString*)tempArr->At(i))->GetString()).CompareTo("stop",TString::kIgnoreCase))
                    triggers[i-1]       = ((TString)((TObjString*)tempArr->At(i))->GetString());
                else
                    i                   = tempArr->GetEntries();
            }
        } else if (tempValue.BeginsWith("deadCellDataCuts",TString::kIgnoreCase)){
            for(Int_t i = 1; i<tempArr->GetEntries() && i < maxSets ; i++){
                if (((TString)((TObjString*)tempArr->At(i))->GetString()).CompareTo("stop",TString::kIgnoreCase))
                    cuts[i-1]       = ((TString)((TObjString*)tempArr->At(i))->GetString());
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
    cout << "nCaloCells:            " << nCaloCells                 << endl;
    cout << "nSets:                 " << nSets                      << endl;
    cout << "nMCSets:               " << nMCSets                    << endl;
    cout << "nTrigger:              " << nTrigger                   << endl;
    cout << "fEnergyFlag:           " << fEnergyFlag.Data()         << endl;
    cout << "mode:                  " << mode                       << endl;
    cout << "addLabelRunlist:       " << addLabelRunlist.Data()     << endl;
    cout << "pathRunLists:          " << pathRunLists.Data()        << endl;
    if (nSets == 0 || !fEnergyFlag.CompareTo("") || mode == -1 || nCaloCells == 0 ){
        cout << "ABORTING: You are missing the nSets, energy, nCaloCells or mode setting, can't continue like that" << endl;
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

    // creating output directory
    gSystem->Exec("mkdir -p "+outputDir);

    fstream fLogRunwiseDeadCells;
    fstream fLogOutput;
    fstream fLogOutputFinal;
    fstream fLogOutputDetailed;
    fstream fLogOutputRunwise;
    std::vector<Int_t>::iterator it;
    std::vector<Int_t> vec[nCaloCells];
    std::vector<Int_t> vecMC[nCaloCells];
    std::vector<Int_t> vecMCNotEnoughNoCold;
    std::vector<Int_t> globalRuns;
    std::map<Int_t,Int_t> mapRuns;
    std::map<Int_t,Int_t> mapRunsReverse;
    std::map<TString,Int_t> mapRunRanges;
    std::vector<Int_t> vecRuns;
    std::vector<TString> vecWhichRuns;
    std::vector<Int_t> vecNCons;
    std::vector<Int_t> vecNCellRuns;
    TString currentRun                          ="";
    TString fileRuns;
    std::vector<TString> runs;

    for(Int_t j=0; j<nSets; j++)
    {
        cout << dataSets[j].Data() << endl;

        //readin of MC cold/dead cells
        if(nMCSets>0){
            runs.clear();
            fileRuns                            = Form("%srunNumbers%s%s.txt",pathRunLists.Data(), mcSets[j].Data(),addLabelRunlist.Data());
            if(!readin(fileRuns, runs, kFALSE)) {
                cout << Form("INFO, no Run Numbers could be found for %s!",mcSets[j].Data()) << "\n looked for list:" <<  fileRuns.Data() << endl;
            }
            TString logFileName         = Form("%s/%s/ClusterQA/%s/Runwise/ColdCellsRunwise-%s.log",mcCuts[j].Data(),fEnergyFlag.Data(),suffix.Data(),mcSets[j].Data());
            fLogRunwiseDeadCells.open(logFileName.Data(), ios::in);
            if(fLogRunwiseDeadCells.good()){
                fLogRunwiseDeadCells.seekg(0L, ios::beg);
                TString fVar;
                while(!fLogRunwiseDeadCells.eof())                {
                    fLogRunwiseDeadCells >> fVar;
                    if(fVar.BeginsWith("Run-")) {
                        TString curr     = fVar;
                        currentRun       = curr.Remove(0,4);
                    } else if(fVar.BeginsWith("NoCold")||fVar.BeginsWith("NotEnough")){
                        it              = find (vecMCNotEnoughNoCold.begin(), vecMCNotEnoughNoCold.end(), currentRun.Atoi());
                        if( it == vecMCNotEnoughNoCold.end() ) vecMCNotEnoughNoCold.push_back(currentRun.Atoi());
                        continue;
                    } else if(fVar.Sizeof()>1) {
                        TObjArray *rNumber  = fVar.Tokenize("-");
                        TObjString* rString = (TObjString*)rNumber->At(0);
                        TString vecString   = rString->GetString();
                        it                  = find (vecMC[vecString.Atoi()].begin(), vecMC[vecString.Atoi()].end(), currentRun.Atoi());
                        if( it == vecMC[vecString.Atoi()].end() ) vecMC[vecString.Atoi()].push_back(currentRun.Atoi());
                    }
                }
            } else {
                cout << "\tCould not open DeadCellsRunwise file for: " << mcSets[j].Data() << endl;
                cout << "tried opening: " << logFileName.Data() << endl;;
            }
            fLogRunwiseDeadCells.close();
        }

        fLogOutput.open(Form("%s/%s.log",outputDir.Data(),dataSets[j].Data()),ios::out);
        fLogOutputDetailed.open(Form("%s/%s-Detailed.log",outputDir.Data(),dataSets[j].Data()),ios::out);
        fLogOutputRunwise.open(Form("%s/%s-Runwise.log",outputDir.Data(),dataSets[j].Data()),ios::out);
        for(Int_t i=0; i<nTrigger; i++){
            runs.clear();
            fileRuns                        = Form("%srunNumbers%s%s.txt", pathRunLists.Data(), dataSets[j].Data(), addLabelRunlist.Data() );
            if(!readin(fileRuns, runs, kFALSE)) {
                cout << Form("INFO, no Run Numbers could be found for %s!",dataSets[j].Data()) << "\n looked for list:" << fileRuns.Data() << endl;
            }

            TString logFileName         = Form("%s/%s/ClusterQA/%s/Runwise/ColdCellsRunwise-%s.log", cuts[i].Data(), fEnergyFlag.Data(), suffix.Data(), dataSets[j].Data());
            fLogRunwiseDeadCells.open(logFileName.Data(), ios::in);
            if(fLogRunwiseDeadCells.good()){
                fLogRunwiseDeadCells.seekg(0L, ios::beg);
                TString fVar;
                while(!fLogRunwiseDeadCells.eof()) {
                    fLogRunwiseDeadCells >> fVar;
                    if(fVar.BeginsWith("Run-")) {
                        TString curr        = fVar;
                        currentRun          = curr.Remove(0,4);
                    } else if(fVar.BeginsWith("NoCold")||fVar.BeginsWith("NotEnough")){
                        continue;
                    }
                    else if(fVar.Sizeof()>1) {
                        TObjArray *rNumber  = fVar.Tokenize("-");
                        TObjString* rString = (TObjString*)rNumber->At(0);
                        TString vecString   = rString->GetString();
                        if(nMCSets>0){
                            it              = find (vecMC[vecString.Atoi()].begin(), vecMC[vecString.Atoi()].end(), currentRun.Atoi());
                            if( !(it == vecMC[vecString.Atoi()].end()) ) continue;

                            it              = find (vecMCNotEnoughNoCold.begin(), vecMCNotEnoughNoCold.end(), currentRun.Atoi());
                            if( !(it == vecMCNotEnoughNoCold.end()) ) continue;
                        }
                        it                  = find (vec[vecString.Atoi()].begin(), vec[vecString.Atoi()].end(), currentRun.Atoi());
                        if( it == vec[vecString.Atoi()].end() ) vec[vecString.Atoi()].push_back(currentRun.Atoi());

                        it                  = find (globalRuns.begin(), globalRuns.end(), currentRun.Atoi());
                        if( it == globalRuns.end() ) globalRuns.push_back(currentRun.Atoi());
                        //cout << "currentRun: " << currentRun << "Vec String: " << vecString << endl;
                    }
                }
            }else{
                cout << "\tCould not open DeadCellsRunwise file for: " << dataSets[j].Data() << endl;
                cout << "tried opening: " << logFileName.Data() << endl;
            }
            fLogRunwiseDeadCells.close();
        }

        selection_sort(globalRuns.begin(),globalRuns.end());
        for(Int_t iRun=0; iRun<(Int_t)globalRuns.size(); iRun++) {
            mapRuns[iRun]=globalRuns.at(iRun);
            mapRunsReverse[globalRuns.at(iRun)]=iRun;
        }

        for(Int_t iCell=0;iCell<nCaloCells;iCell++){
            selection_sort(vec[iCell].begin(), vec[iCell].end());
            //if(vec[iCell].size()>0)cout << ", " << iCell << " - " << vec[iCell].size() << endl;
            if(vec[iCell].size()>0){
                fLogOutputDetailed << "CellID:" << iCell+offSetBadChannel << ",cold/notFiredInRuns:";
                for(Int_t j=0; j<(Int_t)vec[iCell].size(); j++){
                  fLogOutputDetailed << vec[iCell].at(j) << ",";
                  fLogOutputRunwise << iCell+offSetBadChannel << "-"<< vec[iCell].at(j) << "-" <<"10" <<endl;
                }
                fLogOutputDetailed << endl;
            }
            vecNCellRuns.push_back((Int_t)vec[iCell].size());
            if(vec[iCell].size()>1){
                Int_t nRun=1;
                Int_t nBegin=0;
                Int_t nEnd=1;
                Int_t nConsecutiveRuns=1;
                do{
                    if((mapRunsReverse[vec[iCell].at(nEnd)]-mapRunsReverse[vec[iCell].at(nEnd-1)]) == 1) {
                        nEnd++;
                        nConsecutiveRuns++;
                    }
                    else{
                        if(nConsecutiveRuns>1){
                            fLogOutput << iCell+offSetBadChannel << " " << vec[iCell].at(nBegin) << " " << vec[iCell].at(nEnd-1) << endl;
                            vecRuns.push_back(iCell);
                            vecWhichRuns.push_back(Form("%i %i",vec[iCell].at(nBegin),vec[iCell].at(nEnd-1)));
                            vecNCons.push_back(nConsecutiveRuns);
                            if(mapRunRanges.find(Form("%i %i",vec[iCell].at(nBegin),vec[iCell].at(nEnd-1)))==mapRunRanges.end()) mapRunRanges[Form("%i %i",vec[iCell].at(nBegin),vec[iCell].at(nEnd-1))]=1;
                            else mapRunRanges[Form("%i %i",vec[iCell].at(nBegin),vec[iCell].at(nEnd-1))]++;
                        }
                        nBegin=nEnd;
                        nEnd+=1;
                        nConsecutiveRuns=1;
                    }
                    nRun++;
                }while(nRun+1<(Int_t)vec[iCell].size());
                if(vec[iCell].size()==2) nEnd--;
                if(nConsecutiveRuns>1){
                    fLogOutput << iCell+offSetBadChannel << " " << vec[iCell].at(nBegin) << " " << vec[iCell].at(nEnd) << endl;
                    vecRuns.push_back(iCell);
                    vecWhichRuns.push_back(Form("%i %i",vec[iCell].at(nBegin),vec[iCell].at(nEnd)));
                    vecNCons.push_back(nConsecutiveRuns);
                    if(mapRunRanges.find(Form("%i %i",vec[iCell].at(nBegin),vec[iCell].at(nEnd)))==mapRunRanges.end()) mapRunRanges[Form("%i %i",vec[iCell].at(nBegin),vec[iCell].at(nEnd))]=1;
                    else mapRunRanges[Form("%i %i",vec[iCell].at(nBegin),vec[iCell].at(nEnd))]++;
                }
            }
        }

        fLogOutputFinal.open(Form("%s/%s-Final.log",outputDir.Data(),dataSets[j].Data()),ios::out);
        Int_t sortedOutFinal=0;

        std::vector<Int_t> checkDoubleCellsThreshold;
        for(Int_t vecCell=0;vecCell<(Int_t)vecRuns.size();vecCell++){
            if( (((Float_t)vecNCellRuns.at(vecRuns.at(vecCell))/(Float_t)globalRuns.size()) > fractionThesh) ){
                if(find (checkDoubleCellsThreshold.begin(), checkDoubleCellsThreshold.end(), vecRuns.at(vecCell)) == checkDoubleCellsThreshold.end()){
                    checkDoubleCellsThreshold.push_back(vecRuns.at(vecCell));
                    fLogOutputFinal << vecRuns.at(vecCell)+offSetBadChannel << " in " << ((Float_t)vecNCellRuns.at(vecRuns.at(vecCell))/(Float_t)globalRuns.size())*100 << "% of selected runs dead/cold" << endl;
                }
            }else if(mapRunRanges[vecWhichRuns.at(vecCell)]>10 || vecNCons.at(vecCell)>4){
                fLogOutputFinal << vecRuns.at(vecCell)+offSetBadChannel << " " << vecWhichRuns.at(vecCell) << endl;
            }else{
                sortedOutFinal++;
            }
        }
        cout << "\t discarded " << sortedOutFinal << " deadcell ranges while generating final list." << endl;
        fLogOutputFinal.close();
        fLogOutputRunwise.close();

        checkDoubleCellsThreshold.clear();
        for(Int_t iCell=0;iCell<nCaloCells;iCell++){vec[iCell].clear(); vecMC[iCell].clear();}
        vecMCNotEnoughNoCold.clear();
        globalRuns.clear();
        mapRuns.clear();
        mapRunsReverse.clear();
        mapRunRanges.clear();
        vecRuns.clear();
        vecWhichRuns.clear();
        vecNCons.clear();
        vecNCellRuns.clear();
        fLogOutput.close();
        fLogOutputDetailed.close();
    }

    cout << "Done with ClusterQA_DeadCellCompare" << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    return;
}

