/*******************************************************************************
 ******  provided by Gamma Conversion Group, PWGGA,                        *****
 ******     Daniel Muehlheim, d.muehlheim@cern.ch                          *****
 *******************************************************************************/

#include "QA.h"

//***********************************************************************************************
//****************** ClusterQA_CellCompareV2 ****************************************************
// - this macro will provide you with plots regarding the reasons cells have been excluded
// - it will display the dead and hot cells as a function of the runnumber with their
//   deviations according to the runwise output of the dead-cell and hot-cell compare macros
//***********************************************************************************************
void ClusterQA_CellCompareV2(TString configFileName  = "configFile.txt", TString suffix = "eps"){
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "ClusterQA_CellCompare" << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

    gROOT->Reset();

    StyleSettingsThesis();
    SetPlotStyle();

    TString outputDir ="CellCompare";

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
    Int_t CellCompareNSets          = 0;
    Int_t nTrigger                  = 0;

    TString fEnergyFlag             = "";
    TString addLabelRunlist         = "";
    TString DataSets[maxSets];
    TString triggers[maxTriggers];
    TString cut;
    Int_t runStart                  = 0; // beginning of run range to be investigated if runRange is 1
    Int_t runEnd                    = 0;   // end of run range to be investigated if runRange is 1
    Int_t runRange                  = 0; // If 1 will output a log file investigating bad cells in provided range.
    TString DataPath                = "";// path for the hot cell compare output
    TString DataPath_cold           = "";// path for the cold cell compare output



    // initialize arrays
    for (Int_t i = 0; i< maxSets; i++){
        DataSets[i]                 = "";
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
        if (tempValue.BeginsWith("CellCompareNSets",TString::kIgnoreCase)){
            CellCompareNSets= ((TString)((TObjString*)tempArr->At(1))->GetString()).Atoi();
        } else if (tempValue.BeginsWith("CellCompareNTrigger",TString::kIgnoreCase)){
            nTrigger        = ((TString)((TObjString*)tempArr->At(1))->GetString()).Atoi();
        } else if (tempValue.BeginsWith("energy",TString::kIgnoreCase)){
            fEnergyFlag     = (TString)((TObjString*)tempArr->At(1))->GetString();
        } else if (tempValue.BeginsWith("cut",TString::kIgnoreCase)){
            cut             = (TString)((TObjString*)tempArr->At(1))->GetString();
        } else if (tempValue.BeginsWith("hotCellDirName",TString::kIgnoreCase)){
            DataPath        = (TString)((TObjString*)tempArr->At(1))->GetString();
        } else if (tempValue.BeginsWith("coldCellDirName",TString::kIgnoreCase)){
            DataPath_cold   = (TString)((TObjString*)tempArr->At(1))->GetString();
        } else if (tempValue.BeginsWith("runRange",TString::kIgnoreCase)){
            runRange        = ((TString)((TObjString*)tempArr->At(1))->GetString()).Atof();
        } else if (tempValue.BeginsWith("runStart",TString::kIgnoreCase)){
            runStart        = ((TString)((TObjString*)tempArr->At(1))->GetString()).Atof();
        } else if (tempValue.BeginsWith("runEnd",TString::kIgnoreCase)){
            runEnd          = ((TString)((TObjString*)tempArr->At(1))->GetString()).Atof();
        } else if (tempValue.BeginsWith("DataSetNames",TString::kIgnoreCase)){
            for(Int_t i = 1; i<tempArr->GetEntries() && i < maxSets ; i++){
                if (((TString)((TObjString*)tempArr->At(i))->GetString()).CompareTo("stop",TString::kIgnoreCase))
                    DataSets[i-1]   = ((TString)((TObjString*)tempArr->At(i))->GetString());
                else
                    i               = tempArr->GetEntries();
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
    cout << "CellCompareNSets:      " << CellCompareNSets           << endl;
    cout << "nTrigger:              " << nTrigger                   << endl;
    cout << "fEnergyFlag:           " << fEnergyFlag.Data()         << endl;
    cout << "runRange:              " << runRange                   << endl;
    cout << "runStart:              " << runStart                   << endl;
    cout << "runEnd:                " << runEnd                     << endl;
    cout << "cut:                   " << cut                        << endl;
    cout << "hotCellDirName:        " << DataPath                   << endl;
    cout << "coldCellDirName:       " << DataPath_cold              << endl;

    // ******************************* Config Error Checking *******************************************************
    if (CellCompareNSets == 0 || !fEnergyFlag.CompareTo("") ){
        cout << "ABORTING: You are missing the CellCompareNSets  energy setting, can't continue like that" << endl;
        return;
    }

    // check to see if a valid run range is provided
    if (runRange == 1 && ((runStart == 0 | runEnd == 0) | (runStart > runEnd))){
        cout << "ABORTING: You must specify a valid run range in order to continue." << endl;
        return;
    }

    // check to make sure hot cell and cold cell directories are specified
    if (DataPath == "" | DataPath_cold == ""){
        cout << "ABORTING: You must specify a valid hot cell directory AND cold cell directory in order to continue." << endl;
        return;
    }

    TString outputPath = Form("%s/%s/ClusterQA/%s/%s",cut.Data(),fEnergyFlag.Data(),suffix.Data(),outputDir.Data());
    gSystem->Exec("mkdir -p "+outputPath);
    fstream out[CellCompareNSets];

    std::vector<Int_t> cellIDsEnergy[CellCompareNSets];
    std::vector<Int_t> cellIDsTime[CellCompareNSets];
    std::vector<Int_t> cellIDsHotCells1D[CellCompareNSets];
    std::vector<Int_t> cellIDsHotCellsTime1D[CellCompareNSets];
    std::vector<Int_t> cellIDsHotCells2D[CellCompareNSets];
    std::vector<Int_t> cellIDsMissing[CellCompareNSets];
    std::vector<Int_t> allCells[CellCompareNSets];

    TString vecStr[7]={"Energy","Time","HotCells1D","HotCellsTime1D","HotCells2D","Missing","allCells"};

    std::vector<Int_t> globalCells;

    TH1D* temp[7];
    for(Int_t i=0; i<CellCompareNSets; i++){
        TString fFile = Form("%s/%s/ClusterQA/ClusterQA_%s.root", cut.Data(),fEnergyFlag.Data(), DataSets[i].Data());
        TFile* File = new TFile(fFile.Data(),"READ");
        if(File->IsZombie()) {cout << "Warning: ROOT file '" << fFile.Data() << "' could not be openend!" << endl;}

        cout << endl;
        cout << "\t\t----------------------------------------------------------------------------" << endl;
        cout << "\t\tDataSet: " << DataSets[i].Data() << ":" << endl;
        cout << "\t\tProcessing file: " << fFile.Data() << endl;
        cout << "\t\t----------------------------------------------------------------------------" << endl;
        cout << endl;

        temp[0] = (TH1D*) File->Get(vecStr[0]);
        temp[1] = (TH1D*) File->Get(vecStr[1]);
        temp[2] = (TH1D*) File->Get(vecStr[2]);
        temp[3] = (TH1D*) File->Get(vecStr[3]);
        temp[4] = (TH1D*) File->Get(vecStr[4]);
        temp[5] = (TH1D*) File->Get(vecStr[5]);
        temp[6] = (TH1D*) File->Get(vecStr[6]);
        for(Int_t j=0; j<7; j++) cout << " - "<< temp[j] << ", ";
        FillVector(temp[0],cellIDsEnergy[i]);
        FillVector(temp[1],cellIDsTime[i]);
        FillVector(temp[2],cellIDsHotCells1D[i]);
        FillVector(temp[3],cellIDsHotCellsTime1D[i]);
        FillVector(temp[4],cellIDsHotCells2D[i]);
        FillVector(temp[5],cellIDsMissing[i]);
        FillVector(temp[6],allCells[i]);
        FillGlobalCellsVector(allCells[i],globalCells);


        // Make Cell ID vs. Run Number Histo for Bad Cells (combines hot + cold cells)
        // Takes as input Runwise Log File from Hot/Cold Cell Compare (uses directory to find them)
        // In addition to histogram will output a log file with bad cells for either a certain run range (runRange =1)

        //---------------------------------------------------------------------------------------------------------------
        fstream fLogInput_Runwise;
        TString logFileName                     = Form("%s/%s-Runwise-Cleaned.log", DataPath.Data(), DataSets[i].Data());
        fLogInput_Runwise.open(logFileName.Data(), ios::in);

        fstream fLogInput_RunwiseCold;
        TString logFileNameCold                 = Form("%s/%s-Runwise-Cleaned.log", DataPath_cold.Data(), DataSets[i].Data());
        fLogInput_RunwiseCold.open(logFileName.Data(), ios::in);

        //If the user has not cleaned the log file
        if (!fLogInput_Runwise.good()){
            cout << "\n\t\t############################################################################" << endl;
            cout << "\t\tThere is no cleaned log file with the name: " <<  logFileName << endl;
            logFileName                     = Form("%s/%s-Runwise.log", DataPath.Data(), DataSets[i].Data());
            fLogInput_Runwise.open(logFileName.Data(), ios::in);
            cout << "\t\tInstead using the uncleaned version with the name: "<< logFileName << endl;
            cout << "\t\t############################################################################" << endl;
        }
        //If the user has not cleaned the log file
        if (!fLogInput_RunwiseCold.good()){
            cout << "\n\t\t############################################################################" << endl;
            cout << "\t\tThere is no cleaned log file with the name: " <<  logFileNameCold << endl;
            logFileNameCold                     = Form("%s/%s-Runwise.log", DataPath_cold.Data(), DataSets[i].Data());
            fLogInput_RunwiseCold.open(logFileNameCold.Data(), ios::in);
            cout << "\t\tInstead using the uncleaned version with the name: "<< logFileNameCold << endl;
            cout << "\t\t############################################################################" << endl;
        }


        // read in the runlist
        std::vector<TString> vecRuns;
        TString fileRuns[CellCompareNSets];
        fileRuns[i]             = Form("DownloadAndDataPrep/runlists/runNumbers%s.txt", DataSets[i].Data());
        cout << "trying to read: " << fileRuns[i].Data() << endl;
        if(!readin(fileRuns[i], vecRuns, kFALSE)) {cout << "ERROR, no Run Numbers could be found! Returning..." << endl; return;}


        Double_t xbins[vecRuns.size()];
        std::vector<Double_t> vecRunsInd;
        std::vector<Double_t> vecCellIDs;
        std::vector<Double_t> uniqueCells;
        std::vector<Double_t> vecSigma;


        // Read in from the Input File
        if (fLogInput_Runwise.good()){
            fLogInput_Runwise.seekg(0L, ios::beg);
            //loop over hot cells
            while(!fLogInput_Runwise.eof()){
                TString fCurrentLine;
                fLogInput_Runwise >> fCurrentLine;
                if(fCurrentLine.Sizeof()>1) {
                    TObjArray *rNumber          = fCurrentLine.Tokenize("-");
                    TObjString* rCellID         = (TObjString*)rNumber->At(0);
                    TObjString* rRunNumber	    = (TObjString*)rNumber->At(1);
                    TObjString* rSigma          = (TObjString*)rNumber->At(2);
                    TString fCellID             = rCellID->GetString();
                    TString fRunNumber          = rRunNumber->GetString();
                    TString fSigma              = rSigma->GetString();
                    //cout << "Cell ID: " << fCellID << " runNumber: " << fRunNumber << endl;
                    vecRunsInd.push_back(fRunNumber.Atoi());
                    vecCellIDs.push_back(fCellID.Atoi());
                    vecSigma.push_back(fSigma.Atof());
                    // check to see if it is a unique cell
                    if (std::find(uniqueCells.begin(), uniqueCells.end(), fCellID.Atof()) != uniqueCells.end()){
                    }
                    else{
                        uniqueCells.push_back(fCellID.Atof());
                        //cout << "The cell is unique " << fCellID << endl;
                    }
                }
            }
            // loop over cold cells
            if (fLogInput_RunwiseCold.good()){
                while(!fLogInput_RunwiseCold.eof()){
                    TString fCurrentLine;
                    fLogInput_RunwiseCold >> fCurrentLine;
                    if(fCurrentLine.Sizeof()>1) {
                        TObjArray *rNumber          = fCurrentLine.Tokenize("-");
                        TObjString* rCellID         = (TObjString*)rNumber->At(0);
                        TObjString* rRunNumber	    = (TObjString*)rNumber->At(1);
                        TObjString* rSigma          = (TObjString*)rNumber->At(2);
                        TString fCellID	            = rCellID->GetString();
                        TString fRunNumber          = rRunNumber->GetString();
                        TString fSigma              = rSigma->GetString();
                        vecRunsInd.push_back(fRunNumber.Atoi());
                        vecCellIDs.push_back(fCellID.Atoi());
                        vecSigma.push_back(fSigma.Atof());
                        // check to see if it is a unique cell
                        if (std::find(uniqueCells.begin(), uniqueCells.end(), fCellID.Atof()) != uniqueCells.end()){
                        }
                        else{
                            uniqueCells.push_back(fCellID.Atof());
                        }
                    }
                }
            }
            else{
                cout << "\nThere is an error opening the  cold cell log file: " << logFileNameCold << endl;
            }
            // Create the binning for the histogram.
            std::sort(uniqueCells.begin(), uniqueCells.end());
            Double_t ybins[uniqueCells.size()];
            for(UInt_t j = 0; j< vecRuns.size(); j++){
                xbins[j] = j;
            }
            for(UInt_t h = 0; h< uniqueCells.size(); h++){
                ybins[h] = h;
            }

            // Create the canvas as a function of # of Cells and # of Runs
            TCanvas* canvas2 = new TCanvas("canvas2","",200,10,vecRuns.size()*50,uniqueCells.size()*20);  // gives the page size
            Double_t leftMargin = 0.1; Double_t rightMargin = 0.15; Double_t topMargin = 0.1; Double_t bottomMargin = 0.1;
            DrawGammaCanvasSettings(canvas2, leftMargin, rightMargin, topMargin, bottomMargin);
            TH2D* RunwiseHist;
            RunwiseHist = new TH2D(Form("Bad Cells Runwise %s", DataSets[i].Data()), Form("Bad Cells Runwise %s", DataSets[i].Data()),vecRuns.size()-1, xbins, uniqueCells.size()-1, ybins );
            for(UInt_t z= 1; z< vecRuns.size();z++){
                RunwiseHist->GetXaxis()->SetBinLabel(z, vecRuns.at(z));
            }

            fstream fLogRunRange;
            TString name;
            if(runRange == 1){
                cout << endl << "Creating Bad Cell log file for Run Range for Run " << runStart << " to Run " << runEnd << endl;
                name = outputPath + "/" +"BadCells_"+ DataSets[i].Data() + "_" + std::to_string(runStart) + "_" + std::to_string(runEnd) +".log";
                cout << "Log File Path: " << name << endl;
                fLogRunRange.open(name, ios::out);
            }
            else{
                cout << endl << "Creating Bad Cell log file" << endl;
                name = outputPath + "/" +"BadCells" + DataSets[i].Data() +".log";
                cout << "Log File Path: " << name << endl;
                fLogRunRange.open(name, ios::out);
            }

            if (uniqueCells.size()>500){
                cout << "----- Histogram labels may not be clear since there are so many cells. Consider cleaning first.-----" << endl;
            }


            for(UInt_t z= 1; z< uniqueCells.size();z++){
                TString s;
                s.Form("%.0f", uniqueCells.at(z));
                RunwiseHist->GetYaxis()->SetBinLabel(z, s);
                RunwiseHist->GetYaxis()->SetLabelSize(0.01);
                if (runRange == 0){
                    fLogRunRange << uniqueCells.at(z) << endl;
                }
            }

            // Loop over runs to populate hist
            std::vector<Double_t> vecOutput;
            for(UInt_t r= 0; r< vecRunsInd.size();r++){
                if(runRange == 1){
                    if( (runStart <= vecRunsInd.at(r) <= runEnd)){
                        if (std::find(vecOutput.begin(), vecOutput.end(), vecCellIDs.at(r)) != vecOutput.end()){
                            // Do nothing if it is already there.
                        }
                        else{
                            vecOutput.push_back(vecCellIDs.at(r));
                            fLogRunRange << vecCellIDs.at(r) << endl;
                        }
                    }
                }
                Int_t binx  = -1;
                Int_t biny  = -1;;
                // Get proper bin
                for(UInt_t d = 0; d < vecRuns.size(); d++){
                    if(vecRuns.at(d).Atof() == vecRunsInd.at(r)){
                        binx = d;
                        break;
                    }
                }
                for(UInt_t d = 0; d < uniqueCells.size(); d++){
                    if( vecCellIDs.at(r) == uniqueCells.at(d)){
                        biny = d;
                        break;
                    }
                }

                if((binx!=-1)&&(biny!=-1)) {
                    RunwiseHist->SetBinContent(binx,biny, vecSigma.at(r));
                } else{
                    cout << "Error filling RunwiseHist!" << endl;
                    return;
                }
                //cout << "binx: " << binx << " biny: " << biny << "with value: " << vecSigma.at(r) << endl;
            }
            fLogRunRange.close();
            RunwiseHist->LabelsOption("v", "X");
            RunwiseHist->GetXaxis()->SetTitle("Run Number");
            RunwiseHist->GetXaxis()->SetTitleOffset(1.4);
            RunwiseHist->GetYaxis()->SetTitle("CellID");
            // Make the bad cells very obvious
            RunwiseHist->GetZaxis()->SetRangeUser(0,20);
            RunwiseHist->Draw("COLZ");
            canvas2->Update();
            canvas2->SaveAs(Form("%s/BadCellCandidates_Runwise_%s.pdf", outputPath.Data(),DataSets[i].Data()));
            canvas2->Clear();
            delete RunwiseHist;

        }
        else{
            cout << "\nThere is an error opening the hot cell log file: " << logFileName << endl;
        }
        fLogInput_Runwise.close();
    } // end loop over CellCompareNSets

    cout.flush();
    cout << endl;
    cout << "vector.size before sort and unique: " << globalCells.size() << "." << endl;
    cout << "sorting vector...";
    selection_sort(globalCells.begin(),globalCells.end());
    cout << "done" << endl;
    vector<Int_t>::iterator it;
    cout << "unique vector...";
    it = unique(globalCells.begin(),globalCells.end());
    cout << "done" << endl;
    cout << "resize vector...";
    globalCells.resize( distance(globalCells.begin(),it) );
    cout << "done" << endl;
    cout << "vector.size after sort and unique: " << globalCells.size() << "." << endl;

    // --------------------------- Begin other bad cell plotting ----------------------------------------------------------------------
    TCanvas* canvas = new TCanvas("canvas","",200,10,1350,1800);  // gives the page size
    Double_t leftMargin = 0.1; Double_t rightMargin = 0.1; Double_t topMargin = 0.06; Double_t bottomMargin = 0.1;
    DrawGammaCanvasSettings(canvas, leftMargin, rightMargin, topMargin, bottomMargin);

    TH2D* hist;

    Int_t atPlotting = 0;
    Int_t iCount = 0;
    Int_t stepSize = 150;
    do{
        atPlotting = (iCount+1)*stepSize;
        if((Int_t)globalCells.size()<=atPlotting) atPlotting = (Int_t)globalCells.size();
        hist = new TH2D("BadCellCandidates","BadCellCandidates",CellCompareNSets,0,CellCompareNSets,atPlotting-iCount*stepSize,0+iCount*stepSize,atPlotting);
        for(Int_t iX=0; iX<CellCompareNSets; iX++) hist->GetXaxis()->SetBinLabel(iX+1,DataSets[iX]);
        for(Int_t iC=iCount*stepSize; iC<atPlotting; iC++){
            hist->GetYaxis()->SetBinLabel(iC+1-iCount*stepSize,Form("%i",globalCells.at(iC)));
            for(Int_t iX=0; iX<CellCompareNSets; iX++){
                std::vector<Int_t>::iterator it;
                it = find (allCells[iX].begin(), allCells[iX].end(), globalCells.at(iC));
                if (it != allCells[iX].end()){
                    hist->SetBinContent(iX+1,iC+1-iCount*stepSize,1);
                    out[iX] << globalCells.at(iC) << endl;
                }
            }
        }

        hist->GetZaxis()->SetRangeUser(0,2);
        DrawAutoGammaHisto2D(hist,
                             "Bad Cell Candidates",
                             "Period",
                             "CellID",
                             "",
                             0,0,0,
                             0,0,0,
                             1,1,0.035,0.035,0.015,0.035);
        canvas->SetLogx(0); canvas->SetLogy(0); canvas->SetLogz(0); canvas->Update();
        canvas->SaveAs(Form("%s/BadCellCandidates_%i.%s", outputPath.Data(),iCount++,suffix.Data()));
        canvas->Clear();
        delete hist;
    }while(atPlotting != (Int_t)globalCells.size());

    for(Int_t i=0; i<CellCompareNSets; i++) out[i].close();

    cout << "Done with ClusterQA_CellCompare" << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    return;
    >>>>>>> fixed indentations
}

void FillVector(TH1* hist, std::vector<Int_t> &vec){
    if(hist){
        for(Int_t i=1; i<hist->GetXaxis()->GetNbins()+1;i++){
            TString st = hist->GetXaxis()->GetBinLabel(i);
            vec.push_back((Int_t)st.Atoi());
        }
    }
    return;
}

void FillGlobalCellsVector(std::vector<Int_t> &vec,std::vector<Int_t> &vec2){
    for(Int_t i=0; i<(Int_t)vec.size();i++){
        vec2.push_back((Int_t)vec.at(i));
    }
    return;
}

