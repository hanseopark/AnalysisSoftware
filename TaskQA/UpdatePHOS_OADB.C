#include <fstream>

//!  PHOS bad channel map creation macro
/*!
* provide by Nicolas Schmidt:
*   With this macro, the 2D histogram bad channel maps can be added to OADB/PHOS/phosBadMap.root.
 */
//**************************************************************************************************
//************************* Function to convert the absolute cellID ********************************
//*************************     to PHOS module, row and column      ********************************
//**************************************************************************************************
void ReturnPHOSModuleAndPositionFromID(Int_t absId, Int_t * mapArray)
{
    // Converts the absolute numbering into the following array
    //  mapArray[0] = PHOS Module number
    //  mapArray[1] = Row number inside a PHOS module
    //  mapArray[2] = Column number inside a PHOS module
    Float_t id                                  = absId ;     //Float conversion for TMath::Ceil()
    Int_t fNPhi                                 = 64;         //Number of crystals along Phi direction
    Int_t fNZ                                   = 56;         //Number of crystals along Z direction
    Int_t fNCristalsInModule                    = fNPhi*fNZ ; //Number of crystals in one module
    Int_t phosmodulenumber                      = (Int_t)TMath:: Ceil( id / fNCristalsInModule ) ;
    id                                          -= ( phosmodulenumber - 1 ) *  fNPhi * fNZ ;
    mapArray[0]                                 = phosmodulenumber ;
    mapArray[1]                                 = (Int_t)TMath::Ceil( id / fNZ )  ;
    mapArray[2]                                 = (Int_t)( id - ( mapArray[1] - 1 ) * fNZ ) ;
}

//**************************************************************************************************
//************************* Function to read the bad channels from *********************************
//*************************          the given input file          *********************************
//**************************************************************************************************
Bool_t readin(TString fileRuns, std::vector<TString> &vec, Bool_t output)
{
    Bool_t ReturnValue=kTRUE;
    if(output) cout << Form("\nReading from %s...", fileRuns.Data()) << endl;
    fstream file;
    TString fVar;
    Int_t totalN=0;
    file.open(fileRuns.Data(), ios::in);
    if(file.good())
    {
        file.seekg(0L, ios::beg);
        if(output) cout << "Read in bad cells: \"";
        while(!file.eof())
        {
            file >> fVar;
            if(fVar.Sizeof()>1)
            {
                if(output) cout << fVar.Data() << ", ";
                vec.push_back(fVar);
                totalN++;
            }
        }
        if(output) cout << "\"" << endl;
    } else {
        ReturnValue=kFALSE;
    }
    file.close();

    if(output) cout << "...done!\n\nIn total " << totalN << " bad cells were read in!\n" << endl;
    return ReturnValue;
}
void updateFile(const char *fileNameOADB,TString arrName, TString fileName,Int_t lowRun, Int_t highRun, Int_t updateExistingMap);
void sortOutput(const char *fileNameOADB="");



/*******************************************************************
 *  NOTE: Update function that removes existing BC maps from old   *
 *          file and adds new BC maps at the end of the file       *
 *  NOTE: The final file is saved and directly renamed as new      *
 *      input file for the next loop (this reduces total file size)*
 *******************************************************************/
void updateFile(const char *fileNameOADB, TString arrName, TString fileName,Int_t lowRun, Int_t highRun, Int_t updateExistingMap){
    printf("\n\nAdding %s maps to OADB file:\n",arrName.Data());
    Bool_t bPrintOutput=kFALSE;
    //Check if lowRun is smaller than highRun
    if (lowRun > highRun){
        Int_t iBufferVariable=lowRun;
        lowRun=highRun;
        highRun=iBufferVariable;
    }

    // number of super modules in PHOS
    const Int_t nSM = 5;

    // load OADB file and make it a temporary container
    TFile *f                                    = TFile::Open(fileNameOADB);
    f->ls();
    AliOADBContainer *con                       =(AliOADBContainer*)f->Get("phosBadMap");
    con->SetName("Old");

    // make the output container
    AliOADBContainer *con2                      = new AliOADBContainer("phosBadMap");

    // List of brand new arrays
    TObjArray arrayAdd(nSM);
    arrayAdd.SetName(arrName.Data());
    TObjArray arrayAddOverheadLow(nSM);
    arrayAddOverheadLow.SetName(Form("%s_OHLow",arrName.Data()));
    TObjArray arrayAddOverheadHigh(nSM);
    arrayAddOverheadHigh.SetName(Form("%s_OHHigh",arrName.Data()));

    // Read in the list of bad channels from the given file
    std::vector<TString> vecBadChannels;
    if(!readin(fileName, vecBadChannels,bPrintOutput)) {
        cout << "ERROR, no Bad Channels could be found! Returning..., tried to find: " << fileName.Data() << endl;
        return;
    }

    // load OADB file to have access to existing maps
    TFile *foadb                                = TFile::Open(fileNameOADB);
    AliOADBContainer *conAP                     = (AliOADBContainer*)foadb->Get("phosBadMap");
    TObjArray *arrayBClow                       = (TObjArray*)conAP->GetObject(lowRun);
    TObjArray *arrayBChigh                      = (TObjArray*)conAP->GetObject(highRun);

    // define bad channel map histogram array
    TH2I * hBadMap[nSM];
    TH2I * hBadMapOverhead[nSM];
    TCanvas* c = new TCanvas("c","",800,800);
    gStyle->SetOptStat(0);

    // initialize histograms in array - either load existing map or create new, empty map
    for(int i=0;i<nSM;i++){
        hBadMap[i]                              = NULL;
        hBadMapOverhead[i]                      = NULL;
        if(updateExistingMap){
            cout << "Trying to initialize " << Form("PHOS_BadMap_mod%d",i) << " from current AliPhysics OADB file!" << endl;
            if(arrayBClow){
                hBadMap[i]                      = (TH2I*)arrayBClow->FindObject(Form("PHOS_BadMap_mod%d",i));
                if(hBadMap[i])
                    hBadMapOverhead[i]          = (TH2I*)hBadMap[i]->Clone(Form("PHOS_BadMap_mod%d",i));
            } else if (arrayBChigh){
                hBadMap[i]                      = (TH2I*)arrayBChigh->FindObject(Form("PHOS_BadMap_mod%d",i));
                if(hBadMap[i])
                    hBadMapOverhead[i]          = (TH2I*)hBadMap[i]->Clone(Form("PHOS_BadMap_mod%d",i));
            }
            if(hBadMap[i]){
              c->cd();
              hBadMap[i]->SetTitle(Form("mod%d_%d_%d_PHOS_BadMap_old",i,lowRun,highRun));
              hBadMap[i]->Draw("colz");
              c->SaveAs(Form("mod%d_%d_%d_PHOS_BadMap_old.pdf",i,lowRun,highRun));
            }
        }
        if(!hBadMap[i]){
            if (i != 0){
                cout << "Creating new map for " << Form("PHOS_BadMap_mod%d!",i) << endl;
                hBadMap[i]                          = new TH2I(Form("PHOS_BadMap_mod%d",i),"",64,0.,64.,56,0.,56.) ;
            }
        }
    }


    TObjArray *tempArr  = fileName.Tokenize(".");
    const char *fileNameCleaned = Form("%s_cleaned.%s",((TString)((TObjString*)tempArr->At(0))->GetString()).Data(), ((TString)((TObjString*)tempArr->At(1))->GetString()).Data());
    cout << "creating cleaned file with only new bad cells: " << fileNameCleaned << endl;
//     fstream cleanedFile;
//     cleanedFile.open(fileNameCleaned), ios::out);

    // fill bad channel map with new channels from input list
    Int_t positionFromID[3];
    for(Int_t j=0; j<(Int_t) vecBadChannels.size(); j++){
        // Get the position (module, row, column) from the cellID and fill it in 2D map
        ReturnPHOSModuleAndPositionFromID(vecBadChannels.at(j).Atoi(),positionFromID);
        if (bPrintOutput){cout << "supmod: " << positionFromID[0] << "  pos1: " << positionFromID[1] << "  pos2: " << positionFromID[2] << endl;}
        if (positionFromID[0] != 0){
            hBadMap[positionFromID[0]]->SetBinContent(positionFromID[1],positionFromID[2],1);
        } else {
            cout<<"PHOS cellIDs and modules start at 1; The 0 case is skipped"<<endl;
        }
        //if(hBadMap[positionFromID[0]]->GetBinContent(positionFromID[1],positionFromID[2]) == 0){
          //cout  << vecBadChannels.at(j).Atoi() <<  endl;
        //}
//         else cout << vecBadChannels.at(j) << " channel was already marked as bad. " << endl;
    }
//     cleanedFile.close();

    // add histograms to a new array
    for (Int_t mod=1;mod<nSM;mod++){
        arrayAdd.AddAt(hBadMap[mod],mod);
        arrayAddOverheadLow.AddAt(hBadMapOverhead[mod],mod);
        arrayAddOverheadHigh.AddAt(hBadMapOverhead[mod],mod);
    }
    for(int i=0;i<nSM;i++){
        if(hBadMap[i]){
          c->cd();
          hBadMap[i]->SetTitle(Form("mod%d_%d_%d_PHOS_BadMap_new",i,lowRun,highRun));
          hBadMap[i]->Draw("colz");
          c->SaveAs(Form("mod%d_%d_%d_PHOS_BadMap_new.pdf",i,lowRun,highRun));
        }
    }
    Int_t replacedContainerLowLimit             = -1;
    Int_t replacedContainerHighLimit            = -1;
    // check old OADB file for existing BC maps in given run range and remove that old BC map
    for(int i=0;i<con->GetNumberOfEntries();i++){
        if (lowRun >= con->LowerLimit(i) && lowRun <= con->UpperLimit(i)){
            printf("\n!!! Not adding index %d for runrange %d--%d as low run limit %d is contained\n", con->GetObjectByIndex(i), con->LowerLimit(i),con->UpperLimit(i),lowRun);
            replacedContainerLowLimit           = con->LowerLimit(i);
            replacedContainerHighLimit          = con->UpperLimit(i);
        } else if (highRun >= con->LowerLimit(i) && highRun <= con->UpperLimit(i)){
            printf("\n!!! Not adding index %d for runrange %d--%d as high run limit %d is contained\n", con->GetObjectByIndex(i), con->LowerLimit(i),con->UpperLimit(i),highRun);
        } else if (lowRun <= con->LowerLimit(i) && highRun >= con->UpperLimit(i)){
              printf("\n!!! Found object #%d for runrange %d--%d as full run range %d--%d is contained\n", con->GetObjectByIndex(i), con->LowerLimit(i),con->UpperLimit(i),lowRun,highRun);
        }else{
            con2->AddDefaultObject(con->GetObjectByIndex(i));
            con2->AppendObject(con->GetObjectByIndex(i),con->LowerLimit(i),con->UpperLimit(i));
        }
    }
    // add new BC maps at the end of the new OADB file (will be sorted later)
    if(replacedContainerLowLimit>0 && lowRun > replacedContainerLowLimit){
        con2->AddDefaultObject(&arrayAddOverheadLow);
        con2->AppendObject(&arrayAddOverheadLow, replacedContainerLowLimit,lowRun-1);
    }
    con2->AddDefaultObject(&arrayAdd);
    con2->AppendObject(&arrayAdd, lowRun,highRun);
    if(replacedContainerHighLimit>0 && highRun < replacedContainerHighLimit){
        con2->AddDefaultObject(&arrayAddOverheadHigh);
        con2->AppendObject(&arrayAddOverheadHigh, highRun+1,replacedContainerHighLimit);
    }

    // temporarilty save BC map file and rename as new input file
    con2->WriteToFile("tempBC.root");
    gSystem->Exec(Form("mv tempBC.root %s",fileNameOADB));

    printf("\nMaps have been successfully added!\n\n",arrName.Data());

}


/*******************************************************************
 *  NOTE: Sorting function to sort the final OADB file             *
 *                  by ascending runnumber                         *
 *******************************************************************/
void sortOutput(const char *fileNameOADB=""){

    TFile *f                                    = TFile::Open(fileNameOADB);
    f->ls();
    AliOADBContainer *con                       =(AliOADBContainer*)f->Get("phosBadMap");
    con->SetName("Old");

    Int_t indexAdd                              = 0;
    Int_t largerthan                            = 0;
    Int_t currentvalue                          = 0;

    gSystem->Exec("rm PHOSBadMapsNEW.root");
    AliOADBContainer *con2                      = new AliOADBContainer("phosBadMap");
    con2->AddDefaultObject(con->GetObjectByIndex(0));
    con2->AppendObject(con->GetObjectByIndex(0),con->LowerLimit(0),con->UpperLimit(0));
    TString strTemp                             = "";

    // sorting magic happens here
    for(int i=1;i<con->GetNumberOfEntries();i++){
        largerthan                              = con2->UpperLimit(con2->GetNumberOfEntries()-1);
        currentvalue                            = -1;
        indexAdd                                = 0;
        for(int j=0;j<con->GetNumberOfEntries();j++){
            if(con->UpperLimit(j)<=largerthan)
                continue;
            if(currentvalue < 0){
                currentvalue                    = con->UpperLimit(j);
                indexAdd                        = j;
            }
            if(con->UpperLimit(j)<currentvalue){
                currentvalue                    = con->UpperLimit(j);
                indexAdd                        = j;
            }
        }
        con2->AddDefaultObject(con->GetObjectByIndex(indexAdd));
        con2->AppendObject(con->GetObjectByIndex(indexAdd),con->LowerLimit(indexAdd),con->UpperLimit(indexAdd));
    }

    printf("\n\n");
    Int_t nentries2                             = con2->GetNumberOfEntries();
    for(int i=0;i<nentries2;i++){
        printf("\n Entry2 --> %d/%d -->",i,nentries2);
        printf("%d -- %d --> obj = %p , %s", con2->LowerLimit(i),con2->UpperLimit(i),con2->GetObjectByIndex(i),con2->GetObjectByIndex(i)->GetName());
    }
    printf("\n\n");

    con2->WriteToFile("PHOSBadMapsNEW.root");
    gSystem->Exec(Form("rm %s",fileNameOADB));

}

/*******************************************************************
 *  NOTE: Main function which needs to be adjusted for new BC maps *
 *******************************************************************/
void UpdatePHOS_OADB(const char *fileNameOADB="$ALICE_PHYSICS/OADB/PHOS/PHOSBadMaps.root")
{
    gSystem->Load("libOADB");
    const char *fileNameOADBtemp                ="PHOSBadMaps_temp.root";
    gSystem->Exec(Form("cp %s %s",fileNameOADB,fileNameOADBtemp));

    // add additional bad channels to an existing map for the given range (last argument == 1)
    //     updateFile(fileNameOADBtemp,"BadChannels15o_V1","LHC15o/LHC15o_badCellList.log",252666,260000,1);
    // add additional bad channels to an existing map for the given range (last argument == 1)
    //     updateFile(fileNameOADBtemp,"BadChannels15o_V1","LHC15o/LHC15o_badCellList.log",252666,260000,1);

    //-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
    //=======================================================================================================================================================================
    //-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // LHC10x pass4 Florian update pp 7TeV
    /*updateFile(fileNameOADBtemp,"detailed bad channels LHC10b - Florian","badChannelListsPHOS/pp7TeV/BadCellsCleaned_LHC10b_pass4.log",114737,117223,0);
    updateFile(fileNameOADBtemp,"detailed bad channels LHC10c - Florian","badChannelListsPHOS/pp7TeV/BadCellsCleaned_LHC10c_pass4.log",118503,121040,0);
    updateFile(fileNameOADBtemp,"detailed bad channels LHC10d - Florian","badChannelListsPHOS/pp7TeV/BadCellsCleaned_LHC10d_pass4.log",121694,126437,0);
    updateFile(fileNameOADBtemp,"detailed bad channels LHC10e - Florian","badChannelListsPHOS/pp7TeV/BadCellsCleaned_LHC10e_pass4.log",127102,130850,0);
    updateFile(fileNameOADBtemp,"detailed bad channels LHC10f - Florian","badChannelListsPHOS/pp7TeV/BadCellsCleaned_LHC10f_pass4.log",133004,135031,0);*/

    //-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
    //=======================================================================================================================================================================
    //-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // LHC15n pass3 Nico update pp 5TeV
    /*updateFile(fileNameOADBtemp,"detailed bad channels LHC15n - Nico","badChannelListsPHOS/LHC15n_pass3_4.log",244340,244410,1);
    updateFile(fileNameOADBtemp,"detailed bad channels LHC15n - Nico","badChannelListsPHOS/LHC15n_pass3_3.log",244411,244452,1);
    updateFile(fileNameOADBtemp,"detailed bad channels LHC15n - Nico","badChannelListsPHOS/LHC15n_pass3_2.log",244453,244479,1);
    updateFile(fileNameOADBtemp,"detailed bad channels LHC15n - Nico","badChannelListsPHOS/LHC15n_pass3_1.log",244480,244628,1);*/

    //-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
    //=======================================================================================================================================================================
    //-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // LHC16qt pass1 Friederike & Andrea & Toon pPb 5TeV
    /*updateFile(fileNameOADBtemp,"additional bad channels LHC16qt - Andrea","badChannelListsPHOS/BadChannelsFinal16q_265309-265335.log",265309,265335,1);
    updateFile(fileNameOADBtemp,"additional bad channels LHC16qt - Andrea","badChannelListsPHOS/BadChannelsFinal16q_265336-265388.log",265336,265388,1);
    updateFile(fileNameOADBtemp,"additional bad channels LHC16qt - Andrea","badChannelListsPHOS/BadChannelsFinal16q_265419-265525.log",265389,265525,1);
    updateFile(fileNameOADBtemp,"additional bad channels LHC16qt - Andrea","badChannelListsPHOS/BadChannelsFinal16t_267163-267166.log",267163,267166,1);

    // LHC17n pass1 Friederike Xe-Xe 5TeV
    updateFile(fileNameOADBtemp,"additional bad channels LHC17n - Friederike","badChannelListsPHOS/LHC17n_pass1_280234_280235_updatefull2.log",280234,280235,1);
    updateFile(fileNameOADBtemp,"additional bad channels LHC17n - Friederike","badChannelListsPHOS/badChannelsPHOSAdditional_LHC17n.txt",280234,280235,1);
    updateFile(fileNameOADBtemp,"additional bad channels LHC17n - Friederike","badChannelListsPHOS/LHC17n_pass1_280234_280235_Additional.log",280234,280235,1);
*/
    //-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
    //=======================================================================================================================================================================
    //-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
    //Mike
    //updateFile(fileNameOADBtemp,"additional bad channels LHC18q - Mike","BadChannelsLHC18q_295581_296623.log",295581,296623,1);
    //updateFile(fileNameOADBtemp,"additional bad channels LHC18r - Mike","BadChannelsLHC18r_296690_297624.log",296690,297624,1);

    //-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
    //=======================================================================================================================================================================
    //-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
    /*//LHC16 13 TeV Periods: e, g, h, i, j, k, l, o, p //full Periods
    updateFile(fileNameOADBtemp, "additional bad channels LHC16e - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC16e_pass1.log", 253478, 253591, 1);

    updateFile(fileNameOADBtemp, "additional bad channels LHC16g - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC16g_pass1.log", 254128, 254332, 1);

    updateFile(fileNameOADBtemp, "additional bad channels LHC16h - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC16h_pass1.log", 254604, 255467, 1);

    updateFile(fileNameOADBtemp, "additional bad channels LHC16i - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC16i_pass1.log", 255539, 255618, 1);

    updateFile(fileNameOADBtemp, "additional bad channels LHC16j - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC16j_pass1.log", 256219, 256418, 1);

    updateFile(fileNameOADBtemp, "additional bad channels LHC16k Part 1 - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC16k_pass1_Part1.log", 256512, 256944, 1);
    updateFile(fileNameOADBtemp, "additional bad channels LHC16k Part 2 - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC16k_pass1_Part2.log", 257011, 257144, 1);
    updateFile(fileNameOADBtemp, "additional bad channels LHC16k Part 3 - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC16k_pass1_Part3.log", 257204, 257433, 1);
    updateFile(fileNameOADBtemp, "additional bad channels LHC16k Part 4 - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC16k_pass1_Part4.log", 257487, 257604, 1);
    updateFile(fileNameOADBtemp, "additional bad channels LHC16k Part 5 - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC16k_pass1_Part5.log", 257605, 258393, 1);
    updateFile(fileNameOADBtemp, "additional bad channels LHC16k Part 6 - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC16k_pass1_Part6.log", 258452, 258537, 1);

    updateFile(fileNameOADBtemp, "additional bad channels LHC16l - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC16l_pass1.log", 258962, 259888, 1);

    updateFile(fileNameOADBtemp, "additional bad channels LHC16o Part 1 - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC16o_pass1.log", 262424, 262428, 1);
    updateFile(fileNameOADBtemp, "additional bad channels LHC16o Part 2 - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC16o_pass1.log", 262705, 264035, 1);

    updateFile(fileNameOADBtemp, "additional bad channels LHC16p - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC16p_pass1.log", 264076, 264347, 1);*/

    //-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
    /*//LHC16 13 TeV Periods: e, g, h, i, j, k, l, o, p //Good Runs (taken from old OADB File) with Parts in Log Name
    //OADB File: 252613-253591
    updateFile(fileNameOADBtemp, "additional bad channels LHC16e - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC16e_pass1.log", 252613, 253591, 1);

    //OADB File: 254128-254332
    updateFile(fileNameOADBtemp, "additional bad channels LHC16g - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC16g_pass1.log", 254128, 254332, 1);

    //OADB File: 254378-255467
    updateFile(fileNameOADBtemp, "additional bad channels LHC16h - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC16h_pass1.log", 254378, 255467, 1);

    //OADB File: 255539-255618
    updateFile(fileNameOADBtemp, "additional bad channels LHC16i - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC16i_pass1.log", 255539, 255618, 1);

    //OADB File: 256203-256420
    updateFile(fileNameOADBtemp, "additional bad channels LHC16j - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC16j_pass1.log", 256203, 256420, 1);

    //OADB File: 256512-257011; 257012-257145; 257146-258537
    updateFile(fileNameOADBtemp, "additional bad channels LHC16k Part 1 - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC16k_pass1_Part1.log", 256512, 257011, 1);
    updateFile(fileNameOADBtemp, "additional bad channels LHC16k Part 2 - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC16k_pass1_Part2.log", 257012, 257145, 1);
    updateFile(fileNameOADBtemp, "additional bad channels LHC16k Part 3 - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC16k_pass1_Part3.log", 257146, 257433, 1);
    updateFile(fileNameOADBtemp, "additional bad channels LHC16k Part 4 - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC16k_pass1_Part4.log", 257487, 257604, 1);
    updateFile(fileNameOADBtemp, "additional bad channels LHC16k Part 5 - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC16k_pass1_Part5.log", 257605, 258393, 1);
    updateFile(fileNameOADBtemp, "additional bad channels LHC16k Part 6 - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC16k_pass1_Part6.log", 258452, 258537, 1);

    //OADB File: 258883-260014
    updateFile(fileNameOADBtemp, "additional bad channels LHC16l - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC16l_pass1.log", 258883, 260014, 1);

    //real runranges are: 16o=>262395-264035; LHC16p=>264076-264347
    //my split of 16o was: 262424-262428 and 262705-264035
    //OADB File: 262395-262635; 262636-264347
    updateFile(fileNameOADBtemp, "additional bad channels LHC16o Part 1 - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC16o_pass1_Part1.log", 262395, 262635, 1);
    updateFile(fileNameOADBtemp, "additional bad channels LHC16o Part 2 - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC16o_pass1_Part2.log", 262636, 264035, 1);

    //OADB File: 262636-264347, where as real LHC16p=>264076-264347
    updateFile(fileNameOADBtemp, "additional bad channels LHC16p - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC16p_pass1.log", 264076, 264347, 1);*/

    //-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
    /*//LHC16 13 TeV Periods: e, g, h, i, j, k, l, o, p //Good Runs (taken from old OADB File)
    //OADB File: 252613-253591
    updateFile(fileNameOADBtemp, "additional bad channels LHC16e - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC16e_pass1.log", 252613, 253591, 1);

    //OADB File: 254128-254332
    updateFile(fileNameOADBtemp, "additional bad channels LHC16g - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC16g_pass1.log", 254128, 254332, 1);

    //OADB File: 254378-255467
    updateFile(fileNameOADBtemp, "additional bad channels LHC16h - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC16h_pass1.log", 254378, 255467, 1);

    //OADB File: 255539-255618
    updateFile(fileNameOADBtemp, "additional bad channels LHC16i - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC16i_pass1.log", 255539, 255618, 1);

    //OADB File: 256203-256420
    updateFile(fileNameOADBtemp, "additional bad channels LHC16j - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC16j_pass1.log", 256203, 256420, 1);

    //OADB File: 256512-257011; 257012-257145; 257146-258537
    updateFile(fileNameOADBtemp, "additional bad channels LHC16k Part 1 - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC16k_pass1_256512_257011.log", 256512, 257011, 1);
    updateFile(fileNameOADBtemp, "additional bad channels LHC16k Part 2 - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC16k_pass1_257012_257145.log", 257012, 257145, 1);
    updateFile(fileNameOADBtemp, "additional bad channels LHC16k Part 3 - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC16k_pass1_257146_257433.log", 257146, 257433, 1);
    updateFile(fileNameOADBtemp, "additional bad channels LHC16k Part 4 - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC16k_pass1_257487_257604.log", 257487, 257604, 1);
    updateFile(fileNameOADBtemp, "additional bad channels LHC16k Part 5 - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC16k_pass1_257605_258393.log", 257605, 258393, 1);
    updateFile(fileNameOADBtemp, "additional bad channels LHC16k Part 6 - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC16k_pass1_258452_258537.log", 258452, 258537, 1);

    //OADB File: 258883-260014
    updateFile(fileNameOADBtemp, "additional bad channels LHC16l - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC16l_pass1.log", 258883, 260014, 1);

    //real runranges are: 16o=>262395-264035; LHC16p=>264076-264347
    //my split of 16o was: 262424-262428 and 262705-264035
    //OADB File: 262395-262635; 262636-264347
    updateFile(fileNameOADBtemp, "additional bad channels LHC16o Part 1 - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC16o_pass1_262395_262635.log", 262395, 262635, 1);
    updateFile(fileNameOADBtemp, "additional bad channels LHC16o Part 2 - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC16o_pass1_262636_264035.log", 262636, 264035, 1);

    //OADB File: 262636-264347, where as real LHC16p=>264076-264347
    updateFile(fileNameOADBtemp, "additional bad channels LHC16p - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC16p_pass1.log", 264076, 264347, 1);*/

    //-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
    //=======================================================================================================================================================================
    //-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
    //LHC17 13 TeV Periods: c, e, f, h, i, j, k, l, m, o, r //full Periods
    /*updateFile(fileNameOADBtemp, "additional bad channels LHC17c - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC17c_pass1.log", 270531, 270667, 1);

    updateFile(fileNameOADBtemp, "additional bad channels LHC17e - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC17e_pass1.log", 270822, 270830, 1);

    updateFile(fileNameOADBtemp, "additional bad channels LHC17f - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC17f_pass1.log", 270854, 270865, 1);

    updateFile(fileNameOADBtemp, "additional bad channels LHC17h - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC17h_pass1.log", 271870, 273103, 1);

    updateFile(fileNameOADBtemp, "additional bad channels LHC17i - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC17i_pass1.log", 273591, 274442, 1);

    updateFile(fileNameOADBtemp, "additional bad channels LHC17j - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC17j_pass1.log", 274593, 274671,  1);

    updateFile(fileNameOADBtemp, "additional bad channels LHC17k Part 1 - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC17k_pass1_274801_275664.log", 274801, 275664, 1);
    updateFile(fileNameOADBtemp, "additional bad channels LHC17k Part 2 - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC17k_pass1_275847_276020.log", 275847, 276020, 1);
    updateFile(fileNameOADBtemp, "additional bad channels LHC17k Part 3 - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC17k_pass1_276040_276508.log", 276040, 276508, 1);

    updateFile(fileNameOADBtemp, "additional bad channels LHC17l - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC17l_pass1.log", 276551, 278216, 1);

    updateFile(fileNameOADBtemp, "additional bad channels LHC17m - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC17m_pass1.log", 278914, 280140, 0);

    updateFile(fileNameOADBtemp, "additional bad channels LHC17o - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC17o_pass1.log", 280282, 281961, 1);

    updateFile(fileNameOADBtemp, "additional bad channels LHC17r - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC17r_pass1.log", 282528, 282704, 1);*/

    //-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
    /*//LHC17 13 TeV Periods: c, e, f, h, i, j, k, l, m, o, r //Good Runs (taken from old OADB File)
    //OADB File: 270531-270667
    updateFile(fileNameOADBtemp, "additional bad channels LHC17c - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC17c_pass1.log", 270531, 270667, 1);

    //OADB File: 270822-270852
    updateFile(fileNameOADBtemp, "additional bad channels LHC17e - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC17e_pass1.log", 270822, 270852, 1);

    //OADB File: 270854-270865
    updateFile(fileNameOADBtemp, "additional bad channels LHC17f - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC17f_pass1.log", 270854, 270865, 1);

    //OADB File: 272151-273103
    updateFile(fileNameOADBtemp, "additional bad channels LHC17h - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC17h_pass1.log", 272151, 273103, 1);

    //OADB File: 273824-274442
    updateFile(fileNameOADBtemp, "additional bad channels LHC17i - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC17i_pass1.log", 273824, 274442, 1);

    //OADB File: 274591-274671
    updateFile(fileNameOADBtemp, "additional bad channels LHC17j - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC17j_pass1.log", 274591, 274671,  1);

    //OADB File: 274736-276508
    updateFile(fileNameOADBtemp, "additional bad channels LHC17k Part 1 - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC17k_pass1_274801_275664.log", 274801, 275664, 1);
    updateFile(fileNameOADBtemp, "additional bad channels LHC17k Part 2 - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC17k_pass1_275847_276020.log", 275847, 276020, 1);
    updateFile(fileNameOADBtemp, "additional bad channels LHC17k Part 3 - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC17k_pass1_276040_276508.log", 276040, 276508, 1);

    //OADB File: 276551-278216
    updateFile(fileNameOADBtemp, "additional bad channels LHC17l - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC17l_pass1.log", 276551, 278216, 1);

    //LHC17m goes from 278818-270140, while PHOS Good goes from 279005-280140
    //OADB File: 279234-280140
    updateFile(fileNameOADBtemp, "additional bad channels LHC17m - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC17m_pass1.log", 279005, 280140, 1);

    //OADB File: 280282-281961
    updateFile(fileNameOADBtemp, "additional bad channels LHC17o - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC17o_pass1.log", 280282, 281961, 1);

    //OADB File: 282504-282704
    updateFile(fileNameOADBtemp, "additional bad channels LHC17r - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC17r_pass1.log", 282504, 282704, 1);*/

    //-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
    //=======================================================================================================================================================================
    //-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
    /*//LHC18Periods: b, d, e, f, g, h, i, j, k, l, m, n, o, p  //Good Runs (taken from old OADB File)
    //OADB File: 284891-285447
    updateFile(fileNameOADBtemp, "additional bad channels LHC18b - Jens", "badChannelListsPHOS/BadCellsByHand_LHC18b_pass1.log", 284891, 285447, 1);

    //OADB File: 285978-286350
    updateFile(fileNameOADBtemp, "additional bad channels LHC18d - Jens", "badChannelListsPHOS/BadCellsByHand_LHC18d_pass1.log", 285978, 286350, 1);

    //OADB File: 286380-286937
    updateFile(fileNameOADBtemp, "additional bad channels LHC18e - Jens", "badChannelListsPHOS/BadCellsByHand_LHC18e_pass1.log", 286380, 286937, 1);

    //OADB File: 287000-287977
    updateFile(fileNameOADBtemp, "additional bad channels LHC18f - Jens", "badChannelListsPHOS/BadCellsByHand_LHC18f_pass1.log", 287000, 287977, 1);

    //OADB File: 288619-288750
    updateFile(fileNameOADBtemp, "additional bad channels LHC18g - Jens", "badChannelListsPHOS/BadCellsByHand_LHC18g_pass1.log", 288619, 288750, 1);

    //OADB File: 288804-288806
    updateFile(fileNameOADBtemp, "additional bad channels LHC18h - Jens", "badChannelListsPHOS/BadCellsByHand_LHC18h_pass1.log", 288804, 288806, 1);

    //OADB File: 288861-288943
    updateFile(fileNameOADBtemp, "additional bad channels LHC18iandj - Jens", "badChannelListsPHOS/BadCellsByHand_LHC18iandj_pass1.log", 288861, 288943, 1);

    //OADB File: 289165-289201
    updateFile(fileNameOADBtemp, "additional bad channels LHC18k - Jens", "badChannelListsPHOS/BadCellsByHand_LHC18k_pass1.log", 289165, 289201, 1);

    //OADB File: 289240-289971
    updateFile(fileNameOADBtemp, "additional bad channels LHC18l - Jens", "badChannelListsPHOS/BadCellsByHand_LHC18l_pass1.log", 289240, 289971, 1);

    //OADB File: 290222-292839
    updateFile(fileNameOADBtemp, "additional bad channels LHC18m - Jens", "badChannelListsPHOS/BadCellsByHand_LHC18m_pass1.log", 290222, 292839, 1);

    //OADB File: 293357-293362
    updateFile(fileNameOADBtemp, "additional bad channels LHC18n - Jens", "badChannelListsPHOS/BadCellsByHand_LHC18n_pass1.log", 293357, 293362, 1);

    //OADB File: 293368-293898
    updateFile(fileNameOADBtemp, "additional bad channels LHC18o - Jens", "badChannelListsPHOS/BadCellsByHand_LHC18o_pass1.log", 293368, 293898, 1);

    //OADB File: 294009-294925
    updateFile(fileNameOADBtemp, "additional bad channels LHC18p - Jens", "badChannelListsPHOS/BadCellsByHand_LHC18p_pass1.log", 294009, 294925, 1);*/

    //-----------------------------------------------------------------------------------------------------------------------------------------------------------------------

    //LHC18Periods: b, d, e, f, g, h, i, j, k, l, m, n, o, p  //Good Runs (taken from old OADB File)
    //OADB File: 284891-285447
    updateFile(fileNameOADBtemp, "additional bad channels LHC18b - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC18b_pass1.log", 284891, 285447, 1);

    //OADB File: 285978-286350
    updateFile(fileNameOADBtemp, "additional bad channels LHC18d - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC18d_pass1.log", 285978, 286350, 1);

    //OADB File: 286380-286937
    updateFile(fileNameOADBtemp, "additional bad channels LHC18e - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC18e_pass1.log", 286380, 286937, 1);

    //OADB File: 287000-287977
    updateFile(fileNameOADBtemp, "additional bad channels LHC18f - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC18f_pass1.log", 287000, 287977, 1);

    //OADB File: 288619-288750
    updateFile(fileNameOADBtemp, "additional bad channels LHC18g - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC18g_pass1.log", 288619, 288750, 1);

    //OADB File: 288804-288806
    updateFile(fileNameOADBtemp, "additional bad channels LHC18h - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC18h_pass1.log", 288804, 288806, 1);

    //OADB File: 288861-288943
    updateFile(fileNameOADBtemp, "additional bad channels LHC18iandj - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC18iandj_pass1.log", 288861, 288943, 1);

    //OADB File: 289165-289201
    updateFile(fileNameOADBtemp, "additional bad channels LHC18k - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC18k_pass1.log", 289165, 289201, 1);

    //OADB File: 289240-289971
    updateFile(fileNameOADBtemp, "additional bad channels LHC18l - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC18l_pass1.log", 289240, 289971, 1);

    //OADB File: 290222-292839
    updateFile(fileNameOADBtemp, "additional bad channels LHC18m - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC18m_pass1.log", 290222, 292839, 1);

    //OADB File: 293357-293362
    updateFile(fileNameOADBtemp, "additional bad channels LHC18n - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC18n_pass1.log", 293357, 293362, 1);

    //OADB File: 293368-293898
    updateFile(fileNameOADBtemp, "additional bad channels LHC18o - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC18o_pass1.log", 293368, 293898, 1);

    //OADB File: 294009-294925
    updateFile(fileNameOADBtemp, "additional bad channels LHC18p - Jens", "badChannelListsPHOS/BadCellsCleaned_LHC18p_pass1.log", 294009, 294925, 1);
    //-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
    //=======================================================================================================================================================================
    //-----------------------------------------------------------------------------------------------------------------------------------------------------------------------


    // add completely new bad channels, replacing existing map in range (last argument == 0)
    //     updateFile(fileNameOADBtemp,"BadChannelsLHC10c_20171213Update","LHC10c/LHC10c_addCells.log",138125,139517,0);

    // the final output will be sorted by runnumber and saved as a different file
    sortOutput(fileNameOADBtemp);

}
