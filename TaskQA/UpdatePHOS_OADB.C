//!  PHOS bad channel map creation macro
/*!
    With this macro, the 2D histogram bad channel maps can be added to OADB/PHOS/phosBadMap.root.
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
    }
    file.close();
    if(output) cout << "...done!\n\nIn total " << totalN << " bad cells were read in!\n" << endl;
    else return kTRUE;
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

    // LHC15n pass3 Nico update pp 5TeV
    updateFile(fileNameOADBtemp,"detailed bad channels LHC15n - Nico","badChannelListsPHOS/LHC15n_pass3_4.log",244340,244410,1);
    updateFile(fileNameOADBtemp,"detailed bad channels LHC15n - Nico","badChannelListsPHOS/LHC15n_pass3_3.log",244411,244452,1);
    updateFile(fileNameOADBtemp,"detailed bad channels LHC15n - Nico","badChannelListsPHOS/LHC15n_pass3_2.log",244453,244479,1);
    updateFile(fileNameOADBtemp,"detailed bad channels LHC15n - Nico","badChannelListsPHOS/LHC15n_pass3_1.log",244480,244628,1);

    // LHC16qt pass1 Friederike & Andrea & Toon pPb 5TeV
    updateFile(fileNameOADBtemp,"additional bad channels LHC16qt - Friederike","badChannelListsPHOS/LHC16q_addCells_265309-265335.log",265309,265335,1);
    updateFile(fileNameOADBtemp,"additional bad channels LHC16qt - Friederike","badChannelListsPHOS/LHC16q_addCells_265336-265388.log",265336,265388,1);
    updateFile(fileNameOADBtemp,"additional bad channels LHC16qt - Friederike","badChannelListsPHOS/LHC16q_addCells_265389-265525.log",265389,265525,1);
    updateFile(fileNameOADBtemp,"additional bad channels LHC16qt - Friederike","badChannelListsPHOS/LHC16t_addCells.log",267163,267166,1);

    // LHC17n pass1 Friederike Xe-Xe 5TeV
    updateFile(fileNameOADBtemp,"additional bad channels LHC17n - Friederike","badChannelListsPHOS/LHC17n_pass1_280234_280235.log",280234,280235,1);

    // add completely new bad channels, replacing existing map in range (last argument == 0)
//     updateFile(fileNameOADBtemp,"BadChannelsLHC10c_20171213Update","LHC10c/LHC10c_addCells.log",138125,139517,0);

    // the final output will be sorted by runnumber and saved as a different file
    sortOutput(fileNameOADBtemp);

}

/*******************************************************************
 *  NOTE: Update function that removes existing BC maps from old   *
 *          file and adds new BC maps at the end of the file       *
 *  NOTE: The final file is saved and directly renamed as new      *
 *      input file for the next loop (this reduces total file size)*
 *******************************************************************/
void updateFile(const char *fileNameOADB,TString arrName, TString fileName,Int_t lowRun, Int_t highRun, Int_t updateExistingMap){
    printf("\n\nAdding %s maps to OADB file:\n",arrName.Data());

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
    if(!readin(fileName, vecBadChannels,kTRUE)) {
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
        }
        if(!hBadMap[i]){
            if (i != 0){
                cout << "Creating new map for " << Form("PHOS_BadMap_mod%d!",i) << endl;
                hBadMap[i]                          = new TH2I(Form("PHOS_BadMap_mod%d",i),"",64,0.,64.,56,0.,56.) ;
            }
        }
    }

    // fill bad channel map with new channels from input list
    Int_t positionFromID[3];
    for(Int_t j=0; j<(Int_t) vecBadChannels.size(); j++){
        // Get the position (module, row, column) from the cellID and fill it in 2D map
        ReturnPHOSModuleAndPositionFromID(vecBadChannels.at(j).Atoi(),positionFromID);
        cout << "supmod: " << positionFromID[0] << "  pos1: " << positionFromID[1] << "  pos2: " << positionFromID[2] << endl;
        hBadMap[positionFromID[0]]->SetBinContent(positionFromID[1],positionFromID[2],1);
    }

    // add histograms to a new array
    for (Int_t mod=1;mod<nSM;mod++){
        arrayAdd.AddAt(hBadMap[mod],mod);
        arrayAddOverheadLow.AddAt(hBadMapOverhead[mod],mod);
        arrayAddOverheadHigh.AddAt(hBadMapOverhead[mod],mod);
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