//!  EMCal bad channel map creation macro
/*!
    With this macro, the 2D histogram bad channel maps can be added to OADB/EMCAL/TPCdEdxRecalib.root.
*/
void updateFile(const char *fileNameOADB,TString arrName, TString fileNameInput,Int_t lowRun, Int_t highRun);
void sortOutput(const char *fileNameOADB="");
TObjArray *CreatePeriodContainer(TObjArray *inputcont);
void rebuildContainer(const char *fileNameOADB="");



/*******************************************************************
 *  NOTE: Main function which needs to be adjusted for new dEdx maps *
 *******************************************************************/
 void UpdatedEdxRecalib_OADB(const char *fileNameOADBAli="$ALICE_DATA/OADB/PWGGA/TPCdEdxRecalibOADB.root", Bool_t addOldMaps = kFALSE)
{
    gSystem->Load("libOADB");

    // create working copy of OADB file from EOS
    const char *fileNameOADB                ="TPCdEdxRecalib_temp.root";
    gSystem->Exec(Form("cp %s %s",fileNameOADBAli,fileNameOADB));

    // update OADB file with dead, bad and warm cells
    // last parameter:
    //      0: write new map for given run range
    //      1: update existing map for given run range
    //      2: same as 1 but also save over/underhead parts (special setting, not for normal use!)

    if(addOldMaps){
    }

//     updateFile(fileNameOADB,"TPCRecalib_LHC13d","TaskQA/OADB_inputs/LHC13d/DeDxMaps_80010113_00200009327000000000000000_0005314140.root",195681, 195873); //
//     updateFile(fileNameOADB,"TPCRecalib_LHC1Xx","TaskQA/OADB_inputs/LHC13d/DeDxMaps_80010113_00200009327000000000000000_0005314140.root",195874, 999999); //
    updateFile(fileNameOADB,"TPCRecalib_LHC13b","TaskQA/OADB_inputs/LHC13b/DeDxMaps_80010113_00200009327000000000000000_000481b1202.root", 195344, 195484);
    updateFile(fileNameOADB,"TPCRecalib_LHC13c","TaskQA/OADB_inputs/LHC13c/DeDxMaps_80010113_00200009327000000000000000_000481b1202.root", 195529, 195678);
    updateFile(fileNameOADB,"TPCRecalib_LHC15n","TaskQA/OADB_inputs/LHC15n/DeDxMaps_00010113_00200009227300000000000000_000481b1202.root",244340, 244628); //
    updateFile(fileNameOADB,"TPCRecalib_LHC16f","TaskQA/OADB_inputs/LHC16f/DeDxMaps_00010113_00200089227300000000000000_0053b1f1202.root", 253659, 253834);
    updateFile(fileNameOADB,"TPCRecalib_LHC16r","TaskQA/OADB_inputs/LHC16r/DeDxMaps_80010113_00200009327000000000000000_000481b1202.root",265588, 266319); //
    updateFile(fileNameOADB,"TPCRecalib_LHC16s","TaskQA/OADB_inputs/LHC16s/DeDxMaps_80010113_00200009327000000000000000_000481b1202.root",266404, 267111); //
    updateFile(fileNameOADB,"TPCRecalib_LHC17g_1","TaskQA/OADB_inputs/LHC17g/DeDxMaps_00010113_00200089227300000000000000_0053b1f1202_1.root", 270882, 270940);
    updateFile(fileNameOADB,"TPCRecalib_LHC17g_2","TaskQA/OADB_inputs/LHC17g/DeDxMaps_00010113_00200089227300000000000000_0053b1f1202_2.root", 270941, 271777);
    updateFile(fileNameOADB,"TPCRecalib_LHC17n","TaskQA/OADB_inputs/LHC17n/DeDxMaps_14810013_00200009327000000000000000_0053b1f1202.root", 280234, 280235); //
    updateFile(fileNameOADB,"TPCRecalib_LHC17pq","TaskQA/OADB_inputs/LHC17pq/DeDxMaps_00010113_00200009327000000000000000_000481b1202.root",282008, 282441); //
    updateFile(fileNameOADB,"TPCRecalib_LHC18c","TaskQA/OADB_inputs/LHC18c/DeDxMaps_00010113_00200089227300000000000000_0053b1f1202.root", 285466, 285958);
    // the final output will be sorted by runnumber
    sortOutput(fileNameOADB);
    // the OADB container is additionally checked for ownership
    rebuildContainer("TPCdEdxRecalibNEW.root");
}


/*******************************************************************
 *  NOTE: Update function that removes existing BC maps from old   *
 *          file and adds new BC maps at the end of the file       *
 *  NOTE: The final file is saved and directly renamed as new      *
 *      input file for the next loop (this reduces total file size)*
 *******************************************************************/
void updateFile(const char *fileNameOADB,TString arrName, TString fileNameInput,Int_t lowRun, Int_t highRun){
    printf("\n\nAdding %s maps to OADB file:\n",arrName.Data());

    const int nRBins=4;

    // load new input file
    TFile *fInput                                    = TFile::Open(fileNameInput);
    if(fInput->IsZombie()){
        cout << "new input file " << fileNameInput.Data() << " not found! returning..." << endl;
        return;
    }
    // load OADB file
    AliOADBContainer *con = NULL;
    TFile *f                                    = new TFile(fileNameOADB);
    if(!fInput->IsZombie()){
        f->ls();
        con                       =(AliOADBContainer*)f->Get("AliTPCdEdxRecalib");
        if(con) con->SetName("Old");
    } else {
        cout << "previous OADB file not found... this is fine if it is the first time running...." << endl;
    }

    // make the output container
    AliOADBContainer *con2                      = new AliOADBContainer("AliTPCdEdxRecalib");

    // List of brand new array
    TObjArray arrayAdd(16);
    arrayAdd.SetName(arrName.Data());

    TH2S *hdEdxRecalibElecRecalib[nRBins];
    TH2S *hdEdxRecalibPosiRecalib[nRBins];

    // initialize histograms in array - either load existing map or create new, empty map
    for(int i=0;i<nRBins;i++){
        hdEdxRecalibElecRecalib[i] = NULL;
        hdEdxRecalibPosiRecalib[i] = NULL;

        hdEdxRecalibElecRecalib[i] = (TH2S*)fInput->Get(Form("Ele_Cl%d_recalib",i));
        hdEdxRecalibElecRecalib[i]->SetName(Form("Ele_Cl%d_recalib",i));
        hdEdxRecalibPosiRecalib[i] = (TH2S*)fInput->Get(Form("Pos_Cl%d_recalib",i));
        hdEdxRecalibPosiRecalib[i]->SetName(Form("Pos_Cl%d_recalib",i));
        if(!hdEdxRecalibElecRecalib[i] || !hdEdxRecalibPosiRecalib[i] ){
            cout << "Input histogram not found for bin " << i << "! returning...." <<  endl;
            return;
        }
    }

    // add histograms to a new array
    for (Int_t i=0;i<nRBins;i++){
        arrayAdd.Add(hdEdxRecalibElecRecalib[i]);
        arrayAdd.Add(hdEdxRecalibPosiRecalib[i]);
    }

    // check old OADB file for existing BC maps in given run range and remove that old BC map
    Int_t replacedContainerLowLimit             = -1;
    Int_t replacedContainerHighLimit            = -1;
    if(con){
        for(int i=0;i<con->GetNumberOfEntries();i++){
            if (lowRun >= con->LowerLimit(i) && lowRun <= con->UpperLimit(i)){
                printf("\n!!! Not adding index %d for runrange %d--%d as low run limit %d is contained\n", i, con->LowerLimit(i),con->UpperLimit(i),lowRun);
                replacedContainerLowLimit           = con->LowerLimit(i);
                replacedContainerHighLimit          = con->UpperLimit(i);
            } else if (highRun >= con->LowerLimit(i) && highRun <= con->UpperLimit(i)){
                printf("\n!!! Not adding index %d for runrange %d--%d as high run limit %d is contained\n", i, con->LowerLimit(i),con->UpperLimit(i),highRun);
                replacedContainerLowLimit           = con->LowerLimit(i);
                replacedContainerHighLimit          = con->UpperLimit(i);
            } else if ((con->LowerLimit(i)==con->UpperLimit(i))&& con->LowerLimit(i) >= lowRun && con->UpperLimit(i)<=highRun){
                printf("\n!!! Not adding index %d for runrange %d--%d as high run limit %d is contained\n", i, con->LowerLimit(i),con->UpperLimit(i),highRun);
            }else{
                con2->AddDefaultObject(con->GetObjectByIndex(i));
                con2->AppendObject(con->GetObjectByIndex(i),con->LowerLimit(i),con->UpperLimit(i));
            }
        }
    } else {
        cout << "old container not found... is this the first time running???" << endl;
    }
    // add new BC maps at the end of the new OADB file (will be sorted later)
    con2->AddDefaultObject(&arrayAdd);
    con2->AppendObject(&arrayAdd, lowRun,highRun);


    // temporarilty save BC map file and rename as new input file
    con2->WriteToFile("tempCalib.root");
    gSystem->Exec(Form("mv tempCalib.root %s",fileNameOADB));

    printf("\n%s Maps have been successfully added!\n\n",arrName.Data());

//     plotBadChannelMapUpdate(hdEdxRecalib,arrName);

//     printf("...and also plotted!!\n");

}

/*******************************************************************
 *  NOTE: Sorting function to sort the final OADB file             *
 *                  by ascending runnumber                         *
 *******************************************************************/
void sortOutput(const char *fileNameOADB){

    TFile *f                                    = TFile::Open(fileNameOADB);
    f->ls();
    AliOADBContainer *con                       =(AliOADBContainer*)f->Get("AliTPCdEdxRecalib");
    con->SetName("Old");

    Int_t indexAdd                              = 0;
    Int_t largerthan                            = 0;
    Int_t currentvalue                          = 0;

    AliOADBContainer *con2                      = new AliOADBContainer("AliTPCdEdxRecalib");
    con2->SetOwner(1);
    // First entry needs to be added before sorting loop
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
        con2->AppendObject(con->GetObjectByIndex(indexAdd),con->LowerLimit(indexAdd),con->UpperLimit(indexAdd));
    }

    printf("\n\n");
    Int_t nentries2                             = con2->GetNumberOfEntries();
    for(int i=0;i<nentries2;i++){
        printf("\n Entry2 --> %d/%d -->",i,nentries2);
        printf("%d -- %d --> obj = %p , %s", con2->LowerLimit(i),con2->UpperLimit(i),con2->GetObjectByIndex(i),con2->GetObjectByIndex(i)->GetName());
    }
    printf("\n\n");

    con2->WriteToFile("TPCdEdxRecalibNEW.root");
}


TObjArray *CreatePeriodContainer(TObjArray *inputcont){
  TObjArray *newcont = new TObjArray(inputcont->GetEntries());
  newcont->SetName(inputcont->GetName());
  for(int i = 0; i < inputcont->GetEntries(); i++){
    newcont->AddAt(inputcont->At(i)->Clone(), i);
  }
  return newcont;
}

/*******************************************************************
 *                                                                 *
 *        NOTE: Function required to fix OADB ownership            *
 *                                                                 *
 *******************************************************************/
void rebuildContainer(const char *fileNameOADB){
  TFile *reader = TFile::Open(fileNameOADB);
  AliOADBContainer *cont = static_cast<AliOADBContainer *>(reader->Get("AliTPCdEdxRecalib"));
  delete reader;

  AliOADBContainer *newcont = new AliOADBContainer("AliTPCdEdxRecalib");
  newcont->SetOwner(1);
  for(int irun = 0; irun < cont->GetNumberOfEntries(); irun++){
    newcont->AppendObject(CreatePeriodContainer(static_cast<TObjArray *>(cont->GetObjArray()->At(irun))), cont->LowerLimit(irun), cont->UpperLimit(irun));
  }

  newcont->WriteToFile("TPCdEdxRecalibOADBfix.root");

  TFile *reReader = TFile::Open("TPCdEdxRecalibOADBfix.root", "READ");
    AliOADBContainer *reCont = static_cast<AliOADBContainer *>(reReader->Get("AliTPCdEdxRecalib"));
    delete reReader;
    delete reCont;

  gSystem->Exec("mv TPCdEdxRecalibOADBfix.root TPCdEdxRecalibOADB.root");
  gSystem->Exec("rm TPCdEdxRecalibNEW.root");
  gSystem->Exec("rm TPCdEdxRecalib_temp.root");

}

