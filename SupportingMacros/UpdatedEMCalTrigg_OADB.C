//!  EMCal bad channel map creation macro
/*!
    With this macro, the 2D histogram bad channel maps can be added to OADB/EMCAL/TPCdEdxRecalib.root.
*/
void updateFile(const char *fileNameOADB,TString arrName, TString fileNameInput,Int_t lowRun, Int_t highRun);
void sortOutput(const char *fileNameOADB="");
TObjArray *CreatePeriodContainer(TObjArray *inputcont);
void rebuildContainer(const char *fileNameOADB="");



/****************************************************************************
 *  NOTE: Main function which needs to be adjusted for new trigger mimicing *
 ***************************************************************************/
 void UpdatedEMCalTrigg_OADB(const char *fileNameOADBAli="$ALICE_DATA/PWGGA/EMCalTriggerMimic.root", Bool_t addOldMaps = kFALSE)
{
    gSystem->Load("libOADB");

    // create working copy of OADB file from EOS
    const char *fileNameOADB                ="EMCalTriggerMimic_temp.root";
    gSystem->Exec(Form("cp %s %s",fileNameOADBAli,fileNameOADB));

    // update OADB file with dead, bad and warm cells
    // last parameter:
    //      0: write new map for given run range
    //      1: update existing map for given run range
    //      2: same as 1 but also save over/underhead parts (special setting, not for normal use!)

    if(addOldMaps){
    }

    // updateFile(fileNameOADB,"EMCalTriggerThresholds","f_out.root",200000, 300000); //
    updateFile(fileNameOADB,"EMCalTriggerThresholds","144871_145287.root",144871, 145287);
    updateFile(fileNameOADB,"EMCalTriggerThresholds","145288_146374.root",145288, 146374);
    updateFile(fileNameOADB,"EMCalTriggerThresholds","146375_146381.root",146375, 146381);
    updateFile(fileNameOADB,"EMCalTriggerThresholds","146382_146501.root",146382, 146501);
    updateFile(fileNameOADB,"EMCalTriggerThresholds","146502_148521.root",146502, 148521);
    updateFile(fileNameOADB,"EMCalTriggerThresholds","148522_150208.root",148522, 150208);
    updateFile(fileNameOADB,"EMCalTriggerThresholds","150209_153055.root",150209, 153055);
    updateFile(fileNameOADB,"EMCalTriggerThresholds","153056_153910.root",153056, 153910);
    updateFile(fileNameOADB,"EMCalTriggerThresholds","153911_153914.root",153911, 153914);
    updateFile(fileNameOADB,"EMCalTriggerThresholds","153915_158134.root",153915, 158134);
    updateFile(fileNameOADB,"EMCalTriggerThresholds","158135_158135.root",158135, 158135);
    updateFile(fileNameOADB,"EMCalTriggerThresholds","158136_158177.root",158136, 158177);
    updateFile(fileNameOADB,"EMCalTriggerThresholds","158178_158181.root",158178, 158181);
    updateFile(fileNameOADB,"EMCalTriggerThresholds","158182_160682.root",158182, 160682);
    updateFile(fileNameOADB,"EMCalTriggerThresholds","160683_160763.root",160683, 160763);
    updateFile(fileNameOADB,"EMCalTriggerThresholds","160764_161138.root",160764, 161138);
    updateFile(fileNameOADB,"EMCalTriggerThresholds","161139_161255.root",161139, 161255);
    updateFile(fileNameOADB,"EMCalTriggerThresholds","161256_161378.root",161256, 161378);
    updateFile(fileNameOADB,"EMCalTriggerThresholds","161379_161456.root",161379, 161456);
    updateFile(fileNameOADB,"EMCalTriggerThresholds","161457_161524.root",161457, 161524);
    updateFile(fileNameOADB,"EMCalTriggerThresholds","161525_161555.root",161525, 161555);
    updateFile(fileNameOADB,"EMCalTriggerThresholds","161556_161557.root",161556, 161557);
    updateFile(fileNameOADB,"EMCalTriggerThresholds","161558_161608.root",161558, 161608);
    updateFile(fileNameOADB,"EMCalTriggerThresholds","161609_161629.root",161609, 161629);
    updateFile(fileNameOADB,"EMCalTriggerThresholds","161630_161723.root",161630, 161723);
    updateFile(fileNameOADB,"EMCalTriggerThresholds","161724_173730.root",161724, 173730);
    updateFile(fileNameOADB,"EMCalTriggerThresholds","173731_177143.root",173731, 177143);
    updateFile(fileNameOADB,"EMCalTriggerThresholds","177144_177146.root",177144, 177146);
    updateFile(fileNameOADB,"EMCalTriggerThresholds","177147_177652.root",177147, 177652);
    updateFile(fileNameOADB,"EMCalTriggerThresholds","177653_177723.root",177653, 177723);
    updateFile(fileNameOADB,"EMCalTriggerThresholds","177724_178326.root",177724, 178326);
    updateFile(fileNameOADB,"EMCalTriggerThresholds","178327_195179.root",178327, 195179);
    updateFile(fileNameOADB,"EMCalTriggerThresholds","195180_197468.root",195180, 197468);
    updateFile(fileNameOADB,"EMCalTriggerThresholds","197469_197691.root",197469, 197691);
    updateFile(fileNameOADB,"EMCalTriggerThresholds","197692_235194.root",197692, 235194);
    updateFile(fileNameOADB,"EMCalTriggerThresholds","235195_244284.root",235195, 244284);
    updateFile(fileNameOADB,"EMCalTriggerThresholds","244285_244628.root",244285, 244628);
    updateFile(fileNameOADB,"EMCalTriggerThresholds","244629_245140.root",244629, 245140);
    updateFile(fileNameOADB,"EMCalTriggerThresholds","245141_246994.root",245141, 246994);
    updateFile(fileNameOADB,"EMCalTriggerThresholds","246995_255538.root",246995, 255538);
    updateFile(fileNameOADB,"EMCalTriggerThresholds","f_out.root",255539, 258882);
    updateFile(fileNameOADB,"EMCalTriggerThresholds","258883_260215.root",258883, 260215);
    updateFile(fileNameOADB,"EMCalTriggerThresholds","f_out.root",260216, 265014);
    updateFile(fileNameOADB,"EMCalTriggerThresholds","265015_265308.root",265015, 265308);
    updateFile(fileNameOADB,"EMCalTriggerThresholds","265309_265588.root",265309, 265588);
    updateFile(fileNameOADB,"EMCalTriggerThresholds","265589_265784.root",265589, 265784);
    updateFile(fileNameOADB,"EMCalTriggerThresholds","265785_266404.root",265785, 266404);
    updateFile(fileNameOADB,"EMCalTriggerThresholds","f_out.root",266405, 267160);
    updateFile(fileNameOADB,"EMCalTriggerThresholds","267161_270530.root",267161, 270530);
    updateFile(fileNameOADB,"EMCalTriggerThresholds","f_out.root",270531, 282007);
    updateFile(fileNameOADB,"EMCalTriggerThresholds","282008_282503.root",282008, 282503);
    updateFile(fileNameOADB,"EMCalTriggerThresholds","f_out.root",282504, 295232);

    // the final output will be sorted by runnumber
    sortOutput(fileNameOADB);
    // the OADB container is additionally checked for ownership
    rebuildContainer("EMCalTriggerMimicNEW.root");
}


/*******************************************************************
 *  NOTE: Update function that removes existing BC maps from old   *
 *          file and adds new BC maps at the end of the file       *
 *  NOTE: The final file is saved and directly renamed as new      *
 *      input file for the next loop (this reduces total file size)*
 *******************************************************************/
void updateFile(const char *fileNameOADB,TString arrName, TString fileNameInput,Int_t lowRun, Int_t highRun){
    printf("\n\nAdding %s maps to OADB file:\n",arrName.Data());

    const int nTrigger=3;

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
        con                       =(AliOADBContainer*)f->Get("AliEMCalTriggerMimic");
        if(con) con->SetName("Old");
    } else {
        cout << "previous OADB file not found... this is fine if it is the first time running...." << endl;
    }

    // make the output container
    AliOADBContainer *con2                      = new AliOADBContainer("AliEMCalTriggerMimic");

    // List of brand new array
    TObjArray arrayAdd(3);
    arrayAdd.SetName(arrName.Data());

    TH1S *hEMCalTriggerMimic[nTrigger];

    // initialize histograms in array - either load existing map or create new, empty map
    TString TriggNames[3] = {"EMCalL0","EMCalL1G2","EMCalL1G1"};
    for(int i=0;i<nTrigger;i++){
        hEMCalTriggerMimic[i] = NULL;

        hEMCalTriggerMimic[i] = (TH1S*)fInput->Get(Form("%s",TriggNames[i].Data()));
        hEMCalTriggerMimic[i]->SetName(Form("%s",TriggNames[i].Data()));
        if(!hEMCalTriggerMimic[i]){
            cout << "Input histogram not found for bin " << i << "! returning...." <<  endl;
            return;
        }
    }

    // add histograms to a new array
    for (Int_t i=0;i<nTrigger;i++){
        arrayAdd.Add(hEMCalTriggerMimic[i]);
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
void sortOutput(const char *fileNameOADB=""){

    TFile *f                                    = TFile::Open(fileNameOADB);
    f->ls();
    AliOADBContainer *con                       =(AliOADBContainer*)f->Get("AliEMCalTriggerMimic");
    con->SetName("Old");

    Int_t indexAdd                              = 0;
    Int_t largerthan                            = 0;
    Int_t currentvalue                          = 0;

    AliOADBContainer *con2                      = new AliOADBContainer("AliEMCalTriggerMimic");
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

    con2->WriteToFile("EMCalTriggerMimicNEW.root");
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
void rebuildContainer(const char *fileNameOADB=""){
  TFile *reader = TFile::Open(fileNameOADB);
  AliOADBContainer *cont = static_cast<AliOADBContainer *>(reader->Get("AliEMCalTriggerMimic"));
  delete reader;

  AliOADBContainer *newcont = new AliOADBContainer("AliEMCalTriggerMimic");
  newcont->SetOwner(1);
  for(int irun = 0; irun < cont->GetNumberOfEntries(); irun++){
    newcont->AppendObject(CreatePeriodContainer(static_cast<TObjArray *>(cont->GetObjArray()->At(irun))), cont->LowerLimit(irun), cont->UpperLimit(irun));
  }

  newcont->WriteToFile("EMCalTriggerMimicfixOADB.root");

  TFile *reReader = TFile::Open("EMCalTriggerMimicfixOADB.root", "READ");
    AliOADBContainer *reCont = static_cast<AliOADBContainer *>(reReader->Get("AliEMCalTriggerMimic"));
    delete reReader;
    delete reCont;

  gSystem->Exec("mv EMCalTriggerMimicfixOADB.root EMCalTriggerMimicOADB.root");
  gSystem->Exec("rm EMCalTriggerMimicNEW.root");
  gSystem->Exec("rm EMCalTriggerMimic_temp.root");

}
