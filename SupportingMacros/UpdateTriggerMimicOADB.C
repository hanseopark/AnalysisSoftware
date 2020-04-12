//!  EMCal bad channel map creation macro
/*!
    With this macro, the 2D histogram bad channel maps can be added to OADB/EMCAL/TPCdEdxRecalib.root.
*/
#include <TFile.h>
#include <TSystem.h>
#include <TString.h>
#include <TObject.h>
#include <TH1S.h>
#include <AliOADBContainer.h>
#include <iostream>

using namespace std;
void updateFile(const char *fileNameOADB,TString arrName, TString fileNameInput,Int_t lowRun, Int_t highRun, Double_t EMCL0, Double_t EMCdL0, Double_t EMCL1high, Double_t EMCdL1high, Double_t EMCL1low, Double_t EMCdL1low, Double_t PHOSL0, Double_t PHOSdL0);
void sortOutput(const char *fileNameOADB="");
TObjArray *CreatePeriodContainer(TObjArray *inputcont);
void rebuildContainer(const char *fileNameOADB="");



/****************************************************************************
 *  NOTE: Main function which needs to be adjusted for new trigger mimicing *
 ***************************************************************************/
 void UpdateTriggerMimicOADB(const char *fileNameOADBAli="$ALICE_DATA/OADB/PWGGA/EMCalTriggerMimicOADB.root", Bool_t addOldMaps = kTRUE)
{
    gSystem->Load("libOADB");

    // create working copy of OADB file from EOS
    const char *fileNameOADB                ="EMCalTriggerMimic_temp.root";
    TString fileNameOADBAliLocal="/home/jens/oadb/OADB/PWGGA/EMCalTriggerMimicOADB.root";
    gSystem->Exec(Form("cp %s %s",fileNameOADBAliLocal.Data(),fileNameOADB));

    // update OADB file with dead, bad and warm cells
    // last parameter:
    //      0: write new map for given run range
    //      1: update existing map for given run range
    //      2: same as 1 but also save over/underhead parts (special setting, not for normal use!)

    if(addOldMaps){
                                                                      //L0,    dL0, L1high,dL1high,  L1low, dL1low, PHOSL0,    dL0
        // updateFile(fileNameOADB,"11a_7_1",  "",   144871,   145287,   2.11,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00);
        // updateFile(fileNameOADB,"11a_7_2",  "",   145288,   146374,   3.43,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00);
        // updateFile(fileNameOADB,"11a_7_3",  "",   146375,   146381,   1.71,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00);
        // updateFile(fileNameOADB,"11a_7_4",  "",   146382,   146501,   2.05,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00);
        // updateFile(fileNameOADB,"11a_2_1",  "",   146502,   148521,   3.43,   0.65,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00);
        // updateFile(fileNameOADB,"11b_1",    "",   148522,   150208,   1.94,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00);
        // updateFile(fileNameOADB,"11b_2",    "",   150209,   153055,   3.39,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00);
        // updateFile(fileNameOADB,"11bc_1",   "",   153056,   153910,   4.01,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00);
        // updateFile(fileNameOADB,"11c_1",    "",   153911,   153914,   5.25,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00);
        // updateFile(fileNameOADB,"11cd_1",   "",   153915,   158134,   5.50,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00);
        // updateFile(fileNameOADB,"11d_1",    "",   158135,   158135,   2.05,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00);
        // updateFile(fileNameOADB,"11d_2",    "",   158136,   158177,   5.50,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00);
        // updateFile(fileNameOADB,"11d_3",    "",   158178,   158181,   2.05,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00);
        // updateFile(fileNameOADB,"11e_1",    "",   158182,   160682,   5.50,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00);
        // updateFile(fileNameOADB,"11efg_1",  "",   160683,   160763,   2.05,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00);
        // updateFile(fileNameOADB,"11efg_2",  "",   160764,   161138,   1.71,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00);
        // updateFile(fileNameOADB,"11efg_3",  "",   161139,   161255,   5.50,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00);
        // updateFile(fileNameOADB,"11efg_4",  "",   161256,   161378,   1.71,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00);
        // updateFile(fileNameOADB,"11efg_5",  "",   161379,   161456,   5.50,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00);
        // updateFile(fileNameOADB,"11efg_6",  "",   161457,   161524,   1.71,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00);
        // updateFile(fileNameOADB,"11efg_7",  "",   161525,   161555,   5.50,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00);
        // updateFile(fileNameOADB,"11efg_8",  "",   161556,   161557,   1.71,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00);
        // updateFile(fileNameOADB,"11efg_9",  "",   161558,   161608,   5.50,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00);
        // updateFile(fileNameOADB,"11efg_10", "",   161609,   161629,   1.71,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00);
        // updateFile(fileNameOADB,"11efg_11", "",   161630,   161723,   5.50,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00);
        // updateFile(fileNameOADB,"11h_1",    "",   161724,   173730,   1.71,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00);
        // updateFile(fileNameOADB,"12a_1",    "",   173731,   177143,   2.01,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00);
        // updateFile(fileNameOADB,"12a_2",    "",   177144,   177146,   1.75,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00);
        // updateFile(fileNameOADB,"12a_3",    "",   177147,   177652,   1.52,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00);
        // updateFile(fileNameOADB,"12b_1",    "",   177653,   177723,   2.01,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00);
        // updateFile(fileNameOADB,"12b_2",    "",   177724,   178326,   1.52,   0.20,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00);
        // updateFile(fileNameOADB,"12c-i_1",  "",   178327,   179795,   1.85,   0.20,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00);
        // updateFile(fileNameOADB,"12c-i_2",  "",   179796,   195179,   1.85,   0.20,   9.50,   1.00,   2000,   0.10,   0.00,   0.00);
        // updateFile(fileNameOADB,"13b-f_1",  "",   195180,   197468,   3.20,   0.10,  11.50,   0.50,   7.20,   0.30,   0.00,   0.00);
        // updateFile(fileNameOADB,"13g_1",    "",   197469,   197692,   1.80,   0.12,   5.50,   0.60,   3.75,   0.25,   0.00,   0.00);
        // updateFile(fileNameOADB,"15i-h_1",  "",   235195,   244284,   1.80,   0.10,   2000,   1.20,   1.80,   0.10,   0.00,   0.00);
        // updateFile(fileNameOADB,"15i-m_1",  "",   244285,   244628,   5.00,   0.10,   5.00,   0.80,   5.00,   0.10,   0.00,   0.00);
        // updateFile(fileNameOADB,"15n_1",    "",   244629,   245140,   5.00,   0.10,   5.00,   0.00,   5.00,   0.10,   0.00,   0.00);
        // updateFile(fileNameOADB,"15o_1",    "",   245141,   255538,   1.00,   0.10,  10.00,   1.00,   2000,   0.10,   0.00,   0.00);
        // updateFile(fileNameOADB,"16i-k_1",  "f_out.root",   255539,   258882,   2.20,   0.20,   8.50,   0.00,   3.75,   0.30,   0.00,   0.00);
        // updateFile(fileNameOADB,"16l_1",    "",   258883,   260215,   2.20,   0.20,   5.50,   1.00,   3.80,   0.30,   0.00,   0.00);
        // updateFile(fileNameOADB,"16m-p_1",  "f_out.root",   260216,   265014,   2.50,   0.10,   8.50,   0.70,   3.75,   0.30,   0.00,   0.00);
        // updateFile(fileNameOADB,"16q_1",    "",   265015,   265308,   2.50,   0.10,   8.80,   0.80,   3.90,   0.30,   0.00,   0.00);
        // updateFile(fileNameOADB,"16q_2",    "",   265309,   265588,   2.50,   0.10,  10.80,   0.70,   6.30,   0.40,   0.00,   0.00);
        // updateFile(fileNameOADB,"16r_1",    "",   265589,   265784,   2.50,   0.10,   7.70,   1.00,   4.80,   0.40,   0.00,   0.00);
        // updateFile(fileNameOADB,"16r_2",    "",   265785,   266404,   2.50,   0.10,   6.90,   0.80,   4.80,   0.40,   0.00,   0.00);
        // updateFile(fileNameOADB,"16s_1",    "",   266405,   267160,   3.50,   0.10,   6.90,   0.60,   4.70,   0.35,   0.00,   0.00);
        // updateFile(fileNameOADB,"16t_1",    "",   267161,   270530,   2.50,   0.10,   7.80,   0.90,   5.30,   0.30,   0.00,   0.00);
        // updateFile(fileNameOADB,"17c-o_1",  "f_out.root",   270531,   282024,   2.50,   0.10,   8.50,   0.70,   3.75,   0.30,   0.00,   0.00);
        // updateFile(fileNameOADB,"17pq_1",   "",   282025,   282503,   2.50,   0.10,   8.80,   1.00,   3.90,   0.30,   0.00,   0.00);
        // updateFile(fileNameOADB,"17r_1",    "f_out.root",   282504,   295231,   2.50,   0.10,   8.50,   0.70,   3.75,   0.30,   0.00,   0.00);
        // updateFile(fileNameOADB,"18bp_1",   "f_out.root",   295232,   295273,   2.50,   0.10,   8.50,   0.70,   3.75,   0.30,   0.00,   0.00);
        // updateFile(fileNameOADB,"18qr_1",   "",   295274,   999999,   1.00,   0.10,   9.50,   1.00,   5.00,   0.40,   0.00,   0.00);
        //=====================================================PHOS-13-TeV=======================================================================
        //-------------------------------------------------------LHC-2016------------------------------------------------------------------------

        //updateFile(fileNameOADB,"16e-h_1",  "",   253478,   255467,   -2.20,   0.20,   -8.50,   0.00,   -3.75,   0.30,   3.45,   0.75);//Dont Use! They are included in 15o
        updateFile(fileNameOADB,"15o_1",    "",   245141,   255538,   -1.00,   0.10,  -10.00,   1.00,   -2000,   0.10,   3.7,   0.6);
        updateFile(fileNameOADB,"16i-k_1",  "",   255539,   258882,   -2.20,   0.20,   -8.50,   0.00,   -3.75,   0.30,   3.7,   0.6);
        updateFile(fileNameOADB,"16l_1",    "",   258883,   260215,   -2.20,   0.20,   -5.50,   1.00,   -3.80,   0.30,   3.7,   0.6);
        updateFile(fileNameOADB,"16m-p_1",  "",   260216,   265014,   -2.50,   0.10,   -8.50,   0.70,   -3.75,   0.30,   3.7,   0.6);
        //-------------------------------------------------------LHC-2017------------------------------------------------------------------------
        updateFile(fileNameOADB,"17c-o_1",  "",   270531,   282024,   -2.50,   0.10,   -8.50,   0.70,   -3.75,   0.30,   3.7,   0.6);
        updateFile(fileNameOADB,"17pq_1",   "",   282025,   282503,   -2.50,   0.10,   -8.80,   1.00,   -3.90,   0.30,   3.7,   0.6);
        updateFile(fileNameOADB,"17r_1",    "",   282504,   295231,   -2.50,   0.10,   -8.50,   0.70,   -3.75,   0.30,   3.7,   0.6);
        //-------------------------------------------------------LHC-2018------------------------------------------------------------------------
        updateFile(fileNameOADB,"18bp_1",   "",   295232,   295273,   -2.50,   0.10,   -8.50,   0.70,   -3.75,   0.30,   3.7,   0.6);
    }


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
void updateFile(const char *fileNameOADB,TString arrName, TString fileNameInput,Int_t lowRun, Int_t highRun, Double_t EMCL0, Double_t EMCdL0, Double_t EMCL1high, Double_t EMCdL1high, Double_t EMCL1low, Double_t EMCdL1low, Double_t PHOSL0, Double_t PHOSdL0){
    printf("\n\nAdding %s maps to OADB file:\n",arrName.Data());
    Bool_t bKeepOldOADB=kTRUE;
    TString arrName_LowOld=Form("%s_additonalRuns_Low", arrName.Data());
    TString arrName_HighOld=Form("%s_additonalRuns_High", arrName.Data());

    const int nTrigger=4;
    int nTriggerStart=0;

    // load new input file
    TFile *fInput = NULL;
    if(fileNameInput.CompareTo("")){
        fInput                                   = TFile::Open(fileNameInput);
        if(fInput->IsZombie()){
            cout << "new input file " << fileNameInput.Data() << " not found! returning..." << endl;
            return;
        }
    }
    // load OADB file
    AliOADBContainer *con = NULL;
    TFile *foadb                                    = new TFile(fileNameOADB);
    if(!foadb->IsZombie()){
        foadb->ls();
        con                       =(AliOADBContainer*)foadb->Get("AliEMCalTriggerMimic");
        if(con) con->SetName("Old");
    } else {
        cout << "previous OADB file not found... this is fine if it is the first time running...." << endl;
    }

    // make the output container
    AliOADBContainer *con2                      = new AliOADBContainer("AliEMCalTriggerMimic");

    // List of brand new array
    TObjArray arrayAdd(nTrigger);
    TObjArray arrayAdd_LowOld(nTrigger);
    TObjArray arrayAdd_HighOld(nTrigger);
    arrayAdd.SetName(arrName.Data());
    arrayAdd_LowOld.SetName(arrName_LowOld.Data());
    arrayAdd_HighOld.SetName(arrName_HighOld.Data());

    TH1S *hEMCalTriggerMimic[nTrigger];
    TH1S *hEMCalTriggerMimic_Old[nTrigger];
    Bool_t doNewTrigger[nTrigger];
    TObjArray* TriggerMimicObject_Old;
    Double_t LowerLimitForNewMimic=-0.00001;

    // initialize histograms in array - either load existing map or create new, empty map
    TString TriggNames[4] = {"EMCalL0","EMCalL1G2","EMCalL1G1", "PHOSL0"};
    Double_t TriggerValuesThresh[4]={EMCL0, EMCL1low, EMCL1high, PHOSL0};
    Double_t TriggerValueSigma[4]={EMCdL0, EMCdL1low, EMCdL1high, PHOSdL0};
    for(int i=nTriggerStart;i<nTrigger;i++){
        hEMCalTriggerMimic[i] = NULL;
        doNewTrigger[i]=kFALSE;
        if( (!fileNameInput.CompareTo(""))||(i==3) ){
            // make new TH1S and set values
            hEMCalTriggerMimic[i] = new TH1S(TriggNames[i].Data(),TriggNames[i].Data(),1,0,1);
            hEMCalTriggerMimic[i]->SetBinContent(1,TriggerValuesThresh[i]*100);
            hEMCalTriggerMimic[i]->SetBinError(1,TriggerValueSigma[i]*100);
            if (TriggerValuesThresh[i]>=LowerLimitForNewMimic){
                doNewTrigger[i]=kTRUE;
            }
        } else {
            // load SM-wise inputs from .root file
            hEMCalTriggerMimic[i] = (TH1S*)fInput->Get(Form("%s",TriggNames[i].Data()));
            hEMCalTriggerMimic[i]->SetName(Form("%s",TriggNames[i].Data()));
            doNewTrigger[i]=kTRUE;
            if(!hEMCalTriggerMimic[i]){
                cout << "Input histogram not found for bin " << i << "! returning...." <<  endl;
                return;
            }
        }
    }

    // add histograms to a new array
    for (Int_t i=nTriggerStart;i<nTrigger;i++){
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
                TriggerMimicObject_Old=(TObjArray*)con->GetObjectByIndex(i);
            } else if (highRun >= con->LowerLimit(i) && highRun <= con->UpperLimit(i)){
                printf("\n!!! Not adding index %d for runrange %d--%d as high run limit %d is contained\n", i, con->LowerLimit(i),con->UpperLimit(i),highRun);
                replacedContainerLowLimit           = con->LowerLimit(i);
                replacedContainerHighLimit          = con->UpperLimit(i);
                TriggerMimicObject_Old=(TObjArray*)con->GetObjectByIndex(i);

            } else if ((con->LowerLimit(i)==con->UpperLimit(i))&& con->LowerLimit(i) >= lowRun && con->UpperLimit(i)<=highRun){
                printf("\n!!! Not adding index %d for runrange %d--%d as high run limit %d is contained\n", i, con->LowerLimit(i),con->UpperLimit(i),highRun);
                TriggerMimicObject_Old=(TObjArray*)con->GetObjectByIndex(i);
            }else{
                con2->AddDefaultObject(con->GetObjectByIndex(i));
                con2->AppendObject(con->GetObjectByIndex(i),con->LowerLimit(i),con->UpperLimit(i));
            }
        }
    } else {
        cout << "old container not found... is this the first time running???" << endl;
    }
    for (Int_t j=nTriggerStart;j<nTrigger;j++){
        hEMCalTriggerMimic_Old[j]=NULL;
        hEMCalTriggerMimic_Old[j]=(TH1S*)TriggerMimicObject_Old->FindObject(TriggNames[j]);
        if (hEMCalTriggerMimic_Old[j]==NULL){
            hEMCalTriggerMimic_Old[j]= new TH1S(Form("%s_Old",TriggNames[j].Data()),Form("%s",TriggNames[j].Data()),1,0,1);
            hEMCalTriggerMimic_Old[j]->SetBinContent(1,0.*100);
            hEMCalTriggerMimic_Old[j]->SetBinError(1,0.*100);
        }
        if (doNewTrigger[j]==kFALSE){
            arrayAdd.AddAt(hEMCalTriggerMimic_Old[j],j);
        }
        arrayAdd_LowOld.AddAt(hEMCalTriggerMimic_Old[j],j);
        arrayAdd_HighOld.AddAt(hEMCalTriggerMimic_Old[j],j);
    }
    // add new BC maps at the end of the new OADB file (will be sorted later)
    con2->AddDefaultObject(&arrayAdd);
    con2->AppendObject(&arrayAdd, lowRun,highRun);
    if (bKeepOldOADB==kTRUE){
        if (lowRun>replacedContainerLowLimit){
            con2->AddDefaultObject(&arrayAdd_LowOld);
            con2->AppendObject(&arrayAdd_LowOld, replacedContainerLowLimit,lowRun-1);
        }
        if (highRun<replacedContainerHighLimit){
            con2->AddDefaultObject(&arrayAdd_HighOld);
            con2->AppendObject(&arrayAdd_HighOld, highRun+1,replacedContainerHighLimit);
        }
    }
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
void rebuildContainer(const char *fileNameOADB){
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
