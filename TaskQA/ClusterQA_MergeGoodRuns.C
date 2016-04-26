#include "QA.h"

void ChangeStrucToStd(TString nameInputFile, TString namefileOutput, TString nameInputList);
void MakeCLog(const char *inputRootFile = "GammaConvCalo.root", const char *InputName = "");

void ClusterQA_MergeGoodRuns(){

	TString folder = "/home/daniel/data/work/photonconv/AnalysisSoftware";

//    TString Tag = "20160125";
//    const Int_t nMergedFiles = 1;
//    TString MergeFile[nMergedFiles] = {"GammaConvCalo_120"};
//    const Int_t nSets = 21;
//    TString DataSets[nSets]={"LHC12a", "LHC12b", "LHC12c", "LHC12d", "LHC12f", /*"LHC12g",*/ "LHC12h", "LHC12i",
//                             "LHC15h1a1", "LHC15h1b", "LHC15h1c", "LHC15h1d", "LHC15h1f", /*"LHC15h1g",*/ "LHC15h1h", "LHC15h1i",
//                             "LHC15h2a", "LHC15h2b", "LHC15h2c", "LHC15h2d", "LHC15h2f", /*"LHC15h2g",*/ "LHC15h2h", "LHC15h2i"
//                            };
//    // merging periods
//    const Int_t nMerge = 3;
//    Int_t nMergeSets[nMerge] = {6,13,20};
//    TString MergeSets[nMerge] = {"LHC12","LHC15h1","LHC15h2"};


//    TString Tag = "20160309";
//    const Int_t nMergedFiles = 1;
//    TString MergeFile[nMergedFiles] = {"GammaConvCalo_122"};
//    const Int_t nSets = 5;
//    TString DataSets[nSets]={"LHC12c-kEMCEGA", "LHC12d-kEMCEGA", "LHC12f-kEMCEGA", /*"LHC12g-kEMCEGA",*/ "LHC12h-kEMCEGA", "LHC12i-kEMCEGA"};
//    // merging periods
//    const Int_t nMerge = 1;
//    Int_t nMergeSets[nMerge] = {4};
//    TString MergeSets[nMerge] = {"LHC12-kEMCEGA"};

    TString Tag = "20160309";
    const Int_t nMergedFiles = 1;
    TString MergeFile[nMergedFiles] = {"GammaConvCalo_121"};
    const Int_t nSets = 6;
    TString DataSets[nSets]={"LHC12b-kEMC7", "LHC12c-kEMC7", "LHC12d-kEMC7", "LHC12f-kEMC7", /*"LHC12g-kEMC7",*/ "LHC12h-kEMC7", "LHC12i-kEMC7"};
    // merging periods
    const Int_t nMerge = 1;
    Int_t nMergeSets[nMerge] = {5};
    TString MergeSets[nMerge] = {"LHC12-kEMC7"};

//    TString Tag = "20160125";
//    const Int_t nMergedFiles = 1;
//    TString MergeFile[nMergedFiles] = {"GammaConvCalo_120"};
//    const Int_t nSets = 3;
//    TString DataSets[nSets]={"LHC12b","LHC15h1b","LHC15h2b"};
//    // merging periods
//    const Int_t nMerge = 0;
//    Int_t nMergeSets[nMerge] = {};
//    TString MergeSets[nMerge] = {};

	std::vector<TString> vecRuns;
	TString fileRuns;
	TString fDataSet;

    TString mergePeriod[3][nMergedFiles];
    for(Int_t iFiles=0; iFiles<nMergedFiles; iFiles++){
      for(Int_t iM=0; iM<nMerge; iM++){
        TString mergeP = Form("%s/DataQA/%s/%s_%s",folder.Data(), Tag.Data(), MergeSets[iM].Data(), MergeFile[iFiles].Data());
        mergePeriod[iM][iFiles] = Form("hadd -f -k %s.root",mergeP.Data());
      }
    }

	for(Int_t iSet=0; iSet<nSets; iSet++)
	{
		vecRuns.clear();
		fDataSet = DataSets[iSet];
		fileRuns = Form("%s/runNumbers%s.txt", folder.Data(), fDataSet.Data());

		cout << "\n------------------------------------------------------" << endl;
		if(!readin(fileRuns, vecRuns)) cout << "\n\n\n**********************ERROR, no Run Numbers could be found!**********************\n\n\n" << endl;
		cout << "------------------------------------------------------" << endl;

        for(Int_t iFiles=0; iFiles<nMergedFiles; iFiles++){
          TString mergeR = Form("%s/DataQA/%s/%s_%s",folder.Data(), Tag.Data(), fDataSet.Data(), MergeFile[iFiles].Data());
          TString mergeRuns = Form("hadd -f -k %s.root",mergeR.Data());
          for(Int_t iRun=0; iRun<(Int_t)vecRuns.size(); iRun++){
              mergeRuns+=Form(" %s/DataQA/%s/%s/%s/%s.root", folder.Data(), Tag.Data(), fDataSet.Data(), (vecRuns.at(iRun)).Data(), MergeFile[iFiles].Data());
          }

          cout << "Merging " << vecRuns.size() << " runs from period" << fDataSet.Data() << "...";
          gSystem->Exec(mergeRuns.Data());
          //cout << mergeRuns.Data() << endl;
          ChangeStrucToStd(Form("%s.root",mergeR.Data()),Form("%s.root",mergeR.Data()),MergeFile[iFiles].Data());
          cout << "done!" << endl;

          Int_t kMerge=kTRUE;
          for(Int_t iMerge=0; iMerge<nMerge; iMerge++){
            if(kMerge && iSet<=nMergeSets[iMerge]){
              mergePeriod[iMerge][iFiles]+=Form(" %s/DataQA/%s/%s_%s.root",folder.Data(), Tag.Data(), fDataSet.Data(), MergeFile[iFiles].Data());
              kMerge=kFALSE;
            }
          }
        }
	}

    cout << "\n------------------------------------------------------" << endl;
    cout << "Merging Periods: " << endl;
    cout << "------------------------------------------------------\n" << endl;

    for(Int_t iFiles=0; iFiles<nMergedFiles; iFiles++){
      for(Int_t iM=0; iM<nMerge; iM++){
        cout << "Merging " << MergeSets[iM].Data() << " ..." << endl;
        gSystem->Exec(mergePeriod[iM][iFiles].Data());
        //cout << mergePeriod[iM][iFiles].Data() << endl;
        cout << "done!" << endl;
      }
    }

	return;
}

void ChangeStrucToStd(TString nameInputFile, TString namefileOutput, TString nameInputList){

   TFile *fileInput = new TFile(nameInputFile.Data());
   cout << fileInput << endl;

   TList *listInput =(TList*)fileInput->Get(nameInputList.Data());
   if (listInput == NULL){
	  return;
   }else listInput->SetOwner();

   TObjArray *rArr = nameInputList.Tokenize("_");
   TObjString* rString = (TObjString*)rArr->At(0);
   TString string = rString->GetString();

   TFile *fileOutput = new TFile(namefileOutput,"RECREATE");
   TList *listOutput =(TList*)fileOutput->Get(string.Data());
   Bool_t kNewList = kFALSE;
   if (!listOutput){
	  kNewList = kTRUE;
	  listOutput = new TList();
      listOutput->SetName(string.Data());
   }

   for(Int_t i = 0; i<listInput->GetSize(); i++){
	  TList *listToSave = (TList*)listInput->At(i);
	  TString dirname = listToSave->GetName();
	  cout<<dirname<<endl;
	  if(listToSave){
		 cout<<"found"<<endl;
		 listOutput->Add(listToSave);
	  }
   }

   listOutput->Write("",TObject::kSingleKey);
   delete listOutput;
   fileOutput->Close();
   delete fileOutput;

   delete listInput;
   fileInput->Close();
   delete fileInput;

   return;
}

void MakeCLog(const char *inputRootFile, const char *InputName){

  TString filename = inputRootFile;
  TFile *file = new TFile(filename.Data());
  if (file->IsZombie()) return;

  fstream outputFile(InputName,ios::out);
  if(!outputFile.is_open()){
	cout<<"Problem opening file"<<endl;
	return;
  }

  file->ls();

  cout<<file<<endl;
  TList *list = (TList*) file->Get("GammaConvCalo");

  for(Int_t i = 0; i<list->GetEntries(); i++){
	 TList *l2 = (TList*) list->At(i);
	 TString dirname = l2->GetName();
	 if(dirname.BeginsWith("Cut") == 1){
		TString CutNumber(dirname(11,dirname.Length()-11));
		outputFile << CutNumber.Data() <<endl;
		cout<<CutNumber<<endl;
	 }
  }
  outputFile.close();

  file->Close();
  delete file;

}
