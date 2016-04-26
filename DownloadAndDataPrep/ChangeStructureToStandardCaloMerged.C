void ChangeStructureToStandardCaloMerged(TString nameInputFile, TString namefileOutput, TString nameInputList){

   TFile *fileInput = NULL;
   fileInput = new TFile(nameInputFile.Data());
   cout << fileInput << endl;
  
   TList *listInput =(TList*)fileInput->Get(nameInputList.Data());
   if (listInput == NULL){ 
      return;
   }   
   TFile *fileOutput = new TFile(namefileOutput,"RECREATE");
   TList *listOutput =(TList*)fileOutput->Get("GammaCaloMerged");
   Bool_t kNewList = kFALSE;
   if (!listOutput){
      kNewList = kTRUE;
      listOutput = new TList();
      listOutput->SetName("GammaCaloMerged");
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
   fileOutput->Close();
}
