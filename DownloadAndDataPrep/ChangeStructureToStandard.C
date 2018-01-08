void ChangeStructureToStandard(TString nameInputFile, TString namefileOutput, TString nameInputList, Int_t mode == 0){

   TFile *fileInput = NULL;
   fileInput = new TFile(nameInputFile.Data());
   cout << fileInput << endl;

   TList *listInput =(TList*)fileInput->Get(nameInputList.Data());
   if (listInput == NULL){
      return;
   }
   TFile *fileOutput = new TFile(namefileOutput,"RECREATE");

   TString nominalMainDir     = "";
   if (mode == 9 || mode == 0)
       nominalMainDir         = "GammaConvV1";
   else if( mode == 1 )
       nominalMainDir         = "GammaConvDalitzV1";
   else if (mode == 2 || mode == 3 || mode == 13)
       nominalMainDir         = "GammaConvCalo";
   else if (mode == 4 || mode == 12 || mode == 5)
       nominalMainDir         = "GammaCalo";
   else if( mode == 6 || mode == 7 )
       nominalMainDir         = "GammaConvDalitzCalo";
   else if (mode == 10 || mode == 11 )
       nominalMainDir         = "GammaCaloMerged";
   else if (mode == 30 )
       nominalMainDir         = "GammaConvV1";
   else if (mode == 40 || mode == 41 || mode == 42 || mode == 43|| mode == 44 || mode == 45 ||
       mode == 46 || mode == 47 || mode == 48 || mode == 49|| mode == 50)
       nominalMainDir         = "GammaConvNeutralMesonPiPlPiMiPiZero";

   TList *listOutput =(TList*)fileOutput->Get(nominalMainDir.Data());
   Bool_t kNewList = kFALSE;
   if (!listOutput){
      kNewList = kTRUE;
      listOutput = new TList();
      listOutput->SetName(nominalMainDir.Data());
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
