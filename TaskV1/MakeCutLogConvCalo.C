// provided by Gamma Conversion Group, PWG4, Kathrin Koch, kkoch@physi.uni-heidelberg.de

#include <Riostream>
#include <fstream>
using namespace std;

void MakeCutLogConvCalo(const char *inputRootFile = "GammaConvCalo.root", const char *InputName){

  TString filename = inputRootFile; 
  TFile *file = new TFile(filename.Data());  
  if (file->IsZombie()) return;
   
  fstream outputFile(InputName,ios::out);
  if(!outputFile.is_open()){
    cout<<"Problem opening file"<<endl;
    return;
  }

  //  Char_t filename_input1[200] = (Form("%s%s",path,input1));	
  TKey *key;
  file->ls();

  cout<<file<<endl;
  //TIter nextkey(l.GetListOfKeys());
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

}
