// provided by Gamma Conversion Group, PWG4, Kathrin Koch, kkoch@physi.uni-heidelberg.de

#include <Riostream>
#include <fstream>
using namespace std;

void MakeCutLogDalitz(const char *inputRootFile = "AnalysisResults.root", const char *InputName,Int_t mode = 1){

  fstream outputFile(InputName,ios::out);
  if(!outputFile.is_open()){
    cout<<"Problem opening file"<<endl;
    return;
  }


  TString filename = inputRootFile;	
  TFile f(filename.Data());  
  
  TString nameOutputDir = "GammaConvDalitzV1";
  
  cout<<"mode: "<<mode<<endl;
  
  if( mode ==  6 || mode == 7 ){
  
    //nameOutputDir = "GammaCaloDalitz";
    nameOutputDir = "GammaConvDalitzCalo";
  
  }

 TKey *key;
 TIter nextkey(f.GetListOfKeys());
 
  while ((key = (TKey*)nextkey())) {

     const char *classname = key->GetClassName();
     TClass *cl = gROOT->GetClass(classname);

     if (!cl) continue;
    
        if (cl->InheritsFrom(TList::Class())){

            TList* dir = (TList*)key->ReadObj();
            TString listName = dir->GetName();
            if(  listName.Contains(nameOutputDir.Data()) ) { 

            TObject * obj;

            TIter next(dir);

            while ( ( obj = next() ) ){

                TString dirname = obj->GetName();
                    if(dirname.BeginsWith("Cut Number") ){
                        TString CutNumber(dirname(11,100));
                        cout<<CutNumber<<endl;
                        outputFile << CutNumber.Data() <<endl;
                    }
            }
           }
        }
  }
 outputFile.close();
}

  // TList *directories = f.GetListOfKeys(); // get the list of directories in the file
  
  // for(Int_t entFile=0;entFile<directories->GetEntries();entFile++){

  //    TObject * o = f.Get(directories->At(entFile)->GetName()); // get the object in the base directory
     
  //   // if(TString(o->IsA()->GetName())=="TDirectoryFile"){ // means that this is a directory (PWG4......)
      
  //   //   TDirectory *pwg4dir =(TDirectory*)o;
 
  //   //   TString baseDirName = pwg4dir->GetName();
      
  //   //   TString reconstructionFlagString = ""; // this is for new scheme where also the flags are coded in numbers in the PWG4.... name

  //   //   if(baseDirName.Length()>33){
  //   // 	reconstructionFlagString = baseDirName(baseDirName.Index("GammaConversion_")+16,8);
  //   //   }
      
  //   //   TList *pwg4list = pwg4dir->GetListOfKeys(); // list of the yeys inside the base directory

  //   //   for(Int_t entHist=0;entHist<pwg4list->GetEntries();entHist++){
  //   // 	TString name = pwg4list->At(entHist)->GetName();
  //   // 	cout<<name<<endl;

  //   // 	// if(name.Contains("container")==0){ // does not try to read the container (get errors if tried)
  //   // 	    TObject * oHist = pwg4dir->Get(pwg4list->At(entHist)->GetName()); // get the object 
	  
  //   // 	    if(TString(oHist->IsA()->GetName())=="TList"){ // check if the object is a TList
	    
  //   // 	    TString listname = oHist->GetName();
  //   // 	    cout<<"Reading: "<<listname.Data()<<endl;
	    
  //   // 	    TString cutString = listname(listname.Index("_")+1,listname.Length()) + "\n";// get the Cut string from the name


  //   // 	    outputFile << cutString.Data();
  //   // 	  }
  //   // 	}
  //   //   }
  //   }
  // }
  // outputFile.close();
// }
