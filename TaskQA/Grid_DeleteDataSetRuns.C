/*******************************************************************************
 ******  provided by Gamma Conversion Group, PWGGA,                        *****
 ******     Daniel Muehlheim, d.muehlheim@cern.ch                          *****
 *******************************************************************************/

void readin(TString fileRuns, std::vector<TString> &vec);

void Grid_DeleteDataSetRuns(TString folder = "/home/daniel/data/work/pcgGit/AnalysisSoftware")
{
  std::vector<TString> vecRuns;

  TString fDataSet="LHC15h2c";
  TString Tag = "20160104";
  const Int_t nFiles = 2;
//  TString fileDelete[nFiles]={
//    "GammaCalo_110.root",
//    "GammaConvCalo_110.root",
//    "GammaConvCalo_120.root"};
  TString fileDelete[nFiles]={
    "AnalysisResults.root",
    "GammaConvV1_70.root"};

  TString fileRuns = Form("%s/runNumbers%s.txt", folder.Data(), fDataSet.Data());
  cout << "\n------------------------------------------------------" << endl;
  if(!readin(fileRuns, vecRuns)) cout << "\n\n\n**********************ERROR, no Run Numbers could be found!**********************\n\n\n" << endl;
  cout << "------------------------------------------------------" << endl;

  for(Int_t i=0; i<(Int_t)vecRuns.size(); i++){
    for(Int_t j=0; j<nFiles; j++){
      TString fPathLocal = Form("%s/DataQA/%s/%s/%s/%s", folder.Data(), Tag.Data(), fDataSet.Data(), vecRuns.at(i).Data(),fileDelete[j].Data());
      cout << "Deleting " << fPathLocal.Data() << "...";
      gSystem->Exec(Form("rm %s",fPathLocal.Data()));
      cout << "done!" << endl;
    }
  }

  return;
}

Bool_t readin(TString fileRuns, std::vector<TString> &vec){
    cout << Form("Reading from %s...", fileRuns.Data()) << endl;
    fstream file;
    TString fVar;
    Int_t totalN=0;
    file.open(fileRuns.Data(), ios::in);
    if(file.good())
    {
        file.seekg(0L, ios::beg);
        cout << "Processing Runs: \"";
        while(!file.eof())
        {
            file >> fVar;
            if(fVar.Sizeof()>1)
            {
                cout << fVar.Data() << ", ";
                vec.push_back(fVar);
                totalN++;
            }
        }
        cout << "\"" << endl;
    }
    file.close();
    cout << "...done!\n\nIn total " << totalN << " Runs will be processed!" << endl;
    if(totalN > 0) return kTRUE;
    else return kFALSE;
}
