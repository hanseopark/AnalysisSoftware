#include <TString.h>
#include <TFile.h>
#include <TObject.h>
#include <TSystem.h>
#include <TGrid.h>

#include <vector>
#include <iostream>
#include <fstream>

Bool_t copyAlien2Local(TString where, TString loc, TString rem);

void Grid_FailedMergeDLRun(TString folder = "/home/daniel/data/work/photonconv/AnalysisSoftware")
{
    cout<<"Connecting to Alien..."<<endl;
    TGrid::Connect("alien://");
    cout<<"==============================="<<endl;
    cout<<"Successfully connected to Alien"<<endl;
    cout<<"==============================="<<endl;

	//LHC11a
//	const Int_t nFiles = 1;
//	TString Tag = "20150720";
//	TString DataSetsFile[nFiles] = {"GammaConvCalo_1.root"};

//	TString DataSet = "LHC11a_p4_wSDD";
//	TString prefix = "/alice/data/2011/LHC11a/000";
//	TString suffix = "/ESDs/pass4_with_SDD/PWGGA/GA_pp/779_20150721-1421/";
//	Int_t nRuns = 1;
//	Int_t run[nRuns] = {146804};
//	Int_t first[nRuns] = {1};
//	Int_t last[nRuns] = {130};

// 8TeV MC
//	const Int_t nFiles = 1;
//	TString Tag = "20150728";
//	TString DataSetsFile[nFiles] = {"GammaConvCalo_120.root"};

//	TString prefix = "/alice/sim/2014/LHC14e2c/";
//	TString suffix = "/PWGGA/GA_pp_MC/800_20150728-0854/";
//	TString DataSet = "LHC14e2c";
//	const Int_t nRuns = 5;
//	Int_t run[nRuns] = {177597,180201,180564,184371,185217};
//	Int_t first[nRuns] = {1,1,1,1,1};
//	Int_t last[nRuns] = {81,21,21,26,54};

//    const Int_t nFiles = 3;
//    TString Tag = "20160104";
//    TString DataSetsFile[nFiles] = {"GammaCalo_110.root", "GammaConvCalo_110.root", "GammaConvCalo_120.root"};

//    TString prefix = "/alice/data/2012/LHC12h/000";
//    TString suffix = "/pass2/PWGGA/GA_pp/1152_20160105-0012/";
//    TString DataSet = "LHC12h";
//    const Int_t nRuns = 2;
//    Int_t run[nRuns] = {190904,192177};
//    Int_t first[nRuns] = {1,1};
//    Int_t last[nRuns] = {29,75};

//    const Int_t nFiles = 2;
//    TString Tag = "20160104";
//    TString DataSetsFile[nFiles] = {"AnalysisResults.root","GammaConvV1_70.root"};

//    TString prefix = "/alice/data/2012/LHC12h/000";
//    TString suffix = "/pass2/PWGGA/GA_pp/1144_20160104-2355/";
//    TString DataSet = "LHC12h";
//    const Int_t nRuns = 3;
//    Int_t run[nRuns] = {189616,190904,192140};
//    Int_t first[nRuns] = {1,1,1};
//    Int_t last[nRuns] = {48,29,27};

//    const Int_t nFiles = 2;
//    TString Tag = "20160104";
//    TString DataSetsFile[nFiles] = {"AnalysisResults.root","GammaConvV1_70.root"};

//    TString prefix = "/alice/data/2012/LHC12a/000";
//    TString suffix = "/pass2/PWGGA/GA_pp/1138_20160104-2353/";
//    TString DataSet = "LHC12a";
//    const Int_t nRuns = 1;
//    Int_t run[nRuns] = {176926};
//    Int_t first[nRuns] = {1};
//    Int_t last[nRuns] = {32};


//        const Int_t nFiles = 3;
//        TString Tag = "20160104";
//        TString DataSetsFile[nFiles] = {"GammaCalo_110.root", "GammaConvCalo_110.root", "GammaConvCalo_120.root"};

//        TString prefix = "/alice/sim/2015/LHC15h1b/";
//        TString suffix = "/PWGGA/GA_pp_MC/1438_20160105-0050/";
//        TString DataSet = "LHC15h1b";
//        const Int_t nRuns = 1;
//        Int_t run[nRuns] = {177864};
//        Int_t first[nRuns] = {1};
//        Int_t last[nRuns] = {33};

//        const Int_t nFiles = 2;
//        TString Tag = "20160104";
//        TString DataSetsFile[nFiles] = {"AnalysisResults.root","GammaConvV1_70.root"};

//        TString prefix = "/alice/sim/2015/LHC15h2i/";
//        TString suffix = "/PWGGA/GA_pp_MC/1436_20160105-0033/";
//        TString DataSet = "LHC15h2i";
//        const Int_t nRuns = 1;
//        Int_t run[nRuns] = {193007};
//        Int_t first[nRuns] = {1};
//        Int_t last[nRuns] = {23};

//        const Int_t nFiles = 2;
//        TString Tag = "20160104";
//        TString DataSetsFile[nFiles] = {"AnalysisResults.root","GammaConvV1_70.root"};

//        TString prefix = "/alice/sim/2015/LHC15h2f/";
//        TString suffix = "/PWGGA/GA_pp_MC/1433_20160105-0033/";
//        TString DataSet = "LHC15h2f";
//        const Int_t nRuns = 1;
//        Int_t run[nRuns] = {187633};
//        Int_t first[nRuns] = {1};
//        Int_t last[nRuns] = {44};

//    const Int_t nFiles = 2;
//    TString Tag = "20160104";
//    TString DataSetsFile[nFiles] = {"AnalysisResults.root","GammaConvV1_70.root"};

//    TString prefix = "/alice/sim/2015/LHC15h2d/";
//    TString suffix = "/PWGGA/GA_pp_MC/1432_20160105-0032/";
//    TString DataSet = "LHC15h2d";
//    const Int_t nRuns = 2;
//    Int_t run[nRuns] = {183935,183933};
//    Int_t first[nRuns] = {1,1};
//    Int_t last[nRuns] = {60,52};

//    const Int_t nFiles = 2;
//    TString Tag = "20160104";
//    TString DataSetsFile[nFiles] = {"AnalysisResults.root","GammaConvV1_70.root"};

//    TString prefix = "/alice/sim/2015/LHC15h1d/";
//    TString suffix = "/PWGGA/GA_pp_MC/1424_20160105-0029/";
//    TString DataSet = "LHC15h1d";
//    const Int_t nRuns = 2;
//    Int_t run[nRuns] = {184134,184137};
//    Int_t first[nRuns] = {1,1};
//    Int_t last[nRuns] = {61,55};

//            const Int_t nFiles = 2;
//            TString Tag = "20160104";
//            TString DataSetsFile[nFiles] = {"AnalysisResults.root","GammaConvV1_70.root"};

//            TString prefix = "/alice/sim/2015/LHC15h1b/";
//            TString suffix = "/PWGGA/GA_pp_MC/1422_20160105-0029/";
//            TString DataSet = "LHC15h1b";
//            const Int_t nRuns = 1;
//            Int_t run[nRuns] = {177592};
//            Int_t first[nRuns] = {1};
//            Int_t last[nRuns] = {36};


    //LHC10

    const Int_t nFiles = 2;
    TString Tag = "20160206";
    TString DataSetsFile[nFiles] = {"GammaCalo_201.root","GammaConvCalo_201.root"};

    TString prefix = "/alice/sim/2014/LHC14j4e/";
    TString suffix = "/PWGGA/GA_pp_MC/1616_20160207-1047/";
    TString DataSet = "LHC14j4e";
    const Int_t nRuns = 2;
    Int_t run[nRuns] = {128855,128498};
    Int_t first[nRuns] = {1,1};
    Int_t last[nRuns] = {50,43};

	for(Int_t iRun=0; iRun<nRuns; iRun++)
	{
		TString PathDataSet = Form("%s%i%s",prefix.Data(),run[iRun],suffix.Data());

		TString fPathGrid;
		TString fPathLocal;
		for(Int_t k=0; k<nFiles; k++)
		{
			TString merge = Form("hadd %s/DataQA/%s/%s/%i/%s", folder.Data(), Tag.Data(), DataSet.Data(), run[iRun], DataSetsFile[k].Data());

			for(Int_t j=first[iRun]; j<=last[iRun]; j++)
			{
					fPathGrid = Form("%s%04i/%s", PathDataSet.Data(), j, (DataSetsFile[k]).Data());
					fPathLocal = Form("%s/DataQA/%s/%s/%i", folder.Data(), Tag.Data(), DataSet.Data(), run[iRun]);
					//gSystem->Exec(Form("mkdir -p %s",fPathLocal.Data()));

					fPathLocal+="/"; fPathLocal+=j; fPathLocal+="_"; fPathLocal+=DataSetsFile[k];

					cout << endl;
					cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
					cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
					cout << "Copying from (grid): " << fPathGrid.Data() << endl;
					cout << "Copying to (local): " << fPathLocal.Data() << endl;
					cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
					cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

					TFile fileCheck(fPathLocal.Data());
					if(!fileCheck.IsZombie()) {cout << "\n\t\t>>>>>>>>>>>>>>>>>>Info: ROOT-File |" << fPathLocal.Data() << "| does already exist! Continue...<<<<<<<<<<<<<<\n" << endl; continue;}

					if(copyAlien2Local("",fPathGrid,fPathLocal)) {
						merge+=Form(" %s",fPathLocal.Data());
						continue;
					}
					else cout << "\n\n\t*******************Err: copyAlien2Local(), check runlist in photonconv rep!*******************\n" << endl;
			}
			cout << "Merging using command: " << merge.Data() << endl;
			gSystem->Exec(merge.Data());
			cout << "done!" << endl;

            cout << "Deleting files..." << endl;
            for(Int_t j=first[iRun]; j<=last[iRun]; j++)
            {
              fPathLocal = Form("%s/DataQA/%s/%s/%i", folder.Data(), Tag.Data(), DataSet.Data(), run[iRun]);
              fPathLocal+="/"; fPathLocal+=j; fPathLocal+="_"; fPathLocal+=DataSetsFile[k];
              cout << "... " << Form("rm %s",fPathLocal.Data()) << endl;
              gSystem->Exec(Form("rm %s",fPathLocal.Data()));
            }
            cout << "done!" << endl;
		}
	}
    return;
}

Bool_t copyAlien2Local(TString where, TString loc, TString rem)
{
   TString sl(Form("alien://%s", loc.Data()));
   TString sr(Form("file://%s", rem.Data()));
   Bool_t ret = TFile::Cp(sl,sr);
   if (!ret) {
      cout << Form(where.Data(), "Failed to copy %s to %s", sl.Data(), sr.Data()) << endl;
   }
   return ret;
}
