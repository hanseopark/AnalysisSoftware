/*******************************************************************************
 ******  provided by Gamma Conversion Group, PWGGA,                        *****
 ******     Daniel Muehlheim, d.muehlheim@cern.ch                          ***** 
 *******************************************************************************/

#include <TString.h>
#include <TFile.h>
#include <TObject.h>
#include <TSystem.h>
#include <TGrid.h>

#include <vector>
#include <iostream>
#include <fstream>

Bool_t copyAlien2Local(TString where, TString loc, TString rem);

void Grid_FailedMergeDLRun(TString folder = "/home/daniel/data/work/pcgGit/AnalysisSoftware")
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
//    const Int_t nFiles = 1;
//    TString Tag = "20160524";
//    TString DataSetsFile[nFiles] = {"GammaConvCalo_101.root"};

//    TString prefix = "/alice/sim/2015/LHC15h2d/";
//    TString suffix = "/AOD178/PWGGA/GA_pp_MC_AOD/132_20160622-2225/";
//    TString DataSet = "LHC15h2d";
//    const Int_t nRuns = 1;
//    Int_t run[nRuns] = {184147};
//    Int_t first[nRuns] = {1};
//    Int_t last[nRuns] = {39};


//    const Int_t nFiles = 1;
//    TString Tag = "20160524";
//    TString DataSetsFile[nFiles] = {"GammaConvCalo_101.root"};

//    TString prefix = "/alice/sim/2015/LHC15h1c/";
//    TString suffix = "/AOD178/PWGGA/GA_pp_MC_AOD/129_20160621-2332/";
//    TString DataSet = "LHC15h1c";
//    const Int_t nRuns = 1;
//    Int_t run[nRuns] = {179571};
//    Int_t first[nRuns] = {1};
//    Int_t last[nRuns] = {68};

//    const Int_t nFiles = 1;
//    TString Tag = "20160524";
//    TString DataSetsFile[nFiles] = {"GammaConvCalo_101.root"};

//    TString prefix = "/alice/sim/2015/LHC15h1h/";
//    TString suffix = "/AOD178/PWGGA/GA_pp_MC_AOD/131_20160622-2225/";
//    TString DataSet = "LHC15h1h";
//    const Int_t nRuns = 1;
//    Int_t run[nRuns] = {189310};
//    Int_t first[nRuns] = {1};
//    Int_t last[nRuns] = {24};

//    const Int_t nFiles = 2;
//    TString Tag = "20161127";
//    TString DataSetsFile[nFiles] = {"GammaConvCalo_100.root","GammaConvCalo_101.root"};

//    TString prefix = "/alice/sim/2016/LHC16c2/";
//    TString suffix = "/PWGGA/GA_pp_MC/2679_20161124-1407/";
//    TString DataSet = "LHC16c2";
//    const Int_t nRuns = 1;
//    TString run[nRuns] = {/*"10/193051","20/192349","6/192349","15/192073","12/187488",*/"7/187488"/*,"1/187488","5/185687","10/184215","11/182692","19/182692","1/193051"*/};
//    Int_t first[nRuns] = {/*39,42,49,47,33,*/43/*,45,22,54,26,26,17*/};
//    Int_t last[nRuns] = {/*43,47,53,52,39,*/47/*,49,27,58,30,31,21*/};
//    const Int_t nFiles = 2;
//    TString Tag = "20160518";
//    TString DataSetsFile[nFiles] = {"GammaConvV1_70.root", "AnalysisResults.root"};

//    TString prefix = "/alice/sim/2015/LHC15h1a1/";
//    TString suffix = "/PWGGA/GA_pp_MC/2093_20160520-0947/";
//    TString DataSet = "LHC15h1a1";
//    const Int_t nRuns = 1;
//    Int_t run[nRuns] = {176749};
//    Int_t first[nRuns] = {3};
//    Int_t last[nRuns] = {178};

//    const Int_t nFiles = 4;
//    TString Tag = "20170407";
//    TString DataSetsFile[nFiles] = {
//      "GammaCalo_119.root","GammaCalo_120.root",
//      "GammaConvCalo_130.root","GammaConvCalo_131.root"
//    };

//    TString prefix = "/alice/sim/2015/LHC15h2d/";
//    TString suffix = "/PWGGA/GA_pp_MC/2855_20170329-2319/";
//    TString DataSet = "LHC15h2d";
//    const Int_t nRuns = 1;
//    TString run[nRuns] = {"184144"};
//    Int_t first[nRuns] = {1};
//    Int_t last[nRuns] = {20};

    const Int_t nFiles = 4;
    TString Tag = "20170407";
    TString DataSetsFile[nFiles] = {
      "GammaCalo_119.root","GammaCalo_120.root",
      "GammaConvCalo_130.root","GammaConvCalo_131.root"
    };

    TString prefix = "/alice/sim/2015/LHC15h1h/";
    TString suffix = "/PWGGA/GA_pp_MC/2850_20170329-2320/";
    TString DataSet = "LHC15h1h";
    const Int_t nRuns = 2;
    TString run[nRuns] = {"192177","192349"};
    Int_t first[nRuns] = {41,102};
    Int_t last[nRuns] = {49,121};


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

//    const Int_t nFiles = 2;
//    TString Tag = "20160206";
//    TString DataSetsFile[nFiles] = {"GammaCalo_201.root","GammaConvCalo_201.root"};

//    TString prefix = "/alice/sim/2014/LHC14j4e/";
//    TString suffix = "/PWGGA/GA_pp_MC/1616_20160207-1047/";
//    TString DataSet = "LHC14j4e";
//    const Int_t nRuns = 2;
//    Int_t run[nRuns] = {128855,128498};
//    Int_t first[nRuns] = {1,1};
//    Int_t last[nRuns] = {50,43};

	for(Int_t iRun=0; iRun<nRuns; iRun++)
	{
        TString PathDataSet = Form("%s%s%s",prefix.Data(),run[iRun].Data(),suffix.Data());

		TString fPathGrid;
		TString fPathLocal;
		for(Int_t k=0; k<nFiles; k++)
		{
            TString merge = Form("hadd %s/DataQA/%s/%s/%s/%s", folder.Data(), Tag.Data(), DataSet.Data(), run[iRun].Data(), DataSetsFile[k].Data());

			for(Int_t j=first[iRun]; j<=last[iRun]; j++)
			{
					fPathGrid = Form("%s%04i/%s", PathDataSet.Data(), j, (DataSetsFile[k]).Data());
                    fPathLocal = Form("%s/DataQA/%s/%s/%s", folder.Data(), Tag.Data(), DataSet.Data(), run[iRun].Data());
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
              fPathLocal = Form("%s/DataQA/%s/%s/%s", folder.Data(), Tag.Data(), DataSet.Data(), run[iRun].Data());
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
