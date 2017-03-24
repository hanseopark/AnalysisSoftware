/*******************************************************************************
 ******  provided by Gamma Conversion Group, PWGGA,                        *****
 ******     Daniel Muehlheim, d.muehlheim@cern.ch                          ***** 
 *******************************************************************************/

#include <iostream>
#include <vector>
#include <fstream>

#include <TFile.h>
#include <TGrid.h>
#include <TSystem.h>
#include <TString.h>
#include <TObject.h>

Bool_t readin(TString fileTxt, std::vector<TString> &vec);
Bool_t copyAlien2Local(TString loc, TString rem);

void Grid_CopyFiles_Runwise(TString folder = "/home/daniel/data/work/pcgGit/AnalysisSoftware", Bool_t isJetJet=kFALSE)
{
    cout<<"Connecting to Alien..."<<endl;
    TGrid::Connect("alien://");
    cout<<"==============================="<<endl;
    cout<<"Successfully connected to Alien"<<endl;
    cout<<"==============================="<<endl;

    const Int_t nFiles = 2;
    TString Tag = "20170324";
    TString DataSetsFile[nFiles] = {"GammaCalo_201.root","GammaConvCalo_201.root"};

    const Int_t nSets = 10;
    TString DataSets[nSets]={"LHC10b_pass4",
                             "LHC10c_pass4",
                             "LHC10d_pass4",
                             "LHC10e_pass4",
                             "LHC10f_pass4",
                             "LHC14j4b",
                             "LHC14j4c",
                             "LHC14j4d",
                             "LHC14j4e",
                             "LHC14j4f"};
    TString DataSetsFolder[nSets]={"LHC10",
                                   "LHC10",
                                   "LHC10",
                                   "LHC10",
                                   "LHC10",
                                   "LHC14j4",
                                   "LHC14j4",
                                   "LHC14j4",
                                   "LHC14j4",
                                   "LHC14j4"};
    TString PrefixDataSets[nSets]={"/alice/data/2010/LHC10b/000",
                                   "/alice/data/2010/LHC10c/000",
                                   "/alice/data/2010/LHC10d/000",
                                   "/alice/data/2010/LHC10e/000",
                                   "/alice/data/2010/LHC10f/000",
                                   "/alice/sim/2014/LHC14j4b/",
                                   "/alice/sim/2014/LHC14j4c/",
                                   "/alice/sim/2014/LHC14j4d/",
                                   "/alice/sim/2014/LHC14j4e/",
                                   "/alice/sim/2014/LHC14j4f/"};
    TString SuffixDataSets[nSets]={"/pass4/PWGGA/GA_pp/2040_20170321-2035/",
                                   "/pass4/PWGGA/GA_pp/2041_20170321-2035/",
                                   "/pass4/PWGGA/GA_pp/2042_20170321-2036/",
                                   "/pass4/PWGGA/GA_pp/2043_20170321-2036/",
                                   "/pass4/PWGGA/GA_pp/2044_20170321-2036/",
                                   "/PWGGA/GA_pp_MC/2816_20170321-2244/",
                                   "/PWGGA/GA_pp_MC/2817_20170321-2224/",
                                   "/PWGGA/GA_pp_MC/2818_20170321-2244/",
                                   "/PWGGA/GA_pp_MC/2819_20170321-2244/",
                                   "/PWGGA/GA_pp_MC/2820_20170321-2246/"
                                  };

	//pp LHC10 7 TeV PhotonQA
//	const Int_t nFiles = 1;
//	TString Tag = "20150922";
//	TString DataSetsFile[nFiles] = {"AnalysisResults.root"};

//	const Int_t nSets = 10;
//	TString DataSets[nSets]={"LHC10b_pass4",
//							 "LHC10c_pass4",
//							 "LHC10d_pass4",
//							 "LHC10e_pass4",
//							 "LHC10f_pass4",
//							 "LHC14j4b",
//							 "LHC14j4c",
//							 "LHC14j4d",
//							 "LHC14j4e",
//							 "LHC14j4f"};
//	TString PrefixDataSets[nSets]={"/alice/data/2010/LHC10b/000",
//								   "/alice/data/2010/LHC10c/000",
//								   "/alice/data/2010/LHC10d/000",
//								   "/alice/data/2010/LHC10e/000",
//								   "/alice/data/2010/LHC10f/000",
//								   "/alice/sim/2014/LHC14j4b/",
//								   "/alice/sim/2014/LHC14j4c/",
//								   "/alice/sim/2014/LHC14j4d/",
//								   "/alice/sim/2014/LHC14j4e/",
//								   "/alice/sim/2014/LHC14j4f/"};
//	TString SuffixDataSets[nSets]={"/pass4/PWGGA/GA_pp/1002_20150922-1900/",
//								   "/pass4/PWGGA/GA_pp/1003_20150922-1825/",
//								   "/pass4/PWGGA/GA_pp/1005_20150922-1901/",
//								   "/pass4/PWGGA/GA_pp/1007_20150922-1906/",
//								   "/pass4/PWGGA/GA_pp/1006_20150922-1902/",
//								   "/PWGGA/GA_pp_MC/1096_20150925-0935/",
//								   "/PWGGA/GA_pp_MC/1087_20150922-1859/",
//								   "/PWGGA/GA_pp_MC/1089_20150922-1856/",
//								   "/PWGGA/GA_pp_MC/1090_20150922-1857/",
//								   "/PWGGA/GA_pp_MC/1095_20150925-0935/"
//								  };



//pp LHC12
//        const Int_t nFiles = 2;
//        TString Tag = "20160518";
//        //TString DataSetsFile[nFiles] = {"GammaConvCalo_101.root","GammaCalo_101.root"};
//        TString DataSetsFile[nFiles] = {"AnalysisResults.root","GammaConvV1_70.root"};

//        const Int_t nSets = 21;
//        TString DataSets[nSets]={
//          "LHC12a", "LHC12b", "LHC12c", "LHC12d", "LHC12f", "LHC12h", "LHC12i",
//          "LHC15h1a1", "LHC15h1b", "LHC15h1c", "LHC15h1d", "LHC15h1f", "LHC15h1h", "LHC15h1i",
//          "LHC15h2a", "LHC15h2b", "LHC15h2c", "LHC15h2d", "LHC15h2f", "LHC15h2h", "LHC15h2i"
//        };
//        TString DataSetsFolder[nSets]={
//          "LHC12", "LHC12", "LHC12", "LHC12", "LHC12", "LHC12", "LHC12",
//          "LHC15h1", "LHC15h1", "LHC15h1", "LHC15h1", "LHC15h1", "LHC15h1", "LHC15h1",
//          "LHC15h2", "LHC15h2", "LHC15h2", "LHC15h2", "LHC15h2", "LHC15h2", "LHC15h2"
//        };
//        TString PrefixDataSets[nSets]={
//              "/alice/data/2012/LHC12a/000",
//              "/alice/data/2012/LHC12b/000",
//              "/alice/data/2012/LHC12c/000",
//              "/alice/data/2012/LHC12d/000",
//              "/alice/data/2012/LHC12f/000",
//              "/alice/data/2012/LHC12h/000",
//              "/alice/data/2012/LHC12i/000",
//              "/alice/sim/2015/LHC15h1a1/",
//              "/alice/sim/2015/LHC15h1b/",
//              "/alice/sim/2015/LHC15h1c/",
//              "/alice/sim/2015/LHC15h1d/",
//              "/alice/sim/2015/LHC15h1f/",
//              "/alice/sim/2015/LHC15h1h/",
//              "/alice/sim/2015/LHC15h1i/",
//              "/alice/sim/2015/LHC15h2a/",
//              "/alice/sim/2015/LHC15h2b/",
//              "/alice/sim/2015/LHC15h2c/",
//              "/alice/sim/2015/LHC15h2d/",
//              "/alice/sim/2015/LHC15h2f/",
//              "/alice/sim/2015/LHC15h2h/",
//              "/alice/sim/2015/LHC15h2i/"};

////        TString SuffixDataSets[nSets]={
////          "/pass2/PWGGA/GA_pp/1578_20160519-0111/",
////          "/pass2/PWGGA/GA_pp/1579_20160519-0111/",
////          "/pass2/PWGGA/GA_pp/1580_20160519-0110/",
////          "/pass2/PWGGA/GA_pp/1581_20160519-0110/",
////          "/pass2/PWGGA/GA_pp/1582_20160519-0110/",
////          "/pass2/PWGGA/GA_pp/1583_20160519-0110/",
////          "/pass2/PWGGA/GA_pp/1584_20160519-0110/",
////                "/PWGGA/GA_pp_MC/2079_20160519-0918/",
////                "/PWGGA/GA_pp_MC/2080_20160519-0918/",
////                "/PWGGA/GA_pp_MC/2081_20160519-0918/",
////                "/PWGGA/GA_pp_MC/2082_20160519-0918/",
////                "/PWGGA/GA_pp_MC/2083_20160519-0918/",
////                "/PWGGA/GA_pp_MC/2084_20160519-0918/",
////                "/PWGGA/GA_pp_MC/2085_20160519-0918/",
////          "/PWGGA/GA_pp_MC/2086_20160519-0917/",
////          "/PWGGA/GA_pp_MC/2087_20160519-0917/",
////          "/PWGGA/GA_pp_MC/2088_20160519-0917/",
////          "/PWGGA/GA_pp_MC/2089_20160519-0917/",
////          "/PWGGA/GA_pp_MC/2090_20160519-0917/",
////          "/PWGGA/GA_pp_MC/2091_20160519-0917/",
////          "/PWGGA/GA_pp_MC/2092_20160519-0917/"
////          };
//        TString SuffixDataSets[nSets]={
//          "/pass2/PWGGA/GA_pp/1586_20160519-2141/",
//          "/pass2/PWGGA/GA_pp/1587_20160519-2144/",
//          "/pass2/PWGGA/GA_pp/1588_20160519-2144/",
//          "/pass2/PWGGA/GA_pp/1589_20160519-2143/",
//          "/pass2/PWGGA/GA_pp/1590_20160519-2143/",
//          "/pass2/PWGGA/GA_pp/1591_20160519-2143/",
//          "/pass2/PWGGA/GA_pp/1592_20160519-2143/",
//                "/PWGGA/GA_pp_MC/2093_20160520-0947/",
//                "/PWGGA/GA_pp_MC/2094_20160519-2336/",
//                "/PWGGA/GA_pp_MC/2095_20160519-2336/",
//                "/PWGGA/GA_pp_MC/2096_20160519-2336/",
//                "/PWGGA/GA_pp_MC/2097_20160519-2335/",
//                "/PWGGA/GA_pp_MC/2098_20160519-2335/",
//                "/PWGGA/GA_pp_MC/2099_20160519-2335/",
//          "/PWGGA/GA_pp_MC/2100_20160519-2335/",
//          "/PWGGA/GA_pp_MC/2101_20160519-2335/",
//          "/PWGGA/GA_pp_MC/2102_20160519-2335/",
//          "/PWGGA/GA_pp_MC/2103_20160519-2335/",
//          "/PWGGA/GA_pp_MC/2104_20160519-2335/",
//          "/PWGGA/GA_pp_MC/2105_20160519-2334/",
//          "/PWGGA/GA_pp_MC/2106_20160519-2334/"
//          };

//                    const Int_t nFiles = 2;
//                    TString Tag = "20160518";
//                    TString DataSetsFile[nFiles] = {"GammaConvCalo_101.root","GammaCalo_101.root"};
//                    const Int_t nSets = 6;
//                    TString DataSets[nSets]={
//                                             "LHC12b-kEMC7",
//                                             "LHC12c-kEMC7",
//                                             "LHC12d-kEMC7",
//                                             "LHC12f-kEMC7",
//                                             "LHC12h-kEMC7",
//                                             "LHC12i-kEMC7"
//                    };
//                TString DataSetsFolder[nSets]={
//                            "LHC12",
//                            "LHC12",
//                            "LHC12",
//                            "LHC12",
//                            "LHC12",
//                            "LHC12"
//                };
//                    TString PrefixDataSets[nSets]={
////                                                   "/alice/data/2012/LHC12a/000",
//                                                   "/alice/data/2012/LHC12b/000",
//                                                   "/alice/data/2012/LHC12c/000",
//                                                   "/alice/data/2012/LHC12d/000",
//                                                   "/alice/data/2012/LHC12f/000",
//                                                   "/alice/data/2012/LHC12h/000",
//                                                   "/alice/data/2012/LHC12i/000"};
//                    TString SuffixDataSets[nSets]={
////                                                 "/pass2/PWGGA/GA_pp/1578_20160519-0111/",
//                                                 "/pass2/PWGGA/GA_pp/1579_20160519-0111/",
//                                                 "/pass2/PWGGA/GA_pp/1580_20160519-0110/",
//                                                 "/pass2/PWGGA/GA_pp/1581_20160519-0110/",
//                                                 "/pass2/PWGGA/GA_pp/1582_20160519-0110/",
//                                                 "/pass2/PWGGA/GA_pp/1583_20160519-0110/",
//                                                 "/pass2/PWGGA/GA_pp/1584_20160519-0110/",
//                                                  };

//            const Int_t nFiles = 1;
//            TString Tag = "20160523";
//            TString DataSetsFile[nFiles] = {"GammaConvCalo_101.root"};

//            const Int_t nSets = 4;
//            TString DataSets[nSets]={
//              "LHC15h1c", "LHC15h1h", "LHC15h2b", "LHC15h2d"
//            };
//            TString DataSetsFolder[nSets]={
//              "LHC15h1", "LHC15h1",
//              "LHC15h2", "LHC15h2"
//            };
//            TString PrefixDataSets[nSets]={
//                  "/alice/sim/2015/LHC15h1c/",
//                  "/alice/sim/2015/LHC15h1h/",
//                  "/alice/sim/2015/LHC15h2b/",
//                  "/alice/sim/2015/LHC15h2d/"};

//            TString SuffixDataSets[nSets]={
//              "/AOD178/PWGGA/GA_pp_MC_AOD/107_20160524-1504/",
//              "/AOD178/PWGGA/GA_pp_MC_AOD/110_20160524-1504/",
//              "/AOD178/PWGGA/GA_pp_MC_AOD/113_20160524-1504/",
//              "/AOD178/PWGGA/GA_pp_MC_AOD/115_20160524-1504/"
//              };

	//pp LHC12
//            const Int_t nFiles = 3;
//            TString Tag = "20160203";
//            TString DataSetsFile[nFiles] = {"GammaConvCalo_130.root","GammaConvCalo_131.root","GammaConvCalo_132.root"};

//            const Int_t nSets = 9;
//            TString DataSets[nSets]={"LHC12a-kINT8",
//                                     "LHC12b-kINT8",
//                                     "LHC12c-kINT8",
//                                     "LHC12d-kINT8",
//                                     "LHC12e-kINT8",
//                                     "LHC12f-kINT8",
//                                     "LHC12g-kINT8",
//                                     "LHC12h-kINT8",
//                                     "LHC12i-kINT8"};
//            TString DataSetsFolder[nSets]={"LHC12",
//                                     "LHC12",
//                                     "LHC12",
//                                     "LHC12",
//                                     "LHC12",
//                                     "LHC12",
//                                     "LHC12",
//                                     "LHC12",
//                                     "LHC12"};
//            TString PrefixDataSets[nSets]={"/alice/data/2012/LHC12a/000",
//                                           "/alice/data/2012/LHC12b/000",
//                                           "/alice/data/2012/LHC12c/000",
//                                           "/alice/data/2012/LHC12d/000",
//                                           "/alice/data/2012/LHC12e/000",
//                                           "/alice/data/2012/LHC12f/000",
//                                           "/alice/data/2012/LHC12g/000",
//                                           "/alice/data/2012/LHC12h/000",
//                                           "/alice/data/2012/LHC12i/000"};
//            TString SuffixDataSets[nSets]={"/pass2/PWGGA/GA_pp/1225_20160203-1808/",
//                                           "/pass2/PWGGA/GA_pp/1226_20160203-1808/",
//                                           "/pass2/PWGGA/GA_pp/1227_20160203-1808/",
//                                           "/pass2/PWGGA/GA_pp/1228_20160203-1808/",
//                                           "/pass2/PWGGA/GA_pp/1229_20160203-1808/",
//                                           "/pass2/PWGGA/GA_pp/1230_20160203-1808/",
//                                           "/pass2/PWGGA/GA_pp/1231_20160203-1808/",
//                                           "/pass2/PWGGA/GA_pp/1232_20160203-1808/",
//                                           "/pass2/PWGGA/GA_pp/1233_20160203-1807/"
//                                          };


//    const Int_t nFiles = 3;
//    TString Tag = "20160401";
//    TString DataSetsFile[nFiles] = {"GammaConvCalo_102.root","GammaConvCalo_105.root","GammaConvCalo_106.root"};

//    const Int_t nSets = 1;
//    TString DataSets[nSets]={"LHC15h2f"};
//    TString DataSetsFolder[nSets]={"LHC15h2"};
//    TString PrefixDataSets[nSets]={"/alice/sim/2015/LHC15h2f/"};
//    TString SuffixDataSets[nSets]={"/PWGGA/GA_pp_MC/1945_20160321-2202/"};


//    const Int_t nFiles = 1;
//    TString Tag = "20160125";
//    TString DataSetsFile[nFiles] = {"GammaConvCalo_120.root"};

//    const Int_t nSets = 16;
//    TString DataSets[nSets]={
//      "LHC15h1a1", "LHC15h1b", "LHC15h1c", "LHC15h1d", "LHC15h1f", "LHC15h1g", "LHC15h1h", "LHC15h1i",
//      "LHC15h2a", "LHC15h2b", "LHC15h2c", "LHC15h2d", "LHC15h2f", "LHC15h2g", "LHC15h2h", "LHC15h2i"
//    };
//    TString DataSetsFolder[nSets]={
//      "LHC15h1", "LHC15h1", "LHC15h1", "LHC15h1", "LHC15h1", "LHC15h1", "LHC15h1", "LHC15h1",
//      "LHC15h2", "LHC15h2", "LHC15h2", "LHC15h2", "LHC15h2", "LHC15h2", "LHC15h2", "LHC15h2"
//    };
//    TString PrefixDataSets[nSets]={
//      "/alice/sim/2015/LHC15h1a1/",
//      "/alice/sim/2015/LHC15h1b/",
//      "/alice/sim/2015/LHC15h1c/",
//      "/alice/sim/2015/LHC15h1d/",
//      "/alice/sim/2015/LHC15h1f/",
//      "/alice/sim/2015/LHC15h1g/",
//      "/alice/sim/2015/LHC15h1h/",
//      "/alice/sim/2015/LHC15h1i/",
//      "/alice/sim/2015/LHC15h2a/",
//      "/alice/sim/2015/LHC15h2b/",
//      "/alice/sim/2015/LHC15h2c/",
//      "/alice/sim/2015/LHC15h2d/",
//      "/alice/sim/2015/LHC15h2f/",
//      "/alice/sim/2015/LHC15h2g/",
//      "/alice/sim/2015/LHC15h2h/",
//      "/alice/sim/2015/LHC15h2i/"};

//    TString SuffixDataSets[nSets]={
//      "/PWGGA/GA_pp_MC/1536_20160126-2114/",
//      "/PWGGA/GA_pp_MC/1537_20160126-2115/",
//      "/PWGGA/GA_pp_MC/1538_20160126-2116/",
//      "/PWGGA/GA_pp_MC/1539_20160126-2116/",
//      "/PWGGA/GA_pp_MC/1540_20160126-2117/",
//      "/PWGGA/GA_pp_MC/1541_20160126-2117/",
//      "/PWGGA/GA_pp_MC/1542_20160126-2118/",
//      "/PWGGA/GA_pp_MC/1543_20160126-2118/",
//      "/PWGGA/GA_pp_MC/1544_20160126-2119/",
//      "/PWGGA/GA_pp_MC/1545_20160126-2119/",
//      "/PWGGA/GA_pp_MC/1546_20160126-2120/",
//      "/PWGGA/GA_pp_MC/1547_20160126-2120/",
//      "/PWGGA/GA_pp_MC/1548_20160126-2121/",
//      "/PWGGA/GA_pp_MC/1549_20160126-2121/",
//      "/PWGGA/GA_pp_MC/1550_20160126-2122/",
//      "/PWGGA/GA_pp_MC/1551_20160126-2122/"
//    };

	//pp LHC12 MC
//			const Int_t nFiles = 2;
//			TString Tag = "20150823";
//			TString DataSetsFile[nFiles] = {"GammaConvCalo_110.root","GammaCalo_110.root"};

//			const Int_t nSets = 2;
//			TString DataSets[nSets]={"LHC14e2a",
//									 "LHC14e2c"};
//			TString PrefixDataSets[nSets]={"/alice/sim/2014/LHC14e2a/",
//										   "/alice/sim/2014/LHC14e2c/"};
//			TString SuffixDataSets[nSets]={"/PWGGA/GA_pp_MC/916_20150824-1151/",
//										   "/PWGGA/GA_pp_MC/918_20150824-1152/"
//										  };

			// pp LHC12 pass2
//				const Int_t nFiles = 2;
//				TString Tag = "20150731";
//				TString DataSetsFile[nFiles] = {"GammaConvCalo_101.root","GammaCalo_101.root"};

//				const Int_t nSets = 1;
//				TString DataSets[nSets]={
//										 "LHC14e2a"};
//				TString PrefixDataSets[nSets]={
//											   "/alice/sim/2014/LHC14e2a/"};
//				TString SuffixDataSets[nSets]={
//											   "/PWGGA/GA_pp_MC/822_20150801-0541/"
//											  };


    //pp LHC11a
//            const Int_t nFiles = 1;
//            TString Tag = "20160226";
//            TString DataSetsFile[nFiles] = {"GammaCalo_1.root"};

//            const Int_t nSets = 3;
//            TString DataSets[nSets]={"LHC11a_p4_wSDD",
//                                     "LHC12f1a",
//                                     "LHC12f1b"
//                                     };
//            TString DataSetsFolder[nSets]={"",
//                                     "",
//                                     ""
//                                    };
//            TString PrefixDataSets[nSets]={"/alice/data/2011/LHC11a/000",
//                                           "/alice/sim/2012/LHC12f1a/",
//                                           "/alice/sim/2012/LHC12f1b/"
//                                          };
//            TString SuffixDataSets[nSets]={"/ESDs/pass4_with_SDD/PWGGA/GA_pp/1293_20160226-1805/",
//                                           "/PWGGA/GA_pp_MC/1716_20160226-1750/",
//                                           "/PWGGA/GA_pp_MC/1717_20160226-1751/"
//                                          };

//	//pp LHC12h
//    const Int_t nFiles = 1;
//    TString Tag = "20161129";
//    TString DataSetsFile[nFiles] = {"AnalysisResults.root"};

//    const Int_t nSets = 1;
//    TString DataSets[nSets]={"LHC12h"
//                            };
//    TString DataSetsFolder[nSets]={""

//                            };
//    TString PrefixDataSets[nSets]={"/alice/data/2012/LHC12h/000"
//                                  };
//    TString SuffixDataSets[nSets]={"/pass2/PWGJE/Jets_EMC_pp/822_20161123-1006/"
//                                  };

    //pPb
//   const Int_t nFiles           = 1;
//   TString Tag                  = "20160307";
//   TString DataSetsFile[nFiles] = {"GammaConvCalo_20.root"};

//   const Int_t nSets            = 10;
//   TString DataSets[nSets]      ={  "LHC13b",
//                                    "LHC13c",
//                                    "LHC13d",
//                                    "LHC13e",
//                                    "LHC13f",
//                                    "LHC13b2_efix1",
//                                    "LHC13b2_efix2",
//                                    "LHC13b2_efix3",
//                                    "LHC13b2_efix4",
//                                    "LHC13e7"};
//   TString DataSetsFolder[nSets]={  "LHC13b",
//                                    "LHC13c",
//                                    "LHC13d",
//                                    "LHC13e",
//                                    "LHC13f",
//                                    "LHC13b2_efix1",
//                                    "LHC13b2_efix2",
//                                    "LHC13b2_efix3",
//                                    "LHC13b2_efix4",
//                                    "LHC13e7"};

//   TString PrefixDataSets[nSets]={  "/alice/data/2013/LHC13b/000",
//                                    "/alice/data/2013/LHC13c/000",
//                                    "/alice/data/2013/LHC13d/000",
//                                    "/alice/data/2013/LHC13e/000",
//                                    "/alice/data/2013/LHC13f/000",
//                                    "/alice/sim/2013/LHC13b2_efix_p1/",
//                                    "/alice/sim/2013/LHC13b2_efix_p2/",
//                                    "/alice/sim/2013/LHC13b2_efix_p3/",
//                                    "/alice/sim/2013/LHC13b2_efix_p4/",
//                                    "/alice/sim/2013/LHC13e7/"};
//   TString SuffixDataSets[nSets]={  "/ESDs/pass3/PWGGA/GA_pPb/536_20160307-2213/",
//                                    "/ESDs/pass2/PWGGA/GA_pPb/537_20160307-2214/",
//                                    "/pass2/PWGGA/GA_pPb/538_20160307-2214/",
//                                    "/pass2/PWGGA/GA_pPb/539_20160307-2215/",
//                                    "/pass2/PWGGA/GA_pPb/540_20160307-2215/",
//                                    "/PWGGA/GA_pPb_MC/721_20160307-2216/",
//                                    "/PWGGA/GA_pPb_MC/722_20160307-2216/",
//                                    "/PWGGA/GA_pPb_MC/723_20160307-2217/",
//                                    "/PWGGA/GA_pPb_MC/724_20160307-2218/",
//                                    "/PWGGA/GA_pPb_MC/725_20160307-2218/"};

//    TString SuffixDataSets[nSets]={"/ESDs/pass3/PWGGA/GA_pPb/284_20141201-1208/",
//                                   "/ESDs/pass2/PWGGA/GA_pPb/285_20141201-1208/",
//                                   "/PWGGA/GA_pPb_MC/394_20141201-1401/",
//                                   "/PWGGA/GA_pPb_MC/395_20141201-1402/",
//                                   "/PWGGA/GA_pPb_MC/396_20141201-1402/",
//                                   "/PWGGA/GA_pPb_MC/397_20141201-1421/"};

    //PbPb 10h
//    const Int_t nFiles = 7;
//    TString Tag = "20161027";
//    TString DataSetsFile[nFiles] = {"GammaCalo_7.root","GammaCalo_9.root","GammaCalo_11.root",
//                                    "GammaConvCalo_7.root","GammaConvCalo_9.root","GammaConvCalo_11.root","GammaConvCalo_14.root"};

//    const Int_t nSets = 4;
//    TString DataSets[nSets]={
//      "LHC10h",
//      "LHC13d2",
//      "LHC13d2",
//      "LHC13d2b"
//    };
//    TString DataSetsFolder[nSets]={
//      "LHC10h",
//      "LHC13d2",
//      "LHC13d2",
//      "LHC13d2b"
//    };
//    TString PrefixDataSets[nSets]={
//      "/alice/data/2010/LHC10h/000",
//      "/alice/sim/2013/LHC13d2/",
//      "/alice/sim/2013/LHC13d2/",
//      "/alice/sim/2013/LHC13d2b/"
//    };
//    TString SuffixDataSets[nSets]={
//      "/ESDs/pass2/PWGGA/GA_PbPb/252_20161006-1815/",
//      "/PWGGA/GA_PbPb_MC/312_20161007-1846/",
//      "/PWGGA/GA_PbPb_MC/313_20161007-1846/",
//      "/PWGGA/GA_PbPb_MC/314_20161007-1239/"
//    };

    //PbPb 11h
//    const Int_t nFiles = 10;
//    TString Tag = "20161027";
//    TString DataSetsFile[nFiles] = {"GammaCalo_7.root","GammaCalo_9.root","GammaCalo_11.root","GammaCalo_35.root",
//                                    "GammaConvCalo_7.root","GammaConvCalo_9.root","GammaConvCalo_11.root","GammaConvCalo_14.root","GammaConvCalo_33.root","GammaConvCalo_35.root"};

//    const Int_t nSets = 6;
//    TString DataSets[nSets]={
//      "LHC11h_EMCAL",
//      "LHC11h_EMCAL",
//      "LHC14a1a_EMCAL",
//      "LHC14a1b_EMCAL",
//      "LHC14a1c_EMCAL",
//      "LHC14a1b_EMCAL"
//    };
//    TString DataSetsFolder[nSets]={
//      "LHC11h",
//      "LHC11h",
//      "LHC14a1a",
//      "LHC14a1b",
//      "LHC14a1c",
//      "LHC14a1b"
//    };
//    TString PrefixDataSets[nSets]={
//      "/alice/data/2011/LHC11h_2/000",
//      "/alice/data/2011/LHC11h_2/000",
//      "/alice/sim/2014/LHC14a1a/",
//      "/alice/sim/2014/LHC14a1b/"
//      "/alice/sim/2014/LHC14a1c/"
//      "/alice/sim/2014/LHC14a1b/"
//    };
//    TString SuffixDataSets[nSets]={
//      "/ESDs/pass2/PWGGA/GA_PbPb/250_20161006-1813/",
//      "/ESDs/pass2/PWGGA/GA_PbPb/251_20161006-1816/",
//      "/PWGGA/GA_PbPb_MC/315_20161007-1250/"
//      "/PWGGA/GA_PbPb_MC/316_20161007-1303/"
//      "/PWGGA/GA_PbPb_MC/317_20161007-1238/"
//      "/PWGGA/GA_PbPb_MC/318_20161007-1311/"
//    };

    std::vector<TString> vecRuns;
    std::vector<TString> vecBins;
    std::vector<TString> vecErrors[nSets];
    TString fDataSet;
    TString fPathGrid;
    TString fPathLocal;
    TString fileTxt;

    Int_t nErr[nSets];

    for(Int_t i=0; i<nSets; i++)
    {
        nErr[i]=0;
        vecErrors[i].clear();
        vecRuns.clear();
        vecBins.clear();
		fDataSet = DataSets[i];
        if(fDataSet.CompareTo("")==0) continue;
        fileTxt = Form("%s/DownloadAndDataPrep/runlists/runNumbers%s.txt", folder.Data(), fDataSet.Data());
        cout << "\n------------------------------------------------------" << endl;
        if(!readin(fileTxt, vecRuns)) cout << "\n\n\n**********************ERROR, no Run Numbers could be found!**********************\n\n\n" << endl;
        cout << "------------------------------------------------------" << endl;
        if(isJetJet){
          fileTxt = Form("%s/DownloadAndDataPrep/binsJetJet%s.txt", folder.Data(), fDataSet.Data());
          cout << "\n------------------------------------------------------" << endl;
          if(!readin(fileTxt, vecBins)) cout << "\n\n\n**********************ERROR, no Jet Jet Bins could be found!**********************\n\n\n" << endl;
          cout << "------------------------------------------------------" << endl;
        }

		Bool_t doNormalFolder = kFALSE;
		if(DataSetsFolder[i].IsNull()) doNormalFolder = kTRUE;

        for(Int_t j=0; j<(Int_t)vecRuns.size(); j++)
        {
			// exclude files from unavailable SEs
            //if(fDataSet.CompareTo("LHC12g")==0 && vecRuns.at(j).CompareTo("188443")==0) continue;

            Int_t nBins = 1;
            if(isJetJet) nBins = (Int_t)vecBins.size();
            for(Int_t b=0; b<nBins; b++)
            {
              for(Int_t k=0; k<nFiles; k++)
              {
                  if(isJetJet){
                    fPathGrid = Form("%s%s/%s%s%s", PrefixDataSets[i].Data(), vecBins.at(b).Data(), vecRuns.at(j).Data(), SuffixDataSets[i].Data(), DataSetsFile[k].Data());
                    if(doNormalFolder) fPathLocal = Form("%s/DataQA/%s/%s/%s/%s", folder.Data(), Tag.Data(), fDataSet.Data(), vecBins.at(b).Data(), vecRuns.at(j).Data());
                    else fPathLocal = Form("%s/DataQA/%s/%s/%s/%s", folder.Data(), Tag.Data(), DataSetsFolder[i].Data(), vecBins.at(b).Data(), vecRuns.at(j).Data());
                    gSystem->Exec(Form("mkdir -p %s",fPathLocal.Data()));
                  }else{
                    fPathGrid = Form("%s%s%s%s", PrefixDataSets[i].Data(), vecRuns.at(j).Data(), SuffixDataSets[i].Data(), DataSetsFile[k].Data());
                    if(doNormalFolder) fPathLocal = Form("%s/DataQA/%s/%s/%s", folder.Data(), Tag.Data(), fDataSet.Data(), vecRuns.at(j).Data());
                    else fPathLocal = Form("%s/DataQA/%s/%s/%s", folder.Data(), Tag.Data(), DataSetsFolder[i].Data(), vecRuns.at(j).Data());
                    gSystem->Exec(Form("mkdir -p %s",fPathLocal.Data()));
                  }

                  fPathLocal+="/"; fPathLocal+=DataSetsFile[k];

                  cout << endl;
                  cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
                  cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
                  cout << "Copying from (grid): " << fPathGrid.Data() << endl;
                  cout << "Copying to (local): " << fPathLocal.Data() << endl;
                  cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
                  cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

                  TFile fileCheck(fPathLocal.Data());
                  if(!fileCheck.IsZombie()) {cout << "\n\t\t>>>>>>>>>>>>>>>>>>Info: ROOT-File |" << fPathLocal.Data() << "| does already exist! Continue...<<<<<<<<<<<<<<\n" << endl; continue;}

                  if(copyAlien2Local(fPathGrid,fPathLocal)) continue;
                  else{
                    cout << "\n\n\t**********************************************************************************************" << endl;
                    cout << "\t**********************************************************************************************" << endl;
                    cout << "\t*******************Err: copyAlien2Local(), check runlist in photonconv rep!*******************" << endl;
                    cout << "\t**********************************************************************************************" << endl;
                    cout << "\t**********************************************************************************************\n" << endl;
                    nErr[i]++;
                    vecErrors[i].push_back(vecRuns.at(j));
                  }
              }
            }
          }

          if(!doNormalFolder) gSystem->Exec(Form("ln -s %s %s",
											   Form("%s/DataQA/%s/%s", folder.Data(), Tag.Data(), DataSetsFolder[i].Data()),
											   Form("%s/DataQA/%s/%s", folder.Data(), Tag.Data(), fDataSet.Data())));
    }

    for(Int_t i=0; i<nSets; i++){
      cout << "DataSet: " << DataSets[i].Data() << ", number of errors: " << nErr[i] << endl;
      cout << "\t\tRuns: ";
      for(Int_t iRuns=0; iRuns < (Int_t)vecErrors[i].size(); iRuns++) cout << vecErrors[i].at(iRuns).Data() << ", ";
      cout << endl;
    }

    return;
}

Bool_t readin(TString fileTxt, std::vector<TString> &vec){
    cout << Form("Reading from %s...", fileTxt.Data()) << endl;
    fstream file;
    TString fVar;
    Int_t totalN=0;
    file.open(fileTxt.Data(), ios::in);
    if(file.good())
    {
        file.seekg(0L, ios::beg);
        if(fileTxt.Contains("binsJetJet")) cout << "Processing Bins: \"";
        else cout << "Processing Runs: \"";
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
    if(fileTxt.Contains("binsJetJet"))  cout << "...done!\n\nIn total " << totalN << " Bins will be processed!" << endl;
    else cout << "...done!\n\nIn total " << totalN << " Runs will be processed!" << endl;
    if(totalN > 0) return kTRUE;
    else return kFALSE;
}

Bool_t copyAlien2Local(TString loc, TString rem){
   TString sl(Form("alien://%s", loc.Data()));
   TString sr(Form("file://%s", rem.Data()));
   Bool_t ret = TFile::Cp(sl,sr);
   if (!ret) cout << Form("ERROR: Failed to copy %s to %s", sl.Data(), sr.Data()) << endl;
   return ret;
}
