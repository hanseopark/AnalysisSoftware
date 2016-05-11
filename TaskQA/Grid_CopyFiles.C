#include <algorithm>
#include <vector>
#include <iterator>
#include <iostream>
#include <fstream>
#include <TObject.h>
#include <TString.h>
#include <TFile.h>
#include <TGrid.h>
#include <TSystem.h>

void ChangeStrucToStd(TString nameInputFile, TString namefileOutput, TString nameInputList);
Bool_t copyAlien2Local(TString loc, TString rem);
Bool_t readin(TString fileRuns, std::vector<TString> &vec){
    cout << Form("Reading from %s...", fileRuns.Data()) << endl;
    std::fstream file;
    TString fVar;
    Int_t totalN=0;
    file.open(fileRuns.Data(), ios::in);
    if(file.good())
    {
        file.seekg(0L, ios::beg);
        while(!file.eof())
        {
            file >> fVar;
            if(fVar.Sizeof()>1)
            {
                vec.push_back(fVar);
                totalN++;
            }
        }
    }
    file.close();
    if(totalN > 0) return kTRUE;
    else return kFALSE;
}

void Grid_CopyFiles(TString system = "pp", TString type = "ESD", TString folder = "/home/daniel/data/work/Grid")
{
    cout<<"Connecting to Alien..."<<endl;
    TGrid::Connect("alien://");
    cout<<"==============================="<<endl;
    cout<<"Successfully connected to Alien"<<endl;
    cout<<"==============================="<<endl;

//    const Int_t nSets = 21;
//    const Int_t nData = 7;
//    TString DataSets[nSets]={
//      "LHC12a", "LHC12b", "LHC12c", "LHC12d", "LHC12f", "LHC12h", "LHC12i",
//      "LHC15h1a1", "LHC15h1b", "LHC15h1c", "LHC15h1d", "LHC15h1f", "LHC15h1h", "LHC15h1i",
//      "LHC15h2a", "LHC15h2b", "LHC15h2c", "LHC15h2d", "LHC15h2f", "LHC15h2h", "LHC15h2i"
//    };

////    TString train = "Legotrain-vAN-20160120-8TeV_NL";
////    Int_t trainRuns[nSets] = {1180,1181,1182,1183,1184,1186,1187,
////                             1487,1488,1489,1490,1491,1493,1494,
////                             1495,1496,1497,1498,1499,1501,1502
////                             };
////    TString train = "Legotrain-vAN-20160124-8TeV_NL_2";
////    Int_t trainRuns[nSets] = {1188,1189,1190,1191,1192,1194,1195,
////                             1503,1504,1505,1506,1507,1509,1510,
////                             1511,1512,1513,1514,1515,1517,1518
////                             };
////    TString train = "Legotrain-vAN-20160206-8TeV_NL_3";
////    Int_t trainRuns[nSets] = {1241,1242,1243,1244,1245,1247,1248,
////                             1596,1597,1598,1599,1600,1602,1603,
////                             1604,1605,1606,1607,1608,1610,1611
////                             };
////    TString train = "Legotrain-vAN-20160206-8TeV_NL_4";
////    Int_t trainRuns[nSets] = {1241,1242,1243,1244,1245,1247,1248,
////                             1623,1624,1625,1626,1627,11629,1630,
////                             1631,1632,1633,1634,1635,1637,1638
////                             };
////    TString train = "Legotrain-vAN-20160206-8TeV_NL_5";
////    Int_t trainRuns[nSets] = {1241,1242,1243,1244,1245,1247,1248,
////                             1662,1663,1664,1665,1666,1668,1669,
////                             1670,1671,1672,1673,1674,1676,1677
////                             };
//    TString train = "Legotrain-vAN-20160215-MultWeight_DistBC";
//    Int_t trainRuns[nSets] = {1280,1281,1282,1283,1284,1286,1287,
//                             1684,1685,1686,1687,1688,1690,1691,
//                             1692,1693,1694,1695,1696,1698,1699
//                             };

//    TString runlist[nSets] = {"merge_runlist_9","merge_runlist_9","merge_runlist_9","merge_runlist_9","merge_runlist_9","merge_runlist_9","merge_runlist_9",
//                            "merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2",
//                            "merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2"
//                          };

////    const Int_t nFiles = 2;
////    TString Files[nFiles] = {"GammaCalo_111","GammaConvCalo_31"};

////    const Int_t nFiles = 4;
////    TString Files[nFiles] = {"GammaCalo_109","GammaCalo_110","GammaConvCalo_110","GammaConvCalo_111"};

//    const Int_t nFiles = 4;
//    TString Files[nFiles] = {"GammaCalo_101","GammaCalo_108","GammaConvCalo_101","GammaConvCalo_112"};

//    const Int_t nMerge = 11;
//    TString strMerge[nMerge]={"LHC12","LHC15h1","LHC15h2","LHC15h","LHC15ha","LHC15hb","LHC15hc","LHC15hd","LHC15hf","LHC15hh","LHC15hi"};
//    std::vector<Int_t> mergeVec[nMerge];
//    std::vector<Int_t>::iterator it;
//    for(Int_t i=0; i<nSets; i++){
//      if(0<=i && i<=6) mergeVec[0].push_back(i);
//      if(7<=i && i<=13) mergeVec[1].push_back(i);
//      if(14<=i && i<=20) mergeVec[2].push_back(i);
//      if(7<=i && i<=20) mergeVec[3].push_back(i);
//      //merge MCs - for example h1a1 and h2a, h1b and h2b...
//      for(Int_t j=0; j<7; j++){
//        if(i==7+j || i==14+j) mergeVec[4+j].push_back(i);
//      }
//    }

//*********************************************************************************************************************************
//*********************************************************************************************************************************
//*********************************************************************************************************************************

//        const Int_t nSets = 21;
//        const Int_t nData = 7;
//        TString DataSets[nSets]={
//          "LHC12a", "LHC12b", "LHC12c", "LHC12d", "LHC12f", "LHC12h", "LHC12i",
//          "LHC15h1a1", "LHC15h1b", "LHC15h1c", "LHC15h1d", "LHC15h1f", "LHC15h1h", "LHC15h1i",
//          "LHC15h2a", "LHC15h2b", "LHC15h2c", "LHC15h2d", "LHC15h2f", "LHC15h2h", "LHC15h2i"
//        };


//        TString train = "Legotrain-vAN-20160320-8TeV_Systematics_Calo";

//        TString runlist[nSets] = {"merge_runlist_9","merge_runlist_9","merge_runlist_9","merge_runlist_9","merge_runlist_9","merge_runlist_9","merge_runlist_9",
//                                "merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2",
//                                "merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2"
//                              };

//        Int_t trainRuns[nSets] = {1394,1395,1396,1397,1398,1399,1400,
//                                 1858,1859,1860,1861,1862,1863,1864,
//                                 1865,1866,1867,1868,1869,1870,1871
//                                 };
//        const Int_t nFiles = 3;
//        TString Files[nFiles] = {"GammaCalo_101","GammaCalo_106","GammaCalo_121"};
//*********************************************************************************************************************************
//        Int_t trainRuns[nSets] = {1403,1404,1405,1406,1407,1408,1409,
//                                 1872,1873,1874,1875,1876,1877,1878,
//                                 1879,1880,1881,1882,1883,1884,1885
//                                 };
//        const Int_t nFiles = 1;
//        TString Files[nFiles] = {"GammaCalo_101"};
//        for(Int_t n=0; n<nSets; n++){DataSets[n]+="_V1";}
//*********************************************************************************************************************************
//        Int_t trainRuns[nSets] = {1410,1411,1412,1413,1414,1415,1416,
//                                 1892,1893,1894,1895,1896,1897,1898,
//                                 1899,1900,1901,1902,1903,1904,1905
//                                 };
//        const Int_t nFiles = 3;
//        TString Files[nFiles] = {"GammaCalo_102","GammaCalo_103","GammaCalo_107"};
//*********************************************************************************************************************************
//        Int_t trainRuns[nSets] = {1417,1418,1419,1420,1421,1422,1423,
//                                 1906,1907,1908,1909,1910,1911,1912,
//                                 1913,1914,1915,1916,1917,1918,1919
//                                 };
//        const Int_t nFiles = 2;
//        TString Files[nFiles] = {"GammaCalo_104","GammaCalo_105"};
//*********************************************************************************************************************************


//        train = "Legotrain-vAN-20160320-8TeV_Systematics_ConvCalo";
        //*********************************************************************************************************************************
//                Int_t trainRuns[nSets] = {1403,1404,1405,1406,1407,1408,1409,
//                                         1872,1873,1874,1875,1876,1877,1878,
//                                         1879,1880,1881,1882,1883,1884,1885
//                                         };
//                const Int_t nFiles = 1;
//                TString Files[nFiles] = {"GammaConvCalo_101"};
//                for(Int_t n=0; n<nSets; n++){DataSets[n]+="_V1";}
        //*********************************************************************************************************************************
//                Int_t trainRuns[nSets] = {1432,1433,1434,1435,1436,1437,1438,
//                                         1920,1921,1922,1923,1924,1925,1926,
//                                         1927,1928,1929,1930,1931,1932,1933
//                                         };
//                const Int_t nFiles = 3;
//                TString Files[nFiles] = {"GammaConvCalo_103","GammaConvCalo_104","GammaConvCalo_109"};
        //*********************************************************************************************************************************
//                Int_t trainRuns[nSets] = {1440,1441,1442,1443,1444,1445,1446,
//                                         1934,1935,1936,1937,1938,1939,1940,
//                                         1941,1942,1943,1944,1945,1946,1947
//                                         };
//                const Int_t nFiles = 3;
//                TString Files[nFiles] = {"GammaConvCalo_102","GammaConvCalo_105","GammaConvCalo_106"};
        //*********************************************************************************************************************************
//                Int_t trainRuns[nSets] = {1447,1448,1449,1450,1451,1452,1453,
//                                         1948,1949,1950,1951,1952,1953,1954,
//                                         1955,1956,1957,1958,1959,1960,1961
//                                         };
//                const Int_t nFiles = 3;
//                TString Files[nFiles] = {"GammaConvCalo_107","GammaConvCalo_108","GammaConvCalo_118"};
        //*********************************************************************************************************************************
//                Int_t trainRuns[nSets] = {1454,1455,1456,1457,1458,1459,1460,
//                                         1962,1963,1964,1965,1966,1967,1968,
//                                         1969,1970,1971,1972,1973,1974,1975
//                                         };
//                const Int_t nFiles = 3;
//                TString Files[nFiles] = {"GammaConvCalo_119","GammaConvCalo_126","GammaConvCalo_127"};
        //*********************************************************************************************************************************
//        const Int_t nMerge = 11;
//        TString strMerge[nMerge]={"LHC12","LHC15h1","LHC15h2","LHC15h","LHC15ha","LHC15hb","LHC15hc","LHC15hd","LHC15hf","LHC15hh","LHC15hi"};
////        for(Int_t n=0; n<nMerge; n++){strMerge[n]+="_V1";}
//        std::vector<Int_t> mergeVec[nMerge];
//        std::vector<Int_t>::iterator it;
//        for(Int_t i=0; i<nSets; i++){
//          if(0<=i && i<=6) mergeVec[0].push_back(i);
//          if(7<=i && i<=13) mergeVec[1].push_back(i);
//          if(14<=i && i<=20) mergeVec[2].push_back(i);
//          if(7<=i && i<=20) mergeVec[3].push_back(i);
//          //merge MCs - for example h1a1 and h2a, h1b and h2b...
//          for(Int_t j=0; j<7; j++){
//            if(i==7+j || i==14+j) mergeVec[4+j].push_back(i);
//          }
//        }

//*********************************************************************************************************************************
//*********************************************************************************************************************************
//*********************************************************************************************************************************

//    const Int_t nSets = 2;
//    const Int_t nData = 1;
//    TString DataSets[nSets]={
//      "LHC12b", "LHC15h1b"
//    };

//    TString train = "Legotrain-vAN-20160309-8TeV-dAODtest";

//    Int_t trainRuns[nSets] = {54,82};
//    TString runlist[nSets] = {"merge_runlist_9","merge_runlist_2"};

//    const Int_t nFiles = 2;
//    TString Files[nFiles] = {"GammaCalo_101","GammaConvCalo_120"};

//    const Int_t nMerge = 0;
//    TString strMerge[nMerge]={};
//    std::vector<Int_t> mergeVec[nMerge];
//    std::vector<Int_t>::iterator it;

//    const Int_t nSets = 14;
//    const Int_t nData = 7;
//    TString DataSets[nSets]={
//      "LHC12d", "LHC12d_A125","LHC12d_A150","LHC12d_A50","LHC12d_A75","LHC12d_S400","LHC12d_S600",
//      "LHC15h1d","LHC15h1d_A125","LHC15h1d_A150","LHC15h1d_A50","LHC15h1d_A75","LHC15h1d_S400","LHC15h1d_S600"
//    };

//    TString train = "Legotrain-vAN-20160401-8TeV-TenderVariations";
//    Int_t trainRuns[nSets] = {1471,1472,1473,1474,1475,1476,1477,1979,1980,1981,1982,1983,1984,1985};

//    TString runlist[nSets] = {"merge_runlist_9","merge_runlist_9","merge_runlist_9","merge_runlist_9","merge_runlist_9","merge_runlist_9","merge_runlist_9"
//                             ,"merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2"
//                             };

//    const Int_t nFiles = 2;
//    TString Files[nFiles] = {"GammaCalo_101","GammaConvCalo_101"};


//    const Int_t nSets = 4;
//    const Int_t nData = 4;
//    TString DataSets[nSets]={
//      "LHC12d_T50", "LHC12d_T100","LHC12d_T200","LHC12d_T300"
//    };

//    TString train = "Legotrain-vAN-20160401-8TeV-TenderVariations";
//    Int_t trainRuns[nSets] = {1483,1484,1485,1486};

//    TString runlist[nSets] = {"merge_runlist_9","merge_runlist_9","merge_runlist_9","merge_runlist_9"
//                             };

//    const Int_t nFiles = 1;
//    TString Files[nFiles] = {"GammaConvCalo_101"};

//    const Int_t nMerge = 0;
//    TString strMerge[nMerge]={};
//    std::vector<Int_t> mergeVec[nMerge];
//    std::vector<Int_t>::iterator it;


//    mergeVec[0].push_back(11); mergeVec[0].push_back(12);
//    mergeVec[1].push_back(13); mergeVec[1].push_back(14);
//    mergeVec[2].push_back(15); mergeVec[2].push_back(16);
//    mergeVec[3].push_back(17); mergeVec[3].push_back(18);
//    mergeVec[4].push_back(19); mergeVec[4].push_back(20);
//    mergeVec[5].push_back(21); mergeVec[5].push_back(22);

    const Int_t nSets = 1;
    const Int_t nData = 0;
    TString DataSets[nSets]={
      "LHC16c2"
    };

    TString train = "Legotrain-vAN-20160413-8TeV-JetJetValid";
    Int_t trainRuns[nSets] = {2009};

    TString runlist[nSets] = {"merge"};

    const Int_t nFiles = 2;
    TString Files[nFiles] = {"GammaCalo_111","GammaConvCalo_31"};

    const Int_t nMerge = 0;
    TString strMerge[nMerge]={};
    std::vector<Int_t> mergeVec[nMerge];
    std::vector<Int_t>::iterator it;
//*********************************************************************************************************************************
//*********************************************************************************************************************************
//*********************************************************************************************************************************

//    const Int_t nSets = 8;
//    const Int_t nData = 2;
//    TString DataSets[nSets]={
//      "LHC11a", "LHC13g", "LHC12f1a", "LHC12f1b", "LHC15g1a", "LHC15g2", "LHC15a3a","LHC15a3a_plus"
//    };

//    TString train = "Legotrain-vAN-20160226-2.76TeV-Calo-QA";
//    Int_t trainRuns[nSets] = {1293,1294,1716,1717,1720,1721,1718,1719};

//    TString runlist[nSets] = {"merge_runlist_3","merge","merge_runlist_1","merge_runlist_1","merge","merge","merge","merge"};

//    const Int_t nFiles = 2;
//    TString Files[nFiles] = {"GammaCalo_1","GammaCalo_60"};

//    const Int_t nMerge = 1;
//    TString strMerge[nMerge]={"LHC15a3a+"};
//    std::vector<Int_t> mergeVec[nMerge];
//    std::vector<Int_t>::iterator it;
//    mergeVec[0].push_back(6); mergeVec[0].push_back(7);

//    const Int_t nSets = 4;
//    const Int_t nData = 1;
//    TString DataSets[nSets]={
//      "LHC13g", "LHC15g2", "LHC15a3a","LHC15a3a_plus"
//    };

//    TString train = "Legotrain-vAN-20160302-2.76TeV-ConvCalo-QA";
//    Int_t trainRuns[nSets] = {1308,1752,1750,1751};

//    TString runlist[nSets] = {"merge","merge","merge","merge"};

//    const Int_t nFiles = 2;
//    TString Files[nFiles] = {"GammaConvCalo_95","GammaConvCalo_96"};

//    const Int_t nMerge = 1;
//    TString strMerge[nMerge]={"LHC15a3a+"};
//    std::vector<Int_t> mergeVec[nMerge];
//    std::vector<Int_t>::iterator it;
//    mergeVec[0].push_back(2); mergeVec[0].push_back(3);
//*********************************************************************************************************************************
//*********************************************************************************************************************************
//*********************************************************************************************************************************

//    const Int_t nSets = 7;
//    const Int_t nData = 7;
//    TString DataSets[nSets]={
//      "LHC12a", "LHC12b", "LHC12c", "LHC12d", "LHC12f", "LHC12h", "LHC12i"
//    };

//    TString train = "Legotrain-vAN-20160301-8TeV-CaloMerged";
//    Int_t trainRuns[nSets] = {1301,1302,1303,1304,1305,1306,1307};

//    TString runlist[nSets] = {"merge_runlist_11","merge_runlist_11","merge_runlist_11","merge_runlist_11","merge_runlist_11","merge_runlist_11","merge_runlist_11"};

//    const Int_t nFiles = 4;
//    TString Files[nFiles] = {"GammaCaloMerged_115","GammaCaloMerged_116","GammaCaloMerged_117","GammaCaloMerged_118"};

//    const Int_t nMerge = 1;
//    TString strMerge[nMerge]={"LHC12"};
//    std::vector<Int_t> mergeVec[nMerge];
//    std::vector<Int_t>::iterator it;
//    for(Int_t i=0; i<nSets; i++){mergeVec[0].push_back(i);}

    //---------------------------------------------------------------------------------------------------

    TString alienFolder;
    TString alienFolder_MC;

    if(type.CompareTo("ESD")==0){
      alienFolder = Form("/alice/cern.ch/user/a/alitrain/PWGGA/GA_%s/",system.Data());
      alienFolder_MC = Form("/alice/cern.ch/user/a/alitrain/PWGGA/GA_%s_MC/",system.Data());
    }else{
      alienFolder = Form("/alice/cern.ch/user/a/alitrain/PWGGA/GA_%s_AOD/",system.Data());
      alienFolder_MC = Form("/alice/cern.ch/user/a/alitrain/PWGGA/GA_%s_MC_AOD/",system.Data());
    }
    gSystem->Exec(Form("alien_ls %s > tempData.log",alienFolder.Data()));
    gSystem->Exec(Form("alien_ls %s > tempMC.log", alienFolder_MC.Data()));

    TString strTrain[nSets];
    std::vector<TString> vecStrTrain;
    std::vector<TString> vecStrTrainMC;
    if(!readin("tempData.log", vecStrTrain)) cout << "\n\n\n**********************ERROR!**********************\n\n\n" << endl;
    if(!readin("tempMC.log", vecStrTrainMC)) cout << "\n\n\n**********************ERROR!**********************\n\n\n" << endl;

    for(Int_t i=0; i<nSets; i++){
      TString temp;
      TString tempRuns;
      if(i<nData){
        for(Int_t j=0; j<(Int_t)vecStrTrain.size(); j++){
          tempRuns = Form("%i",trainRuns[i]);
          temp = vecStrTrain.at(j);
          if(temp.BeginsWith(tempRuns)){
            strTrain[i] = temp;
            break;
          }
        }
      }else{
        for(Int_t j=0; j<(Int_t)vecStrTrainMC.size(); j++){
          tempRuns = Form("%i",trainRuns[i]);
          temp = vecStrTrainMC.at(j);
          if(temp.BeginsWith(tempRuns)){
            strTrain[i] = temp;
            break;
          }
        }
      }
    }

    Int_t nErr[nSets];
    TString fPathGrid;
    TString fPathLocal;

    fPathLocal = Form("%s/%s", folder.Data(), train.Data());
    gSystem->Exec(Form("mkdir -p %s",fPathLocal.Data()));

    TString mergePeriod[nMerge][nFiles];
    for(Int_t iFiles=0; iFiles<nFiles; iFiles++){
      for(Int_t iM=0; iM<nMerge; iM++){
        TString mergeP = Form("%s/%s_%s.root", fPathLocal.Data(), strMerge[iM].Data(), Files[iFiles].Data());
        mergePeriod[iM][iFiles] = Form("hadd -f -k %s",mergeP.Data());
      }
    }

    for(Int_t i=0; i<nSets; i++)
    {
        nErr[i]=0;
        for(Int_t k=0; k<nFiles; k++)
        {
          if(i<nData) fPathGrid = Form("%s%s/%s/%s.root", alienFolder.Data(), strTrain[i].Data(), runlist[i].Data(), Files[k].Data());
          else fPathGrid = Form("%s%s/%s/%s.root", alienFolder_MC.Data(), strTrain[i].Data(), runlist[i].Data(), Files[k].Data());

          TString fPathTemp = Form("%s/%s", fPathLocal.Data(), strTrain[i].Data());
          gSystem->Exec(Form("mkdir -p %s",fPathTemp.Data()));

          fPathTemp+=Form("/%s_",DataSets[i].Data());
          fPathTemp+=Files[k].Data();
          fPathTemp+=".root";

          cout << endl;
          cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
          cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
          cout << "Copying from (grid): " << fPathGrid.Data() << endl;
          cout << "Copying to (local): " << fPathTemp.Data() << endl;
          cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
          cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

          TFile fileCheck(fPathTemp.Data());
          if(!fileCheck.IsZombie()) {
            cout << "\n\t\t>>>>>>>>>>>>>>>>>>Info: ROOT-File |" << fPathTemp.Data() << "| does already exist! Continue...<<<<<<<<<<<<<<\n" << endl;
            for(Int_t iM=0; iM<nMerge; iM++){
              it = find( mergeVec[iM].begin(), mergeVec[iM].end(), i);
              if( it!=mergeVec[iM].end() ) mergePeriod[iM][k] += Form(" %s",fPathTemp.Data());
            }
            gSystem->Exec(Form("ln -s %s %s/%s_%s.root",fPathTemp.Data(), fPathLocal.Data(), DataSets[i].Data(), Files[k].Data()));
            continue;
          }

          if(copyAlien2Local(fPathGrid,fPathTemp)){
            ChangeStrucToStd(fPathTemp.Data(),fPathTemp.Data(),Files[k].Data());
            for(Int_t iM=0; iM<nMerge; iM++){
              it = find( mergeVec[iM].begin(), mergeVec[iM].end(), i);
              if( it!=mergeVec[iM].end() ) mergePeriod[iM][k] += Form(" %s",fPathTemp.Data());
            }
            gSystem->Exec(Form("ln -s %s %s/%s_%s.root",fPathTemp.Data(), fPathLocal.Data(), DataSets[i].Data(), Files[k].Data()));
            continue;
          }
          else{
            cout << "\n\n\t**********************************************************************************************" << endl;
            cout << "\t**********************************************************************************************" << endl;
            cout << "\t*******************Err: copyAlien2Local()!****************************************************" << endl;
            cout << "\t**********************************************************************************************" << endl;
            cout << "\t**********************************************************************************************\n" << endl;
            nErr[i]++;
          }
        }
    }

    cout << "\n------------------------------------------------------" << endl;
    cout << "Merging: " << endl;
    cout << "------------------------------------------------------\n" << endl;

    for(Int_t iFiles=0; iFiles<nFiles; iFiles++){
      for(Int_t iM=0; iM<nMerge; iM++){
        cout << "Merging " << strMerge[iM].Data() << " ..." << endl;
        gSystem->Exec(mergePeriod[iM][iFiles].Data());
        cout << "done!" << endl;
      }
    }

    for(Int_t i=0; i<nSets; i++){
      cout << "DataSet: " << DataSets[i].Data() << ", number of errors: " << nErr[i] << endl;
    }

    gSystem->Exec("rm tempData.log");
    gSystem->Exec("rm tempMC.log");

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

Bool_t copyAlien2Local(TString loc, TString rem){
   TString sl(Form("alien://%s", loc.Data()));
   TString sr(Form("file://%s", rem.Data()));
   Bool_t ret = TFile::Cp(sl,sr);
   if (!ret) cout << Form("ERROR: Failed to copy %s to %s", sl.Data(), sr.Data()) << endl;
   return ret;
}
