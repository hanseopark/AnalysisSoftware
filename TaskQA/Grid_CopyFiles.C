/*******************************************************************************
 ******  provided by Gamma Conversion Group, PWGGA,                        *****
 ******     Daniel Muehlheim, d.muehlheim@cern.ch                          ***** 
 *******************************************************************************/

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

//    TString train = "Legotrain-vAN-20160120-8TeV_NL";
//    Int_t trainRuns[nSets] = {1180,1181,1182,1183,1184,1186,1187,
//                             1487,1488,1489,1490,1491,1493,1494,
//                             1495,1496,1497,1498,1499,1501,1502
//                             };
//    TString train = "Legotrain-vAN-20160124-8TeV_NL_2";
//    Int_t trainRuns[nSets] = {1188,1189,1190,1191,1192,1194,1195,
//                             1503,1504,1505,1506,1507,1509,1510,
//                             1511,1512,1513,1514,1515,1517,1518
//                             };
//    TString train = "Legotrain-vAN-20160206-8TeV_NL_3";
//    Int_t trainRuns[nSets] = {1241,1242,1243,1244,1245,1247,1248,
//                             1596,1597,1598,1599,1600,1602,1603,
//                             1604,1605,1606,1607,1608,1610,1611
//                             };
//    TString train = "Legotrain-vAN-20160206-8TeV_NL_4";
//    Int_t trainRuns[nSets] = {1241,1242,1243,1244,1245,1247,1248,
//                             1623,1624,1625,1626,1627,11629,1630,
//                             1631,1632,1633,1634,1635,1637,1638
//                             };
//    TString train = "Legotrain-vAN-20160206-8TeV_NL_5";
//    Int_t trainRuns[nSets] = {1241,1242,1243,1244,1245,1247,1248,
//                             1662,1663,1664,1665,1666,1668,1669,
//                             1670,1671,1672,1673,1674,1676,1677
//                             };
//    TString train = "Legotrain-vAN-20160215-MultWeight_DistBC";
//    Int_t trainRuns[nSets] = {1280,1281,1282,1283,1284,1286,1287,
//                             1684,1685,1686,1687,1688,1690,1691,
//                             1692,1693,1694,1695,1696,1698,1699
//                             };
//    TString train = "Legotrain-vAN-20160518-8TeV-QA";
//    Int_t trainRuns[nSets] = {1578,1579,1580,1581,1582,1583,1584,
//                             2079,2080,2081,2082,2083,2084,2085,
//                             2086,2087,2088,2089,2090,2091,2092
//                             };
//    TString train = "Legotrain-vAN-20160522-8TeV-Weighting";
//    Int_t trainRuns[nSets] = {1578,1579,1580,1581,1582,1583,1584,
//                             2107,2108,2109,2110,2111,2112,2113,
//                             2114,2115,2116,2117,2118,2119,2120
//                             };
//    TString train = "Legotrain-vAN-20160527-8TeV-validAOD";
//    Int_t trainRuns[nSets] = {66,67,68,69,70,71,72,
//                             105,106,107,108,109,110,111,
//                             112,113,114,115,116,117,118
//                             };

//    TString runlist[nSets] = {"merge_runlist_9","merge_runlist_9","merge_runlist_9","merge_runlist_9","merge_runlist_9","merge_runlist_9","merge_runlist_9",
//                            "merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2",
//                            "merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2"
//                          };

//    const Int_t nFiles = 2;
//    TString Files[nFiles] = {"GammaCalo_111","GammaConvCalo_31"};

//    const Int_t nFiles = 4;
//    TString Files[nFiles] = {"GammaCalo_109","GammaCalo_110","GammaConvCalo_110","GammaConvCalo_111"};

//    const Int_t nFiles = 4;
//    TString Files[nFiles] = {"GammaCalo_101","GammaCalo_108","GammaConvCalo_101","GammaConvCalo_112"};

//    const Int_t nFiles = 4;
//    TString Files[nFiles] = {"GammaCalo_101","GammaConvV1_78","GammaConvCalo_101","GammaCaloMerged_115"};

//    const Int_t nFiles = 2;
//    TString Files[nFiles] = {"GammaCalo_101","GammaConvCalo_101"};

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


//    const Int_t nSets = 11;
//    const Int_t nData = 11;
//    TString DataSets[nSets]={
//      "LHC12b-kEMC7", "LHC12c-kEMC7", "LHC12d-kEMC7", "LHC12f-kEMC7", "LHC12h-kEMC7", "LHC12i-kEMC7",
//      "LHC12c-kEMCEGA", "LHC12d-kEMCEGA", "LHC12f-kEMCEGA", "LHC12h-kEMCEGA", "LHC12i-kEMCEGA"
//    };

//    TString train = "Legotrain-vAN-20160518-8TeV-QA-Trigger";
//    Int_t trainRuns[nSets] = {1579,1580,1581,1582,1583,1584,
//                              1580,1581,1582,1583,1584
//                             };

//    TString runlist[nSets] = {"merge_runlist_10","merge_runlist_10","merge_runlist_10","merge_runlist_10","merge_runlist_10","merge_runlist_10",
//                            "merge_runlist_11","merge_runlist_11","merge_runlist_11","merge_runlist_11","merge_runlist_11"
//                          };

//    const Int_t nFiles = 4;
//    TString Files[nFiles] = {"GammaCalo_101","GammaConvV1_78","GammaConvCalo_101","GammaCaloMerged_115"};

//    const Int_t nMerge = 2;
//    TString strMerge[nMerge]={"LHC12-kEMC7","LHC12-kEMCEGA"};
//    std::vector<Int_t> mergeVec[nMerge];
//    std::vector<Int_t>::iterator it;
//    for(Int_t i=0; i<nSets; i++){
//      if(0<=i && i<=5) mergeVec[0].push_back(i);
//      if(6<=i && i<=10) mergeVec[1].push_back(i);
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

//    const Int_t nSets = 14;
//    const Int_t nData = 0;
//    TString DataSets[nSets]={
//      "LHC15h1a1", "LHC15h1b", "LHC15h1c", "LHC15h1d", "LHC15h1f", "LHC15h1h", "LHC15h1i",
//      "LHC15h2a", "LHC15h2b", "LHC15h2c", "LHC15h2d", "LHC15h2f", "LHC15h2h", "LHC15h2i"
//    };


//    TString train = "Legotrain-vAN-20160527-8TeV-validAOD";

//    TString runlist[nSets] = {
//                              "merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2",
//                              "merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2"
//                             };

//    Int_t trainRuns[nSets] = {
//                              155,156,157,158,159,160,161,
//                              162,163,164,165,166,167,168
//                             };
//    const Int_t nFiles = 3;
//    TString Files[nFiles] = {"GammaCalo_101","GammaConvCalo_101","GammaConvV1_78"};

//    const Int_t nMerge = 3;
//    TString strMerge[nMerge]={"LHC15h1","LHC15h2","LHC15h"};
//    std::vector<Int_t> mergeVec[nMerge];
//    std::vector<Int_t>::iterator it;
//    for(Int_t i=0; i<nSets; i++){
//      if(0<=i && i<=6) mergeVec[0].push_back(i);
//      if(7<=i && i<=13) mergeVec[1].push_back(i);
//      if(0<=i && i<=13) mergeVec[2].push_back(i);
//    }


//        const Int_t nSets = 1;
//        const Int_t nData = 0;
//        TString DataSets[nSets]={
//          "LHC15h1b"
//        };


//        TString train = "Legotrain-vAN-20161003-8TeV-validAOD_matching";

//        TString runlist[nSets] = {
//                                  "merge_runlist_2"
//                                 };

//        Int_t trainRuns[nSets] = {
//                                  212
//                                 };
//        const Int_t nFiles = 2;
//        TString Files[nFiles] = {"GammaCalo_101","GammaConvCalo_101"};

//        const Int_t nMerge = 0;
//        TString strMerge[nMerge]={};
//        std::vector<Int_t> mergeVec[nMerge];
//        std::vector<Int_t>::iterator it;
//        for(Int_t i=0; i<nSets; i++){

//        }

//*********************************************************************************************************************************
//*********************************************************************************************************************************
//*********************************************************************************************************************************

//    const Int_t nSets = 10;
//    const Int_t nData = 5;
//    TString DataSets[nSets]={
//      "LHC10b", "LHC10c", "LHC10d", "LHC10e", "LHC10f",
//      "LHC14j4b", "LHC14j4c", "LHC14j4d", "LHC14j4e", "LHC14j4f",
//    };


//    TString train = "Legotrain-vAN-20160803-7TeV_Omega";

//    TString runlist[nSets] = {
//      "merge_runlist_3","merge_runlist_3","merge_runlist_3","merge_runlist_3","merge_runlist_3",
//      "merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2"
//    };

//    Int_t trainRuns[nSets] = {
//      1745,1746,1747,1748,1749,
//      2368,2369,2370,2371,2372
//    };
//    const Int_t nFiles = 15;
//    TString Files[nFiles] = {"GammaConvNeutralMesonPiPlPiMiPiZero_0_9","GammaConvNeutralMesonPiPlPiMiPiZero_0_20",
//                             "GammaConvNeutralMesonPiPlPiMiPiZero_0_21", "GammaConvNeutralMesonPiPlPiMiPiZero_0_22",
//                             "GammaConvNeutralMesonPiPlPiMiPiZero_0_23", "GammaConvNeutralMesonPiPlPiMiPiZero_0_24",
//                             "GammaConvNeutralMesonPiPlPiMiPiZero_0_25", "GammaConvNeutralMesonPiPlPiMiPiZero_0_26",
//                             "GammaConvNeutralMesonPiPlPiMiPiZero_0_27",
//                             "OmegaToPiZeroGamma_1","OmegaToPiZeroGamma_101",
//                             "OmegaToPiZeroGamma_201","OmegaToPiZeroGamma_301",
//                             "OmegaToPiZeroGamma_401","OmegaToPiZeroGamma_501"};

//    const Int_t nMerge = 2;
//    TString strMerge[nMerge]={"LHC10","LHC14j4"};
//    std::vector<Int_t> mergeVec[nMerge];
//    std::vector<Int_t>::iterator it;
//    for(Int_t i=0; i<nSets; i++){
//      if(0<=i && i<=4) mergeVec[0].push_back(i);
//      if(5<=i && i<=9) mergeVec[1].push_back(i);
//    }

//    const Int_t nSets = 10;
//    const Int_t nData = 5;
//    TString DataSets[nSets]={
//      "LHC10b", "LHC10c", "LHC10d", "LHC10e", "LHC10f",
//      "LHC14j4b", "LHC14j4c", "LHC14j4d", "LHC14j4e", "LHC14j4f",
//    };

//    TString train = "Legotrain-vAN-20160816-7TeV_Omega";

//    TString runlist[nSets] = {
//      "merge_runlist_3","merge_runlist_3","merge_runlist_3","merge_runlist_3","merge_runlist_3",
//      "merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2"
//    };

//    Int_t trainRuns[nSets] = {
//      1800,1801,1802,1803,1804,
//      2451,2452,2453,2454,2455
//    };
//    const Int_t nFiles = 10;
//    TString Files[nFiles] = {"GammaConvNeutralMesonPiPlPiMiPiZero_0_9","GammaConvNeutralMesonPiPlPiMiPiZero_0_20",
//                             "GammaConvNeutralMesonPiPlPiMiPiZero_0_21", "GammaConvNeutralMesonPiPlPiMiPiZero_0_22",
//                             "GammaConvNeutralMesonPiPlPiMiPiZero_0_23", "GammaConvNeutralMesonPiPlPiMiPiZero_0_24",
//                             "GammaConvNeutralMesonPiPlPiMiPiZero_0_25", "GammaConvNeutralMesonPiPlPiMiPiZero_2_9",
//                             "GammaConvNeutralMesonPiPlPiMiPiZero_2_24", "GammaConvNeutralMesonPiPlPiMiPiZero_2_25"};

//    const Int_t nMerge = 2;
//    TString strMerge[nMerge]={"LHC10","LHC14j4"};
//    std::vector<Int_t> mergeVec[nMerge];
//    std::vector<Int_t>::iterator it;
//    for(Int_t i=0; i<nSets; i++){
//      if(0<=i && i<=4) mergeVec[0].push_back(i);
//      if(5<=i && i<=9) mergeVec[1].push_back(i);
//    }

//*********************************************************************************************************************************
//*********************************************************************************************************************************
//*********************************************************************************************************************************

//    const Int_t nSets = 6;
//    const Int_t nData = 6;
//    TString DataSets[nSets]={
//      "LHC10b", "LHC10c", "LHC10d", "LHC10e", "LHC10f", "LHC10c_900GeV"
//    };

//    TString train = "Legotrain-vAN-20161219-7TeV-std_sys_omega";

//    TString runlist[nSets] = {
//      "merge_runlist_3","merge_runlist_3","merge_runlist_3","merge_runlist_3","merge_runlist_3","merge_runlist_4"
//    };

//    Int_t trainRuns[nSets] = {
//      1980,1985,1986,1987,1988,1985
//    };
//    const Int_t nFiles = 31;
//    TString Files[nFiles] = {
//      "GammaCalo_201","GammaCalo_202","GammaCalo_203","GammaCalo_204","GammaCalo_205",
//      "GammaCalo_206","GammaCalo_207","GammaCalo_208","GammaCalo_209","GammaCalo_210",
//      "GammaCalo_211","GammaCalo_212",
//      "GammaConvCalo_201","GammaConvCalo_202","GammaConvCalo_203","GammaConvCalo_204","GammaConvCalo_205",
//      "GammaConvCalo_206","GammaConvCalo_207","GammaConvCalo_209","GammaConvCalo_210",
//      "GammaConvCalo_212","GammaConvCalo_213","GammaConvCalo_214","GammaConvCalo_215","GammaConvCalo_216",
//      "GammaConvCalo_217","GammaConvCalo_218","GammaConvCalo_219","GammaConvCalo_220","GammaConvV1_160"
//    };

//    const Int_t nMerge = 2;
//    TString strMerge[nMerge]={"LHC10","LHC10_900GeV"};
//    std::vector<Int_t> mergeVec[nMerge];
//    std::vector<Int_t>::iterator it;
//    for(Int_t i=0; i<nSets; i++){
//      if(0<=i && i<=4) mergeVec[0].push_back(i);
//      if(i==5) mergeVec[1].push_back(i);
//    }

//*********************************************************************************************************************************
//*********************************************************************************************************************************
//*********************************************************************************************************************************

//    const Int_t nSets = 6;
//    const Int_t nData = 6;
//    TString DataSets[nSets]={
//      "LHC10b", "LHC10c", "LHC10d", "LHC10e", "LHC10f", "LHC10c_900GeV"
//    };

//    TString train = "Legotrain-vAN-20161219-7TeV-std_sys_omega";

//    TString runlist[nSets] = {
//      "merge_runlist_3","merge_runlist_3","merge_runlist_3","merge_runlist_3","merge_runlist_3","merge_runlist_4"
//    };

//    Int_t trainRuns[nSets] = {
//      1968,1981,1982,1983,1984,1981
//    };
//    const Int_t nFiles = 6;
//    TString Files[nFiles] = {
//      "GammaConvNeutralMesonPiPlPiMiPiZero_0_9","GammaConvNeutralMesonPiPlPiMiPiZero_0_22",
//      "GammaConvNeutralMesonPiPlPiMiPiZero_1_9","GammaConvNeutralMesonPiPlPiMiPiZero_1_22",
//      "GammaConvNeutralMesonPiPlPiMiPiZero_2_9","GammaConvNeutralMesonPiPlPiMiPiZero_2_22"
//    };

//    const Int_t nMerge = 2;
//    TString strMerge[nMerge]={"LHC10","LHC10_900GeV"};
//    std::vector<Int_t> mergeVec[nMerge];
//    std::vector<Int_t>::iterator it;
//    for(Int_t i=0; i<nSets; i++){
//      if(0<=i && i<=4) mergeVec[0].push_back(i);
//      if(i==5) mergeVec[1].push_back(i);
//    }

//*********************************************************************************************************************************
//*********************************************************************************************************************************
//*********************************************************************************************************************************

//    const Int_t nSets = 6;
//    const Int_t nData = 0;
//    TString DataSets[nSets]={
//      "LHC14j4b", "LHC14j4c", "LHC14j4d", "LHC14j4e", "LHC14j4f", "LHC14j4c_900GeV",
//    };

//    TString train = "Legotrain-vAN-20161219-7TeV-std_sys_omega";

//    TString runlist[nSets] = {
//      "merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_4"
//    };

//    Int_t trainRuns[nSets] = {
//      2740,2741,2742,2743,2744,2741
//    };
//    const Int_t nFiles = 14;
//    TString Files[nFiles] = {
//      "GammaCalo_201","GammaCalo_202","GammaCalo_203","GammaCalo_204",
//      "GammaCalo_205","GammaCalo_206","GammaCalo_207","GammaCalo_208",
//      "GammaConvNeutralMesonPiPlPiMiPiZero_0_9","GammaConvNeutralMesonPiPlPiMiPiZero_0_24",
//      "GammaConvNeutralMesonPiPlPiMiPiZero_1_9","GammaConvNeutralMesonPiPlPiMiPiZero_1_22",
//      "GammaConvNeutralMesonPiPlPiMiPiZero_2_9","GammaConvNeutralMesonPiPlPiMiPiZero_2_22"
//    };

//    const Int_t nMerge = 2;
//    TString strMerge[nMerge]={"LHC14j4","LHC14j4_900GeV"};
//    std::vector<Int_t> mergeVec[nMerge];
//    std::vector<Int_t>::iterator it;
//    for(Int_t i=0; i<nSets; i++){
//      if(0<=i && i<=4) mergeVec[0].push_back(i);
//      if(i==5) mergeVec[1].push_back(i);
//    }

//    const Int_t nSets = 6;
//    const Int_t nData = 0;
//    TString DataSets[nSets]={
//      "LHC14j4b", "LHC14j4c", "LHC14j4d", "LHC14j4e", "LHC14j4f", "LHC14j4c_900GeV",
//    };

//    TString train = "Legotrain-vAN-20161219-7TeV-std_sys_omega";

//    TString runlist[nSets] = {
//      "merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_4"
//    };

//    Int_t trainRuns[nSets] = {
//      2745,2746,2747,2748,2749,2746
//    };
//    const Int_t nFiles = 9;
//    TString Files[nFiles] = {
//      "GammaCalo_209","GammaCalo_210","GammaCalo_211","GammaCalo_221",
//      "GammaConvCalo_202","GammaConvCalo_203","GammaConvCalo_204","GammaConvCalo_205",
//      "GammaConvCalo_206"
//    };

//    const Int_t nMerge = 2;
//    TString strMerge[nMerge]={"LHC14j4","LHC14j4_900GeV"};
//    std::vector<Int_t> mergeVec[nMerge];
//    std::vector<Int_t>::iterator it;
//    for(Int_t i=0; i<nSets; i++){
//      if(0<=i && i<=4) mergeVec[0].push_back(i);
//      if(i==5) mergeVec[1].push_back(i);
//    }

//    const Int_t nSets = 6;
//    const Int_t nData = 0;
//    TString DataSets[nSets]={
//      "LHC14j4b", "LHC14j4c", "LHC14j4d", "LHC14j4e", "LHC14j4f", "LHC14j4c_900GeV",
//    };

//    TString train = "Legotrain-vAN-20161219-7TeV-std_sys_omega";

//    TString runlist[nSets] = {
//      "merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_4"
//    };

//    Int_t trainRuns[nSets] = {
//      2750,2751,2752,2753,2754,2751
//    };
//    const Int_t nFiles = 15;
//    TString Files[nFiles] = {
//      "GammaCalo_200","GammaConvCalo_200","GammaConvCalo_201","GammaConvCalo_209",
//      "GammaConvCalo_210","GammaConvCalo_212","GammaConvCalo_213","GammaConvCalo_214",
//      "GammaConvCalo_215","GammaConvCalo_216","GammaConvCalo_217","GammaConvCalo_218",
//      "GammaConvCalo_219","GammaConvCalo_220","GammaConvV1_93"
//    };

//    const Int_t nMerge = 2;
//    TString strMerge[nMerge]={"LHC14j4","LHC14j4_900GeV"};
//    std::vector<Int_t> mergeVec[nMerge];
//    std::vector<Int_t>::iterator it;
//    for(Int_t i=0; i<nSets; i++){
//      if(0<=i && i<=4) mergeVec[0].push_back(i);
//      if(i==5) mergeVec[1].push_back(i);
//    }

    //*********************************************************************************************************************************
    //*********************************************************************************************************************************
    //*********************************************************************************************************************************

//    const Int_t nSets = 12;
//    const Int_t nData = 6;
//    TString DataSets[nSets]={
//      "LHC10b", "LHC10c", "LHC10d", "LHC10e", "LHC10f", "LHC10c_900GeV",
//      "LHC14j4b", "LHC14j4c", "LHC14j4d", "LHC14j4e", "LHC14j4f", "LHC14j4c_900GeV"
//    };

//    TString train = "Legotrain-vAN-2017026-7TeV-IterateECalib";

//    TString runlist[nSets] = {
//      "merge_runlist_3","merge_runlist_3","merge_runlist_3","merge_runlist_3","merge_runlist_3","merge_runlist_4",
//      "merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_4"
//    };

//    Int_t trainRuns[nSets] = {
//      2196,2197,2198,2199,2200,2197,
//      3081,3082,3083,3084,3085,3082
//    };
//    const Int_t nFiles = 4;
//    TString Files[nFiles] = {
//      "GammaConvCalo_209","GammaCalo_209","GammaConvCalo_210","GammaCalo_210"
//    };

//    const Int_t nMerge = 4;
//    TString strMerge[nMerge]={"LHC10","LHC10_900GeV","LHC14j4","LHC14j4_900GeV"};
//    std::vector<Int_t> mergeVec[nMerge];
//    std::vector<Int_t>::iterator it;
//    for(Int_t i=0; i<nSets; i++){
//      if(0<=i && i<=4) mergeVec[0].push_back(i);
//      if(i==5) mergeVec[1].push_back(i);
//      if(6<=i && i<=10) mergeVec[2].push_back(i);
//      if(i==11) mergeVec[3].push_back(i);
//    }

    //*********************************************************************************************************************************
    //*********************************************************************************************************************************
    //*********************************************************************************************************************************

//    const Int_t nSets = 6;
//    const Int_t nData = 0;
//    TString DataSets[nSets]={
//      "LHC14j4b", "LHC14j4c", "LHC14j4d", "LHC14j4e", "LHC14j4f", "LHC14j4c_900GeV"
//    };

//    TString train = "Legotrain-vAN-20170329-7TeV-bugfixForMC";

//    TString runlist[nSets] = {
//      "merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_4"
//    };

//    Int_t trainRuns[nSets] = {
//      2831,2832,2875,2834,2876,2832
//    };
//    const Int_t nFiles = 6;
//    TString Files[nFiles] = {
//      "GammaCalo_200","GammaCalo_201", "GammaConvV1_92", "GammaConvV1_93",
//      "GammaConvCalo_200","GammaConvCalo_201"
//    };

//    const Int_t nMerge = 2;
//    TString strMerge[nMerge]={"LHC14j4","LHC14j4_900GeV"};
//    std::vector<Int_t> mergeVec[nMerge];
//    std::vector<Int_t>::iterator it;
//    for(Int_t i=0; i<nSets; i++){
//      if(0<=i && i<=4) mergeVec[0].push_back(i);
//      if(i==5) mergeVec[1].push_back(i);
//    }

    //*********************************************************************************************************************************
    //*********************************************************************************************************************************
    //*********************************************************************************************************************************

//    const Int_t nSets = 4;
//    const Int_t nData = 4;
//    TString DataSets[nSets]={
//      "LHC10c_100", "LHC10c_900GeV_100",
//      "LHC10c_200", "LHC10c_900GeV_200"
//    };

//    TString train = "Legotrain-vAN-20170317-7TeV-sys-time";

//    TString runlist[nSets] = {
//      "merge_runlist_3","merge_runlist_4",
//      "merge_runlist_3","merge_runlist_4"
//    };

//    Int_t trainRuns[nSets] = {
//      2037,2038,2037,2038
//    };
//    const Int_t nFiles = 2;
//    TString Files[nFiles] = {
//      "GammaCalo_201",
//      "GammaConvCalo_201"
//    };

//    const Int_t nMerge = 4;
//    TString strMerge[nMerge]={"LHC10_100","LHC10_900GeV_100","LHC10_200","LHC10_900GeV_200"};
//    std::vector<Int_t> mergeVec[nMerge];
//    std::vector<Int_t>::iterator it;
//    for(Int_t i=0; i<nSets; i++){
//      if(i==0) mergeVec[0].push_back(i);
//      if(i==1) mergeVec[1].push_back(i);
//      if(i==2) mergeVec[2].push_back(i);
//      if(i==3) mergeVec[3].push_back(i);
//    }

//*********************************************************************************************************************************
//*********************************************************************************************************************************
//*********************************************************************************************************************************

 //convcalo systematics trigger EMC7
//    const Int_t nSets = 6;
//    const Int_t nData = 6;
//    TString DataSets[nSets]={
//      "LHC12b-kEMCEGA", "LHC12c-kEMCEGA", "LHC12d-kEMCEGA", "LHC12f-kEMCEGA", "LHC12h-kEMCEGA", "LHC12i-kEMCEGA"
//    };

//    TString train = "Legotrain-vAN-20160617-8TeV-Systematics_ConvCalo_EMC7";
//    Int_t trainRuns[nSets] = {
//                              1648,1649,1650,1651,1652,1653
//                             };

//    TString runlist[nSets] = {
//                              "merge_runlist_11","merge_runlist_11","merge_runlist_11","merge_runlist_11","merge_runlist_11","merge_runlist_11"
//                             };

//    const Int_t nFiles = 16;
//    TString Files[nFiles] = { "GammaConvCalo_132", "GammaConvCalo_133", "GammaConvCalo_134", "GammaConvCalo_135", "GammaConvCalo_136",
//                              "GammaConvCalo_139", "GammaConvCalo_140", "GammaConvCalo_142", "GammaConvCalo_143", "GammaConvCalo_144",
//                              "GammaConvCalo_145", "GammaConvCalo_146", "GammaConvCalo_147", "GammaConvCalo_148", "GammaConvCalo_149",
//                              "GammaConvCalo_150"
//                            };


//    const Int_t nMerge = 1;
//    TString strMerge[nMerge]={"LHC12-kEMCEGA"};
//    std::vector<Int_t> mergeVec[nMerge];
//    std::vector<Int_t>::iterator it;
//    for(Int_t i=0; i<nSets; i++){
//      if(0<=i && i<=5) mergeVec[0].push_back(i);
//    }

//*********************************************************************************************************************************
//*********************************************************************************************************************************
//*********************************************************************************************************************************

// convcalo systematics trigger EGA
//        const Int_t nSets = 5;
//        const Int_t nData = 5;
//        TString DataSets[nSets]={
//          "LHC12c-kEMCEGA", "LHC12d-kEMCEGA", "LHC12f-kEMCEGA", "LHC12h-kEMCEGA", "LHC12i-kEMCEGA"
//        };

//        TString train = "Legotrain-vAN-20160617-8TeV-Systematics_ConvCalo_EGA";
//        Int_t trainRuns[nSets] = {
//                                  1649,1650,1651,1652,1653
//                                 };

//        TString runlist[nSets] = {
//                                  "merge_runlist_11","merge_runlist_11","merge_runlist_11","merge_runlist_11","merge_runlist_11"
//                                 };

//        const Int_t nFiles = 16;
//        TString Files[nFiles] = { "GammaConvCalo_162", "GammaConvCalo_163", "GammaConvCalo_164", "GammaConvCalo_165", "GammaConvCalo_166",
//                                  "GammaConvCalo_169", "GammaConvCalo_170", "GammaConvCalo_172", "GammaConvCalo_173", "GammaConvCalo_174",
//                                  "GammaConvCalo_175", "GammaConvCalo_176", "GammaConvCalo_177", "GammaConvCalo_178", "GammaConvCalo_179",
//                                  "GammaConvCalo_180"
//                                };


//        const Int_t nMerge = 1;
//        TString strMerge[nMerge]={"LHC12-kEMCEGA"};
//        std::vector<Int_t> mergeVec[nMerge];
//        std::vector<Int_t>::iterator it;
//        for(Int_t i=0; i<nSets; i++){
//          if(0<=i && i<=4) mergeVec[0].push_back(i);
//        }

//*********************************************************************************************************************************
//*********************************************************************************************************************************
//*********************************************************************************************************************************

     //calo systematics trigger EMC7
//        const Int_t nSets = 6;
//        const Int_t nData = 6;
//        TString DataSets[nSets]={
//          "LHC12b-kEMCEGA", "LHC12c-kEMCEGA", "LHC12d-kEMCEGA", "LHC12f-kEMCEGA", "LHC12h-kEMCEGA", "LHC12i-kEMCEGA"
//        };

//        TString train = "Legotrain-vAN-20160617-8TeV-Systematics_Calo_EMC7";
//        Int_t trainRuns[nSets] = {
//                                  1613,1614,1615,1616,1617,1618
//                                 };

//        TString runlist[nSets] = {
//                                  "merge_runlist_11","merge_runlist_11","merge_runlist_11","merge_runlist_11","merge_runlist_11","merge_runlist_11"
//                                 };

//        const Int_t nFiles = 9;
//        TString Files[nFiles] = { "GammaCalo_121", "GammaCalo_122", "GammaCalo_123", "GammaCalo_124", "GammaCalo_125",
//                                  "GammaCalo_126", "GammaCalo_127", "GammaCalo_128", "GammaCalo_129"
//                                };


//        const Int_t nMerge = 1;
//        TString strMerge[nMerge]={"LHC12-kEMCEGA"};
//        std::vector<Int_t> mergeVec[nMerge];
//        std::vector<Int_t>::iterator it;
//        for(Int_t i=0; i<nSets; i++){
//          if(0<=i && i<=5) mergeVec[0].push_back(i);
//        }

    //calo systematics trigger EGA
//        const Int_t nSets = 5;
//        const Int_t nData = 5;
//        TString DataSets[nSets]={
//          "LHC12c-kEMCEGA", "LHC12d-kEMCEGA", "LHC12f-kEMCEGA", "LHC12h-kEMCEGA", "LHC12i-kEMCEGA"
//        };

//        TString train = "Legotrain-vAN-20160617-8TeV-Systematics_Calo_EGA";
//        Int_t trainRuns[nSets] = {
//          1614,1615,1616,1617,1618
//        };

//        TString runlist[nSets] = {
//          "merge_runlist_11","merge_runlist_11","merge_runlist_11","merge_runlist_11","merge_runlist_11"
//        };

//        const Int_t nFiles = 9;
//        TString Files[nFiles] = { "GammaCalo_141",
//                                  "GammaCalo_142", "GammaCalo_143", "GammaCalo_144", "GammaCalo_145", "GammaCalo_146",
//                                  "GammaCalo_147", "GammaCalo_148", "GammaCalo_149"
//                                };


//        const Int_t nMerge = 1;
//        TString strMerge[nMerge]={"LHC12-kEMCEGA"};
//        std::vector<Int_t> mergeVec[nMerge];
//        std::vector<Int_t>::iterator it;
//        for(Int_t i=0; i<nSets; i++){
//          if(0<=i && i<=4) mergeVec[0].push_back(i);
//        }

//*********************************************************************************************************************************
//*********************************************************************************************************************************
//*********************************************************************************************************************************

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

//*********************************************************************************************************************************
//*********************************************************************************************************************************
//*********************************************************************************************************************************

//    const Int_t nSets = 9;
//    const Int_t nData = 9;
//    TString DataSets[nSets]={
//      "LHC12a_0501", "LHC12b_0501","LHC12c_0501",
//      "LHC12a_0408", "LHC12b_0408","LHC12c_0408",
//      "LHC12a_0518", "LHC12b_0518","LHC12c_0518"
//    };

//    TString train = "Legotrain-vAN-20160518-8TeV-DebuggingPhysicsSelec";
//    Int_t trainRuns[nSets] = {1569,1570,1571,1572,1573,1574,1575,1576,1577};

//    TString runlist[nSets] = {"merge_runlist_9","merge_runlist_9","merge_runlist_9",
//                              "merge_runlist_9","merge_runlist_9","merge_runlist_9",
//                              "merge_runlist_9","merge_runlist_9","merge_runlist_9"
//                             };

//    const Int_t nFiles = 2;
//    TString Files[nFiles] = {"GammaConvCalo_101","GammaCalo_101"};

//    const Int_t nMerge = 0;
//    TString strMerge[nMerge]={};
//    std::vector<Int_t> mergeVec[nMerge];
//    std::vector<Int_t>::iterator it;

//*********************************************************************************************************************************
//*********************************************************************************************************************************
//*********************************************************************************************************************************

//    const Int_t nSets = 16;
//    const Int_t nData = 16;
//    TString DataSets[nSets]={
//      "LHC11a", "LHC11a_A125","LHC11a_A150","LHC11a_A50","LHC11a_A75","LHC11a_S100_A50","LHC11a_S400","LHC11a_S600",
//      "LHC13g", "LHC13g_A125","LHC13g_A150","LHC13g_A50","LHC13g_A75","LHC13g_S100_A50","LHC13g_S400","LHC13g_S600"
//    };

//    TString train = "Legotrain-vAN-20160511-2.76TeV-TenderVariations";
//    Int_t trainRuns[nSets] = {1548,1549,1550,1551,1552,1553,1554,1555,1556,1557,1558,1559,1560,1561,1562,1563
//                             };

//    TString runlist[nSets] = {"merge_runlist_3","merge_runlist_3","merge_runlist_3","merge_runlist_3","merge_runlist_3","merge_runlist_3","merge_runlist_3","merge_runlist_3"
//                             ,"merge","merge","merge","merge","merge","merge","merge","merge"
//                             };

//    const Int_t nFiles = 2;
//    TString Files[nFiles] = {"GammaConvCalo_1","GammaConvCalo_40"};

//    const Int_t nMerge = 0;
//    TString strMerge[nMerge]={};
//    std::vector<Int_t> mergeVec[nMerge];
//    std::vector<Int_t>::iterator it;

//*********************************************************************************************************************************
//*********************************************************************************************************************************
//*********************************************************************************************************************************

    const Int_t nSets = 18;
    const Int_t nData = 18;
    TString DataSets[nSets]={
      "LHC12a", "LHC12b", "LHC12c", "LHC12d", "LHC12f", "LHC12h", "LHC12i",
      "LHC12b-kEMC7", "LHC12c-kEMC7", "LHC12d-kEMC7", "LHC12f-kEMC7", "LHC12h-kEMC7", "LHC12i-kEMC7",
      "LHC12c-kEMCEGA", "LHC12d-kEMCEGA", "LHC12f-kEMCEGA", "LHC12h-kEMCEGA", "LHC12i-kEMCEGA"
    };

    TString train = "Legotrain-vAN-20170727-8TeV-merged";
    Int_t trainRuns[nSets] = {
                              2206,2207,2208,2209,2210,2211,2212,
                              2207,2208,2209,2210,2211,2212,
                              2208,2209,2210,2211,2212
                             };

    TString runlist[nSets] = {
                              "merge_runlist_9","merge_runlist_9","merge_runlist_9","merge_runlist_9","merge_runlist_9","merge_runlist_9","merge_runlist_9",
                              "merge_runlist_10","merge_runlist_10","merge_runlist_10","merge_runlist_10","merge_runlist_10","merge_runlist_10",
                              "merge_runlist_11","merge_runlist_11","merge_runlist_11","merge_runlist_11","merge_runlist_11"
                             };

    const Int_t nFiles = 3;
    TString Files[nFiles] = {
      "GammaCaloMerged_114","GammaCaloMerged_116","GammaCaloMerged_119"
                            };

    const Int_t nMerge = 3;
    TString strMerge[nMerge]={"LHC12", "LHC12-kEMC7", "LHC12-kEMCEGA"};
    std::vector<Int_t> mergeVec[nMerge];
    std::vector<Int_t>::iterator it;
    for(Int_t i=0; i<nSets; i++){
      if(0<=i && i<=6) mergeVec[0].push_back(i);
      if(7<=i && i<=12) mergeVec[1].push_back(i);
      if(13<=i && i<=17) mergeVec[2].push_back(i);
    }

//*********************************************************************************************************************************
//*********************************************************************************************************************************
//*********************************************************************************************************************************

//    const Int_t nSets = 18;
//    const Int_t nData = 18;
//    TString DataSets[nSets]={
//      "LHC12a", "LHC12b", "LHC12c", "LHC12d", "LHC12f", "LHC12h", "LHC12i",
//      "LHC12b-kEMC7", "LHC12c-kEMC7", "LHC12d-kEMC7", "LHC12f-kEMC7", "LHC12h-kEMC7", "LHC12i-kEMC7",
//      "LHC12c-kEMCEGA", "LHC12d-kEMCEGA", "LHC12f-kEMCEGA", "LHC12h-kEMCEGA", "LHC12i-kEMCEGA"
//    };

//    TString train = "Legotrain-vAN-20170619-8TeV-EMCal_pastFuture";
//    Int_t trainRuns[nSets] = {
//                              2119,2120,2121,2122,2123,2124,2125,
//                              2120,2121,2122,2123,2124,2125,
//                              2121,2122,2123,2124,2125
//                             };

//    TString runlist[nSets] = {
//                              "merge_runlist_9","merge_runlist_9","merge_runlist_9","merge_runlist_9","merge_runlist_9","merge_runlist_9","merge_runlist_9",
//                              "merge_runlist_10","merge_runlist_10","merge_runlist_10","merge_runlist_10","merge_runlist_10","merge_runlist_10",
//                              "merge_runlist_11","merge_runlist_11","merge_runlist_11","merge_runlist_11","merge_runlist_11"
//                             };

//    const Int_t nFiles = 20;
//    TString Files[nFiles] = {
//      "GammaCalo_101", "GammaCalo_700", "GammaCalo_701", "GammaCalo_702","GammaCalo_703", "GammaCalo_704", "GammaCalo_705",
//      "GammaCalo_706", "GammaCalo_707", "GammaCalo_708",
//      "GammaConvCalo_101", "GammaConvCalo_700", "GammaConvCalo_701", "GammaConvCalo_702","GammaConvCalo_703", "GammaConvCalo_704", "GammaConvCalo_705",
//      "GammaConvCalo_706", "GammaConvCalo_707", "GammaConvCalo_708",
//                            };

//    const Int_t nMerge = 3;
//    TString strMerge[nMerge]={"LHC12", "LHC12-kEMC7", "LHC12-kEMCEGA"};
//    std::vector<Int_t> mergeVec[nMerge];
//    std::vector<Int_t>::iterator it;
//    for(Int_t i=0; i<nSets; i++){
//      if(0<=i && i<=6) mergeVec[0].push_back(i);
//      if(7<=i && i<=12) mergeVec[1].push_back(i);
//      if(13<=i && i<=17) mergeVec[2].push_back(i);
//    }

//*********************************************************************************************************************************
//*********************************************************************************************************************************
//*********************************************************************************************************************************

//    const Int_t nSets = 14;
//    const Int_t nData = 0;
//    TString DataSets[nSets]={
//      "LHC15h1a1", "LHC15h1b", "LHC15h1c", "LHC15h1d", "LHC15h1f", "LHC15h1h", "LHC15h1i",
//      "LHC15h2a", "LHC15h2b", "LHC15h2c", "LHC15h2d", "LHC15h2f", "LHC15h2h", "LHC15h2i"
//    };

//    TString train = "Legotrain-vAN-20170601-8TeV-EMCal_openAngleStudies";
//    Int_t trainRuns[nSets] = {
//      2955,2956,2957,2958,2959,2960,2961,
//      2962,2963,2964,2965,2966,2967,2968
//    };

//    TString runlist[nSets] = {
//      "merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2",
//      "merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2"
//    };

//    const Int_t nFiles = 5;
//    TString Files[nFiles] = {
//      "GammaCalo_114","GammaCalo_115", "GammaCalo_116", "GammaCalo_117", "GammaConvCalo_131"
//    };

//    const Int_t nMerge = 3;
//    TString strMerge[nMerge]={ "LHC15h","LHC15h1","LHC15h2"};
//    std::vector<Int_t> mergeVec[nMerge];
//    std::vector<Int_t>::iterator it;
//    for(Int_t i=0; i<nSets; i++){
//      if(0<=i && i<=13) mergeVec[0].push_back(i);
//      if(0<=i && i<=6) mergeVec[1].push_back(i);
//      if(7<=i && i<=13) mergeVec[2].push_back(i);
//    }

//*********************************************************************************************************************************
//*********************************************************************************************************************************
//*********************************************************************************************************************************

//    const Int_t nSets = 14;
//    const Int_t nData = 0;
//    TString DataSets[nSets]={
//      "LHC15h1a1", "LHC15h1b", "LHC15h1c", "LHC15h1d", "LHC15h1f", "LHC15h1h", "LHC15h1i",
//      "LHC15h2a", "LHC15h2b", "LHC15h2c", "LHC15h2d", "LHC15h2f", "LHC15h2h", "LHC15h2i"
//    };

//    TString train = "Legotrain-vAN-20161219-8TeV-Variations_Tree_Calo";
//    Int_t trainRuns[nSets] = {
//      2755,2756,2757,2758,2759,2760,2761,
//      2762,2763,2764,2765,2766,2767,2768
//    };

//    TString runlist[nSets] = {
//      "merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2",
//      "merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2"
//    };

//    const Int_t nFiles = 1;
//    TString Files[nFiles] = {
//      "GammaCalo_119"
//    };

//    const Int_t nMerge = 1;
//    TString strMerge[nMerge]={ "LHC15h"};
//    std::vector<Int_t> mergeVec[nMerge];
//    std::vector<Int_t>::iterator it;
//    for(Int_t i=0; i<nSets; i++){
//      if(0<=i && i<=13) mergeVec[0].push_back(i);
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

//        TString train = "Legotrain-vAN-20161027-8TeV-stdMB_trMatch";
//        Int_t trainRuns[nSets] = {1887,1888,1889,1890,1891,1892,1893,
//                                  2576,2577,2578,2579,2580,2581,2582,
//                                  2583,2584,2585,2586,2587,2588,2589
//                                 };

//        TString runlist[nSets] = {"merge_runlist_9","merge_runlist_9","merge_runlist_9","merge_runlist_9","merge_runlist_9","merge_runlist_9","merge_runlist_9",
//                                  "merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2",
//                                  "merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2","merge_runlist_2"
//                                 };

//        const Int_t nFiles = 4;
//        TString Files[nFiles] = { "GammaConvCalo_101", "GammaConvCalo_108", "GammaConvCalo_110", "GammaConvV1_89"};

//        const Int_t nMerge = 4;
//        TString strMerge[nMerge]={"LHC12", "LHC15h1", "LHC15h2", "LHC15h"};
//        std::vector<Int_t> mergeVec[nMerge];
//        std::vector<Int_t>::iterator it;
//        for(Int_t i=0; i<nSets; i++){
//          if(0<=i && i<=6) mergeVec[0].push_back(i);
//          if(7<=i && i<=13) mergeVec[1].push_back(i);
//          if(14<=i && i<=20) mergeVec[2].push_back(i);
//          if(7<=i && i<=20) mergeVec[3].push_back(i);
//        }

    //*********************************************************************************************************************************
    //*********************************************************************************************************************************
    //*********************************************************************************************************************************

//    const Int_t nSets = 4;
//    const Int_t nData = 1;
//    TString DataSets[nSets]={
//      "LHC10h","LHC13d2","LHC13d2","LHC13d2b"
//    };

//    TString train = "Legotrain-vAN-20161021-2.76TeV-PbPbQA";
//    Int_t trainRuns[nSets] = {252, 312, 313, 314};

//    TString runlist[nSets] = {"merge", "merge", "merge", "merge"
//                             };

//    const Int_t nFiles = 7;
//    TString Files[nFiles] = {"GammaCalo_7","GammaCalo_9","GammaCalo_11",
//                             "GammaConvCalo_7","GammaConvCalo_9","GammaConvCalo_11",
//                             "GammaConvCalo_14"};

//    const Int_t nMerge = 0;
//    TString strMerge[nMerge]={};
//    std::vector<Int_t> mergeVec[nMerge];
//    std::vector<Int_t>::iterator it;
//    for(Int_t i=0; i<nSets; i++){

//    }

//    const Int_t nSets = 6;
//    const Int_t nData = 2;
//    TString DataSets[nSets]={
//      "LHC11h","LHC11h","LHC14a1a","LHC14a1b","LHC14a1c", "LHC14a1b"
//    };

//    TString train = "Legotrain-vAN-20161021-2.76TeV-PbPbQA_11h";
//    Int_t trainRuns[nSets] = {250, 251, 315, 316, 317, 318};

//    TString runlist[nSets] = {"merge_runlist_7", "merge_runlist_7", "merge_runlist_7", "merge_runlist_7", "merge_runlist_7", "merge_runlist_7"
//                             };

//    const Int_t nFiles = 10;
//    TString Files[nFiles] = {"GammaCalo_7","GammaCalo_9","GammaCalo_11", "GammaCalo_35",
//                             "GammaConvCalo_7","GammaConvCalo_9","GammaConvCalo_11", "GammaConvCalo_35",
//                             "GammaConvCalo_14", "GammaConvCalo_33"};

//    const Int_t nMerge = 0;
//    TString strMerge[nMerge]={};
//    std::vector<Int_t> mergeVec[nMerge];
//    std::vector<Int_t>::iterator it;
//    for(Int_t i=0; i<nSets; i++){

//    }

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
          if(Files[k].BeginsWith("OmegaToPiZeroGamma_")){
            if(i<nData) Files[k].Append("_0");
            else Files[k].Append("_1");
          }
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
            if(Files[k].BeginsWith("OmegaToPiZeroGamma_")){Files[k].Resize(Files[k].Length()-2);}
            continue;
          }

          if(copyAlien2Local(fPathGrid,fPathTemp)){
            ChangeStrucToStd(fPathTemp.Data(),fPathTemp.Data(),Files[k].Data());
            for(Int_t iM=0; iM<nMerge; iM++){
              it = find( mergeVec[iM].begin(), mergeVec[iM].end(), i);
              if( it!=mergeVec[iM].end() ) mergePeriod[iM][k] += Form(" %s",fPathTemp.Data());
            }
            gSystem->Exec(Form("ln -s %s %s/%s_%s.root",fPathTemp.Data(), fPathLocal.Data(), DataSets[i].Data(), Files[k].Data()));
            if(Files[k].BeginsWith("OmegaToPiZeroGamma_")){Files[k].Resize(Files[k].Length()-2);}
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
