/*******************************************************************************
 ******  provided by Gamma Conversion Group, PWGGA,                        *****
 ******     Daniel Muehlheim, d.muehlheim@cern.ch                          ***** 
 *******************************************************************************/

#include "QA.h"
#include "EventQA_Runwise.C"
#include "PhotonQA_Runwise.C"
#include "ClusterQA_Runwise.C"

void QA_Runwise(
                TString select          = "LHC11a",         // set selected
                Bool_t doEventQA        = kFALSE,           // switch on EventQA
                Bool_t doPhotonQA       = kFALSE,           // switch on PCM-PhotonQA
                Bool_t doClusterQA      = kFALSE,           // switch on ClusterQA
                Bool_t doMergedQA       = kFALSE,           // switch on merged ClusterQA
                Int_t mode              = 2,                // standard mode selector
                Int_t cutNr             = -1,               // if -1: you have to choose number at runtime
                Int_t doExtQA           = 2,                // 0: switched off, 1: normal extQA, 2: with Cell level plots
                TString suffix          = "eps"             // output format of plots
){
    //**************************************************************************************************************
    //******************* global settings **************************************************************************
    //**************************************************************************************************************

    TString folderRunlists          = "DownloadAndDataPrep/runlists";

    const Int_t maxSets             = 20;
    TString DataSets[maxSets];
    TString plotDataSets[maxSets];

    Int_t nSets                     = 0;
    Int_t nData                     = 0;
    TString fEnergyFlag;
    TString filePath, fileName;
    TString filePathPhoton          = "";
    Size_t markerSize               = 1.25;

    Bool_t doEquidistantXaxis       = kFALSE;
    Bool_t doTrigger                = kTRUE;
    Bool_t doHistsForEverySet       = kTRUE;
    Bool_t useDataRunListForMC      = kFALSE;
    Bool_t addSubFolder             = kFALSE;

    //**************************************************************************************************************
    //*************************** choose which data set to process *************************************************
    //**************************************************************************************************************
    if(select.CompareTo("LHC13g-11")==0){
        //LHC13g
        doEquidistantXaxis  = kTRUE;
        fEnergyFlag         = "2.76TeV";
        filePath            = "DataQA/20150601";
        fileName            = "GammaConvCalo_11.root";
        nSets               = 1;
        nData               = 1;
        DataSets[0]         = "LHC13g";
        plotDataSets[0]     = "LHC13g";
    }
    //**************************************************************************************************************
    else if(select.CompareTo("LHC13g-kEMCEG1")==0||select.CompareTo("LHC13g-kEMCEG2")==0){
        //LHC13g
        doEquidistantXaxis  = kTRUE;
        fEnergyFlag         = "2.76TeV";
        filePath            = "DataQA/20160302";
        fileName            = "GammaConvCalo_95.root";
        nSets               = 1;
        nData               = 1;
        if(select.CompareTo("LHC13g-kEMCEG1")==0){
            DataSets[0]     = "LHC13g-kEMCEG1";
            cutNr           = 0;
        } else if(select.CompareTo("LHC13g-kEMCEG2")==0){
            DataSets[0]     = "LHC13g-kEMCEG2";
            cutNr           = 1;
        }
        plotDataSets[0]     = "LHC13g";
    }
    //**************************************************************************************************************
    else if(select.CompareTo("LHC13g-kEMC7")==0){
        //LHC13g
        cutNr               = 1;
        doEquidistantXaxis  = kTRUE;
        fEnergyFlag         = "2.76TeV";
        filePath            = "DataQA/20160302";
        fileName            = "GammaConvCalo_96.root";
        nSets               = 1;
        nData               = 1;
        DataSets[0]         = "LHC13g-kEMC7";
        plotDataSets[0]     = "LHC13g";
    }
    //**************************************************************************************************************
    else if(select.CompareTo("LHC13g")==0){
        //LHC13g
        cutNr               = 0;
        doEquidistantXaxis  = kTRUE;
        fEnergyFlag         = "2.76TeV";
        filePath            = "DataQA/20160302";
        fileName            = "GammaConvCalo_96.root";
        nSets               = 2;
        nData               = 1;
        DataSets[0]         = "LHC13g";
        DataSets[1]         = "LHC15g2";
        plotDataSets[0]     = "LHC13g";
        plotDataSets[1]     = "Pythia8";
    }
    //**************************************************************************************************************
    else if(select.CompareTo("LHC13g-Calo")==0){
        //LHC13g
        mode                = 4;
        cutNr               = 0;
        doEquidistantXaxis  = kTRUE;
        fEnergyFlag         = "2.76TeV";
        filePath            = "DataQA/20160226";
        fileName            = "GammaCalo_60.root";
        nSets               = 2;
        nData               = 1;
        DataSets[0]         = "LHC13g";
        DataSets[1]         = "LHC15g2";
        plotDataSets[0]     = "LHC13g";
        plotDataSets[1]     = "Pythia8";
    }
    //**************************************************************************************************************
    else if(select.CompareTo("LHC11a-Calo")==0){
        //LHC11a + MC
        mode                = 4;
        cutNr               = 0;
        doEquidistantXaxis  = kTRUE;
        fEnergyFlag         = "2.76TeV";
        filePath            = "DataQA/20160226";
        fileName            = "GammaCalo_1.root";
        nSets               = 3;
        nData               = 1;
        DataSets[0]         = "LHC11a_p4_wSDD";
        DataSets[1]         = "LHC12f1a";
        DataSets[2]         = "LHC12f1b";
        plotDataSets[0]     = "LHC11a";
        plotDataSets[1]     = "Pythia8";
        plotDataSets[2]     = "Phojet";
    }
    //**************************************************************************************************************
    else if(select.CompareTo("LHC11a-Calo-kEMC1")==0){
        //LHC11a kEMC
        mode                = 4;
        cutNr               = 1;
        doEquidistantXaxis  = kTRUE;
        fEnergyFlag         = "2.76TeV";
        filePath            = "DataQA/20160226";
        fileName            = "GammaCalo_1.root";
        nSets               = 1;
        nData               = 1;
        DataSets[0]         = "LHC11a_p4_wSDD-kEMC1";
        plotDataSets[0]     = "LHC11a";
    }
    //**************************************************************************************************************
    else if(select.CompareTo("LHC13g-Calo-kEMCEG1")==0||select.CompareTo("LHC13g-Calo-kEMCEG2")==0){
        //LHC13g
        mode                = 4;
        doEquidistantXaxis  = kTRUE;
        fEnergyFlag         = "2.76TeV";
        filePath            = "DataQA/20160226";
        fileName            = "GammaCalo_60.root";
        nSets               = 1;
        nData               = 1;
        if(select.CompareTo("LHC13g-Calo-kEMCEG1")==0){
            DataSets[0]     = "LHC13g-kEMCEG1";
            cutNr           = 3;
        } else if(select.CompareTo("LHC13g-Calo-kEMCEG2")==0){
            DataSets[0]     = "LHC13g-kEMCEG2";
            cutNr           = 4;
        }
        plotDataSets[0]= "LHC13g";
    }
    //**************************************************************************************************************
    else if(select.CompareTo("LHC13g-Calo-kEMC7")==0){
        //LHC13g
        mode                = 4;
        cutNr               = 2;
        doEquidistantXaxis  = kTRUE;
        fEnergyFlag = "2.76TeV";
        filePath = "DataQA/20160226";
        fileName = "GammaCalo_60.root";
        nSets = 1;
        nData = 1;
        DataSets[0]= "LHC13g-kEMC7";
        plotDataSets[0]= "LHC13g";
    }
    //**************************************************************************************************************
    else if(select.CompareTo("LHC11a-kEMC1")==0){
        //LHC11a kEMC
        cutNr               = 1;
        doEquidistantXaxis  = kTRUE;
        fEnergyFlag         = "2.76TeV";
        filePath            = "DataQA/20151019";
        fileName            = "GammaConvCalo_1.root";
        nSets               = 1;
        nData               = 1;
        DataSets[0]         = "LHC11a_p4_wSDD-kEMC1";
        plotDataSets[0]     = "LHC11a";
    }
    //**************************************************************************************************************
    else if(select.CompareTo("LHC11a")==0){
        //LHC11a + MC
        cutNr               = 0;
        doEquidistantXaxis  = kTRUE;
        fEnergyFlag         = "2.76TeV";
        filePath            = "DataQA/20151019";
        fileName            = "GammaConvCalo_1.root";
        nSets               = 3;
        nData               = 1;
        DataSets[0]         = "LHC11a_p4_wSDD";
        DataSets[1]         = "LHC12f1a";
        DataSets[2]         = "LHC12f1b";
        plotDataSets[0]     = "LHC11a";
        plotDataSets[1]     = "Pythia8";
        plotDataSets[2]     = "Phojet";
        markerSize          = 1.5;
    }
    //**************************************************************************************************************
    else if(select.CompareTo("LHC11a+AddSig")==0){
        //LHC11a + MC + AddSignals
        cutNr               = 0;
        doEquidistantXaxis  = kTRUE;
        fEnergyFlag         = "2.76TeV";
        filePath            = "DataQA/20150825";
        fileName            = "GammaConvCalo_1.root";
        nSets               = 4;
        nData               = 1;
        DataSets[0]         = "LHC11a_p4_wSDD";
        DataSets[1]         = "LHC12f1a";
        DataSets[2]         = "LHC12f1b";
        DataSets[3]         = "LHC12i3";
        plotDataSets[0]     = "LHC11a";
        plotDataSets[1]     = "Pythia8";
        plotDataSets[2]     = "Phojet"; 
        plotDataSets[3]     = "Pythia8-AddedSignals";
    }
    //**************************************************************************************************************
    else if(select.CompareTo("LHC11aMerged")==0){
        //LHC11a + JJ MC for merged ana
        folderRunlists      = "DownloadAndDataPrep/runlists";
        doEquidistantXaxis  = kTRUE;
        fEnergyFlag         = "2.76TeV";
        filePath            = "DataQA";
        fileName            = "GammaCaloMerged_1.root";
        nSets               = 2;
        nData               = 1;
        DataSets[0]         = "LHC11a_pass4_wSDD";
        DataSets[1]         = "LHC15g1a";
        plotDataSets[0]     = "LHC11a";
        plotDataSets[1]     = "Pythia8-JJ";
    }
    //**************************************************************************************************************
    else if(select.CompareTo("LHC13gMerged")==0){
        //LHC11a + JJ MC for merged ana
        folderRunlists      = "DownloadAndDataPrep/runlists";
        doEquidistantXaxis  = kTRUE;
        fEnergyFlag         = "2.76TeV";
        filePath            = "DataQA";
        fileName            = "GammaCaloMerged_41.root";
        nSets               = 2;
        nData               = 1;
        DataSets[0]         = "LHC13g_pass1";
        DataSets[1]         = "LHC15a3a";
        plotDataSets[0]     = "LHC13g";
        plotDataSets[1]     = "Pythia8-JJ";
    }
    //**************************************************************************************************************
    else if(select.BeginsWith("LHC10") && select.Length()==6){
        //LHC10x
        cutNr               = 0;
        addSubFolder        = kTRUE;
        doEquidistantXaxis  = kTRUE;
        fEnergyFlag         = "7TeV";
        filePath            = "DataQA/20160206";
        fileName            = "GammaConvCalo_201.root";
        nSets               = 2;
        nData               = 1;
        TString temp        = select;
        temp                = temp.Remove(0,temp.Length()-1);
        plotDataSets[0]     = select;
        plotDataSets[1]     = "Pythia8";
        select              += "_pass4";
        DataSets[0]         =select;
        DataSets[1]         =Form("LHC14j4%s",temp.Data());
    }
    //**************************************************************************************************************
    else if(select.CompareTo("LHC10")==0){
        //LHC10
        cutNr                   = 0;
        doHistsForEverySet      = kFALSE;
        fEnergyFlag             = "7TeV";
        filePath                = "DataQA/20160206";
        fileName                = "GammaConvCalo_201.root";
        nSets                   = 6;
        nData                   = 5;
        TString dummyData1[6]   = {"LHC10b_pass4", "LHC10c_pass4", "LHC10d_pass4", "LHC10e_pass4", "LHC10f_pass4", "LHC14j4"};
        TString dummyData2[6]   = {"LHC10b", "LHC10c", "LHC10d", "LHC10e", "LHC10f", "Pythia8"};
        for (Int_t j = 0; j< nSets; j++){
            DataSets[j]     = dummyData1[j];
            plotDataSets[j] = dummyData2[j];
        }
    }
    //**************************************************************************************************************
    else if(select.CompareTo("LHC12")==0){
        //LHC12
        doHistsForEverySet  = kFALSE;
        fEnergyFlag         = "8TeV";
        filePath            = "DataQA/20160125";
        filePathPhoton      = "DataQA/20160104";
        fileName            = "GammaConvCalo_120.root";
        cutNr               = 0;
        nSets               = 9;
        nData               = 7;
        TString dummyData1[9]= {"LHC12a", "LHC12b", "LHC12c", "LHC12d", "LHC12f", "LHC12h", "LHC12i", "LHC15h1", "LHC15h2"};
        TString dummyData2[9]= {"LHC12a", "LHC12b", "LHC12c", "LHC12d", "LHC12f", "LHC12h", "LHC12i", "Pythia8", "Phojet"};
        for (Int_t j = 0; j< nSets; j++){
            DataSets[j]     = dummyData1[j];
            plotDataSets[j] = dummyData2[j];
        }    
    }
    //**************************************************************************************************************
    else if(select.CompareTo("LHC12-single")==0){
        //LHC12
        doHistsForEverySet  = kFALSE;
        fEnergyFlag         = "8TeV";
        filePath            = "DataQA/20160518";
        filePathPhoton      = "DataQA/20160518";
        fileName            = "GammaConvCalo_101.root";
        nSets               = 1;
        nData               = 1;
        DataSets[0]         = "LHC12";
        plotDataSets[0]     = "LHC12";
    }
    //**************************************************************************************************************
    else if(select.CompareTo("LHC12-Calo-single")==0){
        //LHC12
        doHistsForEverySet  = kFALSE;
        fEnergyFlag         = "8TeV";
        filePath            = "DataQA/20160518";
        filePathPhoton      = "DataQA/20160518";
        fileName            = "GammaCalo_101.root";
        nSets               = 1;
        nData               = 1;
        DataSets[0]         = "LHC12";
        plotDataSets[0]     = "LHC12";
    }
    //**************************************************************************************************************
    else if(select.CompareTo("LHC12-Data")==0){
        //LHC12, only Data
        doHistsForEverySet  = kFALSE;
        fEnergyFlag         = "8TeV";
        filePath            = "DataQA/20160125";
        filePathPhoton      = "DataQA/20160104";
        fileName            = "GammaConvCalo_120.root";
        nSets               = 7;
        nData               = 7;
        TString dummyData[7]= {"LHC12a", "LHC12b", "LHC12c", "LHC12d", "LHC12f", "LHC12h", "LHC12i"};
        for (Int_t j = 0; j< nSets; j++){
            DataSets[j]     = dummyData[j];
            plotDataSets[j] = dummyData[j];
        }    
    }
    //**************************************************************************************************************
    else if(select.CompareTo("LHC12-kINT8")==0){
        //LHC12, kINT8
        cutNr                   = 0;
        doHistsForEverySet      = kFALSE;
        fEnergyFlag             = "8TeV";
        filePath                = "DataQA/20160203";
        filePathPhoton          = "";
        fileName                = "GammaConvCalo_130.root";
        nSets                   = 9;
        nData                   = 9;
        TString dummyData1[9]   = {"LHC12a-kINT8", "LHC12b-kINT8", "LHC12c-kINT8", "LHC12d-kINT8", "LHC12e-kINT8", "LHC12f-kINT8", "LHC12g-kINT8", "LHC12h-kINT8", "LHC12i-kINT8"};
        TString dummyData2[9]   = {"LHC12a", "LHC12b", "LHC12c", "LHC12d", "LHC12e", "LHC12f", "LHC12g", "LHC12h", "LHC12i"};
        for (Int_t j = 0; j< nSets; j++){
            DataSets[j]     = dummyData1[j];
            plotDataSets[j] = dummyData2[j];
        }    
    }
    //**************************************************************************************************************
    else if(select.BeginsWith("LHC12") && select.Length()==6){
        //LHC12x
        useDataRunListForMC = kTRUE;
        addSubFolder        = kTRUE;
        doEquidistantXaxis  = kTRUE;
        fEnergyFlag         = "8TeV";
        filePath            = "DataQA/20151002";
        filePathPhoton      = "DataQA/20151018";
        fileName            = "GammaConvCalo_120.root";
        nSets               = 3;
        nData               = 1;
        DataSets[0]         = select;
        DataSets[1]         = "LHC14e2a";
        DataSets[2]         = "LHC14e2c";
        plotDataSets[0]     = select;
        plotDataSets[1]     = "Pythia8";
        plotDataSets[2]     = "Phojet";
    }
    //**************************************************************************************************************
    else if(select.BeginsWith("LHC12") && select.Length()==9 && select.EndsWith("-p2")){
        //LHC12x-p2
        addSubFolder        = kTRUE;
        doEquidistantXaxis  = kTRUE;
        cutNr               = 0;
        fEnergyFlag         = "8TeV";
        filePath            = "DataQA/20160125";
        filePathPhoton      = "DataQA/20160104";
        fileName            = "GammaConvCalo_120.root";
        select.Remove(select.Length()-3,select.Length());
        nSets               = 3;
        nData               = 1;
        TString tmpPeriod   = select;
        TString tmpA        = "";
        tmpPeriod.Remove(0,tmpPeriod.Length()-1);
        if(tmpPeriod.CompareTo("a")==0) 
            tmpA            = "1";
        DataSets[0]         = select;
        DataSets[1]         = Form("LHC15h1%s%s",tmpPeriod.Data(),tmpA.Data()); 
        DataSets[2]         = Form("LHC15h2%s",tmpPeriod.Data());
        plotDataSets[0]     = select;
        plotDataSets[1]     = "Pythia8";
        plotDataSets[2]     = "Phojet";
    }
    //**************************************************************************************************************
    else if(select.BeginsWith("LHC12") && select.EndsWith("-Data") && select.Length()==11){
        //LHC12x - only Data
        select.Remove(select.Length()-5,5);
        addSubFolder        =kTRUE;
        doEquidistantXaxis  = kTRUE;
        cutNr               = 0;
        fEnergyFlag         = "8TeV";
        filePath            = "DataQA/20160125";
        filePathPhoton      = "DataQA/20160104";
        fileName            = "GammaConvCalo_120.root";
        nSets               = 1;
        nData               = 1;
        DataSets[0]         = select;
        plotDataSets[0]     = select;
    }
    //**************************************************************************************************************
    else if(select.CompareTo("LHC12-MC")==0){
        //LHC12 MC
        fEnergyFlag         = "8TeV";
        filePath            = "DataQA/20160125";
        filePathPhoton      = "DataQA/20160104";
        fileName            = "GammaConvCalo_120.root";
        nSets               = 2;
        nData               = 0;
        DataSets[0]         = "LHC15h1";
        DataSets[1]         = "LHC15h2";
        plotDataSets[0]     = "Pythia8";
        plotDataSets[1]     = "Phojet";
    }
    //**************************************************************************************************************
    else if(select.BeginsWith("LHC12-kEMC7")){
        //LHC12 Trigger
        doHistsForEverySet  = kFALSE;
        cutNr               = 1;
        fEnergyFlag         = "8TeV";
        filePath            = "DataQA/20160518";
        fileName            = "GammaConvCalo_101.root";
        if(select.Contains("-Calo-")) fileName = "GammaCalo_101.root";
        if(select.EndsWith("-single")){
            nSets           = 1;
            nData           = 1;
            DataSets[0]     = "LHC12-kEMC7";
            plotDataSets[0] = "LHC12";
        } else{
            nSets                   = 6;
            nData                   = 6;
            TString dummyData1[6]   = {"LHC12b-kEMC7", "LHC12c-kEMC7", "LHC12d-kEMC7", "LHC12f-kEMC7", "LHC12h-kEMC7", "LHC12i-kEMC7"};
            TString dummyData2[6]   = {"LHC12b", "LHC12c", "LHC12d", "LHC12f", "LHC12h", "LHC12i"};
            for (Int_t j = 0; j< nSets; j++){
                DataSets[j]     = dummyData1[j];
                plotDataSets[j] = dummyData2[j];
            }    
        }
    }
    //**************************************************************************************************************
    else if(select.BeginsWith("LHC12-kEMCEGA")){
        //LHC12 Trigger
        doHistsForEverySet  = kFALSE;
        cutNr               = 2;
        fEnergyFlag         = "8TeV";
        filePath            = "DataQA/20160518";
        fileName            = "GammaConvCalo_101.root";
        if(select.Contains("-Calo-")) fileName = "GammaCalo_101.root";
        if(select.EndsWith("-single")){
            nSets           = 1;
            nData           = 1;
            DataSets[0]     = "LHC12-kEMCEGA";
            plotDataSets[0] = "LHC12";
        } else{
            nSets                   = 5;
            nData                   = 5;
            TString dummyData1[5]   = {"LHC12c-kEMCEGA", "LHC12d-kEMCEGA", "LHC12f-kEMCEGA", "LHC12h-kEMCEGA", "LHC12i-kEMCEGA"};
            TString dummyData2[5]   = {"LHC12c", "LHC12d", "LHC12f", "LHC12h", "LHC12i"};
            for (Int_t j = 0; j< nSets; j++){
                DataSets[j]         = dummyData1[j];
                plotDataSets[j]     = dummyData2[j];
            }    
        }
    }
    //**************************************************************************************************************
    else if(!select.BeginsWith("LHC12-kEMC") && select.BeginsWith("LHC12") && select.Contains("-kEMC")){
        //LHC12x - Trigger
        TString number      = "";
        if(select.EndsWith("-kEMC")) 
            number          = "121";
        if(select.EndsWith("-kEMCEGA")) 
            number          = "122";
        doEquidistantXaxis  = kTRUE;
        addSubFolder        = kTRUE;
        cutNr               = 0;
        fEnergyFlag         = "8TeV";
        filePath            = "DataQA/20160125";
        fileName            = Form("GammaConvCalo_%s.root",number.Data());
        nSets               = 1;
        nData               = 1;
        DataSets[0]         = select;
        TString selPlot     = select;
        selPlot.Resize(6);
        plotDataSets[0]=selPlot;
    }
    //**************************************************************************************************************
    else if(select.CompareTo("LHC12-kEMC8")==0){
        //LHC12 Trigger
        doHistsForEverySet  = kFALSE;
        cutNr               = 0;
        fEnergyFlag         = "8TeV";
        filePath            = "DataQA/20160203";
        fileName            = "GammaConvCalo_131.root";
        nSets               = 9;
        nData               = 9;
        TString dummyData1[9]   = {"LHC12a-kEMC8", "LHC12b-kEMC8", "LHC12c-kEMC8", "LHC12d-kEMC8", "LHC12e-kEMC8", "LHC12f-kEMC8", "LHC12g-kEMC8", "LHC12h-kEMC8", "LHC12i-kEMC8"};
        TString dummyData2[9]   = {"LHC12a", "LHC12b", "LHC12c", "LHC12d", "LHC12e", "LHC12f", "LHC12g", "LHC12h", "LHC12i"};
        for (Int_t j = 0; j< nSets; j++){
            DataSets[j]     = dummyData1[j];
            plotDataSets[j] = dummyData2[j];
        }    
    }
    //**************************************************************************************************************
    else if(select.CompareTo("LHC12-kEMC8EGA")==0){
        //LHC12 Trigger
        doHistsForEverySet  = kFALSE;
        cutNr               = 0;
        fEnergyFlag         = "8TeV";
        filePath            = "DataQA/20160203";
        fileName            = "GammaConvCalo_132.root";
        nSets               = 9;
        nData               = 9;
        TString dummyData1[9]   = {"LHC12a-kEMC8EGA", "LHC12b-kEMC8EGA", "LHC12c-kEMC8EGA", "LHC12d-kEMC8EGA", "LHC12e-kEMC8EGA", "LHC12f-kEMC8EGA", "LHC12g-kEMC8EGA", "LHC12h-kEMC8EGA", "LHC12i-kEMC8EGA"};
        TString dummyData2[9]   = {"LHC12a", "LHC12b", "LHC12c", "LHC12d", "LHC12e", "LHC12f", "LHC12g", "LHC12h", "LHC12i"};
        for (Int_t j = 0; j< nSets; j++){
            DataSets[j]     = dummyData1[j];
            plotDataSets[j] = dummyData2[j];
        }    
    }
    //**************************************************************************************************************
    else if(select.CompareTo("LHC13pPb-Data")==0){
        //LHC12
        doHistsForEverySet  = kTRUE;
        fEnergyFlag         = "pPb_5.023TeV";
        filePath            = "DataQA/20160307";
        filePathPhoton      = "DataQA/20160307";
        fileName            = "GammaConvCalo_20.root";
        cutNr               = 0;
        nSets               = 5;
        nData               = 5;
        TString dummyData1[5]   = {"LHC13b", "LHC13c", "LHC13d", "LHC13e", "LHC13f"};
        TString dummyData2[5]   = {"LHC13b", "LHC13c", "LHC13d", "LHC13e", "LHC13f"};
        for (Int_t j = 0; j< nSets; j++){
            DataSets[j]     = dummyData1[j];
            plotDataSets[j] = dummyData2[j];
        }    
        
    }
    //**************************************************************************************************************
    else if(select.CompareTo("LHC13pPb-MBData")==0){
        //pPb MB data
        doEquidistantXaxis  = kTRUE;
        doHistsForEverySet  = kFALSE;
        fEnergyFlag         = "pPb_5.023TeV";
        filePath            = "DataQA/20160307";
        filePathPhoton      = "DataQA/20160307";
        fileName            = "GammaConvCalo_20.root";
        cutNr               = 0;
        nSets               = 2;
        nData               = 2;
        DataSets[0]         = "LHC13b";
        DataSets[1]         = "LHC13c";
        plotDataSets[0]     = "LHC13b";
        plotDataSets[1]     = "LHC13c";
    }
    //**************************************************************************************************************
    else if(select.CompareTo("LHC13pPb-MB")==0){
        //pPb MB
        doEquidistantXaxis  = kTRUE;
        doHistsForEverySet  = kFALSE;
        fEnergyFlag         = "pPb_5.023TeV";
        filePath            = "DataQA/20160307";
        filePathPhoton      = "DataQA/20160307";
        fileName            = "GammaConvCalo_20.root";
        cutNr               = 0;
        nSets               = 4;
        nData               = 2;
        TString dummyData1[4]   = {"LHC13b", "LHC13c", "LHC13b2_efix", "LHC13e7"};
        TString dummyData2[4]   = {"LHC13b", "LHC13c", "DPMJET", "HIJING"};
        for (Int_t j = 0; j< nSets; j++){
            DataSets[j]     = dummyData1[j];
            plotDataSets[j] = dummyData2[j];
        }
    }
    //**************************************************************************************************************
    else if(select.CompareTo("LHC15f")==0){
	  mode=0;
	  fEnergyFlag = "13TeV";
	  nSets = 2;
	  nData = 1;
	  filePath= "/home/meike/analysis/data/GridOutput/GammaConv/pp/";
	  filePathPhoton = "/home/meike/analysis/data/GridOutput/PhotonQA/pp/";
	  fileName = "GammaConvV1.root";
	  DataSets[0]="LHC15f_ESD"; DataSets[1]="LHC15g3a3_ESD";
	  plotDataSets[0]="LHC15f_ESD"; plotDataSets[1]="LHC15g3a3_ESD";
	  markerSize=1.5;
	  doEquidistantXaxis=kTRUE;
    }
    //**************************************************************************************************************
    else if(select.CompareTo("LHC15h")==0){
	  cutNr=0;
	  mode=0;
	  fEnergyFlag = "13TeV";
	  nSets = 1;
	  nData = 1;
	  filePath= "/home/meike/analysis/data/GridOutput/GammaConv/pp/";
	  filePathPhoton = "/home/meike/analysis/data/GridOutput/PhotonQA/pp/";
	  fileName = "GammaConvV1.root";
	  DataSets[0]="LHC15h";
	  plotDataSets[0]="LHC15h";
	  markerSize=3;
	  doEquidistantXaxis=kTRUE;
    }
    //**************************************************************************************************************
    else if(select.CompareTo("LHC15i")==0){
	  cutNr=0;
	  mode=0;
	  fEnergyFlag = "13TeV";
	  nSets = 1;
	  nData = 1;
	  filePath= "/home/meike/analysis/data/GridOutput/GammaConv/pp/";
	  filePathPhoton = "/home/meike/analysis/data/GridOutput/PhotonQA/pp/";
	  fileName = "GammaConvV1.root";
	  DataSets[0]="LHC15i";
	  plotDataSets[0]="LHC15i";
	  markerSize=3;
	  doEquidistantXaxis=kTRUE;
    }
    //**************************************************************************************************************
    else if(select.CompareTo("LHC15g")==0){
	  cutNr=0;
	  mode=0;
	  fEnergyFlag = "13TeV";
	  nSets = 1;
	  nData = 1;
	  filePath= "/home/meike/analysis/data/GridOutput/GammaConv/pp/";
	  filePathPhoton = "/home/meike/analysis/data/GridOutput/PhotonQA/pp/";
	  fileName = "GammaConvV1.root";
	  DataSets[0]="LHC15g";
	  plotDataSets[0]="LHC15g";
	  markerSize=3;
	  doEquidistantXaxis=kTRUE;
    }
    //**************************************************************************************************************
    else if(select.CompareTo("LHC11h")==0){
      // PbPb 2.76TeV
      filePath= "/home/admin1/leardini/GridOutput/PbPb/Legotrain-vAN-20160801-1";
      filePathPhoton = "/home/admin1/leardini/GridOutput/PbPb/Legotrain-vAN-20160801-1";
      fileName = "GammaConvV1_226.root";
      mode=0;
      cutNr=2;
      suffix="pdf";
      fEnergyFlag = "PbPb_2.76TeV";
      nSets               = 4;
      nData               = 1;
      TString dummyDataSets[4]   = {"LHC11h", "LHC11h", "LHC14a1a", "LHC14a1b"};
      TString dummyPlotSets[4]   = {"0-10% LHC11h", "20-50% LHC11h", "0-10% LHC14a1a", "20-50% LHC14a1b"};
      for (Int_t j = 0; j< nSets; j++){
          DataSets[j]     = dummyDataSets[j];
          plotDataSets[j] = dummyPlotSets[j];
      }
      markerSize=1.5;
      doEquidistantXaxis=kTRUE;
    }
    //**************************************************************************************************************
    else if(select.CompareTo("LHC15o")==0){
      // PbPb 5.02TeV
      mode=0;
      fEnergyFlag = "PbPb_5.02TeV";
      nSets = 2;
      nData = 1;
      cutNr = 4;  // 0-100%
      filePath= "/home/meike/analysis/data/GridOutput/GammaConv/PbPb/";
      filePathPhoton = "";
      fileName = "GammaConvV1_246.root";
      DataSets[0]="LHC15o"; DataSets[1]="LHC16h4";
      plotDataSets[0]="LHC15o 0-100%"; plotDataSets[1]="LHC16h4 0-100%";
      doEquidistantXaxis=kTRUE;
    }
    //**************************************************************************************************************
    else{
        cout << "No valid selection! Returning..." << endl;
        return;
    }
    
    //**************************************************************************************************************
    //******************************  Starting individual QA macros ***********************************************
    //**************************************************************************************************************
    if(doEventQA) EventQA_Runwise(  nSets, nData, fEnergyFlag, filePath, fileName, DataSets, plotDataSets, mode, cutNr,
                                    doExtQA,doEquidistantXaxis, doTrigger, doHistsForEverySet, addSubFolder, useDataRunListForMC, markerSize, suffix, folderRunlists);
    if(doPhotonQA){
        TString path                = filePath;
        if(!filePathPhoton.IsNull()) 
            path                    = filePathPhoton;
        PhotonQA_Runwise(   nSets, nData, fEnergyFlag, path, fileName, DataSets, plotDataSets, mode, cutNr,
                            doExtQA, doEquidistantXaxis, doTrigger, doHistsForEverySet, addSubFolder, useDataRunListForMC, markerSize, suffix, folderRunlists);
    }
    if(doClusterQA) ClusterQA_Runwise(  nSets, nData, fEnergyFlag, filePath, fileName, DataSets, plotDataSets, mode, cutNr,
                                        doExtQA, doEquidistantXaxis, doTrigger, doHistsForEverySet, addSubFolder, useDataRunListForMC, markerSize, suffix, folderRunlists);
    if(doMergedQA) ClusterQA_Runwise(  nSets, nData, fEnergyFlag, filePath, fileName, DataSets, plotDataSets, mode, cutNr,
                                       1, doEquidistantXaxis, doTrigger, doHistsForEverySet, addSubFolder, useDataRunListForMC, markerSize, suffix, folderRunlists, kTRUE);
    return;
}
