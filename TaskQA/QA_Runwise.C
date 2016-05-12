#include "QA.h"
#include "EventQA_Runwise.C"
#include "PhotonQA_Runwise.C"
#include "ClusterQA_Runwise.C"

void QA_Runwise(
				TString select = "LHC11a",
				Bool_t doEventQA = kFALSE,
				Bool_t doPhotonQA = kFALSE,
				Bool_t doClusterQA = kFALSE,
				Int_t mode = 2,
				Int_t cutNr = -1,				// if -1: you have to choose number at runtime
				Int_t doExtQA = 2,				// 0: switched off, 1: normal extQA, 2: with Cell level plots
				TString suffix = "eps"
){
    TString folderRunlists = "DownloadAndDataPrep/runlists";

	const Int_t maxSets = 20;
	TString DataSets[maxSets];
	TString plotDataSets[maxSets];

	Int_t nSets = 0;
	Int_t nData = 0;
	TString fEnergyFlag;
	TString filePath, fileName;
	TString filePathPhoton = "";
    Size_t markerSize = 1.25;

	Bool_t doEquidistantXaxis = kFALSE;
	Bool_t doTrigger = kTRUE;
	Bool_t doHistsForEverySet = kTRUE;
	Bool_t useDataRunListForMC = kFALSE;
	Bool_t addSubFolder = kFALSE;

//choose which data set to process
//**************************************************************************************************************
	if(select.CompareTo("LHC13g-11")==0){
		//LHC13g
		doEquidistantXaxis=kTRUE;
		fEnergyFlag = "2.76TeV";
		filePath = "DataQA/20150601";
		fileName = "GammaConvCalo_11.root";
		nSets = 1;
		nData = 1;
		DataSets[0]="LHC13g";
		plotDataSets[0]="LHC13g";
	}
//**************************************************************************************************************
	else if(select.CompareTo("LHC13g-kEMCEG1")==0||select.CompareTo("LHC13g-kEMCEG2")==0){
		//LHC13g
		doEquidistantXaxis=kTRUE;
		fEnergyFlag = "2.76TeV";
        filePath = "DataQA/20160302";
		fileName = "GammaConvCalo_95.root";
		nSets = 1;
		nData = 1;
		if(select.CompareTo("LHC13g-kEMCEG1")==0){
			DataSets[0]="LHC13g-kEMCEG1";
			cutNr=0;
		}else if(select.CompareTo("LHC13g-kEMCEG2")==0){
			DataSets[0]="LHC13g-kEMCEG2";
			cutNr=1;
		}
		plotDataSets[0]="LHC13g";
	}
//**************************************************************************************************************
	else if(select.CompareTo("LHC13g-kEMC7")==0){
		//LHC13g
		cutNr=1;
		doEquidistantXaxis=kTRUE;
		fEnergyFlag = "2.76TeV";
        filePath = "DataQA/20160302";
		fileName = "GammaConvCalo_96.root";
		nSets = 1;
		nData = 1;
		DataSets[0]="LHC13g-kEMC7";
		plotDataSets[0]="LHC13g";
	}
//**************************************************************************************************************
	else if(select.CompareTo("LHC13g")==0){
		//LHC13g
		cutNr=0;
		doEquidistantXaxis=kTRUE;
		fEnergyFlag = "2.76TeV";
        filePath = "DataQA/20160302";
		fileName = "GammaConvCalo_96.root";
		nSets = 2;
		nData = 1;
		DataSets[0]="LHC13g";DataSets[1]="LHC15g2";
		plotDataSets[0]="LHC13g";plotDataSets[1]="Pythia8";
	}
    //**************************************************************************************************************
        else if(select.CompareTo("LHC13g-Calo")==0){
            //LHC13g
            mode = 4;
            cutNr = 0;
            doEquidistantXaxis=kTRUE;
            fEnergyFlag = "2.76TeV";
            filePath = "DataQA/20160226";
            fileName = "GammaCalo_60.root";
            nSets = 2;
            nData = 1;
            DataSets[0]="LHC13g";DataSets[1]="LHC15g2";
            plotDataSets[0]="LHC13g";plotDataSets[1]="Pythia8";
        }
    //**************************************************************************************************************
        else if(select.CompareTo("LHC11a-Calo")==0){
            //LHC11a + MC
            mode = 4;
            cutNr = 0;
            doEquidistantXaxis=kTRUE;
            fEnergyFlag = "2.76TeV";
            filePath = "DataQA/20160226";
            fileName = "GammaCalo_1.root";
            nSets = 3;
            nData = 1;
            DataSets[0]="LHC11a_p4_wSDD";DataSets[1]="LHC12f1a";DataSets[2]="LHC12f1b";
            plotDataSets[0]="LHC11a";plotDataSets[1]="Pythia8";plotDataSets[2]="Phojet";
        }
    //**************************************************************************************************************
        else if(select.CompareTo("LHC11a-Calo-kEMC1")==0){
            //LHC11a kEMC
            mode = 4;
            cutNr = 1;
            doEquidistantXaxis=kTRUE;
            fEnergyFlag = "2.76TeV";
            filePath = "DataQA/20160226";
            fileName = "GammaCalo_1.root";
            nSets = 1;
            nData = 1;
            DataSets[0]="LHC11a_p4_wSDD-kEMC1";
            plotDataSets[0]="LHC11a";
        }
    //**************************************************************************************************************
        else if(select.CompareTo("LHC13g-Calo-kEMCEG1")==0||select.CompareTo("LHC13g-Calo-kEMCEG2")==0){
            //LHC13g
            mode = 4;
            doEquidistantXaxis=kTRUE;
            fEnergyFlag = "2.76TeV";
            filePath = "DataQA/20160226";
            fileName = "GammaCalo_60.root";
            nSets = 1;
            nData = 1;
            if(select.CompareTo("LHC13g-Calo-kEMCEG1")==0){
                DataSets[0]="LHC13g-kEMCEG1";
                cutNr = 3;
            }else if(select.CompareTo("LHC13g-Calo-kEMCEG2")==0){
                DataSets[0]="LHC13g-kEMCEG2";
                cutNr = 4;
            }
            plotDataSets[0]="LHC13g";
        }
    //**************************************************************************************************************
        else if(select.CompareTo("LHC13g-Calo-kEMC7")==0){
            //LHC13g
            mode = 4;
            cutNr = 2;
            doEquidistantXaxis=kTRUE;
            fEnergyFlag = "2.76TeV";
            filePath = "DataQA/20160226";
            fileName = "GammaCalo_60.root";
            nSets = 1;
            nData = 1;
            DataSets[0]="LHC13g-kEMC7";
            plotDataSets[0]="LHC13g";
        }
//**************************************************************************************************************
	else if(select.CompareTo("LHC11a-kEMC1")==0){
		//LHC11a kEMC
		cutNr=1;
		doEquidistantXaxis=kTRUE;
		fEnergyFlag = "2.76TeV";
		filePath = "DataQA/20151019";
		fileName = "GammaConvCalo_1.root";
		nSets = 1;
		nData = 1;
		DataSets[0]="LHC11a_p4_wSDD-kEMC1";
		plotDataSets[0]="LHC11a";
	}
//**************************************************************************************************************
	else if(select.CompareTo("LHC11a")==0){
		//LHC11a + MC
		cutNr=0;
		doEquidistantXaxis=kTRUE;
		fEnergyFlag = "2.76TeV";
		filePath = "DataQA/20151019";
		fileName = "GammaConvCalo_1.root";
		nSets = 3;
		nData = 1;
		DataSets[0]="LHC11a_p4_wSDD";DataSets[1]="LHC12f1a";DataSets[2]="LHC12f1b";
		plotDataSets[0]="LHC11a";plotDataSets[1]="Pythia8";plotDataSets[2]="Phojet";
        markerSize = 1.5;
	}
//**************************************************************************************************************
	else if(select.CompareTo("LHC11a+AddSig")==0){
		//LHC11a + MC + AddSignals
		cutNr=0;
		doEquidistantXaxis=kTRUE;
		fEnergyFlag = "2.76TeV";
		filePath = "DataQA/20150825";
		fileName = "GammaConvCalo_1.root";
		nSets = 4;
		nData = 1;
		DataSets[0]="LHC11a_p4_wSDD";DataSets[1]="LHC12f1a";DataSets[2]="LHC12f1b";DataSets[3]="LHC12i3";
		plotDataSets[0]="LHC11a";plotDataSets[1]="Pythia8";plotDataSets[2]="Phojet";plotDataSets[3]="Pythia8-AddedSignals";
	}
//**************************************************************************************************************
	else if(select.BeginsWith("LHC10") && select.Length()==6){
		//LHC10x
		cutNr=0;
		addSubFolder=kTRUE;
		doEquidistantXaxis=kTRUE;
		fEnergyFlag="7TeV";
        filePath = "DataQA/20160206";
		fileName = "GammaConvCalo_201.root";
		nSets = 2;
		nData = 1;
		TString temp = select;
		temp = temp.Remove(0,temp.Length()-1);
		plotDataSets[0]=select;plotDataSets[1]="Pythia8";
		select+="_pass4";
		DataSets[0]=select;DataSets[1]=Form("LHC14j4%s",temp.Data());
	}
//**************************************************************************************************************
	else if(select.CompareTo("LHC10")==0){
		//LHC10
		cutNr=0;
		doHistsForEverySet = kFALSE;
		fEnergyFlag="7TeV";
        filePath = "DataQA/20160206";
		fileName = "GammaConvCalo_201.root";
        nSets = 6;
        nData = 5;
		DataSets[0]="LHC10b_pass4";DataSets[1]="LHC10c_pass4";DataSets[2]="LHC10d_pass4";DataSets[3]="LHC10e_pass4";DataSets[4]="LHC10f_pass4";
		DataSets[5]="LHC14j4";
		plotDataSets[0]="LHC10b";plotDataSets[1]="LHC10c";plotDataSets[2]="LHC10d";plotDataSets[3]="LHC10e";plotDataSets[4]="LHC10f";
		plotDataSets[5]="Pythia8";
	}
//**************************************************************************************************************
	else if(select.CompareTo("LHC12")==0){
		//LHC12
		doHistsForEverySet = kFALSE;
		fEnergyFlag="8TeV";
        filePath = "DataQA/20160125";
        filePathPhoton = "DataQA/20160104";
		fileName = "GammaConvCalo_120.root";
        cutNr=0;
        nSets = 9;
        nData = 7;
		DataSets[0]="LHC12a";DataSets[1]="LHC12b";DataSets[2]="LHC12c";DataSets[3]="LHC12d";DataSets[4]="LHC12f";
        /*DataSets[5]="LHC12g";*/DataSets[5]="LHC12h";DataSets[6]="LHC12i";DataSets[7]="LHC15h1";DataSets[8]="LHC15h2";
		plotDataSets[0]="LHC12a";plotDataSets[1]="LHC12b";plotDataSets[2]="LHC12c";plotDataSets[3]="LHC12d";plotDataSets[4]="LHC12f";
        /*plotDataSets[5]="LHC12g";*/plotDataSets[5]="LHC12h";plotDataSets[6]="LHC12i";plotDataSets[7]="Pythia8";plotDataSets[8]="Phojet";
	}
//**************************************************************************************************************
	else if(select.CompareTo("LHC12-single")==0){
		//LHC12
		doHistsForEverySet = kFALSE;
		fEnergyFlag="8TeV";
        filePath = "DataQA/20160125";
		//filePath = "DataQA/20150728";
        filePathPhoton = "DataQA/20160104";
		fileName = "GammaConvCalo_120.root";
		nSets = 1;
		nData = 1;
		DataSets[0]="LHC12";
		plotDataSets[0]="LHC12";
	}
//**************************************************************************************************************
	else if(select.CompareTo("LHC12-Data")==0){
		//LHC12, only Data
		doHistsForEverySet = kFALSE;
		fEnergyFlag="8TeV";
        filePath = "DataQA/20160125";
        filePathPhoton = "DataQA/20160104";
		fileName = "GammaConvCalo_120.root";
        nSets = 7;
        nData = 7;
		DataSets[0]="LHC12a";DataSets[1]="LHC12b";DataSets[2]="LHC12c";DataSets[3]="LHC12d";DataSets[4]="LHC12f";
        /*DataSets[5]="LHC12g";*/DataSets[5]="LHC12h";DataSets[6]="LHC12i";
		plotDataSets[0]="LHC12a";plotDataSets[1]="LHC12b";plotDataSets[2]="LHC12c";plotDataSets[3]="LHC12d";plotDataSets[4]="LHC12f";
        /*plotDataSets[5]="LHC12g";*/plotDataSets[5]="LHC12h";plotDataSets[6]="LHC12i";
	}
//**************************************************************************************************************
    else if(select.CompareTo("LHC12-kINT8")==0){
      //LHC12, kINT8
      cutNr=0;
      doHistsForEverySet = kFALSE;
      fEnergyFlag="8TeV";
      filePath = "DataQA/20160203";
      filePathPhoton = "";
      fileName = "GammaConvCalo_130.root";
      nSets = 9;
      nData = 9;
      DataSets[0]="LHC12a-kINT8";DataSets[1]="LHC12b-kINT8";DataSets[2]="LHC12c-kINT8";DataSets[3]="LHC12d-kINT8";DataSets[4]="LHC12e-kINT8";
      DataSets[5]="LHC12f-kINT8";DataSets[6]="LHC12g-kINT8";DataSets[7]="LHC12h-kINT8";DataSets[8]="LHC12i-kINT8";
      plotDataSets[0]="LHC12a";plotDataSets[1]="LHC12b";plotDataSets[2]="LHC12c";plotDataSets[3]="LHC12d";plotDataSets[4]="LHC12e";
      plotDataSets[5]="LHC12f";plotDataSets[6]="LHC12g";plotDataSets[7]="LHC12h";plotDataSets[8]="LHC12i";
    }
//**************************************************************************************************************
	else if(select.BeginsWith("LHC12") && select.Length()==6){
		//LHC12x
		useDataRunListForMC=kTRUE;
		addSubFolder=kTRUE;
		doEquidistantXaxis=kTRUE;
		fEnergyFlag="8TeV";
		filePath = "DataQA/20151002";
		filePathPhoton = "DataQA/20151018";
		//filePath = "DataQA/20150728";
		fileName = "GammaConvCalo_120.root";
		nSets = 3;
		nData = 1;
		DataSets[0]=select;DataSets[1]="LHC14e2a";DataSets[2]="LHC14e2c";
		plotDataSets[0]=select;plotDataSets[1]="Pythia8";plotDataSets[2]="Phojet";
	}
//**************************************************************************************************************
	else if(select.BeginsWith("LHC12") && select.Length()==9 && select.EndsWith("-p2")){
		//LHC12x-p2
        //useDataRunListForMC=kTRUE;
		addSubFolder=kTRUE;
        doEquidistantXaxis=kTRUE;
		cutNr=0;
		fEnergyFlag="8TeV";
        filePath = "DataQA/20160125";
        filePathPhoton = "DataQA/20160104";
		fileName = "GammaConvCalo_120.root";
		select.Remove(select.Length()-3,select.Length());
        nSets = 3;
        nData = 1;
        TString tmpPeriod = select;
        TString tmpA = "";
        tmpPeriod.Remove(0,tmpPeriod.Length()-1);
        if(tmpPeriod.CompareTo("a")==0) tmpA="1";
        DataSets[0]=select;DataSets[1]=Form("LHC15h1%s%s",tmpPeriod.Data(),tmpA.Data()); DataSets[2]=Form("LHC15h2%s",tmpPeriod.Data());
		plotDataSets[0]=select;plotDataSets[1]="Pythia8";plotDataSets[2]="Phojet";
    }
//**************************************************************************************************************
	else if(select.BeginsWith("LHC12") && select.EndsWith("-Data") && select.Length()==11){
		//LHC12x - only Data
		select.Remove(select.Length()-5,5);
        addSubFolder=kTRUE;
		doEquidistantXaxis=kTRUE;
        cutNr=0;
		fEnergyFlag="8TeV";
        filePath = "DataQA/20160125";
        filePathPhoton = "DataQA/20160104";
		fileName = "GammaConvCalo_120.root";
		nSets = 1;
		nData = 1;
		DataSets[0]=select;
		plotDataSets[0]=select;
	}
//**************************************************************************************************************
	else if(select.CompareTo("LHC12-MC")==0){
		//LHC12 MC
		fEnergyFlag="8TeV";
        filePath = "DataQA/20160125";
        filePathPhoton = "DataQA/20160104";
		fileName = "GammaConvCalo_120.root";
		nSets = 2;
		nData = 0;
        DataSets[0]="LHC15h1";DataSets[1]="LHC15h2";
		plotDataSets[0]="Pythia8";plotDataSets[1]="Phojet";
	}
//**************************************************************************************************************
	else if(select.BeginsWith("LHC12-kEMC7")){
		//LHC12 Trigger
		doHistsForEverySet = kFALSE;
		cutNr=0;
		fEnergyFlag="8TeV";
        filePath = "DataQA/20160125";
		fileName = "GammaConvCalo_121.root";
		if(select.EndsWith("-single"))
		{
			nSets = 1;
			nData = 1;
			DataSets[0]="LHC12-kEMC7";
			plotDataSets[0]="LHC12";
		}else{
            nSets = 6;
            nData = 6;
            DataSets[0]="LHC12b-kEMC7";DataSets[1]="LHC12c-kEMC7";DataSets[2]="LHC12d-kEMC7";
            DataSets[3]="LHC12f-kEMC7";/*DataSets[4]="LHC12g-kEMC7";*/DataSets[4]="LHC12h-kEMC7";DataSets[5]="LHC12i-kEMC7";
            plotDataSets[0]="LHC12b";plotDataSets[1]="LHC12c";plotDataSets[2]="LHC12d";
            plotDataSets[3]="LHC12f";/*plotDataSets[4]="LHC12g";*/plotDataSets[4]="LHC12h";plotDataSets[5]="LHC12i";
		}
	}
//**************************************************************************************************************
	else if(select.BeginsWith("LHC12-kEMCEGA")){
		//LHC12 Trigger
		doHistsForEverySet = kFALSE;
		cutNr=0;
		fEnergyFlag="8TeV";
        filePath = "DataQA/20160125";
		fileName = "GammaConvCalo_122.root";
		if(select.EndsWith("-single"))
		{
			nSets = 1;
			nData = 1;
			DataSets[0]="LHC12-kEMCEGA";
			plotDataSets[0]="LHC12";
		}else{
            nSets = 5;
            nData = 5;
            DataSets[0]="LHC12c-kEMCEGA";DataSets[1]="LHC12d-kEMCEGA";DataSets[2]="LHC12f-kEMCEGA";
            /*DataSets[3]="LHC12g-kEMCEGA";*/DataSets[3]="LHC12h-kEMCEGA";DataSets[4]="LHC12i-kEMCEGA";
            plotDataSets[0]="LHC12c";plotDataSets[1]="LHC12d";plotDataSets[2]="LHC12f";
            /*plotDataSets[3]="LHC12g";*/plotDataSets[3]="LHC12h";plotDataSets[4]="LHC12i";
		}
	}
//**************************************************************************************************************
	else if(!select.BeginsWith("LHC12-kEMC") && select.BeginsWith("LHC12") && select.Contains("-kEMC")){
		//LHC12x - Trigger
		TString number = "";
		if(select.EndsWith("-kEMC7")) number = "121";
		if(select.EndsWith("-kEMCEGA")) number = "122";
		doEquidistantXaxis=kTRUE;
        addSubFolder=kTRUE;
		cutNr=0;
		fEnergyFlag="8TeV";
        filePath = "DataQA/20160125";
		fileName = Form("GammaConvCalo_%s.root",number.Data());
		nSets = 1;
		nData = 1;
		DataSets[0]=select;
		TString selPlot = select;
		selPlot.Resize(6);
		plotDataSets[0]=selPlot;
	}
//**************************************************************************************************************
    else if(select.CompareTo("LHC12-kEMC8")==0){
      //LHC12 Trigger
      doHistsForEverySet = kFALSE;
      cutNr=0;
      fEnergyFlag="8TeV";
      filePath = "DataQA/20160203";
      fileName = "GammaConvCalo_131.root";
      nSets = 9;
      nData = 9;
      DataSets[0]="LHC12a-kEMC8";DataSets[1]="LHC12b-kEMC8";DataSets[2]="LHC12c-kEMC8";DataSets[3]="LHC12d-kEMC8";DataSets[4]="LHC12e-kEMC8";
      DataSets[5]="LHC12f-kEMC8";DataSets[6]="LHC12g-kEMC8";DataSets[7]="LHC12h-kEMC8";DataSets[8]="LHC12i-kEMC8";
      plotDataSets[0]="LHC12a";plotDataSets[1]="LHC12b";plotDataSets[2]="LHC12c";plotDataSets[3]="LHC12d";plotDataSets[4]="LHC12e";
      plotDataSets[5]="LHC12f";plotDataSets[6]="LHC12g";plotDataSets[7]="LHC12h";plotDataSets[8]="LHC12i";
    }
//**************************************************************************************************************
    else if(select.CompareTo("LHC12-kEMC8EGA")==0){
      //LHC12 Trigger
      doHistsForEverySet = kFALSE;
      cutNr=0;
      fEnergyFlag="8TeV";
      filePath = "DataQA/20160203";
      fileName = "GammaConvCalo_132.root";
      nSets = 9;
      nData = 9;
      DataSets[0]="LHC12a-kEMC8EGA";DataSets[1]="LHC12b-kEMC8EGA";DataSets[2]="LHC12c-kEMC8EGA";DataSets[3]="LHC12d-kEMC8EGA";DataSets[4]="LHC12e-kEMC8EGA";
      DataSets[5]="LHC12f-kEMC8EGA";DataSets[6]="LHC12g-kEMC8EGA";DataSets[7]="LHC12h-kEMC8EGA";DataSets[8]="LHC12i-kEMC8EGA";
      plotDataSets[0]="LHC12a";plotDataSets[1]="LHC12b";plotDataSets[2]="LHC12c";plotDataSets[3]="LHC12d";plotDataSets[4]="LHC12e";
      plotDataSets[5]="LHC12f";plotDataSets[6]="LHC12g";plotDataSets[7]="LHC12h";plotDataSets[8]="LHC12i";
    }
    //**************************************************************************************************************
    else if(select.CompareTo("LHC13pPb-Data")==0){
        //LHC12
        doHistsForEverySet = kTRUE;
        fEnergyFlag="pPb_5.023TeV";
        filePath = "DataQA/20160307";
        filePathPhoton = "DataQA/20160307";
        fileName = "GammaConvCalo_20.root";
        cutNr=0;
        nSets = 5;
        nData = 5;
        DataSets[0]="LHC13b";DataSets[1]="LHC13c";DataSets[2]="LHC13d";DataSets[3]="LHC13e";DataSets[4]="LHC13f";
        plotDataSets[0]="LHC13b";plotDataSets[1]="LHC13c";plotDataSets[2]="LHC13d";plotDataSets[3]="LHC13e";plotDataSets[4]="LHC13f";
        
    }
    //**************************************************************************************************************
    else if(select.CompareTo("LHC13pPb-MBData")==0){
        //LHC12
        doEquidistantXaxis=kTRUE;
        doHistsForEverySet = kFALSE;
        fEnergyFlag="pPb_5.023TeV";
        filePath = "DataQA/20160307";
        filePathPhoton = "DataQA/20160307";
        fileName = "GammaConvCalo_20.root";
        cutNr=0;
        nSets = 2;
        nData = 2;
        DataSets[0]="LHC13b";DataSets[1]="LHC13c";
        plotDataSets[0]="LHC13b";plotDataSets[1]="LHC13c";
    }
    //**************************************************************************************************************
    else if(select.CompareTo("LHC13pPb-MB")==0){
        //pPb MB
        doEquidistantXaxis=kTRUE;
        doHistsForEverySet = kFALSE;
        fEnergyFlag="pPb_5.023TeV";
        filePath = "DataQA/20160307";
        filePathPhoton = "DataQA/20160307";
        fileName = "GammaConvCalo_20.root";
        cutNr=0;
        nSets = 4;
        nData = 2;
        DataSets[0]="LHC13b";DataSets[1]="LHC13c";DataSets[2]="LHC13b2_efix";DataSets[3]="LHC13e7";
        plotDataSets[0]="LHC13b";plotDataSets[1]="LHC13c";plotDataSets[2]="DPMJET";plotDataSets[3]="HIJING";
    }

    //**************************************************************************************************************
	else{
		cout << "No valid selection! Returning..." << endl;
		return;
	}
//**************************************************************************************************************
	if(doEventQA) EventQA_Runwise(nSets,nData,fEnergyFlag,filePath,fileName,DataSets,plotDataSets,mode,cutNr,
                                  doExtQA,doEquidistantXaxis,doTrigger,doHistsForEverySet,addSubFolder,useDataRunListForMC,markerSize,suffix,folderRunlists);
	if(doPhotonQA){
		TString path = filePath;
		if(!filePathPhoton.IsNull()) path = filePathPhoton;
		PhotonQA_Runwise(nSets,nData,fEnergyFlag,path,fileName,DataSets,plotDataSets,mode,cutNr,
                                    doExtQA,doEquidistantXaxis,doTrigger,doHistsForEverySet,addSubFolder,useDataRunListForMC,markerSize,suffix,folderRunlists);
	}
	if(doClusterQA) ClusterQA_Runwise(nSets,nData,fEnergyFlag,filePath,fileName,DataSets,plotDataSets,mode,cutNr,
                                      doExtQA,doEquidistantXaxis,doTrigger,doHistsForEverySet,addSubFolder,useDataRunListForMC,markerSize,suffix,folderRunlists);
	return;
}
