#include "QA.h"

void ClusterQA_DeadCellCompare(TString suffix = "eps", Int_t mode = 2)
{
	cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
	cout << "ClusterQA_DeadCellCompare" << endl;
	cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

	gROOT->Reset();

	StyleSettingsThesis();
	SetPlotStyle();

//	const Int_t nCaloCells = 11520;
	TString outputDir ="ClusterQA_DeadCellCompare";

// LHC10, with MC comparison
    const Int_t nCaloCells = 4608;
    const Int_t nSets = 5;
    const Int_t nTrigger = 1;
    const Int_t nMCSets = 5;
    outputDir+="/LHC10";
    TString fEnergyFlag = "7TeV";
    TString DataSets[nSets]={"LHC10b_pass4","LHC10c_pass4","LHC10d_pass4","LHC10e_pass4","LHC10f_pass4"};
    TString MCSets[nMCSets]={"LHC14j4b","LHC14j4c","LHC14j4d","LHC14j4e","LHC14j4f"};
    TString MCcut[nMCSets] = {"00000113_00200009327000008250400000_1111100013032230000_0163103100000010",
                              "00000113_00200009327000008250400000_1111100013032230000_0163103100000010",
                              "00000113_00200009327000008250400000_1111100013032230000_0163103100000010",
                              "00000113_00200009327000008250400000_1111100013032230000_0163103100000010",
                              "00000113_00200009327000008250400000_1111100013032230000_0163103100000010"};
    TString Triggers[nTrigger]={""};
    TString cut[nTrigger] = {"00000113_00200009327000008250400000_1111100013032230000_0163103100000010"};


// LHC12 only MB, with MC comparison
//	const Int_t nSets = 8;
//	const Int_t nTrigger = 1;
//	const Int_t nMCSets = 2;
//	outputDir+="/LHC12";
//	TString fEnergyFlag = "8TeV";
//	TString DataSets[nSets]={"LHC12a","LHC12b","LHC12c","LHC12d","LHC12f","LHC12g","LHC12h","LHC12i"};
//	TString MCSets[nMCSets]={"LHC14e2a","LHC14e2c"};
//	TString MCcut[nMCSets] = {"0000011_00200009327000008250400000_11111063032230000_0163103100000000",
//							  "0000011_00200009327000008250400000_11111063032230000_0163103100000000"};
//	TString Triggers[nTrigger]={""};
//	TString cut[nTrigger] = {"0000011_00200009327000008250400000_11111063032230000_0163103100000000"};

//	const Int_t nSets = 9;
//	const Int_t nTrigger = 4;
//	const Int_t nMCSets = 0;
//	TString fEnergyFlag = "8TeV";
//	TString DataSets[nSets]={"LHC12a","LHC12b","LHC12c","LHC12d","LHC12e","LHC12f","LHC12g","LHC12h","LHC12i"};
//	TString Triggers[nTrigger]={"","-kEMC7","-kEMCEGA","-kEMCEJE"};
//	TString cut[nTrigger] = {"0000011_00200009327000008250400000_10000063032230000_0163103100000000",
//							 "0005211_00200009327000008250400000_10000063032230000_0163103100000000",
//							 "0008111_00200009327000008250400000_10000063032230000_0163103100000000",
//							 "0009111_00200009327000008250400000_10000063032230000_0163103100000000"
//							};

//	const Int_t nSets = 1;
//	const Int_t nTrigger = 2;
//	const Int_t nMCSets = 0;
//	TString fEnergyFlag = "2.76TeV";
//	TString DataSets[nSets]={"LHC11a_p4_wSDD"};
//	TString Triggers[nTrigger]={"","-kEMC1"};
//	TString cut[nTrigger] = {"0000311_00200009327000008250400000_11111053032230000_0163103100000000",
//							 "0005111_00200009327000008250400000_11111053032230000_0163103100000000"
//							};

//	const Int_t nSets = 1;
//	const Int_t nTrigger = 4;
//	const Int_t nMCSets = 0;
//	TString fEnergyFlag = "2.76TeV";
//	TString DataSets[nSets]={"LHC13g"};
//	TString Triggers[nTrigger]={"","-kEMC7","-kEMCEG1","-kEMCEG2"};
//	TString cut[nTrigger] = {"0000011_00200009327000008250400000_11111063032230000_0163103100000000",
//							 "0005211_00200009327000008250400000_11111063032230000_0163103100000000",
//							 "0008311_00200009327000008250400000_11111063032230000_0163103100000000",
//							 "0008511_00200009327000008250400000_11111063032230000_0163103100000000"
//							};

	gSystem->Exec("mkdir -p "+outputDir);

	fstream fLogRunwiseDeadCells;
	fstream fLogOutput;
	fstream fLogOutputFinal;
	fstream fLogOutputDetailed;
	std::vector<Int_t>::iterator it;
	std::vector<Int_t> vec[nCaloCells];
	std::vector<Int_t> vecMC[nCaloCells];
	std::vector<Int_t> vecMCNotEnoughNoCold;
	std::vector<Int_t> globalRuns;
	std::map<Int_t,Int_t> mapRuns;
	std::map<Int_t,Int_t> mapRunsReverse;
	std::map<TString,Int_t> mapRunRanges;
	std::vector<Int_t> vecRuns;
	std::vector<TString> vecWhichRuns;
	std::vector<Int_t> vecNCons;
	std::vector<Int_t> vecNCellRuns;
	TString currentRun="";
	TString fileRuns;
	std::vector<TString> runs;

	for(Int_t j=0; j<nSets; j++)
	{
		cout << DataSets[j].Data() << endl;

		//readin of MC cold/dead cells
		if(nMCSets>0){
          runs.clear();
          fileRuns = Form("runNumbers%s.txt",MCSets[j].Data());
          if(!readin(fileRuns, runs, kFALSE)) {cout << Form("INFO, no Run Numbers could be found for %s!",MCSets[j].Data()) << endl;}

          fLogRunwiseDeadCells.open(Form("%s/%s/ClusterQA/%s/Runwise/%s/ColdCellsRunwise-%s.log",MCcut[j].Data(),fEnergyFlag.Data(),suffix.Data(),DataSets[j].Data(),MCSets[j].Data()), ios::in);
          if(fLogRunwiseDeadCells.good())
          {
            fLogRunwiseDeadCells.seekg(0L, ios::beg);
            TString fVar;
            while(!fLogRunwiseDeadCells.eof())
            {
              fLogRunwiseDeadCells >> fVar;
              if(fVar.BeginsWith("Run-")) {
                TString curr = fVar;
                currentRun = curr.Remove(0,4);
              }
              else if(fVar.BeginsWith("NoCold")||fVar.BeginsWith("NotEnough")){
                it = find (vecMCNotEnoughNoCold.begin(), vecMCNotEnoughNoCold.end(), currentRun.Atoi());
                if( it == vecMCNotEnoughNoCold.end() ) vecMCNotEnoughNoCold.push_back(currentRun.Atoi());
                continue;
              }
              else if(fVar.Sizeof()>1) {
                TObjArray *rNumber = fVar.Tokenize("-");
                TObjString* rString = (TObjString*)rNumber->At(0);
                TString vecString = rString->GetString();
                it = find (vecMC[vecString.Atoi()].begin(), vecMC[vecString.Atoi()].end(), currentRun.Atoi());
                if( it == vecMC[vecString.Atoi()].end() ) vecMC[vecString.Atoi()].push_back(currentRun.Atoi());
              }
            }
          }else{
            cout << "\tCould not open DeadCellsRunwise file for: " << MCSets[j].Data() << endl;
          }
          fLogRunwiseDeadCells.close();
		}

		fLogOutput.open(Form("%s/%s.log",outputDir.Data(),DataSets[j].Data()),ios::out);
		fLogOutputDetailed.open(Form("%s/%s-Detailed.log",outputDir.Data(),DataSets[j].Data()),ios::out);
		for(Int_t i=0; i<nTrigger; i++)
		{
//			cout << i << " - '" << Triggers[i].Data() << "'" << endl;
			runs.clear();
			fileRuns = Form("runNumbers%s%s.txt", DataSets[j].Data(),Triggers[i].Data());
			if(!readin(fileRuns, runs, kFALSE)) {cout << Form("INFO, no Run Numbers could be found for %s%s!",DataSets[j].Data(),Triggers[i].Data()) << endl;}

            fLogRunwiseDeadCells.open(Form("%s/%s/ClusterQA/%s/Runwise/%s/ColdCellsRunwise-%s%s.log",cut[i].Data(),fEnergyFlag.Data(),suffix.Data(),DataSets[j].Data(),DataSets[j].Data(),Triggers[i].Data()), ios::in);
			if(fLogRunwiseDeadCells.good())
			{
				fLogRunwiseDeadCells.seekg(0L, ios::beg);
				TString fVar;
				while(!fLogRunwiseDeadCells.eof())
				{
					fLogRunwiseDeadCells >> fVar;
					if(fVar.BeginsWith("Run-")) {
						TString curr = fVar;
						currentRun = curr.Remove(0,4);
					}
					else if(fVar.BeginsWith("NoCold")||fVar.BeginsWith("NotEnough")){
						continue;
					}
					else if(fVar.Sizeof()>1) {
						TObjArray *rNumber = fVar.Tokenize("-");
						TObjString* rString = (TObjString*)rNumber->At(0);
						TString vecString = rString->GetString();
						if(nMCSets>0){
							it = find (vecMC[vecString.Atoi()].begin(), vecMC[vecString.Atoi()].end(), currentRun.Atoi());
							if( !(it == vecMC[vecString.Atoi()].end()) ) continue;

							it = find (vecMCNotEnoughNoCold.begin(), vecMCNotEnoughNoCold.end(), currentRun.Atoi());
							if( !(it == vecMCNotEnoughNoCold.end()) ) continue;
						}
						it = find (vec[vecString.Atoi()].begin(), vec[vecString.Atoi()].end(), currentRun.Atoi());
						if( it == vec[vecString.Atoi()].end() ) vec[vecString.Atoi()].push_back(currentRun.Atoi());

						it = find (globalRuns.begin(), globalRuns.end(), currentRun.Atoi());
						if( it == globalRuns.end() ) globalRuns.push_back(currentRun.Atoi());
					}
				}
			}else{
				cout << "\tCould not open DeadCellsRunwise file for: " << DataSets[j].Data() << Triggers[i].Data() << endl;
			}
			fLogRunwiseDeadCells.close();
//			TString fFileRunwise = Form("%s/%s/ClusterQA/%s%s_ClusterQARunwise_pp.root", cut[i].Data(),fEnergyFlag.Data(),DataSets[j].Data(), Triggers[i].Data());
//			cout << "Opening file: " << fFileRunwise.Data() << endl;
//			TFile* FileRunwise = new TFile(fFileRunwise.Data(),"READ");
//			if(FileRunwise->IsZombie()) {cout << "INFO: ROOT file '" << fFileRunwise.Data() << "' could not be openend, continue..." << endl; continue;}
//			FileRunwise->cd("MissingCells");
//			for(Int_t iRun=0; iRun<(Int_t)runs.size(); iRun++)
//			{
//				Int_t currentRun=((TString)runs.at(iRun)).Atoi();
//				TH2D* temp;
//				gDirectory->GetObject(Form("%i-LHC14e2a", currentRun),temp);
//				cout << iRun << " " << runs.size() << " " << Form("%i-LHC14e2a", currentRun) << " " << temp << endl;
//				if(temp){
//					for(Int_t iCell=6; iCell<(4*nCaloCells)+1; iCell+=4)
//					{
//						if(temp->GetBinContent(iCell)>0){
//							it = find (vec[((iCell-6)/4)].begin(), vec[((iCell-6)/4)].end(), currentRun);
//							if( it == vec[((iCell-6)/4)].end() ) vec[((iCell-6)/4)].push_back(currentRun);
//							it = find (globalRuns.begin(), globalRuns.end(), currentRun);
//							if( it == globalRuns.end() ) globalRuns.push_back(currentRun);
//						}
//					}
//				}else cout << "ERROR: could not open: " << Form("%i-LHC14e2a",currentRun) << endl;
//			}
//			FileRunwise->Close();
//			delete FileRunwise;
		}

		selection_sort(globalRuns.begin(),globalRuns.end());
		for(Int_t iRun=0; iRun<(Int_t)globalRuns.size(); iRun++) {
			mapRuns[iRun]=globalRuns.at(iRun);
			mapRunsReverse[globalRuns.at(iRun)]=iRun;
		}

		for(Int_t iCell=0;iCell<nCaloCells;iCell++){
			selection_sort(vec[iCell].begin(), vec[iCell].end());
			//if(vec[iCell].size()>0)cout << ", " << iCell << " - " << vec[iCell].size() << endl;
			if(vec[iCell].size()>0){
				fLogOutputDetailed << "CellID:" << iCell << ",cold/notFiredInRuns:";
				for(Int_t j=0; j<(Int_t)vec[iCell].size(); j++){fLogOutputDetailed << vec[iCell].at(j) << ",";}
				fLogOutputDetailed << endl;
			}
			vecNCellRuns.push_back((Int_t)vec[iCell].size());
			if(vec[iCell].size()>1){
				Int_t nRun=1;
				Int_t nBegin=0;
				Int_t nEnd=1;
				Int_t nConsecutiveRuns=1;
				do{
					if((mapRunsReverse[vec[iCell].at(nEnd)]-mapRunsReverse[vec[iCell].at(nEnd-1)]) == 1) {
						nEnd++;
						nConsecutiveRuns++;
					}
					else{
						if(nConsecutiveRuns>1){
							fLogOutput << iCell << " " << vec[iCell].at(nBegin) << " " << vec[iCell].at(nEnd-1) << endl;
							vecRuns.push_back(iCell);
							vecWhichRuns.push_back(Form("%i %i",vec[iCell].at(nBegin),vec[iCell].at(nEnd-1)));
							vecNCons.push_back(nConsecutiveRuns);
							if(mapRunRanges.find(Form("%i %i",vec[iCell].at(nBegin),vec[iCell].at(nEnd-1)))==mapRunRanges.end()) mapRunRanges[Form("%i %i",vec[iCell].at(nBegin),vec[iCell].at(nEnd-1))]=1;
							else mapRunRanges[Form("%i %i",vec[iCell].at(nBegin),vec[iCell].at(nEnd-1))]++;
						}
						nBegin=nEnd;
						nEnd+=1;
						nConsecutiveRuns=1;
					}
					nRun++;
					//cout << nRun << " " << vec[iCell].size() << " " << nBegin << " " << nEnd << " " << nConsecutiveRuns << " " << endl;
				}while(nRun+1<(Int_t)vec[iCell].size());
				if(vec[iCell].size()==2) nEnd--;
				if(nConsecutiveRuns>1){
					fLogOutput << iCell << " " << vec[iCell].at(nBegin) << " " << vec[iCell].at(nEnd) << endl;
					vecRuns.push_back(iCell);
					vecWhichRuns.push_back(Form("%i %i",vec[iCell].at(nBegin),vec[iCell].at(nEnd)));
					vecNCons.push_back(nConsecutiveRuns);
					if(mapRunRanges.find(Form("%i %i",vec[iCell].at(nBegin),vec[iCell].at(nEnd)))==mapRunRanges.end()) mapRunRanges[Form("%i %i",vec[iCell].at(nBegin),vec[iCell].at(nEnd))]=1;
					else mapRunRanges[Form("%i %i",vec[iCell].at(nBegin),vec[iCell].at(nEnd))]++;
				}
			}
		}

		fLogOutputFinal.open(Form("%s/%s-Final.log",outputDir.Data(),DataSets[j].Data()),ios::out);
		Int_t sortedOutFinal=0;
		Double_t fractionThesh = 0.5;
		if(DataSets[j].CompareTo("LHC12h")==0) fractionThesh = 0.8;
        if(DataSets[j].Contains("LHC10") && DataSets[j].Contains("_pass4") && DataSets[j].Length()==12) fractionThesh = 0.6;
		std::vector<Int_t> checkDoubleCellsThreshold;
		for(Int_t vecCell=0;vecCell<(Int_t)vecRuns.size();vecCell++){
//			cout << vecNCellRuns.at(vecRuns.at(vecCell)) << " " << globalRuns.size() << endl;
//			cout << mapRunRanges[vecWhichRuns.at(vecCell)] << " " << vecNCons.at(vecCell) << " " << (((Float_t)vecNCellRuns.at(vecRuns.at(vecCell))/(Float_t)globalRuns.size())) << endl;
			if( (((Float_t)vecNCellRuns.at(vecRuns.at(vecCell))/(Float_t)globalRuns.size()) > fractionThesh) ){
				if(find (checkDoubleCellsThreshold.begin(), checkDoubleCellsThreshold.end(), vecRuns.at(vecCell)) == checkDoubleCellsThreshold.end()){
					checkDoubleCellsThreshold.push_back(vecRuns.at(vecCell));
					fLogOutputFinal << vecRuns.at(vecCell) << " in " << ((Float_t)vecNCellRuns.at(vecRuns.at(vecCell))/(Float_t)globalRuns.size())*100 << "% of selected runs dead/cold" << endl;
				}
			}else if(mapRunRanges[vecWhichRuns.at(vecCell)]>10 || vecNCons.at(vecCell)>4){
				fLogOutputFinal << vecRuns.at(vecCell) << " " << vecWhichRuns.at(vecCell) << endl;
			}else{
				sortedOutFinal++;
//				cout << vecCell << ", " << mapRunRanges[vecWhichRuns.at(vecCell)] << ", " << vecNCons.at(vecCell) << ", " << (((Float_t)vecNCellRuns.at(vecRuns.at(vecCell))/(Float_t)runs.size())) << endl;
			}
		}
		cout << "\t discarded " << sortedOutFinal << " deadcell ranges while generating final list." << endl;
		fLogOutputFinal.close();

		checkDoubleCellsThreshold.clear();
		for(Int_t iCell=0;iCell<nCaloCells;iCell++){vec[iCell].clear(); vecMC[iCell].clear();}
		vecMCNotEnoughNoCold.clear();
		globalRuns.clear();
		mapRuns.clear();
		mapRunsReverse.clear();
		mapRunRanges.clear();
		vecRuns.clear();
		vecWhichRuns.clear();
		vecNCons.clear();
		vecNCellRuns.clear();
		fLogOutput.close();
		fLogOutputDetailed.close();
	}

	cout << "Done with ClusterQA_DeadCellCompare" << endl;
	cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
	return;
}

