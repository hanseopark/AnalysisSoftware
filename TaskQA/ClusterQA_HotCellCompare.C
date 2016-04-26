#include "QA.h"

void ClusterQA_HotCellCompare(TString suffix = "eps", Int_t mode = 2)
{
	cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
	cout << "ClusterQA_HotCellCompare" << endl;
	cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

	gROOT->Reset();

	StyleSettingsThesis();
	SetPlotStyle();

	TString outputDir ="ClusterQA_HotCellCompare";

    const Int_t nSets = 5;
    const Int_t nTrigger = 1;
    TString fEnergyFlag = "7TeV";
    outputDir+="/LHC10";
    TString DataSets[nSets]={"LHC10b_pass4","LHC10c_pass4","LHC10d_pass4","LHC10e_pass4","LHC10f_pass4"};
    TString Triggers[nTrigger]={""};
    TString cut[nTrigger] = {"00000113_00200009327000008250400000_1111100013032230000_0163103100000010"};


//	const Int_t nSets = 8;
//	const Int_t nTrigger = 1;
//	TString fEnergyFlag = "8TeV";
//	outputDir+="/LHC12";
//	TString DataSets[nSets]={"LHC12a","LHC12b","LHC12c","LHC12d","LHC12f","LHC12g","LHC12h","LHC12i"};
//	TString Triggers[nTrigger]={""};
//	TString cut[nTrigger] = {"0000011_00200009327000008250400000_11111063032230000_0163103100000000"};

//	const Int_t nSets = 9;
//	const Int_t nTrigger = 4;
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
//	TString fEnergyFlag = "2.76TeV";
//	TString DataSets[nSets]={"LHC11a_p4_wSDD"};
//	TString Triggers[nTrigger]={"","-kEMC1"};
//	TString cut[nTrigger] = {"0000311_00200009327000008250400000_11111053032230000_0163103100000000",
//							 "0005111_00200009327000008250400000_11111053032230000_0163103100000000"
//							};

//	const Int_t nSets = 1;
//	const Int_t nTrigger = 4;
//	TString fEnergyFlag = "2.76TeV";
//	TString DataSets[nSets]={"LHC13g"};
//	TString Triggers[nTrigger]={"","-kEMC7","-kEMCEG1","-kEMCEG2"};
//	TString cut[nTrigger] = {"0000011_00200009327000008250400000_11111063032230000_0163103100000000",
//							 "0005211_00200009327000008250400000_11111063032230000_0163103100000000",
//							 "0008311_00200009327000008250400000_11111063032230000_0163103100000000",
//							 "0008511_00200009327000008250400000_11111063032230000_0163103100000000"
//							};

	gSystem->Exec("mkdir -p "+outputDir);

//	std::vector<Int_t> globalRuns;
//	std::map<Int_t,Int_t> runToNumber;
//	cout << "readin";
//	for(Int_t i=0; i<nTrigger; i++)
//	{
//		TString fileRuns = Form("runNumbersLHC12%s.txt", Triggers[i].Data());
//		if(!readin(fileRuns, globalRuns, kFALSE, kFALSE)) {cout << "ERROR, no Run Numbers could be found! Returning..." << endl; return;}
//	}
//	cout << "done!" << endl;
//	cout << "sort";
//	selection_sort(globalRuns.begin(), globalRuns.end());
//	cout << "done!" << endl;
//	vector<Int_t>::iterator it;
//	cout << "unique";
//	it = unique(globalRuns.begin(),globalRuns.end());
//	globalRuns.resize( distance(globalRuns.begin(),it) );
//	cout << "done!" << endl;

//	for(Int_t i=0; i<(Int_t)globalRuns.size(); i++){
//		runToNumber[globalRuns.at(i)]=i;
//	}

	fstream fLogRunwiseBadCells;
	fstream fLogOutput;
	fstream fLogOutputDetailed;
	std::map <Int_t,Int_t> ma;
	std::map <Int_t,Float_t> maNFired;
	std::map <Int_t,Float_t> maNTimes;
	std::vector<TString> uniqueCellIDRun;
	std::vector<TString>::iterator it;
	std::vector<Int_t> vec;
	TString currentRun="";

	for(Int_t j=0; j<nSets; j++)
	{
		cout << DataSets[j].Data() << endl;
		fLogOutput.open(Form("%s/%s.log",outputDir.Data(),DataSets[j].Data()),ios::out);
		fLogOutputDetailed.open(Form("%s/%s-Detailed.log",outputDir.Data(),DataSets[j].Data()),ios::out);
		for(Int_t i=0; i<nTrigger; i++)
		{
            fLogRunwiseBadCells.open(Form("%s/%s/ClusterQA/%s/Runwise/%s/HotCellsRunwise-%s%s.log",cut[i].Data(),fEnergyFlag.Data(),suffix.Data(),DataSets[j].Data(),DataSets[j].Data(),Triggers[i].Data()), ios::in);
			if(fLogRunwiseBadCells.good())
			{
				fLogRunwiseBadCells.seekg(0L, ios::beg);
				TString fVar;
				while(!fLogRunwiseBadCells.eof())
				{
					fLogRunwiseBadCells >> fVar;
					if(fVar.BeginsWith("Run-")) {
						TString curr = fVar;
						currentRun = curr.Remove(0,4);
					}
					else if(fVar.BeginsWith("NoNoisy")||fVar.BeginsWith("NotEnough")){
						continue;
					}
					else if(fVar.Sizeof()>1) {
						TObjArray *rNumber = fVar.Tokenize("-");
						TObjString* rString = (TObjString*)rNumber->At(0);
						TObjString* rStringNFired = (TObjString*)rNumber->At(1);
						TObjString* rStringNTimes = (TObjString*)rNumber->At(2);
						TString vecString = rString->GetString();
						TString vecStringNFired = rStringNFired->GetString();
						TString vecStringNTimes = rStringNTimes->GetString();
						vecStringNFired.Remove(0,14);
						vecStringNTimes.Remove(0,12);
						TString currentRunVecString = vecString;
						currentRunVecString+="-";
						currentRunVecString+=currentRun;
						it = find (uniqueCellIDRun.begin(), uniqueCellIDRun.end(), currentRunVecString);
						if( it == uniqueCellIDRun.end() ){
							vec.push_back(vecString.Atoi());
							uniqueCellIDRun.push_back(currentRunVecString);
							if( ma.find(vecString.Atoi()) != ma.end() ){
								if(currentRun.Atoi()<ma[vecString.Atoi()]) ma[vecString.Atoi()]=currentRun.Atoi();
							}
							else{
								ma[vecString.Atoi()] = currentRun.Atoi();
							}
						}
						if( maNFired.find(vecString.Atoi()) != maNFired.end() ){
							maNFired[vecString.Atoi()] += vecStringNFired.Atof();
							maNTimes[vecString.Atoi()] += vecStringNTimes.Atof();
						}
						else{
							maNFired[vecString.Atoi()] = vecStringNFired.Atof();
							maNTimes[vecString.Atoi()] = vecStringNTimes.Atof();
						}
					}
				}
			}else{
				cout << "\tCould not open HotCellsRunwise file for: " << DataSets[j].Data() << Triggers[i].Data() << endl;
			}

//			if((Int_t)ma[j].size()==0){
//				cout << "No Bad Cells found for: " << DataSets[j].Data() << Triggers[i].Data() << endl;
//			}
			fLogRunwiseBadCells.close();
		}

		selection_sort(vec.begin(), vec.end());

		Int_t beforeCellID=vec.at(0);
		Int_t currentCellID=0;
		Int_t nFired=1;
		Int_t threshNFired = 1;
		Int_t threshNTotalFired = 50;
		if(DataSets[j].CompareTo("LHC13g")==0) threshNFired=0;
		if(DataSets[j].CompareTo("LHC13g")==0) threshNTotalFired=10;
		for(Int_t iVec=1; iVec<(Int_t)vec.size(); iVec++)
		{
			//cout << "iVec: " << iVec << endl;
			currentCellID=vec.at(iVec);
			if(beforeCellID!=currentCellID){
				//cout << beforeCellID << " " << nFired << " " << ma[beforeCellID] << endl;
				if( (nFired>threshNFired && maNFired[beforeCellID]>threshNTotalFired) || (maNFired[beforeCellID]>threshNTotalFired && (maNFired[beforeCellID]/maNTimes[beforeCellID]>5)) ) {
					fLogOutput << beforeCellID << " " << ma[beforeCellID] << endl;
					fLogOutputDetailed << "CellID: " << beforeCellID << ", occurred for first time in run: " << ma[beforeCellID] << ", in total noisy in - " << nFired << " - differentRuns, where it occurred in clusters " << maNFired[beforeCellID] << " times! That is " << ((Float_t)maNFired[beforeCellID])/maNTimes[beforeCellID] << " times than the average of neighboring cells!" << endl;
				}else{
					cout << "discard " << beforeCellID << ", occurred " << nFired << " times and total EFrac: " << maNFired[beforeCellID] << " which is " << (Float_t)maNFired[beforeCellID]/maNTimes[beforeCellID] << " times than the average of neighboring cells!" << endl;
				}
				beforeCellID=currentCellID;
				nFired=1;
			}else nFired++;
		}
		//cout << beforeCellID << " " << nFired << " " << ma[beforeCellID] << endl;
		if(nFired>threshNFired && maNFired[beforeCellID]>threshNTotalFired) {
			fLogOutput << beforeCellID << " " << ma[beforeCellID] << endl;
			fLogOutputDetailed << "CellID: " << beforeCellID << ", occurred for first time in run: " << ma[beforeCellID] << ", in total noisy in - " << nFired << " - differentRuns, where it occurred in clusters " << maNFired[beforeCellID] << " times! That is " << ((Float_t)maNFired[beforeCellID])/maNTimes[beforeCellID] << " times than the average of neighboring cells!" << endl;
		}

		gSystem->Exec(Form("touch %s/%s-SortedByRun.log",outputDir.Data(),DataSets[j].Data()));
		gSystem->Exec(Form("sort -k 2 %s/%s.log > %s/%s-SortedByRun.log",outputDir.Data(),DataSets[j].Data(),outputDir.Data(),DataSets[j].Data()));

		vec.clear();
		ma.clear();
		uniqueCellIDRun.clear();
		fLogOutput.close();
		fLogOutputDetailed.close();
	}

	cout << "Done with ClusterQA_HotCellCompare" << endl;
	cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
	return;
}

