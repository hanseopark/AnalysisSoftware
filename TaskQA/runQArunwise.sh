#/bin/bash

# $1=outDIR;

fileName0010_Data=PhotonQA23March_0005314141_5010001_02200009247600003190000000_Data.root;
cutNumber0010_Data=5010001_02200009247600003190000000;
cutNumberQA0010_Data=0005314141_5010001_02200009247600003190000000;

fileName1020_Data=PhotonQA23March_0005314141_5120001_02200009247600003190000000_Data.root;
cutNumber1020_Data=5120001_02200009247600003190000000;
cutNumberQA1020_Data=0005314141_5120001_02200009247600003190000000;

fileName2040_Data=PhotonQA23March_0005314141_5240001_02200009247600003190000000_Data.root;
cutNumber2040_Data=5240001_02200009247600003190000000;
cutNumberQA2040_Data=0005314141_5240001_02200009247600003190000000;

fileName4060_Data=PhotonQA23March_0005314141_5460001_02200009247600003190000000_Data.root;
cutNumber4060_Data=5460001_02200009247600003190000000;
cutNumberQA4060_Data=0005314141_5460001_02200009247600003190000000;

fileName0010_MC=PhotonQA23March_0005314141_5010001_02200009247600003190000000_MC.root;
cutNumber0010_MC=5010001_02200009247600003190000000;
cutNumberQA0010_MC=0005314141_5010001_02200009247600003190000000;

fileName1020_MC=PhotonQA23March_0005314141_5120001_02200009247600003190000000_MC.root;
cutNumber1020_MC=5120001_02200009247600003190000000;
cutNumberQA1020_MC=0005314141_5120001_02200009247600003190000000;

fileName2040_MC=PhotonQA23March_0005314141_5240001_02200009247600003190000000_MC.root;
cutNumber2040_MC=5240001_02200009247600003190000000;
cutNumberQA2040_MC=0005314141_5240001_02200009247600003190000000;

fileName4060_MC=PhotonQA23March_0005314141_5460001_02200009247600003190000000_MC.root;
cutNumber4060_MC=5460001_02200009247600003190000000;
cutNumberQA4060_MC=0005314141_5460001_02200009247600003190000000;

filedEdxDATA0010=dEdxRunwise_0010_Data.txt;
GammaTPCSectorsDATA0010=GammaTPCSectors_0010_Data.txt;
GammaTPCSectorsErrorDATA0010=GammaTPCSectorsErrors_0010_Data.txt;
filedEdxMC0010=dEdxRunwise_0010_MC.txt;
GammaTPCSectorsMC0010=GammaTPCSectors_0010_MC.txt;
GammaTPCSectorsErrorMC0010=GammaTPCSectorsErrors_0010_MC.txt;

filedEdxDATA1020=dEdxRunwise_1020_Data.txt;
GammaTPCSectorsDATA1020=GammaTPCSectors_1020_Data.txt;
GammaTPCSectorsErrorDATA1020=GammaTPCSectorsErrors_1020_Data.txt;
filedEdxMC1020=dEdxRunwise_1020_MC.txt;
GammaTPCSectorsMC1020=GammaTPCSectors_1020_MC.txt;
GammaTPCSectorsErrorMC1020=GammaTPCSectorsErrors_1020_MC.txt;

filedEdxDATA2040=dEdxRunwise_2040_Data.txt;
GammaTPCSectorsDATA2040=GammaTPCSectors_2040_Data.txt;
GammaTPCSectorsErrorDATA2040=GammaTPCSectorsErrors_2040_Data.txt;
filedEdxMC2040=dEdxRunwise_2040_MC.txt;
GammaTPCSectorsMC2040=GammaTPCSectors_2040_MC.txt;
GammaTPCSectorsErrorMC2040=GammaTPCSectorsErrors_2040_MC.txt;

filedEdxDATA4060=dEdxRunwise_4060_Data.txt;
GammaTPCSectorsDATA4060=GammaTPCSectors_4060_Data.txt;
GammaTPCSectorsErrorDATA4060=GammaTPCSectorsErrors_4060_Data.txt;
filedEdxMC4060=dEdxRunwise_4060_MC.txt;
GammaTPCSectorsMC4060=GammaTPCSectors_4060_MC.txt;
GammaTPCSectorsErrorMC4060=GammaTPCSectorsErrors_4060_MC.txt;


fileDIRData1=/home/admin1/leardini/GridOutput/PbPb/Legotrain-vAN-20150310/GA_PbPb_AOD-114_20150311-1642;  
fileDIRData2=/home/admin1/leardini/GridOutput/PbPb/Legotrain-vAN-20150310/GA_PbPb_AOD-115_20150311-1642;

fileDIRMC1=/home/admin1/leardini/GridOutput/PbPb/Legotrain-vAN-20150326/GA_PbPb_MC_AOD-317_20150327-1132;
fileDIRMC2=/home/admin1/leardini/GridOutput/PbPb/Legotrain-vAN-20150326/GA_PbPb_MC_AOD-318_20150327-1120;

# fileDIRData=/home/admin1/leardini/GridOutput/PbPb/Legotrain-vAN-20141104/GA_PbPb_AOD-77_20141105-0956;      
# fileDIRMC1=/home/admin1/leardini/GridOutput/PbPb/Legotrain-vAN-20141104/GA_PbPb_MC_AOD-266_20141105-1002;
# fileDIRMC2=/home/admin1/leardini/GridOutput/PbPb/Legotrain-vAN-20141104/GA_PbPb_MC_AOD-268_20141105-1003;

# fileDIRData=/home/admin1/leardini/GridOutput/PbPb/Legotrain-vAN-20141104/GA_PbPb_AOD-77_20141105-0956;      
# fileDIRMC1=/home/admin1/leardini/GridOutput/PbPb/Legotrain-vAN-20141104/GA_PbPb_MC_AOD-266_20141105-1002;
# fileDIRMC2=/home/admin1/leardini/GridOutput/PbPb/Legotrain-vAN-20141104/GA_PbPb_MC_AOD-268_20141105-1003;
# fileDIRData=/home/admin1/leardini/GridOutput/PbPb/Legotrain-vAN-20140609/GA_PbPb_AOD-62_20140610-1406;      
# fileDIRMC1=/home/admin1/leardini/GridOutput/PbPb/Legotrain-vAN-20140609/GA_PbPb_MC_AOD-245_20140610-1403;
# fileDIRMC2=/home/admin1/leardini/GridOutput/PbPb/Legotrain-vAN-20140609/GA_PbPb_MC_AOD-246_20140610-1403;
# outDIR=/home/admin1/leardini/Results/$1; 
outDIR=/home/admin1/leardini/Results/QA26Marzo/QArunwiseEcut13;
mkdir -p $outDIR;


#if [ $1 = "AnalysisQA" ]; then
	runNumbers=`cat listOfRunsLHC11h.txt` 
	for runNumber in $runNumbers; do
	   echo "###########################################"  $runNumber  "for data ####################################################"
 	   echo "0-10%"
	   root -l -b -q -x ExtractRunwiseQA.C++\(\"Data\"\,\"$fileName0010_Data\"\,\"$cutNumber0010_Data\"\,\"$cutNumberQA0010_Data\"\,\"$fileDIRData1\"\,\"$outDIR\"\,$runNumber\,\"pdf\"\,\"0010\"\)
 	   echo "10-20%"
 	   root -l -b -q -x ExtractRunwiseQA.C++\(\"Data\"\,\"$fileName1020_Data\"\,\"$cutNumber1020_Data\"\,\"$cutNumberQA1020_Data\"\,\"$fileDIRData2\"\,\"$outDIR\"\,$runNumber\,\"pdf\"\,\"1020\"\)
 	   echo "20-40%"
 	   root -l -b -q -x ExtractRunwiseQA.C++\(\"Data\"\,\"$fileName2040_Data\"\,\"$cutNumber2040_Data\"\,\"$cutNumberQA2040_Data\"\,\"$fileDIRData2\"\,\"$outDIR\"\,$runNumber\,\"pdf\"\,\"2040\"\)
  	   echo "40-60%"
  	   root -l -b -q -x ExtractRunwiseQA.C++\(\"Data\"\,\"$fileName4060_Data\"\,\"$cutNumber4060_Data\"\,\"$cutNumberQA4060_Data\"\,\"$fileDIRData2\"\,\"$outDIR\"\,$runNumber\,\"pdf\"\,\"4060\"\)

  	   echo "#############################################"  $runNumber  "for MC ##################################################"
 	   echo "0-10%"
   	   root -l -b -q -x ExtractRunwiseQA.C++\(\"MC\"\,\"$fileName0010_MC\"\,\"$cutNumber0010_MC\"\,\"$cutNumberQA0010_MC\"\,\"$fileDIRMC1\"\,\"$outDIR\"\,$runNumber\,\"pdf\"\,\"0010\"\)
 	   echo "10-20%"
 	   root -l -b -q -x ExtractRunwiseQA.C++\(\"MC\"\,\"$fileName1020_MC\"\,\"$cutNumber1020_MC\"\,\"$cutNumberQA1020_MC\"\,\"$fileDIRMC2\"\,\"$outDIR\"\,$runNumber\,\"pdf\"\,\"1020\"\)
 	   echo "20-40%"
 	   root -l -b -q -x ExtractRunwiseQA.C++\(\"MC\"\,\"$fileName2040_MC\"\,\"$cutNumber2040_MC\"\,\"$cutNumberQA2040_MC\"\,\"$fileDIRMC2\"\,\"$outDIR\"\,$runNumber\,\"pdf\"\,\"2040\"\)
 	   echo "40-60%"
 	   root -l -b -q -x ExtractRunwiseQA.C++\(\"MC\"\,\"$fileName4060_MC\"\,\"$cutNumber4060_MC\"\,\"$cutNumberQA4060_MC\"\,\"$fileDIRMC2\"\,\"$outDIR\"\,$runNumber\,\"pdf\"\,\"4060\"\)
	done
# 
# #elif [ $1 = "Plotting" ]; then
	echo "Here start the plotting part"
	root -l -b -q -x  PutToGraphRunwiseQA.C++\(\"Data\"\,\"$outDIR/$filedEdxDATA0010\"\,\"$outDIR/$GammaTPCSectorsDATA0010\"\,\"$outDIR/$GammaTPCSectorsErrorDATA0010\"\,\"0010\"\,\"$outDIR\"\,\"pdf\"\)
	root -l -b -q -x  PutToGraphRunwiseQA.C++\(\"Data\"\,\"$outDIR/$filedEdxDATA1020\"\,\"$outDIR/$GammaTPCSectorsDATA1020\"\,\"$outDIR/$GammaTPCSectorsErrorDATA1020\"\,\"1020\"\,\"$outDIR\"\,\"pdf\"\)
	root -l -b -q -x  PutToGraphRunwiseQA.C++\(\"Data\"\,\"$outDIR/$filedEdxDATA2040\"\,\"$outDIR/$GammaTPCSectorsDATA2040\"\,\"$outDIR/$GammaTPCSectorsErrorDATA2040\"\,\"2040\"\,\"$outDIR\"\,\"pdf\"\)
	root -l -b -q -x  PutToGraphRunwiseQA.C++\(\"Data\"\,\"$outDIR/$filedEdxDATA4060\"\,\"$outDIR/$GammaTPCSectorsDATA4060\"\,\"$outDIR/$GammaTPCSectorsErrorDATA4060\"\,\"4060\"\,\"$outDIR\"\,\"pdf\"\)
	
	root -l -b -q -x  PutToGraphRunwiseQA.C++\(\"MC\"\,\"$outDIR/$filedEdxMC0010\"\,\"$outDIR/$GammaTPCSectorsMC0010\"\,\"$outDIR/$GammaTPCSectorsErrorMC0010\"\,\"0010\"\,\"$outDIR\"\,\"pdf\"\)
	root -l -b -q -x  PutToGraphRunwiseQA.C++\(\"MC\"\,\"$outDIR/$filedEdxMC1020\"\,\"$outDIR/$GammaTPCSectorsMC1020\"\,\"$outDIR/$GammaTPCSectorsErrorMC1020\"\,\"1020\"\,\"$outDIR\"\,\"pdf\"\)
	root -l -b -q -x  PutToGraphRunwiseQA.C++\(\"MC\"\,\"$outDIR/$filedEdxMC2040\"\,\"$outDIR/$GammaTPCSectorsMC2040\"\,\"$outDIR/$GammaTPCSectorsErrorMC2040\"\,\"2040\"\,\"$outDIR\"\,\"pdf\"\)
	root -l -b -q -x  PutToGraphRunwiseQA.C++\(\"MC\"\,\"$outDIR/$filedEdxMC4060\"\,\"$outDIR/$GammaTPCSectorsMC4060\"\,\"$outDIR/$GammaTPCSectorsErrorMC4060\"\,\"4060\"\,\"$outDIR\"\,\"pdf\"\)
#fi

echo "Here start the drawing part"
root -l -b -q -x PutGraphTogetherQA.C++\(kFALSE\,\"$outDIR\"\)
root -l -b -q -x PutGraphTogetherQA.C++\(kTRUE\,\"$outDIR\"\)