#!/usr/local/bin/perl


#$trainDateVersionData="2011-10-28_1917.9959_GammaConv";
#$trainDateVersionDataE="2011-10-28_1917.9959_GammaConv";
#$trainDateVersionMC="2011-10-28_1921.9959_GammaConv";

#$trainDateVersionData="2011-11-03_1924.10120";
#$trainDateVersionDataE="2011-11-03_1924.10120";
#$trainDateVersionMC="2011-11-03_1943.10120";


#$trainDateVersionData="2011-11-08_0002.10238";
#$trainDateVersionDataE="2011-11-08_0002.10238";
#$trainDateVersionMC="2011-11-08_0016.10237";

#$trainDateVersionData="2011-11-17_2005.10417";
#$trainDateVersionDataE="2011-11-17_2005.10417";
#$trainDateVersionMC="2011-11-17_2007.10418";

# $trainDateVersionData="2011-11-24_1627.10535";
# $trainDateVersionDataE="2011-11-24_1627.10535";
# $trainDateVersionMC="2011-11-24_1634.10535";

# $trainDateVersionData="2011-12-02_2219.10647";
# $trainDateVersionDataE="2011-12-02_2219.10647";
# $trainDateVersionMC="2011-12-02_2214.10647";
# $trainDateVersionMCE="2011-12-02_1035.10647";

# central train: Dalitz + Direct Photon
# $trainDateVersionData="2011-12-06_2126.10765";
# $trainDateVersionDataE="2011-12-06_2126.10765";
# $trainDateVersionMC="2011-12-06_2141.10764";
# $trainDateVersionMCE="2011-12-06_2141.10764";


# central train: Dalitz + Martins DirectPhoton cut, Qt cut based on ESDV0s
# $trainDateVersionData="2011-12-13_2044.10865";
# $trainDateVersionDataE="2011-12-13_2044.10865";
# $trainDateVersionMC="2011-12-13_2049.10864";
# $trainDateVersionMCE="2011-12-13_2049.10864";

# central train: Material + Martins DirectPhoton cut, Qt cut based on ESDV0s, new MC's and additional Resolution
# $trainDateVersionData="2011-12-21_0114.11049";
# $trainDateVersionDataE="2011-12-21_0114.11049";
# $trainDateVersionMC="2011-12-20_2210.11047";
# $trainDateVersionMCE="2011-12-20_2210.11047";

# private train: Material Cut + Martins DirectPhoton in eta 0.9, 0.9 < eta < 1.4, new MC's and additional Resolution, only LHC10e, LHC10c for data, had interaction
# $trainDateVersionData="2011-12-22_0938.11051";
# $trainDateVersionDataE="2011-12-22_0938.11051";
# $trainDateVersionMC="2011-12-22_0950.11051";
# $trainDateVersionMCE="2011-12-22_0950.11051";

# private train: Material Cuts in eta 0.9, 0.9 < eta < 1.4, offline, new MC's and additional Resolution, only LHC10e, LHC10c for data, had interaction
# $trainDateVersionData="2011-12-22_1915.11051";
# $trainDateVersionDataE="2011-12-22_1915.11051";
# $trainDateVersionMC="2011-12-22_1914.11051";
# $trainDateVersionMCE="2011-12-22_1914.11051";

# central train: 3* Material 0.9+ extended reach
# $trainDateVersionData="2012-01-06_0014.11113";
# $trainDateVersionDataE="2012-01-06_0014.11113";
# $trainDateVersionMC="2012-01-06_0033.11114";
# $trainDateVersionMCE="2012-01-06_0033.11114";

# central train: pile up rejection tested again
#$trainDateVersionData="2012-05-11_0026.13613";
#$trainDateVersionDataE="2012-05-11_0026.13613";
#$trainDateVersionMC="2012-05-12_1740.13635";
#$trainDateVersionMCE="2012-05-12_1740.13635";

# private train Dalitz new
#$trainDateVersionData="2012-08-27_1955.15418";
#$trainDateVersionDataE="2012-08-27_1955.15418";
#$trainDateVersionMC="2012-08-24_1046.15418";
#$trainDateVersionMCE="2012-08-24_1046.15418";

# private train MC production only Material budget Offline V0 finder
#$trainDateVersionData="2012-12-02_1712.16824";
#$trainDateVersionDataE="2012-12-02_1712.16824";
#$trainDateVersionMC="2012-12-02_1712.16824";
#$trainDateVersionMCE="2012-12-02_1712.16824";

# private train MC production only Material budget Offline V0 finder
# $trainDateVersionData="2012-12-02_1735.16824";
# $trainDateVersionDataE="2012-12-02_1735.16824";
# $trainDateVersionMC="2012-12-02_1735.16824";
# $trainDateVersionMCE="2012-12-02_1735.16824";

# central train x-checks new Software
# $trainDateVersionData="2013-02-18_2217.17676";
# $trainDateVersionDataE="2013-02-18_2217.17676";
# $trainDateVersionMC="2013-02-18_2234.17676";
# $trainDateVersionMCE="2013-02-18_2234.17676";

# private train x-checks new Software V0AND/V0OR
# $trainDateVersionData="2013-02-21_1932.17711";
# $trainDateVersionDataE="2013-02-21_1932.17711";
# $trainDateVersionMC="2013-02-21_1939.17711";
# $trainDateVersionMCE="2013-02-21_1939.17711";

# private train standard cuts 7 TeV cut, 7 TeV and 2.76TeV dataset
# 0000011002093663003800000_01631031009
# 0002011002093663003800000_01631031009
# 0000011002083663000200000_02631031009
# 0000012002093663003800000_01631031009
# 0002012002093663003800000_01631031009
# 0000012002083663000200000_02631031009
# $trainDateVersionData="2013-04-16_1039.18096";
# $trainDateVersionDataE="2013-04-16_1039.18096";
# $trainDateVersionMC="2013-04-13_2012.18096"; 
# $trainDateVersionMCE="2013-04-13_2012.18096"; 

# private train offline V0finder standard cuts 7 TeV cut, 7 TeV and 2.76TeV dataset
# 0000011102093663003800000_01631031009
# 0002011102093663003800000_01631031009
# 0000011102083663000200000_02631031009
# 0000012102093663003800000_01631031009
# 0002012102093663003800000_01631031009
# 0000012102083663000200000_02631031009
# $trainDateVersionData="2013-04-16_1225.18096";
# $trainDateVersionMC="2013-04-12_1251.18096"; 

#private train: test of dca of photons & mesons
$trainDateVersionData="2013-08-07_1742.19409";
$trainDateVersionMC="2013-08-07_1748.19409"; 


# system("echo $PWD");
$outputMergedBase="/home/fredi/Photon/Results/ppDCATests/ppRootFiles/$trainDateVersionData";
system("echo $outputMergedBase");

system("mkdir -p $outputMergedBase/mergedBC");
system("mkdir -p $outputMergedBase/mergedDE");
system("mkdir -p $outputMergedBase/mergedB");
system("mkdir -p $outputMergedBase/mergedC");
system("mkdir -p $outputMergedBase/mergedE");
system("mkdir -p $outputMergedBase/mergedD");
system("mkdir -p $outputMergedBase/mergedALL");


# Merging data 7TeV
$energy="7TeV";
system("ln -fs $outputMergedBase/GammaConvV1-$trainDateVersionData-$energy-LHC10b.root $outputMergedBase/mergedB/GammaConvV1Data.root");
system("ln -fs $outputMergedBase/GammaConvV1-$trainDateVersionData-$energy-LHC10c.root $outputMergedBase/mergedC/GammaConvV1Data.root");
system("ln -fs $outputMergedBase/GammaConvV1-$trainDateVersionData-$energy-LHC10d.root $outputMergedBase/mergedD/GammaConvV1Data.root");
system("ln -fs $outputMergedBase/GammaConvV1-$trainDateVersionData-$energy-LHC10e.root $outputMergedBase/mergedE/GammaConvV1Data.root");
system("ln -fs $outputMergedBase/GammaConvV1-$trainDateVersionMC-$energy-LHC10d1d2d4d4a.root $outputMergedBase/mergedBC/GammaConvV1MC.root");

system("ln -fs $outputMergedBase/GammaConvV1-$trainDateVersionData-$energy-LHC10de.root $outputMergedBase/mergedDE/GammaConvV1Data.root");
system("ln -fs $outputMergedBase/GammaConvV1-$trainDateVersionMC-$energy-LHC10f6f6ae20e21.root $outputMergedBase/mergedDE/GammaConvV1MC.root");

system("ln -fs $outputMergedBase/GammaConvV1-$trainDateVersionData-$energy-LHC10bcde.root $outputMergedBase/mergedALL/GammaConvV1Data.root");
system("ln -fs $outputMergedBase/GammaConvV1-$trainDateVersionMC-$energy-LHC10d1d2d4d4af6f6ae20e21.root $outputMergedBase/mergedALL/GammaConvV1MC.root");

