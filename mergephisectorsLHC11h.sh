#! /bin/bash

BASEDIR=/home/admin1/leardini/GridOutput/PbPb
TRAINDIR=Legotrain-vAN-20160211-1;
OUTPUTDIR=$BASEDIR/$TRAINDIR;
OUTPUTDIRmerged=$BASEDIR/$TRAINDIR/mergeFilesTightCuts;
mkdir -p $OUTPUTDIRmerged

# OUTPUTDIRold=$BASEDIR/Legotrain-vAN-20150507;
# OUTPUTDIRFlat=$OUTPUTDIRold/fileCentFlat;  #$BASEDIR/$TRAINDIR/fileCentFlat;

# LHC11hData=169_20150612-1147; #142_20150510-1013; #122_20150417-0930; # 121_20150417-0928
# LHC11hDataPhi=170_20150612-1147; #143_20150510-1051; #124_20150417-0930;  #123_20150417-0929; # 
# LHC14a1awithEta=365_20150610-1352; #333_20150503-1157; #and with no phi cut
# LHC14a1bwithEta=367_20150610-1352; #335_20150503-1057;
# LHC14a1awithPi0=366_20150610-1408; # 334_20150503-1221; #and with phi cut
# LHC14a1bwithPi0=368_20150610-1402; #336_20150503-1057; 


LHC11hData=260_20160211-1824; 
LHC11hDataPhi=261_20160211-1825; 

LHC14a1awithEta=632_20160211-1832; 
LHC14a1bwithEta=634_20160211-1833;
LHC14a1awithPi0=633_20160211-1832;
LHC14a1bwithPi0=635_20160211-1833 ;



#-> for cross check of weighting effect:
# LHC11hData=250_20160121-1452; 
# LHC11hDataPhi=251_20160121-1452; 

# LHC14a1awithEta=627_20160121-1450; 
# LHC14a1bwithEta=628_20160121-1450;
# LHC14a1awithPi0=626_20160121-1450;
# LHC14a1bwithPi0=625_20160121-1450;


#std slow, double rej, mult not cut
# LHC11hData=241_20151209-1800; #cut 226, 230, 235, 239
# LHC11hDataPhi=242_20151209-1704; 

#first syst check
# LHC14a1awithEta=606_20151219-1059; 
# LHC14a1bwithEta=608_20151219-1059;
# LHC14a1awithPi0=607_20151219-1059;
# LHC14a1bwithPi0=609_20151219-1059;

#second syst check
# LHC14a1awithEta=611_20151220-1102; 
# LHC14a1bwithEta=613_20151220-1103;
# LHC14a1awithPi0=617_20160103-1107;
# LHC14a1bwithPi0=614_20151220-1103;

#third syst check
# LHC11hData=244_20160105-1530; #cut 235 nuovo
# LHC11hDataPhi=245_20160105-1530; 
# 
# LHC14a1awithEta=618_20160104-1012; 
# LHC14a1bwithEta=622_20160104-1212;
# LHC14a1awithPi0=619_20160104-1013;
# LHC14a1bwithPi0=621_20160104-1213;

#Chi2 30. and psi pair 1D (not linked)
# LHC14a1awithEta=602_20151209-1801; 
# LHC14a1bwithEta=604_20151209-1819;
# LHC14a1awithPi0=603_20151209-1711;
# LHC14a1bwithPi0=605_20151209-1818;


#std slow, double rej, mult not cut
# LHC14a1awithEta=598_20151210-0914; 
# LHC14a1bwithEta=600_20151210-1012;
# LHC14a1awithPi0=599_20151210-0914;
# LHC14a1bwithPi0=601_20151210-1033;



#========= std with V0 event cut (slow train) =========
# LHC11hData=238_20151126-1238; 
# LHC11hDataPhi=239_20151126-1041; 
# 
# LHC14a1awithEta=591_20151126-1347; 
# LHC14a1bwithEta=593_20151126-1042;
# LHC14a1awithPi0=592_20151126-1042;
# LHC14a1bwithPi0=594_20151126-1042;

#========= std with no PsiPair with Chi2 cut ========= 
# LHC11hData=235_20151119-1719;
# LHC11hDataPhi=236_20151119-1719;
# 
# LHC14a1awithEta=580_20151119-1745;
# LHC14a1bwithEta=582_20151119-1746;
# LHC14a1awithPi0=581_20151119-1757;
# LHC14a1bwithPi0=583_20151119-1747;

#========= std with with PsiPair no Chi2 cut ========= 
# LHC11hData=235_20151119-1719;
# LHC11hDataPhi=236_20151119-1719;

# LHC14a1awithEta=585_20151120-1018;
# LHC14a1bwithEta=587_20151120-1019;
# LHC14a1awithPi0=586_20151120-1018;
# LHC14a1bwithPi0=588_20151120-1019;


#========= std with no PsiPair/Chi2 cut ========= -> wrong! too much with both off
# LHC11hData=231_20151109-1015;
# LHC11hDataPhi=232_20151109-1116;
# 
# LHC14a1awithEta=571_20151109-1037;
# LHC14a1bwithEta=573_20151109-1037;
# LHC14a1awithPi0=572_20151109-1037;
# LHC14a1bwithPi0=574_20151109-1038;


#========= std with new V0cut and double rejection =========
#226, 227, 228, 229
# LHC11hData=226_20151102-1458;
# LHC11hDataPhi=227_20151102-1458;

# LHC14a1awithEta=563_20151102-1415;
# LHC14a1bwithEta=565_20151102-1416;
# LHC14a1awithPi0=564_20151102-1416;
# LHC14a1bwithPi0=566_20151102-1416;

#standard (data file is same)
# LHC14a1awithEta=567_20151103-1029;
# LHC14a1bwithEta=569_20151103-1030;
# LHC14a1awithPi0=568_20151103-1030;
# LHC14a1bwithPi0=570_20151103-1030;

#=== open dEdx cut ====
# 	LHC11hData=192_20150727-1243;
# 	LHC11hDataPhi=193_20150727-1243;
# 	LHC14a1awithEta=492_20150727-1059;
# 	LHC14a1bwithEta=494_20150727-1100;
# 	LHC14a1awithPi0=493_20150727-1059;
# 	LHC14a1bwithPi0=495_20150727-1100;

#=== R > 35cm ====
# 	LHC11hData=192_20150727-1243;
# 	LHC11hDataPhi=193_20150727-1243;
# 	LHC14a1awithEta=496_20150727-1100;
# 	LHC14a1bwithEta=498_20150727-1101;
# 	LHC14a1awithPi0=497_20150727-1101;
# 	LHC14a1bwithPi0=499_20150727-1101;

#=== re-do cut 146 (new tag) ====
# 	LHC11hData=194_20150729-1034;
# 	LHC11hDataPhi=195_20150729-1034;
# 	LHC14a1awithEta=500_20150729-0907;
# 	LHC14a1bwithEta=502_20150729-0908;
# 	LHC14a1awithPi0=501_20150729-0908;
# 	LHC14a1bwithPi0=503_20150729-0908;


#=========== xrootd fix =========================
# LHC11hData=215_20150921-1117;
# LHC11hDataPhi=216_20150921-1116;
# 
# LHC14a1awithEta=527_20150921-1117;
# LHC14a1bwithEta=529_20150921-1117;
# LHC14a1awithPi0=528_20150921-1132;
# LHC14a1bwithPi0=530_20150921-1117;

# LHC11hData=213_20150907-1315; #std 
# LHC11hDataPhi=214_20150907-1325; #std
# 
# LHC14a1awithEta=521_20150907-1325;  #std
# LHC14a1bwithEta=524_20150907-1341; #std
# LHC14a1awithPi0=522_20150907-1340; #std
# LHC14a1bwithPi0=525_20150907-1340; #std


#======================================================== systematics

#======= standard ========
# LHC11hData=169_20150612-1147;  
# LHC11hDataPhi=170_20150612-1147; 
#LHC14a1awithEta=365_20150610-1352;  
#LHC14a1bwithEta=367_20150610-1352; 
#LHC14a1awithPi0=366_20150610-1408; 
#LHC14a1bwithPi0=368_20150610-1402; 


#special - phi cut variation =========
# 	LHC11hData=169_20150612-1147; # std for no phi cut
# 	LHC11hDataPhi=187_20150630-1308; #186, 188

#======= phi 2.0 - 4.0 ========
# 	LHC14a1awithEta=457_20150711-1550;  #186, 187eta, 189eta
# 	LHC14a1bwithEta=459_20150711-1600; 

#======= phi ========
# 	LHC14a1awithPi0=458_20150711-1558; #188, 187pi0, 189pi0
# 	LHC14a1bwithPi0=460_20150711-1601; 

#======= eta 0.65 ========
##data cut is with std analysis
# LHC14a1awithEta=369_20150612-1010;  #70, 71eta, 73eta
# LHC14a1bwithEta=371_20150612-1008; 
# LHC14a1awithPi0=466_20150719-1912; #wrong: 370_20150612-1009; #72, 71pi0, 73pi0
# LHC14a1bwithPi0=467_20150719-1912; #wrong 372_20150612-1007; 

#================================================================
# LHC11hData=171_20150614-1506; 	 #74, 82, 86
# LHC11hDataPhi=172_20150614-1508; #76, 84, 88
#======= eta 0.75 ========
# LHC14a1awithEta=373_20150614-1513;  #74, 75eta, 77eta
# LHC14a1bwithEta=375_20150614-1514; 
# LHC14a1awithPi0=374_20150614-1513; #76, 75pi0, 77pi0
# LHC14a1bwithPi0=376_20150614-1515; 
#======= single pt 0.075 ========
# LHC14a1awithEta=377_20150615-0837;  #82, 83eta, 85eta
# LHC14a1bwithEta=379_20150615-0838; 
# LHC14a1awithPi0=378_20150615-0837; # 84, 83pi0, 85pi0
# LHC14a1bwithPi0=380_20150615-0838; 
#======= single pt 0.1 ========	
# LHC14a1awithEta=381_20150616-1151; #86, 87eta, 89eta  
# LHC14a1bwithEta=383_20150616-1152; 
# LHC14a1awithPi0=382_20150616-1151; #88, 87pi0, 89pi0 
# LHC14a1bwithPi0=384_20150616-1153; 

#================================================================
# LHC11hData=190_20150720-1040; #better stat 173_20150615-1331; 	 #90, 94, 98 
# LHC11hDataPhi=191_20150720-1418; #better stat 174_20150615-1332; #92, 96, 100
#======= TPC cls 0.7 ========
# LHC14a1awithEta=385_20150617-0811; #90, 91eta, 93eta  
# LHC14a1bwithEta=387_20150617-0812; 
# LHC14a1awithPi0=386_20150617-0811; #92, 91pi0, 93pi 
# LHC14a1bwithPi0=388_20150617-0813; 
#======= TPC cls 0.35 ========
# LHC14a1awithEta=389_20150619-2053;  #94, 95eta, 97eta
# LHC14a1bwithEta=391_20150619-2054; 
# LHC14a1awithPi0=390_20150619-2054;  #96, 95pi, 97pi
# LHC14a1bwithPi0=392_20150619-2056; 
#======= edEdx -4, +5 ========
# LHC14a1awithEta=393_20150624-1059; #98,  99eta, 101eta
# LHC14a1bwithEta=395_20150624-1100; 
# LHC14a1awithPi0=394_20150624-1059; #100, 99pi,  101pi
# LHC14a1bwithPi0=396_20150624-1100; 

#================================================================
# LHC11hData=175_20150615-1939;     #162 not flat, 102
# LHC11hDataPhi=176_20150615-1940;  #164 not flat, 104
#======= edEdx -2.5, +4 ========
# LHC14a1awithEta=397_20150624-1312;  #102, 103eta, 105eta
# LHC14a1bwithEta=399_20150624-1313; 
# LHC14a1awithPi0=398_20150624-1312;  #104, 103pi, 105pi
# LHC14a1bwithPi0=400_20150624-1314; 

#================================================================
# LHC11hData=177_20150619-2051;    #106, 110, 114
# LHC11hDataPhi=189_20150720-1037; #better stat 178_20150619-2051; #108, 112, 116
#======= pi dEdx 2., +1 ========
# LHC14a1awithEta=401_20150626-1142;  #106, 107eta, 109eta
# LHC14a1bwithEta=403_20150626-1144; 
# LHC14a1awithPi0=402_20150626-1144;  #108, 107pi0, 109pi
# LHC14a1bwithPi0=404_20150626-1145; 
#======= pi dEdx 2.5, -- ========
# LHC14a1awithEta=405_20150630-1018;  #110, 111eta, 113eta
# LHC14a1bwithEta=407_20150630-1019; 
# LHC14a1awithPi0=406_20150630-1019;  #112, 111pi,  113pi
# LHC14a1bwithPi0=408_20150630-1020; 
#======= max P pion ========
# LHC14a1awithEta=409_20150630-1313;  #114, 115ta, 117eta
# LHC14a1bwithEta=411_20150630-1314; 
# LHC14a1awithPi0=410_20150630-1313; # 116, 115pi, 117pi
# LHC14a1bwithPi0=412_20150630-1314; 

#================================================================
# LHC11hData=179_20150630-1028;    #118, 122, 150
# LHC11hDataPhi=180_20150630-1029; #120, 124, 152
#======= TOF PID -3, +5 ========
# LHC14a1awithEta=468_20150719-1927; #wrong: 413_20150630-1316; #118, 119eta, 121eta
# LHC14a1bwithEta=470_20150719-1939  #wrong: 415_20150630-1318;
# LHC14a1awithPi0=469_20150719-1917; #worng: 414_20150630-1316; #120, 119pi,  121pi
# LHC14a1bwithPi0=471_20150719-2149; #wrong: 416_20150630-1319; 
#======= TOF PID -2, +3 ======== 
# LHC14a1awithEta=488_20150723-1044; #wrong 417_20150630-1453;   #122, 123eta, 125eta
# LHC14a1bwithEta=490_20150723-1042; #wrong 419_20150630-1517; 
# LHC14a1awithPi0=489_20150723-1044; #wrong 418_20150630-1514;   #124, 123pi,  125pi
# LHC14a1bwithPi0=491_20150721-1656; #wrong 420_20150630-1455; 
#======= cos P.A. -1 ========
# LHC14a1awithEta=421_20150630-1510;  #150, 151eta, 153eta
# LHC14a1bwithEta=423_20150630-1511; 
# LHC14a1awithPi0=422_20150630-1511;  #152, 151pi, 153pi 
# LHC14a1bwithPi0=424_20150630-1512; 

#================================================================ ----------------> data redone, something off; different qt cut too!
# LHC11hData=196_20150731-0835;    #126, 130, 134  
# LHC11hDataPhi=197_20150731-0835; #128, 132, 136
#======= qT 0.03 2D ========
# LHC14a1awithEta=425_20150702-1632; #126, 127eta, 129eta  
# LHC14a1bwithEta=427_20150702-1632; 
# LHC14a1awithPi0=426_20150702-1632; #128, 127pi,  129pi
# LHC14a1bwithPi0=428_20150702-1633; 
# #======= qT 0.06 2D ========
#new qt cut!!!
# LHC14a1awithEta=504_20150730-2013;#130, 131eta, 133eta
# LHC14a1bwithEta=506_20150730-2014;
# LHC14a1awithPi0=505_20150730-2013;#132, 131pi, 133pi
# LHC14a1bwithPi0=507_20150730-2014; 
#======= chi2 50. ========
# LHC14a1awithEta=476_20150720-1740 #wrong 433_20150702-1418;  #134, 135eta, 137eta
# LHC14a1bwithEta=478_20150721-1026 #wrong 435_20150702-1420; 
# LHC14a1awithPi0=477_20150720-1815 #wrong 434_20150702-1419;  #136, 135pi, 137pi
# LHC14a1bwithPi0=479_20150720-1740 #wrong 436_20150702-1422; 

# # LHC11hData=181_20150630-1032;    #126, 130, 134  
# # LHC11hDataPhi=182_20150630-1032; #128, 132, 136
# #======= qT 0.03 2D ========
# # LHC14a1awithEta=425_20150702-1632; #126, 127eta, 129eta  
# # LHC14a1bwithEta=427_20150702-1632; 
# # LHC14a1awithPi0=426_20150702-1632; #128, 127pi,  129pi
# # LHC14a1bwithPi0=428_20150702-1633; 
# #======= qT 0.07 1D ========							
# # LHC14a1awithEta=472_20150719-2149 #wrong 429_20150702-1633;  #130, 131eta, 133eta
# # LHC14a1bwithEta=474_20150720-1024 #wrong 431_20150702-1633; 
# # LHC14a1awithPi0=473_20150720-1019 #wrong 429_20150702-1633;  #132, 131pi,  133pi
# # LHC14a1bwithPi0=475_20150720-1030 #wrong 432_20150702-1633; 
# #======= chi2 50. ========
# # LHC14a1awithEta=476_20150720-1740 #wrong 433_20150702-1418;  #134, 135eta, 137eta
# # LHC14a1bwithEta=478_20150721-1026 #wrong 435_20150702-1420; 
# # LHC14a1awithPi0=477_20150720-1815 #wrong 434_20150702-1419;  #136, 135pi, 137pi
# # LHC14a1bwithPi0=479_20150720-1740 #wrong 436_20150702-1422; 

#================================================================
# LHC11hData=194_20150729-1034; #138, 142, 146 
# LHC11hDataPhi=195_20150729-1034; #140, 144, 148


#wrong?
# LHC11hData=183_20150630-1258;    #138, 142, 146   
# LHC11hDataPhi=184_20150630-1258; #140, 144, 148
#======= chi2 20. ========
# LHC14a1awithEta=480_20150721-1027; #wrong 437_20150702-1527;  #138, 139eta, 141eta
# LHC14a1bwithEta=482_20150721-1631 #wrong 439_20150702-1528; 
# LHC14a1awithPi0=481_20150721-1101; #wrong 438_20150702-1527;  #140, 139pi, 141pi
# LHC14a1bwithPi0=483_20150722-2156 #wrong 440_20150702-1528; 
#======= psi 0.05 2D ========
# LHC14a1awithEta=484_20150723-1042; # wrong 441_20150702-1530;  #142, 143eta, 145eta
# LHC14a1bwithEta=486_20150723-1043; #wrong 443_20150702-1531; 
# LHC14a1awithPi0=485_20150723-1043; # wrong 442_20150702-1530;  #144, 143pi, 145pi
# LHC14a1bwithPi0=487_20150723-1043; #wrong 444_20150702-1531; 
#======= psi 0.2 2D ======== -------------------------> redone, see above
# LHC14a1awithEta=445_20150702-1531;  #146, 147eta, 149eta
# LHC14a1bwithEta=447_20150702-1625; 
# LHC14a1awithPi0=446_20150702-1532;  #148, 147pi, 149pi
# LHC14a1bwithPi0=448_20150702-1626; 

#================================================================
# LHC11hData=185_20150630-1307;    #154, 158   
# LHC11hDataPhi=186_20150630-1307; #156, 160
#======= alpha meson 0.75 ========
# LHC14a1awithEta=449_20150702-1626;  #154, 155eta, 157eta
# LHC14a1bwithEta=451_20150702-1627; 
# LHC14a1awithPi0=450_20150702-1626;  #156, 155pi,  157pi
# LHC14a1bwithPi0=452_20150702-1627; 
#======= alpha meson 1.0 ========
# LHC14a1awithEta=453_20150702-1627;  #158, 159eta, 161eta
# LHC14a1bwithEta=455_20150702-1628; 
# LHC14a1awithPi0=454_20150702-1627;  #160, 159pi, 161pi
# LHC14a1bwithPi0=456_20150702-1628; 

#================================================================
#=================
		
		
standardData=313;
standardPhiData=315;
standard=313;
standardPhi=315; 
added=314;
addedPhi=316;

# standardData=226;
# standardPhiData=228;
# standard=226;
# standardPhi=228; 
# added=227;
# addedPhi=229;

		
###################################### LHC14a1a
echo "For traingconfig data:" $standardData
cat $OUTPUTDIR/GA_PbPb_AOD-$LHC11hData/merge_runlist_7/CutSelection_LHC11h_$standardData.log

hadd -f $OUTPUTDIRmerged/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1a_$standard.root $OUTPUTDIR/GA_PbPb_MC_AOD-$LHC14a1awithEta/merge_runlist_6/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1a_$standard.root $OUTPUTDIR/GA_PbPb_MC_AOD-$LHC14a1awithPi0/merge_runlist_7/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1a_$standardPhi.root

hadd -f $OUTPUTDIRmerged/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1a_Eta_$added.root $OUTPUTDIR/GA_PbPb_MC_AOD-$LHC14a1awithEta/merge_runlist_6/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1a_$added.root $OUTPUTDIR/GA_PbPb_MC_AOD-$LHC14a1awithEta/merge_runlist_7/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1a_$addedPhi.root

hadd -f $OUTPUTDIRmerged/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1a_Pi0_$added.root $OUTPUTDIR/GA_PbPb_MC_AOD-$LHC14a1awithPi0/merge_runlist_6/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1a_$added.root $OUTPUTDIR/GA_PbPb_MC_AOD-$LHC14a1awithPi0/merge_runlist_7/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1a_$addedPhi.root


####################################### LHC14a1b
hadd -f $OUTPUTDIRmerged/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1b_$standard.root $OUTPUTDIR/GA_PbPb_MC_AOD-$LHC14a1bwithEta/merge_runlist_6/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1b_$standard.root $OUTPUTDIR/GA_PbPb_MC_AOD-$LHC14a1bwithPi0/merge_runlist_7/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1b_$standardPhi.root

hadd -f $OUTPUTDIRmerged/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1b_Eta_$added.root $OUTPUTDIR/GA_PbPb_MC_AOD-$LHC14a1bwithEta/merge_runlist_6/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1b_$added.root $OUTPUTDIR/GA_PbPb_MC_AOD-$LHC14a1bwithEta/merge_runlist_7/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1b_$addedPhi.root

hadd -f $OUTPUTDIRmerged/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1b_Pi0_$added.root $OUTPUTDIR/GA_PbPb_MC_AOD-$LHC14a1bwithPi0/merge_runlist_6/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1b_$added.root $OUTPUTDIR/GA_PbPb_MC_AOD-$LHC14a1bwithPi0/merge_runlist_7/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1b_$addedPhi.root


##########################################   DATA   ###########################################   

hadd -f $OUTPUTDIRmerged/GammaConvV1_GA_PbPb_LHC11h-pass2_$standardData.root $OUTPUTDIR/GA_PbPb_AOD-$LHC11hData/merge_runlist_7/GammaConvV1_GA_PbPb_LHC11h-pass2_$standardData.root $OUTPUTDIR/GA_PbPb_AOD-$LHC11hDataPhi/merge/GammaConvV1_GA_PbPb_LHC11h-pass2_$standardPhiData.root

