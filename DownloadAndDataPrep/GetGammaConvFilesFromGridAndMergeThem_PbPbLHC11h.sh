# /bin/bash

# copies files from grid
# creates directory
# changes internal structure
# merges files according to the pPb needs

BASEDIR=/home/admin1/leardini/GridOutput/PbPb
mkdir -p $BASEDIR
TRAINDIR=Legotrain-vAN-20160623-1;


#photon cut studies
LHC11hData=332_20160706-1144; #---> list 7
LHC11hDataWithPhi=332_20160706-1144; #---> list 8

LHC14a1a=800_20160706-1427;
LHC14a1b=801_20160706-1541;
LHC14a1aWithPhi=800_20160706-1427;
LHC14a1bWithPhi=801_20160706-1541;


# #std on ESD with tag of 22nd may
# LHC11hData=238_20160628-1326; #---> list 7
# LHC11hDataWithPhi=238_20160628-1326; #---> list 8
#
# LHC14a1a=290_20160624-1346;
# LHC14a1b=289_20160624-1035;
# LHC14a1aWithPhi=290_20160624-1346;
# LHC14a1bWithPhi=289_20160624-1035;


############################################
# # #std on ESD with tag of 22nd may
# LHC11hData=228_20160522-1730; #---> list 7
# LHC11hDataWithPhi=228_20160522-1730; #---> list 8
# #
# LHC14a1a=278_20160522-1813;
# LHC14a1b=280_20160522-1731;
# LHC14a1aWithPhi=279_20160522-1801;
# LHC14a1bWithPhi=281_20160522-1730;
############################################

#Multiplicity ESDs################
# LHC11hData=223_20160519-0949;
# LHC11hData=227_20160520-1403;
# LHC11hData=229_20160523-1112;
##################################

# LHC11hData=301_20160425-1621; # 74 (spT 0.075), 78 (spT 0.1), 82 (TPCcls 0.7)
# LHC11hDataWithPhi=302_20160425-1622; # 76 (spT 0.075), 80 (spT 0.1), 84 (TPCcls 0.7)

# # 74, 75, 76, 77 (spT 0.075)
# LHC14a1a=707_20160425-1625;
# LHC14a1b=709_20160425-1626;
# LHC14a1aWithPhi=708_20160425-1625;
# LHC14a1bWithPhi=710_20160425-1626;

# # 78, 79, 80, 81 (spT 0.1)
# LHC14a1a=711_20160425-1911;
# LHC14a1b=719_20160429-1049;
# LHC14a1aWithPhi=712_20160425-1911;
# LHC14a1bWithPhi=714_20160425-1913;

# # 82, 83, 84, 85 (TPCcls 0.7)
# LHC14a1a=715_20160426-0942;
# LHC14a1b=717_20160426-0943;
# LHC14a1aWithPhi=716_20160426-0943;
# LHC14a1bWithPhi=718_20160426-0943;


# LHC11hData=303_20160429-1052; # 86 (TPCcls 0.35), 90 (edEdx-4,5), 94 (edEdx-2.5,4)
# LHC11hDataWithPhi=304_20160429-1052; # 88 (TPCcls 0.35), 92 (edEdx-4,5), 96 (edEdx-2.5,4)

# # 86, 87, 88, 89 (TPCcls 0.35)
# LHC14a1a=720_20160429-1104;
# LHC14a1b=722_20160429-1105;
# LHC14a1aWithPhi=721_20160429-1105;
# LHC14a1bWithPhi=737_20160503-1142;

# # 90, 91, 92, 93 (edEdx-4,5)
# LHC14a1a=723_20160429-1106;
# LHC14a1b=725_20160429-1107;
# LHC14a1aWithPhi=724_20160429-1106;
# LHC14a1bWithPhi=726_20160429-1107;

# # 94, 95, 96, 97 (edEdx-2.5,4)
# LHC14a1a=727_20160429-1108;
# LHC14a1b=729_20160429-1108;
# LHC14a1aWithPhi=728_20160429-1108;
# LHC14a1bWithPhi=730_20160429-1109;


# LHC11hData=310_20160502-1751; # 98 (pdEdx2s1s), 102 (pdEdx3s), 106 (pdEdx3s1sp)
# LHC11hDataWithPhi=311_20160502-1752; # 100 (pdEdx2s1s), 104 (pdEdx3s), 108 (pdEdx3s1sp)

# # 98, 99, 100, 101 (pdEdx2s1s) --> 98 is in 784
# LHC14a1a=738_20160511-1233;
# LHC14a1b=740_20160511-1234;
# LHC14a1aWithPhi=739_20160511-1234;
# LHC14a1bWithPhi=741_20160511-1235;

# # 102, 103, 104, 105 (pdEdx3s)
# LHC14a1a=742_20160511-1240;
# LHC14a1b=744_20160511-1242;
# LHC14a1aWithPhi=743_20160511-1241;
# LHC14a1bWithPhi=745_20160511-1242;

# # 106, 107, 108, 109 (pdEdx3s1sp)
# LHC14a1a=746_20160511-1244;
# LHC14a1b=748_20160511-1246;
# LHC14a1aWithPhi=747_20160511-1245;
# LHC14a1bWithPhi=749_20160511-1247;


# LHC11hData=312_20160502-1819; # 118 (qt0.03), 122 (qt0.06), 126 (chi50psi0.1)
# LHC11hDataWithPhi=313_20160502-1854; # 120 (qt0.03), 124 (qt0.06), 128 (chi50psi0.1)

# # 118, 119, 120, 121 (qt0.03)
# LHC14a1a=750_20160517-1015;
# LHC14a1b=752_20160517-1015;
# LHC14a1aWithPhi=751_20160511-1250;
# LHC14a1bWithPhi=753_20160511-1251;

# # 122, 123, 124, 125 (qt0.06)
# LHC14a1a=754_20160511-1253;
# LHC14a1b=756_20160511-1254;
# LHC14a1aWithPhi=755_20160511-1253;
# LHC14a1bWithPhi=757_20160511-1254;

# # 126, 127, 128, 129 (chi50psi0.1)
# LHC14a1a=758_20160511-1255;
# LHC14a1b=760_20160511-1258;
# LHC14a1aWithPhi=759_20160511-1256;
# LHC14a1bWithPhi=761_20160511-1258;



# LHC11hData=314_20160502-1819; # 130 (chi20psi0.1), 134 (chi30psi0.05), 138 (chi30psi0.2)
# LHC11hDataWithPhi=315_20160502-1820; # 132 (chi20psi0.1), 136 (chi30psi0.05), 140 (chi30psi0.2)

# # 130, 131, 132, 133 (chi20psi0.1)
# LHC14a1a=762_20160511-1411;
# LHC14a1b=764_20160511-1412;
# LHC14a1aWithPhi=763_20160511-1412;
# LHC14a1bWithPhi=786_20160523-1541;

# # 134, 135, 136, 137 (chi30psi0.05)
# LHC14a1a=766_20160516-1841;
# LHC14a1b=768_20160516-1843;
# LHC14a1aWithPhi=767_20160516-1842;
# LHC14a1bWithPhi=769_20160516-1843;

# # 138, 139, 140, 141 (chi30psi0.2)
# LHC14a1a=770_20160516-1844;
# LHC14a1b=772_20160516-1846;
# LHC14a1aWithPhi=771_20160516-1845;
# LHC14a1bWithPhi=773_20160516-1846;


# LHC11hData=316_20160503-1050; # 146 (alpha0.75), 150 (alpha0.1)
# LHC11hDataWithPhi=317_20160503-1050; # 148 (alpha0.75), 152 (alpha0.1)

# # 146, 147, 148, 149 (alpha0.75)
# LHC14a1a=774_20160516-1847;
# LHC14a1b=776_20160516-1848;
# LHC14a1aWithPhi=775_20160516-1848;
# LHC14a1bWithPhi=777_20160516-1849;

# # 150, 151, 152, 153 (alpha0.1)
# LHC14a1a=778_20160516-1850;
# LHC14a1b=780_20160516-1851;
# LHC14a1aWithPhi=779_20160516-1850;
# LHC14a1bWithPhi=781_20160516-1851;


# LHC11hDataWithPhi=318_20160502-1821; # 154 (phi2-4), 156 (phi2.4-3.6)

# # 154, 155 (phi2-4) and 156, 157 (phi2.4-3.6)
# LHC14a1a=782_20160516-1852;
# LHC14a1b=784_20160527-1402;
# LHC14a1aWithPhi=785_20160528-1443;
# LHC14a1bWithPhi=783_20160516-1853;


############################################
#std RP after TOF fix (+ AOD check)
# LHC11hData=306_20160429-1406;
# LHC11hDataWithPhi=307_20160429-1328;
#ESD check
# LHC11hData=213_20160429-1124;

#std
# LHC14a1a=731_20160429-1109;
# LHC14a1b=733_20160429-1110;
# LHC14a1aWithPhi=732_20160429-1518;
# LHC14a1bWithPhi=734_20160429-1518;
#AOD check
# LHC14a1a=735_20160429-1517;
# LHC14a1b=736_20160429-1517;
#ESD check
# LHC14a1a=269_20160429-1802;
# LHC14a1b=270_20160429-1122;

#DCA analysis
# LHC11hData=210_20160422-1702;
# LHC14a1b=268_20160429-1041;


OUTPUTDIR=$BASEDIR/$TRAINDIR
if [ $1 = "AODdata" ]; then
   TRAINPATHData=GA_PbPb_AOD
#     OUTPUTDIR_LHC11h=$BASEDIR/$TRAINDIR/$TRAINPATHData-$LHC11hData/merge_runlist_7
#     mkdir -p $OUTPUTDIR_LHC11h
#     alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHData/$LHC11hData/merge_runlist_7/Gam* file:$OUTPUTDIR_LHC11h/
#
#     OUTPUTDIR_LHC11hWithPhi=$BASEDIR/$TRAINDIR/$TRAINPATHData-$LHC11hDataWithPhi/merge
#     mkdir -p $OUTPUTDIR_LHC11hWithPhi
#     alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHData/$LHC11hDataWithPhi/merge/Gam* file:$OUTPUTDIR_LHC11hWithPhi/

else
   TRAINPATHData=GA_PbPb
    OUTPUTDIR_LHC11h=$BASEDIR/$TRAINDIR/$TRAINPATHData-$LHC11hData/merge_runlist_7/Stage_1/
    mkdir -p $OUTPUTDIR_LHC11h
    alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHData/$LHC11hData/merge_runlist_7/Stage_1/* file:$OUTPUTDIR_LHC11h/

    OUTPUTDIR_LHC11hWithPhi=$BASEDIR/$TRAINDIR/$TRAINPATHData-$LHC11hDataWithPhi/merge_runlist_8/Stage_1/
    mkdir -p $OUTPUTDIR_LHC11hWithPhi
    alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHData/$LHC11hDataWithPhi/merge_runlist_8/Stage_1/* file:$OUTPUTDIR_LHC11hWithPhi/

fi

# # for the standard analysis:
# OUTPUTDIR_LHC11h=$BASEDIR/$TRAINDIR/$TRAINPATHData-$LHC11hData/merge_runlist_1
# mkdir -p $OUTPUTDIR_LHC11h
# alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHData/$LHC11hData/merge_runlist_1/Gam* file:$OUTPUTDIR_LHC11h/


# #for the DCA analysis (runwise or stageoutput wise):
# if [ $1 = "DCAdata" ]; then
#
#   OUTPUTDIR_LHC11h=$BASEDIR/$TRAINDIR/$TRAINPATHData-$LHC11hData/merge_runlist_1
# #   stageOutputs=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHData/$LHC11hData/merge_runlist_1/Stage_1`
# # #   for stageOutput in $stageOutputs; do
# # #     mkdir -p $OUTPUTDIR_LHC11h/Stage_1/$stageOutput
#   runs=`cat lhc11hforDCA.txt`
#   for run in $runs; do
#     mkdir -p $OUTPUTDIR_LHC11h/$run
#     alien_cp alien:/alice/data/2011/LHC11h_2/000$run/ESDs/pass2/PWGGA/$TRAINPATHData/$LHC11hData/Gamm* file:$OUTPUTDIR_LHC11h/$run
#   done;
# fi


if [ $2 = "AODmc" ]; then
   TRAINPATHMC=GA_PbPb_MC_AOD

#     OUTPUTDIR_LHC14a1a=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1a/merge_runlist_6
#     OUTPUTDIR_LHC14a1b=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1b/merge_runlist_6
#     OUTPUTDIR_LHC14a1aWithPhi=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1aWithPhi/merge_runlist_7
#     OUTPUTDIR_LHC14a1bWithPhi=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1bWithPhi/merge_runlist_7

#     OUTPUTDIR_LHC14a1a=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1a/merge_runlist_7
#     OUTPUTDIR_LHC14a1b=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1b/merge_runlist_7
#     OUTPUTDIR_LHC14a1aWithPhi=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1aWithPhi/merge_runlist_6
#     OUTPUTDIR_LHC14a1bWithPhi=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1bWithPhi/merge_runlist_6

#     mkdir -p $OUTPUTDIR_LHC14a1a
#     mkdir -p $OUTPUTDIR_LHC14a1b
#     mkdir -p $OUTPUTDIR_LHC14a1aWithPhi
#     mkdir -p $OUTPUTDIR_LHC14a1bWithPhi

#     alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1a/merge_runlist_6/Gam* file:$OUTPUTDIR_LHC14a1a/
#     alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1b/merge_runlist_6/Gam* file:$OUTPUTDIR_LHC14a1b/
#     alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1aWithPhi/merge_runlist_7/Gam* file:$OUTPUTDIR_LHC14a1aWithPhi/
#     alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1bWithPhi/merge_runlist_7/Gam* file:$OUTPUTDIR_LHC14a1bWithPhi/

#     alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1a/merge_runlist_7/Gam* file:$OUTPUTDIR_LHC14a1a/
#     alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1b/merge_runlist_7/Gam* file:$OUTPUTDIR_LHC14a1b/
#     alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1aWithPhi/merge_runlist_6/Gam* file:$OUTPUTDIR_LHC14a1aWithPhi/
#     alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1bWithPhi/merge_runlist_6/Gam* file:$OUTPUTDIR_LHC14a1bWithPhi/

else
   TRAINPATHMC=GA_PbPb_MC

    OUTPUTDIR_LHC14a1a=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1a/merge_runlist_7
    OUTPUTDIR_LHC14a1b=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1b/merge_runlist_7
    OUTPUTDIR_LHC14a1aWithPhi=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1aWithPhi/merge_runlist_8
    OUTPUTDIR_LHC14a1bWithPhi=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1bWithPhi/merge_runlist_8

#     OUTPUTDIR_LHC14a1a=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1a/merge_runlist_8
#     OUTPUTDIR_LHC14a1b=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1b/merge_runlist_8
#     OUTPUTDIR_LHC14a1aWithPhi=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1aWithPhi/merge_runlist_7
#     OUTPUTDIR_LHC14a1bWithPhi=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1bWithPhi/merge_runlist_7

    mkdir -p $OUTPUTDIR_LHC14a1a
    mkdir -p $OUTPUTDIR_LHC14a1b
    mkdir -p $OUTPUTDIR_LHC14a1aWithPhi
    mkdir -p $OUTPUTDIR_LHC14a1bWithPhi

    alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1a/merge_runlist_7/Gam* file:$OUTPUTDIR_LHC14a1a/
    alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1b/merge_runlist_7/Gam* file:$OUTPUTDIR_LHC14a1b/
    alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1aWithPhi/merge_runlist_8/Gam* file:$OUTPUTDIR_LHC14a1aWithPhi/
    alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1bWithPhi/merge_runlist_8/Gam* file:$OUTPUTDIR_LHC14a1bWithPhi/

#     alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1a/merge_runlist_8/Gam* file:$OUTPUTDIR_LHC14a1a/
#     alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1b/merge_runlist_8/Gam* file:$OUTPUTDIR_LHC14a1b/
#     alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1aWithPhi/merge_runlist_7/Gam* file:$OUTPUTDIR_LHC14a1aWithPhi/
#     alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1bWithPhi/merge_runlist_7/Gam* file:$OUTPUTDIR_LHC14a1bWithPhi/

fi

# # OUTPUTDIR_LHC14a1a=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1a/merge_runlist_1
# # alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1a/merge_runlist_1/Gam* file:$OUTPUTDIR_LHC14a1a/
# # OUTPUTDIR_LHC14a1b=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1b/merge_runlist_1
# # alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1b/merge_runlist_1/Gam* file:$OUTPUTDIR_LHC14a1b/




# if [ $2 = "DCAmc" ]; then
#
#   OUTPUTDIR_LHC14a1b=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1b/merge_runlist_1
#   stageOutputs=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1b/merge_runlist_1/Stage_1`
#   for stageOutput in $stageOutputs; do
#     mkdir -p $OUTPUTDIR_LHC14a1b/Stage_1/$stageOutput
#     alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1b/merge_runlist_1/Stage_1/$stageOutput/Gamm* file:$OUTPUTDIR_LHC14a1b/Stage_1/$stageOutput/
#   done;
#   runs=`cat lhc11hforDCA.txt`
#   for run in $runs; do
#     mkdir -p $OUTPUTDIR_LHC14a1b/$run
#     alien_cp alien:/alice/sim/2014/LHC14a1b/$run/PWGGA/$TRAINPATHMC/$LHC14a1b/Gamm* file:$OUTPUTDIR_LHC14a1b/$run
#   done;

# fi



# # # OUTPUTDIR_LHC14a1c=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1c
# # # mkdir -p $OUTPUTDIR_LHC14a1c
# # # alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1c/merge/Gam* file:$OUTPUTDIR_LHC14a1c/
# # #


if [ $2 = "AODmc" ]; then
#################################### normal selection cut ######################################################
   ls $OUTPUTDIR_LHC14a1a/GammaConvV1_*.root > fileLHC14a1a.txt
   fileNumbers=`cat fileLHC14a1a.txt`
   for fileName in $fileNumbers; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 10| cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC14a1a/GammaConvV1_$number.root\"\,\"$OUTPUTDIR_LHC14a1a/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1a_$number.root\"\,\"GammaConvV1_$number\"\)
      root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR_LHC14a1a/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1a_$number.root\"\,\"$OUTPUTDIR_LHC14a1a/CutSelection_LHC14a1a_AOD_$number.log\"\)

   done;

   ls $OUTPUTDIR_LHC14a1b/GammaConvV1_*.root > fileLHC14a1b.txt
   fileNumbersb=`cat fileLHC14a1b.txt`
   for fileName in $fileNumbersb; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC14a1b/GammaConvV1_$number.root\"\,\"$OUTPUTDIR_LHC14a1b/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1b_$number.root\"\,\"GammaConvV1_$number\"\)
      root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR_LHC14a1b/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1b_$number.root\"\,\"$OUTPUTDIR_LHC14a1b/CutSelection_LHC14a1b_AOD_$number.log\"\)
   done;

#    ls $OUTPUTDIR_LHC14a1c/GammaConvV1_*.root > fileLHC14a1c.txt
#    fileNumbersb=`cat fileLHC14a1c.txt`
#    for fileName in $fileNumbersb; do
#       echo $fileName
#       number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
#       echo $number
#       root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC14a1c/GammaConvV1_$number.root\"\,\"$OUTPUTDIR_LHC14a1c/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1c_$number.root\"\,\"GammaConvV1_$number\"\)
#       root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR_LHC14a1c/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1c_$number.root\"\,\"$OUTPUTDIR_LHC14a1c/CutSelection_LHC14a1c_AOD_$number.log\"\)
#    done;


######################################## with phi cut ######################################################
   ls $OUTPUTDIR_LHC14a1aWithPhi/GammaConvV1_*.root > fileLHC14a1aWithPhi.txt
   fileNumbers=`cat fileLHC14a1aWithPhi.txt`
   for fileName in $fileNumbers; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 10| cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC14a1aWithPhi/GammaConvV1_$number.root\"\,\"$OUTPUTDIR_LHC14a1aWithPhi/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1a_$number.root\"\,\"GammaConvV1_$number\"\)
      root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR_LHC14a1aWithPhi/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1a_$number.root\"\,\"$OUTPUTDIR_LHC14a1aWithPhi/CutSelection_LHC14a1a_AOD_$number.log\"\)

   done;

   ls $OUTPUTDIR_LHC14a1bWithPhi/GammaConvV1_*.root > fileLHC14a1bWithPhi.txt
   fileNumbersb=`cat fileLHC14a1bWithPhi.txt`
   for fileName in $fileNumbersb; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC14a1bWithPhi/GammaConvV1_$number.root\"\,\"$OUTPUTDIR_LHC14a1bWithPhi/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1b_$number.root\"\,\"GammaConvV1_$number\"\)
      root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR_LHC14a1bWithPhi/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1b_$number.root\"\,\"$OUTPUTDIR_LHC14a1bWithPhi/CutSelection_LHC14a1b_AOD_$number.log\"\)
   done;


# elif [ $2 = "DCAmc" ]; then
#     echo $OUTPUTDIR_LHC14a1b
#     runs=`cat lhc11hforDCA.txt`
#     for run in $runs; do
#       ls $OUTPUTDIR_LHC14a1b/$run/GammaConvV1_*.root > fileLHC11hDCA.txt
#       fileNumbersDCAmc=`cat fileLHC11hDCA.txt`
#       for fileName in $fileNumbersDCAmc; do
#           echo $fileName
#           number=`echo $fileName  | cut -d "/" -f 11 | cut -d "_" -f 2 | cut -d "." -f1`
#           echo $number
# #           root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC14a1b/$run/GammaConvV1_$number.root\"\,\"$OUTPUTDIR_LHC14a1b/$run/GammaConvV1_GA_PbPb_MC_LHC14a1b_$number.root\"\,\"GammaConvV1_$number\"\)
# #           root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR_LHC14a1b/$run/GammaConvV1_GA_PbPb_MC_LHC14a1b_$number.root\"\,\"$OUTPUTDIR_LHC14a1b/$run/CutSelection_LHC14a1b_$number.log\"\)
#       done;
#     done;
#
#     counter=0;
#     number=40;
#     mkdir -p $OUTPUTDIR_LHC14a1b/merged/
#     echo $number
#     for run in $runs; do
#        echo "run number ---> " $run
#        mv $OUTPUTDIR_LHC14a1b/$run/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1b_$number.root $OUTPUTDIR_LHC14a1b/$run/GammaConvV1_GA_PbPb_MC_LHC14a1b_$number.root
#         if [ -f $OUTPUTDIR_LHC14a1b/$run/GammaConvV1_GA_PbPb_MC_LHC14a1b_$number.root ]; then
#           if [ $counter = 0 ]; then
#             cp $OUTPUTDIR_LHC14a1b/$run/GammaConvV1_GA_PbPb_MC_LHC14a1b_$number.root $OUTPUTDIR_LHC14a1b/intermediate.root
#             counter=$(($counter+1));
#             echo $counter;
#           else
#             hadd -f $OUTPUTDIR_LHC14a1b/GammaConvV1_GA_PbPb_MC_LHC14a1b_merged_$number.root $OUTPUTDIR_LHC14a1b/intermediate.root $OUTPUTDIR_LHC14a1b/$run/GammaConvV1_GA_PbPb_MC_LHC14a1b_$number.root
#             mv $OUTPUTDIR_LHC14a1b/GammaConvV1_GA_PbPb_MC_LHC14a1b_merged_$number.root $OUTPUTDIR_LHC14a1b/intermediate.root
#           fi
#         else
#           echo "file not there for run " $run
#         fi
#     done;
#     mv $OUTPUTDIR_LHC14a1b/intermediate.root  $OUTPUTDIR_LHC14a1b/merged/GammaConvV1_GA_PbPb_MC_LHC14a1b_FinalMerge_$number.root


else

   ls $OUTPUTDIR_LHC14a1a/GammaConvV1_*.root > fileLHC14a1a.txt
   fileNumbers=`cat fileLHC14a1a.txt`
   for fileName in $fileNumbers; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC14a1a/GammaConvV1_$number.root\"\,\"$OUTPUTDIR_LHC14a1a/GammaConvV1_GA_PbPb_MC_LHC14a1a_$number.root\"\,\"GammaConvV1_$number\"\)
      root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR_LHC14a1a/GammaConvV1_GA_PbPb_MC_LHC14a1a_$number.root\"\,\"$OUTPUTDIR_LHC14a1a/CutSelection_LHC14a1a_$number.log\"\)

   done;

   ls $OUTPUTDIR_LHC14a1b/GammaConvV1_*.root > fileLHC14a1b.txt
   fileNumbersb=`cat fileLHC14a1b.txt`
   for fileName in $fileNumbersb; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC14a1b/GammaConvV1_$number.root\"\,\"$OUTPUTDIR_LHC14a1b/GammaConvV1_GA_PbPb_MC_LHC14a1b_$number.root\"\,\"GammaConvV1_$number\"\)
      root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR_LHC14a1b/GammaConvV1_GA_PbPb_MC_LHC14a1b_$number.root\"\,\"$OUTPUTDIR_LHC14a1b/CutSelection_LHC14a1b_$number.log\"\)
   done;

# ######################################## with phi cut ######################################################
   ls $OUTPUTDIR_LHC14a1aWithPhi/GammaConvV1_*.root > fileLHC14a1aWithPhi.txt
   fileNumbers=`cat fileLHC14a1aWithPhi.txt`
   for fileName in $fileNumbers; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 10| cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC14a1aWithPhi/GammaConvV1_$number.root\"\,\"$OUTPUTDIR_LHC14a1aWithPhi/GammaConvV1_GA_PbPb_MC_LHC14a1a_$number.root\"\,\"GammaConvV1_$number\"\)
      root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR_LHC14a1aWithPhi/GammaConvV1_GA_PbPb_MC_LHC14a1a_$number.root\"\,\"$OUTPUTDIR_LHC14a1aWithPhi/CutSelection_LHC14a1a_$number.log\"\)

   done;

   ls $OUTPUTDIR_LHC14a1bWithPhi/GammaConvV1_*.root > fileLHC14a1bWithPhi.txt
   fileNumbersb=`cat fileLHC14a1bWithPhi.txt`
   for fileName in $fileNumbersb; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC14a1bWithPhi/GammaConvV1_$number.root\"\,\"$OUTPUTDIR_LHC14a1bWithPhi/GammaConvV1_GA_PbPb_MC_LHC14a1b_$number.root\"\,\"GammaConvV1_$number\"\)
      root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR_LHC14a1bWithPhi/GammaConvV1_GA_PbPb_MC_LHC14a1b_$number.root\"\,\"$OUTPUTDIR_LHC14a1bWithPhi/CutSelection_LHC14a1b_$number.log\"\)
   done;

fi



#
# if [ $1 = "AODdata" ]; then
# ################################## normal selection cut ######################################################
#
#    ls $OUTPUTDIR_LHC11h/GammaConvV1_*.root > fileLHC11h.txt
#    fileNumbersData=`cat fileLHC11h.txt`
#    for fileName in $fileNumbersData; do
#       echo $fileName
#       number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
#       echo $number
#       root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC11h/GammaConvV1_$number.root\"\,\"$OUTPUTDIR_LHC11h/GammaConvV1_GA_PbPb_LHC11h-pass2_$number.root\"\,\"GammaConvV1_$number\"\)
#       root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR_LHC11h/GammaConvV1_GA_PbPb_LHC11h-pass2_$number.root\"\,\"$OUTPUTDIR_LHC11h/CutSelection_LHC11h_$number.log\"\)
#    done;
#
# ######################################## with phi cut ######################################################
#    ls $OUTPUTDIR_LHC11hWithPhi/GammaConvV1_*.root > fileLHC11hWithPhi.txt
#    fileNumbersData=`cat fileLHC11hWithPhi.txt`
#    for fileName in $fileNumbersData; do
#       echo $fileName
#       number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
#       echo $number
#       root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC11hWithPhi/GammaConvV1_$number.root\"\,\"$OUTPUTDIR_LHC11hWithPhi/GammaConvV1_GA_PbPb_LHC11h-pass2_$number.root\"\,\"GammaConvV1_$number\"\)
#       root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR_LHC11hWithPhi/GammaConvV1_GA_PbPb_LHC11h-pass2_$number.root\"\,\"$OUTPUTDIR_LHC11hWithPhi/CutSelection_LHC11h_$number.log\"\)
#    done;
#
#
# # # elif [ $1 = "DCAdata" ]; then
# # #     runs=`cat lhc11hforDCA.txt`
# # #     for run in $runs; do
# # #       ls $OUTPUTDIR_LHC11h/$run/GammaConvV1_*.root > fileLHC11hDCA.txt
# # #       fileNumbersData=`cat fileLHC11hDCA.txt`
# # #       for fileName in $fileNumbersData; do
# # #           echo $fileName
# # #           number=`echo $fileName  | cut -d "/" -f 11 | cut -d "_" -f 2 | cut -d "." -f1`
# # #           echo $number
# # #           root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC11h/$run/GammaConvV1_$number.root\"\,\"$OUTPUTDIR_LHC11h/$run/GammaConvV1_GA_PbPb_LHC11h-pass2_$number.root\"\,\"GammaConvV1_$number\"\)
# # #           root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR_LHC11h/$run/GammaConvV1_GA_PbPb_LHC11h-pass2_$number.root\"\,\"$OUTPUTDIR_LHC11h/$run/CutSelection_LHC11h_$number.log\"\)
# # #       done;
# # #     done;
# # #     counter=0;
# # #     number=40;
# # #     echo $number
# # #     mkdir -p $OUTPUTDIR_LHC11h/merged
# # #     for run in $runs; do
# # #        echo "run number ---> " $run
# # #        if [ $counter = 0 ]; then
# # #           cp $OUTPUTDIR_LHC11h/$run/GammaConvV1_GA_PbPb_LHC11h-pass2_$number.root $OUTPUTDIR_LHC11h/intermediate.root
# # #           counter=$(($counter+1));
# # #           echo $counter;
# # #        else
# # #           hadd -f $OUTPUTDIR_LHC11h/GammaConvV1_GA_PbPb_LHC11h-pass2_merged_$number.root $OUTPUTDIR_LHC11h/intermediate.root $OUTPUTDIR_LHC11h/$run/GammaConvV1_GA_PbPb_LHC11h-pass2_$number.root
# # #           mv $OUTPUTDIR_LHC11h/GammaConvV1_GA_PbPb_LHC11h-pass2_merged_$number.root $OUTPUTDIR_LHC11h/intermediate.root
# # #        fi
# # #     done;
# # #     mv $OUTPUTDIR_LHC11h/intermediate.root  $OUTPUTDIR_LHC11h/merged/GammaConvV1_GA_PbPb_LHC11h-pass2_FinalMerge_$number.root
#
#
# else
#
#    ls $OUTPUTDIR_LHC11h/GammaConvV1_*.root > fileLHC11h.txt
#    fileNumbersData=`cat fileLHC11h.txt`
#    for fileName in $fileNumbersData; do
#       echo $fileName
#       number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
#       echo $number
#       root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC11h/GammaConvV1_$number.root\"\,\"$OUTPUTDIR_LHC11h/GammaConvV1_GA_PbPb_LHC11h-pass2_$number.root\"\,\"GammaConvV1_$number\"\)
#       root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR_LHC11h/GammaConvV1_GA_PbPb_LHC11h-pass2_$number.root\"\,\"$OUTPUTDIR_LHC11h/CutSelection_LHC11h_$number.log\"\)
#    done;
#
# ######################################## with phi cut ######################################################
#    ls $OUTPUTDIR_LHC11hWithPhi/GammaConvV1_*.root > fileLHC11hWithPhi.txt
#    fileNumbersData=`cat fileLHC11hWithPhi.txt`
#    for fileName in $fileNumbersData; do
#       echo $fileName
#       number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
#       echo $number
#       root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC11hWithPhi/GammaConvV1_$number.root\"\,\"$OUTPUTDIR_LHC11hWithPhi/GammaConvV1_GA_PbPb_LHC11h-pass2_$number.root\"\,\"GammaConvV1_$number\"\)
#       root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR_LHC11hWithPhi/GammaConvV1_GA_PbPb_LHC11h-pass2_$number.root\"\,\"$OUTPUTDIR_LHC11hWithPhi/CutSelection_LHC11h_$number.log\"\)
#    done;
#
# fi


