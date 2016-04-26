#! /bin/bash

# # Yield Cutstudies
#echo -e "" > CutSelection.log
#bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersPi0MBYield


# # Eta Cutstudies # #
# echo -e "6010001_00200009247602008250400000_0152501500000000\n6010001_03200009247602008250400000_0152301500000000\n6010001_04200009247602008250400000_0152201500000000" > CutSelection.log
# cat CutSelection.log
# echo -e "Eta0005\nNo\n0\nYes\nPbPb_2.76TeV\nNo\nNo\nYes" > answersFileSysErr
# bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersFileSysErr


# # SinglePt Cutstudies # #
 echo -e "60100013_00200009247602008250404000_0152501500000000\n6010001_00200049247602008250400000_0152501500000000\n6010001_00200019247602008250400000_0152501500000000" > CutSelection.log
 echo -e "SinglePt0005\nNo\n0\nYes\nPbPb_2.76TeV\nNo\nNo\nYes" > answersFileSysErr
 bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersFileSysErr


# # TPC Cluster Cutstudies # #
 echo -e "60100013_00200009247602008250404000_0152501500000000\n6010001_00200006247602008250400000_0152501500000000\n6010001_00200008247602008250400000_0152501500000000" > CutSelection.log
 echo -e "TPCCluster0005\nNo\n0\nYes\nPbPb_2.76TeV\nNo\nNo\nYes" > answersFileSysErr
 bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersFileSysErr


# # dEdxE Cutstudies # # 
 echo -e "60100013_00200009247602008250404000_0152501500000000\n6010001_00200009347602008250400000_0152501500000000\n6010001_00200009647602008250400000_0152501500000000" > CutSelection.log
 echo -e "dEdxE0005\nNo\n0\nYes\nPbPb_2.76TeV\nNo\nNo\nYes" > answersFileSysErr
 bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersFileSysErr

 
# # dEdxPi Cutstudies # #
 echo -e "60100013_00200009247602008250404000_0152501500000000\n6010001_00200009287602008250400000_0152501500000000\n6010001_00200009237002008250400000_0152501500000000\n6010001_00200009245402008250400000_0152501500000000" > CutSelection.log
 echo -e "dEdxPi0005\nNo\n0\nYes\nPbPb_2.76TeV\nNo\nNo\nYes" > answersFileSysErr
 bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersFileSysErr

 
# pdEdxPi Cutstudies # #
#  echo -e "6010001_00200009247602008250400000_0152501500000000\n6010001_00200009245402008250400000_0152501500000000" > CutSelection.log
#  echo -e "pdEdxPi0005\nNo\n0\nYes\nPbPb_2.76TeV\nNo\nNo\nYes" > answersFileSysErr
# bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersFileSysErr
 

 # cosPA Cutstudies # #
#  echo -e "6010001_00200009247602008250400000_0152501500000000\n6010001_00200009247602008250000000_0152501500000000" > CutSelection.log
#  echo -e "CosPoint0005\nNo\n0\nYes\nPbPb_2.76TeV\nNo\nNo\nYes" > answersFileSysErr
# bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersFileSysErr

 
# TOF Cutstudies # #
#  echo -e "6010001_00200009247602008250400000_0152501500000000\n6010001_00200009247603008250400000_0152501500000000\n6010001_00200009247604008250400000_0152501500000000" > CutSelection.log
#  echo -e "TOF0005\nNo\n0\nYes\nPbPb_2.76TeV\nNo\nNo\nYes" > answersFileSysErr
# bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersFileSysErr
 

# # Qt Cutstudies # #
 echo -e "60100013_00200009247602008250404000_0152501500000000\n6010001_00200009247602009250400000_0152501500000000\n6010001_00200009247602002250400000_0152501500000000" > CutSelection.log
 echo -e "Qt0005\nNo\n0\nYes\nPbPb_2.76TeV\nNo\nNo\nYes" > answersFileSysErr
 bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersFileSysErr

 
# # Chi2 & PsiPair Cutstudies # #
 echo -e "60100013_00200009247602008250404000_0152501500000000\n6010001_00200009247602008150400000_0152501500000000\n6010001_00200009247602008850400000_0152501500000000\n6010001_00200009247602008280400000_0152501500000000\n6010001_00200009247602008260400000_0152501500000000" > CutSelection.log
 echo -e "Chi20005\nNo\n0\nYes\nPbPb_2.76TeV\nNo\nNo\nYes" > answersFileSysErr
 bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersFileSysErr

 
# Alpha Cutstudies # #
echo -e "60100013_00200009247602008250404000_0152501500000000\n6010001_00200009247602008250400000_0152505500000000\n6010001_00200009247602008250400000_0152503500000000" > CutSelection.log
echo -e "Alpha0005\nNo\n0\nYes\nPbPb_2.76TeV\nNo\nNo\nYes" > answersFileSysErr
bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersFileSysErr


# Phi Cutstudies # #
echo -e "60100013_00200009247602008250404000_0152501500000000\n6010001_00215509247602008250400000_0152501500000000\n6010001_00217709247602008250400000_0152501500000000" > CutSelection.log
echo -e "ConvPhi0005\nNo\n0\nYes\nPbPb_2.76TeV\nNo\nNo\nYes" > answersFileSysErr
bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersFileSysErr

# # # BG Cutstudies # #
# echo -e "6010001_00200009247602008250400000_0152501500000000\n0000011_00200009227302008250400000_0252103500000000" > CutSelection.log
# echo -e "BG0005\nNo\n0\nYes\nPbPb_2.76TeV\nNo\nNo\nYes" > answersFileSysErr
# bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersFileSysErr



#! /bin/bash

# # Yield Cutstudies
#echo -e "" > CutSelection.log
#bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersPi0MBYield


# # Eta Cutstudies # #
# echo -e "5010001_00200009247602008250400000_0152501500000000\n5010001_03200009247602008250400000_0152301500000000\n5010001_04200009247602008250400000_0152201500000000" > CutSelection.log
# cat CutSelection.log
# echo -e "Eta0010\nNo\n0\nYes\nPbPb_2.76TeV\nNo\nNo\nYes" > answersFileSysErr
# bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersFileSysErr


# # SinglePt Cutstudies # #
 echo -e "50100013_00200009247602008250404000_0152501500000000\n5010001_00200049247602008250400000_0152501500000000\n5010001_00200019247602008250400000_0152501500000000" > CutSelection.log
 echo -e "SinglePt0010\nNo\n0\nYes\nPbPb_2.76TeV\nNo\nNo\nYes" > answersFileSysErr
 bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersFileSysErr


# # TPC Cluster Cutstudies # #
 echo -e "50100013_00200009247602008250404000_0152501500000000\n5010001_00200006247602008250400000_0152501500000000\n5010001_00200008247602008250400000_0152501500000000" > CutSelection.log
 echo -e "TPCCluster0010\nNo\n0\nYes\nPbPb_2.76TeV\nNo\nNo\nYes" > answersFileSysErr
 bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersFileSysErr


# # dEdxE Cutstudies # # 
 echo -e "50100013_00200009247602008250404000_0152501500000000\n5010001_00200009347602008250400000_0152501500000000\n5010001_00200009647602008250400000_0152501500000000" > CutSelection.log
 echo -e "dEdxE0010\nNo\n0\nYes\nPbPb_2.76TeV\nNo\nNo\nYes" > answersFileSysErr
 bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersFileSysErr

 
# # dEdxPi Cutstudies # #
 echo -e "50100013_00200009247602008250404000_0152501500000000\n5010001_00200009287602008250400000_0152501500000000\n5010001_00200009237002008250400000_0152501500000000\n5010001_00200009245402008250400000_0152501500000000" > CutSelection.log
 echo -e "dEdxPi0010\nNo\n0\nYes\nPbPb_2.76TeV\nNo\nNo\nYes" > answersFileSysErr
 bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersFileSysErr

 
# pdEdxPi Cutstudies # #
#  echo -e "5010001_00200009247602008250400000_0152501500000000\n5010001_00200009245402008250400000_0152501500000000" > CutSelection.log
#  echo -e "pdEdxPi0010\nNo\n0\nYes\nPbPb_2.76TeV\nNo\nNo\nYes" > answersFileSysErr
# bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersFileSysErr
 

 # cosPA Cutstudies # #
#  echo -e "5010001_00200009247602008250400000_0152501500000000\n5010001_00200009247602008250000000_0152501500000000" > CutSelection.log
#  echo -e "CosPoint0010\nNo\n0\nYes\nPbPb_2.76TeV\nNo\nNo\nYes" > answersFileSysErr
# bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersFileSysErr

 
# TOF Cutstudies # #
#  echo -e "5010001_00200009247602008250400000_0152501500000000\n5010001_00200009247603008250400000_0152501500000000\n5010001_00200009247604008250400000_0152501500000000" > CutSelection.log
#  echo -e "TOF0010\nNo\n0\nYes\nPbPb_2.76TeV\nNo\nNo\nYes" > answersFileSysErr
# bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersFileSysErr
 

# # Qt Cutstudies # #
 echo -e "50100013_00200009247602008250404000_0152501500000000\n5010001_00200009247602009250400000_0152501500000000\n5010001_00200009247602002250400000_0152501500000000" > CutSelection.log
 echo -e "Qt0010\nNo\n0\nYes\nPbPb_2.76TeV\nNo\nNo\nYes" > answersFileSysErr
 bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersFileSysErr

 
# # Chi2 & PsiPair Cutstudies # #
 echo -e "50100013_00200009247602008250404000_0152501500000000\n5010001_00200009247602008150400000_0152501500000000\n5010001_00200009247602008850400000_0152501500000000\n5010001_00200009247602008280400000_0152501500000000\n5010001_00200009247602008260400000_0152501500000000" > CutSelection.log
 echo -e "Chi20010\nNo\n0\nYes\nPbPb_2.76TeV\nNo\nNo\nYes" > answersFileSysErr
 bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersFileSysErr

 
# Alpha Cutstudies # #
echo -e "50100013_00200009247602008250404000_0152501500000000\n5010001_00200009247602008250400000_0152505500000000\n5010001_00200009247602008250400000_0152503500000000" > CutSelection.log
echo -e "Alpha0010\nNo\n0\nYes\nPbPb_2.76TeV\nNo\nNo\nYes" > answersFileSysErr
bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersFileSysErr


# Phi Cutstudies # #
echo -e "50100013_00200009247602008250404000_0152501500000000\n5010001_00215509247602008250400000_0152501500000000\n5010001_00217709247602008250400000_0152501500000000" > CutSelection.log
echo -e "ConvPhi0010\nNo\n0\nYes\nPbPb_2.76TeV\nNo\nNo\nYes" > answersFileSysErr
bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersFileSysErr

# # # BG Cutstudies # #
# echo -e "5010001_00200009247602008250400000_0152501500000000\n0000011_00200009227302008250400000_0252103500000000" > CutSelection.log
# echo -e "BG0010\nNo\n0\nYes\nPbPb_2.76TeV\nNo\nNo\nYes" > answersFileSysErr
# bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersFileSysErr


#! /bin/bash

# # Yield Cutstudies
#echo -e "" > CutSelection.log
#bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersPi0MBYield


# # Eta Cutstudies # #
# echo -e "6120001_00200009247602008250400000_0152501500000000\n6120001_03200009247602008250400000_0152301500000000\n6120001_04200009247602008250400000_0152201500000000" > CutSelection.log
# cat CutSelection.log
# echo -e "Eta0510\nNo\n0\nYes\nPbPb_2.76TeV\nNo\nNo\nYes" > answersFileSysErr
# bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersFileSysErr


# # SinglePt Cutstudies # #
 echo -e "61200013_00200009247602008250404000_0152501500000000\n6120001_00200049247602008250400000_0152501500000000\n6120001_00200019247602008250400000_0152501500000000" > CutSelection.log
 echo -e "SinglePt0510\nNo\n0\nYes\nPbPb_2.76TeV\nNo\nNo\nYes" > answersFileSysErr
 bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersFileSysErr


# # TPC Cluster Cutstudies # #
 echo -e "61200013_00200009247602008250404000_0152501500000000\n6120001_00200006247602008250400000_0152501500000000\n6120001_00200008247602008250400000_0152501500000000" > CutSelection.log
 echo -e "TPCCluster0510\nNo\n0\nYes\nPbPb_2.76TeV\nNo\nNo\nYes" > answersFileSysErr
 bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersFileSysErr


# # dEdxE Cutstudies # # 
 echo -e "61200013_00200009247602008250404000_0152501500000000\n6120001_00200009347602008250400000_0152501500000000\n6120001_00200009647602008250400000_0152501500000000" > CutSelection.log
 echo -e "dEdxE0510\nNo\n0\nYes\nPbPb_2.76TeV\nNo\nNo\nYes" > answersFileSysErr
 bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersFileSysErr

 
# # dEdxPi Cutstudies # #
 echo -e "61200013_00200009247602008250404000_0152501500000000\n6120001_00200009287602008250400000_0152501500000000\n6120001_00200009237002008250400000_0152501500000000\n6120001_00200009245402008250400000_0152501500000000" > CutSelection.log
 echo -e "dEdxPi0510\nNo\n0\nYes\nPbPb_2.76TeV\nNo\nNo\nYes" > answersFileSysErr
 bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersFileSysErr

 
# pdEdxPi Cutstudies # #
#  echo -e "6120001_00200009247602008250400000_0152501500000000\n6120001_00200009245402008250400000_0152501500000000" > CutSelection.log
#  echo -e "pdEdxPi0510\nNo\n0\nYes\nPbPb_2.76TeV\nNo\nNo\nYes" > answersFileSysErr
# bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersFileSysErr
 

 # cosPA Cutstudies # #
#  echo -e "6120001_00200009247602008250400000_0152501500000000\n6120001_00200009247602008250000000_0152501500000000" > CutSelection.log
#  echo -e "CosPoint0510\nNo\n0\nYes\nPbPb_2.76TeV\nNo\nNo\nYes" > answersFileSysErr
# bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersFileSysErr

 
# TOF Cutstudies # #
#  echo -e "6120001_00200009247602008250400000_0152501500000000\n6120001_00200009247603008250400000_0152501500000000\n6120001_00200009247604008250400000_0152501500000000" > CutSelection.log
#  echo -e "TOF0510\nNo\n0\nYes\nPbPb_2.76TeV\nNo\nNo\nYes" > answersFileSysErr
# bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersFileSysErr
 

# # Qt Cutstudies # #
 echo -e "61200013_00200009247602008250404000_0152501500000000\n6120001_00200009247602009250400000_0152501500000000\n6120001_00200009247602002250400000_0152501500000000" > CutSelection.log
 echo -e "Qt0510\nNo\n0\nYes\nPbPb_2.76TeV\nNo\nNo\nYes" > answersFileSysErr
 bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersFileSysErr

 
# # Chi2 & PsiPair Cutstudies # #
 echo -e "61200013_00200009247602008250404000_0152501500000000\n6120001_00200009247602008150400000_0152501500000000\n6120001_00200009247602008850400000_0152501500000000\n6120001_00200009247602008280400000_0152501500000000\n6120001_00200009247602008260400000_0152501500000000" > CutSelection.log
 echo -e "Chi20510\nNo\n0\nYes\nPbPb_2.76TeV\nNo\nNo\nYes" > answersFileSysErr
 bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersFileSysErr

 
# Alpha Cutstudies # #
echo -e "61200013_00200009247602008250404000_0152501500000000\n6120001_00200009247602008250400000_0152505500000000\n6120001_00200009247602008250400000_0152503500000000" > CutSelection.log
echo -e "Alpha0510\nNo\n0\nYes\nPbPb_2.76TeV\nNo\nNo\nYes" > answersFileSysErr
bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersFileSysErr


# Phi Cutstudies # #
echo -e "61200013_00200009247602008250404000_0152501500000000\n6120001_00215509247602008250400000_0152501500000000\n6120001_00217709247602008250400000_0152501500000000" > CutSelection.log
echo -e "ConvPhi0510\nNo\n0\nYes\nPbPb_2.76TeV\nNo\nNo\nYes" > answersFileSysErr
bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersFileSysErr

# # # BG Cutstudies # #
# echo -e "6120001_00200009247602008250400000_0152501500000000\n0000011_00200009227302008250400000_0252103500000000" > CutSelection.log
# echo -e "BG0510\nNo\n0\nYes\nPbPb_2.76TeV\nNo\nNo\nYes" > answersFileSysErr
# bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersFileSysErr


#! /bin/bash

# # Yield Cutstudies
#echo -e "" > CutSelection.log
#bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersPi0MBYield


# # Eta Cutstudies # #
# echo -e "5240001_00200009247602008250400000_0152501500000000\n5240001_03200009247602008250400000_0152301500000000\n5240001_04200009247602008250400000_0152201500000000" > CutSelection.log
# cat CutSelection.log
# echo -e "Eta2040\nNo\n0\nYes\nPbPb_2.76TeV\nNo\nNo\nYes" > answersFileSysErr
# bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersFileSysErr


# # SinglePt Cutstudies # #
 echo -e "52400013_00200009247602008250404000_0152501500000000\n5240001_00200049247602008250400000_0152501500000000\n5240001_00200019247602008250400000_0152501500000000" > CutSelection.log
 echo -e "SinglePt2040\nNo\n0\nYes\nPbPb_2.76TeV\nNo\nNo\nYes" > answersFileSysErr
 bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersFileSysErr


# # TPC Cluster Cutstudies # #
 echo -e "52400013_00200009247602008250404000_0152501500000000\n5240001_00200006247602008250400000_0152501500000000\n5240001_00200008247602008250400000_0152501500000000" > CutSelection.log
 echo -e "TPCCluster2040\nNo\n0\nYes\nPbPb_2.76TeV\nNo\nNo\nYes" > answersFileSysErr
 bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersFileSysErr


# # dEdxE Cutstudies # # 
 echo -e "52400013_00200009247602008250404000_0152501500000000\n5240001_00200009347602008250400000_0152501500000000\n5240001_00200009647602008250400000_0152501500000000" > CutSelection.log
 echo -e "dEdxE2040\nNo\n0\nYes\nPbPb_2.76TeV\nNo\nNo\nYes" > answersFileSysErr
 bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersFileSysErr

 
# # dEdxPi Cutstudies # #
 echo -e "52400013_00200009247602008250404000_0152501500000000\n5240001_00200009287602008250400000_0152501500000000\n5240001_00200009237002008250400000_0152501500000000\n5240001_00200009245402008250400000_0152501500000000" > CutSelection.log
 echo -e "dEdxPi2040\nNo\n0\nYes\nPbPb_2.76TeV\nNo\nNo\nYes" > answersFileSysErr
 bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersFileSysErr

 
# pdEdxPi Cutstudies # #
#  echo -e "5240001_00200009247602008250400000_0152501500000000\n5240001_00200009245402008250400000_0152501500000000" > CutSelection.log
#  echo -e "pdEdxPi2040\nNo\n0\nYes\nPbPb_2.76TeV\nNo\nNo\nYes" > answersFileSysErr
# bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersFileSysErr
 

 # cosPA Cutstudies # #
#  echo -e "5240001_00200009247602008250400000_0152501500000000\n5240001_00200009247602008250000000_0152501500000000" > CutSelection.log
#  echo -e "CosPoint2040\nNo\n0\nYes\nPbPb_2.76TeV\nNo\nNo\nYes" > answersFileSysErr
# bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersFileSysErr

 
# TOF Cutstudies # #
#  echo -e "5240001_00200009247602008250400000_0152501500000000\n5240001_00200009247603008250400000_0152501500000000\n5240001_00200009247604008250400000_0152501500000000" > CutSelection.log
#  echo -e "TOF2040\nNo\n0\nYes\nPbPb_2.76TeV\nNo\nNo\nYes" > answersFileSysErr
# bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersFileSysErr
 

# # Qt Cutstudies # #
 echo -e "52400013_00200009247602008250404000_0152501500000000\n5240001_00200009247602009250400000_0152501500000000\n5240001_00200009247602002250400000_0152501500000000" > CutSelection.log
 echo -e "Qt2040\nNo\n0\nYes\nPbPb_2.76TeV\nNo\nNo\nYes" > answersFileSysErr
 bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersFileSysErr

 
# # Chi2 & PsiPair Cutstudies # #
 echo -e "52400013_00200009247602008250404000_0152501500000000\n5240001_00200009247602008150400000_0152501500000000\n5240001_00200009247602008850400000_0152501500000000\n5240001_00200009247602008280400000_0152501500000000\n5240001_00200009247602008260400000_0152501500000000" > CutSelection.log
 echo -e "Chi22040\nNo\n0\nYes\nPbPb_2.76TeV\nNo\nNo\nYes" > answersFileSysErr
 bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersFileSysErr

 
# Alpha Cutstudies # #
echo -e "52400013_00200009247602008250404000_0152501500000000\n5240001_00200009247602008250400000_0152505500000000\n5240001_00200009247602008250400000_0152503500000000" > CutSelection.log
echo -e "Alpha2040\nNo\n0\nYes\nPbPb_2.76TeV\nNo\nNo\nYes" > answersFileSysErr
bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersFileSysErr


# Phi Cutstudies # #
echo -e "52400013_00200009247602008250404000_0152501500000000\n5240001_00215509247602008250400000_0152501500000000\n5240001_00217709247602008250400000_0152501500000000" > CutSelection.log
echo -e "ConvPhi2040\nNo\n0\nYes\nPbPb_2.76TeV\nNo\nNo\nYes" > answersFileSysErr
bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersFileSysErr

# # # BG Cutstudies # #
# echo -e "5240001_00200009247602008250400000_0152501500000000\n0000011_00200009227302008250400000_0252103500000000" > CutSelection.log
# echo -e "BG2040\nNo\n0\nYes\nPbPb_2.76TeV\nNo\nNo\nYes" > answersFileSysErr
# bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersFileSysErr


#! /bin/bash

# # Yield Cutstudies
#echo -e "" > CutSelection.log
#bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersPi0MBYield


# # Eta Cutstudies # #
# echo -e "5250001_00200009247602008250400000_0152501500000000\n5250001_03200009247602008250400000_0152301500000000\n5250001_04200009247602008250400000_0152201500000000" > CutSelection.log
# cat CutSelection.log
# echo -e "Eta2050\nNo\n0\nYes\nPbPb_2.76TeV\nNo\nNo\nYes" > answersFileSysErr
# bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersFileSysErr


# # SinglePt Cutstudies # #
 echo -e "52500013_00200009247602008250404000_0152501500000000\n5250001_00200049247602008250400000_0152501500000000\n5250001_00200019247602008250400000_0152501500000000" > CutSelection.log
 echo -e "SinglePt2050\nNo\n0\nYes\nPbPb_2.76TeV\nNo\nNo\nYes" > answersFileSysErr
 bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersFileSysErr


# # TPC Cluster Cutstudies # #
 echo -e "52500013_00200009247602008250404000_0152501500000000\n5250001_00200006247602008250400000_0152501500000000\n5250001_00200008247602008250400000_0152501500000000" > CutSelection.log
 echo -e "TPCCluster2050\nNo\n0\nYes\nPbPb_2.76TeV\nNo\nNo\nYes" > answersFileSysErr
 bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersFileSysErr


# # dEdxE Cutstudies # # 
 echo -e "52500013_00200009247602008250404000_0152501500000000\n5250001_00200009347602008250400000_0152501500000000\n5250001_00200009647602008250400000_0152501500000000" > CutSelection.log
 echo -e "dEdxE2050\nNo\n0\nYes\nPbPb_2.76TeV\nNo\nNo\nYes" > answersFileSysErr
 bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersFileSysErr

 
# # dEdxPi Cutstudies # #
 echo -e "52500013_00200009247602008250404000_0152501500000000\n5250001_00200009287602008250400000_0152501500000000\n5250001_00200009237002008250400000_0152501500000000\n5250001_00200009245402008250400000_0152501500000000" > CutSelection.log
 echo -e "dEdxPi2050\nNo\n0\nYes\nPbPb_2.76TeV\nNo\nNo\nYes" > answersFileSysErr
 bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersFileSysErr

 
# pdEdxPi Cutstudies # #
#  echo -e "5250001_00200009247602008250400000_0152501500000000\n5250001_00200009245402008250400000_0152501500000000" > CutSelection.log
#  echo -e "pdEdxPi2050\nNo\n0\nYes\nPbPb_2.76TeV\nNo\nNo\nYes" > answersFileSysErr
# bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersFileSysErr
 

 # cosPA Cutstudies # #
#  echo -e "5250001_00200009247602008250400000_0152501500000000\n5250001_00200009247602008250000000_0152501500000000" > CutSelection.log
#  echo -e "CosPoint2050\nNo\n0\nYes\nPbPb_2.76TeV\nNo\nNo\nYes" > answersFileSysErr
# bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersFileSysErr

 
# TOF Cutstudies # #
#  echo -e "5250001_00200009247602008250400000_0152501500000000\n5250001_00200009247603008250400000_0152501500000000\n5250001_00200009247604008250400000_0152501500000000" > CutSelection.log
#  echo -e "TOF2050\nNo\n0\nYes\nPbPb_2.76TeV\nNo\nNo\nYes" > answersFileSysErr
# bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersFileSysErr
 

# # Qt Cutstudies # #
 echo -e "52500013_00200009247602008250404000_0152501500000000\n5250001_00200009247602009250400000_0152501500000000\n5250001_00200009247602002250400000_0152501500000000" > CutSelection.log
 echo -e "Qt2050\nNo\n0\nYes\nPbPb_2.76TeV\nNo\nNo\nYes" > answersFileSysErr
 bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersFileSysErr

 
# # Chi2 & PsiPair Cutstudies # #
 echo -e "52500013_00200009247602008250404000_0152501500000000\n5250001_00200009247602008150400000_0152501500000000\n5250001_00200009247602008850400000_0152501500000000\n5250001_00200009247602008280400000_0152501500000000\n5250001_00200009247602008260400000_0152501500000000" > CutSelection.log
 echo -e "Chi22050\nNo\n0\nYes\nPbPb_2.76TeV\nNo\nNo\nYes" > answersFileSysErr
 bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersFileSysErr

 
# Alpha Cutstudies # #
echo -e "52500013_00200009247602008250404000_0152501500000000\n5250001_00200009247602008250400000_0152505500000000\n5250001_00200009247602008250400000_0152503500000000" > CutSelection.log
echo -e "Alpha2050\nNo\n0\nYes\nPbPb_2.76TeV\nNo\nNo\nYes" > answersFileSysErr
bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersFileSysErr


# Phi Cutstudies # #
echo -e "52500013_00200009247602008250404000_0152501500000000\n5250001_00215509247602008250400000_0152501500000000\n5250001_00217709247602008250400000_0152501500000000" > CutSelection.log
echo -e "ConvPhi2050\nNo\n0\nYes\nPbPb_2.76TeV\nNo\nNo\nYes" > answersFileSysErr
bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersFileSysErr

# # # BG Cutstudies # #
# echo -e "5250001_00200009247602008250400000_0152501500000000\n0000011_00200009227302008250400000_0252103500000000" > CutSelection.log
# echo -e "BG2050\nNo\n0\nYes\nPbPb_2.76TeV\nNo\nNo\nYes" > answersFileSysErr
# bash start_FullMesonAnalysis_TaskV2.sh -dgammaOff dummyFile pdf < answersFileSysErr


