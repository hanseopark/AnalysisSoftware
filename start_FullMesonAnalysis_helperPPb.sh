#! /bin/bash

######################################################################################################
########                                  pPb Setups                                          ########
######################################################################################################

#######################################################################################################
# pPb 5.023TeV
#######################################################################################################

function GiveBinningpPb5()
{
     if [ $DOPI0 -eq 1 ] || [ $DOGAMMA -eq 1 ]; then
        echo "How many p_T bins do you want to use for Pi0? "
        echo "  PCM: 27 (7 GeV), 28 (8 GeV), 29 (10 GeV), 30 (12 GeV), 31 (14 GeV)"
        echo "  PCM-EMC, EMC: 32 (16 GeV), 33 (18 GeV), 34 (20 GeV), 35 (22 GeV), 36 (26 GeV), 37 (30 GeV)";
        echo "  EMC trigger bins reach to 46."
        echo "  gamma dir: 21 (8 GeV), 22 (10 GeV), 23 (14 GeV), PCM-EMC 33"
        read answer
        BINSPTPI0=$answer
        CORRECTPI0=1
        echo "--> You have chosen " $answer " pt bins for Pi0";
        BINSPTGAMMA=$BINSPTPI0
    else
        CORRECTPI0=1
    fi
    if [ $DOETA -eq 1 ] || [ $DOPI0INETABINS -eq 1 ]; then
        echo "How many p_T bins do you want to use for the Eta meson?"
        echo "  PCM: 12 (4 GeV), 14 (6 GeV), 15 (8 GeV), 16 (10 GeV)"
        echo "  PCM-EMC, EMC: 17 (12 GeV), 18 (14 GeV), 19 (16 GeV), 20 (20 GeV), 21 (25 GeV), 22 (30 GeV)";
        echo "  for calo triggers bins reach to 30."
        read answer
        BINSPTETA=$answer
        CORRECTETA=1
        echo "--> You have chosen " $answer " pt bins for Eta";
    fi
    if [ $DOETAPRIME -eq 1 ]; then
        echo "How many p_T bins do you want to use for the Eta' meson? "
        echo "  PCM: 8 (10 GeV)";
        echo "  PCM-EMC: 8 (10 GeV)";
        echo "  EMC: 8 (10 GeV)";
        read answer
        BINSPTETAPRIME=$answer
        CORRECTETAPRIME=1
        echo "--> You have chosen " $answer " pt bins for Eta";
    fi
    if [ $DOGAMMA -eq 1 ]; then
        if [ $MODE == 2 ] || [ $MODE == 13 ] ; then
            echo "How many p_T bins do you want to use for direct photon?  ";
            read answer
            BINSPTGAMMA=$answer
            CORRECTPI0=1
            echo "--> You have chosen " $answer " pt bins for dir gamma";
            directphoton="Gamma"
        fi
    fi
}

#######################################################################################################
# pPb 8TeV
#######################################################################################################

function GiveBinningpPb8()
{
    if [ $DOPI0 -eq 1 ] || [ $DOGAMMA -eq 1 ]; then
        echo "How many p_T bins do you want to use for Pi0? "
        echo "  PCM: 27 (7 GeV), 28 (8 GeV), 29 (10 GeV), 30 (12 GeV), 31 (14 GeV)"
        echo "  PCM-EMC, EMC: 32 (16 GeV), 33 (18 GeV), 34 (20 GeV), 35 (22 GeV), 36 (26 GeV), 37 (30 GeV)";
        echo "  EMC trigger bins reach to 46."
        echo "  gamma dir: 21 (8 GeV), 22 (10 GeV), 23 (14 GeV), PCM-EMC 33"
        read answer
        BINSPTPI0=$answer
        CORRECTPI0=1
        echo "--> You have chosen " $answer " pt bins for Pi0";
        BINSPTGAMMA=$BINSPTPI0
    else
        CORRECTPI0=1
    fi
    if [ $DOETA -eq 1 ] || [ $DOPI0INETABINS -eq 1 ]; then
        echo "How many p_T bins do you want to use for the Eta meson?"
        echo "  PCM: 12 (4 GeV), 14 (6 GeV), 15 (8 GeV), 16 (10 GeV)"
        echo "  PCM-EMC, EMC: 17 (12 GeV), 18 (14 GeV), 19 (16 GeV), 20 (20 GeV), 21 (25 GeV), 22 (30 GeV)";
        echo "  for calo triggers bins reach to 30."
        read answer
        BINSPTETA=$answer
        CORRECTETA=1
        echo "--> You have chosen " $answer " pt bins for Eta";
    fi
    if [ $DOETAPRIME -eq 1 ]; then
        echo "How many p_T bins do you want to use for the Eta' meson? "
        echo "  PCM: 8 (10 GeV)";
        echo "  PCM-EMC: 8 (10 GeV)";
        echo "  EMC: 8 (10 GeV)";
        read answer
        BINSPTETAPRIME=$answer
        CORRECTETAPRIME=1
        echo "--> You have chosen " $answer " pt bins for Eta";
    fi
    if [ $DOGAMMA -eq 1 ]; then
        if [ $MODE == 2 ] || [ $MODE == 13 ] ; then
            echo "How many p_T bins do you want to use for direct photon?  ";
            read answer
            BINSPTGAMMA=$answer
            CORRECTPI0=1
            echo "--> You have chosen " $answer " pt bins for dir gamma";
            directphoton="Gamma"
        fi
    fi
}
