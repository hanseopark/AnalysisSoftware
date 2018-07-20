#! /bin/bash

######################################################################################################
########                               heavy ion Setups                                       ########
######################################################################################################

#######################################################################################################
# PbPb 2.76TeV
#######################################################################################################

function GiveBinningPbPb2760()
{
    DOETAPRIME=0
    if [ $DOPI0 -eq 1 ] ; then
        echo -e "How many p_T bins do you want to use for Pi0?"
        echo -e "\t PCM, EMC, PCM-EMC pi0: 19 (10GeV), 20(12GeV), 21(14GeV), 22(16GeV), 23(18GeV), 24(20GeV), 25(25GeV), 26(30GeV)";
        echo -e "\t PCM dir gamma: 14(6GeV), 17(14GeV), 18(14GeV), 18(11GeV), 19 (20GeV), 22(14GeV)";
        read answer
        if [ $answer -ge 14 ] && [ $answer -lt 27 ] ; then
            echo "--> $answer Bins";
            CORRECTPI0=1
            BINSPTPI0=$answer
        else
            echo "--> Pi0 binning was not set correctly. Please try again.";
            CORRECTPI0=0
        fi
        BINSPTGAMMA=$BINSPTPI0
    fi
    if [ $DOETA -eq 1 ] || [ $DOPI0INETABINS -eq 1 ]; then
        echo "How many p_T bins do you want to use for Eta? 8(8GeV), 9(10GeV), 10(12GeV), 11(15GeV), 12(20GeV), 13(25GeV), 14(30GeV)";
        read answer
        if [ $answer -ge 8 ] && [ $answer -lt 15 ] ; then
        if [ $answer = 8 ]; then
            echo "--> $answer Bins";
            CORRECTETA=1
            BINSPTETA=$answer
        else
            echo "--> Eta Binning was not set correctly. Please try again.";
            CORRECTETA=0
        fi
    fi
}

#######################################################################################################
# PbPb 5.02TeV
#######################################################################################################
function GiveBinningPbPb5TeV()
{
    DOETAPRIME=0
    if [ $DOPI0 -eq 1 ] ; then
        echo "How many p_T bins do you want to use for Pi0? up to 24(20GeV)";
        read answer
        if [ "$answer" -le "24" ]; then
            CORRECTPI0=1
            BINSPTPI0=$answer
        else
            echo "--> Pi0 binning was not set correctly. Please try again.";
            CORRECTPI0=0
        fi
        BINSPTGAMMA=$BINSPTPI0
    fi
    if [ $DOETA -eq 1 ] || [ $DOPI0INETABINS -eq 1 ]; then
        echo "How many p_T bins do you want to use for Eta? up to 22 (30GeV)";
        read answer
        if [ "$answer" -le "22" ]; then
            CORRECTETA=1
            BINSPTETA=$answer
        else
            echo "Eta Binning was not set correctly. Please try again.";
            CORRECTETA=0
        fi
    fi
}


#######################################################################################################
# XeXe 5.44TeV
#######################################################################################################
function GiveBinningXeXe5440GeV()
{
    DOETAPRIME=0
    if [ $DOPI0 -eq 1 ] || [ $DOGAMMA -eq 1 ]; then
        echo "How many p_T bins do you want to use for Pi0? up to 25";
        read answer
        if [ "$answer" -le "25" ]; then
            CORRECTPI0=1
            BINSPTPI0=$answer
        else
            echo "--> Pi0 binning was not set correctly. Please try again.";
            CORRECTPI0=0
        fi
        BINSPTGAMMA=$BINSPTPI0
    fi
    if [ $DOETA -eq 1 ] || [ $DOPI0INETABINS -eq 1 ]; then
        echo "How many p_T bins do you want to use for Eta? up to 9";
        read answer
        if [ "$answer" -le "9" ]; then
            CORRECTETA=1
            BINSPTETA=$answer
        else
            echo "--> Pi0 binning was not set correctly. Please try again.";
            CORRECTETA=0
        fi
    fi
}
