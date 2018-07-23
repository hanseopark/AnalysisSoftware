#! /bin/bash

######################################################################################################
########                                      PP Setups                                       ########
######################################################################################################


#######################################################################################################
# PP 0.9TeV
#######################################################################################################
function GiveBinning900GeV()
{
    DOETAPRIME=0
    if [ $DOPI0 -eq 1 ] || [ $DOGAMMA -eq 1 ]; then
        echo "How many p_T bins do you want to use? 10(3GeV), 11(4GeV)";
        read answer
        if [ $answer -ge 7 ] && [ $answer -lt 13 ] ; then
            echo "--> $answer Bins ";
            CORRECTPI0=1
            BINSPTPI0=$answer
        else
            echo "--> Pi0 binning was not set correctly. Please try again.";
        fi
        BINSPTGAMMA=$BINSPTPI0
    fi

    if [ $DOETA -eq 1 ] || [ $DOPI0INETABINS -eq 1 ]; then
        echo "How many p_T bins do you want to use for the Eta meson? 2 (1.8 GeV), 3 (3 GeV), only EMC related: 4 (5 GeV)"
        read answer
        if [ $answer -ge 2 ] && [ $answer -lt 4 ]; then
            echo "--> $answer Bins ";
            CORRECTETA=1
            BINSPTETA=$answer
        else
            echo "--> Eta binning was not set correctly. Please try again.";
            CORRECTETA=0
        fi
    fi

}


#######################################################################################################
# PP 2.76TeV
#######################################################################################################

function GiveBinning2760GeV()
{
    DOETAPRIME=0
    if [ $DOPI0 -eq 1 ] || [ $DOGAMMA -eq 1 ]; then
        echo "How many p_T bins do you want to use?"
        echo "  PCM: 17 (6GeV), 18 (8GeV), 19 (10GeV)"
        echo "  PCM-EMC: 20 (12 GeV), 21 (15 GeV), 22 (20 GeV), 23 (25 GeV), 24 (30 GeV)";
        echo "  EMC or triggers: 21 (10GeV),  24 (14GeV), 25 (16GeV), 26 (20 GeV), 27 (24 GeV), 28 (28GeV) 29 (30GeV)";
        echo "  mEMC: 25 (16 GeV), 26 (18 GeV), 27 (22 GeV), 28 (26 GeV), 29 (30 GeV), 30 (35 GeV), 31 (40 GeV), 32 (50 GeV)";
        read answer
        if [ $MODE -eq 10 ] || [ $MODE -eq 11 ]; then
            if [ $answer -ge 25 ] && [ $answer -lt 33 ] ; then
                echo "--> $answer Bins";
                CORRECTPI0=1
                BINSPTPI0=$answer
            else
                echo "--> Pi0 binning was not set correctly. Please try again.";
            fi
            DOETA=0;
            DOPI0INETABINS=0;
            CORRECTETA=1;
        else
            if [ $answer -ge 17 ] && [ $answer -lt 30 ] ; then
                echo "--> $answer Bins";
                CORRECTPI0=1
                BINSPTPI0=$answer
            else
                echo "--> Pi0 binning was not set correctly. Please try again.";
            fi
        fi
        BINSPTGAMMA=$BINSPTPI0
    fi

    if [ $DOETA -eq 1 ] || [ $DOPI0INETABINS -eq 1 ]; then
        echo "How many p_T bins do you want to use for the Eta meson?"
        echo "  PCM: 6 (4. GeV), 7 (6 GeV)"
        echo "  PCM-EMC & EMC: 8 (8 GeV), 9 (10 GeV), 10 (12 GeV), 11 (16 GeV), 12 (20 GeV), 13 (25 GeV)"
        read answer
        if [ $answer -ge 6 ] && [ $answer -lt 14 ] ; then
            echo "--> $answer Bins";
            CORRECTETA=1
            BINSPTETA=$answer
        else
            echo "--> Eta binning was not set correctly. Please try again.";
            CORRECTETA=0
        fi
    else
        CORRECTETA=1
    fi
    FILEWEIGHTGAMMA="FitsPaperPP2760GeV_2017_07_10_FrediV2Clusterizer_Gamma_2017_12_18.root";
    echo "Setting weighting file to: " $FILEWEIGHTGAMMA
}


#######################################################################################################
# PP 5TeV
#######################################################################################################


#################################################################
# return number of bins for pi0 and eta for 5 TeV pp running 2015
#################################################################
function GiveBinning5TeV()
{
    if [ $DOPI0 -eq 1 ] || [ $DOGAMMA -eq 1 ]; then
        echo "How many p_T bins do you want to use for the Pi0? "
        echo "  PCM pi0: 12+13+14 & 19 (XX GeV/c) - 24 (XX GeV/c).. 32,33,34 (XX GeV/c), 42 (16 GeV/c)";
        echo "  PCM gamma dir: 19(10.0GeV), 20(12.0GeV), 21(15.0GeV)";
        read answer
        if [ $answer -lt 90 ]; then
            echo "--> $answer bins";
            CORRECTPI0=1
            BINSPTPI0=$answer
        else
            echo "--> Pi0 binning was not set correctly. please try again.";
            CORRECTPI0=0
        fi
        BINSPTGAMMA=$BINSPTPI0
    fi

    if [ $DOETA -eq 1 ] || [ $DOPI0INETABINS -eq 1 ]; then
        echo "How many p_T bins do you want to use for the Eta? 7(XX GeV/c) - 13(XX GeV/c).. 19,20,21,22 (XX GeV/c)";
        read answer
        if [ $answer -lt 40 ]; then
            echo "--> $answer bins";
            CORRECTETA=1
            BINSPTETA=$answer
        else
            echo "--> Eta binning was not set correctly. please try again.";
            CORRECTETA=0
        fi
    fi
}

#################################################################
# return number of bins for pi0 and eta for 5 TeV pp running 2017
#################################################################
function GiveBinning5TeV2017()
{
    if [ $DOPI0 -eq 1 ] || [ $DOGAMMA -eq 1 ]; then
        echo "How many p_T bins do you want to use for the Pi0? "
        echo "  PCM pi0: 26 (as Hikari, 14 GeV/c), 57 (14 GeV/c), 80 (finer, 14 GeV/c)";
        echo "  PCM gamma dir: 19(10.0GeV), 20(12.0GeV), 21(15.0GeV)";
        read answer
        if [ $answer -gt 4 ] && [ $answer -lt 90 ]; then
            echo "--> $answer bins";
            CORRECTPI0=1
            BINSPTPI0=$answer
        else
            echo "--> Pi0 binning was not set correctly. Please try again.";
            CORRECTPI0=0
        fi
        BINSPTGAMMA=$BINSPTPI0
    fi

    if [ $DOETA -eq 1 ] || [ $DOPI0INETABINS -eq 1 ]; then
        echo "How many p_T bins do you want to use for the Eta? 8 (as Hikari 12 GeV/c), 20 (12 GeV/c)";
        read answer
        if [ $answer -gt 4 ] && [ $answer -lt 34 ]; then
            echo "--> $answer bins";
            CORRECTETA=1
            BINSPTETA=$answer
        else
            echo "--> Eta binning was not set correctly. Please try again.";
            CORRECTETA=0
        fi
    fi
}


#######################################################################################################
# PP 7TeV
#######################################################################################################

function GiveBinning7TeV()
{
    DOETAPRIME=0
    if [ $DOPI0 -eq 1 ] || [ $DOGAMMA -eq 1 ]; then
        echo "How many p_T bins do you want to use for the Pi0? "
        echo "  PCM, EMC, PCM-EMC pi0: 36(7GeV/c), 37(8GeV/c), 38(10GeV/c), 39(12GeV/c), 40 (16GeV/c), 41 (20GeV/c), 42 (25GeV/c)";
        echo "  PCM, EMC, PCM-EMC dir gamma: 19(8GeV), 20(10GeV), 21(12GeV), 22(16GeV) 23 (20GeV)"
        read answer
        if [ $answer -ge 19 ] && [ $answer -lt 45 ]; then
            echo "--> $answer bins ";
            CORRECTPI0=1
            BINSPTPI0=$answer
        else
            echo "--> Pi0 binning was not set correctly. please try again.";
            CORRECTPI0=0
        fi
        BINSPTGAMMA=$BINSPTPI0
    fi

    if [ $DOETA -eq 1 ] || [ $DOPI0INETABINS -eq 1 ]; then
        echo "How many p_T bins do you want to use for the Eta? 9(3GeV/c), 11(4GeV/c), 12(5GeV/c), 13(6GeV/c), 14(8GeV/c) 15(10GeV/c)";
        read answer
        if [ $answer -ge 10 ] && [ $answer -lt 19 ] ; then
            echo "--> $answer bins";
            CORRECTETA=1
            BINSPTETA=$answer
        else
            echo "--> Eta binning was not set correctly. please try again.";
            CORRECTETA=0
        fi
    fi
}


#######################################################################################################
# PP 8TeV
#######################################################################################################

#################################################################
# return number of bins for pi0 and eta for 8 TeV pp
#################################################################
function GiveBinning8TeV()
{
    if [ $DOPI0 -eq 1 ] || [ $DOGAMMA -eq 1 ]; then
        echo "How many p_T bins do you want to use for the Pi0? "
        echo "  PCM, PCM-EMC: 23(7GeV), 24(8GeV), 25(10GeV), 26(12GeV) 27 (16GeV) 28 (25GeV)"
        echo "  EMC:  33 (20 GeV)"
        echo "  mEMC: 40 (26GeV), 48 (50GeV), 60 (100GeV)";
        echo "  PCM, EMC, PCM-EMC dir gamma: 19(8GeV), 20(10GeV), 21(12GeV), 22(16GeV) 23 (20GeV)"
        read answer
        if [ $answer -lt 61 ]; then
        echo "--> $answer bins chosen for Pi0...";
            CORRECTPI0=1
            BINSPTPI0=$answer
        else
            echo "--> Pi0 binning was not set correctly. Please try again.";
            CORRECTPI0=0
        fi
    fi

    if [ $DOETA -eq 1 ] || [ $DOPI0INETABINS -eq 1 ]; then
        echo "How many p_T bins do you want to use for the Eta? "
        echo "  PCM: 10(4GeV), 12(6GeV), 13(8GeV), 14(10GeV), 15(12GeV)";
        read answer
        if [ $answer -ge 10 ] && [ $answer -lt 27 ]; then
            echo "--> $answer Bins";
            CORRECTETA=1
            BINSPTETA=$answer
       else
            echo "--> Eta binning was not set correctly. Please try again.";
            CORRECTETA=0
       fi
       BINSPTGAMMA=$BINSPTPI0
    fi
    if [ $DOETAPRIME -eq 1 ]; then
        echo "How many p_T bins do you want to use for the EtaPrime? 10(4GeV), 12(6GeV), 13(8GeV), 14(10GeV), 15(12GeV)";
        read answer
        if [ $answer -ge 10 -a $answer -le 26 ]; then
            CORRECTETAPrime=1
            BINSPTETAPRIME=$answer
        else
            echo "EtaPrime Binning was not set correctly. Please try again.";
            CORRECTETAPrime=0
        fi
        BINSPTGAMMA=$BINSPTPI0
    fi
}



#######################################################################################################
# PP 13TeV
#######################################################################################################

#################################################################
# return number of bins for pi0,eta,eta' for 13 TeV pp running
#################################################################
function GiveBinning13TeV()
{
    if [ $DOPI0 -eq 1 ] || [ $DOGAMMA -eq 1 ]; then
        echo "How many p_T bins do you want to use for the Pi0?"
        echo "  for PCM: Max. 17 for 15f, 20 for 15fhi, 17 for low B 17g ";
        echo "  for other modes: max 37";
        read BINSPTPI0
        CORRECTPI0=1
        echo "--> You have chosen $BINSPTPI0 bins";
        BINSPTGAMMA=$BINSPTPI0
    fi

    if [ $DOETA -eq 1 ] || [ $DOPI0INETABINS -eq 1 ]; then
        echo "How many p_T bins do you want to use for the Eta meson?"
        echo "  for PCM: Max. 7 for 15f, 13 for 15fhi, 4 for low B 15g";
        echo "  for other modes: max 24";
        read BINSPTETA
        CORRECTETA=1
        echo "--> You have chosen $BINSPTETA bins";

    fi
    if [ $DOETAPRIME -eq 1 ]; then
        echo "How many p_T bins do you want to use for the EtaPrime meson?"
        echo "  for PCM: Max. 7 for 15f, 13 for 15fhi, 4 for low B 15g";
        echo "  fFor other modes: max 24";
        read BINSPTETAPRIME
        CORRECTETAPrime=1
        echo "--> You have chosen $BINSPTETAPRIME bins";
    fi
}
