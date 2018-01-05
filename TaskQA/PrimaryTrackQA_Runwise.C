/*******************************************************************************
 ******  provided by Gamma Conversion Group, PWGGA,                        *****
 ******     Daniel Muehlheim, d.muehlheim@cern.ch                          *****
 ******     Friederike Bock, fbock@cern.ch                                 *****
 *******************************************************************************/
#include "QA.h"
//**************************************************************************************************************
//***************************** Main routine *******************************************************************
//**************************************************************************************************************
void PrimaryTrackQA_Runwise(
        Int_t nSetsIn,                          // number of data sets to be analysed
        TString fEnergyFlag,                    // energy flag
        TString* DataSets,                      // technical names of data sets for output
        TString* plotDataSets,                  // labels of data sets in plots
        TString* pathDataSets,                  // path for data sets
        Int_t mode              = 2,            // standard mode for analysis
        Int_t cutNr             = -1,           // if -1: you have to choose number at runtime
        Int_t doExtQA           = 2,            // 0: switched off, 1: normal extQA, 2: with Cell level plots
        TString suffix          = "eps",        // output format of plots
        TString labelData       = "Data",       // Label for data
        Bool_t addSubfolder     = kFALSE        // flag to enable subdirectory creation for primary cut
        )
{
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "PrimaryTrackQA_Runwise" << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "Done with PrimaryTrackQA" << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    return;
}
