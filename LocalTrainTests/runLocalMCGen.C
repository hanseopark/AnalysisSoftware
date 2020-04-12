#if !defined(__CINT__) || defined(__CLING__)
#include "AliMCEventHandler.h"
#include "AliESDInputHandler.h"
#include "AliAODInputHandler.h"
#include "AliDummyHandler.h"
#include "AliAnalysisAlien.h"
#include "AliAnalysisManager.h"

R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
// #include <PWG/EMCAL/macros/AddTaskAodSkim.C>
// #include <PWG/EMCAL/macros/AddTaskEsdSkim.C>
#include <PWGGA/GammaConv/macros/AddTask_GammaMCStudies.C>
#include <PWGGA/GammaConv/macros/Add_MCGenPythia8_TuneX.C>



#endif
#include "TGrid.h"


//______________________________________________________________________________
void runLocalMCGen(
    TString         runMode       = "PQ2HC",
    UInt_t          numEvents     = 50,
    TString         generator     = "MCGenPythia8TuneX",
    Float_t         e_cms         = 5020,
    Int_t           tune          = 14
)
{
    // since we will compile a class, tell root where to look for headers
    #if !defined (__CINT__) || defined (__CLING__)
        gInterpreter->ProcessLine(".include $ROOTSYS/include");
        gInterpreter->ProcessLine(".include $ALICE_ROOT/include");
        //gInterpreter->ProcessLine(".L /home/florianjonas/tools/pythia8240/lib/libpythia8.a");
    #else
        gROOT->ProcessLine(".include $ROOTSYS/include");
        gROOT->ProcessLine(".include $ALICE_ROOT/include");
        //gROOT->ProcessLine(".L $/lib/libAliPythia8.so");

    #endif

    // Create analysis manager
    AliAnalysisManager* mgr                     = new AliAnalysisManager("LocalMCGenTaskRunning");
    AliDummyHandler*    dumH        = new AliDummyHandler();

    // create blank ESD event
    AliESDEvent *esdE               = new AliESDEvent();
    esdE->CreateStdContent();
    AliESDVertex *vtx               = new AliESDVertex(0.,0.,100);
    vtx->SetName("VertexTracks");
    vtx->SetTitle("VertexTracks");
    esdE->SetPrimaryVertexTracks(vtx);
    if(esdE->GetPrimaryVertex()) Printf("vtx set");
    dumH->SetEvent(esdE);
    mgr->SetInputEventHandler(dumH);

    // greate MC input Handler
    AliMCGenHandler* mcInputHandler = new AliMCGenHandler();
    mgr->SetMCtruthEventHandler(mcInputHandler);

    // Create Generator

    AliGenerator* gener             = 0x0;

    if(generator.Contains("MCGenPythia8TuneX")){
	    #if !defined (__CINT__) || defined (__CLING__)
		    gener=reinterpret_cast<AliGenerator*>(
		    gInterpreter->ExecuteMacro("$ALICE_PHYSICS/PWGGA/GammaConv/macros/Add_MCGenPythia8_TuneX.C  ( 5020.,14, kTRUE, 1., 0 )"));
	    #else
		    gROOT->LoadMacro("$ALICE_PHYSICS/PWGGA/GammaConv/macros/Add_MCGenPythia8_TuneX.C");
		    gener = Add_MCGenPythia8_TuneX.C  ( 5020.,14, kTRUE, 1., 0 );
	    #endif
    } else{
        cout << " No generator specified! returning ..." << endl;
    }

    mcInputHandler->SetGenerator(gener);
    mcInputHandler->SetSeedMode(2); // check what this does

    // -----------------------------------------
    //               Gamma MC Studies
    // -----------------------------------------
    if(runMode.Contains("MS")){
        #if !defined (__CINT__) || defined (__CLING__)
            AliAnalysisTask *taskGammaMCStudies=reinterpret_cast<AliAnalysisTask*>(
                    gInterpreter->ExecuteMacro(Form("$ALICE_PHYSICS/PWGGA/GammaConv/macros/AddTask_GammaMCStudies.C ()")));
        #else
            gROOT->LoadMacro("$ALICE_PHYSICS/PWGGA/GammaConv/macros/AddTask_GammaMCStudies.C");
            AliAnalysisTask *taskGammaMCStudies = AddTask_GammaMCStudies();
        #endif
    }


    mgr->SetUseProgressBar(1, 10);
    if (!mgr->InitAnalysis()) return;
    mgr->PrintStatus();
    mgr->EventLoop(numEvents);

   
}

