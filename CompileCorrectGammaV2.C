#include <Riostream.h>
#include <fstream>
#include "TMath.h"
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TPaveLabel.h>
#include <TSystem.h>
#include <TFrame.h>
#include <TStyle.h>
#include <TString.h>
#include "TGaxis.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TF1.h"
#include "TVirtualFitter.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TDatabasePDG.h"
#include "TMinuit.h"
#include "TLatex.h"
#include "TASImage.h"
#include "TPostScript.h"
#include "TGraphErrors.h"
#include "TArrow.h"
#include "TMarker.h"
#include "TGraphAsymmErrors.h" 


void CompileCorrectGammaV2(){
    
     
     //*****************************************************************************************
    // Find out which user is running the code and set RooUnfold directory accordingly
    // RooUnfold is a separate piece of software which needs to be compiled on your system
    // Please download it from https://github.com/skluth/RooUnfold and compile it!
    // Afterwards you have to add the location of the software here for your system
    //*****************************************************************************************
    TString homedirectory = gSystem->HomeDirectory();
	cout << "***************************************" << endl;
	cout << "HomeDirectory: " << homedirectory.Data() << endl;
	cout << "***************************************" << endl;
    if (homedirectory.CompareTo("/home/admin1") == 0){   
        gSystem->AddIncludePath("-I/home/admin1/leardini/photonconv/AnalysisSoftware/RooUnfold/RooUnfold-master/RooUnfold/src"); 
        gSystem->Load("/home/admin1/leardini/photonconv/AnalysisSoftware/RooUnfold/RooUnfold-master/libRooUnfold");
    } else if (homedirectory.CompareTo("/home/fbock") == 0){   
        gSystem->AddIncludePath("-I/home/fbock/Photon/Software/photonconv/AnalysisSoftware/RooUnfold/src"); 
		gSystem->Load("/home/fbock/Photon/Software/photonconv/AnalysisSoftware/RooUnfold/libRooUnfold");
	} else if (homedirectory.CompareTo("/home/daniel") == 0){
		gSystem->AddIncludePath("-I/home/daniel/data/work/photonconv/AnalysisSoftware/RooUnfold/src");
		gSystem->Load("/home/daniel/data/work/photonconv/AnalysisSoftware/RooUnfold/libRooUnfold");
    } else if (homedirectory.CompareTo("/Users/lucasaltenkaemper") == 0) {
        gSystem->AddIncludePath("-I/Volumes/MacintoshHD/Documents/Studium/Master/Thesis/photonconv/AnalysisSoftware/RooUnfold/src");
        gSystem->Load("/Volumes/MacintoshHD/Documents/Studium/Master/Thesis/photonconv/AnalysisSoftware/RooUnfold/libRooUnfold");
    } else {
        cout << "You have not defined where RooUnfold can be found on your system! This macro can't run without it!" << endl;
        return;
    }

    gROOT->LoadMacro("TaskV1/CorrectGammaV2.C+");

    return;
}
