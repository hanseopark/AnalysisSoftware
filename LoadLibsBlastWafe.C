void LoadLibs() {

  gSystem->Load("libCore.so");  
  gSystem->Load("libGeom.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libVMC");
  gSystem->Load("libTree");
  gSystem->Load("libProof");
  gSystem->Load("libMatrix");
  gSystem->Load("libMinuit");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libOADB");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libTENDER");
  gSystem->Load("libCORRFW");
  //  gSystem->Load("libPWG0base");
  gSystem->Load("libMinuit");
  //  gSystem->Load("libPWG2spectra");
  gSystem->Load("libPWGTools");
  gROOT->LoadMacro("$ALICE_ROOT/PWGLF/SPECTRA/PiKaPr/COMBINED/FitParticle.C+");
  // gSystem->Load("libPWG0dep");
  // gSystem->Load("libPWG0selectors");
  // gSystem->Load("libPWG1");
  
  //   gSystem->Load("libSTEER") ;  
  //   gSystem->Load("libSTEERBase") ;  
  //   gSystem->Load("libESD") ;  
  //   gSystem->Load("libAOD") ;  
  //   gSystem->Load("libANALYSIS") ;  
  //   gSystem->Load("libANALYSISalice") ;
  //   gSystem->Load("libTENDER") ;
  //   gSystem->Load("libCORRFW") ;
  

}
