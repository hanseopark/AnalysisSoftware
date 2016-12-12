void run_EMCalAccMap()
{
  gSystem->AddIncludePath("-I/home/daniel/alice/sw/ubuntu1404_x86-64/AliRoot/latest-ali-master/include");
  gSystem->AddIncludePath("-I/home/daniel/alice/sw/ubuntu1404_x86-64/AliPhysics/latest-ali-master/include");
  gROOT->ProcessLine(".L PrepareEMCalAcceptanceMap.C+");

  gROOT->ProcessLine("PrepareEMCalAcceptanceMap(\"SM-0\")");
  gROOT->ProcessLine("PrepareEMCalAcceptanceMap(\"SM-1\")");
  gROOT->ProcessLine("PrepareEMCalAcceptanceMap(\"SM-2\")");
  gROOT->ProcessLine("PrepareEMCalAcceptanceMap(\"SM-3\")");
  gROOT->ProcessLine("PrepareEMCalAcceptanceMap(\"SM-4\")");
  gROOT->ProcessLine("PrepareEMCalAcceptanceMap(\"SM-5\")");
  gROOT->ProcessLine("PrepareEMCalAcceptanceMap(\"SM-6\")");
  gROOT->ProcessLine("PrepareEMCalAcceptanceMap(\"SM-7\")");
  gROOT->ProcessLine("PrepareEMCalAcceptanceMap(\"SM-8\")");
  gROOT->ProcessLine("PrepareEMCalAcceptanceMap(\"SM-9\")");
  gROOT->ProcessLine("PrepareEMCalAcceptanceMap(\"mod_acc-EMCAL\")");
  gROOT->ProcessLine("PrepareEMCalAcceptanceMap(\"mod_acc\")");

  return;
}
