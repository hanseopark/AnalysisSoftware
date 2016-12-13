void run_EMCalAccMap()
{
  gSystem->AddIncludePath("-I/home/daniel/alice/sw/ubuntu1404_x86-64/AliRoot/latest-ali-master/include");
  gSystem->AddIncludePath("-I/home/daniel/alice/sw/ubuntu1404_x86-64/AliPhysics/latest-ali-master/include");
  gROOT->ProcessLine(".L PrepareEMCalAcceptanceMap.C+");

  gROOT->ProcessLine("PrepareEMCalAcceptanceMap(\"SM0\")");
  gROOT->ProcessLine("PrepareEMCalAcceptanceMap(\"SM1\")");
  gROOT->ProcessLine("PrepareEMCalAcceptanceMap(\"SM2\")");
  gROOT->ProcessLine("PrepareEMCalAcceptanceMap(\"SM3\")");
  gROOT->ProcessLine("PrepareEMCalAcceptanceMap(\"SM4\")");
  gROOT->ProcessLine("PrepareEMCalAcceptanceMap(\"SM5\")");
  gROOT->ProcessLine("PrepareEMCalAcceptanceMap(\"SM6\")");
  gROOT->ProcessLine("PrepareEMCalAcceptanceMap(\"SM7\")");
  gROOT->ProcessLine("PrepareEMCalAcceptanceMap(\"SM8\")");
  gROOT->ProcessLine("PrepareEMCalAcceptanceMap(\"SM9\")");
  gROOT->ProcessLine("PrepareEMCalAcceptanceMap(\"ModAccEMCAL\")");
  gROOT->ProcessLine("PrepareEMCalAcceptanceMap(\"ModAcc12\")");
  gROOT->ProcessLine("PrepareEMCalAcceptanceMap(\"ModAcc21\")");
  gROOT->ProcessLine("PrepareEMCalAcceptanceMap(\"ModAcc\")");
  gROOT->ProcessLine("PrepareEMCalAcceptanceMap(\"ModAcc2\")");

  return;
}
