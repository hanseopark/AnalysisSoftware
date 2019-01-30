// R__ADD_INCLUDE_PATH($ALICE_ROOT)
// R__ADD_INCLUDE_PATH($ALICE_PHYSICS)

void run_EMCalAccMap()
{
  // gSystem->AddIncludePath("-I/home/nschmidt/alice/sw/ubuntu1604_x86-64/AliRoot/ali-master-1/include");
  // gSystem->AddIncludePath("-I/home/nschmidt/alice/sw/ubuntu1604_x86-64/AliPhysics/ali-master-1/include");
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
  gROOT->ProcessLine("PrepareEMCalAcceptanceMap(\"SM10\")");
  gROOT->ProcessLine("PrepareEMCalAcceptanceMap(\"SM11\")");
  gROOT->ProcessLine("PrepareEMCalAcceptanceMap(\"SM12\")");
  gROOT->ProcessLine("PrepareEMCalAcceptanceMap(\"SM13\")");
  gROOT->ProcessLine("PrepareEMCalAcceptanceMap(\"SM14\")");
  gROOT->ProcessLine("PrepareEMCalAcceptanceMap(\"SM15\")");
  gROOT->ProcessLine("PrepareEMCalAcceptanceMap(\"SM16\")");
  gROOT->ProcessLine("PrepareEMCalAcceptanceMap(\"SM17\")");
  gROOT->ProcessLine("PrepareEMCalAcceptanceMap(\"SM18\")");
  gROOT->ProcessLine("PrepareEMCalAcceptanceMap(\"SM19\")");
  gROOT->ProcessLine("PrepareEMCalAcceptanceMap(\"ModAccEMCAL\")");
  gROOT->ProcessLine("PrepareEMCalAcceptanceMap(\"ModAcc12\")");
  gROOT->ProcessLine("PrepareEMCalAcceptanceMap(\"ModAcc21\")");
  gROOT->ProcessLine("PrepareEMCalAcceptanceMap(\"ModAcc\")");
  gROOT->ProcessLine("PrepareEMCalAcceptanceMap(\"ModAcc2\")");

  return;
}
