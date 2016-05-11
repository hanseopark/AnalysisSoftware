
void ClusterQA_ReadEMCalBadCellsRun()
{
  gSystem->AddIncludePath("-I/home/daniel/alice/aliroot/master/inst/include");
  gSystem->AddIncludePath("-I/home/daniel/alice/alphysics/master/inst/include");
  gROOT->ProcessLine(".L TaskQA/ClusterQA_ReadEMCalBadCells.C+");
  gROOT->ProcessLine("ClusterQA_ReadEMCalBadCells()");
  return;
}
