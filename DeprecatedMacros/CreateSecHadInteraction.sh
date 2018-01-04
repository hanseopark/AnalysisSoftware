root -l -q -x -b TaskV1/BuildHistogramsForSecHadInteractions.C\+\+\(\"$2/fbock_SecHadIntTree.root\"\,\"$1\"\,0\,0\,\"MappingDetailedSecHadInt_LHC11aPass4\"\)
root -l -q -x -b TaskV1/BuildHistogramsForSecHadInteractions.C\+\+\(\"$3/fbock_SecHadIntTree.root\"\,\"$1\"\,1\,0\,\"MappingDetailedSecHadInt_LHC12f1a_Pythia8\"\)
root -l -q -x -b TaskV1/BuildHistogramsForSecHadInteractions.C\+\+\(\"$4/fbock_SecHadIntTree.root\"\,\"$1\"\,1\,0\,\"MappingDetailedSecHadInt_LHC12f1b_Phojet\"\)

root -l -q -x -b TaskV1/BuildHistogramsForSecHadInteractions.C\+\+\(\"$3/fbock_SecHadIntTree.root\"\,\"$1\"\,1\,0\,\"MappingDetailedSecHadInt_LHC12f1ab\"\)
root -l -q -x -b TaskV1/BuildHistogramsForSecHadInteractions.C\+\+\(\"$4/fbock_SecHadIntTree.root\"\,\"$1\"\,1\,1\,\"MappingDetailedSecHadInt_LHC12f1ab\"\)

root -l -q -x -b TaskV1/MappingMaterialHadronicAdv.C\+\+\(\"MappingDetailedSecHadInt_LHC11aPass4_Data.root\"\,\"MappingDetailedSecHadInt_LHC12f1ab_MC.root\"\,\"$1\",\"$5\",\"2.76TeV\",\"merged\"\,\"\LHC11a\"\)
root -l -q -x -b TaskV1/MappingMaterialHadronicAdv.C\+\+\(\"MappingDetailedSecHadInt_LHC11aPass4_Data.root\"\,\"MappingDetailedSecHadInt_LHC12f1a_Pythia8_MC.root\"\,\"$1\",\"$5\",\"2.76TeV\",\"Pythia8\"\,\"\LHC11a\"\)
root -l -q -x -b TaskV1/MappingMaterialHadronicAdv.C\+\+\(\"MappingDetailedSecHadInt_LHC11aPass4_Data.root\"\,\"MappingDetailedSecHadInt_LHC12f1b_Phojet_MC.root\"\,\"$1\",\"$5\",\"2.76TeV\",\"Phojet\"\,\"\LHC11a\"\)