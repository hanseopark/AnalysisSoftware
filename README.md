# Analysis framework of the ALICE PCM group

The basic structure of this framework is as following:

For the neutral pion/eta analysis the following macros need to be run in this order:
-> TaskV1/ExtractSignalV2.C for data and MC
---> TaskV1/TaskV1/CompareMesonQuantities.C is optional and compare the output of the previous macro between data and MC
-> TaskV1/CorrectSignalV2.C
-> TaskV1/ProduceFinalResultsV2.C (pp w/o triggers), TaskV1/ProduceFinalResultspPb.C (pPb, w/o triggers), TaskV1/ProduceFinalResultsPatchedTriggers.C (pp/pPb, w triggers), TaskV1/CalculateRAA.C (PbPb)

The steering of the first 3 is done by 
-> start_FullMesonAnalysis_TaskV2.sh
