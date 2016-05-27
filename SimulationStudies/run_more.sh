#!/bin/bash

## regular: plot MB, no pt cut
##narrow bins
root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","Pi0",1,1,1,1,0,0)'
root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","Pi0",1,1,1,1,0,2)'
root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","Pi0",1,1,1,1,0,3)'
root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","Pi0",1,1,1,1,0,4)'
root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","Pi0",1,1,1,1,0,5)'
root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","Pi0",0,1,1,1,0,0)'

root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","Eta",1,1,1,1,0,0)'
root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","Eta",1,1,1,1,0,2)'
root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","Eta",1,1,1,1,0,3)'
root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","Eta",1,1,1,1,0,4)'
root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","Eta",1,1,1,1,0,5)'
root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","Eta",0,1,1,1,0,0)'
#
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","EtaPrim",1,1,1,1,0,0)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","EtaPrim",1,1,1,1,0,2)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","EtaPrim",1,1,1,1,0,3)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","EtaPrim",1,1,1,1,0,4)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","EtaPrim",1,1,1,1,0,5)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","EtaPrim",0,1,1,1,0,0)'
#
##wide bins
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-Wide/PythiaAnalysisResults",4,"Pythia8","Monash","Pi0",1,1,0,1,0,0)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-Wide/PythiaAnalysisResults",4,"Pythia8","Monash","Pi0",1,1,0,1,0,2)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-Wide/PythiaAnalysisResults",4,"Pythia8","Monash","Pi0",1,1,0,1,0,3)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-Wide/PythiaAnalysisResults",4,"Pythia8","Monash","Pi0",1,1,0,1,0,4)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-Wide/PythiaAnalysisResults",4,"Pythia8","Monash","Pi0",1,1,0,1,0,5)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-Wide/PythiaAnalysisResults",4,"Pythia8","Monash","Pi0",0,1,0,1,0,0)'
#
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-Wide/PythiaAnalysisResults",4,"Pythia8","Monash","Eta",1,1,0,1,0,0)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-Wide/PythiaAnalysisResults",4,"Pythia8","Monash","Eta",1,1,0,1,0,2)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-Wide/PythiaAnalysisResults",4,"Pythia8","Monash","Eta",1,1,0,1,0,3)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-Wide/PythiaAnalysisResults",4,"Pythia8","Monash","Eta",1,1,0,1,0,4)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-Wide/PythiaAnalysisResults",4,"Pythia8","Monash","Eta",1,1,0,1,0,5)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-Wide/PythiaAnalysisResults",4,"Pythia8","Monash","Eta",0,1,0,1,0,0)'
#
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-Wide/PythiaAnalysisResults",4,"Pythia8","Monash","EtaPrim",1,1,0,1,0,0)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-Wide/PythiaAnalysisResults",4,"Pythia8","Monash","EtaPrim",1,1,0,1,0,2)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-Wide/PythiaAnalysisResults",4,"Pythia8","Monash","EtaPrim",1,1,0,1,0,3)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-Wide/PythiaAnalysisResults",4,"Pythia8","Monash","EtaPrim",1,1,0,1,0,4)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-Wide/PythiaAnalysisResults",4,"Pythia8","Monash","EtaPrim",1,1,0,1,0,5)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-Wide/PythiaAnalysisResults",4,"Pythia8","Monash","EtaPrim",0,1,0,1,0,0)'
#
## don't plot MB, no pt cut
##narrow bins
root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","Pi0",1,1,1,0,0,0)'
root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","Pi0",1,1,1,0,0,2)'
root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","Pi0",1,1,1,0,0,3)'
root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","Pi0",1,1,1,0,0,4)'
root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","Pi0",1,1,1,0,0,5)'
root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","Pi0",0,1,1,0,0,0)'

root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","Eta",1,1,1,0,0,0)'
root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","Eta",1,1,1,0,0,2)'
root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","Eta",1,1,1,0,0,3)'
root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","Eta",1,1,1,0,0,4)'
root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","Eta",1,1,1,0,0,5)'
root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","Eta",0,1,1,0,0,0)'
#
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","EtaPrim",1,1,1,0,0,0)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","EtaPrim",1,1,1,0,0,2)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","EtaPrim",1,1,1,0,0,3)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","EtaPrim",1,1,1,0,0,4)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","EtaPrim",1,1,1,0,0,5)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","EtaPrim",0,1,1,0,0,0)'
#
##wide bins
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-Wide/PythiaAnalysisResults",4,"Pythia8","Monash","Pi0",1,1,0,0,0,0)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-Wide/PythiaAnalysisResults",4,"Pythia8","Monash","Pi0",1,1,0,0,0,2)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-Wide/PythiaAnalysisResults",4,"Pythia8","Monash","Pi0",1,1,0,0,0,3)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-Wide/PythiaAnalysisResults",4,"Pythia8","Monash","Pi0",1,1,0,0,0,4)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-Wide/PythiaAnalysisResults",4,"Pythia8","Monash","Pi0",1,1,0,0,0,5)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-Wide/PythiaAnalysisResults",4,"Pythia8","Monash","Pi0",0,1,0,0,0,0)'
#
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-Wide/PythiaAnalysisResults",4,"Pythia8","Monash","Eta",1,1,0,0,0,0)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-Wide/PythiaAnalysisResults",4,"Pythia8","Monash","Eta",1,1,0,0,0,2)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-Wide/PythiaAnalysisResults",4,"Pythia8","Monash","Eta",1,1,0,0,0,3)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-Wide/PythiaAnalysisResults",4,"Pythia8","Monash","Eta",1,1,0,0,0,4)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-Wide/PythiaAnalysisResults",4,"Pythia8","Monash","Eta",1,1,0,0,0,5)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-Wide/PythiaAnalysisResults",4,"Pythia8","Monash","Eta",0,1,0,0,0,0)'
#
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-Wide/PythiaAnalysisResults",4,"Pythia8","Monash","EtaPrim",1,1,0,0,0,0)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-Wide/PythiaAnalysisResults",4,"Pythia8","Monash","EtaPrim",1,1,0,0,0,2)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-Wide/PythiaAnalysisResults",4,"Pythia8","Monash","EtaPrim",1,1,0,0,0,3)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-Wide/PythiaAnalysisResults",4,"Pythia8","Monash","EtaPrim",1,1,0,0,0,4)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-Wide/PythiaAnalysisResults",4,"Pythia8","Monash","EtaPrim",1,1,0,0,0,5)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-Wide/PythiaAnalysisResults",4,"Pythia8","Monash","EtaPrim",0,1,0,0,0,0)'
#
# don't plot MB, pt cut
#narrow bins
root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","Pi0",1,1,1,0,1.2,0)'
root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","Pi0",1,1,1,0,1.2,2)'
root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","Pi0",1,1,1,0,1.2,3)'
root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","Pi0",1,1,1,0,1.2,4)'
root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","Pi0",1,1,1,0,1.2,5)'
root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","Pi0",0,1,1,0,1.2,0)'

root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","Eta",1,1,1,0,1.2,0)'
root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","Eta",1,1,1,0,1.2,2)'
root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","Eta",1,1,1,0,1.2,3)'
root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","Eta",1,1,1,0,1.2,4)'
root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","Eta",1,1,1,0,1.2,5)'
root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","Eta",0,1,1,0,1.2,0)'

#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","EtaPrim",1,1,1,0,1.2,0)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","EtaPrim",1,1,1,0,1.2,2)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","EtaPrim",1,1,1,0,1.2,3)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","EtaPrim",1,1,1,0,1.2,4)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","EtaPrim",1,1,1,0,1.2,5)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","EtaPrim",0,1,1,0,1.2,0)'

#wide bins
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-WideFull/PythiaAnalysisResultsFull",4,"Pythia8","Monash","Pi0",1,1,0,0,1.2,0)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-WideFull/PythiaAnalysisResultsFull",4,"Pythia8","Monash","Pi0",1,1,0,0,1.2,2)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-WideFull/PythiaAnalysisResultsFull",4,"Pythia8","Monash","Pi0",1,1,0,0,1.2,3)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-WideFull/PythiaAnalysisResultsFull",4,"Pythia8","Monash","Pi0",1,1,0,0,1.2,4)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-WideFull/PythiaAnalysisResultsFull",4,"Pythia8","Monash","Pi0",1,1,0,0,1.2,5)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-WideFull/PythiaAnalysisResultsFull",4,"Pythia8","Monash","Pi0",0,1,0,0,1.2,0)'
#
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-WideFull/PythiaAnalysisResultsFull",4,"Pythia8","Monash","Eta",1,1,0,0,1.2,0)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-WideFull/PythiaAnalysisResultsFull",4,"Pythia8","Monash","Eta",1,1,0,0,1.2,2)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-WideFull/PythiaAnalysisResultsFull",4,"Pythia8","Monash","Eta",1,1,0,0,1.2,3)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-WideFull/PythiaAnalysisResultsFull",4,"Pythia8","Monash","Eta",1,1,0,0,1.2,4)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-WideFull/PythiaAnalysisResultsFull",4,"Pythia8","Monash","Eta",1,1,0,0,1.2,5)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-WideFull/PythiaAnalysisResultsFull",4,"Pythia8","Monash","Eta",0,1,0,0,1.2,0)'
#
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-Wide/PythiaAnalysisResults",4,"Pythia8","Monash","EtaPrim",1,1,0,0,1.2,0)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-Wide/PythiaAnalysisResults",4,"Pythia8","Monash","EtaPrim",1,1,0,0,1.2,2)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-Wide/PythiaAnalysisResults",4,"Pythia8","Monash","EtaPrim",1,1,0,0,1.2,3)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-Wide/PythiaAnalysisResults",4,"Pythia8","Monash","EtaPrim",1,1,0,0,1.2,4)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-Wide/PythiaAnalysisResults",4,"Pythia8","Monash","EtaPrim",1,1,0,0,1.2,5)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-Wide/PythiaAnalysisResults",4,"Pythia8","Monash","EtaPrim",0,1,0,0,1.2,0)'

## don't plot MB, 2nd pt cut
##narrow bins
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","Pi0",1,1,1,0,0.5,0)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","Pi0",1,1,1,0,0.5,2)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","Pi0",1,1,1,0,0.5,3)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","Pi0",1,1,1,0,0.5,4)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","Pi0",1,1,1,0,0.5,5)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","Pi0",0,1,1,0,0.5,0)'
#
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","Eta",1,1,1,0,0.5,0)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","Eta",1,1,1,0,0.5,2)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","Eta",1,1,1,0,0.5,3)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","Eta",1,1,1,0,0.5,4)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","Eta",1,1,1,0,0.5,5)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","Eta",0,1,1,0,0.5,0)'
#
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","EtaPrim",1,1,1,0,0.5,0)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","EtaPrim",1,1,1,0,0.5,2)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","EtaPrim",1,1,1,0,0.5,3)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","EtaPrim",1,1,1,0,0.5,4)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","EtaPrim",1,1,1,0,0.5,5)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",9,"Pythia8","Monash","EtaPrim",0,1,1,0,0.5,0)'
#
##wide bins
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-Wide/PythiaAnalysisResults",4,"Pythia8","Monash","Pi0",1,1,0,0,0.5,0)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-Wide/PythiaAnalysisResults",4,"Pythia8","Monash","Pi0",1,1,0,0,0.5,2)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-Wide/PythiaAnalysisResults",4,"Pythia8","Monash","Pi0",1,1,0,0,0.5,3)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-Wide/PythiaAnalysisResults",4,"Pythia8","Monash","Pi0",1,1,0,0,0.5,4)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-Wide/PythiaAnalysisResults",4,"Pythia8","Monash","Pi0",1,1,0,0,0.5,5)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-Wide/PythiaAnalysisResults",4,"Pythia8","Monash","Pi0",0,1,0,0,0.5,0)'
#
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-Wide/PythiaAnalysisResults",4,"Pythia8","Monash","Eta",1,1,0,0,0.5,0)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-Wide/PythiaAnalysisResults",4,"Pythia8","Monash","Eta",1,1,0,0,0.5,2)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-Wide/PythiaAnalysisResults",4,"Pythia8","Monash","Eta",1,1,0,0,0.5,3)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-Wide/PythiaAnalysisResults",4,"Pythia8","Monash","Eta",1,1,0,0,0.5,4)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-Wide/PythiaAnalysisResults",4,"Pythia8","Monash","Eta",1,1,0,0,0.5,5)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-Wide/PythiaAnalysisResults",4,"Pythia8","Monash","Eta",0,1,0,0,0.5,0)'
#
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-Wide/PythiaAnalysisResults",4,"Pythia8","Monash","EtaPrim",1,1,0,0,0.5,0)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-Wide/PythiaAnalysisResults",4,"Pythia8","Monash","EtaPrim",1,1,0,0,0.5,2)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-Wide/PythiaAnalysisResults",4,"Pythia8","Monash","EtaPrim",1,1,0,0,0.5,3)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-Wide/PythiaAnalysisResults",4,"Pythia8","Monash","EtaPrim",1,1,0,0,0.5,4)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-Wide/PythiaAnalysisResults",4,"Pythia8","Monash","EtaPrim",1,1,0,0,0.5,5)'
#root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia8-Wide/PythiaAnalysisResults",4,"Pythia8","Monash","EtaPrim",0,1,0,0,0.5,0)'
