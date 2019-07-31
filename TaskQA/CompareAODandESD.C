#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include "TLatex.h"

#include "./CommonHeaders/PlottingMeson.h"
#include "./CommonHeaders/PlottingGammaConversionHistos.h"
#include "./CommonHeaders/PlottingGammaConversionAdditional.h"
#include "./CommonHeaders/ExtractSignalBinning.h"
#include "./CommonHeaders/FittingGammaConversion.h"
#include "./CommonHeaders/ConversionFunctions.h"
#include "./CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"


// Macro to compare the histograms contained in and AOD and ESD.
// Note that this macro only works with ROOT6 and that both root files
// have to have the same folder structure.
//
// author: Florian Jonas (florian.jonas@cern.ch)


void CompareAODandESD(TString nameFileESDInput                = "/media/florianjonas/dataslave/data/alice/LocalTests/AnalysisSoftware/LocalTrainTests/pp_5TeV/LHC18j2_cent_woSDD/OESD/GammaConvNeutralMesonPiPlPiMiNeutralMeson_1_0_111.root",
                      TString nameFileAODInput                = "/media/florianjonas/dataslave/data/alice/LocalTests/AnalysisSoftware/LocalTrainTests/pp_5TeV/LHC18j2_cent_woSDD/OAOD/GammaConvNeutralMesonPiPlPiMiNeutralMeson_1_0_111.root",
                      Bool_t  compareTH2                      = kTRUE,   // if activated TH2 are compared by plotting them next to each other
                      Bool_t  isMC                            = kFALSE,  // set to kTRUE if both files are MC files
                      TString fileType                        = "eps"){

    StyleSettingsThesis();
    SetPlotStyle();

    // load input file
    TFile* fileESDInput                   = new TFile(nameFileESDInput);
    TFile* fileAODInput                   = new TFile(nameFileAODInput);

    // find main dir
    TString nameMainDir;
    TKey *key;
    TIter next(fileESDInput->GetListOfKeys());

    // get first main dir found
    while ( (key=(TKey*)next()) ){
        cout << Form("-> found TopDir: %s",key->GetName());
        nameMainDir = key->GetName();
        if(nameMainDir.CompareTo("")!=0){
            break; // take first key found
        }
    }
    if(nameMainDir.CompareTo("")!=0){
        cout << "Succesfully detected " << nameMainDir <<" as MainDir!" << endl;
    } else{
        printf("ERROR: MainDir not found!");
        return;
    }

    TList *TopDirESD =(TList*)fileESDInput->Get(nameMainDir.Data());
    TList *TopDirAOD =(TList*)fileAODInput->Get(nameMainDir.Data());

    if(TopDirESD == NULL){
        cout<<"ERROR: TopDir ESD not Found"<<endl;
        return;
    }
    if(TopDirAOD == NULL){
        cout<<"ERROR: TopDir AOD not Found"<<endl;
        cout<<"ERROR: Trying to find different top dir"<<endl;
        return;
    }

    TList *ESDContainer;
    TList *AODContainer;

    if(TopDirAOD->GetEntries() != TopDirESD->GetEntries()){
        cout<<"ERROR: Number of entries in top dir of AOD and ESD not identical!"<<endl;
        return;
    }

    for(Int_t i = 0; i<TopDirESD->GetEntries(); i++){
        TList *SubDirESD = (TList*) TopDirESD->At(i);
        TList *SubDirAOD = (TList*) TopDirAOD->At(i);
        TString dirname = SubDirESD->GetName();

        TString CutNumber(dirname(11,dirname.Length()-11));

        TString SubDirName = SubDirESD->GetName();
        if(SubDirName.Contains("Cut Number")){
            for(Int_t c = 0; c<SubDirESD->GetEntries();c++){
                ESDContainer             = (TList*) SubDirESD->At(c);
                AODContainer             = (TList*) SubDirAOD->At(c);

                TString ContainerName = ESDContainer->GetName();
                SubDirName.ReplaceAll(" ","_");
                ContainerName.ReplaceAll(" ","_");

                TString outputDir;

                if(isMC){
                    outputDir       = Form("AODandESDComparison/MC_%s/%s/%s",CutNumber.Data(),SubDirName.Data(),ContainerName.Data());
                } else{
                    outputDir       = Form("AODandESDComparison/%s/%s/%s",CutNumber.Data(),SubDirName.Data(),ContainerName.Data());
                }
                gSystem->Exec("mkdir -p "+outputDir);

                // plot all ESD histos
                // Create Canvas
                TCanvas* canvasTH1 = new TCanvas("canvasTH1","",10,10,1000,1000);  // gives the page size
                TCanvas* canvasTH2 = new TCanvas("canvasTH2","",10,10,1000,1000);  // gives the page size


                for(Int_t hist = 0; hist < ESDContainer->GetEntries();hist++){
                    TString histName = ESDContainer->At(hist)->GetName();
                    cout << "Reading ESD histo " << histName.Data()<< endl;
                    if( histName.CompareTo(AODContainer->At(hist)->GetName()) != 0){
                        cout << "Error names not identical" << endl;
                        printf("ESDName = %s \t AODName = %s \n",histName.Data(),AODContainer->At(hist)->GetName());
                        return;
                    }

                    canvasTH1->SetLogz();
                    canvasTH2->SetLogz();
                    canvasTH1->SetLeftMargin(0.13);
                    canvasTH2->SetLeftMargin(0.13);
                    canvasTH1->SetRightMargin(0.12);
                    canvasTH2->SetRightMargin(0.12);
                    canvasTH1->SetTopMargin(0.05);
                    canvasTH2->SetTopMargin(0.05);
                    canvasTH1->SetBottomMargin(0.1);
                    canvasTH2->SetBottomMargin(0.1);

                    canvasTH2->Divide(2,1);

                    TString ESDClassName(ESDContainer->At(hist)->ClassName());
                    TString AODClassName(AODContainer->At(hist)->ClassName());

                    if(ESDClassName.CompareTo(AODClassName) != 0){
                        printf("ERROR: classes are not identical \n");
                        break;
                    }

                    if  (ESDClassName.Contains("TH1")){

                        TH1F* ESDHistTH1 = (TH1F*) ESDContainer->At(hist)->Clone();
                        TH1F* AODHistTH1 = (TH1F*) AODContainer->At(hist)->Clone();
                        AODHistTH1->SetLineColor(kRed);

                        auto rp = new TRatioPlot(ESDHistTH1,AODHistTH1,"divsym");

                        canvasTH1->cd();
                        rp->Draw("");
                        auto* legend = new TLegend(0.8, 0.8, .9, .9);
                        legend->AddEntry(ESDHistTH1, "ESD", "l");
                        legend->AddEntry(AODHistTH1, "AOD", "l");
                        legend->Draw();
                        canvasTH1->Update();
                        canvasTH1->Print(Form("%s/%s.%s",outputDir.Data(),histName.Data(),fileType.Data()),fileType.Data());
                        canvasTH1->Clear();
                    } else if (ESDClassName.Contains("TH2")){

                        if(!compareTH2){
                            printf("TH2 processing disabled, continuing ->\n");
                            continue;
                        }
                        TH2F* ESDHistTH2 = (TH2F*) ESDContainer->At(hist)->Clone();
                        TH2F* AODHistTH2 = (TH2F*) AODContainer->At(hist)->Clone();

                        canvasTH2->cd(1);
                        TText *tESD = new TText(.8,.9,"ESD");
                        tESD->SetNDC(kTRUE);
                        ESDHistTH2->Draw("colz");
                        tESD->SetTextAlign(22);
                        tESD->SetTextColor(kBlack);
                        tESD->SetTextFont(43);
                        tESD->SetTextSize(30);
                        tESD->Draw();

                        canvasTH2->cd(2);
                        TText *tAOD = new TText(.8,.9,"AOD");
                        tAOD->SetNDC(kTRUE);
                        AODHistTH2->Draw("colz");
                        tAOD->SetTextAlign(22);
                        tAOD->SetTextColor(kBlack);
                        tAOD->SetTextFont(43);
                        tAOD->SetTextSize(30);
                        tAOD->Draw();
                        canvasTH2->Update();
                        canvasTH2->Print(Form("%s/%s.%s",outputDir.Data(),histName.Data(),fileType.Data()),fileType.Data());
                        canvasTH2->Clear();
                    } else {
                        // TF1 is for now ignored.
                        printf("Error: unknwon object type %s -> continuing ...\n",ESDClassName.Data());
                        continue;
                    }
                }

            }
        }

        TString outputDir;

        cout << "Plotting container " << SubDirName.Data();
        // Plot preselection histos

        if(isMC){
            outputDir       = Form("AODandESDComparison/MC_%s/%s",CutNumber.Data(),SubDirName.Data());
        } else{
            outputDir       = Form("AODandESDComparison/%s/%s",CutNumber.Data(),SubDirName.Data());
        }
        gSystem->Exec("mkdir -p "+outputDir);

        for(Int_t hist = 0; hist < SubDirESD->GetEntries();hist++){
                    TString histName = SubDirESD->At(hist)->GetName();
                    cout << "Reading ESD histo " << histName.Data()<< endl;
                    if( histName.CompareTo(SubDirAOD->At(hist)->GetName()) != 0){
                        cout << "Error names not identical" << endl;
                        printf("ESDName = %s \t AODName = %s \n",histName.Data(),SubDirAOD->At(hist)->GetName());
                        return;
                    }

                    canvasTH1->SetLogz();
                    canvasTH2->SetLogz();
                    canvasTH1->SetLeftMargin(0.13);
                    canvasTH2->SetLeftMargin(0.13);
                    canvasTH1->SetRightMargin(0.12);
                    canvasTH2->SetRightMargin(0.12);
                    canvasTH1->SetTopMargin(0.05);
                    canvasTH2->SetTopMargin(0.05);
                    canvasTH1->SetBottomMargin(0.1);
                    canvasTH2->SetBottomMargin(0.1);

                    canvasTH2->Divide(2,1);

                    TString ESDClassName(SubDirESD->At(hist)->ClassName());
                    TString AODClassName(SubDirAOD->At(hist)->ClassName());

                    if(ESDClassName.CompareTo(AODClassName) != 0){
                        printf("ERROR: classes are not identical \n");
                        break;
                    }

                    if  (ESDClassName.Contains("TH1")){

                        ESDHistTH1 = (TH1F*) SubDirESD->At(hist)->Clone();
                        AODHistTH1 = (TH1F*) SubDirAOD->At(hist)->Clone();
                        AODHistTH1->SetLineColor(kRed);

                        rp = new TRatioPlot(ESDHistTH1,AODHistTH1,"divsym");

                        canvasTH1->cd();
                        rp->Draw("");
                        legend = new TLegend(0.8, 0.8, .9, .9);
                        legend->AddEntry(ESDHistTH1, "ESD", "l");
                        legend->AddEntry(AODHistTH1, "AOD", "l");
                        legend->Draw();
                        canvasTH1->Update();
                        canvasTH1->Print(Form("%s/%s.%s",outputDir.Data(),histName.Data(),fileType.Data()),fileType.Data());
                        canvasTH1->Clear();
                    } else if (ESDClassName.Contains("TH2")){

                        if(!compareTH2){
                            printf("TH2 processing disabled, continuing ->\n");
                            continue;
                        }
                        ESDHistTH2 = (TH2F*) SubDirESD->At(hist)->Clone();
                        AODHistTH2 = (TH2F*) SubDirAOD->At(hist)->Clone();

                        canvasTH2->cd(1);
                        tESD = new TText(.8,.9,"ESD");
                        tESD->SetNDC(kTRUE);
                        ESDHistTH2->Draw("colz");
                        tESD->SetTextAlign(22);
                        tESD->SetTextColor(kBlack);
                        tESD->SetTextFont(43);
                        tESD->SetTextSize(30);
                        tESD->Draw();

                        canvasTH2->cd(2);
                        tAOD = new TText(.8,.9,"AOD");
                        tAOD->SetNDC(kTRUE);
                        AODHistTH2->Draw("colz");
                        tAOD->SetTextAlign(22);
                        tAOD->SetTextColor(kBlack);
                        tAOD->SetTextFont(43);
                        tAOD->SetTextSize(30);
                        tAOD->Draw();
                        canvasTH2->Update();
                        canvasTH2->Print(Form("%s/%s.%s",outputDir.Data(),histName.Data(),fileType.Data()),fileType.Data());
                        canvasTH2->Clear();
                    } else {
                        // TF1 is for now ignored.
                        printf("Error: unknwon object type %s -> continuing ...\n",ESDClassName.Data());
                        continue;
                    }
            }
        

    }

}
