#include <TH2F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TString.h>
#include <TLatex.h>
#include <TMarker.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <TStyle.h>
#include <algorithm>

// --- ANALYSIS system ---
#include <AliEMCALGeometry.h>


// usage:
// root -q -b -x -l PlotBadMap_EMC.C+
// first argument specifies if only EMCal ("EMCal"), DCal ("DCal") or EMCal + DCal ("") is plotted
// text files containing the channels are required (each channel in a new line)
// three input lists of bad channels can be chosen: bad cells, dead cells and additional bad cells
//bool CompareMaps has to be activated if one wants to compare two bad channel maps. In this case, SFileBadCells has to contain all bad + dead channels from map 1 and SFileDeadCells from map 2

void ReadCellsFromFile(TString myTextFile, std::vector<int> *vMyBadCells){

    std::ifstream file(myTextFile);	                        /// Input file
    std::string str;
    while (std::getline(file, str))
    {
        int i = 0;
        std::stringstream ss(str);
        while (ss >> i)
        {
            vMyBadCells->push_back(i);
            if (ss.peek() == ',' || ss.peek() == ' ') ss.ignore();
        }

    }

  }


void PlotBadMap_EMC(TString SDetector = "", TString SPeriod = "LHC16x", TString SCollSystem = "pp, #sqrt{#it{s}} = 13 TeV", bool CompareMaps = false){
//"EMCal" -> EMCal only, "DCal" -> DCal only, "" -> EMCal + DCal

    AliEMCALGeometry* geom;                                   //get the ALICE EMCal geometry
    geom = AliEMCALGeometry::GetInstance("EMCAL_COMPLETE12SMV1_DCAL_8SM");
    TString SPathToTextFile = ".";                          // path to txt files
    TString SFileBadCells = Form("badcells_%s.txt", SPeriod.Data());
    TString SFileDeadCells = Form("deadcells_%s.txt", SPeriod.Data());
    TString SFileAddCells = Form("AdditionalBadCells_%s.txt", SPeriod.Data());

    std::vector<Int_t> vBadCells;                             // for storing the bad cells
    std::vector<Int_t> vDeadCells;                            // for storing the dead cells
    std::vector<Int_t> vAddBadCells;                          // for storing additional bad cells
    std::vector<Float_t> vAllCells;                           // for storing indormation about all cells


    Int_t nSupMod, nModule, nIphi, nIeta, iphi, ieta;       // Variables needed for the detector geometry
    Double_t eta, phi, etamin, etamax, phimin, phimax;
    Int_t xoffset = 0;
    Int_t yoffset = 0;

    Int_t absIdmax = geom->GetNCells();                     // Number of cells
    Int_t iStartCell = 0;                                   // needed for DCal (first cell in DCal is 12289)
    Int_t iPhiDCalStart = 0;                                // needed for DCal (position in Phi differs from EMCal (offset in phi needed))

    TH2D* hfullbadmap;                                      // 2-dim histogram for plotting the badmap


    // Initialize EMCal / DCal or EMCal+DCal geometry
    if(!SDetector.CompareTo("EMCal")){
        absIdmax = 12288;
        hfullbadmap = new TH2D("BadmapEMCal","BadmapEMCal", 96, 0., 96, 134, 0. , 134. );
    } else if(!SDetector.CompareTo("DCal")){
        iStartCell = 12288;
        iPhiDCalStart = 134;
        hfullbadmap = new TH2D("BadmapDCal","BadmapDCal", 96, 0., 96, 84, 0. , 84. );
    } else {                                                // EMCal + DCal
        hfullbadmap = new TH2D("BadmapEMCalDCal","BadmapEMCalDCal", 96, 0., 96, 218, 0. , 218. );
        SDetector="EMCal, DCal";
    }

    // Get The Bad and Dead Cells from the txt files (if the file does not exist, print something)
    if (std::fstream{ Form("%s/%s",SPathToTextFile.Data(), SFileBadCells.Data()) }) ReadCellsFromFile(Form("%s/%s",SPathToTextFile.Data(), SFileBadCells.Data()), &vBadCells);
    else std::cout << "no bad cell file was found..." << std::endl;
    if (std::fstream{ Form("%s/%s",SPathToTextFile.Data(), SFileDeadCells.Data()) }) ReadCellsFromFile(Form("%s/%s",SPathToTextFile.Data(), SFileDeadCells.Data()), &vDeadCells);
    else std::cout << "no dead cell file was found..." << std::endl;
    if (std::fstream{ Form("%s/%s",SPathToTextFile.Data(), SFileAddCells.Data()) }) ReadCellsFromFile(Form("%s/%s",SPathToTextFile.Data(), SFileAddCells.Data()), &vAddBadCells);
    else std::cout << "no additional bad cell file was found..." << std::endl;

    /// sort cells into bad, dead, add. bad, good
    for(int i = iStartCell; i < absIdmax; i++){
        // Just one map with good, bad, dead and if needed additional bad cells
        if(!CompareMaps){
            if(std::find(vBadCells.begin(), vBadCells.end(), i) != vBadCells.end()){
                vAllCells.push_back(3.5);
            } else if (std::find(vDeadCells.begin(), vDeadCells.end(), i) != vDeadCells.end()){
                vAllCells.push_back(1.5);
            } else if (std::find(vAddBadCells.begin(), vAddBadCells.end(), i) != vAddBadCells.end()){
                vAllCells.push_back(2.5);
            // Cell is good
            } else {
                vAllCells.push_back(0.5);
            }
        // Else statment if comparison between 2 maps is needed
        } else {
            // Cell is only bad in map 1
            if(std::binary_search(vBadCells.begin(), vBadCells.end(), i) && !(std::binary_search(vDeadCells.begin(), vDeadCells.end(), i))){
                vAllCells.push_back(3.5);
            // Cell is only bad in map 2
            } else if (std::binary_search(vDeadCells.begin(), vDeadCells.end(), i) && !(std::binary_search(vBadCells.begin(), vBadCells.end(), i))){
                vAllCells.push_back(1.5);
            // Cell is bad in map 1 and in map 2
            } else if (std::binary_search(vDeadCells.begin(), vDeadCells.end(), i)  && std::binary_search(vBadCells.begin(), vBadCells.end(), i)){
                vAllCells.push_back(2.5);
            // Cell is good
            } else {
                vAllCells.push_back(0.5);
            }
        }
    }


    // Fill the TH2 histogram with the values from vAllCells
    for(int absId = iStartCell; absId < absIdmax; absId++)
    {
        geom->GetCellIndex(absId, nSupMod, nModule, nIphi, nIeta);
        geom->GetCellPhiEtaIndexInSModule(nSupMod, nModule, nIphi, nIeta, iphi, ieta);

        if(nSupMod<20){
            xoffset = nSupMod % 2 * 49 ;
            yoffset = nSupMod / 2 * 25;

            if(nSupMod >11 && nSupMod < 18){
                xoffset = nSupMod % 2 * 65 ;
              }
            if(nSupMod > 11){
                yoffset = nSupMod / 2 * 25 - 16;
              }

            hfullbadmap->SetBinContent(1 +ieta  + xoffset,1 + iphi  + yoffset - iPhiDCalStart,vAllCells[absId - iStartCell]);
          }
    }


    TCanvas* CanBadmap = new TCanvas("","",1200,1200);
    CanBadmap->cd();
    CanBadmap->SetTopMargin(0.2);

    // Define a custom palette for plotting
    Int_t palette[4];
    palette[0] = kGreen - 7;
    palette[1] = kGray + 1;
    palette[2] = kBlack;
    palette[3] = kRed;
    gStyle->SetPalette(4,palette);

    hfullbadmap->SetStats(0);
    hfullbadmap->GetYaxis()->SetTitleOffset(1.3);
    hfullbadmap->GetXaxis()->SetTitleOffset(1.);
    hfullbadmap->GetYaxis()->SetTitleFont(43);
    hfullbadmap->GetYaxis()->SetTitleSize(36);
    hfullbadmap->GetYaxis()->SetLabelFont(43);
    hfullbadmap->GetYaxis()->SetLabelSize(36);
    hfullbadmap->GetXaxis()->SetTitleFont(43);
    hfullbadmap->GetXaxis()->SetTitleSize(36);
    hfullbadmap->GetXaxis()->SetLabelFont(43);
    hfullbadmap->GetXaxis()->SetLabelSize(36);
    hfullbadmap->SetTitle("");
    hfullbadmap->GetXaxis()->SetTitle("cell column (#it{#eta})");
    hfullbadmap->GetYaxis()->SetTitle("cell row (#it{#phi})");
    hfullbadmap->Draw("col");

    // Marker for the legend
    TMarker *markerbad = new TMarker(0,0,21);
    markerbad->SetMarkerColor(kRed);
    markerbad->SetMarkerSize(3.5);
    TMarker *markerdead = new TMarker(0,0,21);
    markerdead->SetMarkerColor(kGray + 1);
    markerdead->SetMarkerSize(3.5);
    TMarker *markergood = new TMarker(0,0,21);
    markergood->SetMarkerColor(kGreen -7);
    markergood->SetMarkerSize(3.5);
    TMarker *markerAdditionalBad = new TMarker(0,0,21);
    markerAdditionalBad->SetMarkerColor(kBlack);
    markerAdditionalBad->SetMarkerSize(3.5);

    // Latex for info: Period, Detecor, collision system & energy
    TLatex InfoLatex;
    InfoLatex.SetTextFont(42);
    InfoLatex.SetTextSize(0.035);
    InfoLatex.SetTextAlign(13);                              //align at top
    InfoLatex.DrawLatexNDC(.15,.96,SCollSystem.Data());
    InfoLatex.DrawLatexNDC(.15,.91,SPeriod.Data());
    InfoLatex.DrawLatexNDC(.15,.86,SDetector.Data());

    //Set the names for the Legend
    TString SnamesForLegend[4] = {"good", "bad", "dead", "add. bad"};
    TString SnamesForLegend_Comp[4] = {"good", "bad (1)", "bad (2)", "bad (1+2)"};
    if(CompareMaps) for(int i = 0; i < 4; ++i) SnamesForLegend[i] = SnamesForLegend_Comp[i];

    TLegend *legBadMap = new TLegend(0.5,0.83, 0.7, 0.97);
    legBadMap->SetBorderSize(0);
    legBadMap->SetTextFont(42);
    legBadMap->SetTextSize(0.035);
    legBadMap->AddEntry(markergood,SnamesForLegend[0].Data(),"p");
    legBadMap->AddEntry(markerAdditionalBad,SnamesForLegend[3].Data(),"p");
    legBadMap->Draw("same");

    TLegend *legBadMap2 = new TLegend(0.75,0.83, 0.97, 0.97);
    legBadMap2->SetBorderSize(0);
    legBadMap2->SetTextFont(42);
    legBadMap2->SetTextSize(0.035);
    legBadMap2->AddEntry(markerbad,SnamesForLegend[1].Data(),"p");
    legBadMap2->AddEntry(markerdead,SnamesForLegend[2].Data(),"p");
    legBadMap2->Draw("same");

    CanBadmap->SaveAs(Form("BadChannelMap_%s.pdf",SPeriod.Data()));
}
