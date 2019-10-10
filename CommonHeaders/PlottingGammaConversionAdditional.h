/***********************************************************************************************
*** provided by Gamma Conversion Group, PWG4, 									******
***	    Friederike Bock, fbock@physi.uni-heidelberg.de ***							******
************************************************************************************************/

#ifndef GAMMACONV_PlottingAdd
#define GAMMACONV_PlottingAdd

    /************************************************************************************************
    ************************************************************************************************
    This header contains the functions to plot additional things, like logos and Lines
    in your histograms
    ************************************************************************************************
    ************************************************************************************************

    The functions are
    - DrawAliceLogo
    - DrawAliceLogoOmega
    - DrawAliceLogo1D
    - DrawAliceLogoPerformance
    - DrawAliceLogoPerformance2D
    - DrawAliceText
    - DrawAliceLogoMonteCarlo
    - DrawAliceLogoPi0MonteCarlo

    - DrawStructure
    - DrawArmenteros
    - DrawdEdxLabels

    - DrawGammaLines
    // */

    #ifndef __CLING__
        #include <Riostream.h>
        #include <TH3.h>
        #include <TLatex.h>
    #endif
    /*********************************************************************************************************
    DrawAliceLogoPi0Performance
    * will draw you the ALICE Logo as well the text "ALICE Performance", "pp @ 7 TeV" and the date
    which you can hand over
    Float_t startTextX, Float_t startTextY, Float_t startPi0TextX, Float_t differenceText,
    Float_t startLogoX, Float_t startLogoY, Float_t widthLogo,
    Float_t textSize,
    TString collisionSystem, Bool_t mcFile
    * float_t startX, float_t startY - give starting Point of Logo
    * float_t widthLogo - gives the width of the logo
    * float_t textHeight - gives you the heigth of the text
    * char * date - date handed over
    **********************************************************************************************************
    **********************************************************************************************************/

    void DrawAliceLogoPi0WithPHOSPerformance(   Float_t startTextX,
                                                Float_t startTextY,
                                                Float_t startPi0TextX,
                                                Float_t differenceText,
                                                Float_t startLogoX,
                                                Float_t startLogoY,
                                                Float_t widthLogo,
                                                Float_t textSize,
                                                Float_t nEvents,
                                                TString collisionSystem,
                                                Bool_t mcFile ,
                                                Bool_t rawData,
                                                Bool_t pi0Label,
                                                Double_t xLengthCanvas,
                                                Double_t yLengthCanvas,
                                                TString date=""){

        TString aliceText   = "ALICE Performance";
        if(rawData)
            aliceText       = "RAW DATA";
        string processText;
        if(pi0Label){
            processText     = "#pi^{0}";
        } else {
            processText     = "#eta";
        }
        string eventsText;
        if(mcFile){
            eventsText      = "MC";
        } else {
            eventsText = "Data";
        }

        differenceText      = textSize + 3/yLengthCanvas;

        Double_t widthLogoPix   = xLengthCanvas*widthLogo;
        Double_t heightLogoPix  = widthLogoPix/0.73447;
        Double_t totalXLogo     = (startLogoX*xLengthCanvas + widthLogoPix)/xLengthCanvas;
        Double_t totalYLogo     = (startLogoY*yLengthCanvas + heightLogoPix)/yLengthCanvas;

        TLatex *alice           = new TLatex(startTextX,(startTextY+(2*differenceText)),Form("%s",aliceText.Data())); // Bo: this was modified
        TLatex *energy          = new TLatex(startTextX,(startTextY+differenceText),collisionSystem.Data()); // Bo: this was modified
        TLatex *process         = new TLatex(startPi0TextX+4.5*textSize, (startTextY+(2*differenceText)), Form("%s #rightarrow #gamma #gamma #rightarrow e^{+}e^{-}  e^{+}e^{-}",processText.c_str()));
        TLatex *process2        = new TLatex(startPi0TextX+4.5*textSize, (startTextY+differenceText), Form("%s #rightarrow #gamma #gamma",processText.c_str()));
        TLatex *process3        = new TLatex(startPi0TextX, (startTextY+differenceText),"PHOS/EMCAL:");
        TLatex *process4        = new TLatex(startPi0TextX, (startTextY+2*differenceText),"PCM:");
        TLatex *events          = 0x00;
        TLatex *latexDate       = 0x00;
        if (nEvents != 0){
            events              = new TLatex(startTextX,startTextY,Form("%s: %2.1e  MinBias events",eventsText.c_str(), nEvents)); // Bo: this was modified
            if (date.CompareTo("") != 0)
                latexDate       = new TLatex(startTextX,startTextY-differenceText,date.Data());
        } else {
            if (date.CompareTo("") != 0)
                latexDate       = new TLatex(startTextX,startTextY,date.Data());
        }

        TPad *myPadLogo = new TPad("myPadLogo", "Pad for ALICE Logo",startLogoX ,startLogoY ,totalXLogo,totalYLogo);

        alice->SetNDC();
        alice->SetTextColor(1);
        alice->SetTextSize(textSize);
        alice->Draw();

        energy->SetNDC();
        energy->SetTextColor(1);
        energy->SetTextSize(textSize);
        energy->Draw();

        if (nEvents != 0){
            events->SetNDC();
            events->SetTextColor(1);
            events->SetTextSize(textSize);
            events->Draw();
        }

        if (date.CompareTo("") != 0){
            latexDate->SetNDC();
            latexDate->SetTextColor(1);
            latexDate->SetTextSize(textSize);
            latexDate->Draw();
        }
        process->SetNDC(kTRUE); // <- use NDC coordinate
        process->SetTextSize(textSize);
        process->Draw();

        process2->SetNDC(kTRUE); // <- use NDC coordinate
        process2->SetTextSize(textSize);
        process2->Draw();

        process3->SetNDC(kTRUE); // <- use NDC coordinate
        process3->SetTextSize(textSize);
        process3->Draw();

        process4->SetNDC(kTRUE); // <- use NDC coordinate
        process4->SetTextSize(textSize);
        process4->Draw();


        //myPadLogo->SetFillStyle(2000); // color to first figure out where is the pad then comment !
        myPadLogo->SetBorderMode(0);
        myPadLogo->SetBorderSize(2);
        myPadLogo->SetFrameBorderMode(0);
        myPadLogo->SetLeftMargin(0.0);
        myPadLogo->SetTopMargin(0.0);
        myPadLogo->SetBottomMargin(0.0);
        myPadLogo->SetRightMargin(0.0);
        TImage *myAliceLogo = TImage::Open("ALICE_logo.eps");
        myPadLogo->Draw();  // to take out for not using a logo.
        myPadLogo->cd();
        myAliceLogo->Draw();

        return;
    }

    // Function added from old PlottingGammaConversionAdditional.h of TaskOmega
    void DrawAliceLogoOmega(Float_t startX, Float_t startY, Float_t widthLogo, Float_t textHeight, Double_t xLengthCanvas, Double_t yLengthCanvas,Double_t entries,TString addText){

        Double_t widthLogoPix = xLengthCanvas*widthLogo;
        Double_t heightLogoPix = widthLogoPix/0.73447;
        Double_t totalXLogo = (startX*xLengthCanvas + widthLogoPix)/xLengthCanvas;
        Double_t totalYLogo = (startY*yLengthCanvas + heightLogoPix)/yLengthCanvas;

        Float_t addTextStartY = startY - textHeight * 3.1;
        TLatex *addTextObj = new TLatex((startX-0.011),addTextStartY,addText.Data()); // Bo: this was modified
        addTextObj->SetNDC();
        addTextObj->SetTextColor(1);
        addTextObj->SetTextFont(62);
        addTextObj->SetTextSize(textHeight);
        addTextObj->SetLineWidth(2);
        addTextObj->Draw("same");
    }

    void DrawAliceLogoPi0WithPHOSOnlyPerformance(Float_t startTextX, Float_t startTextY, Float_t startPi0TextX, Float_t differenceText,
                            Float_t startLogoX, Float_t startLogoY, Float_t widthLogo,
                            Float_t textSize, Float_t nEvents,
                            TString collisionSystem, Bool_t mcFile , Bool_t rawData, Bool_t pi0Label, Double_t xLengthCanvas, Double_t yLengthCanvas, TString date=""){

        TString aliceText = "ALICE Performance";
        if(rawData) aliceText = "RAW DATA";
        string processText;
        if(pi0Label){processText = "#pi^{0}";}
        else { processText = "#eta";}
        string eventsText;
        if(mcFile){ eventsText = "MC";}
        else { eventsText = "Data";}

        differenceText = textSize + 3/yLengthCanvas;

        Double_t widthLogoPix = xLengthCanvas*widthLogo;
        Double_t heightLogoPix= widthLogoPix/0.73447;
        Double_t totalXLogo = (startLogoX*xLengthCanvas + widthLogoPix)/xLengthCanvas;
        Double_t totalYLogo = (startLogoY*yLengthCanvas + heightLogoPix)/yLengthCanvas;

        TLatex *alice = new TLatex(startTextX,(startTextY+(2*differenceText)),Form("%s",aliceText.Data())); // Bo: this was modified
        TLatex *energy = new TLatex(startTextX,(startTextY+differenceText),collisionSystem.Data()); // Bo: this was modified
        TLatex *process = new TLatex(startPi0TextX+2.5*textSize, (startTextY+(2*differenceText)), Form("%s #rightarrow #gamma #gamma #rightarrow e^{+}e^{-}  e^{+}e^{-}",processText.c_str()));
        TLatex *process2 = new TLatex(startPi0TextX+2.5*textSize, (startTextY+differenceText), Form("%s #rightarrow #gamma #gamma",processText.c_str()));
        TLatex *process3 = new TLatex(startPi0TextX, (startTextY+differenceText),"PHOS:");
        TLatex *process4 = new TLatex(startPi0TextX, (startTextY+2*differenceText),"PCM:");
        TLatex *events = 0x00;
            TLatex *latexDate = 0x00;
        if (nEvents != 0){
            events = new TLatex(startTextX,startTextY,Form("%s: %2.1e  MinBias events",eventsText.c_str(), nEvents)); // Bo: this was modified
            if (date.CompareTo("") != 0) latexDate = new TLatex(startTextX,startTextY-differenceText,date.Data());
        } else {
            if (date.CompareTo("") != 0) latexDate = new TLatex(startTextX,startTextY,date.Data());
        }

        TPad *myPadLogo = new TPad("myPadLogo", "Pad for ALICE Logo",startLogoX ,startLogoY ,totalXLogo,totalYLogo);

        alice->SetNDC();
        alice->SetTextColor(1);
        alice->SetTextSize(textSize);
        alice->Draw();

        energy->SetNDC();
        energy->SetTextColor(1);
        energy->SetTextSize(textSize);
        energy->Draw();

        if (nEvents != 0){
            events->SetNDC();
            events->SetTextColor(1);
            events->SetTextSize(textSize);
            events->Draw();
        }

        if (date.CompareTo("") != 0){
            latexDate->SetNDC();
            latexDate->SetTextColor(1);
            latexDate->SetTextSize(textSize);
            latexDate->Draw();
        }
        process->SetNDC(kTRUE); // <- use NDC coordinate
        process->SetTextSize(textSize);
        process->Draw();

        process2->SetNDC(kTRUE); // <- use NDC coordinate
        process2->SetTextSize(textSize);
        process2->Draw();

        process3->SetNDC(kTRUE); // <- use NDC coordinate
        process3->SetTextSize(textSize);
        process3->Draw();

        process4->SetNDC(kTRUE); // <- use NDC coordinate
        process4->SetTextSize(textSize);
        process4->Draw();


        //myPadLogo->SetFillStyle(2000); // color to first figure out where is the pad then comment !
            myPadLogo->SetBorderMode(0);
            myPadLogo->SetBorderSize(2);
            myPadLogo->SetFrameBorderMode(0);
            myPadLogo->SetLeftMargin(0.0);
            myPadLogo->SetTopMargin(0.0);
            myPadLogo->SetBottomMargin(0.0);
            myPadLogo->SetRightMargin(0.0);
            TImage *myAliceLogo = TImage::Open("ALICE_logo.eps");
            myPadLogo->Draw();  // to take out for not using a logo.
            myPadLogo->cd();
            myAliceLogo->Draw();

        return;
    }

    void DrawAliceLogoPi0WithPHOSPrelim(Float_t startTextX, Float_t startTextY, Float_t startPi0TextX, Float_t differenceText,
                                    Float_t startLogoX, Float_t startLogoY, Float_t widthLogo,
                                    Float_t textSize, Float_t nEvents,
                                    TString collisionSystem, Bool_t mcFile , Bool_t rawData, Bool_t pi0Label, Double_t xLengthCanvas, Double_t yLengthCanvas){
        cout <<  "ALICE work in progress" << endl;

        string aliceText = "ALICE Preliminary";
        if(rawData) aliceText = "RAW DATA";
        string processText;
        if(pi0Label){processText = "#pi^{0}";}
        else { processText = "#eta";}
        string eventsText;
        if(mcFile){ eventsText = "MC";}
        else { eventsText = "Data";}

        differenceText = textSize + 3/yLengthCanvas;

        Double_t widthLogoPix = xLengthCanvas*widthLogo;
        Double_t heightLogoPix= widthLogoPix/0.73447;
        Double_t totalXLogo = (startLogoX*xLengthCanvas + widthLogoPix)/xLengthCanvas;
        Double_t totalYLogo = (startLogoY*yLengthCanvas + heightLogoPix)/yLengthCanvas;

        TLatex *alice = new TLatex(startTextX,(startTextY+(2*differenceText)),Form("%s",aliceText.c_str())); // Bo: this was modified
        TLatex *energy = new TLatex(startTextX,(startTextY+differenceText),collisionSystem.Data()); // Bo: this was modified
        TLatex *process = new TLatex(startPi0TextX+3.5*textSize, (startTextY+(2*differenceText)), Form("%s #rightarrow #gamma #gamma #rightarrow e^{+}e^{-}  e^{+}e^{-}",processText.c_str()));
        TLatex *process2 = new TLatex(startPi0TextX+3.5*textSize, (startTextY+differenceText), Form("%s #rightarrow #gamma #gamma",processText.c_str()));
            //	TLatex *process3 = new TLatex(startPi0TextX, (startTextY+differenceText),"PHOS/EMCAL:");
        TLatex *process3 = new TLatex(startPi0TextX, (startTextY+differenceText),"PHOS:");
        TLatex *process4 = new TLatex(startPi0TextX, (startTextY+2*differenceText),"PCM:");
        TLatex *events = 0x00;
        if (nEvents != 0){
            events = new TLatex(startTextX,startTextY,Form("%s: %2.1e  MinBias events",eventsText.c_str(), nEvents)); // Bo: this was modified
        }
        TPad *myPadLogo = new TPad("myPadLogo", "Pad for ALICE Logo",startLogoX ,startLogoY ,totalXLogo,totalYLogo);

        alice->SetNDC();
        alice->SetTextColor(1);
        alice->SetTextSize(textSize);
        alice->Draw();

        energy->SetNDC();
        energy->SetTextColor(1);
        energy->SetTextSize(textSize);
        energy->Draw();

        if (nEvents != 0){
            events->SetNDC();
            events->SetTextColor(1);
            events->SetTextSize(textSize);
            events->Draw();
        }

        process->SetNDC(kTRUE); // <- use NDC coordinate
        process->SetTextSize(textSize);
        process->Draw();

        process2->SetNDC(kTRUE); // <- use NDC coordinate
        process2->SetTextSize(textSize);
        process2->Draw();

        process3->SetNDC(kTRUE); // <- use NDC coordinate
        process3->SetTextSize(textSize);
        process3->Draw();

        process4->SetNDC(kTRUE); // <- use NDC coordinate
        process4->SetTextSize(textSize);
        process4->Draw();


        //myPadLogo->SetFillStyle(2000); // color to first figure out where is the pad then comment !
        myPadLogo->SetBorderMode(0);
        myPadLogo->SetBorderSize(2);
        myPadLogo->SetFrameBorderMode(0);
        myPadLogo->SetLeftMargin(0.0);
        myPadLogo->SetTopMargin(0.0);
        myPadLogo->SetBottomMargin(0.0);
        myPadLogo->SetRightMargin(0.0);
        TImage *myAliceLogo = TImage::Open("ALICE_logo.eps");
        myPadLogo->Draw();  // to take out for not using a logo.
        myPadLogo->cd();
        myAliceLogo->Draw();

    }

    void DrawAliceLogoPi0WithPHOSFinal(Float_t startTextX, Float_t startTextY, Float_t textSize,
                                    TString collisionSystem, Bool_t pi0Label){
        cout <<  "ALICE work in progress" << endl;

        TString aliceText = "ALICE";
        TString processText;
        if(pi0Label){processText = "#pi^{0}";}
        else { processText = "#eta";}

        Float_t differenceText = 1.15*textSize;

        TLatex *alice = new TLatex(startTextX,(startTextY+differenceText),Form("%s",aliceText.Data())); // Bo: this was modified
        TLatex *energy = new TLatex(startTextX,(startTextY),Form("%s, %s", processText.Data(), collisionSystem.Data())); // Bo: this was modified

        alice->SetNDC();
        alice->SetTextColor(1);
        alice->SetTextSize(textSize);
        alice->Draw();

        energy->SetNDC();
        energy->SetTextColor(1);
        energy->SetTextSize(textSize);
        energy->Draw();
    }


    void DrawAliceLogoAllMesonsWithPHOSPrelim(Float_t startTextX, Float_t startTextY, Float_t startPi0TextX, Float_t differenceText,
                                    Float_t startLogoX, Float_t startLogoY, Float_t widthLogo,
                                    Float_t textSize, Float_t nEvents,
                                    TString collisionSystem, Bool_t mcFile , Bool_t rawData, Bool_t pi0Label, Double_t xLengthCanvas, Double_t yLengthCanvas){
        cout <<  "ALICE work in progress" << endl;

        string aliceText = "ALICE Preliminary";
        if(rawData) aliceText = "RAW DATA";
        string processTextPCM = "#pi^{0}, #eta";
        string processTextPHOS = "#pi^{0}, #eta";
        string eventsText;
        if(mcFile){ eventsText = "MC";}
        else { eventsText = "Data";}

        differenceText = textSize + 3/yLengthCanvas;

        Double_t widthLogoPix = xLengthCanvas*widthLogo;
        Double_t heightLogoPix= widthLogoPix/0.73447;
        Double_t totalXLogo = (startLogoX*xLengthCanvas + widthLogoPix)/xLengthCanvas;
        Double_t totalYLogo = (startLogoY*yLengthCanvas + heightLogoPix)/yLengthCanvas;

        TLatex *alice = new TLatex(startTextX,(startTextY+(2*differenceText)),Form("%s",aliceText.c_str())); // Bo: this was modified
        TLatex *energy = new TLatex(startTextX,(startTextY+differenceText),collisionSystem.Data()); // Bo: this was modified
        TLatex *process = new TLatex(startPi0TextX+3.5*textSize, (startTextY+(2*differenceText)), Form("%s #rightarrow #gamma #gamma #rightarrow e^{+}e^{-}  e^{+}e^{-}",processTextPCM.c_str()));
        TLatex *process2 = new TLatex(startPi0TextX+3.5*textSize, (startTextY+differenceText), Form("%s #rightarrow #gamma #gamma",processTextPHOS.c_str()));
        TLatex *process5 = new TLatex(startPi0TextX+3.5*textSize, (startTextY), "#omega #rightarrow #pi^{+}#pi^{-}#pi^{0}");

            //	TLatex *process3 = new TLatex(startPi0TextX, (startTextY+differenceText),"PHOS/EMCAL:");
        TLatex *process3 = new TLatex(startPi0TextX, (startTextY+differenceText),"PHOS:");
        TLatex *process4 = new TLatex(startPi0TextX, (startTextY+2*differenceText),"PCM:");
        TLatex *events = 0x00;
        if (nEvents != 0){
            events = new TLatex(startTextX,startTextY,Form("%s: %2.1e  MinBias events",eventsText.c_str(), nEvents)); // Bo: this was modified
        }
        TPad *myPadLogo = new TPad("myPadLogo", "Pad for ALICE Logo",startLogoX ,startLogoY ,totalXLogo,totalYLogo);

        alice->SetNDC();
        alice->SetTextColor(1);
        alice->SetTextSize(textSize);
        alice->Draw();

        energy->SetNDC();
        energy->SetTextColor(1);
        energy->SetTextSize(textSize);
        energy->Draw();

        if (nEvents != 0){
            events->SetNDC();
            events->SetTextColor(1);
            events->SetTextSize(textSize);
            events->Draw();
        }

        process->SetNDC(kTRUE); // <- use NDC coordinate
        process->SetTextSize(textSize);
        process->Draw();

        process2->SetNDC(kTRUE); // <- use NDC coordinate
        process2->SetTextSize(textSize);
        process2->Draw();

        process5->SetNDC(kTRUE); // <- use NDC coordinate
        process5->SetTextSize(textSize);
        process5->Draw();

        process3->SetNDC(kTRUE); // <- use NDC coordinate
        process3->SetTextSize(textSize);
        process3->Draw();

        process4->SetNDC(kTRUE); // <- use NDC coordinate
        process4->SetTextSize(textSize);
        process4->Draw();


        //myPadLogo->SetFillStyle(2000); // color to first figure out where is the pad then comment !
        myPadLogo->SetBorderMode(0);
        myPadLogo->SetBorderSize(2);
        myPadLogo->SetFrameBorderMode(0);
        myPadLogo->SetLeftMargin(0.0);
        myPadLogo->SetTopMargin(0.0);
        myPadLogo->SetBottomMargin(0.0);
        myPadLogo->SetRightMargin(0.0);
        TImage *myAliceLogo = TImage::Open("ALICE_logo.eps");
        myPadLogo->Draw();  // to take out for not using a logo.
        myPadLogo->cd();
        myAliceLogo->Draw();


    if(pi0Label) pi0Label = kFALSE;

    }


    void DrawAliceLogoPi0WithPHOSPi0EtaPrelim(Float_t startTextX, Float_t startTextY, Float_t startPi0TextX, Float_t differenceText,
                                Float_t startLogoX, Float_t startLogoY, Float_t widthLogo,
                                Float_t textSize, Float_t nEvents,
                                TString collisionSystem, Bool_t mcFile , Bool_t rawData, Bool_t pi0Label, Double_t xLengthCanvas, Double_t yLengthCanvas){

        if(pi0Label > 0)
            {}


        cout <<  "ALICE work in progress" << endl;

        string aliceText = "ALICE Preliminary";
        if(rawData) aliceText = "RAW DATA";
        string processText;
        processText = "#pi^{0} & #eta";
        string eventsText;
        if(mcFile){ eventsText = "MC";}
        else { eventsText = "Data";}

        Double_t widthLogoPix = xLengthCanvas*widthLogo;
        Double_t heightLogoPix= widthLogoPix/0.73447;
        Double_t totalXLogo = (startLogoX*xLengthCanvas + widthLogoPix)/xLengthCanvas;
        Double_t totalYLogo = (startLogoY*yLengthCanvas + heightLogoPix)/yLengthCanvas;

        differenceText = textSize + 3/yLengthCanvas;

        TLatex *alice = new TLatex(startTextX,(startTextY+(2*differenceText)),Form("%s",aliceText.c_str())); // Bo: this was modified
        TLatex *energy = new TLatex(startTextX,(startTextY+differenceText),collisionSystem.Data()); // Bo: this was modified
        TLatex *process = new TLatex(startPi0TextX+3.5*textSize, (startTextY+(2*differenceText)), Form("%s #rightarrow #gamma #gamma #rightarrow e^{+}e^{-}  e^{+}e^{-}",processText.c_str()));
        TLatex *process2 = new TLatex(startPi0TextX+3.5*textSize, (startTextY+differenceText), Form("%s #rightarrow #gamma #gamma",processText.c_str()));

        //	TLatex *process3 = new TLatex(startPi0TextX, (startTextY+differenceText),"PHOS/EMCAL:");
        TLatex *process3 = new TLatex(startPi0TextX, (startTextY+differenceText),"PHOS:");
        TLatex *process4 = new TLatex(startPi0TextX, (startTextY+2*differenceText),"PCM:");
        TLatex *events = 0x00;
        if (nEvents != 0){
            events = new TLatex(startTextX,startTextY,Form("%s: %2.1e  MinBias events",eventsText.c_str(), nEvents)); // Bo: this was modified
        }
        TPad *myPadLogo = new TPad("myPadLogo", "Pad for ALICE Logo",startLogoX ,startLogoY ,totalXLogo,totalYLogo);

        alice->SetNDC();
        alice->SetTextColor(1);
        alice->SetTextSize(textSize);
        alice->Draw();

        energy->SetNDC();
        energy->SetTextColor(1);
        energy->SetTextSize(textSize);
        energy->Draw();

        if (nEvents != 0){
            events->SetNDC();
            events->SetTextColor(1);
            events->SetTextSize(textSize);
            events->Draw();
        }

        process->SetNDC(kTRUE); // <- use NDC coordinate
        process->SetTextSize(textSize);
        process->Draw();

        process2->SetNDC(kTRUE); // <- use NDC coordinate
        process2->SetTextSize(textSize);
        process2->Draw();

        process3->SetNDC(kTRUE); // <- use NDC coordinate
        process3->SetTextSize(textSize);
        process3->Draw();

        process4->SetNDC(kTRUE); // <- use NDC coordinate
        process4->SetTextSize(textSize);
        process4->Draw();


        //myPadLogo->SetFillStyle(2000); // color to first figure out where is the pad then comment !
        myPadLogo->SetBorderMode(0);
        myPadLogo->SetBorderSize(2);
        myPadLogo->SetFrameBorderMode(0);
        myPadLogo->SetLeftMargin(0.0);
        myPadLogo->SetTopMargin(0.0);
        myPadLogo->SetBottomMargin(0.0);
        myPadLogo->SetRightMargin(0.0);
        TImage *myAliceLogo = TImage::Open("ALICE_logo.eps");
        myPadLogo->Draw();  // to take out for not using a logo.
        myPadLogo->cd();
        myAliceLogo->Draw();

    }

    void DrawAliceLogoPi0EtaPrelim(Float_t startTextX, Float_t startTextY, Float_t startPi0TextX, Float_t differenceText,
                                Float_t startLogoX, Float_t startLogoY, Float_t widthLogo,
                                Float_t textSize, Float_t nEvents,
                                TString collisionSystem, Bool_t mcFile , Bool_t rawData, Bool_t pi0Label, Double_t xLengthCanvas, Double_t yLengthCanvas){

    if(pi0Label){}

        string aliceText = "ALICE Preliminary";
        if(rawData) aliceText = "RAW DATA";
        string processText;
        processText = "#pi^{0} & #eta";
        string eventsText;
        if(mcFile){ eventsText = "MC";}
        else { eventsText = "Data";}

        Double_t widthLogoPix = xLengthCanvas*widthLogo;
        Double_t heightLogoPix= widthLogoPix/0.73447;
        Double_t totalXLogo = (startLogoX*xLengthCanvas + widthLogoPix)/xLengthCanvas;
        Double_t totalYLogo = (startLogoY*yLengthCanvas + heightLogoPix)/yLengthCanvas;

        differenceText = textSize + 3/yLengthCanvas;

        TLatex *alice = new TLatex(startTextX,(startTextY+(2*differenceText)),Form("%s",aliceText.c_str())); // Bo: this was modified
        TLatex *energy = new TLatex(startTextX,(startTextY+differenceText),collisionSystem.Data()); // Bo: this was modified
        TLatex *process = new TLatex(startPi0TextX, (startTextY+(2*differenceText)), Form("PCM: %s #rightarrow #gamma #gamma #rightarrow e^{+}e^{-}  e^{+}e^{-}",processText.c_str()));
        TLatex *events = 0x00;
        if (nEvents != 0){
            events = new TLatex(startTextX,startTextY,Form("%s: %2.1e  MinBias events",eventsText.c_str(), nEvents)); // Bo: this was modified
        }
        TPad *myPadLogo = new TPad("myPadLogo", "Pad for ALICE Logo",startLogoX ,startLogoY ,totalXLogo,totalYLogo);

        alice->SetNDC();
        alice->SetTextColor(1);
        alice->SetTextSize(textSize);
        alice->Draw();

        energy->SetNDC();
        energy->SetTextColor(1);
        energy->SetTextSize(textSize);
        energy->Draw();

        if (nEvents != 0){
            events->SetNDC();
            events->SetTextColor(1);
            events->SetTextSize(textSize);
            events->Draw();
        }

        process->SetNDC(kTRUE); // <- use NDC coordinate
        process->SetTextSize(textSize);
        process->Draw();

        //myPadLogo->SetFillStyle(2000); // color to first figure out where is the pad then comment !
        myPadLogo->SetBorderMode(0);
        myPadLogo->SetBorderSize(2);
        myPadLogo->SetFrameBorderMode(0);
        myPadLogo->SetLeftMargin(0.0);
        myPadLogo->SetTopMargin(0.0);
        myPadLogo->SetBottomMargin(0.0);
        myPadLogo->SetRightMargin(0.0);
        TImage *myAliceLogo = TImage::Open("ALICE_logo.eps");
        myPadLogo->Draw();  // to take out for not using a logo.
        myPadLogo->cd();
        myAliceLogo->Draw();

    }


    void DrawAliceLogoPi0Performance(Float_t startTextX, Float_t startTextY, Float_t startPi0TextX, Float_t differenceText,
                            Float_t startLogoX, Float_t startLogoY, Float_t widthLogo,
                            Float_t textSize, Float_t nEvents,
                            TString collisionSystem, Bool_t mcFile , Bool_t rawData, Bool_t pi0Label, Double_t xLengthCanvas, Double_t yLengthCanvas, TString date="",
                            TString centralityLabel = "MinBias", Bool_t dalitz =kFALSE){
        TString aliceText = "ALICE Performance";
        if(rawData) aliceText = "RAW DATA";
        string processText;
        string eventsText;
        if(pi0Label){processText = "#pi^{0}";
        } else { processText = "#eta";}
        if(mcFile){ eventsText = "MC";
        } else { eventsText = "Data";}

        TString decayChainText;
        if(dalitz){ decayChainText = Form("%s #rightarrow #gamma #gamma^{_{#scale[1.2]{*}}} #rightarrow e^{+}e^{-}  e^{+}e^{-}",processText.c_str());
        } else { decayChainText = Form("%s #rightarrow #gamma #gamma #rightarrow e^{+}e^{-}  e^{+}e^{-}",processText.c_str());}

        Double_t widthLogoPix = xLengthCanvas*widthLogo;
        Double_t heightLogoPix= widthLogoPix/0.73447;
        Double_t totalXLogo = (startLogoX*xLengthCanvas + widthLogoPix)/xLengthCanvas;
        Double_t totalYLogo = (startLogoY*yLengthCanvas + heightLogoPix)/yLengthCanvas;
        differenceText = textSize + 3/yLengthCanvas;

        TLatex *alice = new TLatex(startTextX,(startTextY+(2*differenceText)),Form("%s",aliceText.Data())); // Bo: this was modified
        TLatex *energy = new TLatex(startTextX,(startTextY+differenceText),collisionSystem.Data()); // Bo: this was modified
        TLatex *process = new TLatex(startPi0TextX, (startTextY+(2*differenceText)), decayChainText.Data());
        TLatex *events = 0x00;
            TLatex *latexDate = 0x00;
        if (nEvents != 0){
                events = new TLatex(startTextX,startTextY,Form("%s: %2.1e events, %s",eventsText.c_str(), nEvents, centralityLabel.Data())); // Bo: this was modified
                if (date.CompareTo("") != 0) latexDate = new TLatex(startTextX,startTextY-differenceText,date.Data());
        } else {
            if (date.CompareTo("") != 0) latexDate = new TLatex(startTextX,startTextY,date.Data());
        }

        TPad *myPadLogo = new TPad("myPadLogo", "Pad for ALICE Logo",startLogoX,startLogoY,totalXLogo,totalYLogo);

        alice->SetNDC();
        alice->SetTextColor(1);
        alice->SetTextSize(textSize);
        alice->Draw();

        energy->SetNDC();
        energy->SetTextColor(1);
        energy->SetTextSize(textSize);
        energy->Draw();

        process->SetNDC(kTRUE); // <- use NDC coordinate
        process->SetTextSize(textSize);
        process->Draw();

        if (nEvents != 0){
            events->SetNDC();
            events->SetTextColor(1);
            events->SetTextSize(textSize);
            events->Draw();
        }
        if (date.CompareTo("") != 0){
            latexDate->SetNDC();
            latexDate->SetTextColor(1);
            latexDate->SetTextSize(textSize);
            latexDate->Draw();
        }
        //myPadLogo->SetFillStyle(2000); // color to first figure out where is the pad then comment !
        myPadLogo->SetBorderMode(0);
        myPadLogo->SetBorderSize(2);
        myPadLogo->SetFrameBorderMode(0);
        myPadLogo->SetLeftMargin(0.0);
        myPadLogo->SetTopMargin(0.0);
        myPadLogo->SetBottomMargin(0.0);
        myPadLogo->SetRightMargin(0.0);
        TImage *myAliceLogo = TImage::Open("ALICE_logo.eps");
        myPadLogo->Draw();  // to take out for not using a logo.
        myPadLogo->cd();
        myAliceLogo->Draw();
    }

    void DrawAliceLogoPi0PerformanceExtract(Float_t startTextX, Float_t startTextY, Float_t startPi0TextX, Float_t differenceText,
                            Float_t startLogoX, Float_t startLogoY, Float_t widthLogo,
                            Float_t textSize, Float_t nEvents,
                            TString collisionSystem, Bool_t mcFile , Bool_t rawData, Bool_t pi0Label, Double_t xLengthCanvas, Double_t yLengthCanvas,TString date= "",Bool_t dalitz= kFALSE,
                            TString detectionChannel = ""){

    if(startPi0TextX){}

        TString aliceText = "ALICE Performance";
        if(rawData) aliceText = "RAW DATA";
        string processText;
        string eventsText;
        if(pi0Label){processText = "#pi^{0}";}
        else { processText = "#eta";}
        if(mcFile){ eventsText = "MC";}
        else { eventsText = "Data";}

        TString decayChainText;
        if(dalitz){ decayChainText = Form("%s #rightarrow #gamma #gamma^{_{#scale[1.2]{*}}}",processText.c_str());
        } else { decayChainText = Form("%s #rightarrow #gamma #gamma",processText.c_str());}


        Double_t widthLogoPix = xLengthCanvas*widthLogo;
        Double_t heightLogoPix= widthLogoPix/0.73447;
        Double_t totalXLogo = (startLogoX*xLengthCanvas + widthLogoPix)/xLengthCanvas;
        Double_t totalYLogo = (startLogoY*yLengthCanvas + heightLogoPix)/yLengthCanvas;

        differenceText = textSize + 3/yLengthCanvas;

        TLatex *alice = new TLatex(startTextX,(startTextY+(2*differenceText)),Form("%s",aliceText.Data())); // Bo: this was modified
        TLatex *latexDate = 0x00;
        if (date.CompareTo("") != 0) latexDate = new TLatex(startTextX,startTextY+differenceText,date.Data());
        TLatex *energy = new TLatex(startTextX,(startTextY),collisionSystem.Data()); // Bo: this was modified
        TLatex *process = new TLatex(startTextX, (startTextY-(1*differenceText)), decayChainText.Data());
        TLatex *latexDetPro = 0x00;
        if (detectionChannel.CompareTo("") != 0) latexDetPro = new TLatex(startTextX,startTextY-(2*differenceText),detectionChannel.Data());
        TLatex *events = 0x00;
        if (nEvents != 0){
            events = new TLatex(startTextX,startTextY-(3*differenceText),Form("%s: %2.1e events",eventsText.c_str(), nEvents)); // Bo: this was modified
        }

    // 	TPad *myPadLogo = new TPad("myPadLogo", "Pad for ALICE Logo",startLogoX,startLogoY,totalXLogo,totalYLogo);

        alice->SetNDC();
        alice->SetTextColor(1);
        alice->SetTextSize(textSize);
        alice->Draw();

        if (date.CompareTo("") != 0){
            latexDate->SetNDC();
            latexDate->SetTextColor(1);
            latexDate->SetTextSize(textSize);
            latexDate->Draw();
        }

        energy->SetNDC();
        energy->SetTextColor(1);
        energy->SetTextSize(textSize);
        energy->Draw();

        process->SetNDC(kTRUE); // <- use NDC coordinate
        process->SetTextColor(1);
        process->SetTextSize(textSize);
        process->Draw();

        if (detectionChannel.CompareTo("") != 0){
            latexDetPro->SetNDC();
            latexDetPro->SetTextColor(1);
            latexDetPro->SetTextSize(textSize);
            latexDetPro->Draw();
        }


        if (nEvents != 0){
            events->SetNDC();
            events->SetTextColor(1);
            events->SetTextSize(textSize);
            events->Draw();
        }
    // 	//myPadLogo->SetFillStyle(2000); // color to first figure out where is the pad then comment !
    // 	myPadLogo->SetBorderMode(0);
    // 	myPadLogo->SetBorderSize(2);
    // 	myPadLogo->SetFrameBorderMode(0);
    // 	myPadLogo->SetLeftMargin(0.0);
    // 	myPadLogo->SetTopMargin(0.0);
    // 	myPadLogo->SetBottomMargin(0.0);
    // 	myPadLogo->SetRightMargin(0.0);
    // 	TImage *myAliceLogo = TImage::Open("ALICE_logo.eps");
    // 	myPadLogo->Draw();  // to take out for not using a logo.
    // 	myPadLogo->cd();
    // 	myAliceLogo->Draw();
    }





    void DrawAliceLogoPi0WorkInProgress(Float_t startTextX, Float_t startTextY, Float_t startPi0TextX, Float_t differenceText,
                                Float_t startLogoX, Float_t startLogoY, Float_t widthLogo,
                                Float_t textSize, Float_t nEvents,
                                TString collisionSystem, Bool_t mcFile , Bool_t rawData, Bool_t pi0Label, Double_t xLengthCanvas, Double_t yLengthCanvas, Bool_t dalitz = kFALSE, TString stringCentrality=""){

        string aliceText = "ALICE";
        if(rawData) aliceText = "RAW DATA";
        string processText;
        string eventsText;
        if(pi0Label){processText = "#pi^{0}";}
        else { processText = "#eta";}
        if(mcFile){ eventsText = "MC";}
        else { eventsText = "Data";}
        TString decayChainText;
        if(dalitz){ decayChainText = Form("%s #rightarrow #gamma #gamma^{_{#scale[1.2]{*}}}",processText.c_str());
        } else { decayChainText = Form("%s #rightarrow #gamma #gamma",processText.c_str());}


        differenceText = textSize + 3/yLengthCanvas;

        Double_t widthLogoPix = xLengthCanvas*widthLogo;
        Double_t heightLogoPix= widthLogoPix/0.73447;
        Double_t totalXLogo = (startLogoX*xLengthCanvas + widthLogoPix)/xLengthCanvas;
        Double_t totalYLogo = (startLogoY*yLengthCanvas + heightLogoPix)/yLengthCanvas;
        TLatex *alice = new TLatex(startTextX,(startTextY+(2*differenceText)),Form("%s",aliceText.c_str())); // Bo: this was modified
        TLatex *energy = new TLatex(startTextX,(startTextY+differenceText),collisionSystem.Data()); // Bo: this was modified
        TLatex *process = new TLatex(startPi0TextX, (startTextY+(2*differenceText)), decayChainText.Data());
        TLatex *events = 0x00;
        if (nEvents != 0 && stringCentrality.CompareTo("") == 0){
            events = new TLatex(startTextX,startTextY,Form("%s: %2.1e  MinBias events",eventsText.c_str(), nEvents)); // Bo: this was modified
        } else if (nEvents != 0 ){
            events = new TLatex(startTextX,startTextY,Form("%s: %2.1e events %s",eventsText.c_str(), nEvents, stringCentrality.Data())); // Bo: this was modified
        }
        TPad *myPadLogo = new TPad("myPadLogo", "Pad for ALICE Logo",startLogoX,startLogoY,totalXLogo,totalYLogo);

        alice->SetNDC();
        alice->SetTextColor(1);
        alice->SetTextSize(textSize);
        alice->Draw();

        energy->SetNDC();
        energy->SetTextColor(1);
        energy->SetTextSize(textSize);
        energy->Draw();

        process->SetNDC(kTRUE); // <- use NDC coordinate
        process->SetTextSize(textSize);
        process->Draw();

        if (nEvents != 0){
            events->SetNDC();
            events->SetTextColor(1);
            events->SetTextSize(textSize);
            events->Draw();
        }
        //myPadLogo->SetFillStyle(2000); // color to first figure out where is the pad then comment !
        myPadLogo->SetBorderMode(0);
        myPadLogo->SetBorderSize(2);
        myPadLogo->SetFrameBorderMode(0);
        myPadLogo->SetLeftMargin(0.0);
        myPadLogo->SetTopMargin(0.0);
        myPadLogo->SetBottomMargin(0.0);
        myPadLogo->SetRightMargin(0.0);
        TImage *myAliceLogo = TImage::Open("ALICE_logo.eps");
        myPadLogo->Draw();  // to take out for not using a logo.
        myPadLogo->cd();
        //myAliceLogo->Draw();
    }



    /*********************************************************************************************************
    DrawAliceLogoPi0Performance
    * will draw you the ALICE Logo as well the text "ALICE Performance", "pp @ 7 TeV" and the date
    which you can hand over
    Float_t startTextX, Float_t startTextY, Float_t startPi0TextX, Float_t differenceText,
    Float_t startLogoX, Float_t startLogoY, Float_t widthLogo,
    Float_t textSize,
    Bool_t collisionSystem, Bool_t mcFile
    * float_t startX, float_t startY - give starting Point of Logo
    * float_t widthLogo - gives the width of the logo
    * float_t textHeight - gives you the heigth of the text
    * char * date - date handed over
    **********************************************************************************************************
    **********************************************************************************************************/

    void DrawAliceLogoCombined(Float_t startTextX, Float_t startTextY, Float_t startPi0TextX, Float_t differenceText,
                            Float_t startLogoX, Float_t startLogoY, Float_t widthLogo,
                            Float_t textSize, Float_t nEvents,
                        TString collisionSystem, Bool_t mcFile , Bool_t rawData, Bool_t pi0Label, Double_t xLengthCanvas, Double_t yLengthCanvas, TString date ="", Bool_t dalitz = kFALSE){

    if(pi0Label){}

        TString aliceText = "ALICE Performance";
        if(rawData) aliceText = "RAW DATA";

        string processText = "#pi^{0} and #eta";
        string eventsText;
        if(mcFile){ eventsText = "MC";}
        else { eventsText = "Data";}

        TString decayChainText;
        if(dalitz){ decayChainText = Form("%s #rightarrow #gamma #gamma^{_{#scale[1.2]{*}}} #rightarrow e^{+}e^{-}  e^{+}e^{-}",processText.c_str());
        } else { decayChainText = Form("%s #rightarrow #gamma #gamma #rightarrow e^{+}e^{-}  e^{+}e^{-}",processText.c_str());}

        Double_t widthLogoPix = xLengthCanvas*widthLogo;
        Double_t heightLogoPix= widthLogoPix/0.73447;
        Double_t totalXLogo = (startLogoX*xLengthCanvas + widthLogoPix)/xLengthCanvas;
        Double_t totalYLogo = (startLogoY*yLengthCanvas + heightLogoPix)/yLengthCanvas;

        differenceText = textSize + 3/yLengthCanvas;

        TLatex *alice = new TLatex(startTextX,(startTextY+(2*differenceText)),Form("%s",aliceText.Data())); // Bo: this was modified
        TLatex *energy = new TLatex(startTextX,(startTextY+differenceText),collisionSystem.Data()); // Bo: this was modified
        TLatex *process = new TLatex(startPi0TextX, (startTextY+(2*differenceText)),decayChainText.Data());
        TLatex *events = 0x00;
        TLatex *latexDate = 0x00;
        if (nEvents != 0){
            events = new TLatex(startTextX,startTextY,Form("%s: %2.1e  MinBias events",eventsText.c_str(), nEvents)); // Bo: this was modified
            if (date.CompareTo("") != 0) latexDate = new TLatex(startTextX,startTextY-differenceText,date.Data());
        } else {
            if (date.CompareTo("") != 0) latexDate = new TLatex(startTextX,startTextY,date.Data());
        }

        TPad *myPadLogo = new TPad("myPadLogo", "Pad for ALICE Logo",startLogoX,startLogoY,totalXLogo,totalYLogo);


        alice->SetNDC();
        alice->SetTextColor(1);
        alice->SetTextSize(textSize);
        alice->Draw();

        energy->SetNDC();
        energy->SetTextColor(1);
        energy->SetTextSize(textSize);
        energy->Draw();

        process->SetNDC(kTRUE); // <- use NDC coordinate
        process->SetTextSize(textSize);
        process->Draw();

        if (nEvents != 0){
            events->SetNDC();
            events->SetTextColor(1);
            events->SetTextSize(textSize);
            events->Draw();
        }
        if (date.CompareTo("") != 0){
            latexDate->SetNDC();
            latexDate->SetTextColor(1);
            latexDate->SetTextSize(textSize);
            latexDate->Draw();
        }
        //myPadLogo->SetFillStyle(2000); // color to first figure out where is the pad then comment !
        myPadLogo->SetBorderMode(0);
        myPadLogo->SetBorderSize(2);
        myPadLogo->SetFrameBorderMode(0);
        myPadLogo->SetLeftMargin(0.0);
        myPadLogo->SetTopMargin(0.0);
        myPadLogo->SetBottomMargin(0.0);
        myPadLogo->SetRightMargin(0.0);
        TImage *myAliceLogo = TImage::Open("ALICE_logo.eps");
        myPadLogo->Draw();  // to take out for not using a logo.
        myPadLogo->cd();
        myAliceLogo->Draw();
    }


    void DrawAliceLogoCombinedPreliminary(Float_t startTextX, Float_t startTextY, Float_t startPi0TextX, Float_t differenceText,
                                    Float_t startLogoX, Float_t startLogoY, Float_t widthLogo,
                                    Float_t textSize, Float_t nEvents,
                                TString collisionSystem, Bool_t mcFile , Bool_t rawData, Bool_t pi0Label, Double_t xLengthCanvas, Double_t yLengthCanvas, Bool_t dalitz = kFALSE){

    if(pi0Label){}


        string aliceText = "ALICE Preliminary";
        if(rawData) aliceText = "RAW DATA";
        string processText = "#pi^{0} and #eta";
        string eventsText;
        if(mcFile){ eventsText = "MC";}
        else { eventsText = "Data";}

        TString decayChainText;
        if(dalitz){ decayChainText = Form("%s #rightarrow #gamma #gamma^{_{#scale[1.2]{*}}} #rightarrow e^{+}e^{-}  e^{+}e^{-}",processText.c_str());
        } else { decayChainText = Form("%s #rightarrow #gamma #gamma (#rightarrow e^{+}e^{-}  e^{+}e^{-})",processText.c_str());}

        Double_t widthLogoPix = xLengthCanvas*widthLogo;
        Double_t heightLogoPix= widthLogoPix/0.73447;
        Double_t totalXLogo = (startLogoX*xLengthCanvas + widthLogoPix)/xLengthCanvas;
        Double_t totalYLogo = (startLogoY*yLengthCanvas + heightLogoPix)/yLengthCanvas;
        differenceText = textSize + 3/yLengthCanvas;

        TLatex *alice = new TLatex(startTextX,(startTextY+(2*differenceText)),Form("%s",aliceText.c_str())); // Bo: this was modified
        TLatex *energy = new TLatex(startTextX,(startTextY+differenceText),collisionSystem.Data()); // Bo: this was modified
        TLatex *process = new TLatex(startPi0TextX, (startTextY+(2*differenceText)), decayChainText.Data());
        TLatex *events =0x00;
        if (nEvents != 0){
            events = new TLatex(startTextX,startTextY,Form("%s: %2.1e  MinBias events",eventsText.c_str(), nEvents)); // Bo: this was modified
        }
        TPad *myPadLogo = new TPad("myPadLogo", "Pad for ALICE Logo",startLogoX,startLogoY,totalXLogo,totalYLogo);

        alice->SetNDC();
        alice->SetTextColor(1);
        alice->SetTextSize(textSize);
        alice->Draw();

        energy->SetNDC();
        energy->SetTextColor(1);
        energy->SetTextSize(textSize);
        energy->Draw();

        process->SetNDC(kTRUE); // <- use NDC coordinate
        process->SetTextSize(textSize);
        process->Draw();

        if (nEvents != 0){
            events->SetNDC();
            events->SetTextColor(1);
            events->SetTextSize(textSize);
            events->Draw();
        }
        //myPadLogo->SetFillStyle(2000); // color to first figure out where is the pad then comment !
        myPadLogo->SetBorderMode(0);
        myPadLogo->SetBorderSize(2);
        myPadLogo->SetFrameBorderMode(0);
        myPadLogo->SetLeftMargin(0.0);
        myPadLogo->SetTopMargin(0.0);
        myPadLogo->SetBottomMargin(0.0);
        myPadLogo->SetRightMargin(0.0);
        TImage *myAliceLogo = TImage::Open("ALICE_logo.eps");
        myPadLogo->Draw();  // to take out for not using a logo.
        myPadLogo->cd();
        myAliceLogo->Draw();
    }


    /*************************************************************************************************
    DrawAliceLogo draws you the Alice logo + "work in progress" and " ALICE Performance"
    be careful you have to set the path of the alice-logo for 	your system
    * float_t startX, float_t startY - give starting Point of Logo
    * float_t widthLogo - gives the width of the logo
    * float_t textHeight - gives you the heigth of the text
    **************************************************************************************************
    **************************************************************************************************/

    void DrawAliceLogo(Float_t startX, Float_t startY, Float_t widthLogo, Float_t textHeight, Double_t xLengthCanvas, Double_t yLengthCanvas){

        Float_t aliceStartY = startY - textHeight * 1.1;
        TLatex *alice = new TLatex(startX-0.0225,aliceStartY,"ALICE Performance"); // Bo: this was modified
        alice->SetNDC();
        alice->SetTextColor(1);
        alice->SetTextFont(42);
        alice->SetTextSize(textHeight);
        alice->SetLineWidth(2);
        alice->Draw("same");


        Double_t widthLogoPix = xLengthCanvas*widthLogo;
        Double_t heightLogoPix = widthLogoPix/0.73447;
        Double_t totalXLogo = (startX*xLengthCanvas + widthLogoPix)/xLengthCanvas;
        Double_t totalYLogo = (startY*yLengthCanvas + heightLogoPix)/yLengthCanvas;


        Float_t wipStartY = startY - textHeight *(1 + 1.1);
        TLatex *wip = new TLatex((startX-0.011),wipStartY,"work in progress"); // Bo: this was modified
        wip->SetNDC();
        wip->SetTextColor(1);
        wip->SetTextFont(51);
        wip->SetTextSize(textHeight);
        wip->SetLineWidth(2);
        wip->Draw("same");

        TPad *myPadLogo = new TPad("myPadLogo", "Pad for ALICE Logo",startX,startY,totalXLogo,totalYLogo);
        //  myPadLogo->SetFillColor(2); // color to first figure out where is the pad then comment !
        myPadLogo->SetBorderMode(0);
        myPadLogo->SetBorderSize(2);
        myPadLogo->SetFrameBorderMode(0);
        myPadLogo->SetLeftMargin(0.0);
        myPadLogo->SetTopMargin(0.0);
        myPadLogo->SetBottomMargin(0.0);
        myPadLogo->SetRightMargin(0.0);
        TImage *myAliceLogo = TImage::Open("ALICE_logo.eps");
        myPadLogo->Draw();  // to take out for not using a logo.
        myPadLogo->cd();
        myAliceLogo->Draw("same");

    }

    /*************************************************************************************************
    DrawAliceLogo1D draws you the Alice logo + "work in progress" and " ALICE Performance"
    be careful you have to set the path of the alice-logo for 	your system
    * float_t startX, float_t startY - give starting Point of Logo
    * float_t widthLogo - gives the width of the logo
    * float_t textHeight - gives you the heigth of the text
    **************************************************************************************************
    **************************************************************************************************/

    void DrawAliceLogo1D(Float_t startX, Float_t startY, Float_t widthLogo, Float_t textHeight, Double_t xLengthCanvas, Double_t yLengthCanvas){

        Float_t aliceStartY = startY - textHeight * 1.1;
        TLatex *alice = new TLatex(startX+0.008,aliceStartY,"ALICE Performance"); // Bo: this was modified
        alice->SetNDC();
        alice->SetTextColor(1);
        alice->SetTextFont(42);
        alice->SetTextSize(textHeight);
        alice->SetLineWidth(2);
        alice->Draw("same");

        Double_t widthLogoPix = xLengthCanvas*widthLogo;
        Double_t heightLogoPix = widthLogoPix/0.73447;
        Double_t totalXLogo = (startX*xLengthCanvas + widthLogoPix)/xLengthCanvas;
        Double_t totalYLogo = (startY*yLengthCanvas + heightLogoPix)/yLengthCanvas;

        Float_t wipStartY = startY - textHeight *(1 + 1.1);
        TLatex *wip = new TLatex((startX+0.022),wipStartY,"work in progress"); // Bo: this was modified
        wip->SetNDC();
        wip->SetTextColor(1);
        wip->SetTextFont(51);
        wip->SetTextSize(textHeight);
        wip->SetLineWidth(2);
        wip->Draw("same");

        TPad *myPadLogo = new TPad("myPadLogo", "Pad for ALICE Logo",startX,startY,totalXLogo,totalYLogo);
        //  myPadLogo->SetFillColor(2); // color to first figure out where is the pad then comment !
        myPadLogo->SetBorderMode(0);
        myPadLogo->SetBorderSize(2);
        myPadLogo->SetFrameBorderMode(0);
        myPadLogo->SetLeftMargin(0.0);
        myPadLogo->SetTopMargin(0.0);
        myPadLogo->SetBottomMargin(0.0);
        myPadLogo->SetRightMargin(0.0);
        TImage *myAliceLogo = TImage::Open("ALICE_logo.eps");
        myPadLogo->Draw();  // to take out for not using a logo.
        myPadLogo->cd();
        myAliceLogo->Draw("same");
    }


    /*********************************************************************************************************
    DrawAliceLogoPerformance
    * will draw you the ALICE Logo as well the text "ALICE Performance", "pp @ 7 TeV" and the date
    which you can hand over
    * float_t startX, float_t startY - give starting Point of Logo
    * float_t widthLogo - gives the width of the logo
    * float_t textHeight - gives you the heigth of the text
    * char * date - date handed over
    **********************************************************************************************************
    **********************************************************************************************************/

    void DrawAliceLogoPerformance(Float_t startX, Float_t startY, Float_t widthLogo, Float_t textHeight, Float_t decrease, TString date, TString collisionSystem, TString textGenerator, TString textPeriod, Double_t xLengthCanvas, Double_t yLengthCanvas){

        TString ALICEPerform = "ALICE Performance";
        Float_t aliceStartY = startY - textHeight * 1.15;
        TLatex *alice = new TLatex((startX-decrease),aliceStartY,ALICEPerform); // Bo: this was modified
        alice->SetNDC();
        alice->SetTextColor(1);
        alice->SetTextFont(62);
        alice->SetTextSize(textHeight);
        alice->SetLineWidth(2);
        alice->Draw("same");
        TLatex *pp7 = NULL;
        if( collisionSystem.CompareTo("PbPb @ #sqrt{#it{s}_{_{NN}}} = 2.76 TeV") == 0){
            pp7 = new TLatex((startX-2*decrease),(aliceStartY-textHeight*1.15),collisionSystem.Data()); // Bo: this was modified
        } else {
            pp7 = new TLatex((startX+2*decrease),(aliceStartY-textHeight*1.15),collisionSystem.Data()); // Bo: this was modified
        }
        pp7->SetNDC();
        pp7->SetTextColor(1);
        pp7->SetTextFont(62);
        pp7->SetTextSize(textHeight);
        pp7->SetLineWidth(2);
        pp7->Draw("same");
        TLatex *today = new TLatex((startX+2*decrease),(aliceStartY-2*textHeight*1.15),date.Data()); // Bo: this was modified
        today->SetNDC();
        today->SetTextColor(1);
        today->SetTextFont(62);
        today->SetTextSize(textHeight);
        today->SetLineWidth(2);
        today->Draw("same");
        if (textGenerator.CompareTo("")!=0 && textPeriod.CompareTo("")!=0){
            TLatex *generator = new TLatex((startX+decrease),(aliceStartY-3*textHeight*1.15),Form("%s   %s",textGenerator.Data(),textPeriod.Data())); // Bo: this was modified
            generator->SetNDC();
            generator->SetTextColor(1);
            generator->SetTextFont(62);
            generator->SetTextSize(textHeight);
            generator->SetLineWidth(2);
            generator->Draw("same");
        } else if (textGenerator.CompareTo("")!=0) {
            TLatex *generator = new TLatex((startX+decrease),(aliceStartY-3*textHeight*1.15),Form("%s",textGenerator.Data())); // Bo: this was modified
            generator->SetNDC();
            generator->SetTextColor(1);
            generator->SetTextFont(62);
            generator->SetTextSize(textHeight);
            generator->SetLineWidth(2);
            generator->Draw("same");
        }
        Double_t widthLogoPix = xLengthCanvas*widthLogo;
        Double_t heightLogoPix = widthLogoPix/0.73447;
        Double_t totalXLogo = (startX*xLengthCanvas + widthLogoPix)/xLengthCanvas;
        Double_t totalYLogo = (startY*yLengthCanvas + heightLogoPix)/yLengthCanvas;

            TPad *myPadLogo = new TPad("myPadLogo", "Pad for ALICE Logo",startX,startY,totalXLogo,totalYLogo);
        //  myPadLogo->SetFillColor(2); // color to first figure out where is the pad then comment !
        myPadLogo->SetBorderMode(0);
        myPadLogo->SetBorderSize(2);
        myPadLogo->SetFrameBorderMode(0);
        myPadLogo->SetLeftMargin(0.0);
        myPadLogo->SetTopMargin(0.0);
        myPadLogo->SetBottomMargin(0.0);
        myPadLogo->SetRightMargin(0.0);
        TImage *myAliceLogo = TImage::Open("ALICE_logo.eps");
        myPadLogo->Draw();  // to take out for not using a logo.
        myPadLogo->cd();
        myAliceLogo->Draw("same");

    }

    void DrawLabelsEvents(Float_t startX, Float_t startY, Float_t textHeight, Float_t decrease,  TString collisionSystem, TString textGenerator, TString textPeriod){

        Float_t aliceStartY = startY - textHeight * 1.15;
        TLatex *pp7 = NULL;
        if( collisionSystem.CompareTo("PbPb @ #sqrt{#it{s}_{_{NN}}} = 2.76 TeV") == 0){
            pp7 = new TLatex((startX-2*decrease),(aliceStartY),collisionSystem.Data()); // Bo: this was modified
        } else {
            pp7 = new TLatex((startX+2*decrease),(aliceStartY),collisionSystem.Data()); // Bo: this was modified
        }
        pp7->SetNDC();
        pp7->SetTextColor(1);
        pp7->SetTextFont(62);
        pp7->SetTextSize(textHeight);
        pp7->SetLineWidth(2);
        pp7->Draw("same");
        if (textGenerator.CompareTo("")!=0 && textPeriod.CompareTo("")!=0){
            TLatex *generator = new TLatex((startX+decrease),(aliceStartY-1*textHeight*1.15),Form("%s   %s",textGenerator.Data(),textPeriod.Data())); // Bo: this was modified
            generator->SetNDC();
            generator->SetTextColor(1);
            generator->SetTextFont(62);
            generator->SetTextSize(textHeight);
            generator->SetLineWidth(2);
            generator->Draw("same");
        } else if (textGenerator.CompareTo("")!=0) {
            TLatex *generator = new TLatex((startX+decrease),(aliceStartY-1*textHeight*1.15),Form("%s",textGenerator.Data())); // Bo: this was modified
            generator->SetNDC();
            generator->SetTextColor(1);
            generator->SetTextFont(62);
            generator->SetTextSize(textHeight);
            generator->SetLineWidth(2);
            generator->Draw("same");
        }
    }

    /*********************************************************************************************************
    DrawAliceLogoPerformance2D
    * will draw you the ALICE Logo as well the text "ALICE Performance", "pp @ 7 TeV" and the date
    which you can hand over
    * float_t startX, float_t startY - give starting Point of Logo
    * float_t widthLogo - gives the width of the logo
    * float_t textHeight - gives you the heigth of the text
    * float_t decrease - gives percentage in the canvas to which the text should be decreased in y
    compared to the logo
    * char * date - date handed over
    **********************************************************************************************************
    **********************************************************************************************************/

    void DrawAliceLogoPerformance2D(Float_t startX, Float_t startY, Float_t widthLogo, Float_t textHeight, Float_t decrease, TString date, TString collisionSystem, TString textGenerator, TString textPeriod, Double_t xLengthCanvas, Double_t yLengthCanvas){


        Double_t widthLogoPix = xLengthCanvas*widthLogo;
        Double_t heightLogoPix = widthLogoPix/0.73447;
        Double_t totalXLogo = (startX*xLengthCanvas + widthLogoPix)/xLengthCanvas;
        Double_t totalYLogo = (startY*yLengthCanvas - heightLogoPix)/yLengthCanvas;

        Float_t aliceStartY = totalYLogo - 0.03;
        TLatex *alice = new TLatex((startX-2*decrease),aliceStartY,"ALICE Performance"); // Bo: this was modified
        alice->SetNDC();
        alice->SetTextColor(1);
        alice->SetTextFont(62);
        alice->SetTextSize(textHeight);
        alice->SetLineWidth(2);
        alice->Draw("same");
        TLatex *pp7 = new TLatex((startX-decrease),(aliceStartY-textHeight*1.1),collisionSystem.Data()); // Bo: this was modified
        pp7->SetNDC();
        pp7->SetTextColor(1);
        pp7->SetTextFont(62);
        pp7->SetTextSize(textHeight);
        pp7->SetLineWidth(2);
        pp7->Draw("same");
        TLatex *today = new TLatex((startX-decrease),(aliceStartY-2*textHeight*1.1),date.Data()); // Bo: this was modified
        today->SetNDC();
        today->SetTextColor(1);
        today->SetTextFont(62);
        today->SetTextSize(textHeight);
        today->SetLineWidth(2);
        today->Draw("same");
        if (textGenerator.CompareTo("")!=0 && textPeriod.CompareTo("")!=0){
            TLatex *generator = new TLatex((startX-2*decrease),(aliceStartY-3*textHeight*1.1),Form("%s   %s",textGenerator.Data(),textPeriod.Data())); // Bo: this was modified
            generator->SetNDC();
            generator->SetTextColor(1);
            generator->SetTextFont(62);
            generator->SetTextSize(textHeight);
            generator->SetLineWidth(2);
            generator->Draw("same");
        } else if (textGenerator.CompareTo("")!=0) {
            TLatex *generator = new TLatex((startX-2*decrease),(aliceStartY-3*textHeight*1.1),Form("%s",textGenerator.Data())); // Bo: this was modified
            generator->SetNDC();
            generator->SetTextColor(1);
            generator->SetTextFont(62);
            generator->SetTextSize(textHeight);
            generator->SetLineWidth(2);
            generator->Draw("same");
        }


        TPad *myPadLogo = new TPad("myPadLogo", "Pad for ALICE Logo",startX , totalYLogo, totalXLogo,startY);
        //  myPadLogo->SetFillColor(2); // color to first figure out where is the pad then comment !
        myPadLogo->SetBorderMode(0);
        myPadLogo->SetBorderSize(2);
        myPadLogo->SetFrameBorderMode(0);
        myPadLogo->SetLeftMargin(0.0);
        myPadLogo->SetTopMargin(0.0);
        myPadLogo->SetBottomMargin(0.0);
        myPadLogo->SetRightMargin(0.0);
        TImage *myAliceLogo = TImage::Open("ALICE_logo.eps");
        myPadLogo->Draw();  // to take out for not using a logo.
        myPadLogo->cd();
        myAliceLogo->Draw("same");

    }


    //***********************************************************************************************************
    //* DrawAliceText draws you the "work in progress" and " ALICE Performance"
    //* float_t startX, float_t startY - give starting Point of Logo
    //* float_t textHeight - gives you the heigth of the text
    //***********************************************************************************************************
    //***********************************************************************************************************

    void DrawAliceText(Float_t startX, Float_t startY, Float_t textHeight){
        Float_t aliceStartY = startY - textHeight * 1.1;
        TLatex *alice = new TLatex(startX,aliceStartY,"ALICE Performance"); // Bo: this was modified
        alice->SetNDC();
        alice->SetTextColor(1);
        alice->SetTextFont(42);
        alice->SetTextSize(textHeight);
        alice->SetLineWidth(2);
        alice->Draw();

        Float_t wipStartY = startY - textHeight *(1 + 1.1);
        TLatex *wip = new TLatex(startX,wipStartY,"work in progress"); // Bo: this was modified
        wip->SetNDC();
        wip->SetTextColor(1);
        wip->SetTextFont(51);
        wip->SetTextSize(textHeight);
        wip->SetLineWidth(2);
        wip->Draw();
    }

    /*****************************************************************************************************
    DrawStructure() draws the structure of the Inner Alice Detectors labeled for a xy - Plot
    ******************************************************************************************************
    ******************************************************************************************************/

    void DrawStructure(){
        TLatex *ssdText = new TLatex(0.16,0.87,"SSD");
        ssdText->SetNDC();
        ssdText->SetTextFont(72);
        ssdText->SetTextSize(0.03);
        ssdText->SetLineWidth(4);
        ssdText->Draw();

        TLatex *sddText = new TLatex(0.14,0.78,"SDD");
        sddText->SetNDC();
        sddText->SetTextFont(72);
        sddText->SetTextSize(0.03);
        sddText->SetLineWidth(4);
        sddText->Draw();

        TLatex *spdText = new TLatex(0.14,0.27,"SPD");
        spdText->SetNDC();
        spdText->SetTextFont(72);
        spdText->SetTextSize(0.03);
        spdText->SetLineWidth(4);
        spdText->Draw();

        TLatex *tpcRText = new TLatex(0.52,0.095,"TPC Rods");
        tpcRText->SetNDC();
        tpcRText->SetTextFont(72);
        tpcRText->SetTextSize(0.03);
        tpcRText->SetLineWidth(4);
        tpcRText->Draw();

        TLatex *tpcIText = new TLatex(0.14,0.18,"TPC inner");
        tpcIText->SetNDC();
        tpcIText->SetTextFont(72);
        tpcIText->SetTextSize(0.03);
        tpcIText->SetLineWidth(4);
        tpcIText->Draw();

        TLatex *tpcIFText = new TLatex(0.14,0.15,"field cage");
        tpcIFText->SetNDC();
        tpcIFText->SetTextFont(72);
        tpcIFText->SetTextSize(0.03);
        tpcIFText->SetLineWidth(4);
        tpcIFText->Draw();

        TLatex *tpcIFVText = new TLatex(0.16,0.12,"vessel");
        tpcIFVText->SetNDC();
        tpcIFVText->SetTextFont(72);
        tpcIFVText->SetTextSize(0.03);
        tpcIFVText->SetLineWidth(4);
        tpcIFVText->Draw();

        TLatex *tpc1IText = new TLatex(0.705,0.18,"TPC inner");
        tpc1IText->SetNDC();
        tpc1IText->SetTextFont(72);
        tpc1IText->SetTextSize(0.03);
        tpc1IText->SetLineWidth(4);
        tpc1IText->Draw();

        TLatex *tpcICText = new TLatex(0.69,0.15,"containment");
        tpcICText->SetNDC();
        tpcICText->SetTextFont(72);
        tpcICText->SetTextSize(0.03);
        tpcICText->SetLineWidth(4);
        tpcICText->Draw();

        TLatex *tpcICVText = new TLatex(0.72,0.12,"vessel");
        tpcICVText->SetNDC();
        tpcICVText->SetTextFont(72);
        tpcICVText->SetTextSize(0.03);
        tpcICVText->SetLineWidth(10);
        tpcICVText->Draw();

        TLatex *tpcGasText = new TLatex(0.79,0.30,"TPC");
        tpcGasText->SetNDC();
        tpcGasText->SetTextFont(72);
        tpcGasText->SetTextSize(0.03);
        tpcGasText->SetLineWidth(10);
        tpcGasText->Draw();

        TLatex *tpcGas2Text = new TLatex(0.77,0.27,"drift gas");
        tpcGas2Text->SetNDC();
        tpcGas2Text->SetTextFont(72);
        tpcGas2Text->SetTextSize(0.03);
        tpcGas2Text->SetLineWidth(10);
        tpcGas2Text->Draw();


        TArrow *arrow = new TArrow(-150.8049,145.,-11.99843,38.629599,0.02,">"); //SSD arrow
        arrow->SetFillColor(1);
        arrow->SetFillStyle(1001);
        arrow->SetLineWidth(2.);
        arrow->Draw();

        TArrow *arrow1 = new TArrow(-160.,105.,-11.99843,25.,0.02,">"); //SDD arrow
        arrow1->SetFillColor(1);
        arrow1->SetFillStyle(1001);
        arrow1->SetLineWidth(2.);
        arrow1->Draw();


        TArrow *arrow2 = new TArrow(-150,-108.,-7.,2.,0.02,">");  //SPD arrow
        arrow2->SetFillColor(1);
        arrow2->SetFillStyle(1001);
        arrow2->SetLineWidth(2.);
        arrow2->Draw();

        TArrow *arrow3 = new TArrow(-105.,-160.,-30.,-75.,0.02,">"); //TPC field cage vessel arrow
        arrow3->SetFillColor(1);
        arrow3->SetFillStyle(1001);
        arrow3->SetLineWidth(2.);
        arrow3->Draw();


        TArrow *arrow4 = new TArrow(50.,-178.,15.,-88.,0.02,">"); //TPC rods arrow
        arrow4->SetFillColor(1);
        arrow4->SetFillStyle(1001);
        arrow4->SetLineWidth(2.);
        arrow4->Draw();

        TArrow *arrow5 = new TArrow(130.,-130.,50.,-38.,0.02,">");// TPC inner constainment vessel arrow
        arrow5->SetFillColor(1);
        arrow5->SetFillStyle(1001);
        arrow5->SetLineWidth(2.);
        arrow5->Draw();


        TArrow *arrow6 = new TArrow(160.,-90.,130.,-70.,0.02,">");// TPC gas
        arrow6->SetFillColor(1);
        arrow6->SetFillStyle(1001);
        arrow6->SetLineWidth(2.);
        arrow6->Draw();

    }

    /*****************************************************************************************************
    DrawStructure() draws the structure of the Inner Alice Detectors labeled for a xy - Plot
    ******************************************************************************************************
    ******************************************************************************************************/

    void DrawStructureNew(){
        TLatex *ssdText = new TLatex(0.16,0.87,"SSD");
        ssdText->SetTextColor(800);
        ssdText->SetNDC();
        ssdText->SetTextFont(72);
        ssdText->SetTextSize(0.03);
        ssdText->SetLineWidth(4);
        ssdText->Draw();

        TLatex *sddText = new TLatex(0.14,0.78,"SDD");
        sddText->SetTextColor(800);
        sddText->SetNDC();
        sddText->SetTextFont(72);
        sddText->SetTextSize(0.03);
        sddText->SetLineWidth(4);
        sddText->Draw();

        TLatex *spdText = new TLatex(0.14,0.27,"SPD");
        spdText->SetTextColor(800);
        spdText->SetNDC();
        spdText->SetTextFont(72);
        spdText->SetTextSize(0.03);
        spdText->SetLineWidth(4);
        spdText->Draw();

        TLatex *tpcRText = new TLatex(0.52,0.095,"TPC Rods");
        tpcRText->SetTextColor(800);
        tpcRText->SetNDC();
        tpcRText->SetTextFont(72);
        tpcRText->SetTextSize(0.03);
        tpcRText->SetLineWidth(4);
        tpcRText->Draw();

        TLatex *tpcIText = new TLatex(0.14,0.18,"TPC inner");
        tpcIText->SetTextColor(800);
        tpcIText->SetNDC();
        tpcIText->SetTextFont(72);
        tpcIText->SetTextSize(0.03);
        tpcIText->SetLineWidth(4);
        tpcIText->Draw();

        TLatex *tpcIFText = new TLatex(0.14,0.15,"field cage");
        tpcIFText->SetTextColor(800);
        tpcIFText->SetNDC();
        tpcIFText->SetTextFont(72);
        tpcIFText->SetTextSize(0.03);
        tpcIFText->SetLineWidth(4);
        tpcIFText->Draw();

        TLatex *tpcIFVText = new TLatex(0.16,0.12,"vessel");
        tpcIFVText->SetTextColor(800);
        tpcIFVText->SetNDC();
        tpcIFVText->SetTextFont(72);
        tpcIFVText->SetTextSize(0.03);
        tpcIFVText->SetLineWidth(4);
        tpcIFVText->Draw();

        TLatex *tpc1IText = new TLatex(0.705,0.18,"TPC inner");
        tpc1IText->SetTextColor(800);
        tpc1IText->SetNDC();
        tpc1IText->SetTextFont(72);
        tpc1IText->SetTextSize(0.03);
        tpc1IText->SetLineWidth(4);
        tpc1IText->Draw();

        TLatex *tpcICText = new TLatex(0.69,0.15,"containment");
        tpcICText->SetTextColor(800);
        tpcICText->SetNDC();
        tpcICText->SetTextFont(72);
        tpcICText->SetTextSize(0.03);
        tpcICText->SetLineWidth(4);
        tpcICText->Draw();

        TLatex *tpcICVText = new TLatex(0.72,0.12,"vessel");
        tpcICVText->SetTextColor(800);
        tpcICVText->SetNDC();
        tpcICVText->SetTextFont(72);
        tpcICVText->SetTextSize(0.03);
        tpcICVText->SetLineWidth(10);
        tpcICVText->Draw();

        TLatex *tpcGasText = new TLatex(0.79,0.30,"TPC");
        tpcGasText->SetTextColor(800);
        tpcGasText->SetNDC();
        tpcGasText->SetTextFont(72);
        tpcGasText->SetTextSize(0.03);
        tpcGasText->SetLineWidth(10);
        tpcGasText->Draw();

        TLatex *tpcGas2Text = new TLatex(0.77,0.27,"drift gas");
        tpcGas2Text->SetTextColor(800);
        tpcGas2Text->SetNDC();
        tpcGas2Text->SetTextFont(72);
        tpcGas2Text->SetTextSize(0.03);
        tpcGas2Text->SetLineWidth(10);
        tpcGas2Text->Draw();


        TArrow *arrow = new TArrow(-82.,90.,-11.99843,38.629599,0.02,">"); //SSD arrow
        arrow->SetLineColor(kBlack);
        arrow->SetFillColor(kBlack);
        arrow->SetFillStyle(1001);
        arrow->SetLineWidth(2.);
        arrow->Draw();

        TArrow *arrow1 = new TArrow(-90.,65.,-11.99843,25.,0.02,">"); //SDD arrow
        arrow1->SetLineColor(kBlack);
        arrow1->SetFillColor(kBlack);
        arrow1->SetFillStyle(1001);
        arrow1->SetLineWidth(2.);
        arrow1->Draw();


        TArrow *arrow2 = new TArrow(-90.,-65.,-7.,2.,0.02,">");  //SPD arrow
        arrow2->SetLineColor(kBlack);
        arrow2->SetFillColor(kBlack);
        arrow2->SetFillStyle(1001);
        arrow2->SetLineWidth(2.);
        arrow2->Draw();

        TArrow *arrow3 = new TArrow(-70.,-100.,-30.,-75.,0.02,">"); //TPC field cage vessel arrow
        arrow3->SetLineColor(kBlack);
        arrow3->SetFillColor(kBlack);
        arrow3->SetFillStyle(1001);
        arrow3->SetLineWidth(2.);
        arrow3->Draw();


        TArrow *arrow4 = new TArrow(20.,-110.,15.,-83.,0.02,">"); //TPC rods arrow
        arrow4->SetLineColor(kBlack);
        arrow4->SetFillColor(kBlack);
        arrow4->SetFillStyle(1001);
        arrow4->SetLineWidth(2.);
        arrow4->Draw();

        TArrow *arrow5 = new TArrow(80.,-85.,50.,-38.,0.02,">");// TPC inner constainment vessel arrow
        arrow5->SetLineColor(kBlack);
        arrow5->SetFillColor(kBlack);
        arrow5->SetFillStyle(1001);
        arrow5->SetLineWidth(2.);
        arrow5->Draw();


        TArrow *arrow6 = new TArrow(100.,-50.,90.,-10.,0.02,">");// TPC gas
        arrow6->SetLineColor(kBlack);
        arrow6->SetFillColor(kBlack);
        arrow6->SetFillStyle(1001);
        arrow6->SetLineWidth(2.);
        arrow6->Draw();

    }


    void DrawStructureZR(){
        TLatex *ssdText = new TLatex(0.14,0.25,"SSD");
        ssdText->SetNDC();
        ssdText->SetTextFont(72);
        ssdText->SetTextSize(0.03);
        ssdText->SetLineWidth(4);
        ssdText->Draw();

        TLatex *sddText = new TLatex(0.14,0.17,"SDD");
        sddText->SetNDC();
        sddText->SetTextFont(72);
        sddText->SetTextSize(0.03);
        sddText->SetLineWidth(4);
        sddText->Draw();

        TLatex *spdText = new TLatex(0.14,0.12,"SPD & Beam pipe");
        spdText->SetNDC();
        spdText->SetTextFont(72);
        spdText->SetTextSize(0.03);
        spdText->SetLineWidth(4);
        spdText->Draw();

        TLatex *tpcRText = new TLatex(0.14,0.53,"TPC Rods");
        tpcRText->SetNDC();
        tpcRText->SetTextFont(72);
        tpcRText->SetTextSize(0.03);
        tpcRText->SetLineWidth(4);
        tpcRText->Draw();

        TLatex *tpcIText = new TLatex(0.14,0.45,"TPC inner field");
        tpcIText->SetNDC();
        tpcIText->SetTextFont(72);
        tpcIText->SetTextSize(0.03);
        tpcIText->SetLineWidth(4);
        tpcIText->Draw();

        TLatex *tpcIFVText = new TLatex(0.14,0.42,"cage vessel");
        tpcIFVText->SetNDC();
        tpcIFVText->SetTextFont(72);
        tpcIFVText->SetTextSize(0.03);
        tpcIFVText->SetLineWidth(4);
        tpcIFVText->Draw();

        TLatex *tpc1IText = new TLatex(0.14,0.35,"TPC inner");
        tpc1IText->SetNDC();
        tpc1IText->SetTextFont(72);
        tpc1IText->SetTextSize(0.03);
        tpc1IText->SetLineWidth(4);
        tpc1IText->Draw();

        TLatex *tpcICText = new TLatex(0.14,0.32,"containment vessel");
        tpcICText->SetNDC();
        tpcICText->SetTextFont(72);
        tpcICText->SetTextSize(0.03);
        tpcICText->SetLineWidth(4);
        tpcICText->Draw();

        TLatex *tpcGasText = new TLatex(0.14,0.7,"TPC");
        tpcGasText->SetNDC();
        tpcGasText->SetTextFont(72);
        tpcGasText->SetTextSize(0.03);
        tpcGasText->SetLineWidth(10);
        tpcGasText->Draw();

        TLatex *tpcGas2Text = new TLatex(0.14,0.67,"drift gas");
        tpcGas2Text->SetNDC();
        tpcGas2Text->SetTextFont(72);
        tpcGas2Text->SetTextSize(0.03);
        tpcGas2Text->SetLineWidth(10);
        tpcGas2Text->Draw();

        TLatex *tpcMebText = new TLatex(0.70,0.54,"TPC central ");
        tpcMebText->SetNDC();
        tpcMebText->SetTextFont(72);
        tpcMebText->SetTextSize(0.03);
        tpcMebText->SetLineWidth(10);
        tpcMebText->Draw();

        TLatex *tpcMeb2Text = new TLatex(0.72,0.51,"electrode");
        tpcMeb2Text->SetNDC();
        tpcMeb2Text->SetTextFont(72);
        tpcMeb2Text->SetTextSize(0.03);
        tpcMeb2Text->SetLineWidth(10);
        tpcMeb2Text->Draw();



        TArrow *arrow = new TArrow(-90.,10.,-10,4.,0.02,">"); //SPD arrow
        arrow->SetFillColor(1);
        arrow->SetFillStyle(1001);
        arrow->SetLineWidth(2.);
        arrow->Draw();

        TArrow *arrow1 = new TArrow(-180.,22.,-27,15.,0.02,">"); //SDD1 arrow
        arrow1->SetFillColor(1);
        arrow1->SetFillStyle(1001);
        arrow1->SetLineWidth(2.);
        arrow1->Draw();


        TArrow *arrow2 = new TArrow(-180,22.,-25.,25.,0.02,">");  //SDD2 arrow
        arrow2->SetFillColor(1);
        arrow2->SetFillStyle(1001);
        arrow2->SetLineWidth(2.);
        arrow2->Draw();

        TArrow *arrow3 = new TArrow(-150.,40.,-35.,40.,0.02,">"); //SSD1 arrow
        arrow3->SetFillColor(1);
        arrow3->SetFillStyle(1001);
        arrow3->SetLineWidth(2.);
        arrow3->Draw();


        TArrow *arrow4 = new TArrow(-150.,40.,-45.,45.,0.02,">"); //SSD2 arrow
        arrow4->SetFillColor(1);
        arrow4->SetFillStyle(1001);
        arrow4->SetLineWidth(2.);
        arrow4->Draw();

        TArrow *arrow5 = new TArrow(-100.,60.,-60.,62.,0.02,">");// TPC inner containment vessel arrow
        arrow5->SetFillColor(1);
        arrow5->SetFillStyle(1001);
        arrow5->SetLineWidth(2.);
        arrow5->Draw();


        TArrow *arrow6 = new TArrow(-120.,80.,-80.,80.,0.02,">");// TPC inner field cage vessel
        arrow6->SetFillColor(1);
        arrow6->SetFillStyle(1001);
        arrow6->SetLineWidth(2.);
        arrow6->Draw();

        TArrow *arrow7 = new TArrow(-150.,100.,-67.,82.5,0.02,">");// TPC rods
        arrow7->SetFillColor(1);
        arrow7->SetFillStyle(1001);
        arrow7->SetLineWidth(2.);
        arrow7->Draw();

        TArrow *arrow8 = new TArrow(-150.,135.,-100.,140.,0.02,">");// TPC gas
        arrow8->SetFillColor(1);
        arrow8->SetFillStyle(1001);
        arrow8->SetLineWidth(2.);
        arrow8->Draw();

        TArrow *arrow9 = new TArrow(130.,100.,2.,110.,0.02,">");// TPC gas
        arrow9->SetFillColor(1);
        arrow9->SetFillStyle(1001);
        arrow9->SetLineWidth(2.);
        arrow9->Draw();


    }

    void DrawStructureZRNew(){
        TLatex *ssdText = new TLatex(0.13,0.45,"SSD");
        ssdText->SetNDC();
        ssdText->SetTextFont(72);
        ssdText->SetTextSize(0.03);
        ssdText->SetLineWidth(4);
        ssdText->Draw();

        TLatex *sddText = new TLatex(0.13,0.24,"SDD");
        sddText->SetNDC();
        sddText->SetTextFont(72);
        sddText->SetTextSize(0.03);
        sddText->SetLineWidth(4);
        sddText->Draw();

        TLatex *spdText = new TLatex(0.13,0.12,"SPD & Beam pipe");
        spdText->SetNDC();
        spdText->SetTextFont(72);
        spdText->SetTextSize(0.03);
        spdText->SetLineWidth(4);
        spdText->Draw();

        TLatex *tpcRText = new TLatex(0.13,0.85,"TPC Rods");
        tpcRText->SetNDC();
        tpcRText->SetTextFont(72);
        tpcRText->SetTextSize(0.03);
        tpcRText->SetLineWidth(4);
        tpcRText->Draw();

        TLatex *tpcIText = new TLatex(0.13,0.74,"TPC inner field");
        tpcIText->SetNDC();
        tpcIText->SetTextFont(72);
        tpcIText->SetTextSize(0.03);
        tpcIText->SetLineWidth(4);
        tpcIText->Draw();

        TLatex *tpcIFVText = new TLatex(0.13,0.71,"cage vessel");
        tpcIFVText->SetNDC();
        tpcIFVText->SetTextFont(72);
        tpcIFVText->SetTextSize(0.03);
        tpcIFVText->SetLineWidth(4);
        tpcIFVText->Draw();

        TLatex *tpc1IText = new TLatex(0.13,0.65,"TPC inner");
        tpc1IText->SetNDC();
        tpc1IText->SetTextFont(72);
        tpc1IText->SetTextSize(0.03);
        tpc1IText->SetLineWidth(4);
        tpc1IText->Draw();

        TLatex *tpcICText = new TLatex(0.13,0.62,"containment");
        tpcICText->SetNDC();
        tpcICText->SetTextFont(72);
        tpcICText->SetTextSize(0.03);
        tpcICText->SetLineWidth(4);
        tpcICText->Draw();

        TLatex *tpcICText2 = new TLatex(0.13,0.59,"vessel");
        tpcICText2->SetNDC();
        tpcICText2->SetTextFont(72);
        tpcICText2->SetTextSize(0.03);
        tpcICText2->SetLineWidth(4);
        tpcICText2->Draw();


        TLatex *tpcGasText = new TLatex(0.81,0.9,"TPC");
        tpcGasText->SetNDC();
        tpcGasText->SetTextFont(72);
        tpcGasText->SetTextSize(0.03);
        tpcGasText->SetLineWidth(10);
        tpcGasText->Draw();

        TLatex *tpcGas2Text = new TLatex(0.76,0.87,"drift gas");
        tpcGas2Text->SetNDC();
        tpcGas2Text->SetTextFont(72);
        tpcGas2Text->SetTextSize(0.03);
        tpcGas2Text->SetLineWidth(10);
        tpcGas2Text->Draw();

        TLatex *tpcMebText = new TLatex(0.73,0.74,"TPC central ");
        tpcMebText->SetNDC();
        tpcMebText->SetTextFont(72);
        tpcMebText->SetTextSize(0.03);
        tpcMebText->SetLineWidth(10);
        tpcMebText->Draw();

        TLatex *tpcMeb2Text = new TLatex(0.76,0.71,"electrode");
        tpcMeb2Text->SetNDC();
        tpcMeb2Text->SetTextFont(72);
        tpcMeb2Text->SetTextSize(0.03);
        tpcMeb2Text->SetLineWidth(10);
        tpcMeb2Text->Draw();



        TArrow *arrow = new TArrow(-55.,5.,-5,5.,0.02,">"); //SPD arrow
        arrow->SetFillColor(1);
        arrow->SetFillStyle(1001);
        arrow->SetLineWidth(2.);
        arrow->Draw();

        TArrow *arrow1 = new TArrow(-120.,19.,-15,15.,0.02,">"); //SDD1 arrow
        arrow1->SetFillColor(1);
        arrow1->SetFillStyle(1001);
        arrow1->SetLineWidth(2.);
        arrow1->Draw();


        TArrow *arrow2 = new TArrow(-120,19.,-25.,25.,0.02,">");  //SDD2 arrow
        arrow2->SetFillColor(1);
        arrow2->SetFillStyle(1001);
        arrow2->SetLineWidth(2.);
        arrow2->Draw();

        TArrow *arrow3 = new TArrow(-120.,42.,-35.,39.,0.02,">"); //SSD1 arrow
        arrow3->SetFillColor(1);
        arrow3->SetFillStyle(1001);
        arrow3->SetLineWidth(2.);
        arrow3->Draw();


        TArrow *arrow4 = new TArrow(-120.,42.,-45.,43.,0.02,">"); //SSD2 arrow
        arrow4->SetFillColor(1);
        arrow4->SetFillStyle(1001);
        arrow4->SetLineWidth(2.);
        arrow4->Draw();

        TArrow *arrow5 = new TArrow(-80.,60.,-60.,60.,0.02,">");// TPC inner containment vessel arrow
        arrow5->SetFillColor(1);
        arrow5->SetFillStyle(1001);
        arrow5->SetLineWidth(2.);
        arrow5->Draw();


        TArrow *arrow6 = new TArrow(-110.,76.,-80.,78.,0.02,">");// TPC inner field cage vessel
        arrow6->SetFillColor(1);
        arrow6->SetFillStyle(1001);
        arrow6->SetLineWidth(2.);
        arrow6->Draw();

        TArrow *arrow7 = new TArrow(-95.,86.,-65.,82.5,0.02,">");// TPC rods
        arrow7->SetFillColor(1);
        arrow7->SetFillStyle(1001);
        arrow7->SetLineWidth(2.);
        arrow7->Draw();

        TArrow *arrow8 = new TArrow(110.,91.,80.,95.,0.02,">");// TPC gas
        arrow8->SetFillColor(1);
        arrow8->SetFillStyle(1001);
        arrow8->SetLineWidth(2.);
        arrow8->Draw();

        TArrow *arrow9 = new TArrow(100.,76.,0.,95.,0.02,">");// TPC gas
        arrow9->SetFillColor(1);
        arrow9->SetFillStyle(1001);
        arrow9->SetLineWidth(2.);
        arrow9->Draw();


    }



    /*************************************************************************************************************
    DrawArmenteros() draws the labels for the particles in the Armenteros plot
    **************************************************************************************************************/

    void DrawArmenteros(){
        TLatex *k0s = new TLatex(0.48,0.73,"K^{0}_{s}");
        k0s->SetNDC();
        k0s->SetTextFont(62);
        k0s->SetTextSize(0.05);
        k0s->SetLineWidth(4);
        k0s->Draw();

        TLatex *lambda = new TLatex(0.72,0.47,"#Lambda");
        lambda->SetNDC();
        lambda->SetTextFont(62);
        lambda->SetTextSize(0.05);
        lambda->SetLineWidth(4);
        lambda->Draw();

        TLatex *lambdabar = new TLatex(0.25,0.48,"#bar{#Lambda}");
        lambdabar->SetNDC();
        lambdabar->SetTextFont(62);
        lambdabar->SetTextSize(0.05);
        lambdabar->SetLineWidth(4);
        lambdabar->Draw();

        TLatex *gamma = new TLatex(0.5,0.18,"#gamma");
        gamma->SetNDC();
        gamma->SetTextFont(62);
        gamma->SetTextSize(0.05);
        gamma->SetLineWidth(4);
        gamma->Draw();

    }


    /************************************************************************************************
    DrawdEdxLabels()
    * Float_t linewidth - will give the final linewidth in the plot

    *************************************************************************************************
    *************************************************************************************************/
    void DrawdEdxLabel(){
        TLatex *muon = new TLatex(0.16,0.6,"#mu"); //text at muon line
        muon->SetNDC();
        muon->SetTextColor(1);
        muon->SetTextFont(62);
        muon->SetTextSize(0.04);
        muon->SetLineWidth(2);

        TLatex *pion = new TLatex(0.4,0.17,"#pi"); //text at pion line
        pion->SetNDC();
        pion->SetTextColor(1);
        pion->SetTextFont(62);
        pion->SetTextSize(0.04);
        pion->SetLineWidth(2);

        TLatex *kaon = new TLatex(0.35,0.8,"K"); //text at kaon line
        kaon->SetNDC();
        kaon->SetTextColor(1);
        kaon->SetTextFont(62);
        kaon->SetTextSize(0.04);
        kaon->SetLineWidth(2);

        TLatex *proton = new TLatex(0.44,0.7,"p"); //text at proton line
        proton->SetNDC();
        proton->SetTextColor(1);
        proton->SetTextFont(62);
        proton->SetTextSize(0.04);
        proton->SetLineWidth(2);

        TLatex *electron = new TLatex(0.58,0.4,"e"); //text at electron line
        electron->SetNDC();
        electron->SetTextColor(1);
        electron->SetTextFont(62);
        electron->SetTextSize(0.04);
        electron->SetLineWidth(2);

        electron->Draw("same");
        muon->Draw("same");
        kaon->Draw("same");
        proton->Draw("same");
        pion->Draw("same");
    }

    /************************************************************************************************
    DrawTOFLabels()
    * Float_t linewidth - will give the final linewidth in the plot

    *************************************************************************************************
    *************************************************************************************************/
    void DrawTOFLabels(){
        TLatex *muon = new TLatex(0.19,0.3,"#mu"); //text at muon line
        muon->SetNDC();
        muon->SetTextColor(1);
        muon->SetTextFont(62);
        muon->SetTextSize(0.04);
        muon->SetLineWidth(2);

        TLatex *pion = new TLatex(0.22,0.3,"#pi"); //text at pion line
        pion->SetNDC();
        pion->SetTextColor(1);
        pion->SetTextFont(62);
        pion->SetTextSize(0.04);
        pion->SetLineWidth(2);

        TLatex *kaon = new TLatex(0.25,0.8,"K"); //text at kaon line
        kaon->SetNDC();
        kaon->SetTextColor(1);
        kaon->SetTextFont(62);
        kaon->SetTextSize(0.04);
        kaon->SetLineWidth(2);

        TLatex *proton = new TLatex(0.37,0.7,"p"); //text at proton line
        proton->SetNDC();
        proton->SetTextColor(1);
        proton->SetTextFont(62);
        proton->SetTextSize(0.04);
        proton->SetLineWidth(2);

        TLatex *electron = new TLatex(0.2,0.18,"e"); //text at electron line
        electron->SetNDC();
        electron->SetTextColor(1);
        electron->SetTextFont(62);
        electron->SetTextSize(0.04);
        electron->SetLineWidth(2);

        electron->Draw("same");
        muon->Draw("same");
        kaon->Draw("same");
        proton->Draw("same");
        pion->Draw("same");
    }

    /* // DrawGammaLines will draw the lines in the histogram for you
    * startX - starting point of drawing in x
    * endX - end point of drawing in x
    * startY -starting point of drawing in y
    * endY - end point of drawing in y
    * linew - line width
    */
    void DrawGammaLines(Float_t startX, Float_t endX,
                    Float_t startY, Float_t endY,
                    Float_t linew, Float_t lineColor = 4, Style_t lineStyle = 1){
        TLine * l1 = new TLine (startX,startY,endX,endY);
        l1->SetLineColor(lineColor);
        l1->SetLineWidth(linew);
        l1->SetLineStyle(lineStyle);
        l1->Draw("same");
    }


    TBox* CreateBoxFromGraph( TGraphAsymmErrors* graph, Double_t xStart, Double_t yStart, Double_t xEnd, Double_t yEnd ) {
        TBox* box = new TBox(xStart ,yStart , xEnd, yEnd);
        box->SetLineColor(graph->GetMarkerColor());
        box->SetLineWidth(graph->GetLineWidth());
        box->SetFillStyle(0);
        return box;
    }

    TBox* CreateBoxFromGraphWithFill( TGraphAsymmErrors* graph, Double_t xStart, Double_t yStart, Double_t xEnd, Double_t yEnd ) {
        TBox* box = new TBox(xStart ,yStart , xEnd, yEnd);
        box->SetLineColor(graph->GetMarkerColor());
        box->SetLineWidth(graph->GetLineWidth());
        box->SetFillStyle(graph->GetFillStyle());
        box->SetFillColor(graph->GetFillColor());
        return box;
    }


    TBox* CreateBoxConv(Color_t colorBox, Double_t xStart, Double_t yStart, Double_t xEnd, Double_t yEnd, Style_t fillStyle = 1001 ) {
        TBox* box = new TBox(xStart ,yStart , xEnd, yEnd);
        box->SetLineColor(colorBox);
        box->SetFillColor(colorBox);
        box->SetFillStyle(fillStyle);
        return box;
    }


    TBox* CreateBoxFromGraph( TGraphErrors* graph, Double_t xStart, Double_t yStart, Double_t xEnd, Double_t yEnd ) {
        TBox* box =  new TBox(xStart ,yStart , xEnd, yEnd);
        box->SetLineColor(graph->GetMarkerColor());
        box->SetLineWidth(graph->GetLineWidth());;
        box->SetFillStyle(0);
        return box;
    }


    TBox* CreateBoxFromGraph( TGraph* graph, Double_t xStart, Double_t yStart, Double_t xEnd, Double_t yEnd ) {
        TBox* box =  new TBox(xStart ,yStart , xEnd, yEnd);
        box->SetLineColor(graph->GetMarkerColor());
        box->SetLineWidth(graph->GetLineWidth());
        box->SetFillStyle(0);
        return box;
    }

    TBox* CreateBoxFromHisto( TH1* histo, Double_t xStart, Double_t yStart, Double_t xEnd, Double_t yEnd ) {
        TBox* box = new TBox(xStart ,yStart , xEnd, yEnd);
        box->SetLineColor(histo->GetMarkerColor());
        box->SetLineWidth(histo->GetLineWidth());
        box->SetFillStyle(0);
        return box;
    }

    TMarker* CreateMarkerFromGraph( TGraphAsymmErrors* graph, Double_t x, Double_t y, Double_t scaleSize) {
        TMarker* marker= new TMarker(x,y, graph->GetMarkerStyle());
        marker->SetMarkerColor(graph->GetMarkerColor());
        marker->SetMarkerSize(graph->GetMarkerSize() * scaleSize);
        return marker;
    }

    TMarker* CreateMarkerFromGraph( TGraphErrors* graph, Double_t x, Double_t y, Double_t scaleSize) {
        TMarker* marker= new TMarker(x,y, graph->GetMarkerStyle());
        marker->SetMarkerColor(graph->GetMarkerColor());
        marker->SetMarkerSize(graph->GetMarkerSize() * scaleSize);
        return marker;
    }

    TMarker* CreateMarkerFromGraph( TGraph* graph, Double_t x, Double_t y, Double_t scaleSize) {
        TMarker* marker= new TMarker(x,y, graph->GetMarkerStyle());
        marker->SetMarkerColor(graph->GetMarkerColor());
        marker->SetMarkerSize(graph->GetMarkerSize() * scaleSize);
        return marker;
    }

    TMarker* CreateMarkerFromHisto( TH1* histo, Double_t x, Double_t y, Double_t scaleSize) {
        TMarker* marker= new TMarker(x,y, histo->GetMarkerStyle());
        marker->SetMarkerColor(histo->GetMarkerColor());
        marker->SetMarkerSize(histo->GetMarkerSize() * scaleSize);
        return marker;
    }

    TLine* CreateLineFromFit(TF1* fit, Double_t xStart, Double_t yStart, Double_t xEnd, Double_t yEnd, Double_t scaleWidth) {
        TLine *line = new TLine(xStart,yStart,xEnd,yEnd);
        line->SetLineColor(fit->GetLineColor());
        line->SetLineWidth(fit->GetLineWidth()*scaleWidth);
        line->SetLineStyle(fit->GetLineStyle());
        return line;
    }

    TLine* CreateLineFromHisto(TH1* histo, Double_t xStart, Double_t yStart, Double_t xEnd, Double_t yEnd, Double_t scaleWidth) {
        TLine * line = new TLine (xStart,yStart,xEnd,yEnd);
        line->SetLineColor(histo->GetLineColor());
        line->SetLineWidth(histo->GetLineWidth()*scaleWidth);
        line->SetLineStyle(histo->GetLineStyle());
        return line;
    }

    TLine* CreateLineFromGraph(TGraph* graph, Double_t xStart, Double_t yStart, Double_t xEnd, Double_t yEnd, Double_t scaleWidth) {
        TLine * line = new TLine (xStart,yStart,xEnd,yEnd);
        line->SetLineColor(graph->GetLineColor());
        line->SetLineWidth(graph->GetLineWidth()*scaleWidth);
        line->SetLineStyle(graph->GetLineStyle());
        return line;
    }

    void drawLatexAdd(TString latextext, Double_t textcolumn, Double_t textrow, Double_t textSizePixel,Bool_t setFont = kFALSE, Bool_t setFont2 = kFALSE, Bool_t alignRight = kFALSE, Color_t textcolor = kBlack){
    TLatex *latexDummy                  = new TLatex(textcolumn ,textrow,latextext);
    SetStyleTLatex( latexDummy, textSizePixel,4);
    if(setFont)
        latexDummy->SetTextFont(62);
    if(setFont2)
        latexDummy->SetTextFont(43);
    if(alignRight)
        latexDummy->SetTextAlign(31);
    latexDummy->SetTextColor(textcolor);
    latexDummy->Draw();
}

    void DrawAliceLogoSimple(Float_t startX, Float_t startY, Float_t widthLogo, Double_t xLengthCanvas, Double_t yLengthCanvas){


        Double_t widthLogoPix = xLengthCanvas*widthLogo;
        Double_t heightLogoPix = widthLogoPix/0.73447;
        Double_t totalXLogo = (startX*xLengthCanvas + widthLogoPix)/xLengthCanvas;
        Double_t totalYLogo = (startY*yLengthCanvas + heightLogoPix)/yLengthCanvas;


    TPad *myPadLogo = new TPad("myPadLogo", "Pad for ALICE Logo",startX,startY, totalXLogo,totalYLogo);
    myPadLogo->SetFillColor(0);
    myPadLogo->SetBorderMode(0);
    myPadLogo->SetBorderSize(2);
    myPadLogo->SetFrameBorderMode(0);
    myPadLogo->SetLeftMargin(0.0);
    myPadLogo->SetTopMargin(0.0);
    myPadLogo->SetBottomMargin(0.0);
    myPadLogo->SetRightMargin(0.0);
    TImage *myAliceLogo = TImage::Open("ALICE_logo.eps");
    myPadLogo->Draw();  // to take out for not using a logo.
    myPadLogo->cd();

    myAliceLogo->Draw();

    }

    void DrawAliceLogoSimplePreliminary(Float_t startX, Float_t startY, Float_t widthLogo, Size_t textSize, Double_t xLengthCanvas, Double_t yLengthCanvas){


        Double_t widthLogoPix = xLengthCanvas*widthLogo;
        Double_t heightLogoPix = widthLogoPix/0.73447;
        Double_t totalXLogo = (startX*xLengthCanvas + widthLogoPix)/xLengthCanvas;
        Double_t totalYLogo = (startY*yLengthCanvas + heightLogoPix)/yLengthCanvas;

        Float_t aliceStartY = startY - textSize- 0.03;
        TLatex *alice = new TLatex(startX,aliceStartY,"PRELIMINARY"); // Bo: this was modified
        alice->SetNDC();
        alice->SetTextColor(1);
        alice->SetTextFont(62);
        alice->SetTextSize(textSize);
        alice->SetLineWidth(2);
        alice->Draw("same");


    TPad *myPadLogo = new TPad("myPadLogo", "Pad for ALICE Logo",startX,startY, totalXLogo,totalYLogo);
    myPadLogo->SetFillColor(0);
    myPadLogo->SetBorderMode(0);
    myPadLogo->SetBorderSize(2);
    myPadLogo->SetFrameBorderMode(0);
    myPadLogo->SetLeftMargin(0.0);
    myPadLogo->SetTopMargin(0.0);
    myPadLogo->SetBottomMargin(0.0);
    myPadLogo->SetRightMargin(0.0);
    TImage *myAliceLogo = TImage::Open("ALICE_logo.eps");
    myPadLogo->Draw();  // to take out for not using a logo.
    myPadLogo->cd();

    myAliceLogo->Draw();

    }


    void PrintLevyFitResults(Float_t startX,
                            Float_t startY,
                            Float_t spacingxVal,
                            Float_t spacingxText,
                            Float_t spacingy,
                            Size_t textsize,
                            Double_t* fitResults){
        cout << fitResults[0] << "\t" << fitResults[1] << "\t" << fitResults[2] << endl;
        cout << fitResults[3] << "\t" << fitResults[4] << "\t" << fitResults[5] << endl;
        cout << fitResults[6] << "\t" << fitResults[7] << "\t" << fitResults[8] << endl;

        TLatex *latexdNdyLabel = new TLatex(startX,startY,"dN/dy:");
        TLatex *latexdNdyValue = new TLatex(startX + spacingxText,startY, Form("%1.2f",fitResults[0]));
        TLatex *latexdNdyStatErr = new TLatex(startX + spacingxText + spacingxVal ,startY, Form("#pm %1.2f_{stat}",fitResults[1]));
        TLatex *latexdNdySysErr = new TLatex(startX + spacingxText + 2*spacingxVal ,startY, Form("#pm %1.2f_{sys}",fitResults[2]));

        TLatex *latexNLabel = new TLatex(startX,startY - spacingy,"n:");
        TLatex *latexNValue = new TLatex(startX + spacingxText,startY - spacingy, Form("%1.2f",fitResults[3]));
        TLatex *latexNStatErr = new TLatex(startX + spacingxText + spacingxVal ,startY -spacingy, Form("#pm %1.2f_{stat}",fitResults[4]));
        TLatex *latexNSysErr = new TLatex(startX + spacingxText + 2*spacingxVal ,startY- spacingy , Form("#pm %1.2f_{sys}",fitResults[5]));

        TLatex *latexTLabel = new TLatex(startX,startY - 2*spacingy,"T_{Levy}:");
        TLatex *latexTValue = new TLatex(startX + spacingxText,startY -2*spacingy , Form("%1.3f",fitResults[6]));
        TLatex *latexTStatErr = new TLatex(startX + spacingxText + spacingxVal ,startY-2*spacingy , Form("#pm %1.3f_{stat}",fitResults[7]));
        TLatex *latexTSysErr = new TLatex(startX + spacingxText + 2*spacingxVal ,startY-2*spacingy , Form("#pm %1.3f_{sys}",fitResults[8]));


        latexdNdyLabel->SetNDC();
        latexdNdyLabel->SetTextColor(1);
        latexdNdyLabel->SetTextSize(textsize);
        latexdNdyLabel->Draw();

        latexdNdyValue->SetNDC();
        latexdNdyValue->SetTextColor(1);
        latexdNdyValue->SetTextSize(textsize);
        latexdNdyValue->Draw();

        latexdNdyStatErr->SetNDC();
        latexdNdyStatErr->SetTextColor(1);
        latexdNdyStatErr->SetTextSize(textsize);
        latexdNdyStatErr->Draw();

        latexdNdySysErr->SetNDC();
        latexdNdySysErr->SetTextColor(1);
        latexdNdySysErr->SetTextSize(textsize);
        latexdNdySysErr->Draw();

        latexNLabel->SetNDC();
        latexNLabel->SetTextColor(1);
        latexNLabel->SetTextSize(textsize);
        latexNLabel->Draw();

        latexNValue->SetNDC();
        latexNValue->SetTextColor(1);
        latexNValue->SetTextSize(textsize);
        latexNValue->Draw();

        latexNStatErr->SetNDC();
        latexNStatErr->SetTextColor(1);
        latexNStatErr->SetTextSize(textsize);
        latexNStatErr->Draw();

        latexNSysErr->SetNDC();
        latexNSysErr->SetTextColor(1);
        latexNSysErr->SetTextSize(textsize);
        latexNSysErr->Draw();

        latexTLabel->SetNDC();
        latexTLabel->SetTextColor(1);
        latexTLabel->SetTextSize(textsize);
        latexTLabel->Draw();

        latexTValue->SetNDC();
        latexTValue->SetTextColor(1);
        latexTValue->SetTextSize(textsize);
        latexTValue->Draw();

        latexTStatErr->SetNDC();
        latexTStatErr->SetTextColor(1);
        latexTStatErr->SetTextSize(textsize);
        latexTStatErr->Draw();

        latexTSysErr->SetNDC();
        latexTSysErr->SetTextColor(1);
        latexTSysErr->SetTextSize(textsize);
        latexTSysErr->Draw();

    }

    void DrawBinShiftingText( Float_t startX,
                            Float_t startY ,
                            Size_t textsize,
                            Double_t differenceText=0.021){  // kk

        TLatex *latexBinShifting1= new TLatex(startX,startY,"Correction in E #frac{d^{3}#sigma}{dp^{3}}");
            latexBinShifting1->SetNDC();
        latexBinShifting1->SetTextColor(1);
        latexBinShifting1->SetTextSize(textsize);
        //latexBinShifting1->SetTextAngle(90);
        latexBinShifting1->Draw();

        TLatex *latexBinShifting2= new TLatex(startX,startY-(differenceText+textsize),"for bin width applied");
            latexBinShifting2->SetNDC();
        latexBinShifting2->SetTextColor(1);
        latexBinShifting2->SetTextSize(textsize);
        //latexBinShifting2->SetTextAngle(90);
        latexBinShifting2->Draw();

        TLatex *latexBinShifting3= new TLatex(startX,startY-2*differenceText-1.5*textsize,"NIM A 355 (1995) 541");
            latexBinShifting3->SetNDC();
        latexBinShifting3->SetTextColor(1);
        latexBinShifting3->SetTextSize(textsize);
        //latexBinShifting3->SetTextAngle(90);
        latexBinShifting3->Draw();
    }


    void DrawPrevPrelimPi0(Float_t startX,
                        Float_t startY ,
                        Size_t textsize){  // kk
        TLatex *latex7000prev1 = new TLatex(startX,startY, "Previous prelimimary #pi^{0} @ 7 TeV:  #sigma_{MB} = 67mb #pm 10% used");
        //	TLatex *latex7000Errorprev1 = new TLatex(startX + spacingx ,startY-1.25*spacingy,"#pm 4.3 mb");

        latex7000prev1->SetNDC();
        latex7000prev1->SetTextColor(1);
        latex7000prev1->SetTextSize(textsize);
        latex7000prev1->Draw();
    }

    void DrawNormalizationErrorText( Float_t startX,
                                    Float_t startY,
                                    Float_t spacingx,
                                    Float_t spacingy,
                                    Size_t textsize,
                                    TString energy){

        TLatex *latexNormErrorText1 = new TLatex(startX,startY,"#sigma_{pp} uncertainty");
        latexNormErrorText1->SetNDC();
        latexNormErrorText1->SetTextColor(1);
        latexNormErrorText1->SetTextSize(textsize);
        latexNormErrorText1->Draw();

        if(energy == "7TeV"){
            TLatex *latex7000 = new TLatex(startX,startY -1.25*spacingy, "7 TeV");
            TLatex *latex7000Error = new TLatex(startX + spacingx ,startY-1.25*spacingy,"#pm 3.5%");  //4.3 mb
            latex7000->SetNDC();
            latex7000->SetTextColor(1);
            latex7000->SetTextSize(textsize);
            latex7000->Draw();

            latex7000Error->SetNDC();
            latex7000Error->SetTextColor(1);
            latex7000Error->SetTextSize(textsize);
            latex7000Error->Draw();
        }


        if(energy == "2.76TeV"){
            TLatex *latex2760 = new TLatex(startX,startY -1.25*spacingy, "2.76 TeV");
            TLatex *latex2760Error = new TLatex(startX + spacingx ,startY-1.25*spacingy,"#pm  7%"); //3.9bm
            latex2760->SetNDC();
            latex2760->SetTextColor(1);
            latex2760->SetTextSize(textsize);
            latex2760->Draw();

            latex2760Error->SetNDC();
            latex2760Error->SetTextColor(1);
            latex2760Error->SetTextSize(textsize);
            latex2760Error->Draw();
        }

        if(energy == "900GeV"){
            TLatex *latex900 = new TLatex(startX,startY -1.25*spacingy, "0.9 TeV");
            TLatex *latex900Error = new TLatex(startX + spacingx,startY-1.25*spacingy,"+ 5.0%, - 3.9%");

            latex900->SetNDC();
            latex900->SetTextColor(1);
            latex900->SetTextSize(textsize);
            latex900->Draw();

            latex900Error->SetNDC();
            latex900Error->SetTextColor(1);
            latex900Error->SetTextSize(textsize);
            latex900Error->Draw();
        }


        if(energy == "all"){
            TLatex *latex7000 = new TLatex(startX,startY -1.25*spacingy, "7 TeV");
            TLatex *latex7000Error = new TLatex(startX + spacingx ,startY-1.25*spacingy,"#pm 3.5%");

            TLatex *latex2760 = new TLatex(startX,startY -2.1*spacingy, "2.76 TeV");
            TLatex *latex2760Error = new TLatex(startX + spacingx ,startY-2.1*spacingy,"#pm  7%");

            TLatex *latex900 = new TLatex(startX,startY -2.95*spacingy, "0.9 TeV");
            TLatex *latex900Error = new TLatex(startX + spacingx ,startY-2.95*spacingy,"+ 5.0%, - 3.9%");

            latex7000->SetNDC();
            latex7000->SetTextColor(1);
            latex7000->SetTextSize(textsize);
            latex7000->Draw();

            latex7000Error->SetNDC();
            latex7000Error->SetTextColor(1);
            latex7000Error->SetTextSize(textsize);
            latex7000Error->Draw();

            latex2760->SetNDC();
            latex2760->SetTextColor(1);
            latex2760->SetTextSize(textsize);
            latex2760->Draw();

            latex2760Error->SetNDC();
            latex2760Error->SetTextColor(1);
            latex2760Error->SetTextSize(textsize);
            latex2760Error->Draw();

            latex900->SetNDC();
            latex900->SetTextColor(1);
            latex900->SetTextSize(textsize);
            latex900->Draw();

            latex900Error->SetNDC();
            latex900Error->SetTextColor(1);
            latex900Error->SetTextSize(textsize);
            latex900Error->Draw();

        }
        if(energy == "No2.76"){
            TLatex *latex7000 = new TLatex(startX,startY -1.25*spacingy, "7 TeV");
            TLatex *latex7000Error = new TLatex(startX + spacingx ,startY-1.25*spacingy,"#pm 3.5%");

            TLatex *latex900 = new TLatex(startX,startY -2.1*spacingy, "0.9 TeV");
            TLatex *latex900Error = new TLatex(startX + spacingx ,startY-2.1*spacingy,"+ 5.0%, - 3.9%");

            latex7000->SetNDC();
            latex7000->SetTextColor(1);
            latex7000->SetTextSize(textsize);
            latex7000->Draw();

            latex7000Error->SetNDC();
            latex7000Error->SetTextColor(1);
            latex7000Error->SetTextSize(textsize);
            latex7000Error->Draw();

            latex900->SetNDC();
            latex900->SetTextColor(1);
            latex900->SetTextSize(textsize);
            latex900->Draw();

            latex900Error->SetNDC();
            latex900Error->SetTextColor(1);
            latex900Error->SetTextSize(textsize);
            latex900Error->Draw();
        }
    }

    void DrawIndividualTextSlicesR ( Double_t* yValues,
                                    Size_t textsize,
                                    TString option ){

    // 	Double_t rValuesLin[9] = {4.,10.,19.,28.,40,48.,64.,90.,110.};
    // 	Double_t rValuesLog[9] = {4.,8.,19.,28.,40,48.,64.,83.,110.};
        Double_t rValuesLin[9] = {4.,10.,19.,28.,40,48.,64.,90.,110.};
        Double_t rValuesLog[9] = {4.,8.,19.,28.,40,48.,64.,83.,110.};

        Double_t rValues[9];
        if (option.CompareTo("lin") == 0){
            for ( Int_t i = 0; i<9 ; i++){
                rValues[i] = rValuesLin[i];
            }
        } else {
            for ( Int_t i = 0; i<9 ; i++){
                rValues[i] = rValuesLog[i];
            }
        }

        TLatex *latex1 = NULL;
        TLatex *latex0 = NULL;
        if (option.CompareTo("lin") == 0){
            latex0 = new TLatex(rValues[0],yValues[0], "Beam Pipe");
            latex1 = new TLatex(rValues[1],yValues[1], "SPD");
        } else {
            latex1 = new TLatex(rValues[1],yValues[1], "Beam Pipe & SPD");
        }
        TLatex *latex2 = new TLatex(rValues[2],yValues[2], "SDD 1^{st} Layer");
        TLatex *latex3 = new TLatex(rValues[3],yValues[3], "SDD 2^{nd} Layer + Support Structures");
        TLatex *latex4 = new TLatex(rValues[4],yValues[4], "SSD 1^{st} Layer");
        TLatex *latex5 = new TLatex(rValues[5],yValues[5], "SSD 2^{nd} Layer");
        TLatex *latex6 = new TLatex(rValues[6],yValues[6], "TPC Inner Containment Vessel");
        TLatex *latex7 = new TLatex(rValues[7],yValues[7], "TPC Inner Field Cage Vessel");
        TLatex *latex8 = new TLatex(rValues[8],yValues[8], "TPC Gas");
        if (option.CompareTo("lin") == 0){
            latex0->SetTextColor(1);
            latex0->SetTextSize(textsize);
            latex0->SetTextAngle(90);
            latex0->Draw();
        }

        latex1->SetTextColor(1);
        latex1->SetTextSize(textsize);
        latex1->SetTextAngle(90);
        latex1->Draw();

        latex2->SetTextColor(1);
        latex2->SetTextSize(textsize);
        latex2->SetTextAngle(90);
        latex2->Draw();

        latex3->SetTextColor(1);
        latex3->SetTextSize(textsize);
        latex3->SetTextAngle(90);
        latex3->Draw();

        latex4->SetTextColor(1);
        latex4->SetTextSize(textsize);
        latex4->SetTextAngle(90);
        latex4->Draw();

        latex5->SetTextColor(1);
        latex5->SetTextSize(textsize);
        latex5->SetTextAngle(90);
        latex5->Draw();

        latex6->SetTextColor(1);
        latex6->SetTextSize(textsize);
        latex6->SetTextAngle(90);
        latex6->Draw();

        latex7->SetTextColor(1);
        latex7->SetTextSize(textsize);
        latex7->SetTextAngle(90);
        latex7->Draw();

        latex8->SetTextColor(1);
        latex8->SetTextSize(textsize);
        latex8->SetTextAngle(90);
        latex8->Draw();

    }


    void DrawCentrality(Double_t x_low,
                        Double_t y_low,
                        TString centrality){
        TText *centr = new TText(x_low,y_low,centrality);
        centr->SetNDC();
        centr->SetTextFont(42);
        centr->SetTextSize(0.04);
        centr->Draw();
    }

    void DrawSystem( Double_t x_low,
                    Double_t y_low,
                    TString energy,
                    TString centrality){

        TString collisionenergy ="";
        if(energy.CompareTo("7TeV") == 0) collisionenergy = "pp, #sqrt{#it{s}} = 7 TeV";
        if(energy.CompareTo("HI") == 0) collisionenergy = centrality+"  Pb-Pb, #sqrt{#it{s}_{_{NN}}} = 2.76 TeV";
        if(energy.CompareTo("900GeV") == 0) collisionenergy = "pp, #sqrt{#it{s}} = 900 GeV";
        if(energy.CompareTo("2760GeV") == 0) collisionenergy = "pp, #sqrt{#it{s}} = 2.76 TeV";
        TLatex *CollisionEnergy = new TLatex(x_low,y_low,collisionenergy);
        CollisionEnergy->SetNDC();
        CollisionEnergy->SetTextColor(kBlack);
        CollisionEnergy->SetTextFont(42);
        CollisionEnergy->SetTextSize(0.05);
        CollisionEnergy->SetLineWidth(2);
        CollisionEnergy->Draw();

        /* TText *centr = new TText(x_low,y_low,centrality); */
        /*   centr->SetNDC(); */
        /*   centr->SetTextFont(42); */
        /*   centr->SetTextSize(0.05); */
        /*   centr->Draw(); */

    }




    void PutProcessLabelAndEnergyOnPlot(    Double_t startTextX,
                                            Double_t startTextY,
                                            Size_t textHeight,
                                            TString fEnergy,
                                            TString fDecayChannel,
                                            TString fDetectionChannel,
                                            Style_t textFont            = 62,
                                            Size_t textHeightRel        = 0.03,
                                            TString fPeriodName         = "",
                                            Color_t textColor           = 1,
                                            Float_t textHeightFac       = 1.25,
                                            Int_t textAlign             = 11
                                    ){

        Double_t differenceText     = textHeight*textHeightFac;
        if (textFont == 63 || textFont == 43) differenceText = textHeightRel*textHeightFac;

        Double_t startYPeriod       = startTextY-2*differenceText;
        Double_t startYProcess      = startTextY-2*differenceText;
        Double_t startYdDetProcess  = startTextY-3*differenceText;
        if (fPeriodName.CompareTo("") != 0 && fPeriodName.CompareTo("No") != 0 ){
            startYProcess         = startTextY-3*differenceText;
            startYdDetProcess     = startTextY-4*differenceText;
        }


        TLatex *energy          = new TLatex(startTextX, (startTextY-differenceText), fEnergy);
        TLatex *process         = new TLatex(startTextX, startYProcess, fDecayChannel);
        TLatex *detprocess      = new TLatex(startTextX, startYdDetProcess, fDetectionChannel);

        energy->SetNDC();
        energy->SetTextColor(textColor);
        energy->SetTextFont(textFont);
        energy->SetTextSize(textHeight);
        energy->SetTextAlign(textAlign);
        energy->Draw();

        process->SetNDC();
        process->SetTextColor(textColor);
        process->SetTextFont(textFont);
        process->SetTextSize(textHeight);
        process->SetTextAlign(textAlign);
        process->Draw();

        detprocess->SetNDC();
        detprocess->SetTextColor(textColor);
        detprocess->SetTextFont(textFont);
        detprocess->SetTextSize(textHeight);
        detprocess->SetTextAlign(textAlign);
        detprocess->Draw();

        if (fPeriodName.CompareTo("") != 0 && fPeriodName.CompareTo("No") != 0 ){
            TLatex *period          = new TLatex(startTextX, startYPeriod, fPeriodName);
            period->SetNDC();
            period->SetTextColor(textColor);
            period->SetTextFont(textFont);
            period->SetTextSize(textHeight);
            period->SetTextAlign(textAlign);
            period->Draw();
        }
    }

    void PutThisThesisLabel(Double_t    startTextX,
                            Double_t    startTextY,
                            Size_t      textHeight,
                            Size_t      textHeightRel   = 0.03,
                            Float_t     textHeightFac   = 1.25,
                            Style_t     textFont        = 62,
                            Color_t     textColor       = 1,
                            Bool_t      alignLeft       = kTRUE) {

        Double_t differenceText = textHeight*textHeightFac;
        if (textFont == 63 || textFont == 43)
            differenceText        = textHeightRel*textHeightFac;

        TString labelString     = "ALICE this thesis";
        TLatex *label           = new TLatex(startTextX, (startTextY-differenceText), labelString);

        if(!alignLeft) label->SetTextAlign(31);
        label->SetNDC();
        label->SetTextColor(textColor);
        label->SetTextFont(textFont);
        label->SetTextSize(textHeight);
        label->Draw();
    }

    void PutALICESimulationLabel(Double_t   startTextX,
                                Double_t   startTextY,
                                Size_t     textHeight,
                                Size_t     textHeightRel   = 0.03,
                                Float_t    textHeightFac   = 1.25,
                                Style_t    textFont        = 62,
                                Color_t    textColor       = 1,
                                Bool_t     alignLeft       = kTRUE) {

        Double_t differenceText = textHeight*textHeightFac;
        if (textFont == 63 || textFont == 43)
            differenceText        = textHeightRel*textHeightFac;

        TString labelString     = "ALICE simulation";
        TLatex *label           = new TLatex(startTextX, (startTextY-differenceText), labelString);

        if(!alignLeft) label->SetTextAlign(31);
        label->SetNDC();
        label->SetTextColor(textColor);
        label->SetTextFont(textFont);
        label->SetTextSize(textHeight);
        label->Draw();
    }

    Color_t GetColorDefaultColor( TString energy,
                                TString generator,
                                TString centrality,
                                Bool_t kBox = kFALSE,
                                Bool_t output = kTRUE){

        if(output) cout << energy.Data() << "\t" << generator.Data() << "\t" << centrality.Data() << "\t" << kBox << endl;
        if (!energy.CompareTo("900GeV")){
            if (!kBox){
                if (!generator.CompareTo("")||!generator.CompareTo("LHC10c_900GeV_pass4")){
                    return kRed+2;
                } else if (!generator.CompareTo("Pythia")||!generator.CompareTo("LHC14j4c_900GeV")){
                    return kRed+4;
                } else if (!generator.CompareTo("Pythia2")){
                    return kRed+3;
                } else if (!generator.CompareTo("PythiaAddSig")){
                    return kRed-4;
                } else if (!generator.CompareTo("Phojet")){
                    return kRed-2;
                }
            } else {
                return kRed -5;
            }
        }

        if (!energy.CompareTo("2.76TeV")){
            if (!kBox){
                if (!generator.CompareTo("")||!generator.CompareTo("LHC13g")||!generator.CompareTo("LHC13g-kEMC7")||
                        !generator.CompareTo("LHC13g-kEMCEG1")||!generator.CompareTo("LHC13g-kEMCEG2")||
                        !generator.CompareTo("LHC11a_p4_wSDD")||!generator.CompareTo("LHC11a_p4_wSDD-kEMC1") ||
                        !generator.CompareTo("LHC11a_pass4_wSDD")||!generator.CompareTo("LHC11a_pass4_wSDD-kEMC1")
                ){
                    return kMagenta+2;
                } else if (!generator.CompareTo("Pythia")||!generator.CompareTo("LHC12f1a")||!generator.CompareTo("LHC15g2")){
                    return kMagenta+4;
                } else if (!generator.CompareTo("Pythia2") || !generator.CompareTo("LHC12i3")){
                    return kMagenta+3;
                } else if (!generator.CompareTo("PythiaAddSig")||!generator.CompareTo("LHC12i3_AddSig")){
                    return kMagenta-4;
                } else if (!generator.CompareTo("Phojet")||!generator.CompareTo("LHC12f1b")){
                    return kMagenta-2;
                } else if (!generator.CompareTo("Pythia8-JJ")||!generator.CompareTo("JJ") || !generator.CompareTo("LHC15g1a") || !generator.CompareTo("LHC15a3a")){
                    return 807;
                }
            } else {
                return kMagenta-5;
            }
        }

        if (!energy.CompareTo("5TeV") || !energy.CompareTo("5.023TeV") || !energy.CompareTo("5.02TeV") || energy.Contains("5TeV2017")){
            if (!generator.CompareTo("")){
                return kOrange+1;
            } else if (!generator.CompareTo("LHC16h8a")) {
                return kBlue;
            } else if (!generator.CompareTo("LHC16h8b")) {
                return kRed;
            } else if ( generator.Contains("Data")) {
                return kGray+2;
            } else if (!generator.CompareTo("LHC16k5a")) {
                return kBlue;
            } else if (!generator.CompareTo("LHC16k5b")) {
                return kRed;
            } else if (!generator.CompareTo("LHC17e2")) {
                return kGreen+2;
            } else if (!generator.CompareTo("LHC18j3") || generator.Contains("Pythia8")) {
                return kRed+2;
            } else if (!generator.CompareTo("LHC17p_fast")){
                return kBlue+1;
            } else if (!generator.CompareTo("LHC17q_fast")){
                return kGreen+2;
            } else if (!generator.CompareTo("LHC17p_wSDD")){
                return kAzure+1;
            } else if (!generator.CompareTo("LHC17q_wSDD")){
                return kSpring+3;
            } else if (!generator.CompareTo("LHC17p_woSDD")){
                return kCyan+1;
            } else if (!generator.CompareTo("LHC17q_woSDD")){
                return kYellow+2;
            } else if (generator.Contains("LHC17pq_fast")){
                return kBlue+1;
            } else if (generator.Contains("LHC17pq_woSDD")){
                return kGreen+2;
            } else if (!generator.CompareTo("LHC17l3b_fast")){
                return kRed+1;
            } else if (!generator.CompareTo("LHC17l4b_fast")){
                return kViolet+2;
            } else if (!generator.CompareTo("LHC18j2_fast")){
                return kViolet+2;
            } else if (!generator.CompareTo("LHC17l3b_wSDD")){
                return kOrange+7;
            } else if (!generator.CompareTo("LHC17l4b_wSDD")){
                return kMagenta+2;
            } else if (!generator.CompareTo("LHC17l3b_woSDD")){
                return kRed-6;
            } else if (!generator.CompareTo("LHC17l4b_woSDD")){
                return kPink+6;
            } else if (!generator.CompareTo("LHC18j2_woSDD")){
                return kPink+6;
            } else if (generator.Contains("LHC17l3b")){
                return kOrange+7;
            } else {
                return kOrange-9;
            }
        }

        if (!energy.CompareTo("7TeV")){
            if (!kBox){
                if (!generator.CompareTo("")){
                    return kBlue+2;
                } else if (!generator.CompareTo("Pythia")){
                    return kBlue+5;
                } else if (!generator.CompareTo("Pythia2")){
                    return kBlue+4;
                } else if (!generator.CompareTo("PythiaAddSig")){
                    return kBlue-3;
                } else if (!generator.CompareTo("Phojet")){
                    return kBlue-1;
                } else if (!generator.CompareTo("LHC10_pass4")){
                    return kCyan+3;
                } else if (!generator.CompareTo("LHC10b_pass4") || !generator.CompareTo("LHC11a")){
                    return kBlack;
                } else if (!generator.CompareTo("LHC10c_pass4") || !generator.CompareTo("LHC11b")){
                    return 633;
                } else if (!generator.CompareTo("LHC10d_pass4") || !generator.CompareTo("LHC11c")){
                    return 807;
                } else if (!generator.CompareTo("LHC10e_pass4") || !generator.CompareTo("LHC11d")){
                    return 800;
                } else if (!generator.CompareTo("LHC10f_pass4")){
                    return 418;
                } else if (!generator.CompareTo("LHC14j4")){
                    return 852;
                } else if (!generator.CompareTo("LHC14j4b") || !generator.CompareTo("LHC14b7")){
                    return kGreen+4;
                } else if (!generator.CompareTo("LHC14j4c")){
                    return 435;
                } else if (!generator.CompareTo("LHC14j4d")){
                    return 601;
                } else if (!generator.CompareTo("LHC14j4e")){
                    return 879;
                } else if (!generator.CompareTo("LHC14j4f")){
                    return 806;
                } else if (!generator.CompareTo("LHC10_pass4_EMC")){
                    return kCyan+3;
                } else if (!generator.CompareTo("LHC10b_pass4_EMC")){
                    return kBlack;
                } else if (!generator.CompareTo("LHC10c_pass4_EMC")){
                    return 633;
                } else if (!generator.CompareTo("LHC10d_pass4_EMC")){
                    return 807;
                } else if (!generator.CompareTo("LHC10e_pass4_EMC")){
                    return 800;
                } else if (!generator.CompareTo("LHC10f_pass4_EMC")){
                    return 418;
                } else if (!generator.CompareTo("LHC14j4_EMC")){
                    return 852;
                } else if (!generator.CompareTo("LHC14j4b_EMC")){
                    return kBlack;
                } else if (!generator.CompareTo("LHC14j4c_EMC")){
                    return 633;
                } else if (!generator.CompareTo("LHC14j4d_EMC")){
                    return 807;
                } else if (!generator.CompareTo("LHC14j4e_EMC")){
                    return 800;
                } else if (!generator.CompareTo("LHC14j4f_EMC")){
                    return 418;
                }
            } else {
                return kBlue-5;
            }
        }
        if (energy.BeginsWith("8TeV")){
            if (!kBox){
                if (!generator.CompareTo("")){
                    return kGreen+2;
                } else if (!generator.CompareTo("Pythia")){
                    return kGreen+4;
                } else if (!generator.CompareTo("Pythia2")){
                    return kGreen+3;
                } else if (!generator.CompareTo("PythiaAddSig")){
                    return kGreen-4;
                } else if (!generator.CompareTo("Phojet")){
                    return kGreen-2;
                } else if (!generator.CompareTo("LHC12") || !generator.CompareTo("LHC12-kEMC7") || !generator.CompareTo("LHC12-kEMCEGA") || !generator.CompareTo("LHC12-kEMCEJE")){
                    return kCyan+3;
                } else if (generator.Contains("LHC12a")){
                    return kBlack;
                } else if (generator.Contains("LHC12b")){
                    return 633;
                } else if (generator.Contains("LHC12c")){
                    return 807;
                } else if (generator.Contains("LHC12d")){
                    return 800;
                } else if (generator.Contains("LHC12e")){
                    return kGreen+4;
                } else if (generator.Contains("LHC12f")){
                    return 418;
                } else if (generator.Contains("LHC12g")){
                    return 435;
                } else if (generator.Contains("LHC12h")){
                    return 601;
                } else if (generator.Contains("LHC12i")){
                    return 879;
                } else if (!generator.CompareTo("LHC14e2a") || generator.Contains("LHC15h1")){
                    return 852;
                } else if (!generator.CompareTo("LHC14e2b")){
                    return 426;
                } else if (!generator.CompareTo("LHC14e2c") || generator.Contains("LHC15h2")){
                    return 806;
                }
            } else {
                return kGreen-5;
            }
        }

        if (!energy.CompareTo("13TeV")||!energy.CompareTo("13TeVLowB")){
            if (!generator.CompareTo("LHC15f") || !generator.CompareTo("LHC15f_ESD") || !generator.CompareTo("LHC15g") || !generator.CompareTo("LHC15h") || !generator.CompareTo("LHC15i") || !generator.CompareTo("LHC15j") || !generator.CompareTo("LHC15k") || !generator.CompareTo("LHC15l")){
                return kBlack;
            } else if (!generator.CompareTo("LHC15f_AOD")) {
                return kOrange+7;
            } else if (!generator.CompareTo("LHC15g3a3") || !generator.CompareTo("LHC15g3a3_ESD") ) {
                return kBlue;
            } else if (!generator.CompareTo("LHC15g3a3_AOD") || generator.Contains("LHC18c_woSDD")) {
                return kMagenta+2;
            } else if (!generator.CompareTo("LHC15g3c3") || !generator.CompareTo("LHC15g3c3_ESD")) {
                return kRed;
            } else if (!generator.CompareTo("LHC15g3c3_AOD") || generator.Contains("LHC18c_wSDD")) {
                return kGreen+2;
            } else if (generator.Contains("LHC17k1") || generator.Contains("LHC17d1") || generator.Contains("LHC17h3")){
                return kGreen+2;
            } else if (generator.Contains("Pythia")) {
                return kCyan+2;
            } else if ((!generator.CompareTo("LHC17h11"))||(!generator.CompareTo("LHC17l5"))||(!generator.CompareTo("LHC17h1"))||(!generator.CompareTo("LHC17k4"))||(!generator.CompareTo("LHC18c13"))||(!generator.CompareTo("LHC18a8"))||(!generator.CompareTo("LHC18a9"))||(!generator.CompareTo("LHC18a1"))||(!generator.CompareTo("LHC18c12"))) {
                return kRed;
            } else if (generator.Contains("Phojet")) {
                return kGreen+3;
            } else if (generator.Contains("AOD")){
                return kRed;
            } else if (generator.Contains("FILTER")) {
                return kGreen;
            } else if (generator.Contains("ESD")){
                return kGreen-2;;
            } else if (generator.Contains("MC")){
                return kCyan+3;
            } else if (generator.Contains("LHC16d")||generator.Contains("LHC17c") ||generator.Contains("LHC18b")) {
                return kGreen-2;
            } else if (generator.Contains("LHC16e")||generator.Contains("LHC17e") ||generator.Contains("LHC18d")) {
                return 601;
            } else if (generator.Contains("LHC16f") ||generator.Contains("LHC18f")) {
                return kRed+1;
            } else if (generator.Contains("LHC16g")||generator.Contains("LHC17f") ||generator.Contains("LHC18g")) {
                return kCyan+3;
            } else if (generator.Contains("LHC16h")||generator.Contains("LHC17h") ||generator.Contains("LHC18h")) {
                return kBlack;
            } else if (generator.Contains("LHC16i")||generator.Contains("LHC17i") ||generator.Contains("LHC18i")) {
                return 633;
            } else if (generator.Contains("LHC16j")||generator.Contains("LHC17j") ||generator.Contains("LHC18j")) {
                return 807;
            } else if (generator.Contains("LHC16k")||generator.Contains("LHC17k") ||generator.Contains("LHC18k")) {
                return 800;
            } else if (generator.Contains("LHC16l")||generator.Contains("LHC17l") ||generator.Contains("LHC18l")) {
                return kGreen+4;
            } else if (generator.Contains("LHC16o")||generator.Contains("LHC17m") ||generator.Contains("LHC18m")) {
                return 418;
            } else if (generator.Contains("LHC16p")||generator.Contains("LHC17o") ||generator.Contains("LHC18n")) {
                return 435;
            } else if (generator.Contains("LHC17g") ||generator.Contains("LHC18o")) {
                return kCyan+3;
            } else if (generator.Contains("LHC17r") ||generator.Contains("LHC18p")) {
                return kMagenta+2;
            } else {
                return kRed;
            }
        }

        if (!energy.CompareTo("pPb_5.023TeV") || energy.CompareTo("pPb_5.023TeVCent") == 0 || !energy.CompareTo("pPb_5.023TeVRun2")){
            if (!kBox){
                if (!generator.CompareTo("")){
                    if (!centrality.CompareTo("0-10%")){
                        return kRed+1;
                    } else if (!centrality.CompareTo("0-20%")){
                        return kRed+1;
                    } else if (!centrality.CompareTo("0-5%")){
                        return kRed+3;
                    } else if (!centrality.CompareTo("5-10%")){
                        return 807;
                    } else if (!centrality.CompareTo("10-20%")){
                        return 800;
                    } else if (!centrality.CompareTo("20-40%")){
                        return kGreen+2;
                    } else if (!centrality.CompareTo("40-60%")){
                        return kCyan+2;
                    } else if (!centrality.CompareTo("60-80%")){
                        return kBlue+1;
                    } else if (!centrality.CompareTo("60-100%")){
                        return kViolet+4;
                    } else {
                        return kViolet+2;
                    }
                } else if (!generator.CompareTo("Hijing")){
                    return kViolet+1;
                } else if (!generator.CompareTo("DPMJET")){
                    return kViolet+6;
                } else if (!generator.CompareTo("LHC13b") || !generator.CompareTo("LHC13b_pass4")){
                    return 633;
                } else if (!generator.CompareTo("LHC13c")|| !generator.CompareTo("LHC13c_pass4")){
                    return 807;
                } else if (generator.Contains("LHC13b2_efix") || !generator.CompareTo("LHC18j5")){
                    return 418;
                } else if (generator.Contains("LHC13e7") || generator.Contains("LHC19a4")){
                    return 601;
                } else if (!generator.CompareTo("HijingAddSig")){
                    return kViolet-3;
                } else if (!generator.CompareTo("LHC13d") || !generator.CompareTo("LHC13d_pass4")){
                    return 800;
                } else if (!generator.CompareTo("LHC13e") || !generator.CompareTo("LHC13e_pass4")){
                    return kGreen+4;
                } else if (!generator.CompareTo("LHC13f") || !generator.CompareTo("LHC13f_pass4")){
                    return 418;
                } else if (generator.Contains("LHC16q_fast")){
                    return 597;
                } else if (generator.Contains("LHC16q_woSSD")){
                    return 419;
                } else if (generator.Contains("LHC16q")){
                    return 633;
                } else if (generator.Contains("LHC16t_fast")){
                    return 807;
                } else if (generator.Contains("LHC16t_woSSD")){
                    return 633;
                } else if (generator.Contains("LHC16t")){
                    return 807;
                } else if (generator.Contains("LHC17f2a") || generator.Contains("EPOS") || generator.Contains("LHC18f3")){
                    return kViolet+1;
                } else if (generator.Contains("LHC17f2b") || generator.Contains("DPMJET")){
                    return kGreen+2;
                } else {
                    if (!centrality.CompareTo("0-10%")){
                        return kRed+3;
                    } else if (!centrality.CompareTo("0-20%")){
                        return kRed+3;
                    } else if (!centrality.CompareTo("0-5%")){
                        return kRed+2;
                    } else if (!centrality.CompareTo("5-10%")){
                        return 809;
                    } else if (!centrality.CompareTo("10-20%")){
                        return 802;
                    } else if (!centrality.CompareTo("20-40%")){
                        return kGreen+4;
                    } else if (!centrality.CompareTo("40-60%")){
                        return kCyan+3;
                    } else if (!centrality.CompareTo("60-80%")){
                        return kBlue+3;
                    } else if (!centrality.CompareTo("60-100%")){
                        return kViolet+3;
                    } else {
                        return kViolet+1;
                    }
                }
            } else {
                if (!centrality.CompareTo("0-10%")){
                    return kRed-6;
                } else if (!centrality.CompareTo("0-20%")){
                    return kRed-5;
                } else if (!centrality.CompareTo("10-20%")){
                    return 800-3;
                } else if (!centrality.CompareTo("20-40%")){
                    return kGreen-5;
                } else if (!centrality.CompareTo("40-60%")){
                    return kCyan-5;
                } else if (!centrality.CompareTo("60-100%")){
                    return kViolet+2;
                } else {
                    return kViolet+6;
                }
            }
        }
        if (energy.Contains("pPb_8TeV") ){
            if (!kBox){
                if (!generator.CompareTo("")){
                    return kViolet+2;
                } else if (!generator.CompareTo("Hijing")){
                    return kViolet+1;
                } else if (!generator.CompareTo("DPMJET")){
                    return kViolet+6;
                } else if (!generator.CompareTo("LHC16r")){
                    return 633;
                } else if (!generator.CompareTo("LHC16s")){
                    return 807;
                } else if (generator.Contains("LHC17f3a") || generator.Contains("LHC18f3b")){
                    return 418;
                } else if (generator.Contains("LHC17f4a") || generator.Contains("LHC18f3c")){
                    return 601;
                } else if (generator.Contains("LHC17f3b")){
                    return kViolet-3;
                } else if (generator.Contains("LHC17f4b")){
                    return kViolet+2;
                } else if (generator.Contains("LHC16r") && generator.Contains("woSDD")){
                    return 633;
                } else if (generator.Contains("LHC16r") && (generator.Contains("wSDD")|| generator.Contains("mcp2"))){
                    return 807;
                } else if (generator.Contains("LHC16r") && (generator.Contains("fast") || generator.Contains("FAST")) ){
                    return 800;
                } else if (generator.Contains("LHC16s") && generator.Contains("woSDD")){
                    return kGreen+2;
                } else if (generator.Contains("LHC16s") && (generator.Contains("wSDD")|| generator.Contains("mcp2"))){
                    return 418;
                } else if (generator.Contains("LHC16s") && (generator.Contains("fast") || generator.Contains("FAST")) ){
                    return 633;
                } else if (generator.Contains("LHC16r")){
                    return 633;
                } else if (generator.Contains("LHC16s")){
                    return 807;
                } else if (generator.Contains("LHC17g8b") || generator.Contains("EPOS")){
                    return kViolet+1;
                } else if (generator.Contains("LHC17g8c") || generator.Contains("DPMJET")){
                    return kGreen+2;
                }
            } else {
                return kViolet+6;
            }
        }

        if (!energy.CompareTo("PbPb_2.76TeV")){
            if (!kBox){
                if (!generator.CompareTo("")){
                    if (!centrality.CompareTo("0-10%")){
                        return kRed+1;
                    } else if (!centrality.CompareTo("0-20%")){
                        return kRed+1;
                    } else if (!centrality.CompareTo("0-40%")){
                        return kMagenta+2;
                    } else if (!centrality.CompareTo("0-5%")){
                        return kRed+1;
                    } else if (!centrality.CompareTo("5-10%")){
                        return 807;
                    } else if (!centrality.CompareTo("10-20%")|| !centrality.CompareTo("10-30%") ){
                        return 800;
                    } else if (!centrality.CompareTo("20-40%")){
                        return kGreen+2;
                    } else if (!centrality.CompareTo("20-50%")){
                        return kCyan+2;
                    } else if (!centrality.CompareTo("40-60%")){
                        return kCyan+2;
                    } else if (!centrality.CompareTo("60-80%") || !centrality.CompareTo("60-92%")){
                        return kBlue+1;
                    } else if (!centrality.CompareTo("40-80%") || !centrality.CompareTo("0-94%")){
                        return kCyan+2;
                    } else if (!centrality.CompareTo("30-100%") ){
                        return kViolet+2;
                    }
                } else {
                    if (!centrality.CompareTo("0-10%")){
                        return kRed+3;
                    } else if (!centrality.CompareTo("0-20%")){
                        return kRed+3;
                    } else if (!centrality.CompareTo("0-40%")){
                        return kMagenta+4;
                    } else if (!centrality.CompareTo("0-5%")){
                        return kRed+2;
                    } else if (!centrality.CompareTo("5-10%")){
                        return 809;
                    } else if (!centrality.CompareTo("10-20%") || !centrality.CompareTo("10-30%") ){
                        return 802;
                    } else if (!centrality.CompareTo("20-40%")){
                        return kGreen+4;
                    } else if (!centrality.CompareTo("20-50%")){
                        return kCyan+3;
                    } else if (!centrality.CompareTo("40-60%") ){
                        return kCyan+3;
                    } else if (!centrality.CompareTo("60-80%") || !centrality.CompareTo("60-92%")){
                        return kBlue+3;
                    } else if (!centrality.CompareTo("40-80%") || !centrality.CompareTo("0-94%")){
                        return kCyan+3;
                    } else if (!centrality.CompareTo("30-100%") ){
                        return kViolet-6;
                    }
                }
            } else {
                if (!generator.CompareTo("")){
                    if (!centrality.CompareTo("0-10%")){
                        return kRed-6;
                    } else if (!centrality.CompareTo("0-20%")){
                        return kRed-5;
                    } else if (!centrality.CompareTo("10-20%")){
                        return 800-3;
                    } else if (!centrality.CompareTo("20-40%")){
                        return kGreen-5;
                    } else if (!centrality.CompareTo("20-50%")){
                        return kCyan-5;
                    } else if (!centrality.CompareTo("40-80%") || !centrality.CompareTo("0-94%")){
                        return kCyan-5;
                    }
                }
            }
        }
        if (!energy.CompareTo("PbPb_5.02TeV")){
            if(!centrality.CompareTo("")){ // generator=periodname, centrality="" (for QA)
                if (!generator.CompareTo("LHC15o") || !generator.CompareTo("LHC15o_red")){
                    return kBlack;
                } else if (!generator.CompareTo("LHC18q") || !generator.CompareTo("LHC18r")){
                    return kBlack;
                } else if (!generator.CompareTo("LHC16g1") || !generator.CompareTo("LHC16h4") || !generator.CompareTo("LHC18e1")){
                    return kGray+2;
                } else if (generator.Contains("LHC18l8")){ //"LHC18l8a4", "LHC18l8b2", "LHC18l8b2_cent2mb", "LHC18l8c2", "LHC18l8c2_mb2cent"
                  if (!generator.CompareTo("LHC18l8a4")){
                    return kGray+2;
                  } else if (!generator.CompareTo("LHC18l8b2")){
                    return kRed+2;
                  } else if (!generator.CompareTo("LHC18l8b2_cent2mb")){
                    return kGreen+2;
                  } else if (!generator.CompareTo("LHC18l8c2")){
                    return kBlue+3;
                  } else if (!generator.CompareTo("LHC18l8c2_mb2cent")){
                    return kMagenta+2;
                  }
                } else if (generator.Contains("a")){
                    return kRed;
                } else if (generator.Contains("b")){
                    return kGreen+2;
                } else if (generator.Contains("c")){
                    return kBlue+3;
                }
            } else if (!generator.CompareTo("")){  // generator="" (data),  centrality defined (for combination)
                if (!centrality.CompareTo("0-10%"))      return kRed+1;
                else if (!centrality.CompareTo("0-20%")) return kRed+1;
                else if (!centrality.CompareTo("0-40%")) return kMagenta+2;
                else if (!centrality.CompareTo("0-5%"))  return kRed+1;
                else if (!centrality.CompareTo("5-10%")) return 807;
                else if (!centrality.CompareTo("10-20%")) return 800;
                else if (!centrality.CompareTo("20-30%") ) return kSpring+5;
                else if (!centrality.CompareTo("10-30%") ) return kYellow-3;
                else if (!centrality.CompareTo("20-40%") || !centrality.CompareTo("30-40%")) return kGreen+2;
                else if (!centrality.CompareTo("20-50%")) return kCyan+2;
                else if (!centrality.CompareTo("30-50%") || !centrality.CompareTo("40-50%")) return kTeal+4;
                else if (!centrality.CompareTo("40-60%") || !centrality.CompareTo("50-60%")) return kCyan+2;
                else if (!centrality.CompareTo("60-70%") ) return kAzure+5;
                else if (!centrality.CompareTo("70-80%") ) return kAzure-1;
                else if (!centrality.CompareTo("60-80%") || !centrality.CompareTo("60-92%")) return kBlue+1;
                else if (!centrality.CompareTo("40-80%") || !centrality.CompareTo("0-94%")) return kCyan+2;
                else if (!centrality.CompareTo("50-90%") )   return kViolet+2;
                else if (!centrality.CompareTo("80-90%") )   return kViolet+4;
                else if (!centrality.CompareTo("30-100%") )  return kViolet+2;
            } else { // generator defined (MC) and centrality defined (for combination)
                if (!centrality.CompareTo("0-10%"))          return kRed+3;
                else if (!centrality.CompareTo("0-20%"))     return kRed+3;
                else if (!centrality.CompareTo("0-40%"))     return kMagenta+4;
                else if (!centrality.CompareTo("0-5%"))      return kRed+2;
                else if (!centrality.CompareTo("5-10%"))     return 809;
                else if (!centrality.CompareTo("10-20%") )   return 802;
                else if (!centrality.CompareTo("20-30%") )   return kSpring+4;
                else if (!centrality.CompareTo("10-30%") )   return kYellow-2;
                else if (!centrality.CompareTo("20-40%") || !centrality.CompareTo("30-40%"))    return kGreen+4;
                else if (!centrality.CompareTo("20-50%"))    return kCyan+3;
                else if (!centrality.CompareTo("30-50%") || !centrality.CompareTo("40-50%"))    return kTeal-7;
                else if (!centrality.CompareTo("40-60%") || !centrality.CompareTo("50-60%"))   return kCyan+3;
                else if (!centrality.CompareTo("60-70%") )   return kAzure+4;
                else if (!centrality.CompareTo("70-80%") )   return kAzure-6;
                else if (!centrality.CompareTo("60-80%") || !centrality.CompareTo("60-92%"))  return kBlue+3;
                else if (!centrality.CompareTo("40-80%") || !centrality.CompareTo("0-94%"))   return kCyan+3;
                else if (!centrality.CompareTo("50-90%") )   return kViolet-6;
                else if (!centrality.CompareTo("80-90%") )   return kViolet+3;
                else if (!centrality.CompareTo("30-100%") )  return kViolet-6;
            }
        }

        if (!energy.CompareTo("XeXe_5.44TeV")){
            if (!generator.CompareTo("LHC17n") || !generator.CompareTo("") ){
                if (!centrality.CompareTo("0-20%")|| !centrality.CompareTo("0-10%")){
                    return kRed+1;
                } else if (!centrality.CompareTo("20-40%")){
                    return kGreen+2;
                } else if (!centrality.CompareTo("0-40%")){
                    return kOrange+7;
                } else if (!centrality.CompareTo("40-80%")){
                    return kBlue+2;
                } else if (!centrality.CompareTo("0-80%")){
                    return kBlack;
                } else if (!centrality.CompareTo("60-90%")){
                  return kBlack;
                } else return kGray+1;
            } else if (!generator.CompareTo("LHC17j7") || !generator.CompareTo("HIJING") || !generator.CompareTo("Hijing") || !generator.CompareTo("LHC18d2")){
                if (!centrality.CompareTo("0-20%") || !centrality.CompareTo("0-10%") ){
                    return kRed-6;
                } else if (!centrality.CompareTo("20-40%")){
                    return kGreen-6;
                } else if (!centrality.CompareTo("0-40%")){
                    return kOrange-3;
                } else if (!centrality.CompareTo("40-80%")){
                    return kBlue-6;
                } else if (!centrality.CompareTo("0-80%")){
                    return kGray+2;
                } else if (!centrality.CompareTo("60-90%")){
                  return kGray+2;
                } else return kGray+2;
            }
        }
        cout << "GetColorDefaultColor(): no valid input parameters given..." << endl;
        return kBlack;
    }

    Style_t GetDefaultMarkerStyle(  TString energy,
                                    TString generator,
                                    TString centrality){
        if (!energy.CompareTo("900GeV")){
            if (!generator.CompareTo("")||!generator.CompareTo("LHC10c_900GeV_pass4")){
                return 21;
            } else {
                return 25;
            }
        }
        if (!energy.CompareTo("2.76TeV")){
            if (!generator.CompareTo("")||!generator.CompareTo("LHC11a_p4_wSDD")|| !generator.CompareTo("LHC11a_pass4_wSDD")){
                return 29;
            } else if(!generator.CompareTo("LHC13g")||!generator.CompareTo("LHC13g-kEMC7")|
                    !generator.CompareTo("LHC13g-kEMCEG1")||!generator.CompareTo("LHC13g-kEMCEG2") || !generator.CompareTo("LHC13g_pass1")){
                return 29;
            } else if(!generator.CompareTo("LHC12f1a")||!generator.CompareTo("LHC12f1b")||!generator.CompareTo("LHC12i3")){
                return 30;
            } else {
                return 30;
            }
        }
        if (!energy.CompareTo("5TeV") || !energy.CompareTo("5.023TeV") || !energy.CompareTo("5.02TeV") || energy.Contains("5TeV2017")){
            if (!generator.CompareTo("") || generator.Contains("Pythia8")){
                return 34;
            } else if(generator.Contains("LHC16h8a")){
                return 30;
            } else if(generator.Contains("LHC16h8b")){
                return 28;
            } else if(generator.Contains("LHC16k5a")){
                return 30;
            } else if(generator.Contains("LHC16k5b")){
                return 28;
            } else if(generator.Contains("LHC17e2") || generator.Contains("Data")){
                return 28;
            } else if(generator.Contains("LHC18j3")){
                return 20;
            } else if(generator.Contains("LHC17p")){
                return 20;
            } else if(generator.Contains("LHC17q")){
                return 20;
            } else if(generator.Contains("LHC17l3b")){
                return 28;
            } else if(generator.Contains("LHC17l4b")){
                return 29;
            } else if(generator.Contains("LHC18j2")){
                return 29;
            } else {
                return 34;
            }
        }
        if (!energy.CompareTo("7TeV")){
            if (!generator.CompareTo("")){
                return 20;
            } else if( ( generator.Contains("LHC10") && generator.Contains("_pass4") ) || generator.Contains("LHC11")){
                return 29;
            } else if(generator.Contains("LHC14j4") || !generator.CompareTo("LHC14b7")){
                return 30;
            } else {
                return 24;
            }
        }
        if (energy.BeginsWith("8TeV")){
            if (!generator.CompareTo("")){
                return 33;
            } else if(generator.CompareTo("LHC12")==0 || (generator.BeginsWith("LHC12")&&generator.Length()==6)
                    || (generator.BeginsWith("LHC12")&&generator.EndsWith("-kEMC7"))
                    || (generator.BeginsWith("LHC12")&&generator.EndsWith("-kEMCEGA"))
                    || (generator.BeginsWith("LHC12")&&generator.EndsWith("-kEMC8"))
                    || (generator.BeginsWith("LHC12")&&generator.EndsWith("-kEMC8EGA"))
                    || (generator.BeginsWith("LHC12")&&generator.EndsWith("-kEMCEJE"))){
                return 29;
            } else if(!generator.CompareTo("LHC14e2a") || generator.Contains("LHC15h1")){
                return 30;
            } else if(!generator.CompareTo("LHC14e2b")){
                return 27;
            } else if(!generator.CompareTo("LHC14e2c") || generator.Contains("LHC15h2")){
                return 28;
            } else {
                return 27;
            }
        }

        if (!energy.CompareTo("13TeV")){
            if (!generator.CompareTo("")) {
                return 33;
            } else if (!generator.CompareTo("LHC15f") || !generator.CompareTo("LHC15f_ESD")) {
                return 24;
            } else if (!generator.CompareTo("LHC15f_AOD")) {
                return 5;
            } else if (!generator.CompareTo("LHC15g3a3") || !generator.CompareTo("LHC15g3a3_ESD")) {
                return 28;
            } else if (!generator.CompareTo("LHC15g3a3_AOD")) {
                return 25;
            } else if (!generator.CompareTo("LHC15g3c3") || !generator.CompareTo("LHC15g3c3_ESD")) {
                return 31;
            } else if (!generator.CompareTo("LHC15g3c3_AOD")) {
                return 27;
            } else if (generator.Contains("Pythia")) {
                return 28;
            } else if (generator.Contains("Phojet")) {
                return 27;
            } else if ((generator.Contains("LHC16d")) || (generator.Contains("LHC16e")) || (generator.Contains("LHC16g")) || (generator.Contains("LHC16h")) || (generator.Contains("LHC16i")) || (generator.Contains("LHC16j")) || (generator.Contains("LHC16k")) || (generator.Contains("LHC16l")) || (generator.Contains("LHC16o")) || (generator.Contains("LHC16p"))) {
                return 24;
            } else if ((generator.Contains("LHC17c")) || (generator.Contains("LHC17e")) || (generator.Contains("LHC17f")) || (generator.Contains("LHC17h")) || (generator.Contains("LHC17i")) ||
                        (generator.Contains("LHC17j")) || (generator.Contains("LHC17k")) || (generator.Contains("LHC17l")) || (generator.Contains("LHC17m")) || (generator.Contains("LHC17o")) ||
                        (generator.Contains("LHC17r"))) {
                return 25;
            } else if ( (generator.Contains("LHC18b")) || (generator.Contains("LHC18d"))|| (generator.Contains("LHC18e")) || (generator.Contains("LHC18f")) || (generator.Contains("LHC18h")) ||
                        (generator.Contains("LHC18i")) || (generator.Contains("LHC18j")) || (generator.Contains("LHC18k")) || (generator.Contains("LHC18l")) || (generator.Contains("LHC18m")) ||
                        (generator.Contains("LHC18n")) || (generator.Contains("LHC18o")) || (generator.Contains("LHC18p"))) {
                return 26;
            } else if ((generator.Contains("LHC16")) || (generator.Contains("LHC17")) || (generator.Contains("LHC18"))) {
                return 24;
            } else if (generator.Contains("MC")) {
                return 25;
            } else {
                return 34;
            }
        }

        if (!energy.CompareTo("13TeVLowB")){
            if (!generator.CompareTo("")) {
                return 33;
            } else if (generator.Contains("LHC18c")) {
                return 24;
            } else if (generator.Contains("Phytia")) {
                return 30;
            } else if (generator.Contains("Phojet")) {
                return 27;
            } else if (generator.Contains("LHC17k1")) {
                return 28;
            } else if (generator.Contains("LHC17d1") || generator.Contains("LHC17h3")) {
                return 28;
            } else if ((generator.Contains("LHC16d")) || (generator.Contains("LHC16e")) || (generator.Contains("LHC16g")) || (generator.Contains("LHC16h")) || (generator.Contains("LHC16i")) || (generator.Contains("LHC16j")) || (generator.Contains("LHC16k")) || (generator.Contains("LHC16l")) || (generator.Contains("LHC16o")) || (generator.Contains("LHC16p"))) {
                return 29;
            } else {
                return 34;
            }
        }

        if (!energy.CompareTo("pPb_5.023TeV") || energy.CompareTo("pPb_5.023TeVCent") == 0 || !energy.CompareTo("pPb_5.023TeVRun2")){
            if (!generator.CompareTo("")){
                if (!centrality.CompareTo("0-10%")){
                  return 20;
                } else if (!centrality.CompareTo("0-20%")){
                  return 20;
                } else if (!centrality.CompareTo("0-40%")){
                  return 20;
                } else if (!centrality.CompareTo("0-5%")){
                  return 20;
                } else if (!centrality.CompareTo("5-10%")){
                  return 21;
                } else if (!centrality.CompareTo("10-20%")){
                  return 29;
                } else if (!centrality.CompareTo("20-40%")){
                  return 33;
                } else if (!centrality.CompareTo("20-50%")){
                  return 34;
                } else if (!centrality.CompareTo("40-60%")){
                  return 20;
                } else if (!centrality.CompareTo("60-80%")){
                  return 21;
                } else if (!centrality.CompareTo("60-100%")){
                  return 34;
                } else {
                  return 33;
                }
            } else if (!generator.CompareTo("LHC13b") || !generator.CompareTo("LHC13b_pass4")){
                return 20;
            } else if (!generator.CompareTo("LHC13c") || !generator.CompareTo("LHC13c_pass4")){
                return 21;
            } else if (!generator.CompareTo("LHC13d") || !generator.CompareTo("LHC13d_pass4")){
                return 20;
            } else if (!generator.CompareTo("LHC13e") || !generator.CompareTo("LHC13e_pass4")){
                return 21;
            } else if (!generator.CompareTo("LHC13f") || !generator.CompareTo("LHC13f_pass4")){
                return 20;
            } else if (generator.Contains("LHC13b2_efix") || !generator.CompareTo("LHC18j5")){
                return 26;
            } else if (generator.Contains("LHC13e7") || generator.Contains("LHC19a4")){
                return 28;
            } else if (generator.Contains("LHC16q")){
                return 20;
            } else if (generator.Contains("LHC16t")){
                return 21;
            } else if (generator.Contains("LHC17f2a") || generator.Contains("LHC18f3")){
                return 24;
            } else if (generator.Contains("LHC17f2b")){
                return 25;
            } else {
              if (!centrality.CompareTo("0-10%")){
                return 24;
              } else if (!centrality.CompareTo("0-20%")){
                return 24;
              } else if (!centrality.CompareTo("0-40%")){
                return 24;
              } else if (!centrality.CompareTo("0-5%")){
                return 24;
              } else if (!centrality.CompareTo("5-10%")){
                return 25;
              } else if (!centrality.CompareTo("10-20%")){
                return 30;
              } else if (!centrality.CompareTo("20-40%")){
                return 27;
              } else if (!centrality.CompareTo("20-50%")){
                return 28;
              } else if (!centrality.CompareTo("40-60%")){
                return 24;
              } else if (!centrality.CompareTo("60-80%")){
                return 25;
              } else if (!centrality.CompareTo("60-100%")){
                return 28;
              } else {
                return 27;
              }
            }
        }
        if (energy.Contains("pPb_8TeV") ){
            if (!generator.CompareTo("")){
                return 33;
            } else if (generator.Contains("LHC16r")){
                return 24;
            } else if (generator.Contains("LHC16s")){
                return 27;
            } else if (generator.Contains("LHC17f3") || generator.Contains("LHC18f3b")){
                return 30;
            } else if (generator.Contains("LHC17f4") || generator.Contains("LHC18f3c")){
                return 26;
            } else if (generator.Contains("LHC17g8b")){
                return 24;
            } else if (generator.Contains("LHC17g8c")){
                return 25;
            } else {
                return 27;
            }
        }

        if (!energy.CompareTo("PbPb_2.76TeV")){
            if (!generator.CompareTo("")){
                if (!centrality.CompareTo("0-10%")){
                    return 20;
                } else if (!centrality.CompareTo("0-20%")){
                    return 20;
                } else if (!centrality.CompareTo("0-40%")){
                    return 20;
                } else if (!centrality.CompareTo("0-5%")){
                    return 20;
                } else if (!centrality.CompareTo("5-10%")){
                    return 21;
                } else if (!centrality.CompareTo("10-20%") || !centrality.CompareTo("10-30%") ){
                    return 29;
                } else if (!centrality.CompareTo("20-40%")){
                    return 33;
                } else if (!centrality.CompareTo("20-50%")){
                    return 34;
                } else if (!centrality.CompareTo("40-60%")){
                    return 20;
                } else if (!centrality.CompareTo("30-100%") ){
                    return 21;
                } else if (!centrality.CompareTo("60-80%") || !centrality.CompareTo("60-92%")){
                    return 21;
                } else if (!centrality.CompareTo("40-80%") || !centrality.CompareTo("0-94%") || !centrality.CompareTo("0-86%") || !centrality.CompareTo("0-92%") || !centrality.CompareTo("0-80%") ){
                    return 34;
                }
            } else {
                if (!centrality.CompareTo("0-10%")){
                    return 24;
                } else if (!centrality.CompareTo("0-20%")){
                    return 24;
                } else if (!centrality.CompareTo("0-40%")){
                    return 24;
                } else if (!centrality.CompareTo("0-5%")){
                    return 24;
                } else if (!centrality.CompareTo("5-10%")){
                    return 25;
                } else if (!centrality.CompareTo("10-20%")|| !centrality.CompareTo("10-30%") ){
                    return 30;
                } else if (!centrality.CompareTo("20-40%")){
                    return 27;
                } else if (!centrality.CompareTo("20-50%")){
                    return 28;
                } else if (!centrality.CompareTo("40-60%")){
                    return 24;
                } else if (!centrality.CompareTo("30-100%") ){
                    return 25;
                } else if (!centrality.CompareTo("60-80%") || !centrality.CompareTo("60-92%")){
                    return 25;
                } else if (!centrality.CompareTo("40-80%") || !centrality.CompareTo("0-94%") || !centrality.CompareTo("0-86%") || !centrality.CompareTo("0-92%") || !centrality.CompareTo("0-80%") ){
                    return 28;
                }
            }
        }
        if (!energy.CompareTo("PbPb_5.02TeV")){

            if (!centrality.CompareTo("")){  // generator=periodname, centrality="" (for QA)
                if (!generator.CompareTo("LHC15o") || !generator.CompareTo("LHC15o_red")){
                    return 20;
                } else if (!generator.CompareTo("LHC18q")){
                    return 20;
                } else if (!generator.CompareTo("LHC18r")){
                    return 21;
                } else if (!generator.CompareTo("LHC16g1") || !generator.CompareTo("LHC16h4") || !generator.CompareTo("LHC18e1")){
                    return 24;
                } else if (generator.Contains("LHC18l8")){ //"LHC18l8a4", "LHC18l8b2", "LHC18l8b2_cent2mb", "LHC18l8c2", "LHC18l8c2_mb2cent"
                  if (!generator.CompareTo("LHC18l8a4")){
                    return 24;
                  } else if (!generator.CompareTo("LHC18l8b2")){
                    return 25;
                  } else if (!generator.CompareTo("LHC18l8b2_cent2mb")){
                    return 27;
                  } else if (!generator.CompareTo("LHC18l8c2")){
                    return 28;
                  } else if (!generator.CompareTo("LHC18l8c2_mb2cent")){
                    return 30;
                  }
                } else if (generator.Contains("a")){
                    return 25;
                } else if (generator.Contains("b")){
                    return 28;
                } else if (generator.Contains("c")){
                    return 27;
                } else
                    return 5;
            } else if (!generator.CompareTo("")){ // generator="" (data),  centrality defined (for combination)
                if (!centrality.CompareTo("0-10%"))         return 20;
                else if (!centrality.CompareTo("0-20%"))  return 20;
                else if (!centrality.CompareTo("0-40%"))  return 20;
                else if (!centrality.CompareTo("0-5%"))   return 20;
                else if (!centrality.CompareTo("5-10%"))  return 21;
                else if (!centrality.CompareTo("10-20%")) return 29;
                else if (!centrality.CompareTo("20-30%") ) return 47;
                else if (!centrality.CompareTo("10-30%") ) return 47;
                else if (!centrality.CompareTo("20-40%") || !centrality.CompareTo("30-40%")) return 33;
                else if (!centrality.CompareTo("40-50%")) return 34;
                else if (!centrality.CompareTo("50-60%")) return 41;
                else if (!centrality.CompareTo("60-70%")) return 49;
                else if (!centrality.CompareTo("70-80%")) return 45;
                else if (!centrality.CompareTo("80-90%")) return 43;
                else if (!centrality.CompareTo("20-50%")) return 34;
                else if (!centrality.CompareTo("30-50%")) return 34;
                else if (!centrality.CompareTo("40-60%")) return 20;
                else if (!centrality.CompareTo("30-100%") ) return 21;
                else if (!centrality.CompareTo("50-90%"))   return 43;
                else if (!centrality.CompareTo("60-80%") || !centrality.CompareTo("60-92%")) return 21;
                else if (!centrality.CompareTo("40-80%") || !centrality.CompareTo("0-94%") || !centrality.CompareTo("0-86%") || !centrality.CompareTo("0-92%") || !centrality.CompareTo("0-80%") ) return 34;
            } else { // generator defined (MC) and centrality defined (for combination)
                if (!centrality.CompareTo("0-10%"))         return 24;
                else if (!centrality.CompareTo("0-20%"))    return 24;
                else if (!centrality.CompareTo("0-40%"))    return 24;
                else if (!centrality.CompareTo("0-5%"))     return 24;
                else if (!centrality.CompareTo("5-10%"))    return 25;
                else if (!centrality.CompareTo("10-20%"))   return 30;
                else if (!centrality.CompareTo("20-30%"))   return 46;
                else if (!centrality.CompareTo("10-30%"))   return 46;
                else if (!centrality.CompareTo("20-40%") || !centrality.CompareTo("30-40%"))   return 27;
                else if (!centrality.CompareTo("40-50%")) return 28;
                else if (!centrality.CompareTo("50-60%")) return 40;
                else if (!centrality.CompareTo("60-70%")) return 28;
                else if (!centrality.CompareTo("70-80%")) return 44;
                else if (!centrality.CompareTo("80-90%")) return 42;
                else if (!centrality.CompareTo("20-50%"))   return 28;
                else if (!centrality.CompareTo("30-50%"))   return 28;
                else if (!centrality.CompareTo("40-60%"))   return 24;
                else if (!centrality.CompareTo("30-100%") ) return 25;
                else if (!centrality.CompareTo("50-90%"))   return 42;
                else if (!centrality.CompareTo("60-80%") || !centrality.CompareTo("60-92%")) return 25;
                else if (!centrality.CompareTo("40-80%") || !centrality.CompareTo("0-94%") || !centrality.CompareTo("0-86%") || !centrality.CompareTo("0-92%") || !centrality.CompareTo("0-80%") )
                    return 28;
            }
        }
        if (!energy.CompareTo("XeXe_5.44TeV")){
            if (!generator.CompareTo("LHC17n") || !generator.CompareTo("")){
                if (!centrality.CompareTo("0-80%")){
                    return kFullSquare;
                } else if (!centrality.CompareTo("0-20%") || !centrality.CompareTo("0-10%")){
                    return kFullCircle;
                } else if (!centrality.CompareTo("20-40%")){
                    return kFullDiamond;
                } else if (!centrality.CompareTo("0-40%")){
                    return kFullCross;
                } else if (!centrality.CompareTo("40-80%")){
                    return kFullStar;
                } else if (!centrality.CompareTo("60-90%")){
                  return kFullStar;
                } else return kFullSquare;
            } else if (!generator.CompareTo("LHC17j7") || !generator.CompareTo("HIJING") || !generator.CompareTo("Hijing") || !generator.CompareTo("LHC18d2") ){
                if (!centrality.CompareTo("0-80%")){
                    return kOpenSquare;
                } else if (!centrality.CompareTo("0-20%") || !centrality.CompareTo("0-10%")){
                    return kOpenCircle;
                } else if (!centrality.CompareTo("20-40%")){
                    return kOpenDiamond;
                } else if (!centrality.CompareTo("0-40%")){
                    return kOpenCross;
                } else if (!centrality.CompareTo("40-80%")){
                    return kOpenStar;
                } else if (!centrality.CompareTo("60-90%")){
                  return kOpenStar;
                } else return kOpenSquare;
            }
        }
        cout << "GetDefaultMarkerStyle(): no valid input parameters given..." << endl;
        return 0;
    }

    Size_t GetDefaultMarkerSize( TString energy,
                                TString generator,
                                TString centrality){

        if (!energy.CompareTo("900GeV")){
            if (!generator.CompareTo("")){
                return 1.8;
            } else {
                return 1.8;
            }
        }

        if (!energy.CompareTo("2.76TeV")){
            if (!generator.CompareTo("")){
                return 2.2;
            } else {
                return 2.2;
            }
        }

        if (!energy.CompareTo("5TeV") || !energy.CompareTo("5.023TeV") || !energy.CompareTo("5.02TeV") || energy.Contains("5TeV2017")){
            if (!generator.CompareTo("")){
                return 2;
            } else {
                return 2.2;
            }
        }

        if (!energy.CompareTo("7TeV")){
            if (!generator.CompareTo("")){
                return 2.2;
            } else {
                return 2.2;
            }
        }

        if (energy.BeginsWith("8TeV")){
            if (!generator.CompareTo("")){
                return 2.2;
            } else {
                return 2.2;
            }
        }

        if (energy.Contains("13TeV")){
            return 2.2;
        }

        if (!energy.CompareTo("pPb_5.023TeV") || energy.CompareTo("pPb_5.023TeVCent") == 0 || !energy.CompareTo("pPb_5.023TeVRun2")){
            if (!generator.CompareTo("")){
              if (!centrality.CompareTo("0-10%")){
                return 2;
              } else if (!centrality.CompareTo("0-20%")){
                return 2;
              } else if (!centrality.CompareTo("0-40%")){
                return 2;
              } else if (!centrality.CompareTo("0-5%")){
                return 2;
              } else if (!centrality.CompareTo("5-10%")){
                return 2;
              } else if (!centrality.CompareTo("10-20%")){
                return 2.5;
              } else if (!centrality.CompareTo("20-40%")){
                return 2.5;
              } else if (!centrality.CompareTo("20-50%")){
                return 2.;
              } else if (!centrality.CompareTo("40-60%")){
                return 2;
              } else if (!centrality.CompareTo("60-80%")){
                return 2;
              } else if (!centrality.CompareTo("60-100%")){
                return 2;
              } else {
                return 2.2;
              }
            } else if (generator.CompareTo("LHC13")==0 || (generator.BeginsWith("LHC13")&&generator.Length()==6)){
                return 2.2;
            } else {
                return 2.2;
            }
        }

        if (!energy.CompareTo("PbPb_2.76TeV")){
            if (!generator.CompareTo("")){
                if (!centrality.CompareTo("0-10%")){
                    return 2;
                } else if (!centrality.CompareTo("0-20%")){
                    return 2;
                } else if (!centrality.CompareTo("0-40%")){
                    return 2;
                } else if (!centrality.CompareTo("0-5%")){
                    return 2;
                } else if (!centrality.CompareTo("5-10%")){
                    return 2;
                } else if (!centrality.CompareTo("10-20%") || !centrality.CompareTo("10-30%") ){
                    return 2.5;
                } else if (!centrality.CompareTo("20-40%")){
                    return 2.5;
                } else if (!centrality.CompareTo("20-50%")){
                    return 2.;
                } else if (!centrality.CompareTo("40-60%") || !centrality.CompareTo("30-100%") ){
                    return 2;
                } else if (!centrality.CompareTo("60-80%") || !centrality.CompareTo("60-92%")){
                    return 2;
                } else if (!centrality.CompareTo("40-80%") || !centrality.CompareTo("0-94%")){
                    return 2;
                } else {
                    return 2;
                }
            }
        }

        if (!energy.CompareTo("PbPb_5.02TeV")){
            if (!generator.CompareTo("")){
                if (!centrality.CompareTo("0-10%")){
                    return 2;
                } else if (!centrality.CompareTo("0-20%")){
                    return 2;
                } else if (!centrality.CompareTo("0-40%")){
                    return 2;
                } else if (!centrality.CompareTo("0-5%")){
                    return 2;
                } else if (!centrality.CompareTo("5-10%")){
                    return 2;
                } else if (!centrality.CompareTo("10-20%") || !centrality.CompareTo("10-30%") ){
                    return 2.5;
                } else if (!centrality.CompareTo("20-40%")){
                    return 2.5;
                } else if (!centrality.CompareTo("20-50%")){
                    return 2.;
                } else if (!centrality.CompareTo("40-60%") || !centrality.CompareTo("30-100%") ){
                    return 2;
                } else if (!centrality.CompareTo("60-80%") || !centrality.CompareTo("60-92%")){
                    return 2;
                } else if (!centrality.CompareTo("40-80%") || !centrality.CompareTo("0-94%")){
                    return 2;
                } else {
                    return 2;
                }

            } else {
              return 2.2;
            }
        }
        if (!energy.CompareTo("XeXe_5.44TeV")){
            if (!generator.CompareTo("")){
                if (!centrality.CompareTo("0-80%")){
                    return 2;
                } else if (!centrality.CompareTo("0-20%")){
                    return 2;
                } else if (!centrality.CompareTo("20-40%")){
                    return 2;
                } else if (!centrality.CompareTo("0-40%")){
                    return 2;
                } else return 2;
            }
        }
        return 0;
    }

    Color_t GetDefaultColorDiffDetectors( TString detector,
                                        Bool_t isMC,
                                        Bool_t isBox = kFALSE,
                                        Bool_t isHighPt = kTRUE){
        if (isMC){ // MC
            if (detector.CompareTo("PCM") == 0){
                return kGray+1;
            } else if (detector.CompareTo("PHOS") == 0){
                return kRed-7;
            } else if (detector.CompareTo("EMCal") == 0 || detector.CompareTo("EMC") == 0){
                return kGreen-6;
            } else if (detector.CompareTo("PCM-EMCal") == 0 || detector.CompareTo("PCM-EMC") == 0){
                return kBlue-6;
            } else if (detector.CompareTo("PCM-PHOS") == 0){
                return kOrange+1;
            } else if (detector.CompareTo("Dalitz") == 0 || detector.CompareTo("PCM-Dalitz") == 0){
                return kViolet-4;
            } else if (detector.CompareTo("EMCal high pT") == 0){
                return kOrange+1;
            } else if (detector.CompareTo("EMCal merged") == 0 || detector.CompareTo("EMC merged") == 0 || detector.CompareTo("mEMC") == 0){
                return kCyan-6;
            }	else {
                return kCyan-6;
            }

        } else { // data
            if (!isBox){
                if (detector.CompareTo("PCM") == 0){
                    if(isHighPt) return kBlack;
                    else  return kGray+1;
                } else if (detector.CompareTo("PHOS") == 0){
                    if(isHighPt) return kRed+1;
                    else return kRed-6;
                } else if (detector.CompareTo("EMCal") == 0 || detector.CompareTo("EMC") == 0){
                    if(isHighPt) return kGreen+2;
                    else return kGreen-6;
                } else if (detector.CompareTo("PCM-EMCal") == 0 || detector.CompareTo("PCM-EMC") == 0){
                    if(isHighPt) return kBlue+1;
                    else return kBlue-6;
                } else if (detector.CompareTo("PCM-PHOS") == 0){
                    if(isHighPt) return 807;
                        else return 806;
                } else if (detector.CompareTo("Dalitz") == 0 || detector.CompareTo("PCM-Dalitz") == 0){
                    if(isHighPt) return kViolet;
                        else return kViolet-4;
                } else if (detector.CompareTo("Comb") == 0){
                    if(isHighPt) return kBlack;
                    else  return kGray+1;
                }  else if (detector.CompareTo("EMCal high pT") == 0){
                    if(isHighPt) return 807;
                        else return 806;
                } else if (detector.CompareTo("EMCal merged") == 0 || detector.CompareTo("EMC merged") == 0 || detector.CompareTo("mEMC") == 0){
                    if(isHighPt) return kCyan+2;
                        else return kCyan-6;
                } else {
                    if(isHighPt) return kCyan+2;
                        else return kCyan-6;
                }
            } else {
                if (detector.CompareTo("PCM") == 0){
                    return kGray+1;
                } else if (detector.CompareTo("PHOS") == 0){
                    return kRed-6;
                } else if (detector.CompareTo("EMCal") == 0 || detector.CompareTo("EMC") == 0){
                    return kGreen-6;
                } else if (detector.CompareTo("PCM-EMCal") == 0 || detector.CompareTo("PCM-EMC") == 0){
                    return kBlue-7;
                } else if (detector.CompareTo("PCM-PHOS") == 0){
                    return 806;
                } else if (detector.CompareTo("Dalitz") == 0 || detector.CompareTo("PCM-Dalitz") == 0){
                    return kViolet-4;
                }  else if (detector.CompareTo("EMCal high pT") == 0){
                    return 806;
                } else if (detector.CompareTo("EMCal merged") == 0 || detector.CompareTo("EMC merged") == 0 || detector.CompareTo("mEMC") == 0){
                    return kCyan-6;
                } else {
                    return kCyan-6;
                }
            }
        }
        return kBlack;
    }



    Style_t GetDefaultMarkerStyleDiffDetectors( TString detector,
                                                Bool_t isMC){
        if (isMC){ // MC
            if (detector.CompareTo("PCM") == 0){
                return 24;
            } else if (detector.CompareTo("PHOS") == 0){
                return 25;
            } else if (detector.CompareTo("EMCal") == 0 || detector.CompareTo("EMC") == 0){
                return 27;
            } else if (detector.CompareTo("PCM-EMCal") == 0 || detector.CompareTo("PCM-EMC") == 0){
                return 28;
            } else if (detector.CompareTo("PCM-PHOS") == 0){
                return 28;
            } else if (detector.CompareTo("Dalitz") == 0 || detector.CompareTo("PCM-Dalitz") == 0){
                return 30;
            }  else if (detector.CompareTo("EMCal high pT") == 0){
                return 28;
            } else if (detector.CompareTo("EMCal merged") == 0 || detector.CompareTo("EMC merged") == 0 || detector.CompareTo("mEMC") == 0){
                return 30;
            } else if (detector.CompareTo("Comb") == 0){
                return 20;
            } else {
                return 30;
            }
        } else { // data
            if (detector.CompareTo("PCM") == 0){
                return 20;
            } else if (detector.CompareTo("PHOS") == 0){
                return 21;
            } else if (detector.CompareTo("EMCal") == 0 || detector.CompareTo("EMC") == 0){
                return 33;
            } else if (detector.CompareTo("PCM-EMCal") == 0 || detector.CompareTo("PCM-EMC") == 0){
                return 34;
            } else if (detector.CompareTo("PCM-PHOS") == 0){
                return 34;
            } else if (detector.CompareTo("Dalitz") == 0 || detector.CompareTo("PCM-Dalitz") == 0){
                return 29;
            }  else if (detector.CompareTo("EMCal high pT") == 0){
                return 34;
            } else if (detector.CompareTo("EMCal merged") == 0 || detector.CompareTo("EMC merged") == 0 || detector.CompareTo("mEMC") == 0){
                return 29;
            } else if (detector.CompareTo("Comb") == 0){
                return 20;
            } else {
                return 29;
            }
        }
        return 1;
    }

    Size_t GetDefaultMarkerSizeDiffDetectors( TString detector,
                                            Bool_t isMC){
        if (isMC){ // MC
            if (detector.CompareTo("PCM") == 0){
                return 2.2;
            } else if (detector.CompareTo("PHOS") == 0){
                return 2.2;
            } else if (detector.CompareTo("EMCal") == 0 || detector.CompareTo("EMC") == 0){
                return 3.3;
            } else if (detector.CompareTo("PCM-EMCal") == 0 || detector.CompareTo("PCM-EMC") == 0){
                return 2.2;
            } else if (detector.CompareTo("PCM-PHOS") == 0){
                return 2.2;
            } else if (detector.CompareTo("Dalitz") == 0 || detector.CompareTo("PCM-Dalitz") == 0){
                return 2.2;
            }  else if (detector.CompareTo("EMCal high pT") == 0){
                return 2.2;
            } else if (detector.CompareTo("EMCal merged") == 0 || detector.CompareTo("EMC merged") == 0 || detector.CompareTo("mEMC") == 0){
                return 2.6;
            } else {
                return 2.2;
            }
        } else { // data
            if (detector.CompareTo("PCM") == 0){
                return 2.2;
            } else if (detector.CompareTo("PHOS") == 0){
                return 2.2;
            } else if (detector.CompareTo("EMCal") == 0 || detector.CompareTo("EMC") == 0){
                return 3.3;
            } else if (detector.CompareTo("PCM-EMCal") == 0 || detector.CompareTo("PCM-EMC") == 0){
                return 2.2;
            } else if (detector.CompareTo("PCM-PHOS") == 0){
                return 2.2;
            } else if (detector.CompareTo("Dalitz") == 0 || detector.CompareTo("PCM-Dalitz") == 0){
                return 2.2;
            }  else if (detector.CompareTo("EMCal high pT") == 0){
                return 2.2;
            } else if (detector.CompareTo("EMCal merged") == 0 || detector.CompareTo("EMC merged") == 0 || detector.CompareTo("mEMC") == 0){
                return 2.6;
            } else {
                return 2.2;
            }
        }
        return 1.;
    }

    void PlotErrorBarAtUpperEdgeOfTGraphAsymErr(TGraphAsymmErrors* graph,
                                                Double_t widthTick = 0.05,
                                                Bool_t isLog = kFALSE,
                                                TGraphAsymmErrors* graphpTbinning = NULL){

        TGraphAsymmErrors* dummy = (TGraphAsymmErrors*)graph->Clone("dummyPlotErrorsAtEdge");
        if(!graphpTbinning){
            for (Int_t i=0; i < dummy->GetN();i++){
                Double_t widthTickUsedUp = widthTick;
                Double_t widthTickUsedDown = widthTick;
                DrawGammaLines(dummy->GetX()[i]-widthTickUsedUp, dummy->GetX()[i]+widthTickUsedDown, dummy->GetY()[i]+dummy->GetEYhigh()[i], dummy->GetY()[i]+dummy->GetEYhigh()[i],
                            dummy->GetLineWidth(), dummy->GetLineColor(),dummy->GetLineStyle());

            }
        } else {
            for (Int_t i=0; i < dummy->GetN();i++){
                Double_t widthTickUsedUp = graphpTbinning->GetEXhigh()[i];
                Double_t widthTickUsedDown = graphpTbinning->GetEXlow()[i];
                DrawGammaLines(dummy->GetX()[i]-widthTickUsedUp, dummy->GetX()[i]+widthTickUsedDown, dummy->GetY()[i]+dummy->GetEYhigh()[i], dummy->GetY()[i]+dummy->GetEYhigh()[i],
                            dummy->GetLineWidth(), dummy->GetLineColor(),dummy->GetLineStyle());

            }
        }
    }
    void PlotErrorBarAtUpperEdgeOfTGraphErr(TGraphErrors* graph,
                                            Double_t widthTick = 0.05,
                                            Bool_t isLog = kFALSE,
                                            TGraphErrors* graphpTbinning = NULL){

        TGraphErrors* dummy = (TGraphErrors*)graph->Clone("dummyPlotErrorsAtEdge");
        if(!graphpTbinning){
            for (Int_t i=0; i < dummy->GetN();i++){
                Double_t widthTickUsedUp = widthTick;
                Double_t widthTickUsedDown = widthTick;
                DrawGammaLines(dummy->GetX()[i]-widthTickUsedUp, dummy->GetX()[i]+widthTickUsedDown, dummy->GetY()[i]+dummy->GetEY()[i], dummy->GetY()[i]+dummy->GetEY()[i],
                               dummy->GetLineWidth(), dummy->GetLineColor(),dummy->GetLineStyle());

            }
        } else {
            for (Int_t i=0; i < dummy->GetN();i++){
                Double_t widthTickUsedUp = graphpTbinning->GetEX()[i];
                Double_t widthTickUsedDown = graphpTbinning->GetEX()[i];
                DrawGammaLines(dummy->GetX()[i]-widthTickUsedUp, dummy->GetX()[i]+widthTickUsedDown, dummy->GetY()[i]+dummy->GetEY()[i], dummy->GetY()[i]+dummy->GetEY()[i],
                               dummy->GetLineWidth(), dummy->GetLineColor(),dummy->GetLineStyle());

            }
        }
    }


    Color_t GetDefaultTriggerColor( TString periodName,
                                    Int_t SpecialTrigger){
        if( periodName.CompareTo("LHC11a") == 0                 ||
            periodName.CompareTo("LHC11a_p4_wSDD") == 0         ||
            periodName.CompareTo("LHC11a_p4_wSDD-kEMC1") == 0   ||
            periodName.CompareTo("LHC11a_pass4_wSDD") == 0
        ){
            switch(SpecialTrigger)
            {
            case 3:
                return kBlack;
            case 51:
                return kAzure;
            default:
                cout << "GetDefaultTriggerMarker(): no valid input parameters '" << SpecialTrigger << "' given for period: " << periodName.Data() << endl;
                return kBlack;
            }
        }else if(periodName.CompareTo("LHC13g") == 0            || periodName.CompareTo("LHC13g-kEMC7") == 0        ||
                periodName.CompareTo("LHC13g-kEMCEG1") == 0    || periodName.CompareTo("LHC13g-kEMCEG2") == 0      ||
                periodName.CompareTo("LHC13g_pass1") == 0 )	{
            switch(SpecialTrigger)
            {
            case 0:
                return kBlack;
            case 10:
                return kBlack;
            case 52:
                return kAzure;
            case 83:
                return kGreen+2;
            case 85:
                return kOrange+2;
            case 93:
                return kRed;
            case 95:
                return kViolet;
            default:
                cout << "GetDefaultTriggerMarker(): no valid input parameters '" << SpecialTrigger << "' given for period: " << periodName.Data() << endl;
                return kBlack;
            }
        }else if(periodName.CompareTo("LHC12")==0 || (periodName.BeginsWith("LHC12")&&periodName.Length()==6)
                || (periodName.BeginsWith("LHC12")&&periodName.EndsWith("-kEMC7"))
                || (periodName.BeginsWith("LHC12")&&periodName.EndsWith("-kEMCEGA"))
                || (periodName.BeginsWith("LHC12")&&periodName.EndsWith("-kEMC8"))
                || (periodName.BeginsWith("LHC12")&&periodName.EndsWith("-kEMC8EGA"))
                || (periodName.BeginsWith("LHC12")&&periodName.EndsWith("-kEMCEJE")))	{
            switch(SpecialTrigger)
            {
            case 0:
                return kBlack;
            case 11:
                return kBlack;
            case 52:
                return kAzure;
            case 53:
                return kAzure;
            case 81:
                return kGreen+2;
            case 82:
                return kGreen+2;
            case 91:
                return kRed;
            default:
                cout << "GetDefaultTriggerMarker(): no valid input parameters '" << SpecialTrigger << "' given for period: " << periodName.Data() << endl;
                return kBlack;
            }
        }

        cout << "GetDefaultTriggerColor(): no valid input parameters given..." << endl;
        return kBlack;
    }

    Style_t GetDefaultTriggerMarker( TString periodName,
                                    Int_t SpecialTrigger){
        if( periodName.CompareTo("LHC11a") == 0                 ||
            periodName.CompareTo("LHC11a_p4_wSDD") == 0         ||
            periodName.CompareTo("LHC11a_p4_wSDD-kEMC1") == 0   ||
            periodName.CompareTo("LHC11a_pass4_wSDD") == 0
        ){
            switch(SpecialTrigger)
            {
            case 3:
                return 20;
            case 51:
                return 21;
            default:
                cout << "GetDefaultTriggerMarker(): no valid input parameters '" << SpecialTrigger << "' given for period: " << periodName.Data() << endl;
                return 29;
            }
        }else if(   periodName.CompareTo("LHC13g") == 0         || periodName.CompareTo("LHC13g-kEMC7") == 0    ||
                    periodName.CompareTo("LHC13g-kEMCEG1") == 0 || periodName.CompareTo("LHC13g-kEMCEG2") == 0  ||
                    periodName.CompareTo("LHC13g_pass1") == 0
                ){
            switch(SpecialTrigger)
            {
            case 0:
                return 20;
            case 10:
                return 20;
            case 52:
                return 21;
            case 83:
                return 22;
            case 85:
                return 23;
            case 93:
                return 24;
            case 95:
                return 25;
            default:
                cout << "GetDefaultTriggerMarker(): no valid input parameters '" << SpecialTrigger << "' given for period: " << periodName.Data() << endl;
                return 29;
            }
        }else if(periodName.CompareTo("LHC12")==0 || (periodName.BeginsWith("LHC12")&&periodName.Length()==6)
                    || (periodName.BeginsWith("LHC12")&&periodName.EndsWith("-kEMC7"))
                    || (periodName.BeginsWith("LHC12")&&periodName.EndsWith("-kEMCEGA"))
                    || (periodName.BeginsWith("LHC12")&&periodName.EndsWith("-kEMC8"))
                    || (periodName.BeginsWith("LHC12")&&periodName.EndsWith("-kEMC8EGA"))
                    || (periodName.BeginsWith("LHC12")&&periodName.EndsWith("-kEMCEJE"))){
            switch(SpecialTrigger)
            {
            case 0:
                return 20;
            case 10:
                return 20;
            case 52:
                return 21;
            case 53:
                return 21;
            case 81:
                return 22;
            case 82:
                return 22;
            case 91:
                return 24;
            default:
                cout << "GetDefaultTriggerMarker(): no valid input parameters '" << SpecialTrigger << "' given for period: " << periodName.Data() << endl;
                return 29;
            }
        }

        cout << "GetDefaultTriggerMarker(): no valid input parameters '" << SpecialTrigger << "' given for period: " << periodName.Data() << endl;
        return 29;
    }

    Color_t GetDefaultTriggerColorName (TString triggerName, Bool_t isShade, TString opt = ""){
        if ( (triggerName.Contains("MB") || triggerName.Contains("INT1") || (opt.BeginsWith("8TeV") && triggerName.Contains("INT7"))) && !isShade)     return kBlack;
        else if ((triggerName.Contains("MB") || triggerName.Contains("INT1") || (opt.BeginsWith("8TeV") && triggerName.Contains("INT7"))) && isShade)  return kGray+1;
        else if (triggerName.Contains("INT7") && !isShade) return kGray+1;
        else if (triggerName.Contains("INT7") && isShade)  return kGray;
        else if (triggerName.Contains("EMC1") && !isShade) return kRed+2;
        else if (triggerName.Contains("EMC1") && isShade)  return kRed-6;
        else if ((triggerName.Contains("EMC7")|| triggerName.Contains("PHI7")) && !isShade) return kBlue+2;
        else if ((triggerName.Contains("EMC7")|| triggerName.Contains("PHI7")) && isShade)  return kBlue-6;
        else if ((triggerName.Contains("EG2") || triggerName.Contains("EGA")) && !isShade)  return kGreen+3;
        else if ((triggerName.Contains("EG2") || triggerName.Contains("EGA")) && isShade)   return kGreen-8;
        else if (triggerName.Contains("EG1") && !isShade)  return kCyan+2;
        else if (triggerName.Contains("EG1") && isShade)   return kCyan-6;
        else if ((triggerName.Contains("EJ2") ) && !isShade)  return kRed+2;
        else if ((triggerName.Contains("EJ2") ) && isShade)   return kRed-6;
        else if (triggerName.Contains("EJ1") && !isShade)  return kRed+2;
        else if (triggerName.Contains("EJ1") && isShade)   return kRed-6;
        else if (triggerName.Contains("PHI7") && !isShade) return kRed+2;
        else if (triggerName.Contains("PHI7") && isShade) return kRed-6;
        else return kBlack;

        return kBlack;
    }

    Marker_t GetDefaultTriggerMarkerStyleName (TString triggerName, Bool_t isShade ){
        if ((triggerName.Contains("MB_NLM1") || triggerName.Contains("INT1_NLM1") )&& !isShade)  return 33;
        else if ((triggerName.Contains("MB") || triggerName.Contains("INT1") )&& !isShade)       return 20;
        else if ((triggerName.Contains("MB_NLM1") || triggerName.Contains("INT1_NLM1") )&& isShade)        return 27;
        else if ((triggerName.Contains("MB") || triggerName.Contains("INT1") )&& isShade)        return 24;
        else if ((triggerName.Contains("INT7_NLM1") || triggerName.Contains("INT7_CF")) && !isShade) return 27;
        else if (triggerName.Contains("INT7") && !isShade)      return 20;
        else if ((triggerName.Contains("INT7_NLM1") || triggerName.Contains("INT7_CF")) && isShade)  return 33;
        else if (triggerName.Contains("INT7") && isShade)       return 24;
        else if (triggerName.Contains("EMC1_NLM1") && !isShade) return 25;
        else if (triggerName.Contains("EMC1") && !isShade) return 21;
        else if (triggerName.Contains("EMC1") && isShade)  return 25;
        else if (triggerName.Contains("EMC7_NLM1") && !isShade) return 28;
        else if ((triggerName.Contains("EMC7") || triggerName.Contains("PHI7")) && !isShade) return 34;
        else if ((triggerName.Contains("EMC7") || triggerName.Contains("PHI7")) && isShade)  return 28;
        else if (triggerName.Contains("EG2_NLM1") && !isShade)  return 30;
        else if ((triggerName.Contains("EG2") || triggerName.Contains("EGA")) && !isShade)  return 29;
        else if ((triggerName.Contains("EG2") || triggerName.Contains("EGA")) && isShade)   return 30;
        else if (triggerName.Contains("EG1_NLM1") && !isShade)  return 27;
        else if (triggerName.Contains("EG1") && !isShade)  return 33;
        else if (triggerName.Contains("EG1") && isShade)   return 27;
        else if (triggerName.Contains("EJ2") && !isShade) return 21;
        else if (triggerName.Contains("EJ2") && isShade)  return 25;
        else if (triggerName.Contains("EJ1") && !isShade) return 33;
        else if (triggerName.Contains("EJ1") && isShade)  return 24;
        else if (triggerName.Contains("PHI7") && !isShade) return 21;
        else if (triggerName.Contains("PHI7") && isShade) return 25;
        else return 20;

        return 20;
    }

    Size_t GetDefaultTriggerMarkerSizeName (TString triggerName, Bool_t isShade ){
        if ((triggerName.Contains("MB_NLM1") || triggerName.Contains("INT1_NLM1") )&& !isShade)  return (Size_t)2;
        else if ((triggerName.Contains("MB") || triggerName.Contains("INT1") )&& !isShade)       return (Size_t)1.5;
        else if ((triggerName.Contains("MB_NLM1") || triggerName.Contains("INT1_NLM1") )&& isShade) return (Size_t)2;
        else if ((triggerName.Contains("MB") || triggerName.Contains("INT1") )&& isShade)        return (Size_t)1.5;
        else if ((triggerName.Contains("INT7_NLM1") || triggerName.Contains("INT7_CF"))&& !isShade) return (Size_t)2;
        else if (triggerName.Contains("INT7") && !isShade)      return (Size_t)1.5;
        else if ((triggerName.Contains("INT7_NLM1") || triggerName.Contains("INT7_CF")) && isShade)  return (Size_t)2;
        else if (triggerName.Contains("INT7") && isShade)       return (Size_t)1.5;
        else if (triggerName.Contains("EMC1_NLM1") && !isShade) return (Size_t)1.5;
        else if (triggerName.Contains("EMC1") && !isShade) return (Size_t)1.5;
        else if (triggerName.Contains("EMC1") && isShade)  return (Size_t)1.5;
        else if (triggerName.Contains("EMC7_NLM1") && !isShade) return (Size_t)2.;
        else if (triggerName.Contains("EMC7") && !isShade) return (Size_t)2.;
        else if (triggerName.Contains("EMC7") && isShade)  return (Size_t)2.;
        else if (triggerName.Contains("EG2_NLM1") && !isShade)  return (Size_t)2;
        else if ((triggerName.Contains("EG2") || triggerName.Contains("EGA")) && !isShade)  return (Size_t)2;
        else if ((triggerName.Contains("EG2") || triggerName.Contains("EGA")) && isShade)   return (Size_t)2;
        else if (triggerName.Contains("EG1_NLM1") && !isShade)  return (Size_t)2;
        else if (triggerName.Contains("EG1") && !isShade)  return (Size_t)2;
        else if (triggerName.Contains("EG1") && isShade)   return (Size_t)2;
        else if (triggerName.Contains("EJ2") && !isShade)  return (Size_t)2;
        else if (triggerName.Contains("EJ2") && isShade)   return (Size_t)2;
        else if (triggerName.Contains("EJ1") && !isShade)  return (Size_t)2;
        else if (triggerName.Contains("EJ1") && isShade)   return (Size_t)2;
        else if (triggerName.Contains("PHI7") && !isShade) return (Size_t)2;
        else if (triggerName.Contains("PHI7") && isShade) return (Size_t)2;
        else return (Size_t)1.5;

        return (Size_t)1.5;
    }




    void DrawMergedClusterLambdaCuts (Int_t nlm = 1, Float_t maxPt = 50.05){
        if (nlm == 1 || nlm == 0 ){
            TF1 *min1 = new TF1("min1","TMath::Exp(2.135-0.245*x)",4.95,14.63);
            min1->SetLineColor(kWhite);
            min1->SetLineWidth(4);
            min1->SetLineStyle(2);
            min1->Draw("same");
            TF1 *min = new TF1("min","TMath::Exp(2.135-0.245*x)",4.95,14.63);
            min->SetLineColor(kRed+2);
            min->SetLineWidth(3);
            min->SetLineStyle(2);
            min->Draw("same");
            TF1 *cnst1 = new TF1("cnst","0.27",13.63,maxPt);
            cnst1->SetLineColor(kWhite);
            cnst1->SetLineWidth(4);
            cnst1->SetLineStyle(2);
            cnst1->Draw("same");
            TF1 *cnst = new TF1("cnst","0.27",13.63,maxPt);
            cnst->SetLineColor(kRed+2);
            cnst->SetLineStyle(2);
            cnst->SetLineWidth(3);
            cnst->Draw("same");
            TF1 *max1 = new TF1("max","TMath::Exp(0.0662-0.0201*x)-0.0955+0.00186*x+21.9/x",9.5,maxPt);
            max1->SetLineColor(kWhite);
            max1->SetLineWidth(4);
            max1->SetLineStyle(2);
            max1->Draw("same");
            TF1 *max = new TF1("max","TMath::Exp(0.0662-0.0201*x)-0.0955+0.00186*x+21.9/x",9.5,maxPt);
            max->SetLineColor(kRed+2);
            max->SetLineWidth(3);
            max->SetLineStyle(2);
            max->Draw("same");
        } else if (nlm == 2){
            TF1 *min = new TF1("min","TMath::Exp(2.135-0.245*x)",4.95,14.63);
            min->SetLineColor(kWhite);
            min->SetLineWidth(2);
            min->SetLineStyle(2);
            min->Draw("same");
            min->SetLineColor(kBlack);
            min->SetLineWidth(1);
            min->SetLineStyle(2);
            min->Draw("same");
            TF1 *cnst = new TF1("cnst","0.27",13.63,maxPt);
            cnst->SetLineColor(kWhite);
            cnst->SetLineWidth(2);
            cnst->SetLineStyle(2);
            cnst->Draw("same");
            cnst->SetLineColor(kBlack);
            cnst->SetLineStyle(2);
            cnst->SetLineWidth(1);
            cnst->Draw("same");
            TF1 *max = new TF1("max","TMath::Exp(0.353-0.0264*x)-0.524+0.00559*x+21.9/x",9.5,maxPt);
            max->SetLineColor(kWhite);
            max->SetLineWidth(2);
            max->SetLineStyle(2);
            max->Draw("same");
            max->SetLineColor(kBlack);
            max->SetLineWidth(1);
            max->SetLineStyle(2);
            max->Draw("same");
        }
    }

    //**********************************************************************************
    //**********************************************************************************
    //**********************************************************************************
    Color_t GetColorSystematics(    TString variationName,
                                    TString mode            = ""
                                ){

        if ( variationName.Contains("BGEstimateIterations") || variationName.CompareTo("PileupDCA") == 0 || variationName.CompareTo("OOBPileupPi0") == 0)
            return kTeal+3;
        if ( variationName.Contains("BGEstimate")  || variationName.CompareTo("Pileup") == 0 || variationName.CompareTo("pileup") == 0 || variationName.CompareTo("OOBPileupGamma") == 0 )
            return kTeal-5;
        if ( variationName.CompareTo("BG") == 0 )
            return 800;
        if ( variationName.Contains("Periods") )
            return kBlue-7;
        if ( variationName.CompareTo("SPD") == 0 )
            return kCyan+4;
        if ( variationName.CompareTo("DistanceBadChannel") == 0 )
            return kGreen-1;
        if ( variationName.CompareTo("Purity") == 0 )
            return kBlue+2;
        if ( variationName.CompareTo("CosPoint") == 0 )
            return kSpring+2;
        if ( variationName.CompareTo("Eta") == 0 )
            return kTeal-7;
        if ( variationName.CompareTo("MCSmearing") == 0 )
            return kAzure-5;
        if ( variationName.CompareTo("PsiPair") == 0 )
            return kRed-2;
        if ( variationName.CompareTo("YieldExtraction") == 0 || variationName.CompareTo("YieldExtraction_pp") == 0 || variationName.CompareTo("IntRange") == 0 )
            return kBlue;
        if ( variationName.CompareTo("dEdxE") == 0 )
            return kRed+1;
        if ( variationName.CompareTo("dEdxPi") == 0 )
            return kOrange+7;
        if ( variationName.CompareTo("TPCCluster") == 0 )
            return kPink+8;
        if ( variationName.CompareTo("SinglePt") == 0 )
            return kGreen+2;
        if ( variationName.CompareTo("Chi2") == 0 || variationName.CompareTo("Chi2PsiPair") == 0 || variationName.CompareTo("PsiPairChi2") == 0 )
            return kYellow+2;
        if ( variationName.CompareTo("Qt") == 0 )
            return kOrange+2;
        if ( variationName.CompareTo("Alpha") == 0 || variationName.CompareTo("Asym") == 0)
            return kBlue+2;
        if ( variationName.CompareTo("Rapidity") == 0 )
            return kPink-3;
        if ( variationName.CompareTo("ConvPhi") == 0 )
            return kCyan-2;
        if ( variationName.CompareTo("DoubleCount") == 0)
            return kGray+2;
        if ( variationName.CompareTo("ClusterMinEnergy") == 0 )
            return kViolet+1;
        if ( variationName.CompareTo("ClusterNCells") == 0 )
            return kAzure+1;
        if ( variationName.CompareTo("ClusterNonLinearity") == 0 )
            return kAzure+4;
        if ( variationName.CompareTo("ClusterEnergyScale") == 0 )
            return kPink+4;
        if ( variationName.CompareTo("ClusterTrackMatching") == 0 || variationName.CompareTo("ClusterTrackMatchingCalo") == 0 )
            return kOrange;
        if ( variationName.CompareTo("ClusterM02") == 0 )
            return kGreen-6;
        if ( variationName.CompareTo("CellTiming") == 0 )
            return kTeal-5;
        if ( variationName.CompareTo("ClusterMaterialTRD") == 0 )
            return kOrange+5;
        if ( variationName.CompareTo("Trigger") == 0 || variationName.CompareTo("PileUp") == 0 )
            return 404;
        if ( variationName.CompareTo("Efficiency") == 0 )
            return kPink-6;
        if ( variationName.CompareTo("YieldExtractionPi0") == 0 || variationName.CompareTo("IntRange") == 0 )
            return 860;
        if ( variationName.CompareTo("ClusterTime") == 0 || variationName.CompareTo("TOF") == 0 || variationName.CompareTo("ClusterTiming") == 0 )
            return kBlue-4;
        if ( variationName.CompareTo("ClusterizationEnergy") == 0 )
            return kRed-4;
        if ( variationName.CompareTo("Periods") == 0 )
            return kGreen-4;
        if ( variationName.CompareTo("InnerMaterial") == 0 || variationName.CompareTo("Material") == 0)
            return kMagenta-2;
        if ( variationName.CompareTo("OpeningAngle") == 0 )
            return kCyan-3;
        if ( variationName.CompareTo("MesonResolution") == 0 )
            return kPink-3;
        if ( variationName.CompareTo("Secondary") == 0 )
            return kViolet-3;
        if ( variationName.CompareTo("CellMinE") == 0 )
            return kSpring-6;
        if ( variationName.CompareTo("Cocktail") == 0 )
            return kRed-6;

        cout << "GetColorSystematics(): no valid input parameters given for '" << variationName.Data() <<  "'..." << endl;
        return kBlack;
    }

    //**********************************************************************************
    //************** Determine marker style for systematics ****************************
    //**********************************************************************************
    Style_t GetMarkerStyleSystematics(  TString variationName,
                                        TString mode            = ""
                                    ){

        if ( variationName.Contains("BGEstimateIterations") || variationName.CompareTo("PileupDCA") == 0 || variationName.CompareTo("OOBPileupPi0") == 0)
            return 25;
        if ( variationName.Contains("BGEstimate")  || variationName.CompareTo("Pileup") == 0 || variationName.CompareTo("pileup") == 0 || variationName.CompareTo("OOBPileupGamma") == 0  )
            return 24;
        if ( variationName.CompareTo("BG") == 0 )
            return 28;
        if ( variationName.Contains("Periods") )
            return 33;
        if ( variationName.CompareTo("SPD") == 0 )
            return 29;
        if ( variationName.CompareTo("DistanceBadChannel") == 0 )
            return 29;
        if ( variationName.CompareTo("Purity") == 0 )
            return 27;
        if ( variationName.CompareTo("CosPoint") == 0 )
            return 24;
        if ( variationName.CompareTo("Eta") == 0 )
            return 30;
        if ( variationName.CompareTo("MCSmearing") == 0 )
            return 21;
        if ( variationName.CompareTo("PsiPair") == 0 )
            return 23;
        if ( variationName.CompareTo("YieldExtraction") == 0 || variationName.CompareTo("YieldExtraction_pp") == 0 )
            return 24;
        if ( variationName.CompareTo("dEdxE") == 0 )
            return 21;
        if ( variationName.CompareTo("dEdxPi") == 0 )
            return 22;
        if ( variationName.CompareTo("TPCCluster") == 0 )
            return 23;
        if ( variationName.CompareTo("SinglePt") == 0 )
            return 20;
        if ( variationName.CompareTo("Chi2") == 0 || variationName.CompareTo("Chi2PsiPair") == 0 || variationName.CompareTo("PsiPairChi2") == 0 )
            return 25;
        if ( variationName.CompareTo("Qt") == 0 )
            return 26;
        if ( variationName.CompareTo("Alpha") == 0 ||  variationName.CompareTo("Asym") == 0)
            return 27;
        if ( variationName.CompareTo("Rapidity") == 0 )
            return 26;
        if ( variationName.CompareTo("ConvPhi") == 0 )
            return 28;
        if ( variationName.CompareTo("DoubleCount") == 0)
            return 29;
        if ( variationName.CompareTo("ClusterMinEnergy") == 0 )
            return 29;
        if ( variationName.CompareTo("ClusterNCells") == 0 )
            return 30;
        if ( variationName.CompareTo("ClusterNonLinearity") == 0 )
            return 31;
        if ( variationName.CompareTo("ClusterEnergyScale") == 0 )
            return 29;
        if ( variationName.CompareTo("ClusterTrackMatching") == 0 || variationName.CompareTo("ClusterTrackMatchingCalo") == 0 )
            return 32;
        if ( variationName.CompareTo("ClusterM02") == 0 )
            return 33;
        if ( variationName.CompareTo("CellTiming") == 0 || variationName.CompareTo("TOF") == 0)
            return 24;
        if ( variationName.CompareTo("ClusterMaterialTRD") == 0 )
            return 21;
        if ( variationName.CompareTo("Trigger") == 0  || variationName.CompareTo("PileUp") == 0)
            return 22;
        if ( variationName.CompareTo("Efficiency") == 0 )
            return 23;
        if ( variationName.CompareTo("YieldExtractionPi0") == 0 || variationName.CompareTo("IntRange") == 0 )
            return 20;
        if ( variationName.CompareTo("ClusterTime") == 0 || variationName.CompareTo("ClusterTiming") == 0)
            return 25;
        if ( variationName.CompareTo("ClusterizationEnergy") == 0 )
            return 24;
        if ( variationName.CompareTo("Periods") == 0 )
            return 21;
        if ( variationName.CompareTo("InnerMaterial") == 0 || variationName.CompareTo("Material") == 0 )
            return 20;
        if ( variationName.CompareTo("OpeningAngle") == 0 )
            return 25;
        if ( variationName.CompareTo("MesonResolution") == 0 )
            return 24;
        if ( variationName.CompareTo("Secondary") == 0 )
            return 30;
        if ( variationName.CompareTo("CellMinE") == 0 )
            return 33;
        if ( variationName.CompareTo("Cocktail") == 0 )
            return 24;

        cout << "GetMarkerStyleSystematics(): no valid input parameters given for '" << variationName.Data() <<  "'..." << endl;
        return 20;
    }

    //**********************************************************************************
    //********************** Determine name for systematics ****************************
    //**********************************************************************************
    TString GetSystematicsName  (  TString variationName = "" ){

        if ( variationName.CompareTo("PileupDCA") == 0 || variationName.Contains("BGEstimateIterations") )
            return "OB pileup B";
        if ( variationName.CompareTo("OOBPileupGamma") == 0)
            return "OB pileup #gamma";
        if ( variationName.CompareTo("OOBPileupPi0") == 0)
            return "OB pileup #pi^{0}";
        if ( variationName.Contains("BGEstimate") )
            return "OB pileup A";
        if ( variationName.CompareTo("Pileup") == 0 || variationName.CompareTo("pileup") == 0 )
            return "IB pileup ";
        if ( variationName.CompareTo("BG") == 0 )
            return "BG method";
        if ( variationName.CompareTo("Purity") == 0 )
            return "#epsilon_{pur} correction";
        if ( variationName.Contains("Periods"))
            return "periods";
        if ( variationName.CompareTo("SPD") == 0 )
            return "IB pileup";
        if ( variationName.CompareTo("DistanceBadChannel") == 0 )
            return "Dist. BC";
        if ( variationName.CompareTo("CosPoint") == 0 )
            return "cos(#Theta_{point})";
        if ( variationName.CompareTo("Eta") == 0 )
            return "#eta_{#gamma, e^{#pm}}";
        if ( variationName.CompareTo("MCSmearing") == 0 )
            return "MC smearing";
        if ( variationName.CompareTo("PsiPair") == 0 )
            return "#psi_{pair} #gamma";
        if ( variationName.CompareTo("YieldExtraction") == 0 || variationName.CompareTo("YieldExtraction_pp") == 0 )
            return "yield extraction";
        if ( variationName.CompareTo("dEdxE") == 0 )
            return "dE/dx e-line";
        if ( variationName.CompareTo("dEdxPi") == 0 )
            return "dE/dx #pi-line";
        if ( variationName.CompareTo("TPCCluster") == 0 )
            return "TPC cluster";
        if ( variationName.CompareTo("SinglePt") == 0 )
            return "single e^{#pm} #it{p}_{T}";
        if ( variationName.CompareTo("Chi2") == 0 || variationName.CompareTo("Chi2PsiPair") == 0 || variationName.CompareTo("PsiPairChi2") == 0 )
            return "2D #chi^{2} #gamma, #psi_{pair} #gamma";
        if ( variationName.CompareTo("Qt") == 0 )
            return "2D q_{T}";
        if ( variationName.CompareTo("Alpha") == 0 )
            return "#alpha meson";
        if ( variationName.CompareTo("Asym") == 0 )
            return "asymmetry";
        if ( variationName.CompareTo("Rapidity") == 0 )
            return "rapidity";
        if ( variationName.CompareTo("ConvPhi") == 0 )
            return "#varphi_{conv}";
        if ( variationName.CompareTo("TOF") == 0 )
            return "TOF";
        if ( variationName.CompareTo("DoubleCount") == 0)
            return "double count";
        if ( variationName.CompareTo("ClusterMinEnergy") == 0 )
            return "min E_{cluster}";
        if ( variationName.CompareTo("ClusterNCells") == 0 )
            return "min # cells";
        if ( variationName.CompareTo("ClusterNonLinearity") == 0 )
            return "clu. energy calibration";
        if ( variationName.CompareTo("ClusterEnergyScale") == 0 )
            return "clu. energy scale";
        if ( variationName.CompareTo("ClusterTrackMatching") == 0 )
            return "V^{0}-tr. match. to cl.";
        if ( variationName.CompareTo("ClusterTrackMatchingCalo") == 0 )
            return "tr. match. to cl.";
        if ( variationName.CompareTo("ClusterM02") == 0 )
            return "#sigma_{long}^{2}";
        if ( variationName.CompareTo("CellTiming") == 0 )
            return "cell time";
        if ( variationName.CompareTo("ClusterMaterialTRD") == 0 )
            return "mat. infront of EMCal";
        if ( variationName.CompareTo("Trigger") == 0 )
            return "trigger normalization";
        if ( variationName.CompareTo("PileUp") == 0 )
            return "SPD pileup";
        if ( variationName.CompareTo("Efficiency") == 0 )
            return "efficiency";
        if ( variationName.CompareTo("YieldExtractionPi0") == 0 || variationName.CompareTo("IntRange") == 0 )
            return "yield extraction #pi^{0}";
        if ( variationName.CompareTo("YieldExtractionEta") == 0 )
            return "yield extraction #eta";
        if ( variationName.CompareTo("ClusterTime") == 0 || variationName.CompareTo("ClusterTiming") == 0)
            return "clu. timing";
        if ( variationName.CompareTo("ClusterizationEnergy") == 0 )
            return "clusterization energy";
        if ( variationName.CompareTo("Periods") == 0 )
            return "periods";
        if ( variationName.CompareTo("InnerMaterial") == 0 )
            return "inner material";
        if ( variationName.CompareTo("Material") == 0 )
            return "material";
        if ( variationName.CompareTo("OpeningAngle") == 0 )
            return "#theta_{#gamma#gamma}";
        if ( variationName.CompareTo("MesonResolution") == 0 )
            return "#pi^{0} resolution";
        if ( variationName.CompareTo("Secondary") == 0 )
            return "secondary corr.";
        if ( variationName.CompareTo("CellMinE") == 0 )
            return "cell E_{agg}";
        if ( variationName.CompareTo("Cocktail") == 0 )
            return "cocktail";

        cout << "GetSystematicsName(): no valid input parameters given..." << endl;
        return "";

    }

    void    PlotLabelsInvMassInPtPlots          (   Double_t startTextX,
                                                    Double_t startTextY,
                                                    Double_t textHeight,
                                                    Double_t differenceText,
                                                    TString textAlice,
                                                    TString dateDummy,
                                                    TString fEnergy,
                                                    TString fDecayChannel,
                                                    TString fDetectionChannel,
                                                    TString textEvents          = "",
                                                    Double_t fNEvents           = 0,
                                                    Int_t textAlign             = 11
    ){

        TLatex *alice           = new TLatex(startTextX, startTextY, Form("%s",textAlice.Data()));
        SetStyleTLatex( alice, textHeight*1.3, 1, 1, 42, kTRUE, textAlign);
        alice->Draw();
        TLatex *latexDate       = new TLatex(startTextX, (startTextY-1.25*differenceText), dateDummy.Data());
        SetStyleTLatex( latexDate, textHeight, 1, 1, 42, kTRUE, textAlign);
        latexDate->Draw();
        TLatex *energy          = new TLatex(startTextX, (startTextY-2.25*differenceText), fEnergy);
        SetStyleTLatex( energy, textHeight*1, 1, 1, 42, kTRUE, textAlign);
        energy->Draw();
        TLatex *process         = new TLatex(startTextX, (startTextY-3.25*differenceText), fDecayChannel);
        SetStyleTLatex( process, textHeight*1, 1, 1, 42, kTRUE, textAlign);
        process->Draw();
        TLatex *detprocess      = new TLatex(startTextX, (startTextY-4.25*differenceText), fDetectionChannel);
        SetStyleTLatex( detprocess, textHeight*1, 1, 1, 42, kTRUE, textAlign);
        detprocess->Draw();
        TLatex *events          = NULL;
        if (textEvents.CompareTo("") != 0 && fNEvents > 0){
            events              = new TLatex(startTextX, (startTextY-5.25*differenceText), Form("%s: %2.1e events",textEvents.Data(), fNEvents));
            SetStyleTLatex( events, textHeight*1, 1, 1, 42, kTRUE, textAlign);
            events->Draw();
        }
    }


#endif
