#include "TLatex.h"

//void StyleSettingsThesis();
void StyleSettings();

class TGradientParFunction {
	
	public:
		
		TGradientParFunction(int ipar, TF1 * f)  :
		fPar(ipar),
		fFunc(f)
		{}
		
		double operator() (double * x, double *) const
		{
			// evaluate gradient vector of functions at point x
			return fFunc->GradientPar(fPar,x);
		}
		
	private:
		
		unsigned int fPar;
		mutable TF1 * fFunc;
};


// void StyleSettings(){
// 	
// 	//gStyle->SetOptTitle(kFALSE);
// 	gStyle->SetOptDate(0);   //show day and time
// 	gStyle->SetOptStat(0);  //show statistic
// 	gStyle->SetPalette(1,0);
// 	gStyle->SetFrameBorderMode(0);
// 	gStyle->SetFrameFillColor(0);
// 	gStyle->SetTitleFillColor(0);
// 	gStyle->SetTextSize(0.5);
// 	gStyle->SetLabelSize(0.03,"xyz");
// 	gStyle->SetLabelOffset(0.002,"xyz");
// 	gStyle->SetTitleFontSize(0.04);
// 	gStyle->SetTitleOffset(1,"y");
// 	gStyle->SetTitleOffset(0.7,"x");
// 	gStyle->SetCanvasColor(0);
// 	gStyle->SetPadTickX(1);
// 	gStyle->SetPadTickY(1);
// 	gStyle->SetLineWidth(0.01);
// 	
// 	gStyle->SetPadTopMargin(0.1);
// 	gStyle->SetPadBottomMargin(0.05);
// 	gStyle->SetPadRightMargin(0.05);
// 	gStyle->SetPadLeftMargin(0.05);
// 	
// 	
// 	//	TGaxis::SetMaxDigits(3);
// 	
// }

/* DrawAutoGammaHisto is function used for styling a histograma of the gamma conversion group with standart settings
* histo1 - first histogram (Data)
* Title - histogram title
* XTitle - X-axis title
* YTitle - Y-axis title
*xMin - minimum Y
*xMax - maximum Y
*/

//MY ADDITION!!!
void DrawSingleHisto( TH1* histo1, TString Title, TString XTitle, TString YTitle, Bool_t is2D) {

	histo1->SetTitle("");

	if(XTitle.Length() > 0){
		histo1->SetXTitle(XTitle.Data());
	}
	if(YTitle.Length() > 0){
		histo1->SetYTitle(YTitle.Data());
	}
	//cout<<"HERE??"<<endl;
	histo1->GetYaxis()->SetLabelSize(0.02);
	histo1->GetYaxis()->SetTitleSize(0.025);
	histo1->GetYaxis()->SetDecimals();
	histo1->GetYaxis()->SetTitleOffset(0.5);
	histo1->GetXaxis()->SetTitleSize(0.025);
	histo1->GetXaxis()->SetLabelSize(0.02);
	histo1->SetMarkerStyle(20);
	histo1->SetMarkerColor(1);
	histo1->SetLineColor(37);
	histo1->SetLineWidth(1);
	histo1->SetMarkerSize(0.2);
	histo1->SetTitleOffset(1.2,"xy");		
	histo1->SetTitleSize(0.05,"xy");		
	histo1->GetYaxis()->SetLabelSize(0.05);
	histo1->GetXaxis()->SetLabelSize(0.05);
	
	if (is2D)
	histo1->DrawCopy("colz");
	else
	histo1->DrawCopy("e,hist");
	/*	
	if(Title.Length() > 0){
		histo1->SetTitle("");
		TLatex *alice = new TLatex(0.1,0.95,Form("%s",Title.Data())); // Bo: this was 
		alice->SetNDC();
		alice->SetTextColor(1);
		alice->SetTextSize(0.062);
		alice->Draw();		
	}*/
}
//END OF MY ADDITION!!!


void DrawGammaHisto( TH1* histo1,
				 TString Title, TString XTitle, TString YTitle,
				 Float_t xMin, Float_t xMax,Int_t bck) {
	
	histo1->GetXaxis()->SetRangeUser(xMin, xMax);
	if(Title.Length() > 0){
		histo1->SetTitle("");
		TLatex *alice = new TLatex(0.1,0.95,Form("%s",Title.Data())); // Bo: this was 
		alice->SetNDC();
		alice->SetTextColor(1);
		alice->SetTextSize(0.062);
		alice->Draw();		
	}
	if (bck != 1 || bck != 2){
      Double_t yMin = 0;
      Double_t yMax = 0;
      for (Int_t i = histo1->GetXaxis()->FindBin(xMin); i < histo1->GetXaxis()->FindBin(xMax); i++){
         if (histo1->GetBinContent(i) < yMin){
            yMin = histo1->GetBinContent(i);
         }
         if (histo1->GetBinContent(i) > yMax){
            yMax = histo1->GetBinContent(i);
         }
      }   
      if (xMin > 0.2)  histo1->GetYaxis()->SetRangeUser(yMin, 1.5*yMax);
         else histo1->GetYaxis()->SetRangeUser(yMin, 1.2*yMax);

   }

	if(XTitle.Length() > 0){
		histo1->SetXTitle(XTitle.Data());
	}
	if(YTitle.Length() > 0){
		histo1->SetYTitle(YTitle.Data());
	}
	histo1->GetYaxis()->SetLabelSize(0.02);
	histo1->GetYaxis()->SetTitleSize(0.025);
	histo1->GetYaxis()->SetDecimals();
	histo1->GetYaxis()->SetTitleOffset(0.5);
	histo1->GetXaxis()->SetTitleSize(0.025);
	histo1->GetXaxis()->SetLabelSize(0.02);
	histo1->SetMarkerStyle(20);
	histo1->SetMarkerColor(1);
	histo1->SetLineColor(1);
	histo1->SetLineWidth(0.5);
	histo1->SetMarkerSize(0.2);
	histo1->SetTitleOffset(1.2,"xy");		
	histo1->SetTitleSize(0.05,"xy");		
	histo1->GetYaxis()->SetLabelSize(0.05);
	histo1->GetXaxis()->SetLabelSize(0.05);
	histo1->GetXaxis()->SetNdivisions(507,kTRUE);
	
	switch (bck)
	{
		default: histo1->DrawCopy("e1,p"); break;
		case 1:
			histo1->SetLineStyle(1);		
			histo1->SetLineColor(4);
			histo1->SetMarkerColor(4);
			histo1->SetMarkerStyle(24);
			histo1->SetLineWidth(0.7);
			histo1->DrawCopy("hist,same"); break;
		case 2:
			histo1->DrawCopy("same"); break;
		case 3:
			histo1->SetLineStyle(1);		
			histo1->SetLineColor(1);
			histo1->SetMarkerColor(4);
			histo1->SetMarkerStyle(24);
			histo1->SetLineWidth(0.7);
			histo1->DrawCopy("hist,same"); break;
	}
	/*
	if(bck==1){
		histo1->SetLineStyle(1);		
		histo1->SetLineColor(4);
		histo1->SetMarkerColor(4);
		histo1->SetMarkerStyle(24);
		histo1->SetLineWidth(0.7);
		histo1->DrawCopy("hist,same");
	}else{
		if(bck==2){
			histo1->DrawCopy("same");
		}else {
			histo1->DrawCopy("e1,p");
		}
	}*/
}

void DrawGammaHistoColored( TH1* histo1,
				 TString Title, TString XTitle, TString YTitle,
				 Float_t xMin, Float_t xMax,Int_t firstTime, Int_t color) {
	
	histo1->GetXaxis()->SetRangeUser(xMin, xMax);
	if(Title.Length() > 0){
		histo1->SetTitle("");
		TLatex *alice = new TLatex(0.1,0.95,Form("%s",Title.Data())); // Bo: this was 
		alice->SetNDC();
		alice->SetTextColor(1);
		alice->SetTextSize(0.062);
		alice->Draw();		
	}
	if (firstTime == 1){
      Double_t yMin = 0;
      Double_t yMax = 0;
      for (Int_t i = histo1->GetXaxis()->FindBin(xMin); i < histo1->GetXaxis()->FindBin(xMax); i++){
         if (histo1->GetBinContent(i) < yMin){
            yMin = histo1->GetBinContent(i);
         }
         if (histo1->GetBinContent(i) > yMax){
            yMax = histo1->GetBinContent(i);
         }
      }   
      if (xMin > 0.2)  histo1->GetYaxis()->SetRangeUser(yMin, 1.5*yMax);
         else histo1->GetYaxis()->SetRangeUser(yMin, 1.2*yMax);
   }
	if(XTitle.Length() > 0){
		histo1->SetXTitle(XTitle.Data());
	}
	if(YTitle.Length() > 0){
		histo1->SetYTitle(YTitle.Data());
	}
	histo1->GetYaxis()->SetLabelSize(0.02);
	histo1->GetYaxis()->SetTitleSize(0.025);
	histo1->GetYaxis()->SetDecimals();
	histo1->GetYaxis()->SetTitleOffset(0.5);
	histo1->GetXaxis()->SetTitleSize(0.025);
	histo1->GetXaxis()->SetLabelSize(0.02);
	histo1->SetMarkerStyle(20);
	histo1->SetMarkerColor(color);
	histo1->SetLineColor(color);
	histo1->SetLineWidth(0.5);
	histo1->SetMarkerSize(0.2);
	histo1->SetTitleOffset(1.2,"xy");		
	histo1->SetTitleSize(0.05,"xy");		
	histo1->GetYaxis()->SetLabelSize(0.05);
	histo1->GetXaxis()->SetLabelSize(0.05);
	histo1->GetXaxis()->SetNdivisions(507,kTRUE);
	if (firstTime == 1) {
		histo1->DrawCopy("hist");
	} else {	
		histo1->DrawCopy("hist,same");	
	}
}

void DrawGammaHistoBigger( TH1* histo1,
				 TString Title, TString XTitle, TString YTitle,
				 Float_t xMin, Float_t xMax,Int_t bck) {
	
	histo1->GetXaxis()->SetRangeUser(xMin, xMax);
   if (bck != 1 || bck != 2){
      Double_t yMin = 0;
      Double_t yMax = 0;
      for (Int_t i = histo1->GetXaxis()->FindBin(xMin); i < histo1->GetXaxis()->FindBin(xMax); i++){
         if (histo1->GetBinContent(i) < yMin){
            yMin = histo1->GetBinContent(i);
         }
         if (histo1->GetBinContent(i) > yMax){
            yMax = histo1->GetBinContent(i);
         }
      }   
      if (xMin > 0.2)  histo1->GetYaxis()->SetRangeUser(yMin, 1.5*yMax);
         else histo1->GetYaxis()->SetRangeUser(yMin, 1.2*yMax);
   }
	if(Title.Length() > 0){
		histo1->SetTitle("");
		TLatex *alice = new TLatex(0.1,0.95,Form("%s",Title.Data())); // Bo: this was 
		alice->SetNDC();
		alice->SetTextColor(1);
		alice->SetTextSize(0.04);
		alice->Draw();		
	}
	if(XTitle.Length() > 0){
		histo1->SetXTitle(XTitle.Data());
	}
	if(YTitle.Length() > 0){
		histo1->SetYTitle(YTitle.Data());
	}
	histo1->GetYaxis()->SetLabelSize(0.04);
	histo1->GetYaxis()->SetTitleSize(0.045);
	histo1->GetYaxis()->SetDecimals();
	histo1->GetYaxis()->SetTitleOffset(1.2);
	histo1->GetXaxis()->SetTitleOffset(1.);
	histo1->GetXaxis()->SetTitleSize(0.045);
	histo1->GetXaxis()->SetLabelSize(0.04);
	histo1->SetMarkerStyle(20);
	histo1->SetMarkerColor(1);
	histo1->SetLineColor(1);
	histo1->SetLineWidth(1);
	histo1->SetMarkerSize(0.8);
/*	histo1->SetTitleOffset(1.2,"xy");		
	histo1->SetTitleSize(0.05,"xy");		*/
/*	histo1->GetYaxis()->SetLabelSize(0.05);
	histo1->GetXaxis()->SetLabelSize(0.05);*/
	histo1->GetXaxis()->SetNdivisions(507,kTRUE);
	if(bck==1){
		histo1->SetLineStyle(1);		
		histo1->SetLineColor(4);
		histo1->SetMarkerColor(4);
		histo1->SetMarkerStyle(24);
		histo1->SetLineWidth(2.);
		histo1->DrawCopy("hist,same");
	}else{
		if(bck==2){
			histo1->DrawCopy("same");
		}else {
			histo1->DrawCopy("e1,p");
		}
	}
}