/**********************************************************************************************
******      Friederike Bock, friederike.bock@cern.ch                                      *****
******		Lucia Leardini, lucia.leardini@cern.ch										  *****
**********************************************************************************************/

#include <Riostream.h>
#include "TMath.h"
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TPaveLabel.h>
#include <TSystem.h>
#include <TFrame.h>
#include <TStyle.h>
#include <TString.h>
#include "TGaxis.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TF1.h"
#include "TVirtualFitter.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TDatabasePDG.h"
#include "TMinuit.h"
#include "TBenchmark.h"
#include "TRandom.h"
#include "TLatex.h"
#include "TASImage.h"
#include "TPostScript.h"
#include "TGraphErrors.h"
#include "TArrow.h"
#include "TGraphAsymmErrors.h" 
#include "TGaxis.h"
#include "TMarker.h"
#include "Math/WrappedTF1.h"
#include "Math/BrentRootFinder.h"
#include "CommonHeaders/PlottingGammaConversionHistos.h"
#include "CommonHeaders/PlottingGammaConversionAdditional.h"
#include "CommonHeaders/FittingGammaConversion.h"
#include "CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "CommonHeaders/ConversionFunctions.h"


extern TRandom*   gRandom;
extern TBenchmark*   gBenchmark;
extern TSystem*   gSystem;
extern TMinuit*   gMinuit;

Double_t		xSectionpp7TeV =			62.22*1e-3;
Double_t		xSectionpp7TeVV0AND =			54.31*1e-3;
Double_t		xSectionpp7TeVErrUp =		2.18;
Double_t		xSectionpp7TeVErrDown =	2.18;
Double_t		xSectionpp900GeV =		47.78*1e-3;
Double_t		xSectionpp900GeVV0AND =		40.06*1e-3;
Double_t		xSectionpp900GeVErrUp =		 2.39;
Double_t		xSectionpp900GeVErrDown =		 1.86;
Double_t		xSectionpp2760GeV = 		55.416*1e-3;
Double_t		xSectionpp2760GeVV0AND = 		47.73*1e-3;
Double_t		xSectionpp2760GeVErr = 	3.9;
// Double_t 		recalcBarn = 			1e12; //NLO in pbarn!!!!
Double_t 		xSection7TeVppINEL = 73.2*1e9;
Double_t 		xSection2760GeVppINEL = 62.8*1e9;
Double_t 		xSection900GeVppINEL = 52.5*1e9;


void get_phi2pi(TGraphErrors* g,TBox** b){
  /*
  TH2D* h=phi_final_over_pix_c2
  int j
  double x[8],v[8],t[8],y[8]
  for(int j=0;j<8;j++) x[j]=h->GetXaxis()->GetBinCenter(j+1)
  for(int j=0;j<8;j++) v[j]=h->GetBinContent(j+1,1)
  for(int j=0;j<8;j++) t[j]=h->GetBinError(j+1,1)
  for(int j=0;j<8;j++) y[j]=h->GetBinContent(j+1,2)
  for(int j=0;j<8;j++) printf("  g->SetPoint(%i,%1.2f,a*%1.9e); g->SetPointError(%i,0.,a*%1.9e);\n  b[%i]=new TBox(%1.2f-dx,a*%1.9e,%1.2f+dx,a*%1.9e);\n",j,x[j],v[j],j,t[j],j,x[j],v[j]-y[j],x[j],v[j]+y[j])
  */

  double dx=0.1,a=1;
  g->SetPoint(0,0.65,a*1.353657541e-02); g->SetPointError(0,0.,a*6.957752378e-04);
  b[0]=new TBox(0.65-dx,a*1.229208816e-02,0.65+dx,a*1.478106265e-02);
  g->SetPoint(1,0.90,a*3.070164564e-02); g->SetPointError(1,0.,a*1.151929435e-03);
  b[1]=new TBox(0.90-dx,a*2.821480778e-02,0.90+dx,a*3.318848349e-02);
  g->SetPoint(2,1.25,a*5.945575846e-02); g->SetPointError(2,0.,a*3.003430892e-03);
  b[2]=new TBox(1.25-dx,a*5.450533182e-02,1.25+dx,a*6.440618510e-02);
  g->SetPoint(3,1.75,a*1.386005779e-01); g->SetPointError(3,0.,a*7.727299101e-03);
  b[3]=new TBox(1.75-dx,a*1.229326293e-01,1.75+dx,a*1.542685265e-01);
  g->SetPoint(4,2.25,a*2.337683977e-01); g->SetPointError(4,0.,a*1.495423881e-02);
  b[4]=new TBox(2.25-dx,a*2.096313435e-01,2.25+dx,a*2.579054519e-01);
  g->SetPoint(5,2.75,a*3.033753113e-01); g->SetPointError(5,0.,a*2.122702618e-02);
  b[5]=new TBox(2.75-dx,a*2.718312198e-01,2.75+dx,a*3.349194027e-01);
  g->SetPoint(6,3.50,a*3.446229720e-01); g->SetPointError(6,0.,a*2.260536484e-02);
  b[6]=new TBox(3.50-dx,a*2.932587157e-01,3.50+dx,a*3.959872283e-01);
  g->SetPoint(7,4.50,a*3.294960583e-01); g->SetPointError(7,0.,a*3.560356503e-02);
  b[7]=new TBox(4.50-dx,a*2.681613542e-01,4.50+dx,a*3.908307624e-01);

  return;
  
  
  
}

void get_phi2pipp7TeV(TGraphErrors* g,TBox** b){
	
  double dx=0.1,a=1;
  g->SetPoint(0,0.45,a*1.119247202e-02); g->SetPointError(0,0.,a*8.994972179e-04);
  b[0]=new TBox(0.40,a*(1.119247202e-02-1.101869365e-03),0.50,a*(1.119247202e-02+1.101869365e-03));
  g->SetPoint(1,0.55,a*1.533208402e-02); g->SetPointError(1,0.,a*6.295941749e-04);
  b[1]=new TBox(0.50,a*(1.533208402e-02-1.581905727e-03),0.60,a*(1.533208402e-02+1.581905727e-03));
  g->SetPoint(2,0.65,a*1.864215424e-02); g->SetPointError(2,0.,a*5.232385277e-04);
  b[2]=new TBox(0.60,a*(1.864215424e-02-1.868848646e-03),0.70,a*(1.864215424e-02+1.868848646e-03));
  g->SetPoint(3,0.75,a*2.452232769e-02); g->SetPointError(3,0.,a*5.496848650e-04);
  b[3]=new TBox(0.70,a*(2.452232769e-02-2.514044412e-03),0.80,a*(2.452232769e-02+2.514044412e-03));
  g->SetPoint(4,0.85,a*2.890036189e-02); g->SetPointError(4,0.,a*6.062008267e-04);
  b[4]=new TBox(0.80,a*(2.890036189e-02-2.956792689e-03),0.90,a*(2.890036189e-02+2.956792689e-03));
  g->SetPoint(5,0.95,a*3.571668561e-02); g->SetPointError(5,0.,a*7.334431535e-04);
  b[5]=new TBox(0.90,a*(3.571668561e-02-3.547199518e-03),1.0,a*(3.571668561e-02+3.547199518e-03));
  g->SetPoint(6,1.05,a*4.189824060e-02); g->SetPointError(6,0.,a*8.999940406e-04);
  b[6]=new TBox(1.0,a*(4.189824060e-02-4.026962770e-03),1.1,a*(4.189824060e-02+4.026962770e-03));
  g->SetPoint(7,1.15,a*4.949514549e-02); g->SetPointError(7,0.,a*1.106071036e-03);
  b[7]=new TBox(1.1,a*(4.949514549e-02-4.988083954e-03),1.2,a*(4.949514549e-02+4.988083954e-03));
  g->SetPoint(8,1.25,a*5.651248406e-02); g->SetPointError(8,0.,a*1.363279269e-03);
  b[8]=new TBox(1.2,a*(5.651248406e-02-6.658769101e-03),1.3,a*(5.651248406e-02+6.658769101e-03));
  g->SetPoint(9,1.35,a*6.307407920e-02); g->SetPointError(9,0.,a*1.625844248e-03);
  b[9]=new TBox(1.3,a*(6.307407920e-02-7.587164846e-03),1.4,a*(6.307407920e-02+7.587164846e-03));
  g->SetPoint(10,1.45,a*6.801247840e-02); g->SetPointError(10,0.,a*1.976728860e-03);
  b[10]=new TBox(1.4,a*(6.801247840e-02-7.668559573e-03),1.5,a*(6.801247840e-02+7.668559573e-03));
  g->SetPoint(11,1.55,a*7.533635873e-02); g->SetPointError(11,0.,a*2.274333731e-03);
  b[11]=new TBox(1.5,a*(7.533635873e-02-9.005238564e-03),1.6,a*(7.533635873e-02+9.005238564e-03));
  g->SetPoint(12,1.65,a*8.174400361e-02); g->SetPointError(12,0.,a*2.601051492e-03);
  b[12]=new TBox(1.6,a*(8.174400361e-02-9.625709494e-03),1.7,a*(8.174400361e-02+9.625709494e-03));
  g->SetPoint(13,1.75,a*8.764237237e-02); g->SetPointError(13,0.,a*2.908019836e-03);
  b[13]=new TBox(1.7,a*(8.764237237e-02-9.771773650e-03),1.8,a*(8.764237237e-02+9.771773650e-03));
  g->SetPoint(14,1.85,a*8.917655620e-02); g->SetPointError(14,0.,a*3.161382295e-03);
  b[14]=new TBox(1.8,a*(8.917655620e-02-9.568485948e-03),1.9,a*(8.917655620e-02+9.568485948e-03));
  g->SetPoint(15,1.95,a*8.857730995e-02); g->SetPointError(15,0.,a*3.300159174e-03);
  b[15]=new TBox(1.9,a*(8.857730995e-02-1.039718102e-02),2.0,a*(8.857730995e-02+1.039718102e-02));
  g->SetPoint(16,2.1,a*9.587081473e-02); g->SetPointError(16,0.,a*2.629077487e-03);
  b[16]=new TBox(2.05,a*(9.587081473e-02-1.005044416e-02),2.15,a*(9.587081473e-02+1.005044416e-02));
  g->SetPoint(17,2.3,a*1.044300342e-01); g->SetPointError(17,0.,a*3.151066585e-03);
  b[17]=new TBox(2.25,a*(1.044300342e-01-1.152118092e-02),2.35,a*(1.044300342e-01+1.152118092e-02));
  g->SetPoint(18,2.5,a*1.050055171e-01); g->SetPointError(18,0.,a*3.487316359e-03);
  b[18]=new TBox(2.45,a*(1.050055171e-01-1.372232580e-02),2.55,a*(1.050055171e-01+1.372232580e-02));
  g->SetPoint(19,2.7,a*1.149071190e-01); g->SetPointError(19,0.,a*4.192586576e-03);
  b[19]=new TBox(2.65,a*(1.149071190e-01-1.509959004e-02),2.75,a*(1.149071190e-01+1.509959004e-02));
  g->SetPoint(20,2.9,a*1.182595988e-01); g->SetPointError(20,0.,a*4.801603508e-03);
  b[20]=new TBox(2.85,a*(1.182595988e-01-1.545271193e-02),2.95,a*(1.182595988e-01+1.545271193e-02));
  g->SetPoint(21,3.25,a*1.260789449e-01); g->SetPointError(21,0.,a*3.941883188e-03);
  b[21]=new TBox(3.2,a*(1.260789449e-01-1.529558821e-02),3.3,a*(1.260789449e-01+1.529558821e-02));
  g->SetPoint(22,3.75,a*1.409705043e-01); g->SetPointError(22,0.,a*5.574422171e-03);
  b[22]=new TBox(3.7,a*(1.409705043e-01-1.769544912e-02),3.8,a*(1.409705043e-01+1.769544912e-02));
  g->SetPoint(23,4.25,a*1.550439739e-01); g->SetPointError(23,0.,a*7.619842933e-03);
  b[23]=new TBox(4.2,a*(1.550439739e-01-1.962113206e-02),4.3,a*(1.550439739e-01+1.962113206e-02));
  g->SetPoint(24,4.75,a*1.637287914e-01); g->SetPointError(24,0.,a*1.059080219e-02);
  b[24]=new TBox(4.7,a*(1.637287914e-01-2.166891568e-02),4.8,a*(1.637287914e-01+2.166891568e-02));
  g->SetPoint(25,5.5,a*1.586523927e-01); g->SetPointError(25,0.,a*1.019251973e-02);
  b[25]=new TBox(5.45,a*(1.586523927e-01-2.542619154e-02),5.55,a*(1.586523927e-01+2.542619154e-02));
  
  return;
}

TGraphAsymmErrors* get_phipp7TeV(){

	// Plot: p8208_d2x1y1
	double p8208_d2x1y1_xval[] = { 0.45, 0.55, 0.65, 0.75, 0.85, 
								   0.95, 1.05, 1.15, 1.25, 1.35, 
								   1.45, 1.55, 1.65, 1.75, 1.85, 
								   1.95, 2.1, 2.3, 2.5, 2.7, 
								   2.9, 3.25, 3.75, 4.25, 4.75,
								   5.5 };
	double p8208_d2x1y1_xerrminus[] = { 0.04999999999999999, 0.050000000000000044, 0.050000000000000044, 0.050000000000000044, 0.04999999999999993, 
										0.04999999999999993, 0.050000000000000044, 0.04999999999999982, 0.050000000000000044, 0.050000000000000044,
										0.050000000000000044, 0.050000000000000044, 0.04999999999999982, 0.050000000000000044, 0.050000000000000044,
										0.050000000000000044, 0.10000000000000009, 0.09999999999999964, 0.10000000000000009, 0.10000000000000009,
										0.10000000000000009, 0.25, 0.25, 0.25, 0.25, 
										0.5 };
	double p8208_d2x1y1_xerrplus[] = { 	0.04999999999999999, 0.04999999999999993, 0.04999999999999993, 0.050000000000000044, 0.050000000000000044, 
										0.050000000000000044, 0.050000000000000044, 0.050000000000000044, 0.050000000000000044, 0.04999999999999982,
										0.050000000000000044, 0.050000000000000044, 0.050000000000000044, 0.050000000000000044, 0.04999999999999982, 
										0.050000000000000044, 0.10000000000000009, 0.10000000000000009, 0.10000000000000009, 0.09999999999999964, 
										0.10000000000000009, 0.25, 0.25, 0.25, 0.25, 
										0.5 };
	double p8208_d2x1y1_yval[] = { 		0.02464, 0.02461, 0.02284, 0.02289, 0.0202, 
										0.01873, 0.01653, 0.01502, 0.01333, 0.0117, 
										0.01003, 0.00903, 0.007908, 0.006922, 0.00577, 
										0.004736, 0.003904, 0.002982, 0.002141, 0.001679, 
										0.001309, 8.073E-4, 4.712E-4, 2.835E-4, 1.751E-4, 
										8.503E-5 };
	double p8208_d2x1y1_yerrminus[] = { 0.00291547594742265, 0.0023482972554597936, 0.0019859506539690254, 0.002015589243868899, 0.0017705366418123065, 
										0.001586190404711868, 0.0014236923825040294, 0.0013800362314084368, 0.0012419339757008018, 0.0011884864324004712, 
										9.646242791885347E-4, 9.965550662156106E-4, 8.821655173492103E-4, 7.333437120477682E-4, 5.956584591861345E-4, 
										5.382434393469185E-4, 3.820013088982811E-4, 3.1072174046886387E-4, 2.490401574043833E-4, 1.7369225659193908E-4, 
										1.3261975720080322E-4, 7.486848469149086E-5, 4.7309724158992934E-5, 2.9924237667817037E-5, 2.0999523804124702E-5, 
										1.2698952712723992E-5 };
	double p8208_d2x1y1_yerrplus[] = { 	0.00291547594742265, 0.0023482972554597936, 0.0019859506539690254, 0.002015589243868899, 0.0017705366418123065, 
										0.001586190404711868, 0.0014236923825040294, 0.0013800362314084368, 0.0012419339757008018, 0.0011884864324004712, 
										9.646242791885347E-4, 9.965550662156106E-4, 8.821655173492103E-4, 7.333437120477682E-4, 5.956584591861345E-4, 
										5.382434393469185E-4, 3.820013088982811E-4, 3.1072174046886387E-4, 2.490401574043833E-4, 1.7369225659193908E-4, 
										1.3261975720080322E-4, 7.486848469149086E-5, 4.7309724158992934E-5, 2.9924237667817037E-5, 2.0999523804124702E-5, 
										1.2698952712723992E-5 };
	double p8208_d2x1y1_ystatminus[] = { 	0.00198, 0.00101, 6.4E-4, 5.1E-4, 4.2E-4, 
											3.8E-4, 3.5E-4, 3.3E-4, 3.2E-4, 3.0E-4,
											2.9E-4, 2.71E-4, 2.5E-4, 2.28E-4, 2.03E-4,
											1.75E-4, 1.06E-4, 8.8E-5, 7.0E-5, 6.0E-5, 
											5.2E-5, 2.52E-5, 1.86E-5, 1.39E-5, 1.13E-5,
											5.45E-6 };
	double p8208_d2x1y1_ystatplus[] = { 	0.00198, 0.00101, 6.4E-4, 5.1E-4, 4.2E-4, 
											3.8E-4, 3.5E-4, 3.3E-4, 3.2E-4, 3.0E-4, 
											2.9E-4, 2.71E-4, 2.5E-4, 2.28E-4, 2.03E-4,
											1.75E-4, 1.06E-4, 8.8E-5, 7.0E-5, 6.0E-5,
											5.2E-5, 2.52E-5, 1.86E-5, 1.39E-5, 1.13E-5, 
											5.45E-6 };
	double p8208_d2x1y1_ysystminus[] = { 	0.00214, 0.00212, 0.00188, 0.00195, 0.00172,
											0.00154, 0.00138, 0.00134, 0.0012, 0.00115, 
											9.2E-4, 9.59E-4, 8.46E-4, 6.97E-4, 5.6E-4, 
											5.09E-4, 3.67E-4, 2.98E-4, 2.39E-4, 1.63E-4, 
											1.22E-4, 7.05E-5, 4.35E-5, 2.65E-5, 1.77E-5, 
											1.147E-5 };
																					
	int p8208_d2x1y1_numpoints = 26;
	for (Int_t i = 0; i < p8208_d2x1y1_numpoints; i++){
		p8208_d2x1y1_yval[i]= p8208_d2x1y1_yval[i]/(2*TMath::Pi()*p8208_d2x1y1_xval[i])*xSectionpp7TeV*recalcBarn;
		p8208_d2x1y1_yerrminus[i]= p8208_d2x1y1_yerrminus[i]/(2*TMath::Pi()*p8208_d2x1y1_xval[i])*xSectionpp7TeV*recalcBarn;
		p8208_d2x1y1_yerrplus[i]= p8208_d2x1y1_yerrplus[i]/(2*TMath::Pi()*p8208_d2x1y1_xval[i])*xSectionpp7TeV*recalcBarn;
		p8208_d2x1y1_ystatminus[i]= p8208_d2x1y1_ystatminus[i]/(2*TMath::Pi()*p8208_d2x1y1_xval[i])*xSectionpp7TeV*recalcBarn;
		p8208_d2x1y1_ystatplus[i]= p8208_d2x1y1_ystatplus[i]/(2*TMath::Pi()*p8208_d2x1y1_xval[i])*xSectionpp7TeV*recalcBarn;
		p8208_d2x1y1_ysystminus[i]= p8208_d2x1y1_ysystminus[i]/(2*TMath::Pi()*p8208_d2x1y1_xval[i])*xSectionpp7TeV*recalcBarn;
		
	}	
	
	TGraphAsymmErrors* p8208_d2x1y1 = new TGraphAsymmErrors(p8208_d2x1y1_numpoints, p8208_d2x1y1_xval, p8208_d2x1y1_yval, p8208_d2x1y1_xerrminus, p8208_d2x1y1_xerrplus, p8208_d2x1y1_yerrminus, p8208_d2x1y1_yerrplus);
	return p8208_d2x1y1;
}


TGraphAsymmErrors* get_phiToPi0pp7TeV(TF1* fitPi0){

	// Plot: p8208_d2x1y1
	double p8208_d2x1y1_xval[] = { 0.45, 0.55, 0.65, 0.75, 0.85, 
								   0.95, 1.05, 1.15, 1.25, 1.35, 
								   1.45, 1.55, 1.65, 1.75, 1.85, 
								   1.95, 2.1, 2.3, 2.5, 2.7, 
								   2.9, 3.25, 3.75, 4.25, 4.75,
								   5.5 };
	double p8208_d2x1y1_xerrminus[] = { 0.04999999999999999, 0.050000000000000044, 0.050000000000000044, 0.050000000000000044, 0.04999999999999993, 
										0.04999999999999993, 0.050000000000000044, 0.04999999999999982, 0.050000000000000044, 0.050000000000000044,
										0.050000000000000044, 0.050000000000000044, 0.04999999999999982, 0.050000000000000044, 0.050000000000000044,
										0.050000000000000044, 0.10000000000000009, 0.09999999999999964, 0.10000000000000009, 0.10000000000000009,
										0.10000000000000009, 0.25, 0.25, 0.25, 0.25, 
										0.5 };
	double p8208_d2x1y1_xerrplus[] = { 	0.04999999999999999, 0.04999999999999993, 0.04999999999999993, 0.050000000000000044, 0.050000000000000044, 
										0.050000000000000044, 0.050000000000000044, 0.050000000000000044, 0.050000000000000044, 0.04999999999999982,
										0.050000000000000044, 0.050000000000000044, 0.050000000000000044, 0.050000000000000044, 0.04999999999999982, 
										0.050000000000000044, 0.10000000000000009, 0.10000000000000009, 0.10000000000000009, 0.09999999999999964, 
										0.10000000000000009, 0.25, 0.25, 0.25, 0.25, 
										0.5 };
	double p8208_d2x1y1_yval[] = { 		0.02464, 0.02461, 0.02284, 0.02289, 0.0202, 
										0.01873, 0.01653, 0.01502, 0.01333, 0.0117, 
										0.01003, 0.00903, 0.007908, 0.006922, 0.00577, 
										0.004736, 0.003904, 0.002982, 0.002141, 0.001679, 
										0.001309, 8.073E-4, 4.712E-4, 2.835E-4, 1.751E-4, 
										8.503E-5 };
	double p8208_d2x1y1_yerrminus[] = { 0.00291547594742265, 0.0023482972554597936, 0.0019859506539690254, 0.002015589243868899, 0.0017705366418123065, 
										0.001586190404711868, 0.0014236923825040294, 0.0013800362314084368, 0.0012419339757008018, 0.0011884864324004712, 
										9.646242791885347E-4, 9.965550662156106E-4, 8.821655173492103E-4, 7.333437120477682E-4, 5.956584591861345E-4, 
										5.382434393469185E-4, 3.820013088982811E-4, 3.1072174046886387E-4, 2.490401574043833E-4, 1.7369225659193908E-4, 
										1.3261975720080322E-4, 7.486848469149086E-5, 4.7309724158992934E-5, 2.9924237667817037E-5, 2.0999523804124702E-5, 
										1.2698952712723992E-5 };
	double p8208_d2x1y1_yerrplus[] = { 	0.00291547594742265, 0.0023482972554597936, 0.0019859506539690254, 0.002015589243868899, 0.0017705366418123065, 
										0.001586190404711868, 0.0014236923825040294, 0.0013800362314084368, 0.0012419339757008018, 0.0011884864324004712, 
										9.646242791885347E-4, 9.965550662156106E-4, 8.821655173492103E-4, 7.333437120477682E-4, 5.956584591861345E-4, 
										5.382434393469185E-4, 3.820013088982811E-4, 3.1072174046886387E-4, 2.490401574043833E-4, 1.7369225659193908E-4, 
										1.3261975720080322E-4, 7.486848469149086E-5, 4.7309724158992934E-5, 2.9924237667817037E-5, 2.0999523804124702E-5, 
										1.2698952712723992E-5 };
	double p8208_d2x1y1_ystatminus[] = { 	0.00198, 0.00101, 6.4E-4, 5.1E-4, 4.2E-4, 
											3.8E-4, 3.5E-4, 3.3E-4, 3.2E-4, 3.0E-4,
											2.9E-4, 2.71E-4, 2.5E-4, 2.28E-4, 2.03E-4,
											1.75E-4, 1.06E-4, 8.8E-5, 7.0E-5, 6.0E-5, 
											5.2E-5, 2.52E-5, 1.86E-5, 1.39E-5, 1.13E-5,
											5.45E-6 };
	double p8208_d2x1y1_ystatplus[] = { 	0.00198, 0.00101, 6.4E-4, 5.1E-4, 4.2E-4, 
											3.8E-4, 3.5E-4, 3.3E-4, 3.2E-4, 3.0E-4, 
											2.9E-4, 2.71E-4, 2.5E-4, 2.28E-4, 2.03E-4,
											1.75E-4, 1.06E-4, 8.8E-5, 7.0E-5, 6.0E-5,
											5.2E-5, 2.52E-5, 1.86E-5, 1.39E-5, 1.13E-5, 
											5.45E-6 };
	double p8208_d2x1y1_ysystminus[] = { 	0.00214, 0.00212, 0.00188, 0.00195, 0.00172,
											0.00154, 0.00138, 0.00134, 0.0012, 0.00115, 
											9.2E-4, 9.59E-4, 8.46E-4, 6.97E-4, 5.6E-4, 
											5.09E-4, 3.67E-4, 2.98E-4, 2.39E-4, 1.63E-4, 
											1.22E-4, 7.05E-5, 4.35E-5, 2.65E-5, 1.77E-5, 
											1.147E-5 };
																					
	int p8208_d2x1y1_numpoints = 26;
	for (Int_t i = 0; i < p8208_d2x1y1_numpoints; i++){
		Double_t integralPi0 = fitPi0->Integral(p8208_d2x1y1_xval[i]-p8208_d2x1y1_xerrminus[i], p8208_d2x1y1_xval[i]+p8208_d2x1y1_xerrplus[i])/(p8208_d2x1y1_xerrminus[i]+p8208_d2x1y1_xerrplus[i]);
		
		p8208_d2x1y1_yval[i]= p8208_d2x1y1_yval[i]/(integralPi0*2*TMath::Pi()*p8208_d2x1y1_xval[i])*xSectionpp7TeV*recalcBarn;
// 		cout << i << "\t" <<p8208_d2x1y1_yval[i] << "\t"<< integralPi0 << endl;
		p8208_d2x1y1_yerrminus[i]= p8208_d2x1y1_yerrminus[i]/(integralPi0*2*TMath::Pi()*p8208_d2x1y1_xval[i])*xSectionpp7TeV*recalcBarn;
		p8208_d2x1y1_yerrplus[i]= p8208_d2x1y1_yerrplus[i]/(integralPi0*2*TMath::Pi()*p8208_d2x1y1_xval[i])*xSectionpp7TeV*recalcBarn;
		p8208_d2x1y1_ystatminus[i]= p8208_d2x1y1_ystatminus[i]/(integralPi0*2*TMath::Pi()*p8208_d2x1y1_xval[i])*xSectionpp7TeV*recalcBarn;
		p8208_d2x1y1_ystatplus[i]= p8208_d2x1y1_ystatplus[i]/(integralPi0*2*TMath::Pi()*p8208_d2x1y1_xval[i])*xSectionpp7TeV*recalcBarn;
		p8208_d2x1y1_ysystminus[i]= p8208_d2x1y1_ysystminus[i]/(integralPi0*2*TMath::Pi()*p8208_d2x1y1_xval[i])*xSectionpp7TeV*recalcBarn;
		
	}	
	
	TGraphAsymmErrors* p8208_d2x1y1 = new TGraphAsymmErrors(p8208_d2x1y1_numpoints, p8208_d2x1y1_xval, p8208_d2x1y1_yval, p8208_d2x1y1_xerrminus, p8208_d2x1y1_xerrplus, p8208_d2x1y1_yerrminus, p8208_d2x1y1_yerrplus);
	return p8208_d2x1y1;
}


void CombineEtaToPi0RatiosLHC11h(TString suffix = "pdf", TString nameFilePbPb = "data_PCMResults_PbPb.root", Bool_t runDrawReweighted = kTRUE, Bool_t checkStatFluc = kTRUE){


	gROOT->Reset();   
	gROOT->SetStyle("Plain");

	TString dateForOutput = ReturnDateStringForOutput();
	TString outputDir = Form("%s/%s/CombineEtaToPi0Ratios",suffix.Data(),dateForOutput.Data());
	gSystem->Exec("mkdir -p "+outputDir);
	gSystem->Exec(Form("cp %s %s/InputFilePCMPbPb.root ",nameFilePbPb.Data(),outputDir.Data() ));

	StyleSettingsThesis();  
	SetPlotStyle();

// 	Color_t  colorCombpPb         	   = kBlack;
	Color_t  colorCombpPb         	   = kGreen+3;

	Color_t  colorCombpPb0020          = kRed+1;
	Color_t  colorCombpPb2050          = kGreen+2;
	Color_t  colorCombpPb4060          = kCyan+2;
	Color_t  colorCombpPb6080          = kBlue+1;
	Color_t  colorCombpPb60100         = kViolet+3;

	Style_t  markerStylepPb  = 34 ;
	Style_t  markerStylepPb0020   = 20 ;
	Style_t  markerStylepPb2050   = 21 ;
	Style_t  markerStylepPb4060   = 29 ;
	Style_t  markerStylepPb6080   = 33 ;
	Style_t  markerStylepPb60100 = 34 ;

	Size_t   markerSizepPb  = 2.;
	Size_t   markerSizepPb0020  = 2.;
	Size_t   markerSizepPb2050  = 2.;
	Size_t   markerSizepPb4060  = 2.5;
	Size_t   markerSizepPb6080  = 2.5;
	Size_t   markerSizepPb60100  = 2.;

	Color_t colorPi0 = kRed+2;
	Color_t colorEta = kMagenta+2;
	Color_t colorOmega = kBlue+2;
	Color_t colorPhi = kGreen+2;
	
	Color_t  color900GeV 			= kRed +2;
	Color_t  color2760GeV 			= kMagenta+2;
	Color_t  color7TeV				= kBlue+2;
	Color_t  color900GeVBox = color900GeV-10;
	Color_t  color2760GeVBox = color2760GeV-10;
	Color_t  color7TeVBox = color7TeV-10;

   Color_t  colorCombPbPb0010          = kRed+1;
   Color_t  colorCombPbPb2050          = kAzure+1;
   
	Style_t  markerStyleSpectrum7TeV 	= 21 ;
	Style_t  markerStyleSpectrum900GeV = 21 ;
	Style_t  markerStyleSpectrum2760GeV = 29 ;
	Style_t  markerStylePbPb0010  = 20 ;
    Style_t  markerStylePbPb2050  = 21 ;
 
	Size_t  markerSizePi0PP7TeV 	= 1.8;
	Size_t  markerSizePi0PP900GeV = 1.8;
	Size_t  markerSizePi0PP2760GeV 	= 2.2;

	TString collisionSystemPbPb0010 = "0-10% Pb-Pb #sqrt{s_{_{NN}}} = 2.76 TeV"	;      
	TString collisionSystemPbPb2050 = "20-50% Pb-Pb #sqrt{s_{_{NN}}} = 2.76 TeV";      
	
	TString collisionSystemPP2760GeV = "pp #sqrt{#it{s}} = 2.76 TeV";		
	TString collisionSystemPP7TeV = "pp #sqrt{#it{s}} = 7 TeV";		
	TString collisionSystemPP900GeV = "pp #sqrt{#it{s}} = 0.9 TeV";		
	TString collisionSystempPb = "p-Pb  #sqrt{s_{_{NN}}} = 5.02 TeV";
		
	
	cout << "PCM Eta for pPb" << endl;
	TFile* filePCMpPb 					= new TFile("data_PCMResults_pPb_20150624_standard_dc4.root");
	TDirectory* fEtaToPi0pPbContainer 	= (TDirectory*) filePCMpPb->GetDirectory("Eta_pPb_5.023TeV_0-100%");
	TH1D* histoPCMEtaToPi0RatiopPb 		= (TH1D*)fEtaToPi0pPbContainer->Get("EtatoPi0Ratio");
	TGraphAsymmErrors* graphPCMEtaToPi0RatioSysErrpPb=    (TGraphAsymmErrors*)fEtaToPi0pPbContainer->Get("EtatoPi0RatioSys"); 

	cout << "PCM Eta for pp" << endl;
	TFile* fileCombinedpp 				= new TFile("FinalResults/CombinedResultsPP_ShiftedX_PaperRAA_16_May_2014.root ");
// 	TFile* fileCombinedpp 				= new TFile("CombinedResultsPaperX_18_Feb_2014.root");
	TGraphAsymmErrors* graphCombEtaToPi0Ratiopp7TeV =         (TGraphAsymmErrors*)fileCombinedpp->Get("graphEtaToPi0Comb7TeVStat");
	TGraphAsymmErrors* graphCombEtaToPi0RatioSysErrpp7TeV=    (TGraphAsymmErrors*)fileCombinedpp->Get("graphEtaToPi0Comb7TeVSys"); 

	TF1* fitEtaToPi0pp7TeVHighPt = new TF1("fitEtaToPi0pp7TeVHighPt","[0]");
	graphCombEtaToPi0Ratiopp7TeV->Fit(fitEtaToPi0pp7TeVHighPt,"NRMEX0+","",3.5,15.);
	cout << WriteParameterToFile(fitEtaToPi0pp7TeVHighPt)<< endl;		

	TGraphAsymmErrors* graphCombPi0SpecCombpp7TeV=    		(TGraphAsymmErrors*)fileCombinedpp->Get("graphInvCrossSectionPi0Comb7TeVStatErr"); 
	cout << "Pi0 spec" << endl;
	graphCombPi0SpecCombpp7TeV->Print();
	TGraphAsymmErrors* graphCombEtaSpecCombpp7TeV=    		(TGraphAsymmErrors*)fileCombinedpp->Get("graphInvCrossSectionEtaComb7TeVStatErr"); 
	cout << "Eta spec" << endl;
	graphCombEtaSpecCombpp7TeV->Print();

	TGraphAsymmErrors* graphCombEtaSpecCombpp7TeVMod=    		(TGraphAsymmErrors*)graphCombEtaSpecCombpp7TeV->Clone("graphInvCrossSectionEtaComb7TeVStatErr"); 

	Double_t* xValuesEta = graphCombEtaSpecCombpp7TeVMod->GetX();
	Int_t counter = 0;
	while (xValuesEta[counter] < 4){
		cout <<xValuesEta[counter] << endl;
		counter++;
	}	
	
	cout << counter << endl;
	Int_t nPointsToBeRemoved = graphCombEtaSpecCombpp7TeVMod->GetN()-counter;
	cout << nPointsToBeRemoved << endl;
	while(nPointsToBeRemoved > 0){
		graphCombEtaSpecCombpp7TeVMod->RemovePoint(graphCombEtaSpecCombpp7TeVMod->GetN()-1);
		nPointsToBeRemoved--;
	}
	graphCombEtaSpecCombpp7TeVMod->Print();
	Double_t* xValuesPi0 = graphCombPi0SpecCombpp7TeV->GetX();
	counter = 0;
	while (xValuesPi0[counter] < 4){
		cout <<xValuesPi0[counter] << endl;
		counter++;
	}	
	cout << counter << endl;
	cout <<graphCombEtaSpecCombpp7TeVMod->GetN()<< endl;
	Double_t xValuePi0 = 0;
	Double_t yValuePi0 = 0;
	Double_t yErrorPi0Up = 0;
	Double_t yErrorPi0Down = 0;
	while (counter < graphCombPi0SpecCombpp7TeV->GetN()){
		graphCombPi0SpecCombpp7TeV->GetPoint(counter,xValuePi0,yValuePi0);
		Double_t yValuePi0Int = yValuePi0*fitEtaToPi0pp7TeVHighPt->GetParameter(0);
		yErrorPi0Down = graphCombPi0SpecCombpp7TeV->GetErrorYlow(counter)*fitEtaToPi0pp7TeVHighPt->GetParameter(0);
		yErrorPi0Up = graphCombPi0SpecCombpp7TeV->GetErrorYhigh(counter)*fitEtaToPi0pp7TeVHighPt->GetParameter(0);
		Int_t newEtaPoint = graphCombEtaSpecCombpp7TeVMod->GetN();
		graphCombEtaSpecCombpp7TeVMod->SetPoint(newEtaPoint,xValuePi0,yValuePi0Int);
		graphCombEtaSpecCombpp7TeVMod->SetPointError(newEtaPoint,graphCombPi0SpecCombpp7TeV->GetErrorXlow(counter),graphCombPi0SpecCombpp7TeV->GetErrorXhigh(counter),yErrorPi0Down,yErrorPi0Up);
		counter++;
	}	
	graphCombEtaSpecCombpp7TeVMod->Print();
	
	TFile* fileOmegapp7TeV = 					new TFile("ExternalInput/PHOS/7TeV/PHOS_pp_omega_7TeV_07082012.root");
	TGraphErrors* graphOmegaToPi0Ratiopp7TeV =         (TGraphErrors*)fileOmegapp7TeV->Get("omega_to_pi0_stat_err");
	TGraphErrors* graphOmegaToPi0RatioSysErrpp7TeV=    (TGraphErrors*)graphOmegaToPi0Ratiopp7TeV->Clone("omega_to_pi0_syst_err"); 
	TGraphAsymmErrors* graphOmegaSpecStatpp7TeV = 	(TGraphAsymmErrors*)fileOmegapp7TeV->Get("graphOmegaStat");
	TGraphAsymmErrors* graphOmegaSpecSyspp7TeV = 		(TGraphAsymmErrors*)fileOmegapp7TeV->Get("graphOmegaSyst");
	TGraphAsymmErrors* graphOmegaSpecCombpp7TeV = 	(TGraphAsymmErrors*)fileOmegapp7TeV->Get("graphOmegaStat");//graphOmegaComb
	cout << "Omega spec" << endl;	
	graphOmegaSpecCombpp7TeV->Print();
	
	graphOmegaToPi0RatioSysErrpp7TeV->SetPointError(0,0.2,0.601407*0.155304);
	graphOmegaToPi0RatioSysErrpp7TeV->SetPointError(1,0.2,0.756806*0.0892788);
	graphOmegaToPi0RatioSysErrpp7TeV->SetPointError(2,0.2,0.802541*0.0939079);
	graphOmegaToPi0RatioSysErrpp7TeV->SetPointError(3,0.2,0.854334*0.168217);
	graphOmegaToPi0RatioSysErrpp7TeV->SetPointError(4,0.2,0.95354*0.133078);
	graphOmegaToPi0RatioSysErrpp7TeV->SetPointError(5,0.2,1.09497*0.175439);
	graphOmegaToPi0RatioSysErrpp7TeV->SetPointError(6,0.2,0.8821*0.298216);
	
	Double_t paramGraph[3] = {1.0e12,7.,0.13};
	TF1* fitInvCrossSectionPi07TeV = FitObject("l","fitInvCrossSectionPi07TeV","Pi0",graphCombPi0SpecCombpp7TeV,0.3,20 ,paramGraph,"QNRMEX0+");
	cout << "Pi0 7TeV ____________________________________________" << endl;
	cout << WriteParameterToFile(fitInvCrossSectionPi07TeV)<< endl;		
	cout << "dN/dy: \t" << fitInvCrossSectionPi07TeV->GetParameter(0)/xSection7TeVppINEL << endl << endl;
	
	
	TF1* fitInvCrossSectionOmega7TeV = FitObject("l","fitInvCrossSectionOmega7TeV","Omega",graphOmegaSpecCombpp7TeV,2.,17 ,paramGraph,"QNRMEX0+");
//  fitInvCrossSectionOmega7TeV->SetParameter(1,fitInvCrossSectionPi07TeV->GetParameter(1));
// 	fitInvCrossSectionOmega7TeV->SetParLimits(1,fitInvCrossSectionPi07TeV->GetParameter(1)*0.97,fitInvCrossSectionPi07TeV->GetParameter(1)*1.03);
// 	graphOmegaSpecCombpp7TeV->Fit(fitInvCrossSectionOmega7TeV,"QNRMEX0+","",2,17.);
	cout << "Omega 7TeV ____________________________________________" << endl;
	cout << WriteParameterToFile(fitInvCrossSectionOmega7TeV)<< endl;		
	cout << "dN/dy: \t" << fitInvCrossSectionOmega7TeV->GetParameter(0)/xSection7TeVppINEL << endl << endl;
	
	paramGraph[1] = 5.; 
	TF1* fitInvCrossSectionEta7TeV = FitObject("l","fitInvCrossSectionEta7TeV","Pi0",graphCombEtaSpecCombpp7TeVMod,0.4,15 ,paramGraph,"QNRMEX0+");
//  	fitInvCrossSectionEta7TeV->FixParameter(1,fitInvCrossSectionPi07TeV->GetParameter(1));
// // 	fitInvCrossSectionEta7TeV->SetParLimits(1,fitInvCrossSectionPi07TeV->GetParameter(1)*0.97,fitInvCrossSectionPi07TeV->GetParameter(1)*1.03);
	graphCombEtaSpecCombpp7TeVMod->Fit(fitInvCrossSectionEta7TeV,"QNRMEX0+","",0.4,15.);
	cout << "Eta 7TeV ____________________________________________" << endl;
	cout << WriteParameterToFile(fitInvCrossSectionEta7TeV)<< endl;		
	cout << "dN/dy: \t" << fitInvCrossSectionEta7TeV->GetParameter(0)/xSection7TeVppINEL << endl << endl;
	
	TGraphAsymmErrors* graphCombPhiSpecCombpp7TeV =get_phipp7TeV();
	graphCombPhiSpecCombpp7TeV->Print();

	TF1* fitInvCrossSectionPhi7TeV = FitObject("l","fitInvCrossSectionPhi7TeV","Phi",graphCombPhiSpecCombpp7TeV,0.45,6 ,paramGraph,"QNRMEX0+");
//  	fitInvCrossSectionPhi7TeV->SetParameter(1,fitInvCrossSectionPi07TeV->GetParameter(1));
// 	fitInvCrossSectionPhi7TeV->SetParLimits(1,fitInvCrossSectionPi07TeV->GetParameter(1)*0.97,fitInvCrossSectionPi07TeV->GetParameter(1)*1.03);
// 	graphCombPhiSpecCombpp7TeV->Fit(fitInvCrossSectionPhi7TeV,"QNRMEX0+","",2,17.);
	cout << "Phi 7TeV ____________________________________________" << endl;
	cout << WriteParameterToFile(fitInvCrossSectionPhi7TeV)<< endl;		
	cout << "dN/dy: \t" << fitInvCrossSectionPhi7TeV->GetParameter(0)/xSection7TeVppINEL << endl << endl;
	
	TGraphAsymmErrors* graphCombPhiPi0RatioCombpp7TeV = get_phiToPi0pp7TeV(fitInvCrossSectionPi07TeV);
	graphCombPhiPi0RatioCombpp7TeV->Print();
	
	
	TF1 *fitEtaToPi0pp7TeV = new TF1("fitEtaToPi0pp7TeV","([0]*([1]-1.)*([1]-2.)/ ([1]*[2]*([1]*[2]+[3]*([1]-2.))) * pow(1.+(sqrt(x*x+[3]*[3])-[3])/([1]*[2]), -[1]))/ ([4]*([5]-1.)*([5]-2.) / ([5]*[6]*([5]*[6]+[7]*([5]-2.)))  * pow(1.+(sqrt(x*x+[7]*[7])-[7])/([5]*[6]), -[5]))");
	fitEtaToPi0pp7TeV->SetParameter(0,fitInvCrossSectionEta7TeV->GetParameter(0));
	fitEtaToPi0pp7TeV->SetParameter(1,fitInvCrossSectionEta7TeV->GetParameter(1));
	fitEtaToPi0pp7TeV->SetParameter(2,fitInvCrossSectionEta7TeV->GetParameter(2));
	fitEtaToPi0pp7TeV->SetParameter(3,TDatabasePDG::Instance()->GetParticle(221)->Mass());
	fitEtaToPi0pp7TeV->SetParameter(4,fitInvCrossSectionPi07TeV->GetParameter(0));
	fitEtaToPi0pp7TeV->SetParameter(5,fitInvCrossSectionPi07TeV->GetParameter(1));
	fitEtaToPi0pp7TeV->SetParameter(6,fitInvCrossSectionPi07TeV->GetParameter(2));
	fitEtaToPi0pp7TeV->SetParameter(7,TDatabasePDG::Instance()->GetParticle(111)->Mass());

	TF1 *fitOmegaToPi0pp7TeV = new TF1("fitOmegaToPi0pp7TeV","([0] / ([1]*[2]*([1]*[2]+[3]*([1]-2.)))  * pow(1.+(sqrt(x*x+[3]*[3])-[3])/([1]*[2]), -[1]))/([4] *([5]-1.)*([5]-2.) / ([5]*[6]*([5]*[6]+[7]*([5]-2.)))  * pow(1.+(sqrt(x*x+[7]*[7])-[7])/([5]*[6]), -[5]))");
	fitOmegaToPi0pp7TeV->SetParameter(0,fitInvCrossSectionOmega7TeV->GetParameter(0));
	fitOmegaToPi0pp7TeV->SetParameter(1,fitInvCrossSectionOmega7TeV->GetParameter(1));
	fitOmegaToPi0pp7TeV->SetParameter(2,fitInvCrossSectionOmega7TeV->GetParameter(2));
	fitOmegaToPi0pp7TeV->SetParameter(3,TDatabasePDG::Instance()->GetParticle(223)->Mass());
	fitOmegaToPi0pp7TeV->SetParameter(4,fitInvCrossSectionPi07TeV->GetParameter(0));
	fitOmegaToPi0pp7TeV->SetParameter(5,fitInvCrossSectionPi07TeV->GetParameter(1));
	fitOmegaToPi0pp7TeV->SetParameter(6,fitInvCrossSectionPi07TeV->GetParameter(2));
	fitOmegaToPi0pp7TeV->SetParameter(7,TDatabasePDG::Instance()->GetParticle(111)->Mass());

	TF1 *fitPhiDivPi0pp7TeV = new TF1("fitPhiDivPi0pp7TeV","([0] *([1]-1.)*([1]-2.) / ([1]*[2]*([1]*[2]+[3]*([1]-2.)))  * pow(1.+(sqrt(x*x+[3]*[3])-[3])/([1]*[2]), -[1]))/([4] *([5]-1.)*([5]-2.) / ([5]*[6]*([5]*[6]+[7]*([5]-2.)))  * pow(1.+(sqrt(x*x+[7]*[7])-[7])/([5]*[6]), -[5]))");
	fitPhiDivPi0pp7TeV->SetParameter(0,fitInvCrossSectionPhi7TeV->GetParameter(0));
	fitPhiDivPi0pp7TeV->SetParameter(1,fitInvCrossSectionPhi7TeV->GetParameter(1));
	fitPhiDivPi0pp7TeV->SetParameter(2,fitInvCrossSectionPhi7TeV->GetParameter(2));
	fitPhiDivPi0pp7TeV->SetParameter(3,TDatabasePDG::Instance()->GetParticle(333)->Mass());
	fitPhiDivPi0pp7TeV->SetParameter(4,fitInvCrossSectionPi07TeV->GetParameter(0));
	fitPhiDivPi0pp7TeV->SetParameter(5,fitInvCrossSectionPi07TeV->GetParameter(1));
	fitPhiDivPi0pp7TeV->SetParameter(6,fitInvCrossSectionPi07TeV->GetParameter(2));
	fitPhiDivPi0pp7TeV->SetParameter(7,TDatabasePDG::Instance()->GetParticle(111)->Mass());

	TF1 *fitEtaToPi0pp7TeVFitted = new TF1("fitEtaToPi0pp7TeVFitted","([0] *([1]-1.)*([1]-2.) / ([1]*[2]*([1]*[2]+[3]*([1]-2.)))  * pow(1.+(sqrt(x*x+[3]*[3])-[3])/([1]*[2]), -[1]))/ ([4] *([1]-1.)*([1]-2.) / ([1]*[5]*([1]*[5]+[6]*([1]-2.)))  * pow(1.+(sqrt(x*x+[6]*[6])-[6])/([1]*[5]), -[1]))");
	fitEtaToPi0pp7TeVFitted->SetParameter(0,fitInvCrossSectionEta7TeV->GetParameter(0));
	fitEtaToPi0pp7TeVFitted->SetParameter(1,fitInvCrossSectionEta7TeV->GetParameter(1));
	fitEtaToPi0pp7TeVFitted->SetParameter(2,fitInvCrossSectionEta7TeV->GetParameter(2));
	fitEtaToPi0pp7TeVFitted->FixParameter(3,TDatabasePDG::Instance()->GetParticle(221)->Mass());
	fitEtaToPi0pp7TeVFitted->SetParameter(4,fitInvCrossSectionPi07TeV->GetParameter(0));
// 	fitEtaToPi0pp7TeVFitted->SetParameter(5,fitInvCrossSectionPi07TeV->GetParameter(1));
	fitEtaToPi0pp7TeVFitted->SetParameter(5,fitInvCrossSectionPi07TeV->GetParameter(2));
	fitEtaToPi0pp7TeVFitted->FixParameter(6,TDatabasePDG::Instance()->GetParticle(111)->Mass());
	graphCombEtaToPi0Ratiopp7TeV->Fit(fitEtaToPi0pp7TeVFitted,"NRMEX0+","",0.4,12.);	
	cout << WriteParameterToFile(fitEtaToPi0pp7TeVFitted)<< endl;		
	
	TF1 *fitOmegaToPi0pp7TeVFitted = new TF1("fitOmegaToPi0pp7TeVFitted","([0] / ( 2 * TMath::Pi())*([1]-1.)*([1]-2.) / ([1]*[2]*([1]*[2]+[3]*([1]-2.)))  * pow(1.+(sqrt(x*x+[3]*[3])-[3])/([1]*[2]), -[1]))/ ([4] *([1]-1.)*([1]-2.) / ([1]*[5]*([1]*[5]+[6]*([1]-2.)))  * pow(1.+(sqrt(x*x+[6]*[6])-[6])/([1]*[5]), -[1]))");
	fitOmegaToPi0pp7TeVFitted->SetParameter(0,fitInvCrossSectionEta7TeV->GetParameter(0));
	fitOmegaToPi0pp7TeVFitted->SetParameter(1,fitInvCrossSectionEta7TeV->GetParameter(1));
	fitOmegaToPi0pp7TeVFitted->SetParameter(2,fitInvCrossSectionEta7TeV->GetParameter(2));
	fitOmegaToPi0pp7TeVFitted->FixParameter(3,TDatabasePDG::Instance()->GetParticle(223)->Mass());
	fitOmegaToPi0pp7TeVFitted->SetParameter(4,fitInvCrossSectionPi07TeV->GetParameter(0));
// 	fitOmegaToPi0pp7TeVFitted->SetParameter(5,fitInvCrossSectionPi07TeV->GetParameter(1));
	fitOmegaToPi0pp7TeVFitted->SetParameter(5,fitInvCrossSectionPi07TeV->GetParameter(2));
	fitOmegaToPi0pp7TeVFitted->FixParameter(6,TDatabasePDG::Instance()->GetParticle(111)->Mass());
	graphOmegaToPi0Ratiopp7TeV->Fit(fitOmegaToPi0pp7TeVFitted,"NRMEX0+","",0.4,12.);	
	cout << WriteParameterToFile(fitOmegaToPi0pp7TeVFitted)<< endl;		
	
	
	TH1D* histoSpecEtaConstructed7TeV = new TH1D("histoSpecEtaConstructed7TeV","histoSpecEtaConstructed7TeV",100,0,40);
	for (Int_t i = 1; i < histoSpecEtaConstructed7TeV->GetNbinsX()+1;i++){
		histoSpecEtaConstructed7TeV->SetBinContent(i, fitEtaToPi0pp7TeVFitted->Eval(histoSpecEtaConstructed7TeV->GetBinCenter(i))*fitInvCrossSectionPi07TeV->Eval(histoSpecEtaConstructed7TeV->GetBinCenter(i)));
	}	

	TH1D* histoSpecOmegaConstructed7TeV = new TH1D("histoSpecOmegaConstructed7TeV","histoSpecEtaConstructed7TeV",100,0,40);
	for (Int_t i = 1; i < histoSpecOmegaConstructed7TeV->GetNbinsX()+1;i++){
		histoSpecOmegaConstructed7TeV->SetBinContent(i, fitOmegaToPi0pp7TeVFitted->Eval(histoSpecOmegaConstructed7TeV->GetBinCenter(i))*fitInvCrossSectionPi07TeV->Eval(histoSpecOmegaConstructed7TeV->GetBinCenter(i)));
	}	
	
	TH1D* ratioOmegaToPhi7TeVConstructed = CalculateHistoRatioToFit(histoSpecOmegaConstructed7TeV,fitInvCrossSectionPhi7TeV);
	TH1D* ratioEtaToPhi7TeVConstructed = CalculateHistoRatioToFit(histoSpecEtaConstructed7TeV,fitInvCrossSectionPhi7TeV);
	
	
	
	cout << "PCM Eta for PbPb" << endl;
	TFile* filePCMPbPb = 					new TFile(nameFilePbPb); 
	TDirectory*	directoryPCMEtaPbPb0010 = 				(TDirectory*)filePCMPbPb->Get("Eta_PbPb_2.76TeV_0-10%"); 
	TH1D* histoPCMEtaToPi0RatioPbPb0010 =            (TH1D*)directoryPCMEtaPbPb0010->Get("EtatoPi0Ratio");
	TGraphAsymmErrors* graphPCMEtaToPi0RatioSysErrPbPb0010=    (TGraphAsymmErrors*)directoryPCMEtaPbPb0010->Get("EtatoPi0RatioSys"); 

	TDirectory*	directoryPCMEtaPbPb2050 = 				(TDirectory*)filePCMPbPb->Get("Eta_PbPb_2.76TeV_20-50%"); 
	TH1D* histoPCMEtaToPi0RatioPbPb2050 =            (TH1D*)directoryPCMEtaPbPb2050->Get("EtatoPi0Ratio");
	TGraphAsymmErrors* graphPCMEtaToPi0RatioSysErrPbPb2050=    (TGraphAsymmErrors*)directoryPCMEtaPbPb2050->Get("EtatoPi0RatioSys"); 
	
	cout << "eta to pi0 syst" << endl;
// 	graphPCMEtaToPi0RatioSysErrPbPb0010->RemovePoint(graphPCMEtaToPi0RatioSysErrPbPb0010->GetN()-1);
// 	graphPCMEtaToPi0RatioSysErrPbPb2050->RemovePoint(graphPCMEtaToPi0RatioSysErrPbPb2050->GetN()-1);

	graphPCMEtaToPi0RatioSysErrPbPb2050->Print();
	graphPCMEtaToPi0RatioSysErrPbPb0010->Print();

	
	TFile* cocktailFileK0sScaled = new TFile("CocktailInput/cocktail_PbPb_0020c_dmtsallis_RatioEta.root");
	TFile* cocktailFileK0sScaled0040 = new TFile("CocktailInput/cocktail_PbPb_0040c_dmtsallis_RatioEta.root");
	TFile* cocktailFileMtScaledEta = new TFile("CocktailInput/cocktail_PbPb_0020c_dmtsallis_MTEta.root");
	TFile* cocktailFile7TeV = new TFile("CocktailInput/cocktail_7TeV_pi0LevyetaLevy.root");
	TDirectory* cocktailDirK0sScaled = (TDirectoryFile*) cocktailFileK0sScaled->Get("cocktail_PbPb_0020c_dmtsallis_RatioEta");
	TDirectory* cocktailDirK0sScaled0040 = (TDirectoryFile*) cocktailFileK0sScaled0040->Get("cocktail_PbPb_0040c_dmtsallis_RatioEta");
	TDirectory* cocktailDirMtScaledEta = (TDirectoryFile*) cocktailFileMtScaledEta->Get("cocktail_PbPb_0020c_dmtsallis_MTEta");
	TDirectory* cocktailDir7TeV = (TDirectoryFile*) cocktailFile7TeV->Get("cocktail_7TeV_pi0LevyetaLevy");
	
	Int_t rebin = 9;
	Double_t fBinsPi07TeVPt[22] =					{0.0, 0.1, 0.2, 0.3, 0.4,
													 0.5, 0.6, 0.8, 1.0, 1.5,  
												     1.9, 2.2, 2.7, 3.2,
													 4.0,  5.0, 6.0, 8.0,
													 10.0, 12.0, 18.0, 22.0};

	
	TH1D* cocktailPi0_K0Scaled = (TH1D* )cocktailDirK0sScaled->Get("ptPi0");
	cocktailPi0_K0Scaled->Sumw2();
	TH1D* cocktailPi0_K0ScaledRebinned = (TH1D* )cocktailPi0_K0Scaled->Rebin(21,"ptPionK0sScaeldRebinned",fBinsPi07TeVPt);
	TH1D* cocktailEta_K0Scaled = (TH1D* )cocktailDirK0sScaled->Get("ptEta");
	cocktailEta_K0Scaled->Sumw2();
	TH1D* cocktailEta_K0ScaledRebinned = (TH1D* )cocktailEta_K0Scaled->Rebin(21,"ptEtaK0sScaledRebinned",fBinsPi07TeVPt);
	TH1D* cocktailEtaToPi0Ratio_K0Scaled = (TH1D* )cocktailEta_K0Scaled->Clone("EtaToPi0Ratio_K0Scaled");
	cocktailEtaToPi0Ratio_K0Scaled->Sumw2();
	cocktailEtaToPi0Ratio_K0Scaled->Divide(cocktailEtaToPi0Ratio_K0Scaled,cocktailPi0_K0Scaled);
	TH1D* cocktailEtaToPi0Ratio_K0ScaledRebinned = (TH1D* )cocktailEta_K0ScaledRebinned->Clone("EtaToPi0Ratio_K0ScaledRebinned");
	cocktailEtaToPi0Ratio_K0ScaledRebinned->Sumw2();
	cocktailEtaToPi0Ratio_K0ScaledRebinned->Divide(cocktailEtaToPi0Ratio_K0ScaledRebinned,cocktailPi0_K0ScaledRebinned);
	
	TH1D* cocktailPi0_MtScaled = (TH1D* )cocktailDirMtScaledEta->Get("ptPi0");
	cocktailPi0_MtScaled->Sumw2();
	TH1D* cocktailPi0_MtScaledRebinned = (TH1D* )cocktailPi0_MtScaled->Rebin(21,"ptPionMTScaledRebinned",fBinsPi07TeVPt);
	TH1D* cocktailEta_MtScaled = (TH1D* )cocktailDirMtScaledEta->Get("ptEta");
	cocktailEta_MtScaled->Sumw2();
	TH1D* cocktailEta_MtScaledRebinned = (TH1D* )cocktailEta_MtScaled->Rebin(21,"ptEtaMTScaledRebinned",fBinsPi07TeVPt);
	TH1D* cocktailOmega_MtScaled = (TH1D* )cocktailDirMtScaledEta->Get("ptOmega");
	cocktailOmega_MtScaled->Sumw2();
	TH1D* cocktailOmega_MtScaledRebinned = (TH1D* )cocktailOmega_MtScaled->Rebin(21,"ptOmegaMTScaledRebinned",fBinsPi07TeVPt);

	TH1D* cocktailEtaToPi0Ratio_MtScaled = (TH1D* )cocktailEta_MtScaled->Clone("EtaToPi0Ratio_MtScaled");
	cocktailEtaToPi0Ratio_MtScaled->Sumw2();
	cocktailEtaToPi0Ratio_MtScaled->Divide(cocktailEtaToPi0Ratio_MtScaled,cocktailPi0_MtScaled);
	TH1D* cocktailEtaToPi0Ratio_MtScaledRebinned = (TH1D* )cocktailEta_MtScaledRebinned->Clone("EtaToPi0Ratio_MtScaledRebinned");
	cocktailEtaToPi0Ratio_MtScaledRebinned->Sumw2();
	cocktailEtaToPi0Ratio_MtScaledRebinned->Divide(cocktailEtaToPi0Ratio_MtScaledRebinned,cocktailPi0_MtScaledRebinned);
	
	TH1D* cocktailOmegaToPi0Ratio_MtScaledRebinned = (TH1D* )cocktailOmega_MtScaledRebinned->Clone("OmegaToPi0Ratio_MtScaledRebinned");
	cocktailOmegaToPi0Ratio_MtScaledRebinned->Sumw2();
	cocktailOmegaToPi0Ratio_MtScaledRebinned->Divide(cocktailOmegaToPi0Ratio_MtScaledRebinned,cocktailPi0_MtScaledRebinned);

	
	TH1D* cocktailPi07TeV = (TH1D* )cocktailDir7TeV->Get("ptPion");
	cocktailPi07TeV->Sumw2();
	TH1D* cocktailPi07TeVRebined = (TH1D* )cocktailPi07TeV->Rebin(21,"ptPionRebinned",fBinsPi07TeVPt);
	cocktailPi07TeV->Rebin(rebin);
	cout << "here"<< endl;
	TH1D* cocktailEta7TeV = (TH1D* )cocktailDir7TeV->Get("ptEta");
	cocktailEta7TeV->Sumw2();
	TH1D* cocktailEta7TeVRebined = (TH1D* )cocktailEta7TeV->Rebin(21,"ptEtaRebinned",fBinsPi07TeVPt);
	cocktailEta7TeV->Rebin(rebin);
	cout << "here"<< endl;
	TH1D* cocktailOmega7TeV = (TH1D* )cocktailDir7TeV->Get("ptOmega");
	cocktailOmega7TeV->Sumw2();
	TH1D* cocktailOmega7TeVRebinned = (TH1D* )cocktailOmega7TeV->Rebin(21,"ptOmegaRebinned",fBinsPi07TeVPt);
	cocktailOmega7TeV->Rebin(rebin);
	
	TH1D* cocktailPhi = (TH1D* )cocktailDirK0sScaled0040->Get("ptPhi");
	cocktailPhi->Sumw2();
	TH1D* cocktailPhiRebinned = (TH1D* )cocktailPhi->Rebin(21,"ptPhiRebinned",fBinsPi07TeVPt);
// 	cocktailPhi->Rebin(rebin);
	TH1D* cocktailOmegaPbPb0040 = (TH1D* )cocktailDirK0sScaled0040->Get("ptOmega");
	cocktailOmegaPbPb0040->Sumw2();
	TH1D* cocktailOmegaPbPb0040Rebinned = (TH1D* )cocktailOmegaPbPb0040->Rebin(21,"ptOmegaRebinned",fBinsPi07TeVPt);

	TH1D* cocktailPi0PbPb0040 = (TH1D* )cocktailDirK0sScaled0040->Get("ptPi0");
	cocktailPi0PbPb0040->Sumw2();
	TH1D* cocktailPi0PbPb0040Rebinned_2 = (TH1D* )cocktailPi0PbPb0040->Rebin(21,"ptPi0Rebinned_2",fBinsPi07TeVPt);
// 	cocktailPi0PbPb0040->Rebin(rebin);

	
	TH1D* cocktailEtaToPi0Ratio7TeV = (TH1D* )cocktailEta7TeV->Clone("EtaToPi0Ratio7TeV");
	cocktailEtaToPi0Ratio7TeV->Sumw2();
	cocktailEtaToPi0Ratio7TeV->Divide(cocktailEtaToPi0Ratio7TeV,cocktailPi07TeV);
	
	TH1D* cocktailEtaToPi0Ratio7TeVRebined = (TH1D* )cocktailEta7TeVRebined->Clone("EtaToPi0Ratio7TeVRebined");
	cocktailEtaToPi0Ratio7TeVRebined->Sumw2();
	cocktailEtaToPi0Ratio7TeVRebined->Divide(cocktailEtaToPi0Ratio7TeVRebined,cocktailPi07TeVRebined);
	
	TH1D* cocktailOmegaToPi0Ratio7TeVRebined = (TH1D* )cocktailOmega7TeVRebinned->Clone("OmegaToPi0Ratio7TeVRebined");
	cocktailOmegaToPi0Ratio7TeVRebined->Sumw2();
	cocktailOmegaToPi0Ratio7TeVRebined->Divide(cocktailOmegaToPi0Ratio7TeVRebined,cocktailPi07TeVRebined);
	
	TH1D* cocktailPhiPi0RatioPbPb0040Rebined = (TH1D* )cocktailPhiRebinned->Clone("PhiPi0RatioPbPb0040Rebined");
	cocktailPhiPi0RatioPbPb0040Rebined->Sumw2();
	cocktailPhiPi0RatioPbPb0040Rebined->Divide(cocktailPhiPi0RatioPbPb0040Rebined,cocktailPi0PbPb0040Rebinned_2);
	TH1D* cocktailOmegaToPi0RatioPbPb0040Rebined = (TH1D* )cocktailOmegaPbPb0040Rebinned->Clone("OmegaToPi0RatioPbPb0040Rebined");
	cocktailOmegaToPi0RatioPbPb0040Rebined->Sumw2();
	cocktailOmegaToPi0RatioPbPb0040Rebined->Divide(cocktailOmegaToPi0RatioPbPb0040Rebined,cocktailPi0PbPb0040Rebinned_2);
	
	
	Width_t  widthLinesBoxes            = 1.4;
	Width_t  widthCommonFit             = 2.;
	Width_t  widthStatErrBars           = 1.5;
	Width_t  widthCommonErrors          = 1.1;
	Width_t  widthCommonSpectrumBoxes         = 0.99;
	if (suffix.CompareTo("eps")==0){
		widthLinesBoxes            = 1.4;
		widthCommonFit             = 2.;
		widthStatErrBars           = 1.5;
		widthCommonErrors          = 1.1;
		widthCommonSpectrumBoxes         = 0.99;
	} else {
		widthLinesBoxes            = 2.3;
		widthCommonFit             = 2.6;
		widthStatErrBars           = 2.6;
		widthCommonErrors          = 2.;
		widthCommonSpectrumBoxes         = 2.3;
	}

	//*********************** Comparison of Measuremnts in ALICE ***************************************************
	TCanvas* canvasSpec7TeV = new TCanvas("canvasSpec7TeV","",200,10,1350,1500);  // gives the page size
	DrawGammaCanvasSettings( canvasSpec7TeV, 0.15, 0.01, 0.015, 0.115);
	canvasSpec7TeV->SetLogy(1);
	
	TH2D *histo2DSpec7TeV = new TH2D("histo2DSpec7TeV", "histo2DSpec7TeV", 20,0.,21.,1000.,1e2,1e12);
	SetStyleHistoTH2ForGraphs(histo2DSpec7TeV, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )", 0.032,0.04, 0.032,0.04, 1,1.55);
// 	histo2DSpec7TeV->GetYaxis()->SetRangeUser(0.,1.02);
	histo2DSpec7TeV->Draw();
	
	DrawGammaSetMarkerTGraphAsym(graphOmegaSpecCombpp7TeV, markerStylepPb, markerSizepPb, colorOmega, colorOmega, widthLinesBoxes, kTRUE);
	graphOmegaSpecCombpp7TeV->Draw("same,pe");
	fitInvCrossSectionOmega7TeV->SetRange(0,21.);
	DrawGammaSetMarkerTF1(fitInvCrossSectionOmega7TeV, 1, widthCommonFit, colorOmega);
// 	histoSpecOmegaConstructed7TeV->Draw("same,hist");
	fitInvCrossSectionOmega7TeV->Draw("same");
	
	DrawGammaSetMarkerTGraphAsym(graphCombPi0SpecCombpp7TeV, markerStylepPb, markerSizepPb, colorPi0, colorPi0, widthLinesBoxes, kTRUE);
	graphCombPi0SpecCombpp7TeV->Draw("same,pe");
	fitInvCrossSectionPi07TeV->SetRange(0,21.);
	DrawGammaSetMarkerTF1(fitInvCrossSectionPi07TeV, 1, widthCommonFit, colorPi0);
	fitInvCrossSectionPi07TeV->Draw("same");
	DrawGammaSetMarkerTGraphAsym(graphCombEtaSpecCombpp7TeV, markerStylepPb, markerSizepPb, colorEta, colorEta, widthLinesBoxes, kTRUE);
	graphCombEtaSpecCombpp7TeV->Draw("same,pe");
// 	histoSpecEtaConstructed7TeV->Draw("same,hist");
	DrawGammaSetMarkerTGraphAsym(graphCombEtaSpecCombpp7TeVMod, markerStylepPb+4, markerSizepPb, colorEta, colorEta, widthLinesBoxes, kTRUE);
	graphCombEtaSpecCombpp7TeVMod->Draw("same,pe");
	
	
	fitInvCrossSectionEta7TeV->SetRange(0,21.);
	DrawGammaSetMarkerTF1(fitInvCrossSectionEta7TeV, 1, widthCommonFit, colorEta);
	fitInvCrossSectionEta7TeV->Draw("same");
	
	DrawGammaSetMarkerTGraphAsym(graphCombPhiSpecCombpp7TeV, markerStylepPb, markerSizepPb, colorPhi, colorPhi, widthLinesBoxes, kTRUE);
	graphCombPhiSpecCombpp7TeV->Draw("same,pe");
	fitInvCrossSectionPhi7TeV->SetRange(0,21.);
	DrawGammaSetMarkerTF1(fitInvCrossSectionPhi7TeV, 1, widthCommonFit, colorPhi);
	fitInvCrossSectionPhi7TeV->Draw("same");
	

	TLegend* legendSpec7TeV = new TLegend(0.65,0.7,0.977,0.93);
	legendSpec7TeV->SetTextSize(0.04);			
	legendSpec7TeV->SetFillColor(0);
	legendSpec7TeV->SetFillStyle(0);
	legendSpec7TeV->SetBorderSize(0);
	legendSpec7TeV->SetNColumns(2);
	legendSpec7TeV->SetMargin(0.2);
	legendSpec7TeV->SetHeader(collisionSystemPP7TeV.Data());
	legendSpec7TeV->AddEntry(graphCombPi0SpecCombpp7TeV,"#pi^{0}","pe");
	legendSpec7TeV->AddEntry(fitInvCrossSectionPi07TeV,"#pi^{0} fit","l");
	legendSpec7TeV->AddEntry(graphCombEtaSpecCombpp7TeV,"#eta","pe");
	legendSpec7TeV->AddEntry(fitInvCrossSectionEta7TeV,"#eta fit","l");
// 	legendSpec7TeV->AddEntry(graphCombEtaSpecCombpp7TeVMod,"#eta modified","pe");
// 	legendSpec7TeV->AddEntry((TObject*)0,"","");
	legendSpec7TeV->AddEntry(graphOmegaSpecCombpp7TeV,"#omega","pe");
	legendSpec7TeV->AddEntry(fitInvCrossSectionOmega7TeV,"#omega fit","l");
	legendSpec7TeV->AddEntry(graphCombPhiSpecCombpp7TeV,"#phi","pe");
	legendSpec7TeV->AddEntry(fitInvCrossSectionPhi7TeV,"#phi fit","l");
	legendSpec7TeV->Draw();

	
	canvasSpec7TeV->Update();
// 	canvasSpec7TeV->SaveAs(Form("%s/SpectraPlusFits_pp7TeV.%s",outputDir.Data(), suffix.Data()));


	TCanvas* canvasRatioToFitsALICE = new TCanvas("canvasRatioToFitsALICE","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings( canvasRatioToFitsALICE, 0.09, 0.01, 0.015, 0.115);
	
	TH2D *histo2DRatioFitToDatapp7TeV = new TH2D("histo2DRatioFitToDatapp7TeV", "histo2DRatioFitToDatapp7TeV", 20,0.,20.,1000.,-0.4,2.2);
	SetStyleHistoTH2ForGraphs(histo2DRatioFitToDatapp7TeV, "#it{p}_{T} (GeV/#it{c})","#eta/#pi^{0}", 0.046,0.058, 0.046,0.058, 0.85,0.74, 510, 510);
	histo2DRatioFitToDatapp7TeV->GetYaxis()->SetRangeUser(0.45,1.45);
	histo2DRatioFitToDatapp7TeV->Draw();
	
	TGraphAsymmErrors* graphRatioPi0ToFitpp7TeV = CalculateGraphErrRatioToFit(graphCombPi0SpecCombpp7TeV,fitInvCrossSectionPi07TeV);
	DrawGammaSetMarkerTGraphAsym(graphRatioPi0ToFitpp7TeV,  markerStylepPb, markerSizepPb, colorPi0, colorPi0, widthLinesBoxes, kTRUE);
	graphRatioPi0ToFitpp7TeV->Draw("same,pe");

	TGraphAsymmErrors* graphRatioEtaToFitpp7TeV = CalculateGraphErrRatioToFit(graphCombEtaSpecCombpp7TeV,fitInvCrossSectionEta7TeV);
	DrawGammaSetMarkerTGraphAsym(graphRatioEtaToFitpp7TeV,  markerStylepPb, markerSizepPb, colorEta, colorEta, widthLinesBoxes, kTRUE);
	graphRatioEtaToFitpp7TeV->Draw("same,pe");

	TGraphAsymmErrors* graphRatioEtaToFitpp7TeVmod = CalculateGraphErrRatioToFit(graphCombEtaSpecCombpp7TeVMod,fitInvCrossSectionEta7TeV);
	DrawGammaSetMarkerTGraphAsym(graphRatioEtaToFitpp7TeVmod,  markerStylepPb+4, markerSizepPb, colorEta, colorEta, widthLinesBoxes, kTRUE);
// 	graphRatioEtaToFitpp7TeVmod->Draw("same,pe");

	TGraphAsymmErrors* graphRatioOmegaToFitpp7TeV = CalculateGraphErrRatioToFit(graphOmegaSpecCombpp7TeV,fitInvCrossSectionOmega7TeV);
	DrawGammaSetMarkerTGraphAsym(graphRatioOmegaToFitpp7TeV,  markerStylepPb, markerSizepPb, colorOmega, colorOmega, widthLinesBoxes, kTRUE);
	graphRatioOmegaToFitpp7TeV->Draw("same,pe");

	TGraphAsymmErrors* graphRatioPhiToFitpp7TeV = CalculateGraphErrRatioToFit(graphCombPhiSpecCombpp7TeV,fitInvCrossSectionPhi7TeV);
	DrawGammaSetMarkerTGraphAsym(graphRatioPhiToFitpp7TeV,  markerStylepPb, markerSizepPb, colorPhi, colorPhi, widthLinesBoxes, kTRUE);
	graphRatioPhiToFitpp7TeV->Draw("same,pe");
// 	graphRatioPhiToFitpp7TeV->Print();
	
	DrawGammaLines(0, 20,1,1,1,kGray,4);
				
// 	TLegend* legendRatioALICE = new TLegend(0.3,0.16,0.977,0.34);
// 	legendRatioALICE->SetTextSize(0.04);			
// 	legendRatioALICE->SetFillColor(0);
// 	legendRatioALICE->SetFillStyle(0);
// 	legendRatioALICE->SetBorderSize(0);
// 	legendRatioALICE->SetMargin(0.1);
// 	legendRatioALICE->AddEntry(graphPCMEtaToPi0RatioSysErrpPb,collisionSystempPb.Data(),"pe");
// 	legendRatioALICE->AddEntry(graphCombEtaToPi0RatioSysErrpp7TeV,collisionSystemPP7TeV.Data(),"pe");
// 	legendRatioALICE->AddEntry(cocktailEtaToPi0Ratio_MtScaledRebinned,"p-Pb #sqrt{#it{s}_{_{NN}}} = 5.02 TeV, #eta from m_{T} scaled #pi^{0}","pl");
// 	legendRatioALICE->Draw();
	
	canvasRatioToFitsALICE->Update();
// 	canvasRatioToFitsALICE->SaveAs(Form("%s/RatioFitspp7TeV.%s",outputDir.Data(), suffix.Data()));
	
	
	
	//*********************** Comparison of Measuremnts in ALICE ***************************************************
	TCanvas* canvasRatioEtaToPi0ALICE = new TCanvas("canvasRatioEtaToPi0ALICE","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings( canvasRatioEtaToPi0ALICE, 0.09, 0.01, 0.015, 0.115);
	
	TH2D *histo2DRatioEtaToPi0ALICE = new TH2D("histo2DRatioEtaToPi0ALICE", "histo2DRatioEtaToPi0ALICE", 20,0.,20.,1000.,-0.4,1.2);
	SetStyleHistoTH2ForGraphs(histo2DRatioEtaToPi0ALICE, "#it{p}_{T} (GeV/#it{c})","#eta/#pi^{0}", 0.046,0.058, 0.046,0.058, 0.85,0.74, 510, 510);
	histo2DRatioEtaToPi0ALICE->GetYaxis()->SetRangeUser(0.,1.02);
	histo2DRatioEtaToPi0ALICE->Draw("copy");
	
	DrawGammaSetMarker(histoPCMEtaToPi0RatioPbPb0010, markerStylePbPb0010, markerSizepPb, colorCombPbPb0010, colorCombPbPb0010);
	DrawGammaSetMarker(histoPCMEtaToPi0RatioPbPb2050, markerStylePbPb2050, markerSizepPb, colorCombPbPb2050, colorCombPbPb2050);
	DrawGammaSetMarker(histoPCMEtaToPi0RatiopPb, markerStylepPb, markerSizepPb, colorCombpPb, colorCombpPb);
	
	DrawGammaSetMarkerTGraphAsym(graphPCMEtaToPi0RatioSysErrpPb, markerStylepPb, markerSizepPb, colorCombpPb, colorCombpPb, widthLinesBoxes, kTRUE);
	DrawGammaSetMarkerTGraphAsym(graphPCMEtaToPi0RatioSysErrPbPb0010, markerStylePbPb0010, markerSizepPb, colorCombPbPb0010, colorCombPbPb0010, widthLinesBoxes, kTRUE);
	DrawGammaSetMarkerTGraphAsym(graphPCMEtaToPi0RatioSysErrPbPb2050, markerStylePbPb2050, markerSizepPb, colorCombPbPb2050, colorCombPbPb2050, widthLinesBoxes, kTRUE);
	

	TGraphAsymmErrors* graphCombEtaToPi0Ratiopp7TeVNoXErrors = (TGraphAsymmErrors*)graphCombEtaToPi0Ratiopp7TeV->Clone();
	ProduceGraphAsymmFixedXErrors(graphCombEtaToPi0Ratiopp7TeVNoXErrors, 0.);
	DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0Ratiopp7TeVNoXErrors, markerStyleSpectrum7TeV, markerSizePi0PP7TeV, color7TeV, color7TeV);
	
	DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0RatioSysErrpp7TeV, markerStyleSpectrum7TeV, markerSizePi0PP7TeV, color7TeV, color7TeV, widthLinesBoxes, kTRUE);
	
	DrawGammaSetMarker(cocktailEtaToPi0Ratio7TeVRebined, 2, 0, kBlue-2, kBlue-2);
	DrawGammaSetMarker(cocktailEtaToPi0Ratio_K0ScaledRebinned, 2, 0, kRed+2, kRed+2);
	DrawGammaSetMarker(cocktailEtaToPi0Ratio_MtScaledRebinned, 2, 0, kGreen+2, kGreen+2);
// 	graphPCMEtaToPi0RatioSysErrpPb->Draw("same,pE2");  
	graphCombEtaToPi0RatioSysErrpp7TeV->Draw("same,pE2");
	
// 	histoPCMEtaToPi0RatiopPb->Draw("same,peX0");
	graphCombEtaToPi0Ratiopp7TeVNoXErrors->Draw("same,pe");
	cocktailEtaToPi0Ratio_MtScaledRebinned->GetXaxis()->SetRangeUser(0.2,20);
	cocktailEtaToPi0Ratio_MtScaledRebinned->Draw("same,hist,c");
	
	TLegend* legendRatioALICE = new TLegend(0.3,0.16,0.977,0.34);
	legendRatioALICE->SetTextSize(0.04);			
	legendRatioALICE->SetFillColor(0);
	legendRatioALICE->SetFillStyle(0);
	legendRatioALICE->SetBorderSize(0);
	legendRatioALICE->SetMargin(0.1);
// 	legendRatioALICE->AddEntry(graphPCMEtaToPi0RatioSysErrpPb,collisionSystempPb.Data(),"pe");
	legendRatioALICE->AddEntry(graphCombEtaToPi0RatioSysErrpp7TeV,collisionSystemPP7TeV.Data(),"pe");
	legendRatioALICE->AddEntry(cocktailEtaToPi0Ratio_MtScaledRebinned,"Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV, #eta from m_{T} scaled #pi^{0}","pl");
	legendRatioALICE->Draw();
	
	canvasRatioEtaToPi0ALICE->Update();
// 	canvasRatioEtaToPi0ALICE->SaveAs(Form("%s/EtaToPi0Ratio_pp7TeVwithmTscaled.%s",outputDir.Data(), suffix.Data()));

	
	TCanvas* canvasRatioEtaToPi0ALICEPbPb = new TCanvas("canvasRatioEtaToPi0ALICEPbPb","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings( canvasRatioEtaToPi0ALICEPbPb, 0.09, 0.01, 0.015, 0.115);
	
	TH2D *histo2DRatioEtaToPi0ALICEPbPb = new TH2D("histo2DRatioEtaToPi0ALICEPbPb", "histo2DRatioEtaToPi0ALICEPbPb", 20,0.1,15.01,1000.,-0.4,1.3);
	SetStyleHistoTH2ForGraphs(histo2DRatioEtaToPi0ALICEPbPb, "#it{p}_{T} (GeV/#it{c})","#eta/#pi^{0}", 0.046,0.058, 0.046,0.058, 0.85,0.74, 510, 510);
	histo2DRatioEtaToPi0ALICEPbPb->GetYaxis()->SetRangeUser(0.,1.32);
	histo2DRatioEtaToPi0ALICEPbPb->Draw("copy");

	graphCombEtaToPi0RatioSysErrpp7TeV->Draw("same,pE2");
	graphCombEtaToPi0Ratiopp7TeVNoXErrors->Draw("same,pe");
// 	cocktailEtaToPi0Ratio7TeVRebined->GetXaxis()->SetRangeUser(0.2,20);
// 	cocktailEtaToPi0Ratio_K0ScaledRebinned->GetXaxis()->SetRangeUser(0.2,20);
// 	cocktailEtaToPi0Ratio_MtScaledRebinned->GetXaxis()->SetRangeUser(0.2,20);
	cocktailEtaToPi0Ratio7TeVRebined->SetLineStyle(5);
	cocktailEtaToPi0Ratio7TeVRebined->SetLineWidth(2.5);
	cocktailEtaToPi0Ratio_K0ScaledRebinned->SetLineStyle(6);
	cocktailEtaToPi0Ratio_K0ScaledRebinned->SetLineWidth(2.5);
	cocktailEtaToPi0Ratio_MtScaledRebinned->SetLineStyle(7);
	cocktailEtaToPi0Ratio_MtScaledRebinned->SetLineWidth(2.5);
	cocktailEtaToPi0Ratio7TeVRebined->Draw("same,hist,c");
	cocktailEtaToPi0Ratio_K0ScaledRebinned->Draw("same,hist,c");
	cocktailEtaToPi0Ratio_MtScaledRebinned->Draw("same,hist,c");
	
	histoPCMEtaToPi0RatioPbPb2050->Draw("same,peX0");
	histoPCMEtaToPi0RatioPbPb0010->Draw("same,peX0");
	histoPCMEtaToPi0RatiopPb->Draw("same,peX0");
	graphPCMEtaToPi0RatioSysErrpPb->Draw("same,pE2");
	graphPCMEtaToPi0RatioSysErrPbPb2050->Draw("same,pE2");  
	graphPCMEtaToPi0RatioSysErrPbPb0010->Draw("same,pE2");

	
	TLegend* legendRatioALICEdata2 = new TLegend(0.12,0.7,0.32,0.96);
	legendRatioALICEdata2->SetTextSize(0.04);			
	legendRatioALICEdata2->SetFillColor(0);
	legendRatioALICEdata2->SetFillStyle(0);
	legendRatioALICEdata2->SetBorderSize(0);
// 	legendRatioALICEdata2->SetMargin(0.1);
	legendRatioALICEdata2->AddEntry(graphCombEtaToPi0RatioSysErrpp7TeV,Form("ALICE, %s",collisionSystemPP7TeV.Data()),"pe");
	legendRatioALICEdata2->AddEntry((TObject*)0,"Phys.Lett. B717 (2012) 162-172","");
	legendRatioALICEdata2->AddEntry(histoPCMEtaToPi0RatioPbPb0010,Form("%s (2011)",collisionSystemPbPb0010.Data()),"pe");
	legendRatioALICEdata2->AddEntry(histoPCMEtaToPi0RatioPbPb2050,Form("%s (2011)",collisionSystemPbPb2050.Data()),"pe");
	legendRatioALICEdata2->AddEntry(histoPCMEtaToPi0RatiopPb,Form("%s, LHC13b+c",collisionSystempPb.Data()),"pe");
	legendRatioALICEdata2->Draw();
	
	TLegend* legendRatioALICE2 = new TLegend(0.28,0.13,0.977,0.27);
	legendRatioALICE2->SetTextSize(0.04);			
	legendRatioALICE2->SetFillColor(0);
	legendRatioALICE2->SetFillStyle(0);
	legendRatioALICE2->SetBorderSize(0);
	legendRatioALICE2->SetMargin(0.1);
	legendRatioALICE2->AddEntry(cocktailEtaToPi0Ratio7TeVRebined,Form("#eta from %s as input",collisionSystemPP7TeV.Data()),"pl");
	legendRatioALICE2->AddEntry(cocktailEtaToPi0Ratio_MtScaledRebinned,"#eta from m_{T} scaled #pi^{0}","pl");
	legendRatioALICE2->AddEntry(cocktailEtaToPi0Ratio_K0ScaledRebinned,"0-20% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV (2010), #eta from K^{0}_{s} scaled ","pl");
	legendRatioALICE2->Draw();
	
	canvasRatioEtaToPi0ALICEPbPb->Update();  //qui
	canvasRatioEtaToPi0ALICEPbPb->SaveAs(Form("%s/EtaToPi0Ratio_DiffSystems_PbPb.%s",outputDir.Data(), suffix.Data()));
	
	if(checkStatFluc){
		TH1D* histoChi2 = (TH1D*)histoPCMEtaToPi0RatioPbPb0010->Clone("histoChi2");
		
		Int_t nA = histoPCMEtaToPi0RatioPbPb0010->GetNbinsX();
		Int_t nB = histoPCMEtaToPi0RatioPbPb2050->GetNbinsX();
		cout << "Number of bins for the Chi2 study: " << nA << endl;
		
		Double_t *Chi2 = new Double_t[nA];
		Double_t *pointA = new Double_t[nA];
		Double_t *pointB = new Double_t[nA];
		Double_t *sigmaChi2 = new Double_t[nA];
		Double_t *sigmaA = new Double_t[nA];
		Double_t *sigmaB = new Double_t[nA];
		Double_t SumChi2;
		if(nA = nB){
			for(Int_t i=2; i<nA; i++){
				pointA[i] = histoPCMEtaToPi0RatioPbPb0010->GetBinContent(i+1);
				pointB[i] = histoPCMEtaToPi0RatioPbPb2050->GetBinContent(i+1);
				sigmaA[i] = histoPCMEtaToPi0RatioPbPb0010->GetBinError(i+1);
				sigmaB[i] = histoPCMEtaToPi0RatioPbPb2050->GetBinError(i+1);
				sigmaChi2[i] = pow(sigmaA[i],2) + pow(sigmaB[i],2);
				
				cout << "For bin " << i << " pointA = " << pointA[i] << " +- " << sigmaA[i] <<  endl;
				cout << "For bin " << i << " pointB = " << pointB[i] << " +- " << sigmaB[i] << endl;
				cout << "For bin " << i << " sigmaChi2 = " << sigmaChi2[i] << endl;

				Chi2[i] = pow( (pointA[i] - pointB[i]), 2) / sigmaChi2[i];
				cout << "Chi2 for bin " << i << " = " << Chi2[i] << endl;
				histoChi2->SetBinContent(i+1,Chi2[i]);
// 				cout << "Chi2/ndf for bin " << i << " = " << Chi2[i]/10 << endl;
				
				SumChi2+=Chi2[i];
				cout << "Probability for bin " << i << " is " << TMath::Prob(SumChi2,i-1) << endl;
			}
		cout << "Probability for is " << TMath::Prob(SumChi2,10) << endl;
		cout << "The total chi2 is: " << SumChi2 << endl;
			
		} else {
			cout << "histos with different binning!" << endl;	
		}
		

		
		
		TCanvas* canvasTestChi2 = new TCanvas("canvasTestChi2","",200,10,1350,900);  // gives the page size
		DrawGammaCanvasSettings( canvasTestChi2, 0.09, 0.01, 0.015, 0.115);
	
// 		TH2D *histo2DTestChi2 = new TH2D("histo2DTestChi2", "histo2DTestChi2", 20,0.,20.,1000.,0.,2.);
// 		SetStyleHistoTH2ForGraphs(histo2DTestChi2, "#it{p}_{T} (GeV/#it{c})","#chi^{2}", 0.046,0.058, 0.046,0.058, 0.85,0.74, 510, 510);
// 		histo2DTestChi2->GetYaxis()->SetRangeUser(0.,20.);
// 		histo2DTestChi2->Draw("copy");
		
		histoChi2->GetXaxis()->SetRangeUser(1.,10.);
		histoChi2->GetYaxis()->SetRangeUser(0.,20.);
		histoChi2->GetYaxis()->SetTitle("#chi^{2}/ndf");
		histoChi2->SetLineColor(kBlue+2);
		histoChi2->SetLineWidth(2);
		histoChi2->Draw("same,hist");
		
		 
		canvasTestChi2->Update();  
		canvasTestChi2->SaveAs(Form("%s/Chi2forEtaToPi0Ratio_PbPb.%s",outputDir.Data(), suffix.Data()));

		TCanvas* c1 = new TCanvas("c1","",200,10,1350,900);  // gives the page size
		DrawGammaCanvasSettings( c1, 0.09, 0.01, 0.015, 0.115);
		
			
		TH2D *histo2DTestChi2 = new TH2D("histo2DTestChi2", "histo2DTestChi2", 20,0.,20.,1000.,0.,2.);
		SetStyleHistoTH2ForGraphs(histo2DTestChi2, "#chi^{2}", "", 0.046,0.058, 0.046,0.058, 0.85,0.74, 510, 510);
		histo2DTestChi2->GetXaxis()->SetRangeUser(0.,20.);
		histo2DTestChi2->GetYaxis()->SetRangeUser(0.,0.3);
		histo2DTestChi2->Draw("copy");
		
		TF1 *fc2 = new TF1("fc2","TMath::GammaDist(x,[0],[1],[2])",0,20.); 
		fc2->FixParameter(0,2.); 
		fc2->FixParameter(1,0.); 
		fc2->FixParameter(2,2.); 
		fc2->Draw("same"); 
		
		TF1 *pdf = new TF1("pdf", "ROOT::Math::chisquared_pdf(x, [0], [1])",0,20);
		pdf->SetParameters(3.5,0.1);
		pdf -> SetLineColor(2);
		pdf -> Draw("same"); 
		
		c1->Update(); 
		c1->SaveAs(Form("%s/c1.%s",outputDir.Data(), suffix.Data()));
		
		if(Chi2) delete Chi2;
		if(pointA) delete pointA;
		if(pointB) delete pointB;
		if(sigmaA) delete sigmaA;
		if(sigmaB) delete sigmaB;
		if(sigmaChi2) delete sigmaChi2;
	}
	
	

	
	
	histo2DRatioEtaToPi0ALICE->Draw("copy");
	graphCombEtaToPi0RatioSysErrpp7TeV->Draw("same,pE2");
	
	graphCombEtaToPi0Ratiopp7TeVNoXErrors->Draw("same,pe");
	DrawGammaSetMarkerTF1(fitEtaToPi0pp7TeVHighPt, 5, widthCommonFit, kBlack);
	fitEtaToPi0pp7TeVHighPt->SetRange(3.5,20);
	fitEtaToPi0pp7TeVHighPt->Draw("same");
	
	fitEtaToPi0pp7TeVFitted->SetRange(0,20);
	DrawGammaSetMarkerTF1(fitEtaToPi0pp7TeVFitted, 1, widthCommonFit, kBlack);
	fitEtaToPi0pp7TeVFitted->Draw("same");
	
	fitEtaToPi0pp7TeV->SetRange(0,20);
	DrawGammaSetMarkerTF1(fitEtaToPi0pp7TeV, 1, widthCommonFit, kGray);
	fitEtaToPi0pp7TeV->SetLineStyle(4);
// 	fitEtaToPi0pp7TeV->Draw("same");
	
	
	DrawGammaSetMarker(cocktailEtaToPi0Ratio7TeVRebined, 2, 0, kBlue-2, kBlue-2);
	DrawGammaSetMarker(cocktailEtaToPi0Ratio_K0ScaledRebinned, 2, 0, kRed+2, kRed+2);
	DrawGammaSetMarker(cocktailEtaToPi0Ratio_MtScaledRebinned, 2, 0, kGreen+2, kGreen+2);
	cocktailEtaToPi0Ratio7TeVRebined->GetXaxis()->SetRangeUser(0.2,20);
	cocktailEtaToPi0Ratio_K0ScaledRebinned->GetXaxis()->SetRangeUser(0.2,20);
	cocktailEtaToPi0Ratio_MtScaledRebinned->GetXaxis()->SetRangeUser(0.2,20);
	cocktailEtaToPi0Ratio7TeVRebined->SetLineStyle(5);
	cocktailEtaToPi0Ratio7TeVRebined->SetLineWidth(2.5);
	cocktailEtaToPi0Ratio_K0ScaledRebinned->SetLineStyle(6);
	cocktailEtaToPi0Ratio_K0ScaledRebinned->SetLineWidth(2.5);
	cocktailEtaToPi0Ratio_MtScaledRebinned->SetLineStyle(7);
	cocktailEtaToPi0Ratio_MtScaledRebinned->SetLineWidth(2.5);
	cocktailEtaToPi0Ratio7TeVRebined->Draw("same,hist,c");
	cocktailEtaToPi0Ratio_K0ScaledRebinned->Draw("same,hist,c");
	cocktailEtaToPi0Ratio_MtScaledRebinned->Draw("same,hist,c");
	
	TLegend* legendRatioALICEdata = new TLegend(0.12,0.82,0.32,0.92);
	legendRatioALICEdata->SetTextSize(0.04);			
	legendRatioALICEdata->SetFillColor(0);
	legendRatioALICEdata->SetFillStyle(0);
	legendRatioALICEdata->SetBorderSize(0);
// 	legendRatioALICEdata->SetMargin(0.1);
	legendRatioALICEdata->AddEntry(graphCombEtaToPi0RatioSysErrpp7TeV,Form("ALICE, %s",collisionSystemPP7TeV.Data()),"pe");
	legendRatioALICEdata->AddEntry((TObject*)0,"Phys.Lett. B717 (2012) 162-172","");
	legendRatioALICEdata->Draw();
	
	TLegend* legendRatioALICE3 = new TLegend(0.25,0.16,0.977,0.34);
	legendRatioALICE3->SetTextSize(0.04);			
	legendRatioALICE3->SetFillColor(0);
	legendRatioALICE3->SetFillStyle(0);
	legendRatioALICE3->SetBorderSize(0);
	legendRatioALICE3->SetMargin(0.06);
	legendRatioALICE3->SetHeader("Cocktail Calculations");
	legendRatioALICE3->AddEntry(cocktailEtaToPi0Ratio7TeVRebined,Form("#eta from %s as input",collisionSystemPP7TeV.Data()),"pl");
	legendRatioALICE3->AddEntry(cocktailEtaToPi0Ratio_MtScaledRebinned,"#eta from m_{T} scaled #pi^{0}","pl");
	legendRatioALICE3->AddEntry(cocktailEtaToPi0Ratio_K0ScaledRebinned,"0-20% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV, #eta from K^{0}_{s} scaled ","pl");
	
	legendRatioALICE3->Draw();
	
	canvasRatioEtaToPi0ALICE->Update();
// 	canvasRatioEtaToPi0ALICE->SaveAs(Form("%s/EtaToPi0Ratio_pp_7TeV_withCocktails.%s",outputDir.Data(), suffix.Data()));

	//*********************** Comparison of Measuremnts in ALICE ***************************************************
	TCanvas* canvasRatioOmegaToPi0ALICE = new TCanvas("canvasRatioOmegaToPi0ALICE","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings( canvasRatioOmegaToPi0ALICE, 0.09, 0.01, 0.015, 0.115);
	
	TH2D *histo2DRatioOmegaToPi0ALICE = new TH2D("histo2DRatioOmegaToPi0ALICE", "histo2DRatioOmegaToPi0ALICE", 200,0.,20.05,1000.,-0.4,1.7);
	SetStyleHistoTH2ForGraphs(histo2DRatioOmegaToPi0ALICE, "#it{p}_{T} (GeV/#it{c})","#omega/#pi^{0}", 0.046,0.058, 0.046,0.058, 0.85,0.74, 510, 510);
	histo2DRatioOmegaToPi0ALICE->GetYaxis()->SetRangeUser(0.,1.55);
	histo2DRatioOmegaToPi0ALICE->Draw();
	
	DrawGammaSetMarkerTGraphErr(graphOmegaToPi0Ratiopp7TeV, markerStyleSpectrum7TeV, markerSizePi0PP7TeV, color7TeV, color7TeV);
	DrawGammaSetMarkerTGraphErr(graphOmegaToPi0RatioSysErrpp7TeV, markerStyleSpectrum7TeV, markerSizePi0PP7TeV, color7TeV, color7TeV, widthLinesBoxes, kTRUE);
	
	
	
	histo2DRatioOmegaToPi0ALICE->Draw();

	TF1* fitOmegaToPi0pp7TeVHighPt = new TF1("fitOmegaToPi0pp7TeVHighPt","[0]");
	graphOmegaToPi0Ratiopp7TeV->Fit(fitOmegaToPi0pp7TeVHighPt,"NRMEX0+","",6,17.);
	cout << WriteParameterToFile(fitOmegaToPi0pp7TeVHighPt)<< endl;		
	fitOmegaToPi0pp7TeVHighPt->SetRange(6,20);

	graphOmegaToPi0RatioSysErrpp7TeV->Draw("same,pE2");	
	graphOmegaToPi0Ratiopp7TeV->Draw("same,pe");
	
	fitOmegaToPi0pp7TeVHighPt->Draw("same");
	
	fitOmegaToPi0pp7TeV->SetRange(0,20);
	DrawGammaSetMarkerTF1(fitOmegaToPi0pp7TeV, 1, widthCommonFit, kGray);
	fitOmegaToPi0pp7TeV->SetLineStyle(4);
// 	fitOmegaToPi0pp7TeV->Draw("same");
	
	fitOmegaToPi0pp7TeVFitted->SetRange(0,20);
	DrawGammaSetMarkerTF1(fitOmegaToPi0pp7TeVFitted, 1, widthCommonFit, kBlack);
// 	fitOmegaToPi0pp7TeVFitted->Draw("same");
	
	
	DrawGammaSetMarker(cocktailOmegaToPi0Ratio_MtScaledRebinned, 2, 0, kGreen+2, kGreen+2);
	cocktailOmegaToPi0Ratio_MtScaledRebinned->GetXaxis()->SetRangeUser(0.2,20);
	cocktailOmegaToPi0Ratio_MtScaledRebinned->SetLineStyle(7);
	cocktailOmegaToPi0Ratio_MtScaledRebinned->SetLineWidth(2.5);
	cocktailOmegaToPi0Ratio_MtScaledRebinned->Draw("same,hist,c");
	
	TLegend* legendOmegaRatioALICEdata = new TLegend(0.12,0.82,0.32,0.92);
	legendOmegaRatioALICEdata->SetTextSize(0.04);			
	legendOmegaRatioALICEdata->SetFillColor(0);
	legendOmegaRatioALICEdata->SetBorderSize(0);
// 	legendOmegaRatioALICEdata->SetMargin(0.1);
	legendOmegaRatioALICEdata->AddEntry(graphOmegaToPi0RatioSysErrpp7TeV,Form("ALICE preliminary, %s",collisionSystemPP7TeV.Data()),"pe");
// 	legendOmegaRatioALICEdata->AddEntry((TObject*)0,"Phys.Lett. B717 (2012) 162-172","");
	legendOmegaRatioALICEdata->Draw();

	TLegend* legendRatioOmegaCocktail = new TLegend(0.25,0.18,0.977,0.30);
	legendRatioOmegaCocktail->SetTextSize(0.04);			
	legendRatioOmegaCocktail->SetFillColor(0);
	legendRatioOmegaCocktail->SetBorderSize(0);
	legendRatioOmegaCocktail->SetMargin(0.06);
	legendRatioOmegaCocktail->SetHeader("Cocktail Calculations");
	legendRatioOmegaCocktail->AddEntry(cocktailOmegaToPi0Ratio_MtScaledRebinned,"#omega from m_{T} scaled #pi^{0}","pl");
	legendRatioOmegaCocktail->Draw();
	
	canvasRatioOmegaToPi0ALICE->Update();
// 	canvasRatioOmegaToPi0ALICE->SaveAs(Form("%s/OmegaToPi0Ratio_pp_7TeV_withCocktails.%s",outputDir.Data(), suffix.Data()));

	
	histo2DRatioOmegaToPi0ALICE->Draw();
	cocktailOmegaToPi0Ratio_MtScaledRebinned->Draw("same,hist,c");
	DrawGammaSetMarker(cocktailOmegaToPi0Ratio7TeVRebined, 2, 0, kBlue+2, kBlue+2);
	cocktailOmegaToPi0Ratio7TeVRebined->Draw("same,hist,c");
	DrawGammaSetMarker(cocktailOmegaToPi0RatioPbPb0040Rebined, 2, 0, kRed+2, kRed+2);
	cocktailOmegaToPi0RatioPbPb0040Rebined->Draw("same,hist,c");
	
	canvasRatioOmegaToPi0ALICE->Update();
// 	canvasRatioOmegaToPi0ALICE->SaveAs(Form("%s/OmegaToPi0Ratio_differntMtscaledVersions.%s",outputDir.Data(), suffix.Data()));

	
	//*********************** Comparison of Measuremnts in ALICE ***************************************************
	TCanvas* canvasRatioPhiPi0ALICE = new TCanvas("canvasRatioPhiPi0ALICE","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings( canvasRatioPhiPi0ALICE, 0.09, 0.01, 0.015, 0.115);
	
	TH2D *histo2DRatioPhiPi0ALICE = new TH2D("histo2DRatioPhiPi0ALICE", "histo2DRatioPhiPi0ALICE", 200,0.,6,1000.,-0.4,1.7);
	SetStyleHistoTH2ForGraphs(histo2DRatioPhiPi0ALICE, "#it{p}_{T} (GeV/#it{c})","#phi/#pi", 0.046,0.058, 0.046,0.058, 0.85,0.74, 510, 510);
	histo2DRatioPhiPi0ALICE->GetYaxis()->SetRangeUser(-0.0,0.55);
	histo2DRatioPhiPi0ALICE->Draw();

	int j;	
	TGraphErrors* graphPhiToPiPbPb0010=new TGraphErrors(8);
	TBox* boxPhiToPiPbPb0010[8];
	get_phi2pi(graphPhiToPiPbPb0010,boxPhiToPiPbPb0010);
	TGraphErrors* graphPhiToPipp7TeV=new TGraphErrors(26);
	TBox* boxPhiToPipp7TeV[26];
	get_phi2pipp7TeV(graphPhiToPipp7TeV,boxPhiToPipp7TeV);

	
	TF1* fitPhiToPi0pp7TeVHighPt = new TF1("fitPhiToPi0pp7TeVHighPt","[0]");
	graphPhiToPipp7TeV->Fit(fitPhiToPi0pp7TeVHighPt,"NRMEX0+","",3,6.);
	cout << WriteParameterToFile(fitPhiToPi0pp7TeVHighPt)<< endl;		
	fitPhiToPi0pp7TeVHighPt->SetRange(3,10);

	
	DrawGammaSetMarkerTGraphErr(graphPhiToPipp7TeV, markerStyleSpectrum7TeV, markerSizePi0PP7TeV, color7TeV, color7TeV);
	for(j=0;j<26;j++){boxPhiToPipp7TeV[j]->SetLineColor(color7TeV); boxPhiToPipp7TeV[j]->SetFillStyle(0);}
	
	
// 	DrawGammaSetMarkerTGraphAsym(graphCombPhiPi0RatioCombpp7TeV, markerStyleSpectrum7TeV, markerSizePi0PP7TeV, color7TeV-4, color7TeV-4);
	
	DrawGammaSetMarkerTGraphErr(graphPhiToPiPbPb0010, markerStylePbPb0010, markerSizePi0PP7TeV, colorCombPbPb0010, colorCombPbPb0010);
	for(j=0;j<8;j++){boxPhiToPiPbPb0010[j]->SetLineColor(colorCombPbPb0010); boxPhiToPiPbPb0010[j]->SetFillStyle(0);}
	
	for(j=0;j<26;j++) boxPhiToPipp7TeV[j]->Draw();
	for(j=0;j<8;j++) boxPhiToPiPbPb0010[j]->Draw();
	
// 	graphCombPhiPi0RatioCombpp7TeV->Draw("pzsame");
	graphPhiToPipp7TeV->Draw("pzsame");
	graphPhiToPiPbPb0010->Draw("pzsame");
	fitPhiToPi0pp7TeVHighPt->Draw("same");
	
// 	fitPhiDivPi0pp7TeV->SetRange(0,20);
// 	DrawGammaSetMarkerTF1(fitPhiDivPi0pp7TeV, 1, widthCommonFit, kBlack);
// 	fitPhiDivPi0pp7TeV->Draw("same");
	
	DrawGammaSetMarker(cocktailPhiPi0RatioPbPb0040Rebined, 2, 0, kGreen+2, kGreen+2);
	cocktailPhiPi0RatioPbPb0040Rebined->GetXaxis()->SetRangeUser(0.2,20);
	cocktailPhiPi0RatioPbPb0040Rebined->SetLineStyle(7);
	cocktailPhiPi0RatioPbPb0040Rebined->SetLineWidth(2.5);
	cocktailPhiPi0RatioPbPb0040Rebined->Draw("same,hist,c");
	
	
	TLegend* legendPhiRatioALICEdata = new TLegend(0.13,0.68,0.48,0.93);
	legendPhiRatioALICEdata->SetTextSize(0.04);			
	legendPhiRatioALICEdata->SetFillColor(0);
	legendPhiRatioALICEdata->SetBorderSize(0);
	legendPhiRatioALICEdata->SetMargin(0.15);
	legendPhiRatioALICEdata->AddEntry(graphPhiToPipp7TeV,"#phi/(0.5 (#pi^{-} + #pi^{+})) ALICE preliminary","pe");
	legendPhiRatioALICEdata->AddEntry((TObject*)0,"pp #sqrt{#it{s}} = 7 TeV","");
	legendPhiRatioALICEdata->AddEntry(graphPhiToPiPbPb0010,"#phi/(0.5 (#pi^{-} + #pi^{+})) ALICE, ","pe");
	legendPhiRatioALICEdata->AddEntry((TObject*)0,"0-10% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV","");
	legendPhiRatioALICEdata->AddEntry((TObject*)0,"arXiv:1404.0495","");
	legendPhiRatioALICEdata->Draw();


	TLegend* legendRatioPhiCocktail = new TLegend(0.53,0.15,0.95,0.27);
	legendRatioPhiCocktail->SetTextSize(0.04);			
	legendRatioPhiCocktail->SetFillStyle(0);
	legendRatioPhiCocktail->SetFillColor(0);
	legendRatioPhiCocktail->SetBorderSize(0);
	legendRatioPhiCocktail->SetMargin(0.1);
	legendRatioPhiCocktail->SetHeader("Cocktail Calculations");
	legendRatioPhiCocktail->AddEntry(cocktailPhiPi0RatioPbPb0040Rebined,"#phi/#pi^{0}, #phi from m_{T} scaled #pi^{0}","pl");
	legendRatioPhiCocktail->Draw();

	
	canvasRatioPhiPi0ALICE->Update();
// 	canvasRatioPhiPi0ALICE->SaveAs(Form("%s/PhiToPi0Ratio_PbPb0010Andpp_withCocktails.%s",outputDir.Data(), suffix.Data()));

	
	
	
	//*********************** Comparison of Measuremnts in ALICE ***************************************************
	TCanvas* canvasRatioEtaToPhiALICE = new TCanvas("canvasRatioEtaToPhiALICE","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings( canvasRatioEtaToPhiALICE, 0.09, 0.01, 0.015, 0.115);
	
	TH2D *histo2DRatioEtaToPhi = new TH2D("histo2DRatioEtaToPhi", "histo2DRatioPhiPi0ALICE", 200,0.,20.,1000.,-0.,15);
	SetStyleHistoTH2ForGraphs(histo2DRatioEtaToPhi, "#it{p}_{T} (GeV/#it{c})","#eta/#phi", 0.046,0.058, 0.046,0.058, 0.85,0.74, 510, 510);
// 	histo2DRatioEtaToPhi->GetYaxis()->SetRangeUser(-0.0,0.55);
	histo2DRatioEtaToPhi->Draw();

	
	DrawGammaSetMarker(ratioEtaToPhi7TeVConstructed, markerStyleSpectrum7TeV, markerSizePi0PP7TeV, color7TeV, color7TeV);
	ratioEtaToPhi7TeVConstructed->Draw("clsame");
	
	canvasRatioEtaToPhiALICE->Update();
// 	canvasRatioEtaToPhiALICE->SaveAs(Form("%s/EtaToPhiRatio_pp7TeV.%s",outputDir.Data(), suffix.Data()));

	//*********************** Comparison of Measuremnts in ALICE ***************************************************
	TCanvas* canvasRatioOmegaToPhiALICE = new TCanvas("canvasRatioOmegaToPhiALICE","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings( canvasRatioOmegaToPhiALICE, 0.09, 0.01, 0.015, 0.115);
	
	TH2D *histo2DRatioOmegaToPhi = new TH2D("histo2DRatioOmegaToPhi", "histo2DRatioOmegaToPhi", 200,0.,20.,1000.,-0.,15);
	SetStyleHistoTH2ForGraphs(histo2DRatioOmegaToPhi, "#it{p}_{T} (GeV/#it{c})","#omega/#phi", 0.046,0.058, 0.046,0.058, 0.85,0.74, 510, 510);
// 	histo2DRatioEtaToPhi->GetYaxis()->SetRangeUser(-0.0,0.55);
	histo2DRatioOmegaToPhi->Draw();

	
	DrawGammaSetMarker(ratioOmegaToPhi7TeVConstructed, markerStyleSpectrum7TeV, markerSizePi0PP7TeV, color7TeV, color7TeV);
	ratioOmegaToPhi7TeVConstructed->Draw("clsame");
	
	canvasRatioOmegaToPhiALICE->Update();
// 	canvasRatioOmegaToPhiALICE->SaveAs(Form("%s/OmegaToPhiRatio_pp7TeV.%s",outputDir.Data(), suffix.Data()));
	
	
	
   TFile fEtatoPi0input("EtaToPi0InputsForCombination.root","UPDATE");
	histoPCMEtaToPi0RatioPbPb0010->Write("histoPCMEtaToPi0RatioPbPb0010");
	histoPCMEtaToPi0RatioPbPb2050->Write("histoPCMEtaToPi0RatioPbPb2050");
	histoPCMEtaToPi0RatiopPb->Write("histoPCMEtaToPi0RatiopPb");
	graphPCMEtaToPi0RatioSysErrpPb->Write("graphPCMEtaToPi0RatioSysErrpPb");
	
	graphCombEtaToPi0RatioSysErrpp7TeV->Write("graphCombEtaToPi0RatioSysErrpp7TeV");
	graphCombEtaToPi0Ratiopp7TeVNoXErrors->Write("graphCombEtaToPi0Ratiopp7TeVNoXErrors");

	cocktailEtaToPi0Ratio7TeVRebined->Write("cocktailEtaToPi0Ratio7TeVRebined");
	cocktailEtaToPi0Ratio_MtScaledRebinned->Write("cocktailEtaToPi0Ratio_MtScaledRebinned");
	cocktailEtaToPi0Ratio_K0ScaledRebinned->Write("cocktailEtaToPi0Ratio_K0ScaledRebinned");

   fEtatoPi0input.Close();

	
	
	
}
