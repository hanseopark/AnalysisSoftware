// Fit Functions
// Just call TF1 *function = FitObject(TObject *Obj_Dummy, Int_t type, TString FunctinName, Double_t xmin = 0., Double_t xmax = 7.,Double_t Parameter = NULL, TString FitOptions = "NQRME+");

#ifndef GAMMACONV_Fitting
#define GAMMACONV_Fitting

    TF1* FitObject( TString type, 
                    TString FunctionName, 
                    TString mesonType       = "Pi0", 
                    TObject *Obj_Dummy      = NULL,  
                    Double_t xmin           = 0., 
                    Double_t xmax           = 10.,
                    Double_t Parameter[]    = NULL, 
                    TString FitOptions      = "IQNRME+", 
                    TString fixingPar       = "",
                    Bool_t limitPar         = kTRUE
                ){

        Double_t mass;

        if (mesonType.CompareTo("Pi0")==0){
            mass = TDatabasePDG::Instance()->GetParticle(111)->Mass();
            cout <<"Meson Mass: "<< mass << endl;
        } else if (mesonType.CompareTo("ChargedPi")==0){
            mass = TDatabasePDG::Instance()->GetParticle(211)->Mass();
            cout <<"Meson Mass: "<< mass << endl;
        } else if (mesonType.CompareTo("Eta")==0){
            mass = TDatabasePDG::Instance()->GetParticle(221)->Mass();
            cout <<"Meson Mass: "<< mass << endl;
        } else if (mesonType.CompareTo("Omega")==0){
            mass = TDatabasePDG::Instance()->GetParticle(223)->Mass();
            cout <<"Meson Mass: "<< mass << endl;
        } else if (mesonType.CompareTo("Phi")==0){
            mass = TDatabasePDG::Instance()->GetParticle(333)->Mass();
            cout <<"Meson Mass: "<< mass << endl;
        } else if (mesonType.CompareTo("K")==0){
            mass = TDatabasePDG::Instance()->GetParticle(321)->Mass();
            cout <<"Meson Mass: "<< mass << endl;
        } else if (mesonType.CompareTo("K0s")==0){
            mass = TDatabasePDG::Instance()->GetParticle(310)->Mass();
            cout <<"Meson Mass: "<< mass << endl;
        } else if (mesonType.CompareTo("P")==0){
            mass = TDatabasePDG::Instance()->GetParticle(2212)->Mass();
            cout <<"Baryon Mass: "<< mass << endl;
        } else if (mesonType.CompareTo("Lambda")==0){
            mass = TDatabasePDG::Instance()->GetParticle(3122)->Mass();
            cout <<"Baryon Mass: "<< mass << endl;
        } else {
            mass = 0;
            cout <<"Meson Mass: "<< mass << endl;
        }

        if(Obj_Dummy == NULL) {
            if(type.BeginsWith("h") || type.BeginsWith("H")){
                TF1 *Hagedorn_Dummy = new TF1("Hagedorn_Dummy","[0]*pow([2]/([2]+x),[1])");
                Hagedorn_Dummy->SetParNames("C_{H}","n","p_{0} (GeV/c)");
                Hagedorn_Dummy->SetParameters(19.,6.8,0.84);
                Hagedorn_Dummy->SetName(FunctionName);
                return Hagedorn_Dummy;
            } else if(type.CompareTo("powPure")==0 || type.CompareTo("PowPure")==0){
                cout <<Form("fitting %s with Pure Powerlaw",FunctionName.Data()) << endl;
                TF1 *PowerLaw_Dummy = new TF1("PowerLawPure_Dummy","[0] * 1/pow(x,[1])");
                PowerLaw_Dummy->SetParNames("A_{pow}","n");
                PowerLaw_Dummy->SetParameters(2.,5.);
                PowerLaw_Dummy->SetName(FunctionName);
                return PowerLaw_Dummy;
            } else if(type.BeginsWith("p") || type.BeginsWith("P")){
                TF1 *PowerLaw_Dummy = new TF1("PowerLaw_Dummy","[0] * 2 / TMath::Pi() * ([1]-1.)*([1]-2.)/pow([1]-3.,2) /x * pow(1+2*x/[2]/([1]-3),-[1])");
                PowerLaw_Dummy->SetParNames("A_{pow}","n","p_{0}");
                PowerLaw_Dummy->SetParameters(2.,5.,0.37);
                PowerLaw_Dummy->SetName(FunctionName);
                return PowerLaw_Dummy;
            } else if(type.BeginsWith("l") || type.BeginsWith("L")){
                TF1 *Levy_Dummy = new TF1("Levy_Dummy",Form("[0] / ( 2 * TMath::Pi())*([1]-1.)*([1]-2.) / ([1]*[2]*([1]*[2]+%.10f*([1]-2.)))  * pow(1.+(sqrt(x*x+%.10f*%.10f)-%.10f)/([1]*[2]), -[1])",mass,mass,mass,mass));
                Levy_Dummy->SetParNames("dN/dy","n","T_{Levy} (GeV/c)");
                Levy_Dummy->SetParameters(2.,5.,0.18); // standard parameter optimize if necessary
                Levy_Dummy->SetName(FunctionName);
                return Levy_Dummy;
            } else if(type.BeginsWith("b") || type.BeginsWith("B")){
                TF1 *Boltzmann_Dummy = new TF1("Boltzmann_Dummy",Form("[0] *sqrt(x*x+%.10f*%.10f)* exp(- sqrt(x*x+%.10f*%.10f)/[1])",mass,mass,mass,mass));
                Boltzmann_Dummy->SetParNames("C_{B}","T_{Boltzm.} (GeV/c)");
                Boltzmann_Dummy->SetParameters(2.5,0.3); // standard parameter optimize if necessary
                Boltzmann_Dummy->SetName(FunctionName);
                return Boltzmann_Dummy;
            } else if(type.BeginsWith("e") || type.BeginsWith("E")){
                TF1 *Exponential_Dummy = new TF1("Exponential_Dummy",Form("[0]/( 2 * TMath::Pi())/([1]*(%.10f+[1]))*exp(-(sqrt(x*x+%.10f*%.10f)-%.10f)/[1])",mass,mass,mass,mass));
                Exponential_Dummy->SetParameters(1.1,0.37); // standard parameter optimize if necessary
                Exponential_Dummy->SetParNames("C_{E}","T_{exp} (GeV/c)");
                Exponential_Dummy->SetName(FunctionName);
                return Exponential_Dummy;
            } else if(type.BeginsWith("m") || type.BeginsWith("M")){
                TF1 *ModPowerLaw_Dummy = new TF1("ModPowerLaw_Dummy","[0]*pow((1 + (x)/[1]),-[2])");
                ModPowerLaw_Dummy->SetParameters(2.,0.37,5.); // standard parameter optimize if necessary
                ModPowerLaw_Dummy->SetParNames("A","p_{0}","n");
                ModPowerLaw_Dummy->SetName(FunctionName);
                return ModPowerLaw_Dummy;
            } else if(type.BeginsWith("6pol") || type.BeginsWith("6POL")){
                TF1 *Pol6_Dummy = new TF1("Pol6_Dummy","[0]+[1]*pow(x,2)+[2]*pow(x,4)+[3]*pow(x,6)");
                Pol6_Dummy->SetParameters(1.,1.,1.,1.); // standard parameter optimize if necessary
                Pol6_Dummy->SetParNames("a","b","c","d");
                Pol6_Dummy->SetName(FunctionName);
                return Pol6_Dummy;
            } else if(type.BeginsWith("doubqcd") || type.BeginsWith("doubqcd")){
                cout <<Form("fitting %s with doubqcd",FunctionName.Data()) << endl;
                TF1 *QCD_Dummy = new TF1("QCD_Dummy","(x<=[5])*[0]*TMath::Power(x,-1*([1]+[2]/(TMath::Power(x,[3])+[4])))+(x>[5])*[6]*TMath::Power(x,-1*([7]+[8]/(TMath::Power(x,[9])+[10])))");
                QCD_Dummy->SetParameters(24,6.7,-6.5,1.,10,24,2,6.7,-6.5,1.,10); // standard parameter optimize if necessary
                QCD_Dummy->SetParNames("a1","b1","c1","d1","e1","pT","a2","b2","c2","d2","e2");
                QCD_Dummy->SetName(FunctionName);
                return QCD_Dummy;
            } else if(type.BeginsWith("qcdtsal") || type.BeginsWith("qcdtsal")){
                cout <<Form("fitting %s with qcdtsal",FunctionName.Data()) << endl;
                TF1 *QCD_Dummy = new TF1("QCD_Dummy","(x<=[5])*[0]*TMath::Power(x,-1*([1]+[2]/(TMath::Power(x,[3])+[4])))+(x>[5])*[6] / ( 2 * TMath::Pi())*([7]-1.)*([7]-2.) / ([7]*[8]*([7]*[8]+%.10f*([7]-2.)))  * pow(1.+(sqrt(x*x+%.10f*%.10f)-%.10f)/([7]*[8]), -[7])");
                QCD_Dummy->SetParameters(24,6.7,-6.5,1.,10,2,2.,5.,0.18); // standard parameter optimize if necessary
                QCD_Dummy->SetParNames("a1","b1","c1","d1","e1","pT","a2","b2","c2");
                QCD_Dummy->SetName(FunctionName);
                return QCD_Dummy;
            } else if(type.BeginsWith("rad") || type.BeginsWith("RAD")){
                TF1 *Rad_Dummy = new TF1("Rad_Dummy",Form("(x<=[3])*x*[0]*([1]-1.)*([1]-2.)/([1]*[2]*([1]*[2]+%.10f*([1]-2.)))*pow(1.+(sqrt(x*x+%.10f*%.10f)-%.10f)/([1]*[2]), -[1]) + (x>[3]) * x *[0] * ([1]-1.)*([1]-2.)/([1]*[2]*([1]*[2]+%.10f*([1]-2.))) * pow([1]*[2]/([3]+[1]*[2]-%.10f),[1]) * pow([3],[4]) * pow( 1./x, [4] )",mass,mass,mass,mass,mass,mass));
                Rad_Dummy->SetLineWidth(1);
                Double_t par0min3 = 150.;
                Double_t par0max3 = 1500.;
                Double_t par1min3 = 4.;
                Double_t par1max3 = 40.;
                Double_t par2min3 = 0.007;
                Double_t par2max3 = .7;
                Double_t par3min3 = 4.36;
                Double_t par3max3 = 4.40;
                Double_t par4min3 = 2.;
                Double_t par4max3 = 18.;
                Rad_Dummy->SetParLimits(0, par0min3, par0max3);
                Rad_Dummy->SetParLimits(1, par1min3, par1max3);
                Rad_Dummy->SetParLimits(2, par2min3, par2max3);
                Rad_Dummy->SetParLimits(3, par3min3, par3max3);
                Rad_Dummy->SetParLimits(4, par4min3, par4max3);
                Rad_Dummy->SetParNames("a","b","c","d");
                Rad_Dummy->SetName(FunctionName);
                return Rad_Dummy;
            } else if(type.BeginsWith("tmpt") || type.BeginsWith("TMPT")){
                // Tsallis Dummy multiplied with pt
                TF1 *Levy_Dummy = new TF1("Tsallis_Dummy",Form("[0] / ( 2 * TMath::Pi())*([1]-1.)*([1]-2.) / ([1]*[2]*([1]*[2]+%.10f*([1]-2.)))  * x* pow(1.+(sqrt(x*x+%.10f*%.10f)-%.10f)/([1]*[2]), -[1])",mass,mass,mass,mass));
                Levy_Dummy->SetParNames("dN/dy","n","T_{Levy} (GeV/c)");
                Levy_Dummy->SetParameters(2.e11,7., 0.137) ; // standard parameter optimize if necessary
                Levy_Dummy->SetName(FunctionName);
                return Levy_Dummy;
            } else if(type.BeginsWith("hmpt") || type.BeginsWith("HMPT")){
                // Hagedorn Dummy multiplied with pt
                TF1 *Hagedorn_Dummy = new TF1("Hagedorn_Dummy","[0]*x*pow([2]/([2]+x),[1])");
                Hagedorn_Dummy->SetParNames("C_{H}","n","p_{0} (GeV/c)");
                Hagedorn_Dummy->SetParameters(1.,7.,0.37);
                Hagedorn_Dummy->SetName(FunctionName);
                return Hagedorn_Dummy;
            } else if(type.BeginsWith("qmpt") || type.BeginsWith("QMPT")){
                TF1 *QCD_Dummy = new TF1("QCD_Dummy","([0]*x*TMath::Power(x,-1*([1]+[2]/(TMath::Power(x,[3])+[4]))))");
                QCD_Dummy->SetParNames("a","b","c","d","e");
                QCD_Dummy->SetParameters(24,6.7,-6.5,1.,10); // standard parameter optimize if necessary
                QCD_Dummy->SetName(FunctionName);
                return QCD_Dummy;
            } else if(type.BeginsWith("oHag") || type.BeginsWith("OHag")){
                cout << "entered"<< endl;
                TF1 *ModPowerLaw_Dummy2 = new TF1("ModHagedorn_Dummy","[0]*pow(exp(-[1]*x-abs([2])*x*x)+x/[3],-[4])");
                ModPowerLaw_Dummy2->SetParNames("a","b","c","d","e");
                ModPowerLaw_Dummy2->SetParameters(30.,0.37,0.07,0.68,6.1);
                ModPowerLaw_Dummy2->SetName(FunctionName);
                return ModPowerLaw_Dummy2;
            } else if(type.BeginsWith("mohag") || type.BeginsWith("MOHAG")){
                cout <<Form("fitting %s with ModHagedorn",FunctionName.Data()) << endl;
                TF1 *ModPowerLaw_Dummy = new TF1("ModHagedorn_Dummy","[0]*x*pow(exp(-[1]*x-abs([2])*x*x)+x/[3],-[4])");
                //TF1 *ModPowerLaw_Dummy = new TF1("ModHagedorn_Dummy","[0]*pow(exp(-[1]*x)+x/[2],-[3])");
                ModPowerLaw_Dummy->SetParameters(450.,0.37,0.2,0.7,8.); // standard parameter optimize if necessary
                ModPowerLaw_Dummy->SetParNames("a","b","c","d","e");
                ModPowerLaw_Dummy->SetName(FunctionName);
                return ModPowerLaw_Dummy;
            } else if(type.BeginsWith("tcmlow") || type.BeginsWith("TCMLOW")){ // Two component model fit [A. Bylinkin and A. Rostovtsev, Phys. Atom. Nucl 75 (2012) 999-1005]
                cout <<Form("fitting %s with low Pt two component model by Bylinkin",FunctionName.Data()) << endl;
                TF1 *TwoCompModel_Dummy = new TF1("twoCompModelLow_Dummy",Form("[0]*exp(-(TMath::Sqrt(x*x+%.10f*%.10f)-%.10f)/[1])",mass,mass,mass));
                if (mesonType.CompareTo("Pi0")==0){
                    TwoCompModel_Dummy->SetParameters(450.,0.3); // standard parameter optimize if necessary
                } else if (mesonType.CompareTo("Eta")==0){    
                    TwoCompModel_Dummy->SetParameters(450.,0.3); // standard parameter optimize if necessary
                }
                TwoCompModel_Dummy->SetParNames("Ae","Te");
                TwoCompModel_Dummy->SetName(FunctionName);
                return TwoCompModel_Dummy;
            } else if(type.BeginsWith("tcmhigh") || type.BeginsWith("TCMHIGH")){ // Two component model fit [A. Bylinkin and A. Rostovtsev, Phys. Atom. Nucl 75 (2012) 999-1005]
                cout <<Form("fitting %s with high pt two component model by Bylinkin",FunctionName.Data()) << endl;
                TF1 *TwoCompModel_Dummy = new TF1("twoCompModelHigh_Dummy","[0]/(TMath::Power(1+x*x/([1]*[1]*[2]),[2]) )");
                if (mesonType.CompareTo("Pi0")==0){
                    TwoCompModel_Dummy->SetParameters(1,0.3,8.); // standard parameter optimize if necessary
                } else if (mesonType.CompareTo("Eta")==0){    
                    TwoCompModel_Dummy->SetParameters(1,0.3,8.); // standard parameter optimize if necessary
                }
                TwoCompModel_Dummy->SetParNames("A","T","n");
                TwoCompModel_Dummy->SetName(FunctionName);
                return TwoCompModel_Dummy;
            } else if(type.BeginsWith("tcm") || type.BeginsWith("TCM")){ // Two component model fit [A. Bylinkin and A. Rostovtsev, Phys. Atom. Nucl 75 (2012) 999-1005]
                cout <<Form("fitting %s with two component model by Bylinkin",FunctionName.Data()) << endl;
                TF1 *TwoCompModel_Dummy = new TF1("twoCompModel_Dummy",Form("[0]*exp(-(TMath::Sqrt(x*x+%.10f*%.10f)-%.10f)/[1]) + [2]/(TMath::Power(1+x*x/([3]*[3]*[4]),[4]) )",mass,mass,mass));
                if (mesonType.CompareTo("Pi0")==0){
                    TwoCompModel_Dummy->SetParameters(450.,0.3,1,0.3,8.); // standard parameter optimize if necessary
                } else if (mesonType.CompareTo("Eta")==0){    
                    TwoCompModel_Dummy->SetParameters(450.,0.3,1,0.3,8.); // standard parameter optimize if necessary
                }
                TwoCompModel_Dummy->SetParNames("Ae","Te","A","T","n");
                TwoCompModel_Dummy->SetName(FunctionName);
                return TwoCompModel_Dummy;
            } else if(type.BeginsWith("tcmpt") || type.BeginsWith("TCMPT")){ // Two component model fit [A. Bylinkin and A. Rostovtsev, Phys. Atom. Nucl 75 (2012) 999-1005]
                cout <<Form("fitting %s with two component model by Bylinkin",FunctionName.Data()) << endl;
                TF1 *TwoCompModel_Dummy = new TF1("twoCompModel_Dummy",Form("x*[0]*exp(-(TMath::Sqrt(x*x+%.10f*%.10f)-%.10f)/[1]) + x*[2]/(TMath::Power(1+x*x/([3]*[3]*[4]),[4]) )",mass,mass,mass));
                if (mesonType.CompareTo("Pi0")==0){
                    TwoCompModel_Dummy->SetParameters(450.,0.3,1,0.3,8.); // standard parameter optimize if necessary
                } else if (mesonType.CompareTo("Eta")==0){    
                    TwoCompModel_Dummy->SetParameters(450.,0.3,1,0.3,8.); // standard parameter optimize if necessary
                }
                TwoCompModel_Dummy->SetParNames("Ae","Te","A","T","n");
                TwoCompModel_Dummy->SetName(FunctionName);
                return TwoCompModel_Dummy;
            }
            
        }

        TF1 *FitFunction = new TF1();
        TString ClassName = Obj_Dummy->ClassName();
        
        if(ClassName.BeginsWith("TH1")){
            TH1D *Obj = (TH1D*)Obj_Dummy;
            if(type.BeginsWith("xqcd") || type.BeginsWith("XQCD")){
                cout <<Form("fitting %s with XQCD",FunctionName.Data()) << endl;
                TF1 *XQCD_Dummy = new TF1("QCD_Dummy","[0]*TMath::Power(x,-1*([1]+[2]/(TMath::Power(x,[3])+[4]*TMath::Power(x,0.5))))");
                FitFunction = (TF1*)XQCD_Dummy->Clone(FunctionName);
                FitFunction->SetRange(xmin, xmax);
                if(Parameter == NULL)FitFunction->SetParameters(24,6.7,-6.5,1.,10,0.5); // standard parameter optimize if necessary
                else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4]);
                FitFunction->SetParNames("a","b","c","d","e");
                Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
            } else if(type.BeginsWith("doubqcd") || type.BeginsWith("QCD")){
                cout <<Form("fitting %s with QCD",FunctionName.Data()) << endl;
                TF1 *QCD_Dummy = new TF1("QCD_Dummy","(x<=[5])*[0]*TMath::Power(x,-1*([1]+[2]/(TMath::Power(x,[3])+[4])))+(x>[5])*[6]*TMath::Power(x,-1*([7]+[8]/(TMath::Power(x,[9])+[10])))");
                FitFunction = (TF1*)QCD_Dummy->Clone(FunctionName);
                FitFunction->SetRange(xmin, xmax);
                if(Parameter == NULL)FitFunction->SetParameters(24,6.7,-6.5,1.,10,24,2,6.7,-6.5,1.,10); // standard parameter optimize if necessary
                else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4],Parameter[5],Parameter[6],Parameter[7],Parameter[8],Parameter[9],Parameter[10]);
                FitFunction->SetParNames("a1","b1","c1","d1","e1","pT","a2","b2","c2","d2","e2");
                Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
            } else if(type.BeginsWith("qcd") || type.BeginsWith("QCD")){
                cout <<Form("fitting %s with QCD",FunctionName.Data()) << endl;
                TF1 *QCD_Dummy = new TF1("QCD_Dummy","[0]*TMath::Power(x,-1*([1]+[2]/(TMath::Power(x,[3])+[4])))");
                FitFunction = (TF1*)QCD_Dummy->Clone(FunctionName);
                FitFunction->SetRange(xmin, xmax);
                if(Parameter == NULL)FitFunction->SetParameters(24,6.7,-6.5,1.,10); // standard parameter optimize if necessary
                else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4]);
                FitFunction->SetParNames("a","b","c","d","e");
                Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
            } else if(type.BeginsWith("h") || type.BeginsWith("H")){
                cout <<Form("fitting %s with Hagedorn",FunctionName.Data()) << endl;
                TF1 *Hagedorn_Dummy = new TF1("Hagedorn_Dummy","[0]*pow([2]/([2]+x),[1])");
                FitFunction = (TF1*)Hagedorn_Dummy->Clone(FunctionName);
                FitFunction->SetRange(xmin, xmax);
                if(Parameter == NULL) FitFunction->SetParameters(19.,6.8,0.84);
                else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2]);
                FitFunction->SetParNames("C_{H}","n","p_{0} (GeV/c)");
                Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
            } else if(type.BeginsWith("p") || type.BeginsWith("P")){
                cout <<Form("fitting %s with Powerlaw",FunctionName.Data()) << endl;
                TF1 *PowerLaw_Dummy = new TF1("PowerLaw_Dummy","[0] * 2 / TMath::Pi() * ([1]-1.)*([1]-2.)/pow([1]-3.,2) /x * pow(1+2*x/[2]/([1]-3),-[1])");
                FitFunction = (TF1*)PowerLaw_Dummy->Clone(FunctionName);
                FitFunction->SetRange(xmin, xmax);
                if(Parameter == NULL)FitFunction->SetParameters(2.,5.,0.37); // standard parameter optimize if necessary
                else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2]);
                FitFunction->SetParNames("A_{pow}","n","p_{0}");
                Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
            } else if(type.BeginsWith("l") || type.BeginsWith("L")){
                cout <<Form("fitting %s with Levy",FunctionName.Data()) << endl;
                TF1 *Levy_Dummy = new TF1("Levy_Dummy",Form("[0] / ( 2 * TMath::Pi())*([1]-1.)*([1]-2.) / ([1]*[2]*([1]*[2]+%.10f*([1]-2.)))  * pow(1.+(sqrt(x*x+%.10f*%.10f)-%.10f)/([1]*[2]), -[1])",mass,mass,mass,mass));
                FitFunction = (TF1*)Levy_Dummy->Clone(FunctionName);
                FitFunction->SetRange(xmin, xmax);
                if(Parameter == NULL)FitFunction->SetParameters(2.,5.,0.18); // standard parameter optimize if necessary
                else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2]);
                FitFunction->SetParNames("dN/dy","n","T_{Levy} (GeV/c)");
                Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
            } else if(type.BeginsWith("b") || type.BeginsWith("B")){
                cout <<Form("fitting %s with Boltzmann",FunctionName.Data()) << endl;
                TF1 *Boltzmann_Dummy =  new TF1("Boltzmann_Dummy",Form("[0] *sqrt(x*x+%.10f*%.10f)* exp(- sqrt(x*x+%.10f*%.10f)/[1])",mass,mass,mass,mass));
                FitFunction = (TF1*)Boltzmann_Dummy->Clone(FunctionName);
                FitFunction->SetRange(xmin, xmax);
                if(Parameter == NULL)FitFunction->SetParameters(2.5,0.3); // standard parameter optimize if necessary
                else FitFunction->SetParameters(Parameter[0],Parameter[1]);
                FitFunction->SetParNames("C_{B}","T_{Boltzm.} (GeV/c)");
                Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
            } else if(type.BeginsWith("e") || type.BeginsWith("E")){
                cout <<Form("fitting %s with Exponential",FunctionName.Data()) << endl;
                TF1 *Exponential_Dummy = new TF1("Exponential_Dummy",Form("[0]/( 2 * TMath::Pi())/([1]*(%.10f+[1]))*exp(-(sqrt(x*x+%.10f*%.10f)-%.10f)/[1])",mass,mass,mass,mass));
                FitFunction = (TF1*)Exponential_Dummy->Clone(FunctionName);
                FitFunction->SetRange(xmin, xmax);
                if(Parameter == NULL)FitFunction->SetParameters(1.1,0.37); // standard parameter optimize if necessary
                else FitFunction->SetParameters(Parameter[0],Parameter[1]);
                FitFunction->SetParNames("C_{E}","T_{exp} (GeV/c)");
                Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
            } else if(type.BeginsWith("m") || type.BeginsWith("M")){
                cout <<Form("fitting %s with ModPowerlaw",FunctionName.Data()) << endl;
                TF1 *ModPowerLaw_Dummy = new TF1("ModPowerLaw_Dummy","[0]*pow((1 + (x)/[1]),-[2])");
                FitFunction = (TF1*)ModPowerLaw_Dummy->Clone(FunctionName);
                FitFunction->SetRange(xmin, xmax);
                if(Parameter == NULL)FitFunction->SetParameters(2.,0.37,5.); // standard parameter optimize if necessary
                else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2]);
                FitFunction->SetParNames("A","p_{0}","n");
                Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
            } else if(type.BeginsWith("oHag") || type.BeginsWith("OHag")){
                cout <<Form("fitting %s with ModHagedorn",FunctionName.Data()) << endl;
                TF1 *ModPowerLaw_Dummy = new TF1("ModHagedorn_Dummy","[0]*pow(exp(-[1]*x-abs([2])*x*x)+x/[3],-[4])");
                FitFunction = (TF1*)ModPowerLaw_Dummy->Clone(FunctionName);
                FitFunction->SetRange(xmin, xmax);
                if(Parameter == NULL)FitFunction->SetParameters(450.,0.37,0.2,0.7,8.); // standard parameter optimize if necessary
                else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4]);
                FitFunction->SetParNames("A (mbGeV^{-2}c^{3})","a [(GeV/c)^{-1}]","b [(GeV/c)^{-1}]","p_{0} (GeV/c)","n");
                Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
            } else if(type.BeginsWith("6pol") || type.BeginsWith("6POL")){
                cout <<Form("fitting %s with Polynom order 6",FunctionName.Data()) << endl;
                TF1 *Pol6_Dummy = new TF1("Pol6_Dummy","[0]+[1]*pow(x,2)+[2]*pow(x,4)+[3]*pow(x,6)");
                FitFunction = (TF1*)Pol6_Dummy->Clone(FunctionName);
                FitFunction->SetRange(xmin, xmax);
                if(Parameter == NULL)FitFunction->SetParameters(1.,1.,1.,1.); // standard parameter optimize if necessary
                else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3]);
                FitFunction->SetParNames("a","b","c","d");
                Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
            } else if(type.BeginsWith("rad") || type.BeginsWith("RAD")){
                cout <<Form("fitting %s with Radoslav Func",FunctionName.Data()) << endl;
                TF1 *Rad_Dummy = new TF1("Rad_Dummy",Form("(x<=[3])*x*[0]*([1]-1.)*([1]-2.)/([1]*[2]*([1]*[2]+%.10f*([1]-2.)))*pow(1.+(sqrt(x*x+%.10f*%.10f)-%.10f)/([1]*[2]), -[1]) + (x>[3]) * x *[0] * ([1]-1.)*([1]-2.)/([1]*[2]*([1]*[2]+%.10f*([1]-2.))) * pow([1]*[2]/([3]+[1]*[2]-%.10f),[1]) * pow([3],[4]) * pow( 1./x, [4] )",mass,mass,mass,mass,mass,mass));
                FitFunction = (TF1*)Rad_Dummy->Clone(FunctionName);
                FitFunction->SetRange(xmin, xmax);
                Double_t par0min3 = 10.;
                Double_t par0max3 = 1500.;
                Double_t par1min3 = 4.;
                Double_t par1max3 = 40.;
                Double_t par2min3 = 0.07;
                Double_t par2max3 = .7;
                Double_t par3min3 = 3.;
                Double_t par3max3 = 5.;
                Double_t par4min3 = 2.;
                Double_t par4max3 = 18.;

                if(Parameter != NULL){
                    par0min3 = Parameter[0];
                    par0max3 = Parameter[1];
                    par1min3 = Parameter[2];
                    par1max3 = Parameter[3];
                    par2min3 = Parameter[4];
                    par2max3 = Parameter[5];
                    par3min3 = Parameter[6];
                    par3max3 = Parameter[7];
                    par4min3 = Parameter[8];
                    par4max3 = Parameter[9];
                }
                FitFunction->SetParLimits(0, par0min3, par0max3);
                FitFunction->SetParLimits(1, par1min3, par1max3);
                FitFunction->SetParLimits(2, par2min3, par2max3);
                FitFunction->SetParLimits(3, par3min3, par3max3);
                FitFunction->SetParLimits(4, par4min3, par4max3);
                FitFunction->SetParNames("a","b","c","d");
                Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
                return FitFunction;
            } else if(type.BeginsWith("tcm") || type.BeginsWith("TCM")){ // Two component model fit [A. Bylinkin and A. Rostovtsev, Phys. Atom. Nucl 75 (2012) 999-1005]
                cout <<Form("fitting %s with two component model by Bylinkin",FunctionName.Data()) << endl;
                TF1 *TwoCompModel_Dummy = new TF1("twoCompModel_Dummy",Form("[0]*exp(-(TMath::Sqrt(x*x+%.10f*%.10f)-%.10f)/[1]) + [2]/(TMath::Power(1+x*x/([3]*[3]*[4]),[4]) )",mass,mass,mass));
                FitFunction = (TF1*)TwoCompModel_Dummy->Clone(FunctionName);
                FitFunction->SetRange(xmin, xmax);
                if (mesonType.CompareTo("Pi0")==0){
                    if(Parameter == NULL)FitFunction->SetParameters(450.,0.3,1,0.3,8.); // standard parameter optimize if necessary
                    else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4]);
                } else if (mesonType.CompareTo("Eta")==0){    
                    if(Parameter == NULL)FitFunction->SetParameters(450.,0.3,1,0.3,8.); // standard parameter optimize if necessary
                    else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4]);
                }
                FitFunction->SetParNames("Ae","Te","A","T","n");
                Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
            }
            return FitFunction;
        }


        if(ClassName.CompareTo("TGraphErrors") == 0){
            TGraphErrors *Obj = (TGraphErrors*)Obj_Dummy;

            if(type.BeginsWith("h") || type.BeginsWith("H")){
                cout <<Form("fitting %s with Hagedorn",FunctionName.Data()) << endl;
                TF1 *Hagedorn_Dummy = new TF1("Hagedorn_Dummy","[0]*pow([2]/([2]+x),[1])");
                FitFunction = (TF1*)Hagedorn_Dummy->Clone(FunctionName);
                FitFunction->SetRange(xmin, xmax);
                if(Parameter == NULL)FitFunction->SetParameters(1.,7.,0.37); // standard parameter optimize if necessary
                else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2]);
                FitFunction->SetParNames("C_{H}","n","p_{0} (GeV/c)");
                Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
            } else if(type.CompareTo("p")==0 || type.CompareTo("P")==0){
                cout <<Form("fitting %s with Powerlaw",FunctionName.Data()) << endl;
                TF1 *PowerLaw_Dummy = new TF1("PowerLaw_Dummy","[0] * 2 / TMath::Pi() * ([1]-1.)*([1]-2.)/pow([1]-3.,2) /x * pow(1+2*x/[2]/([1]-3),-[1])");
                FitFunction = (TF1*)PowerLaw_Dummy->Clone(FunctionName);
                FitFunction->SetRange(xmin, xmax);
                if(Parameter == NULL)FitFunction->SetParameters(2.,0.37,5.); // standard parameter optimize if necessary
                else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2]);
                FitFunction->SetParNames("A","p_{0}","n");
                Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
            } else if(type.BeginsWith("l") || type.BeginsWith("L")){
                cout <<Form("fitting %s with Levy",FunctionName.Data()) << endl;
                TF1 *Levy_Dummy = new TF1("Levy_Dummy",Form("[0]/(2.*TMath::Pi())*([1]-1.)*([1]-2.)/([1]*[2]*([1]*[2]+%.10f*([1]-2.)))*pow(1.+(sqrt(x*x+%.10f*%.10f)-%.10f)/([1]*[2]), -[1])",mass,mass,mass,mass));
                FitFunction = (TF1*)Levy_Dummy->Clone(FunctionName);
                FitFunction->SetRange(xmin, xmax);
                if(Parameter == NULL)FitFunction->SetParameters(2.,5.,0.18); // standard parameter optimize if necessary
                else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2]);
                if (fixingPar.CompareTo("fixn") == 0){
                    cout << "entered fixing" << endl;
                    FitFunction->FixParameter(1, Parameter[1]);
                }
                FitFunction->SetParNames("dN/dy","n","T_{Levy} (GeV/c)");
                Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
            } else if(type.BeginsWith("qcd") || type.BeginsWith("QCD")){
                cout <<Form("fitting %s with QCD",FunctionName.Data()) << endl;
                TF1 *QCD_Dummy = new TF1("QCD_Dummy","[0]*TMath::Power(x,-1*([1]+[2]/(TMath::Power(x,[3])+[4])))");
                FitFunction = (TF1*)QCD_Dummy->Clone(FunctionName);
                FitFunction->SetRange(xmin, xmax);
                if(Parameter == NULL)FitFunction->SetParameters(24,6.7,-6.5,1.,10); // standard parameter optimize if necessary
                else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4]);
                FitFunction->SetParNames("a","b","c","d","e");
                Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
            } else if(type.BeginsWith("oHag") || type.BeginsWith("OHag")){
                cout <<Form("fitting %s with ModHagedorn",FunctionName.Data()) << endl;
                TF1 *ModPowerLaw_Dummy = new TF1("ModHagedorn_Dummy","[0]*pow(exp(-[1]*x-abs([2])*x*x)+x/[3],-[4])");
                FitFunction = (TF1*)ModPowerLaw_Dummy->Clone(FunctionName);
                FitFunction->SetRange(xmin, xmax);
                if(Parameter == NULL)FitFunction->SetParameters(450.,0.37,0.2,0.7,8.); // standard parameter optimize if necessary
                else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4]);
                FitFunction->SetParNames("A (mbGeV^{-2}c^{3})","a [(GeV/c)^{-1}]","b [(GeV/c)^{-1}]","p_{0} (GeV/c)","n");
                Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
            } else if(type.BeginsWith("b") || type.BeginsWith("B")){
                cout <<Form("fitting %s with Boltzmann",FunctionName.Data()) << endl;
                TF1 *Boltzmann_Dummy = new TF1("Boltzmann_Dummy",Form("[0] *sqrt(x*x+%.10f*%.10f)* exp(- sqrt(x*x+%.10f*%.10f)/[1])",mass,mass,mass,mass));
                FitFunction = (TF1*)Boltzmann_Dummy->Clone(FunctionName);
                FitFunction->SetRange(xmin, xmax);
                if(Parameter == NULL)FitFunction->SetParameters(2.5,0.3); // standard parameter optimize if necessary
                else FitFunction->SetParameters(Parameter[0],Parameter[1]);
                FitFunction->SetParNames("C_{B}","T_{Boltzm.} (GeV/c)");
                Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
            } else if(type.BeginsWith("e") || type.BeginsWith("E")){
                cout <<Form("fitting %s with Exponential",FunctionName.Data()) << endl;
                TF1 *Exponential_Dummy = new TF1("Exponential_Dummy",Form("[0]/( 2 * TMath::Pi())/([1]*(%.10f+[1]))*exp(-(sqrt(x*x+%.10f*%.10f)-%.10f)/[1])",mass,mass,mass,mass));
                FitFunction = (TF1*)Exponential_Dummy->Clone(FunctionName);
                FitFunction->SetRange(xmin, xmax);
                if(Parameter == NULL)FitFunction->SetParameters(1.1,0.37); // standard parameter optimize if necessary
                else FitFunction->SetParameters(Parameter[0],Parameter[1]);
                FitFunction->SetParNames("C_{E}","T_{exp} (GeV/c)");
                Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
            } else if(type.BeginsWith("m") || type.BeginsWith("M")){
                cout <<Form("fitting %s with ModPowerlaw",FunctionName.Data()) << endl;
                TF1 *ModPowerLaw_Dummy = new TF1("ModPowerLaw_Dummy","[0]*pow((1 + (x)/[1]),-[2])");
                FitFunction = (TF1*)ModPowerLaw_Dummy->Clone(FunctionName);
                FitFunction->SetRange(xmin, xmax);
                if(Parameter == NULL)FitFunction->SetParameters(2.,0.37,5.); // standard parameter optimize if necessary
                else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2]);
                FitFunction->SetParNames("A","p_{0}","n");
                Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
            } else if(type.BeginsWith("6pol") || type.BeginsWith("6POL")){
                cout <<Form("fitting %s with Polynom order 6",FunctionName.Data()) << endl;
                TF1 *Pol6_Dummy = new TF1("Pol6_Dummy","[0]*pow(x,4)+[1]*pow(x,2)+[2]*x+[3]");
                FitFunction = (TF1*)Pol6_Dummy->Clone(FunctionName);
                FitFunction->SetRange(xmin, xmax);
                if(Parameter == NULL)FitFunction->SetParameters(1E-8,1E-3,0.1,0.); // standard parameter optimize if necessary
                else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3]);
                FitFunction->SetParNames("a","b","c","d");
                Obj->Fit(FitFunction,"RMNEX0+","",xmin,xmax);
            } else if(type.BeginsWith("tcm") || type.BeginsWith("TCM")){ // Two component model fit [A. Bylinkin and A. Rostovtsev, Phys. Atom. Nucl 75 (2012) 999-1005]
                cout <<Form("fitting %s with two component model by Bylinkin",FunctionName.Data()) << endl;
                TF1 *TwoCompModel_Dummy = new TF1("twoCompModel_Dummy",Form("[0]*exp(-(TMath::Sqrt(x*x+%.10f*%.10f)-%.10f)/[1]) + [2]/(TMath::Power(1+x*x/([3]*[3]*[4]),[4]) )",mass,mass,mass));
                FitFunction = (TF1*)TwoCompModel_Dummy->Clone(FunctionName);
                FitFunction->SetRange(xmin, xmax);
                if (mesonType.CompareTo("Pi0")==0){
                    if(Parameter == NULL){
                        FitFunction->SetParameters(450.,0.3,1,0.3,8.); // standard parameter optimize if necessary
                    } else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4]);
                } else if (mesonType.CompareTo("Eta")==0){
                    if(Parameter == NULL){
                        FitFunction->SetParameters(450.,0.3,1,0.3,8.); // standard parameter optimize if necessary
                    } else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4]);
                }
                FitFunction->SetParNames("Ae","Te","A","T","n");
                Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
            } else if(type.BeginsWith("tcmpt") || type.BeginsWith("TCMPT")){ // Two component model fit [A. Bylinkin and A. Rostovtsev, Phys. Atom. Nucl 75 (2012) 999-1005]
                cout <<Form("fitting %s with two component model by Bylinkin multiplied by pT",FunctionName.Data()) << endl;
                TF1 *TwoCompModel_Dummy = new TF1("twoCompModel_Dummy",Form(" x*[0]*exp(-(TMath::Sqrt(x*x+%.10f*%.10f)-%.10f)/[1]) + x*[2]/(TMath::Power(1+x*x/([3]*[3]*[4]),[4]) )",mass,mass,mass));
                FitFunction = (TF1*)TwoCompModel_Dummy->Clone(FunctionName);
                FitFunction->SetRange(xmin, xmax);
                if (mesonType.CompareTo("Pi0")==0){
                    if(Parameter == NULL){
                        FitFunction->SetParameters(450.,0.3,1,0.3,8.); // standard parameter optimize if necessary
                    } else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4]);
                } else if (mesonType.CompareTo("Eta")==0){
                    if(Parameter == NULL){
                        FitFunction->SetParameters(450.,0.3,1,0.3,8.); // standard parameter optimize if necessary
                    } else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4]);
                }
                FitFunction->SetParNames("Ae","Te","A","T","n");
                Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
            } else if(type.BeginsWith("rad") || type.BeginsWith("RAD")){
                cout <<Form("fitting %s with Radoslav Func",FunctionName.Data()) << endl;
                TF1 *Rad_Dummy = new TF1("Rad_Dummy",Form("(x<=[3])*x*[0]*([1]-1.)*([1]-2.)/([1]*[2]*([1]*[2]+%.10f*([1]-2.)))*pow(1.+(sqrt(x*x+%.10f*%.10f)-%.10f)/([1]*[2]), -[1]) + (x>[3]) * x *[0] * ([1]-1.)*([1]-2.)/([1]*[2]*([1]*[2]+%.10f*([1]-2.))) * pow([1]*[2]/([3]+[1]*[2]-%.10f),[1]) * pow([3],[4]) * pow( 1./x, [4] )",mass,mass,mass,mass,mass,mass));
                FitFunction = (TF1*)Rad_Dummy->Clone(FunctionName);
                FitFunction->SetRange(xmin, xmax);
                Double_t par0min3 = 150.;
                Double_t par0max3 = 1500.;
                Double_t par1min3 = 4.;
                Double_t par1max3 = 40.;
                Double_t par2min3 = 0.07;
                Double_t par2max3 = .7;
                Double_t par3min3 = 4.36;
                Double_t par3max3 = 4.40;
                Double_t par4min3 = 2.;
                Double_t par4max3 = 18.;
                if(Parameter != NULL){
                    par0min3 = Parameter[0];
                    par0max3 = Parameter[1];
                    par1min3 = Parameter[2];
                    par1max3 = Parameter[3];
                    par2min3 = Parameter[4];
                    par2max3 = Parameter[5];
                    par3min3 = Parameter[6];
                    par3max3 = Parameter[7];
                    par4min3 = Parameter[8];
                    par4max3 = Parameter[9];
                }
                FitFunction->SetParLimits(0, par0min3, par0max3);
                FitFunction->SetParLimits(1, par1min3, par1max3);
                FitFunction->SetParLimits(2, par2min3, par2max3);
                FitFunction->SetParLimits(3, par3min3, par3max3);
                FitFunction->SetParLimits(4, par4min3, par4max3);
                FitFunction->SetParNames("a","b","c","d");
                Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
                return FitFunction;
            }
        }

        if(ClassName.CompareTo("TGraphAsymmErrors") == 0){
            TGraphAsymmErrors *Obj = (TGraphAsymmErrors*)Obj_Dummy;
            if(type.CompareTo("h")==0 || type.CompareTo("H")==0){
                cout <<Form("fitting %s with Hagedorn",FunctionName.Data()) << endl;
                TF1 *Hagedorn_Dummy = new TF1("Hagedorn_Dummy","[0]*pow([2]/([2]+x),[1])");
                FitFunction = (TF1*)Hagedorn_Dummy->Clone(FunctionName);
                FitFunction->SetRange(xmin, xmax);
                if(Parameter == NULL) FitFunction->SetParameters(19.,6.8,0.84); // standard parameter optimize if necessary
                else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2]);
                FitFunction->SetParNames("C_{H}","n","p_{0} (GeV/c)");
                Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
            } else if(type.CompareTo("powPure")==0 || type.CompareTo("PowPure")==0){
                cout <<Form("fitting %s with Pure Powerlaw",FunctionName.Data()) << endl;
                TF1 *PowerLaw_Dummy = new TF1("PowerLawPure_Dummy","[0] * 1/pow(x,[1])");
                FitFunction = (TF1*)PowerLaw_Dummy->Clone(FunctionName);
                FitFunction->SetRange(xmin, xmax);
                if(Parameter == NULL)FitFunction->SetParameters(2.,5.); // standard parameter optimize if necessary
                else FitFunction->SetParameters(Parameter[0],Parameter[1]);
                FitFunction->SetParNames("A_{pow}","n");
                Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
            } else if(type.CompareTo("p")==0 || type.CompareTo("P")==0){
                cout <<Form("fitting %s with Powerlaw",FunctionName.Data()) << endl;
                TF1 *PowerLaw_Dummy = new TF1("PowerLaw_Dummy","[0] * 2 / TMath::Pi() * ([1]-1.)*([1]-2.)/pow([1]-3.,2) /x * pow(1+2*x/[2]/([1]-3),-[1])");
                FitFunction = (TF1*)PowerLaw_Dummy->Clone(FunctionName);
                FitFunction->SetRange(xmin, xmax);
                if(Parameter == NULL)FitFunction->SetParameters(2.,5.,0.37); // standard parameter optimize if necessary
                else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2]);
                FitFunction->SetParNames("A_{pow}","n","p_{0}");
                Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
            } else if(type.CompareTo("l") ==0|| type.CompareTo("L")==0){
                cout <<Form("fitting %s with Levy",FunctionName.Data()) << endl;
                TF1 *Levy_Dummy = new TF1("Levy_Dummy",Form("[0] / ( 2 * TMath::Pi())*([1]-1.)*([1]-2.) / ([1]*[2]*([1]*[2]+%.10f*([1]-2.)))  * pow(1.+(sqrt(x*x+%.10f*%.10f)-%.10f)/([1]*[2]), -[1])",mass,mass,mass,mass));
                FitFunction = (TF1*)Levy_Dummy->Clone(FunctionName);
                FitFunction->SetRange(xmin, xmax);
                if(Parameter == NULL)FitFunction->SetParameters(2.,5.,0.18); // standard parameter optimize if necessary
                else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2]);
                if (fixingPar.CompareTo("fixn") == 0){
                    cout << "entered fixing" << endl;
                    FitFunction->FixParameter(1, Parameter[1]);
                }
                FitFunction->SetParNames("dN/dy","n","T_{Levy} (GeV/c)");
                Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
            } else if(type.BeginsWith("qcd") || type.BeginsWith("QCD")){
                cout <<Form("fitting %s with QCD",FunctionName.Data()) << endl;
                TF1 *QCD_Dummy = new TF1("QCD_Dummy","[0]*TMath::Power(x,-1*([1]+[2]/(TMath::Power(x,[3])+[4])))");
                FitFunction = (TF1*)QCD_Dummy->Clone(FunctionName);
                FitFunction->SetRange(xmin, xmax);
                if(Parameter == NULL)FitFunction->SetParameters(24,6.7,-6.5,1.,10); // standard parameter optimize if necessary
                else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4]);
                FitFunction->SetParNames("a","b","c","d","e");
                Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
            } else if(type.BeginsWith("oHag") || type.BeginsWith("OHag")){
                cout <<Form("fitting %s with ModHagedorn",FunctionName.Data()) << endl;
                TF1 *ModPowerLaw_Dummy = new TF1("ModHagedorn_Dummy","[0]*pow(exp(-[1]*x-abs([2])*x*x)+x/[3],-[4])");
                FitFunction = (TF1*)ModPowerLaw_Dummy->Clone(FunctionName);
                FitFunction->SetRange(xmin, xmax);
                if(Parameter == NULL)FitFunction->SetParameters(450.,0.37,0.2,0.7,8.); // standard parameter optimize if necessary
                else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4]);
                FitFunction->SetParNames("A (mbGeV^{-2}c^{3})","a [(GeV/c)^{-1}]","b [(GeV/c)^{-1}]","p_{0} (GeV/c)","n");
                Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
            } else if(type.CompareTo("b")==0 || type.CompareTo("B")==0){
                cout <<Form("fitting %s with Boltzmann",FunctionName.Data()) << endl;
                TF1 *Boltzmann_Dummy =  new TF1("Boltzmann_Dummy",Form("[0] *sqrt(x*x+%.10f*%.10f)* exp(- sqrt(x*x+%.10f*%.10f)/[1])",mass,mass,mass,mass));
                FitFunction = (TF1*)Boltzmann_Dummy->Clone(FunctionName);
                FitFunction->SetRange(xmin, xmax);
                if(Parameter == NULL)FitFunction->SetParameters(2.5,0.3); // standard parameter optimize if necessary
                else FitFunction->SetParameters(Parameter[0],Parameter[1]);
                FitFunction->SetParNames("C_{B}","T_{Boltzm.} (GeV/c)");
                Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
            } else if(type.CompareTo("e")==0 || type.CompareTo("E")==0){
                cout <<Form("fitting %s with Exponential",FunctionName.Data()) << endl;
                TF1 *Exponential_Dummy = new TF1("Exponential_Dummy",Form("[0]/( 2 * TMath::Pi())/([1]*(%.10f+[1]))*exp(-(sqrt(x*x+%.10f*%.10f)-%.10f)/[1])",mass,mass,mass,mass));
                FitFunction = (TF1*)Exponential_Dummy->Clone(FunctionName);
                FitFunction->SetRange(xmin, xmax);
                if(Parameter == NULL)FitFunction->SetParameters(1.1,0.37); // standard parameter optimize if necessary
                else FitFunction->SetParameters(Parameter[0],Parameter[1]);
                FitFunction->SetParNames("C_{E}","T_{exp} (GeV/c)");
                Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
            } else if(type.CompareTo("m")==0 || type.CompareTo("M")==0){
                cout <<Form("fitting %s with ModPowerlaw",FunctionName.Data()) << endl;
                TF1 *ModPowerLaw_Dummy = new TF1("ModPowerLaw_Dummy","[0]*pow((1 + x/[1]),-[2])");
                FitFunction = (TF1*)ModPowerLaw_Dummy->Clone(FunctionName);
                FitFunction->SetRange(xmin, xmax);
                if(Parameter == NULL)FitFunction->SetParameters(2.,0.37,5.); // standard parameter optimize if necessary
                else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2]);
                FitFunction->SetParNames("A","p_{0}","n");
                Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
            } else if(type.CompareTo("tmpt") ==0|| type.CompareTo("TMPT")==0){
                cout <<Form("fitting %s with Levy",FunctionName.Data()) << endl;
                TF1 *Levy_Dummy = new TF1("Tsallis_Dummy",Form("[0] / ( 2 * TMath::Pi())*([1]-1.)*([1]-2.) / ([1]*[2]*([1]*[2]+%.10f*([1]-2.)))  * x* pow(1.+(sqrt(x*x+%.10f*%.10f)-%.10f)/([1]*[2]), -[1])",mass,mass,mass,mass));
                FitFunction = (TF1*)Levy_Dummy->Clone(FunctionName);
                FitFunction->SetRange(xmin, xmax);
                if(Parameter == NULL)FitFunction->SetParameters(2.,5.,0.18); // standard parameter optimize if necessary
                else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2]);
                if (fixingPar.CompareTo("fixn") == 0){
                    cout << "entered fixing" << endl;
                    FitFunction->FixParameter(1, Parameter[1]);
                }
                FitFunction->SetParNames("dN/dy","n","T_{Levy} (GeV/c)");
                Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
            } else if(type.BeginsWith("tcmmod") || type.BeginsWith("TCMMOD")){ // Two component model fit [A. Bylinkin and A. Rostovtsev, Phys. Atom. Nucl 75 (2012) 999-1005]
                cout <<Form("fitting %s with orginal two component model by Bylinkin",FunctionName.Data()) << endl;
                TF1 *TwoCompModel_Dummy = new TF1("twoCompModel_Dummy",Form("[0]*exp(-(TMath::Sqrt(x*x+%.10f*%.10f)-%.10f)/[1]) + [2]/(TMath::Power(1+x*x/([3]*[3]*[4]),[4]) )",mass,mass,mass));
                FitFunction = (TF1*)TwoCompModel_Dummy->Clone(FunctionName);
                FitFunction->SetRange(xmin, xmax);
                if (mesonType.CompareTo("Pi0")==0){
                    if(Parameter == NULL)FitFunction->SetParameters(450.,0.3,1,0.3,8.); // standard parameter optimize if necessary
                    else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4]);
                } else if (mesonType.CompareTo("Eta")==0){
                    if(Parameter == NULL)FitFunction->SetParameters(450.,0.3,1,0.3,8.); //TGraphAsymmErrors
                    else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4]);
                }
                FitFunction->SetParNames("Ae","Te","A","T","n");
                if(mesonType.CompareTo("Pi0")==0 && limitPar){
                    FitFunction->SetParLimits(0,0,10000);
                    FitFunction->SetParLimits(2,0,10000);
                } else {
                    if(FunctionName.Contains("0010") && limitPar){
                        FitFunction->SetParLimits(0,0,100);
                        FitFunction->SetParLimits(2,0,100);   
                    } else if(FunctionName.Contains("2050") && limitPar){
                        FitFunction->SetParLimits(0,0,100);
                        FitFunction->SetParLimits(2,0,100);   
                    }
                }
                if (limitPar){
                    FitFunction->SetParLimits(1,0,10000);
                    FitFunction->SetParLimits(3,0,10000);
                    FitFunction->SetParLimits(4,0,10);
                }    
                Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
            } else if(type.BeginsWith("tcmpt") || type.BeginsWith("TCMPT")){ // Two component model fit [A. Bylinkin and A. Rostovtsev, Phys. Atom. Nucl 75 (2012) 999-1005]
                cout <<Form("fitting %s with two component model by Bylinkin multiplied by pT",FunctionName.Data()) << endl;
                TF1 *TwoCompModel_Dummy = new TF1("twoCompModel_Dummy",Form(" x*[0]*exp(-(TMath::Sqrt(x*x+%.10f*%.10f)-%.10f)/[1]) + x*[2]/(TMath::Power(1+x*x/([3]*[3]*[4]),[4]) )",mass,mass,mass));
                FitFunction = (TF1*)TwoCompModel_Dummy->Clone(FunctionName);
                FitFunction->SetRange(xmin, xmax);
                if (mesonType.CompareTo("Pi0")==0){
                    if(Parameter == NULL)FitFunction->SetParameters(450.,0.3,1,0.3,8.); // standard parameter optimize if necessary
                    else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4]);
                } else if (mesonType.CompareTo("Eta")==0){
                    if(Parameter == NULL)FitFunction->SetParameters(450.,0.3,1,0.3,8.); //(450.,0.01,1,0.3,.5); //
                    else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4]);
                }
                FitFunction->SetParNames("Ae","Te","A","T","n");
                if(mesonType.CompareTo("Pi0")==0 && limitPar){
                    FitFunction->SetParLimits(0,0,10000);
                    FitFunction->SetParLimits(2,0,10000);
                } else if (limitPar){
                    FitFunction->SetParLimits(0,0,100);
                    FitFunction->SetParLimits(2,0,100);
                }
                if (limitPar){
                    FitFunction->SetParLimits(1,0,10000);
                    FitFunction->SetParLimits(3,0,10000);
                    FitFunction->SetParLimits(4,0,10);
                }
                Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
            } else if(type.BeginsWith("tcm") || type.BeginsWith("TCM")){ // Two component model fit [A. Bylinkin and A. Rostovtsev, Phys. Atom. Nucl 75 (2012) 999-1005]
                cout <<Form("fitting %s with two component model by Bylinkin",FunctionName.Data()) << endl;
                TF1 *TwoCompModel_Dummy = new TF1("twoCompModel_Dummy",Form("[0]*exp(-(TMath::Sqrt(x*x+%.10f*%.10f)-%.10f)/[1]) + [2]/(TMath::Power(1+x*x/([3]*[3]*[4]),[4]) )",mass,mass,mass));
                FitFunction = (TF1*)TwoCompModel_Dummy->Clone(FunctionName);
                FitFunction->SetRange(xmin, xmax);
                if (mesonType.CompareTo("Pi0")==0){
                    if(Parameter == NULL)FitFunction->SetParameters(450.,0.3,1,0.3,8.); // standard parameter optimize if necessary
                    else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4]);
                } else if (mesonType.CompareTo("Eta")==0){
                    if(Parameter == NULL)FitFunction->SetParameters(450.,0.3,1,0.3,8.); //TGraphAsymmErrors
                    else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4]);
                }
                FitFunction->SetParNames("Ae","Te","A","T","n");
                if(mesonType.CompareTo("Pi0")==0 && limitPar){
                    if(FunctionName.Contains("0010")){
                    FitFunction->SetParLimits(0,1.,1e3);
    //                   FitFunction->SetParLimits(2,1.,1e3);
    //                   FitFunction->FixParameter(0,160.);
                    FitFunction->FixParameter(2,840.);
                    } else if(FunctionName.Contains("2050")){
                    FitFunction->SetParLimits(0,1.,1e2);
    //                   FitFunction->SetParLimits(2,1.,1e2);
    //                   FitFunction->FixParameter(0,30.);
                    FitFunction->FixParameter(2,80.);
                    } else {
                    FitFunction->SetParLimits(0,0,10000);
                    FitFunction->SetParLimits(2,0,10000);
                    }
                } else if(mesonType.CompareTo("Eta")==0 && limitPar){
                    if(FunctionName.Contains("0010")){
                        FitFunction->SetParLimits(0,0.01,1e2);
    //                     FitFunction->SetParLimits(2,1.,1e2);
    //                     FitFunction->FixParameter(0,15.);
                        FitFunction->FixParameter(2,100.);
                    } else if(FunctionName.Contains("2050")){
                        FitFunction->SetParLimits(0,0.01,1e2);
    //                     FitFunction->SetParLimits(2,1.,1e2);
    //                     FitFunction->FixParameter(0,5.);
                        FitFunction->FixParameter(2,2.);
                    }
                } else if(limitPar){
                    FitFunction->SetParLimits(1,0,10000);
                    FitFunction->SetParLimits(3,0,10000);
                    FitFunction->SetParLimits(4,0,10);
                }    
                Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
            } else if(type.BeginsWith("Lowtcm") || type.BeginsWith("LOWTCM")){ // Two component model fit [A. Bylinkin and A. Rostovtsev, Phys. Atom. Nucl 75 (2012) 999-1005]
                cout <<Form("fitting %s with low pt component of Bylinkin",FunctionName.Data()) << endl;
                TF1 *TwoCompModel_Dummy = new TF1("twoCompModel_Dummy",Form("[0]*exp(-(TMath::Sqrt(x*x+%.10f*%.10f)-%.10f)/[1])",mass,mass,mass));
                FitFunction = (TF1*)TwoCompModel_Dummy->Clone(FunctionName);
                FitFunction->SetRange(xmin, xmax);
                if (mesonType.CompareTo("Pi0")==0){
                    if(Parameter == NULL){
                        FitFunction->SetParameters(450.,0.3); // standard parameter optimize if necessary
                    } else FitFunction->SetParameters(Parameter[0],Parameter[1]);
                } else if (mesonType.CompareTo("Eta")==0){
                    if(Parameter == NULL){
                        FitFunction->SetParameters(450.,0.3); // standard parameter optimize if necessary
                    } else FitFunction->SetParameters(Parameter[0],Parameter[1]);
                }
                FitFunction->SetParNames("Ae","Te");
                Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
            } else if(type.BeginsWith("Hightcm") || type.BeginsWith("HIGHTCM")){ // Two component model fit [A. Bylinkin and A. Rostovtsev, Phys. Atom. Nucl 75 (2012) 999-1005]
                cout <<Form("fitting %s with high pt component of Bylinkin",FunctionName.Data()) << endl;
                TF1 *TwoCompModel_Dummy = new TF1("twoCompModel_Dummy",Form("[0]/(TMath::Power(1+x*x/([1]*[1]*[2]),[2]) )"));
                FitFunction = (TF1*)TwoCompModel_Dummy->Clone(FunctionName);
                FitFunction->SetRange(xmin, xmax);
                if (mesonType.CompareTo("Pi0")==0){
                    if(Parameter == NULL){
                        FitFunction->SetParameters(1,0.3,8.); // standard parameter optimize if necessary
                    } else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2]);
                } else if (mesonType.CompareTo("Eta")==0){
                    if(Parameter == NULL){
                        FitFunction->SetParameters(1,0.3,8.); // standard parameter optimize if necessary
                    } else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2]);
                }
                FitFunction->SetParNames("A","T","n");
                Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
            } else if(type.CompareTo("rad") ==0|| type.CompareTo("RAD")==0){
                cout <<Form("fitting %s with Radoslav Func",FunctionName.Data()) << endl;
                TF1 *Rad_Dummy = new TF1("Rad_Dummy",Form("(x<=[3])*x*[0]*([1]-1.)*([1]-2.)/([1]*[2]*([1]*[2]+%.10f*([1]-2.)))*pow(1.+(sqrt(x*x+%.10f*%.10f)-%.10f)/([1]*[2]), -[1]) + (x>[3]) * x *[0] * ([1]-1.)*([1]-2.)/([1]*[2]*([1]*[2]+%.10f*([1]-2.))) * pow([1]*[2]/([3]+[1]*[2]-%.10f),[1]) * pow([3],[4]) * pow( 1./x, [4] )",mass,mass,mass,mass,mass,mass));
                FitFunction = (TF1*)Rad_Dummy->Clone(FunctionName);
                FitFunction->SetRange(xmin, xmax);
                Double_t par0min3 = 150.;
                Double_t par0max3 = 1500.;
                Double_t par1min3 = 4.;
                Double_t par1max3 = 40.;
                Double_t par2min3 = 0.07;
                Double_t par2max3 = .7;
                Double_t par3min3 = 4.36;
                Double_t par3max3 = 4.40;
                Double_t par4min3 = 2.;
                Double_t par4max3 = 18.;
                if(Parameter != NULL){
                    par0min3 = Parameter[0];
                    par0max3 = Parameter[1];
                    par1min3 = Parameter[2];
                    par1max3 = Parameter[3];
                    par2min3 = Parameter[4];
                    par2max3 = Parameter[5];
                    par3min3 = Parameter[6];
                    par3max3 = Parameter[7];
                    par4min3 = Parameter[8];
                    par4max3 = Parameter[9];
                }
                FitFunction->SetParLimits(0, par0min3, par0max3);
                FitFunction->SetParLimits(1, par1min3, par1max3);
                FitFunction->SetParLimits(2, par2min3, par2max3);
                FitFunction->SetParLimits(3, par3min3, par3max3);
                FitFunction->SetParLimits(4, par4min3, par4max3);
                FitFunction->SetParNames("a","b","c","d");
                Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
                return FitFunction;
            } else if(type.CompareTo("ptredt")==0 || type.CompareTo("ptredt")==0){
                cout <<Form("fitting %s with Tsallis (*Pt)",FunctionName.Data()) << endl;
                TF1 *Levy_Dummy = new TF1("Levy_Dummy",Form("[0] / ( 2 * TMath::Pi())*([1]-1.)*([1]-2.) / ([1]*[2]*([1]*[2]+%.10f*([1]-2.)))  * pow(1.+(sqrt(x*x+%.10f*%.10f)-%.10f)/([1]*[2]), -[1])",mass,mass,mass,mass)); 
                FitFunction = (TF1*)Levy_Dummy->Clone(FunctionName);
                FitFunction->SetRange(xmin, xmax);
                if(Parameter == NULL)FitFunction->SetParameters(2.,5.,0.18); // standard parameter optimize if necessary
                else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2]);
                if (fixingPar.CompareTo("fixn") == 0){
                    cout << "entered fixing" << endl;
                    FitFunction->FixParameter(1, Parameter[1]);
                }
                FitFunction->SetParNames("dN/dy","n","T_{Levy} (GeV/c)");
                Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
            }
        
        }

        if(ClassName.CompareTo("TGraph") == 0){
            TGraph *Obj = (TGraph*)Obj_Dummy;
            if(type.BeginsWith("h") || type.BeginsWith("H")){
                cout <<Form("fitting %s with Hagedorn",FunctionName.Data()) << endl;
                TF1 *Hagedorn_Dummy = new TF1("Hagedorn_Dummy","[0]*pow([2]/([2]+x),[1])");
                FitFunction = (TF1*)Hagedorn_Dummy->Clone(FunctionName);
                FitFunction->SetRange(xmin, xmax);
                if(Parameter == NULL) FitFunction->SetParameters(19.,6.8,0.84); // standard parameter optimize if necessary
                else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2]);
                FitFunction->SetParNames("C_{H}","n","p_{0} (GeV/c)");
                Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
            } else if(type.BeginsWith("p") || type.BeginsWith("P")){
                cout <<Form("fitting %s with Powerlaw",FunctionName.Data()) << endl;
                TF1 *PowerLaw_Dummy = new TF1("PowerLaw_Dummy","[0] * 2 / TMath::Pi() * ([1]-1.)*([1]-2.)/pow([1]-3.,2) /x * pow(1+2*x/[2]/([1]-3),-[1])");
                FitFunction = (TF1*)PowerLaw_Dummy->Clone(FunctionName);
                FitFunction->SetRange(xmin, xmax);
                if(Parameter == NULL)FitFunction->SetParameters(2.,5.,0.37); // standard parameter optimize if necessary
                else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2]);
                FitFunction->SetParNames("A_{pow}","n","p_{0}");
                Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
            } else if(type.BeginsWith("l") || type.BeginsWith("L")){
                cout <<Form("fitting %s with Levy",FunctionName.Data()) << endl;
                TF1 *Levy_Dummy = new TF1("Levy_Dummy",Form("[0] / ( 2 * TMath::Pi())*([1]-1.)*([1]-2.) / ([1]*[2]*([1]*[2]+%.10f*([1]-2.)))  * pow(1.+(sqrt(x*x+%.10f*%.10f)-%.10f)/([1]*[2]), -[1])",mass,mass,mass,mass));
                FitFunction = (TF1*)Levy_Dummy->Clone(FunctionName);
                FitFunction->SetRange(xmin, xmax);
                if(Parameter == NULL)FitFunction->SetParameters(2.,5.,0.18); // standard parameter optimize if necessary
                else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2]);
                if (fixingPar.CompareTo("fixn") == 0){
                    cout << "entered fixing" << endl;
                    FitFunction->FixParameter(1, Parameter[1]);
                }
                FitFunction->SetParNames("dN/dy","n","T_{Levy} (GeV/c)");
                Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
            } else if(type.BeginsWith("qcd") || type.BeginsWith("QCD")){
                cout <<Form("fitting %s with QCD",FunctionName.Data()) << endl;
                TF1 *QCD_Dummy = new TF1("QCD_Dummy","[0]*TMath::Power(x,-1*([1]+[2]/(TMath::Power(x,[3])+[4])))");
                FitFunction = (TF1*)QCD_Dummy->Clone(FunctionName);
                FitFunction->SetRange(xmin, xmax);
                if(Parameter == NULL)FitFunction->SetParameters(24,6.7,-6.5,1.,10); // standard parameter optimize if necessary
                else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4]);
                FitFunction->SetParNames("a","b","c","d","e");
                Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
            } else if(type.BeginsWith("oHag") || type.BeginsWith("OHag")){
                cout <<Form("fitting %s with ModHagedorn",FunctionName.Data()) << endl;
                TF1 *ModPowerLaw_Dummy = new TF1("ModHagedorn_Dummy","[0]*pow(exp(-[1]*x-abs([2])*x*x)+x/[3],-[4])");
                FitFunction = (TF1*)ModPowerLaw_Dummy->Clone(FunctionName);
                FitFunction->SetRange(xmin, xmax);
                if(Parameter == NULL)FitFunction->SetParameters(450.,0.37,0.2,0.7,8.); // standard parameter optimize if necessary
                else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4]);
                FitFunction->SetParNames("A (mbGeV^{-2}c^{3})","a [(GeV/c)^{-1}]","b [(GeV/c)^{-1}]","p_{0} (GeV/c)","n");
                Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
            } else if(type.BeginsWith("b") || type.BeginsWith("B")){
                cout <<Form("fitting %s with Boltzmann",FunctionName.Data()) << endl;
                TF1 *Boltzmann_Dummy =  new TF1("Boltzmann_Dummy",Form("[0] *sqrt(x*x+%.10f*%.10f)* exp(- sqrt(x*x+%.10f*%.10f)/[1])",mass,mass,mass,mass));
                FitFunction = (TF1*)Boltzmann_Dummy->Clone(FunctionName);
                FitFunction->SetRange(xmin, xmax);
                if(Parameter == NULL)FitFunction->SetParameters(2.5,0.3); // standard parameter optimize if necessary
                else FitFunction->SetParameters(Parameter[0],Parameter[1]);
                FitFunction->SetParNames("C_{B}","T_{Boltzm.} (GeV/c)");
                Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
            } else if(type.BeginsWith("e") || type.BeginsWith("E")){
                cout <<Form("fitting %s with Exponential",FunctionName.Data()) << endl;
                TF1 *Exponential_Dummy = new TF1("Exponential_Dummy",Form("[0]/( 2 * TMath::Pi())/([1]*(%.10f+[1]))*exp(-(sqrt(x*x+%.10f*%.10f)-%.10f)/[1])",mass,mass,mass,mass));
                FitFunction = (TF1*)Exponential_Dummy->Clone(FunctionName);
                FitFunction->SetRange(xmin, xmax);
                if(Parameter == NULL)FitFunction->SetParameters(1.1,0.37); // standard parameter optimize if necessary
                else FitFunction->SetParameters(Parameter[0],Parameter[1]);
                FitFunction->SetParNames("C_{E}","T_{exp} (GeV/c)");
                Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
            } else if(type.BeginsWith("m") || type.BeginsWith("M")){
                cout <<Form("fitting %s with ModPowerlaw",FunctionName.Data()) << endl;
                TF1 *ModPowerLaw_Dummy = new TF1("ModPowerLaw_Dummy","[0]*pow((1 + (x)/[1]),-[2])");
                FitFunction = (TF1*)ModPowerLaw_Dummy->Clone(FunctionName);
                FitFunction->SetRange(xmin, xmax);
                if(Parameter == NULL)FitFunction->SetParameters(2.,0.37,5.); // standard parameter optimize if necessary
                else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2]);
                FitFunction->SetParNames("A","p_{0}","n");
                Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
            } else if(type.BeginsWith("tcm") || type.BeginsWith("TCM")){ // Two component model fit [A. Bylinkin and A. Rostovtsev, Phys. Atom. Nucl 75 (2012) 999-1005]
                cout <<Form("fitting %s with two component model by Bylinkin",FunctionName.Data()) << endl;
                TF1 *TwoCompModel_Dummy = new TF1("twoCompModel_Dummy",Form("[0]*exp(-(TMath::Sqrt(x*x+%.10f*%.10f)-%.10f)/[1]) + [2]/(TMath::Power(1+x*x/([3]*[3]*[4]),[4]) )",mass,mass,mass));
                FitFunction = (TF1*)TwoCompModel_Dummy->Clone(FunctionName);
                FitFunction->SetRange(xmin, xmax);
                if (mesonType.CompareTo("Pi0")==0){
                    if(Parameter == NULL)FitFunction->SetParameters(450.,0.3,1,0.3,8.); // standard parameter optimize if necessary
                    else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4]);
                } else if (mesonType.CompareTo("Eta")==0){ 
                    if(Parameter == NULL)FitFunction->SetParameters(450.,0.3,1,0.3,8.); // 450.,0.01,1,0.3,.5)
                    else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4]);
                }
                if (limitPar){
                    FitFunction->SetParLimits(1,0,10000);
                    FitFunction->SetParLimits(3,0,10000);
                    FitFunction->SetParLimits(4,0,10);
                }    
                FitFunction->SetParNames("Ae","Te","A","T","n");
                Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
            } else if(type.BeginsWith("tcmpt") || type.BeginsWith("TCMPT")){ // Two component model fit [A. Bylinkin and A. Rostovtsev, Phys. Atom. Nucl 75 (2012) 999-1005]
                cout <<Form("fitting %s with two component model by Bylinkin multiplied by pT",FunctionName.Data()) << endl;
                TF1 *TwoCompModel_Dummy = new TF1("twoCompModel_Dummy",Form(" x*[0]*exp(-(TMath::Sqrt(x*x+%.10f*%.10f)-%.10f)/[1]) + x*[2]/(TMath::Power(1+x*x/([3]*[3]*[4]),[4]) )",mass,mass,mass));
                FitFunction = (TF1*)TwoCompModel_Dummy->Clone(FunctionName);
                FitFunction->SetRange(xmin, xmax);
                if (mesonType.CompareTo("Pi0")==0){
                    if(Parameter == NULL)FitFunction->SetParameters(450.,0.3,1,0.3,8.); // standard parameter optimize if necessary
                    else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4]);
                } else if (mesonType.CompareTo("Eta")==0){
                    if(Parameter == NULL)FitFunction->SetParameters(450.,0.3,1,0.3,8.); 
                    else FitFunction->SetParameters(Parameter[0],Parameter[1],Parameter[2],Parameter[3],Parameter[4]);
                }
                if (limitPar){
                    FitFunction->SetParLimits(1,0,10000);
                    FitFunction->SetParLimits(3,0,10000);
                    FitFunction->SetParLimits(4,0,10);
                }    
                FitFunction->SetParNames("Ae","Te","A","T","n");
                Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
            } else if(type.BeginsWith("rad") || type.BeginsWith("RAD")){
                cout <<Form("fitting %s with Radoslav Func",FunctionName.Data()) << endl;
                TF1 *Rad_Dummy = new TF1("Rad_Dummy",Form("(x<=[3])*x*[0]*([1]-1.)*([1]-2.)/([1]*[2]*([1]*[2]+%.10f*([1]-2.)))*pow(1.+(sqrt(x*x+%.10f*%.10f)-%.10f)/([1]*[2]), -[1]) + (x>[3]) * x *[0] * ([1]-1.)*([1]-2.)/([1]*[2]*([1]*[2]+%.10f*([1]-2.))) * pow([1]*[2]/([3]+[1]*[2]-%.10f),[1]) * pow([3],[4]) * pow( 1./x, [4] )",mass,mass,mass,mass,mass,mass));
                FitFunction = (TF1*)Rad_Dummy->Clone(FunctionName);
                FitFunction->SetRange(xmin, xmax);
                Double_t par0min3 = 150.;
                Double_t par0max3 = 1500.;
                Double_t par1min3 = 4.;
                Double_t par1max3 = 40.;
                Double_t par2min3 = 0.07;
                Double_t par2max3 = .7;
                Double_t par3min3 = 4.36;
                Double_t par3max3 = 4.40;
                Double_t par4min3 = 2.;
                Double_t par4max3 = 18.;
                if(Parameter != NULL){
                    par0min3 = Parameter[0];
                    par0max3 = Parameter[1];
                    par1min3 = Parameter[2];
                    par1max3 = Parameter[3];
                    par2min3 = Parameter[4];
                    par2max3 = Parameter[5];
                    par3min3 = Parameter[6];
                    par3max3 = Parameter[7];
                    par4min3 = Parameter[8];
                    par4max3 = Parameter[9];
                }
                FitFunction->SetParLimits(0, par0min3, par0max3);
                FitFunction->SetParLimits(1, par1min3, par1max3);
                FitFunction->SetParLimits(2, par2min3, par2max3);
                FitFunction->SetParLimits(3, par3min3, par3max3);
                FitFunction->SetParLimits(4, par4min3, par4max3);
                FitFunction->SetParNames("a","b","c","d");
                Obj->Fit(FitFunction,FitOptions,"",xmin,xmax);
                return FitFunction;
            }
        }
        return FitFunction;
    }

    void SetParametersLimitsForFit (TF1* fit, Int_t nParams, Double_t parameters[] = NULL){
        for (Int_t i = 0; i < nParams; i++){
            fit->SetParLimits(i, parameters[i*2], parameters[i*2+1]);
            fit->SetParameter(i, parameters[i+nParams*2]);
        }
    }

    TString WriteParameterToFile(TF1 *Funktion){

        TString name = Funktion->GetName();
        Int_t NumberParam = Funktion->GetNpar();
        TString  paramName[20];
        Double_t param[20];
        Double_t paramErr[20];
        Double_t chiSquare = Funktion->GetChisquare();
        Double_t NDF = Funktion->GetNDF();

        TString  paramLine[20];

        for( Int_t i = 0; i < NumberParam; i++){
            paramName[i] = Funktion->GetParName(i);
            param[i] = Funktion->GetParameter(i);
            paramErr[i] = Funktion->GetParError(i);
            paramLine[i] = Form("%s \t %.10f \t %.10f",paramName[i].Data(),param[i],paramErr[i]);
        }

        TString stringOutput = "\n\n"+name+"\n";
        for( Int_t i = 0; i < NumberParam; i++){
            stringOutput = stringOutput+paramLine[i]+"\n";
        }
        stringOutput = stringOutput + Form("Chi2 \t %.10f \t ndf \t %.10f \n",chiSquare,NDF);
        if (NDF!=0) stringOutput = stringOutput + Form("Chi2/ndf \t %.10f \n",chiSquare/NDF);


        return stringOutput;
    }

    TString WriteParameterToFileLatexTable(TF1 *Funktion, Bool_t isRad = 0){

        TString name = Funktion->GetName();
        Int_t NumberParam = Funktion->GetNpar();
        TString  paramName[20];
        Double_t param[20];
        Double_t paramErr[20];
        Double_t chiSquare = Funktion->GetChisquare();
        Double_t NDF = Funktion->GetNDF();

        TString  paramLine;
        TString  paramLine0;
        for( Int_t i = 0; i < NumberParam; i++){
            paramName[i] = Funktion->GetParName(i);
            param[i] = Funktion->GetParameter(i);
            paramErr[i] = Funktion->GetParError(i);
            if (isRad && i == 2) {
                param[i] = param[i]*1000;
                paramErr[i] = paramErr[i]*1000;
            }
            paramLine0 = paramLine0 + Form("%s \t ",paramName[i].Data());
            paramLine = paramLine + Form("$%.2f +- %.2f$ &\t",param[i],paramErr[i]);
        }

        TString stringOutput = "\n\n"+name+"\n" + paramLine0 + "\n" + paramLine ;
        stringOutput = stringOutput + Form("\n Chi2 \t %.10f \t ndf \t %.10f \n",chiSquare,NDF);
        if (NDF!=0) stringOutput = stringOutput + Form("Chi2/ndf \t %.10f \n",chiSquare/NDF);
        return stringOutput;
    }


    TH1D *FitToHist(TF1* FitFunction, TH1 *FittedHist, TString FitHistName){

        TH1D *FitHist = new TH1D();
        FitHist->Sumw2();
        FitHist = (TH1D*) FittedHist->Clone(FitHistName);
        Double_t binning = FittedHist->GetNbinsX();
        for(Int_t bin = 1; bin<= binning; bin++){
            //      FitHist->SetBinError(bin,FittedHist->GetBinError(bin)); //--->To be done
            FitHist->SetBinError(bin,0.0000000001);
            FitHist->SetBinContent(bin,FitFunction->Eval(FittedHist->GetBinCenter(bin)));
        }

        return FitHist;
    }

    TF1* FitTH1DRecursivelyGaussian (TH1D* histo, Double_t precision, Double_t fitRangeMin, Double_t fitRangeMax ){
        TF1 *f0 = new TF1("f0", "gaus", fitRangeMin,fitRangeMax);
        histo->Fit(f0,"0RMEQ");
        Double_t rp = f0->GetParameter(2);
        Double_t mp = f0->GetParameter(1);
        Double_t ymin = mp -(rp * 2.);
        Double_t ymax = mp + (rp * 2.);
        Double_t deviation = 100;
        Int_t counter = 0;
        TF1* f1 = new TF1 ("f1", "gaus", ymin, ymax);
        while(deviation > precision && counter < 100){
            f1->SetRange(ymin,ymax);
            histo->Fit(f1,"0RMEQ");
            Double_t rp2 = f1->GetParameter(2);
            if (rp2>rp){ deviation = rp2-rp;} 
                else {deviation = rp -rp2 ;}
            rp = rp2 ;
            mp = f1->GetParameter(1);
            ymin = mp -(rp * 2);
            ymax = mp +(rp * 2);
            counter++;
        }
        delete f0;
        return f1;
    }


    void ResolutionFitting( TH2* histoRebinned, 
                            TH1D** histoArray, 
                            TF1** fitArray,
                            Int_t nBins, 
                            TH1F* histoMean, 
                            TH1F* histoSigma, 
                            TString fitName, 
                            Double_t rangeMin,
                            Double_t rangeMax, 
                            Double_t precision, 
                            TString defaultName){
        
        for (Int_t i = 1; i < nBins + 1 ; i ++){
            fitArray[i-1] = 0x00;
            histoArray[i-1] = histoRebinned->ProjectionY(Form("%s_%i",defaultName.Data(),i-1), i, i);
            if (histoArray[i-1]->GetEntries() > 0){
                TF1 *f0;
                if (fitName.CompareTo("gaus(0)+gaus(3)")==0){
                    f0 = new TF1("f0", "gaus" , rangeMin ,rangeMax);
                } else {
                    f0 = new TF1("f0", fitName.Data() , rangeMin ,rangeMax);
                }
                if (fitName.CompareTo("gaus")==0){
                    f0->SetParameter(0,1.5*histoArray[i-1]->GetMaximum());
                    f0->SetParameter(0,0.);
                }
                histoArray[i-1]->Fit(f0,"0RMEQ");
                Double_t rp = f0->GetParameter(2);
                Double_t mp = f0->GetParameter(1);
                Double_t ymin = mp -(rp * 2.);
                Double_t ymax = mp + (rp * 2.);
                Double_t deviation = 100;
                Int_t counter = 0;
                fitArray[i-1]  = new TF1 ("f1", fitName.Data(), ymin, ymax);
                if (fitName.CompareTo("gaus")==0){
                    fitArray[i-1]->SetParameter(0,1.5*histoArray[i-1]->GetMaximum());
                    fitArray[i-1]->SetParameter(0,0.);
                } else if (fitName.CompareTo("gaus(0)+gaus(3)")==0){
                    fitArray[i-1]->SetRange(rangeMin ,rangeMax);
                    fitArray[i-1]->SetParameter(0,histoArray[i-1]->GetMaximum());
                    fitArray[i-1]->SetParameter(1,0);
                    fitArray[i-1]->SetParameter(2,0.5);
                    fitArray[i-1]->SetParameter(3,0.5*histoArray[i-1]->GetMaximum());
                    fitArray[i-1]->SetParameter(4,f0->GetParameter(1));
                    fitArray[i-1]->SetParameter(5,f0->GetParameter(2));
                }

                if (fitName.CompareTo("gaus(0)+gaus(3)")!=0){
                    while(deviation > precision && counter < 100){
                        fitArray[i-1]->SetRange(ymin,ymax);
                        
                        histoArray[i-1]->Fit(fitArray[i-1],"0RMEQ");
                        Double_t rp2 = fitArray[i-1]->GetParameter(2);
                        if (rp2>rp){ deviation = rp2-rp;} 
                            else {deviation = rp -rp2 ;}
                        rp = rp2 ;
                        mp = fitArray[i-1]->GetParameter(1);
                        ymin = mp -(rp * 2);
                        ymax = mp +(rp * 2);
                        counter++;
                    }
                } else {
                    histoArray[i-1]->Fit(fitArray[i-1],"0RMEQ");
                }
                Double_t mean;
                Double_t meanE;
                Double_t sigma;
                Double_t sigmaE;
                
                mean = fitArray[i-1]->GetParameter(1);
                meanE = fitArray[i-1]->GetParError(1);
                sigma = fitArray[i-1]->GetParameter(2);
                sigmaE = fitArray[i-1]->GetParError(2);
                
                histoMean->SetBinContent(i, mean);
                histoMean->SetBinError(i, meanE);
                histoSigma->SetBinContent(i, sigma);
                histoSigma->SetBinError(i, sigmaE);
                delete f0;
            } else {
                histoMean->SetBinContent(i, 0);
                histoMean->SetBinError(i, 0);
                histoSigma->SetBinContent(i, 0);
                histoSigma->SetBinError(i, 0);
            }
        } //end of fitting 
    }

    void ResolutionFittingNormalized( TH2* histoRebinned, TH1D** histoArray, TF1** fitArray , TH1D* histoNorm, Int_t nBins, TH1F* histoMean, TH1F* histoSigma, TString fitName, Double_t rangeMin, Double_t rangeMax, Double_t precision, TString defaultName){
        for (Int_t i = 1; i < nBins + 1 ; i ++){
            fitArray[i-1] = 0x00;
            histoArray[i-1] = histoRebinned->ProjectionY(Form("%s_%i",defaultName.Data(),i-1), i, i);
            histoArray[i-1]->Sumw2();
            histoArray[i-1]->Divide(histoNorm);
            if (histoArray[i-1]->GetEntries() > 0){
                TF1 *f0;
                if (fitName.CompareTo("gaus(0)+gaus(3)")==0){
                    f0 = new TF1("f0", "gaus" , rangeMin ,rangeMax);
                } else if (fitName.CompareTo("gaus(0)+gaus(3)")==0){
                    f0 = new TF1("GaussExpLinear","(x>=[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(exp(-0.5*((x-[1])/[2])^2)-1)))+(x<[1])*([0]*exp(-0.5*((x-[1])/[2])^2))",rangeMin,rangeMax);
                } else {
                    f0 = new TF1("f0", fitName.Data() , rangeMin ,rangeMax);
                }
                        
                if (fitName.CompareTo("gaus")==0){
                    f0->SetParameter(0,1.5*histoArray[i-1]->GetMaximum());
                    f0->SetParameter(0,0.);
                }
                histoArray[i-1]->Fit(f0,"0RMEQF");
                Double_t rp = f0->GetParameter(2);
                Double_t mp = f0->GetParameter(1);
                Double_t ymin = mp -(rp * 2.);
                Double_t ymax = mp + (rp * 2.);
                Double_t deviation = 100;
                Int_t counter = 0;
                fitArray[i-1]  = new TF1 ("f1", fitName.Data(), ymin, ymax);
                if (fitName.CompareTo("gaus")==0){
                    fitArray[i-1]->SetParameter(0,1.5*histoArray[i-1]->GetMaximum());
                    fitArray[i-1]->SetParameter(0,0.);
                } else if (fitName.CompareTo("gaus(0)+gaus(3)")==0){
                    fitArray[i-1]->SetRange(rangeMin ,rangeMax);
                    fitArray[i-1]->SetParameter(0,histoArray[i-1]->GetMaximum());
                    fitArray[i-1]->SetParameter(1,0);
                    fitArray[i-1]->SetParameter(2,0.5);
                    fitArray[i-1]->SetParameter(3,0.5*histoArray[i-1]->GetMaximum());
                    fitArray[i-1]->SetParameter(4,f0->GetParameter(1));
                    fitArray[i-1]->SetParameter(5,f0->GetParameter(2));
                }

                if (fitName.CompareTo("gaus(0)+gaus(3)")!=0){
                    while(deviation > precision && counter < 100){
                        fitArray[i-1]->SetRange(ymin,ymax);
                        
                        histoArray[i-1]->Fit(fitArray[i-1],"0RMEQ");
                        Double_t rp2 = fitArray[i-1]->GetParameter(2);
                        if (rp2>rp){ deviation = rp2-rp;} 
                            else {deviation = rp -rp2 ;}
                        rp = rp2 ;
                        mp = fitArray[i-1]->GetParameter(1);
                        ymin = mp -(rp * 2);
                        ymax = mp +(rp * 2);
                        counter++;
                    }
                } else {
                    histoArray[i-1]->Fit(fitArray[i-1],"0RMEQ");
                }
                Double_t mean;
                Double_t meanE;
                Double_t sigma;
                Double_t sigmaE;
                
                mean = fitArray[i-1]->GetParameter(1);
                meanE = fitArray[i-1]->GetParError(1);
                sigma = fitArray[i-1]->GetParameter(2);
                sigmaE = fitArray[i-1]->GetParError(2);
                
                histoMean->SetBinContent(i, mean);
                histoMean->SetBinError(i, meanE);
                histoSigma->SetBinContent(i, sigma);
                histoSigma->SetBinError(i, sigmaE);
                delete f0;
            } else {
                histoMean->SetBinContent(i, 0);
                histoMean->SetBinError(i, 0);
                histoSigma->SetBinContent(i, 0);
                histoSigma->SetBinError(i, 0);
            }
        } //end of fitting 
    }


    void ResolutionFittingRebined( TH2* histo, 
                                TH1D** histoArray, 
                                TF1** fitArray ,  
                                Int_t nBins,
                                Double_t* ptRanges,
                                Int_t* rebin, 
                                TH1F* histoMean, 
                                TH1F* histoSigma, 
                                TString fitName, 
                                Double_t rangeMin, 
                                Double_t rangeMax, 
                                Double_t precision, 
                                TString defaultName){
        for (Int_t i = 1; i < nBins + 1 ; i ++){
            cout << "Bin \t" << i-1 << "\t pt range: "<< ptRanges[i-1]<< " - " << ptRanges[i] <<  endl;
            fitArray[i-1] = 0x00;
            histoArray[i-1] = histo->ProjectionY(Form("%s_%i",defaultName.Data(),i-1), histo->GetXaxis()->FindBin(ptRanges[i-1]+0.001), histo->GetXaxis()->FindBin(ptRanges[i]-0.001));
            histoArray[i-1]->Sumw2();
            histoArray[i-1]->Rebin(rebin[i-1]);
            if (histoArray[i-1]->GetEntries() > 0){
                TF1 *f0;
                if (fitName.CompareTo("gaus(0)+gaus(3)")==0){
                    f0 = new TF1("f0", "gaus" , rangeMin ,rangeMax);
                } else if (fitName.CompareTo("gaus(0)+gaus(3)")==0){
                    f0 = new TF1("GaussExpLinear","(x>=[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(exp(-0.5*((x-[1])/[2])^2)-1)))+(x<[1])*([0]*exp(-0.5*((x-[1])/[2])^2))",rangeMin,rangeMax);
                } else {
                    f0 = new TF1("f0", fitName.Data() , rangeMin ,rangeMax);
                }
                        
                if (fitName.CompareTo("gaus")==0){
                    f0->SetParameter(0,1.5*histoArray[i-1]->GetMaximum());
                    f0->SetParameter(1,0.);
                    f0->SetParameter(2,0.02);
                }
                histoArray[i-1]->Fit(f0,"0RMEQF");
                Double_t rp = f0->GetParameter(2);
                Double_t mp = f0->GetParameter(1);
                Double_t ymin = mp -(rp * 2.);
                Double_t ymax = mp + (rp * 2.);
                Double_t deviation = 100;
                Int_t counter = 0;
                fitArray[i-1]  = new TF1 ("f1", fitName.Data(), ymin, ymax);
                if (fitName.CompareTo("gaus")==0){
                    fitArray[i-1]->SetParameter(0,1.5*histoArray[i-1]->GetMaximum());
                    fitArray[i-1]->SetParameter(0,0.);
                } else if (fitName.CompareTo("gaus(0)+gaus(3)")==0){
                    fitArray[i-1]->SetRange(rangeMin ,rangeMax);
                    fitArray[i-1]->SetParameter(0,histoArray[i-1]->GetMaximum());
                    fitArray[i-1]->SetParameter(1,0);
                    fitArray[i-1]->SetParameter(2,0.5);
                    fitArray[i-1]->SetParameter(3,0.5*histoArray[i-1]->GetMaximum());
                    fitArray[i-1]->SetParameter(4,f0->GetParameter(1));
                    fitArray[i-1]->SetParameter(5,f0->GetParameter(2));
                }

                if (fitName.CompareTo("gaus(0)+gaus(3)")!=0){
                    while(deviation > precision && counter < 100){
                        fitArray[i-1]->SetRange(ymin,ymax);
                        
                        histoArray[i-1]->Fit(fitArray[i-1],"0RMEQ");
                        Double_t rp2 = fitArray[i-1]->GetParameter(2);
                        if (rp2>rp){ deviation = rp2-rp;} 
                            else {deviation = rp -rp2 ;}
                        rp = rp2 ;
                        mp = fitArray[i-1]->GetParameter(1);
                        ymin = mp -(rp * 2);
                        ymax = mp +(rp * 2);
                        counter++;
                    }
                } else {
                    histoArray[i-1]->Fit(fitArray[i-1],"0RMEQ");
                }
                Double_t mean;
                Double_t meanE;
                Double_t sigma;
                Double_t sigmaE;
                
                mean = fitArray[i-1]->GetParameter(1);
                meanE = fitArray[i-1]->GetParError(1);
                sigma = fitArray[i-1]->GetParameter(2);
                sigmaE = fitArray[i-1]->GetParError(2);
                cout << "Mean: " << mean << " +- " << meanE << ", sigma: " << sigma << " +- " << sigmaE << endl;
                histoMean->SetBinContent(i, mean);
                histoMean->SetBinError(i, meanE);
                histoSigma->SetBinContent(i, sigma);
                histoSigma->SetBinError(i, sigmaE);
                delete f0;
            } else {
                histoMean->SetBinContent(i, 0);
                histoMean->SetBinError(i, 0);
                histoSigma->SetBinContent(i, 0);
                histoSigma->SetBinError(i, 0);
            }
        } //end of fitting 
    }

    void ResolutionFittingArray(TH1D** histoArray, TF1** fitArray , Int_t nBins, TH1F* histoMean, TH1F* histoSigma, TString fitName, Double_t rangeMin, Double_t rangeMax, Double_t precision, Double_t meanInput =0., Double_t sigmaInput=0.5, Double_t offsetMean = 0, Double_t fractionSecondGauss=0.5){
        for (Int_t i = 1; i < nBins + 1 ; i ++){
            fitArray[i-1] = 0x00;
            if (histoArray[i-1]->GetEntries() > 0){
                TF1 *f0;
                if (fitName.CompareTo("gaus(0)+gaus(3)")==0){
                    f0 = new TF1("f0", "gaus" , rangeMin ,rangeMax);
                } else {
                    f0 = new TF1("f0", fitName.Data() , rangeMin ,rangeMax);
                }
                if (fitName.CompareTo("gaus")==0){
                    f0->SetParameter(0,1.5*histoArray[i-1]->GetMaximum());
                    f0->SetParameter(0,meanInput);
                }
                histoArray[i-1]->Fit(f0,"0RMEQ");
                Double_t rp = f0->GetParameter(2);
                Double_t mp = f0->GetParameter(1);
                Double_t ymin = mp -(rp * 2.);
                Double_t ymax = mp + (rp * 2.);
                Double_t deviation = 100;
                Int_t counter = 0;
                fitArray[i-1]  = new TF1 ("f1", fitName.Data(), ymin, ymax);
                if (fitName.CompareTo("gaus")==0){
                    fitArray[i-1]->SetParameter(0,1.5*histoArray[i-1]->GetMaximum());
                    fitArray[i-1]->SetParameter(0,0.);
                } else if (fitName.CompareTo("gaus(0)+gaus(3)")==0){
                    fitArray[i-1]->SetRange(rangeMin ,rangeMax);
                    fitArray[i-1]->SetParameter(0,histoArray[i-1]->GetMaximum());
                    fitArray[i-1]->SetParameter(1,meanInput);
                    fitArray[i-1]->SetParameter(2,sigmaInput);
                    fitArray[i-1]->SetParameter(3,fractionSecondGauss*histoArray[i-1]->GetMaximum());
                    fitArray[i-1]->SetParameter(4,histoArray[i-1]->GetMean()+offsetMean);
                    fitArray[i-1]->SetParameter(5,f0->GetParameter(2));
                }

                if (fitName.CompareTo("gaus(0)+gaus(3)")!=0){
                    while(deviation > precision && counter < 100){
                        fitArray[i-1]->SetRange(ymin,ymax);
                        
                        histoArray[i-1]->Fit(fitArray[i-1],"0RMEQ");
                        Double_t rp2 = fitArray[i-1]->GetParameter(2);
                        if (rp2>rp){ deviation = rp2-rp;} 
                            else {deviation = rp -rp2 ;}
                        rp = rp2 ;
                        mp = fitArray[i-1]->GetParameter(1);
                        ymin = mp -(rp * 2);
                        ymax = mp +(rp * 2);
                        counter++;
                    }
                } else {
                    histoArray[i-1]->Fit(fitArray[i-1],"0RMEQ");
                }
                Double_t mean;
                Double_t meanE;
                Double_t sigma;
                Double_t sigmaE;
                
                mean = fitArray[i-1]->GetParameter(1);
                meanE = fitArray[i-1]->GetParError(1);
                sigma = fitArray[i-1]->GetParameter(2);
                sigmaE = fitArray[i-1]->GetParError(2);
                
                fitArray[i-1]->SetRange(rangeMin, rangeMax);
                
                histoMean->SetBinContent(i, mean);
                histoMean->SetBinError(i, meanE);
                histoSigma->SetBinContent(i, sigma);
                histoSigma->SetBinError(i, sigmaE);
                delete f0;
            } else {
                histoMean->SetBinContent(i, 0);
                histoMean->SetBinError(i, 0);
                histoSigma->SetBinContent(i, 0);
                histoSigma->SetBinError(i, 0);
            }
        } //end of fitting 
    }

    void ResolutionFittingArrayTH1D(TH1D** histoArray, TF1** fitArray , Int_t nBins, TH1D* histoMean, TH1D* histoSigma, TString fitName, Double_t rangeMin, Double_t rangeMax, Double_t precision, Double_t meanInput =0., Double_t sigmaInput=0.5, Double_t offsetMean = 0, Double_t fractionSecondGauss=0.5){
        for (Int_t i = 1; i < nBins + 1 ; i ++){
            fitArray[i-1] = 0x00;
            if (histoArray[i-1]->GetEntries() > 0){
                TF1 *f0;
                if (fitName.CompareTo("gaus(0)+gaus(3)")==0){
                    f0 = new TF1("f0", "gaus" , rangeMin ,rangeMax);
                } else {
                    f0 = new TF1("f0", fitName.Data() , rangeMin ,rangeMax);
                }
                if (fitName.CompareTo("gaus")==0){
                    f0->SetParameter(0,1.5*histoArray[i-1]->GetMaximum());
                    f0->SetParameter(0,meanInput);
                }
                histoArray[i-1]->Fit(f0,"0RMEQ");
                Double_t rp = f0->GetParameter(2);
                Double_t mp = f0->GetParameter(1);
                Double_t ymin = mp -(rp * 2.);
                Double_t ymax = mp + (rp * 2.);
                Double_t deviation = 100;
                Int_t counter = 0;
                fitArray[i-1]  = new TF1 ("f1", fitName.Data(), ymin, ymax);
                if (fitName.CompareTo("gaus")==0){
                    fitArray[i-1]->SetParameter(0,1.5*histoArray[i-1]->GetMaximum());
                    fitArray[i-1]->SetParameter(0,0.);
                } else if (fitName.CompareTo("gaus(0)+gaus(3)")==0){
                    fitArray[i-1]->SetRange(rangeMin ,rangeMax);
                    fitArray[i-1]->SetParameter(0,histoArray[i-1]->GetMaximum());
                    fitArray[i-1]->SetParameter(1,meanInput);
                    fitArray[i-1]->SetParameter(2,sigmaInput);
                    fitArray[i-1]->SetParameter(3,fractionSecondGauss*histoArray[i-1]->GetMaximum());
                    fitArray[i-1]->SetParameter(4,histoArray[i-1]->GetMean()+offsetMean);
                    fitArray[i-1]->SetParameter(5,f0->GetParameter(2));
                }

                if (fitName.CompareTo("gaus(0)+gaus(3)")!=0){
                    while(deviation > precision && counter < 100){
                        fitArray[i-1]->SetRange(ymin,ymax);
                        
                        histoArray[i-1]->Fit(fitArray[i-1],"0RMEQ");
                        Double_t rp2 = fitArray[i-1]->GetParameter(2);
                        if (rp2>rp){ deviation = rp2-rp;} 
                            else {deviation = rp -rp2 ;}
                        rp = rp2 ;
                        mp = fitArray[i-1]->GetParameter(1);
                        ymin = mp -(rp * 2);
                        ymax = mp +(rp * 2);
                        counter++;
                    }
                } else {
                    histoArray[i-1]->Fit(fitArray[i-1],"0RMEQ");
                }
                Double_t mean;
                Double_t meanE;
                Double_t sigma;
                Double_t sigmaE;
                
                mean = fitArray[i-1]->GetParameter(1);
                meanE = fitArray[i-1]->GetParError(1);
                sigma = fitArray[i-1]->GetParameter(2);
                sigmaE = fitArray[i-1]->GetParError(2);
                
                fitArray[i-1]->SetRange(rangeMin, rangeMax);
                
                histoMean->SetBinContent(i, mean);
                histoMean->SetBinError(i, meanE);
                histoSigma->SetBinContent(i, sigma);
                histoSigma->SetBinError(i, sigmaE);
                delete f0;
            } else {
                histoMean->SetBinContent(i, 0);
                histoMean->SetBinError(i, 0);
                histoSigma->SetBinContent(i, 0);
                histoSigma->SetBinError(i, 0);
            }
        } //end of fitting 
    }

    void ResolutionFittingTH1D( TH2* histoRebinned, TH1D** histoArray, TF1** fitArray , Int_t nBins, TH1D* histoMean, TH1D* histoSigma, TString fitName, Double_t rangeMin, Double_t rangeMax, Double_t precision, TString defaultName){
        for (Int_t i = 1; i < nBins + 1 ; i ++){
            fitArray[i-1] = 0x00;
            histoArray[i-1] = histoRebinned->ProjectionY(Form("%s_%i",defaultName.Data(),i-1), i, i);
            if (histoArray[i-1]->GetEntries() > 0){
                TF1 *f0;
                if (fitName.CompareTo("gaus(0)+gaus(3)")==0){
                    f0 = new TF1("f0", "gaus" , rangeMin ,rangeMax);
                } else {
                    f0 = new TF1("f0", fitName.Data() , rangeMin ,rangeMax);
                }
                if (fitName.CompareTo("gaus")==0){
                    f0->SetParameter(0,1.5*histoArray[i-1]->GetMaximum());
                    f0->SetParameter(0,0.);
                }
                histoArray[i-1]->Fit(f0,"0RMEQ");
                Double_t rp = f0->GetParameter(2);
                Double_t mp = f0->GetParameter(1);
                Double_t ymin = mp -(rp * 2.);
                Double_t ymax = mp + (rp * 2.);
                Double_t deviation = 100;
                Int_t counter = 0;
                fitArray[i-1]  = new TF1 ("f1", fitName.Data(), ymin, ymax);
                if (fitName.CompareTo("gaus")==0){
                    fitArray[i-1]->SetParameter(0,1.5*histoArray[i-1]->GetMaximum());
                    fitArray[i-1]->SetParameter(0,0.);
                } else if (fitName.CompareTo("gaus(0)+gaus(3)")==0){
                    fitArray[i-1]->SetRange(rangeMin ,rangeMax);
                    fitArray[i-1]->SetParameter(0,histoArray[i-1]->GetMaximum());
                    fitArray[i-1]->SetParameter(1,0);
                    fitArray[i-1]->SetParameter(2,0.5);
                    fitArray[i-1]->SetParameter(3,0.5*histoArray[i-1]->GetMaximum());
                    fitArray[i-1]->SetParameter(4,f0->GetParameter(1));
                    fitArray[i-1]->SetParameter(5,f0->GetParameter(2));
                }

                if (fitName.CompareTo("gaus(0)+gaus(3)")!=0){
                    while(deviation > precision && counter < 100){
                        fitArray[i-1]->SetRange(ymin,ymax);
                        
                        histoArray[i-1]->Fit(fitArray[i-1],"0RMEQ");
                        Double_t rp2 = fitArray[i-1]->GetParameter(2);
                        if (rp2>rp){ deviation = rp2-rp;} 
                            else {deviation = rp -rp2 ;}
                        rp = rp2 ;
                        mp = fitArray[i-1]->GetParameter(1);
                        ymin = mp -(rp * 2);
                        ymax = mp +(rp * 2);
                        counter++;
                    }
                } else {
                    histoArray[i-1]->Fit(fitArray[i-1],"0RMEQ");
                }
                Double_t mean;
                Double_t meanE;
                Double_t sigma;
                Double_t sigmaE;
                
                mean = fitArray[i-1]->GetParameter(1);
                meanE = fitArray[i-1]->GetParError(1);
                sigma = fitArray[i-1]->GetParameter(2);
                sigmaE = fitArray[i-1]->GetParError(2);
                
                histoMean->SetBinContent(i, mean);
                histoMean->SetBinError(i, meanE);
                histoSigma->SetBinContent(i, sigma);
                histoSigma->SetBinError(i, sigmaE);
                delete f0;
            } else {
                histoMean->SetBinContent(i, 0);
                histoMean->SetBinError(i, 0);
                histoSigma->SetBinContent(i, 0);
                histoSigma->SetBinError(i, 0);
            }
        } //end of fitting 
    }

    //*****************************************************************/
    //* BOLTZMANN
    //*****************************************************************/

    Double_t
    Boltzmann_Func(const Double_t *x, const Double_t *p)
    {
    /* dN/dpt */
    
    Double_t pt = x[0];
    Double_t mass = p[0];
    Double_t mt = TMath::Sqrt(pt * pt + mass * mass);
    Double_t T = p[1];
    Double_t norm = p[2];

    return pt * norm * mt * TMath::Exp(-mt / T);
    }

    TF1 *
    Boltzmann(const Char_t *name, Double_t mass, Double_t T = 0.1, Double_t norm = 1.)
    {
    
    TF1 *fBoltzmann = new TF1(name, Boltzmann_Func, 0., 10., 3);
    fBoltzmann->SetParameters(mass, T, norm);
    fBoltzmann->SetParNames("mass", "T", "norm");
    fBoltzmann->FixParameter(0, mass);
    return fBoltzmann;
    }

    /*****************************************************************/
    /* LEVY-TSALLIS */
    /*****************************************************************/

    Double_t
    LevyTsallis_Func(const Double_t *x, const Double_t *p)
    {
    /* dN/dpt */
    
    Double_t pt = x[0];
    Double_t mass = p[0];
    Double_t mt = TMath::Sqrt(pt * pt + mass * mass);
    Double_t n = p[1];
    Double_t C = p[2];
    Double_t norm = p[3];

    Double_t part1 = (n - 1.) * (n - 2);
    Double_t part2 = n * C * (n * C + mass * (n - 2.));
    Double_t part3 = part1 / part2;
    Double_t part4 = 1. + (mt - mass) / n / C;
    Double_t part5 = TMath::Power(part4, -n);
    return pt * norm * part3 * part5;
    }

    TF1 *
    LevyTsallis(const Char_t *name, Double_t mass, Double_t n = 5., Double_t C = 0.1, Double_t norm = 1.)
    {
    
    TF1 *fLevyTsallis = new TF1(name, LevyTsallis_Func, 0., 10., 4);
    fLevyTsallis->SetParameters(mass, n, C, norm);
    fLevyTsallis->SetParNames("mass", "n", "C", "norm");
    fLevyTsallis->FixParameter(0, mass);
    return fLevyTsallis;
    }

    /*****************************************************************/
    /* BOLTZMANN-GIBBS BLAST-WAVE */
    /*****************************************************************/

    static TF1 *fBGBlastWave_Integrand = NULL;
    Double_t
    BGBlastWave_Integrand(const Double_t *x, const Double_t *p)
    {
    
    /* 
        x[0] -> r (radius)
        p[0] -> mT (transverse mass)
        p[1] -> pT (transverse momentum)
        p[2] -> beta_max (surface velocity)
        p[3] -> T (freezout temperature)
        p[4] -> n (velocity profile)
    */
    
    Double_t r = x[0];
    Double_t mt = p[0];
    Double_t pt = p[1];
    Double_t beta_max = p[2];
    Double_t temp_1 = 1. / p[3];
    Double_t n = p[4];

    Double_t beta = beta_max * TMath::Power(r, n);
    if (beta > 0.9999999999999999) beta = 0.9999999999999999;
    Double_t rho = TMath::ATanH(beta);
    Double_t argI0 = pt * TMath::SinH(rho) * temp_1;
    if (argI0 > 700.) argI0 = 700.;
    Double_t argK1 = mt * TMath::CosH(rho) * temp_1;
    //  if (argI0 > 100 || argI0 < -100)
    //    printf("r=%f, pt=%f, beta_max=%f, temp=%f, n=%f, mt=%f, beta=%f, rho=%f, argI0=%f, argK1=%f\n", r, pt, beta_max, 1. / temp_1, n, mt, beta, rho, argI0, argK1);
    return r * mt * TMath::BesselI0(argI0) * TMath::BesselK1(argK1);
    
    }

    Double_t
    BGBlastWave_Func(const Double_t *x, const Double_t *p)
    {
    /* dN/dpt */
    
    Double_t pt = x[0];
    Double_t mass = p[0];
    Double_t mt = TMath::Sqrt(pt * pt + mass * mass);
    Double_t beta_max = p[1];
    Double_t temp = p[2];
    Double_t n = p[3];
    Double_t norm = p[4];
    
    if (!fBGBlastWave_Integrand)
        fBGBlastWave_Integrand = new TF1("fBGBlastWave_Integrand", BGBlastWave_Integrand, 0., 1., 5);
    fBGBlastWave_Integrand->SetParameters(mt, pt, beta_max, temp, n);
    Double_t integral = fBGBlastWave_Integrand->Integral(0., 1.);
    return norm * pt * integral;
    }

    Double_t
    BGBlastWave_Func_OneOverPt(const Double_t *x, const Double_t *p)
    {
    /* 1/pt dN/dpt */
    
    Double_t pt = x[0];
    Double_t mass = p[0];
    Double_t mt = TMath::Sqrt(pt * pt + mass * mass);
    Double_t beta_max = p[1];
    Double_t temp = p[2];
    Double_t n = p[3];
    Double_t norm = p[4];
    
    if (!fBGBlastWave_Integrand)
        fBGBlastWave_Integrand = new TF1("fBGBlastWave_Integrand", BGBlastWave_Integrand, 0., 1., 5);
    fBGBlastWave_Integrand->SetParameters(mt, pt, beta_max, temp, n);
    Double_t integral = fBGBlastWave_Integrand->Integral(0., 1.);

    return norm * integral;
    }


    TF1 *
    BGBlastWave(const Char_t *name, Double_t mass, Double_t xmin = 0, Double_t xmax = 10,Double_t beta_max = 0.9, Double_t temp = 0.1, Double_t n = 1., Double_t norm = 1.e6)
    {
    
    TF1 *fBGBlastWave = new TF1(name, BGBlastWave_Func, xmin, xmax, 5);
    fBGBlastWave->SetParameters(mass, beta_max, temp, n, norm);
    fBGBlastWave->SetParNames("mass", "beta_max", "T", "n", "norm");
    fBGBlastWave->FixParameter(0, mass);
    fBGBlastWave->SetParLimits(1, 0.01, 0.99);
    fBGBlastWave->SetParLimits(2, 0.01, 1.);
    fBGBlastWave->SetParLimits(3, 0.01, 10.);
    return fBGBlastWave;
    }

    TF1 * BGBlastWave_OneOverPT(const Char_t *name, Double_t mass, Double_t xmin = 0, Double_t xmax = 1, Double_t beta_max = 0.9, Double_t temp = 0.1, Double_t n = 1., Double_t norm = 1.e6)
    {
    
    TF1 *fBGBlastWave = new TF1(name, BGBlastWave_Func_OneOverPt, xmin, xmax, 5);
    fBGBlastWave->SetParameters(mass, beta_max, temp, n, norm);
    fBGBlastWave->SetParNames("mass", "beta_max", "T", "n", "norm");
    fBGBlastWave->FixParameter(0, mass);
    fBGBlastWave->SetParLimits(1, 0.01, 0.99);
    fBGBlastWave->SetParLimits(2, 0.01, 1.);
    fBGBlastWave->SetParLimits(3, 0.01, 10.);
    return fBGBlastWave;
    }

    /*****************************************************************/
    /* TSALLIS BLAST-WAVE */
    /*****************************************************************/

    static TF1 *fTsallisBlastWave_Integrand_r = NULL;
    Double_t
    TsallisBlastWave_Integrand_r(const Double_t *x, const Double_t *p)
    {
    /* 
        x[0] -> r (radius)
        p[0] -> mT (transverse mass)
        p[1] -> pT (transverse momentum)
        p[2] -> beta_max (surface velocity)
        p[3] -> T (freezout temperature)
        p[4] -> n (velocity profile)
        p[5] -> q
        p[6] -> y (rapidity)
        p[7] -> phi (azimuthal angle)
    */
    
    Double_t r = x[0];
    Double_t mt = p[0];
    Double_t pt = p[1];
    Double_t beta_max = p[2];
    Double_t temp_1 = 1. / p[3];
    Double_t n = p[4];
    Double_t q = p[5];
    Double_t y = p[6];
    Double_t phi = p[7];

    if (q <= 1.) return r;

    Double_t beta = beta_max * TMath::Power(r, n);
    Double_t rho = TMath::ATanH(beta);
    
    Double_t part1 = mt * TMath::CosH(y) * TMath::CosH(rho);
    Double_t part2 = pt * TMath::SinH(rho) * TMath::Cos(phi);
    Double_t part3 = part1 - part2;
    Double_t part4 = 1 + (q - 1.) * temp_1 * part3;
    Double_t expo = -1. / (q - 1.);
    //  printf("part1=%f, part2=%f, part3=%f, part4=%f, expo=%f\n", part1, part2, part3, part4, expo);
    Double_t part5 = TMath::Power(part4, expo);

    return r * part5;
    }

    static TF1 *fTsallisBlastWave_Integrand_phi = NULL;
    Double_t
    TsallisBlastWave_Integrand_phi(const Double_t *x, const Double_t *p)
    {
    /* 
        x[0] -> phi (azimuthal angle)
    */
    if(p){};
    Double_t phi = x[0];
    fTsallisBlastWave_Integrand_r->SetParameter(7, phi);
    Double_t integral = fTsallisBlastWave_Integrand_r->Integral(0., 1.);
    return integral;
    }

    static TF1 *fTsallisBlastWave_Integrand_y = NULL;
    Double_t
    TsallisBlastWave_Integrand_y(const Double_t *x, const Double_t *p)
    {
    /* 
        x[0] -> y (rapidity)
    */
    if(p){};
    Double_t y = x[0];
    fTsallisBlastWave_Integrand_r->SetParameter(6, y);
    Double_t integral = fTsallisBlastWave_Integrand_phi->Integral(-TMath::Pi(), TMath::Pi());
    return TMath::CosH(y) * integral;
    }

    Double_t
    TsallisBlastWave_Func(const Double_t *x, const Double_t *p)
    {
    /* dN/dpt */
    
    Double_t pt = x[0];
    Double_t mass = p[0];
    Double_t mt = TMath::Sqrt(pt * pt + mass * mass);
    Double_t beta_max = p[1];
    Double_t temp = p[2];
    Double_t n = p[3];
    Double_t q = p[4];
    Double_t norm = p[5];

    if (!fTsallisBlastWave_Integrand_r)
        fTsallisBlastWave_Integrand_r = new TF1("fTsallisBlastWave_Integrand_r", TsallisBlastWave_Integrand_r, 0., 1., 8);
    if (!fTsallisBlastWave_Integrand_phi)
        fTsallisBlastWave_Integrand_phi = new TF1("fTsallisBlastWave_Integrand_phi", TsallisBlastWave_Integrand_phi, -TMath::Pi(), TMath::Pi(), 0);
    if (!fTsallisBlastWave_Integrand_y)
        fTsallisBlastWave_Integrand_y = new TF1("fTsallisBlastWave_Integrand_y", TsallisBlastWave_Integrand_y, -0.5, 0.5, 0);

    fTsallisBlastWave_Integrand_r->SetParameters(mt, pt, beta_max, temp, n, q, 0., 0.);
    Double_t integral = fTsallisBlastWave_Integrand_y->Integral(-0.5, 0.5);
    return norm * pt * integral;
    }

    TF1 *
    TsallisBlastWave(const Char_t *name, Double_t mass, Double_t xmin =0., Double_t xmax = 10.,  Double_t beta_max = 0.9, Double_t temp = 0.1, Double_t n = 1., Double_t q = 2., Double_t norm = 1.e6)
    {
    
    TF1 *fTsallisBlastWave = new TF1(name, TsallisBlastWave_Func, xmin, xmax, 6);
    fTsallisBlastWave->SetParameters(mass, beta_max, temp, n, q, norm);
    fTsallisBlastWave->SetParNames("mass", "beta_max", "T", "n", "q", "norm");
    fTsallisBlastWave->FixParameter(0, mass);
    fTsallisBlastWave->SetParLimits(1, 0.01, 0.99);
    fTsallisBlastWave->SetParLimits(2, 0.01, 1.);
    fTsallisBlastWave->SetParLimits(3, 0.1, 10.);
    fTsallisBlastWave->SetParLimits(4, 1., 10.);
    return fTsallisBlastWave;
    }

    /*****************************************************************/
    /*****************************************************************/
    /*****************************************************************/


    TF1 *
    BGBlastWave_SingleFit(TH1 *h, Double_t mass, Option_t *opt = "", Double_t xmin = 0., Double_t xmax = 10.)
    {

    TF1 *f = BGBlastWave(Form("fBGBW_%s", h->GetName()), mass, xmin, xmax);
    //   h->Fit(f);
    //   h->Fit(f);
    h->Fit(f, opt);
    return f;
    
    }

    TF1 *
    TsallisBlastWave_SingleFit(TH1 *h, Double_t mass, Option_t *opt = "", Double_t xmin = 0., Double_t xmax = 10.)
    {

    TF1 *f = TsallisBlastWave(Form("fTsBW_%s", h->GetName()), mass, xmin, xmax);
    //   h->Fit(f);
    //   h->Fit(f);
    h->Fit(f, opt);
    return f;
    
    }


    TF1 *
    BGBlastWave_SingleFit_OneOverPT(TObject *h, Double_t mass, Option_t *opt = "", Double_t xmin = 0., Double_t xmax = 10.)
    {
    TF1 *f = BGBlastWave_OneOverPT(Form("fBGBW_%s", h->GetName()), mass,xmin,xmax);
    TString ClassName = h->ClassName();
    if(ClassName.BeginsWith("TH1")){
        TH1 *Obj = (TH1*)h;
    //    h->Fit(f);
    //    h->Fit(f);
        Obj->Fit(f, opt);
    }
    if(ClassName.BeginsWith("TGraphErrors")){
        TGraphErrors *Obj = (TGraphErrors*)h;
    //    h->Fit(f);
    //    h->Fit(f);
        Obj->Fit(f, opt);
    }
    if(ClassName.BeginsWith("TGraphAsymmErrors")){
        TGraphAsymmErrors *Obj = (TGraphAsymmErrors*)h;
    //    h->Fit(f);
    //    h->Fit(f);
        Obj->Fit(f, opt);
    }
    return f;
    
    }

    void GetYieldAndMean(TH1 *h, TF1 *f, Double_t &yield, Double_t &yielderr, Double_t &yielderrcorr, Double_t &mean, Double_t &meanerr, Double_t &meanerrcorr, Double_t min, Double_t max, Double_t *partyield, Double_t *partyielderr, Double_t *partyielderrcorr)
    {
    
    /* find lowest edge in histo */
    Int_t binlo = 0;
    Double_t lo = 0;
    for (Int_t ibin = 1; ibin < h->GetNbinsX() + 1; ibin++) {
        if (h->GetBinContent(ibin) != 0.) {
        binlo = ibin;
        lo = h->GetBinLowEdge(ibin);
        break;
        }
    }
    
    /* find highest edge in histo */
    Int_t binhi = 0;
    Double_t hi = 0;
    for (Int_t ibin = h->GetNbinsX(); ibin > 0; ibin--) {
        if (h->GetBinContent(ibin) != 0.) {
        binhi = ibin + 1;
        hi = h->GetBinLowEdge(ibin + 1);
        break;
        }
    }
    
    /* integrate the data */
    Double_t cont, err, width, cent, integral_data = 0., integralerr_data = 0., integralerrcorr_data = 0., meanintegral_data = 0., meanintegralerr_data = 0., meanintegralerrcorr_data = 0.;
    for (Int_t ibin = binlo; ibin < binhi; ibin++) {
        cent = h->GetBinCenter(ibin);
        width = h->GetBinWidth(ibin);
        cont = h->GetBinContent(ibin);
        err = h->GetBinError(ibin);
        /* check we didn't get an empty bin in between */
        if (cont != 0. && err != 0.) {
        /* all right, use data */
        integral_data += cont * width;
        integralerr_data += err * err * width * width;
        integralerrcorr_data += err * width;
        meanintegral_data += cont * width * cent;
        meanintegralerr_data += err * err * width * width * cent * cent;
        meanintegralerrcorr_data += err * width * cent;
        }
        else {
        /* missing data-point, complain and use function */
        printf("WARNING: missing data-point at %f\n", cent);
        printf("         using function as a patch\n");
        integral_data += f->Integral(h->GetBinLowEdge(ibin), h->GetBinLowEdge(ibin+1));
        integralerr_data += f->IntegralError(h->GetBinLowEdge(ibin), h->GetBinLowEdge(ibin+1), 0, 0, 1.e-6);
        meanintegral_data += f->Mean(h->GetBinLowEdge(ibin), h->GetBinLowEdge(ibin+1)) * f->Integral(h->GetBinLowEdge(ibin), h->GetBinLowEdge(ibin+1));
        meanintegralerr_data += f->Mean(h->GetBinLowEdge(ibin), h->GetBinLowEdge(ibin+1)) * f->IntegralError(h->GetBinLowEdge(ibin), h->GetBinLowEdge(ibin+1), 0, 0, 1.e-6);
        }
    }
    integralerr_data = TMath::Sqrt(integralerr_data);
    meanintegralerr_data = TMath::Sqrt(meanintegralerr_data);
    
    /* integrate below the data */
    Double_t integral_lo = min < lo ? f->Integral(min, lo) : 0.;
    Double_t integralerr_lo = min < lo ? f->IntegralError(min, lo, 0, 0, 1.e-6) : 0.;
    Double_t meanintegral_lo = min < lo ? f->Mean(min, lo) * integral_lo : 0.;
    Double_t meanintegralerr_lo = min < lo ? f->Mean(min, lo) * integralerr_lo : 0.;
    
    /* integrate above the data */
    Double_t integral_hi = max > hi ? f->Integral(hi, max) : 0.;
    Double_t integralerr_hi = max > hi ? f->IntegralError(hi, max, 0, 0, 1.e-6) : 0.;
    Double_t meanintegral_hi = max > hi ? f->Mean(hi, max) * integral_hi : 0.;
    Double_t meanintegralerr_hi = max > hi ? f->Mean(hi, max) * integralerr_hi : 0.;

    /* compute integrated yield */
    yield = integral_data + integral_lo + integral_hi;
    yielderr = TMath::Sqrt(integralerr_data * integralerr_data + 
            integralerr_lo * integralerr_lo + 
            integralerr_hi * integralerr_hi);
    yielderrcorr = TMath::Sqrt(integralerrcorr_data * integralerrcorr_data + 
                integralerr_lo * integralerr_lo + 
                integralerr_hi * integralerr_hi);
    
    /* compute integrated mean */
    mean = (meanintegral_data + meanintegral_lo + meanintegral_hi) / yield;
    meanerr = TMath::Sqrt(meanintegralerr_data * meanintegralerr_data + 
            meanintegralerr_lo * meanintegralerr_lo + 
            meanintegralerr_hi * meanintegralerr_hi) / yield;
    meanerrcorr = TMath::Sqrt(meanintegralerrcorr_data * meanintegralerrcorr_data + 
                meanintegralerr_lo * meanintegralerr_lo + 
                meanintegralerr_hi * meanintegralerr_hi) / yield;

    /* set partial yields */
    partyield[0] = integral_data;
    partyielderr[0] = integralerr_data;
    partyielderrcorr[0] = integralerrcorr_data;
    partyield[1] = integral_lo;
    partyielderr[1] = integralerr_lo;
    partyielderrcorr[1] = integralerr_lo;
    partyield[2] = integral_hi;
    partyielderr[2] = integralerr_hi;
    partyielderrcorr[2] = integralerr_hi;

    }

    /*****************************************************************/

    void IntegratedProduction(TH1 *h, Double_t mass, Option_t *opt = "", Double_t xmin = 0., Double_t xmax = 10.)
    {

    Double_t yield, yielderr, yielderrcorr, mean, meanerr, meanerrcorr, partyield[3], partyielderr[3], partyielderrcorr[3];
    TF1 *f = BGBlastWave_SingleFit(h, mass, opt,xmin,xmax);
    GetYieldAndMean(h, f, yield, yielderr, yielderrcorr, mean, meanerr, meanerrcorr, 0., 10., partyield, partyielderr, partyielderrcorr);

    //  Double_t fint = f->Integral(0.,10.);
    //  Double_t finte = f->IntegralError(0.,10.);
    //  Double_t fmean = f->Mean(0., 10.);

    printf("dN/dy        = %f +- %f (%f)\n", yield, yielderr, yielderrcorr);
    printf("<pt>         = %f +- %f (%f)\n", mean, meanerr, meanerrcorr);
    printf("dN/dy (data) = %f +- %f (%f)\n", partyield[0], partyielderr[0], partyielderrcorr[0]);
    printf("dN/dy (low)  = %f +- %f (%f)\n", partyield[1], partyielderr[1], partyielderrcorr[1]);
    printf("dN/dy (high) = %f +- %f (%f)\n", partyield[2], partyielderr[2], partyielderrcorr[2]);
        
    //  printf("----\n");
    //  printf("dN/dy (func) = %f +- %f\n", fint, finte);
    //  printf("<pT> (func)  = %f +- %f\n", fmean, 0.);
    
    //  TH1 *hr = (TH1 *)h->Clone("hr");
    //  hr->Divide(f);
    //  new TCanvas;
    //  hr->Draw();

    //  TProfile *p = new TProfile("p", "", 100, 0., 10.);
    //  gROOT->LoadMacro("HistoUtils.C");
    //  HistoUtils_Function2Profile(f, p);
    //  p->Draw();
    }


    /*****************************************************************/

    Double_t 
    y2eta(Double_t pt, Double_t mass, Double_t y){
    Double_t mt = TMath::Sqrt(mass * mass + pt * pt);
    return TMath::ASinH(mt / pt * TMath::SinH(y));
    }
    Double_t 
    eta2y(Double_t pt, Double_t mass, Double_t eta){
    Double_t mt = TMath::Sqrt(mass * mass + pt * pt);
    return TMath::ASinH(pt / mt * TMath::SinH(eta));
    }



    void GetFitParameter(TString type, TString centrality, Double_t* params){
    if(!centrality.CompareTo("pp")){
        if(!type.CompareTo("dmtsal")){
            params[0] = 2.1064;
            params[1] = 0.134723;
            params[2] = 5.15457;
            params[3] = 262.57;
            params[4] = 2667.16;
            params[5] = 4.79596;
            params[6] = -2.77055;
        }
    }
    if(!centrality.CompareTo("0-40%")){
        if(!type.CompareTo("oHag")){
            params[0] = 2.05231e+04;
            params[1] = -2.63361e-01;
            params[2] = 6.33381e-01;
            params[3] = 4.62521e-01;
            params[4] = 6.74947e+00;
        }
        if(!type.CompareTo("qcd")){
            params[0] = 16.9746;
            params[1] = 5.92734;
            params[2] = -21.6339;
            params[3] = 2.56968;
            params[4] = 7.64226;
        }
        if(!type.CompareTo("dmtsal")){
            params[0] = 33.376;
            params[1] = 1.57542;
            params[2] = 6.73397;
            params[3] = 33.5075;
            params[4] = -8.48103;
            params[5] = 3.04902;
            params[6] = -1.40645;
        }
    }
    else if(!centrality.CompareTo("0-5%")){
        if(!type.CompareTo("oHag")){
            params[0] = 153111;
            params[1] = -0.567233;
            params[2] = 0.648865;
            params[3] = 0.416557;
            params[4] = 7.15563;

        }
        if(!type.CompareTo("qcd")){
            params[0] = 2.71598e+01;
            params[1] = 5.94784e+00;
            params[2] = -3.36471e+01;
            params[3] = 2.87856e+00;
            params[4] = 1.23700e+01;
        }
        if(!type.CompareTo("rad")){
            params[0] = 89.3562;
            params[1] = 19.9729;
            params[2] = 0.211878;
            params[3] = 3.14135;
            params[4] = 7.80294;
        }
    }
    else if(!centrality.CompareTo("5-10%")){
        if(!type.CompareTo("oHag")){
            params[0] = 1.61223e+05;
            params[1] = -6.21140e-01;
            params[2] = 6.23594e-01;
            params[3] = 3.92189e-01;
            params[4] = 6.99367e+00;
        }
        if(!type.CompareTo("qcd")){
            params[0] = 2.27531e+01;
            params[1] = 5.93873e+00;
            params[2] = -4.08455e+01;
            params[3] = 2.82786e+00;
            params[4] = 1.65429e+01;
        }
        if(!type.CompareTo("rad")){
            params[0] = 85.5714;
            params[1] = 14.3427;
            params[2] = 0.186864;
            params[3] = 3.8672;
            params[4] = 7.81067;
        }
        if(!type.CompareTo("xqcd")){
            params[0] = 22.5395;
            params[1] = 5.88241;
            params[2] = -2.53228;
            params[3] = 0.0322852;
            params[4] = 0.0640679;
            params[5] = 2.84297;
        }
    }
    else if(!centrality.CompareTo("0-20%")){
        if(!type.CompareTo("oHag")){
            params[0] = 6.76330e+04;
            params[1] = -5.37960e-01;
            params[2] = 7.21961e-01;
            params[3] = 4.20035e-01;
            params[4] = 6.88750e+00;
        }
        if(!type.CompareTo("qcd")){
            params[0] = 21.5687;
            params[1] = 5.99981;
            params[2] = -20.0468;
            params[3] = 2.37839;
            params[4] = 6.21999;
        }
        if(!type.CompareTo("rad")){
            params[0] = 87.8356;
            params[1] = 15.3679;
            params[2] = 0.186791;
            params[3] = 2.45483;
            params[4] = 7.84368;
        }
        if(!type.CompareTo("xqcd")){
            params[0] = 21.0163;
            params[1] = 5.76098;
            params[2] = -174.133;
            params[3] = 4.20146;
            params[4] = 62.0668;
            params[5] = 0.605968;
        }
        if(!type.CompareTo("dmtsal")){
            params[0] = 122.169;
            params[1] = 0.916029;
            params[2] = 6.02524;
            params[3] = 27.0573;
            params[4] = -5.84535;
            params[5] = 3.9012;
            params[6] = -2.68122;
        }
    }
    else if(!centrality.CompareTo("0-10%")){
        if(!type.CompareTo("oHag")){
            params[0] = 30463.8;
            params[1] = -0.185472;
            params[2] = 0.488276;
            params[3] = 0.506966;
            params[4] = 7.06264;
        }
        if(!type.CompareTo("qcd")){
            params[0] = 24.8867;
            params[1] = 5.9897;
            params[2] = -36.6162;
            params[3] = 2.83522;
            params[4] = 13.6085;
        }
        if(!type.CompareTo("xqcd")){
            params[0] = 24.2976;
            params[1] = 5.82612;
            params[2] = -435.844;
            params[3] = 4.69103;
            params[4] = 142.667;
            params[5] = 0.738401;
        }
        if(!type.CompareTo("rad")){
            params[0] = 79.2539;
            params[1] = 19.0808;
            params[2] = 0.211468;
            params[3] = 3.27727;
            params[4] = 7.79216;
        }
        if(!type.CompareTo("dmtsal")){
            params[0] = 1549.6;
            params[1] = 0.419372;
            params[2] = 5.5094;
            params[3] = 18.8853;
            params[4] = -3.3233;
            params[5] = 4.5521;
            params[6] = -3.47764;
        }
    }
    else if(!centrality.CompareTo("10-20%")){
        if(!type.CompareTo("oHag")){
            params[0] = 814.152;
            params[1] = 0.549949;
            params[2] = 0.314121;
            params[3] = 0.757416;
            params[4] = 6.7398;
        }
        if(!type.CompareTo("qcd")){
            params[0] = 19.4279;
            params[1] = 6.0361;
            params[2] = -10.7205;
            params[3] = 1.90016;
            params[4] = 2.65789;
        }
        if(!type.CompareTo("xqcd")){
            params[0] = 22.5395;
            params[1] = 5.88241;
            params[2] = -2.53228;
            params[3] = 0.0322852;
            params[4] = 0.0640679;
            params[5] = 2.84297;
        }
        if(!type.CompareTo("dmtsal")){
            params[0] = 2.09114;
            params[1] = 2.7482; 
            params[2] = 6.911;
            params[3] = 51.1383;
            params[4] = -10.6896;
            params[5] = 1.30818;
            params[6] = -1.59137;
        }
    }
    else if(!centrality.CompareTo("20-40%")){
        if(!type.CompareTo("oHag")){
            params[0] = 600.761;
            params[1] = 0.615125;
            params[2] = -0.238326;
            params[3] = 0.724804;
            params[4] = 6.64027;
        }
        if(!type.CompareTo("qcd")){
            params[0] = 11.8906;
            params[1] = 5.9101;
            params[2] = -5.0689;
            params[3] = 2.02887;
            params[4] = 3.4672;
        }
        if(!type.CompareTo("rad")){
            params[0] = 6.89306e+01;
            params[1] = 1.15234e+01;
            params[2] = 1.48406e-01;
            params[3] = 2.59509e+00;
            params[4] = 7.56185e+00;
        }
        if(!type.CompareTo("xqcd")){
            params[0] = 11.4355;
            params[1] = 5.72001;
            params[2] = -667.524;
            params[3] = 4.83608;
            params[4] = 254.033;
            params[5] = 0.925496;
            /* params[0] = 11.8869; */
            /* params[1] = 5.60128; */
            /* params[2] = -24570.1; */
            /* params[3] = 7.79897; */
            /* params[4] = 9435.37; */
            /* params[5] = 1.17752; */
        }
        if(!type.CompareTo("dmtsal")){
            params[0] = 3.14272;
            params[1] = 2.8254;
            params[2] = 7.47365;
            params[3] = 39.5217;
            params[4] = -9.03567;
            params[5] = 0.978345;
            params[6] = -0.821876;
        }
    }
    else if(!centrality.CompareTo("20-50%")){
        if(!type.CompareTo("oHag")){
            params[0] = 600.761;
            params[1] = 0.615125;
            params[2] = -0.238326;
            params[3] = 0.724804;
            params[4] = 6.64027;
        }
        if(!type.CompareTo("qcd")){
            params[0] = 11.8906;
            params[1] = 5.9101;
            params[2] = -5.0689;
            params[3] = 2.02887;
            params[4] = 3.4672;
        }
        if(!type.CompareTo("rad")){
            params[0] = 6.89306e+01;
            params[1] = 1.15234e+01;
            params[2] = 1.48406e-01;
            params[3] = 2.59509e+00;
            params[4] = 7.56185e+00;
        }
        if(!type.CompareTo("xqcd")){
            params[0] = 11.4355;
            params[1] = 5.72001;
            params[2] = -667.524;
            params[3] = 4.83608;
            params[4] = 254.033;
            params[5] = 0.925496;
            /* params[0] = 11.8869; */
            /* params[1] = 5.60128; */
            /* params[2] = -24570.1; */
            /* params[3] = 7.79897; */
            /* params[4] = 9435.37; */
            /* params[5] = 1.17752; */
        }
        if(!type.CompareTo("dmtsal")){
            params[0] = 3.14272;
            params[1] = 2.8254;
            params[2] = 7.47365;
            params[3] = 39.5217;
            params[4] = -9.03567;
            params[5] = 0.978345;
            params[6] = -0.821876;
        }
    }
    else if(!centrality.CompareTo("40-80%")){
        if(!type.CompareTo("oHag")){
            params[0] = 1.30488e+03;
            params[1] = 1.22506e-01;
            params[2] = 4.29101e-01;
            params[3] = 4.75866e-01;
            params[4] = 6.29019e+00;
        }
        if(!type.CompareTo("qcd")){
            params[0] = 2.61438;
            params[1] = 5.77234;
            params[2] = -7.78985;
            params[3] = 1.70584;
            params[4] = 2.36665;
        }
        if(!type.CompareTo("xqcd")){
            params[0] = 2.47731;
            params[1] = 5.95872;
            params[2] = -4.1513;
            params[3] = 1.1465;
            params[4] = 0.963824;
            params[5] =  -1.50422;


    /* params[0] = 2.58789; */
    /*       params[1] = 6.3124; */
    /*       params[2] = -29333.9; */
    /*       params[3] = 1.10274; */
    /*       params[4] = 9993.26; */
    /*       params[5] = 0.644336; */
        }
        if(!type.CompareTo("dmtsal")){
            /* params[0] = 37.1902; */
            /* params[1] = 0.222623; */
            /* params[2] = 5.03864; */
            /* params[3] = 138.652; */
            /* params[4] = 289.536; */
            /* params[5] = 3.36717; */
            /* params[6] = -2.33195; */

            params[0] = 14.4026;
            params[1] = 0.404077;
            params[2] = 5.09444;
            params[3] = 120.794;
            params[4] = 29.8208;
            params[5] = 3.64673;
            params[6] = -2.20963;
        }
    }
    else if(!centrality.CompareTo("60-80%")){
        if(!type.CompareTo("oHag")){
            params[0] = 2.56880e+01;
            params[1] = 8.03456e-01;
            params[2] = -7.88919e-16;
            params[3] = 8.27464e-01;
            params[4] = 6.42907e+00;
        }
        if(!type.CompareTo("qcd")){
            params[0] = 0.988512;
            params[1] = 7.74394;
            params[2] = -0.533419;
            params[3] = 0.0559894;
            params[4] = -0.887449;
        }
    }
}
#endif