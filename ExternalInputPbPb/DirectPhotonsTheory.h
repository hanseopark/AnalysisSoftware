/*
 *    DirectPhotonsTheory
 *
 * Created by Sami Räsänen on 9th February 2014
 * Contact: sami.s.rasanen@jyu.fi
 *
 * This class contains theory results from the following papers:
 *
 * Thermal photons and thermal photon v2 from ideal hydro, no fluctuations, no viscosity:
 *    Holopainen, Räsänen and Eskola, Phys.Rev. C84 (2011) 064903
 *
 * Thermal photons from event-by-event ideal hydro, i.e. fluctuation included but no viscosity
 *    Chatterjee, Holopainen, Renk and Eskola, Phys. Rev. C85 (2012) 064910
 *
 * pQCD photons (prompt + fragmentation) using impact parameter dependent nPDF's, EPS09s.
 * (For details of EPS09s, see Helenius, Eskola, Honkanen and Salgado: JHEP 1207 (2012) 073)
 * Calculation itself is presented in paper
 *    Helenius, Eskola and Paukkunen, JHEP 1305 (2013) 030
 *
 * and the results are discussed together with thermal production in
 *    Chatterjee, Holopainen, Helenius, Renk and Eskola, Phys. Rev. C88 (2013) 034901
 *
 * Thermal photon v2 and v3 from event-by-event ideal hydro
 *    Chatterjee, Srivastava and Renk, arXiv:1401.7464 [hep-ph]
 *
 * Note that v2, presented in paper above, was first discussed in PRC88 (2913) 034901.
 *
 * More detailed comments on calculations below
 *
 */

#ifndef DIRECTPHOTONSTHEORY_H
#define DIRECTPHOTONSTHEORY_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>

class DirectPhotonsTheory{

public:
	// Constructor - all data points in the .cxx file
	DirectPhotonsTheory(void);
	
	// Destructor
	~DirectPhotonsTheory(void){;}
	
   /* =========================
	* H. Holopainen, S.S. Räsänen and K.J. Eskola
	* Phys.Rev. C84 (2011) 064903
    *
    * TGraph's gThermal(...) give invariant yield of thermal photons,
    * E dN/d^3p, at y = 0 as a function of pT. TGraph's gv2(...)
    * give the pT dependence of the v2. Extension spesifies the
    * the centrality class of the results.
    *
    * "Basic" 2+1 -dimensional ideal hydro evolution, no fluctuations
    * in the initial state, nor the viscosity in the evolution
    *
    * The default calculation is done with resummed photon
    * emission rates in plasma by Arnold, Moore and Yaffe (AMY),
    * emission rates by in hadron gas by Turbide, Rapp and Gale (TRG).
    * Equation of state is a thermodynamically consistent parametrization
    * to lattice results (EoS L) by Huovinen and Petreczky. In the case
    * of EoS L, the default switching temperature between rates was
    * chosen to Tswitch = 170 MeV, i.e. below 100 % hadron gas and above
    * 100 % plasma. The default results are without "extension" in
    * the TGraph names.
    *
    * To estimate the spread in theoretical results, the default
    * results are compared to case where switching temperature is
    * raised to 200 MeV (extension "Tc200") and to case where, instead
    * of EoS L, one uses a traditional bag equation with 1st order phase
    * transition at Tc = 165 MeV (EoS Bag), in which case in the
    * mixed phase the plasma-hadron gas ratio is given by Gibbs criteria.
    * The calculation with EoS Bag is done with "old" hadron gas emission
    * rates (R92) by Kapusta, et al.
    *
    * The spectra resulting from the different combinations of rate and EoS
    * are inbetween the extremes given by EoS L with Tswitch = 200 + TRG and
    * EoS Bag + R92.
    *
    * Although EoS bag + R92 gives the largest v2 of photons, the main
    * results of this paper are from the default calculation
    *
    * The results for 00-20 % and 20-40 % were in the paper,
    * predictions before the data came out. 00 - 40 % is the
    * same calculation done for the preliminary data
    */
    
    // 00 - 20 %
	TGraph *gThermalHolopainen00to20;
	TGraph *gThermalHolopainen00to20Tc200;
	TGraph *gThermalHolopainen00to20Bag;
	TGraph *gv2Holopainen00to20;
	TGraph *gv2Holopainen00to20Tc200;
	TGraph *gv2Holopainen00to20Bag;
    // 20 - 40 %
	TGraph *gThermalHolopainen20to40;
	TGraph *gThermalHolopainen20to40Tc200;
	TGraph *gThermalHolopainen20to40Bag;
	TGraph *gv2Holopainen20to40;
	TGraph *gv2Holopainen20to40Tc200;
	TGraph *gv2Holopainen20to40Bag;
    // 00 - 40 %
    TGraph *gThermalHolopainen00to40;
	TGraph *gThermalHolopainen00to40Tc200;
	TGraph *gThermalHolopainen00to40Bag;
	TGraph *gv2Holopainen00to40;
	TGraph *gv2Holopainen00to40Tc200;
	TGraph *gv2Holopainen00to40Bag;
	
   /* =========================
    * Chatterjee, Holopainen, Renk and Eskola
    *    Phys. Rev. C85 (2012) 064910
    *
    * TGraph's gThermal(...) consist invariant spectra for thermal photons,
    * now for 2+1 -dimensional event-by-event hydro evolution includes
    * fluctuations at the initial state, but no viscosity.
    *
    * The calculation implements the same lattice Equation of State
    * (EoS L), plasma emission rate AMY and hadron gas rate TRG as
    * descibed above.
    *
    * 0 - 20 %, sigma=0.4 fm (MC Glauber), tau_0=0.14 fm
    *
    * Note: the points are significantly more sparse as compared to above, 
    * because numerics require significantly more CPU -time. Hence
    * recommend to plot as marker + line.
    *
    */
    
    TGraph *gThermalChatterjee00to20;
    TGraph *gThermalChatterjee20to40;
    TGraph *gThermalChatterjee40to60;
    TGraph *gThermalChatterjee00to40;
    TGraph *gpQCDChatterjee00to40;
    TGraph *gDirectChatterjee00to40;
    
   /* =========================
    * Chatterjee, Srivastava and Renk
    *    arXiv:1401.7464 [hep-ph]
    *
    * TGraph's gv2(...) and gv3(...) give the thermal photon v2 and v3.
    * Following extension are used:
    *    PPFIC = wrt Participant Plane in case of Fluctuating Initial Conditions
    *    RPFIC = wrt Reaction Plane in case of Fluctuating Initial Conditions
    *    SIC   = Smooth Initial Conditions
    *
    * The hydro setup used is the same as in paper above (PRC85(2012)064910),
    * i.e. event-by-event by no viscosity.
    *
    * Note: the paper discussed also results with more recent NLO -emission 
    * rates at plasma presented in JHEP 1305 (2013) 010. So far these NLO rates 
    * are included into SIC -case, others are with AMY rate.
    *
    * Again, points sparse.
    *
    */
    
    TGraph *gv2ChatterjeePPFIC00to40;
    TGraph *gv2ChatterjeeRPFIC00to40;
    TGraph *gv2ChatterjeeSIC00to40;
    TGraph *gv3ChatterjeePPFIC00to40;
    
   /* =========================
    * Helenius, Eskola, Paukkunen
    *      JHEP 1305 (2013) 030
    *
    * (See also
    *   Chatterjee, Holopainen, Helenius, Renk and Eskola, Phys. Rev. C88 (2013) 034901)
    *
    * pQCD photons (=prompt + fragmentation) using
    * impact parameter dependent EPS09s parametrizations
    * for nuclear effects.
    *
    * Description:
    * NLO pQCD calculations for invariant differential yield of inclusive pQCD photons (=prompt + fragmentation)
    * production in Pb+Pb collisions with \sqrt{s}=2.76 TeV at y=0, calculated with INCNLO.
    * Based on the calculations published in JHEP 1305 (2013) 030 which are applied also
    * in Phys.Rev. C88 (2013) 034901 (Fig. 10).
    *
    * PDFs: CTEQ6.6m with EPS09s nuclear modifications
    * FFs (for fragmentation contribution): BFG set II
    * Scale choices: \mu_fact = \mu_ren = \mu_frag = p_T
    * sigma_inel^pp = 64 mb
    *
    * Centrality class = 0-20%, <N_bin> = 1255.82
    * Centrality class = 20-40%, <N_bin> = 430.999
    * Centrality class = 40-80%, <N_bin> = 66.1789
    *
    * Upper and lower limits here represent here uncertanties
    * originating from EPS09 error sets
    *
    */
    
    // Central PDF set
    TGraph *gpQCDHelenius00to20;
    TGraph *gpQCDHelenius20to40;
    TGraph *gpQCDHelenius40to80;
    // EPS09 error sets - high limit
    TGraph *gpQCDHeleniusHigh00to20;
    TGraph *gpQCDHeleniusHigh20to40;
    TGraph *gpQCDHeleniusHigh40to80;
    // EPS09 error sets - low limit
    TGraph *gpQCDHeleniusLow00to20;
    TGraph *gpQCDHeleniusLow20to40;
    TGraph *gpQCDHeleniusLow40to80;
    
};

#endif