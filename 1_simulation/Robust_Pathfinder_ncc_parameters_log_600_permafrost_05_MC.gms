** Copyright IIASA (International Institute for Applied Systems Analysis), 2019-2023
** Contributor(s): Côme Cheritel (come.cheritel@psemail.eu), Armon Rezai (armon.rezai@wu.ac.at), Thomas Gasser (gasser@iiasa.ac.at), Artem Baklanov and Mickael Obersteiner
** This software is a computer program whose purpose is to estimate robust mitigation pathways under a large number of global climate (and environmental) change scenarios.
** This software is governed by the CeCILL license under French law and abiding by the rules of distribution of free software.  You can use, modify and/ or redistribute the software under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at the following URL "http://www.cecill.info". 
** As a counterpart to the access to the source code and rights to copy, modify and redistribute granted by the license, users are provided only with a limited warranty and the software's author, the holder of the economic rights, and the successive licensors have only limited liability. 
** In this respect, the user's attention is drawn to the risks associated with loading, using, modifying and/or developing or reproducing the software by the user in light of its specific status of free software, that may mean that it is complicated to manipulate, and that also therefore means that it is reserved for developers and experienced professionals having in-depth computer knowledge. Users are therefore encouraged to load and test the software's suitability as regards their requirements in conditions enabling the security of their systems and/or data to be ensured and,  more generally, to use and operate it in the same conditions as regards security. 
** The fact that you are presently reading this means that you have had knowledge of the CeCILL license and that you accept its terms.

** This code provides results for:
**  - the Monte Carlo simulation with a constraint limiting permafrost thawing to 50%

** With socio-economic parameters, technological constraints, and damage functions as in Hänsel, M. C., Drupp, M. A., Johansson, D. J., Nesje, F., Azar, C., Freeman, M. C., ... & Sterner, T. (2020). Climate economics support for the UN climate targets. Nature Climate Change, 10(8), 781-789.

** This parameter sets which model variant should be run: ifmod=0 unmodified DICE, ifmod=1 for restriction on abatement and climate damages from Hansel et al (2020).
Parameter ifmod; ifmod = 1;

Option nlp =conopt4;
Option dnlp=conopt4;
Option Reslim =  86400;

$onecho > conopt4.opt
rtredg 3E-13
$offecho

** Defining the different states of the world
set m /1*600 /;
** Defining the different time periods ranging from 2015 to 2510
set t /1*100 /;

variable dt; dt.up = 5;

display t, m, dt.l;

** Defining different boxes of the ocean and permafrost box models
set j /1*5  /,
  pfl /1*3  /;
  
** Loading first set of parameters
set col1;
parameter p1(m,col1);
parameter config(m), phi(m), T2x(m), THs(m), THd(m), th(m), eheat(m), aOHC(m), Lthx(m), lgla_0(m), Ggla_1(m), Ggla_3(m), ggla(m), lgis_0(m), Lgis_1(m), Lgis_3(m), Lais_smb(m), lais_0(m), aais(m), k_toc(m), vgx(m), ggx(m), To(m), bdic(m), gdic(m), npp0(m), vfire(m), vharv(m), vmort(m), vstab(m), vrh1(m), vrh23(m), apass(m), bnpp(m), anpp(m), gnpp(m), bfire(m), gfire(m), brh(m), grh(m), ka(m), ga(m), k_tth(m), Cfr0(m), CO2pi(m), Lgla(m), Lais(m), tgla(m), tgis(m), tais(m), adic(m), aoc_1(m), aoc_2(m), aoc_3(m), aoc_4(m), aoc_5(m), toc_1(m), toc_2(m), toc_3(m), toc_4(m), toc_5(m), vrh3(m), aLST(m), grt1(m), grt2(m), krt(m), amin(m), vthaw(m), vfroz(m), ath_1(m), ath_2(m), ath_3(m), tth_1(m), tth_2(m), tth_3(m), aCO2(m), k_pH(m), T2x0(m) ;


$onecho>par_m.tab
dset=m   rng=Par_v1!a2:a1001     rdim=1
set=col1 rng=Par_v1!b1:bz1               cdim=1
par=p1   rng=par_v1!a1           rdim=1  cdim=1
$offecho
$call gdxxrw i=par_v1_05A_MC.xlsx o=par_m.gdx @par_m.tab

$gdxin par_m.gdx
$load col1, p1
$gdxin


phi(m) = p1(m, "phi");
T2x(m) = p1(m, "T2x");
THs(m) = p1(m, "THs");
THd(m) = p1(m, "THd");
th(m) = p1(m, "th");
eheat(m) = p1(m, "eheat");
aOHC(m) = p1(m, "aOHC");
Lthx(m)= p1(m, "Lthx");
lgla_0(m) = p1(m, "lgla0");
Ggla_1(m) = p1(m, "Ggla1");
Ggla_3(m) = p1(m, "Ggla3");
ggla(m) = p1(m, "ggla");
lgis_0(m) = p1(m, "lgis0");
Lgis_1(m) = p1(m, "Lgis1");
Lgis_3(m) = p1(m, "Lgis3");
Lais_smb(m) = p1(m, "Lais_smb");
lais_0(m) = p1(m, "lais0");
aais(m)= p1(m, "aais");
k_toc(m) = p1(m, "k_toc");
vgx(m) = p1(m, "vgx");
ggx(m) = p1(m, "ggx");
To(m) = p1(m, "To");
bdic(m) = p1(m, "bdic");
gdic(m) = p1(m, "gdic");
npp0(m) = p1(m, "npp0");
vfire(m) = p1(m, "vfire");
vharv(m) = p1(m, "vharv");
vmort(m) = p1(m, "vmort");
vstab(m) = p1(m, "vstab");
vrh1(m) = p1(m, "vrh1");
vrh23(m) = p1(m, "vrh23");
apass(m) = p1(m, "apass");
bnpp(m) = p1(m, "bnpp");
anpp(m) = p1(m, "anpp");
gnpp(m) = p1(m, "gnpp");
bfire(m) = p1(m, "bfire");
gfire(m) = p1(m, "gfire");
brh(m) = p1(m, "brh");
grh(m) = p1(m, "grh");
ka(m) = p1(m, "ka");
ga(m) = p1(m, "ga");
k_tth(m) = p1(m, "k_tth");
Cfr0(m) = p1(m, "Cfr0");
CO2pi(m) = p1(m, "CO2pi");
Lgla(m) = p1(m, "Lgla");
Lais(m) = p1(m, "Lais");
tgla(m) = p1(m, "tgla");
tgis(m) = p1(m, "tgis");
tais(m) = p1(m, "tais");
adic(m) = p1(m, "adic");
aoc_1(m) = p1(m, "aoc_1");
aoc_2(m) = p1(m, "aoc_2");
aoc_3(m) = p1(m, "aoc_3");
aoc_4(m) = p1(m, "aoc_4");
aoc_5(m) = p1(m, "aoc_5");
toc_1(m) = p1(m, "toc_1");
toc_2(m) = p1(m, "toc_2");
toc_3(m) = p1(m, "toc_3");
toc_4(m) = p1(m, "toc_4");
toc_5(m) = p1(m, "toc_5");
vrh3(m) = p1(m, "vrh3");
aLST(m) = p1(m, "aLST");
grt1(m) = p1(m, "grt1");
grt2(m) = p1(m, "grt2");
krt(m) = p1(m, "krt");
amin(m) = p1(m, "amin");
vthaw(m) = p1(m, "vthaw");
vfroz(m) = p1(m, "vfroz");
ath_1(m) = p1(m, "ath_1");
ath_2(m) = p1(m, "ath_2");
ath_3(m) = p1(m, "ath_3");
tth_1(m) = p1(m, "tth_1");
tth_2(m) = p1(m, "tth_2");
tth_3(m) = p1(m, "tth_3");
aCO2(m) = p1(m, "aCO2");
k_pH(m) = p1(m, "k_pH");
T2x0(m) = p1(m, "T2x0");

display phi,T2x,THs,THd,th,eheat,aOHC,Lthx,lgla_0,Ggla_1,Ggla_3,ggla,lgis_0,Lgis_1,Lgis_3,Lais_smb,lais_0,aais,k_toc,vgx,ggx,To,bdic,gdic,npp0,vfire,vharv,vmort,vstab,vrh1,vrh23,apass,bnpp,anpp,gnpp,bfire,gfire,brh,grh,ka,ga,k_tth,Cfr0,CO2pi,Lgla,Lais,tgla,tgis,tais,adic,aoc_1,aoc_2,aoc_3,aoc_4,aoc_5,toc_1,toc_2,toc_3,toc_4,toc_5,vrh3,aLST,grt1,grt2,krt,amin,vthaw,vfroz,ath_1,ath_2,ath_3,tth_1,tth_2,tth_3,aCO2,k_pH,T2x0;



** Loading second set of parameters
set col2;
parameter p2(m,col2);

parameter T_0, d_T_0, CO2_0, d_CO2_0, Td_0, Hgla_0, Hgis_0, Hais_smb_0, Hais_0, Co_1_0, Co_2_0, Co_3_0, Co_4_0, Co_5_0, Cd_0, Cv_0, Cs1_0, Cs2_0, Cs3_0, a_0, Cth_1_0, Cth_2_0, Cth_3_0 ;


$onecho>ini_m.tab
dset=m   rng=ini_2015_v1!a1:a1001 rdim=1
set=col2 rng=ini_2015_v1!b1:x1           cdim=1
par=p2   rng=ini_2015_v1!a1       rdim=1  cdim=1
$offecho
$call gdxxrw i=ini_2015_v1_05A_MC.xlsx o=ini_m.gdx @ini_m.tab

$gdxin ini_m.gdx
$load col2, p2
$gdxin


T_0(m) = p2(m,"T");
d_T_0(m) = p2(m,"d_T");
CO2_0(m) = p2(m,"CO2");
d_CO2_0(m) = p2(m,"d_CO2");
Td_0(m) = p2(m,"Td");
Hgla_0(m) = p2(m,"Hgla");
Hgis_0(m) = p2(m,"Hgis");
Hais_smb_0(m) = p2(m,"Hais_smb");
Hais_0(m) = p2(m,"Hais");
Co_1_0(m) = p2(m,"Co_1");
Co_2_0(m) = p2(m,"Co_2");
Co_3_0(m) = p2(m,"Co_3");
Co_4_0(m) = p2(m,"Co_4");
Co_5_0(m) = p2(m,"Co_5");
Cd_0(m) = p2(m,"Cd");
Cv_0(m) = p2(m,"Cv");
Cs1_0(m) = p2(m,"Cs1");
Cs2_0(m) = p2(m,"Cs2");
Cs3_0(m) = p2(m,"Cs3");
a_0(m) = p2(m,"a");
Cth_1_0(m) = p2(m,"Cth_1");
Cth_2_0(m) = p2(m,"Cth_2");
Cth_3_0(m) = p2(m,"Cth_3");


display T_0, d_T_0, CO2_0, d_CO2_0, Td_0, Hgla_0, Hgis_0, Hais_smb_0, Hais_0, Co_1_0, Co_2_0, Co_3_0, Co_4_0, Co_5_0, Cd_0, Cv_0, Cs1_0, Cs2_0, Cs3_0, a_0, Cth_1_0, Cth_2_0, Cth_3_0 ;


** Loading third and fourth sets of parameters
parameter p3(m,t), p4(m,t) ;

parameter Erfx(m,t), Eluc(m,t) ;



$onecho>ERFx_for_dice.tab
dset=m   rng=ERFx_for_dice_2!a2:a601 rdim=1
dset=t  rng=ERFx_for_dice_2!b1:cw1           cdim=1
par=p3   rng=ERFx_for_dice_2!a1:cw601       rdim=1  cdim=1
$offecho
$call gdxxrw i=ERFx_for_dice_2_05A_MC.xlsx o=ERFx_for_dice.gdx @ERFx_for_dice.tab

$gdxin ERFx_for_dice.gdx
$load p3
$gdxin

display p3;

Erfx(m,t) = p3(m,t);

display Erfx;


$onecho>Eluc_for_dice.tab
dset=m   rng=Eluc_for_dice_2!a2:a601 rdim=1
dset=t  rng=Eluc_for_dice_2!b1:cw1           cdim=1
par=p4   rng=Eluc_for_dice_2!a1:cw601       rdim=1  cdim=1
$offecho
$call gdxxrw i=Eluc_for_dice_2_05A_MC.xlsx o=Eluc_for_dice.gdx @Eluc_for_dice.tab

$gdxin Eluc_for_dice.gdx
$load p4
$gdxin

display p4;

Eluc(m,t) = p4(m,t);

display Eluc;

parameter aoc(m,j), toc(m,j);
aoc(m,"1") = aoc_1(m); aoc(m,"2") = aoc_2(m); aoc(m,"3") = aoc_3(m); aoc(m,"4") = aoc_4(m); aoc(m,"5") = aoc_5(m);
toc(m,"1") = toc_1(m); toc(m,"2") = toc_2(m); toc(m,"3") = toc_3(m); toc(m,"4") = toc_4(m); toc(m,"5") = toc_5(m);

parameter ath(m,pfl),tth(m,pfl);
ath(m,"1") = ath_1(m); ath(m,"2") = ath_2(m); ath(m,"3") = ath_3(m);
tth(m,"1") = tth_1(m); tth(m,"2") = tth_2(m); tth(m,"3") = tth_3(m);

parameters
** Availability of fossil fuels
        fosslim  Maximum cumulative extraction fossil fuels (GtC)  /6000/
**Time Step
        tstep    Years per Period                                    /5/
** Preferences
        elasmu   Elasticity of marginal utility of consumption     /  1.45 /
        prstp    Initial rate of social time preference per year   / .015  /
** Population and technology
        gama     Capital elasticity in production function        /.300    /
        pop0     Initial world population 2015 (millions)         /7403    /
        popadj   Growth rate to calibrate to 2050 pop projection  /0.134   /
        popasym  Asymptotic population (millions)                 /11500   /
        dk       Depreciation rate on capital (per year)          /.100    /
        q0       Initial world gross output 2015 (trill 2010 USD) /105.5   /
        kD0      Initial capital value 2015 (trill 2010 USD)      /223     /
        a0       Initial level of total factor productivity       /5.115    /
        ga0      Initial growth rate for TFP per 5 years          /0.076   /
        dela     Decline rate of TFP per 5 years                  /0.005   /
** Emissions parameters
        gsigma1  Initial growth of sigma (per year)                   /-0.0152 /
        dsig     Decline rate of decarbonization (per period)         /-0.001  /
        miu0     Initial emissions control rate for base case 2015    /.03     /
** These are for declaration and are defined later
        sig0     Carbon intensity 2010 (kgCO2 per output 2005 USD 2010)
** Climate damage parameters
        a10       Initial damage intercept                         /0       /
        a20       Initial damage quadratic term
        a1D       Damage intercept                                 /0       /
        a2D       Damage quadratic term                            /0.00236 /
        a3D       Damage exponent                                  /2.00    /
** Abatement cost
        expcost2  Exponent of control cost function               / 2.6  /
        pback     Cost of backstop 2010$ per tCO2 2015            / 550  /
        gback     Initial cost decline backstop cost per period   / .025 /
        limmiu    Upper limit on control rate after 2150          / 1.2 /
        tnopol    Period before which no emissions controls base  / 45   /
        cprice0   Initial base carbon price (2010$ per tCO2)      / 2    /
        gcprice   Growth rate of base carbon price per year       /.02   /
** newly added parameters
        eff0      Industrial emissions 2015 (GtC per year)        / 9.4    /
        chi       loss aversion parameter in meta model           / 0  /
        a2dmod    damage parameter for modified verion            / 0.007438 /
        elasmumod   Elasticity of marginal utility of consumption     /  1.000001 /
        prstpmod   Initial rate of social time preference per year   / .005  /
;


PARAMETERS
        l(t)          Level of population and labor
        al(t)         Level of total factor productivity
        sigma(m,t)    CO2-equivalent-emissions output ratio
        rr(t)         Average utility social discount rate
        gaD(t)         Growth rate of productivity from
        forcoth(t)    Exogenous forcing for other greenhouse gases
        gl(t)         Growth rate of labor
        gcost1        Growth of cost factor
        gsig(t)       Change in sigma (cumulative improvement of energy efficiency)
        cost1(t)      Adjusted cost for backstop
        partfract(t)  Fraction of emissions in control regime
        lam           Climate model parameter
        gfacpop(t)    Growth factor population
        pbacktime(t)  Backstop price
        optlrsav      Optimal long-run savings rate used for transversality
        scc(m,t)      Social cost of carbon
        cprice(m,t)   carbon price
        cpricebase(t) Carbon price in base case
        photel(t)     Carbon Price under no damages (Hotelling rent condition);

* Further definitions of parameters
        sig0 = 3.666 * Eff0 /(q0*(1-miu0));
        l("1") = pop0;
        loop(t, l(t+1)=l(t););
        loop(t, l(t+1)=l(t)*(popasym/L(t))**popadj ;);

        gaD(t)=ga0*exp(-dela*5*((t.val-1)));
        al("1") = a0; loop(t, al(t+1)=al(t)/((1-gaD(t))););
        gsig("1")=gsigma1; loop(t,gsig(t+1)=gsig(t)*((1+dsig)**tstep) ;);
        sigma(m,"1")=sig0;   loop(t,sigma(m,t+1)=(sigma(m,t)*exp(gsig(t)*tstep)););
        pbacktime(t)=pback*(1-gback)**(t.val-1);
        cost1(t) = pbacktime(t)*(Sum(m,sigma(m,t))/Card(m))/expcost2/1000;

parameter tflip, miulim(t), etree(m,t);
tflip = 30;

** Setting socio-economic, damage and technological parameters as in Hansel et al. (2020)
If((ifmod eq 1),
  a2D = 0.007438;
  tflip = 7;
  elasmu = elasmumod;
  prstp = prstpmod;
  display a2D, tflip;
);

miulim(t) = limmiu; miulim(t)$(t.val<tflip) = 1; miulim("1") = miu0 ;
etree(m,t)= eluc(m,t);

** Equation of the climate model
** This is a translation of the pathfinder model model in GAMS.
** Please cite :
** Bossy, T., Gasser, T., and Ciais, P.: Pathfinder v1.0.1: a Bayesian-inferred simple carbon–climate model to explore climate change scenarios, Geosci. Model Dev., 15, 8831–8868, https://doi.org/10.5194/gmd-15-8831-2022, 2022. 
** https://doi.org/10.5194/gmd-15-8831-2022
** Note that the variable t(m,t) was changed to Tatm(m,t) because t is already the time set.
** The translation is based on an implicit Euler method


variable Tatm(m,t), CO2(m,t) ;

variables Eco2(m,t), pH(m,t), CO2_dot(m,t);

variables ERF(m,t), Td(m,t), RFco2(m,t), T_dot(m,t), Td_dot(m,t);

variables Fland(m,t), Cs1(m,t), Cs2(m,t), Cs3(m,t), Cs(m,t), Cv(m,t), NPP(m,t), Efire(m,t), Eharv(m,t), Fmort(m,t), RH1(m,t), RH2(m,t), RH3(m,t), Fstab(m,t), r_rh(m,t), Fpass(m,t)  ;

variables Cv_dot(m,t), Cs1_dot(m,t), Cs2_dot(m,t), Cs3_dot(m,t);

variables Co_i(m,t, j), Co(m,t), Cd(m,t), pCO2(m,t), Focean(m,t), Co_i_dot(m,t,j), Cd_dot(m,t), dic(m,t), pdic(m,t);

variables a(m,t), abar(m,t), r_rt(m,t),   Cth_i(m,t, pfl), Epf(m,t) , expr(m,t) , a_dot(m,t), Cfr_dot(m,t), Cth_i_dot(m,t, pfl);

variables Hthx(m,t), Hgis(m,t), Hais(m,t), Hais_smb(m,t), Hgla(m, t), Htot(m,t), OHC(m,t), OHC_dot(m,t), Hthx_dot(m,t), Hgis_dot(m,t), Hais_dot(m,t), Hais_smb_dot(m,t), Hgla_dot(m, t), Htot_dot(m,t);



** Atmosperic module

equations Eco2_eq(m,t), pH_eq(m,t), CO2_tr(m,t) ;

model atmos /all/ ;

pH_eq(m,t)..     pH(m,t)         =E= k_pH(m) * (8.5541 - 0.00173 * CO2(m,t) + 1.3264E-6 * CO2(m,t)**2 - 4.4943E-10 * CO2(m,t)**3) ;
Eco2_eq(m,t)..   aCO2(m)*CO2_dot(m,t) =E= (Eco2(m,t) + Epf(m,t) - Fland(m,t) - Focean(m,t) ) ;

CO2_tr(m,t+1)..  CO2(m,t+1)      =E= CO2(m,t) + dt * CO2_dot(m,t+1) ;


** Climate module

equations RFco2_eq(m,t), ERF_eq(m,t), T_dot_eq(m,t), Td_dot_eq(m,t), Tatm_tr(m,t), Td_tr(m,t) ;

model clim /all - atmos/ ;

RFco2_eq(m,t)..  RFco2(m,t)      =E= phi(m) * log(CO2(m,t)/CO2pi(m)) ;
ERF_eq(m,t)..    ERF(m,t)        =E= RFco2(m,t) + ERFx(m,t) ;

T_dot_eq(m,t)..  THs(m) * T_dot(m,t)     =E= ERF(m,t) - phi(m)*log(2)*Tatm(m,t)/T2x(m) - eheat(m)*th(m)*(Tatm(m,t)-Td(m,t)) ;
Td_dot_eq(m,t).. THd(m) * Td_dot(m,t)    =E= th(m)*(Tatm(m,t)-Td(m,t)) ;

Tatm_tr(m,t+1).. Tatm(m,t+1)     =E= Tatm(m,t) + dt * T_dot(m,t+1) ;
Td_tr(m,t+1)..   Td(m,t+1)       =E= Td(m,t)   + dt * Td_dot(m,t+1) ;


** Land module

equations NPP_eq(m,t), r_rh_eq(m,t), Efire_eq(m,t), Eharv_eq(m,t), Fmort_eq(m,t), RH1_eq(m,t), Fstab_eq(m,t), RH2_eq(m,t), Fpass_eq(m,t), RH3_eq(m,t), Fland_eq(m,t), Cv_dot_eq(m,t), Cs1_dot_eq(m,t), Cs2_dot_eq(m,t), Cs3_dot_eq(m,t), Cv_tr(m,t), Cs1_tr(m,t), Cs2_tr(m,t), Cs3_tr(m,t) ;

model land /all - atmos - clim/ ;

NPP_eq(m,t)..    NPP(m,t)        =E= npp0(m)*(1. + bnpp(m) / anpp(m) * (1. - (CO2(m,t) / CO2pi(m)+ 0.00001)**(-anpp(m)))) *(1. + gnpp(m) * Tatm(m,t)) ;

r_rh_eq(m,t)..   r_rh(m,t)       =E= (1. + brh(m) * (Cs1(m,t) / (Cs1(m,t) + Cs2(m,t) + Cs3(m,t)+ 0.00001) * (1. + vstab(m) / vrh23(m)) - 1.)) * exp(grh(m) * Tatm(m,t)) ;
Efire_eq(m,t)..  Efire(m,t)      =E= vfire(m) *Cv(m,t) *(1.+ bfire(m)*(CO2(m,t)/CO2pi(m) - 1.))*(1. + gfire(m)* Tatm(m,t)) ;
Eharv_eq(m,t)..  Eharv(m,t)      =E= vharv(m)*Cv(m,t) ;
Fmort_eq(m,t)..  Fmort(m,t)      =E= vmort(m)*Cv(m,t) ;
RH1_eq(m,t)..    RH1(m,t)        =E= vrh1(m) * r_rh(m,t) * Cs1(m,t) ;
Fstab_eq(m,t)..  Fstab(m,t)       =E= vstab(m) * r_rh(m,t) * Cs1(m,t) ;
RH2_eq(m,t)..    RH2(m,t)        =E= (vrh23(m) - vrh3(m) * apass(m)) / (1. - apass(m)) * r_rh(m,t) * Cs2(m,t) ;
Fpass_eq(m,t)..  Fpass(m,t)      =E= vrh3(m) * apass(m) / (1. - apass(m)) *  r_rh(m,t) * Cs2(m,t) ;
RH3_eq(m,t)..    RH3(m,t)        =E= vrh3(m) * r_rh(m,t)  * Cs3(m,t) ;
Fland_eq(m,t)..  Fland(m,t)      =E= NPP(m,t) - Efire(m,t) - Eharv(m,t) - RH1(m,t) - RH2(m,t) - RH3(m,t)  ;
Cv_dot_eq(m,t).. Cv_dot(m,t)     =E= NPP(m,t) - Efire(m,t) - Eharv(m,t) - Fmort(m,t) ;
Cs1_dot_eq(m,t)..Cs1_dot(m,t)    =E= Fmort(m,t) - Fstab(m,t)- RH1(m,t) ;
Cs2_dot_eq(m,t)..Cs2_dot(m,t)    =E= Fstab(m,t) - Fpass(m,t) - RH2(m,t) ;
Cs3_dot_eq(m,t)..Cs3_dot(m,t)    =E= Fpass(m,t) - RH3(m,t) ;

Cv_tr(m,t+1)..   Cv(m,t+1)       =E= Cv(m,t)   + dt * Cv_dot(m,t+1) ;
Cs1_tr(m,t+1)..  Cs1(m,t+1)      =E= Cs1(m,t)  + dt * Cs1_dot(m,t+1) ;
Cs2_tr(m,t+1)..  Cs2(m,t+1)      =E= Cs2(m,t)  + dt * Cs2_dot(m,t+1) ;
Cs3_tr(m,t+1)..  Cs3(m,t+1)      =E= Cs3(m,t)  + dt * Cs3_dot(m,t+1) ;


** Ocean module

equations Co_eq(m,t), dic_eq(m,t), pdic_eq(m,t), pCO2_eq(m,t), Focean_eq(m,t), Co_ind_eq(m,t, j), Co_i_tr(m,t, j) ;

model ocean /all - atmos - clim - land/ ;

Co_eq(m,t)..     Co(m,t)         =E= sum(j, Co_i(m,t, j) );
dic_eq(m,t)..    dic(m,t)        =E= adic(m) / bdic(m) * Co(m,t) ;
pdic_eq(m,t)..   pdic(m,t)       =E= ((1.5568 - 1.3993E-2 * To(m)) * dic(m,t)
                                    + (7.4706 - 0.20207   * To(m)) * 1E-3  * Power(dic(m,t),2)
                                    - (1.2748 - 0.12015   * To(m)) * 1E-5  * Power(dic(m,t),3)
                                    + (2.4491 - 0.12639   * To(m)) * 1E-7  * Power(dic(m,t),4)
                                    - (1.5768 - 0.15326   * To(m)) * 1E-10 * Power(dic(m,t),5) ) ;
pCO2_eq(m,t)..   pCO2(m,t)       =E= (pdic(m,t) + CO2pi(m)) * exp(gdic(m) * Tatm(m,t)) ;
Focean_eq(m,t).. Focean(m,t)     =E= vgx(m) * (1. + ggx(m) * Tatm(m,t)) *(CO2(m,t) - pCO2(m,t)) ;
Co_ind_eq(m,t,j)..Co_i_dot(m,t,j)=E= - Co_i(m,t, j) /(k_toc(m)*toc(m,j)) + aoc(m,j) * Focean(m,t)  ;

Co_i_tr(m,t+1,j)..Co_i(m,t+1,j)   =E= Co_i(m,t,j) + dt * Co_i_dot(m,t+1,j) ;


** Permafrost module

equations r_rt_eq(m,t), abar_eq(m,t), Epf_eq(m,t) , expr_eq(m,t), a_dot_eq(m,t), Cth_i_dot_eq(m,t, pfl), a_tr(m,t), Cth_i_tr(m,t,pfl) ;

model perma /all - atmos - clim - land - ocean/ ;

r_rt_eq(m,t)..   r_rt(m,t)       =E= exp(krt(m) * grt1(m) * aLST(m) * Tatm(m,t) - krt(m) * grt2(m) * Power(aLST(m) * Tatm(m,t),2) );
expr_eq(m,t)..   expr(m,t)       =E= exp(- ka(m) * ga(m) *  aLST(m) * Tatm(m,t)) ;
abar_eq(m,t)..   abar(m,t)       =E= -amin(m) + (1. + amin(m)) / (1. + ((1 + 1./amin(m))**ka(m) - 1.) * expr(m,t))**(1./ka(m)) ;
Epf_eq(m,t)..    Epf(m,t)        =E= sum(pfl, Cth_i(m,t, pfl)/tth(m,pfl) ) * r_rt(m,t)/k_tth(m) ;
a_dot_eq(m,t)..  a_dot(m,t)      =E= (abar(m,t) - a(m,t)) * (vfroz(m) + ( vthaw(m) - vfroz(m) ) * sigmoid(500 * (abar(m,t) - a(m,t)) ) );
Cth_i_dot_eq(m,t, pfl)..
             Cth_i_dot(m,t, pfl) =E= ath(m,pfl) * a_dot(m,t) * Cfr0(m) -  Cth_i(m,t, pfl)/ tth(m,pfl) * r_rt(m,t)/ k_tth(m) ;

a_tr(m,t+1)..    a(m,t+1)        =E= a(m,t)   + dt * a_dot(m,t+1) ;
Cth_i_tr(m,t+1,pfl)..
              Cth_i(m,t+1,pfl)   =E= Cth_i(m,t,pfl) + dt * Cth_i_dot(m,t+1,pfl) ;

model IIASA /all/ ;


** Sea level rise module

equations OHC_eq(m,t), OHC_dot_eq(m,t), Hthx_eq(m,t), Hthx_dot_eq(m,t), Hgis_dot_eq(m, t), Hais_smb_dot_eq(m,t), Hais_dot_eq(m,t), Hgla_dot_eq(m,t), Hgis_tr(m, t), Hais_smb_tr(m,t), Hais_tr(m,t), Hgla_tr(m,t), Htot_dot_eq(m,t), Htot_eq(m,t);

OHC_eq(m,t)..      OHC(m,t)      =E= aOHC(m) * (THs(m) * Tatm(m,t) + THd(m) * Td(m,t)) ;
OHC_dot_eq(m,t)..  OHC_dot(m,t)  =E= aOHC(m) * (THs(m) * T_dot(m,t) + THd(m) * Td_dot(m,t)) ;

Hthx_eq(m,t)..     Hthx(m,t)     =E= Lthx(m) * OHC(m,t) ;
Hthx_dot_eq(m,t).. Hthx_dot(m,t) =E= Lthx(m) * OHC_dot(m,t) ;

Hgis_dot_eq(m,t)..     Hgis_dot(m,t)     =E= lgis_0(m) + 1/tgis(m) * (Lgis_1(m)*Tatm(m,t) + Lgis_3(m)* power(Tatm(m,t), 3) - Hgis(m,t)) ;
Hais_smb_dot_eq(m,t).. Hais_smb_dot(m,t) =E= - Lais_smb(m) * Tatm(m,t) ;
Hais_dot_eq(m,t)..     Hais_dot(m,t)     =E= Hais_smb_dot(m,t) + lais_0(m) + (Lais(m) * Tatm(m,t) - ( Hais(m,t) - Hais_smb(m,t) ) )/ tais(m) * (1 + aais(m) * (Hais(m,t) - Hais_smb(m,t)) ) ;
Hgla_dot_eq(m,t)..     Hgla_dot(m,t)     =E= lgla_0(m) + exp( ggla(m) * Tatm(m,t) )/tgla(m) * ( Lgla(m) * (1- exp( - Ggla_1(m)*Tatm(m,t) - Ggla_3(m)* power(Tatm(m,t), 3) ) ) - Hgla(m,t) ) ;


Hgis_tr(m,t+1)..      Hgis(m,t+1)     =E= Hgis(m,t)    + dt * Hgis_dot(m,t+1) ;
Hais_smb_tr(m,t+1)..  Hais_smb(m,t+1) =E= Hais_smb(m,t)+ dt * Hais_smb_dot(m,t+1) ;
Hais_tr(m, t+1)..     Hais(m,t+1)     =E= Hais(m,t)    + dt * Hais_dot(m,t+1) ;
Hgla_tr(m,t+1)..      Hgla(m,t+1)     =E= Hgla(m,t)    + dt * Hgla_dot(m,t+1) ;


Htot_dot_eq(m,t).. Htot_dot(m,t) =E= Hthx_dot(m,t) + Hgis_dot(m,t) + Hais_dot(m,t) +Hgla_dot(m,t) ;
Htot_eq(m,t) ..    Htot(m,t)     =E= Hthx(m,t) + Hgis(m,t) + Hais(m,t) + Hgla(m,t) ;


model IIASAplus /all/ ;



** essential DICE equations

VARIABLES
        MIU(m,t)          Emission control rate GHGs
        EIND(m,t)         Industrial emissions (GtCO2 per year)
        C(m,t)            consumption (trillions 2005 US dollars per year)
        KD(m,t)            Capital stock (trillions 2005 US dollars)
        I(m,t)            Investment (trillions 2005 USD per year)
        Y(m,t)            Gross world product net of abatement and damages (trillions 2005 USD per year)
        CCA(m,t)          Cumulative industrial carbon emissions (GTC)
        UTILITY(m)        Welfare function
        METAUTILITY       Meta-Welfare function
;

NONNEGATIVE VARIABLES  MIU, Y, C, KD, I, UTILITY;
VARIABLES a2D_var;

EQUATIONS
        EEQ(m,t)           Emissions equation
        EINDEQ(m,t)        Industrial emissions
        CCACCA(m,t)        Cumulative carbon emissions
        YY(m,t)            Output net equation
        CC(m,t)            Consumption equation
        KK(m,t)            Capital balance equation
        UTIL(m)            Objective function      ;

 eeq(m,t)..             ECO2(m,t)        =E= EIND(m,t) * 1/3.666 + Etree(m,t) ;
 eindeq(m,t)..          EIND(m,t)        =E= sigma(m,t) * (al(t)*(L(t)/1000)**(1-GAMA))*(KD(m,t)**GAMA) * (1-(MIU(m,t)));
 ccacca(m,t+1)..        CCA(m,t+1)       =E= CCA(m,t)+ EIND(m,t)*5/3.666;
 yy(m,t)..              Y(m,t)           =E= (al(t)*(L(t)/1000)**(1-GAMA))*(KD(m,t)**GAMA) * ( 1- ( a1D*TATM(m,t) + a2D_var*Power(TATM(m,t),a3D) ) - (cost1(t) * (MIU(m,t)**expcost2)) );
 cc(m,t)..              C(m,t)           =E= Y(m,t) - I(m,t);
 kk(m,t+1)..            KD(m,t+1)        =E= (1-dk)**tstep * KD(m,t) + tstep * I(m,t);
 util(m)..              UTILITY(m)       =E= SUM(T,  1/( (1+prstp)**( tstep * (ord(T)-1) ) ) * L(T)* log( tstep * C(m,t)*1000/L(T) ) ) ;


model DICE / All - IIASAplus / ;


Equations Miu_Ramp1(m,t), Miu_Ramp2(m,t);
Miu_Ramp1(m,t+1)$(t.val LT tflip+1).. EInd(m,t) - EInd(m,t-1)=G= -10;
Miu_Ramp2(m,t+1)$(t.val GE tflip+1).. miu(m,t)               =L= miu(m,t-1) * 1.1 ;
model Ramp / Miu_Ramp1 + Miu_Ramp2 /;
model DICERamp / DICE, Ramp /;

Tatm.fx(m,"1") = T_0(m);
CO2.fx(m,"1")  = CO2_0(m);
Td.fx(m,"1")   = Td_0(m);

Hgla.fx(m,"1")     = Hgla_0(m);
Hgis.fx(m,"1")     = Hgis_0(m);
Hais.fx(m,"1")     = Hais_0(m);
Hais_smb.fx(m,"1") = Hais_smb_0(m);


Co_i.fx(m,"1","1") = Co_1_0(m);
Co_i.fx(m,"1","2") = Co_2_0(m);
Co_i.fx(m,"1","3") = Co_3_0(m);
Co_i.fx(m,"1","4") = Co_4_0(m);
Co_i.fx(m,"1","5") = Co_5_0(m);


Cd.fx(m,"1") = Cd_0(m);
Cv.fx(m,"1") = Cv_0(m);

Cs1.fx(m,"1") = Cs1_0(m);
Cs2.fx(m,"1") = Cs2_0(m);
Cs3.fx(m,"1") = Cs3_0(m);

a.fx(m,"1") = a_0(m);

Cth_i.fx(m,"1","1") = Cth_1_0(m);
Cth_i.fx(m,"1","2") = Cth_2_0(m);
Cth_i.fx(m,"1","3") = Cth_3_0(m);


** Bounds on drivers and states

ECO2.lo(m,t) = -inf; ECO2.up(m,t) = inf;
ECO2.lo(m,t) = -inf; ECO2.up(m,t) = inf;

If((ifmod eq 1), miu.fx(m,t) = 1  ;
else             miu.fx(m,t) = 1  ;
);
miu.fx(m,"1") = miu0;

CCA.lo(m,t)       = -inf;
CCA.up(m,t)       = fosslim;

KD.LO(m,t)        = 1;
C.LO(m,t)         = 2;


** Initial values

Tatm.fx(m,"1") = T_0(m);
CO2.fx(m,"1") = CO2_0(m);

Co_i.fx(m,"1","1") = Co_1_0(m);
Co_i.fx(m,"1","2") = Co_2_0(m);
Co_i.fx(m,"1","3") = Co_3_0(m);
Co_i.fx(m,"1","4") = Co_4_0(m);
Co_i.fx(m,"1","5") = Co_5_0(m);

Cv.fx(m,"1") = Cv_0(m);
Cs1.fx(m,"1") = Cs1_0(m);
Cs2.fx(m,"1") = Cs2_0(m);
Cs3.fx(m,"1") = Cs3_0(m);

a.fx(m,"1") = a_0(m);
Cth_i.fx(m,"1","1") = Cth_1_0(m);
Cth_i.fx(m,"1","2") = Cth_2_0(m);
Cth_i.fx(m,"1","3") = Cth_3_0(m);

Td.fx(m,"1")   = Td_0(m);
Hgla.fx(m,"1")     = Hgla_0(m);
Hgis.fx(m,"1")     = Hgis_0(m);
Hais.fx(m,"1")     = Hais_0(m);
Hais_smb.fx(m,"1") = Hais_smb_0(m);

CCA.FX(m,"1")    = 400;
KD.FX(m,"1")     = kD0;


** Defining meta-utility over several states of the world

equations MetaUTIL, unimiueq(m,t), uniceq(m,t); variable unimiu(t), MetaUTILITY; positive variable unic(t);
scalar scale1; scale1 = 1;
Metautil..       MetaUTILITY     =E= 1/Card(m) * ( Sum(m, (utility(m)/ scale1)**(1-chi) )  )**(1/(1-chi));
unimiueq(m,T)..                     miu(m,t)        =E= unimiu(t);
uniceq(m,t)$(t.val LT Card(T))..    C(m,t)          =E= Unic(t) ;
equation FOCc(T);

model DICEuni    / DICE + MetaUtil + unimiueq + Uniceq / ;
model DICEuniFOC / DICE + MetaUtil + unimiueq + Uniceq + FOCc / ;
model META / IIASAplus + DICEuni / ;
model MonteCarlo / Meta - unimiueq - uniceq / ;

model METAramp / Meta + Ramp / ;
model MonteCarloRamp / MonteCarlo + Ramp / ;


META.optfile = 1 ;
METAramp.optfile = 1 ;

parameter prob_below, solstat;

parameter lambdaC, lambdaUniC, lambdaE, lambdaUniE, focunic, focmiu1, focmiu2, focunimiu, foclambda1, foclambda, mutil;


** Definition of equations for probabilistic targets
equation MetautilAux, MaxTempEq, MaxTempConEq, MaxAEq, MaxAConEq, MinpHEq, MinpHConEq;
variable MetaUTILITYAux, MaxTempCon, MaxTemp, MaxACon, MaxA, MinpHCon, MinpH;
parameter prob, alpha, targ_temp, targ_pH, targ_a, pen_temp, pen_a, pen_pH, prob_targ ;
alpha = 50; targ_temp = 2.0; targ_a = 0.5;
pen_temp = 0;
pen_a = 0; pen_pH = 0; prob_targ = 1;

MetautilAux..    MetaUTILITYAux     =E= MetaUtility - pen_Temp * Power(MaxTempCon - prob_targ, 2)
                                                    - pen_pH   * Power(MinpHCon   - prob_targ, 2)
                                                    - pen_A    * Power(MaxAcon    - prob_targ, 2)    ;

** Max and min are set implicitely as, in the case of max a value above all values. Then the optimization code is going to find the smallest value above all values i.e. the max
MaxTempEq(m,t).. MaxTemp(m) =G= Tatm(m,t) ;
MaxAEq(m,t)..    MaxA(m)    =G= a(m,t) ;
MinPhEq(m,t)..   MinPH(m)   =L= pH(m,t) ;
MaxTempConEq..   MaxTempCon =E= Sum(m, Sigmoid( 100* (MaxTemp(m) - Targ_Temp ) )) / Card(m) ;
MaxAConEq..      MaxACon    =E= Sum(m, Sigmoid( 100* (MaxA(m)    - Targ_A    ) )) / Card(m) ;
MinPhConEq..     MinpHCon   =E= 1-Sum(m,Sigmoid(250* (MinpH(m)   - Targ_pH(m)) )) / Card(m) ;


model METAmajor / META, MaxTempEq, MaxTempConEq / ;
model METAmajorramp / METAmajor + Ramp / ;

model METAmajorA / META, MaxAEq, MaxAConEq / ;
model METAmajorAramp / METAmajorA + Ramp / ;

model METAmajorpH / META, MinpHEq, MinpHConEq / ;
model METAmajorpHramp / METAmajorpH + Ramp / ;

model METAplus / META, MetautilAux, MaxTempEq, MaxTempConEq, MaxAEq, MaxAConEq, MinpHEq, MinpHConEq / ;
model METAplusramp / METAplus + Ramp / ;

METAmajor.optfile = 1 ;  METAmajorramp.optfile = 1 ;  METAplus.optfile = 1 ;  METAplusramp.optfile = 1 ;
METAmajorA.optfile = 1 ; METAmajorAramp.optfile = 1 ; METAmajorpH.optfile = 1 ; METAmajorpHramp.optfile = 1 ;


** Setting intermediate parameters to load a result file to run the model
parameter Tatm_param(m,t), CO2_param(m,t) ;

parameter ECO2_param(m,t), ph_param(m,t), CO2_dot_param(m,t);

parameter ERF_param(m,t), Td_param(m,t), RFco2_param(m,t), T_dot_param(m,t), Td_dot_param(m,t);

parameter Fland_param(m,t), Cs1_param(m,t), Cs2_param(m,t), Cs3_param(m,t), Cs_param(m,t), Cv_param(m,t), NPP_param(m,t), Efire_param(m,t), Eharv_param(m,t), Fmort_param(m,t), RH1_param(m,t), RH2_param(m,t), RH3_param(m,t), Fstab_param(m,t), r_rh_param(m,t), Fpass_param(m,t)  ;

parameter Cv_dot_param(m,t), Cs1_dot_param(m,t), Cs2_dot_param(m,t), Cs3_dot_param(m,t);

parameter Co_i_param(m,t, j), Co_param(m,t), Cd_param(m,t), pCO2_param(m,t), Focean_param(m,t), Co_i_dot_param(m,t,j), Cd_dot_param(m,t), dic_param(m,t), pdic_param(m,t);

parameter a_param(m,t), abar_param(m,t), r_rt_param(m,t), Cth_i_param(m,t, pfl), Epf_param(m,t), expr_param(m,t), a_dot_param(m,t), Cfr_dot_param(m,t), Cth_i_dot_param(m,t, pfl);

parameter Hthx_param(m,t), Hgis_param(m,t), Hais_param(m,t), Hais_smb_param(m,t), Hgla_param(m, t), Htot_param(m,t), OHC_param(m,t), OHC_dot_param(m,t), Hthx_dot_param(m,t), Hgis_dot_param(m,t), Hais_dot_param(m,t), Hais_smb_dot_param(m,t), Hgla_dot_param(m, t), Htot_dot_param(m,t);



parameter MIU_param(m,t), EIND_param(m,t), C_param(m,t), KD_param(m,t), I_param(m,t), Y_param(m,t), CCA_param(m,t), UTILITY_param(m), METAUTILITY_param;

parameter a2D_var_param;

parameter MetaUTILITYAux_param, MaxTempCon_param, MaxTemp_param(m), MaxACon_param, MaxA_param(m), MinpHCon_param, MinpH_param(m);


** Loading previous result file
$gdxin results2deg95.gdx
$load Tatm_param =  Tatm.l
$load CO2_param =  CO2.l

$load ECO2_param =  ECO2.l
$load pH_param =  pH.l
$load CO2_dot_param =  CO2_dot.l

$load ERF_param =  ERF.l
$load Td_param =  Td.l
$load RFco2_param =  RFco2.l
$load T_dot_param =  T_dot.l
$load Td_dot_param =  Td_dot.l

$load Fland_param =  Fland.l
$load Cs1_param =  Cs1.l
$load Cs2_param =  Cs2.l
$load Cs3_param =  Cs3.l
$load Cs_param =  Cs.l
$load Cv_param =  Cv.l
$load NPP_param =  NPP.l
$load Efire_param =  Efire.l
$load Eharv_param =  Eharv.l
$load Fmort_param =  Fmort.l
$load RH1_param =  RH1.l
$load RH2_param =  RH2.l
$load RH3_param =  RH3.l
$load Fstab_param =  Fstab.l
$load r_rh_param =  r_rh.l
$load Fpass_param =  Fpass.l

$load Cv_dot_param =  Cv_dot.l
$load Cs1_dot_param =  Cs1_dot.l
$load Cs2_dot_param =  Cs2_dot.l
$load Cs3_dot_param =  Cs3_dot.l

$load Co_i_param =  Co_i.l
$load Co_param =  Co.l
$load Cd_param =  Cd.l
$load pCO2_param =  pCO2.l
$load Focean_param =  Focean.l
$load Co_i_dot_param =  Co_i_dot.l
$load Cd_dot_param =  Cd_dot.l
$load dic_param =  dic.l
$load pdic_param =  pdic.l

$load a_param = a.l
$load abar_param = abar.l
$load r_rt_param = r_rt.l
$load Cth_i_param = Cth_i.l
$load Epf_param = Epf.l
$load expr_param = expr.l
$load a_dot_param = a_dot.l
$load Cfr_dot_param = Cfr_dot.l
$load Cth_i_dot_param = Cth_i_dot.l


$load Hthx_param = Hthx.l
$load Hgis_param = Hgis.l
$load Hais_param = Hais.l
$load Hais_smb_param = Hais_smb.l
$load Hgla_param = Hgla.l
$load Htot_param = Htot.l
$load OHC_param = OHC.l
$load OHC_dot_param = OHC_dot.l
$load Hthx_dot_param = Hthx_dot.l
$load Hgis_dot_param = Hgis_dot.l
$load Hais_dot_param = Hais_dot.l
$load Hais_smb_dot_param = Hais_smb_dot.l
$load Hgla_dot_param = Hgla_dot.l
$load Htot_dot_param = Htot_dot.l


$load MIU_param = MIU.l
$load EIND_param = EIND.l
$load C_param = C.l
$load KD_param = KD.l
$load I_param = I.l
$load Y_param = Y.l
$load CCA_param = CCA.l
$load UTILITY_param = UTILITY.l
$load METAUTILITY_param = METAUTILITY.l

$load a2D_var_param = a2D_var.l

$load MetaUTILITYAux_param = MetaUTILITYAux.l
$load MaxTempCon_param = MaxTempCon.l
$load MaxTemp_param = MaxTemp.l
$load MaxACon_param = MaxACon.l
$load MaxA_param = MaxA.l
$load MinpHCon_param = MinpHCon.l
$load MinpH_param = MinpH.l

$gdxin

** Affect result values to the model variables as a first guess
Tatm.l(m,t) = Tatm_param(m,t);
CO2.l(m,t) = CO2_param(m,t);

ECO2.l(m,t) = ECO2_param(m,t);
pH.l(m,t) = pH_param(m,t);
CO2_dot.l(m,t)= CO2_dot_param(m,t);

ERF.l(m,t) = ERF_param(m,t);
Td.l(m,t) = Td_param(m,t);
RFco2.l(m,t) = RFco2_param(m,t);
T_dot.l(m,t) = T_dot_param(m,t);
Td_dot.l(m,t) = Td_dot_param(m,t);

Fland.l(m,t) = Fland_param(m,t);
Cs1.l(m,t) = Cs1_param(m,t);
Cs2.l(m,t) = Cs2_param(m,t);
Cs3.l(m,t) = Cs3_param(m,t);
Cs.l(m,t) = Cs_param(m,t);
Cv.l(m,t) = Cv_param(m,t);
NPP.l(m,t) = NPP_param(m,t);
Efire.l(m,t) = Efire_param(m,t);
Eharv.l(m,t) = Eharv_param(m,t);
Fmort.l(m,t) = Fmort_param(m,t);
RH1.l(m,t) = RH1_param(m,t);
RH2.l(m,t) = RH2_param(m,t);
RH3.l(m,t) = RH3_param(m,t);
Fstab.l(m,t) = Fstab_param(m,t);
r_rh.l(m,t) = r_rh_param(m,t);
Fpass.l(m,t) = Fpass_param(m,t)  ;

Cv_dot.l(m,t) = Cv_dot_param(m,t);
Cs1_dot.l(m,t) = Cs1_dot_param(m,t);
Cs2_dot.l(m,t) = Cs2_dot_param(m,t);
Cs3_dot.l(m,t) = Cs3_dot_param(m,t);

Co_i.l(m,t, j) = Co_i_param(m,t, j);
Co.l(m,t) = Co_param(m,t);
Cd.l(m,t) = Cd_param(m,t);
pCO2.l(m,t) = pCO2_param(m,t);
Focean.l(m,t) = Focean_param(m,t);
Co_i_dot.l(m,t,j) = Co_i_dot_param(m,t,j);
Cd_dot.l(m,t) = Cd_dot_param(m,t);
dic.l(m,t) = dic_param(m,t);
pdic.l(m,t) = pdic_param(m,t);

a.l(m,t) = a_param(m,t);
abar.l(m,t) = abar_param(m,t);
r_rt.l(m,t) = r_rt_param(m,t);
Cth_i.l(m,t, pfl) = Cth_i_param(m,t, pfl);
Epf.l(m,t) = Epf_param(m,t);
expr.l(m,t) = expr_param(m,t);
a_dot.l(m,t) = a_dot_param(m,t);
Cfr_dot.l(m,t) = Cfr_dot_param(m,t);
Cth_i_dot.l(m,t, pfl) = Cth_i_dot_param(m,t, pfl);

Hthx.l(m,t) = Hthx_param(m,t);
Hgis.l(m,t) = Hgis_param(m,t);
Hais.l(m,t) = Hais_param(m,t);
Hais_smb.l(m,t) = Hais_smb_param(m,t);
Hgla.l(m,t) = Hgla_param(m,t);
Htot.l(m,t) = Htot_param(m,t);
OHC.l(m,t) = OHC_param(m,t);
OHC_dot.l(m,t) = OHC_dot_param(m,t);
Hthx_dot.l(m,t) = Hthx_dot_param(m,t);
Hgis_dot.l(m,t) = Hgis_dot_param(m,t);
Hais_dot.l(m,t) = Hais_dot_param(m,t);
Hais_smb_dot.l(m,t) = Hais_smb_dot_param(m,t);
Hgla_dot.l(m,t) = Hgla_dot_param(m,t);
Htot_dot.l(m,t) = Htot_dot_param(m,t);

MIU.l(m,t) = MIU_param(m,t);
EIND.l(m,t) = EIND_param(m,t);
C.l(m,t) = C_param(m,t);
KD.l(m,t) = KD_param(m,t);
I.l(m,t) = I_param(m,t);
Y.l(m,t) = Y_param(m,t);
CCA.l(m,t) = CCA_param(m,t);
UTILITY.l(m) = UTILITY_param(m);
METAUTILITY.l = METAUTILITY_param;

a2D_var.l = a2D_var_param;

MetaUTILITYAux.l = MetaUTILITYAux_param;
MaxTempCon.l = MaxTempCon_param;
MaxTemp.l(m) = MaxTemp_param(m);
MaxACon.l = MaxACon_param;
MaxA.l(m) = MaxA_param(m);
MinpHCon.l = MinpHCon_param;
MinpH.l(m) = MinpH_param(m);

execute_loadpoint "results2deg95.gdx"  ;

dt.fx=5;

miu.lo(m,t)=0;
miu.up(m,t)=miulim(t);
miu.fx(m,"1")=0.03;

MaxTempCon.up = inf;
MaxaCon.up    = inf;
MinpHCon.up   = inf;
targ_temp = inf;
targ_A    = inf;
targ_pH(m)= inf;
pen_temp = 0;
pen_A    = 0;
pen_pH   = 0;


a2D_var.fx = 0;



** Setting up model in one state only (i.e. without m dimension) for MonteCarlo
parameters config1, phi1, T2x1, THs1, THd1, th1, eheat1, aOHC1, Lthx1, lgla_01, Ggla_11, Ggla_31, ggla1, lgis_01, Lgis_11, Lgis_31, Lais_smb1, lais_01, aais1, k_toc1, vgx1, ggx1, To1, bdic1, gdic1, npp01, vfire1, vharv1, vmort1, vstab1, vrh11, vrh231, apass1, bnpp1, anpp1, gnpp1, bfire1, gfire1, brh1, grh1, ka1, ga1, k_tth1, Cfr01, CO2pi1, Lgla1, Lais1, tgla1, tgis1, tais1, adic1, aoc_11, aoc_21, aoc_31, aoc_41, aoc_51, toc_11, toc_21, toc_31, toc_41, toc_51, vrh31, aLST1, grt11, grt21, krt1, amin1, vthaw1, vfroz1, ath_11, ath_21, ath_31, tth_11, tth_21, tth_31, aCO21, k_pH1, T2x01 ;
parameters aoc1(j), toc1(j);
parameters ath1(pfl), tth1(pfl);
parameters ERFx1(t), etree1(t), sigma1(t);

variables Eco21(t), pH1(t), CO2_dot1(t), ERF1(t), Td1(t), RFco21(t), T_dot1(t), Td_dot1(t) ;
variables Fland1(t), Cs11(t), Cs21(t), Cs31(t), Cv1(t), NPP1(t), Efire1(t), Eharv1(t), Fmort1(t), RH11(t), RH21(t), RH31(t), Fstab1(t), r_rh1(t), Fpass1(t)  ;
variables Cv_dot1(t), Cs1_dot1(t), Cs2_dot1(t), Cs3_dot1(t);
variables Co_i1(t, j), Co1(t), Cd1(t), pCO21(t), Focean1(t), Co_i_dot1(t,j), Cd_dot1(t), dic1(t), pdic1(t);
variables a1(t), abar1(t), r_rt1(t),   Cth_i1(t, pfl), Epf1(t), expr1(t), a_dot1(t), Cfr_dot1(t), Cth_i_dot1(t, pfl);
variables OHC1(t), OHC_dot1(t), Hthx1(t), Hthx_dot1(t), Htot1(t), Htot_dot1(t), Hgis1(t), Hgis_dot1(t), Hais1(t), Hais_dot1(t), Hais_smb1(t), Hais_smb_dot1(t), Hgla1(t), Hgla_dot1(t)  ;

variables Tatm1(t), CO21(t) ;

equations Eco2_eq1(t), pH_eq1(t), CO2_tr1(t);
equations RFco2_eq1(t), ERF_eq1(t), T_dot_eq1(t), Td_dot_eq1(t), Tatm_tr1(t), Td_tr1(t);
equations NPP_eq1(t), r_rh_eq1(t), Efire_eq1(t), Eharv_eq1(t), Fmort_eq1(t), RH1_eq1(t), Fstab_eq1(t), RH2_eq1(t), Fpass_eq1(t), RH3_eq1(t), Fland_eq1(t), Cv_dot_eq1(t), Cs1_dot_eq1(t), Cs2_dot_eq1(t), Cs3_dot_eq1(t), Cv_tr1(t), Cs1_tr1(t), Cs2_tr1(t), Cs3_tr1(t);
equations Co_eq1(t), dic_eq1(t), pdic_eq1(t), pCO2_eq1(t), Focean_eq1(t),Co_ind_eq1(t, j), Co_i_tr1(t, j);
equations r_rt_eq1(t), abar_eq1(t), Epf_eq1(t) , expr_eq1(t), a_dot_eq1(t), Cth_i_dot_eq1(t, pfl), a_tr1(t), Cth_i_tr1(t,pfl);
equations OHC_eq1(t), OHC_dot_eq1(t), Hthx_eq1(t), Hthx_dot_eq1(t), Htot_eq1(t), Htot_dot_eq1(t), Hgis_dot_eq1(t), Hgis_tr1(t), Hais_dot_eq1(t), Hais_tr1(t), Hais_smb_dot_eq1(t), Hais_smb_tr1(t), Hgla_dot_eq1(t), Hgla_tr1(t);

** Atmospheric model
pH_eq1(t)..     pH1(t)         =E= k_pH1 * (8.5541 - 0.00173 * CO21(t) + 1.3264E-6 * CO21(t)**2 - 4.4943E-10 * CO21(t)**3) ;
Eco2_eq1(t)..   aCO21*CO2_dot1(t) =E= (Eco21(t) + Epf1(t) - Fland1(t) - Focean1(t) ) ;
CO2_tr1(t+1)..  CO21(t+1)      =E= CO21(t) + dt * CO2_dot1(t+1) ;

** Climate response model
RFco2_eq1(t)..  RFco21(t)      =E= phi1 * log(CO21(t)/CO2pi1) ;
ERF_eq1(t)..    ERF1(t)        =E= RFco21(t) + ERFx1(t) ;
T_dot_eq1(t)..  THs1 * T_dot1(t)     =E= ERF1(t) - phi1*log(2)*Tatm1(t)/T2x1 - eheat1*th1*(Tatm1(t)-Td1(t)) ;
Td_dot_eq1(t).. THd1 * Td_dot1(t)    =E= th1*(Tatm1(t)-Td1(t)) ;
Tatm_tr1(t+1).. Tatm1(t+1)     =E= Tatm1(t) + dt * T_dot1(t+1) ;
Td_tr1(t+1)..   Td1(t+1)       =E= Td1(t)   + dt * Td_dot1(t+1) ;

** Land model
NPP_eq1(t)..    NPP1(t)        =E= npp01*(1. + bnpp1 / anpp1 * (1. - (CO21(t) / CO2pi1+ 0.00001)**(-anpp1))) *(1. + gnpp1 * Tatm1(t)) ;
r_rh_eq1(t)..   r_rh1(t)       =E= (1. + brh1 * (Cs11(t) / (Cs11(t) + Cs21(t) + Cs31(t)+ 0.00001) * (1. + vstab1 / vrh231) - 1.)) * exp(grh1 * Tatm1(t)) ;
Efire_eq1(t)..  Efire1(t)      =E= vfire1 *Cv1(t) *(1.+ bfire1*(CO21(t)/CO2pi1 - 1.))*(1. + gfire1* Tatm1(t)) ;
Eharv_eq1(t)..  Eharv1(t)      =E= vharv1*Cv1(t) ;
Fmort_eq1(t)..  Fmort1(t)      =E= vmort1*Cv1(t) ;
RH1_eq1(t)..    RH11(t)        =E= vrh11 * r_rh1(t) * Cs11(t) ;
Fstab_eq1(t)..  Fstab1(t)      =E= vstab1 * r_rh1(t) * Cs11(t) ;
RH2_eq1(t)..    RH21(t)        =E= (vrh231 - vrh31 * apass1) / (1. - apass1) * r_rh1(t) * Cs21(t) ;
Fpass_eq1(t)..  Fpass1(t)      =E= vrh31 * apass1 / (1. - apass1) *  r_rh1(t) * Cs21(t) ;
RH3_eq1(t)..    RH31(t)        =E= vrh31 * r_rh1(t)  * Cs31(t) ;
Fland_eq1(t)..  Fland1(t)      =E= NPP1(t) - Efire1(t) - Eharv1(t) - RH11(t) - RH21(t) - RH31(t)  ;
Cv_dot_eq1(t).. Cv_dot1(t)     =E= NPP1(t) - Efire1(t) - Eharv1(t) - Fmort1(t) ;
Cs1_dot_eq1(t)..Cs1_dot1(t)    =E= Fmort1(t) - Fstab1(t)- RH11(t) ;
Cs2_dot_eq1(t)..Cs2_dot1(t)    =E= Fstab1(t) - Fpass1(t) - RH21(t) ;
Cs3_dot_eq1(t)..Cs3_dot1(t)    =E= Fpass1(t) - RH31(t) ;
Cv_tr1(t+1)..   Cv1(t+1)       =E= Cv1(t)   + dt * Cv_dot1(t+1) ;
Cs1_tr1(t+1)..  Cs11(t+1)      =E= Cs11(t)  + dt * Cs1_dot1(t+1) ;
Cs2_tr1(t+1)..  Cs21(t+1)      =E= Cs21(t)  + dt * Cs2_dot1(t+1) ;
Cs3_tr1(t+1)..  Cs31(t+1)      =E= Cs31(t)  + dt * Cs3_dot1(t+1) ;

** Ocean model
Co_eq1(t)..     Co1(t)         =E= sum(j, Co_i1(t, j) );
dic_eq1(t)..    dic1(t)        =E= adic1 / bdic1 * Co1(t) ;
pdic_eq1(t)..   pdic1(t)       =E= ((1.5568 - 1.3993E-2 * To1) * dic1(t)
                                    + (7.4706 - 0.20207   * To1) * 1E-3  * Power(dic1(t),2)
                                    - (1.2748 - 0.12015   * To1) * 1E-5  * Power(dic1(t),3)
                                    + (2.4491 - 0.12639   * To1) * 1E-7  * Power(dic1(t),4)
                                    - (1.5768 - 0.15326   * To1) * 1E-10 * Power(dic1(t),5) ) ;
pCO2_eq1(t)..   pCO21(t)       =E= (pdic1(t) + CO2pi1) * exp(gdic1 * Tatm1(t)) ;
Focean_eq1(t).. Focean1(t)     =E= vgx1 * (1. + ggx1 * Tatm1(t)) *(CO21(t) - pCO21(t)) ;
Co_ind_eq1(t,j)..Co_i_dot1(t,j)=E= - Co_i1(t, j) / (k_toc1*toc1(j)) + aoc1(j) * Focean1(t)  ;
Co_i_tr1(t+1,j)..Co_i1(t+1,j)   =E= Co_i1(t,j) + dt * Co_i_dot1(t+1,j) ;

** Permafrost model
r_rt_eq1(t)..   r_rt1(t)       =E= exp(krt1 * grt11 * aLST1 * Tatm1(t) - krt1 * grt21 * Power(aLST1 * Tatm1(t),2) );
expr_eq1(t)..   expr1(t)       =E= exp(- ka1 * ga1 *  aLST1 * Tatm1(t)) ;
abar_eq1(t)..   abar1(t)       =E= -amin1 + (1. + amin1) / (1. + ((1 + 1./amin1)**ka1 - 1.) * expr1(t))**(1./ka1) ;
Epf_eq1(t)..    Epf1(t)        =E= sum(pfl, Cth_i1(t, pfl)/tth1(pfl) ) * r_rt1(t)/k_tth1 ;
a_dot_eq1(t)..  a_dot1(t)      =E= (abar1(t) - a1(t)) * (vfroz1 + ( vthaw1 - vfroz1 ) * sigmoid(500 * (abar1(t) - a1(t)) ) );
Cth_i_dot_eq1(t, pfl)..
             Cth_i_dot1(t, pfl) =E= ath1(pfl) * a_dot1(t) * Cfr01 -  Cth_i1(t, pfl)/ tth1(pfl) * r_rt1(t)/ k_tth1 ;
a_tr1(t+1)..    a1(t+1)        =E= a1(t)   + dt * a_dot1(t+1) ;
Cth_i_tr1(t+1,pfl)..
              Cth_i1(t+1,pfl)   =E= Cth_i1(t,pfl) + dt * Cth_i_dot1(t+1,pfl) ;

** Sea-level rise and Ocean heat content model

OHC_eq1(t)..      OHC1(t)      =E= aOHC1 * (THs1 * Tatm1(t) + THd1 * Td1(t)) ;
OHC_dot_eq1(t)..  OHC_dot1(t)  =E= aOHC1 * (THs1 * T_dot1(t) + THd1 * Td_dot1(t)) ;

Hthx_eq1(t)..     Hthx1(t)     =E= Lthx1 * OHC1(t) ;
Hthx_dot_eq1(t).. Hthx_dot1(t) =E= Lthx1 * OHC_dot1(t) ;

Hgis_dot_eq1(t)..     Hgis_dot1(t)     =E= lgis_01 + 1/tgis1 * (Lgis_11*Tatm1(t) + Lgis_31* power(Tatm1(t),3) - Hgis1(t)) ;
Hais_smb_dot_eq1(t).. Hais_smb_dot1(t) =E= - Lais_smb1 * Tatm1(t) ;
Hais_dot_eq1(t)..     Hais_dot1(t)     =E= Hais_smb_dot1(t) + lais_01 + (Lais1 * Tatm1(t) - ( Hais1(t) - Hais_smb1(t) ) )/ tais1 * (1 + aais1 * (Hais1(t) - Hais_smb1(t)) ) ;
Hgla_dot_eq1(t)..     Hgla_dot1(t)     =E= lgla_01 + exp( ggla1 * Tatm1(t) )/tgla1 * ( Lgla1 * (1- exp( - Ggla_11*Tatm1(t) - Ggla_31*power(Tatm1(t),3) ) ) - Hgla1(t) ) ;


Hgis_tr1(t+1)..      Hgis1(t+1)     =E= Hgis1(t)    + dt * Hgis_dot1(t+1) ;
Hais_smb_tr1(t+1)..  Hais_smb1(t+1) =E= Hais_smb1(t)+ dt * Hais_smb_dot1(t+1) ;
Hais_tr1(t+1)..       Hais1(t+1)     =E= Hais1(t)    + dt * Hais_dot1(t+1) ;
Hgla_tr1(t+1)..       Hgla1(t+1)     =E= Hgla1(t)    + dt * Hgla_dot1(t+1) ;


Htot_dot_eq1(t).. Htot_dot1(t) =E= Hthx_dot1(t) + Hgis_dot1(t) + Hais_dot1(t) +Hgla_dot1(t) ;
Htot_eq1(t) ..    Htot1(t)     =E= Hthx1(t) + Hgis1(t) + Hais1(t) + Hgla1(t) ;

VARIABLES MIU1(t), EIND1(t), C1(t), KD1(t), I1(t), Y1(t), CCA1(t), UTILITY1, a2d_var1;
NONNEGATIVE VARIABLES  MIU1, Y1, C1, KD1, I1;
EQUATIONS EEQ1(t), EINDEQ1(t), CCACCA1(t), YY1(t), CC1(t), KK1(t), UTIL1;

eeq1(t)..             ECO21(t)        =E= EIND1(t) * 1/3.666 + etree1(t);
eindeq1(t)..          EIND1(t)        =E= sigma1(t) * (al(t)*(L(t)/1000)**(1-GAMA))*(KD1(t)**GAMA) * (1-(MIU1(t)));
ccacca1(t+1)..        CCA1(t+1)       =E= CCA1(t)+ EIND1(t)*5/3.666;
yy1(t)..              Y1(t)           =E= (al(t)*(L(t)/1000)**(1-GAMA))*(KD1(t)**GAMA) * (1- (a1D*TATM1(t) + a2D_var1*Power(TATM1(t),a3D) )   - (cost1(t) * (MIU1(t)**expcost2)) );
cc1(t)..              C1(t)           =E= Y1(t) - I1(t);
kk1(t+1)..            KD1(t+1)        =E= (1-dk)**tstep * KD1(t) + tstep * I1(t);
util1..               UTILITY1        =E=  SUM(T,  1/( (1+prstp)**( tstep * (ord(T)-1) ) ) * L(T)* log(  tstep * C1(t)*1000/L(T) ) ) ;

model IIASAsingle / Eco2_eq1, pH_eq1, CO2_tr1, RFco2_eq1, ERF_eq1, T_dot_eq1, Td_dot_eq1, Tatm_tr1, Td_tr1, NPP_eq1, r_rh_eq1, Efire_eq1, Eharv_eq1, Fmort_eq1,  RH1_eq1, Fstab_eq1, RH2_eq1, Fpass_eq1, RH3_eq1, Fland_eq1, Cv_dot_eq1, Cs1_dot_eq1, Cs2_dot_eq1, Cs3_dot_eq1, Cv_tr1, Cs1_tr1, Cs2_tr1, Cs3_tr1, Co_eq1, dic_eq1, pdic_eq1, pCO2_eq1, Focean_eq1,Co_ind_eq1, Co_i_tr1, r_rt_eq1, abar_eq1, Epf_eq1, expr_eq1, a_dot_eq1, Cth_i_dot_eq1, a_tr1, Cth_i_tr1 , OHC_eq1, OHC_dot_eq1, Hthx_eq1, Hthx_dot_eq1, Htot_eq1, Htot_dot_eq1, Hgis_dot_eq1, Hgis_tr1, Hais_dot_eq1, Hais_tr1, Hais_smb_dot_eq1, Hais_smb_tr1, Hgla_dot_eq1, Hgla_tr1 / ;
model DICEsingle  / EEQ1, EINDEQ1, CCACCA1, YY1, CC1, KK1, UTIL1 / ;
model Nizza       / IIASAsingle, DICEsingle / ;

equation obj_eq1; variable obj1;
obj_eq1.. obj1 =E= dt + 0*tatm1("1") ;
model Init / IIASAsingle + DICEsingle + obj_eq1/ ;

equation Miu_Ramp11(t), Miu_Ramp12(t);
Miu_Ramp11(t+1)$(t.val LT tflip+1).. EInd1(t) - EInd1(t-1)=G= -10;
Miu_Ramp12(t+1)$(t.val GE tflip+1).. miu1(t)              =L= miu1(t-1) * 1.1 ;
model Ramp1 / Miu_Ramp11, Miu_Ramp12 / ;
model NizzaRamp / Nizza, Ramp1 / ;

equation deg_eq1, utility_2deg_eq1; variable lim2deg, Utility_2deg1; scalar pen2deg, tar2deg;
deg_eq1(t)..          Tatm1(t)        =L= lim2deg ;
utility_2deg_eq1..    Utility_2deg1   =E= Utility1 - pen2deg* Power(lim2deg - tar2deg, 2);
model Nizza2deg / Nizza + deg_eq1 + utility_2deg_eq1/ ;
model Nizza2degRamp / Nizza2deg, Ramp1 / ;

equation A_eq1, utility_A_eq1; variable limA, Utility_A1; scalar penA1, tarA1;
A_eq1(t)..            a1(t)        =L= limA ;
utility_A_eq1..       Utility_A1   =E= Utility1 - penA1* Power(limA - tarA1, 2);
model NizzaA / Nizza + A_eq1 + utility_A_eq1/ ;
model NizzaAramp / NizzaA, Ramp1 / ;


** bounds (state-independent)
CO21.lo(t) = 0.5 * CO2pi("1");
Tatm1.lo(t) = -10;
Tatm1.up(t) = 20;
Td1.lo(t) = -inf;

dt.fx = 5;

miu1.lo(t) = 0;
miu1.up(t) = miulim(t);
miu1.fx("1")= miu0;

c1.lo(t) = 2;
c1.up(t) = inf;

cca1.lo(t) = -inf;
cca1.up(t) = fosslim;

parameter scc1, cprice1, optsol;

miu1.lo(t) = 0;
miu1.up(t) = miulim(t);
miu1.fx("1")= miu0;

c1.lo(t) = 2;
c1.up(t) = inf;

cca1.lo(t) = -inf;
cca1.up(t) = fosslim;

a2D_var1.lo = a2d;

CO21.lo(t) = 0.5 * CO2pi("1");
Tatm1.lo(t) = -10; Tatm1.up(t) = inf;
Td1.lo(t) = -inf;


** setting penalty paramers for thawing of permafrost above 50%
penA1=10**9; tarA1 = 0.5 ;

a2d_var1.fx = 0;

** Looping for the Monte Carlo simulation
loop(m,
phi1 = phi(m) ; T2x1 = T2x(m) ; THs1 = THs(m) ; THd1 = THd(m) ; th1 = th(m) ; eheat1 = eheat(m) ; aOHC1 = aOHC(m) ; Lthx1 = Lthx(m) ; lgla_01 = lgla_0(m) ; Ggla_11 = Ggla_1(m) ; Ggla_31 = Ggla_3(m) ; ggla1 = ggla(m) ; lgis_01 = lgis_0(m) ; Lgis_11 = Lgis_1(m) ; Lgis_31 = Lgis_3(m) ; Lais_smb1 = Lais_smb(m) ; lais_01 = lais_0(m) ; aais1 = aais(m) ; k_toc1 = k_toc(m) ; vgx1 = vgx(m) ; ggx1 = ggx(m) ; To1 = To(m) ; bdic1 = bdic(m) ; gdic1 = gdic(m) ; npp01 = npp0(m) ; vfire1 = vfire(m) ; vharv1 = vharv(m) ; vmort1 = vmort(m) ; vstab1 = vstab(m) ; vrh11 = vrh1(m) ; vrh231 = vrh23(m) ; apass1 = apass(m) ; bnpp1 = bnpp(m) ; anpp1 = anpp(m) ; gnpp1 = gnpp(m) ; bfire1 = bfire(m) ; gfire1 = gfire(m) ; brh1 = brh(m) ; grh1 = grh(m) ; ka1 = ka(m) ; ga1 = ga(m) ; k_tth1 = k_tth(m) ; Cfr01 = Cfr0(m) ; CO2pi1 = CO2pi(m) ; Lgla1 = Lgla(m) ; Lais1 = Lais(m) ;  tgla1 = tgla(m) ; tgis1 = tgis(m) ; tais1 = tais(m) ; adic1 = adic(m) ; aoc_11 = aoc_1(m) ; aoc_21 = aoc_2(m) ; aoc_31 = aoc_3(m) ; aoc_41 = aoc_4(m) ; aoc_51 = aoc_5(m) ; toc_11 = toc_1(m) ; toc_21 = toc_2(m) ; toc_31 = toc_3(m) ; toc_41 = toc_4(m) ; toc_51 = toc_5(m) ; vrh31 = vrh3(m) ; aLST1 = aLST(m) ; grt11 = grt1(m) ; grt21 = grt2(m) ; krt1 = krt(m) ; amin1 = amin(m) ; vthaw1 = vthaw(m) ; vfroz1 = vfroz(m) ; ath_11 = ath_1(m) ; ath_21 = ath_2(m) ; ath_31 = ath_3(m) ; tth_11 = tth_1(m) ; tth_21 = tth_2(m) ; tth_31 = tth_3(m) ; aCO21 = aCO2(m) ; k_pH1 = k_pH(m) ; T2x01 = T2x0(m) ;

ERFx1(t) = ERFx(m,t); etree1(t) = Eluc(m,t); sigma1(t) = sigma(m,t);

aoc1(j) = aoc(m,j); toc1(j) = toc(m,j);
ath1(pfl)= ath(m,pfl); tth1(pfl) = tth(m,pfl);

** lower bounds (state-dependent)
Efire1.lo(t) = 0;
Cv1.lo(t) = 0.01*Cv_0(m);
Cs11.lo(t) = 0.01* Cs1_0(m);
Cs21.lo(t) = 0.01* Cs2_0(m);
Cs31.lo(t) = 0.01* Cs3_0(m);
a1.lo(t) = -amin1 ; a1.up(t) = 1;
abar1.lo(t) = -amin1 ; abar1.up(t) = 1;

KD1.lo(t)= 1;

** initial conditions
CCA1.FX("1")    = 400;
KD1.FX("1")     = kD0;

Tatm1.fx("1") = T_0(m);
CO21.fx("1") = CO2_0(m);

Co_i1.fx("1","1") = Co_1_0(m);
Co_i1.fx("1","2") = Co_2_0(m);
Co_i1.fx("1","3") = Co_3_0(m);
Co_i1.fx("1","4") = Co_4_0(m);
Co_i1.fx("1","5") = Co_5_0(m);

Cv1.fx("1") = Cv_0(m);
Cd1.fx("1") = Cd_0(m);

Cs11.fx("1") = Cs1_0(m);
Cs21.fx("1") = Cs2_0(m);
Cs31.fx("1") = Cs3_0(m);

a1.fx("1") = a_0(m);
Cth_i1.fx("1","1") = Cth_1_0(m);
Cth_i1.fx("1","2") = Cth_2_0(m);
Cth_i1.fx("1","3") = Cth_3_0(m);

Td1.fx("1") = Td_0(m);

Hgla1.fx("1")     = Hgla_0(m);
Hgis1.fx("1")     = Hgis_0(m);
Hais1.fx("1")     = Hais_0(m);
Hais_smb1.fx("1") = Hais_smb_0(m);

If((ifmod eq 1), solve NizzaARamp maxmizing utility_A1 using nlp;
else             solve NizzaA maxmizing utility_A1 using nlp;
);

scc1(t)   = -1000*eeq1.m(t)/(cc1.m(t)+0.00001);
CPRICE1(t)= 3.666* 1000/sigma1(t) * cost1(t) * expcost2 * (MIU1.l(t))**(expcost2-1);
display scc1, cprice1;

Eco2.l(m,t) = Eco21.l(t) ; pH.l(m,t) = pH1.l(t) ; CO2_dot.l(m,t) = CO2_dot1.l(t) ; ERF.l(m,t) = ERF1.l(t) ; Td.l(m,t) = Td1.l(t) ; RFco2.l(m,t) = RFco21.l(t) ; T_dot.l(m,t) = T_dot1.l(t) ; Td_dot.l(m,t) = Td_dot1.l(t) ; Fland.l(m,t) = Fland1.l(t) ; Cs1.l(m,t) = Cs11.l(t) ; Cs2.l(m,t) = Cs21.l(t) ; Cs3.l(m,t) = Cs31.l(t) ; Cv.l(m,t) = Cv1.l(t) ; NPP.l(m,t) = NPP1.l(t) ; Efire.l(m,t) = Efire1.l(t) ; Eharv.l(m,t) = Eharv1.l(t) ; Fmort.l(m,t) = Fmort1.l(t) ; RH1.l(m,t) = RH11.l(t) ; RH2.l(m,t) = RH21.l(t) ; RH3.l(m,t) = RH31.l(t) ; Fstab.l(m,t) = Fstab1.l(t) ; r_rh.l(m,t) = r_rh1.l(t) ; Fpass.l(m,t) = Fpass1.l(t) ; Cv_dot.l(m,t) = Cv_dot1.l(t) ; Cs1_dot.l(m,t) = Cs1_dot1.l(t) ; Cs2_dot.l(m,t) = Cs2_dot1.l(t) ; Cs3_dot.l(m,t) = Cs3_dot1.l(t) ; Co_i.l(m,t,j) = Co_i1.l(t,j) ; Co.l(m,t) = Co1.l(t) ; Cd.l(m,t) = Cd1.l(t) ; pCO2.l(m,t) = pCO21.l(t) ; Focean.l(m,t) = Focean1.l(t) ; Co_i_dot.l(m,t,j) = Co_i_dot1.l(t,j) ; dic.l(m,t) = dic1.l(t) ; pdic.l(m,t) = pdic1.l(t) ; a.l(m,t) = a1.l(t) ; abar.l(m,t) = abar1.l(t) ; r_rt.l(m,t) = r_rt1.l(t) ; Cth_i.l(m,t,pfl) = Cth_i1.l(t,pfl) ; Epf.l(m,t) = Epf1.l(t) ; expr.l(m,t) = expr1.l(t) ; a_dot.l(m,t) = a_dot1.l(t) ; Cth_i_dot.l(m,t,pfl) = Cth_i_dot1.l(t,pfl) ;
OHC.l(m,t) = OHC1.l(t) ; OHC_dot.l(m,t) = OHC_dot1.l(t) ;
Hthx.l(m,t) = Hthx1.l(t) ; Hthx_dot.l(m,t) = Hthx_dot1.l(t) ; Htot.l(m,t) = Htot1.l(t) ; Htot_dot.l(m,t) = Htot_dot1.l(t) ; Hgis.l(m,t) = Hgis1.l(t) ; Hgis_dot.l(m,t) = Hgis_dot1.l(t) ; Hais.l(m,t) = Hais1.l(t) ; Hais_dot.l(m,t) = Hais_dot1.l(t) ; Hais_smb.l(m,t) = Hais_smb1.l(t) ; Hais_smb_dot.l(m,t) = Hais_smb_dot1.l(t) ; Hgla.l(m,t) = Hgla1.l(t) ; Hgla_dot.l(m,t) = Hgla_dot1.l(t)  ;
Tatm.l(m,t) = Tatm1.l(t) ; CO2.l(m,t) = CO21.l(t) ;

Miu.l(m,t) = MIU1.l(t) ; EIND.l(m,t) = EIND1.l(t) ; C.l(m,t) = C1.l(t); KD.l(m,t) = KD1.l(t); I.l(m,t) = I1.l(t); Y.l(m,t) = Y1.l(t); CCA.l(m,t) = CCA1.l(t); Utility.l(m) = UTILITY1.l;
Cprice(m,t) = Cprice1(t);
optsol(m) = NizzaARamp.modelstat$(ifmod eq 1) + nizzaA.modelstat$(ifmod eq 0) ;


display phi1, T2x1, THs1, THd1, th1, eheat1, aOHC1, Lthx1, lgla_01, Ggla_11, Ggla_31, ggla1, lgis_01, Lgis_11, Lgis_31, Lais_smb1, lais_01, aais1, k_toc1, vgx1, ggx1, To1, bdic1, gdic1, npp01, vfire1, vharv1, vmort1, vstab1, vrh11, vrh231, apass1, bnpp1, anpp1, gnpp1, bfire1, gfire1, brh1, grh1, ka1, ga1, k_tth1, Cfr01, CO2pi1, Lgla1, Lais1, tgla1, tgis1, tais1, adic1, aoc_11, aoc_21, aoc_31, aoc_41, aoc_51, toc_11, toc_21, toc_31, toc_41, toc_51, vrh31, aLST1, grt11, grt21, krt1, amin1, vthaw1, vfroz1, ath_11, ath_21, ath_31, tth_11, tth_21, tth_31, aCO21, k_pH1, T2x01 ;
display aoc1, toc1;
display ath1, tth1;

);

File NiceSA05 /NiceSA05.txt/;
NiceSA05.pc=5;
NiceSA05.pw=32767;
Put NiceSA05;
Put / "Time"; Loop (t, put t.val::0);
Put / "CPRICE(t) = 3.70 *  1000/sigma1(t) * cost1(t) * expcost2 * (MIU(t))**(expcost2-1)"; Put "expcost2"; Put expcost2::5;
Put / "pbacktime(t)"; Loop(t, put pbacktime(t):: 5 ) ;
Put / "Monte Carlo"; Put "Damage Parameter"; Put a2d_var1.l::5;
Put / "Mitigation";
Loop(m, Put / m.val; Loop (t, put miu.l(m,t)::5 ); );
Put / "Cprice";
Loop(m, Put / m.val; Loop (t, put Cprice(m,t)::5 ); );
Put / "Temperature";
Loop(m, Put / m.val; Loop (t, put tatm.l(m,t)::5 ); );
Put / "Consumption";
Loop(m, Put / m.val; Loop (t, put C.l(m,t)::5 ); );
Put / / "Solution Status > 2  --> not optimal";
Loop(m, Put / m.val; Put optsol(m)::0 );
putclose;

pen2deg = 0; tar2deg = 0 ;


a2d_var.fx = 0;
miu.fx(m,t) = miu.l(m,t);
c.fx(m,t) = c.l(m,t);

MetaUtility.l =  ( Sum(m, (utility.l(m)/ scale1)**(1-chi) )  )**(1/(1-chi)) ;
C.lo(m,t)$(t.val EQ Card(t)) =  0.5 * c.l(m,t); C.up(m,t)$(t.val EQ Card(t)) = inf;

If((ifmod eq 1), solve MonteCarloRamp maximize MetaUTILITY using NLP;
else             solve MonteCarlo maximize MetaUTILITY using NLP;
);

parameter Amax, Amaxmax;
Amax(m) = smax(t, A.l(m,t));
Amaxmax = smax(m, Amax(m));
display Amax, Amaxmax;

A.up(m,t) = max(0.5, Amax(m));
A.fx(m,"1") = A_0(m);

c.lo(m,t) = 0.5 * c.l(m,T);
c.up(m,t) = inf;

miu.lo(m,t) = 0;
miu.up(m,t) = miulim(t);
miu.fx(m,"1")= miu0;


If((ifmod eq 1), solve MonteCarloRamp maximize MetaUTILITY using NLP; solstat = MonteCarloRamp.modelstat ;
else             solve MonteCarlo maximize MetaUTILITY using NLP; solstat = MonteCarlo.modelstat ;
);
If((ifmod eq 1), solve MonteCarloRamp maximize MetaUTILITY using NLP;   solstat = MonteCarloRamp.modelstat ;
else             solve MonteCarlo maximize MetaUTILITY using NLP;     solstat = MonteCarlo.modelstat ;
);
If((ifmod eq 1), solve MonteCarloRamp maximize MetaUTILITY using NLP;   solstat = MonteCarloRamp.modelstat ;
else             solve MonteCarlo maximize MetaUTILITY using NLP;     solstat = MonteCarlo.modelstat ;
);
If((ifmod eq 1), solve MonteCarloRamp maximize MetaUTILITY using NLP; solstat = MonteCarloRamp.modelstat ;
else             solve MonteCarlo maximize MetaUTILITY using NLP; solstat = MonteCarlo.modelstat ;
);
If((ifmod eq 1), solve MonteCarloRamp maximize MetaUTILITY using NLP;   solstat = MonteCarloRamp.modelstat ;
else             solve MonteCarlo maximize MetaUTILITY using NLP;     solstat = MonteCarlo.modelstat ;
);
If((ifmod eq 1), solve MonteCarloRamp maximize MetaUTILITY using NLP;   solstat = MonteCarloRamp.modelstat ;
else             solve MonteCarlo maximize MetaUTILITY using NLP;     solstat = MonteCarlo.modelstat ;
);


parameter lambdaC, lambdaUniC, lambdaE, lambdaUniE, focunic, focmiu1, focmiu2, focunimiu, foclambda1, foclambda, mutil;

cprice(m,t) = 3.666* 1000/sigma(m,t) * cost1(t) * expcost2 * (MIU.l(m,t))**(expcost2-1);
scc(m,t)    = -1000*eindeq.m(m,t) / (cc.m(m,t) + 10**(-15) )  ;
mutil(m,t)  = Card(m) * ( (tstep * C.l(m,t)*1000 /L(T)) ) ** ( - elasmu ) ;

lambdaC(m,t) = 1/5 * CC.m(m,t)     * 1/ (1+prstp)**( - tstep * (t.val-1) ) ;
lambdaE(m,t) = 1/5 * EIndeq.m(m,t) * 1/ (1+prstp)**( - tstep * (t.val-1) ) ;
focunic(t) = sum(m, lambdac(m,t) ) - 1/Card(m) * sum(m, 1000* ( (tstep * C.l(m,t)*1000 /L(T)) ) ** ( - elasmu ) ) ;

focmiu1(m,t) = tstep * lambdac(m,t) * ( cost1(t) * expcost2 * miu.l(m,t)**(expcost2-1) ) * (al(t)*(L(t)/1000)**(1-GAMA))*(KD.l(m,t)**GAMA) ;
focmiu2(m,t) = sigma(m,t) * 5/1000 * scc(m,t) * lambdac(m,t) * (al(t)*(L(t)/1000)**(1-GAMA))*(KD.l(m,t)**GAMA) ;
focunimiu(t) = sum(m, focmiu1(m,t)) - sum(m, focmiu2(m,t));

foclambda1(m,t)$(t.val LT Card(T)) =  1/(1+prstp)**   tstep   * ( (1-DK)**tstep + tstep * GAMA * 1/ KD.l(m,T+1) * (1 - ( a1D*TATM.l(m,t+1) + a2D*Power(TATM.l(m,t+1),a3D) ) - (  cost1(T+1) * MIU.l(m,T+1)**expcost2   ) - (1 - MIU.l(m,T+1)) * (  expcost2 * cost1(t+1) * MIU.l(m,T+1)**(expcost2-1) )  ) * (al(t+1)*(L(t+1)/1000)**(1-GAMA))*(KD.l(m,T+1)**GAMA)  ) * lambdac(m,t+1);
foclambda(m,t) = foclambda1(m,t) - lambdac(m,t) ;

execute_unload "resultsMC05A.gdx"
execute 'gdxxrw.exe resultsMC05A.gdx o=resultsMC05A.xlsx text="Damage Parameter" rng=stats!a3 var=a2d_var.l rng=stats!b3 text="Solution status" rng=stats!a4 par=solstat rng=stats!b4 text="modified DICE (yes if 1)" rng=stats!a5 par=ifmod rng=stats!b5 text="population" rng=stats!a6 par=l rng=stats!b6 ' ;
execute 'gdxxrw.exe resultsMC05A.gdx o=resultsMC05A.xlsx var=Eind.l rng=Eind!a1 text="Eind" rng=Eind!A1 par=Etree rng=Etree!a1 text="Etree" rng=Etree!A1 var=ECO2.l rng=ECO2!a1 text="ECO2" rng=ECO2!A1 var=Fland.l rng=Fland!a1 text="Fland" rng=Fland!A1  var=Focean.l rng=Focean!a1 text="Focean" rng=Focean!A1 var=Epf.l rng=Epf!a1 text="Epf" rng=Epf!A1 var=Co2.l rng=Co2!a1 text="Co2" rng=Co2!A1 var=Co2_dot.l rng=Co2_dot!a1 text="Co2_dot" rng=Co2_dot!A1';
execute 'gdxxrw.exe resultsMC05A.gdx o=resultsMC05A.xlsx var=RFCO2.l rng=RFCO2!a1 text="RFCO2" rng=RFCO2!A1 par=ERFx rng=ERFx!a1 text="ERFx" rng=ERFx!A1 var=ERF.l rng=ERF!a1 text="ERF" rng=ERF!A1';
execute 'gdxxrw.exe resultsMC05A.gdx o=resultsMC05A.xlsx var=miu.l rng=miu!a1 text="miu" rng=miu!A1  var=C.l rng=C!a1 text="C" rng=C!A1  var=Y.l rng=Y!a1 text="Y" rng=Y!A1  var=Kd.l rng=Kd!a1 text="Kd" rng=Kd!A1 var=Y.l rng=Y!a1 text="Y" rng=Y!A1';
execute 'gdxxrw.exe resultsMC05A.gdx o=resultsMC05A.xlsx var=Tatm.l rng=Tatm!a1 text="Tatm" rng=Tatm!A1  var=pH.l rng=pH!a1 text="pH" rng=pH!A1 var=a.l rng=a!a1 text="a" rng=a!A1 var=Td.l rng=Td!a1 text="Td" rng=Td!A1 var=T_dot.l rng=T_dot!a1 text="T_dot" rng=T_dot!A1 var=Td_dot.l rng=Td_dot!a1 text="Td_dot" rng=Td_dot!A1' ;
execute 'gdxxrw.exe resultsMC05A.gdx o=resultsMC05A.xlsx var=OHC.l rng=OHC!a1 text="OHC" rng=OHC!A1 var=OHC_dot.l rng=OHC_dot!a1 text="OHC_dot" rng=OHC_dot!A1 var=Hthx_dot.l rng=Hthx_dot!a1 text="Hthx_dot" rng=Hthx_dot!A1 var=Hgla_dot.l rng=Hgla_dot!a1 text="Hgla_dot" rng=Hgla_dot!A1 var=Hgis_dot.l rng=Hgis_dot!a1 text="Hgis_dot" rng=Hgis_dot!A1 var=Hais_dot.l rng=Hais_dot!a1 text="Hais_dot" rng=Hais_dot!A1 var=Hais_smb_dot.l rng=Hais_smb_dot!a1 text="Hais_smb_dot" rng=Hais_smb_dot!A1 var=Htot_dot.l rng=Htot_dot!a1 text="Htot_dot" rng=Htot_dot!A1 var=Hthx.l rng=Hthx!a1 text="Hthx" rng=Hthx!A1 var=Hgla.l rng=Hgla!a1 text="Hgla" rng=Hgla!A1 var=Hgis.l rng=Hgis!a1 text="Hgis" rng=Hgis!A1 var=Hais.l rng=Hais!a1 text="Hais" rng=Hais!A1 var=Hais_smb.l rng=Hais_smb!a1 text="Hais_smb" rng=Hais_smb!A1 var=Htot.l rng=Htot!a1 text="Htot" rng=Htot!A1' ;
execute 'gdxxrw.exe resultsMC05A.gdx o=resultsMC05A.xlsx par=lambdaC rng=lambdaC!a1 text="lambdaC" rng=lambdaC!A1  par=lambdaE rng=lambdaE!a1 text="lambdaE" rng=lambdaE!A1  par=Cprice rng=Cprice!a1 text="Cprice" rng=Cprice!A1  par=SCC rng=SCC!a1 text="false SCC" rng=SCC!A1  par=FOCuniC rng=FOCuniC!a1 text="FOCuniC" rng=FOCuniC!A1  par=FOCuniMIU rng=FOCuniMIU!a1 text="FOCuniMIU" rng=FOCuniMIU!A1  par=FOClambda rng=FOClambda!a1 text="FOClambda" rng=FOClambda!A1 par=mutil rng=marg_util!a1 text="marginal utility" rng=marg_util!A1' ;
execute 'gdxxrw.exe resultsMC05A.gdx o=resultsMC05A.xlsx var=utility.l rng=utility!a1 text="Utility" rng=utility!A1  var=Metautility.l rng=Metautility!a2 text="Metautility" rng=Metautility!A1' ;


$exit

