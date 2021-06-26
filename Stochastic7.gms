$ontext
This is the beta version of DICE-2016R. The major changes are outlined in Nordhaus,
"Revisiting the social cost of carbon: Estimates from the DICE-2016R model,"
September 30, 2016," available from the author.

Version is DICE-2016R-091916ap.gms
$offtext

$title        DICE-2016R September 2016 (DICE-2016R-091216a.gms)

$set timehorizon 50

$set timesteplength 5

sets        t  Time periods (5 years per period)                                            /1*%timehorizon%/
            st1  Time step when the first catastrophy can occour -> changes tcre            /0,10,20,30,40/
            st2  Time step when the second catastrophy can occour -> changes emissions      /0,10,20,30,40/
;


* This set is one of the most important. We have a tuple with (current time, time of catastrophy 1, time of catastrophy 2)
* The current time is forced to be greater than the time of the catastrophy, so that we don't know that a catastrophy will occour
* until it happens. In this way we describe in an unique way the nodes of a tree where there is a branching at each time step when
* a catastrophy can occour. If the catastrophy occours its effects will remain forever.
sets as(t,st1,st2) active states
;

* This is the part where we define the as set that will used in the equations. Check the .gdx to better understand its behaviour 

parameter os(t,st1,st2) off-set pointer
;

os(t,st1,st2)$(t.val < st1.val or t.val < st2.val) = min(-st1.val,-st2.val) ;


as(t,st1,st2) = yes$(os(t,st1,st2)=0);

*-----------------------------------------------------------------------------------




PARAMETERS
** Availability of fossil fuels
        fosslim  Maximum cumulative extraction fossil fuels (GtC)  /6000/
**Time Step
        tstep    Years per Period                                  /%timesteplength%/
** If optimal control
        ifopt    Indicator where optimized is 1 and base is 0      /1/
** montecarlo control
        montecarlo Indicator where montecarlo is 1 and base is 0   /0/

** Preferences
        elasmu   Elasticity of marginal utility of consumption     /1.45 /
        prstp    Initial rate of social time preference per year   /.015  /
** Population and technology
        gama     Capital elasticity in production function        /.300    /
        pop0     Initial world population 2015 (millions)         /7403    /
        popadj   Growth rate to calibrate to 2050 pop projection  /0.134   /
        popasym  Asymptotic population (millions)                 /11500   /
        dk       Depreciation rate on capital (per year)          /.100    /
        q0       Initial world gross output 2015 (trill 2010 USD) /105.5   /
        k0       Initial capital value 2015 (trill 2010 USD)      /223     /
        a0       Initial level of total factor productivity       /5.115   /
        ga0      Initial growth rate for TFP per 5 years          /0.076   /
        dela     Decline rate of TFP per 5 years                  /0.005   /
** Emissions parameters
        gsigma1  Initial growth of sigma (per year)                   /-0.0152 /
        dsig     Decline rate of decarbonization (per period)         /-0.001  /
        eland0   Carbon emissions from land 2015 (GtCO2 per year)     / 2.6    /
        deland   Decline rate of land emissions (per period)          / .115   /
        e0       Industrial emissions 2015 (GtCO2 per year)           /35.85    /
        miu0     Initial emissions control rate for base case 2015    /.03     /
        sig0     Carbon intensity 2010 (kgCO2 per output 2005 USD 2010)
        
** Climate model parameters
* original
        tcre     transient climate response (oC per TtC:range .8-2.4)  / 1.6  /        
        fex0     2015 forcings of non-CO2 GHG (Wm-2)              / 0.5  /
        fex1     2100 forcings of non-CO2 GHG (Wm-2)              / 1.0  /
        tatm0    Initial atmospheric temp change (C from 1900)    /0.85  /
* effects of tipping points 
        tcre_cat catastrophy has higher climate response              /2.4/
        eind_cat Coefficient of industrial emission in catastrophy    /1.27/
        
** Climate damage parameters
* original, are not used anymore
        a10       Initial damage intercept                         /0       /
        a20       Initial damage quadratic term
        a1        Damage intercept                                 /0       /
        a2        Damage quadratic term                            /0.00236 /
        a3        Damage exponent                                  /2.00    /
* new damage function (NEW)
        w1        Damage Weitzman first coefficient                /20.64/     
        w2        Damage Weitzman second coefficient               /6.081/
        w3        Damage Weitzman third coefficient                /6.754/
        
        c1_a      Damage increment in tipping 1                    /1.2/
        c1_b      Damage increment in tipping 2                    /1.2/  
        c2        Damage increment in both tipping                 /1.5/
                        
        
** Abatement cost
        expcost2  Exponent of control cost function               / 2.6  /
        pback     Cost of backstop 2010$ per tCO2 2015            / 550  /
        gback     Initial cost decline backstop cost per period   / .025 /
        limmiu    Upper limit on control rate after 2150          / 1.2 /
        tnopol    Period before which no emissions controls base  / 45   /
        cprice0   Initial base carbon price (2010$ per tCO2)      / 0    /
        gcprice   Growth rate of base carbon price per year       /.02   /
        
** Probability parameters 
        beta      Coefficient of tanh for hazards rates             /0.5/
        tref      Traslation of tanh for hazards rates             /4/
        

** Scaling and inessential parameters
* Note that these are unnecessary for the calculations
* They ensure that MU of first period's consumption =1 and PV cons = PV utilty
        scale1      Multiplicative scaling coefficient           /0.0302455265681763 /
        scale2      Additive scaling coefficient                 /-10993.704/
        
;

* Program control variables
sets     tfirst(t), tlast(t), tearly(t), tlate(t);

PARAMETERS
        l(t)          Level of population and labor
        al(t)         Level of total factor productivity
        sigma(t)      CO2-equivalent-emissions output ratio
        rr(t)         Average utility social discount rate
        ga(t)         Growth rate of productivity from
        forcoth(t)    Exogenous forcing for other greenhouse gases
        gl(t)         Growth rate of labor
        gcost1        Growth of cost factor
        gsig(t)       Change in sigma (cumulative improvement of energy efficiency)
        etree(t)      Emissions from deforestation
        cumetree(t)   Cumulative from land
        cost1(t)      Adjusted cost for backstop
        lam           Climate model parameter
        gfacpop(t)    Growth factor population
        pbacktime(t)  Backstop price
        optlrsav      Optimal long-run savings rate used for transversality
        scc(t)        Social cost of carbon
        cpricebase(t) Carbon price in base case
        photel(t)     Carbon Price under no damages (Hotelling rent condition)
;
* Program control definitions
        tfirst(t) = yes$(t.val eq 1);
        tlast(t)  = yes$(t.val eq card(t));
* Further definitions of parameters
        a20 = a2;
        sig0 = e0/(q0*(1-miu0));
        l("1") = pop0;
        loop(t, l(t+1)=l(t););
        loop(t, l(t+1)=l(t)*(popasym/L(t))**popadj ;);
    display l;
        ga(t)=ga0*exp(-dela*5*((t.val-1)));
        al("1") = a0; loop(t, al(t+1)=al(t)/((1-ga(t))););
        gsig("1")=gsigma1; loop(t,gsig(t+1)=gsig(t)*((1+dsig)**tstep) ;);
        sigma("1")=sig0;   loop(t,sigma(t+1)=(sigma(t)*exp(gsig(t)*tstep)););

        pbacktime(t)=pback*(1-gback)**(t.val-1);
        cost1(t) = pbacktime(t)*sigma(t)/expcost2/1000;

        etree(t) = eland0*(1-deland)**(t.val-1);
        cumetree("1")= 100; loop(t,cumetree(t+1)=cumetree(t)+etree(t)*(5/3.666););

        rr(t) = 1/((1+prstp)**(tstep*(t.val-1)));
        forcoth(t) = fex0+ (1/17)*(fex1-fex0)*(t.val-1)$(t.val lt 18)+ (fex1-fex0)$(t.val ge 18);
        optlrsav = (dk + .004)/(dk + .004*elasmu + prstp)*gama;

*Base Case Carbon Price
        cpricebase(t)= cprice0*(1+gcprice)**(tstep*(t.val-1));

*
* Variables are now defined over the triple (t,st1,st2) but in reality only the ones over the set as(t,st1,st2) are meaningful
* we are forced to define them over a bigger (static) set since GAMS doesn't allow definitions over dynamic sets like as(t,st1,st2)
* we will fix this by defining equations only for as(t,st1,st2) so that we have equations at each node of the tree
* 

VARIABLES
        MIU(t,st1,st2)          Emission control rate GHGs
        TATM(t,st1,st2)         Increase temperature of atmosphere (degrees C from 1900)
        E(t,st1,st2)            Total CO2 emissions (GtCO2 per year)
        EIND(t,st1,st2)         Industrial emissions (GtCO2 per year)
        C(t,st1,st2)            Consumption (trillions 2005 US dollars per year)
        K(t,st1,st2)            Capital stock (trillions 2005 US dollars)
        CPC(t,st1,st2)          Per capita consumption (thousands 2005 USD per year)
        I(t,st1,st2)            Investment (trillions 2005 USD per year)
        S(t,st1,st2)            Gross savings rate as fraction of gross world product
        RI(t,st1,st2)           Real interest rate (per annum)
        Y(t,st1,st2)            Gross world product net of abatement and damages (trillions 2005 USD per year)
        YGROSS(t,st1,st2)       Gross world product GROSS of abatement and damages (trillions 2005 USD per year)
        YNET(t,st1,st2)         Output net of damages equation (trillions 2005 USD per year)
        DAMAGES(t,st1,st2)      Damages (trillions 2005 USD per year)
        DAMFRAC(t,st1,st2)      Damages as fraction of gross output
        ABATECOST(t,st1,st2)    Cost of emissions reductions  (trillions 2005 USD per year)
        MCABATE(t,st1,st2)      Marginal cost of abatement (2005$ per ton CO2)
        CCA(t,st1,st2)          Cumulative carbon emissions (GTC)
        CCATOT(t,st1,st2)       Total carbon emissions (GtC)
        PERIODU(t,st1,st2)      One period utility function
        CPRICE(t,st1,st2)       Carbon price (2005$ per ton of CO2)
        CEMUTOTPER(t,st1,st2)   Period utility
        UTILITY(st1,st2)        Welfare function for each possible scenario
        STOCH_UTILITY           Expected value of Welfare function
        
*Probabilities
        PROB(st1,st2)           Probability of a scenario
        H1(st1,st2)             Probability of the first catastrophy in a node 
        H2(st1,st2)             Probability of the second catastrophy in a node
        
;
         

NONNEGATIVE VARIABLES  MIU, TATM, Y, YGROSS, C, K, I;

EQUATIONS
*Emissions and Damages
        EEQ(t,st1,st2)           Emissions equation
        EINDEQ(t,st1,st2)        Industrial emissions
        CCACCA(t,st1,st2)        Cumulative industrial carbon emissions
        CCATOTEQ(t,st1,st2)        Cumulative total carbon emissions
        DAMFRACEQ(t,st1,st2)     Equation for damage fraction
        DAMEQ(t,st1,st2)         Damage equation
        ABATEEQ(t,st1,st2)       Cost of emissions reductions equation
        MCABATEEQ(t,st1,st2)     Equation for MC abatement
        CARBPRICEEQ(t,st1,st2)   Carbon price equation from abatement

*Climate and carbon cycle
        TATMEQ(t,st1,st2)        Temperature-climate equation for atmosphere

*Economic variables
        YGROSSEQ(t,st1,st2)      Output gross equation
        YNETEQ(t,st1,st2)        Output net of damages equation
        YY(t,st1,st2)            Output net equation
        CC(t,st1,st2)            Consumption equation
        CPCE(t,st1,st2)          Per capita consumption definition
        SEQ(t,st1,st2)           Savings rate equation
        KK(t,st1,st2)            Capital balance equation
        RIEQ(t,st1,st2)          Interest rate equation

* Utility
        CEMUTOTPEREQ(t,st1,st2)  Period utility
        PERIODUEQ(t,st1,st2)     Instantaneous utility function equation
        UTIL(st1,st2)            Objective function
        STOCH_UTIL
        
* Probability
        HEQ_1(st1,st2)
        HEQ_2(st1,st2)
        PROBEQ(st1,st2)          Function for the probs
        PROB0(st1,st2)           Ensures the scenario with no cat
        PROB1(st1,st2)
        PROB2(st1,st2)
;

** Equations of the model
* all the equations are defined on the dinamic set as(t,st1,st2)

alias(st1,index1,n_index1);
alias(st2,index2,n_index2);

*Emissions and Damages
 eeq(as(t,st1,st2))..             E(t,st1,st2)           =E= EIND(t,st1,st2) + etree(t);
 
* changes when we incour in tipping 2
 eindeq(as(t,st1,st2))..          EIND(t,st1,st2)        =E= (sigma(t) * YGROSS(t,st1,st2) * (1-(MIU(t,st1,st2))))$(st2.val eq 0)+
                                                             (sigma(t) * YGROSS(t,st1,st2) * (1-(MIU(t,st1,st2))))*eind_cat$(st2.val ne 0)
 ;
 

* For the equations that link t+1 with t, we need to better specify which node is the "previous" one, since nodes like (10,10,0)
* follow not (9,10,0) but (9,0,0) instead. The time needs to be always greater than the catastrophy time
*
 ccacca(as(t+1,st1,st2))..        CCA(t+1,st1,st2)       =E= (CCA(t,st1,st2)+ EIND(t,st1,st2)*tstep/3.666)$(t.val ge st1.val and t.val ge st2.val) +
                                                             (CCA(t,'0',st2)+ EIND(t,'0',st2)*tstep/3.666)$(t.val lt st1.val and t.val ge st2.val) +
                                                             (CCA(t,st1,'0')+ EIND(t,st1,'0')*tstep/3.666)$(t.val ge st1.val and t.val lt st2.val) +
                                                             (CCA(t,'0','0')+ EIND(t,'0','0')*tstep/3.666)$(t.val lt st1.val and t.val lt st2.val)
 

;
 


 ccatoteq(as(t,st1,st2))..        CCATOT(t,st1,st2)      =E= CCA(t,st1,st2)+cumetree(t);

* Here we defined new damage functions for the scenarios were a catastrophy has occoured 

* damfraceq(as(t,st1,st2)) ..      DAMFRAC(t,st1,st2)     =E= (((a1*TATM(t,st1,st2))+(a2*TATM(t,st1,st2)**a3))*1.3)$(st1.val ne 0 and st2.val eq 0) +
*                                                            (((a1*TATM(t,st1,st2))+(a2*TATM(t,st1,st2)**a3))*1.4)$(st1.val eq 0 and st2.val ne 0 ) +
*                                                            (((a1*TATM(t,st1,st2))+(a2*TATM(t,st1,st2)**a3))*2)$(st1.val ne 0 and st2.val ne 0) +
*                                                            (((a1*TATM(t,st1,st2))+(a2*TATM(t,st1,st2)**a3)))$(st1.val eq 0 and st2.val eq 0);
 

 damfraceq(as(t,st1,st2))..       DAMFRAC(t,st1,st2)      =E=  ((1-1/(1+(TATM(t,st1,st2)/w1)**2+(TATM(t,st1,st2)/w2)**w3))*c1_a*ifopt)$(st1.val ne 0 and st2.val eq 0) +
                                                            ((1-1/(1+(TATM(t,st1,st2)/w1)**2+(TATM(t,st1,st2)/w2)**w3))*c1_b*ifopt)$(st1.val eq 0 and st2.val ne 0 ) +
                                                            ((1-1/(1+(TATM(t,st1,st2)/w1)**2+(TATM(t,st1,st2)/w2)**w3))*c2*ifopt)$(st1.val ne 0 and st2.val ne 0) +
                                                            ((1-1/(1+(TATM(t,st1,st2)/w1)**2+(TATM(t,st1,st2)/w2)**w3))*1*ifopt)$(st1.val eq 0 and st2.val eq 0) ;
 
*---------------------------------------------------------------------------------
                                           

 dameq(as(t,st1,st2))..           DAMAGES(t,st1,st2)     =E= YGROSS(t,st1,st2) * DAMFRAC(t,st1,st2);
 abateeq(as(t,st1,st2))..         ABATECOST(t,st1,st2)   =E= YGROSS(t,st1,st2) * cost1(t) * (MIU(t,st1,st2)**expcost2);
 mcabateeq(as(t,st1,st2))..       MCABATE(t,st1,st2)     =E= pbacktime(t) * MIU(t,st1,st2)**(expcost2-1);
 carbpriceeq(as(t,st1,st2))..     CPRICE(t,st1,st2)      =E= pbacktime(t) * (MIU(t,st1,st2))**(expcost2-1);

*Climate

* sensitivity changes when we incour in tipping 1
 tatmeq(as(t+1,st1,st2))..        TATM(t+1,st1,st2)      =E=  (CCATOT(t+1,st1,st2)*tcre/1000 + 0.75*forcoth(t+1))$(st1.val eq 0) +
                                                              (CCATOT(t+1,st1,st2)*tcre_cat/1000 + 0.75*forcoth(t+1))$(st1.val ne 0)
 ;

*Economic variables
 ygrosseq(as(t,st1,st2))..        YGROSS(t,st1,st2)      =E= (al(t)*(L(t)/1000)**(1-GAMA))*(K(t,st1,st2)**GAMA);
 yneteq(as(t,st1,st2))..          YNET(t,st1,st2)        =E= YGROSS(t,st1,st2)*(1-damfrac(t,st1,st2));
 yy(as(t,st1,st2))..              Y(t,st1,st2)           =E= YNET(t,st1,st2) - ABATECOST(t,st1,st2);
 cc(as(t,st1,st2))..              C(t,st1,st2)           =E= Y(t,st1,st2) - I(t,st1,st2);
 cpce(as(t,st1,st2))..            CPC(t,st1,st2)         =E= 1000 * C(t,st1,st2) / L(t);
 seq(as(t,st1,st2))..             I(t,st1,st2)           =E= S(t,st1,st2) * Y(t,st1,st2);
 

 kk(as(t+1,st1,st2))..            K(t+1,st1,st2)         =L= ((1-dk)**tstep * K(t,st1,st2) + tstep * I(t,st1,st2))$(t.val ge st1.val and t.val ge st2.val) +
                                                             ((1-dk)**tstep * K(t,'0',st2) + tstep * I(t,'0',st2))$(t.val lt st1.val and t.val ge st2.val) +
                                                            ((1-dk)**tstep * K(t,st1,'0') + tstep * I(t,st1,'0'))$(t.val ge st1.val and t.val lt st2.val) +
                                                            ((1-dk)**tstep * K(t,'0','0') + tstep * I(t,'0','0'))$(t.val lt st1.val and t.val lt st2.val)
;

 rieq(as(t+1,st1,st2))..           CPC(t+1,st1,st2) =E= (((RI(t,st1,st2) + 1)/(1+prstp))**(tstep/elasmu)*CPC(t,st1,st2))$(t.val ge st1.val and t.val ge st2.val) +
                                                        (((RI(t,'0',st2) + 1)/(1+prstp))**(tstep/elasmu)*CPC(t,'0',st2))$(t.val lt st1.val and t.val ge st2.val) +
                                                        (((RI(t,st1,'0') + 1)/(1+prstp))**(tstep/elasmu)*CPC(t,st1,'0'))$(t.val ge st1.val and t.val lt st2.val) +
                                                        (((RI(t,'0','0') + 1)/(1+prstp))**(tstep/elasmu)*CPC(t,'0','0'))$(t.val lt st1.val and t.val lt st2.val)
 

;
 


*Utility
 cemutotpereq(as(t,st1,st2))..    CEMUTOTPER(t,st1,st2)  =E= PERIODU(t,st1,st2) * L(t) * rr(t);
 periodueq(as(t,st1,st2))..       PERIODU(t,st1,st2)     =E= ((C(t,st1,st2)*1000/L(T))**(1-elasmu)-1)/(1-elasmu)-1;
 

* For the utility we need to sum all the nodes of a scenario

 util(st1,st2)..               UTILITY(st1,st2)       =E= tstep * scale1 * (sum(t $ (t.val ge st1.val and t.val ge st2.val),  CEMUTOTPER(t,st1,st2)) +
                                                                            sum(t $ (t.val lt st1.val and t.val lt st2.val),  CEMUTOTPER(t,'0','0')) +
                                                                            sum(t $ (t.val ge st1.val and t.val lt st2.val),  CEMUTOTPER(t,st1,'0')) +
                                                                            sum(t $ (t.val lt st1.val and t.val ge st2.val),  CEMUTOTPER(t,'0',st2))) + scale2
 
;


*Stocastic things

alias(t,tt);

* hazard rate and probabilities definitions

heq_1(st1,st2)..            H1(st1,st2) =E= sum(tt $ (tt.val eq st1.val),(tanh((beta*(TATM(tt,'0','0')-tref)))+1)/2)$(st1.val le st2.val) +
                                            sum(tt $ (tt.val eq st1.val),(tanh((beta*(TATM(tt,'0',st2)-tref)))+1)/2)$(st1.val gt st2.val);
heq_2(st1,st2)..            H2(st1,st2) =E= sum(tt $ (tt.val eq st2.val),(tanh((beta*(TATM(tt,'0','0')-tref)))+1)/2)$(st2.val le st1.val) +
                                            sum(tt $ (tt.val eq st2.val),(tanh((beta*(TATM(tt,st1,'0')-tref)))+1)/2)$(st2.val gt st1.val);

probeq(st1,st2)$(st1.val ne 0 and st2.val ne 0)..              PROB(st1,st2) =E= prod((index1,index2) $((index1.val lt st1.val) and (index1.val eq index2.val)),1-H1(index1,index2)-H2(index1,index2)+H1(index1,index2)*H2(index1,index2))*
                                                                                        prod(index2 $ ((index2.val ge st1.val) and (index2.val lt st2.val)),1-H2(st1,index2))*
                                                                                        sum(n_index2 $(n_index2.val eq st1.val),H1(st1,n_index2))*
                                                                                        H2(st1,st2)$(st1.val lt st2.val) +
                                                                                        prod((index1,index2) $((index1.val lt st2.val) and (index1.val eq index2.val)),1-H1(index1,index2)-H2(index1,index2)+H1(index1,index2)*H2(index1,index2))*
                                                                                        prod(index1 $ ((index1.val ge st2.val) and (index1.val lt st1.val)),1-H1(index1,st2))*
                                                                                        sum(n_index1 $(n_index1.val eq st2.val),H2(n_index1,st2))*
                                                                                        H1(st1,st2)$(st1.val gt st2.val) +
                                                                                        prod((index1,index2) $((index1.val lt st2.val) and (index1.val eq index2.val)),1-H1(index1,index2)-H2(index1,index2)+H1(index1,index2)*H2(index1,index2))*
                                                                                        H1(st1,st2)*H2(st1,st2)$(st1.val eq st2.val);


prob1(st1,st2)$((st1.val eq 0 and st2.val ne 0))..                 PROB(st1,st2) =E= prod((index1,index2) $((index1.val lt st2.val) and (index1.val eq index2.val)),1-H1(index1,index2)-H2(index1,index2)+H1(index1,index2)*H2(index1,index2))*
                                                                                     prod(index1 $ (index1.val ge st2.val),1-H1(index1,st2))*
                                                                                     H2(st1,st2);

prob2(st1,st2)$((st1.val ne 0 and st2.val eq 0))..                 PROB(st1,st2) =E= prod((index1,index2) $((index1.val lt st1.val) and (index1.val eq index2.val)),1-H1(index1,index2)-H2(index1,index2)+H1(index1,index2)*H2(index1,index2))*
                                                                                     prod(index2 $ (index2.val ge st1.val),1-H2(st1,index2))*
                                                                                     H1(st1,st2); 

prob0(st1,st2)$(st1.val eq 0 and st2.val eq 0)..                PROB(st1,st2) =E= 1- sum((index1,index2) $(not(index1.val eq 0 and index2.val eq 0)), PROB(index1,index2));



* Expected Utility

stoch_util..              STOCH_UTILITY =E= sum((st1,st2),UTILITY(st1,st2)*PROB(st1,st2));

 
* GIVE SOME INERTIA

equations miu_control_up(t,st1,st2)
          miu_control_down(t,st1,st2)
;

miu_control_up(as(t+1,st1,st2)).. miu(t+1,st1,st2) =L= 0.1 + miu(t,st1,st2)$(t.val ge st1.val and t.val ge st2.val) +
                                                          miu(t,'0',st2)$(t.val lt st1.val and t.val ge st2.val) +
                                                          miu(t,st1,'0')$(t.val ge st1.val and t.val lt st2.val) +
                                                          miu(t,'0','0')$(t.val lt st1.val and t.val lt st2.val)
;                                                          


miu_control_down(as(t+1,st1,st2)).. miu(t+1,st1,st2) =G= -0.075 + miu(t,st1,st2)$(t.val ge st1.val and t.val ge st2.val) +
                                                          miu(t,'0',st2)$(t.val lt st1.val and t.val ge st2.val) +
                                                          miu(t,st1,'0')$(t.val ge st1.val and t.val lt st2.val) +
                                                          miu(t,'0','0')$(t.val lt st1.val and t.val lt st2.val)
;                                                          


* just to check that our computations are right

variable SOMMA_PROB;
equation somma_probeq;

somma_probeq.. SOMMA_PROB =E= sum((st1,st2),PROB(st1,st2));

*---------------------------------

*Resource limit
CCA.up(t,st1,st2)       = fosslim;

* Control rate limits
MIU.up(t,st1,st2)            = limmiu;
MIU.up(t,st1,st2)$(t.val<10) = 1;


**  Upper and lower bounds for stability
K.LO(t,st1,st2)         = 1;
C.LO(t,st1,st2)         = 2;
CPC.LO(t,st1,st2)       = .01;
TATM.UP(t,st1,st2)      = 12;

* Control variables
set lag10(t) ;
lag10(t) =  yes$(t.val gt card(t)-10);
S.FX(lag10(t),st1,st2) = optlrsav;

* Initial conditions
CCA.FX(tfirst,st1,st2)    = 400;
K.FX(tfirst,st1,st2)      = k0;
TATM.FX(tfirst,st1,st2)   = tatm0;


** Solution options
option iterlim = 99900;
option reslim = 99999;
option solprint = on;
option limrow = 0;
option limcol = 0;
model  CO2 /all/;

* For base run, this subroutine calculates Hotelling rents
* Carbon price is maximum of Hotelling rent or baseline price
* The cprice equation is different from 2013R. Not sure what went wrong.
If (ifopt eq 0,
       a2 = 0;
       solve CO2 maximizing STOCH_UTILITY using nlp;
*       photel(t)=cprice.l(t);
*       a2 = a20;
*      cprice.up(t)$(t.val<tnopol+1) = max(photel(t),cpricebase(t));
);



miu.fx('1',st1,st2)$(ifopt=1) = miu0;



           
If ( montecarlo eq 0,
  
     solve co2 maximizing stoch_utility using nlp;
 )



* Here you can change the number of iterations for montecarlo
set iter iterations / i1*i100 /;


parameters report(*,*,*,*,*) report file
           report2(*,*,*,*);
                
If ( montecarlo eq 1,

    loop(iter,
* change the parameter of interest
        eind_cat = 1.1 + ord(iter)*0.003;

            

* solve model
            solve co2 maximizing stoch_utility using nlp;
* save values that ar eof interest
* if the variable has 3 indices use report, if it has only 2 use report2
            report('TATM',t,st1,st2,iter) = TATM.l(t,st1,st2);
            report('E',t,st1,st2,iter)    = E.l(t,st1,st2);
            report('DAMFRAC',t,st1,st2,iter)    = DAMFRAC.l(t,st1,st2);
            report('CCATOT',t,st1,st2,iter)    = CCATOT.l(t,st1,st2);
            report('MIU',t,st1,st2,iter)    = MIU.l(t,st1,st2);
            report2('H1',st1,st2,iter)    = H1.l(st1,st2);
            report2('H2',st1,st2,iter)    = H2.l(st1,st2);
        );
* unload report parameter to gdx
            execute_unload 'rep.gdx', report , report2;

)




** POST-SOLVE
* Calculate social cost of carbon and other variables
*scc(t) = -1000*eeq.m(t)/(.00001+cc.m(t));
