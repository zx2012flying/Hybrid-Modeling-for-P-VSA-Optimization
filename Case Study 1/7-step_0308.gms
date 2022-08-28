Set

comp                     'propane PA: pass, propene PE: adsorb'                          /PA, PE/

t                        'time discretization'                                           /0 * 70/
t_start(t)               'starting point'                                                /0/
t_pre(t)                 'time discretization for pressurization'                        /1 * 10/
t_ads(t)                 'time discretization for adsorption'                            /11 * 20/
t_equ(t)                 'time discretization for equalization'                          /21 * 30/
t_rin(t)                 'time discretization for rinse'                                 /31 * 40/
t_des(t)                 'time discretization for desorption'                            /41 * 50/
t_eqn(t)                 'time discretization for equalization'                          /51 * 60/
t_ref(t)                 'time discretization for rinse reflux'                          /61 * 70/

z                        'column axis discretization'                                    /0 * 21/
z_start(z)               'axis starting point'                                           /1/
z_mid(z)                 'axis middle section'                                           /2 * 19/
z_end(z)                 'axis ending point'                                             /20/
z_left(z)                'left ghost'                                                    /0/
z_axis(z)                'axis point'                                                    /1 * 20/
z_right(z)               'right ghost'                                                   /21/;

Parameters

u_ini(t, z)              'velocity initial value:                m/s'
P_ini(t, z)              'pressure initial value:                Pa'
yE_ini(t, z)             'PE composition initial value:          /'
yA_ini(t, z)             'PE composition initial value:          /';

$CALL 'del 7_step_PSA_0308.gdx'
$CALL GDXXRW I=7_step_PSA_0308.xlsx O=7_step_PSA_0308.gdx Index=In_index!A1

$GDXIN 7_step_PSA_0308.gdx
$LOAD u_ini, P_ini, yE_ini, yA_ini

Parameters

yE_feed                  'input PE composition'                                          /0.85/
yA_feed                  'input PA composition'
yE_rinse                 'rinse input PE composition'                                    /0.99/
yA_rinse                 'rinse input PA composition'
P_input                  'input pressure from FCC unit           Pa'                     /202650/
P_atm                    'atmosphere pressure                    Pa'                     /101325/
T_feed                   'feed temperature                       K'                      /323/
T_amb                    'ambient temperature                    K'                      /298/

deltaH_PA                'adsorption of heat                     J/mol'                  /15605/
deltaH_PE                'adsorption of heat                     J/mol'                  /28265/
MW_PA                    'molecular weight of PA:                kg/mol'                 /0.0441/
MW_PE                    'molecular weight of PA:                kg/mol'                 /0.04208/
kA                       'mass transfer coefficient PA:          1/s'                    /7.56e-6/
kE                       'mass transfer coefficient PE:          1/s'                    /2.24e-2/

Length                   'length of column pilot plant           m'                      /10/
Diameter                 'diameter of column pilot plant         m'
R_inner                  'inner radius                           m'
R_out                    'outer radius                           m'
a_W                      'cross section area                     m2'
CSA                      'cross section area                     m2'
Porosity_bed             'porosity of column bed                 /'                      /0.45/
Nporosity                '1 - Porosity_bed                       '
rho_s                    'particle density of solid adsorbent:   kg/m3'                  /1210/
viscosity                'viscosity of gas:                      kg/m/s'                 /8e-6/
Radius_SP                'solid particle radius                  m'                      /8e-4/
Cp_g                     'heat capacity of gases                 J/mol/K'                /119/
Cp_s                     'heat capacity of solid                 J/kg/K'                 /920/
Cp_a                     'heat capacity of air                   J/mol/K'                /30/
K_z                      'thermal conductivity of gases          J/m/s/K'                /0.0903/

K_W                      'thermal diffuusivity of wall:          J/m/s/K'                /16/
rho_W                    'density of wall                        kg/m3'                  /8238/
Cp_W                     'heat capacity of wall                  J/kg/K'                 /500/
h_W_in                   'heat transfer coefficient inner        W/m2/K'                 /8.6/
h_W_out                  'heat transfer coefficient outer        W/m2/K'                 /2.5/

pi                       'pi'                                                            /3.1416/
R                        'gas constant                           kg*m2/k/mol/s2'         /8.3145/
Eff                      '1 / efficiency of vacuum and compressor'                       /1.333/
IC                       'isentropic coeff for PE k = 1.15,      k/k-1'                  /7.67/

Nubz                     'number of z interval+1'
Dz                       'delta z: m'
DzS                      'square of delta z'
KDzS                     'K_z / Dz / Dz'
t_total                  'total duration                                 60'
t_equ_equ                'duration between equ and equ'
t_rin_rin                'duration between rinse and rinse reflux'

EM                       'parameter Ergun equation               kg/m3/s '
FM                       'parameter Ergun equation               1/m     '
M                        'parameter in mass balance                      '
ERTI                     'Eff * R * Temp * IC                             0.026'
CPRT                     'CSA * Porosity_bed / (R * Temp)'
EB_1                     'parameter energy balance'
EB_2                     'parameter energy balance'
EB_3                     'parameter energy balance'
CPR                      'parameter energy balance Cp_g * Porosity_bed / R'
CPDR                     'parameter energy balance Cp_g * Porosity_bed / Dz / R'
THR                      'parameter energy balance 2 * h_W_in / R_inner'
NPCC                     'parameter energy balance Nporosity * (Cp_g - Cp_a) * Tc(t, z)';

yA_feed = 1 - yE_feed;
yA_rinse = 1 - yE_rinse;

Diameter = Length / 5;
R_inner = Diameter / 2;
R_out = R_inner + 0.06;
a_W = pi * (R_out * R_out - R_inner * R_inner);
CSA = pi * Diameter * Diameter / 4;
Nporosity = 1 - Porosity_bed;

Nubz = card(z_axis)+1;
Dz = Length / card(z_axis);
DzS = Dz * Dz;
KDzS = K_z / Dz / Dz;
t_total = card(t) - 1;
t_equ_equ = card(t_equ) + card(t_rin) + card(t_des);
t_rin_rin = card(t_rin) + card(t_des) + card(t_eqn);

EM = 150 * viscosity * (1 - Porosity_bed) * (1 - Porosity_bed) / (4 * Radius_SP * Radius_SP * Porosity_bed * Porosity_bed * Porosity_bed);
FM = 1.75 * (1 - Porosity_bed) / (2 * Radius_SP * Porosity_bed * Porosity_bed * Porosity_bed);
M = (1 - Porosity_bed) * R * rho_s / Porosity_bed;
ERTI = Eff * IC * R * T_feed / 1e6;
CPRT = CSA * Porosity_bed / R;
EB_1 = K_W / (rho_W * Cp_W);
EB_2 = 2 * pi * h_W_in * R_inner / (Cp_W * rho_W * a_W);
EB_3 = 2 * pi * h_W_out * R_out / (Cp_W * rho_W * a_W);
CPR = Cp_g * Porosity_bed / R;
CPDR = Cp_g * Porosity_bed / Dz / R;
THR = 2 * h_W_in / R_inner;
NPCC = Nporosity * (Cp_g - Cp_a);

Parameters

Qsat_1_PA                'isotherm parameter: mol/kg'                                    /1.7527/
Qsat_2_PA                'isotherm parameter: mol/kg'                                    /0/
b0_1_PA                  'isotherm parameter'                                            /-5.4727/
b0_2_PA                  'isotherm parameter'                                            /0/
H_1_PA                   'heat of adsorption / R: J/mol'                                 /1952/
H_2_PA                   'heat of adsorption / R: J/mol'                                 /0/

Qsat_1_PE                'isotherm parameter: mol/kg'                                    /1.1866/
Qsat_2_PE                'isotherm parameter: mol/kg'                                    /0.7656/
b0_1_PE                  'isotherm parameter'                                            /-5.4059/
b0_2_PE                  'isotherm parameter'                                            /-14.603/
H_1_PE                   'heat of adsorption / R: J/mol'                                 /2452.3/
H_2_PE                   'heat of adsorption / R: J/mol'                                 /6977/;

Positive variables

Dt_pre                   'time interval at pressurization                second'
Dt_ads                   'time interval at adsorption                    second'
Dt_equ                   'time interval at equialization one             second'
Dt_rin                   'time interval at resin                         second'
Dt_des                   'time interval at desorption                    second'
Dt_eqn                   'time interval at equialization two             second'
Dt_ref                   'time interval at rinse reflux                  second'
P_high                   'high pressure for adsorption                   Pa'
P_low                    'low pressure for desorption                    Pa'
P_end                    'adsorption rinse end pressure                  Pa'

yA(t, z)                 'gas molar fraction of component PA             /'
qA(t, z)                 'solid loading of PA                            mol/kg'
qA_Star(t, z)            'equilibrium loading of PA                      mol/kg'
yE(t, z)                 'gas molar fraction of component PE             /'
qE(t, z)                 'solid loading of PE                            mol/kg'
qE_Star(t, z)            'equilibrium loading of PE                      mol/kg'
P(t, z)                  'pressure in the column                         Pa'

b_1_PA(t, z)             'isotherm variables                             Pa-1'
b_2_PA(t, z)             'isotherm variables                             Pa-1'
b_1_PE(t, z)             'isotherm variables                             Pa-1'
b_2_PE(t, z)             'isotherm variables                             Pa-1'
Tc(t, z)                 'temperature in the column                      K'
Tw(t, z)                 'wall temperature                               K'
EB(t, z)                 'demoninator in energy balance                  '

PA_pre_feed              'PA feed at pre: / R / Temp                     mol'
PA_ads_feed              'PA feed at ads: / R / Temp                     mol'
PA_ads_out               'PA output at pre: / R / Temp                   mol'
PA_rin_feed              'PA feed at rin: / R / Temp                     mol'
PA_des_exhaust           'PA exhaust at des: / R / Temp                  mol'

PE_pre_feed              'PE feed at pre: / R / Temp                     mol'
PE_ads_feed              'PE feed at ads: / R / Temp                     mol'
PE_ads_out               'PE output at pre: / R / Temp                   mol'
PE_rin_feed              'PE feed at rin: / R / Temp                     mol'
PE_des_exhaust           'PE exhaust at des: / R / Temp                  mol'

Equ_out                  'equalization one output                        mol'
yE_equ_out               'PE equalizaion one output                      mol'
Eqn_feed                 'equalization two feed                          mol'
Rinse_out                'rinse output                                   mol'
yE_rin_out               'PE molar fraction at rinse output              /'
Reflux_feed              'rinse reflux feed                              mol'
;

Variables

u(t, z)                  'gas velocity in the column                     m/s'

PurityA                  'purity of PA                                   /'
PurityE                  'purity of PE                                   /'
RecoveryA                'recovery of PA                                 /'
RecoveryE                'recovery of PE                                 /'
PE_product               '1% PA + 99% PE product                         mol'
Mass_balanceA            'overall mass balances of PA: input - output    mol'
Mass_balanceE            'overall mass balances of PE: input - output    mol'
Productivity

W_12_feed                'vacuum pump workload at pre+ads feed           MJ'
W_23_out                 'vacuum pump workload at ads+rin out            MJ'
W_3_feed                 'vacuum pump workload at rin feed               MJ'
W_4_out                  'vacuum pump workload at des out                MJ'

OBJ                      'objective function                             MJ/mol';

Dt_pre.lo = 1;
Dt_pre.up = 10;
Dt_pre.l = 8;

Dt_ads.lo = 1;
Dt_ads.up = 20;
Dt_ads.l = 8;

Dt_rin.lo = 1;
Dt_rin.up = 10;
Dt_rin.l = 1;

Dt_equ.lo = 0.1;
Dt_equ.up = 4;
Dt_equ.l = 1;

Dt_des.lo = 1;
Dt_des.up = 20;
Dt_des.l = 18;

Dt_eqn.lo = 0.001;
Dt_eqn.up = 500;
Dt_eqn.l = 1;

Dt_ref.lo = 0.1;
Dt_ref.up = 500;
Dt_ref.l = 1;

P_high.lo = P_atm;
P_high.up = 4e5;
P_high.l = 1.8e5;

P_low.lo = 1e3;
P_low.up = P_atm;
P_low.l = 4e4;

P_end.lo = P_atm;
P_end.up = 4e5;
P_end.l = 1.1e5;

yA.lo(t, z) = 0;
yA.up(t, z) = 1;
yA.l(t, z) = yA_ini(t, z);

yE.lo(t, z) = 0;
yE.up(t, z) = 1;
yE.l(t, z) = yE_ini(t, z);

qA.lo(t, z) = 0;
qA.up(t, z) = 3;
qA.l(t, z) = 0.1;

qE.lo(t, z) = 0;
qE.up(t, z) = 5;
qE.l(t, z) = 2;

qA_Star.lo(t, z) = 0;
qA_Star.up(t, z) = 3;
qA_Star.l(t, z) = 0.2;

qE_Star.lo(t, z) = 0;
qE_Star.up(t, z) = 5;
qE_Star.l(t, z) = 2;

P.lo(t, z) = 1e3;
P.up(t, z) = 5e5;
P.l(t, z) = P_ini(t, z);

u.lo(t, z) = -40;
u.up(t, z) = 40;
u.l(t, z) = u_ini(t, z);

b_1_PA.lo(t, z) = 0.5;
b_1_PA.up(t, z) = 3;
b_1_PA.l(t, z) = 1.75;

b_2_PA.lo(t, z) = 0;
b_2_PA.up(t, z) = 0.1;
b_2_PA.l(t, z) = 0;

b_1_PE.lo(t, z) = 1;
b_1_PE.up(t, z) = 20;
b_1_PE.l(t, z) = 9;

b_2_PE.lo(t, z) = 100;
b_2_PE.up(t, z) = 5000;
b_2_PE.l(t, z) = 1000;

Tc.lo(t, z) = 303;
Tc.up(t, z) = 353;
Tc.l(t, z) = 318;

Tw.lo(t, z) = 283;
Tw.up(t, z) = 323;
Tw.l(t, z) = 310;

EB.lo(t, z) = 1e5;
EB.up(t, z) = 1e6;
EB.l(t, z) = 6e5;

PA_pre_feed.lo = 0.1;
PA_pre_feed.up = 1e3;
PA_pre_feed.l = 500;

PA_ads_feed.lo = 0.1;
PA_ads_feed.up = 2e3;
PA_ads_feed.l = 500;

PA_ads_out.lo = 0.1;
PA_ads_out.up = 2e3;
PA_ads_out.l = 500;

PA_rin_feed.lo = 0.1;
PA_rin_feed.up = 1e3;
PA_rin_feed.l = 500;

PA_des_exhaust.lo = 0.1;
PA_des_exhaust.up = 1e3;
PA_des_exhaust.l = 1000;

PE_pre_feed.lo = 0.1;
PE_pre_feed.up = 1e4;
PE_pre_feed.l = 1000;

PE_ads_feed.lo = 0.1;
PE_ads_feed.up = 1e4;
PE_ads_feed.l = 1000;

PE_ads_out.lo = 1;
PE_ads_out.up = 1e4;
PE_ads_out.l = 1000;

PE_rin_feed.lo = 0.1;
PE_rin_feed.up = 1e4;
PE_rin_feed.l = 500;

PE_des_exhaust.lo = 0.1;
PE_des_exhaust.up = 1e4;
PE_des_exhaust.l = 1000;

PE_product.l = 1000;

Mass_balanceA.l = 0;
Mass_balanceE.l = 0;

W_12_feed.l = 1.2;
W_23_out.l = 0.5;
W_3_feed.l = 0.8;
W_4_out.l = 1;
Productivity.l = 50;

OBJ.l = 0.015;

Equ_out.l = 100;

yE_equ_out.lo = 0.1;
yE_equ_out.up = 1;
yE_equ_out.l = 0.85;

Eqn_feed.l = 100;
Rinse_out.l = 100;

yE_rin_out.lo = 0.5;
yE_rin_out.up = 1;
yE_rin_out.l = 0.85;

Reflux_feed.l = 100;

Equations

LGEq00, LGEq01, LGEq02, LGEq03, LGEq04, LGEq05, LGEq06, LGEq07
ISEq00, ISEq01, ISEq02, ISEq03, ISEq04, ISEq05
PREq00, PREq01, PREq02, PREq04, PREq05, PREq06, PREq07, PREq08, PREq10, PREq11, PREq12, PREq13, PREq14, PREq16, PREq17, PREq18
ADEq00, ADEq01, ADEq02, ADEq04, ADEq05, ADEq06, ADEq07, ADEq09, ADEq10, ADEq11, ADEq12, ADEq13, ADEq14, ADEq16, ADEq17, ADEq18
EQEq00, EQEq01, EQEq02, EQEq04, EQEq05, EQEq06, EQEq07, EQEq09, EQEq10, EQEq11, EQEq12, EQEq13, EQEq14, EQEq16, EQEq17, EQEq18
REEq00, REEq01, REEq02, REEq04, REEq05, REEq06, REEq07, REEq09, REEq10, REEq11, REEq12, REEq13, REEq14, REEq16, REEq17, REEq18
DEEq00, DEEq01, DEEq02, DEEq03, DEEq05, DEEq06, DEEq07, DEEq08, DEEq10, DEEq11, DEEq12, DEEq13, DEEq14, DEEq16, DEEq17, DEEq18
EAEq00, EAEq01, EAEq02, EAEq04, EAEq05, EAEq06, EAEq07, EAEq08, EAEq10, EAEq11, EAEq12, EAEq13, EAEq14, EAEq16, EAEq17, EAEq18
RREq00, RREq01, RREq02, RREq04, RREq05, RREq06, RREq07, RREq08, RREq10, RREq11, RREq12, RREq13, RREq14, RREq16, RREq17, RREq18
CSEq00, CSEq01, CSEq02, CSEq03, CSEq04, CSEq05
CCEq00, CCEq01, CCEq02, CCEq03, CCEq04, CCEq05, CCEq06, CCEq07
MBEq00, MBEq01, MBEq02, MBEq03, MBEq04, MBEq05, MBEq06, MBEq07, MBEq08, MBEq09
SSEq00, SSEq01, SSEq02, SSEq03, SSEq04, SSEq05, SSEq06, SSEq07, SSEq08, SSEq09, SSEq10
ECEq00, ECEq01, ECEq02, ECEq03, ECEq04
;

*****************************************
*          logical constraint           *
*****************************************
LGEq00(t, z)                             ..      yA(t, z)                           =e= 1 - yE(t, z)                                ;
LGEq01(t, z)$(t_ads(t) and z_axis(z))    ..      u(t, z)                            =g= 0                                           ;
LGEq02(t, z)$(t_equ(t) and z_axis(z))    ..      u(t, z)                            =g= 0                                           ;
LGEq03(t, z)$(t_rin(t) and z_axis(z))    ..      u(t, z)                            =g= 0                                           ;
LGEq04                                   ..      P_high                             =g= 1.05 * P_end                                ;
LGEq05                                   ..      P_end                              =g= 1.05 * P_low                                ;
LGEq06(t, z)$(z_left(z))                 ..      Tw(t, z)                           =e= T_amb                                       ;
LGEq07(t, z)$(z_right(z))                ..      Tw(t, z)                           =e= T_amb                                       ;
*****************************************
*      Dual-site Langmuir isotherm      *
*****************************************
ISEq00(t, z)$z_axis(z)                   ..      b_1_PA(t, z)                       =e= exp(b0_1_PA + H_1_PA / Tc(t, z))            ;
ISEq01(t, z)$z_axis(z)                   ..      b_2_PA(t, z)                       =e= b0_2_PA                                     ;
ISEq02(t, z)$z_axis(z)                   ..      b_1_PE(t, z)                       =e= exp(b0_1_PE + H_1_PE / Tc(t, z))            ;
ISEq03(t, z)$z_axis(z)                   ..      b_2_PE(t, z)                       =e= exp(b0_2_PE + H_2_PE / Tc(t, z))            ;
ISEq04(t, z)$z_axis(z)                   ..      qA_Star(t, z)                      =e= Qsat_1_PA * b_1_PA(t, z) * P(t, z) * yA(t, z) / (1e5 + b_1_PA(t, z) * P(t, z) * yA(t, z) + b_1_PE(t, z) * P(t, z) * yE(t, z)) + Qsat_2_PA * b_2_PA(t, z) * P(t, z) * yA(t, z) / (1e5 + b_2_PA(t, z) * P(t, z) * yA(t, z) + b_2_PE(t, z) * P(t, z) * yE(t, z));
ISEq05(t, z)$z_axis(z)                   ..      qE_Star(t, z)                      =e= Qsat_1_PE * b_1_PE(t, z) * P(t, z) * yE(t, z) / (1e5 + b_1_PA(t, z) * P(t, z) * yA(t, z) + b_1_PE(t, z) * P(t, z) * yE(t, z)) + Qsat_2_PE * b_2_PE(t, z) * P(t, z) * yE(t, z) / (1e5 + b_2_PA(t, z) * P(t, z) * yA(t, z) + b_2_PE(t, z) * P(t, z) * yE(t, z));
*****************************************
*         Pressurization step           *
*****************************************
PREq00(t, z)$(t_pre(t) and z_left(z))    ..      yE(t, z)                           =e= yE_feed                                    ;
PREq01(t, z)$(t_pre(t) and z_left(z))    ..      P(t, z)                            =e= P_high                                     ;
PREq02(t, z)$(t_pre(t) and z_left(z))    ..      Tc(t, z)                           =e= T_feed                                     ;

PREq04(t, z)$(t_pre(t) and z_right(z))   ..      yE(t, z)                           =e= yE(t, z - 1)                               ;
PREq05(t, z)$(t_pre(t) and z_right(z))   ..      P(t, z)                            =e= P(t, z - 1)                                ;
PREq06(t, z)$(t_pre(t) and z_right(z))   ..      u(t, z)                            =e= u(t, z - 1)                                ;
PREq07(t, z)$(t_pre(t) and z_right(z))   ..      u(t, z)                            =e= 0                                          ;
PREq08(t, z)$(t_pre(t) and z_right(z))   ..      Tc(t, z)                           =e= Tc(t, z - 1)                               ;

PREq10(t, z)$(t_pre(t) and z_axis(z))    ..      qA(t, z) - qA(t-1, z)              =e= Dt_pre * kA * (qA_Star(t, z) - qA(t, z))   ;
PREq11(t, z)$(t_pre(t) and z_axis(z))    ..      qE(t, z) - qE(t-1, z)              =e= Dt_pre * kE * (qE_Star(t, z) - qE(t, z))   ;

PREq12(t, z)$(t_pre(t) and z_axis(z))    ..      yE(t, z) - yE(t-1, z)              =e= - Dt_pre * u(t, z) * (yE(t, z) - yE(t, z-1)) / Dz - M * Tc(t, z) * ( (qE(t, z) - qE(t-1, z)) - yE(t, z) * (qE(t, z) - qE(t-1, z) + qA(t, z) - qA(t-1, z)) ) / P(t, z)                                                                                                                                      ;
PREq13(t, z)$(t_pre(t) and z_axis(z))    ..      (P(t, z) - P(t-1, z)) / Dt_pre     =e= - (P(t, z) * u(t, z) - P(t, z-1) * u(t, z-1)) / Dz - M * Tc(t, z) * ( qE(t,z) - qE(t-1, z) + qA(t,z) - qA(t-1, z) ) / Dt_pre + P(t, z) * (Tc(t, z) - Tc(t-1, z)) / (Tc(t, z) * Dt_pre) + P(t, z) * u(t, z) * (Tc(t, z) - Tc(t, z-1)) / (Tc(t, z) * Dz)                                                     ;
PREq14(t, z)$(t_pre(t) and z_axis(z))    ..      (P(t, z) - P(t, z-1)) / Dz         =e= - EM * u(t, z-1) - FM * u(t, z-1) * ((u(t, z-1)*u(t, z-1))**0.5) * P(t, z-1) * MW_PE / (R * Tc(t, z))                                                                                                                                                                                                      ;

PREq16(t, z)$(t_pre(t) and z_axis(z))    ..      EB(t, z)                           =e= Nporosity * ( rho_s * Cp_s + Cp_a * (qA(t, z) + qE(t, z)) ) + CPR * P(t, z) / Tc(t, z)                                                                                                                                                                                                                     ;
PREq17(t, z)$(t_pre(t) and z_axis(z))    ..      (Tc(t, z) - Tc(t-1, z)) * EB(t, z) =e= - CPDR * Dt_pre * P(t, z) * u(t, z) * (Tc(t, z) - Tc(t, z-1)) / Tc(t, z) + NPCC * Tc(t, z) * (qA(t, z) - qA(t-1, z) + qE(t, z) - qE(t-1, z)) + Nporosity * (deltaH_PA * qA(t, z) - deltaH_PA * qA(t-1, z) + deltaH_PE * qE(t, z) - deltaH_PE * qE(t-1, z)) - THR * Dt_pre * (Tc(t, z) - Tw(t, z))          ;
PREq18(t, z)$(t_pre(t) and z_axis(z))    ..      (Tw(t, z) - Tw(t-1, z)) / Dt_pre   =e= EB_1 * (Tw(t, z+1) + Tw(t, z-1) - 2 * Tw(t, z)) / Dz / Dz + EB_2 * (Tc(t, z) - Tw(t, z)) - EB_3 * (Tw(t, z) - T_amb)                                                                                                                                                                                       ;
*****************************************
*           adsorption step             *
*****************************************
ADEq00(t, z)$(t_ads(t) and z_left(z))    ..      yE(t, z)                           =e= yE_feed                                    ;
ADEq01(t, z)$(t_ads(t) and z_left(z))    ..      P(t, z)                            =e= P_high                                     ;
ADEq02(t, z)$(t_ads(t) and z_left(z))    ..      Tc(t, z)                           =e= T_feed                                     ;

ADEq04(t, z)$(t_ads(t) and z_right(z))   ..      yE(t, z)                           =e= yE(t, z - 1)                               ;
ADEq05(t, z)$(t_ads(t) and z_right(z))   ..      P(t, z)                            =e= P_end                                      ;
ADEq06(t, z)$(t_ads(t) and z_right(z))   ..      u(t, z)                            =e= u(t, z - 1)                                ;
ADEq07(t, z)$(t_ads(t) and z_right(z))   ..      Tc(t, z)                           =e= Tc(t, z - 1)                               ;

ADEq09(t, z)$(t_ads(t) and z_axis(z))    ..      qA(t, z) - qA(t-1, z)              =e= Dt_ads * kA * (qA_Star(t, z) - qA(t, z))   ;
ADEq10(t, z)$(t_ads(t) and z_axis(z))    ..      qE(t, z) - qE(t-1, z)              =e= Dt_ads * kE * (qE_Star(t, z) - qE(t, z))   ;

ADEq11(t, z)$(t_ads(t) and z_axis(z))    ..      yE(t, z) - yE(t-1, z)              =e= - Dt_ads * u(t, z) * (yE(t, z) - yE(t, z-1)) / Dz - M * Tc(t, z) * ( (qE(t, z) - qE(t-1, z)) - yE(t, z) * (qE(t, z) - qE(t-1, z) + qA(t, z) - qA(t-1, z)) ) / (P(t, z) + 1e-10)                                                                                                                            ;
ADEq12(t, z)$(t_ads(t) and z_axis(z))    ..      (P(t, z) - P(t-1, z)) / Dt_ads     =e= - (P(t, z) * u(t, z) - P(t, z-1) * u(t, z-1)) / Dz - M * Tc(t, z) * ( qE(t,z) - qE(t-1, z) + qA(t,z) - qA(t-1, z) ) / Dt_ads + P(t, z) * (Tc(t, z) - Tc(t-1, z)) / (1e-10 + Tc(t, z) * Dt_ads) + P(t, z) * u(t, z) * (Tc(t, z) - Tc(t, z-1)) / (1e-10 + Tc(t, z) * Dz)                                     ;
ADEq13(t, z)$(t_ads(t) and z_axis(z))    ..      (P(t, z) - P(t, z-1)) / Dz         =e= - EM * u(t, z-1) - FM * u(t, z-1) * ((u(t, z-1)*u(t, z-1))**0.5) * P(t, z-1) * MW_PE / (R * Tc(t, z))                                                                                                                                                                                                      ;
ADEq14(t, z)$(t_ads(t) and z_right(z))   ..      (P(t, z) - P(t, z-1)) / Dz         =e= - EM * u(t, z-1) - FM * u(t, z-1) * ((u(t, z-1)*u(t, z-1))**0.5) * P(t, z-1) * MW_PE / (R * Tc(t, z))                                                                                                                                                                                                      ;

ADEq16(t, z)$(t_ads(t) and z_axis(z))    ..      EB(t, z)                           =e= Nporosity * ( rho_s * Cp_s + Cp_a * (qA(t, z) + qE(t, z)) ) + CPR * P(t, z) / Tc(t, z)                                                                                                                                                                                                                     ;
ADEq17(t, z)$(t_ads(t) and z_axis(z))    ..      (Tc(t, z) - Tc(t-1, z)) * EB(t, z) =e= - CPDR * Dt_ads * P(t, z) * u(t, z) * (Tc(t, z) - Tc(t, z-1)) / Tc(t, z) + NPCC * Tc(t, z) * (qA(t, z) - qA(t-1, z) + qE(t, z) - qE(t-1, z)) + Nporosity * (deltaH_PA * qA(t, z) - deltaH_PA * qA(t-1, z) + deltaH_PE * qE(t, z) - deltaH_PE * qE(t-1, z)) - THR * Dt_ads * (Tc(t, z) - Tw(t, z))          ;
ADEq18(t, z)$(t_ads(t) and z_axis(z))    ..      (Tw(t, z) - Tw(t-1, z)) / Dt_ads   =e= EB_1 * (Tw(t, z+1) + Tw(t, z-1) - 2 * Tw(t, z)) / Dz / Dz + EB_2 * (Tc(t, z) - Tw(t, z)) - EB_3 * (Tw(t, z) - T_amb)                                                                                                                                                                                       ;
*****************************************
*       Equalization step one           *
*****************************************
EQEq00(t, z)$(t_equ(t) and z_left(z))    ..      yE(t, z)                           =e= yE(t, z + 1)                               ;
EQEq01(t, z)$(t_equ(t) and z_left(z))    ..      P(t, z)                            =e= P(t, z + 1)                                ;
EQEq02(t, z)$(t_equ(t) and z_left(z))    ..      Tc(t, z)                           =e= Tc(t, z + 1)                               ;

EQEq04(t, z)$(t_equ(t) and z_right(z))   ..      yE(t, z)                           =e= yE(t, z - 1)                               ;
EQEq05(t, z)$(t_equ(t) and z_right(z))   ..      P(t, z)                            =e= P_end                                      ;
EQEq06(t, z)$(t_equ(t) and z_right(z))   ..      u(t, z)                            =e= u(t, z - 1)                                ;
EQEq07(t, z)$(t_equ(t) and z_right(z))   ..      Tc(t, z)                           =e= Tc(t, z - 1)                               ;

EQEq09(t, z)$(t_equ(t) and z_axis(z))    ..      qA(t, z) - qA(t-1, z)              =e= Dt_equ * kA * (qA_Star(t, z) - qA(t, z))   ;
EQEq10(t, z)$(t_equ(t) and z_axis(z))    ..      qE(t, z) - qE(t-1, z)              =e= Dt_equ * kE * (qE_Star(t, z) - qE(t, z))   ;

EQEq11(t, z)$(t_equ(t) and z_axis(z))    ..      yE(t, z) - yE(t-1, z)              =e= - Dt_equ * u(t, z) * (yE(t, z) - yE(t, z-1)) / Dz - M * Tc(t, z) * ( (qE(t, z) - qE(t-1, z)) - yE(t, z) * (qE(t, z) - qE(t-1, z) + qA(t, z) - qA(t-1, z)) ) / (P(t, z) + 1e-10)                                                                                                                            ;
EQEq12(t, z)$(t_equ(t) and z_axis(z))    ..      (P(t, z) - P(t-1, z)) / Dt_equ     =e= - (P(t, z) * u(t, z) - P(t, z-1) * u(t, z-1)) / Dz - M * Tc(t, z) * ( qE(t,z) - qE(t-1, z) + qA(t,z) - qA(t-1, z) ) / Dt_equ + P(t, z) * (Tc(t, z) - Tc(t-1, z)) / (1e-10 + Tc(t, z) * Dt_equ) + P(t, z) * u(t, z) * (Tc(t, z) - Tc(t, z-1)) / (1e-10 + Tc(t, z) * Dz)                                     ;
EQEq13(t, z)$(t_equ(t) and z_axis(z))    ..      (P(t, z) - P(t, z-1)) / Dz         =e= - EM * u(t, z-1) - FM * u(t, z-1) * ((u(t, z-1)*u(t, z-1))**0.5) * P(t, z-1) * MW_PE / (R * Tc(t, z))                                                                                                                                                                                                      ;
EQEq14(t, z)$(t_equ(t) and z_right(z))   ..      (P(t, z) - P(t, z-1)) / Dz         =e= - EM * u(t, z-1) - FM * u(t, z-1) * ((u(t, z-1)*u(t, z-1))**0.5) * P(t, z-1) * MW_PE / (R * Tc(t, z))                                                                                                                                                                                                      ;

EQEq16(t, z)$(t_equ(t) and z_axis(z))    ..      EB(t, z)                           =e= Nporosity * ( rho_s * Cp_s + Cp_a * (qA(t, z) + qE(t, z))) + CPR * P(t, z) / Tc(t, z)                                                                                                                                                                                                                      ;
EQEq17(t, z)$(t_equ(t) and z_axis(z))    ..      (Tc(t, z) - Tc(t-1, z)) * EB(t, z) =e= - CPDR * Dt_equ * P(t, z) * u(t, z) * (Tc(t, z) - Tc(t, z-1)) / Tc(t, z) + NPCC * Tc(t, z) * (qA(t, z) - qA(t-1, z) + qE(t, z) - qE(t-1, z)) + Nporosity * (deltaH_PA * qA(t, z) - deltaH_PA * qA(t-1, z) + deltaH_PE * qE(t, z) - deltaH_PE * qE(t-1, z)) - THR * Dt_equ * (Tc(t, z) - Tw(t, z))          ;
EQEq18(t, z)$(t_equ(t) and z_axis(z))    ..      (Tw(t, z) - Tw(t-1, z)) / Dt_equ   =e= EB_1 * (Tw(t, z+1) + Tw(t, z-1) - 2 * Tw(t, z)) / Dz / Dz + EB_2 * (Tc(t, z) - Tw(t, z)) - EB_3 * (Tw(t, z) - T_amb)                                                                                                                                                                                       ;
*****************************************
*              Rinse step               *
*****************************************
REEq00(t, z)$(t_rin(t) and z_left(z))    ..      yE(t, z)                           =e= yE_rinse                                   ;
REEq01(t, z)$(t_rin(t) and z_left(z))    ..      P(t, z)                            =e= P_high                                     ;
REEq02(t, z)$(t_rin(t) and z_left(z))    ..      Tc(t, z)                           =e= T_feed                                     ;

REEq04(t, z)$(t_rin(t) and z_right(z))   ..      yE(t, z)                           =e= yE(t, z - 1)                               ;
REEq05(t, z)$(t_rin(t) and z_right(z))   ..      P(t, z)                            =e= P_end                                      ;
REEq06(t, z)$(t_rin(t) and z_right(z))   ..      u(t, z)                            =e= u(t, z - 1)                                ;
REEq07(t, z)$(t_rin(t) and z_right(z))   ..      Tc(t, z)                           =e= Tc(t, z - 1)                               ;

REEq09(t, z)$(t_rin(t) and z_axis(z))    ..      qA(t, z) - qA(t-1, z)              =e= Dt_rin * kA * (qA_Star(t, z) - qA(t, z))   ;
REEq10(t, z)$(t_rin(t) and z_axis(z))    ..      qE(t, z) - qE(t-1, z)              =e= Dt_rin * kE * (qE_Star(t, z) - qE(t, z))   ;

REEq11(t, z)$(t_rin(t) and z_axis(z))    ..      yE(t, z) - yE(t-1, z)              =e= - Dt_rin * u(t, z) * (yE(t, z) - yE(t, z-1)) / Dz - M * Tc(t, z) * ( (qE(t, z) - qE(t-1, z)) - yE(t, z) * (qE(t, z) - qE(t-1, z) + qA(t, z) - qA(t-1, z)) ) / (P(t, z) + 1e-10)                                                                                                                            ;
REEq12(t, z)$(t_rin(t) and z_axis(z))    ..      (P(t, z) - P(t-1, z)) / Dt_rin     =e= - (P(t, z) * u(t, z) - P(t, z-1) * u(t, z-1)) / Dz - M * Tc(t, z) * ( qE(t,z) - qE(t-1, z) + qA(t,z) - qA(t-1, z) ) / Dt_rin + P(t, z) * (Tc(t, z) - Tc(t-1, z)) / (1e-10 + Tc(t, z) * Dt_rin) + P(t, z) * u(t, z) * (Tc(t, z) - Tc(t, z-1)) / (1e-10 + Tc(t, z) * Dz)                                     ;
REEq13(t, z)$(t_rin(t) and z_axis(z))    ..      (P(t, z) - P(t, z-1)) / Dz         =e= - EM * u(t, z-1) - FM * u(t, z-1) * ((u(t, z-1)*u(t, z-1))**0.5) * P(t, z-1) * MW_PE / (R * Tc(t, z))                                                                                                                                                                                                      ;
REEq14(t, z)$(t_rin(t) and z_right(z))   ..      (P(t, z) - P(t, z-1)) / Dz         =e= - EM * u(t, z-1) - FM * u(t, z-1) * ((u(t, z-1)*u(t, z-1))**0.5) * P(t, z-1) * MW_PE / (R * Tc(t, z))                                                                                                                                                                                                      ;

REEq16(t, z)$(t_rin(t) and z_axis(z))    ..      EB(t, z)                           =e= Nporosity * ( rho_s * Cp_s + Cp_a * (qA(t, z) + qE(t, z)) ) + CPR * P(t, z) / Tc(t, z)                                                                                                                                                                                                                     ;
REEq17(t, z)$(t_rin(t) and z_axis(z))    ..      (Tc(t, z) - Tc(t-1, z)) * EB(t, z) =e= - CPDR * Dt_rin * P(t, z) * u(t, z) * (Tc(t, z) - Tc(t, z-1)) / Tc(t, z) + NPCC * Tc(t, z) * (qA(t, z) - qA(t-1, z) + qE(t, z) - qE(t-1, z)) + Nporosity * (deltaH_PA * qA(t, z) - deltaH_PA * qA(t-1, z) + deltaH_PE * qE(t, z) - deltaH_PE * qE(t-1, z)) - THR * Dt_rin * (Tc(t, z) - Tw(t, z))          ;
REEq18(t, z)$(t_rin(t) and z_axis(z))    ..      (Tw(t, z) - Tw(t-1, z)) / Dt_rin   =e= EB_1 * (Tw(t, z+1) + Tw(t, z-1) - 2 * Tw(t, z)) / Dz / Dz + EB_2 * (Tc(t, z) - Tw(t, z)) - EB_3 * (Tw(t, z) - T_amb)                                                                                                                                                                                       ;
*****************************************
*            desorption step            *
*****************************************
DEEq00(t, z)$(t_des(t) and z_left(z))    ..      yE(t, z)                           =e= yE(t, z + 1)                               ;
DEEq01(t, z)$(t_des(t) and z_left(z))    ..      P(t, z)                            =e= P_low                                      ;
DEEq02(t, z)$(t_des(t) and z_left(z))    ..      u(t, z)                            =e= u(t, z + 1)                                ;
DEEq03(t, z)$(t_des(t) and z_left(z))    ..      Tc(t, z)                           =e= Tc(t, z + 1)                               ;

DEEq05(t, z)$(t_des(t) and z_right(z))   ..      yE(t, z)                           =e= yE(t, z - 1)                               ;
DEEq06(t, z)$(t_des(t) and z_right(z))   ..      P(t, z)                            =e= P(t, z - 1)                                ;
DEEq07(t, z)$(t_des(t) and z_right(z))   ..      u(t, z)                            =e= 0                                          ;
DEEq08(t, z)$(t_des(t) and z_right(z))   ..      Tc(t, z)                           =e= Tc(t, z - 1)                               ;

DEEq10(t, z)$(t_des(t) and z_axis(z))    ..      qA(t, z) - qA(t-1, z)              =e= Dt_des * kA * (qA_Star(t, z) - qA(t, z))   ;
DEEq11(t, z)$(t_des(t) and z_axis(z))    ..      qE(t, z) - qE(t-1, z)              =e= Dt_des * kE * (qE_Star(t, z) - qE(t, z))   ;

DEEq12(t, z)$(t_des(t) and z_axis(z))    ..      yE(t, z) - yE(t-1, z)              =e= - Dt_des * u(t, z) * (yE(t, z+1) - yE(t, z)) / Dz - M * Tc(t, z) * ( (qE(t, z) - qE(t-1, z)) - yE(t, z) * (qE(t, z) - qE(t-1, z) + qA(t, z) - qA(t-1, z)) ) / (P(t, z) + 1e-10)                                                                                                                            ;
DEEq13(t, z)$(t_des(t) and z_axis(z))    ..      (P(t, z) - P(t-1, z)) / Dt_des     =e= - (P(t, z+1) * u(t, z+1) - P(t, z) * u(t, z)) / Dz - M * Tc(t, z) * ( qE(t,z) - qE(t-1, z) + qA(t,z) - qA(t-1, z) ) / Dt_des + P(t, z) * (Tc(t, z) - Tc(t-1, z)) / (1e-10 + Tc(t, z) * Dt_des) + P(t, z) * u(t, z) * (Tc(t, z+1) - Tc(t, z)) / (1e-10 + Tc(t, z) * Dz)                                     ;
DEEq14(t, z)$(t_des(t) and z_axis(z))    ..      (P(t, z) - P(t, z-1)) / Dz         =e= - EM * u(t, z) - FM * u(t, z) * ((u(t, z)*u(t, z))**0.5) * P(t, z) * MW_PE / (R * Tc(t, z))                                                                                                                                                                                                                ;

DEEq16(t, z)$(t_des(t) and z_axis(z))    ..      EB(t, z)                           =e= Nporosity * ( rho_s * Cp_s + Cp_a * (qA(t, z) + qE(t, z)) ) + CPR * P(t, z) / Tc(t, z)                                                                                                                                                                                                                     ;
DEEq17(t, z)$(t_des(t) and z_axis(z))    ..      (Tc(t, z) - Tc(t-1, z)) * EB(t, z) =e= - CPDR * Dt_des * P(t, z) * u(t, z) * (Tc(t, z+1) - Tc(t, z)) / Tc(t, z) + NPCC * Tc(t, z) * (qA(t, z) - qA(t-1, z) + qE(t, z) - qE(t-1, z)) + Nporosity * (deltaH_PA * qA(t, z) - deltaH_PA * qA(t-1, z) + deltaH_PE * qE(t, z) - deltaH_PE * qE(t-1, z)) - THR * Dt_des * (Tc(t, z) - Tw(t, z))          ;
DEEq18(t, z)$(t_des(t) and z_axis(z))    ..      (Tw(t, z) - Tw(t-1, z)) / Dt_des   =e= EB_1 * (Tw(t, z+1) + Tw(t, z-1) - 2 * Tw(t, z)) / Dz / Dz + EB_2 * (Tc(t, z) - Tw(t, z)) - EB_3 * (Tw(t, z) - T_amb)                                                                                                                                                                                       ;
*****************************************
*       Equalization step two           *
*****************************************
EAEq00(t, z)$(t_eqn(t) and z_left(z))    ..      yE(t, z)                           =e= yE_equ_out                                 ;
EAEq01(t, z)$(t_eqn(t) and z_left(z))    ..      P(t, z)                            =e= P_end                                      ;
EAEq02(t, z)$(t_eqn(t) and z_left(z))    ..      Tc(t, z)                           =e= T_feed                                     ;

EAEq04(t, z)$(t_eqn(t) and z_right(z))   ..      yE(t, z)                           =e= yE(t, z - 1)                               ;
EAEq05(t, z)$(t_eqn(t) and z_right(z))   ..      P(t, z)                            =e= P(t, z - 1)                                ;
EAEq06(t, z)$(t_eqn(t) and z_right(z))   ..      u(t, z)                            =e= u(t, z - 1)                                ;
EAEq07(t, z)$(t_eqn(t) and z_right(z))   ..      u(t, z)                            =e= 0                                          ;
EAEq08(t, z)$(t_eqn(t) and z_right(z))   ..      Tc(t, z)                           =e= Tc(t, z - 1)                               ;

EAEq10(t, z)$(t_eqn(t) and z_axis(z))    ..      qA(t, z) - qA(t-1, z)              =e= Dt_eqn * kA * (qA_Star(t, z) - qA(t, z))   ;
EAEq11(t, z)$(t_eqn(t) and z_axis(z))    ..      qE(t, z) - qE(t-1, z)              =e= Dt_eqn * kE * (qE_Star(t, z) - qE(t, z))   ;

EAEq12(t, z)$(t_eqn(t) and z_axis(z))    ..      yE(t, z) - yE(t-1, z)              =e= - Dt_eqn * u(t, z) * (yE(t, z) - yE(t, z-1)) / Dz - M * Tc(t, z) * ( (qE(t, z) - qE(t-1, z)) - yE(t, z) * (qE(t, z) - qE(t-1, z) + qA(t, z) - qA(t-1, z)) ) / (P(t, z) + 1e-10)                                                                                                                            ;
EAEq13(t, z)$(t_eqn(t) and z_axis(z))    ..      (P(t, z) - P(t-1, z)) / Dt_eqn     =e= - (P(t, z) * u(t, z) - P(t, z-1) * u(t, z-1)) / Dz - M * Tc(t, z) * ( qE(t,z) - qE(t-1, z) + qA(t,z) - qA(t-1, z) ) / Dt_eqn + P(t, z) * (Tc(t, z) - Tc(t-1, z)) / (1e-10 + Tc(t, z) * Dt_eqn) + P(t, z) * u(t, z) * (Tc(t, z) - Tc(t, z-1)) / (1e-10 + Tc(t, z) * Dz)                                     ;
EAEq14(t, z)$(t_eqn(t) and z_axis(z))    ..      (P(t, z) - P(t, z-1)) / Dz         =e= - EM * u(t, z-1) - FM * u(t, z-1) * ((u(t, z-1)*u(t, z-1))**0.5) * P(t, z-1) * MW_PE / (R * Tc(t, z))                                                                                                                                                                                                      ;

EAEq16(t, z)$(t_eqn(t) and z_axis(z))    ..      EB(t, z)                           =e= Nporosity * ( rho_s * Cp_s + Cp_a * (qA(t, z) + qE(t, z)) ) + CPR * P(t, z) / Tc(t, z)                                                                                                                                                                                                                     ;
EAEq17(t, z)$(t_eqn(t) and z_axis(z))    ..      (Tc(t, z) - Tc(t-1, z)) * EB(t, z) =e= - CPDR * Dt_eqn * P(t, z) * u(t, z) * (Tc(t, z) - Tc(t, z-1)) / Tc(t, z) + NPCC * Tc(t, z) * (qA(t, z) - qA(t-1, z) + qE(t, z) - qE(t-1, z)) + Nporosity * (deltaH_PA * qA(t, z) - deltaH_PA * qA(t-1, z) + deltaH_PE * qE(t, z) - deltaH_PE * qE(t-1, z)) - THR * Dt_eqn * (Tc(t, z) - Tw(t, z))          ;
EAEq18(t, z)$(t_eqn(t) and z_axis(z))    ..      (Tw(t, z) - Tw(t-1, z)) / Dt_eqn   =e= EB_1 * (Tw(t, z+1) + Tw(t, z-1) - 2 * Tw(t, z)) / Dz / Dz + EB_2 * (Tc(t, z) - Tw(t, z)) - EB_3 * (Tw(t, z) - T_amb)                                                                                                                                                                                       ;
*****************************************
*             Rinse reflux              *
*****************************************
RREq00(t, z)$(t_ref(t) and z_left(z))    ..      yE(t, z)                           =e= yE_rin_out                                 ;
RREq01(t, z)$(t_ref(t) and z_left(z))    ..      P(t, z)                            =e= P_end                                      ;
RREq02(t, z)$(t_ref(t) and z_left(z))    ..      Tc(t, z)                           =e= T_feed                                     ;

RREq04(t, z)$(t_ref(t) and z_right(z))   ..      yE(t, z)                           =e= yE(t, z - 1)                               ;
RREq05(t, z)$(t_ref(t) and z_right(z))   ..      P(t, z)                            =e= P(t, z - 1)                                ;
RREq06(t, z)$(t_ref(t) and z_right(z))   ..      u(t, z)                            =e= u(t, z - 1)                                ;
RREq07(t, z)$(t_ref(t) and z_right(z))   ..      u(t, z)                            =e= 0                                          ;
RREq08(t, z)$(t_ref(t) and z_right(z))   ..      Tc(t, z)                           =e= Tc(t, z - 1)                               ;

RREq10(t, z)$(t_ref(t) and z_axis(z))    ..      qA(t, z) - qA(t-1, z)              =e= Dt_ref * kA * (qA_Star(t, z) - qA(t, z))   ;
RREq11(t, z)$(t_ref(t) and z_axis(z))    ..      qE(t, z) - qE(t-1, z)              =e= Dt_ref * kE * (qE_Star(t, z) - qE(t, z))   ;

RREq12(t, z)$(t_ref(t) and z_axis(z))    ..      yE(t, z) - yE(t-1, z)              =e= - Dt_ref * u(t, z) * (yE(t, z) - yE(t, z-1)) / Dz - M * Tc(t, z) * ( (qE(t, z) - qE(t-1, z)) - yE(t, z) * (qE(t, z) - qE(t-1, z) + qA(t, z) - qA(t-1, z)) ) / (P(t, z) + 1e-10)                                                                                                                            ;
RREq13(t, z)$(t_ref(t) and z_axis(z))    ..      (P(t, z) - P(t-1, z)) / Dt_ref     =e= - (P(t, z) * u(t, z) - P(t, z-1) * u(t, z-1)) / Dz - M * Tc(t, z) * ( qE(t,z) - qE(t-1, z) + qA(t,z) - qA(t-1, z) ) / Dt_ref + P(t, z) * (Tc(t, z) - Tc(t-1, z)) / (1e-10 + Tc(t, z) * Dt_ref) + P(t, z) * u(t, z) * (Tc(t, z) - Tc(t, z-1)) / (1e-10 + Tc(t, z) * Dz)                                     ;
RREq14(t, z)$(t_ref(t) and z_axis(z))    ..      (P(t, z) - P(t, z-1)) / Dz         =e= - EM * u(t, z-1) - FM * u(t, z-1) * ((u(t, z-1)*u(t, z-1))**0.5) * P(t, z-1) * MW_PE / (R * Tc(t, z))                                                                                                                                                                                                      ;

RREq16(t, z)$(t_ref(t) and z_axis(z))    ..      EB(t, z)                           =e= Nporosity * ( rho_s * Cp_s + Cp_a * (qA(t, z) + qE(t, z)) ) + CPR * P(t, z) / Tc(t, z)                                                                                                                                                                                                                     ;
RREq17(t, z)$(t_ref(t) and z_axis(z))    ..      (Tc(t, z) - Tc(t-1, z)) * EB(t, z) =e= - CPDR * Dt_ref * P(t, z) * u(t, z) * (Tc(t, z) - Tc(t, z-1)) / Tc(t, z) + NPCC * Tc(t, z) * (qA(t, z) - qA(t-1, z) + qE(t, z) - qE(t-1, z)) + Nporosity * (deltaH_PA * qA(t, z) - deltaH_PA * qA(t-1, z) + deltaH_PE * qE(t, z) - deltaH_PE * qE(t-1, z)) - THR * Dt_ref * (Tc(t, z) - Tw(t, z))          ;
RREq18(t, z)$(t_ref(t) and z_axis(z))    ..      (Tw(t, z) - Tw(t-1, z)) / Dt_ref   =e= EB_1 * (Tw(t, z+1) + Tw(t, z-1) - 2 * Tw(t, z)) / Dz / Dz + EB_2 * (Tc(t, z) - Tw(t, z)) - EB_3 * (Tw(t, z) - T_amb)                                                                                                                                                                                       ;
*****************************************
*       cyclic steady state: CSS        *
*****************************************
CSEq00(t, z)$(t_start(t) and z_axis(z))  ..      yE(t, z) - yE(t+t_total, z)      =e= 0                                          ;
CSEq01(t, z)$(t_start(t) and z_axis(z))  ..      P(t, z)  - P(t+t_total, z)       =e= 0                                          ;
CSEq02(t, z)$(t_start(t) and z_axis(z))  ..      qE(t, z) - qE(t+t_total, z)      =e= 0                                          ;
CSEq03(t, z)$(t_start(t) and z_axis(z))  ..      qA(t, z) - qA(t+t_total, z)      =e= 0                                          ;
CSEq04(t, z)$(t_start(t) and z_axis(z))  ..      Tc(t, z) - Tc(t+t_total, z)      =e= 0                                          ;
CSEq05(t, z)$(t_start(t) and z_axis(z))  ..      Tw(t, z) - Tw(t+t_total, z)      =e= 0                                          ;
*****************************************
*              connections              *
*****************************************
CCEq00(z)$z_right(z)                     ..      Equ_out         =e= P_end * sum(t_equ, u(t_equ, z) / Tc(t_equ, z)) * Dt_equ * CPRT                                              ;
CCEq01(z)$z_right(z)                     ..      yE_equ_out      =e= sum(t_equ, u(t_equ, z) * yE(t_equ, z) / Tc(t_equ, z)) / (1e-10 + sum(t_equ, u(t_equ, z) / Tc(t_equ, z)))    ;
CCEq02(z)$z_left(z)                      ..      Eqn_feed        =e= P_end * sum(t_eqn, u(t_eqn, z) / Tc(t_eqn, z)) * Dt_eqn * CPRT                                              ;
CCEq03                                   ..      Eqn_feed        =e= Equ_out                                                                                                     ;
CCEq04(z)$z_right(z)                     ..      Rinse_out       =e= P_end * sum(t_rin, u(t_rin, z) / Tc(t_rin, z)) * Dt_rin * CPRT                                              ;
CCEq05(z)$z_right(z)                     ..      yE_rin_out      =e= sum(t_rin, u(t_rin, z) * yE(t_rin, z) / Tc(t_rin, z)) / (1e-10 + sum(t_rin, u(t_rin, z) / Tc(t_rin, z)))    ;
CCEq06(z)$z_left(z)                      ..      Reflux_feed     =e= P_end * sum(t_ref, u(t_ref, z) / Tc(t_ref, z)) * Dt_ref * CPRT                                              ;
CCEq07                                   ..      Rinse_out       =e= Reflux_feed                                                                                                 ;
*****************************************
*              mass balances            *
*****************************************
MBEq00(z)$z_left(z)                      ..      PA_pre_feed     =e= P_high * sum(t_pre, u(t_pre, z)) * yA_feed * Dt_pre * CPRT / T_feed                                         ;
MBEq01(z)$z_left(z)                      ..      PA_ads_feed     =e= P_high * sum(t_ads, u(t_ads, z)) * yA_feed * Dt_ads * CPRT / T_feed                                         ;
MBEq02(z)$z_right(z)                     ..      PA_ads_out      =e= sum(t_ads, P(t_ads, z) * u(t_ads, z) * yA(t_ads, z) / Tc(t_ads, z)) * Dt_ads * CPRT                         ;
MBEq03(z)$z_left(z)                      ..      PA_rin_feed     =e= P_high * sum(t_rin, u(t_rin, z)) * yA_rinse * Dt_rin * CPRT / T_feed                                        ;
MBEq04(z)$z_left(z)                      ..      PA_des_exhaust  =e= - sum(t_des, P(t_des, z) * u(t_des, z) * yA(t_des, z) / Tc(t_des, z)) * Dt_des * CPRT                       ;
MBEq05(z)$z_left(z)                      ..      PE_pre_feed     =e= P_high * sum(t_pre, u(t_pre, z)) * yE_feed * Dt_pre * CPRT / T_feed                                         ;
MBEq06(z)$z_left(z)                      ..      PE_ads_feed     =e= P_high * sum(t_ads, u(t_ads, z)) * yE_feed * Dt_ads * CPRT / T_feed                                         ;
MBEq07(z)$z_right(z)                     ..      PE_ads_out      =e= sum(t_ads, P(t_ads, z) * u(t_ads, z) * yE(t_ads, z) / Tc(t_ads, z)) * Dt_ads * CPRT                         ;
MBEq08(z)$z_left(z)                      ..      PE_rin_feed     =e= P_high * sum(t_rin, u(t_rin, z)) * yE_rinse * Dt_rin * CPRT / T_feed                                        ;
MBEq09(z)$z_left(z)                      ..      PE_des_exhaust  =e= - sum(t_des, P(t_des, z) * u(t_des, z) * yE(t_des, z) / Tc(t_des, z)) * Dt_des * CPRT                       ;
*****************************************
*       separation specifications       *
*****************************************
SSEq00                                   ..      PurityE         =e= PE_des_exhaust / (PA_des_exhaust + PE_des_exhaust)                                                  ;
SSEq01                                   ..      PurityA         =e= PA_ads_out / (PA_ads_out + PE_ads_out + 1E-10)                                                      ;
SSEq02                                   ..      RecoveryE       =e= (PE_des_exhaust - PE_rin_feed) / (PE_pre_feed + PE_ads_feed + 1E-10)                                ;
SSEq03                                   ..      RecoveryA       =e= PA_ads_out / (PA_pre_feed + PA_ads_feed + 1E-10)                                                    ;
SSEq04                                   ..      PE_product      =e= PE_des_exhaust - PE_rin_feed                                                                        ;
SSEq05                                   ..      Mass_balanceA   =e= (PA_pre_feed + PA_ads_feed + PA_rin_feed) - (PA_ads_out + PA_des_exhaust)                           ;
SSEq06                                   ..      Mass_balanceE   =e= (PE_pre_feed + PE_ads_feed + PE_rin_feed) - (PE_ads_out + PE_des_exhaust)                           ;
SSEq07                                   ..      Productivity    =e= PE_product / (Dt_pre + Dt_ads + Dt_equ + Dt_rin + Dt_des + Dt_eqn + Dt_ref)                         ;
SSEq08                                   ..      PurityE         =g= 0.99                                                                                                ;
SSEq09                                   ..      RecoveryE       =g= 0.5                                                                                                 ;
SSEq10                                   ..      Productivity    =g= 50;
*****************************************
*           energy consumption          *
*****************************************
ECEq00                                   ..      W_12_feed       =e= 0                                                                                                   ;

ECEq01                                   ..      W_23_out        =e= 0                                                                                                   ;

ECEq02                                   ..      W_3_feed        =e= ERTI * (PA_rin_feed + PE_rin_feed) * ((P_high / P_low)**(1/IC) - 1)                                 ;
ECEq03                                   ..      W_4_out         =e= ERTI * (PA_des_exhaust + PE_des_exhaust) * ((P_atm / P_low)**(1/IC) - 1)                            ;

ECEq04                                   ..      OBJ             =e= (W_12_feed + W_23_out + W_3_feed + W_4_out) / (1e-10 + PE_product)                                  ;

***********************************************************************************

Model PSA  /all/;

option NLP = CONOPT4;
option reslim = 36000;
option domlim = 500000;

solve PSA using NLP minimizing OBJ;

option OBJ:6;
Display OBJ.l

$ontext
OBJ is MJ/mol PE
OBJ/3.6 is kWh/mol PE
OBJ/(3.6*0.04208) = OBJ/0.151488 kWh/kg PE
OBJ / 0.000151488    kWh/ton PE
$offtext
