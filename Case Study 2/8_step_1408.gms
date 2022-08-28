$ontext
referring to Simulation and optimization of VPSA system based on pseudo transient continuation method. AlChE J. 2022:e17729.

created by Xiang Zhang, zhangx@mpi-magdeburg.mpg.de, 24.06.2022
$offtext

Set

comp                     'N2 adsorbed, O2 pass, Ar pass'                         /PN, PO/

t                        'time discretization'                                   /0 * 80/
t_start(t)               'starting step'                                         /0/
t_step(t)                'time for all steps'                                    /1 * 80/
t_pre(t)                 'time discretization for pressurization'                /1 * 10/
t_ads(t)                 'time discretization for adsorption one'                /11 * 20/
t_adp(t)                 'time discretization for adsorption two'                /21 * 30/
t_equ(t)                 'time discretization for equilization one'              /31 * 40/
t_des(t)                 'time discretization for desorption one'                /41 * 50/
t_dep(t)                 'time discretization for desorption two'                /51 * 60/
t_pur(t)                 'time discretization for purge'                         /61 * 70/
t_eqn(t)                 'time discretization for equilization two'              /71 * 80/

z                        'column axis discretization'                            /0 * 11/
z_start(z)               'axis starting point'                                   /1/
z_mid(z)                 'axis middle section'                                   /2 * 9/
z_end(z)                 'axis ending point'                                     /10/
z_left(z)                'left ghost'                                            /0/
z_axis(z)                'axis point'                                            /1 * 10/
z_right(z)               'right ghost'                                           /11/
z_ar(z)                  'axis + right ghost'                                    /1 * 11/
z_sr(z)                  'second + right ghost'                                  /2 * 11/;

Parameters

u_ini(t, z)              'velocity initial value:                m/s'
P_ini(t, z)              'pressure initial value:                Pa'
yN_ini(t, z)             'PE composition initial value:          /'
yO_ini(t, z)             'PE composition initial value:          /';

$CALL 'del 8_step_PSA_1408.gdx'
$CALL GDXXRW I=8_step_PSA_1408.xlsx O=8_step_PSA_1408.gdx Index=In_index!A1

$GDXIN 8_step_PSA_1408.gdx
$LOAD u_ini, P_ini, yN_ini, yO_ini

Parameters

yN_feed                  'input N2 composition'                                  /0.78/
yO_feed                  'input O2 composition'
MW_PN                    'molecular weight of N2                 kg/mol'         /0.028/
MW_PO                    'molecular weight of O2                 kg/mol'         /0.032/
kN                       'mass transfer coefficient N2           1/s'            /1/
kO                       'mass transfer coefficient CO2          1/s'            /2.5/
deltaH_PN                'adsorption heat of N2                  J/mol'          /23430/
deltaH_PO                'adsorption heat of O2                  J/mol'          /13220/

T_feed                   'feed temperature                       K'              /298/
T_amb                    'ambient temperature                    K'              /298/
P_atm                    'ambient pressure                       Pa'             /101325/
P_2atm                   '2 * ambinent'                                          /202652/

Length                   'bed length                             m'              /2/
Diameter                 'bed diameter                           m'              /3.5/
R_inner                  'inner radius                           m'
CSA                      'cross section area                     m2'
Porosity_bed             'porosity of column bed                 /'              /0.36/
Nporosity                '1 - porosity_bed                       /'
viscosity                'viscosity of gas                       kg/m/s'         /8e-6/
Radius_SP                'solid particle radius                  m'              /8.5e-4/
rho_s                    'particle density of solid              kg/m3'          /630/
Cp_g                     'heat capacity of gases                 J/mol/K'        /30/
Cp_s                     'heat capacity of solid                 J/kg/K'         /1210/
Cp_a                     'heat capacity of air                   J/mol/K'        /28.85/
h_W_in                   'heat transfer coefficient              W/m2/K'         /0.3/

pi                       'constant'                                              /3.1416/
R                        'gas constant'                                          /8.3145/
Nubz                     'number of z interval+1'
Dz                       'delta z: m'
DzS                      'square of delta z'
t_total                  'total duration'
t_adp_pur                'interval between adsorption two and purge'             /40/
t_equ_eqn                'equ and eqn'                                           /40/

EM                       'parameter Ergun equation               kg/m3/s '
FM                       'parameter Ergun equation               1/m     '
M                        'parameter in mass balance                      '
CPR                      'energy balance parameter'
CPDR                     'energy balance parameter'
THR                      'parameter energy balance 2 * h_W_in / R_inner'
NPCC                     'parameter energy balance Nporosity * (Cp_g - Cp_a) * Tc(t, z)'
CPRT                     'CSA * Porosity_bed / R';
;

yO_feed = 1 - yN_feed;
Nporosity = 1 - Porosity_bed;
R_inner = Diameter / 2;
CSA = pi * Diameter * Diameter / 4;

Nubz = card(z_axis)+1;
Dz = Length / card(z_axis);
DzS = Dz * Dz;
t_total = card(t) - 1;

EM = 150 * viscosity * (1 - Porosity_bed) * (1 - Porosity_bed) / (4 * Radius_SP * Radius_SP * Porosity_bed * Porosity_bed * Porosity_bed);
FM = 1.75 * (1 - Porosity_bed) / (2 * Radius_SP * Porosity_bed * Porosity_bed * Porosity_bed);
M = (1 - Porosity_bed) * R * rho_s / Porosity_bed;
CPR = Cp_g * Porosity_bed / R;
CPDR = Cp_g * Porosity_bed / Dz / R;
THR = 2 * h_W_in / R_inner;
NPCC = Nporosity * (Cp_g - Cp_a);
CPRT = CSA * Porosity_bed / R;

Parameters

Q_PN                     'N2 isotherm top log                    mol/kg/Pa'      /-9.5518/
b_PN                     'N2 isotherm top                        K'              /2910/
Q_PO                     'O2 isotherm top log                    mol/kg/Pa'      /-7.2845/
b_PO                     'O2 isotherm top                        K'              /1567/

Q_PN_b                   'N2 isotherm bot log                    mol/kg/Pa'      /-5.9666/
b_PN_b                   'N2 isotherm bot                        K'              /1612/
Q_PO_b                   'O2 isotherm bot log                    mol/kg/Pa'      /-5.3763/
b_PO_b                   'O2 isotherm bot                        K'              /441.3/
;

Positive Variables

Dt_pre                   'time interval at pressurization        s'
Dt_ads                   'time interval at adsorption one        s'
Dt_adp                   'time interval at adsorption two        s'
Dt_equ                   'time interval at equilization one      s'
Dt_des                   'time interval at desorption one        s'
Dt_dep                   'time interval at desorption two        s'
Dt_pur                   'time interval at purge                 s'
Dt_eqn                   'time interval at equilization two      s'
P_high                   'adsorption pressure                    Pa'
P_end                    'adsorption output pressure             Pa'
P_end_AP                 'adsorption output pressure             Pa'
P_end_PU                 'adsorption output pressure             Pa'
P_end_EQ                 'adsorption output pressure             Pa'
P_end_EA                 'adsorption output pressure             Pa'
P_low                    'desorption pressure                    Pa'

yN(t, z)                 'gas molar fraction of N2               /'
qN(t, z)                 'solid loading of N2                    mol/kg'
b_1_PN(t, z)             'variable in isotherm'
b_2_PN(t, z)             'variable in isotherm'
qN_Star(t, z)            'equilibrium loading of N2              mol/kg'
yO(t, z)                 'gas molar fraction of O2               /'
qO(t, z)                 'solid loading of O2                    mol/kg'
b_1_PO(t, z)             'variable in isotherm'
b_2_PO(t, z)             'variable in isotherm'
qO_Star(t, z)            'equilibrium loading of O2              mol/kg'

P(t, z)                  'pressure                               Pa'
EB(t, z)                 'denominator                            '
Tc(t, z)                 'temperature                            K'

Adp_out                  'adsorption two output                  mol/s'
yN_adp_out               'adsorption two output of N2 fraction   /'
yO_adp_out               'adsorption two output of O2 fraction   /'
Pur_feed                 'purge feed                             mol/s'
Equ_out                  'equalization one output                mol/s'
yN_equ_out               'equalization output of N2 fraction     /'
yO_equ_out               'equalization output of O2 fraction     /'
Eqn_feed                 'equalization two feed                  mol/s'

PO_ads_out               'adsorption one output of O2            mol/s'
PO_adp_out               'adsorption two output of O2            mol/s'
PO_pur_feed              'purge feed of O2                       mol/s'
P_ads_out                'total adsorption one output            mol/s'
P_adp_out                'total adsorption two output            mol/s';

Variables

u(t, z)                  'gas velocity in the column             m/s'

PurityO                  'purity of O2                           /'
ProductivityO            'productivity of O2                     mol/s';

Dt_pre.lo = 0.1;
Dt_pre.up = 1;
Dt_pre.l = 0.5;

Dt_ads.lo = 0.1;
Dt_ads.up = 1;
Dt_ads.l = 0.8;

Dt_adp.lo = 0.3;
Dt_adp.up = 1;
Dt_adp.l = 0.3;

Dt_equ.lo = 0.1;
Dt_equ.up = 1;
Dt_equ.l = 0.1;


Dt_des.lo = 0.1;
Dt_des.up = 1;
Dt_des.l = 0.5;

Dt_dep.lo = 0.1;
Dt_dep.up = 1;
Dt_dep.l = 0.8;

Dt_pur.lo = 0.3;
Dt_pur.up = 1;
Dt_pur.l = 0.3;

Dt_eqn.lo = 0.1;
Dt_eqn.up = 1;
Dt_eqn.l = 0.1;

P_high.lo = 1.2e5;
P_high.up = 1.5e5;
P_high.l = 1.3e5;

P_end.lo = 1e4;
P_end.up = 1.5e5;
P_end.l = 1.2e5;

P_end_AP.l = 1e5;
P_end_PU.l = 1e5;
P_end_EQ.l = 1e5;
P_end_EA.l = 1e5;

P_low.lo = 3.5e4;
P_low.up = 8e4;
P_low.l = 3.5e4;

P.lo(t, z) = 1e4;
P.up(t, z) = 2e5;
P.l(t, z) = P_ini(t, z);

yN.lo(t, z) = 0;
yN.up(t, z) = 1;
yN.l(t, z) = yN_ini(t, z);

qN.lo(t, z) = 0;
qN.up(t, z) = 3;
qN.l(t, z) = 0.5;

b_1_PN.lo(t, z) = 0.5;
b_1_PN.up(t, z) = 2.5;
b_1_PN.l(t, z) = 1.2;

b_2_PN.lo(t, z) = 0.4;
b_2_PN.up(t, z) = 0.9;
b_2_PN.l(t, z) = 0.57;

qN_Star.lo(t, z) = 0;
qN_Star.up(t, z) = 3;
qN_Star.l(t, z) = 1;

yO.lo(t, z) = 0;
yO.up(t, z) = 1;
yO.l(t, z) = yO_ini(t, z);

qO.lo(t, z) = 0;
qO.up(t, z) = 0.5;
qO.l(t, z) = 0.05;

b_1_PO.lo(t, z) = 0.05;
b_1_PO.up(t, z) = 0.3;
b_1_PO.l(t, z) = 0.12;

b_2_PO.lo(t, z) = 0.01;
b_2_PO.up(t, z) = 0.03;
b_2_PO.l(t, z) = 0.02;

qO_Star.lo(t, z) = 0;
qO_Star.up(t, z) = 0.5;
qO_Star.l(t, z) = 0.2;

EB.lo(t, z) = 4.6e5;
EB.up(t, z) = 5e5;
EB.l(t, z) = 4.88e5;

Tc.lo(t, z) = 293;
Tc.up(t, z) = 303;
Tc.l(t, z) = 298;

Adp_out.l = 5;
Pur_feed.l = 5;
Equ_out.l = 5;
Eqn_feed.l = 5;
PO_ads_out.l = 5;
PO_adp_out.l = 5;
PO_pur_feed.l = 5;
P_ads_out.l = 5;
P_adp_out.l = 5;

PurityO.l = 0.8;

u.lo(t, z) = -10;
u.up(t, z) = 10;
u.l(t, z) = u_ini(t, z);

ProductivityO.l = 100;


Equations

LGEq00, LGEq01, LGEq02, LGEq03, LGEq04, LGEq05, LGEq06, LGEq07, LGEq08, LGEq09, LGEq10, LGEq11, LGEq12, LGEq13, LGEq14
ISEq00, ISEq01, ISEq02, ISEq03, ISEq04, ISEq05
PREq00, PREq01, PREq02, PREq03, PREq04, PREq05, PREq06, PREq07, PREq08, PREq09, PREq10, PREq11, PREq12, PREq13
ADEq00, ADEq01, ADEq02, ADEq03, ADEq04, ADEq05, ADEq06, ADEq07, ADEq08, ADEq09, ADEq10, ADEq11, ADEq12
APEq00, APEq01, APEq02, APEq03, APEq04, APEq05, APEq06, APEq07, APEq08, APEq09, APEq10, APEq11, APEq12
EQEq00, EQEq01, EQEq02, EQEq03, EQEq04, EQEq05, EQEq06, EQEq07, EQEq08, EQEq09, EQEq10, EQEq11, EQEq12
DEEq00, DEEq01, DEEq02, DEEq03, DEEq04, DEEq05, DEEq06, DEEq07, DEEq08, DEEq09, DEEq10, DEEq11, DEEq12, DEEq13
DPEq00, DPEq01, DPEq02, DPEq03, DPEq04, DPEq05, DPEq06, DPEq07, DPEq08, DPEq09, DPEq10, DPEq11, DPEq12, DPEq13
PUEq00, PUEq01, PUEq02, PUEq03, PUEq04, PUEq05, PUEq06, PUEq07, PUEq08, PUEq09, PUEq10, PUEq11, PUEq12
EAEq00, EAEq01, EAEq02, EAEq03, EAEq04, EAEq05, EAEq06, EAEq07, EAEq08, EAEq09, EAEq10, EAEq11, EAEq12, EAEq13
CSEq00, CSEq01, CSEq02, CSEq03, CSEq04
CCEq00, CCEq01, CCEq02, CCEq03, CCEq04, CCEq05
MBEq00, MBEq01, MBEq02, MBEq03, MBEq04
SSEq00, SSEq01, SSEq02
;

*****************************************
*          logical constraint           *
*****************************************
LGEq00(t, z)                             ..      yO(t, z)                         =e= 1 - yN(t, z)                                                                                       ;
LGEq01(t, z)$z_axis(z)                   ..      EB(t, z)                         =e= Nporosity * ( rho_s * Cp_s + Cp_a * (qN(t, z) + qO(t, z)) ) + CPR * P(t, z) / Tc(t, z)             ;
LGEq02                                   ..      P_high                           =g= 1.05 * P_end                                                                                       ;
LGEq03                                   ..      P_end                            =g= 1.05 * P_low                                                                                       ;
LGEq04                                   ..      Dt_pre                           =e= Dt_des                                                                                             ;
LGEq05                                   ..      Dt_ads                           =e= Dt_dep                                                                                             ;
LGEq06                                   ..      Dt_adp                           =e= Dt_pur                                                                                             ;
LGEq07                                   ..      Dt_equ                           =e= Dt_eqn                                                                                             ;
LGEq08                                   ..      P_end                            =g= P_end_AP                                                                                       ;
LGEq09                                   ..      P_end_AP                         =g= 1.05 * P_end_PU                                                                                       ;
LGEq10                                   ..      P_end_EQ                         =g= 1.05 * P_end_EA                                                                                       ;
LGEq11(t, z)$(t_adp(t) and z_right(z))   ..      u(t, z)                          =g= 0;
LGEq12(t, z)$(t_equ(t) and z_right(z))   ..      u(t, z)                          =g= 0;
LGEq13(t, z)$(t_pur(t) and z_right(z))   ..      u(t, z)                          =l= 0;
LGEq14(t, z)$(t_eqn(t) and z_right(z))   ..      u(t, z)                          =l= 0;
*****************************************
*                isotherm               *
*****************************************
ISEq00(t, z)$z_axis(z)                   ..      b_1_PN(t, z)                     =e= exp(Q_PN + b_PN / Tc(t, z))                ;
ISEq01(t, z)$z_axis(z)                   ..      b_1_PO(t, z)                     =e= exp(Q_PO + b_PO / Tc(t, z))                ;
ISEq02(t, z)$z_axis(z)                   ..      b_2_PN(t, z)                     =e= exp(Q_PN_b + b_PN_b / Tc(t, z))            ;
ISEq03(t, z)$z_axis(z)                   ..      b_2_PO(t, z)                     =e= exp(Q_PO_b + b_PO_b / Tc(t, z))            ;
ISEq04(t, z)$z_axis(z)                   ..      qN_Star(t, z)                    =e= b_1_PN(t, z) * P(t, z) * yN(t, z) / (1e5 + b_2_PN(t, z) * P(t, z) * yN(t, z) + b_2_PO(t, z) * P(t, z) * yO(t, z))  ;
ISEq05(t, z)$z_axis(z)                   ..      qO_Star(t, z)                    =e= b_1_PO(t, z) * P(t, z) * yO(t, z) / (1e5 + b_2_PN(t, z) * P(t, z) * yN(t, z) + b_2_PO(t, z) * P(t, z) * yO(t, z))  ;
*****************************************
*         Pressurization step           *
*****************************************
PREq00(t, z)$(t_pre(t) and z_left(z))    ..      yN(t, z)                         =e= yN_feed                                    ;
PREq01(t, z)$(t_pre(t) and z_left(z))    ..      P(t, z)                          =e= P_high                                     ;
PREq02(t, z)$(t_pre(t) and z_left(z))    ..      Tc(t, z)                         =e= T_feed                                     ;

PREq03(t, z)$(t_pre(t) and z_right(z))   ..      yN(t, z)                         =e= yN(t, z - 1)                               ;
PREq04(t, z)$(t_pre(t) and z_right(z))   ..      P(t, z)                          =e= P(t, z - 1)                                ;
PREq05(t, z)$(t_pre(t) and z_right(z))   ..      u(t, z)                          =e= u(t, z - 1)                                ;
PREq06(t, z)$(t_pre(t) and z_right(z))   ..      u(t, z)                          =e= 0                                          ;
PREq07(t, z)$(t_pre(t) and z_right(z))   ..      Tc(t, z)                         =e= Tc(t, z - 1)                               ;

PREq08(t, z)$(t_pre(t) and z_axis(z))    ..      qN(t, z) - qN(t-1, z)            =e= Dt_pre * kN * (qN_Star(t, z) - qN(t, z))   ;
PREq09(t, z)$(t_pre(t) and z_axis(z))    ..      qO(t, z) - qO(t-1, z)            =e= Dt_pre * kO * (qO_Star(t, z) - qO(t, z))   ;

PREq10(t, z)$(t_pre(t) and z_axis(z))    ..      yN(t, z) - yN(t-1, z)            =e= - Dt_pre * u(t, z) * (yN(t, z) - yN(t, z-1)) / Dz - M * Tc(t, z) * ( (qN(t, z) - qN(t-1, z)) - yN(t, z) * (qN(t, z) - qN(t-1, z) + qO(t, z) - qO(t-1, z)) ) / P(t, z)                                                                                                                              ;
PREq11(t, z)$(t_pre(t) and z_axis(z))    ..      (P(t, z) - P(t-1, z)) / Dt_pre   =e= - (P(t, z) * u(t, z) - P(t, z-1) * u(t, z-1)) / Dz - M * Tc(t, z) * ( qN(t, z) - qN(t-1, z) + qO(t, z) - qO(t-1, z) ) / Dt_pre + P(t, z) * (Tc(t, z) - Tc(t-1, z)) / (Tc(t, z) * Dt_pre) + P(t, z) * u(t, z) * (Tc(t, z) - Tc(t, z-1)) / (Tc(t, z) * Dz)                                           ;
PREq12(t, z)$(t_pre(t) and z_axis(z))    ..      (P(t, z) - P(t, z-1)) / Dz       =e= - EM * u(t, z-1) - FM * u(t, z-1) * ((u(t, z-1)*u(t, z-1) + 1e-10)**0.5) * P(t, z-1) * MW_PN / (R * Tc(t, z))                                                                                                                                                                                      ;

PREq13(t, z)$(t_pre(t) and z_axis(z))    ..      (Tc(t, z) - Tc(t-1, z)) * EB(t, z) =e= - CPDR * Dt_pre * P(t, z) * u(t, z) * (Tc(t, z) - Tc(t, z-1)) / Tc(t, z) + NPCC * Tc(t, z) * (qN(t, z) - qN(t-1, z) + qO(t, z) - qO(t-1, z)) + Nporosity * (deltaH_PN * qN(t, z) - deltaH_PN * qN(t-1, z) + deltaH_PO * qO(t, z) - deltaH_PO * qO(t-1, z)) - THR * Dt_pre * (Tc(t, z) - T_amb)   ;
*****************************************
*           adsorption step 1           *
*****************************************
ADEq00(t, z)$(t_ads(t) and z_left(z))    ..      yN(t, z)                         =e= yN_feed                                    ;
ADEq01(t, z)$(t_ads(t) and z_left(z))    ..      P(t, z)                          =e= P_high                                     ;
ADEq02(t, z)$(t_ads(t) and z_left(z))    ..      Tc(t, z)                         =e= T_feed                                     ;

ADEq03(t, z)$(t_ads(t) and z_right(z))   ..      yN(t, z)                         =e= yN(t, z - 1)                               ;
ADEq04(t, z)$(t_ads(t) and z_right(z))   ..      P(t, z)                          =e= P_end                                      ;
ADEq05(t, z)$(t_ads(t) and z_right(z))   ..      u(t, z)                          =e= u(t, z - 1)                                ;
ADEq06(t, z)$(t_ads(t) and z_right(z))   ..      Tc(t, z)                         =e= Tc(t, z - 1)                               ;

ADEq07(t, z)$(t_ads(t) and z_axis(z))    ..      qN(t, z) - qN(t-1, z)            =e= Dt_ads * kN * (qN_Star(t, z) - qN(t, z))   ;
ADEq08(t, z)$(t_ads(t) and z_axis(z))    ..      qO(t, z) - qO(t-1, z)            =e= Dt_ads * kO * (qO_Star(t, z) - qO(t, z))   ;

ADEq09(t, z)$(t_ads(t) and z_axis(z))    ..      yN(t, z) - yN(t-1, z)            =e= - Dt_ads * u(t, z) * (yN(t, z) - yN(t, z-1)) / Dz - M * Tc(t, z) * ( (qN(t, z) - qN(t-1, z)) - yN(t, z) * (qN(t, z) - qN(t-1, z) + qO(t, z) - qO(t-1, z) ) ) / P(t, z)                                                                                                                             ;
ADEq10(t, z)$(t_ads(t) and z_axis(z))    ..      (P(t, z) - P(t-1, z)) / Dt_ads   =e= - (P(t, z) * u(t, z) - P(t, z-1) * u(t, z-1)) / Dz - M * Tc(t, z) * ( qN(t, z) - qN(t-1, z) + qO(t, z) - qO(t-1, z) ) / Dt_ads + P(t, z) * (Tc(t, z) - Tc(t-1, z)) / (Tc(t, z) * Dt_ads) + P(t, z) * u(t, z) * (Tc(t, z) - Tc(t, z-1)) / (Tc(t, z) * Dz)                                           ;
ADEq11(t, z)$(t_ads(t) and z_ar(z))      ..      (P(t, z) - P(t, z-1)) / Dz       =e= - EM * u(t, z-1) - FM * u(t, z-1) * ((u(t, z-1)*u(t, z-1) + 1e-10)**0.5) * P(t, z-1) * MW_PN / (R * Tc(t, z))                                                                                                                                                                                      ;

ADEq12(t, z)$(t_ads(t) and z_axis(z))    ..      (Tc(t, z) - Tc(t-1, z)) * EB(t, z) =e= - CPDR * Dt_ads * P(t, z) * u(t, z) * (Tc(t, z) - Tc(t, z-1)) / Tc(t, z) + NPCC * Tc(t, z) * (qN(t, z) - qN(t-1, z) + qO(t, z) - qO(t-1, z)) + Nporosity * (deltaH_PN * qN(t, z) - deltaH_PN * qN(t-1, z) + deltaH_PO * qO(t, z) - deltaH_PO * qO(t-1, z)) - THR * Dt_ads * (Tc(t, z) - T_amb)   ;
*****************************************
*           adsorption step 2           *
*****************************************
APEq00(t, z)$(t_adp(t) and z_left(z))    ..      yN(t, z)                         =e= yN_feed                                    ;
APEq01(t, z)$(t_adp(t) and z_left(z))    ..      P(t, z)                          =e= P_high                                     ;
APEq02(t, z)$(t_adp(t) and z_left(z))    ..      Tc(t, z)                         =e= T_feed                                     ;

APEq03(t, z)$(t_adp(t) and z_right(z))   ..      yN(t, z)                         =e= yN(t, z - 1)                               ;
APEq04(t, z)$(t_adp(t) and z_right(z))   ..      P(t, z)                          =e= P_end_AP                                   ;
APEq05(t, z)$(t_adp(t) and z_right(z))   ..      u(t, z)                          =e= u(t, z - 1)                                ;
APEq06(t, z)$(t_adp(t) and z_right(z))   ..      Tc(t, z)                         =e= Tc(t, z - 1)                               ;

APEq07(t, z)$(t_adp(t) and z_axis(z))    ..      qN(t, z) - qN(t-1, z)            =e= Dt_adp * kN * (qN_Star(t, z) - qN(t, z))   ;
APEq08(t, z)$(t_adp(t) and z_axis(z))    ..      qO(t, z) - qO(t-1, z)            =e= Dt_adp * kO * (qO_Star(t, z) - qO(t, z))   ;

APEq09(t, z)$(t_adp(t) and z_axis(z))    ..      yN(t, z) - yN(t-1, z)            =e= - Dt_adp * u(t, z) * (yN(t, z) - yN(t, z-1)) / Dz - M * Tc(t, z) * ( (qN(t, z) - qN(t-1, z)) - yN(t, z) * (qN(t, z) - qN(t-1, z) + qO(t, z) - qO(t-1, z) ) ) / P(t, z)                                                                                                                             ;
APEq10(t, z)$(t_adp(t) and z_axis(z))    ..      (P(t, z) - P(t-1, z)) / Dt_adp   =e= - (P(t, z) * u(t, z) - P(t, z-1) * u(t, z-1)) / Dz - M * Tc(t, z) * ( qN(t, z) - qN(t-1, z) + qO(t, z) - qO(t-1, z) ) / Dt_adp + P(t, z) * (Tc(t, z) - Tc(t-1, z)) / (Tc(t, z) * Dt_adp) + P(t, z) * u(t, z) * (Tc(t, z) - Tc(t, z-1)) / (Tc(t, z) * Dz)                                           ;
APEq11(t, z)$(t_adp(t) and z_ar(z))      ..      (P(t, z) - P(t, z-1)) / Dz       =e= - EM * u(t, z-1) - FM * u(t, z-1) * ((u(t, z-1)*u(t, z-1) + 1e-10)**0.5) * P(t, z-1) * MW_PN / (R * Tc(t, z))                                                                                                                                                                                      ;

APEq12(t, z)$(t_adp(t) and z_axis(z))    ..      (Tc(t, z) - Tc(t-1, z)) * EB(t, z) =e= - CPDR * Dt_adp * P(t, z) * u(t, z) * (Tc(t, z) - Tc(t, z-1)) / Tc(t, z) + NPCC * Tc(t, z) * (qN(t, z) - qN(t-1, z) + qO(t, z) - qO(t-1, z)) + Nporosity * (deltaH_PN * qN(t, z) - deltaH_PN * qN(t-1, z) + deltaH_PO * qO(t, z) - deltaH_PO * qO(t-1, z)) - THR * Dt_adp * (Tc(t, z) - T_amb)   ;
*****************************************
*       Equalization step one           *
*****************************************
EQEq00(t, z)$(t_equ(t) and z_left(z))    ..      yN(t, z)                         =e= yN(t, z + 1)                               ;
EQEq01(t, z)$(t_equ(t) and z_left(z))    ..      P(t, z)                          =e= P(t, z + 1)                                ;
EQEq02(t, z)$(t_equ(t) and z_left(z))    ..      Tc(t, z)                         =e= Tc(t, z + 1)                               ;

EQEq03(t, z)$(t_equ(t) and z_right(z))   ..      yN(t, z)                         =e= yN(t, z - 1)                               ;
EQEq04(t, z)$(t_equ(t) and z_right(z))   ..      P(t, z)                          =e= P_end_EQ                                   ;
EQEq05(t, z)$(t_equ(t) and z_right(z))   ..      u(t, z)                          =e= u(t, z - 1)                                ;
EQEq06(t, z)$(t_equ(t) and z_right(z))   ..      Tc(t, z)                         =e= Tc(t, z - 1)                               ;

EQEq07(t, z)$(t_equ(t) and z_axis(z))    ..      qN(t, z) - qN(t-1, z)            =e= Dt_equ * kN * (qN_Star(t, z) - qN(t, z))   ;
EQEq08(t, z)$(t_equ(t) and z_axis(z))    ..      qO(t, z) - qO(t-1, z)            =e= Dt_equ * kO * (qO_Star(t, z) - qO(t, z))   ;

EQEq09(t, z)$(t_equ(t) and z_axis(z))    ..      yN(t, z) - yN(t-1, z)            =e= - Dt_equ * u(t, z) * (yN(t, z) - yN(t, z-1)) / Dz - M * Tc(t, z) * ( (qN(t, z) - qN(t-1, z)) - yN(t, z) * (qN(t, z) - qN(t-1, z) + qO(t, z) - qO(t-1, z) ) ) / P(t, z)                                                                                                                             ;
EQEq10(t, z)$(t_equ(t) and z_axis(z))    ..      (P(t, z) - P(t-1, z)) / Dt_equ   =e= - (P(t, z) * u(t, z) - P(t, z-1) * u(t, z-1)) / Dz - M * Tc(t, z) * ( qN(t, z) - qN(t-1, z) + qO(t, z) - qO(t-1, z) ) / Dt_equ + P(t, z) * (Tc(t, z) - Tc(t-1, z)) / (Tc(t, z) * Dt_equ) + P(t, z) * u(t, z) * (Tc(t, z) - Tc(t, z-1)) / (Tc(t, z) * Dz)                                           ;
EQEq11(t, z)$(t_equ(t) and z_ar(z))      ..      (P(t, z) - P(t, z-1)) / Dz       =e= - EM * u(t, z-1) - FM * u(t, z-1) * ((u(t, z-1)*u(t, z-1) + 1e-10)**0.5) * P(t, z-1) * MW_PN / (R * Tc(t, z))                                                                                                                                                                                      ;

EQEq12(t, z)$(t_equ(t) and z_axis(z))    ..      (Tc(t, z) - Tc(t-1, z)) * EB(t, z) =e= - CPDR * Dt_equ * P(t, z) * u(t, z) * (Tc(t, z) - Tc(t, z-1)) / Tc(t, z) + NPCC * Tc(t, z) * (qN(t, z) - qN(t-1, z) + qO(t, z) - qO(t-1, z)) + Nporosity * (deltaH_PN * qN(t, z) - deltaH_PN * qN(t-1, z) + deltaH_PO * qO(t, z) - deltaH_PO * qO(t-1, z)) - THR * Dt_equ * (Tc(t, z) - T_amb)   ;
*****************************************
*          desorption step one          *
*****************************************
DEEq00(t, z)$(t_des(t) and z_left(z))    ..      yN(t, z)                         =e= yN(t, z + 1)                               ;
DEEq01(t, z)$(t_des(t) and z_left(z))    ..      P(t, z)                          =e= P_low                                      ;
DEEq02(t, z)$(t_des(t) and z_left(z))    ..      u(t, z)                          =e= u(t, z + 1)                                ;
DEEq03(t, z)$(t_des(t) and z_left(z))    ..      Tc(t, z)                         =e= Tc(t, z + 1)                               ;

DEEq04(t, z)$(t_des(t) and z_right(z))   ..      yN(t, z)                         =e= yN(t, z - 1)                               ;
DEEq05(t, z)$(t_des(t) and z_right(z))   ..      P(t, z)                          =e= P(t, z - 1)                                ;
DEEq06(t, z)$(t_des(t) and z_right(z))   ..      u(t, z)                          =e= 0                                          ;
DEEq07(t, z)$(t_des(t) and z_right(z))   ..      Tc(t, z)                         =e= Tc(t, z - 1)                               ;

DEEq08(t, z)$(t_des(t) and z_axis(z))    ..      qN(t, z) - qN(t-1, z)            =e= Dt_des * kN * (qN_Star(t, z) - qN(t, z))   ;
DEEq09(t, z)$(t_des(t) and z_axis(z))    ..      qO(t, z) - qO(t-1, z)            =e= Dt_des * kO * (qO_Star(t, z) - qO(t, z))   ;

DEEq10(t, z)$(t_des(t) and z_axis(z))    ..      yN(t, z) - yN(t-1, z)            =e= - Dt_des * u(t, z) * (yN(t, z+1) - yN(t, z)) / Dz - M * Tc(t, z) * ( (qN(t, z) - qN(t-1, z)) - yN(t, z) * (qN(t, z) - qN(t-1, z) + qO(t, z) - qO(t-1, z) ) ) / P(t, z)                                                                                                                             ;
DEEq11(t, z)$(t_des(t) and z_axis(z))    ..      (P(t, z) - P(t-1, z)) / Dt_des   =e= - (P(t, z+1) * u(t, z+1) - P(t, z) * u(t, z)) / Dz - M * Tc(t, z) * ( qN(t, z) - qN(t-1, z) + qO(t, z) - qO(t-1, z) ) / Dt_des + P(t, z) * (Tc(t, z) - Tc(t-1, z)) / (Tc(t, z) * Dt_des) + P(t, z) * u(t, z) * (Tc(t, z+1) - Tc(t, z)) / (Tc(t, z) * Dz)                                           ;
DEEq12(t, z)$(t_des(t) and z_axis(z))    ..      (P(t, z) - P(t, z-1)) / Dz       =e= - EM * u(t, z) - FM * u(t, z) * ((u(t, z)*u(t, z) + 1e-10)**0.5) * P(t, z) * MW_PN / (R * Tc(t, z))                                                                                                                                                                                                ;

DEEq13(t, z)$(t_des(t) and z_axis(z))    ..      (Tc(t, z) - Tc(t-1, z)) * EB(t, z) =e= - CPDR * Dt_des * P(t, z) * u(t, z) * (Tc(t, z+1) - Tc(t, z)) / Tc(t, z) + NPCC * Tc(t, z) * (qN(t, z) - qN(t-1, z) + qO(t, z) - qO(t-1, z)) + Nporosity * (deltaH_PN * qN(t, z) - deltaH_PN * qN(t-1, z) + deltaH_PO * qO(t, z) - deltaH_PO * qO(t-1, z)) - THR * Dt_des * (Tc(t, z) - T_amb)   ;
*****************************************
*          desorption step two          *
*****************************************
DPEq00(t, z)$(t_dep(t) and z_left(z))    ..      yN(t, z)                         =e= yN(t, z + 1)                               ;
DPEq01(t, z)$(t_dep(t) and z_left(z))    ..      P(t, z)                          =e= P_low                                      ;
DPEq02(t, z)$(t_dep(t) and z_left(z))    ..      u(t, z)                          =e= u(t, z + 1)                                ;
DPEq03(t, z)$(t_dep(t) and z_left(z))    ..      Tc(t, z)                         =e= Tc(t, z + 1)                               ;

DPEq04(t, z)$(t_dep(t) and z_right(z))   ..      yN(t, z)                         =e= yN(t, z - 1)                               ;
DPEq05(t, z)$(t_dep(t) and z_right(z))   ..      P(t, z)                          =e= P(t, z - 1)                                ;
DPEq06(t, z)$(t_dep(t) and z_right(z))   ..      u(t, z)                          =e= 0                                          ;
DPEq07(t, z)$(t_dep(t) and z_right(z))   ..      Tc(t, z)                         =e= Tc(t, z - 1)                               ;

DPEq08(t, z)$(t_dep(t) and z_axis(z))    ..      qN(t, z) - qN(t-1, z)            =e= Dt_dep * kN * (qN_Star(t, z) - qN(t, z))   ;
DPEq09(t, z)$(t_dep(t) and z_axis(z))    ..      qO(t, z) - qO(t-1, z)            =e= Dt_dep * kO * (qO_Star(t, z) - qO(t, z))   ;

DPEq10(t, z)$(t_dep(t) and z_axis(z))    ..      yN(t, z) - yN(t-1, z)            =e= - Dt_dep * u(t, z) * (yN(t, z+1) - yN(t, z)) / Dz - M * Tc(t, z) * ( (qN(t, z) - qN(t-1, z)) - yN(t, z) * (qN(t, z) - qN(t-1, z) + qO(t, z) - qO(t-1, z) ) ) / P(t, z)                                                                                                                             ;
DPEq11(t, z)$(t_dep(t) and z_axis(z))    ..      (P(t, z) - P(t-1, z)) / Dt_dep   =e= - (P(t, z+1) * u(t, z+1) - P(t, z) * u(t, z)) / Dz - M * Tc(t, z) * ( qN(t, z) - qN(t-1, z) + qO(t, z) - qO(t-1, z) ) / Dt_dep + P(t, z) * (Tc(t, z) - Tc(t-1, z)) / (Tc(t, z) * Dt_dep) + P(t, z) * u(t, z) * (Tc(t, z+1) - Tc(t, z)) / (Tc(t, z) * Dz)                                           ;
DPEq12(t, z)$(t_dep(t) and z_axis(z))    ..      (P(t, z) - P(t, z-1)) / Dz       =e= - EM * u(t, z) - FM * u(t, z) * ((u(t, z)*u(t, z) + 1e-10)**0.5) * P(t, z) * MW_PN / (R * Tc(t, z))                                                                                                                                                                                                ;

DPEq13(t, z)$(t_dep(t) and z_axis(z))    ..      (Tc(t, z) - Tc(t-1, z)) * EB(t, z) =e= - CPDR * Dt_dep * P(t, z) * u(t, z) * (Tc(t, z+1) - Tc(t, z)) / Tc(t, z) + NPCC * Tc(t, z) * (qN(t, z) - qN(t-1, z) + qO(t, z) - qO(t-1, z)) + Nporosity * (deltaH_PN * qN(t, z) - deltaH_PN * qN(t-1, z) + deltaH_PO * qO(t, z) - deltaH_PO * qO(t-1, z)) - THR * Dt_dep * (Tc(t, z) - T_amb)   ;
*****************************************
*             purge step                *
*****************************************
PUEq00(t, z)$(t_pur(t) and z_left(z))    ..      yN(t, z)                         =e= yN(t, z + 1)                               ;
PUEq01(t, z)$(t_pur(t) and z_left(z))    ..      P(t, z)                          =e= P_low                                      ;
PUEq02(t, z)$(t_pur(t) and z_left(z))    ..      u(t, z)                          =e= u(t, z + 1)                                ;
PUEq03(t, z)$(t_pur(t) and z_left(z))    ..      Tc(t, z)                         =e= Tc(t, z + 1)                               ;

PUEq04(t, z)$(t_pur(t) and z_right(z))   ..      yN(t, z)                         =e= yN(t - t_adp_pur, z)                       ;
PUEq05(t, z)$(t_pur(t) and z_right(z))   ..      P(t, z)                          =e= P_end_PU                                   ;
PUEq06(t, z)$(t_pur(t) and z_right(z))   ..      Tc(t, z)                         =e= Tc(t - t_adp_pur, z)                       ;

PUEq07(t, z)$(t_pur(t) and z_axis(z))    ..      qN(t, z) - qN(t-1, z)            =e= Dt_pur * kN * (qN_Star(t, z) - qN(t, z))   ;
PUEq08(t, z)$(t_pur(t) and z_axis(z))    ..      qO(t, z) - qO(t-1, z)            =e= Dt_pur * kO * (qO_Star(t, z) - qO(t, z))   ;

PUEq09(t, z)$(t_pur(t) and z_axis(z))    ..      yN(t, z) - yN(t-1, z)            =e= - Dt_pur * u(t, z) * (yN(t, z+1) - yN(t, z)) / Dz - M * Tc(t, z) * ( (qN(t, z) - qN(t-1, z)) - yN(t, z) * (qN(t, z) - qN(t-1, z) + qO(t, z) - qO(t-1, z) ) ) / P(t, z)                                                                                                                             ;

PUEq10(t, z)$(t_pur(t) and z_axis(z))    ..      (P(t, z) - P(t-1, z)) / Dt_pur   =e= - (P(t, z+1) * u(t, z+1) - P(t, z) * u(t, z)) / Dz - M * Tc(t, z) * ( qN(t, z) - qN(t-1, z) + qO(t, z) - qO(t-1, z) ) / Dt_pur + P(t, z) * (Tc(t, z) - Tc(t-1, z)) / (Tc(t, z) * Dt_pur) + P(t, z) * u(t, z) * (Tc(t, z+1) - Tc(t, z)) / (Tc(t, z) * Dz)                                           ;
PUEq11(t, z)$(t_pur(t) and z_ar(z))      ..      (P(t, z) - P(t, z-1)) / Dz       =e= - EM * u(t, z) - FM * u(t, z) * ((u(t, z)*u(t, z) + 1e-10)**0.5) * P(t, z) * MW_PN / (R * Tc(t, z))                                                                                                                                                                                                ;

PUEq12(t, z)$(t_pur(t) and z_axis(z))    ..      (Tc(t, z) - Tc(t-1, z)) * EB(t, z) =e= - CPDR * Dt_pur * P(t, z) * u(t, z) * (Tc(t, z+1) - Tc(t, z)) / Tc(t, z) + NPCC * Tc(t, z) * (qN(t, z) - qN(t-1, z) + qO(t, z) - qO(t-1, z)) + Nporosity * (deltaH_PN * qN(t, z) - deltaH_PN * qN(t-1, z) + deltaH_PO * qO(t, z) - deltaH_PO * qO(t-1, z)) - THR * Dt_pur * (Tc(t, z) - T_amb)   ;
*****************************************
*       Equalization step two           *
*****************************************
EAEq00(t, z)$(t_eqn(t) and z_left(z))    ..      yN(t, z)                         =e= yN(t, z + 1)                               ;
EAEq01(t, z)$(t_eqn(t) and z_left(z))    ..      P(t, z)                          =e= P(t, z + 1)                                ;
EAEq02(t, z)$(t_eqn(t) and z_left(z))    ..      u(t, z)                          =e= 0                                          ;
EAEq03(t, z)$(t_eqn(t) and z_left(z))    ..      u(t, z)                          =e= u(t, z + 1)                                ;
EAEq04(t, z)$(t_eqn(t) and z_left(z))    ..      Tc(t, z)                         =e= Tc(t, z + 1)                               ;

EAEq05(t, z)$(t_eqn(t) and z_right(z))   ..      yN(t, z)                         =e= yN(t - t_equ_eqn, z)                       ;
EAEq06(t, z)$(t_eqn(t) and z_right(z))   ..      P(t, z)                          =e= P_end_EA                                   ;
EAEq07(t, z)$(t_eqn(t) and z_right(z))   ..      Tc(t, z)                         =e= Tc(t - t_equ_eqn, z)                       ;

EAEq08(t, z)$(t_eqn(t) and z_axis(z))    ..      qN(t, z) - qN(t-1, z)            =e= Dt_eqn * kN * (qN_Star(t, z) - qN(t, z))   ;
EAEq09(t, z)$(t_eqn(t) and z_axis(z))    ..      qO(t, z) - qO(t-1, z)            =e= Dt_eqn * kO * (qO_Star(t, z) - qO(t, z))   ;

EAEq10(t, z)$(t_eqn(t) and z_axis(z))    ..      yN(t, z) - yN(t-1, z)            =e= - Dt_eqn * u(t, z) * (yN(t, z+1) - yN(t, z)) / Dz - M * Tc(t, z) * ( (qN(t, z) - qN(t-1, z)) - yN(t, z) * (qN(t, z) - qN(t-1, z) + qO(t, z) - qO(t-1, z) ) ) / P(t, z)                                                                                                                             ;

EAEq11(t, z)$(t_eqn(t) and z_axis(z))    ..      (P(t, z) - P(t-1, z)) / Dt_eqn   =e= - (P(t, z+1) * u(t, z+1) - P(t, z) * u(t, z)) / Dz - M * Tc(t, z) * ( qN(t, z) - qN(t-1, z) + qO(t, z) - qO(t-1, z) ) / Dt_eqn + P(t, z) * (Tc(t, z) - Tc(t-1, z)) / (Tc(t, z) * Dt_eqn) + P(t, z) * u(t, z) * (Tc(t, z+1) - Tc(t, z)) / (Tc(t, z) * Dz)                                           ;
EAEq12(t, z)$(t_eqn(t) and z_sr(z))      ..      (P(t, z) - P(t, z-1)) / Dz       =e= - EM * u(t, z) - FM * u(t, z) * ((u(t, z)*u(t, z) + 1e-10)**0.5) * P(t, z) * MW_PN / (R * Tc(t, z))                                                                                                                                                                                                ;

EAEq13(t, z)$(t_eqn(t) and z_axis(z))    ..      (Tc(t, z) - Tc(t-1, z)) * EB(t, z) =e= - CPDR * Dt_eqn * P(t, z) * u(t, z) * (Tc(t, z+1) - Tc(t, z)) / Tc(t, z) + NPCC * Tc(t, z) * (qN(t, z) - qN(t-1, z) + qO(t, z) - qO(t-1, z)) + Nporosity * (deltaH_PN * qN(t, z) - deltaH_PN * qN(t-1, z) + deltaH_PO * qO(t, z) - deltaH_PO * qO(t-1, z)) - THR * Dt_eqn * (Tc(t, z) - T_amb)   ;
*****************************************
*       cyclic steady state: CSS        *
*****************************************
CSEq00(t, z)$(t_start(t) and z_axis(z))  ..      yN(t, z) - yN(t+t_total, z)      =e= 0                                          ;
CSEq01(t, z)$(t_start(t) and z_axis(z))  ..      P(t, z)  - P(t+t_total, z)       =e= 0                                          ;
CSEq02(t, z)$(t_start(t) and z_axis(z))  ..      qN(t, z) - qN(t+t_total, z)      =e= 0                                          ;
CSEq03(t, z)$(t_start(t) and z_axis(z))  ..      qO(t, z) - qO(t+t_total, z)      =e= 0                                          ;
CSEq04(t, z)$(t_start(t) and z_axis(z))  ..      Tc(t, z) - Tc(t+t_total, z)      =e= 0                                          ;
*****************************************
*              connections              *
*****************************************
CCEq00(z)$z_right(z)                     ..      Adp_out                          =e= P_end_AP * sum(t_adp, u(t_adp, z) / Tc(t_adp, z)) * Dt_adp                                                                 ;
CCEq01(z)$z_right(z)                     ..      Pur_feed                         =e= - P_end_PU * sum(t_pur, u(t_pur, z) / Tc(t_pur, z)) * Dt_pur                                                               ;
CCEq02                                   ..      Adp_out                          =e= Pur_feed                                                                                                                   ;
CCEq03(z)$z_right(z)                     ..      Equ_out                          =e= P_end_EQ * sum(t_equ, u(t_equ, z) / Tc(t_equ, z)) * Dt_equ                                                                 ;
CCEq04(z)$z_right(z)                     ..      Eqn_feed                         =e= - P_end_EA * sum(t_eqn, u(t_eqn, z) / Tc(t_eqn, z)) * Dt_eqn                                                               ;
CCEq05                                   ..      Equ_out                          =e= Eqn_feed                                                                                                                   ;
*****************************************
*              mass balances            *
*****************************************
MBEq00(z)$z_right(z)                     ..      PO_ads_out                       =e= P_end * sum(t_ads, u(t_ads, z) * yO(t_ads, z) / Tc(t_ads, z)) * Dt_ads * CPRT                      ;
MBEq01(z)$z_right(z)                     ..      PO_adp_out                       =e= P_end_EQ * sum(t_adp, u(t_adp, z) * yO(t_adp, z) / Tc(t_adp, z)) * Dt_adp * CPRT                   ;
MBEq02(z)$z_right(z)                     ..      PO_pur_feed                      =e= - P_end_PU * sum(t_pur, u(t_pur, z) * yO(t_pur, z) / Tc(t_pur, z)) * Dt_pur * CPRT                 ;
MBEq03(z)$z_right(z)                     ..      P_ads_out                        =e= P_end * sum(t_ads, u(t_ads, z) / Tc(t_ads, z)) * Dt_ads * CPRT                                     ;
MBEq04(z)$z_right(z)                     ..      P_adp_out                        =e= P_end_EQ * sum(t_adp, u(t_adp, z) / Tc(t_adp, z)) * Dt_adp * CPRT                                  ;
*****************************************
*       separation specifications       *
*****************************************
SSEq00                                   ..      PurityO                          =e= PO_ads_out / P_ads_out                                                                             ;
SSEq01                                   ..      ProductivityO                    =e= - 0.5 * (PO_ads_out + PO_adp_out - PO_pur_feed) / (Dt_pre + Dt_ads + Dt_adp + Dt_equ)              ;
SSEq02                                   ..      PurityO                          =g= 0.8                                                                                                ;

***********************************************************************************

Model PSA  /all/;

option NLP = CONOPT4;
option reslim = 36000;
option domlim = 500000;

solve PSA using NLP minimizing ProductivityO;

$ontext
 - ProductivityO / 10: 20.3766 mol/s
Adsorbent mass: 3.1416*1.75*1.75*2*630*(1-0.36) = 7758.5 kg
1 mol/s/kg = 80640 Nm3/h/ton

result = - result * 80640 / 10 / 7758.5
$offtext
