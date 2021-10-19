function v = configure_fun(data)

%% Assign options

pathway_opt = {'C3'}; % Default is C3 pathway

% N.B., subsequent variables all have units (_u), and a description (_d)

%% Assign environmental measurements 

Q = data.Qin./1e6;
Q_u = 'mol PAR m-2 s-1';
Q_d = 'PPFD';

T = data.Tin;
T_u = 'degrees C';
T_d = 'Leaf temperature';

P = data.Pin;
P_u = 'bar';
P_d = 'Total pressure';

C_m = data.Cin./1e6;
C_m_u = 'bar';
C_m_d = 'Partial pressure of CO2';

O_m = data.Oin;
O_m_u = 'bar';
O_m_d = 'Partial pressure of O2';

%% Assign values for physiological variables

Abs = 0.85; 
Abs_u = 'mol absorbed mol-1 incident';
Abs_d = 'Total leaf absorptance to PAR';

abs_frac = 0;
abs_frac_u = 'dimensionless';
abs_frac_d = 'Bundle sheath fraction of total leaf absorptance to PAR';

a2_m_frac = 0.52;
a2_m_frac = 'dimensionless';
a2_m_frac = 'PS II fraction of mesophyll absorptance to PAR';

a2_s_frac = 0;
a2_s_frac = 'dimensionless';
a2_s_frac = 'PS II fraction of bundle sheath absorptance';

vq_frac = 0;
vq_frac_u = 'dimensionless';
vq_frac_d = 'Bundle sheath fraction of Cyt b6f density';

vc_frac = 0;
vc_frac_u = 'dimensionless';
vc_frac_d = 'Bundle sheath fraction of Rubisco density';

CB6F = (350./300)./1e6; 
CB6F_u = 'mol sites m-2';
CB6F_d = 'Cyt b6f density';

RUB = (100./3.6)./1e6;
RUB_u = 'mol sites m-2';
RUB_d = 'Rubisco density';

Rdsc = 0.01; 
Rdsc_u = 'dimensionless';
Rdsc_d = 'Scalar for mitochondrial (dark) respiration';

Vpmax = 0; 
Vpmax_u = 'mol CO2 m-2 s-1';
Vpmax_d = 'Maximum activity of PEP carboxylase';

gbs = 0.003;
gbs_u = 'mol CO2 m-2 s-1 bar-1';
gbs_d = 'Bundle sheath conductance to CO2';

gbso = gbs.*0.047;
gbso_u = 'mol O2 m-2 s-1 bar-1';
gbso_d = 'Bundle sheath conductance to O2';

%% Assign values for photochemical constants

Kf = 0.05e09;
Kf_u = 's-1';
Kf_d = 'Rate constant for fluoresence at PSII and PSI';

Kd = 0.55e09; 
Kd_u = 's-1';
Kd_d = 'Rate constant for constitutive heat loss at PSII and PSI';

Kp1 = 14.5e09; 
Kp1_u = 's-1';
Kp1_d = 'Rate constant for photochemistry at PSI';

Kn1 = 14.5e09; 
Kn1_u = 's-1';
Kn1_d = 'Rate constant for regulated heat loss at PSI';

Kp2 = 4.5e09; 
Kp2_u = 's-1';
Kp2_d = 'Rate constant for photochemistry at PSII';

Ku2 = 0e09; 
Ku2_u = 's-1';
Ku2_d = 'Rate constant for exciton sharing at PSII';

%% Assign values for biochemical constants 

kq = 300; 
kq_u = 'mol e-1 mol sites-1 s-1';
kq_d = 'Cyt b6f kcat for PQH2';

nl = 0.75; 
nl_u = 'ATP/e-';
nl_d = 'ATP per e- in linear flow';

nc = 1.00;
nc_u = 'ATP/e-';
nc_d = 'ATP per e- in cyclic flow';

kc = 3.6; 
kc_u = 'mol CO2 mol sites-1 s-1';
kc_d = 'Rubisco kcat for CO2';

ko = 3.6.*0.27; 
ko_u = 'mol O2 mol sites-1 s-1';
ko_d = 'Rubisco kcat for O2';

Kc = 260./1e6; 
Kc_u = 'bar';
Kc_d = 'Rubisco Km for CO2';

Ko = 179000./1e6;
Ko_u = 'bar';
Ko_d = 'Rubisco Km for O2';

Kp = 80./1e6;
Kp_u = 'bar';
Kp_d = 'PEPc Km for CO2';

%% Transfer functions for fluorescence
%   N.B., values of eps1 > 0 mix PSI F with PSII F in PAM calculations

eps1 = 0; 
eps1_u = 'mol PSI F to detector mol-1 PSI F emitted';
eps1_d = 'PS I transfer function';

eps2 = 1; 
eps2_u = 'mol PSII F to detector mol-1 PSII F emitted';
eps2_d = 'PS II transfer function';

%% Assign functions to blank slots in 'v'    

c3c4_solve_cc = 'NAN';
c3c4_solve_cj = 'NAN';
c3c4_solve_jc = 'NAN';
c3c4_solve_jj = 'NAN';
c4_solve_cc = 'NAN';
c4_solve_cj = 'NAN';
c4_solve_jc = 'NAN';
c4_solve_jj = 'NAN';
    
%% Create structure to hold variables

clearvars data;
v = workspace2struct_fun();

end