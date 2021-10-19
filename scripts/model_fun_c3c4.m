function m = model_fun_c3c4(v)

%% Update input variables

% Rename measured environmental variables to drive simulations
[pathway_opt, Q, T, P, O_m, C_m, ...
Abs, abs_frac, a2_m_frac, a2_s_frac, eps1, eps2, ...
CB6F, vq_frac, RUB, Rdsc, vc_frac, Vpmax, gbs, gbso, ...
Kf, Kd, Kp1, Kn1, Kp2, Ku2, kq, nl, nc, kc, ko, Kc, Ko, Kp, ...
c3c4_solve_cc, c3c4_solve_cj, c3c4_solve_jc, c3c4_solve_jj, ...
c4_solve_cc, c4_solve_cj, c4_solve_jc, c4_solve_jj] = loadvars_fun(v);

%% Calculate derived variables

% Absorptance of mesophyll (M) vs. bundle sheath (BS), and PS I vs. PS II 
Abs_m = Abs.*(1-abs_frac);  % M total, mol PPFD abs mol-1 PPFD incident
Abs_s = Abs.*(abs_frac);    % BS total, mol PPFD abs mol-1 PPFD incident
a2_m = Abs_m.*(a2_m_frac);  % PSII, mol PPFD abs PS2 mol-1 PPFD incident
a1_m = Abs_m - a2_m;        % PSI, mol PPFD abs PS1 mol-1 PPFD incident
a2_s = Abs_s.*(a2_s_frac);  % PSII, mol PPFD abs PS2 mol-1 PPFD incident
a1_s = Abs_s - a2_s;        % PSI, mol PPFD abs PS1 mol-1 PPFD incident

%% Temperature scaling

% Parameters at 25C are scaled to leaf temperature using three functions:
%
% (i) Activation only functions: 
%       P = P.*exp(Ha./R.*(1./Tref_K - 1./T_K));
% (ii) Activation and deactivation functions:
%       P = P.*exp(Ha./R.*(1./Tref_K - 1./T_K)).*...
%           (1+exp((Tref_K.*En-Hd)./(R.*Tref_K)))./...
%           (1+exp((T_K.*En-Hd)./(R.*T_K)))
% (iii) Deactivation only functions:
%       P = P.*(1+exp((Tref_K.*En-Hd)./(R.*Tref_K)))./...
%           (1+exp((T_K.*En-Hd)./(R.*T_K)))
%
% N.B., the parameterization of these functions is modified from:
%
% Farquhar GD, von Caemmerer S, Berry JA (1980a) A biochemical model of 
% photosynthetic CO2 assimilation in leaves of C3 species. Planta 149:78–90 
% https://doi.org/10.1007/bf00386231
%
% von Caemmerer S, Farquhar GD, Berry JA (2009) Biochemical Model 
% of C3 Photosynthesis. In: Laisk A, Nedbal L, Govindjee (eds) 
% Photosynthesis in silico: Understanding complexity from molecules 
% to ecosystems. Springer Netherlands, Dordrecht, pp 209–230 
% https://doi.org/10.1007/978-1-4020-9237-4_9

% Shared terms
T_K = T + repmat(273.15,length(T),1); % Leaf temperature, K
Tref_K = 25 + 273.15; % Reference temperature, K
R = 0.008314; % Ideal gas constant, kJ mol-1 K-1

% Cytochrome b6f complex
Ha = 37; % Activation energy, kJ mol-1 
kq = kq.*exp(Ha./R.*(1./Tref_K - 1./T_K)); % Cyt b6f kcat for PQH2, s-1
Vqmax = CB6F.*kq; % Max Cyt b6f activity, mol e- m-2 s-1
Vqmax_m = Vqmax.*(1-vq_frac); % Max Cyt b6f activity in M, mol e- m-2 s-1
Vqmax_s = Vqmax.*(vq_frac); % Max Cyt b6f activity in BS, mol e- m-2 s-1

% Coupling to proton circuit
En = 0.710; % Entropy, kJ K-1 mol-1
Hd = 220; % Enthalpy, kJ mol-1
nl = nl.*(1+exp((Tref_K.*En-Hd)./(R.*Tref_K)))./...
    (1+exp((T_K.*En-Hd)./(R.*T_K))); % ATP per e- in linear flow
nc = nc.*(1+exp((Tref_K.*En-Hd)./(R.*Tref_K)))./...
    (1+exp((T_K.*En-Hd)./(R.*T_K))); % ATP per e- in cyclic flow

% Mitochondrial respiration
Rd = RUB.*kc.*Rdsc; % Respiration, mol CO2 m-2 s-1
Ha = 66; % Activation energy, kJ mol-1 
Rd = Rd.*exp(Ha./R.*(1./Tref_K - 1./T_K)); % Resp., mol CO2 m-2 s-1
Rd_m = Rd.*(1-abs_frac); % Respiration in M, mol CO2 m-2 s-1
Rd_s = Rd.*(abs_frac); % Respiration in BS, mol CO2 m-2 s-1

% Rubisco
S = (kc./Kc).*(Ko./ko); % Rubisco specificity for CO2/O2, dimensionless
Ha = 23; % Activation energy, kJ mol-1 
S = 1./(1./S.*(exp(Ha./R.*(1./Tref_K - 1./T_K))));
Ha = 59; % Activation energy, kJ mol-1 
Kc = Kc.*exp(Ha./R.*(1./Tref_K - 1./T_K)); % Rubisco Km for CO2, bar
Ha = 36; % Activation energy, kJ mol-1 
Ko = Ko.*exp(Ha./R.*(1./Tref_K - 1./T_K)); % Rubisco Km for O2, bar
Ha = 58; % Activation energy, kJ mol-1 
kc = kc.*exp(Ha./R.*(1./Tref_K - 1./T_K)); % Rubisco kcat for CO2, s-1
Vcmax = RUB.*kc; % Max Rubisco activity, mol CO2 m-2 s-1
Vcmax_m = Vcmax.*(1-vc_frac); % Max Rubisco activity in M, mol CO2 m-2 s-1
Vcmax_s = Vcmax.*(vc_frac); % Max Rubisco activity in BS, mol CO2 m-2 s-1

% PEP carboxylase
Ha = 58; % Activation energy (nb, scaled to match kc, above), kJ mol-1 
Vpmax = Vpmax.*exp(Ha./R.*(1./Tref_K - 1./T_K));
Vpmax_m = Vpmax; % Max PEPC activity in mesophyll, mol CO2 m-2 s-1

%% Calculate limiting rates for gas-exchange and electron transport

% Symbol, description, units:
%   JP700,  PS I electron transport rate, mol e- m-2 s-1
%   JP680,  PS II electron transport rate, mol e- m-2 s-1
%   Vc,     Rubisco carboxylation rate, mol CO2 m-2 s-1
%   Vo,     Rubisco oxygenation rate, mol O2 m-2 s-1
%   Vg,     Glycine shuttling rate, mol C m-2 s-1
%   Vp,     Malate shuttling rate, mol C m-2 s-1
%   L,      Leak rate, mol CO2 m-2 s-1
%   Ag,     Gross CO2 assimilation rate, mol CO2 m-2 s-1
%   An,     Net CO2 assimilation rate, mol CO2 m-2 s-1
%
% Subscripts:
%   _m,     in the mesophyll
%   _s,     in the bundle sheath
%   _j,     light-lighted state
%   _jj,    light-limited states in mesophyll and bundle sheath
%   _jc,    light-limited mesophyll and light-saturated bundle sheath
%   _c,     light-saturated state
%   _cc,    light-saturated states in mesophyll and bundle sheath
%   _cj,    light-saturated state mesophyll and light-limited bundle sheath
%   _a,     actual state
%
% ------------------------------------------------------------------- %
%   1. MESOPHYLL
% ------------------------------------------------------------------- %

if Vpmax_m == 0    
if strcmp(pathway_opt,'C3') == 1 || strcmp(pathway_opt,'Type-I-C3-C4') == 1
        
% Cytochrome b6f-limited rates (_mj), OECO-21: Eqn. 7a, 7b, 11a 
JP700_mj = Q.*Vqmax_m./(Q+Vqmax_m./(a1_m.*(Kp1./(Kp1 + Kd + Kf))));
JP680_mj = JP700_mj./(1-(nl./nc)+(3+7.*O_m./(2.*S.*C_m))./...
    ((4+4.*O_m./(S.*C_m)).*nc));
Vc_mj = JP680_mj./(4.*(1+O_m./(S.*C_m)));
Vo_mj = Vc_mj.*O_m./(S.*C_m);
Ag_mj = Vc_mj - Vo_mj./2;
An_mj = Ag_mj - Rd_m;
Vp_mj = 0;

% Rubisco-limited rates (_mc), OECO-21: Eqn. 7a, 7b, 17a
Vc_mc = C_m.*Vcmax_m./(C_m + Kc.*(1+O_m./Ko));
Vo_mc = Vc_mc.*O_m./(S.*C_m);
Ag_mc = Vc_mc - Vo_mc./2;
An_mc = Ag_mc - Rd_m;
JP680_mc = Ag_mc.*4.*(1+O_m./(S.*C_m))./(1-O_m./(2.*S.*C_m)); 
JP700_mc = JP680_mc.*(1-(nl./nc)+(3+7.*O_m./(2.*S.*C_m))./...
    ((4+4.*O_m./(S.*C_m)).*nc));
Vp_mc = 0;

% Rates of glycine shuttling (_mj,_mc), OECO-21: Eqn. 12a, 18a    
if Vcmax_s == 0
    Vg_mj = 0;
    Vg_mc = 0;
else
    Vg_mj = Vo_mj./2;
    Vg_mc = Vo_mc./2;
end
        
end
end

if Vpmax_m > 0  
if strcmp(pathway_opt,'NADP-ME-C4') == 1

% Cytochrome b6f-limited rates (_mj), OECO-21: Eqn. 9a, 11a, 12b 
JP700_mj = Q.*Vqmax_m./(Q+Vqmax_m./(a1_m.*(Kp1./(Kp1 + Kd + Kf))));
JP680_mj = JP700_mj./(1 - (nl./nc) + 1./nc); 
Vp_mj = JP680_mj./2; 

% PEPC-limited rates (_mc), OECO-21: Eqn. 9a, 12b, 18b
Vp_mc = Vpmax_m .*C_m./(Kp + C_m);
JP680_mc = Vp_mc.*2;  
JP700_mc = JP680_mc.*(1 - (nl./nc) + 1./nc);  

% Other terms, OECO-21: Eqn. 9b
Vc_mj = 0;
Vo_mj = 0;
Vg_mj = 0;
Ag_mj = 0;
An_mj = - Rd_m;
Vc_mc = 0;
Vo_mc = 0;
Vg_mc = 0;
Ag_mc = 0;
An_mc = - Rd_m;

end
end

% ------------------------------------------------------------------- %
%   2. BUNDLE SHEATH
% ------------------------------------------------------------------- %

% ------------------------------------------------------------------- %
%   2a. C3 solution
% ------------------------------------------------------------------- %

if Vcmax_s == 0
    
JP700_sjj = 0;
JP680_sjj = 0;
JP700_sjc = 0;
JP680_sjc = 0;
JP700_scc = 0;
JP680_scc = 0;
JP700_scj = 0;
JP680_scj = 0;
JP700_sj = 0;
JP680_sj = 0;
JP700_sc = 0;
JP680_sc = 0;
JP700_sa = 0;
JP680_sa = 0;
which_JP700_sj = 0;
which_JP700_sc = 0;
Ag_sa = 0;
An_sa = 0;
C_sa = 0;
O_sa = 0;
L_C_sa = 0;

end

% ------------------------------------------------------------------- %
%   2b. C3-C4 solutions
% ------------------------------------------------------------------- %

if mean(Vcmax_s) > 0 && mean(Vpmax_m) == 0
    
% Cyt b6f-limited M & BS (_sjj),
%   OECO-21: Eqn. 8a, 8b, 11b, 12a, 13-16
JP700_sjj = Q.*Vqmax_s./(Q+Vqmax_s./(a1_s.*(Kp1./(Kp1 + Kd + Kf))));
An_sjj = c3c4_solve_jj(C_m,JP680_mj,JP700_sjj,O_m,P,Rd_s,S,gbs,gbso,nc,nl);
C_sjj = C_m + ((JP680_mj.*O_m./(2.*S.*C_m)./...
    (4.*(1+O_m./(S.*C_m)))) - An_sjj).*P./gbs;
O_sjj = O_m + (An_sjj).*P./gbso;
Ag_sjj = An_sjj + Rd_s;
JP680_sjj = Ag_sjj.*4.*(1+O_sjj./(S.*C_sjj))./(1-O_sjj./(2.*S.*C_sjj)); 
Vc_sjj = JP680_sjj./(4.*(1+O_sjj./(S.*C_sjj)));
Vo_sjj = Vc_sjj.*O_sjj./(S.*C_sjj);

% Cyt b6f-limited M & Rubiso-limited BS (_sjc), 
%   OECO-21: Eqn. 8a, 8b, 12a, 13, 15a, 15b, 17b, 21a, 21b
An_sjc = c3c4_solve_jc(C_m,JP680_mj,Kc,Ko,O_m,P,Rd_s,S,Vcmax_s,gbs,gbso);
C_sjc = C_m + ((JP680_mj.*O_m./(2.*S.*C_m)./...
    (4.*(1+O_m./(S.*C_m)))) - An_sjc).*P./gbs;
O_sjc = O_m + (An_sjc).*P./gbso;
Ag_sjc = An_sjc + Rd_s;
JP680_sjc = Ag_sjc.*4.*(1+O_sjc./(S.*C_sjc))./(1-O_sjc./(2.*S.*C_sjc)); 
JP700_sjc = JP680_sjc.*(1-(nl./nc)+(3+7.*O_sjc./(2.*S.*C_sjc))./...
    ((4+4.*O_sjc./(S.*C_sjc)).*nc));
Vc_sjc = JP680_sjc./(4.*(1+O_sjc./(S.*C_sjc)));
Vo_sjc = Vc_sjc.*O_sjc./(S.*C_sjc);  

% Rubisco-limited M & Cyt b6f-limited BS (_scj),
%   OECO-21: Eqn. 8a, 8b, 11b, 13, 15a, 15b, 18a, 22a, 22b
JP700_scj = JP700_sjj;
An_scj = c3c4_solve_cj(C_m,JP700_scj,Kc,Ko,O_m,P,Rd_s,S,Vcmax_m,gbs,gbso,nc,nl);
C_scj = C_m + ((Vcmax_m.*O_m./(2.*S)./...
    (C_m + Kc.*(1+O_m./Ko))) - An_scj).*P./gbs;
O_scj = O_m + (An_scj).*P./gbso;
Ag_scj = An_scj + Rd_s;
JP680_scj = Ag_scj.*4.*(1+O_scj./(S.*C_scj))./(1-O_scj./(2.*S.*C_scj)); 
Vc_scj = JP680_scj./(4.*(1+O_scj./(S.*C_scj)));
Vo_scj = Vc_scj.*O_scj./(S.*C_scj);

% Rubisco-limited M & BS (_scc),
%   OECO-21: Eqn. 8a, 8b, 13, 15a, 15b, 17b, 18a, 19, 20
An_scc = c3c4_solve_cc(C_m,Kc,Ko,O_m,P,Rd_s,S,Vcmax_m,Vcmax_s,gbs,gbso);
C_scc = C_m + ((Vcmax_m.*O_m./(2.*S)./(C_m + Kc.*(1+O_m./Ko))) - An_scc).*P./gbs;
O_scc = O_m + (An_scc).*P./gbso;
Ag_scc = An_scc + Rd_s;
JP680_scc = Ag_scc.*4.*(1+O_scc./(S.*C_scc))./(1-O_scc./(2.*S.*C_scc)); 
JP700_scc = JP680_scc.*(1-(nl./nc)+(3+7.*O_scc./(2.*S.*C_scc))./...
    ((4+4.*O_scc./(S.*C_scc)).*nc));
Vc_scc = JP680_scc./(4.*(1+O_scc./(S.*C_scc)));
Vo_scc = Vc_scc.*O_scc./(S.*C_scc);
        
end

% ------------------------------------------------------------------- %
%   2c. C4 solutions
% ------------------------------------------------------------------- %

if Vpmax_m > 0

% ------------------------------------------------------------------- %
%   NADP-ME type C4 solutions
% ------------------------------------------------------------------- %

if strcmp(pathway_opt,'NADP-ME-C4') == 1
    
% Cyt b6f-limited M & BS (_sjj),
%   OECO-21: Eqn. 10a, 10b, 11b, 12b, 13-16
JP700_sjj = Q.*Vqmax_s./(Q+Vqmax_s./(a1_s.*(Kp1./(Kp1 + Kd + Kf))));
An_sjj = c4_solve_jj(C_m,JP700_mj,JP700_sjj,O_m,P,Rd_s,S,gbs,gbso,nc,nl);
C_sjj = C_m + (Vp_mj - An_sjj).*P./gbs;
O_sjj = O_m + (An_sjj - Vp_mj./2).*P./gbso;
Ag_sjj = An_sjj + Rd_s;
JP680_sjj = Ag_sjj.*4.*...
    (1+O_sjj./(S.*C_sjj))./(1-O_sjj./(2.*S.*C_sjj)) - 2.*Vp_mj; 

% Rubiso-limited BS & Cyt b6f-limited M (_sjc), 
%   OECO-21: Eqn. 10a, 10b, 12b, 13, 15a, 15b, 17b, 21a, 21b
An_sjc = c4_solve_jc(C_m,JP700_mj,Kc,Ko,O_m,P,Rd_s,S,Vcmax_s,gbs,gbso,nc,nl); 
C_sjc = C_m + (Vp_mj - An_sjc).*P./gbs;
O_sjc = O_m + (An_sjc - Vp_mj./2).*P./gbso;
Ag_sjc = An_sjc + Rd_s;
JP680_sjc = Ag_sjc.*4.*(1+O_sjc./(S.*C_sjc))./...
    (1-O_sjc./(2.*S.*C_sjc)) - 2.*Vp_mj; 
JP700_sjc = JP680_sjc.*(1-(nl./nc)+(3+7.*O_sjc./(2.*S.*C_sjc))./...
    ((4+4.*O_sjc./(S.*C_sjc)).*nc)) + (3+7.*O_sjc./(2.*S.*C_sjc))./...
    ((4+4.*O_sjc./(S.*C_sjc)).*nc).*(Vp_mj./2);

% Cyt b6f-limited BS & PEPC-limited M (_scj),
%   OECO-21: Eqn. 10a, 10b, 11b, 13, 15a, 15b, 18b, 22a, 22b
JP700_scj = Q.*Vqmax_s./(Q+Vqmax_s./(a1_s.*(Kp1./(Kp1 + Kd + Kf))));
An_scj = c4_solve_cj(C_m,JP700_scj,Kp,O_m,P,Rd_s,S,Vpmax_m,gbs,gbso,nc,nl);
C_scj = C_m + (Vp_mc - An_scj).*P./gbs;
O_scj = O_m + (An_scj - Vp_mc./2).*P./gbso; 
Ag_scj = An_scj + Rd_s;
JP680_scj = Ag_scj.*4.*(1+O_scj./(S.*C_scj))./...
    (1-O_scj./(2.*S.*C_scj)) - 2.*Vp_mc; 

% Rubisco-limited BS & PEPC-limited M (_scc),
%   OECO-21: Eqn. 10a, 10b, 13, 15a, 15b, 17b, 18b, 19, 20
An_scc = c4_solve_cc(C_m,Kc,Ko,Kp,O_m,P,Rd_s,S,Vcmax_s,Vpmax_m,gbs,gbso); 
C_scc = C_m + (Vp_mc - An_scc).*P./gbs;
O_scc = O_m + (An_scc - Vp_mc./2).*P./gbso;
Ag_scc = An_scc + Rd_s;
JP680_scc = (Ag_scc.*4.*(1+O_scc./(S.*C_scc))./...
    (1-O_scc./(2.*S.*C_scc))) - 2.*Vp_mc; 
JP700_scc = JP680_scc.*(1-(nl./nc)+(3+7.*O_scc./(2.*S.*C_scc))./...
    ((4+4.*O_scc./(S.*C_scc)).*nc)) + (3+7.*O_scc./(2.*S.*C_scc))./...
    ((4+4.*O_scc./(S.*C_scc)).*nc).*(Vp_mc./2);
    
end

end

%% Select min of light-limited and light-saturated rates
% Use min(A,[],2) for column vector with minimum value of each row

% 1. Actual state of mesophyll, OECO-21: Eqn. 23

    [JP700_ma,which_JP700_ma] = min([JP700_mj,JP700_mc],[],2);
    JP680_ma = JP680_mj.*logical(which_JP700_ma == 1) + ...
               JP680_mc.*logical(which_JP700_ma == 2);
    Vg_ma = Vg_mj.*logical(which_JP700_ma == 1) + ...
            Vg_mc.*logical(which_JP700_ma == 2);
    Vp_ma = Vp_mj.*logical(which_JP700_ma == 1) + ...
            Vp_mc.*logical(which_JP700_ma == 2);
    An_ma = An_mj.*logical(which_JP700_ma == 1) + ...
            An_mc.*logical(which_JP700_ma == 2); 
    Ag_ma = An_ma + Rd_m;          

% 2. Potential states of bundle sheath

if Vcmax_s > 0 
    
    % Mesophyll Cyt b6f-limited, OECO-21: Eqn. 24a, 24b, 24c
    [JP700_sj,which_JP700_sj] = min([JP700_sjj,JP700_sjc],[],2);
    JP680_sj = JP680_sjj.*logical(which_JP700_sj == 1) + ...
               JP680_sjc.*logical(which_JP700_sj == 2);
    C_sj = C_sjj.*logical(which_JP700_sj == 1) + ...
           C_sjc.*logical(which_JP700_sj == 2);
    O_sj = O_sjj.*logical(which_JP700_sj == 1) + ...
           O_sjc.*logical(which_JP700_sj == 2);
    An_sj = An_sjj.*logical(which_JP700_sj == 1) + ...
            An_sjc.*logical(which_JP700_sj == 2);
        
    % Mesophyll Rubisco-limited, OECO-21: Eqn. 24a, 24b, 24c
    [JP700_sc,which_JP700_sc] = min([JP700_scj,JP700_scc],[],2);
    JP680_sc = JP680_scj.*logical(which_JP700_sc == 1) + ...
               JP680_scc.*logical(which_JP700_sc == 2);
    C_sc = C_scj.*logical(which_JP700_sc == 1) + ...
           C_scc.*logical(which_JP700_sc == 2);
    O_sc = O_scj.*logical(which_JP700_sc == 1) + ...
           O_scc.*logical(which_JP700_sc == 2);
    An_sc = An_scj.*logical(which_JP700_sc == 1) + ...
            An_scc.*logical(which_JP700_sc == 2);
       
% 3. Actual state of bundle sheath

    JP700_sa = JP700_sj.*logical(which_JP700_ma == 1) + ...
               JP700_sc.*logical(which_JP700_ma == 2);
    JP680_sa = JP680_sj.*logical(which_JP700_ma == 1) + ...
               JP680_sc.*logical(which_JP700_ma == 2);
    C_sa = C_sj.*logical(which_JP700_ma == 1) + ...
           C_sc.*logical(which_JP700_ma == 2);
    O_sa = O_sj.*logical(which_JP700_ma == 1) + ...
           O_sc.*logical(which_JP700_ma == 2);
    An_sa = An_sj.*logical(which_JP700_ma == 1) + ...
            An_sc.*logical(which_JP700_ma == 2);
    L_C_sa = gbs./P.*(C_sa - C_m);    
    Ag_sa = An_sa + Rd_s;
    
end

% 4. Totals for leaf, OECO-21: Eqn. 25
    JP700_a = JP700_ma + JP700_sa;
    JP680_a = JP680_ma + JP680_sa;
    An_a = An_ma + An_sa;
    Ag_a = An_a + Rd;
    
%% Derive fluorescence parameters from gas-exchange and electron transport
%
% The calculations below follow the approach in PRES-21:
%
%   Johnson, JE, Berry, JA (2021) The role of Cytochrome b6f in
%   the control of steady-state photosynthesis: a conceptual
%   and quantitative model. Photosynthesis Research 148: 101–136.

% Primary fluorescence parameters for mesophyll
CB6F_m = CB6F.*(1-vq_frac);
RUB_m = RUB.*(1-vc_frac);
CB6F_ma = JP700_mj./kq;                 % PRES-21 Eqns. 21, 30a, 34
phi1P_ma = JP700_ma./(Q.*a1_m);         % PRES-21 Eqn. 20
q1_ma = phi1P_ma.*((Kp1+Kd+Kf)./Kp1);   % PRES-21 Eqn. 19a
phi2P_ma = JP680_ma./(Q.*a2_m);         % PRES-21 Eqn. 26  
q2_ma = 1-CB6F_ma./CB6F_m;              % PRES-21 Eqn. 28, 34

% N.B., rearrange PRES-21 Eqn. 25a to solve for Kn2_a
Kn2_ma = ((Kp2.^2.*phi2P_ma.^2-2.*Kp2.^2.*phi2P_ma.*q2_ma+...
   Kp2.^2.*q2_ma.^2-4.*Kp2.*Ku2.*phi2P_ma.^2.*q2_ma+...
   2.*Kp2.*Ku2.*phi2P_ma.^2+2.*Kp2.*Ku2.*phi2P_ma.*q2_ma+...
   Ku2.^2.*phi2P_ma.^2).^(1./2)-Kp2.*phi2P_ma+Ku2.*phi2P_ma+...
   Kp2.*q2_ma)./(2.*phi2P_ma)-Kf-Ku2-Kd;

% Primary fluorescence parameters for bundle sheath
if Vpmax_m == 0

CB6F_s = 0;
RUB_s = 0;
CB6F_sa = 0;
phi1P_sa = 0;
q1_sa = 0;
phi2P_sa = 0;
q2_sa = 0;
Kn2_sa = 0;
    
end

if Vcmax_s > 0

CB6F_s = CB6F.*(vq_frac); 
RUB_s = RUB.*(vc_frac); 
CB6F_sa = JP700_sjj./kq;                % PRES-21 Eqns. 21, 30a, 34
phi1P_sa = JP700_sa./(Q.*a1_s);         % PRES-21 Eqn. 20
q1_sa = phi1P_sa.*((Kp1+Kd+Kf)./Kp1);   % PRES-21 Eqn. 19a
phi2P_sa = JP680_sa./(Q.*a2_s);         % PRES-21 Eqn. 26  
q2_sa = 1-CB6F_sa./CB6F_s;              % PRES-21 Eqn. 28, 34

% N.B., rearrange PRES-21 Eqn. 25a to solve for Kn2_a
Kn2_sa = ((Kp2.^2.*phi2P_sa.^2-2.*Kp2.^2.*phi2P_sa.*q2_sa+...
   Kp2.^2.*q2_sa.^2-4.*Kp2.*Ku2.*phi2P_sa.^2.*q2_sa+...
   2.*Kp2.*Ku2.*phi2P_sa.^2+2.*Kp2.*Ku2.*phi2P_sa.*q2_sa+...
   Ku2.^2.*phi2P_sa.^2).^(1./2)-Kp2.*phi2P_sa+Ku2.*phi2P_sa+...
   Kp2.*q2_sa)./(2.*phi2P_sa)-Kf-Ku2-Kd;

end

% Derived fluorescence parameters -- 'True values'

% For Photosystem II, PRES-21 Eqns. 23a-23e and 25a-25d
phi2p_ma = (q2_ma).*Kp2./(Kp2+Kn2_ma+Kd+Kf+Ku2);
phi2n_ma = (q2_ma).*Kn2_ma./(Kp2+Kn2_ma+Kd+Kf+Ku2)+...
    (1-q2_ma).*Kn2_ma./(Kn2_ma+Kd+Kf+Ku2);
phi2d_ma = (q2_ma).*Kd./(Kp2+Kn2_ma+Kd+Kf+Ku2)+...
    (1-q2_ma).*Kd./(Kn2_ma+Kd+Kf+Ku2);
phi2f_ma = (q2_ma).*Kf./(Kp2+Kn2_ma+Kd+Kf+Ku2)+...
    (1-q2_ma).*Kf./(Kn2_ma+Kd+Kf+Ku2);
phi2u_ma = (q2_ma).*Ku2./(Kp2+Kn2_ma+Kd+Kf+Ku2)+...
    (1-q2_ma).*Ku2./(Kn2_ma+Kd+Kf+Ku2);
phi2P_ma = phi2p_ma./(1-phi2u_ma);
phi2N_ma = phi2n_ma./(1-phi2u_ma);
phi2D_ma = phi2d_ma./(1-phi2u_ma);
phi2F_ma = phi2f_ma./(1-phi2u_ma);

% For Photosystem I, PRES-21 Eqns. 19a-19d
phi1P_ma = (q1_ma).*Kp1./(Kp1+Kd+Kf);
phi1N_ma = (1-q1_ma).*Kn1./(Kn1+Kd+Kf);
phi1D_ma = (q1_ma).*Kd./(Kp1+Kd+Kf)+(1-q1_ma).*Kd./(Kn1+Kd+Kf);
phi1F_ma = (q1_ma).*Kf./(Kp1+Kd+Kf)+(1-q1_ma).*Kf./(Kn1+Kd+Kf);

% For Photosystem II, PRES-21 Eqns. 23a-23e and 25a-25d
phi2p_sa = (q2_sa).*Kp2./(Kp2+Kn2_sa+Kd+Kf+Ku2);
phi2n_sa = (q2_sa).*Kn2_sa./(Kp2+Kn2_sa+Kd+Kf+Ku2)+ ...
    (1-q2_sa).*Kn2_sa./(Kn2_sa+Kd+Kf+Ku2);
phi2d_sa = (q2_sa).*Kd./(Kp2+Kn2_sa+Kd+Kf+Ku2)+ ...
    (1-q2_sa).*Kd./(Kn2_sa+Kd+Kf+Ku2);
phi2f_sa = (q2_sa).*Kf./(Kp2+Kn2_sa+Kd+Kf+Ku2)+ ...
    (1-q2_sa).*Kf./(Kn2_sa+Kd+Kf+Ku2);
phi2u_sa = (q2_sa).*Ku2./(Kp2+Kn2_sa+Kd+Kf+Ku2)+...
    (1-q2_sa).*Ku2./(Kn2_sa+Kd+Kf+Ku2);
phi2P_sa = phi2p_sa./(1-phi2u_sa);
phi2N_sa = phi2n_sa./(1-phi2u_sa);
phi2D_sa = phi2d_sa./(1-phi2u_sa);
phi2F_sa = phi2f_sa./(1-phi2u_sa);

% For Photosystem I, PRES-21 Eqns. 19a-19d
phi1P_sa = (q1_sa).*Kp1./(Kp1+Kd+Kf);
phi1N_sa = (1-q1_sa).*Kn1./(Kn1+Kd+Kf);
phi1D_sa = (q1_sa).*Kd./(Kp1+Kd+Kf)+(1-q1_sa).*Kd./(Kn1+Kd+Kf);
phi1F_sa = (q1_sa).*Kf./(Kp1+Kd+Kf)+(1-q1_sa).*Kf./(Kn1+Kd+Kf);

% Derived fluorescence parameters -- 'Observed values'  

% PAM measured fluorescence levels for mesophyll, PRES-21 Eqns. 38-42
Fs_ma = a2_m.*phi2F_ma.*eps2+a1_m.*phi1F_ma.*eps1;
Fm_ma = a2_m.*Kf./(Kd+Kf).*eps2+a1_m.*Kf./(Kn1+Kd+Kf).*eps1;
Fo_ma = a2_m.*Kf./(Kp2+Kd+Kf).*eps2+a1_m.*Kf./(Kp1+Kd+Kf).*eps1;
Fmp_ma = a2_m.*Kf./(Kn2_ma+Kd+Kf).*eps2+a1_m.*Kf./(Kn1+Kd+Kf).*eps1;
Fop_ma = a2_m.*Kf./(Kp2+Kn2_ma+Kd+Kf).*eps2+a1_m.*Kf./(Kp1+Kd+Kf).*eps1;

% PAM measured fluorescence levels for bundle sheath, PRES-21 Eqns. 38-42
Fs_sa = a2_s.*phi2F_sa.*eps2+a1_s.*phi1F_sa.*eps1;
Fm_sa = a2_s.*Kf./(Kd+Kf).*eps2+a1_s.*Kf./(Kn1+Kd+Kf).*eps1;
Fo_sa = a2_s.*Kf./(Kp2+Kd+Kf).*eps2+a1_s.*Kf./(Kp1+Kd+Kf).*eps1;
Fmp_sa = a2_s.*Kf./(Kn2_sa+Kd+Kf).*eps2+a1_s.*Kf./(Kn1+Kd+Kf).*eps1;
Fop_sa = a2_s.*Kf./(Kp2+Kn2_sa+Kd+Kf).*eps2+a1_s.*Kf./(Kp1+Kd+Kf).*eps1;

% Sum mesophyll and bundle sheath
Fs_a = Fs_ma+Fs_sa;
Fm_a = Fm_ma+Fm_sa;
Fo_a = Fo_ma+Fo_sa;
Fmp_a = Fmp_ma+Fmp_sa;
Fop_a = Fop_ma+Fop_sa;

% PAM indices derived from fluorescence levels
PAM1_a = (Fmp_a-Fs_a)./(Fmp_a-Fop_a);                   % qP
PAM2_a = (Fmp_a-Fs_a).*Fop_a./((Fmp_a-Fop_a).*Fs_a);    % qL
PAM3_a = Fm_a./Fmp_a-1;                                 % NPQ
PAM4_a = 1-Fs_a./Fmp_a;                                 % PhiP
PAM5_a = Fs_a.*(1./Fmp_a-1./Fm_a);                      % PhiN
PAM6_a = Fs_a./Fm_a;                                    % PhiD + PhiF

%% Clean up and save results

% Remove model input variables from workspace
clearvars   v ...
            Abs CB6F RUB Rdsc ...
            Kf Kd Kp1 Kn1 Kp2 Ku2 ...
            kq nl nc kc ko Kc Ko ...
            c3c4_solve_cc c3c4_solve_cj c3c4_solve_jc c3c4_solve_jj ...
            c4_solve_cc c4_solve_cj c4_solve_jc c4_solve_jj...
            eps1 eps2;

% Create structure to hold all remaining model output
m = workspace2struct_fun();

end

