function ss = symsolver_c3c4_fun(pathway_opt)

%% Clean up working environment
clearvars -except currdir workdir datadir resultsdir pathway_opt v; 
warning('off'); 

%% Define symbols
syms Ag_sc An_sc Ag_sj An_sj C_s C_m gbs gbso JP700_sj JP700_mj ...
    JP680_sj JP680_mj Kc Ko Kp Kp1 Kp2 Kd Kf Ku2 kc ko kq nl nc ...
    O_s O_m P phi2P_ma phi2p_ma phi2u_ma phi2P_sa phi2p_sa ...
    phi2u_sa Rd_s S Vg_c Vp_c Vg_j Vp_j Vc_mc Vo_mc Vc_sc ...
    Vo_sc Vc_mj Vo_mj Vc_sj Vo_sj Vpmax_m Vcmax_s Vcmax_m;

%% Solutions for Type I C3-C4 net CO2 assimilation in bundle sheath (As)      

if strcmp(pathway_opt,'Type-I-C3-C4') == 1

%% (1) Cyt b6f-limitation M & BS ---------------------------------------- %
eq1 = Ag_sj - Rd_s - An_sj == 0;
eq2 = subs(eq1, Ag_sj, Vc_sj - Vo_sj./2);
eq3 = solve(eq2, An_sj);
eq4 = eq3 - An_sj == 0;
eq5 = subs(eq4, Vo_sj, Vc_sj.*O_s./(S.*C_s));
eq6 = solve(eq5, An_sj);
eq7 = eq6 - An_sj == 0;
eq8 = subs(eq7, Vc_sj, JP680_sj./(4.*(1+O_s./(S.*C_s))));
eq9 = solve(eq8, An_sj);
eq10 = eq9 - An_sj == 0;
eq11 = subs(eq10, JP680_sj, JP700_sj./...
    (1-(nl./nc)+(3+3.5.*O_s./(S.*C_s))./((4+4.*O_s./(S.*C_s)).*nc)));
eq12 = solve(eq11, An_sj);
eq13 = eq12 - An_sj == 0;
eq14 = subs(eq13, C_s, C_m + (Vg_j - An_sj).*P./gbs);
eq15 = solve(eq14, An_sj,'PrincipalValue',true);
eq16 = eq15 - An_sj == 0;
eq17 = subs(eq16, Vg_j, JP680_mj.*O_m./(2.*S.*C_m)./(4.*(1+O_m./(S.*C_m))));
eq18 = solve(eq17, An_sj);
eq19 = eq18 - An_sj == 0;
eq20 = subs(eq19, O_s, O_m + (An_sj).*P./gbso);
c3c4_soln_jj = solve(eq20, An_sj); 
clear eq1 eq2 eq3 eq4 eq5 eq6 eq7 eq8 eq9 eq10...
      eq11 eq12 eq13 eq14 eq15 eq16 eq17 eq18 eq19 eq20;
% @(C_m,JP680_mj,JP700_sj,O_m,P,Rd_s,S,gbs,gbso,nc,nl)

%% (2) Cyt b6f-limited M, Rubiso-limited BS ----------------------------- %        
eq1 = Ag_sc - Rd_s - An_sc == 0;
eq2 = subs(eq1, Ag_sc, Vc_sc - Vo_sc./2);
eq3 = solve(eq2, An_sc);
eq4 = eq3 - An_sc == 0;
eq5 = subs(eq4, Vo_sc, Vc_sc.*O_s./(S.*C_s));
eq6 = solve(eq5, An_sc);
eq7 = eq6 - An_sc == 0;
eq8 = subs(eq7, Vc_sc, C_s.*Vcmax_s./(C_s + Kc.*(1+O_s./Ko)));
eq9 = solve(eq8, An_sc);
eq10 = eq9 - An_sc == 0;
eq11 = subs(eq10, C_s, C_m + (Vg_j - An_sc).*P./gbs);
eq12 = solve(eq11, An_sc,'PrincipalValue',true);
eq13 = eq12 - An_sc == 0;
eq14 = subs(eq13, Vg_j, JP680_mj.*O_m./(2.*S.*C_m)./(4.*(1+O_m./(S.*C_m))));
eq15 = solve(eq14, An_sc);
eq16 = eq15 - An_sc == 0;
eq17 = subs(eq16, O_s, O_m + (An_sc).*P./gbso);
c3c4_soln_jc = solve(eq17, An_sc); 
clear eq1 eq2 eq3 eq4 eq5 eq6 eq7 eq8 eq9 eq10...
      eq11 eq12 eq13 eq14 eq15 eq16 eq17;
% @(C_m,JP680_mj,Kc,Ko,O_m,P,Rd_s,S,Vcmax_s,gbs,gbso)

%% (3) Rubisco-limitation M & Cyt b6f-limitation BS --------------------- %  
eq1 = Ag_sj - Rd_s - An_sj == 0;
eq2 = subs(eq1, Ag_sj, Vc_sj - Vo_sj./2);
eq3 = solve(eq2, An_sj);
eq4 = eq3 - An_sj == 0;
eq5 = subs(eq4, Vo_sj, Vc_sj.*O_s./(S.*C_s));
eq6 = solve(eq5, An_sj);
eq7 = eq6 - An_sj == 0;
eq8 = subs(eq7, Vc_sj, JP680_sj./(4.*(1+O_s./(S.*C_s))));
eq9 = solve(eq8, An_sj);
eq10 = eq9 - An_sj == 0;
eq11 = subs(eq10, JP680_sj, JP700_sj./...
    (1-(nl./nc)+(3+3.5.*O_s./(S.*C_s))./((4+4.*O_s./(S.*C_s)).*nc)));
eq12 = solve(eq11, An_sj);
eq13 = eq12 - An_sj == 0;
eq14 = subs(eq13, C_s, C_m + (Vg_c - An_sj).*P./gbs);
eq15 = solve(eq14, An_sj,'PrincipalValue',true);
eq16 = eq15 - An_sj == 0;
eq17 = subs(eq16, Vg_c, Vcmax_m.*O_m./(2.*S)./(C_m + Kc.*(1+O_m./Ko)));
eq18 = solve(eq17, An_sj);
eq19 = eq18 - An_sj == 0;
eq20 = subs(eq19, O_s, O_m + (An_sj).*P./gbso);
c3c4_soln_cj = solve(eq20, An_sj); 
clear eq1 eq2 eq3 eq4 eq5 eq6 eq7 eq8 eq9 eq10 ...
      eq11 eq12 eq13 eq14 eq15 eq16 eq17 eq18 eq19 eq20;
% @(C_m,JP700_sj,Kc,Ko,O_m,P,Rd_s,S,Vcmax_m,gbs,gbso,nc,nl)

%% (4) Rubisco-limited M & BS ------------------------------------------- %
eq1 = Ag_sc - Rd_s - An_sc == 0;
eq2 = subs(eq1, Ag_sc, Vc_sc - Vo_sc./2);
eq3 = solve(eq2, An_sc);
eq4 = eq3 - An_sc == 0;
eq5 = subs(eq4, Vo_sc, Vc_sc.*2.*O_s./(2.*S.*C_s));
eq6 = solve(eq5, An_sc);
eq7 = eq6 - An_sc == 0;
eq8 = subs(eq7, Vc_sc, C_s.*Vcmax_s./(C_s + Kc.*(1+O_s./Ko)));
eq9 = solve(eq8, An_sc);
eq10 = eq9 - An_sc == 0;
eq11 = subs(eq10, C_s, C_m + (Vg_c - An_sc).*P./gbs);
eq12 = solve(eq11, An_sc,'PrincipalValue',true);
eq13 = eq12 - An_sc == 0;
eq14 = subs(eq13, Vg_c, Vcmax_m.*O_m./(2.*S)./(C_m + Kc.*(1+O_m./Ko)));
eq15 = solve(eq14, An_sc);
eq16 = eq15 - An_sc == 0;
eq17 = subs(eq16, O_s, O_m + (An_sc).*P./gbso);
c3c4_soln_cc = solve(eq17, An_sc); 
clear eq1 eq2 eq3 eq4 eq5 eq6 eq7 eq8 eq9 eq10...
      eq11 eq12 eq13 eq14 eq15 eq16 eq17;
% @(C_m,Kc,Ko,O_m,P,Rd_s,S,Vcmax_m,Vcmax_s,gbs,gbso)

% Export Type I C3-C4 functions:
c3c4_solve_jj = matlabFunction(simplify(c3c4_soln_jj(1))); 
c3c4_solve_jc = matlabFunction(simplify(c3c4_soln_jc(1))); 
c3c4_solve_cj = matlabFunction(simplify(c3c4_soln_cj(2)));
c3c4_solve_cc = matlabFunction(simplify(c3c4_soln_cc(2)));

% Clean up workspace
symObj = syms;
cellfun(@clear,symObj);
clear symObj;
warning('on');

clearvars -except currdir workdir datadir resultsdir pathway_opt...
    c3c4_solve_jj c3c4_solve_jc c3c4_solve_cj c3c4_solve_cc;

% Create structure to hold all remaining model output
ss = workspace2struct_fun();

end

%% Solutions for NADP-ME C4 net CO2 assimilation in bundle sheath (As) 

if strcmp(pathway_opt,'NADP-ME-C4') == 1

%% (1) Cyt b6f-limitation M & BS ---------------------------------------- %
eq1 = Ag_sj - Rd_s - An_sj == 0;
eq2 = subs(eq1, Ag_sj, Vc_sj - Vo_sj./2);
eq3 = solve(eq2, An_sj);
eq4 = eq3 - An_sj == 0;
eq5 = subs(eq4, Vo_sj, Vc_sj.*O_s./(S.*C_s));
eq6 = solve(eq5, An_sj);
eq7 = eq6 - An_sj == 0;
eq8 = subs(eq7, Vc_sj, (JP680_sj+2.*Vp_j)./(4.*(1+O_s./(S.*C_s)))); 
eq9 = solve(eq8, An_sj);
eq10 = eq9 - An_sj == 0;
eq11 = subs(eq10, JP680_sj, (JP700_sj - (Vp_j./2).*...
    (3+3.5.*O_s./(S.*C_s))./((4+4.*O_s./(S.*C_s)).*nc))./...
    (1-(nl./nc)+(3+3.5.*O_s./(S.*C_s))./((4+4.*O_s./(S.*C_s)).*nc)));
eq12 = simplify(solve(eq11, An_sj));
eq13 = eq12 - An_sj == 0;
eq14 = subs(eq13, C_s, C_m + (Vp_j - An_sj).*P./gbs);
eq15 = subs(eq14, O_s, O_m + (An_sj - Vp_j./2).*P./gbso);
eq16 = subs(eq15, Vp_j, JP700_mj./(2.*(1 -(nl./nc)+1./nc)));
c4_soln_jj = solve(eq16, An_sj,'MaxDegree',3);
clear eq1 eq2 eq3 eq4 eq5 eq6 eq7 eq8 eq9 eq10... 
      eq11 eq12 eq13 eq14 eq15 eq16;
% @(C_m,JP700_mj,JP700_sj,O_m,P,Rd_s,S,gbs,gbso,nc,nl)

%% (2) Cyt b6f-limited M, Rubiso-limited BS ----------------------------- %
eq1 = Ag_sc - Rd_s - An_sc == 0;
eq2 = subs(eq1, Ag_sc, Vc_sc - Vo_sc./2);
eq3 = solve(eq2, An_sc);
eq4 = eq3 - An_sc == 0;
eq5 = subs(eq4, Vo_sc, Vc_sc.*O_s./(S.*C_s));
eq6 = solve(eq5, An_sc);
eq7 = eq6 - An_sc == 0;
eq8 = subs(eq7, Vc_sc, C_s.*Vcmax_s./(C_s + Kc.*(1+O_s./Ko)));
eq9 = solve(eq8, An_sc);
eq10 = eq9 - An_sc == 0;
eq11 = subs(eq10, C_s, C_m + (Vp_j - An_sc).*P./gbs);
eq12 = solve(eq11, An_sc,'PrincipalValue',true);
eq13 = eq12 - An_sc == 0;
eq14 = simplify(solve(eq13, An_sc));
eq15 = eq14 - An_sc == 0;
eq16 = subs(eq15, O_s, O_m + (An_sc - Vp_j./2).*P./gbso); 
eq17 = simplify(solve(eq16, An_sc));
eq18 = eq17(2) - An_sc == 0;
eq19 = subs(eq18, Vp_j, JP700_mj./(2.*(1 -(nl./nc)+1./nc))); 
c4_soln_jc = solve(eq19, An_sc); 
clear eq1 eq2 eq3 eq4 eq5 eq6 eq7 eq8 eq9 eq10... 
      eq11 eq12 eq13 eq14 eq15 eq16 eq17 eq18 eq19;
% @(C_m,JP700_mj,Kc,Ko,O_m,P,Rd_s,S,Vcmax_s,gbs,gbso,nc,nl)

%% (3) Rubisco-limitation M & Cyt b6f-limitation BS --------------------- %
eq1 = Ag_sj - Rd_s - An_sj == 0;
eq2 = subs(eq1, Ag_sj, Vc_sj - Vo_sj./2);
eq3 = solve(eq2, An_sj);
eq4 = eq3 - An_sj == 0;
eq5 = subs(eq4, Vo_sj, Vc_sj.*O_s./(S.*C_s));
eq6 = simplify(solve(eq5, An_sj));
eq7 = eq6 - An_sj == 0;
eq8 = subs(eq7, Vc_sj, (JP680_sj+2.*Vp_c)./(4.*(1+O_s./(S.*C_s))));
eq9 = subs(eq8, JP680_sj, (JP700_sj - (Vp_c./2).*...
    (3+3.5.*O_s./(S.*C_s))./((4+4.*O_s./(S.*C_s)).*nc))./...
    (1-(nl./nc)+(3+3.5.*O_s./(S.*C_s))./((4+4.*O_s./(S.*C_s)).*nc)));
eq10 = simplify(solve(eq9, An_sj));
eq11 = eq10 - An_sj == 0;
eq12 = subs(eq11, C_s, C_m + (Vp_c - An_sj).*P./gbs);
eq13 = subs(eq12, O_s, O_m + (An_sj - Vp_c./2).*P./gbso);  
eq14 = subs(eq13, Vp_c, Vpmax_m.*C_m./(Kp + C_m));
c4_soln_cj = solve(eq14, An_sj, 'MaxDegree',3);
clear eq1 eq2 eq3 eq4 eq5 eq6 eq7 eq8 eq9 eq10... 
      eq11 eq12 eq13 eq14 eq15 eq16;
% @(C_m,JP700_sj,Kp,O_m,P,Rd_s,S,Vpmax_m,gbs,gbso,nc,nl)

%% (4) Rubisco-limited M & BS ------------------------------------------- %
eq1 = Ag_sc - Rd_s - An_sc == 0;
eq2 = subs(eq1, Ag_sc, Vc_sc - Vo_sc./2);
eq3 = solve(eq2, An_sc);
eq4 = eq3 - An_sc == 0;
eq5 = subs(eq4, Vo_sc, Vc_sc.*O_s./(S.*C_s));
eq6 = solve(eq5, An_sc);
eq7 = eq6 - An_sc == 0;
eq8 = subs(eq7, Vc_sc, C_s.*Vcmax_s./(C_s + Kc.*(1+O_s./Ko)));
eq9 = solve(eq8, An_sc);
eq10 = eq9 - An_sc == 0;
eq11 = subs(eq10, C_s, C_m + (Vp_c - An_sc).*P./gbs);
eq12 = simplify(solve(eq11, An_sc,'PrincipalValue',true));
eq13 = eq12 - An_sc == 0;
eq14 = subs(eq13, O_s, O_m + (An_sc - Vp_c./2).*P./gbso);
eq15 = solve(eq14, An_sc);
eq16 = eq15(2) - An_sc == 0;
eq17 = subs(eq16, Vp_c, Vpmax_m.*C_m./(Kp + C_m));
c4_soln_cc = solve(eq17, An_sc);
clear eq1 eq2 eq3 eq4 eq5 eq6 eq7 eq8 eq9 eq10...
      eq11 eq12 eq13 eq14 eq15 eq16 eq17;
% @(C_m,Kc,Ko,Kp,O_m,P,Rd_s,S,Vcmax_s,Vpmax_m,gbs,gbso)

%% Export NADP-ME C4 functions:
c4_solve_jj = matlabFunction(real(c4_soln_jj(2))); 
c4_solve_jc = matlabFunction(simplify(c4_soln_jc(1)));
c4_solve_cj = matlabFunction(real(c4_soln_cj(3))); 
c4_solve_cc = matlabFunction(simplify(c4_soln_cc(1)));

% Clean up workspace
symObj = syms;
cellfun(@clear,symObj);
clear symObj;
warning('on');

clearvars -except currdir workdir datadir resultsdir pathway_opt...
    c4_solve_jj c4_solve_jc c4_solve_cj c4_solve_cc;

% Create structure to hold all remaining model output
ss = workspace2struct_fun();

end

end