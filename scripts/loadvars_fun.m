function [pathway_opt, Q, T, P, O_m, C_m, ...
Abs, abs_frac, a2_m_frac, a2_s_frac, eps1, eps2, ...
CB6F, vq_frac, RUB, Rdsc, vc_frac, Vpmax, gbs, gbso, ...
Kf, Kd, Kp1, Kn1, Kp2, Ku2, kq, nl, nc, kc, ko, Kc, Ko, Kp, ...
c3c4_solve_cc, c3c4_solve_cj, c3c4_solve_jc, c3c4_solve_jj, ...
c4_solve_cc, c4_solve_cj, c4_solve_jc, c4_solve_jj] = loadvars_fun(v)

% Rename default values of parameters passed in from 'v' structure

pathway_opt = v.pathway_opt;
Q = v.Q;
T = v.T;
P = v.P;
O_m = v.O_m;
C_m = v.C_m;
Abs = v.Abs;
abs_frac = v.abs_frac;
a2_m_frac = v.a2_m_frac;
a2_s_frac = v.a2_s_frac;
vq_frac = v.vq_frac;
vc_frac = v.vc_frac;
CB6F = v.CB6F;
RUB = v.RUB;
Rdsc = v.Rdsc;
Vpmax = v.Vpmax;
gbs = v.gbs;
gbso = v.gbso;
Kf = v.Kf;
Kd = v.Kd;
Kp1 = v.Kp1;
Kp2 = v.Kp2;
Kn1 = v.Kn1; 
Ku2 = v.Ku2;
kq = v.kq;
nl = v.nl;
nc = v.nc;
kc = v.kc;
ko = v.ko;
Kc = v.Kc;
Ko = v.Ko;
Kp = v.Kp;
eps1 = v.eps1;
eps2 = v.eps2;
c3c4_solve_cc = v.c3c4_solve_cc;
c3c4_solve_cj = v.c3c4_solve_cj;
c3c4_solve_jc = v.c3c4_solve_jc;
c3c4_solve_jj = v.c3c4_solve_jj;
c4_solve_cc = v.c4_solve_cc;
c4_solve_cj = v.c4_solve_cj;
c4_solve_jc = v.c4_solve_jc;
c4_solve_jj = v.c4_solve_jj;

end