%% Header
%
%  Citation:
%  Johnson, JE, Field, CB, Berry, JA (2021) The limiting factors 
%  and regulatory processes that control the environmental responses of C3, 
%  C3-C4 intermediate, and C4 photosynthesis. Oecologia, 
%  https://doi.org/10.1007/s00442-021-05062-y
%
%  Description:
%   Example simulation of light response with forward model
%
%  Compatibile with:
%   MATLAB R2020b
%
%  Last revised: 
%   2021-10-18
%
%% Set up environment

% Clean up working environment
clear all; % variables
close all; % figures
clc; % command window
tic; % start timer

% Set up subdirectories
currdir = pwd;
workdir = fullfile(currdir,'/scripts');
resultsdir = fullfile(currdir,'/outputs');

%% Select inputs

% Choose pathway for simulation
cd(workdir)
%pathway_opt = 'C3';
%pathway_opt = 'Type-I-C3-C4'; 
pathway_opt = 'NADP-ME-C4';

% Create output directory
cd(resultsdir);
outputname = strcat('Example1-Light-',pathway_opt);
mkdir(outputname);
cd(workdir);

% Specify environmental conditions
n = 100;                                   % Steps in vector
data.Qin = transpose(linspace(1,2400,n));   % PAR, umol PPFD m-2 s-1
data.Tin = repmat(25,n,1);                  % Leaf temperature, C
data.Cin = repmat(250,n,1);                 % Mesophyll CO2, ubar
data.Oin = repmat(0.209,n,1);               % Atmospheric O2, bar
data.Pin = repmat(1,n,1);                   % Total pressure, bar

%% Configure model simulations

% Set all parameters
v = configure_fun(data);
v.Ku2 = 2e09;

% Select subset of parameters to fit

if strcmp(pathway_opt,'C3') == 1
    % Case 1: Pure C3  
    v.CB6F = 175./v.kq.*1e-06; 
    v.RUB = 50./v.kc.*1e-06; 
    v.abs_frac = 0.0;
    v.vq_frac = 0.0;
    v.vc_frac = 0.0;
    v.a2_s_frac = 0.0;
    v.a2_m_frac = 0.52;
end

if strcmp(pathway_opt,'Type-I-C3-C4') == 1
    % Case 2: Type I C3-C4 w/ "neutral" shuttle 
    ss = symsolver_c3c4_fun(pathway_opt);
    v.CB6F = 175./v.kq.*1e-06; 
    v.RUB = 50./v.kc.*1e-06; 
    v.abs_frac = 0.05;
    v.vq_frac = 0.05;
    v.vc_frac = 0.1;
    v.a2_s_frac = 0.52;
    v.a2_m_frac = 0.52;
end
     
if strcmp(pathway_opt,'NADP-ME-C4') == 1
    % Case 3: NADP-ME-C4 w/ "neutral" shuttle  
    ss = symsolver_c3c4_fun(pathway_opt);
    v.CB6F = 175./v.kq.*1e-06; 
    v.RUB = 30./v.kc.*1e-06;
    v.Vpmax = v.RUB.*v.kc.*2;
    v.abs_frac = 0.4;
    v.vq_frac = 0.4;
    v.vc_frac = 1;
    v.a2_s_frac = 0.47;
    v.a2_m_frac = 0.47;
end

%% Run simulation and visualize results

% % Update objective function (optional)
v.Model_id = {'model_fun_c3c4'};
v.pathway_opt = pathway_opt;

% Assign functions to slots in 'v'
if strcmp(pathway_opt,'Type-I-C3-C4') == 1
v.c3c4_solve_cc = ss.c3c4_solve_cc;
v.c3c4_solve_cj = ss.c3c4_solve_cj;
v.c3c4_solve_jc = ss.c3c4_solve_jc;
v.c3c4_solve_jj = ss.c3c4_solve_jj;
end
if strcmp(pathway_opt,'NADP-ME-C4') == 1
v.c4_solve_cc = ss.c4_solve_cc;
v.c4_solve_cj = ss.c4_solve_cj;
v.c4_solve_jc = ss.c4_solve_jc;
v.c4_solve_jj = ss.c4_solve_jj;
end

% Run simulation
m = model_fun_c3c4(v);

 % Plot results
[s] = plotter_forward_fun_c3c4(outputname,v,m);

%% Save results
cd(workdir)

% Save figures in graphic format (.png)
set(s.figure1,'PaperPositionMode','auto');
set(s.figure1,'PaperUnits','inches','PaperPosition',[0 0 12 12])
print(s.figure1,fullfile(resultsdir,outputname,...
    strcat(outputname,'-figure.png')),'-dpng','-r300');

% Save results in Matlab file format (.mat)
save(fullfile(resultsdir,outputname,...
    strcat(outputname,'-modelinputs.mat')),'v');
save(fullfile(resultsdir,outputname,...
    strcat(outputname,'-modeloutputs.mat')),'m'); 

% Return to project directory
cd(currdir);

toc;
