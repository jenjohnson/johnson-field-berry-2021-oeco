function [s] = plotter_forward_fun_c3c4(outputname,v,m)

% Dock figures
set(0,'DefaultFigureWindowStyle','docked');

%% Make figure

if range(m.Q)>0
    x = m.Q.*1e6;
    xlab = 'PAR (umol PPFD m-2 s-1)';
end

if range(m.C_m)>0
    x = m.C_m.*1e6;
    xlab = 'Cm (ubar CO2)';
end

if range(m.T)>0
    x = m.T;
    xlab = 'Tleaf (C)';
end

if range(m.Q)>0 && range(m.T)>0
    x = linspace(1,1440,1440)./60;
    xlab = 'Time (hours)';
end
    
figure1 = figure;
xpos = 0.85;
ypos = 0.9;
 
subplot(5,4,1);
plot(x, m.JP680_ma.*1e6,'-b')
hold on
plot(x, m.JP680_sa.*1e6,'-r')
hold off
xlim([min(x) max(x)]);
ylim([0 250]);
ylabel('PS2 ETR (umol e- m-2 s-1)');
legend({'Mesophyll','Bundle sheath'});
xlim_curr = get(gca,'xlim');
ylim_curr = get(gca,'ylim');
text(xlim_curr(1) + (xlim_curr(2)-xlim_curr(1)).*xpos,...
     ylim_curr(1) + (ylim_curr(2)-ylim_curr(1)).*ypos,...
     '(a)')

subplot(5,4,2)
plot(x, m.phi2P_ma.*100,'-b')
hold on
plot(x, m.phi2P_sa.*100,'-r')
hold off
xlim([min(x) max(x)]);
ylim([0 100]);
ylabel('Phi2P (%)');
xlim_curr = get(gca,'xlim');
ylim_curr = get(gca,'ylim');
text(xlim_curr(1) + (xlim_curr(2)-xlim_curr(1)).*xpos,...
     ylim_curr(1) + (ylim_curr(2)-ylim_curr(1)).*ypos,...
     '(b)')
 
subplot(5,4,3)
plot(x, m.phi2N_ma.*100,'-b')
hold on
plot(x, m.phi2N_sa.*100,'-r')
hold off
xlim([min(x) max(x)]);
ylim([0 100]);
ylabel('Phi2N (%)');
xlim_curr = get(gca,'xlim');
ylim_curr = get(gca,'ylim');
text(xlim_curr(1) + (xlim_curr(2)-xlim_curr(1)).*xpos,...
     ylim_curr(1) + (ylim_curr(2)-ylim_curr(1)).*ypos,...
     '(c)')
 
subplot(5,4,4)
plot(x, (m.phi2D_ma + m.phi2F_ma).*100,'-b')
hold on
plot(x, (m.phi2D_sa + m.phi2F_sa).*100,'-r')
hold off
xlim([min(x) max(x)]);
ylim([0 100]);
ylabel('Phi2DF (%)');
xlim_curr = get(gca,'xlim');
ylim_curr = get(gca,'ylim');
text(xlim_curr(1) + (xlim_curr(2)-xlim_curr(1)).*xpos,...
     ylim_curr(1) + (ylim_curr(2)-ylim_curr(1)).*ypos,...
     '(d)')
 
subplot(5,4,5)
plot(x, m.JP700_ma.*1e6,'-b')
hold on
plot(x, m.JP700_sa.*1e6,'-r')
hold off
xlim([min(x) max(x)]);
ylim([0 250]);
ylabel('PS1 ETR (umol e- m-2 s-1)');
xlim_curr = get(gca,'xlim');
ylim_curr = get(gca,'ylim');
text(xlim_curr(1) + (xlim_curr(2)-xlim_curr(1)).*xpos,...
     ylim_curr(1) + (ylim_curr(2)-ylim_curr(1)).*ypos,...
     '(e)')
 
subplot(5,4,6)
plot(x, m.phi1P_ma.*100,'-b')
hold on
plot(x, m.phi1P_sa.*100,'-r')
hold off
xlim([min(x) max(x)]);
ylim([0 100]);
ylabel('Phi1P (%)');
xlim_curr = get(gca,'xlim');
ylim_curr = get(gca,'ylim');
text(xlim_curr(1) + (xlim_curr(2)-xlim_curr(1)).*xpos,...
     ylim_curr(1) + (ylim_curr(2)-ylim_curr(1)).*ypos,...
     '(f)')
 
subplot(5,4,7)
plot(x, m.phi1N_ma.*100,'-b')
hold on
plot(x, m.phi1N_sa.*100,'-r')
hold off
xlim([min(x) max(x)]);
ylim([0 100]);
ylabel('Phi1N (%)');
xlim_curr = get(gca,'xlim');
ylim_curr = get(gca,'ylim');
text(xlim_curr(1) + (xlim_curr(2)-xlim_curr(1)).*xpos,...
     ylim_curr(1) + (ylim_curr(2)-ylim_curr(1)).*ypos,...
     '(g)')
 
subplot(5,4,8)
plot(x, (m.phi1D_ma + m.phi1F_ma).*100,'-b')
hold on
plot(x, (m.phi1D_sa + m.phi1F_sa).*100,'-r')
hold off
xlim([min(x) max(x)]);
ylim([0 100]);
ylabel('Phi1DF (%)');
xlim_curr = get(gca,'xlim');
ylim_curr = get(gca,'ylim');
text(xlim_curr(1) + (xlim_curr(2)-xlim_curr(1)).*xpos,...
     ylim_curr(1) + (ylim_curr(2)-ylim_curr(1)).*ypos,...
     '(h)')
 
subplot(5,4,9)
plot(x, m.JP700_ma.*1e6,'-b')
hold on
plot(x, m.JP700_sa.*1e6,'-r')
hold off
xlim([min(x) max(x)]);
ylim([0 250]);
ylabel('Cyt b6f ETR (umol e- m-2 s-1)');
xlim_curr = get(gca,'xlim');
ylim_curr = get(gca,'ylim');
text(xlim_curr(1) + (xlim_curr(2)-xlim_curr(1)).*xpos,...
     ylim_curr(1) + (ylim_curr(2)-ylim_curr(1)).*ypos,...
     '(i)')
 
subplot(5,4,10) 
plot(x, (1-m.q2_ma).*100,'-b')
hold on
plot(x, (1-m.q2_sa).*100,'-r')
hold off
xlim([min(x) max(x)]);
ylim([0 100]);
ylabel('PQH2/[PQ+PQH2] (%)');
xlim_curr = get(gca,'xlim');
ylim_curr = get(gca,'ylim');
text(xlim_curr(1) + (xlim_curr(2)-xlim_curr(1)).*xpos,...
     ylim_curr(1) + (ylim_curr(2)-ylim_curr(1)).*ypos,...
     '(j)')
 
subplot(5,4,11)
plot(x, m.JP700_ma./m.JP700_mj.*(m.Vqmax_m./m.CB6F_m),'-b')
hold on
plot(x, m.JP700_sa./...
    (m.JP700_sjj.*logical(m.which_JP700_ma == 1) + ...
     m.JP700_scj.*logical(m.which_JP700_ma == 2))...
     .*(m.Vqmax_s./m.CB6F_s),'-r')
hold off
xlim([min(x) max(x)]);
ylim([0 600]);
ylabel('Cyt b6f rate constant (s-1)');
xlim_curr = get(gca,'xlim');
ylim_curr = get(gca,'ylim');
text(xlim_curr(1) + (xlim_curr(2)-xlim_curr(1)).*xpos,...
     ylim_curr(1) + (ylim_curr(2)-ylim_curr(1)).*ypos,...
     '(k)')
 
subplot(5,4,12) 
plot(x, (1-m.JP700_ma./m.JP700_mj).*100,'-b')
hold on
plot(x, (1-m.JP700_sa./...
    (m.JP700_sjj.*logical(m.which_JP700_ma == 1) + ...
     m.JP700_scj.*logical(m.which_JP700_ma == 2))).*100,'-r')
hold off
xlim([min(x) max(x)]);
ylim([0 100]);
ylabel('Source downregulation (%)');
xlim_curr = get(gca,'xlim');
ylim_curr = get(gca,'ylim');
text(xlim_curr(1) + (xlim_curr(2)-xlim_curr(1)).*xpos,...
     ylim_curr(1) + (ylim_curr(2)-ylim_curr(1)).*ypos,...
     '(l)')
 
subplot(5,4,13)
plot(x, m.An_ma.*1e6,'-b')
hold on
plot(x, m.An_sa.*1e6,'-r')
hold off
xlim([min(x) max(x)]);
ylim([0 30]);
ylabel('An (umol CO2 m-2 s-1)');
xlim_curr = get(gca,'xlim');
ylim_curr = get(gca,'ylim');
text(xlim_curr(1) + (xlim_curr(2)-xlim_curr(1)).*xpos,...
     ylim_curr(1) + (ylim_curr(2)-ylim_curr(1)).*ypos,...
     '(m)')
 
subplot(5,4,14)
plot(x, m.C_m./m.O_m.*1000,'-b')
hold on
plot(x, m.C_sa./m.O_sa.*1000,'-r')
hold off
xlim([min(x) max(x)]);
ylim([0 30]);
ylabel('CO2:O2 (mbar bar-1)');
xlim_curr = get(gca,'xlim');
ylim_curr = get(gca,'ylim');
text(xlim_curr(1) + (xlim_curr(2)-xlim_curr(1)).*xpos,...
     ylim_curr(1) + (ylim_curr(2)-ylim_curr(1)).*ypos,...
     '(n)')
 
subplot(5,4,15) 
plot(x, m.JP700_ma./m.JP700_mc.*m.RUB_m.*1e6,'-b')
hold on
plot(x, m.JP700_sa./...
    (m.JP700_sjc.*logical(m.which_JP700_ma == 1) + ...
     m.JP700_scc.*logical(m.which_JP700_ma == 2))...
     .*m.RUB_s.*1e6,'-r')
hold off
xlim([min(x) max(x)]);
ylim([0 15]);
ylabel('Rub active sites (umol m-2)');
xlim_curr = get(gca,'xlim');
ylim_curr = get(gca,'ylim');
text(xlim_curr(1) + (xlim_curr(2)-xlim_curr(1)).*xpos,...
     ylim_curr(1) + (ylim_curr(2)-ylim_curr(1)).*ypos,...
     '(o)')
 
subplot(5,4,16) 
plot(x, (1-m.JP700_ma./m.JP700_mc).*100,'-b')
hold on
plot(x, (1-m.JP700_sa./...
    (m.JP700_sjc.*logical(m.which_JP700_ma == 1) + ...
     m.JP700_scc.*logical(m.which_JP700_ma == 2))).*100,'-r')
hold off
xlim([min(x) max(x)]);
ylim([0 100]);
ylabel('Sink downregulation (%)');
xlim_curr = get(gca,'xlim');
ylim_curr = get(gca,'ylim');
text(xlim_curr(1) + (xlim_curr(2)-xlim_curr(1)).*xpos,...
     ylim_curr(1) + (ylim_curr(2)-ylim_curr(1)).*ypos,...
     '(p)')

subplot(5,4,17) 
plot(x, m.which_JP700_ma,'-b')
hold on
plot(x, (m.which_JP700_sj.*logical(m.which_JP700_ma == 1) + ... % _jj v _jc
         m.which_JP700_sc.*logical(m.which_JP700_ma == 2)),'-r') % _cj v _cc
hold off
xlim([min(x) max(x)]);
ylim([0 3]);
xlabel(xlab);
ylabel('Limiting states');
xlim_curr = get(gca,'xlim');
ylim_curr = get(gca,'ylim');
text(xlim_curr(1) + (xlim_curr(2)-xlim_curr(1)).*0.5,0.66,...
    '1: Light-limited','HorizontalAlignment','center');
text(xlim_curr(1) + (xlim_curr(2)-xlim_curr(1)).*0.5,0.33,...
    '2: Light-saturated','HorizontalAlignment','center');
text(xlim_curr(1) + (xlim_curr(2)-xlim_curr(1)).*xpos,...
     ylim_curr(1) + (ylim_curr(2)-ylim_curr(1)).*ypos,...
     '(q)')
 
subplot(5,4,18)
plot(x, m.Kn2_ma./1e9,'-b')
hold on
plot(x, m.Kn2_sa./1e9,'-r')
hold off
xlim([min(x) max(x)]);
ylim([0 6]);
xlabel(xlab);
ylabel('Kn2 (ns-1)');
xlim_curr = get(gca,'xlim');
ylim_curr = get(gca,'ylim');
text(xlim_curr(1) + (xlim_curr(2)-xlim_curr(1)).*xpos,...
     ylim_curr(1) + (ylim_curr(2)-ylim_curr(1)).*ypos,...
     '(r)')    

subplot(5,4,19)
plot(x, (m.JP700_ma-m.JP680_ma)./m.JP700_ma.*100,'-b')
hold on
plot(x, (m.JP700_sa-m.JP680_sa)./m.JP700_sa.*100,'-r')
hold off
xlim([min(x) max(x)]);
ylim([0 100]);
xlabel(xlab);
ylabel('CEF1 fraction (%)');
xlim_curr = get(gca,'xlim');
ylim_curr = get(gca,'ylim');
text(xlim_curr(1) + (xlim_curr(2)-xlim_curr(1)).*xpos,...
     ylim_curr(1) + (ylim_curr(2)-ylim_curr(1)).*ypos,...
     '(s)')

subplot(5,4,20) 
if strcmp(v.pathway_opt,'C3')
    plot(x,zeros(length(x),1),'-r')
else
plot(x, (m.L_C_sa)./(m.Vp_ma+m.Vg_ma).*100,'-r');
end
hold off
xlim([min(x) max(x)]);
ylim([0 100]);
xlabel(xlab);
ylabel('L/(Vp+Vg) (%)');
xlim_curr = get(gca,'xlim');
ylim_curr = get(gca,'ylim');
text(xlim_curr(1) + (xlim_curr(2)-xlim_curr(1)).*xpos,...
     ylim_curr(1) + (ylim_curr(2)-ylim_curr(1)).*ypos,...
     '(t)')

% % Add title to plots
%sgtitle(figure1,outputname);
 
% Remove model input variables from workspace
clearvars -except figure1;

% Create structure to hold all remaining model output
s = workspace2struct_fun();

end