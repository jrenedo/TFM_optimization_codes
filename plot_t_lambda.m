%% With the results from the study of F_opt and E_opt:

clc
clear all
close all

set(0,'DefaultTextFontsize',8,'DefaultTextFontname','Times New Roman','DefaultAxesFontname','Times New Roman')
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');


plotheight=8*1.5;
plotwidth=14*1.5;

subplotsx=2;
subplotsy=2;   
leftedge=1;
rightedge=0.4;   
topedge=0.4;
bottomedge=1;
spacex=1;
spacey=1.3;
fontsize=10; 
linewidth_1 = 1;

%% 



n_dm_g = 4;


load('results_cores_tlambda.mat')
vec_t_lambda=[0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6];
dm_g=0.8;

[N_cores, N_t_lambda] = meshgrid(vec_cores, vec_t_lambda);
Mat_K_opt = 0*N_cores;
Mat_K_B = 0*N_cores;


f=figure('visible','on');
%clf(f);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
set(gcf,'Units','centimeters','Position',[5 5 plotwidth plotheight])
% contourf(N_cores,Kb,100*N_c_gaps)
% contourf(N_cores,Kb,N_c_gaps,'ShowText','on')
%contourf(N_cores,Kb,N_t_lambda,'ShowText','on','LineWidth',2)
surf(N_cores,N_t_lambda,Kb)
% xlim([0 80])
% ylim([0 1])

xlabel('$N_{cores}$','Interpreter','latex')
ylabel('$t/\tau$','Interpreter','latex')
zlabel('$K_B$','Interpreter','latex')
colorbar
view(-50,-40)
title('$d_m/g=0.8$','Interpreter','latex')

%print('K_B_tlambda','-dpng','-r300')
print('K_B_tlambda','-depsc')

f=figure('visible','on');
%clf(f);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
set(gcf,'Units','centimeters','Position',[5 5 plotwidth plotheight])
% contourf(N_cores,Kb,100*N_c_gaps)
% contourf(N_cores,Kb,N_c_gaps,'ShowText','on')
%contourf(N_cores,Kb,N_t_lambda,'ShowText','on','LineWidth',2)
surf(N_cores,N_t_lambda,K_opt)
% xlim([0 80])
% ylim([0 1])

xlabel('$N_{cores}$','Interpreter','latex')
ylabel('$t/\tau$','Interpreter','latex')
zlabel('$f_{opt}$','Interpreter','latex')
colorbar
view(50,20)
title('$d_m/g=0.8$','Interpreter','latex')

%print('f_opt_tlambda','-dpng','-r300')
print('f_opt_tlambda','-depsc')


return
