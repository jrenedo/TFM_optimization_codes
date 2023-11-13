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


load('results_cores_a2')
vec_gaps_pc=[0.1 0.5 0.75 1 1.5 2]

[N_cores, N_c_gaps] = meshgrid(vec_cores, vec_gaps_pc);
Mat_K_opt = 0*N_cores;
Mat_K_B = 0*N_cores;


%figure
% contourf(N_cores,Kb,100*N_c_gaps)
% contourf(N_cores,Kb,N_c_gaps,'ShowText','on')
%contourf(N_cores,Kb,N_c_gaps,'ShowText','on','LineWidth',linewidth_1)
% xlim([0 80])
% ylim([0 1])

% return

%% Plot F_opt:
F_opt_max = 15;
%setting the Matlab figure
f=figure('visible','on');
%clf(f);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
set(gcf,'Units','centimeters','Position',[5 5 plotwidth plotheight])

load('results_cores_a2')
vec_gaps_pc=[0.1 0.5 0.75 1 1.5 2]
v = [0.1, 0.2, 0.5, 1, 1.5, 2];
%v = [0.2, 0.5, 1, 1.5, 2];


[N_cores, N_c_gaps] = meshgrid(vec_cores, vec_gaps_pc);
Mat_K_opt = 0*N_cores;
Mat_K_B = 0*N_cores;


grey=[0.4,0.4,0.4];

subplot(2,2,1)
% contourf(N_cores,K_opt,N_c_gaps,'ShowText','on')
% contour(N_cores,K_opt,N_c_gaps,'ShowText','on','LineWidth',2)
[C,h] =  contourf(N_cores,K_opt,N_c_gaps,'ShowText','on','LineWidth',linewidth_1);
clabel(C,h,v)
%clabel(C)

%xlabel('$N_{cores}$','Interpreter','latex')
ylabel('$f_{opt}$','Interpreter','latex')

title('$d_m/g=0.6$','Interpreter','latex')
ylim([0 F_opt_max])

load('results_cores_b2')
vec_gaps_pc=[0.1 0.5 0.75 1 1.5 2]

[N_cores, N_c_gaps] = meshgrid(vec_cores, vec_gaps_pc);
Mat_K_opt = 0*N_cores;
Mat_K_B = 0*N_cores;

subplot(2,2,2)
% contourf(N_cores,K_opt,N_c_gaps,'ShowText','on')
%contourf(N_cores,K_opt,N_c_gaps,'ShowText','on','LineWidth',2)
[C,h] =  contourf(N_cores,K_opt,N_c_gaps,'ShowText','on','LineWidth',linewidth_1);
clabel(C,h,v)
%clabel(C)
%xlabel('$N_{cores}$','Interpreter','latex')
%ylabel('$f_{opt}$','Interpreter','latex')

title('$d_m/g=0.7$','Interpreter','latex')
ylim([0 F_opt_max])

load('results_cores_c2')
vec_gaps_pc=[0.1 0.5 0.75 1 1.5 2]

[N_cores, N_c_gaps] = meshgrid(vec_cores, vec_gaps_pc);
Mat_K_opt = 0*N_cores;
Mat_K_B = 0*N_cores;

subplot(2,2,3)
% contourf(N_cores,K_opt,N_c_gaps,'ShowText','on')
% contour(N_cores,K_opt,N_c_gaps,'ShowText','on','LineWidth',2)
[C,h] =  contourf(N_cores,K_opt,N_c_gaps,'ShowText','on','LineWidth',linewidth_1);
clabel(C,h,v)
%clabel(C)
xlabel('$N_{cores}$','Interpreter','latex')
ylabel('$f_{opt}$','Interpreter','latex')

title('$d_m/g=0.8$','Interpreter','latex')
ylim([0 F_opt_max])

load('results_cores_d2')
vec_gaps_pc=[0.1 0.5 0.75 1 1.5 2]

[N_cores, N_c_gaps] = meshgrid(vec_cores, vec_gaps_pc);
Mat_K_opt = 0*N_cores;
Mat_K_B = 0*N_cores;

subplot(2,2,4)
% contourf(N_cores,K_opt,N_c_gaps,'ShowText','on')
% contour(N_cores,K_opt,N_c_gaps,'ShowText','on','LineWidth',2)
[C,h] =  contourf(N_cores,K_opt,N_c_gaps,'ShowText','on','LineWidth',linewidth_1);
clabel(C,h,v)
%clabel(C)
xlabel('$N_{cores}$','Interpreter','latex')
%ylabel('$f_{opt}$','Interpreter','latex')

title('$d_m/g=0.9$','Interpreter','latex')
ylim([0 F_opt_max])
%print('f_opt','-dpng','-r300')
print('f_opt','-depsc')

%% Plot K_B:
K_B_max = 0.8;
%setting the Matlab figure
f=figure('visible','on');
%clf(f);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
set(gcf,'Units','centimeters','Position',[5 5 plotwidth plotheight])

load('results_cores_a2')
vec_gaps_pc=[0.1 0.5 0.75 1 1.5 2]
v = [0.1, 0.2, 0.5, 1, 1.5, 2];
%v = [0.2, 0.5, 1, 1.5, 2];


[N_cores, N_c_gaps] = meshgrid(vec_cores, vec_gaps_pc);
Mat_K_opt = 0*N_cores;
Mat_K_B = 0*N_cores;


grey=[0.4,0.4,0.4];

subplot(2,2,1)
% contourf(N_cores,K_opt,N_c_gaps,'ShowText','on')
% contour(N_cores,K_opt,N_c_gaps,'ShowText','on','LineWidth',2)
[C,h] =  contourf(N_cores,Kb,N_c_gaps,'ShowText','on','LineWidth',linewidth_1);
clabel(C,h,v)
%clabel(C)

%xlabel('$N_{cores}$','Interpreter','latex')
ylabel('$K_B$','Interpreter','latex')

title('$d_m/g=0.6$','Interpreter','latex')
ylim([0 K_B_max])

load('results_cores_b2')
vec_gaps_pc=[0.1 0.5 0.75 1 1.5 2]

[N_cores, N_c_gaps] = meshgrid(vec_cores, vec_gaps_pc);
Mat_K_opt = 0*N_cores;
Mat_K_B = 0*N_cores;

subplot(2,2,2)
% contourf(N_cores,K_opt,N_c_gaps,'ShowText','on')
%contourf(N_cores,K_opt,N_c_gaps,'ShowText','on','LineWidth',2)
[C,h] =  contourf(N_cores,Kb,N_c_gaps,'ShowText','on','LineWidth',linewidth_1);
clabel(C,h,v)
%clabel(C)
%xlabel('$N_{cores}$','Interpreter','latex')
%ylabel('$f_{opt}$','Interpreter','latex')

title('$d_m/g=0.7$','Interpreter','latex')
ylim([0 K_B_max])

load('results_cores_c2')
vec_gaps_pc=[0.1 0.5 0.75 1 1.5 2]

[N_cores, N_c_gaps] = meshgrid(vec_cores, vec_gaps_pc);
Mat_K_opt = 0*N_cores;
Mat_K_B = 0*N_cores;

subplot(2,2,3)
% contourf(N_cores,K_opt,N_c_gaps,'ShowText','on')
% contour(N_cores,K_opt,N_c_gaps,'ShowText','on','LineWidth',2)
[C,h] =  contourf(N_cores,Kb,N_c_gaps,'ShowText','on','LineWidth',linewidth_1);
clabel(C,h,v)
%clabel(C)
xlabel('$N_{cores}$','Interpreter','latex')
ylabel('$K_B$','Interpreter','latex')

title('$d_m/g=0.8$','Interpreter','latex')
ylim([0 K_B_max])

load('results_cores_d2')
vec_gaps_pc=[0.1 0.5 0.75 1 1.5 2]

[N_cores, N_c_gaps] = meshgrid(vec_cores, vec_gaps_pc);
Mat_K_opt = 0*N_cores;
Mat_K_B = 0*N_cores;

subplot(2,2,4)
% contourf(N_cores,K_opt,N_c_gaps,'ShowText','on')
% contour(N_cores,K_opt,N_c_gaps,'ShowText','on','LineWidth',2)
[C,h] =  contourf(N_cores,Kb,N_c_gaps,'ShowText','on','LineWidth',linewidth_1);
clabel(C,h,v)
%clabel(C)
xlabel('$N_{cores}$','Interpreter','latex')
%ylabel('$f_{opt}$','Interpreter','latex')

title('$d_m/g=0.9$','Interpreter','latex')
ylim([0 K_B_max])
%print('K_B_plot','-dpng','-r300')
print('K_B_plot','-depsc')

return

%% Plot Kb:
KB_max = 0.8;
%setting the Matlab figure
f=figure('visible','on');
clf(f);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
set(gcf,'Units','centimeters','Position',[5 5 plotwidth plotheight])

load('results_cores_a2')
vec_gaps_pc=[0.1 0.5 0.75 1 1.5 2]

[N_cores, N_c_gaps] = meshgrid(vec_cores, vec_gaps_pc);
Mat_K_opt = 0*N_cores;
Mat_K_B = 0*N_cores;


grey=[0.4,0.4,0.4];

ax=axes('position',sub_pos{1,2},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');

% contourf(N_cores,Kb,N_c_gaps,'ShowText','on')
contour(N_cores,Kb,N_c_gaps,'ShowText','on','LineWidth',2)
xlabel('$N_{cores}$','Interpreter','latex')
ylabel('$K_B$','Interpreter','latex')

title('$d_m/g=0.6$','Interpreter','latex')
ylim([0 KB_max])

load('results_cores_b2')
vec_gaps_pc=[0.1 0.5 0.75 1 1.5 2]

[N_cores, N_c_gaps] = meshgrid(vec_cores, vec_gaps_pc);
Mat_K_opt = 0*N_cores;
Mat_K_B = 0*N_cores;


ax=axes('position',sub_pos{2,2},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
% contourf(N_cores,Kb,N_c_gaps,'ShowText','on')
contour(N_cores,Kb,N_c_gaps,'ShowText','on','LineWidth',2)
xlabel('$N_{cores}$','Interpreter','latex')
ylabel('$K_B$','Interpreter','latex')

title('$d_m/g=0.7$','Interpreter','latex')
ylim([0 KB_max])

load('results_cores_c2')
vec_gaps_pc=[0.1 0.5 0.75 1 1.5 2]

[N_cores, N_c_gaps] = meshgrid(vec_cores, vec_gaps_pc);
Mat_K_opt = 0*N_cores;
Mat_K_B = 0*N_cores;

ax=axes('position',sub_pos{1,1},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
% contourf(N_cores,Kb,N_c_gaps,'ShowText','on')
contour(N_cores,Kb,N_c_gaps,'ShowText','on','LineWidth',2)
xlabel('$N_{cores}$','Interpreter','latex')
ylabel('$K_B$','Interpreter','latex')

title('$d_m/g=0.8$','Interpreter','latex')
ylim([0 KB_max])

load('results_cores_d2')
vec_gaps_pc=[0.1 0.5 0.75 1 1.5 2]

[N_cores, N_c_gaps] = meshgrid(vec_cores, vec_gaps_pc);
Mat_K_opt = 0*N_cores;
Mat_K_B = 0*N_cores;

ax=axes('position',sub_pos{2,1},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
% contourf(N_cores,Kb,N_c_gaps,'ShowText','on')
contour(N_cores,Kb,N_c_gaps,'ShowText','on','LineWidth',2)
xlabel('$N_{cores}$','Interpreter','latex')
ylabel('$K_B$','Interpreter','latex')

title('$d_m/g=0.9$','Interpreter','latex')
ylim([0 KB_max])


return

%% Plot K_B:

%setting the Matlab figure
f=figure('visible','on');
clf(f);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
set(gcf,'Units','centimeters','Position',[5 5 plotwidth plotheight])

load('results_cores_a2')


grey=[0.4,0.4,0.4];

ax=axes('position',sub_pos{1,2},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
contourf(N_cores,Kb,N_c_gaps,'ShowText','on')
xlabel('$N_{cores}$','Interpreter','latex')

ylabel('$K_B$','Interpreter','latex')

title('$d_m/g=0.6$','Interpreter','latex')


load('results_cores_b2')


grey=[0.4,0.4,0.4];

ax=axes('position',sub_pos{2,2},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
contourf(N_cores,Kb,N_c_gaps,'ShowText','on')
xlabel('$N_{cores}$','Interpreter','latex')
ylabel('$K_B$','Interpreter','latex')

title('$d_m/g=0.7$','Interpreter','latex')

load('results_cores_c2')


grey=[0.4,0.4,0.4];

ax=axes('position',sub_pos{1,1},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
contourf(N_cores,Kb,N_c_gaps,'ShowText','on')
xlabel('$N_{cores}$','Interpreter','latex')
ylabel('$K_B$','Interpreter','latex')

title('$d_m/g=0.8$','Interpreter','latex')

load('results_cores_d2')


grey=[0.4,0.4,0.4];

figure
hold on
plot(vec_cores,Kb(1,:),'k-','LineWidth',2)


xlabel('$N_{cores}$','Interpreter','latex')
ylabel('$K_B$','Interpreter','latex')

title('$d_m/g=0.9$','Interpreter','latex')

% size of the figure

% legend('c_{gap}=0.5% of D','c_{gap}=0.75% of D','c_{gap}=1% of D','c_{gap}=1.5% of D','c_{gap}=2% of D')



%% Plot KE_opt:

%setting the Matlab figure
f=figure('visible','on');
clf(f);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
set(gcf,'Units','centimeters','Position',[5 5 plotwidth plotheight])

load('results_cores_a')


grey=[0.4,0.4,0.4];

ax=axes('position',sub_pos{1,2},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
hold on
plot(vec_cores,K_gamma(1,:),'k-','LineWidth',2)
plot(vec_cores,K_gamma(2,:),'k--','LineWidth',2)
plot(vec_cores,K_gamma(3,:),'k-.','LineWidth',2)
plot(vec_cores,K_gamma(4,:),'k:','LineWidth',2)
plot(vec_cores,K_gamma(5,:),'LineWidth',2,'Color',grey)

xlabel('$N_{cores}$','Interpreter','latex')
ylabel('$f_e$','Interpreter','latex')

title('$d_m/g=0.6$','Interpreter','latex')

load('results_cores_b')


grey=[0.4,0.4,0.4];

ax=axes('position',sub_pos{2,2},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
hold on
plot(vec_cores,K_gamma(1,:),'k-','LineWidth',2)
plot(vec_cores,K_gamma(2,:),'k--','LineWidth',2)
plot(vec_cores,K_gamma(3,:),'k-.','LineWidth',2)
plot(vec_cores,K_gamma(4,:),'k:','LineWidth',2)
plot(vec_cores,K_gamma(5,:),'LineWidth',2,'Color',grey)

xlabel('$N_{cores}$','Interpreter','latex')
ylabel('$f_e$','Interpreter','latex')

title('$d_m/g=0.7$','Interpreter','latex')

load('results_cores_c')


grey=[0.4,0.4,0.4];

ax=axes('position',sub_pos{1,1},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
hold on
plot(vec_cores,K_gamma(1,:),'k-','LineWidth',2)
plot(vec_cores,K_gamma(2,:),'k--','LineWidth',2)
plot(vec_cores,K_gamma(3,:),'k-.','LineWidth',2)
plot(vec_cores,K_gamma(4,:),'k:','LineWidth',2)
plot(vec_cores,K_gamma(5,:),'LineWidth',2,'Color',grey)

xlabel('$N_{cores}$','Interpreter','latex')
ylabel('$f_e$','Interpreter','latex')

title('$d_m/g=0.8$','Interpreter','latex')

load('results_cores_d')


grey=[0.4,0.4,0.4];

ax=axes('position',sub_pos{2,1},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
hold on
plot(vec_cores,K_gamma(1,:),'k-','LineWidth',2)
plot(vec_cores,K_gamma(2,:),'k--','LineWidth',2)
plot(vec_cores,K_gamma(3,:),'k-.','LineWidth',2)
plot(vec_cores,K_gamma(4,:),'k:','LineWidth',2)
plot(vec_cores,K_gamma(5,:),'LineWidth',2,'Color',grey)

xlabel('$N_{cores}$','Interpreter','latex')
ylabel('$f_e$','Interpreter','latex')

title('$d_m/g=0.9$','Interpreter','latex')

% size of the figure

% legend('c_{gap}=0.5% of D','c_{gap}=0.75% of D','c_{gap}=1% of D','c_{gap}=1.5% of D','c_{gap}=2% of D')




%% Plot K_e:

%setting the Matlab figure
f=figure('visible','on');
clf(f);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
set(gcf,'Units','centimeters','Position',[5 5 plotwidth plotheight])

load('results_cores_a')


grey=[0.4,0.4,0.4];

ax=axes('position',sub_pos{1,2},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
hold on
plot(vec_cores,Ke(1,:),'k-','LineWidth',2)
plot(vec_cores,Ke(2,:),'k--','LineWidth',2)
plot(vec_cores,Ke(3,:),'k-.','LineWidth',2)
plot(vec_cores,Ke(4,:),'k:','LineWidth',2)
plot(vec_cores,Ke(5,:),'LineWidth',2,'Color',grey)

xlabel('$N_{cores}$','Interpreter','latex')
ylabel('$K_E$','Interpreter','latex')

title('$d_m/g=0.6$','Interpreter','latex')

load('results_cores_b')


grey=[0.4,0.4,0.4];

ax=axes('position',sub_pos{2,2},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
hold on
plot(vec_cores,Ke(1,:),'k-','LineWidth',2)
plot(vec_cores,Ke(2,:),'k--','LineWidth',2)
plot(vec_cores,Ke(3,:),'k-.','LineWidth',2)
plot(vec_cores,Ke(4,:),'k:','LineWidth',2)
plot(vec_cores,Ke(5,:),'LineWidth',2,'Color',grey)

xlabel('$N_{cores}$','Interpreter','latex')
ylabel('$K_E$','Interpreter','latex')

title('$d_m/g=0.6$','Interpreter','latex')

load('results_cores_c')


grey=[0.4,0.4,0.4];

ax=axes('position',sub_pos{1,1},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
hold on
plot(vec_cores,Ke(1,:),'k-','LineWidth',2)
plot(vec_cores,Ke(2,:),'k--','LineWidth',2)
plot(vec_cores,Ke(3,:),'k-.','LineWidth',2)
plot(vec_cores,Ke(4,:),'k:','LineWidth',2)
plot(vec_cores,Ke(5,:),'LineWidth',2,'Color',grey)

xlabel('$N_{cores}$','Interpreter','latex')
ylabel('$K_E$','Interpreter','latex')

title('$d_m/g=0.6$','Interpreter','latex')

load('results_cores_d')


grey=[0.4,0.4,0.4];

ax=axes('position',sub_pos{2,1},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
hold on
plot(vec_cores,Ke(1,:),'k-','LineWidth',2)
plot(vec_cores,Ke(2,:),'k--','LineWidth',2)
plot(vec_cores,Ke(3,:),'k-.','LineWidth',2)
plot(vec_cores,Ke(4,:),'k:','LineWidth',2)
plot(vec_cores,Ke(5,:),'LineWidth',2,'Color',grey)

xlabel('$N_{cores}$','Interpreter','latex')
ylabel('$K_E$','Interpreter','latex')

title('$d_m/g=0.6$','Interpreter','latex')

% size of the figure

% legend('c_{gap}=0.5% of D','c_{gap}=0.75% of D','c_{gap}=1% of D','c_{gap}=1.5% of D','c_{gap}=2% of D')


