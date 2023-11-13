%% With the results from the study of F_opt and E_opt:

clc
clear all

set(0,'DefaultTextFontsize',8,'DefaultTextFontname','Times New Roman','DefaultAxesFontname','Times New Roman')

plotheight=10;
plotwidth=18;

subplotsx=2;
subplotsy=2;   
leftedge=1;
rightedge=0.4;   
topedge=0.4;
bottomedge=1;
spacex=1;
spacey=1.3;
fontsize=8;  

sub_pos=subplot_pos(plotwidth,plotheight,leftedge,rightedge,bottomedge,topedge,subplotsx,subplotsy,spacex,spacey);


%% Plot F_opt:

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
plot(vec_cores,K_opt(1,:),'k-','LineWidth',2)
plot(vec_cores,K_opt(2,:),'k--','LineWidth',2)
plot(vec_cores,K_opt(3,:),'k-.','LineWidth',2)
plot(vec_cores,K_opt(4,:),'k:','LineWidth',2)
plot(vec_cores,K_opt(5,:),'LineWidth',2,'Color',grey)

xlabel('$N_{cores}$','Interpreter','latex')
ylabel('$f_{opt}$','Interpreter','latex')

title('$d_m/g=0.6$','Interpreter','latex')

ylim([0 12])


load('results_cores_b')




ax=axes('position',sub_pos{2,2},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
hold on
plot(vec_cores,K_opt(1,:),'k-','LineWidth',2)
plot(vec_cores,K_opt(2,:),'k--','LineWidth',2)
plot(vec_cores,K_opt(3,:),'k-.','LineWidth',2)
plot(vec_cores,K_opt(4,:),'k:','LineWidth',2)
plot(vec_cores,K_opt(5,:),'LineWidth',2,'Color',grey)

xlabel('$N_{cores}$','Interpreter','latex')
ylabel('$f_{opt}$','Interpreter','latex')

title('$d_m/g=0.7$','Interpreter','latex')


load('results_cores_c')




ax=axes('position',sub_pos{1,1},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
hold on
plot(vec_cores,K_opt(1,:),'k-','LineWidth',2)
plot(vec_cores,K_opt(2,:),'k--','LineWidth',2)
plot(vec_cores,K_opt(3,:),'k-.','LineWidth',2)
plot(vec_cores,K_opt(4,:),'k:','LineWidth',2)
plot(vec_cores,K_opt(5,:),'LineWidth',2,'Color',grey)

xlabel('$N_{cores}$','Interpreter','latex')
ylabel('$f_{opt}$','Interpreter','latex')

title('$d_m/g=0.8$','Interpreter','latex')


load('results_cores_d')




figure
hold on
plot(vec_cores,K_opt(1,:),'k-','LineWidth',2)


xlabel('$N_{cores}$','Interpreter','latex')
ylabel('$f_{opt}$','Interpreter','latex')

title('$d_m/g=0.9$','Interpreter','latex')


% size of the figure

% legend('c_{gap}=0.5% of D','c_{gap}=0.75% of D','c_{gap}=1% of D','c_{gap}=1.5% of D','c_{gap}=2% of D')




%% Plot K_B:

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
plot(vec_cores,Kb(1,:),'k-','LineWidth',2)
plot(vec_cores,Kb(2,:),'k--','LineWidth',2)
plot(vec_cores,Kb(3,:),'k-.','LineWidth',2)
plot(vec_cores,Kb(4,:),'k:','LineWidth',2)
plot(vec_cores,Kb(5,:),'LineWidth',2,'Color',grey)

xlabel('$N_{cores}$','Interpreter','latex')
ylabel('$K_B$','Interpreter','latex')

title('$d_m/g=0.6$','Interpreter','latex')


load('results_cores_b')


grey=[0.4,0.4,0.4];

ax=axes('position',sub_pos{2,2},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
hold on
plot(vec_cores,Kb(1,:),'k-','LineWidth',2)
plot(vec_cores,Kb(2,:),'k--','LineWidth',2)
plot(vec_cores,Kb(3,:),'k-.','LineWidth',2)
plot(vec_cores,Kb(4,:),'k:','LineWidth',2)
plot(vec_cores,Kb(5,:),'LineWidth',2,'Color',grey)

xlabel('$N_{cores}$','Interpreter','latex')
ylabel('$K_B$','Interpreter','latex')

title('$d_m/g=0.7$','Interpreter','latex')

load('results_cores_c')


grey=[0.4,0.4,0.4];

ax=axes('position',sub_pos{1,1},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
hold on
plot(vec_cores,Kb(1,:),'k-','LineWidth',2)
plot(vec_cores,Kb(2,:),'k--','LineWidth',2)
plot(vec_cores,Kb(3,:),'k-.','LineWidth',2)
plot(vec_cores,Kb(4,:),'k:','LineWidth',2)
plot(vec_cores,Kb(5,:),'LineWidth',2,'Color',grey)

xlabel('$N_{cores}$','Interpreter','latex')
ylabel('$K_B$','Interpreter','latex')

title('$d_m/g=0.8$','Interpreter','latex')

load('results_cores_d')


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


