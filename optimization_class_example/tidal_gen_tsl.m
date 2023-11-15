%% Optimisation of the shape of a C-core:
% You have a Bs_max which comes from the level of flux from the PM, see my
% own notes for more details. The constraints of the C-core are
% geometrical and the air-gap is given. Specifications of the Casimere
% tidal generator from TSL

clc
clear all

Dyoke=422; %[mm]

N_phase=3;
N_rows=1;

Length=150; % [mm]
hs=45; % [mm]
ws=Length/N_phase/N_rows*0.8; % [mm]


dm=3; % [mm]
c_gap=2; % [mm] Grauers
g=c_gap+dm; % [mm]

D_ccores=Dyoke-2*dm-2*c_gap; %[mm]

J=5; % [A/mm^2]

% 150 rpm more less 
n_rpm=150; % [rpm]
n_rad=n_rpm/60*2*pi; % [rad/s]

% rho_steel=7650; % kg/m^3

rho_steel=7.650/1000; % kg/cm^3
rho_pm=7.3/1000; % kg/cm^3
rho_cu=8.96/1000; % kg/cm^3

cost_pm=45; % [euro/kg]
cost_steel=3.5; % [euro/kg]
cost_cu=12.4; % [euro/kg]

slot_fill=0.35;



%% Calculate the torque for a given t, s, lambda...


mu_0=4*pi*10^-7;

N=5;
res_B=500;

Mag=1.12/(4*pi*10^-7);

t_lambda=0.3;

%% Study 1
% variations of the number of C-cores


% 
count=1
vec_cores=[65];
for N_cores=vec_cores
    N_cores
    
    
    matrix=torque_optimised( hs,ws,dm,c_gap,D_ccores,N_cores,t_lambda,J,Mag, N, res_B,N_phase,N_rows);
    KB_s1_rect(count)=matrix(1,1);

    Torque_s1(count)=matrix(1,2); % [kNm]
%     Torque_s1(count)
%     Torque_s1(count)*n_rad
    
    vol_pm(count)=matrix(4,1)/1000; % [cm^3]
    vol_cores(count)=matrix(4,2)/1000; % [cm^3]
    vol_cu(count)=matrix(5,1)/1000; % [cm^3]
    
    mass_pm(count)=rho_pm*vol_pm(count); % [kg]
    mass_cores(count)=rho_steel*vol_cores(count); % [kg]
    mass_cu(count)=slot_fill*rho_cu*vol_cu(count); % [kg]
    
    tot_cost_pm(count)=cost_pm*mass_pm(count); % [euros]
    tot_cost_cores(count)=cost_steel*mass_cores(count); % [euros]
    tot_cost_cu(count)=cost_cu*mass_cu(count); % [euros]
    
    count=count+1;
end
figure
plot(vec_cores,Torque_s1,'o')
xlabel('C-cores')
ylabel('Torque [kNm]')

figure
plot(vec_cores,Torque_s1*n_rad,'o')
xlabel('C-cores')
ylabel('Power @ 150rpm [kW]')

Torque_s1*n_rad
N_cores
Torque_s1
n_rpm
J
N_phase
N_rows
K_B=matrix(1,1)

mass_pm
mass_cu
mass_cores

tot_cost_pm
tot_cost_cu
tot_cost_cores

Din=D_ccores-2*hs
D_out=430
c_gap
dm
N_phase*N_rows*ws
hs
ws
R_g=0.5*(D_ccores+g);
t=0.3*2*R_g*pi/N_cores

Vol_tot=pi*N_phase*N_rows*ws/4*(D_out^2-Din^2)/1000; % [cm^3]

Torque_s1/(mass_pm+mass_cu+mass_cores)*1000 % [Nm/kg]
Torque_s1/(mass_pm)*1000 % [Nm/kg]
Torque_s1/(tot_cost_pm+tot_cost_cu+tot_cost_cores)*1000 % [Nm/€]
Torque_s1/Vol_tot*(100^3) % [kNm/m^3]

A_eff=2*N_phase*N_rows*2*pi*R_g*matrix(3,2)/100; %[cm^2]

sigma=Torque_s1/A_eff*(100^2)/R_g*1000 % [kN/m^2]

% figure
% plot(vec_cores,tot_cost_pm+tot_cost_cores+tot_cost_cu,'o')
% xlabel('C-cores')
% ylabel('cost - euros')
% 
% 
% figure
% hold on
% plot(vec_cores,(mass_pm+mass_cores+mass_cu))
% plot(vec_cores,mass_pm)
% plot(vec_cores,mass_cores)
% plot(vec_cores,mass_cu)
% xlabel('C-cores')
% ylabel('mass [kg]')
% 
% legend('total','PM','steel','copper')


%% Study 2
% change magnet thickness, (more than 10mm the improvement of torque is
% small)

% clear Torque_s1
% N_cores=65
% % 
% count=1
% vec_dm=[0.1 0.5 1 2 3 4];
% for temp_dm=vec_dm
%     temp_dm
    
%     matrix=torque_optimised( hs,ws,temp_dm,c_gap,D_ccores,N_cores,t_lambda,J,Mag, N, res_B,N_phase,N_rows);
%     KB_s1_rect(count)=matrix(1,1);

%     Torque_s1(count)=matrix(1,2); % [kNm]
%     KB_s1_rect(count)
%     Torque_s1(count)
%     vol_pm(count)=matrix(4,1)/1000; % [cm^3]
%     vol_cores(count)=matrix(4,2)/1000; % [cm^3]
%     vol_cu(count)=matrix(5,1)/1000; % [cm^3]
    
%     mass_pm(count)=rho_pm*vol_pm(count); % [kg]
%     mass_cores(count)=rho_steel*vol_cores(count); % [kg]
%     mass_cu(count)=slot_fill*rho_cu*vol_cu(count); % [kg]
    
%     tot_cost_pm(count)=cost_pm*mass_pm(count); % [euros]
%     tot_cost_cores(count)=cost_steel*mass_cores(count); % [euros]
%     tot_cost_cu(count)=cost_cu*mass_cu(count); % [euros]
    
%     count=count+1;
% end
% % figure
% % plot(vec_dm,Torque_s1)
% % xlabel('Magnet thickness [mm]')
% % ylabel('Torque [kNm]')

% figure
% plot(vec_dm,Torque_s1*n_rad)
% xlabel('Magnet thickness [mm]')
% ylabel('Power @ 150rpm [kW]')

% figure
% hold on
% plot(vec_dm,(tot_cost_pm+tot_cost_cores+tot_cost_cu)/1000)
% plot(vec_dm,tot_cost_pm/1000,'g')
% plot(vec_dm,tot_cost_cores/1000,'c')
% plot(vec_dm,tot_cost_cu/1000,'r')
% xlabel('Magnet thickness [mm]')
% ylabel('cost - thousands euros')

% legend('total','PM','steel','copper')

% figure
% hold on
% plot(vec_dm,(mass_pm+mass_cores+mass_cu))
% plot(vec_dm,mass_pm,'g')
% plot(vec_dm,mass_cores,'c')
% plot(vec_dm,mass_cu,'r')
% xlabel('Magnet thickness [mm]')
% ylabel('mass [kg]')

% legend('total','PM','steel','copper')

% %% Study 3
% % change t_lambda

% % clear Torque_s1
% N_cores=65

% count=1
% vec_t_lambda=[0.28 0.3 0.32];
% for temp_t_lambda=vec_t_lambda
%     temp_t_lambda
%     matrix=torque_optimised( hs,ws,dm,c_gap,D_ccores,N_cores,temp_t_lambda,J,Mag, N, res_B,N_phase,N_rows)
%     KB_s1_rect(count)=matrix(1,1);

%     Torque_s1(count)=matrix(1,2); % [kNm]
%     KB_s1_rect(count)
%     Torque_s1(count)
    
%     count=count+1;
% end

% % figure
% % plot(vec_t_lambda,Torque_s1)
% % xlabel('t/lambda')
% % ylabel('Torque [kNm]')

% figure
% plot(vec_t_lambda,Torque_s1*n_rad)
% xlabel('t/lambda')
% ylabel('Power @ 150rpm [kW]')

