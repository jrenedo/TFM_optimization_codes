function [ output_args ] = torque_optimised( hs,ws,dm,c_gap,D_ccores,N_cores,t_lambda,J,Mag, N, res_B,N_phase,N_rows )
% This function optimises the 
%   Detailed explanation goes here


%% Calculate the torque for a given t, s, lambda...


g=c_gap+dm; % [mm]
R_g=D_ccores/2+g/2;
mu_0=4*pi*10^-7;




Fm=Mag*dm*10^-3;



%% Do the optimisation of the C-core leg width, wm

%syms x 
%fun = symfun(-(ws-2*x)*(hs-1.1*x)*x, x)

fun = @(x)(-(ws-2*x(1))*(hs-1.1*x(1))*x(2));

x0=[7 7 hs ws g dm J];

A=[];
b=[];
Aeq=[];
beq=[];
lb=[0 0 hs ws g dm J];
ub=[min([hs ws/2]) min([hs ws/2]) hs ws g dm J];


% you need to update in the constraints the geometrical values:

x_opt = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,@constraints_core_ref2);

wm=x_opt(2)
wl=x_opt(1)
hc=hs-1.1*wl
wc=ws-2*wl



F=J*hc*wc/2; % [A]
Bs=mu_0*F/g*1000 % [T]

%% Calculate the torque:
lambda=2*pi*R_g/N_cores;
depth=hc+1.1*wl;
t=t_lambda*lambda;
s=lambda-t;

KB_rect=KB_VRPM_rect( s,g,t,depth,dm,N,res_B ); 
Torque_1=N_rows*N_phase*2*N_cores*KB_rect*Bs*Fm*(wm*10^-3)*2*(R_g*10^-3)/1000; % [kNm]

%% Calculate volumes


vol_pm=2*pi*(R_g+g/2-dm/2)*wm*dm*2*N_rows*N_phase; % [mm^3]
vol_cores=N_rows*N_phase*N_cores*(hs*ws-hc*wc)*t; % [mm^3]
vol_yoke=N_rows*N_phase*ws*pi*((R_g+g/2+wl)^2-(R_g+g/2)^2); % [mm^3]
vol_cu=2*pi*(R_g-g/2-hc/2)*wc*hc*N_rows*N_phase; % [mm^3]

vol_tot=N_rows*N_phase*pi*(R_g+g/2+wl)^2*ws; % [mm^3]

output_args=[KB_rect Torque_1; hc wc; wl wm; vol_pm vol_cores+vol_yoke; vol_cu vol_tot];


end

