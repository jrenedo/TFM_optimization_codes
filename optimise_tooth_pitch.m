%   This code studies the optimal number of C-cores as a function of the normalised parameters.
%   It is based in conformal mapping and uses the
%   SC-toolbox by Driscoll and a logarithmic conformal transformation.
%   Author: Jaime Renedo Anglada, University of Southampton
%   s: slot width
%   g: airgap
%   t: tooth width
%   d: slot depth
%   R_coreback: radius of the coreback of the machine
%   dm: magnet width

% with: N=11; res_B=500; it takes around 14 hours

%% Add the path of the library and define standar variables:

clc
clear all

path(path,'/Users/jaime/Documents/work/sc-toolbox-3.1.3')

% Choose a function of the clearance gap f(D). Later a dm/g, with this we can optimise the number of C-cores
t_lambda=0.4;



N=11;
res_B=100;

vec_cores=2:2:80;
D=430*10^-3;
c_gap=1/100*D;
vec_t_lambda=[0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6];
dm_g=0.8;

% The time approx is: length(vec_cores)*length(vec_gaps) minutes.
% Total takes around 5 hours.

tot_count=1;

for t_lambda = vec_t_lambda
    count=1;
    g=c_gap/(1-dm_g);
    dm=dm_g*g;

    for N_cores=vec_cores
        count

        lambda=pi*D/N_cores
        t=t_lambda*lambda;
        s=lambda-t;
        d=3*lambda;
        temp=K_tot_VRPM_rect( s,g,t,d,dm,N,res_B );
        Kb(tot_count,count)=temp(1);
        K_opt(tot_count,count)=Kb(tot_count,count)*N_cores*dm/g;
        Ke(tot_count,count)=temp(2);
        Ke_opt(tot_count,count)=Ke(tot_count,count)*dm/g;
        K_gamma(tot_count,count)=temp(3);
        count=count+1;

    end
    tot_count=tot_count+1;
end

save('results_cores_tlambda')


grey=[0.4,0.4,0.4];

%% results plot

figure
hold on
plot(vec_cores,K_opt(1,:),'LineWidth',2)
plot(vec_cores,K_opt(2,:),'LineWidth',2)
plot(vec_cores,K_opt(3,:),'LineWidth',2)
plot(vec_cores,K_opt(4,:),'LineWidth',2)
plot(vec_cores,K_opt(5,:),'LineWidth',2)
plot(vec_cores,K_opt(6,:),'LineWidth',2)
plot(vec_cores,K_opt(8,:),'LineWidth',2)
plot(vec_cores,K_opt(9,:),'LineWidth',2)
plot(vec_cores,K_opt(10,:),'LineWidth',2)
plot(vec_cores,K_opt(11,:),'LineWidth',2)

xlabel('N_{cores}')
ylabel('F_{opt}')

title('d_m/g=0.9')

