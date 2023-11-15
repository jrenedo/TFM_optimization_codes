function [c,ceq] = constraints_core_ref2(x)

hs=x(3); % [mm]
ws=x(4); % [mm]

dm=x(6); % [mm]

g=x(5); % [mm]
mu_0=4*pi*10^-7;
J=x(7); % [A/mm^2]

Bcore_max=1.4;



wm=x(2);
wl=x(1);
hc=hs-1.1*wl;
wc=ws-2*wl;


F=J*hc*wc/2; % [A]
Bs=mu_0*F/g*1000; % [T]
Bcore=Bs*wm/wl;
Bcore_mag=(1.12*dm/g)*(wm/wl);

c = [((Bcore_mag+Bcore)-Bcore_max); wm-1.5*wl;2*(wm-wl)-0.8*(wc); wl-wm];     % Compute nonlinear inequalities at x.
ceq=[];

end