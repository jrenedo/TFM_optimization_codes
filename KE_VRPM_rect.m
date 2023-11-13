function [ K_E ] = KE_VRPM_rect( s,g,t,d,dm,N,res_B )

%   KE_VRPM: the function estimates the flux factor K_E to calculate the
%   torque in VRPM machines. It is based in conformal mapping and uses the
%   SC-toolbox by Driscoll and a logarithmic conformal transformation.
%   Author: Jaime Renedo Anglada, University of Southampton
%   s: slot width
%   g: airgap
%   t: tooth width
%   d: slot depth
%   R_coreback: radius of the coreback of the machine
%   dm: magnet width

%% Add the path of the library and define standar variables:


lambda=s+t;



%% Step 1
% Doing the map of the w-plane and calculating the constants.



g_new=g;


V=g_new;

% Geometric operations for the map:
% Generate a polygon with the geometry of the problem:
v=[d*i t/2+d*i t/2 t/2+s t/2+s+d*i t+s+d*i t+s+(d+g_new)*i (d+g_new)*i];
v1=v(1);
v2=v(2);
v3=v(3);
v4=v(4);
v5=v(5);
v6=v(6);
v7=v(7);
v8=v(8);

p=polygon(v);

% Indicates the right angles in the Canonical Domain:
alpha=[0.5 1 1 1 1 0.5 0.5 0.5];

% Remember that acording to this criteria X is the plane with the toothed
% member and w is the plane with the canonical rectangle.

% Define the Canonical Domain:

f=crrectmap(p,alpha);

% Vertices of the canonical rectangle:
vc1=evalinv(f,v1);
vc6=evalinv(f,v6);
vc7=evalinv(f,v7);
vc8=evalinv(f,v8);

h_canonical=abs(vc1-vc8); % height of the canonical rectangle.

B_w=V/h_canonical;

%% Iteration to obtain K_Bi:

vec_KB=zeros(1,N);
l_ag_line=length(0:(lambda/res_B):lambda);

for h=1:N
    Br=zeros(l_ag_line,1);
    
    % Geometric description of gamma_h:
    delta_h_real=h*dm/N;
%     delta_h=R_mg*log(R_coreback/(R_coreback-delta_h_real));
%     R_h=R_coreback-delta_h;
    l_h=delta_h_real;
    %l_h=R_mg*log(R_coreback/R_h);
    
    % Generate the line gamma_h (in X):
    ag_line=(d+g_new-l_h)*i+(0:(lambda/res_B):lambda);
    x=real(ag_line);
    
    for count=2:(l_ag_line-1)
        dif=evaldiff(f,evalinv(f,ag_line(count)));
        B_temp=B_w/conj(dif);
        
        Br(count)=real(B_temp);
        
    end
    Br(1)=Br(2);
    
    count=l_ag_line;
    Br(count)=Br(count-1);

    

    
    % Using trapz fuction:
    alpha=trapz(x,Br)/lambda;
    Bm=alpha;
    
    
    F1=Br.*transpose(cos(2*1*pi*x./lambda));
    F3=Br.*transpose(cos(2*3*pi*x./lambda));
    F5=Br.*transpose(cos(2*5*pi*x./lambda));
    
    gamma_1=trapz(x,F1)*2/(lambda*Bm);
    gamma_3=trapz(x,F3)*2/(lambda*Bm);
    gamma_5=trapz(x,F5)*2/(lambda*Bm);
    
    
    vec_KB(h)=gamma_1;
    vec_KB(h);
    
    
end


K_E=mean(vec_KB);


end

