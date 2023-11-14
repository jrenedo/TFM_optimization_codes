classdef airgapVRPM_no_screen_curvature < handle
% airgap_screen: represents the airgap of the VRPM machine with
% superconducting screens to clock the magnetic field. Considering the
% effect of curvature.
% Jaime Renedo Anglada, University of Southampton

    properties
        s_self
        t_self
        g_self
        d_self
        gnew_self
        Rg_self
        n_cores
        
        polygon_self
        map_self
        h_canonical
        
        coeffs_pf
    end
    methods 
        %% Definition of the object:
        
        function obj = airgapVRPM_no_screen_curvature (input_s, input_t, input_g, input_d, input_Rg,n_cores)
            obj.s_self = input_s;
            obj.t_self = input_t;
            obj.g_self = input_g;
            obj.Rg_self = input_Rg;
            obj.d_self = input_d;
            obj.n_cores=n_cores;
        end
        
        %% Building the conformal map from the SC toolbox
        
        function obj = build_map (obj)
            s=obj.s_self;
            t=obj.t_self;
            g=obj.g_self;
            Rg=obj.Rg_self;
            d=obj.d_self;
            
            gnew=Rg*log((Rg+g/2)/(Rg-g/2));
            obj.gnew_self=gnew;
            
            % Polygon for SC toolbox
            path(path,'C:\local\Documents Jaime\PhD\MATLAB\Driscoll SA toolbox\sc')
            % Geometric operations for the map:
            % Generate a polygon with the geometry of the problem:
            v=[d*i t/2+d*i t/2 t/2+s t/2+s+d*i t+s+d*i t+s+(d+gnew)*i (d+gnew)*i];
            v1=v(1);
            v2=v(2);
            v3=v(3);
            v4=v(4);
            v5=v(5);
            v6=v(6);
            v7=v(7);
            v8=v(8);
            
            p=polygon(v);
            obj.polygon_self=p;
            
            % Indicates the right angles in the Canonical Domain:
            alpha=[0.5 1 1 1 1 0.5 0.5 0.5];
            
            % Remember that acording to this criteria X is the plane with the toothed
            % member and w is the plane with the canonical rectangle.
            
            % Define the Canonical Domain:
            
            f=crrectmap(p,alpha);
            obj.map_self=f;
            
            % Vertices of the canonical rectangle:
            vc1=evalinv(f,v1);
            vc6=evalinv(f,v6);
            vc7=evalinv(f,v7);
            vc8=evalinv(f,v8);
            
            obj.h_canonical=abs(vc1-vc8); % height of the canonical rectangle.
            
        end
        %% Magnetic field waveform
        % delta is the distance to the coreback, V is the MMF and n_points
        % is the number of points evaluated. In case we want to plot it.
        
        function result = B_func (obj,delta,V,n_points)
            
            s=obj.s_self;
            t=obj.t_self;
            g=obj.g_self;
            d=obj.d_self;
            gnew=obj.gnew_self;
            Rg=obj.Rg_self;
            
            f=obj.map_self;
            
            result=zeros(3,2*n_points-1);
            
            
            B_w=V/obj.h_canonical;
            
            l_ag_line=n_points;
            Br=zeros(l_ag_line,1);
            Bt=zeros(l_ag_line,1);
            
            Ri=Rg+g/2-delta;
            delta_SC=Rg*log((Rg+g/2)/(Rg+g/2-delta)); %distance in the w plane
            
            ag_line=(gnew-delta_SC+d)*i+(s+t)/2*(0:n_points)/n_points;
            x=real(ag_line);
            
            for count=2:(l_ag_line-1)
                dif=evaldiff(obj.map_self,evalinv(f,ag_line(count)));
                B_temp=B_w/conj(dif)*Rg/Ri;
                
                Br(count)=real(B_temp);
                Bt(count)=imag(B_temp);
                
            end
            Br(1)=Br(2);
            Bt(1)=0;
            
            count=l_ag_line;
            Br(count)=Br(count-1);
            Bt(count)=0;
            
            for count=1:n_points
                result(1,count)=Br(count);
                result(2,count)=Bt(count);
                result(3,count)=x(count);
                
                if count<n_points
                    result(1,2*n_points-count)=Br(count);
                    result(2,2*n_points-count)=-Bt(count);
                    result(3,2*n_points-count)=s+t-x(count);
                end
                
            end
            
        end
        
        
        %% Evaluation of lambda(theta,r) at one point
        % delta is the distance to the coreback, V is the MMF and n_points
        % is the number of points evaluated. In case we want to plot it.
        
        function result = eval_lambda (obj,theta_in,r_in)
            
            s=obj.s_self;
            t=obj.t_self;
            g=obj.g_self;
            d=obj.d_self;
            gnew=obj.gnew_self;
            Rg=obj.Rg_self;
            
            f=obj.map_self;
            
            V=g;
            
            B_w=V/obj.h_canonical;
            
            
            
            Ri=r_in;
            delta=Rg+g/2-Ri;
            delta_SC=Rg*log((Rg+g/2)/(Rg+g/2-delta)); %distance in the w plane
            
            
            
            point=(gnew-delta_SC+d)*i+(s+t)/(2*pi)*theta_in;
            
            dif=evaldiff(obj.map_self,evalinv(f,point));
            B_temp=B_w/conj(dif)*Rg/Ri;
            
            
            
            result=B_temp;
            
        end
        
                        
        %% Evaluation of lambda(theta,r) at one point
        % delta is the distance to the coreback, V is the MMF and n_points
        % is the number of points evaluated. In case we want to plot it.
        
        function result = gen_vec_points(obj,N_points,dm)
            g=obj.g_self;
            Rg=obj.Rg_self;
            
            for k=1:N_points
                result(k,:)=[(0.01+0.98*rand(1,1))*2*pi Rg+g/2-(0.01+0.98*rand(1,1))*dm ];
            end
            
            
        end        
        
        
        %% Evaluation of lambda(theta,r) at one point
        % delta is the distance to the coreback, V is the MMF and n_points
        % is the number of points evaluated. In case we want to plot it.
        
        function result = eval_coeffs_lambdar (obj,vec_points,N_harmonics,N_poly)
            
            s=obj.s_self;
            t=obj.t_self;
            g=obj.g_self;
            d=obj.d_self;
            gnew=obj.gnew_self;
            Rg=obj.Rg_self;
            
            f=obj.map_self;
            
            result=zeros(N_harmonics+1,N_poly+1);
            lambda_sol=zeros(1,length(vec_points(:,1)));
            delta_temp=zeros(1,length(vec_points(:,1)));
            
            for count1=1:length(vec_points(:,1))
%                 Rg+g/2-vec_points(count1,2)
                lambda_sol(count1)=real(obj.eval_lambda(vec_points(count1,1),vec_points(count1,2)))*vec_points(count1,2)/Rg;
                delta_temp(count1)=Rg+g/2-vec_points(count1,2);
                theta_temp(count1)=vec_points(count1,1);
            end
            
            Delta=zeros(length(vec_points(:,1)),1+N_harmonics*(N_poly+1));
            
            Delta(:,1)=Delta(:,1)+1;
            
            for count1=1:length(vec_points(:,1))
                control_count=0;
                control_count2=1;
                for count2=2:(N_harmonics*(N_poly+1)+1)
                    Delta(count1,count2)=delta_temp(count1)^control_count*cos(control_count2*theta_temp(count1));
                    if control_count==N_poly
                        control_count=0;
                        control_count2=control_count2+1;
                    else
                        control_count=control_count+1;
                    end
                end
            end
            
%             result=Delta;
%             Delta
            
            Gamma=(inv(transpose(Delta)*Delta)*transpose(Delta))*transpose(lambda_sol);
            
            result(1,1)=Gamma(1);
            count_aux=2;
            for count1=2:(N_harmonics+1)
                
                for count2=1:(N_poly+1)
                    result(count1,count2)=Gamma(count_aux)/Gamma(1);
                    count_aux=count_aux+1;
                end
            end
            
            
            obj.coeffs_pf=result;
            
            
            
        end
        
        
         %% Evaluation of lambda(theta,r) at one point
        % delta is the distance to the coreback, V is the MMF and n_points
        % is the number of points evaluated. In case we want to plot it.
        
        function result = pflux_linkage (obj,dm)
            
            coeffs_func=obj.coeffs_pf;
            n_cores=obj.n_cores;
            
            resolution=25;
            
            N_terms=length(coeffs_func(1,:));
            N_harm=length(coeffs_func(:,1));
            
            coeffs_func(1,1)
            
            coeffs_mod=coeffs_func((2:N_harm),:);
            
            count2=1;
            for delta_c=0:(dm/resolution):dm
                for h=1:N_terms
                    temp_delta(h,1)=delta_c^(h-1);
                end
                
                gamma_matrix(count2,:)=coeffs_mod*temp_delta;
                
                count2=count2+1;
            end
%             gamma_matrix
            gamma_av1=mean(gamma_matrix(:,1))
            gamma_av3=mean(gamma_matrix(:,3))
            gamma_av5=mean(gamma_matrix(:,5))
            gamma_av7=mean(gamma_matrix(:,7))
            
            mu_0=4*pi*10^-7;
            Rg=obj.Rg_self
            L=0.036/2;
            Mag=1.05/(4*pi*10^-7);
            n=20;
            N_wind=230;
            g=obj.g_self
            
            K=4*mu_0*Rg*L/g/n_cores; % we need to check this
            
            beta=0:(2*pi/100):(2*pi);
            M_func=K*coeffs_func(1,1)*(gamma_av1*cos(beta)-gamma_av3*cos(3*beta)/3+gamma_av5*cos(5*beta)/5-gamma_av7*cos(7*beta)/7);
            result=M_func;
            
%             Phi_func=n*M_func*Mag*dm*10^-3*N_wind;
%             
%             plot(beta,Phi_func)
%             for beta=0:(2*pi/100):(2*pi)
%                 
%                 
%             end
            
        end
        
        
        
        %% Expression of the magnetic field as Fourier Series:
        
        function fourier_coeffs = fourier_series(obj,delta,V,n_points)
            temp_matrix=obj.B_func (delta,V,n_points);
            lambda=obj.s_self+obj.t_self; % because we defined x of length lambda
            
            x=temp_matrix(3,:);
%             x(length(x))
%             lambda
            Br=temp_matrix(1,:);
            Bt=temp_matrix(2,:);
            
            
            % Using trapz fuction:
            alpha=trapz(x,Br)/lambda;
            Bm=alpha;
            
            
            F1=Br.*(cos(2*1*pi*x./lambda));
            F2=Br.*(cos(2*2*pi*x./lambda));
            F3=Br.*(cos(2*3*pi*x./lambda));
            F4=Br.*(cos(2*4*pi*x./lambda));
            F5=Br.*(cos(2*5*pi*x./lambda));
            
            gamma_1=trapz(x,F1)*2/(lambda*Bm);
            gamma_2=trapz(x,F2)*2/(lambda*Bm);
            gamma_3=trapz(x,F3)*2/(lambda*Bm);
            gamma_4=trapz(x,F4)*2/(lambda*Bm);
            gamma_5=trapz(x,F5)*2/(lambda*Bm);
            
            fourier_coeffs=[Bm gamma_1 gamma_2 gamma_3 gamma_4 gamma_5];
        end
        
        
        %% KB calculation:
        % n_int: the number of points for the integration layers
        
        function Kb_out = Kb_calc_square (obj,dm)
            coeffs_func=obj.coeffs_pf;
            n_cores=obj.n_cores;
            
            resolution=50;
            
            N_terms=length(coeffs_func(1,:));
            N_harm=length(coeffs_func(:,1));
            
            coeffs_func(1,1);
            
            coeffs_mod=coeffs_func((2:N_harm),:);
            
            count2=1;
            for delta_c=0:(dm/resolution):dm
                for h=1:N_terms
                    temp_delta(h,1)=delta_c^(h-1);
                end
                
                gamma_matrix(count2,:)=coeffs_mod*temp_delta;
                
                count2=count2+1;
            end
%             gamma_matrix
            gamma_av1=mean(gamma_matrix(:,1));
            gamma_av3=mean(gamma_matrix(:,3));
            gamma_av5=mean(gamma_matrix(:,5));
%             gamma_av7=mean(gamma_matrix(:,7));
            gamma_av7=0;
            
%             for h=1:N_terms
%                 dm_int(h,1)=dm^(h)/h;
%             end
            
%             int_gamma=coeffs_mod*dm_int
            
            
%             Kb_out1=4/pi*coeffs_func(1,1)/dm*(int_gamma(1)-int_gamma(3)/3+int_gamma(5)/7-int_gamma(7)/7)
            
            
%             gamma_av1=trapz(0:(dm/resolution):dm,gamma_matrix(:,1))/dm;
%             gamma_av3=trapz(0:(dm/resolution):dm,gamma_matrix(:,3))/dm;
%             gamma_av5=trapz(0:(dm/resolution):dm,gamma_matrix(:,5))/dm;
%             gamma_av7=trapz(0:(dm/resolution):dm,gamma_matrix(:,7))/dm;
            
            Kb_out=4/pi*coeffs_func(1,1)*(gamma_av1-gamma_av3/3+gamma_av5/5-gamma_av7/7)
            
        end
        
        
              
                
        %% KB calculation:
        % n_int: the number of points for the integration layers
        
        function Kb_out = Kb_calc_sine (obj,dm)
            coeffs_func=obj.coeffs_pf;
            n_cores=obj.n_cores;
            
            resolution=50;
            
            N_terms=length(coeffs_func(1,:));
            N_harm=length(coeffs_func(:,1));
            
            coeffs_func(1,1);
            
            coeffs_mod=coeffs_func((2:N_harm),:);
            
            count2=1;
            for delta_c=0:(dm/resolution):dm
                for h=1:N_terms
                    temp_delta(h,1)=delta_c^(h-1);
                end
                
                gamma_matrix(count2,:)=coeffs_mod*temp_delta;
                
                count2=count2+1;
            end
%             gamma_matrix
            gamma_av1=mean(gamma_matrix(:,1));
            gamma_av3=mean(gamma_matrix(:,3));
            gamma_av5=mean(gamma_matrix(:,5));
%             gamma_av7=mean(gamma_matrix(:,7));
            gamma_av7=0;
            
%             for h=1:N_terms
%                 dm_int(h,1)=dm^(h)/h;
%             end
            
%             int_gamma=coeffs_mod*dm_int
            
            
%             Kb_out1=4/pi*coeffs_func(1,1)/dm*(int_gamma(1)-int_gamma(3)/3+int_gamma(5)/7-int_gamma(7)/7)
            
            
%             gamma_av1=trapz(0:(dm/resolution):dm,gamma_matrix(:,1))/dm;
%             gamma_av3=trapz(0:(dm/resolution):dm,gamma_matrix(:,3))/dm;
%             gamma_av5=trapz(0:(dm/resolution):dm,gamma_matrix(:,5))/dm;
%             gamma_av7=trapz(0:(dm/resolution):dm,gamma_matrix(:,7))/dm;
            
            Kb_out=coeffs_func(1,1)*gamma_av1*sqrt(2);
            
        end
        
        
              
        %% Plot the magnetic field:
        % n_int: the number of points for the integration layers
        
        function result = plot_field (obj,n_pointsu,n_pointsv,Bs)
            
            g=obj.g_self;
            s=obj.s_self;
            t=obj.t_self;
            d=obj.d_self;
            Rg=obj.Rg_self;
            R_g=Rg;
            lambda=s+t;
            f=obj.map_self;
            n_cores=obj.n_cores;
            V=g;
            
            B_w=V/obj.h_canonical;
            theta_lambda=lambda/Rg;
            theta_s=s/Rg;
            theta_t=t/Rg;
            
            
            Z_coord=zeros(2*n_pointsu,n_pointsv);
            B_com=zeros(2*n_pointsu,n_pointsv);
            Br=zeros(2*n_pointsu,n_pointsv);
            Bt=zeros(2*n_pointsu,n_pointsv);
            
            distance_g=(0:(2*n_pointsu))*g*2/n_pointsu;
            %distance_g_w=Rg*log((Rg+g/2)/(Rg+g/2-distance_g))
            theta_vec=(0:n_pointsv)*theta_lambda/n_pointsv;
            
            for k=1:(2*n_pointsu)
                for h=1:(n_pointsv)
                    r_temp=Rg+g/2-distance_g(k);
                    Z_coord(k,h)=r_temp*exp(i*theta_vec(h));
                    d_w=Rg*log((Rg+g/2)/r_temp);
%                     Z_coord(k,h)
                    w_temp=i*(d+g-d_w)+theta_vec(h)*Rg;
                    if r_temp>(Rg-g/2) 
                        if k==1
                            r_temp=Rg-distance_g(2)+g/2;
                            d_w=Rg*log((Rg+g/2)/r_temp);
                            w_temp=i*(d+g-d_w)+theta_vec(h)*Rg+lambda/10000;
                            dif=evaldiff(obj.map_self,evalinv(f,w_temp));
                            B_temp=B_w/conj(dif)*Rg/r_temp;
                            B_com(k,h)=B_temp;
                        elseif h==1
                            w_temp=i*(d+g-d_w)+theta_vec(2)*Rg;
                            B_temp=B_w/conj(dif)*Rg/r_temp;
                            B_com(k,h)=B_temp;
                        elseif h==n_pointsv
                            w_temp=i*(d+g-d_w)+theta_vec(n_pointsv-1)*Rg;
                            B_temp=B_w/conj(dif)*Rg/r_temp;
                            B_com(k,h)=B_temp;
                        
                        else
%                             w_temp
                            dif=evaldiff(obj.map_self,evalinv(f,w_temp));
                            B_temp=B_w/conj(dif)*Rg/r_temp;
                            B_com(k,h)=B_temp;
                        end
                        
                    elseif theta_vec(h)>(theta_t/2) && theta_vec(h)<(theta_t/2+theta_s)
                        dif=evaldiff(obj.map_self,evalinv(f,w_temp));
                        B_temp=B_w/conj(dif)*Rg/r_temp;
                        B_com(k,h)=B_temp;
%                         B_com(k,h)=0;
                        
                    else
                        B_com(k,h)=NaN;
                        
                    end
                    
                end
            end
            X=real(Z_coord);
            Y=imag(Z_coord);
            
           
            
            
            F_z=B_com*Bs;
            angle_step=2*pi/n_cores;
            
            grey = [0.4,0.4,0.4];
            

            
            figure
            hold on
            
            
        
        
                
            for count=1:n_cores
                temp_angle=count*angle_step;
                temp_number=cos(temp_angle)+i*sin(temp_angle);
                Z_temp=Z_coord*temp_number;
                z2=[R_g-g/2+t/2*i R_g-g/2-t/2*i R_g-g/2-d/2-t/2*i R_g-g/2-d/2+t/2*i]*temp_number;
                x_temp=real(z2);
                y_temp=imag(z2);
                X=real(Z_temp);
                Y=imag(Z_temp);
                contourf(X,Y,real(F_z),100,'LineStyle','none')
                fill(x_temp,y_temp,grey)
                
            end
            for count=1:n_cores
                temp_angle=count*angle_step;
                temp_number=cos(temp_angle)+i*sin(temp_angle);
                
                z2=[R_g-g/2+t/2*i R_g-g/2-t/2*i R_g-g/2-d/2-t/2*i R_g-g/2-d/2+t/2*i]*temp_number;
                x_temp=real(z2);
                y_temp=imag(z2);
                
                fill(x_temp,y_temp,grey)
            end
            t = (0:0.001:1)'*2*pi;
            x1 = (R_g+3*g/2)*sin(t);
            y1 = (R_g+3*g/2)*cos(t);
            
            y2 = (R_g+g/2)*cos(t);
            x2 = (R_g+g/2)*sin(t);
            
            x3=[x1; x2];
            y3=[y1; y2];
            grey = [0.4,0.4,0.4];
            fill(x3,y3,grey)
            colorbar
            h = colorbar;
            colorbar
            h = colorbar;
            
            ylabel2 = get(h,'YTickLabel');
            
            mm = repmat(' T',size(ylabel2,1),1);
            
            ylabel2 = [ylabel2 mm];
            
            set(h,'YTickLabel',ylabel2);
            h = colorbar;
            
            ylabel2 = get(h,'YTickLabel');
            
            mm = repmat(' T',size(ylabel2,1),1);
            
            ylabel2 = [ylabel2 mm];
            
            set(h,'YTickLabel',ylabel2);
            caxis([0 Bs])
            xlim([0 Rg+2*g])
            ylim([0 Rg+2*g])
            axis square
            ylabel('y (mm)','Interpreter','latex');
            xlabel('x (mm)','Interpreter','latex');
            title('Radial component, $B_r$ (T)','Interpreter','latex')
%             t = (0:0.001:1)'*2*pi;
%             x = (R_g+3*g/2)*sin(t);
%             y = (R_g+3*g/2)*cos(t);
%             fill(x,y,grey)
%             clear x y
%             x = (R_g+g/2)*sin(t);
%             y = (R_g+g/2)*cos(t);
%             fill(x,y,'w')
            clear z2
            figure
            hold on
            for count=1:n_cores
                temp_angle=count*angle_step;
                temp_number=cos(temp_angle)+i*sin(temp_angle);
                Z_temp=Z_coord*temp_number;
                z2=[R_g-g/2+t/2*i R_g-g/2-t/2*i R_g-g/2-d/2-t/2*i R_g-g/2-d/2+t/2*i]*temp_number;
                x_temp=real(z2);
                y_temp=imag(z2);
                X=real(Z_temp);
                Y=imag(Z_temp);
                contourf(X,Y,abs(imag(F_z)),100,'LineStyle','none')
                fill(x_temp,y_temp,grey)
            end
            for count=1:n_cores
                temp_angle=count*angle_step;
                temp_number=cos(temp_angle)+i*sin(temp_angle);
                
                z2=[R_g-g/2+t/2*i R_g-g/2-t/2*i R_g-g/2-d/2-t/2*i R_g-g/2-d/2+t/2*i]*temp_number;
                x_temp=real(z2);
                y_temp=imag(z2);
                
                fill(x_temp,y_temp,grey)
            end
            t = (0:0.001:1)'*2*pi;
            x1 = (R_g+3*g/2)*sin(t);
            y1 = (R_g+3*g/2)*cos(t);
            
            y2 = (R_g+g/2)*cos(t);
            x2 = (R_g+g/2)*sin(t);
            
            x3=[x1; x2];
            y3=[y1; y2];
            grey = [0.4,0.4,0.4];
            fill(x3,y3,grey)
            colorbar
            h = colorbar;
            
            ylabel2 = get(h,'YTickLabel');
            
            mm = repmat(' T',size(ylabel2,1),1);
            
            ylabel2 = [ylabel2 mm];
            
            set(h,'YTickLabel',ylabel2);
            h = colorbar;
            
            ylabel2 = get(h,'YTickLabel');
            
            mm = repmat(' T',size(ylabel2,1),1);
            
            ylabel2 = [ylabel2 mm];
            
            set(h,'YTickLabel',ylabel2);
            caxis([0 Bs])
            xlim([0 Rg+2*g])
            ylim([0 Rg+2*g])
            axis square
            ylabel('y (mm)','Interpreter','latex')
            xlabel('x (mm)','Interpreter','latex')
            title('Tangential component','Interpreter','latex')
            
            
%             figure
%             surf(X,Y,real(F_z))
%             
%             shading interp;
%             colormap(jet(128))
            
%             figure
%             %VectorField2d([real(F_z), imag(F_z)], X,Y);
%             quiver(X,Y,real(F_z), imag(F_z))
            
%            figure
%             %VectorField2d([real(F_z), imag(F_z)], X,Y);
%             contourf(X,Y,real(F_z),100,'LineStyle','none')
%             
%             figure
%             %VectorField2d([real(F_z), imag(F_z)], X,Y);
%             contourf(X,Y,imag(F_z),100,'LineStyle','none')
            
            
        end
        
              
        %% Plot the magnetic field contour:
        % n_int: the number of points for the integration layers
        
        function result = plot_field_contour (obj,n_pointsu,n_pointsv,Bs)
            
            g=obj.g_self;
            s=obj.s_self;
            t=obj.t_self;
            d=obj.d_self;
            Rg=obj.Rg_self;
            lambda=s+t;
            f=obj.map_self;
            n_cores=obj.n_cores;
            V=g;
            
            B_w=V/obj.h_canonical;
            theta_lambda=lambda/Rg;
            theta_s=s/Rg;
            theta_t=t/Rg;
            
            
            Z_coord=zeros(2*n_pointsu,n_pointsv);
            B_com=zeros(2*n_pointsu,n_pointsv);
            Br=zeros(2*n_pointsu,n_pointsv);
            Bt=zeros(2*n_pointsu,n_pointsv);
            
            distance_g=(0:(2*n_pointsu))*g*2/n_pointsu;
            %distance_g_w=Rg*log((Rg+g/2)/(Rg+g/2-distance_g))
            theta_vec=(0:n_pointsv)*theta_lambda/n_pointsv;
            
            for k=1:(2*n_pointsu)
                for h=1:(n_pointsv)
                    r_temp=Rg+g/2-distance_g(k);
                    Z_coord(k,h)=r_temp*exp(i*theta_vec(h));
                    d_w=Rg*log((Rg+g/2)/r_temp);
%                     Z_coord(k,h)
                    w_temp=i*(d+g-d_w)+theta_vec(h)*Rg;
                    if r_temp>(Rg-g/2) 
                        if k==1
                            r_temp=Rg-distance_g(2)+g/2;
                            d_w=Rg*log((Rg+g/2)/r_temp);
                            w_temp=i*(d+g-d_w)+theta_vec(h)*Rg+lambda/10000;
                            dif=evaldiff(obj.map_self,evalinv(f,w_temp));
                            B_temp=B_w/conj(dif)*Rg/r_temp;
                            B_com(k,h)=B_temp;
                        elseif h==1
                            w_temp=i*(d+g-d_w)+theta_vec(2)*Rg;
                            B_temp=B_w/conj(dif)*Rg/r_temp;
                            B_com(k,h)=B_temp;
                        elseif h==n_pointsv
                            w_temp=i*(d+g-d_w)+theta_vec(n_pointsv-1)*Rg;
                            B_temp=B_w/conj(dif)*Rg/r_temp;
                            B_com(k,h)=B_temp;
                        
                        else
%                             w_temp
                            dif=evaldiff(obj.map_self,evalinv(f,w_temp));
                            B_temp=B_w/conj(dif)*Rg/r_temp;
                            B_com(k,h)=B_temp;
                        end
                        
                    elseif theta_vec(h)>(theta_t/2) && theta_vec(h)<(theta_t/2+theta_s)
                        dif=evaldiff(obj.map_self,evalinv(f,w_temp));
                        B_temp=B_w/conj(dif)*Rg/r_temp;
                        B_com(k,h)=B_temp;
%                         B_com(k,h)=0;
                        
                    else
                        B_com(k,h)=NaN;
                        
                    end
                    
                end
            end
            X=real(Z_coord);
            Y=imag(Z_coord);
            
           
            
            
            F_z=B_com*Bs;
            angle_step=2*pi/n_cores;
            
            grey = [0.4,0.4,0.4];

            
            figure
            hold on
            
            
        
        
                
            for count=1:n_cores
                temp_angle=count*angle_step;
                temp_number=cos(temp_angle)+i*sin(temp_angle);
                Z_temp=Z_coord*temp_number;
                X=real(Z_temp);
                Y=imag(Z_temp);
                contour(X,Y,real(F_z),20)
            end
            xlim([-(Rg+g) Rg+g])
            ylim([-(Rg+g) Rg+g])
            axis square
            colorbar
            caxis([0 Bs])
            ylabel('y (mm)')
            xlabel('x (mm)')
            title('Radial component')
            
            figure
            hold on
            for count=1:n_cores
                temp_angle=count*angle_step;
                temp_number=cos(temp_angle)+i*sin(temp_angle);
                Z_temp=Z_coord*temp_number;
                X=real(Z_temp);
                Y=imag(Z_temp);
                contour(X,Y,abs(imag(F_z)),20)
            end
            xlim([-(Rg+g) Rg+g])
            ylim([-(Rg+g) Rg+g])
            axis square
            colorbar
            caxis([0 Bs])
            ylabel('y (mm)')
            xlabel('x (mm)')
            title('Tangential component')
            
            
%             figure
%             surf(X,Y,real(F_z))
%             
%             shading interp;
%             colormap(jet(128))
            
%             figure
%             %VectorField2d([real(F_z), imag(F_z)], X,Y);
%             quiver(X,Y,real(F_z), imag(F_z))
            
%            figure
%             %VectorField2d([real(F_z), imag(F_z)], X,Y);
%             contourf(X,Y,real(F_z),100,'LineStyle','none')
%             
%             figure
%             %VectorField2d([real(F_z), imag(F_z)], X,Y);
%             contourf(X,Y,imag(F_z),100,'LineStyle','none')
            
            
        end
        
    end
    
        
        
  
    
    
%     methods(Static)
%         function d = distanceBetweenShapes(shape1,shape2)
%             xDist = abs(shape1.centerX - shape2.centerX);
%             yDist = abs(shape1.centerY - shape2.centerY);
%             d = sqrt(xDist^2 + yDist^2);
%         end
%     end
end 