classdef airgapVRPM_no_screen < handle
% airgap_screen: represents the airgap of the VRPM machine with
% superconducting screens to clock the magnetic field. Not considering the
% effect of curvature.
% Jaime Renedo Anglada, University of Southampton

    properties
        s_self
        t_self
        g_self
        d_self
        
        polygon_self
        map_self
        h_canonical
    end
    methods 
        %% Definition of the object:
        
        function obj = airgapVRPM_no_screen (input_s, input_t, input_g, input_d)
            obj.s_self = input_s;
            obj.t_self = input_t;
            obj.g_self = input_g;
            obj.d_self = input_d;
        end
        
        %% Building the conformal map from the SC toolbox
        
        function obj = build_map (obj)
            s=obj.s_self;
            t=obj.t_self;
            g=obj.g_self;
            d=obj.d_self;
            
            % Polygon for SC toolbox
            path(path,'C:\local\Documents Jaime\PhD\MATLAB\Driscoll SA toolbox\sc')
            % Geometric operations for the map:
            % Generate a polygon with the geometry of the problem:
            v=[d*i t/2+d*i t/2 t/2+s t/2+s+d*i t+s+d*i t+s+(d+g)*i (d+g)*i];
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
            
            f=obj.map_self;
            
            result=zeros(3,2*n_points-1);
            
            
            B_w=V/obj.h_canonical;
            
            l_ag_line=n_points;
            Br=zeros(l_ag_line,1);
            Bt=zeros(l_ag_line,1);
            
            ag_line=(g-delta+d)*i+(s+t)/2*(0:n_points)/n_points;
            x=real(ag_line);
            
            for count=2:(l_ag_line-1)
                dif=evaldiff(obj.map_self,evalinv(f,ag_line(count)));
                B_temp=B_w/conj(dif);
                
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
        
        %% Expression of the magnetic field as Fourier Series:
        
        function fourier_coeffs = fourier_series(obj,delta,V,n_points)
            temp_matrix=obj.B_func (delta,V,n_points);
            lambda=obj.s_self+obj.t_self;
            
            x=temp_matrix(3,:);
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
        
        function Kb_out = Kb_calc (obj,n_int,n_points,dm)
            KBi=zeros(1,n_int+1);
            g=obj.g_self;
            for k=1:n_int
                delta=dm/n_int*k;
                vec_temp=obj.fourier_series(delta,g,n_points);
                KBi(k+1)=4/pi*vec_temp(1)*(vec_temp(2)-vec_temp(4)/3+vec_temp(6)/5);
                
            end
            KBi(1)=KBi(2);
            
            Kb_out=mean(KBi);
            
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