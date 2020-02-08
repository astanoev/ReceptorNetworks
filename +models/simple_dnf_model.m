classdef simple_dnf_model < models.model
    properties
        % from ode file
        par = struct(... %2.5 - bistability, 2.957 - memory, 3.5 - rev. bistability, 4.3 - monostability
            'g1',2.957,...
            'a1',0.0017,'a2',0.3,'a3',1.0,...
            'b1',36.0558,'k21',0.5,...
            'kR',0.8,'k1',0.01,...
            'kon',0.003,'Kd',5.56,'koff',0);
    end
    
    methods
        function obj = simple_dnf_model()
            obj.par.koff=obj.par.Kd*obj.par.kon;
            obj.set_initial_cond();
            obj.labels = {'Ra', 'Pdnf_a', 'LRa'};
        end
        
        function obj = set_initial_cond(obj)
            obj.init_conds = [0.01, 0.4, 0.0];
            obj.traj_cols = [219,23,222]./255;%[1.0, 0.34, 0.34];
            obj.marker_cols = [0,1,0;1, 0, 0; 0, 0, 1];
        end
                
        function nc = get_nullclines(obj, Ra, Lt)
            nc = zeros(length(Ra),2);
            vals = obj.quasi_steady_state(Ra, Lt);
            Pdnf_a_qss = vals(1,:); LRa = vals(2,:);
            Ri = (1 - Ra - LRa);
            nc(:,1) = Ri.*(obj.par.a1.*Ri +obj.par.a2.*Ra +obj.par.a3.*LRa)./(Ra.*obj.par.g1); % Ra-nc
            nc(:,2) = Pdnf_a_qss;
        end
        
        function [fr, br] = get_rate_balance_rates(obj, Ra, Lt)
            qsss = obj.quasi_steady_state(Ra, Lt);
            Pdnf_a_qss = qsss(1,:); LRa = qsss(2,:);
            Ri = (1 - Ra - LRa);
            fr = obj.par.kR*Ri.*(obj.par.a1.*Ri +obj.par.a2.*Ra +obj.par.a3.*LRa); 
            br = obj.par.kR*obj.par.g1*Pdnf_a_qss.*Ra;
        end
        
        function vals = quasi_steady_state(obj, Ra, Lt)
            LRa = Lt./(Lt+obj.par.Kd).*ones(size(Ra)); % without depletion
            if isinf(Lt); LRa = 1; end
            Pdnf_a_qss = 1./(obj.par.k21 + 1 +obj.par.b1*(Ra+LRa));
            vals = [Pdnf_a_qss; LRa];
        end
        
        function [dydt] = df_model(obj, t, y, experiment, qss)
            Lt = experiment.get_input(t);
            
            dydt = zeros(size(y));
            Ra = y(1);
                        
            if nargin>=6
                qsss = obj.quasi_steady_state(Ra, Lt);
                vals = qss.*qsss + (1-qss).*y(2:end);
                Pdnf_a = vals(1); LRa = vals(2);
            else
                Pdnf_a = y(2);
                LRa = y(3);
            end
            Pdnf_i = 1-Pdnf_a;
            Ri = 1 -Ra -LRa;
            
            dydt(1) = obj.par.kR*(Ri*(obj.par.a1*Ri +obj.par.a2*Ra +obj.par.a3*LRa) -obj.par.g1*Pdnf_a*Ra) - obj.par.kon*Ra*Lt +0.5*obj.par.koff*LRa; %Ra
            dydt(2) = obj.par.k1*(Pdnf_i -obj.par.k21*Pdnf_a -obj.par.b1*(Ra+LRa)*Pdnf_a); %Pdnf_a
            dydt(3) = obj.par.kon*(Ra +Ri)*Lt -obj.par.koff*LRa; %LRa
            if nargin>=5; dydt = dydt.*[1;1-qss]; end
        end
    end
end