classdef model_2comp < models.model
    properties
        par = struct(...
            'dyn_rec',1,'Akt_t',2.5,...
            'Rt',1.3499,...
            ...% calculated from Rt-krec 2-par bifurcation
            'Rt_asymp',1.086035,'Rt_crit',1.3499,...
            'g1',3,'g2',3,'g3',0.001,...
            ...% set to 1s for redundancy with gamma pars
            'Pdnf_t',1.0,'Pnf_t',1.0,'Pnr_t',1.0,...
            'a1',0.0017,'a2',0.3,'a3',1.0,...
            'b1',36.0558,'k21',0.5,...
            'k1',0.01,...
            'kon',0.003,'Kd',5.56,'koff',0,...
            'kint',0.01,'krec',0.021,...
            'ktraf',2,'knpint',0.2,'kdeg',0.2);
    end
        
    methods
        function obj = model_2comp()
            obj.par.koff=obj.par.Kd*obj.par.kon;
            obj.labels = {'Ra', 'Ri', 'Pdnf_a', 'Rendo_a', 'Rendo_i', 'LRa', 'LRendo_a'};
            obj.init_conds = [0.05, 0.795, 0.8, 0., 0.155, 0., 0.];
        end
        
        function nc = get_nullclines(obj, Ra, Lt)
            nc = zeros(length(Ra),2);
            vals = obj.quasi_steady_state(Ra, Lt);
            Ri = vals(1,:); Pdnf_a = vals(2,:); LRa = vals(5,:);
            nc(:,1) = (obj.par.Rt.*Ri.*(obj.par.a1.*Ri +obj.par.a2.*Ra +obj.par.a3.*LRa) -obj.par.g3.*obj.par.Pnr_t.*Ra -obj.par.ktraf.*obj.par.kint.*Ra -obj.par.kon.*Ra.*Lt + 0.5*obj.par.koff.*LRa)./(obj.par.g1.*obj.par.Pdnf_t.*Ra); %Ra-nc
            nc(:,2) = Pdnf_a;
        end
        
        function [fr, br] = get_rate_balance_rates(obj, Ra, Lt)
            vals = obj.quasi_steady_state(Ra, Lt);
            Ri = vals(1,:); Pdnf_a = vals(2,:); LRa = vals(5,:);
            fr = obj.par.Rt*Ri.*(obj.par.a1.*Ri +obj.par.a2.*Ra +obj.par.a3.*LRa) +0.5*obj.par.koff.*LRa;
            br = obj.par.g1.*obj.par.Pdnf_t.*Pdnf_a.*Ra +obj.par.g3.*obj.par.Pnr_t.*Ra +obj.par.ktraf*obj.par.kint*Ra +obj.par.kon.*Ra.*Lt;
        end
        
        function vals = quasi_steady_state(obj, Ra, Lt)
            LRa = Lt/(Lt+obj.par.Kd).*ones(size(Ra)); 
            if isinf(Lt); LRa = 1; end
            
            Pdnf_a_qss = 1./(1 +obj.par.k21 +obj.par.b1*obj.par.Rt*(Ra +LRa));
            Rendo_a_qss = obj.par.ktraf*obj.par.kint*Ra./(obj.par.g2*obj.par.Pnf_t);
            Rendo_i_qss = (obj.par.ktraf*obj.par.kint*Ra +obj.par.ktraf*obj.par.knpint*obj.par.kint*(1 -Ra -LRa -obj.par.ktraf*obj.par.kint*Ra./(obj.par.g2*obj.par.Pnf_t)))./(obj.par.ktraf*obj.par.krec +obj.par.ktraf*obj.par.knpint*obj.par.kint);
            Ri_qss = 1 -Ra -Rendo_a_qss -Rendo_i_qss -LRa;

            vals = [Ri_qss; Pdnf_a_qss; Rendo_a_qss; Rendo_i_qss; LRa; zeros(size(Ra))];
        end
        
        function [dydt] = df_model(obj, t, y, experiment, qss)
            Lt = experiment.get_input(t);
            
            dydt = zeros(7,1);
            Ra = y(1);
            
            if nargin>=5
                qsss = obj.quasi_steady_state(Ra, Lt);
                vals = qss.*qsss' + (1-qss).*y(2:end);
                Ri = vals(1); Pdnf_a = vals(2);
                Rendo_a = vals(3); Rendo_i = vals(4);
                LRa = vals(5); LRendo_i = vals(6);
            else
                Ri = y(2);
                Pdnf_a = y(3);
                Rendo_a = y(4); Rendo_i = y(5);
                LRa = y(6); LRendo_i = y(7);
            end
            Pdnf_i = 1-Pdnf_a;
            
            Akt_a = 1;
            if obj.par.dyn_rec == 1 % if feedback coupling is present
                x = obj.par.Rt -LRendo_i*obj.par.Rt;
                if x < obj.par.Rt_asymp
                    Akt_a = obj.par.Akt_t;
                else
                    Akt_a = (obj.par.Rt_crit-obj.par.Rt_asymp)./(x-obj.par.Rt_asymp);
                    Akt_a = min(Akt_a, obj.par.Akt_t);
                    if Akt_a<1 % do not compensate if in bistable or preactivation region
                        Akt_a = 1;
                    end
                end
            end

            dydt(1) = obj.par.Rt*Ri*(obj.par.a1*Ri +obj.par.a2*Ra +obj.par.a3*LRa) -obj.par.g1*obj.par.Pdnf_t*Pdnf_a*Ra -obj.par.g3*obj.par.Pnr_t*Ra -obj.par.ktraf*obj.par.kint*Ra -obj.par.kon*Ra*Lt + 0.5*obj.par.koff*LRa; %Ra
            dydt(2) = -obj.par.Rt*Ri*(obj.par.a1*Ri +obj.par.a2*Ra +obj.par.a3*LRa) +obj.par.g1*obj.par.Pdnf_t*Pdnf_a*Ra +obj.par.g3*obj.par.Pnr_t*Ra -obj.par.ktraf*obj.par.knpint*obj.par.kint*Ri +Akt_a*obj.par.ktraf*obj.par.krec*Rendo_i -obj.par.kon*Ri*Lt + 0.5*obj.par.koff*LRa; %Ri
            dydt(3) = obj.par.k1*(Pdnf_i -obj.par.k21*Pdnf_a -obj.par.b1*obj.par.Rt*(Ra +LRa)*Pdnf_a); %Pdnf_a
            dydt(4) = obj.par.ktraf*obj.par.kint*Ra -obj.par.g2*obj.par.Pnf_t*Rendo_a; %Rendo_a
            dydt(5) = obj.par.g2*obj.par.Pnf_t*Rendo_a - Akt_a*obj.par.ktraf*obj.par.krec*Rendo_i +obj.par.ktraf*obj.par.knpint*obj.par.kint*Ri; %Rendo_i
            dydt(6) = -obj.par.ktraf*obj.par.kdeg*obj.par.kint*LRa + obj.par.kon*(Ra +Ri)*Lt -obj.par.koff*LRa; %LRa
            dydt(7) = obj.par.ktraf*obj.par.kdeg*obj.par.kint*LRa; %LRendo_i
        end
        
    end
end