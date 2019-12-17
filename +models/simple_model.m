classdef simple_model < models.model
    properties
        % from ode file
        par = struct(... %2.0 - bistability, 2.43 - memory, 3.6 - monostability
            'g1',2.43,'g2',0.061,'g3',0.001,...%'g1',3.5,'g2',0.0,'g3',0.0,...
            'a1',0.0017,'a2',0.3,'a3',1.0,...
            'b1',36.0558,'k21',0.5,'b2',1.60248,'k34',0.22,...
            'EGFRt',0.8,'k1',0.01,'k4',2.25,'eps',0.001,...
            'kon',0.003,'Kd',5.56,'koff',0);
    end
    
    properties
        rand_num = randn(4,10000000)
    end
    
    methods
        function obj = simple_model()
            %obj@models.model();
            obj.par.koff=obj.par.Kd*obj.par.kon;
            obj.set_initial_cond();
            obj.labels = {'EGFRpt', 'PTPRGat', 'EGF\_EGFRpt', 'PTPN2at'};
        end
        
        function obj = set_initial_cond(obj)
            obj.init_conds = [0.01, 0.4, 0.0, 0.1];
            %obj.traj_cols = [0.3, 0.8, 0.4; 0.9, 0.1, 0.5];
            obj.traj_cols = [1.0, 0.34, 0.34];
            %obj.traj_cols = [0.39, 0.39, 1.0];
            %obj.marker_cols = [0, 1, 0; 1, 0, 0];%[1, 0.6, 0; 0.235, 0.67, 0.9];
            obj.marker_cols = [1, 0, 0; 0, 0, 1];
            %obj.marker_cols = [0, 0, 1];
        end
                
        function nc = get_nullclines(obj, EGFRpt, EGFfree)
            nc = zeros(length(EGFRpt),2);
            vals = obj.quasi_steady_state(EGFRpt, EGFfree);
            PTPRGat_qss = vals(1,:); EGF_EGFRpt = vals(2,:); PTPN2at_qss = vals(3,:);
            EGFRnpt = (1 - EGFRpt - EGF_EGFRpt);
            nc(:,1) = (EGFRnpt.*(obj.par.a1.*EGFRnpt +obj.par.a2.*EGFRpt +obj.par.a3.*EGF_EGFRpt) - obj.par.g3.*EGFRpt - obj.par.g2.*PTPN2at_qss.*EGFRpt)./(EGFRpt.*obj.par.g1); % EGFRpt-nc
            nc(:,2) = PTPRGat_qss;
        end
        
        function [fr, br] = get_rate_balance_rates(obj, EGFRpt, EGFfree)
            qsss = obj.quasi_steady_state(EGFRpt, EGFfree);
            PTPRGat_qss = qsss(1,:); EGF_EGFRpt = qsss(2,:); PTPN2at_qss = qsss(3,:);
            EGFRnpt = (1 - EGFRpt - EGF_EGFRpt);
            fr = obj.par.EGFRt*EGFRnpt.*(obj.par.a1.*EGFRnpt +obj.par.a2.*EGFRpt +obj.par.a3.*EGF_EGFRpt); 
            br = obj.par.EGFRt*(obj.par.g1*PTPRGat_qss.*EGFRpt +obj.par.g2*PTPN2at_qss.*EGFRpt +obj.par.g3*EGFRpt);
        end
        
        function vals = quasi_steady_state(obj, EGFRpt, EGFfree)
            %EGF_EGFRpt = (obj.EGFRt+EGFfree+obj.Kd-sqrt((obj.EGFRt+EGFfree+obj.Kd).^2-4*obj.EGFRt*EGFfree))/(2*obj.EGFRt)*ones(size(EGFRpt));
            EGF_EGFRpt = EGFfree./(EGFfree+obj.par.Kd).*ones(size(EGFRpt)); % without depletion
            if isinf(EGFfree); EGF_EGFRpt = 1; end
            PTPRGat_qss = 1./(obj.par.k21 + 1 +obj.par.b1*(EGFRpt+EGF_EGFRpt));
            PTPN2at_qss = 1 - 1./(obj.par.k34 + 1 +obj.par.b2*(EGFRpt+EGF_EGFRpt));
            vals = [PTPRGat_qss; EGF_EGFRpt; PTPN2at_qss];
        end
        
        function [dydt] = df_model(obj, t, y, experiment, qss, noise)
            %t_vec, EGF_vec
            %[~, index] = min(abs(t_vec-tt));
            %EGFfree = EGF_vec(index);
            EGFfree = experiment.get_input(t);
            
            dydt = zeros(size(y));
            EGFRpt = y(1);
                        
            if nargin>=6
                qsss = obj.quasi_steady_state(EGFRpt, EGFfree);
                vals = qss.*qsss + (1-qss).*y(2:end);
                PTPRGat = vals(1); EGF_EGFRpt = vals(2); PTPN2at = vals(3);
            else
                PTPRGat = y(2);
                EGF_EGFRpt = y(3);
                PTPN2at = y(4);
            end
            PTPRGit = 1-PTPRGat;
            PTPN2it = 1-PTPN2at;
            EGFRnpt = 1 -EGFRpt -EGF_EGFRpt;
            
            dydt(1) = obj.par.EGFRt*(EGFRnpt*(obj.par.a1*EGFRnpt +obj.par.a2*EGFRpt +obj.par.a3*EGF_EGFRpt) -obj.par.g1*PTPRGat*EGFRpt -obj.par.g2*PTPN2at*EGFRpt -obj.par.g3*EGFRpt) - obj.par.kon*EGFRpt*EGFfree +0.5*obj.par.koff*EGF_EGFRpt; %EGFRpt
            %dydt(2) = obj.par.k1*(PTPRGit -obj.par.k21*PTPRGat -obj.par.b1*obj.par.EGFRt*(EGFRpt+EGF_EGFRpt)*PTPRGat); %PTPRGat
            dydt(2) = obj.par.k1*(PTPRGit -obj.par.k21*PTPRGat -obj.par.b1*(EGFRpt+EGF_EGFRpt)*PTPRGat); %PTPRGat
            dydt(3) = obj.par.kon*(EGFRpt +EGFRnpt)*EGFfree -obj.par.koff*EGF_EGFRpt; %EGF_EGFRpt
            %dydt(4) = obj.par.k4*(obj.par.k34*PTPN2it -PTPN2at +obj.par.b2*obj.par.EGFRt*(EGFRpt+EGF_EGFRpt)*PTPN2it); %PTPN2a
            dydt(4) = obj.par.k4*(obj.par.k34*PTPN2it -PTPN2at +obj.par.b2*(EGFRpt+EGF_EGFRpt)*PTPN2it); %PTPN2a
            if nargin>=6; dydt = dydt.*[1;1-qss]; end
            
            if nargin>=7 && noise
%                 persistent temp_ind
%                 if isempty(temp_ind) temp_ind = 1; else temp_ind = temp_ind+1; end;
%                 dydt = dydt + 0.000002.*(1-[0;qss]).*y.*obj.rand_num(:,temp_ind);
                dydt = dydt + 0.02.*(1-[0;qss]).*y.*randn(4,1);
            end
        end
    end
end