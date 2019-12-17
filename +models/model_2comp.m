classdef model_2comp < models.model
    properties
        par = struct(...
            'Akt',1,'EGFRt',1.3499,'PTPN2t',1.0,...
            'g1',3,'g2',3,'g3',0.001,...
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
            obj.labels = {'EGFRpt', 'EGFRnpt', 'PTPRGat', 'EGFRendo_{pt}', 'EGFRendo_{npt}', 'EGF\_EGFRpt', 'EGF\_EGFRendo_{npt}', 'Akt_a'};
            obj.init_conds = [0.22, 0.625, 0.1, 0., 0.155, 0., 0.];
        end
        
        function nc = get_nullclines(obj, EGFRpt, EGFfree)
            nc = zeros(length(EGFRpt),2);
            vals = obj.quasi_steady_state(EGFRpt, EGFfree);
            EGFRnpt = vals(1,:); PTPRGat = vals(2,:); EGF_EGFRpt = vals(5,:);
            nc(:,1) = (obj.par.EGFRt.*EGFRnpt.*(obj.par.a1.*EGFRnpt +obj.par.a2.*EGFRpt +obj.par.a3.*EGF_EGFRpt) -obj.par.g3.*EGFRpt -obj.par.ktraf.*obj.par.kint.*EGFRpt -obj.par.kon.*EGFRpt.*EGFfree + 0.5*obj.par.koff.*EGF_EGFRpt)./(obj.par.g1.*EGFRpt); %EGFRpt-nc
            nc(:,2) = PTPRGat;
        end
        
        function [fr, br] = get_rate_balance_rates(obj, EGFRpt, EGFfree)
            vals = obj.quasi_steady_state(EGFRpt, EGFfree);
            EGFRnpt = vals(1,:); PTPRGat = vals(2,:); EGF_EGFRpt = vals(5,:);
            fr = obj.par.EGFRt*EGFRnpt.*(obj.par.a1.*EGFRnpt +obj.par.a2.*EGFRpt +obj.par.a3.*EGF_EGFRpt) +0.5*obj.par.koff.*EGF_EGFRpt;
            br = obj.par.g1.*PTPRGat.*EGFRpt +obj.par.g3*EGFRpt +obj.par.ktraf*obj.par.kint*EGFRpt +obj.par.kon.*EGFRpt.*EGFfree;
        end
        
        function vals = quasi_steady_state(obj, EGFRpt, EGFfree)
            %EGF_EGFRpt = (obj.EGFRt+EGFfree+obj.Kd-sqrt((obj.EGFRt+EGFfree+obj.Kd).^2-4*obj.EGFRt*EGFfree))/(2*obj.EGFRt)*ones(size(EGFRpt));
            EGF_EGFRpt = EGFfree/(EGFfree+obj.par.Kd).*ones(size(EGFRpt)); % without depletion
            if isinf(EGFfree); EGF_EGFRpt = 1; end
            
            PTPRGat_qss = 1./(1 +obj.par.k21 +obj.par.b1*obj.par.EGFRt*(EGFRpt +EGF_EGFRpt));
            EGFRendo_pt_qss = obj.par.ktraf*obj.par.kint*EGFRpt./(obj.par.g2*obj.par.PTPN2t);
            EGFRendo_npt_qss = (obj.par.ktraf*obj.par.kint*EGFRpt +obj.par.ktraf*obj.par.knpint*obj.par.kint*(1 -EGFRpt -EGF_EGFRpt -obj.par.ktraf*obj.par.kint*EGFRpt./(obj.par.g2*obj.par.PTPN2t)))./(obj.par.ktraf*obj.par.krec +obj.par.ktraf*obj.par.knpint*obj.par.kint);
            EGFRnpt_qss = 1 -EGFRpt -EGFRendo_pt_qss -EGFRendo_npt_qss -EGF_EGFRpt;

            vals = [EGFRnpt_qss; PTPRGat_qss; EGFRendo_pt_qss; EGFRendo_npt_qss; EGF_EGFRpt; zeros(size(EGFRpt))];
        end
        
        function [dydt] = df_model(obj, t, y, experiment, qss)
            EGFfree = experiment.get_input(t);
            
            dydt = zeros(7,1);
            EGFRpt = y(1);
            
            if nargin>=5
                qsss = obj.quasi_steady_state(EGFRpt, EGFfree);
                vals = qss.*qsss' + (1-qss).*y(2:end);
                EGFRnpt = vals(1); PTPRGat = vals(2);
                EGFRendo_pt = vals(3); EGFRendo_npt = vals(4);
                EGF_EGFRpt = vals(5); EGF_EGFRendo_npt = vals(6);
            else
                EGFRnpt = y(2);
                PTPRGat = y(3);
                EGFRendo_pt = y(4); EGFRendo_npt = y(5);
                EGF_EGFRpt = y(6); EGF_EGFRendo_npt = y(7);
            end
            PTPRGit = 1-PTPRGat;
            PTPN2a = obj.par.PTPN2t;

            Akt_a = 1;
            if obj.par.Akt == 1 % if feedback coupling is present
                x = obj.par.EGFRt-EGF_EGFRendo_npt*obj.par.EGFRt;
                EGFRt_asymp = 1.086035;
                Akt_a = (obj.par.EGFRt-EGFRt_asymp)./(x-EGFRt_asymp);
                Akt_t = 2.5;
                Akt_a = min(Akt_a, Akt_t);
                if Akt_a<0; Akt_a = Akt_t; end
            end

            dydt(1) = obj.par.EGFRt*EGFRnpt*(obj.par.a1*EGFRnpt +obj.par.a2*EGFRpt +obj.par.a3*EGF_EGFRpt) -obj.par.g1*PTPRGat*EGFRpt -obj.par.g3*EGFRpt -obj.par.ktraf*obj.par.kint*EGFRpt -obj.par.kon*EGFRpt*EGFfree + 0.5*obj.par.koff*EGF_EGFRpt; %EGFRpt
            dydt(2) = -obj.par.EGFRt*EGFRnpt*(obj.par.a1*EGFRnpt +obj.par.a2*EGFRpt +obj.par.a3*EGF_EGFRpt) +obj.par.g1*PTPRGat*EGFRpt +obj.par.g3*EGFRpt -obj.par.ktraf*obj.par.knpint*obj.par.kint*EGFRnpt +Akt_a*obj.par.ktraf*obj.par.krec*EGFRendo_npt -obj.par.kon*EGFRnpt*EGFfree + 0.5*obj.par.koff*EGF_EGFRpt; %EGFRnpt
            dydt(3) = obj.par.k1*(PTPRGit -obj.par.k21*PTPRGat -obj.par.b1*obj.par.EGFRt*(EGFRpt +EGF_EGFRpt)*PTPRGat); %PTPRGat
            dydt(4) = obj.par.ktraf*obj.par.kint*EGFRpt -obj.par.g2*PTPN2a*EGFRendo_pt; %EGFRendo_pt
            dydt(5) = obj.par.g2*PTPN2a*EGFRendo_pt - Akt_a*obj.par.ktraf*obj.par.krec*EGFRendo_npt +obj.par.ktraf*obj.par.knpint*obj.par.kint*EGFRnpt; %EGFRendo_npt
            dydt(6) = -obj.par.ktraf*obj.par.kdeg*obj.par.kint*EGF_EGFRpt + obj.par.kon*(EGFRpt +EGFRnpt)*EGFfree -obj.par.koff*EGF_EGFRpt; %EGF_EGFRpt
            dydt(7) = obj.par.ktraf*obj.par.kdeg*obj.par.kint*EGF_EGFRpt; %EGF_EGFRendo_npt
        end
        
    end
end