classdef model < handle
    
    properties (Access = private)
        par = struct(...
        'Akt',1,'EGFRt',1.0,'PTPN2t',1.8,...
        'a1',0.0017,'a2',0.3,'a3',1.0,...
        'b1',36.0558,'b2',1.60248,'k1',0.1,...
        'g1',2.43,'g2',0.061,'g3',0.001,...
        'k21',0.5,'k34',0.22,'k4',2.25,'eps',0.001,... % k21 is k in Aneta's ode file
        'kon',0.003,'Kd',5.56,'koff',5.56*0.003,...
        'kint',0.005,'krec',0.005,'ktraf',2,'knpint',0.2,'kdeg',1);
    end
    
    properties
        init_conds
        labels
        traj_cols
        marker_cols
    end
    
    methods(Static)
        function obj = loadobj(s)
            if isstruct(s)
                obj = model(s);
            else
                obj = s;
            end
        end
    end

    methods
        % constructor
        function obj = model(model)
            if nargin>0
                fn = fieldnames(model);
                for i=1:length(fn)
                    disp('asd');
                    obj.(fn{i}) = model.(fn{i});
                end
            end
        end
        
        function savepar(obj, name)
            if nargin>1
                save(strcat('mat/models/',name),obj.par);
            else
                save(strcat('mat/models/',obj.model_name),obj.par);
            end
        end
        
        function loadpar(obj, name)
            if nargin>1
                obj.par = load(strcat('mat/models/',name));
            else
                obj.par = load(strcat('mat/models/',obj.model_name));
            end
        end
        
        function name = model_name(obj)
            fullmodel = strsplit(class(obj),'.');
            name = fullmodel(end);
        end
        
        function str = get_parameters_string(obj, EGF)
            str = 'quasi_potential_landscape';
            par_fields = fields(obj.par);
            for i=1:length(par_fields)
                par_name = par_fields{i};
                par_val = obj.par.(par_name);
                str = strcat(str,'_',par_name,'_',num2str(par_val));
            end
        end
                
        function obj = set_initial_cond(obj, init_conds)
            if nargin < 2; init_conds = [0., 1., 0.5, 0., 0., 0., 0., 0.]; end
            obj.init_conds = init_conds;
            obj.traj_cols = [0.3, 0.8, 0.4; 0.9, 0.1, 0.5];
            obj.marker_cols = [0, 1, 0; 1, 0, 0];%[1, 0.6, 0; 0.235, 0.67, 0.9];
        end
        
        function egf_free = tune_egf_free(obj, target_egf_egfr_p)
            egf_free = target_egf_egfr_p.*obj.par.Kd./(1-target_egf_egfr_p);
        end
        
        function ind = get_index(obj, label)
            ind = find(strncmpi(obj.labels,label,length(label)));
        end
        
        function q = set_state_space(obj, fig, q, EGFfree)
            ax = fig.CurrentAxes;
            xy_grid_spacing = 0.05;
            x_lim = 0.8; y_lim = 0.8;
            x_array = 0 : xy_grid_spacing : x_lim;
            y_array = 0 : xy_grid_spacing : y_lim;
            [A,B]=meshgrid(x_array,y_array);
            dA = zeros(size(A));
            dB = zeros(size(B));
            exp = experiments.experiment(obj);
            exp.input = EGFfree;
            RG_index = obj.get_index('PTPRGat');
            for i=1:size(A,1)
                for j=1:size(A,2)
                    % calculate quasi-steady-state values of all the other
                    % variable, given current state, except PTPRG
                    qsss = obj.quasi_steady_state(A(i,j), EGFfree); % set the other variables in a quasi-steady-state
                    y = [A(i,j);qsss(1:end)];
                    y(RG_index) = B(i,j);
                    qss = ones(size(qsss));
                    qss(RG_index-1) = 0;
                    % calculate derivative for each point in state space
                    %dydt = obj.df_model(0, y, 0, EGFfree, qss);
                    dydt = obj.df_model(0, y, exp, qss);
                    dA(i,j) = dydt(1);
                    dB(i,j) = dydt(RG_index);
                end
            end
            dAB = sqrt(dA.^2+dB.^2);
            if isempty(q)
                q = quiver(ax,B,A,(dB./dAB.^0.75),(dA./dAB.^0.75));hold on;
            else
                set(q,'udata',(dB./dAB.^0.75),'vdata',(dA./dAB.^0.75));
            end
            uistack(q,'bottom');
            plotting().plot_state_space(fig, obj, EGFfree, dAB, q, x_lim, y_lim);
        end
        
        function p = plot_nullclines(obj, EGFRpt, EGFfree, p)
            nc = obj.get_nullclines(EGFRpt, EGFfree);
            if isempty(p)
                p = [];
                p(1) = plot(nc(:,1), EGFRpt, 'k', 'LineWidth',1);
                p(2) = plot(nc(:,2), EGFRpt, 'k', 'LineWidth',1);
            else
                set(p(1),'XData',nc(:,1));
                set(p(2),'XData',nc(:,2));
            end
        end
        
        function p_sep = plot_separatrix(obj, EGFRpt, EGFfree, p_sep)
            % plot separatrix when an unstable steady state is present:
            % revert the vectors in phase space and move according to them
            % (using Euler's method),
            [EGFRpt_ss,PTPRGat_ss,~,stable] = obj.get_steady_states(EGFRpt, EGFfree);
            if ~isempty(find(stable==0, 1))
                exp = experiments.experiment(obj);
                exp.input = EGFfree;
                tsteps = 1000;
                egfr = zeros(1,tsteps);
                ptprg = zeros(1,tsteps);
                for step=[1,-1] % for two opposite perturbations
                    ind = tsteps/2+step+1*(step==-1); % midpoint in array
                    egfr(ind) = EGFRpt_ss(stable==0)+step*1e-3; % perturb
                    ptprg(ind) = PTPRGat_ss(stable==0)+step*1e-3;
                    dt = 1; % dt of Euler's method
                    % iterate till 2/(tsteps-1)
                    for i=ind:step:(2*(step==-1)+(tsteps-1)*(step==1))
                        % same procedure as for state space drawing
                        qsss = obj.quasi_steady_state(egfr(i), EGFfree);
                        y = [egfr(i);qsss(1:end)];
                        y(2) = ptprg(i);
                        qss = ones(size(qsss));
                        qss(1) = 0;
                        dydt = obj.df_model(0, y, exp, qss);
                        % flip the derivative sign from '+' to '-' and
                        % update the next state
                        egfr(i+step) = egfr(i) - dt*dydt(1);
                        ptprg(i+step) = ptprg(i) - dt*dydt(2);
                    end
                end
                if isempty(p_sep)
                    p_sep = plot(ptprg,egfr,'--','LineWidth',1,'Color',[0.57,0.57,0.57]);
                else
                    set(p_sep,'XData',ptprg);
                    set(p_sep,'YData',egfr);
                end
            else
                if ~isempty(p_sep)
                    set(p_sep,'XData',[]);
                    set(p_sep,'YData',[]);
                end
            end
        end
        
        function [EGFRpt_ss,PTPRGat_ss,PTPN2at_ss,stable] = get_steady_states(obj, EGFRpt, EGFfree)
            % get all steady state using rate-balance plot, to check also
            % for stability
            [fr,br] = obj.get_rate_balance_rates(EGFRpt, EGFfree);
            % intersections of FR and BR gives the steady states (x-axis)
            [EGFRpt_ss,~] = external_tools.intersections(EGFRpt,fr,EGFRpt,br,1);
            EGFRpt_ss = EGFRpt_ss';
            stable = ones(size(EGFRpt_ss));
            for i=1:length(EGFRpt_ss)
                [~,ind] = min(abs(EGFRpt-EGFRpt_ss(i))); % closest element of EGFRpt
                if ind>1 && fr(ind-1)-br(ind-1)<0; % if left of s.s. BR>FR
                    stable(i) = 0;
                end
                if ind<length(EGFRpt) && fr(ind+1)-br(ind+1)>0; % if right of s.s. FR>BR
                    stable(i) = 0;
                end
            end
            qsss = obj.quasi_steady_state(EGFRpt_ss, EGFfree);
            PTPRGat_ss = qsss(1,:);
            PTPN2at_ss = qsss(3,:);
            %nc = obj.get_nullclines(EGFRpt, EGFfree);
            %[y,x] = intersections(EGFRpt,nc(:,1),EGFRpt,nc(:,2),1);
        end
        
        function [ss_s,ss_u] = plot_steady_states(obj, EGFRpt, EGFfree, ss_s, ss_u)
            [EGFRpt_ss,PTPRGat_ss,~,stable] = obj.get_steady_states(EGFRpt, EGFfree);
            if isempty(ss_s)
                ss_s = plot(PTPRGat_ss(stable==1),EGFRpt_ss(stable==1),'o','markerfacecolor',[1, 0.6, 0],'markersize',8);
            else
                set(ss_s,'XData',PTPRGat_ss(stable==1));
                set(ss_s,'YData',EGFRpt_ss(stable==1));
            end
            if isempty(ss_u)
                ss_u = plot(PTPRGat_ss(stable==0),EGFRpt_ss(stable==0),'o','markerfacecolor',[0.235, 0.67, 0.9],'markersize',8);
            else
                set(ss_u,'XData',PTPRGat_ss(stable==0));
                set(ss_u,'YData',EGFRpt_ss(stable==0));
            end
        end
        
        function [ss_s,ss_u] = plot_steady_states_3d(obj, EGFRpt, EGFfree, ss_s, ss_u)
            [EGFRpt_ss,PTPRGat_ss,PTPN2at_ss,stable] = obj.get_steady_states(EGFRpt, EGFfree);
            if isempty(ss_s)
                ss_s = plot3(PTPRGat_ss(stable==1),PTPN2at_ss(stable==1),EGFRpt_ss(stable==1),'o','markerfacecolor',[1, 0.6, 0],'markersize',8);
            else
                set(ss_s,'XData',PTPRGat_ss(stable==1));
                set(ss_s,'YData',PTPN2at_ss(stable==1));
                set(ss_s,'ZData',EGFRpt_ss(stable==1));
            end
            if isempty(ss_u)
                ss_u = plot3(PTPRGat_ss(stable==0),PTPN2at_ss(stable==0),EGFRpt_ss(stable==0),'o','markerfacecolor',[0.235, 0.67, 0.9],'markersize',8);
            else
                set(ss_u,'XData',PTPRGat_ss(stable==0));
                set(ss_u,'YData',PTPN2at_ss(stable==0));
                set(ss_u,'ZData',EGFRpt_ss(stable==0));
            end
        end
        
        function nc = get_nullclines(obj, EGFRpt, EGFfree)
            nc = zeros(length(EGFRpt),2);
            vals = quasi_steady_state(obj, EGFRpt, EGFfree);
            PTPRGat_qss = vals(1,:); EGF_EGFRpt = vals(2,:); PTPN2at_qss = vals(3,:);
            EGFRnpt = (1 - EGFRpt - EGF_EGFRpt);
            nc(:,1) = (EGFRnpt.*(obj.par.a1.*EGFRnpt +obj.par.a2.*EGFRpt +obj.par.a3.*EGF_EGFRpt) - obj.par.g3.*EGFRpt - obj.par.g2.*PTPN2at_qss.*EGFRpt)./(EGFRpt.*obj.par.g1); % EGFRpt-nc
            nc(:,2) = PTPRGat_qss;
        end
        
        function p = plot_rate_balance(obj, EGFRpt, EGFfree, p)
            [fr,br] = obj.get_rate_balance_rates(EGFRpt, EGFfree);
            if isempty(p)
                p = [];
                hold on;
                p(1) = plot(EGFRpt, fr, 'b', 'LineWidth',2);
                p(2) = plot(EGFRpt, br, 'r', 'LineWidth',2);
            else
                set(p(1),'YData',fr);
                set(p(2),'YData',br);
            end
        end
        
        function [fr,br] = get_rate_balance_rates(obj, EGFRpt, EGFfree)
            [PTPRGat_qss, EGF_EGFRpt, PTPN2at_qss] = deal(obj.quasi_steady_state(EGFRpt, EGFfree));
            EGFRnpt = (1 - EGFRpt - EGF_EGFRpt);
            fr = obj.par.EGFRt*EGFRnpt.*(obj.par.a1.*EGFRnpt +obj.par.a2.*EGFRpt +obj.par.a3.*EGF_EGFRpt); % EGFRpt-nc
            br = obj.par.EGFRt*(obj.par.g1*PTPRGat_qss.*EGFRpt +obj.par.g2*PTPN2at_qss.*EGFRpt +obj.par.g3*EGFRpt);
        end
        
        function vals = quasi_steady_state(obj, EGFRpt, EGFfree)
            %EGF_EGFRpt = (obj.EGFRt+EGFfree+obj.Kd-sqrt((obj.EGFRt+EGFfree+obj.Kd).^2-4*obj.EGFRt*EGFfree))/(2*obj.EGFRt)*ones(size(EGFRpt));
            EGF_EGFRpt = EGFfree/(EGFfree+obj.par.Kd).*ones(size(EGFRpt)); % without depletion
            PTPRGat_qss = 1./(obj.par.k21 + 1 +obj.par.b1*(EGFRpt+EGF_EGFRpt));
            PTPN2at_qss = 1 - 1./(obj.par.k34 + 1 +obj.par.b2*(EGFRpt+EGF_EGFRpt));
            vals = [PTPRGat_qss; EGF_EGFRpt; PTPN2at_qss];
        end
        
        function y = find_steady_state(obj, EGFRpt_init, input, return_saddle)
            options = optimoptions('fsolve','Display','none','TolFun',1e-12,'TolX',1e-12,'Algorithm','levenberg-marquardt');
            EGFRpt_init_init = EGFRpt_init;
            has_minus = 1;
            while has_minus>0 && EGFRpt_init>=0
                init_cond = [EGFRpt_init, obj.quasi_steady_state(EGFRpt_init, input)'];
                has_minus = sum(init_cond<0);
                EGFRpt_init = EGFRpt_init - 0.05;
            end
            func = @(y) obj.df_model(0, y, 0, input);
            
%             options = odeset('RelTol',1e-12,'AbsTol',1e-12);
%             t = 0:.1:500;
%             func = @(tt, y) obj.df_model(tt, y, t, input*ones(size(t)))';
%             sol = ode45(func, [0 500], init_cond, options);
%             y = sol.y(:,end)';
%             return;
            [y,~,exitflag,~,J] = fsolve(func,init_cond,options);
            if (exitflag<=0) || sum(real(eig(J))>1e-12)>0
            %if (exitflag~=1 && exitflag~=3)% || sum(real(eig(J))>=0)>0
                if return_saddle; y=input; return; end
                EGFRpt_init = 1-EGFRpt_init;
                has_minus = 1;
                while has_minus>0 && EGFRpt_init>=0
                    init_cond = [EGFRpt_init, obj.quasi_steady_state(EGFRpt_init, input)'];
                    has_minus = sum(init_cond<0);
                    EGFRpt_init = EGFRpt_init - 0.05;
                end
                %init_cond = [EGFRpt_init, obj.quasi_steady_state(EGFRpt_init, input)'];
                if has_minus==0
                    [y,~,exitflag,~,J] = fsolve(func,init_cond,options);
                else
                    EGFRpt_init = EGFRpt_init_init;
                    init_cond = [EGFRpt_init, obj.quasi_steady_state(EGFRpt_init, input)'];
                end
                if (exitflag<=0) || sum(real(eig(J))>1e-12)>0
                    options = odeset('RelTol',1e-12,'AbsTol',1e-12);
                    t = 0:.1:500;
                    func = @(tt, y) obj.df_model(tt, y, t, input*ones(size(t)));
                    sol = ode45(func, [0 500], init_cond, options);
                    y = sol.y(:,end)';
%                     assert(norm(y(2:end)-obj.quasi_steady_state(y(1), input)')<1e-12, 'quasi-steady-state not matching final steady state!');
                end
            end
        end
        
        function [] = eventfun(t,y,x)
            
        end
        
        function [dydt] = df_model(obj, tt, y, t_vec, EGF_vec)
            [~, index] = min(abs(t_vec-tt));
            EGFfree = EGF_vec(index);
            
            dydt = zeros(8,1);
            
            EGFRpt = y(1); EGFRnpt = y(2);
            PTPRGat = y(3); PTPRGit = 1-PTPRGat;
            EGFRendo_pt = y(4); EGFRendo_npt = y(5);
            EGF_EGFRpt = y(6); EGF_EGFRendo_pt = y(7);
            PTPN2a = obj.par.PTPN2t;
            %EGFRnpt = 1 -EGFRpt -EGFRendo_pt -EGFRendo_npt -EGF_EGFRpt -EGF_EGFRendo_pt;
            
            a = 3.252e+09; b = -21.68; c = 44.27; d = -2.765;
            x = obj.par.EGFRt-EGF_EGFRendo_pt*obj.par.EGFRt;
            Akt_a = a.*exp(b.*x)+c.*exp(d.*x);
            if obj.Akt==0; Akt_ef = 1; else; Akt_ef = Akt_a; end

            dydt(1) = obj.par.EGFRt*EGFRnpt*(obj.par.a1*EGFRnpt +obj.par.a2*EGFRpt +obj.par.a3*EGF_EGFRpt) -obj.par.g1*PTPRGat*EGFRpt -obj.par.g3*EGFRpt -obj.par.ktraf*obj.par.kint*EGFRpt -obj.par.kon*EGFRpt*EGFfree + 0.5*obj.par.koff*EGF_EGFRpt; %EGFRpt
            dydt(2) = -obj.par.EGFRt*EGFRnpt*(obj.par.a1*EGFRnpt +obj.par.a2*EGFRpt +obj.par.a3*EGF_EGFRpt) +obj.par.g1*PTPRGat*EGFRpt +obj.par.g3*EGFRpt -obj.par.ktraf*obj.par.knpint*obj.par.kint*EGFRnpt +Akt_ef*obj.par.ktraf*obj.par.krec*EGFRendo_npt -obj.par.kon*EGFRnpt*EGFfree + 0.5*obj.par.koff*EGF_EGFRpt; %EGFRnpt
            dydt(3) = obj.par.k1*(PTPRGit -obj.par.k21*PTPRGat -obj.par.b1*obj.par.EGFRt*(EGFRpt +EGF_EGFRpt)*PTPRGat); %PTPRGat
            dydt(4) = obj.par.ktraf*obj.par.kint*EGFRpt -obj.par.g2*PTPN2a*EGFRendo_pt; %EGFRendo_pt
            dydt(5) = obj.par.g2*PTPN2a*EGFRendo_pt - Akt_ef*obj.par.ktraf*obj.par.krec*EGFRendo_npt +obj.par.ktraf*obj.par.knpint*obj.par.kint*EGFRnpt; %EGFRendo_npt
            dydt(6) = -obj.par.ktraf*obj.par.kdeg*obj.par.kint*EGF_EGFRpt + obj.par.kon*(EGFRpt +EGFRnpt)*EGFfree -obj.par.koff*EGF_EGFRpt; %EGF_EGFRpt
            dydt(7) = obj.par.ktraf*obj.par.kdeg*obj.par.kint*EGF_EGFRpt; %EGF_EGFRendo_npt
            %dydt(8) = kf*freq*(Akt_t-Akt_a) -kb*Akt_a; %Akt_a
        end
    end
end