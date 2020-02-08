classdef model < handle
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
        
        function str = get_parameters_string(obj)
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
            obj.marker_cols = [0, 1, 0; 1, 0, 0];
        end
        
        function Lt = tune_Lt(obj, LRa)
            Lt = LRa.*obj.par.Kd./(1-LRa);
        end
        
        function ind = get_index(obj, label)
            ind = find(strncmpi(obj.labels,label,length(label)));
        end
        
        function q = set_state_space(obj, ax, q, Lt)
            xy_grid_spacing = 0.05;
            x_lim = 0.8; y_lim = 0.8;
            x_array = 0 : xy_grid_spacing : x_lim;
            y_array = 0 : xy_grid_spacing : y_lim;
            [A,B]=meshgrid(x_array,y_array);
            dA = zeros(size(A));
            dB = zeros(size(B));
            sce = experiments.simple_convergence_experiment();
            sce.set_up_input(Lt);
            Pdnf_a_index = obj.get_index('Pdnf_a');
            
            qss = ones(numel(obj.labels)-1,1);
            qss(Pdnf_a_index-1) = 0;
            for j=1:size(A,2)
                % calculate quasi-steady-state values of all the other
                % variable, given current state, except Pdnf_a
                qsss = obj.quasi_steady_state(A(1,j), Lt); % set the other variables in a quasi-steady-state
                y = [A(1,j);qsss(1:end)];
                for i=1:size(A,1)
                    y(Pdnf_a_index) = B(i,j);
                    % calculate derivative for each point in state space
                    dydt = obj.df_model(0, y, sce, qss);
                    dA(i,j) = dydt(1);
                    dB(i,j) = dydt(Pdnf_a_index);
                end
            end
            dAB = sqrt(dA.^2+dB.^2);
            if isempty(q)
                q = quiver(ax,B,A,(dB./dAB.^0.75),(dA./dAB.^0.75));hold on;
                uistack(q,'bottom');
                plotting().plot_state_space(ax, dAB, q, x_lim, y_lim);
            else
                set(q,'udata',(dB./dAB.^0.75),'vdata',(dA./dAB.^0.75));
            end
        end
        
        function p = plot_nullclines(obj, ax, Ra, Lt, p)
            nc = obj.get_nullclines(Ra, Lt);
            if isempty(p)
                p = [];
                p(1) = plot(ax, nc(:,1), Ra, 'k', 'LineWidth',1);
                p(2) = plot(ax, nc(:,2), Ra, 'k', 'LineWidth',1);
            else
                set(p(1),'XData',nc(:,1));
                set(p(2),'XData',nc(:,2));
            end
        end
        
        function p_sep = plot_separatrix(obj, ax, Ra, Lt, p_sep)
            % plot separatrix when an unstable steady state is present:
            % revert the vectors in phase space and move according to them
            % (using Euler's method),
            [Ra_ss,Pdnf_a_ss,stable] = obj.get_steady_states(Ra, Lt);
            unst_inxs = find(stable==0);
            if ~isempty(unst_inxs)
                sce = experiments.simple_convergence_experiment();
                sce.set_up_input(Lt);
                tsteps = 100;
                R_sep = zeros(tsteps,numel(unst_inxs));
                Pdnf_sep = zeros(tsteps,numel(unst_inxs));
                for j = 1:numel(unst_inxs) % in case there are multiple unstable ss
                    for step=[1,-1] % for two opposite perturbations
                        ind = tsteps/2+step+1*(step==-1); % midpoint in array
                        R_sep(ind,j) = Ra_ss(unst_inxs(j))+step*1e-3; % perturb
                        Pdnf_sep(ind,j) = Pdnf_a_ss(unst_inxs(j))+step*1e-3;
                        dt = 1; % dt of Euler's method
                        % iterate till 2/(tsteps-1)
                        for i=ind:step:(2*(step==-1)+(tsteps-1)*(step==1))
                            % same procedure as for state space drawing
                            qsss = obj.quasi_steady_state(R_sep(i), Lt);
                            y = [R_sep(i,j);qsss(1:end)];
                            y(2) = Pdnf_sep(i,j);
                            qss = ones(size(qsss));
                            qss(1) = 0;
                            dydt = obj.df_model(0, y, sce, qss);
                            % flip the derivative sign from '+' to '-' and
                            % update the next state
                            R_sep(i+step,j) = R_sep(i,j) - dt*dydt(1);
                            Pdnf_sep(i+step,j) = Pdnf_sep(i,j) - dt*dydt(2);
                        end
                    end
                end
                if isempty(p_sep)
                    p_sep = plot(ax,Pdnf_sep,R_sep,'--','LineWidth',1,'Color',[0.57,0.57,0.57]);
                else
                    set(p_sep,'XData',Pdnf_sep);
                    set(p_sep,'YData',R_sep);
                end
            else
                if ~isempty(p_sep)
                    set(p_sep,'XData',[]);
                    set(p_sep,'YData',[]);
                end
            end
        end
        
        function [Ra_ss,Pdnf_a_ss,stable] = get_steady_states(obj, Ra, Lt)
            % get all steady state using rate-balance plot, to check also
            % for stability
            [fr,br] = obj.get_rate_balance_rates(Ra, Lt);
            % intersections of FR and BR gives the steady states (x-axis)
            [Ra_ss,~] = external_tools.intersections(Ra,fr,Ra,br,1);
            Ra_ss = Ra_ss';
            tol = 1e-8;
            for i=1:length(Ra_ss)
                [fr_ss_i,br_ss_i] = obj.get_rate_balance_rates(Ra_ss(i), Lt);
                step_move = 1e-6;
                % iteratively refine the numerical solution until the
                % difference between fr and br is smaller than a tolerance
                while abs(fr_ss_i-br_ss_i)>tol
                    Ra_ss_lr = [Ra_ss(i)-step_move,Ra_ss(i)+step_move];
                    [fr_ss_ii,br_ss_ii] = obj.get_rate_balance_rates(Ra_ss_lr, Lt);
                    [fb_diff, i_lr] = min(abs(fr_ss_ii-br_ss_ii));
                    if fb_diff<abs(fr_ss_i-br_ss_i) % update Ra_ss(i)
                        Ra_ss(i) = Ra_ss_lr(i_lr);
                        fr_ss_i = fr_ss_ii(i_lr);
                        br_ss_i = br_ss_ii(i_lr);
                    else % reduce moving step
                        step_move = 1/3*step_move;
                    end
                end
            end
            [fr_ss_l,br_ss_l] = obj.get_rate_balance_rates(Ra_ss-1e-6, Lt);
            [fr_ss_r,br_ss_r] = obj.get_rate_balance_rates(Ra_ss+1e-6, Lt);
            rate_diff_ss_l = (fr_ss_l-br_ss_l)>0;
            rate_diff_ss_r = (fr_ss_r-br_ss_r)<0;
            stable = rate_diff_ss_l .* rate_diff_ss_r;
            stable(rate_diff_ss_l~=rate_diff_ss_r) = nan;
            qsss = obj.quasi_steady_state(Ra_ss, Lt);
            Pdnf_a_ss = qsss(1,:);
        end
        
        function [ss_s,ss_u] = plot_steady_states(obj, ax, Ra, Lt, ss_s, ss_u)
            [Ra_ss,Pdnf_a_ss,stable] = obj.get_steady_states(Ra, Lt);
            if isempty(ss_s)
                ss_s = plot(ax, Pdnf_a_ss(stable==1),Ra_ss(stable==1),'o','markerfacecolor',[1, 0.6, 0],'markersize',8);
            else
                set(ss_s,'XData',Pdnf_a_ss(stable==1));
                set(ss_s,'YData',Ra_ss(stable==1));
            end
            if isempty(ss_u)
                ss_u = plot(ax, Pdnf_a_ss(stable==0),Ra_ss(stable==0),'o','markerfacecolor',[0.235, 0.67, 0.9],'markersize',8);
            else
                set(ss_u,'XData',Pdnf_a_ss(stable==0));
                set(ss_u,'YData',Ra_ss(stable==0));
            end
        end
        
        function nc = get_nullclines(obj, Ra, Lt)
        end
        
        function p = plot_rate_balance(obj, ax, Ra, Lt, p)
            [fr,br] = obj.get_rate_balance_rates(Ra, Lt);
            if isempty(p)
                p = [];
                hold on;
                p(1) = plot(ax, Ra, fr, 'b', 'LineWidth',2);
                p(2) = plot(ax, Ra, br, 'r', 'LineWidth',2);
            else
                set(p(1),'YData',fr);
                set(p(2),'YData',br);
            end
        end
        
        function [fr,br] = get_rate_balance_rates(obj)
        end
        
        function vals = quasi_steady_state(obj)
        end
        
        function [dydt] = df_model(obj, t, y)
        end
    end
end