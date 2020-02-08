classdef model_simulation < handle
    properties
        input_vec;
        state;
        time;
        experiment;
        model;
        sol;
    end
    
    methods
        function obj = model_simulation(model, experiment)
            if nargin>=1; obj.set_model(model); end
            if nargin>=2; obj.set_experiment(experiment); end
        end
        
        function set_experiment(obj, experiment)
            obj.experiment = experiment;
            obj.experiment.set_up_model(obj.model);
            obj.input_vec = obj.experiment.set_up_input(); 
        end
        
        function set_model(obj, model)
            obj.model = model;
        end
        
        function [active_fraction, active_patches] = calc_active_fraction(obj)
            Ra = obj.state(strcmp(obj.model.labels,'Ra'),:);
            included_frames = (obj.time>obj.experiment.pre_stimulus_min*60) & (obj.time<obj.experiment.t_total_sec - obj.experiment.post_stimulus_min*60);
            active_times = find(included_frames,1)-1 + find(Ra(included_frames)>=0.44);
            patch_times = find(included_frames,1)-1 + find(Ra(included_frames)>=0.025);
            total_times = find(included_frames,1)-1 + find(Ra(included_frames)>=0.0);
            active = sum(obj.time(active_times)-obj.time(active_times-1));
            patch = nnz((patch_times - [0,patch_times(1:end-1)+1])>0);
            total = sum(obj.time(total_times)-obj.time(total_times-1));
            active_fraction = active/total;
            active_patches = patch;
        end
        
        function plot_fraction_active(obj, varargin)
            p = inputParser;
            addOptional(p,'ax',[]);
            addOptional(p,'color',[]);
            addOptional(p,'minor_format',false);
            addOptional(p,'plot_pulses',true);
            addOptional(p,'plot_input',true);
            parse(p,varargin{:});
            ax = p.Results.ax;
            if isempty(ax)
                fig = figure();
                ax = axes('Parent',fig,'Position',[0.15,0.15,0.75,0.75]);
                hold(ax,'on');
            else
                fig = ax.Parent;
            end
            Ra = obj.state(strcmp(obj.model.labels,'Ra'),:);
            Ri = obj.state(strcmp(obj.model.labels,'Ri'),:);
            LRa = obj.state(strcmp(obj.model.labels,'LRa'),:);
            if isempty(LRa); LRa = zeros(size(Ra)); end
            % for 2-comp model
            if isempty(Ri); Rpm = ones(size(Ra)); else; Rpm = Ra+Ri+LRa; end
            time_min = obj.experiment.correct_time(obj.time);
            if p.Results.plot_input
                plot(ax, time_min,LRa./Rpm,'LineWidth',3,'Color',[0.75,0.75,0.75]);
            end
            if isempty(p.Results.color)
                plot(ax, time_min,(Ra+LRa)./Rpm,'LineWidth',3);
            else
                plot(ax, time_min,(Ra+LRa)./Rpm,'LineWidth',3,'Color',p.Results.color);
            end
            plotting().plot_fraction_active(fig, ax, obj.experiment.time_min, p.Results.minor_format, obj);
            if p.Results.plot_pulses
                obj.plot_pulses(ax);
            end
        end
        
        function p = plot_pulses(obj, ax)
            if nargin<2; ax = gca; end
            yVal = ylim(ax);
            y1 = [yVal(1),yVal(2),yVal(2),yVal(1)];
            p = [];
            n_patches = 0;
            j = 1;
            while j < length(obj.input_vec)
                ii = find(obj.input_vec(j:end)>0,1);
                if isempty(ii); break; end
                i = j + ii;
                jj = find(obj.input_vec(i:end)==0,1);
                if isempty(jj); jj = length(obj.input_vec)-i; end
                j = i + jj;
                n_patches = n_patches +1;
                p(n_patches) = patch(ax, obj.experiment.time_min([i,i,j,j]),y1,'y','FaceColor',[221,187,24]./255,'EdgeColor','none','FaceAlpha',39/255); %#ok<AGROW>
            end
            % set the patches in the background (front of ui children array)
            set(ax,'children',circshift(get(ax,'children'),-n_patches));
        end
        
        function plot_state_space(obj, varname1, varname2, ax)
            ind1 = obj.model.get_index(varname1);
            ind2 = obj.model.get_index(varname2);
            figure();hold on;
            plot(ax, obj.state(ind1,:)',obj.state(ind2,:)','LineWidth',3);
            xlabel(varname1); ylabel(varname2);
        end
        
        function simulate(obj)
            [obj.time, obj.state] = obj.time_profile(obj.model.init_conds);
        end
        
        function stochastic_simulation(obj)
            ftb = 0;
            rng('shuffle');
            X_std = 0.05;
            F = @(t,y) obj.model.df_model(t, y, obj.experiment);
            G = @(t,X) X_std.*X.*(1-X);
            dt = 0.01;
            n_periods = obj.experiment.time(end)/dt;
            if ftb && license('test', 'financial_toolbox')
                sde_solve = sde(F, G, 'StartState', obj.model.init_conds');
                [obj.state, obj.time] = sde_solve.simByEuler(n_periods,'DeltaTime', dt);
                obj.state = obj.state';
                obj.time = obj.time';
            else % if financial toolbox is missing
                [obj.state, obj.time] = obj.sde_solver_rk(F, G, obj.model.init_conds', dt, n_periods);
            end
            obj.time = obj.experiment.correct_time(obj.time);
        end
        
        function [y_m, time] = sde_solver_rk(obj, F, G, init_state, dt, n_periods) %#ok<INUSL>
            % Runge-Kutta SDE method
            time = 0:dt:((n_periods-1)*dt);
            y_m = zeros(length(init_state), n_periods);
            y_m(:,1) = init_state;
            dW = sqrt(dt)*randn(size(y_m));
            for i = 2:n_periods
                t = (i-1) * dt;
                y = y_m(:,i-1);
                y_cap = y + F(t,y)*dt + G(t, y).*dt^0.5;
                y_m(:,i) = y + F(t,y)*dt + G(t, y).*dW(:,i) + 0.5.*(G(t,y_cap)-G(t,y)).*(dW(:,i).^2-dt).*dt^(-0.5);
            end
        end
        
        function [T, Y] = time_profile(obj, init_cond)
            options = odeset('RelTol',1e-8,'AbsTol',1e-10);
            obj.sol = ode45(@(tt,y) obj.model.df_model(tt, y, obj.experiment), [obj.experiment.time(1) obj.experiment.time(end)], init_cond, options);
            T = obj.sol.x;
            Y = obj.sol.y;
        end
        
        function obj = animate(obj)
            loop = 1;
            skipframes = 100;
            traj_color = obj.model.traj_cols(1,:);
            traj_memory = 10000;
            % colormap for the trajectory trail
            traj_trail = uint8(255*[fliplr(parula(traj_memory)'); linspace(0,0.7,traj_memory)]);
            
            fig = figure();
            ax1 = subplot(1,2,1,'Parent',fig);
            box(ax1,'on');
            hold(ax1,'on');
            q = obj.model.set_state_space(ax1, [], 0);
            
            set(fig, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1.0, 0.96]);
            ax2 = subplot(1,2,2,'Parent',fig);
            set(ax2, 'XLim', [-5,35], 'YLim', [0,0.8]);
            set(ax2,'YTick',0:.2:1);
            linkaxes([ax1,ax2],'y');
            linkprop([ax1 ax2],'YLimMode');
            set(ax2, 'FontSize', 20);
            xlabel(ax2, 'Time (min)', 'fontsize', 25);
            ylabel(ax2, 'R_a', 'fontsize', 25);
            box(ax2,'on');
            hold(ax2,'on');
            obj.realign_axes(ax2, ax1); % resize ax2 upon figure size change
            addlistener(ax1, 'Position', 'PostSet', @(src,evn) obj.realign_axes(ax2, ax1));
            set(fig, 'SizeChangedFcn',@(src,evn) obj.realign_axes(ax2, ax1));
            
            LRa = obj.state(obj.model.get_index('LRa'),:);
            Pdnf_a = obj.state(obj.model.get_index('Pdnf_a'),:);
            Ra = obj.state(obj.model.get_index('Ra'),:); 
            p_pulse = obj.plot_pulses(ax2);
            set(p_pulse,'HitTest','off'); % set all elements but axis as 'non-clickable'
            p_LRa = plot(ax2, obj.time, LRa,'LineWidth',2.5,'Color',[0.75,0.75,0.75]);
            set(p_LRa,'HitTest','off');
            %trajectory_state_space = animatedline(ax1,'LineWidth',2,'Color',obj.model.traj_cols(1,:));
            time_profile = animatedline(ax2,'LineWidth',2.5,'Color',traj_color);
            set(time_profile,'HitTest','off');
            
            p_nc = [];
            p_ss_s = [];
            p_ss_u = [];
            p_sep = [];
            Ra_state_space = 0:0.001:1;
            i_updated = 1; % record latest state-space update
            i = 1;
            set(ax2,'ButtonDownFcn',@(e,s) change_time_stamp());
            while true % animate in loop
                if ~isvalid(fig); break; end % exit if figure is closed in the meantime
                if (i==1) || (i>1 && abs(LRa(i)-LRa(i_updated))>1e-6) % update state space if input (LRa) is changed
                    Lt = obj.model.tune_Lt(LRa(i)); % because state-space works with Lt as input
                    q = obj.model.set_state_space(ax1, q, Lt);
                    p_nc = obj.model.plot_nullclines(ax1, Ra_state_space, Lt, p_nc);
                    p_sep = obj.model.plot_separatrix(ax1, Ra_state_space, Lt, p_sep);
                    if i==1 % reset trajctories and points
                        p_ss = plot(ax1, repelem(Pdnf_a(1),traj_memory), repelem(Ra(1),traj_memory),'-r','LineWidth',2.5); drawnow;
                        p_ss.Color(4)=0.7; drawnow;
                        set(p_ss.Edge, 'ColorBinding','interpolated', 'ColorData', traj_trail); drawnow;
                        % add current trajectory markers in state space
                        p = plot(ax1, Pdnf_a(i), Ra(i),'o','MarkerSize',8,'MarkerFaceColor', traj_color);
                    end
                    [p_ss_s, p_ss_u] = obj.model.plot_steady_states(ax1, Ra_state_space, Lt, p_ss_s, p_ss_u);
                    i_updated = i; % record change
                end
                % update trajectories
                %addpoints(trajectory_state_space, Pdnf_a(max(1,i-skipframes):i), Ra(max(1,i-skipframes):i));
                p_ss.XData(max(1,traj_memory-i+1):end) = Pdnf_a(max(1,i-traj_memory+1):i);
                p_ss.YData(max(1,traj_memory-i+1):end) = Ra(max(1,i-traj_memory+1):i);
                set(p,'XData',Pdnf_a(i));
                set(p,'YData',Ra(i));
                addpoints(time_profile, obj.time(max(1,i-skipframes):i), Ra(max(1,i-skipframes):i));
                drawnow;
                i = i + skipframes;
                if i>length(obj.time)
                    if loop
                        i = 1;
                        clearpoints(time_profile);
                        delete(p_ss);
                        delete(p);
                    else
                        break; %#ok<UNRCH>
                    end
                end
            end
            drawnow;
            
            function change_time_stamp()
                pt = get(ax2,'CurrentPoint');
                xlims = get(ax2,'xlim');
                ylims = get(ax2,'ylim');
                t = pt(1,1); yv = pt(1,2);
                if t>xlims(2) || t<xlims(1) || yv>ylims(2) || yv<ylims(1); return; end
                [~,i] = min(abs(obj.time-t));
                clearpoints(time_profile);
                addpoints(time_profile, obj.time(1:i), Ra(1:i));
            end
        end
        
        function realign_axes(obj, ax2, ax1) %#ok<INUSL>
            drawnow;
            pos1 = external_tools.plotboxpos(ax1); 
            pos2 = external_tools.plotboxpos(ax2); 
            set(ax2,'Units', 'Normalized', 'Position',[pos2(1),pos1(2),pos2(3),pos1(4)]); drawnow;
        end
    end
end