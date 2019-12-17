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
            EGFRpt = obj.state(strcmp(obj.model.labels,'EGFRpt'),:);
            included_frames = (obj.time>obj.experiment.pre_stimulus_min*60) & (obj.time<obj.experiment.t_total_sec - obj.experiment.post_stimulus_min*60);
            active_times = find(included_frames,1)-1 + find(EGFRpt(included_frames)>=0.44);
            patch_times = find(included_frames,1)-1 + find(EGFRpt(included_frames)>=0.025);
            total_times = find(included_frames,1)-1 + find(EGFRpt(included_frames)>=0.0);
            active = sum(obj.time(active_times)-obj.time(active_times-1));
            patch = nnz((patch_times - [0,patch_times(1:end-1)+1])>0);
            total = sum(obj.time(total_times)-obj.time(total_times-1));
            active_fraction = active/total;
            active_patches = patch;
        end
        
        function plot_fraction_phosphorylated(obj, varargin)
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
            EGFRpt = obj.state(strcmp(obj.model.labels,'EGFRpt'),:);
            EGFRnpt = obj.state(strcmp(obj.model.labels,'EGFRnpt'),:);
            EGF_EGFRpt = obj.state(strcmp(obj.model.labels,'EGF\_EGFRpt'),:);
            EGF_EGFRnpt = obj.state(strcmp(obj.model.labels,'EGF\_EGFRnpt'),:);
            if isempty(EGF_EGFRpt); EGF_EGFRpt = zeros(size(EGFRpt)); end
            if isempty(EGF_EGFRnpt); EGF_EGFRnpt = zeros(size(EGF_EGFRpt)); end
            if isempty(EGFRnpt); EGFRpm = ones(size(EGFRpt)); else; EGFRpm = EGFRpt+EGFRnpt+EGF_EGFRpt+EGF_EGFRnpt; end
            time_min = obj.experiment.correct_time(obj.time);
            if p.Results.plot_input
                plot(ax, time_min,(EGF_EGFRpt+EGF_EGFRnpt)./EGFRpm,'LineWidth',3,'Color',[0.75,0.75,0.75]);
            end
            if isempty(p.Results.color)
                plot(ax, time_min,(EGFRpt+EGF_EGFRpt)./EGFRpm,'LineWidth',3);
            else
                plot(ax, time_min,(EGFRpt+EGF_EGFRpt)./EGFRpm,'LineWidth',3,'Color',p.Results.color);
            end
            plotting().plot_fraction_phosphorylated(fig, ax, obj.experiment.time_min, p.Results.minor_format, obj);
            if p.Results.plot_pulses
                obj.plot_pulses(ax);
            end
        end
        
        function plot_fraction_phosphorylated_old(obj, plot_input)
            if nargin<2; plot_input = true; end
            EGFRpt = obj.state(strcmp(obj.model.labels,'EGFRpt'),:);%./obj.model.EGFRt;
            % carefull with this, if it is derived then you need to extract
            % it that way, and not use the ode variable, which is appprox.
            EGFRnpt = obj.state(strcmp(obj.model.labels,'EGFRnpt'),:);%./obj.model.EGFRt;
            EGF_EGFRpt = obj.state(strcmp(obj.model.labels,'EGF\_EGFRpt'),:);%./obj.model.EGFRt;
            EGF_EGFRnpt = obj.state(strcmp(obj.model.labels,'EGF\_EGFRnpt'),:);%./obj.model.EGFRt;
            if isempty(EGF_EGFRnpt); EGF_EGFRnpt = zeros(size(EGF_EGFRpt)); end
            if isempty(EGFRnpt); EGFRpm = ones(size(EGFRpt)); else; EGFRpm = EGFRpt+EGFRnpt+EGF_EGFRpt+EGF_EGFRnpt; end
            EGFRendo_pt = obj.state(strcmp(obj.model.labels,'EGFRendo_{pt}'),:);
            EGFRendo_npt = obj.state(strcmp(obj.model.labels,'EGFRendo_{npt}'),:);
            EGF_EGFRendo_npt = obj.state(strcmp(obj.model.labels,'EGF\_EGFRendo_{npt}'),:);
            %EGFRt = EGFRpt+EGFRnpt+EGF_EGFRpt+EGF_EGFRnpt+EGFRendo_pt+EGFRendo_npt+EGF_EGFRendo_npt;
            fig = figure; hold on;
            % correction of 1+0.5 is needed: 1 - for 0-indexing, 0.5 for
            % setting the input to the one of the closest integer time
            % point, hence it is split at .5
            if plot_input==true
                subplot(4,1,1);
                plot(obj.experiment.time_min, obj.input_vec, 'LineWidth',3, 'Color',[0.15,0.75,0.15]);
                plotting().plot_fraction_phosphorylated(fig, obj, obj.experiment.time_min);
                set(gca,'xticklabel',[]);
                set(gca,'yticklabel',[]);
                ylabel('EGF');
                xlabel('');
                subplot(4,1,2:4); hold on;
            end
            time_min = obj.experiment.correct_time(obj.time);
            plot(time_min,(EGF_EGFRpt+EGF_EGFRnpt)./EGFRpm,'LineWidth',3,'Color',[0.75,0.75,0.75]);
            %area(t,(EGF_EGFRpt+EGF_EGFRnpt)./EGFRpm);%,'LineWidth',3);
            plot(time_min,(EGFRpt+EGF_EGFRpt)./EGFRpm,'LineWidth',3);
            %plot(t,obj.EGF_free,'LineWidth',3,'Color',[0.15,0.75,0.15]);
            %plot(t,EGFRpm,'LineWidth',3);
            %plot(t,EGFRt,'LineWidth',3);
            plotting().plot_fraction_phosphorylated(fig, obj, obj.experiment.time_min);
            obj.plot_pulses();
            three_d = false;
            if three_d
                figure(50);
                view([49,34]);
                grid on;
                hold on;
                set(gca, 'XLim', [0,0.8], 'YLim', [0,0.8], 'ZLim', [0,0.8]);
                set(gca, 'FontSize', 20);
                xlabel('PTPRG activity', 'fontsize', 25);
                ylabel('PTPN2 activity', 'fontsize', 25);
                zlabel('EGFR phosphorylation', 'fontsize', 25);
                x = obj.state(strcmpi(obj.model.labels,'ptprgat'),:);
                y = obj.state(strcmpi(obj.model.labels,'ptpn2at'),:);
                z = obj.state(strcmpi(obj.model.labels,'egfrpt'),:);
                plot3(x,y,z);
                A = 0.0566; B = 0.0189; C = -0.0615; D = -0.00006;
                t = (D-B.*x-C.*y-A.*z)./(A^2+B^2+C^2);
                x2 = x + 2*t.*B;
                y2 = y + 2*t.*C;
                z2 = z + 2*t.*A;
                plot3(x2,y2,z2);
            end
        end
        
        function plot_full_state(obj)
            figure();hold on;
            t_sec = obj.time-obj.experiment.pre_stimulus+1.5;
            t_min = t_sec./60;
            h1 = plot(t_min,obj.state','LineWidth',3);
            legend(obj.model.labels);
            obj.plot_pulses();
        end
        
        function plot_pulses(obj, ax)
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
                %p(n_patches) = patch(ax, obj.experiment.time_min([i,i,j,j]),y1,'y','FaceColor',[50,50,50]./255,'EdgeColor','none','FaceAlpha',39/255);
                p(n_patches) = patch(ax, obj.experiment.time_min([i,i,j,j]),y1,'y','FaceColor',[221,187,24]./255,'EdgeColor','none','FaceAlpha',39/255);
            end
            % set the patches in the background (front of ui children array)
            set(ax,'children',circshift(get(ax,'children'),-n_patches));
        end
        
        function plot_state_space(obj, varname1, varname2)
            ind1 = obj.model.get_index(varname1);
            ind2 = obj.model.get_index(varname2);
            figure();hold on;
            plot(obj.state(ind1,:)',obj.state(ind2,:)','LineWidth',3);
            xlabel(varname1); ylabel(varname2);
        end
        
        function simulate(obj)
            [obj.time, obj.state] = obj.time_profile(obj.model.init_conds);
        end
        
        function stochastic_simulation(obj)
            %F = @(t,X) 0.1 * X;
            rng('shuffle');
            F = @(tt,y) obj.model.df_model(tt, y, obj.experiment.time, obj.input_vec);
            G = @(t,X) 0.05 .* X .* (1-X);
            dt = 0.2;
            t_max = 60*60;
            nPeriods = t_max/dt;
            sde_solve = sde(F, G, 'StartState',obj.model.init_conds');
            [obj.state, obj.time] = sde_solve.simByEuler(nPeriods,'DeltaTime', dt);
            obj.state = obj.state';
            obj.time = obj.time';
            figure; hold on;
            plot(obj.experiment.time, obj.input_vec);
            plot(obj.time,obj.state);
        end
        
        function [T, Y] = time_profile(obj, init_cond)
            if isempty(obj.sol) || true
                %options = odeset();
                options = odeset('RelTol',1e-8,'AbsTol',1e-10);%,'Events',@obj.eventfun);
                %options = odeset('RelTol',1e-3,'AbsTol',1e-4);
                noise = 0;
                qss = [0;0;1];
                T = [];
                Y = [];
                while true
                    obj.sol = ode45(@(tt,y) obj.model.df_model(tt, y, obj.experiment), [obj.experiment.time(1) obj.experiment.time(end)], init_cond, options); %, qss, noise), [0 tmax], init_cond, options);%[obj.experiment.time(1) obj.experiment.time(end)], init_cond, options); %, qss, noise), [0 tmax], init_cond, options);
                    %T = obj.experiment.time;
                    %Y = deval(obj.sol, obj.experiment.time);
                    T = [T, obj.sol.x];
                    Y = [Y, obj.sol.y];
                    break;
%                     if isempty(obj.sol.ie)
%                         break;
%                     else
%                         t_init = T(end);
%                         init_cond = Y(:,end)';
%                         inxs = find(init_cond(1:3:end)>=obj.model.par.vp);
%                         init_cond((1:3:end)*inxs) = obj.model.par.c;
%                         init_cond((2:3:end)*inxs) = init_cond((2:3:end)*inxs) + obj.model.par.d;
%                     end
                end
%                 qsss = obj.model.quasi_steady_state(obj.sol.y(1,:), EGF_free);
%                 obj.sol.y(find([0;qss]),:) = qsss(find(qss),:);
            end
            %Y = deval(obj.sol,t);
%             Y = obj.sol.y;
%             t = obj.sol.x;
        end
        
        function [position,isterm,dir] = eventfun(obj,t,y)
            position = 0;%double(sum(double(y(1:3:end)>=obj.model.par.vp))>0);
            isterm = true;
            dir = 0;  %or -1, doesn't matter
        end
        
        function obj = animate(obj, savefig, time_profile, skipframes)
            if nargin<2; savefig = 0; end
            if nargin<3; time_profile = 1; end
            if nargin<4; skipframes = 1; end
            
            dt = datetime('now');
            dt.TimeZone = 'America/Toronto';
            date = datestr(dt,'ddmmyyyy');
            folder = fullfile('C:/Users/',getenv('USERNAME'),'/Documents/Python Images/PTPs/temp/animations/');
            
            % set max time and time vector for eval
            tmax = obj.experiment.timemax;
            fr_per_sec = obj.experiment.fr_per_sec;
            t = 0:1/fr_per_sec:(tmax-0.25);
            
            % for saving pdf of certain time points
            times_snap = [25,75,200,400];
            
            fig = figure();
            s = load('mat/2xsubplots.mat'); % load print preview for double plot
            setprinttemplate(fig,s.template);
            
            if time_profile; figure(fig); ax1 = subplot(1,2,1); end; % make double plot if time profile is plotted
            hold on;
            q = obj.model.set_state_space(fig, [], 0);
            %obj.model.plot_nullclines(0:.001:1, 0, []);
            set(fig, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.6, 0.96]);
            if time_profile
                set(fig, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1.0, 0.96]);
                set(0,'CurrentFigure',fig); ax2 = subplot(1,2,2); hold on;
                set(ax2, 'XLim', [-5,30], 'YLim', [0,0.8]);
                set(ax2,'YTick',0:.2:1);
                linkaxes([ax1,ax2],'y');
                linkprop([ax1,ax2],'y'); 
                linkprop([ax1 ax2],'YLimMode');
                set(gca, 'FontSize', 20);
                xlabel('Time (min)', 'fontsize', 25);
                %ylabel('EGFRp', 'fontsize',20);
                %ylabel('R_a', 'fontsize',20);
                ylabel('EGFR phosphorylation', 'fontsize', 25);
                box on;
                for i=1:5
                    pause(0.1);
                    pos1 = plotboxpos(ax1);
                    pos2 = plotboxpos(ax2);
                    set(ax2,'Units', 'Normalized', 'Position',[pos2(1),pos1(2),pos2(3),pos1(4)]);
                end
                %setprinttemplate(fig,s.template);
            end
            
            % initialize the animation if obj.time is empty, i.e. there is no previous animation
            if isempty(obj.time)
                trajectory_state_spaces = [];
                time_profiles = [];
                %obj.init_conds = obj.init_conds.*repmat(ones(),size(obj.init_conds,1),1);
                obj.set_initial_cond();
                obj.egf_free_t = obj.set_dose_temporal(t, fr_per_sec);
                
                obj.time = zeros(size(obj.init_conds,1),length(t));
                obj.egfrp = zeros(size(obj.init_conds,1),length(t));
                obj.ptprga = zeros(size(obj.init_conds,1),length(t));
                obj.ptpn2a = zeros(size(obj.init_conds,1),length(t));
                obj.egf_egfrpt = zeros(size(obj.init_conds,1),length(t));
                % iterate over initial conditions
                for ind = 1:size(obj.init_conds,1)
                    init_cond = obj.init_conds(ind,:);
                    [obj.time(ind,:), Y] = obj.time_profile(t, init_cond, obj.egf_free_t);
                    obj.populate_state(Y, ind);
                    close(fig);
                    return;
                    if time_profile; figure(fig); subplot(1,2,1); end;
                    trajectory_state_spaces = [trajectory_state_spaces; animatedline('LineWidth',2,'Color',obj.traj_cols(ind,:))];
                    if time_profile;
                        set(0,'CurrentFigure',fig); subplot(1,2,2);
                        time_profiles = [time_profiles; animatedline('LineWidth',2,'Color',obj.traj_cols(ind,:))];
                    end;
                end
            else % load previously calculated trajectories
                trajectory_state_spaces = [];
                time_profiles = [];
                egf_egfrpt = obj.state(obj.model.get_index('egf_egfrpt'),:);
                ptprga = obj.state(obj.model.get_index('ptprga'),:);
                egfrp = obj.state(obj.model.get_index('egfrp'),:); 
%                 if time_profile;
%                     set(0,'CurrentFigure',fig); axes(ax2); hold on;
%                     plot(obj.time,egf_egfrpt,'LineWidth',2,'Color',[0.75,0.75,0.75]);
%                 end
                for ind=1:size(obj.time,1)
                    if time_profile; figure(fig); axes(ax1); end;
                    trajectory_state_spaces = [trajectory_state_spaces; animatedline('LineWidth',2,'Color',obj.model.traj_cols(ind,:))];
                    if time_profile; figure(fig); axes(ax2); end;
                    time_profiles = [time_profiles; animatedline('LineWidth',2,'Color',obj.model.traj_cols(ind,:))];
                end
            end
            
            if savefig==1
                setprinttemplate(fig,s.template);
                saveas(fig,fullfile(folder,strcat('EGFRt=',num2str(obj.model.par.EGFRt),'_init.png')),'png');
            end
            
            % plot and save calculated trajectories
            if time_profile; % plot egf_egfr in time with gray
                set(0,'CurrentFigure',fig); axes(ax2); hold on;
                ttt = (obj.time-obj.experiment.pre_stimulus+1.5)./60;% convert to mins
                %plot(ttt,egf_egfrpt,'LineWidth',2,'Color',[0.75,0.75,0.75]);
                set(gca,'children',flipud(get(gca,'children')));
                obj.plot_pulses();
                set(0,'CurrentFigure',fig); axes(ax1);
            end
            p = []; % add current trajectory markers in state space
            p_nc = [];
            p_ss_s = [];
            p_ss_u = [];
            p_sep = [];
            for ind = 1:size(obj.model.init_conds,1)
                p(ind) = plot(ptprga(1), egfrp(1),'o','MarkerFaceColor', obj.model.marker_cols(ind,:));
            end
            i_updated = 1;
            for i=1:size(obj.time,2) % animate
%                 [~,ii] = min(abs(obj.experiment.time-obj.time(i)));
%                 if i>1; [~,ii1] = min(abs(obj.experiment.time-obj.time(i-1))); end
                if (i==1) || (i>1 && abs(egf_egfrpt(i)-egf_egfrpt(i_updated))>1e-8 && mod(i,skipframes)==1) % update state space if input (egf_egfr) is changed
                %if (i==1) || (i>1 && ii>1 && abs(obj.input_vec(ii)-obj.input_vec(ii1))>1e-4) % update state space if input (egf_egfr) is changed
                    egfrpt = 0:0.001:1;
                    %if time_profile; figure(fig); axes(ax1); end;
                    %q = obj.model.set_state_space(fig, q, egf_egfrpt(i));
                    egf_free = obj.model.tune_egf_free(egf_egfrpt(i));
                    q = obj.model.set_state_space(fig, q, egf_free);
                    %q = obj.model.set_state_space(fig, q, obj.input_vec(ii));
                    %p_nc = obj.model.plot_nullclines(egfrpt, egf_egfrpt(i), p_nc);
                    p_nc = obj.model.plot_nullclines(egfrpt, egf_free, p_nc);
                    %p_nc = obj.model.plot_nullclines(egfrpt, obj.input_vec(ii), p_nc);
%                     [p_ss_s, p_ss_u] = obj.model.plot_steady_states(egfrpt, obj.input_vec(ii), p_ss_s, p_ss_u);
%                     p_sep = obj.model.plot_separatrix(egfrpt, obj.input_vec(ii), p_sep);
                    [p_ss_s, p_ss_u] = obj.model.plot_steady_states(egfrpt, egf_free, p_ss_s, p_ss_u);
                    p_sep = obj.model.plot_separatrix(egfrpt, egf_free, p_sep);
                    if i==1
                        cld = get(gca,'children');
                        als = [];
                        for j=1:length(cld)
                            if strcmp(cld(j).Type,'animatedline')==1
                                als = [als,j-1,j];
                            end
                        end
                        set(gca,'children',[cld(als);cld(setdiff(1:length(cld),als))]);% cld(5);cld(1:4);cld(6:end)]);
                    end
                    %drawnow;
                end
                
                for ind=size(obj.time,1):-1:1 % iterate over trajectories (order is reversed just because of color)
                    %if time_profile; figure(fig); axes(ax1); end;
                    addpoints(trajectory_state_spaces(ind), ptprga(ind,i), egfrp(ind,i));
                    set(p(ind),'XData',ptprga(ind,i));
                    set(p(ind),'YData',egfrp(ind,i));
                    if time_profile
                        %set(0,'CurrentFigure',fig); axes(ax2);
                        ttt = (obj.time(ind,i)-obj.experiment.pre_stimulus)./60;% convert to mins
                        addpoints(time_profiles(ind), ttt, egfrp(ind,i));%+egf_egfrpt(i));
                    end
                end
                
                if mod(i,skipframes)==1
                    drawnow; i_updated = i;
                    if savefig == 1 %&& ismember(i/fr_per_sec,times_snap)
                        setprinttemplate(fig,s.template);
                        %saveas(fig,fullfile(folder,strcat('t=',num2str(i/fr_per_sec),'_',sprintf('%05d',
                        %round(1+i/skipframes)),'.png')),'png');
                        export_fig(fig,fullfile(folder,strcat('img',sprintf('%05d', round(i/skipframes)),'.png')));
                    end
                end
            end
            if savefig==1
                setprinttemplate(fig,s.template);
                saveas(fig,fullfile(folder,strcat('EGFRt=',num2str(obj.model.par.EGFRt),'_final.png')),'png');
            end
        end
        
        function obj = animate_3d(obj, savefig, time_profile, skipframes)
            if nargin<2; savefig = 0; end
            if nargin<3; time_profile = 1; end
            if nargin<4; skipframes = 1; end
            
            dt = datetime('now');
            dt.TimeZone = 'America/Toronto';
            date = datestr(dt,'ddmmyyyy');
            folder = fullfile('C:/Users/',getenv('USERNAME'),'/Documents/Python Images/PTPs/temp/animations/');
            
            % set max time and time vector for eval
            tmax = obj.experiment.timemax;
            fr_per_sec = obj.experiment.fr_per_sec;
            t = 0:1/fr_per_sec:(tmax-0.25);
            
            % for saving pdf of certain time points
            times_snap = [25,75,200,400];
            
            fig = figure();
            s = load('mat/2xsubplots.mat'); % load print preview for double plot
            setprinttemplate(fig,s.template);
            
            if time_profile; figure(fig); ax1 = subplot(1,2,1); end; % make double plot if time profile is plotted
            view([49,34]);
            grid on;
            hold on;
            set(ax1, 'XLim', [0,0.8], 'YLim', [0,0.8], 'ZLim', [0,0.8]);
            set(ax1, 'FontSize', 20);
            xlabel('PTPRG activity', 'fontsize', 25);
            ylabel('PTPN2 activity', 'fontsize', 25);
            zlabel('EGFR phosphorylation', 'fontsize', 25);
            
            %q = obj.model.set_state_space(fig, [], 0);
            %obj.model.plot_nullclines(0:.001:1, 0, []);
            %set(fig, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.6, 0.96]);
            if time_profile
                %set(fig, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1.0, 0.96]);
                set(0,'CurrentFigure',fig); ax2 = subplot(1,2,2); hold on;
                set(ax2, 'XLim', [-5,60], 'YLim', [0,0.8]);
                set(ax2,'YTick',0:.2:1);
                set(ax2,'XTick',0:10:60);
                %linkaxes([ax1,ax2],'y');
                %linkprop([ax1,ax2],'y'); 
                linkprop([ax1 ax2],'YLimMode');
                set(ax2, 'FontSize', 20);
                xlabel('Time (min)', 'fontsize', 25);
                %ylabel('EGFRp', 'fontsize',20);
                %ylabel('R_a', 'fontsize',20);
                ylabel('EGFR phosphorylation', 'fontsize', 25);
                box on;
                for i=1:5
                    pause(0.1);
                    pos1 = plotboxpos(ax1);
                    pos2 = plotboxpos(ax2);
                    %set(ax2,'Units', 'Normalized', 'Position',[pos2(1),pos1(2),pos2(3),pos1(4)]);
                    set(ax2,'Units', 'Normalized', 'Position',[0.6,pos1(2),0.3,pos1(4)]);
                end
                set(fig, 'Units', 'Normalized', 'OuterPosition', [0, 0.1, 0.92, 0.7]);
                %setprinttemplate(fig,s.template);
            end
            
            % initialize the animation if obj.time is empty, i.e. there is no previous animation
            if isempty(obj.time)
                trajectory_state_spaces = [];
                time_profiles = [];
                %obj.init_conds = obj.init_conds.*repmat(ones(),size(obj.init_conds,1),1);
                obj.set_initial_cond();
                obj.egf_free_t = obj.set_dose_temporal(t, fr_per_sec);
                
                obj.time = zeros(size(obj.init_conds,1),length(t));
                obj.egfrp = zeros(size(obj.init_conds,1),length(t));
                obj.ptprga = zeros(size(obj.init_conds,1),length(t));
                obj.ptpn2a = zeros(size(obj.init_conds,1),length(t));
                obj.egf_egfrpt = zeros(size(obj.init_conds,1),length(t));
                % iterate over initial conditions
                for ind = 1:size(obj.init_conds,1)
                    init_cond = obj.init_conds(ind,:);
                    [obj.time(ind,:), Y] = obj.time_profile(t, init_cond, obj.egf_free_t);
                    obj.populate_state(Y, ind);
                    close(fig);
                    return;
                    if time_profile; figure(fig); subplot(1,2,1); end;
                    trajectory_state_spaces = [trajectory_state_spaces; animatedline('LineWidth',2,'Color',obj.traj_cols(ind,:))];
                    if time_profile;
                        set(0,'CurrentFigure',fig); subplot(1,2,2);
                        time_profiles = [time_profiles; animatedline('LineWidth',2,'Color',obj.traj_cols(ind,:))];
                    end;
                end
            else % load previously calculated trajectories
                trajectory_state_spaces = [];
                time_profiles = [];
                egf_egfrpt = obj.state(obj.model.get_index('egf_egfrpt'),:);
                ptprga = obj.state(obj.model.get_index('ptprga'),:);
                egfrp = obj.state(obj.model.get_index('egfrp'),:);
                ptpn2a = obj.state(obj.model.get_index('ptpn2a'),:);
%                 if time_profile;
%                     set(0,'CurrentFigure',fig); axes(ax2); hold on;
%                     plot(obj.time,egf_egfrpt,'LineWidth',2,'Color',[0.75,0.75,0.75]);
%                 end
                for ind=1:size(obj.time,1)
                    if time_profile; figure(fig); axes(ax1); end;
                    trajectory_state_spaces = [trajectory_state_spaces; animatedline('LineWidth',2,'Color',obj.model.traj_cols(ind,:))];
                    if time_profile; figure(fig); axes(ax2); end;
                    time_profiles = [time_profiles; animatedline('LineWidth',2,'Color',obj.model.traj_cols(ind,:))];
                end
            end
            
            if savefig==1
                setprinttemplate(fig,s.template);
                saveas(fig,fullfile(folder,strcat('EGFRt=',num2str(obj.model.par.EGFRt),'_init.png')),'png');
            end
            
            % plot and save calculated trajectories
            if time_profile; % plot egf_egfr in time with gray
                set(0,'CurrentFigure',fig); axes(ax2); hold on;
                ttt = (obj.time-obj.experiment.pre_stimulus+1.5)./60;% convert to mins
                %plot(ttt,egf_egfrpt,'LineWidth',2,'Color',[0.75,0.75,0.75]);
                set(gca,'children',flipud(get(gca,'children')));
                obj.plot_pulses();
                set(0,'CurrentFigure',fig); axes(ax1);
            end
            p = []; % add current trajectory markers in state space
            p_nc = [];
            p_ss_s = [];
            p_ss_u = [];
            p_sep = [];
            for ind = 1:size(obj.model.init_conds,1)
                p(ind) = plot3(ptprga(1), ptpn2a(1), egfrp(1),'o','MarkerFaceColor', obj.model.marker_cols(ind,:),'MarkerSize',10);
            end
            i_updated = 1;
            for i=1:size(obj.time,2) % animate
%                 [~,ii] = min(abs(obj.experiment.time-obj.time(i)));
%                 if i>1; [~,ii1] = min(abs(obj.experiment.time-obj.time(i-1))); end
                if (i==1) || (i>1 && abs(egf_egfrpt(i)-egf_egfrpt(i_updated))>1e-8 && mod(i,skipframes)==1) % update state space if input (egf_egfr) is changed
                %if (i==1) || (i>1 && ii>1 && abs(obj.input_vec(ii)-obj.input_vec(ii1))>1e-4) % update state space if input (egf_egfr) is changed
                    egfrpt = 0:0.001:1;
                    egf_free = obj.model.tune_egf_free(egf_egfrpt(i));
                    %q = obj.model.set_state_space(fig, q, egf_free);
                    %p_nc = obj.model.plot_nullclines(egfrpt, egf_free, p_nc);
                    [p_ss_s, p_ss_u] = obj.model.plot_steady_states_3d(egfrpt, egf_free, p_ss_s, p_ss_u);
                    %p_sep = obj.model.plot_separatrix(egfrpt, egf_free, p_sep);
                    if i==1
                        cld = get(gca,'children');
                        als = [];
                        for j=1:length(cld)
                            if strcmp(cld(j).Type,'animatedline')==1
                                als = [als,j-1,j];
                            end
                        end
                        set(gca,'children',[cld(als);cld(setdiff(1:length(cld),als))]);% cld(5);cld(1:4);cld(6:end)]);
                    end
                    %drawnow;
                end
                
                for ind=size(obj.time,1):-1:1 % iterate over trajectories (order is reversed just because of color)
                    %if time_profile; figure(fig); axes(ax1); end;
                    addpoints(trajectory_state_spaces(ind), ptprga(ind,i), ptpn2a(ind,i), egfrp(ind,i));
                    set(p(ind),'XData',ptprga(ind,i));
                    set(p(ind),'YData',ptpn2a(ind,i));
                    set(p(ind),'ZData',egfrp(ind,i));
                    if time_profile
                        %set(0,'CurrentFigure',fig); axes(ax2);
                        ttt = (obj.time(ind,i)-obj.experiment.pre_stimulus)./60;% convert to mins
                        addpoints(time_profiles(ind), ttt, egfrp(ind,i));%+egf_egfrpt(i));
                    end
                end
                
                if mod(i,skipframes)==1
                    drawnow; i_updated = i;
                    if savefig == 1 %&& ismember(i/fr_per_sec,times_snap)
                        setprinttemplate(fig,s.template);
                        %saveas(fig,fullfile(folder,strcat('t=',num2str(i/fr_per_sec),'_',sprintf('%05d',
                        %round(1+i/skipframes)),'.png')),'png');
                        export_fig(fig,fullfile(folder,strcat('img',sprintf('%05d', round(i/skipframes)),'.png')));
                    end
                end
            end
            if savefig==1
                setprinttemplate(fig,s.template);
                saveas(fig,fullfile(folder,strcat('EGFRt=',num2str(obj.model.par.EGFRt),'_final.png')),'png');
            end
        end
    end
end