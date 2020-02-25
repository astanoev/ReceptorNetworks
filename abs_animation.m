classdef abs_animation < handle
    %ABS_ANIMATION Class that animates streams of agent_based_simulation
    %objects
    %   
    
    properties 
        R_state_cols = [[255,225,225]./255;[248,102,101]./255;[0,0.8,0]];
        P_dnf_state_cols = [[223,223,253]./255;[104,108,249]./255];
        marker_size = 100;
        square_size = 3.5; % um - width and length
        folder_data = fullfile('data','agent_based_simulations');
        filename_format = 'agent_based_ics_hl=%d_g1=%g_rep=%d_part=%d.mat';
        t_mem = 50;
        % sa - 0-w/o saving, 1-save as separate pngs, 2-save directly into video
        save_animation = 0;
        newly_generated = false;
        t_current;
        % animation starting frame (dashed line)
        t_start = 1;
        t_stop = Inf;
        combined = false;
        paused = false;
        plot_conc_ratio_bkg = false;
        plot_R_active = true;
        plot_R_inactive = true;
        plot_P_dnf_active = false;
        plot_P_dnf_inactive = false;
        
        g1
        ics_hl
        rep
        n_parts = 1;        
        model = models.agent_based_model;
        
        time_steps;
        R_positions;
        P_dnf_positions;
        R_states;
        P_dnf_states;
        R_activation_events;
        anim_fig;
        anim_ax;
        anim_btn_pause_resume;
        anim_slider;
        time_fig;
        time_ax;
        time_ax_time_stamp;
        conc_ratio_fig; 
        conc_ratio_ax; 
        conc_ratio_p;
        R_plot;
        P_dnf_plot;
    end
    
    methods
        function obj = abs_animation(varargin)
            p = inputParser;
            addOptional(p, 'abs_obj', []);
            addOptional(p, 'ics_hl', 1);
            addOptional(p, 'g1', 4.95);
            addOptional(p, 'rep', 1);
            addOptional(p, 'n_parts', 1);
            parse(p, varargin{:});
            if isempty(p.Results.abs_obj)
                obj.ics_hl = p.Results.ics_hl;
                obj.g1 = p.Results.g1;
                obj.rep = p.Results.rep;
                obj.n_parts = p.Results.n_parts;
                if ~obj.load_files()
                    disp('unsuccessful loading - generating new simulation!');
                    obj.simulate();
                    if ~obj.load_files()
                        disp('error loading - stopping program!');
                        return;
                    end
                    obj.newly_generated = true;
                end
            else
                abs_obj = p.Results.abs_obj;
                obj.ics_hl = abs_obj.ics_hl;
                obj.n_parts = abs_obj.n_parts;
                obj.g1 = abs_obj.model.g1*(abs_obj.model.interaction_radius^2*pi);
                obj.R_states = single(abs_obj.R_states);
                obj.P_dnf_states = single(abs_obj.P_dnf_states);
                obj.R_activation_events = abs_obj.R_activation_events.get();
                if isempty(abs_obj.R_positions)
                    abs_obj.calculate_positions();
                end
                obj.R_positions = abs_obj.R_positions;
                obj.P_dnf_positions = abs_obj.P_dnf_positions;
                obj.model = abs_obj.model;
                obj.time_steps = abs_obj.model.time_steps;
            end
            if ~isempty(dbstack(1)); return; end
            % animate if called from editor
            obj.plot_all;
        end
        
        function success = load_files(obj)
            str = sprintf('loaded %d/%d',0,obj.n_parts);
            fprintf(str);
            for part = 1:obj.n_parts
                % load separate data files to piece together animation
                filename = sprintf(obj.filename_format, obj.ics_hl, obj.g1, obj.rep, part);
                if ~exist(fullfile(obj.folder_data, filename),'file')
                    success = 0;
                    fprintf('...missing files\n');
                    return;
                end
                s = load(fullfile(obj.folder_data, filename));
                if part==1
                    obj.R_states = s.abs_obj.R_states;
                    obj.P_dnf_states = s.abs_obj.P_dnf_states;
                    if isempty(s.abs_obj.R_positions)
                        try
                            s.abs_obj.calculate_positions();
                        catch
                            % older versions of matlab do not have the
                            % random number generator that was used for the
                            % single-molecule positions
                            success = 0;
                            fprintf('...threefry rng missing (introduced from Matlab R2018b)\n');
                            return;
                        end
                    end
                    obj.R_positions = s.abs_obj.R_positions;
                    obj.P_dnf_positions = s.abs_obj.P_dnf_positions;
                    obj.R_activation_events = s.abs_obj.R_activation_events.get();
                    obj.model = s.abs_obj.model;
                    if isempty(obj.model.Rt)
                        obj.model.Rt = size(obj.R_positions,2);
                        obj.model.P_dnf_t = size(obj.P_dnf_positions,2);
                    end
                    obj.time_steps = s.abs_obj.model.time_steps;
                    obj.square_size = s.abs_obj.model.square_size;
                else
                    obj.R_states = cat(1,obj.R_states,s.abs_obj.R_states);
                    obj.P_dnf_states = cat(1,obj.P_dnf_states,s.abs_obj.P_dnf_states);
                    if isempty(s.abs_obj.R_positions)
                        s.abs_obj.calculate_positions();
                    end
                    obj.R_positions = cat(1,obj.R_positions,s.abs_obj.R_positions);
                    obj.P_dnf_positions = cat(1,obj.P_dnf_positions,s.abs_obj.P_dnf_positions);
                    R_act = s.abs_obj.R_activation_events.get();
                    R_act(:, 1) = R_act(:, 1) + obj.time_steps;
                    obj.R_activation_events = cat(1,obj.R_activation_events, R_act);
                    obj.time_steps = obj.time_steps + s.abs_obj.time_steps;
                end
                fprintf(repmat('\b', 1, length(str)));
                str = sprintf('%d/%d', part, obj.n_parts);
                fprintf(str);
            end
            fprintf('...done\n');
            success = 1;
        end
        
        function simulate(obj)
            for part = 1:obj.n_parts
                filename = sprintf(obj.filename_format, obj.ics_hl, obj.g1, obj.rep, part);
                if exist(fullfile(obj.folder_data,filename),'file')
                    try
                        s = load(fullfile(obj.folder_data,filename));
                        abs_mini_obj = abs_mini(s.abs_obj); % for next part
                        continue;
                    catch
                        fprintf('replacing previous file..');
                    end
                end
                abs_obj = agent_based_simulation(part,obj.n_parts,'ics_hl',obj.ics_hl);
                if part==1
                    abs_obj.model.g1 = obj.g1./(abs_obj.model.interaction_radius^2*pi);
                    obj.model = abs_obj.model;
                else
                    abs_obj.set_continuation_variables(abs_mini_obj);
                end
                abs_obj.simulation();
                abs_animation.parsave(fullfile(obj.folder_data,filename), abs_obj);
                abs_mini_obj = abs_mini(abs_obj);
            end
        end
        
        function ret = plot_state(obj,t)
            ret = 0;
            if ~isvalid(obj.anim_fig); ret=1; return; end
            if obj.plot_R_active || obj.plot_R_inactive
                % update positions
                obj.R_plot.XData = obj.R_positions(t,:,1);
                % exclude markers of active or inactive molecules
                if ~obj.plot_R_active
                    obj.R_plot.XData(obj.R_states(t,:,1)==0) = nan;
                end
                if ~obj.plot_R_inactive
                    obj.R_plot.XData(obj.R_states(t,:,1)==1) = nan;
                end
                obj.R_plot.YData = obj.R_positions(t,:,2);
                % update colors according to the current state
                obj.R_plot.CData = single(squeeze(obj.R_states(t,:,:)))*obj.R_state_cols;
                size_vec = obj.marker_size*ones(size(obj.R_plot.SizeData));
                change_vec = zeros(size(obj.R_plot.SizeData));
                % increase size when new activation event happens to R mol.
                if t>1
                    ind_start = find(obj.R_activation_events(:,1)>=max(t-obj.t_mem, 2), 1, 'first');
                    ind_end = find(obj.R_activation_events(:,1)<=t, 1, 'last');
                    changed = obj.R_activation_events(ind_start:ind_end,3);
                    added_size = 300-(300/obj.t_mem)*(t-obj.R_activation_events(ind_start:ind_end,1));
                    change_vec(changed) = added_size;
                end
                obj.R_plot.SizeData = size_vec + change_vec;
            end
            if obj.plot_P_dnf_active || obj.plot_P_dnf_inactive
                obj.P_dnf_plot.XData = obj.P_dnf_positions(t,:,1);
                if ~obj.plot_P_dnf_active
                    obj.P_dnf_plot.XData(obj.P_dnf_states(t,:,1)==0) = nan;
                end
                if ~obj.plot_P_dnf_inactive
                    obj.P_dnf_plot.XData(obj.P_dnf_states(t,:,1)==1) = nan;
                end
                obj.P_dnf_plot.YData = obj.P_dnf_positions(t,:,2);
                obj.P_dnf_plot.CData = single(squeeze(obj.P_dnf_states(t,:,:)))*obj.P_dnf_state_cols;
            end
            title(obj.anim_ax, sprintf('t=%.3fs', (t-1)*obj.model.delta_t));
            if isvalid(obj.time_fig) % update dashed line
                set(obj.time_ax_time_stamp,'XData',[(t-1)*obj.model.delta_t,(t-1)*obj.model.delta_t],'YData',[0,1]);
            end
            if ~isempty(obj.conc_ratio_fig)
                if isvalid(obj.conc_ratio_fig)
                    obj.plot_conc_ratio();
                end
            end
        end
        
        function init_animation(obj)
            if ~ishandle(obj.time_fig)
                obj.combined = false;
            end
            if obj.combined
                obj.anim_fig = obj.time_fig;
                obj.anim_ax = subplot(1,2,1,'Parent', obj.anim_fig);
            else
                obj.anim_fig = figure('OuterPosition', [400, 300, 750, 700]);
                obj.anim_ax = axes('Parent', obj.anim_fig);
            end
            
            hold(obj.anim_ax,'on');
            if obj.combined
                % length of slider+button
                ax_len = 0.5;
            else
                ax_len = 0.9;
            end
            obj.anim_btn_pause_resume = uicontrol('Parent',obj.anim_fig,'Style', 'push', 'String', '||', 'Units', 'normal', 'Position', [0.1*ax_len 0.02 0.1*ax_len 0.05],'CallBack', @(source,event) pause_resume());
            obj.paused = false;
            obj.anim_slider = uicontrol('Parent',obj.anim_fig,'Style','slider','Units', 'normal', 'Position', [0.2*ax_len 0.02 0.8*ax_len 0.05],...
                'Min',1,'Max',obj.time_steps,'Value',1,'SliderStep',[1/(obj.time_steps-1),0.001]);
            
            if obj.plot_conc_ratio_bkg
                obj.plot_conc_ratio;
            end
            
            if obj.plot_R_active || obj.plot_R_inactive
                % if both are to be hidden do not generate plot
                obj.R_plot = scatter(obj.anim_ax,obj.R_positions(1,:,1),obj.R_positions(1,:,2),obj.marker_size*ones(size(obj.R_positions,2),1),...
                    single(squeeze(obj.R_states(1,:,:)))*obj.R_state_cols,'filled','MarkerEdgeColor',[1,0,0],'MarkerFaceAlpha',0.8);
                if ~obj.plot_R_active
                    obj.R_plot.XData(obj.R_states(1,:,1)==0) = nan;
                end
                if ~obj.plot_R_inactive
                    obj.R_plot.XData(obj.R_states(1,:,1)==1) = nan;
                end
            end
            if obj.plot_P_dnf_active || obj.plot_P_dnf_inactive
                % if both are to be hidden do not generate plot
                obj.P_dnf_plot = scatter(obj.anim_ax,obj.P_dnf_positions(1,:,1),obj.P_dnf_positions(1,:,2),obj.marker_size*ones(size(obj.P_dnf_positions,2),1),...
                    single(squeeze(obj.P_dnf_states(1,:,:)))*obj.P_dnf_state_cols,'filled','s','MarkerEdgeColor',[0,0.5,0.8],'MarkerFaceAlpha',0.8);
                if ~obj.plot_P_dnf_active
                    obj.P_dnf_plot.XData(obj.P_dnf_states(1,:,1)==0) = nan;
                end
                if ~obj.plot_P_dnf_inactive
                    obj.P_dnf_plot.XData(obj.P_dnf_states(1,:,1)==1) = nan;
                end
            end
            addlistener(obj.anim_slider,'Value','PostSet',@(s,e) slide());
            
            axis(obj.anim_ax,'equal');
            xlim(obj.anim_ax,[0,obj.square_size]);
            ylim(obj.anim_ax,[0,obj.square_size]);
            set(obj.anim_ax,'xtick',[]);
            set(obj.anim_ax,'ytick',[]);
            box(obj.anim_ax,'on');
            set(obj.anim_ax,'fontsize',15);
            drawnow;
            if isvalid(obj.time_fig)
                figure(obj.time_fig);
                set(obj.time_fig,'paperpositionmode','auto');
                set(obj.time_ax_time_stamp,'XData',[0,0]);
                obj.realign_axes(obj.time_ax, obj.anim_ax); % resize time_ax upon figure size change
                addlistener(obj.anim_ax, 'Position', 'PostSet', @(src,evn) obj.realign_axes(obj.time_ax, obj.anim_ax));
                set(obj.time_fig, 'SizeChangedFcn',@(src,evn) obj.realign_axes(obj.time_ax, obj.anim_ax));
            end
            drawnow;
            
            function pause_resume()
                % button callback
                if obj.paused
                    obj.paused = false;
                    obj.anim_btn_pause_resume.String = '||';
                    t = round(get(obj.anim_slider,'Value'));
                    obj.loop(t);
                else
                    obj.paused = true;
                    obj.anim_btn_pause_resume.String = '>';
                end
            end
            
            function slide()
                % slider update
                t = round(get(obj.anim_slider,'Value'));
                obj.t_current = t;
                obj.plot_state(t);
                drawnow;
            end
        end
                
        function loop(obj, t_start)
            obj.t_current = t_start;
            br = false;
            [i, video] = obj.init_save_animation();            
            while true
                if obj.t_current>obj.time_steps; break; end
                if obj.t_current>obj.t_stop; break; end
                if obj.paused
                    br = true;
                    break; 
                end
                if isvalid(obj.anim_slider)
                    % main figure update
                    set(obj.anim_slider,'Value',obj.t_current); 
                else % if figure is closed
                    br = true;
                    break
                end
                obj.t_current = obj.t_current+1; % progress time
                [i, video] = obj.update_save_animation(i, video);
            end
            obj.final_save_animation(video);
            if br == false
                obj.paused = true;
                obj.anim_btn_pause_resume.String = '>';
            end
        end
        
        function [i, video] = init_save_animation(obj)
            i = 0; video = [];
            % create folder for animation if it doesn't exist
            if obj.save_animation > 0 && ~exist(fullfile(obj.folder_data,'animations'),'dir')
                mkdir(fullfile(obj.folder_data,'animations'));
            end
            if obj.save_animation == 1
                i = 1; % track of file number
                set(obj.time_fig, 'Color', 'w');
            elseif obj.save_animation == 2
                %% Initialize video
                video = VideoWriter(fullfile(obj.folder_data,'animations','animation')); %open video file
                video.FrameRate = 30;  %can adjust this, 5 - 10 works well for me
                video.Quality = 100;
                open(video);
            end
        end
        
        function [i, video] = update_save_animation(obj, i, video)
            if obj.save_animation == 1
                skipframes = 10;
                % save every 10th frame (for speed purposes)
                if mod(i,skipframes) == 0
                    saveas(obj.time_fig, fullfile(obj.folder_data,'animations',strcat('img',sprintf('%05d', round(i/skipframes)),'.png')));
                end
                i = i+1;
            elseif obj.save_animation == 2
                frame = getframe(obj.anim_fig); %get frame
                writeVideo(video, frame);
            end
        end
        
        function final_save_animation(obj, video)
            if obj.save_animation == 2; close(video); end
        end
        
        function animation(obj)
            obj.init_animation();
            obj.loop(obj.t_start);
            waitfor(obj.anim_fig);
        end
        
        function plot_temporal(obj)
            obj.time_fig = figure();
            if obj.combined
                obj.time_ax = subplot(1,2,2,'Parent',obj.time_fig);
            else
                obj.time_ax = axes('Parent',obj.time_fig);
            end
            hold(obj.time_ax,'on'); box(obj.time_ax,'on');
            set(obj.time_fig, 'OuterPosition', [163,170,1702,715]);
            x = 0:obj.model.delta_t:(obj.time_steps-1)*obj.model.delta_t;
            xlim(obj.time_ax, [min(x),max(x)]);
            ylim(obj.time_ax, [0,1]);
            xlabel(obj.time_ax, 'Time (s)');
            ylabel(obj.time_ax, 'R_a');
            set(obj.time_ax,'fontsize',20);
            h1 = plot(obj.time_ax, x, squeeze(sum(obj.R_states(:,:,2),2))./obj.model.Rt);
            set(h1, {'color'}, num2cell(obj.R_state_cols(2,:),2));
            set(h1, 'HitTest', 'off'); % so that it does not intercept change_time_stamp()
            h2 = plot(obj.time_ax, x, squeeze(sum(obj.P_dnf_states(:,:,2),2))./obj.model.P_dnf_t);
            set(h2, {'color'}, num2cell(obj.P_dnf_state_cols(2,:),2));
            set(h2, 'HitTest', 'off');
            obj.time_ax_time_stamp = plot(obj.time_ax,[0,0],[0,1],'k--');
            set(obj.time_ax,'ButtonDownFcn',@(e,s) change_time_stamp());
            
            function change_time_stamp()
                % axes click callback - update time
                pt = get(obj.time_ax,'CurrentPoint');
                t = round(pt(1,1)/obj.model.delta_t);
                if isvalid(obj.anim_slider)
                    if obj.paused
                        set(obj.anim_slider,'Value',t);
                    else
                        obj.t_current = t;
                    end
                else
                    set(obj.time_ax_time_stamp,'XData',[(t-1)*obj.model.delta_t,(t-1)*obj.model.delta_t],'YData',[0,1]);
                end
            end
        end
        
        function realign_axes(obj, ax2, ax1) %#ok<INUSL>
            drawnow;
            pos1 = external_tools.plotboxpos(ax1); 
            pos2 = external_tools.plotboxpos(ax2); 
            set(ax2,'Units', 'Normalized', 'Position',[pos2(1),pos1(2),pos2(3),pos1(4)]); drawnow;
        end
        
        function plot_conc_ratio(obj)
            % plot as background the local P_dnf_t/R_t ratio (R0)
            time_average = 30;
            n_bins = 12;
            edges_bins = linspace(0,obj.square_size,n_bins+1);
            Rh = histcounts2(squeeze(median(obj.R_positions(max(obj.t_current-time_average+1,1):obj.t_current,:,1),1)),squeeze(median(obj.R_positions(max(obj.t_current-time_average+1,1):obj.t_current,:,2),1)),edges_bins, edges_bins);
            Ph = histcounts2(squeeze(median(obj.P_dnf_positions(max(obj.t_current-time_average+1,1):obj.t_current,:,1),1)),squeeze(median(obj.P_dnf_positions(max(obj.t_current-time_average+1,1):obj.t_current,:,2),1)),edges_bins, edges_bins);
            % calculate local ratio, and compare to main figure for SN points
            ratio = obj.model.g1*(obj.model.interaction_radius^2*pi)*(Ph./Rh)/(obj.model.P_dnf_t/obj.model.Rt);
            ratio = ratio';
            ratio(n_bins+1,n_bins+1) = 0;
            if isempty(obj.conc_ratio_fig) || ~isvalid(obj.conc_ratio_fig)
                obj.conc_ratio_fig = obj.anim_fig;
                obj.conc_ratio_ax = obj.anim_ax;
                hold(obj.conc_ratio_ax, 'on');
                [X,Y] = meshgrid(edges_bins, edges_bins);
                set(obj.conc_ratio_ax, 'xtick',[]); set(obj.conc_ratio_ax, 'ytick', []);
                % color the plot with local ratios
                obj.conc_ratio_p = pcolor(obj.conc_ratio_ax, X, Y, ratio);
                set(obj.conc_ratio_p, 'EdgeColor', 'none');
                set(obj.conc_ratio_p, 'facealpha', 0.7);
                box(obj.conc_ratio_ax, 'on');
                colormap(obj.conc_ratio_ax, flipud(plotting().rgb));
                % normalize to SN-points
                caxis(obj.conc_ratio_ax, [3.8, 4.8]);
                c = colorbar(obj.conc_ratio_ax, 'FontSize',20);
                c.Label.FontSize = 20;
                c.Label.Interpreter = 'latex';
                c.Label.String = 'local $$\tilde{\gamma}_{DNF}P_{DNF,T}/R_T\propto 1/R_0$$';
            else
                obj.conc_ratio_p.CData = ratio;
            end
            drawnow;
        end
        
        function plot_all(obj)
            obj.t_current = 1;
            obj.combined = true;
            % combine temporal plot with single-molecule animation
            obj.plot_temporal;
            obj.animation;
        end
    end
    
    methods (Static)        
        function parsave(fname, abs_obj)  %#ok<INUSD>
            % for saving data while in parallel mode
            save(fname, 'abs_obj');
        end
    end
end