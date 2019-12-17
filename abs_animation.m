classdef abs_animation < handle
    %ABS_ANIMATION Class that animates streams of agent_based_simulation
    %objects
    %   
    
    properties 
        EGFR_state_cols = [[255,225,225]./255;[248,102,101]./255;[0,0.8,0]];
        PTPRG_state_cols = [[223,223,253]./255;[104,108,249]./255];
        marker_size = 100;
        square_size = 3.5; % um - width and length
        folder_data = 'data'; %'\\\\billy.storage.mpi-dortmund.mpg.de\\abt2\\group\\agkoseska\\PTP_Theory_paper\\Analysis\\single_molecule\\data\\'
        filename_format = 'agent_based_ics_hl=%d_g1=%g_rep=%d_part=%d.mat';
        t_mem = 50;
        save_animation = 0;
        t_stop = Inf;
        plot_bif = false;
        g1
        ics_hl
        
        model = models.agent_based_model;
        time_steps;
        EGFR_positions;
        PTPRG_positions;
        EGFR_states;
        PTPRG_states;
        EGFR_phosphorylations;
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
        conc_ratio_bif_fig;
        conc_ratio_bif_ax;
        conc_ratio_bif_sc;
        EGFR_plot;
        PTPRG_plot;
        paused = false;
        t_current;
        gpu;
        plot_conc_ratio_bkg = false;
        plot_egfr_active = true;
        plot_egfr_inactive = false;
        plot_rg_active = true;
        plot_rg_inactive = false;
    end
    
    methods
        function obj = abs_animation(abs_obj)
            if nargin<1
                date = 'abs'; %date = '06082019_v30_v31';
                obj.ics_hl = 1;
                obj.g1 = 6;
                rep = 7;
                n_parts = 1;
                if obj.load_files(date, obj.ics_hl, obj.g1, rep, n_parts)
                    obj.plot_all;
                end
            else
                obj.EGFR_states = single(abs_obj.EGFR_states);
                obj.PTPRG_states = single(abs_obj.PTPRG_states);
                obj.EGFR_phosphorylations = abs_obj.EGFR_phosphorylations.get();
                if isempty(abs_obj.EGFR_positions)
                    abs_obj.calculate_positions();
                end
                obj.EGFR_positions = abs_obj.EGFR_positions;
                obj.PTPRG_positions = abs_obj.PTPRG_positions;
                obj.model = abs_obj.model;
                obj.time_steps = abs_obj.model.time_steps;
            end
        end
        
        function success = load_files(obj, date, ics_hl, g1, rep, n_parts)
            str = sprintf('loaded %d/%d',0,n_parts);
            fprintf(str);
            for part = 1:n_parts
                filename = sprintf(obj.filename_format,ics_hl, g1, rep, part);
                if ~exist(fullfile(obj.folder_data, date, filename),'file')
                    success = 0;
                    fprintf('...missing files\n');
                    return;
                end
                s = load(fullfile(obj.folder_data, date, filename));
                if part==1
                    obj.EGFR_states = s.abs_obj.EGFR_states;
                    obj.PTPRG_states = s.abs_obj.PTPRG_states;
                    if isempty(s.abs_obj.EGFR_positions)
                        s.abs_obj.calculate_positions();
                    end
                    obj.EGFR_positions = s.abs_obj.EGFR_positions;
                    obj.PTPRG_positions = s.abs_obj.PTPRG_positions;
                    obj.EGFR_phosphorylations = s.abs_obj.EGFR_phosphorylations.get();
                    obj.model = s.abs_obj.model;
                    if isempty(obj.model.EGFRt)
                        obj.model.EGFRt = size(obj.EGFR_positions,2);
                        obj.model.PTPRGt = size(obj.PTPRG_positions,2);
                    end
                    obj.time_steps = s.abs_obj.model.time_steps;
                    obj.square_size = s.abs_obj.model.square_size;
                else
                    obj.EGFR_states = cat(1,obj.EGFR_states,s.abs_obj.EGFR_states);
                    obj.PTPRG_states = cat(1,obj.PTPRG_states,s.abs_obj.PTPRG_states);
                    if isempty(s.abs_obj.EGFR_positions)
                        s.abs_obj.calculate_positions();
                    end
                    obj.EGFR_positions = cat(1,obj.EGFR_positions,s.abs_obj.EGFR_positions);
                    obj.PTPRG_positions = cat(1,obj.PTPRG_positions,s.abs_obj.PTPRG_positions);
                    EGFRph = s.abs_obj.EGFR_phosphorylations.get();
                    EGFRph(:, 1) = EGFRph(:, 1) + obj.time_steps;
                    obj.EGFR_phosphorylations = cat(1,obj.EGFR_phosphorylations, EGFRph);
                    obj.time_steps = obj.time_steps + s.abs_obj.time_steps;
                end
                fprintf(repmat('\b', 1, length(str)));
                str = sprintf('%d/%d',part,n_parts);
                fprintf(str);
            end
            fprintf('...done\n');
            success = 1;
        end
        
        function ret = plot_state(obj,t)
            ret = 0;
            if ~isvalid(obj.anim_fig); ret=1; return; end
            obj.EGFR_plot.XData = obj.EGFR_positions(t,:,1);
            if ~obj.plot_egfr_active
                obj.EGFR_plot.XData(obj.EGFR_states(t,:,1)==1) = obj.EGFR_plot.XData(obj.EGFR_states(t,:,1)==1) + 5;
            end
            if ~obj.plot_egfr_inactive
                obj.EGFR_plot.XData(obj.EGFR_states(t,:,1)==0) = obj.EGFR_plot.XData(obj.EGFR_states(t,:,1)==0) + 5;
            end
            obj.EGFR_plot.YData = obj.EGFR_positions(t,:,2);
            obj.EGFR_plot.CData = single(squeeze(obj.EGFR_states(t,:,:)))*obj.EGFR_state_cols;
            size_vec = obj.marker_size*ones(size(obj.EGFR_plot.SizeData));
            change_vec = zeros(size(obj.EGFR_plot.SizeData));
            if t>1
                ind_start = find(obj.EGFR_phosphorylations(:,1)>=max(t-obj.t_mem, 2), 1, 'first');
                ind_end = find(obj.EGFR_phosphorylations(:,1)<=t, 1, 'last');
                changed = obj.EGFR_phosphorylations(ind_start:ind_end,3);
                added_size = 300-(300/obj.t_mem)*(t-obj.EGFR_phosphorylations(ind_start:ind_end,1));
                change_vec(changed) = added_size;
            end
            obj.EGFR_plot.SizeData = size_vec + change_vec;
            obj.PTPRG_plot.XData = obj.PTPRG_positions(t,:,1);
            if ~obj.plot_rg_active
                obj.PTPRG_plot.XData(obj.PTPRG_states(t,:,1)==1) = obj.PTPRG_plot.XData(obj.PTPRG_states(t,:,1)==1) + 5;
            end
            if ~obj.plot_rg_inactive
                obj.PTPRG_plot.XData(obj.PTPRG_states(t,:,1)==0) = obj.PTPRG_plot.XData(obj.PTPRG_states(t,:,1)==0) + 5;
            end
            obj.PTPRG_plot.YData = obj.PTPRG_positions(t,:,2);
            obj.PTPRG_plot.CData = single(squeeze(obj.PTPRG_states(t,:,:)))*obj.PTPRG_state_cols;
            title(obj.anim_ax, sprintf('t=%.3fs', (t-1)*obj.model.delta_t));
            if isvalid(obj.time_fig)
                set(obj.time_ax_time_stamp,'XData',[(t-1)*obj.model.delta_t,(t-1)*obj.model.delta_t],'YData',[0,1]);
            end
            if ~isempty(obj.conc_ratio_fig)
                if isvalid(obj.conc_ratio_fig)
                    obj.plot_conc_ratio();
                end
            end
        end
        
        function ret = plot_state_double_size(obj,t)
            ret = 0;
            if ~isvalid(obj.anim_fig); ret=1; return; end
            obj.EGFR_plot.XData = [obj.EGFR_positions(t,:,1),obj.EGFR_positions(t,:,1),obj.square_size+obj.EGFR_positions(t,:,1),obj.square_size+obj.EGFR_positions(t,:,1)];
            obj.EGFR_plot.YData = [obj.EGFR_positions(t,:,2),obj.square_size+obj.EGFR_positions(t,:,2),obj.EGFR_positions(t,:,2),obj.square_size+obj.EGFR_positions(t,:,2)];
            obj.EGFR_plot.CData = [squeeze(obj.EGFR_states(t,:,:));squeeze(obj.EGFR_states(t,:,:));squeeze(obj.EGFR_states(t,:,:));squeeze(obj.EGFR_states(t,:,:))]*obj.EGFR_state_cols;
            size_vec = obj.marker_size*ones(size(obj.EGFR_positions,2),1);
            change_vec = zeros(size(obj.EGFR_positions,2),1);
            if t>1
                ind_start = find(obj.EGFR_phosphorylations(:,1)==max(t-obj.t_mem, 2), 1, 'first');
                ind_end = find(obj.EGFR_phosphorylations(:,1)==t, 1, 'last');
                changed = obj.EGFR_phosphorylations(ind_start:ind_end,3);
                added_size = 1000-(1000/obj.t_mem)*(t-obj.EGFR_phosphorylations(ind_start:ind_end,1));
                change_vec(changed) = added_size;
            end
            obj.EGFR_plot.SizeData = [size_vec + change_vec;size_vec + change_vec;size_vec + change_vec;size_vec + change_vec];
            obj.PTPRG_plot.XData = [obj.PTPRG_positions(t,:,1),obj.PTPRG_positions(t,:,1),obj.square_size+obj.PTPRG_positions(t,:,1),obj.square_size+obj.PTPRG_positions(t,:,1)];
            obj.PTPRG_plot.YData = [obj.PTPRG_positions(t,:,2),obj.square_size+obj.PTPRG_positions(t,:,2),obj.PTPRG_positions(t,:,2),obj.square_size+obj.PTPRG_positions(t,:,2)];
            obj.PTPRG_plot.CData = [squeeze(obj.PTPRG_states(t,:,:));squeeze(obj.PTPRG_states(t,:,:));squeeze(obj.PTPRG_states(t,:,:));squeeze(obj.PTPRG_states(t,:,:))]*obj.PTPRG_state_cols;
            title(obj.anim_ax, sprintf('t=%.3fs', (t-1)*obj.model.delta_t));
            if isvalid(obj.time_fig)
                set(obj.time_ax_time_stamp,'XData',[(t-1)*obj.model.delta_t,(t-1)*obj.model.delta_t],'YData',[0,1]);
            end
        end
        
        function init_animation(obj)
            %obj.anim_fig = figure('visible','off'); hold on;
            %obj.anim_fig = figure(); hold on;
            obj.anim_fig = obj.time_fig;
            %set(obj.anim_fig, 'OuterPosition', [400, 300, 750, 700]);
            %obj.anim_ax = gca;
            obj.anim_ax = subplot(1,2,1); hold on; 
            %figure(obj.conc_ratio_fig);
            %axes(obj.conc_ratio_ax);
            ax_len = 0.5;
            obj.anim_btn_pause_resume = uicontrol('Style', 'push', 'String', '||', 'Units', 'normal', 'Position', [0.1*ax_len 0.02 0.1*ax_len 0.05],'CallBack', @(source,event) pause_resume());
            obj.paused = false;
            obj.anim_slider = uicontrol('Style','slider','Units', 'normal', 'Position', [0.2*ax_len 0.02 0.8*ax_len 0.05],...
                'Min',1,'Max',obj.time_steps,'Value',1,'SliderStep',[1/(obj.time_steps-1),0.001]);%1/(obj.time_steps-1)]);
            obj.EGFR_plot = scatter(obj.EGFR_positions(1,:,1),obj.EGFR_positions(1,:,2),obj.marker_size*ones(size(obj.EGFR_positions,2),1),...
                single(squeeze(obj.EGFR_states(1,:,:)))*obj.EGFR_state_cols,'filled','MarkerEdgeColor',[1,0,0],'MarkerFaceAlpha',0.8);
            obj.PTPRG_plot = scatter(obj.PTPRG_positions(1,:,1),obj.PTPRG_positions(1,:,2),obj.marker_size*ones(size(obj.PTPRG_positions,2),1),...
                single(squeeze(obj.PTPRG_states(1,:,:)))*obj.PTPRG_state_cols,'filled','s','MarkerEdgeColor',[0,0.5,0.8],'MarkerFaceAlpha',0.8);
            addlistener(obj.anim_slider,'Value','PostSet',@(s,e) slide);
            
            axis equal;
            xlim([0,obj.square_size]);
            ylim([0,obj.square_size]);
            set(gca,'xtick',[]);
            set(gca,'ytick',[]);
            box on;
            if isvalid(obj.time_fig)
                figure(obj.time_fig);
                set(obj.time_ax_time_stamp,'XData',[0,0]);
            end
            drawnow;
            
            function pause_resume()
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
                t = round(get(obj.anim_slider,'Value'));
                obj.t_current = t;
                obj.plot_state(t);
                drawnow;
            end
        end
        
        function init_animation_double_size(obj)
            obj.anim_fig = figure();hold on;
            set(obj.anim_fig, 'OuterPosition', [220, 50, 1110, 1050]);
            obj.anim_ax = gca;
            %obj.anim_slider = uicontrol('Style','slider','Position',[50 10 500 10],'Min',1,'Max',obj.time_steps,'Value',1);
            obj.anim_slider = uicontrol('Style','slider','Units', 'normal', 'Position', [0.1 0.02 0.8 0.05],...
                'Min',1,'Max',obj.time_steps,'Value',1,'SliderStep',[obj.model.delta_t, obj.model.delta_t]);
            addlistener(obj.anim_slider,'Value','PostSet',@(s,e) slide);
            plot([obj.square_size,obj.square_size],[0,2*obj.square_size],'k--');
            plot([0,2*obj.square_size],[obj.square_size,obj.square_size],'k--');
            x = [obj.EGFR_positions(1,:,1),obj.EGFR_positions(1,:,1),obj.square_size+obj.EGFR_positions(1,:,1),obj.square_size+obj.EGFR_positions(1,:,1)];
            y = [obj.EGFR_positions(1,:,2),obj.square_size+obj.EGFR_positions(1,:,2),obj.EGFR_positions(1,:,2),obj.square_size+obj.EGFR_positions(1,:,2)];
            states = [squeeze(obj.EGFR_states(1,:,:));squeeze(obj.EGFR_states(1,:,:));squeeze(obj.EGFR_states(1,:,:));squeeze(obj.EGFR_states(1,:,:))];
            obj.EGFR_plot = scatter(x, y, obj.marker_size*ones(size(x,2),1),...
                states*obj.EGFR_state_cols,'filled','MarkerEdgeColor',[1,0,0],'MarkerFaceAlpha',0.6);
            x = [obj.PTPRG_positions(1,:,1),obj.PTPRG_positions(1,:,1),obj.square_size+obj.PTPRG_positions(1,:,1),obj.square_size+obj.PTPRG_positions(1,:,1)];
            y = [obj.PTPRG_positions(1,:,2),obj.square_size+obj.PTPRG_positions(1,:,2),obj.PTPRG_positions(1,:,2),obj.square_size+obj.PTPRG_positions(1,:,2)];
            states = [squeeze(obj.PTPRG_states(1,:,:));squeeze(obj.PTPRG_states(1,:,:));squeeze(obj.PTPRG_states(1,:,:));squeeze(obj.PTPRG_states(1,:,:))];
            obj.PTPRG_plot = scatter(x, y, obj.marker_size*ones(size(x,2),1),...
                states*obj.PTPRG_state_cols,'filled','^','MarkerEdgeColor',[0,0.8,0.8],'MarkerFaceAlpha',0.3);
            axis equal;
            xlim([0,2*obj.square_size]);
            ylim([0,2*obj.square_size]);
            set(gca,'xtick',[]);
            set(gca,'ytick',[]);
            box on;
            if isvalid(obj.time_fig)
                figure(obj.time_fig);
                obj.time_ax_time_stamp = plot([0,0],[0,1],'k--');
            end
            drawnow;
            
            function slide()
                t = round(get(obj.anim_slider,'Value'));
                obj.plot_state_double_size(t);
                drawnow;
            end
        end
        
        function loop(obj, t_start)
            obj.t_current = t_start;
            br = false;
            if obj.save_animation == 1
                folder = fullfile('C:/Users/',getenv('USERNAME'),'/Documents/Python Images/PTPs/temp/animations/');
                i = 1;
                set(obj.time_fig, 'Color', 'w');
                %set(obj.time_fig, 'visible', 'off');
            elseif obj.save_animation == 2
                %% Initialize video
                myVideo = VideoWriter('myVideoFile'); %open video file
                myVideo.FrameRate = 30;  %can adjust this, 5 - 10 works well for me
                myVideo.Quality = 100;
                open(myVideo);
            end
            
            while true
                if obj.t_current>obj.time_steps; break; end
                if obj.t_current>obj.t_stop; break; end
                if obj.paused
                    br = true;
                    break; 
                end
                if isvalid(obj.anim_slider)
                    set(obj.anim_slider,'Value',obj.t_current); 
                else % if figure is closed
                    br = true;
                    break
                end
                obj.t_current = obj.t_current+1;
                if obj.save_animation == 1
                    skipframes = 10;
                    if mod(i,skipframes) == 0
                        export_fig(obj.time_fig,fullfile(folder,strcat('img',sprintf('%05d', round(i/skipframes)),'.png')));
                    end
                    i = i+1;
                elseif obj.save_animation == 2
                    %frame = getframe(obj.anim_fig); %get frame
                    frame = getframe(obj.conc_ratio_fig);
                    writeVideo(myVideo, frame);
                end
                %pause(0.01);
            end
            if obj.save_animation == 2; close(myVideo); end
            if br == false
                obj.paused = true;
                obj.anim_btn_pause_resume.String = '>';
            end
        end
        
        function animation(obj)
            obj.init_animation();
            obj.loop(11670);
        end
        
        function plot_temporal(obj)
            obj.time_fig = figure();
            obj.time_ax = subplot(1,2,2);
            %obj.time_ax = axes('Parent', obj.time_fig);
            hold on; box on;
            %set(obj.time_fig, 'OuterPosition', [1150, 300, 750, 700]);
            set(obj.time_fig, 'OuterPosition', [163,170,1702,715]);
            x = 0:obj.model.delta_t:(obj.model.time_steps-1)*obj.model.delta_t;
            xlim(obj.time_ax,[min(x),max(x)]);
            ylim(obj.time_ax,[0,1]);
            xlabel('Time (s)');
            ylabel('R_a');
            set(obj.time_ax,'fontsize',20);
            h1 = plot(obj.time_ax, x, squeeze(sum(obj.EGFR_states(:,:,2),2))./obj.model.EGFRt);
            set(h1, {'color'}, num2cell(obj.EGFR_state_cols(2,:),2));
            h2 = plot(obj.time_ax, x, squeeze(sum(obj.PTPRG_states(:,:,2),2))./obj.model.PTPRGt);
            set(h2, {'color'}, num2cell(obj.PTPRG_state_cols(2,:),2));
%             h3 = plot(obj.time_ax, x,(1-obj.model.k1/(obj.model.k1+obj.model.k2))*ones(size(x)),'--');
%             set(h3, {'color'}, num2cell(obj.PTPRG_state_cols(1,:),2));
%             h4 = plot(obj.time_ax, x,(obj.model.k1/(obj.model.k1+obj.model.k2))*ones(size(x)),'--');
%             set(h4, {'color'}, num2cell(obj.PTPRG_state_cols(2,:),2));
            obj.time_ax_time_stamp = plot(obj.time_ax,[0,0],[0,1],'k--');
            set(obj.time_ax,'ButtonDownFcn',@(e,s) change_time_stamp());
            
            function change_time_stamp()
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
        
        function frac_div = calc_frac_div(obj)
            time_average = 100;
            n_bins = 12;
            edges_bins = linspace(0,5,n_bins+1);
            frac_div = zeros(obj.time_steps,1);
            Ehh = zeros(obj.time_steps,n_bins*n_bins);
            Phh = zeros(obj.time_steps,n_bins*n_bins);
            parfor t=1:obj.time_steps
                [Eh,~,~,binx,biny] = histcounts2(squeeze(median(obj.EGFR_positions(max(t-time_average+1,1):t,:,1),1)),squeeze(median(obj.EGFR_positions(max(t-time_average+1,1):t,:,2),1)),edges_bins, edges_bins);
                Ph = histcounts2(squeeze(median(obj.PTPRG_positions(max(t-time_average+1,1):t,:,1),1)),squeeze(median(obj.PTPRG_positions(max(t-time_average+1,1):t,:,2),1)),edges_bins, edges_bins);
                ratio = Ph./Eh;
                ratio_log = log2(ratio)';
                frac_div(t) = sum(ratio_log(:)<log2(obj.model.PTPRGt/obj.model.EGFRt))/length(ratio_log(:));
                Ehh(t,:) = Eh(:);
                Phh(t,:) = Ph(:);
            end
            x = 0:obj.model.delta_t:(obj.time_steps-1)*obj.model.delta_t;
            figure(); plot(x,frac_div);
            xlim([min(x),max(x)]);
            figure(); Eh_h = histogram(Ehh(:)); hold on;
            Ph_h = histogram(Phh(:)); hold on;
            Eu = mean(Ehh(:))
            Es = std(Ehh(:))
            Pu = mean(Phh(:))
            Ps = std(Phh(:))
            x = 0:.01:10;
            y_E = exp(-0.5*(x-Eu).^2./Es^2);
            y_E = max(Eh_h.BinCounts)/max(y_E).*y_E;
            y_P = exp(-0.5*(x-Pu).^2./Ps^2);
            y_P = max(Ph_h.BinCounts)/max(y_P).*y_P;
            plot(x,y_E);
            plot(x,y_P);
            figure(); pe_h = histogram(pe(:)); hold on;
            figure(); plot(x,y_P./y_E);
        end
        
        function plot_dist_distr(obj)
            epi = find(obj.EGFR_phosphorylations(:,2)>=0);
            dists = zeros(length(epi),1);
            EGFRp = obj.EGFR_phosphorylations(epi,:);
            EGFRpos = obj.EGFR_positions;
            parfor i=1:length(epi)
                t = EGFRp(i,1);
                i1 = EGFRp(i,2);
                i2 = EGFRp(i,3);
                dists(i) = sqrt(sum((EGFRpos(t,i1,:)-EGFRpos(t,i2,:)).^2));
            end
            figure; histogram(dists);
        end
        
        function plot_conc_ratio(obj)
            time_average = 30;
            n_bins = 12;
            edges_bins = linspace(0,obj.square_size,n_bins+1);
            [Eh,~,~,binx,biny] = histcounts2(squeeze(median(obj.EGFR_positions(max(obj.t_current-time_average+1,1):obj.t_current,:,1),1)),squeeze(median(obj.EGFR_positions(max(obj.t_current-time_average+1,1):obj.t_current,:,2),1)),edges_bins, edges_bins);
            Ph = histcounts2(squeeze(median(obj.PTPRG_positions(max(obj.t_current-time_average+1,1):obj.t_current,:,1),1)),squeeze(median(obj.PTPRG_positions(max(obj.t_current-time_average+1,1):obj.t_current,:,2),1)),edges_bins, edges_bins);
            ratio = obj.model.g1*(obj.model.interaction_radius^2*pi)*(Ph./Eh)/(obj.model.PTPRGt/obj.model.EGFRt);
            ratio_log = ratio';%log2(ratio)';
            ratio_log(n_bins+1,n_bins+1) = 0;
            if obj.plot_bif
                inxs = sub2ind(size(ratio),binx,biny);
                Es = zeros(n_bins*n_bins,1);
                for i=1:n_bins*n_bins
                    Es(i) = median(sum(obj.EGFR_states(max(obj.t_current-time_average,1):obj.t_current,inxs==i,2),2))/length(find(inxs==i));
                end
            end
            if isempty(obj.conc_ratio_fig) || ~isvalid(obj.conc_ratio_fig)
                [X,Y] = meshgrid(edges_bins, edges_bins);
                obj.conc_ratio_fig = obj.time_fig; %figure(); 
                figure(obj.conc_ratio_fig);
                %obj.conc_ratio_ax = axes('Parent', obj.conc_ratio_fig );
                obj.conc_ratio_ax = subplot(1,2,1);
                hold on; 
                set(gca,'xtick',[]); set(gca,'ytick',[]);
                %set(obj.conc_ratio_fig, 'OuterPosition', [1900, 300, 750, 700]);
                obj.conc_ratio_p = pcolor(X,Y,ratio_log);
                set(obj.conc_ratio_p, 'EdgeColor', 'none');
                set(obj.conc_ratio_p, 'facealpha', 0.7);
                box on;
                rgb = [ ...
                    94    79   162
                    50   136   189
                    102   194   165
                    171   221   164
                    230   245   152
                    255   255   191
                    254   224   139
                    253   174    97
                    244   109    67
                    213    62    79
                    158     1    66  ] / 255;
                colormap(flipud(rgb));
                %caxis([0.5344, 0.66]);
                %caxis([4.3, 5.0]);
                caxis([3.8, 4.8]);
                c = colorbar('FontSize',20);
                c.Label.FontSize = 20;
                c.Label.Interpreter = 'latex';
                c.Label.String = 'local $$\tilde{\gamma}_{DNF}P_{DNF,T}/R_T\propto 1/R_0$$';
                
                if obj.plot_bif
                    obj.conc_ratio_bif_fig = figure(); 
                    obj.conc_ratio_bif_ax = axes('Parent', obj.conc_ratio_bif_fig);
                    obj.conc_ratio_bif_sc = scatter(obj.conc_ratio_bif_ax, ratio(:), Es);
                    xlim([0,obj.square_size]); ylim([0,1]);
                end
            else
                obj.conc_ratio_p.CData = ratio_log;
                if obj.plot_bif
                    obj.conc_ratio_bif_sc.XData = ratio(:);
                    obj.conc_ratio_bif_sc.YData = Es;
                end
            end
            drawnow;
        end
        
        function plot_all(obj, varargin)
            p = inputParser;
            addOptional(p,'plot_rg',false);
            addOptional(p,'plot_conc_ratio',false);
            parse(p,varargin{:});
            close all
            obj.t_current = 1;
            obj.plot_temporal;
            if obj.plot_conc_ratio_bkg
                obj.plot_conc_ratio;
            end
            obj.animation;
        end
    end
end