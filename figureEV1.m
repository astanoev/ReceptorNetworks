classdef figureEV1
    properties
        fontsize = 15;
        fig_size = [1000, 500];
        colors = [[128,0,0];[212,0,0];[255,42,42];[255,128,128];[255,213,213]]./255;
        betas = [1.5*1e-4, 5.5*1e-4, 1e-3, 3.5*1e-3, 1e-2];
        n_pulses = 12;
    end
    
    methods
        function obj = figureEV1()
            obj.figureEV1a();
            obj.figureEV1b();
            obj.figureEV1c();
        end
        
        function figureEV1a(obj)
            fig = figure('Name','Figure EV1a');
            ax = axes('Parent',fig);
            hold(ax,'on');
            
            model_sdm = models.simple_dnf_model();
            mpe = experiments.multi_pulse_experiment(1);
            mpe.post_stimulus_min = 80;
            mpe.set_pulse_interval_min(20);
            ms = model_simulation(model_sdm, mpe);
            ms.simulate();
            ms.plot_fraction_active('ax', ax, 'color', [219,23,222]./255, 'minor_format', true, 'plot_pulses', true, 'plot_input', false);
            
            Ra = ms.state(strcmp(model_sdm.labels,'Ra'),:);
            LRa = ms.state(strcmp(model_sdm.labels,'LRa'),:);
            Ra_ss = max(Ra+LRa); % get steady-state value
            I = max(mpe.input); % get input amplitude from experiment
            model = models.exponential_decay_model();
            p = zeros(1,numel(obj.betas));
            for j=1:numel(obj.betas)
                model.par.beta = obj.betas(j);
                model.set_Rt_ss(Ra_ss, I); % set to match the steady_state
                ms.set_model(model);
                ms.simulate();
                time_min = ms.experiment.correct_time(ms.time);
                Ra = ms.state(strcmp(ms.model.labels,'Ra'),:);
                exponent = floor(log10(model.par.beta));
                basic = model.par.beta*10^(-exponent);
                p(j) = plot(ax, time_min, Ra, '--', 'LineWidth', 3, 'Color', obj.colors(j,:), 'DisplayName', sprintf('$\\beta=%.1f \\times 10^{%d}$',basic,exponent));
            end
            
            legend(ax, p, 'Interpreter','latex');
            
            obj.finish_plot_ax(ax, 'Time (min)', 'Response');
            obj.finish_plot_fig(fig);
        end
        
        function figureEV1b(obj)
            if ~exist('data\active_betas.mat','file')
                % generate distributions if previously-generated files do
                % not exist
                obj.distributions_betas_stream_pulses();
            end
            s = load('data\active_betas.mat');
            fig = figure('Name','Figure EV1b');
            ax = axes('Parent', fig, 'Position', [0.18,0.18,0.75,0.75]);
            hold(ax,'on'); box(ax,'on');
            bin_size = 0.025; b_min = 0; b_max = 1;
            for j = 1:numel(s.betas)
                exponent = floor(log10(s.betas(j)));
                basic = s.betas(j)*10^(-exponent);
                histogram(ax, s.active_fraction(:,j), (b_min-bin_size/2):bin_size:(b_max+bin_size/2), 'DisplayName', sprintf('$\\beta=%.1f \\times 10^{%d}$',basic,exponent), 'FaceColor', obj.colors(j,:), 'FaceAlpha', 0.8);
            end
            legend(ax,'Interpreter','latex');
            obj.finish_plot_ax(ax, 'Duration of receptor activity (Fraction of time)', 'Frequency');
            obj.finish_plot_fig(fig);
        end
        
        function figureEV1c(obj)
            if ~exist('data\active_betas.mat','file')
                % generate distributions if previously-generated files do
                % not exist
                obj.distributions_betas_stream_pulses();
            end
            s = load('data\active_betas.mat');
            fig = figure('Name','Figure EV1c');
            ax = axes('Parent', fig, 'Position', [0.18,0.18,0.75,0.75]);
            hold(ax,'on'); box(ax,'on');
            bin_size = 1; b_min = 1; b_max = obj.n_pulses;
            for j = 1:numel(s.betas)
                exponent = floor(log10(s.betas(j)));
                basic = s.betas(j)*10^(-exponent);
                histogram(ax, s.active_patches(:,j), (b_min-bin_size/2):bin_size:(b_max+bin_size/2), 'DisplayName',sprintf('$\\beta=%.1f \\times 10^{%d}$',basic,exponent), 'FaceColor', obj.colors(j,:), 'FaceAlpha', 0.8);
            end
            xlim([0,13]);
            legend(ax,'Interpreter','latex');
            obj.finish_plot_ax(ax, {'Number of disjoint intervals','of receptor activity'}, 'Frequency');
            obj.finish_plot_fig(fig);
        end
        
        function distributions_betas_stream_pulses(obj)
            n_reps = 50; % should be increased to 1000 for better statistics
            
            active_fraction = zeros(n_reps,numel(obj.betas)); active_patches = zeros(n_reps,numel(obj.betas));
            
            lambda = 1.5; % pulses per hour
            t_per_pulse_min = 60/lambda;
            t_pre_stimulus_min = 2;
            t_post_stimulus_min = 2.5;
            t_total_min = obj.n_pulses*t_per_pulse_min +t_pre_stimulus_min +t_post_stimulus_min;
                        
            % values precalculated from criticality high steady state, as in Figure EV1a
            Ra_ss = 0.8152;
            I = 0.9812;
            
            betas = obj.betas; %#ok<*PROP> % for parfor loop
            n_pulses = obj.n_pulses;
                        
            parfor i=1:n_reps
                try task = getCurrentTask(); catch; task = []; end
                if ~isempty(task); tstr = num2str(task.ID); else; tstr = 'single'; end
                disp(['lab ',tstr,': ',num2str(i)]);
                model = models.exponential_decay_model();
                spe = experiments.stream_pulse_experiment(n_pulses, lambda, t_total_min);
                spe.pulse_amplitude = I;
                ms = model_simulation(model, spe);
                active_fraction_i = zeros(1,numel(betas));
                active_patches_i = zeros(1,numel(betas));
                for j=1:numel(betas)
                    model.par.beta = betas(j);
                    model.set_Rt_ss(Ra_ss, I); % set to match the steady_state
                    ms.set_model(model);
                    ms.simulate();
                    [active_fraction_i(j), active_patches_i(j)] = ms.calc_active_fraction();
                end
                active_fraction(i,:) = active_fraction_i;  %#ok<PFOUS>
                active_patches(i,:) = active_patches_i; %#ok<PFOUS>
            end
            save('data\active_betas.mat','active_fraction','active_patches','betas');
        end
        
        function finish_plot_ax(obj, ax, xlab, ylab)
            xlabel(ax, xlab, 'fontsize', obj.fontsize);
            ylabel(ax, ylab, 'fontsize', obj.fontsize);
            set(ax, 'FontSize', obj.fontsize);
            box on;
        end
        
        function finish_plot_fig(obj, fig)
            pix_SS = get(0,'screensize');
            fig_pos = get(fig,'Position');
            fig_pos = min([fig_pos(1:2);pix_SS(3:4)-obj.fig_size-100],[],1);
            set(fig, 'Position', [fig_pos(1), fig_pos(2), obj.fig_size(1), obj.fig_size(2)]);
        end
    end
end