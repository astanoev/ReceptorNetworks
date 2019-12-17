classdef figure3supp
    properties
        fontsize = 15;
        fig_size = [1000, 500];
    end
    
    methods
        function obj = figure3supp()
            obj.figure3supp_ab();
            obj.figure3supp_c();
            obj.figure3supp_d();
        end
        
        function figure3supp_ab(obj)
            fig = figure('Name','Exponential decay models');
            ax1 = subplot(1,3,1,'Parent',fig);
            ax1.ActivePositionProperty = 'position';
            hold(ax1,'on');
            colors = [1,0,0; [219,23,222]./255; 0,0,1];
            
            model_sdm = models.simple_dnf_model();
            mpe = experiments.multi_pulse_experiment(1);
            mpe.set_pulse_interval_min(20);
            ms = model_simulation(model_sdm, mpe);
            ms.simulate();
            ms.plot_fraction_phosphorylated('ax', ax1, 'color', colors(2,:), 'minor_format', true, 'plot_pulses', true, 'plot_input', false);
            
            EGFRpt = ms.state(strcmp(model_sdm.labels,'EGFRpt'),:);
            EGF_EGFRpt = ms.state(strcmp(model_sdm.labels,'EGF\_EGFRpt'),:);
            X_ss = max(EGFRpt+EGF_EGFRpt); % get steady-state value
            I = max(mpe.input); % get input amplitude from experiment
            
            model_edm_fast = models.exponential_decay_model();
            model_edm_fast.set_Xt_ss(X_ss, I); % set to match the steady_state
            ms = model_simulation(model_edm_fast, mpe);
            ms.simulate();
            ms.plot_fraction_phosphorylated('ax', ax1, 'color', colors(1,:), 'minor_format', true, 'plot_pulses', false, 'plot_input', false);
            
            model_edm_slow = models.exponential_decay_model();
            model_edm_slow.par.beta = 5.5*1e-4;
            model_edm_slow.set_Xt_ss(X_ss, I);
            ms = model_simulation(model_edm_slow, mpe);
            ms.simulate();
            ms.plot_fraction_phosphorylated('ax', ax1, 'color', colors(3,:), 'minor_format', true, 'plot_pulses', false, 'plot_input', false);
            obj.finish_plot_ax(ax1, 'Time (min)', 'Response');
            
            ax2 = subplot(1,3,2:3,'Parent',fig);
            ax2.ActivePositionProperty = 'position';
            hold(ax2,'on');
            
            n_pulses = 6;
            lambda = 1.5; % pulses per hour
            t_per_pulse_min = 60/lambda;
            t_pre_stimulus_min = 2;
            t_post_stimulus_min = 2.5;
            t_total_min = n_pulses*t_per_pulse_min +t_pre_stimulus_min +t_post_stimulus_min;
            sppe = experiments.stream_poisson_pulse_experiment(n_pulses, lambda, t_total_min);
            ms = model_simulation(model_edm_fast, sppe);
            models_cell = {model_edm_fast; model_sdm; model_edm_slow};
            for i = 1:3
                ms.set_model(models_cell{i});
                ms.simulate();
                ms.plot_fraction_phosphorylated('ax', ax2, 'color', colors(i,:), 'minor_format', true, 'plot_pulses', i==1, 'plot_input', i==1);
            end
            obj.finish_plot_ax(ax2, 'Time (min)', 'Response');
            obj.finish_plot_fig(fig);
        end
        
        function figure3supp_c(obj)
            if ~exist('data\active_betas.mat','file')
                % generate distributions if previously-generated files do
                % not exist
                obj.distributions_betas_stream_pulses();
            end
            s = load('data\active_betas.mat');
            fig = figure('Name','Supp Figure 3C');
            ax = axes('Parent', fig, 'Position', [0.18,0.18,0.75,0.75]);
            hold(ax,'on'); box(ax,'on');
            for j = 1:numel(s.betas)
                exponent = floor(log10(s.betas(j)));
                basic = s.betas(j)*10^(-exponent);
                histogram(ax, s.active_fraction(:,j), 0:.025:1, 'DisplayName',sprintf('$\\beta=%.1f \\times 10^{%d}$',basic,exponent));
            end
            legend(ax,'Interpreter','latex');
            obj.finish_plot_ax(ax, 'Duration of receptor activity (Fraction of time)', 'Frequency');
            obj.finish_plot_fig(fig);
        end
        
        function figure3supp_d(obj)
            if ~exist('data\active_betas.mat','file')
                % generate distributions if previously-generated files do
                % not exist
                obj.distributions_betas_stream_pulses();
            end
            s = load('data\active_betas.mat');
            fig = figure('Name','Supp Figure 3D');
            ax = axes('Parent', fig, 'Position', [0.18,0.18,0.75,0.75]);
            hold(ax,'on'); box(ax,'on');
            for j = 1:numel(s.betas)
                exponent = floor(log10(s.betas(j)));
                basic = s.betas(j)*10^(-exponent);
                histogram(ax, s.active_patches(:,j), 1:12, 'DisplayName',sprintf('$\\beta=%.1f \\times 10^{%d}$',basic,exponent));
            end
            xlim([0,13]);
            legend(ax,'Interpreter','latex');
            obj.finish_plot_ax(ax, {'Number of disjoint intervals','of receptor activity'}, 'Frequency');
            obj.finish_plot_fig(fig);
        end
        
        function distributions_betas_stream_pulses(obj)
            n_reps = 1000; % should be increased to 1000 for better statistics
            
            betas = [1.5*1e-4, 5.5*1e-4, 1e-3, 3.5*1e-3, 1e-2];
            
            active_fraction = zeros(n_reps,numel(betas)); active_patches = zeros(n_reps,numel(betas));
            
            n_pulses = 12;
            lambda = 1.5; % pulses per hour
            t_per_pulse_min = 60/lambda;
            t_pre_stimulus_min = 2;
            t_post_stimulus_min = 2.5;
            t_total_min = n_pulses*t_per_pulse_min +t_pre_stimulus_min +t_post_stimulus_min;
                        
            X_ss = 0.8152;
            I = 0.9812;
                        
            try task = getCurrentTask(); catch; task = []; end
            if ~isempty(task); tstr = num2str(task.ID); else; tstr = 'single'; end
            parfor i=1:n_reps
                disp(['lab ',tstr,': ',num2str(i)]);
                model = models.exponential_decay_model();
                sppe = experiments.stream_poisson_pulse_experiment(n_pulses, lambda, t_total_min);
                ms = model_simulation(model, sppe);
                active_fraction_i = zeros(1,numel(betas));
                active_patches_i = zeros(1,numel(betas));
                for j=1:numel(betas)
                    model.par.beta = betas(j);
                    model.set_Xt_ss(X_ss, I); % set to match the steady_state
                    ms.set_model(model);
                    ms.simulate();
                    [active_fraction_i(j), active_patches_i(j)] = ms.calc_active_fraction();
                end
                active_fraction(i,:) = active_fraction_i; 
                active_patches(i,:) = active_patches_i;
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