classdef figure3
    properties
        fontsize = 20;
        fig_size = [700, 650];
    end
    
    methods
        function obj = figure3()
            obj.figure3a();
            obj.figure3b();
            obj.figure3c();
            obj.figure3d();
            obj.figure3e();
        end
        
        function figure3a(obj)
            fig = figure('Name','Figure 3A');
            ax = axes('Parent',fig,'Position',[0.15,0.15,0.75,0.75]);
            hold(ax,'on');
            
            px = plot_xpp();
            regimes = fields(px.gamma_dnf);
            model = models.simple_dnf_model;
            
            for i = 1:numel(regimes)
                % set regime
                gamma_dnf = px.gamma_dnf.(regimes{i});
                model.par.g1 = gamma_dnf; 
                % plot simulation with double-pulse
                obj.figure3_double_pulse(model, ax, px.colors.(regimes{i}), i==1);
            end
            obj.finish_plot_ax(ax, 'Time (min)', 'Response');
            obj.finish_plot_fig(fig);
        end
        
        function figure3b(obj)
            if ~exist('data\active_fractions.mat','file')
                % generate distributions if previously-generated files do
                % not exist
                obj.distributions_stream_pulses();
            end
            s = load('data\active_fractions.mat');
            fig = figure('Name','Figure 3B');
            ax = axes('Parent', fig, 'Position', [0.18,0.18,0.75,0.75]);
            hold(ax,'on'); box(ax,'on');
            px = plot_xpp();
            histogram(ax, s.active_fraction_mono, 0:.025:1, 'FaceColor', px.colors.monostable);
            histogram(ax, s.active_fraction_revbst, 0:.025:1, 'FaceColor', px.colors.rev_bistable);
            histogram(ax, s.active_fraction_mem, 0:.025:1, 'FaceColor', px.colors.criticality);
            histogram(ax, s.active_fraction_bst, 0:.025:1, 'FaceColor', px.colors.bistable);
            obj.finish_plot_ax(ax, 'Duration of receptor activity (Fraction of time)', 'Frequency');
            obj.finish_plot_fig(fig);
        end
        
        function figure3c(obj)
            if ~exist('data\active_patches.mat','file')
                % generate distributions if previously-generated files do
                % not exist
                obj.distributions_stream_pulses();
            end
            s = load('data\active_patches.mat');
            fig = figure('Name','Figure 3C');
            ax = axes('Parent', fig, 'Position', [0.18,0.18,0.75,0.75]);
            hold(ax,'on'); box(ax,'on');
            px = plot_xpp();
            histogram(ax, s.active_patches_mono, 1:12, 'FaceColor', px.colors.monostable);
            histogram(ax, s.active_patches_revbst, 1:12, 'FaceColor', px.colors.rev_bistable);
            histogram(ax, s.active_patches_mem, 1:12, 'FaceColor', px.colors.criticality);
            histogram(ax, s.active_patches_bst, 1:12, 'FaceColor', px.colors.bistable);
            xlim([0,13]);
            obj.finish_plot_ax(ax, {'Number of disjoint intervals','of receptor activity'}, 'Frequency');
            obj.finish_plot_fig(fig);
        end
        
        function figure3d(obj)
            fig = figure('Name','Figure 3D');
            ax1 = subplot(3,1,1,'Parent',fig);
            hold(ax1,'on'); box(ax1,'on');
            xlabel(ax1,'Time (min)'); ylabel(ax1,'Response');
            ax2 = subplot(3,1,2,'Parent',fig);
            hold(ax2,'on'); box(ax2,'on');
            xlabel(ax2,'Time (min)'); ylabel(ax2,'Response');
            ax3 = subplot(3,1,3,'Parent',fig);
            hold(ax3,'on'); box(ax3,'on');
            xlabel(ax3,'Time (min)'); ylabel(ax3,'Response');
            px = plot_xpp();
            regimes = fields(px.gamma_dnf);
            model = models.simple_dnf_model;
            % only plot criticality - comment out below if all of them
            % should be plotted
            regimes = regimes(2);
            
            n_pulses = 12;
            lambda = 1.5; % pulses per hour
            t_per_pulse_min = 60/lambda;
            t_pre_stimulus_min = 2;
            t_post_stimulus_min = 2.5;
            t_total_min = n_pulses*t_per_pulse_min +t_pre_stimulus_min +t_post_stimulus_min;
            sppe = experiments.stream_poisson_pulse_experiment(n_pulses, lambda, t_total_min);
            ms = model_simulation(model, sppe);
            for i = 1:numel(regimes)
                % set regime
                gamma_dnf = px.gamma_dnf.(regimes{i});
                model.par.g1 = gamma_dnf;
                ms.set_model(model);
                ms.simulate();
                ms.plot_fraction_phosphorylated('ax', ax2, 'color', px.colors.(regimes{i}), 'minor_format', true, 'plot_pulses', i==1, 'plot_input', i==1);
                drawnow;
                
                mpe = experiments.multi_pulse_experiment(n_pulses);
                mpe.set_pulse_interval_min(40);
                ms = model_simulation(model, mpe);
                ms.simulate();
                ms.plot_fraction_phosphorylated('ax', ax1, 'color', px.colors.(regimes{i}), 'minor_format', true, 'plot_pulses', i==1, 'plot_input', i==1);
                drawnow;
                
                mpe = experiments.multi_pulse_experiment(n_pulses);
                mpe.post_stimulus_min = mpe.post_stimulus_min + (40-22.5)*n_pulses;
                mpe.set_pulse_interval_min(22.5);
                ms = model_simulation(model, mpe);
                ms.simulate();
                ms.plot_fraction_phosphorylated('ax', ax3, 'color', px.colors.(regimes{i}), 'minor_format', true, 'plot_pulses', i==1, 'plot_input', i==1);
                drawnow;
            end
            set(fig,'Position',[250, 250, 1100, 875]);
        end
        
        function figure3e(obj)
            fig = figure('Name','Figure 3E');
            ax = axes('Parent',fig,'Position',[0.15,0.15,0.75,0.75]);
            hold on;
            px = plot_xpp();
            px.plot_response_amplitude(fig);
            obj.finish_plot_ax(ax, 'P_{DNF,T}/R_T', 'Response amplitude');
            obj.finish_plot_fig(fig);
        end
        
        function figure3_double_pulse(obj, model, ax, color, with_pulses_input)
            mpe = experiments.multi_pulse_experiment(2);
            mpe.set_pulse_interval_min(20);
            ms = model_simulation(model, mpe);
            ms.simulate();
            ms.plot_fraction_phosphorylated('ax', ax, 'color', color, 'plot_pulses', with_pulses_input, 'plot_input', with_pulses_input);
        end
        
        function distributions_stream_pulses(obj)
            n_reps = 12; % should be increased to 1000 for better statistics
            active_fraction_mem = zeros(n_reps,1); active_patches_mem = zeros(n_reps,1);
            active_fraction_revbst = zeros(n_reps,1); active_patches_revbst = zeros(n_reps,1);
            active_fraction_mono = zeros(n_reps,1); active_patches_mono = zeros(n_reps,1);
            active_fraction_bst = zeros(n_reps,1); active_patches_bst = zeros(n_reps,1);
            
            n_pulses = 12;
            lambda = 1.5; % pulses per hour
            t_per_pulse_min = 60/lambda;
            t_pre_stimulus_min = 2;
            t_post_stimulus_min = 2.5;
            t_total_min = n_pulses*t_per_pulse_min +t_pre_stimulus_min +t_post_stimulus_min;
            
            px = plot_xpp();
            gamma_dnf = px.gamma_dnf;
            
            try task = getCurrentTask(); catch; task = []; end
            if ~isempty(task); tstr = num2str(task.ID); else; tstr = 'single'; end
            parfor i=1:n_reps
                disp(['lab ',tstr,': ',num2str(i)]);
                mod = models.simple_dnf_model;
                mod.par.g1 = gamma_dnf.('criticality');
                ms = model_simulation(mod,experiments.stream_poisson_pulse_experiment(n_pulses, lambda, t_total_min));
                ms.simulate();
                [active_fraction_mem(i), active_patches_mem(i)] = ms.calc_active_fraction();
                mod.par.g1 = gamma_dnf.('rev_bistable');
                ms.set_model(mod)
                ms.simulate();
                [active_fraction_revbst(i), active_patches_revbst(i)] = ms.calc_active_fraction();
                mod.par.g1 = gamma_dnf.('monostable');
                ms.set_model(mod)
                ms.simulate();
                [active_fraction_mono(i), active_patches_mono(i)] = ms.calc_active_fraction();
                mod.par.g1 = gamma_dnf.('bistable');
                ms.set_model(mod)
                ms.simulate();
                [active_fraction_bst(i), active_patches_bst(i)] = ms.calc_active_fraction();
            end
            save('data\active_fractions.mat','active_fraction_mem','active_fraction_revbst','active_fraction_mono','active_fraction_bst');
            save('data\active_patches.mat','active_patches_mem','active_patches_revbst','active_patches_mono','active_patches_bst');
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