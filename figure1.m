classdef figure1
    properties
        fontsize = 20;
        fig_size = [650, 650];
    end
    
    methods
        function obj = figure1()
            obj.figure1b();
            obj.figure1c();
            obj.figure1d();
            obj.figure1e();
        end
        
        function figure1b(obj)
            fig = figure('Name','Figure 1B');
            ax = axes('Parent',fig,'Position',[0.15,0.15,0.75,0.75]);
            hold on;
            px = plot_xpp();
            px.plot_bifurcation(fig);
            obj.finish_plot(fig, ax, 'P_{DNF,T}/R_T', 'R_a', [0,5], [0,1]);
        end
        
        function figure1c(obj)
            fig = figure('Name','Figure 1C');
            ax = axes('Parent',fig,'Position',[0.15,0.15,0.75,0.75]);
            hold on;
            px = plot_xpp();
            px.plot_two_parameter_bifurcation(fig);
            obj.finish_plot(fig, ax, 'P_{DNF,T}/R_T', 'Input', [0,5], [0,0.06]);
        end
        
        function figure1d(obj)
            fig = figure('Name','Figure 1D');
            ax = axes('Parent',fig,'Position',[0.15,0.15,0.75,0.75]);
            hold on;
            px = plot_xpp();
            regimes = fields(px.gamma_dnf);
            for i = numel(regimes):-1:1
                px.plot_dose_response_bifurcation(regimes{i}, fig);
            end
            obj.finish_plot(fig, ax, 'Input', 'Response', [0,0.06], [0,1.0]);
        end
        
        function figure1e(obj)
            fig = figure('Name','Figure 1E');
            ax = axes('Parent',fig);
            hold on;
            % create model to simulate
            model = models.simple_dnf_model;
            px = plot_xpp();
            regimes = fields(px.gamma_dnf);
            for i = numel(regimes):-1:1
                gamma_dnf = px.gamma_dnf.(regimes{i});
                model.par.g1 = gamma_dnf; % set regime
                % define simulation with multistep experiment type
                ms = model_simulation(model,experiments.multistep_experiment(5));
                % integrate ode
                ms.simulate();
                % convert time to minutes
                time_min = ms.experiment.correct_time(ms.time);
                % get profiles for LRa and Ra
                LRa = ms.state(strcmp(ms.model.labels,'LRa'),:);
                Ra = ms.state(strcmp(ms.model.labels,'Ra'),:);
                if i==numel(regimes)
                    % plot the LRa input with the first plot
                    plot(ax, time_min, LRa, 'LineWidth', 3, 'Color', [0.75,0.75,0.75]);
                end
                % plot response
                plot(ax, time_min, Ra+LRa, 'LineWidth', 3, 'Color', px.colors.(regimes{i}));
            end
            obj.finish_plot(fig, ax, 'Time (min)', 'Response', [time_min(1),time_min(end)], [0,1.0]);
        end
        
        function finish_plot(obj, fig, ax, xlab, ylab, xlim, ylim)
            set(ax, 'XLim', xlim, 'YLim', ylim);
            xlabel(ax, xlab, 'fontsize', obj.fontsize);
            ylabel(ax, ylab, 'fontsize', obj.fontsize);
            set(ax, 'FontSize', obj.fontsize);
            box on;
            pix_SS = get(0,'screensize');
            fig_pos = get(fig,'Position');
            fig_pos = min([fig_pos(1:2);pix_SS(3:4)-obj.fig_size-100],[],1);
            set(fig, 'Position', [fig_pos(1), fig_pos(2), obj.fig_size(1), obj.fig_size(2)]);
            drawnow;
        end
    end
end