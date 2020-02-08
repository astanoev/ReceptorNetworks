classdef figure5
    properties
        fontsize = 20;
        fig_size = [700, 500];
    end
    
    methods
        function obj = figure5()
            obj.figure5b();
            obj.figure5d();
        end
        
        function figure5b(obj)
            fig = figure('Name','Figure 5B');
            ax = axes('Parent',fig,'Position',[0.15,0.15,0.75,0.75]);
            hold(ax,'on'); box(ax,'on');
            px = plot_xpp();
            model = models.simple_dnf_model;
            obj.figure5_pulses(model, ax, px.colors.criticality, false);
            ax.Children(1).LineStyle = '--';
            model = models.model_2comp();
            model.par.dyn_rec = 0;
            obj.figure5_pulses(model, ax, px.colors.monostable, true);
            xlim(ax, [-3, 123]);
            obj.finish_plot_fig(fig);
        end
        
        function figure5d(obj)
            fig = figure('Name','Figure 5D');
            ax = axes('Parent',fig,'Position',[0.15,0.15,0.75,0.75]);
            hold(ax,'on'); box(ax,'on');
            model = models.model_2comp();
            model.par.dyn_rec = 1;
            obj.figure5_pulses(model, ax, [217, 83, 25]./255, true);
            xlim(ax, [-3, 123]);
            obj.finish_plot_fig(fig);
        end
 
        function figure5_pulses(obj, model, ax, color, with_pulses_input) %#ok<INUSL>
            mpe = experiments.multi_pulse_experiment(4);
            mpe.pulse_amplitude = 0.2926;
            mpe.set_pulse_interval_min(30);
            mpe.set_up_input();
            ms = model_simulation(model, mpe);
            ms.simulate();
            ms.plot_fraction_active('ax', ax, 'color', color, 'plot_pulses', with_pulses_input, 'plot_input', with_pulses_input);
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