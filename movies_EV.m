classdef movies_EV
    properties
        fontsize = 10;
        fig_size = [1000, 1000];
    end
    
    methods
        function obj = movies_EV()
            obj.movie_EV1();
%             obj.movie_EV2();
%             obj.movie_EV3();
%             obj.movie_EV4();
%             obj.movie_EV5();
        end
        
        function movie_EV1(obj)
            model = models.simple_dnf_model;
            mpe = experiments.multi_pulse_experiment(1);
            ms = model_simulation(model, mpe);
            ms.stochastic_simulation();
            ms.animate();
        end
        
        function movie_EV2(obj)
            obj.ics_hl = 1;
            obj.g1 = 4.55;
            rep = 7;
            n_parts = 1;
        end
        
        function movie_EV3(obj)
            obj.ics_hl = 1;
            obj.g1 = 4.95;
            rep = 7;
            n_parts = 1;
        end
        
        function movie_EV4(obj)
            obj.ics_hl = 1;
            obj.g1 = 4.95;
            rep = 7;
            n_parts = 1;
            plot_conc_ratio_bkg = true;
        end
        
        function movie_EV5(obj)
            obj.ics_hl = 1;
            obj.g1 = 6;
            rep = 7;
            n_parts = 1;
        end
        
        function finish_plot_ax(obj, ax, xlab, ylab, xlim, ylim, regime)
            set(ax, 'XLim', xlim, 'YLim', ylim);
            xlabel(ax, xlab, 'fontsize', obj.fontsize);
            ylabel(ax, ylab, 'fontsize', obj.fontsize);
            set(ax, 'FontSize', obj.fontsize);
            box on;
            title(ax, regime);
        end
        
        function finish_plot_fig(obj, fig)
            pix_SS = get(0,'screensize');
            fig_pos = get(fig,'Position');
            fig_pos = min([fig_pos(1:2);pix_SS(3:4)-obj.fig_size-100],[],1);
            set(fig, 'Position', [fig_pos(1), fig_pos(2), obj.fig_size(1), obj.fig_size(2)]);
        end
    end
end