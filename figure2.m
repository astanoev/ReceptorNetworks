classdef figure2
    properties
        fontsize = 10;
        fig_size = [1000, 1000];
    end
    
    methods
        function obj = figure2()
            obj.figure2_main();
        end
        
        function figure2_main(obj)
            % three figures - for pulses simulations, state space and
            % quasy-potential landscape
            % all regimes and conditions are plotted together in each
            % figure
            
            fig_p = figure('Name','Figure 2 - pulses');
            fig_ss = figure('Name','Figure 2 - state spaces');
            %fig_qpl = figure('Name','Figure 2 - quasy-potential landscapes');
            
            px = plot_xpp();
            regimes = fields(px.gamma_dnf);
            % order: monostable, bistable, criticality
            regimes = regimes([4,1,2]);
            model = models.simple_dnf_model;
            
            for i = 1:3
                % set regime
                gamma_dnf = px.gamma_dnf.(regimes{i});
                model.par.g1 = gamma_dnf; 
                
                % plot simulation with pulses
                ax_p(i) = subplot(2,2,i,'Parent',fig_p);
                ax_p(i).ActivePositionProperty = 'position';
                obj.figure2_pulses(model, ax_p(i), regimes{i}, px.colors.(regimes{i}));
                
                % plot state space
                ax_ss(i) = subplot(2,2,i,'Parent',fig_ss);
                obj.figure2_state_space(model, ax_ss(i), 0, regimes{i});
                
                %plot quasy-potential landscape
                %qpl = quasi_potential_landscape(model, 0);
            end
            
            % plot state space with input signal (producing LRa=0.15)
            gamma_dnf = px.gamma_dnf.(regimes{3});
            model.par.g1 = gamma_dnf;
            target_LRa = 0.15;
            input = model.tune_egf_free(target_LRa);
            ax_ss(4) = subplot(2,2,4,'Parent',fig_ss);
            obj.figure2_state_space(model, ax_ss(4), input, 'with input');
            
            % plot quasy-potential landscape with input signal (producing LRa=0.15)
            %qpl = quasi_potential_landscape(model, target_LRa);
            
            obj.finish_plot_fig(fig_p);
            obj.finish_plot_fig(fig_ss);
            %obj.finish_plot_fig(fig_qpl);
        end
        
        function figure2_pulses(obj, model, ax, regime, color)
            % plot single-pulse simulations in main axes
            ms = model_simulation(model, experiments.multi_pulse_experiment(1));
            ms.simulate();
            ms.plot_fraction_phosphorylated('ax', ax, 'color', color);
            % plot double-pulse simulation as inset
            ax_p_pos = ax.Position;
            ax_p2_pos = [ax_p_pos(1:2)+0.57*ax_p_pos(3:4), 0.4*ax_p_pos(3:4)];
            ax_p2 = axes('Parent',ax.Parent,'Position', ax_p2_pos);
            mpe = experiments.multi_pulse_experiment(2);
            mpe.set_pulse_interval_min(30);
            ms = model_simulation(model,mpe);
            ms.simulate();
            ms.plot_fraction_phosphorylated('ax', ax_p2, 'color', color, 'minor_format', true);
            title(ax, regime);
        end
        
        function figure2_state_space(obj, model, ax, input, regime)
            hold(ax,'on');
            model.set_state_space(ax.Parent, [], input);
            model.plot_nullclines(0:.0001:1, input, []);
            model.plot_separatrix(0:.0001:1, input, []);
            model.plot_steady_states(0:.0001:1, input, [], []);
            if strcmp(regime,'criticality') && input==0
                % approximate and plot 'ghost' attractor
                model.par.g1 = model.par.g1-0.01;
                [EGFRpt_ss,PTPRGat_ss,~,stable] = model.get_steady_states(0:.0001:1, 0);
                EGFRpt_gh = mean(EGFRpt_ss(EGFRpt_ss>=EGFRpt_ss(stable==0)));
                PTPRGat_gh = mean(PTPRGat_ss(EGFRpt_ss>=EGFRpt_ss(stable==0)));
                plot(ax, PTPRGat_gh, EGFRpt_gh, 'o','markerfacecolor','none','markeredgecolor',[1, 0.6, 0],'markersize',8);
            end
            obj.finish_plot_ax(ax, 'P_{DNF,a}', 'R_a', [0,0.8], [0,0.8], regime);
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