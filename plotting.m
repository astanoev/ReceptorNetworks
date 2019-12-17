classdef plotting
    properties
        savefig = 0;
    end
    
    methods
        function obj = plotting(savefig)
            if nargin>0; obj.savefig = savefig; end
        end
    end
    
    methods
        function plot_fraction_phosphorylated(obj, fig, ax, t, minor_format, ms_obj)
            xlim(ax, [min(t),max(t)+1e-5]);
            ylim(ax, [0,1]);
            box(ax,'on');
            set(ax,'fontname','Arial');
            if ~minor_format
                xlabel(ax, 'Time (min)');
                ylabel(ax, 'Response');
                set(ax,'XTick',0:15:max(t),'YTick',0:.2:1);
                set(ax,'fontsize',20);
            end
            if obj.savefig==1
                s = load('mat/centered.mat');
                setprinttemplate(fig,s.template);
                set(fig,'Position',[150   150   650   600]);
                filename = sprintf('figure_g1=%.5f',ms_obj.model.par.g1);
                obj.save_figure(filename, fig); 
            end
        end
        
        function plot_state_space(obj, fig, model, EGFfree, dAB, q, x_lim, y_lim)
            ax = fig.CurrentAxes;
            mags = log(dAB);
            currentColormap = colormap(hsv);
            currentColormap = flipud(currentColormap);
            [~, ~, ind] = histcounts(mags, size(currentColormap, 1));
            ind = max(ind,20); % substitute circular red-violet colors in colormap
            cmap = uint8(ind2rgb(ind(:), currentColormap) * 255);
            cmap(:,:,4) = 255;
            cmap = permute(repmat(cmap, [1 3 1]), [2 1 3]);
            set(q, 'AutoScaleFactor', 0.75, 'MaxHeadSize',5,'LineWidth',1.5);
            set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(cmap(1:3,:,:), [], 4).');   %'
            set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(cmap(1:2,:,:), [], 4).');
            
            axis equal;
            set(ax,'XTick',0:.2:1,'YTick',0:.2:1);
            xlim(ax,[0,y_lim]); ylim(ax,[0,x_lim]);
            set(ax, 'FontSize', 20, 'fontname', 'Arial');
            %xlabel('PTPRGa'); ylabel('EGFRp');
            %xlabel('P_{DNF,a}'); ylabel('R_a');
            xlabel(ax,'PTP activity', 'FontSize', 25); ylabel(ax,'EGFR phosphorylation', 'FontSize', 25);
            box on;
            
            %s = load('mat/png.mat');
            %setprinttemplate(fig,s.template);
            %set(fig, 'Position', [488, 162, 650, 577]);
            
            filename = sprintf('figure_ss_g1=%.5f_egf=%.3f', model.par.g1, EGFfree);
            if obj.savefig==1; obj.save_figure(filename, fig); end
        end
        
        function save_figure(obj, filename, fig)
            dt = datetime('now','TimeZone','America/New_York');
            date = datestr(dt,'ddmmyyyy');
            folder = fullfile('C:/Users/',getenv('USERNAME'),'/Dropbox/Projects/PTP-2 project/Paper/Figures/Fig3/',date,'/');
            if ~exist(folder, 'dir'); mkdir(folder); end
            saveas(fig,fullfile(folder,strcat(filename,'.svg')),'svg');
            saveas(fig,fullfile(folder,strcat(filename,'.fig')),'fig');
        end
    end
end