classdef plot_xpp
    %PLOT_XPP 
    %   For plotting profiles generated from xpp 
    %   .dat files are provided already, but can be generated using the .ode
    %   file (Ra-Pdnfa model) and xppaut
    
    properties
        colors = struct('bistable',[128,179,36]./255,'criticality',[219,23,222]./255,'rev_bistable',[222,163,71]./255,'monostable',[10,168,222]./255);
        gamma_dnf = struct('bistable',2.5,'criticality',2.957,'rev_bistable',3.5,'monostable',4.3);
    end
    
    methods
        function obj = plot_xpp()
            
        end
        
        function plot_bifurcation(obj, fig)
            ax = fig.CurrentAxes;
            folder = 'data';
            filename = 'bifurcation_over_gamma_dnf_with_LRa=0.dat';
            M = dlmread(fullfile(folder,filename),' ');
            x = M(:,1); % gamma_dnf (prop to P_dnft/R_t)
            y = M(:,2); % Ra response
            x_branch = gradient(x)>0;
            idx_unst = find(x_branch==0);
            sn1 = idx_unst(end);
            sn2 = idx_unst(1);
            lw = 3;
            patch(ax,'Faces',1:4,'Vertices',[x(sn1),0;x(sn1),1;x(sn2),1;x(sn2),0],'EdgeColor','none','FaceColor',obj.colors.bistable,'FaceAlpha',.35);
            %patch(ax,'Faces',1:4,'Vertices',[x(sn2),0;x(sn2),1;5,1;5,0],'EdgeColor','none','FaceColor',[obj.colors.criticality;obj.colors.criticality;obj.colors.monostable;obj.colors.monostable],'FaceAlpha',.35);
            vtx = [x(sn2),0;x(sn2),1;x(sn2)+0.1*(5-x(sn2)),1;x(sn2)+0.2*(5-x(sn2)),1;5,1;5,0;x(sn2)+0.2*(5-x(sn2)),0;x(sn2)+0.1*(5-x(sn2)),0];
            cc = obj.colors.criticality;
            cm = obj.colors.monostable;
            col = [cc;cc;cc;cm;cm;cm;cm;cc];
            patch(ax,'Faces',1:8,'Vertices',vtx,'EdgeColor','none','FaceVertexCData',col,'FaceColor','interp','FaceAlpha',.35);
            plot(ax, x(1:sn2), y(1:sn2),'k-','LineWidth',lw);
            plot(ax, x(sn2:sn1), y(sn2:sn1),'k--','LineWidth',lw);
            plot(ax, x(sn1:end), y(sn1:end),'k-','LineWidth',lw);
            text(ax, x(sn1), 2*y(sn1),'SN_1','HorizontalAlignment','right','fontsize',20);
            text(ax, x(sn2), y(sn2),'SN_2','fontsize',20);
        end
        
        function plot_response_amplitude(obj, fig)
            folder = 'data';
            filename1 = 'bifurcation_over_gamma_dnf_with_LRa=0.15.dat';
            filename2 = 'bifurcation_over_gamma_dnf_with_LRa=0.dat';
            M = dlmread(fullfile(folder,filename2),' ');
            xq=0:0.001:5;
            x = M(:,1);
            y = M(:,2);
            [~,inxs] = sort(y,'descend');
            x = x(inxs);
            y = y(inxs);
            y = y((x>=0)&(x<=5));
            x = x((x>=0)&(x<=5));
            x_branch = gradient(x)>0;
            if ~isempty(x(x_branch==0))
                idx = find(x_branch==0);
                y = [y(1:idx(1)-1);y(x>x(idx(1)-1))];
                x = [x(1:idx(1)-1);x(x>x(idx(1)-1))];
                yq1 = interp1(x,y,xq);
            end
            M = dlmread(fullfile(folder,filename1),' ');
            x2 = M(:,1);
            y2 = M(:,2);
            [~,inxs] = sort(y2,'descend');
            x2 = x2(inxs);
            y2 = y2(inxs);
            y2 = y2((x2>=0)&(x2<=5));
            x2 = x2((x2>=0)&(x2<=5));
            yq2 = 0.15+interp1(x2,y2,xq);
            
            ax = fig.CurrentAxes;
            plot(ax, xq, yq2-yq1, 'k-', 'LineWidth',3);
            xlim([0, 5]); ylim([0, 1]);
        end
        
        function plot_two_parameter_bifurcation(obj, fig)
            ax = fig.CurrentAxes;
            folder = 'data';
            filename = 'bifurcation_LRa_vs_gamma_dnf.dat';
            M = dlmread(fullfile(folder,filename),' ');
            x = M(:,1); % gamma_dnf (prop to P_dnft/R_t)
            y = M(:,2); % LRa input
            % plot green bistable patch
            patch(ax,'Faces',1:numel(x),'Vertices',[x,y],'EdgeColor','none','FaceColor',obj.colors.bistable,'FaceAlpha',.33);
            plot(ax, x, y,'-','LineWidth',3,'color',obj.colors.bistable);
            % plot parameter organization points (vertical lines) for regimes
            regimes = fields(obj.gamma_dnf);
            for i = 1:numel(regimes)
                plot(ax, [obj.gamma_dnf.(regimes{i}), obj.gamma_dnf.(regimes{i})], [0, 1],'--','LineWidth',2,'color',obj.colors.(regimes{i}));
            end
            % find and plot cusp point
            [~,idx_cusp] = max(x);
            plot(x(idx_cusp),y(idx_cusp),'k*','MarkerSize',20,'LineWidth',3);
        end
        
        function plot_dose_response_bifurcation(obj, regime, fig)
            ax = fig.CurrentAxes;
            folder = 'data';
            filename = ['bifurcation_over_LRa_with_gamma_dnf=',num2str(obj.gamma_dnf.(regime)),'.dat'];
            M = dlmread(fullfile(folder,filename),' ');
            x = M(:,1); % LRa input
            y = M(:,2); % Ra response
            y = x+y; % system response = Ra+LRa
            % distinguish between stable and unstable parts of the profile
            x_branch = gradient(x)>0;
            lw = 3;
            if ~isempty(x(x_branch==0)) % if there is bistable region
                idx_unst = find((x_branch==0)&(x>=0));
                % find saddle node indices in the unstable branch
                idx_sn2_unst = idx_unst(1);
                idx_sn1_unst = idx_unst(end)+1;
                % find first point in upper branch (relevant for bistable regime)
                idx_st_up_1 = find((x_branch==1)&(y>=y(idx_sn1_unst)&(x>=0)),1,'first')-1;
                % plot stable and unstable branches
                plot(ax, x(1:idx_sn2_unst), y(1:idx_sn2_unst), '-', 'LineWidth',lw, 'Color',obj.colors.(regime));
                plot(ax, x(idx_sn2_unst:idx_sn1_unst), y(idx_sn2_unst:idx_sn1_unst), ':', 'LineWidth',lw, 'Color',obj.colors.(regime));
                plot(ax, x(idx_st_up_1:end), y(idx_st_up_1:end), '-','LineWidth',lw, 'Color',obj.colors.(regime));
                % find saddle node indices in the stable branches:
                % sn2 in upper branch
                idx_sn2_st = find(x(idx_st_up_1:end)>x(idx_sn2_unst),1,'first')+idx_st_up_1-1;
                % sn1 in lower branch
                idx_sn1_st = find(x(1:idx_sn2_unst-1)<x(idx_sn1_unst),2,'last');
                % draw arrows for switch-like activation/deactivation
                ar_on = annotation('arrow');
                set(ar_on,'parent',ax,'position',[x(idx_sn2_unst), y(idx_sn2_unst), 0, y(idx_sn2_st)-y(idx_sn2_unst)],'linestyle','--','linewidth',lw,'headwidth',5*lw,'headlength',5*lw,'color',obj.colors.(regime));
                if ~isempty(idx_sn1_st) % if sn1>=0
                    ar_off = annotation('arrow');
                    idx_sn1_st = idx_sn1_st(1);
                    set(ar_off,'parent',ax,'position',[x(idx_sn1_unst), y(idx_sn1_unst), 0, y(idx_sn1_st)-y(idx_sn1_unst)],'linestyle','--','linewidth',lw,'headwidth',5*lw,'headlength',5*lw,'color',obj.colors.(regime));
                end
            else % if monostable simply plot profile
                plot(ax, x, y,'-','LineWidth',lw,'color',obj.colors.(regime));
            end
        end
    end
end

