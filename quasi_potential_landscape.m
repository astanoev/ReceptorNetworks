classdef quasi_potential_landscape
    %QUASI_POTENTIAL_LANDSCAPE takes model as an input and plots the
    %quasi-potential landscape
    
    properties
        model
    end
    
    properties (Access=private)
        xy_grid_spacing = 0.02;
        x_array = 0:0.02:1; 
        y_array = 0:0.02:1; 
        x_index;
        y_index;
        x_path;
        y_path;
        v_path;
        input = 0;
        time_sec = 150;
        time_sec_max = 5000;
        n_time_step_samples = 250;
        save_qpl_paths = true;
        zmax = 0.01;
    end
    
    methods
        function obj = quasi_potential_landscape(model, input)
            if nargin<1
                obj.model = models.simple_dnf_model;
            else
                obj.model = model;
            end
            if nargin<2
                obj.input = 0;
            else
                obj.input = input;
            end
            obj.x_index = find(strcmp(obj.model.labels,'Pdnf_a'));
            obj.y_index = find(strcmp(obj.model.labels,'Ra'));
            [obj.x_path, obj.y_path, obj.v_path] = obj.calc_landscape();
            
            if ~isempty(dbstack(1)); return; end
            obj.plot_landscape();
        end
        
        function [x_path,y_path,v_path] = calc_landscape(obj)
            % load phase-space trajectory data for landscape visualization
            str = obj.get_parameters_string(obj.input);
            if exist(fullfile('data','quasi_potential_landscape',strcat(str,'.mat')),'file')
                s = load(fullfile('data','quasi_potential_landscape',strcat(str,'.mat')));
                [x_path,y_path,v_path] = deal(s.x_path,s.y_path,s.v_path);
                return;
            end
            
            [x_path, y_path, v_path] = obj.calc_paths();
            
            [v_path, attractors_id_x_y_v, path_attractor_ids] = obj.assign_level_attractors(x_path, y_path, v_path);
            
            [x_path, y_path, v_path, attractor_graph, sepx_matrix_attr1_attr2_node1_node2, path_attractor_ids] = obj.match_basins_at_separatrix(x_path, y_path, v_path, attractors_id_x_y_v, path_attractor_ids);
            
            [v_path, ~] = obj.stitch_basins(v_path, attractors_id_x_y_v, path_attractor_ids, attractor_graph, sepx_matrix_attr1_attr2_node1_node2);

            if obj.save_qpl_paths; save(fullfile('data','quasi_potential_landscape',strcat(str,'.mat')),'x_path','y_path','v_path'); end
        end
    
        function [x_path, y_path, v_path] = calc_paths(obj)
            %% calculate phase-space trajectories and potentials starting 
            % from points on the grid as initial conditions
            sce = experiments.simple_convergence_experiment();
            sce.set_up_input(obj.input);
            
            % initialize "path" variable matrices
            x_path = nan(length(obj.x_array),length(obj.y_array),obj.n_time_step_samples);  %x-coord. of phase-space traj.
            y_path = nan(length(obj.x_array),length(obj.y_array),obj.n_time_step_samples);  %y-coord. of phase-space traj.
            v_path = nan(length(obj.x_array),length(obj.y_array),obj.n_time_step_samples);  %calculated potential along traj.
            % to skip areas in phase-space that are calculated by prev calcs
            grid_visited = zeros(length(obj.x_array),length(obj.y_array));
            % spiral outwards to inwards
            sp = spiral(length(obj.x_array));
            [~,idx] = sort(sp(:),'descend');
            [R,C] = ind2sub([length(obj.x_array),length(obj.x_array)],idx);

            % run trajectories from seeds from the grid and save relative potential values
            % Loop over x-y grid along the spiral
            fprintf('0.00\n');
            perc_rc = 0;
            for rc =1:length(R)
                i = R(rc); j = C(rc);
                x0 = obj.x_array(i); y0 = obj.y_array(j);
                % do not skip border grid points, even if 'visited'
                if grid_visited(i,j)==1 && rc>2*(length(obj.x_array)+length(obj.y_array))-2; continue; end
                % update timer when percentage changes
                if floor(100*rc/length(R))>perc_rc
                    perc_rc = floor(100*rc/length(R));
                    for i_erase=1:5; fprintf('\b'); end
                    fprintf('%.2f\n',rc/length(R));
                end
                
                [x_path(i,j,:), y_path(i,j,:), v_path(i,j,:), sol] = obj.calc_prune_traj(x0, y0, sce);
                
                % exclude phase-space points as initial conditions where
                % any trajectory has already passed (more comp. and data efficient)
                around_point = obj.xy_grid_spacing/2;
                bin = floor(([sol.y(obj.x_index,:);sol.y(obj.y_index,:)]'+around_point)/(2*around_point))+1;
                bin = unique(bin,'rows');
                for k = 1:size(bin,1)
                    if grid_visited(bin(k,1),bin(k,2))==1; continue; end
                    grid_visited(bin(k,1),bin(k,2)) = 1;
                    % set end-point v_path, important for later
                    x_path(bin(k,1),bin(k,2),end) = x_path(i,j,end);
                    y_path(bin(k,1),bin(k,2),end) = y_path(i,j,end);
                    v_path(bin(k,1),bin(k,2),end) = v_path(i,j,end);
                end
            end
            for i_erase=1:5; fprintf('\b'); end % erase timer
        end
        
        function [x_path, y_path, v_path, sol] = calc_prune_traj(obj, x0, y0, sce)
            %% calculate trajectory from initial conditions (x0,y0)
            % and prune the data points to only keep the most 'diverse' of them
            % sce = simple_convergence_experiment
            
            init_conds = [y0; obj.model.quasi_steady_state(y0, obj.input); 0]; % last element is zero for quasi_potential
            init_conds([obj.x_index,obj.y_index]) = [x0,y0];
            % calculate with quasy-steady state assumptions for all vars except
            % x and y vars
            qss = ones(length(init_conds)-1,1);
            qss([obj.x_index,obj.y_index]) = [0,0];
            qss(obj.y_index) = [];
            
            t_final = obj.time_sec;
            sol = ode45(@(t,z) obj.deriv(t,z,sce,qss), [0,obj.time_sec], init_conds, odeset('RelTol',1e-6,'AbsTol',1e-8));
            sol_der = obj.deriv(sol.x(:,end),sol.y(:,end),sce,qss);
            convergence = norm(sol_der);
            % continue the solver if not converged to steady state
            while convergence>1e-8 && t_final < obj.time_sec_max
                t_final = t_final + obj.time_sec;
                sol = odextend(sol, [], t_final);
                sol_der = obj.deriv(sol.x(:,end),sol.y(:,end),sce,qss);
                convergence = norm(sol_der);
            end

            % down-size iteratively if too many sample points along trajectory
            IA = 1:size(sol.y,2);
            tol = 0.0001;
            while length(IA)>obj.n_time_step_samples 
                [~,IA] = uniquetol(cumsum(sqrt(sum((sol.y(:,2:end)-sol.y(:,1:end-1)).^2,1))),tol,'OutputAllIndices',false);
                tol = tol*2;
            end
            Z = nan(size(sol.y,1),obj.n_time_step_samples);
            Z(:,1:length(IA)) = sol.y(:,IA);
            Z(:,end) = sol.y(:,end);

            x_path = Z(obj.x_index,:);
            y_path = Z(obj.y_index,:);
            v_path = Z(end,:);
        end
        
        function [v_path, attractors_id_x_y_v, path_attractor_ids] = assign_level_attractors(obj, x_path, y_path, v_path)
            %% assign trajectories to a converging attractor point and update potential
            % values of each trajectory so that they are all level at the
            % steady-state point within the basin
            attractors_id_x_y_v = []; % assign array to keep track of attractors and their coordinates; and pot.
            path_attractor_ids = ones(length(obj.x_array),length(obj.y_array),1);   %attractor id for given path (to denote basin of attraction)
            for i = 1:length(obj.x_array)
                for j = 1:length(obj.y_array)
                    x_ss = x_path(i,j,end);
                    y_ss = y_path(i,j,end);
                    v_ss = v_path(i,j,end);
                    attractor_id = obj.get_attractor_id(x_ss, y_ss, attractors_id_x_y_v);
                    if attractor_id == -1
                        attractor_id = size(attractors_id_x_y_v,1)+1;
                        current_att_num_x_y_v = [attractor_id, x_ss, y_ss, v_ss];
                        attractors_id_x_y_v = [attractors_id_x_y_v; current_att_num_x_y_v]; %#ok<AGROW>
                    else % update path to match the attractor potential
                        v_path(i,j, :) = v_path(i,j, :) + (attractors_id_x_y_v(attractor_id, 4) - v_ss);
                    end
                    path_attractor_ids(i,j) = attractor_id;
                end
            end
        end
        
        function [x_path, y_path, v_path, attractor_graph, sepx_matrix_attr1_attr2_node1_node2, path_attractor_ids] = match_basins_at_separatrix(obj, x_path, y_path, v_path, attractors_id_x_y_v, path_attractor_ids)
            %% find neighboring points in the grid whose trajectories converge to
            % different attractors - separatrix pairs - used for matching relative
            % attractor potentials
            sce = experiments.simple_convergence_experiment();
            sce.set_up_input(obj.input);
            qss = ones(length(obj.model.init_conds),1);
            qss([obj.x_index,obj.y_index]) = [0,0];
            qss(obj.y_index) = [];
            % neighborhood graph between basins of attractors
            attractor_graph = zeros(size(attractors_id_x_y_v,1),size(attractors_id_x_y_v,1));
            % assign array to keep track of pairs of points accross separatrices
            sepx_matrix_attr1_attr2_node1_node2 = []; 
            i = 1;
            while i<=length(obj.x_array)
                j = 1;
                while j<=length(obj.y_array)
                    if i>1 && path_attractor_ids(i,j) ~= path_attractor_ids(i-1,j)
                        [x_path, y_path, v_path, path_attractor_ids, z1, ~] = obj.check_separatrix_point(i, j, sce, x_path, y_path, v_path, attractors_id_x_y_v, path_attractor_ids, qss);
                        [x_path, y_path, v_path, path_attractor_ids, z2, point_misassigned] = obj.check_separatrix_point(i-1, j, sce, x_path, y_path, v_path, attractors_id_x_y_v, path_attractor_ids, qss);
                        % return to previous i-point if misassigned to wrong attractor
                        if point_misassigned; i = i-1; continue; end
                        % connect attractors as neighboring in the graph
                        attractor_graph(path_attractor_ids(i-1,j), path_attractor_ids(i,j)) = 1;
                        attractor_graph(path_attractor_ids(i,j), path_attractor_ids(i-1,j)) = 1;
                        % calculate cosine of angle between diverging trajectories
                        cos_theta = dot(z1(1:2),z2(1:2))/(norm(z1(1:2))*norm(z2(1:2)));
                        curr_sepx = [path_attractor_ids(i-1,j), path_attractor_ids(i,j), i-1, j, i, j, cos_theta];
                        curr_sepx_reversed = [path_attractor_ids(i,j), path_attractor_ids(i-1,j), i, j, i-1, j, cos_theta];
                        % record separatrix pair in matrix
                        sepx_matrix_attr1_attr2_node1_node2 = [sepx_matrix_attr1_attr2_node1_node2; curr_sepx; curr_sepx_reversed]; %#ok<AGROW>
                    end
                    if j>1 && path_attractor_ids(i,j) ~= path_attractor_ids(i,j-1)
                        [x_path, y_path, v_path, path_attractor_ids, z1, ~] = obj.check_separatrix_point(i, j, sce, x_path, y_path, v_path, attractors_id_x_y_v, path_attractor_ids, qss);
                        [x_path, y_path, v_path, path_attractor_ids, z2, point_misassigned] = obj.check_separatrix_point(i, j-1, sce, x_path, y_path, v_path, attractors_id_x_y_v, path_attractor_ids, qss);
                        % return to previous i-point if misassigned to wrong attractor
                        if point_misassigned; j = j-1; continue; end
                        
                        attractor_graph(path_attractor_ids(i,j-1), path_attractor_ids(i,j)) = 1;
                        attractor_graph(path_attractor_ids(i,j), path_attractor_ids(i,j-1)) = 1;
                        cos_theta = dot(z1(1:2),z2(1:2))/(norm(z1(1:2))*norm(z2(1:2)));
                        curr_sepx = [path_attractor_ids(i,j-1), path_attractor_ids(i,j), i, j-1, i, j, cos_theta];
                        curr_sepx_reversed = [path_attractor_ids(i,j), path_attractor_ids(i,j-1), i, j, i, j-1, cos_theta];
                        sepx_matrix_attr1_attr2_node1_node2 = [sepx_matrix_attr1_attr2_node1_node2; curr_sepx; curr_sepx_reversed]; %#ok<AGROW>
                    end
                    j = j+1;
                end
                i = i+1;
            end
        end
        
        function [x_path, y_path, v_path, path_attractor_ids, z1, point_misassigned] = check_separatrix_point(obj, i, j, sce, x_path, y_path, v_path, attractors_id_x_y_v, path_attractor_ids, qss)
            %% check if separatrix point has not been calculated or misassigned
            % update attractor assignment and v-path accordingly
            point_misassigned = false;
            x = obj.x_array(i); y = obj.y_array(j);
            if isnan(x_path(i,j,1))
                [x_path(i,j,:), y_path(i,j,:), v_path(i,j,:), ~] = obj.calc_prune_traj(x, y, sce);
                if obj.get_attractor_id(x_path(i,j,end), y_path(i,j,end), attractors_id_x_y_v) ~= path_attractor_ids(i,j)
                    path_attractor_ids(i,j) = obj.get_attractor_id(x_path(i,j,end), y_path(i,j,end), attractors_id_x_y_v);
                    point_misassigned = true;
                end
                v_path(i,j,:) = v_path(i,j,:) + (attractors_id_x_y_v(path_attractor_ids(i,j), 4) - v_path(i,j,end));
            end
            init_conds = [y;obj.model.quasi_steady_state(y, obj.input);0];
            init_conds([obj.x_index,obj.y_index]) = [x,y];
            z1 = obj.deriv(0,init_conds,sce,qss);
        end
        
        function [v_path, attractors_id_x_y_v] = stitch_basins(obj, v_path, attractors_id_x_y_v, path_attractor_ids, attractor_graph, sepx_matrix_attr1_attr2_node1_node2) %#ok<INUSL>
            %% 'stitch' the neighboring attractors using the separatrix pairs
            % use queue to process attractor graph - breadth-first search
            queued = zeros(size(attractors_id_x_y_v,1),1);
            processed = zeros(size(attractors_id_x_y_v,1),1);
            queued(1) = 1;
            queue = (1);
            while(~isempty(queue))
                attr_id = queue(1);
                ref_attrs = find((attractor_graph(:,attr_id)==1).*(processed==1));
                if ~isempty(ref_attrs)
                    % find separatrix pairs between ref_attr and attr
                    idxs = find((sepx_matrix_attr1_attr2_node1_node2(:,1)==attr_id).*ismember(sepx_matrix_attr1_attr2_node1_node2(:,2),ref_attrs));
                    v_path_ref_attr_basin = zeros(length(idxs),1);
                    v_path_attr_basin = zeros(length(idxs),1);
                    weights = zeros(length(idxs),1);
                    for i=1:length(idxs)
                        v_path_ref_attr_basin(i) = v_path(sepx_matrix_attr1_attr2_node1_node2(idxs(i),5),sepx_matrix_attr1_attr2_node1_node2(idxs(i),6));
                        v_path_attr_basin(i) = v_path(sepx_matrix_attr1_attr2_node1_node2(idxs(i),3),sepx_matrix_attr1_attr2_node1_node2(idxs(i),4));
                        weights(i) = .5*(1-sepx_matrix_attr1_attr2_node1_node2(idxs(i),7));
                    end
                    % give more weight to the pairs of points in state
                    % space that point away from each other
                    % (theta->180deg)
                    k = sum(weights.*(v_path_ref_attr_basin-v_path_attr_basin))/sum(weights);
                    attr_path_idxs = find(path_attractor_ids==attr_id);
                    [ii,jj] = ind2sub(size(path_attractor_ids),attr_path_idxs);
                    for i=1:length(ii)
                        v_path(ii(i),jj(i),:) = v_path(ii(i),jj(i),:) + k;
                    end
                    attractors_id_x_y_v(attr_id,4) = attractors_id_x_y_v(attr_id,4) + k;
                end
                processed(attr_id) = 1;
                queue(1) = [];
                next_attrs = find((attractor_graph(:,attr_id)==1).*(queued==0));
                if ~isempty(next_attrs)
                    for i=1:length(next_attrs)
                        queue = [queue;next_attrs(i)]; %#ok<AGROW>
                        queued(next_attrs(i)) = 1;
                    end
                end
            end
            
            v_path = v_path - min(min(min(v_path)));
        end
        
        function plot_landscape(obj, fig, ax, zmax)
            if nargin<4; zmax = obj.zmax; end
            xxx = obj.x_path(:);
            yyy = obj.y_path(:);
            vvv = obj.v_path(:);
            xxx = xxx(~isnan(vvv));
            yyy = yyy(~isnan(vvv));
            vvv = vvv(~isnan(vvv));
            % remove duplicates, so receive warning from griddata only when v-values are not matched properly
            [~,ids] = uniquetol([xxx,yyy,vvv],1e-6,'ByRows',true,'OutputAllIndices',false);
            
            [Xq,Yq] = meshgrid([0:0.01:0.1,0.15:0.05:1.0]);
            Vq = griddata(xxx(ids),yyy(ids),vvv(ids),Xq,Yq);
                
            %% -- Plot 3d potential landscape --
            if nargin <2; fig = figure(); end
            if nargin <3; ax = axes('Parent',fig); end
            hold(ax,'on');
            
            mesh(ax, Xq, Yq, Vq, 'FaceColor','flat', 'EdgeColor',0.55*ones(1,3), 'FaceAlpha', 0.75, 'LineWidth', 0.1);
            plot3(ax, Xq(1,:), Yq(1,:), Vq(1,:),'Color',0.5*ones(1,3),'linewidth',1);
            plot3(ax, Xq(:,1), Yq(:,1), Vq(:,1),'Color',0.5*ones(1,3),'linewidth',1);
            plot3(ax, Xq(:,end), Yq(:,end), Vq(:,end),'Color',0.5*ones(1,3),'linewidth',1);
            plot3(ax, Xq(end,:), Yq(end,:), Vq(end,:),'Color',0.5*ones(1,3),'linewidth',1);
            colormap(ax, pink);
            caxis(ax, [0 zmax]);
            grid(ax, 'on');
            box(ax, 'on');
            view(ax, -45, 20);
            xlabel(ax, 'P_{DNF,a}', 'fontsize',20);
            ylabel(ax, 'R_a', 'fontsize',20);
            zlabel(ax, 'Quasi-Potential', 'fontsize',20);
            xlim(ax, [0,1]);
            ylim(ax, [0,1]);
            zlim(ax, [0,zmax]);
            set(ax,'fontname','Arial');
            set(ax, 'FontSize', 20);
            set(fig,'renderer','painters');
            %set(fig,'Position',[115,55,1132,676]);

            return;
        end
            
        function str = get_parameters_string(obj, Lt)
            str = 'quasi_potential_landscape';
            cl_model = class(obj.model); % get class handle
            model_name = strsplit(cl_model,'.');
            model_name = model_name{2};
            str = strcat(str,'_',model_name);
            cl_constr = str2func(cl_model); % get constructor handle
            model_original = cl_constr(); % create new object from constructor
            par_fields = fields(obj.model.par);
            for i=1:length(par_fields)
                par_name = par_fields{i};
                par_val = obj.model.par.(par_name);
                par_val_orig = model_original.par.(par_name);
                if par_val ~= par_val_orig
                    str = strcat(str,'_',par_name,'_',num2str(par_val));
                end
            end
            if nargin>1
                if Lt>0
                    strEGF = strcat('Lt_',num2str(Lt,3));
                    str = strcat(str,'_',strEGF);
                end
            end
        end

        function attractor_id = get_attractor_id(obj, x_curr, y_curr, attractors_id_x_y_v) %#ok<INUSL>
            tol = 1e-3;
            for attractor_id = 1 : size(attractors_id_x_y_v,1)
                x_att = attractors_id_x_y_v(attractor_id,2);
                y_att = attractors_id_x_y_v(attractor_id,3);
                if norm([x_att-x_curr,y_att-y_curr]) < tol % this attractor has been identified before
                    return;
                end
            end
            attractor_id = -1;
        end

        function dydt = deriv(obj, t, y, exp, qss)
            dydt = obj.model.df_model(t, y(1:end-1), exp, qss);
            quasi_potential = -dydt(obj.x_index).*dydt(obj.x_index) -dydt(obj.y_index).*dydt(obj.y_index);
            dydt = [dydt; quasi_potential];
        end
    end
end

