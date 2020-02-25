classdef agent_based_simulation < handle %matlab.mixin.Copyable
    %AGENT_BASED_SIMULATION

    properties (Constant)
    end
    
    properties
        time_steps = 20001; % 2 seconds
        marker_size = 100;
        square_size = 3.5; % um - width and length
        model = models.agent_based_model();
        R_states
        P_dnf_states
        R_activation_events
        R_positions_init
        P_dnf_positions_init
        R_pos_randsubstream = 2;
        R_pos_randsubstream_state;
        P_dnf_pos_randsubstream = 3;
        P_dnf_pos_randsubstream_state;
        rng_init;
        rng_pos;
        rng_state;
        rng_final;
        part_id = 1;
        n_parts = 1;
        ics_hl = 0;
        gpu;
    end
    
    properties (Transient)
        R_positions
        P_dnf_positions
        geo_net
    end
    
    methods
        function obj = agent_based_simulation(varargin)
            %if ~isempty(dbstack(1)); return; end % do not simulate directly if called from somewhere else
            if nargin==0; return; end
            p = inputParser;
            addOptional(p,'part_id',1);
            addOptional(p,'n_parts',1);
            addOptional(p,'ics_hl',0);
            addOptional(p,'model',models.agent_based_model(obj.time_steps, obj.square_size));
            addParameter(p,'rng_init',rng('shuffle','combRecursive'));%RandStream('mlfg6331_64','Seed','shuffle'));%
            addParameter(p,'gpu',false);
            parse(p,varargin{:});

            obj.R_activation_events = utils.expandable_array(int32(zeros(0,3)));
            obj.part_id = p.Results.part_id;
            obj.n_parts = p.Results.n_parts;
            obj.ics_hl = p.Results.ics_hl;
            obj.model = p.Results.model;
            obj.rng_init = p.Results.rng_init;
            obj.gpu = p.Results.gpu;
        end
        
        function initialize_simulation(obj, varargin)
            p = inputParser;
            addRequired(p,'obj');
            
            if obj.gpu
                obj.R_states = zeros(obj.time_steps,obj.model.Rt,3,'gpuArray');
                obj.P_dnf_states = zeros(obj.time_steps,obj.model.P_dnf_t,2,'gpuArray');
                obj.R_positions = zeros(obj.time_steps,obj.model.Rt,2,'gpuArray');
                obj.P_dnf_positions = zeros(obj.time_steps,obj.model.P_dnf_t,2,'gpuArray');
            else
                obj.R_states = zeros(obj.time_steps,obj.model.Rt,3,'int8');
                obj.P_dnf_states = zeros(obj.time_steps,obj.model.P_dnf_t,2,'int8');
                obj.R_positions = zeros(obj.time_steps,obj.model.Rt,2);
                obj.P_dnf_positions = zeros(obj.time_steps,obj.model.P_dnf_t,2);
            end
            obj.R_activation_events = utils.expandable_array(int32(zeros(0,3)));
                        
            % different schemes for initial conditions
            % most often used were 1 (high Ra) and 4 (low Ra)
            if obj.ics_hl==0
                R_init = [obj.model.Rt,0,0];%[inactive,active,ligand-bound]
                P_dnf_init = [0,obj.model.P_dnf_t];%[i,a]
            elseif obj.ics_hl==1
                R_init = [0,obj.model.Rt,0];%[inactive,active,ligand-bound]
                P_dnf_init = [obj.model.P_dnf_t,0];%[i,a]
            elseif obj.ics_hl==2
                R_init = [int32(0.99*obj.model.Rt),int32(0.01*obj.model.Rt),0];%[inactive,active,ligand-bound]
                P_dnf_init = [0,obj.model.P_dnf_t];%[i,a]
             elseif obj.ics_hl==3
                R_init = [0,obj.model.Rt,0];%[inactive,active,ligand-bound]
                P_dnf_init = [0,obj.model.P_dnf_t];%[i,a]
            elseif obj.ics_hl == 4
                R_init = [obj.model.Rt,0,0];%[inactive,active,ligand-bound]
                P_dnf_i = round(1-obj.model.k1/(obj.model.k1+obj.model.k2));
                P_dnf_init = [P_dnf_i,obj.model.P_dnf_t-P_dnf_i];%[i,a]
            elseif obj.ics_hl == 5 || obj.ics_hl == 6
                n_groups = 5;
                n_per_group = 4;
                n_ligand_bound = n_groups*n_per_group;
                R_init = [obj.model.Rt-n_ligand_bound,0,n_ligand_bound];%[inactive,active,ligand-bound]
                P_dnf_i = round(1-obj.model.k1/(obj.model.k1+obj.model.k2));
                P_dnf_init = [P_dnf_i,obj.model.P_dnf_t-P_dnf_i];%[i,a]
            elseif obj.ics_hl == 7
                n_groups = 10;
                n_per_group = 4;
                n_ligand_bound = n_groups*n_per_group;
                R_init = [obj.model.Rt-n_ligand_bound,0,n_ligand_bound];%[inactive,active,ligand-bound]
                P_dnf_i = round(1-obj.model.k1/(obj.model.k1+obj.model.k2));
                P_dnf_init = [P_dnf_i,obj.model.P_dnf_t-P_dnf_i];%[i,a]
            elseif obj.ics_hl == 8
                n_groups = 5;
                n_per_group = 8;
                n_ligand_bound = n_groups*n_per_group;
                R_init = [obj.model.Rt-n_ligand_bound,0,n_ligand_bound];%[inactive,active,ligand-bound]
                P_dnf_i = round(1-obj.model.k1/(obj.model.k1+obj.model.k2));
                P_dnf_init = [P_dnf_i,obj.model.P_dnf_t-P_dnf_i];%[i,a]
            end
            
            init_ind = 1;
            for i=1:length(R_init)
                obj.R_states(1,init_ind:(init_ind+R_init(i)-1),i) = 1;
                init_ind = init_ind + R_init(i);
            end
            init_ind = 1;
            for i=1:length(P_dnf_init)
                obj.P_dnf_states(1,init_ind:(init_ind+P_dnf_init(i)-1),i) = 1;
                init_ind = init_ind + P_dnf_init(i);
            end
                        
            act = find(squeeze(obj.R_states(1,:,2))==1);
            if ~isempty(act)
                ar = ones(length(act),1)*[1,0,0]; ar(:,3) = act;
                obj.R_activation_events.extend(ar);
            end
            
            obj.R_positions_init = obj.square_size.*rand(1,obj.model.Rt,2);
            obj.P_dnf_positions_init = obj.square_size.*rand(1,obj.model.P_dnf_t,2);
            
            if ismember(obj.ics_hl,[6,7,8])
                inxs = find(obj.R_states(1,:,3)==1);
                pos_groups = obj.square_size.*rand(n_groups,2);
                angle = 2*pi/n_per_group;
                radius = 2*obj.model.interaction_radius;
                for j=1:n_groups
                    obj.R_positions_init(1,inxs(j:n_groups:end),:) = repmat(pos_groups(j,:),n_per_group,1) + radius.*[cos((1:n_per_group)'.*angle),sin((1:n_per_group)'.*angle)];
                end
            end
        end
        
        function set_continuation_variables(obj,abs_mini)
            obj.model = abs_mini.model;
            obj.R_states(1,:,:) = abs_mini.R_states;
            obj.P_dnf_states(1,:,:) = abs_mini.P_dnf_states;
            obj.R_positions_init = abs_mini.R_positions;
            obj.P_dnf_positions_init = abs_mini.P_dnf_positions;
            obj.rng_init = abs_mini.rng_final;
        end
        
        function calculate_positions(obj)
            rng(obj.rng_pos);
            obj.R_positions = calc_pos_mol(obj.R_positions_init, squeeze(obj.R_states(1,:,:)));
            obj.P_dnf_positions = calc_pos_mol(obj.P_dnf_positions_init);
            function pos = calc_pos_mol(pos_init, st_init)
                mol_total = size(pos_init,2);
                std = sqrt(4*obj.model.diffusion*obj.model.delta_t);
                if nargin>1 % if R ligand-bound positions are to be kept constant
                    pos = zeros(obj.model.time_steps, mol_total, 2);
                    pos(:,st_init(:,3)==0,:) = std.*randn(obj.model.time_steps,nnz(st_init(:,3)==0),2);
                else
                    pos = std.*randn(obj.model.time_steps,mol_total,2);
                end
                pos = cumsum(pos,1);
                pos = pos + repmat(pos_init,obj.model.time_steps,1,1);
                pos = mod(pos,obj.model.square_size);
            end
        end
                
        function update_positions(obj, t)
            std = sqrt(4*obj.model.diffusion*obj.model.delta_t); 
            obj.R_positions(t,:,:) = squeeze(obj.R_positions(t-1,:,:)) + std.*randn(obj.model.Rt,2).*repmat((squeeze(obj.R_states(t-1,:,3))==0)',1,2);
            % periodic boundary conditions
            obj.R_positions(t,:,:) = squeeze(obj.R_positions(t,:,:)) + obj.square_size*ones(obj.model.Rt,2);
            obj.R_positions(t,:,:) = squeeze(obj.R_positions(t,:,:)) - obj.square_size*floor(squeeze(obj.R_positions(t,:,:))./obj.square_size);
            
            obj.P_dnf_positions(t,:,:) = squeeze(obj.P_dnf_positions(t-1,:,:)) + std.*randn(obj.model.P_dnf_t,2);
            obj.P_dnf_positions(t,:,:) = squeeze(obj.P_dnf_positions(t,:,:)) + obj.square_size*ones(obj.model.P_dnf_t,2);
            obj.P_dnf_positions(t,:,:) = squeeze(obj.P_dnf_positions(t,:,:)) - obj.square_size*floor(squeeze(obj.P_dnf_positions(t,:,:))./obj.square_size);
        end
                
        function update_states(obj, t)
            Rpos = squeeze(obj.R_positions(t,:,:));
            P_dnf_pos = squeeze(obj.P_dnf_positions(t,:,:));
            Rst = squeeze(obj.R_states(t-1,:,:));
            P_dnf_st = squeeze(obj.P_dnf_states(t-1,:,:));
            % find potential interactions between an (in)active and inactive molecule
            R_R_pot_int = obj.geo_net.potential_interactions(Rpos(Rst(:,1)==1,:),Rpos(:,:));
            % Rst(:,1)==0 includes both Ra and LRa
            % find potential interactions between an Ra and Pdnf_a
            R_P_dnf_pot_int = obj.geo_net.potential_interactions(Rpos(Rst(:,1)==0,:),P_dnf_pos(P_dnf_st(:,1)==0,:));
            Ri = find(Rst(:,1)==1);
            Ra = find(Rst(:,1)==0); 
            P_dnf_a = find(P_dnf_st(:,1)==0);
            
            % copy prev. states and only change upon new event
            obj.R_states(t,:,:) = Rst;
            % pairs of inactive and active Rs
            [inact1,act1] = find(R_R_pot_int);
            inact1 = Ri(inact1);
            diff_links = inact1~=act1; % to remove self-links
            inact1 = inact1(diff_links(:));
            act1 = act1(diff_links(:));
            sample_R_R = rand(numel(inact1),1);
            % set reaction prob. according to rates
            R_act_prob = single(Rst(act1,:))*[obj.model.a1;obj.model.a2;obj.model.a3].*obj.model.delta_t;
            % sample for successful reactions
            R_R_int_pre = R_act_prob>sample_R_R;
            R_R_int = inact1(R_R_int_pre);
            % update activations
            obj.R_states(t,R_R_int,:) = repmat([0,1,0],numel(R_R_int),1);
            
            sample_R_P_dnf = rand(nnz(R_P_dnf_pot_int),2);
            % pairs of active Ras and active Pdnfs
            [act2,act3] = find(R_P_dnf_pot_int);
            act2 = Ra(act2(:)); 
            % set reaction prob
            R_inact_prob = obj.model.g1*obj.model.delta_t;
            % sample
            R_P_dnf_int_pre = R_inact_prob>sample_R_P_dnf(:,1);
            R_P_dnf_int = act2(R_P_dnf_int_pre);
            R_P_dnf_int(Rst(R_P_dnf_int,3)==1) = []; % remove changes to LRa
            % update inactivations to Ra
            obj.R_states(t,R_P_dnf_int,:) = repmat([1,0,0],numel(R_P_dnf_int),1);
            % update R-activation events (for R0)
            obj.record_R_transitions(t,R_R_int,act1(R_R_int_pre),R_P_dnf_int);
            % spontaneous transitions to Pdnf
            obj.P_dnf_states(t,:,:) = P_dnf_st(:,:);
            P_dnf_spont_change = rand(obj.model.P_dnf_t,1);
            P_dnf_spont_changed = P_dnf_spont_change<=single(P_dnf_st(:,:))*[obj.model.k1;obj.model.k2].*obj.model.delta_t;
            obj.P_dnf_states(t,P_dnf_spont_changed,:) = abs(1-P_dnf_st(P_dnf_spont_changed,:));
            % inactivations to Pdnf_a
            P_dnf_inact_prob = obj.model.b1*obj.model.delta_t;
            P_dnf_R_int_pre = P_dnf_inact_prob>sample_R_P_dnf(:,2);
            act3 = P_dnf_a(act3(:));
            P_dnf_R_int = act3(P_dnf_R_int_pre(~P_dnf_spont_changed(act3(P_dnf_R_int_pre))));
            obj.P_dnf_states(t,P_dnf_R_int,:) = repmat([1,0],numel(P_dnf_R_int),1);
        end
        
        function record_R_transitions(obj, t, receivers, donors, minus_ones)
            % Record state transition details (time_of_(in)activation,
            % (-1)donor_id, receiver_id)
            
            [uvals, ~, uidx] = unique(receivers);
            if numel(uvals)==numel(receivers)
                obj.R_activation_events.extend([repmat(t,numel(receivers),1),donors,receivers]);
            else
                cnts = accumarray(uidx, 1);
                act_array = zeros(numel(uvals),3);
                act_array(:,1) = t;
                act_array(cnts==1,2) = donors(cnts==1);
                act_array(:,3) = uvals;
                ct2 = find(cnts>1);
                for i=1:numel(ct2) % choose random donor from the successful interactions
                    act_array(ct2(i),2) = donors(randsample(find(uvals(uidx)==uvals(ct2(i))),1));
                end
                obj.R_activation_events.extend(act_array);
            end
            obj.R_activation_events.extend([repmat([t,-1],numel(minus_ones),1),minus_ones]);
        end
                
        function simulation(obj)
            rng(obj.rng_init);
            % do not initialize when continuing simulation (if part>1)
            if obj.part_id==1; obj.initialize_simulation(); end
                        
            obj.geo_net = utils.geometric_network(obj.model.interaction_radius, obj.square_size);
            
            if isempty(obj.rng_pos); obj.rng_pos = rng; end
            obj.calculate_positions();
            
            if isempty(obj.rng_state); obj.rng_state = rng; end
            for t=2:obj.time_steps
                % no-need anymore since they are generated initially
                %obj.update_positions(t);
                obj.update_states(t);
            end
            obj.rng_final = rng;
        end
    end
end