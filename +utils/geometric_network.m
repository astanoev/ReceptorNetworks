classdef geometric_network < handle
    %build a geometric network from points with uniform (x,y) distribution
    %on a 2-d plane, by connecting pairs within a given distance
        
    properties
        neigh
        interaction_radius
        square_size
    end
    
    methods
        function obj = geometric_network(interaction_radius, square_size)
            obj.interaction_radius = interaction_radius;
            obj.square_size = square_size;
            obj.neigh = obj.neighborhood_mat();
        end
        
        function pot_int = potential_interactions(obj, pos1, pos2)
            %% hybrid algorithm for finding pairs of points within a given distance
            % uses three different algorithms that operate optimally in
            % different conditions, then chooses the best one given the
            % size of pos1 and pos2 (based on empirical data):
            % pdist is good for smaller sets, it doesn't depend on density
            % kdtree performs good all around and for different densities
            % bin scales really really good, so it's best for large sets
            % and large radiuses
            if size(pos1,1)==0 || size(pos2,1)==0; pot_int = sparse(size(pos1,1),size(pos2,1)); return; end
            if nargin<=2
                if size(pos1,1)<=300
                    pot_int = obj.pdist_method(pos1);
                elseif obj.interaction_radius/obj.square_size<=0.004
                    pot_int = obj.kdtree_method(pos1);
                else
                    pot_int = obj.bin_method(pos1);
                end
            else
                if (size(pos1,1)<=400 && size(pos2,1)<200) || (size(pos1,1)<200 && size(pos2,1)<=400)
                    pot_int = obj.pdist_method(pos1, pos2);
                elseif obj.interaction_radius/obj.square_size<=0.004
                    pot_int = obj.kdtree_method(pos1, pos2);
                else
                    pot_int = obj.bin_method(pos1, pos2);
                end
            end
        end
        
        function pot_int = pdist_method(obj, pos1, pos2)
            %% pdist is good for smaller sets, it doesn't depend on density
            if nargin<=2 % intra-distance
                dists_x = pdist(pos1(:,1));
                % have to do this way because of the periodic boundary conditions
                dists_x = min(dists_x,obj.square_size-dists_x);
                dists_y = pdist(pos1(:,2));
                dists_y = min(dists_y,obj.square_size-dists_y);
                dists = sqrt(dists_x.^2+dists_y.^2);
                pot_int = squareform_sp(dists<=obj.interaction_radius);
            else % inter-distance
                dists_x = pdist2(pos1(:,1),pos2(:,1));
                dists_x = min(dists_x,obj.square_size-dists_x);
                dists_y = pdist2(pos1(:,2),pos2(:,2));
                dists_y = min(dists_y,obj.square_size-dists_y);
                dists = sqrt(dists_x.^2+dists_y.^2);
                pot_int = sparse(dists<=obj.interaction_radius);
            end
        end
        
        function pot_int = kdtree_method(obj, pos1, pos2)
            %% kdtree performs good all around and for different densities
            intra = nargin<=2;
            if intra
                pos2 = pos1;
            end
            % periodic boundary conditions
            % add the points near the border to the other side and diag.
            inxs1 = find(pos2(:,1)<=obj.interaction_radius);
            inxs2 = find(pos2(:,1)>=obj.square_size-obj.interaction_radius);
            inxs3 = find(pos2(:,2)<=obj.interaction_radius);
            inxs4 = find(pos2(:,2)>=obj.square_size-obj.interaction_radius);
            inxs5 = find(pos2(:,1)<=obj.interaction_radius & pos2(:,2)<=obj.interaction_radius);
            inxs6 = find(pos2(:,1)<=obj.interaction_radius & pos2(:,2)>=obj.square_size-obj.interaction_radius);
            inxs7 = find(pos2(:,1)>=obj.square_size-obj.interaction_radius & pos2(:,2)<=obj.interaction_radius);
            inxs8 = find(pos2(:,1)>=obj.square_size-obj.interaction_radius & pos2(:,2)>=obj.square_size-obj.interaction_radius);
            inxs = [(1:size(pos2,1))';inxs1;inxs2;inxs3;inxs4;inxs5;inxs6;inxs7;inxs8];
            pos22 = [pos2;[pos2(inxs1,1)+obj.square_size,pos2(inxs1,2)];[pos2(inxs2,1)-obj.square_size,pos2(inxs2,2)];...
                [pos2(inxs3,1),pos2(inxs3,2)+obj.square_size];[pos2(inxs4,1),pos2(inxs4,2)-obj.square_size];...
                [pos2(inxs5,1)+obj.square_size,pos2(inxs5,2)+obj.square_size];[pos2(inxs6,1)+obj.square_size,pos2(inxs6,2)-obj.square_size];...
                [pos2(inxs7,1)-obj.square_size,pos2(inxs7,2)+obj.square_size];[pos2(inxs8,1)-obj.square_size,pos2(inxs8,2)-obj.square_size]];
            
            Idx = rangesearch(pos22,pos1,obj.interaction_radius,'BucketSize',200);
            Idx_ct = cellfun('length',Idx);
            links_inxs = find(Idx_ct>intra); % points with neighbors
            if isempty(links_inxs)
                pot_int = sparse(size(pos1,1),size(pos2,1));
                return;
            end
            l = repelem(links_inxs,Idx_ct(links_inxs)-intra);
            r = inxs(horzcat(Idx{links_inxs}));
            if intra && ~isempty(l) % remove self-links
                erase = [1;1+cumsum(Idx_ct(links_inxs(1:end-1)))];
                r(erase) = [];
            end
            pot_int = sparse(l,r,1,size(pos1,1),size(pos2,1));
        end
        
        function pot_int = bin_method(obj, pos1, pos2)
            %% bin scales really good, so it's best for large sets
            n_bins_dim = floor(obj.square_size/obj.interaction_radius); % number of bins
            bins1 = ceil(pos1./obj.interaction_radius); % allocate points to bins
            bins1(bins1>n_bins_dim) = n_bins_dim; % if uneven bins, merge last two
            inxs1 = sub2ind([n_bins_dim,n_bins_dim],bins1(:,1),bins1(:,2)); 
            % N*n_el matrix - points to bins
            A = sparse(1:size(pos1,1),inxs1,1,size(pos1,1),n_bins_dim*n_bins_dim);
            if nargin<=2
                pos2 = pos1;
                B = A;
            else 
                % allocation for second set of points
                bins2 = ceil(pos2./obj.interaction_radius);
                bins2(bins2>n_bins_dim) = n_bins_dim;
                inxs2 = sub2ind([n_bins_dim,n_bins_dim],bins2(:,1),bins2(:,2));
                B = sparse(1:size(pos2,1),inxs2,1,size(pos2,1),n_bins_dim*n_bins_dim);
            end
            % neigh matrix contains connections between neighboring
            % bins - 9 in total inc. self link
            % A*neigh*B' will connect points within neighboring areas
            pot_con = A*obj.neigh*B';
            if nargin<=2; pot_con = pot_con-speye(size(A,1)); end
            [i,j] = find(pot_con);
            dists = abs(pos1(i,:)-pos2(j,:));
            dists = sqrt(sum(min(dists,obj.square_size-dists).^2,2));
            % filter the remaining by distance
            inxs = find(dists<=obj.interaction_radius);
            pot_int = sparse(i(inxs),j(inxs),1,size(A,1),size(B,1));
        end
        
        function N = neighborhood_mat(obj)
            %% builds network of neighboring bins on the 2-d plane
            n_bins_dim = floor(obj.square_size/obj.interaction_radius);
            n_bins = n_bins_dim*n_bins_dim;
            B = reshape(1:n_bins,n_bins_dim,n_bins_dim);
            D = repmat((1:n_bins)',9,1);
            E = zeros(size(D));
            ct = 1;
            % covers periodic boundary conditions too
            for i=-1:1
                for j=-1:1
                    C = circshift(circshift(B,i,1),j,2);
                    E((1+(ct-1)*n_bins):(ct*n_bins)) = C(:);
                    ct = ct+1;
                end
            end
            N = sparse(D,E,1,n_bins,n_bins);
        end
    end
end

