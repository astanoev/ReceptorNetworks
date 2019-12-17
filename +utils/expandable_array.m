classdef expandable_array < handle
    %EXPANDABLE_ARRAY - An array that automatically expands as more data is
    %appended. This is done at every 25% increase in size of the currently
    %stored array (matrix).
    
    properties
        exp_array;
        increase_percentage = 1.25;
        along_dimension = 1; % 1 - add rows to matrix, 2 - add columns
        current_size = 0;
    end
    
    methods
        function obj = expandable_array(varargin)
            p = inputParser;
            addOptional(p,'exp_array',zeros(0,1));
            addOptional(p,'increase_percentage',1.25);
            addOptional(p,'along_dimension',1);
            parse(p,varargin{:});
            
            obj.exp_array = p.Results.exp_array;
            obj.increase_percentage = p.Results.increase_percentage;
            obj.along_dimension = p.Results.along_dimension;
            obj.current_size = size(obj.exp_array, obj.along_dimension);
            % expand array initially, +1 added for case where array has
            % zero rows/columns
            if obj.along_dimension==1 
                obj.exp_array(ceil((obj.current_size+1)*obj.increase_percentage),end) = 0;
            else
                obj.exp_array(end,ceil((obj.current_size+1)*obj.increase_percentage)) = 0;
            end
        end
        
        function array_ret = get(obj)
            if obj.along_dimension == 1
                array_ret = obj.exp_array(1:obj.current_size,:);
            else
                array_ret = obj.exp_array(:,1:obj.current_size);
            end
        end
        
        function extend(obj, array_ext)
            if isempty(array_ext); return; end
            c_len = size(obj.exp_array,obj.along_dimension);
            cur_len = size(obj.exp_array,obj.along_dimension);
            len = size(array_ext,obj.along_dimension);
            new_size = obj.current_size + len;
            while new_size>cur_len; cur_len = cur_len*obj.increase_percentage; end
            if cur_len>c_len
                if obj.along_dimension==1; obj.exp_array(ceil(cur_len),end) = 0; else; obj.exp_array(end,ceil(cur_len)) = 0; end
            end
            if obj.along_dimension==1
                obj.exp_array((obj.current_size+1):new_size,:) = array_ext;
            else
                obj.exp_array(:,(obj.current_size+1):new_size) = array_ext;
            end
            obj.current_size = new_size;
        end
    end
end

