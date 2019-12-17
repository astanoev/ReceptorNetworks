classdef util_methods < handle
    properties
    end
    
    methods
        function obj = util_methods()
        end
    end
    
    methods (Static)
        function parsave(folder, fname, vars, varnames) 
            % for saving data while in parallel mode
            try
                s = struct();
                for i=1:length(vars)
                    s.(varnames{i}) = vars{i};
                end
                if ~exist(folder,'dir'); mkdir(folder); end
                save(fullfile(folder,fname), '-struct', 's');
            catch
                disp('error saving!');
            end
        end
    end
end

