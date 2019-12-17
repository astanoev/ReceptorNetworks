classdef exponential_decay_model < models.model
    properties
        par = struct('beta',0.0035,'alpha',0.023,'Xt',0.87);
    end
    
    properties
    end
    
    methods
        function obj = exponential_decay_model()
            %obj@models.model();
            obj.set_initial_cond();
            obj.labels = {'EGFRpt'};
        end
        
        function obj = set_initial_cond(obj)
            obj.init_conds = [0.01];
            %obj.traj_cols = [0.3, 0.8, 0.4; 0.9, 0.1, 0.5];
            obj.traj_cols = [1.0, 0.34, 0.34];
            %obj.traj_cols = [0.39, 0.39, 1.0];
            %obj.marker_cols = [0, 1, 0; 1, 0, 0];%[1, 0.6, 0; 0.235, 0.67, 0.9];
            obj.marker_cols = [1, 0, 0; 0, 0, 1];
            %obj.marker_cols = [0, 0, 1];
        end
        
        function set_Xt_ss(obj, X_ss, I)
            obj.par.Xt = X_ss*(obj.par.beta + I*obj.par.alpha)/(I*obj.par.alpha);
        end
                        
        function [dydt] = df_model(obj, t, y, experiment)
            %t_vec, EGF_vec
            %[~, index] = min(abs(t_vec-tt));
            %EGFfree = EGF_vec(index);
            input = experiment.get_input(t);
            
            dydt = input.*obj.par.alpha*(obj.par.Xt-y) -obj.par.beta*y;
        end
    end
end