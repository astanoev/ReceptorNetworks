classdef exponential_decay_model < models.model
    properties
        par = struct('beta',0.0035,'alpha',0.023,'Rt',0.87);
    end
    
    properties
    end
    
    methods
        function obj = exponential_decay_model()
            obj.set_initial_cond();
            obj.labels = {'Ra'};
        end
        
        function obj = set_initial_cond(obj)
            obj.init_conds = [0.01];
            obj.traj_cols = [1.0, 0.34, 0.34];
            obj.marker_cols = [1, 0, 0; 0, 0, 1];
        end
        
        function set_Rt_ss(obj, Ra_ss, I)
            obj.par.Rt = Ra_ss*(obj.par.beta + I*obj.par.alpha)/(I*obj.par.alpha);
        end
                        
        function [dydt] = df_model(obj, t, y, experiment)
            input = experiment.get_input(t);
            dydt = input.*obj.par.alpha*(obj.par.Rt-y) -obj.par.beta*y;
        end
    end
end