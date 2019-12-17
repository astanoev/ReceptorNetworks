classdef simple_convergence_experiment < experiments.experiment
    %SIMPLE_CONVERGENCE_EXPERIMENT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
    end
    
    properties %(Access=private)
        tmax
        fr_per_tp
        t
    end
    
    methods
        function obj = simple_convergence_experiment(model, tmax)
            obj@experiments.experiment();
            if nargin >=1; obj.set_up_model(model); end
            if nargin >=2; obj.tmax = tmax; end
            obj.set_up_time();
        end
        
        function set_up_model(obj, model)
            obj.model = model;
        end
        
        function set_up_time(obj)
            if isempty(obj.tmax); obj.tmax = 5000; end
            obj.fr_per_tp = 1;
            obj.t = 0:1/obj.fr_per_tp:(obj.tmax-1/obj.fr_per_tp);
        end
        
        function egf_free_t = set_up_input(obj)
            target_egf_egfr_p = 0;%0.1;
            pulse_amp = obj.model.tune_egf_free(target_egf_egfr_p); %pulse_amp = 1.5;
            egf_free_t = target_egf_egfr_p*ones(size(obj.t));
        end
    end
end

