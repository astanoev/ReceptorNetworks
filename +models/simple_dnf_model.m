classdef simple_dnf_model < models.simple_model
    properties
        
    end
    
    methods
        function obj = simple_dnf_model()
            %obj@models.simple_model();
            obj.par.g1 = 2.957; %2.5 - bistability, 2.957 - memory, 3.5 - rev. bistability, 4.3 - monostability
            obj.par.g2 = 0;
            obj.par.g3 = 0;
        end
    end
end