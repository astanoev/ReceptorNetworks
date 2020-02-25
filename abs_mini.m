classdef abs_mini
    % reduced information from agent_based_simulation
    
    properties
        model;
        R_states;
        P_dnf_states;
        R_positions;
        P_dnf_positions;
        rng_final;
    end
    
    methods
        function obj = abs_mini(abs)
            obj.model = abs.model;
            obj.R_states = abs.R_states(end,:,:);
            obj.P_dnf_states = abs.P_dnf_states(end,:,:);
            if isempty(abs.R_positions)
                abs.calculate_positions();
            end
            obj.R_positions = abs.R_positions(end,:,:);
            obj.P_dnf_positions = abs.P_dnf_positions(end,:,:);
            obj.rng_final = abs.rng_final;
        end
    end 
end

