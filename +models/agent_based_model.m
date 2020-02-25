classdef agent_based_model < handle
    %AGENT_BASED_MODEL A model that is analogous to the deterministic
    %reaction-diffusion PDE model, using the microscopic rate constants,
    %which then extends to the ODE model with the macroscopic rate
    %constants
    
    properties
        k1=5; % arbitrary scaling, preserving reaction prob. k1*dt<=1
        k2=0.5;
        % Second-order reaction rate constants are normalized to the cross
        % section area where the particles interact
        % In essence, microscopic rates can be related to macroscopic rates
        % with k_mic = k_mac/(1-k_mac/(4*pi*D_k))
        % and k_mac = k_mic/(1+k_mic/(4*pi*D_k))
        % maximum 2nd order parameter value should be 1/obj.model.delta_t*(obj.model.interaction_radius^2*pi)
        % as the probability for interaction surpasses one if higher, so
        % delta_t should be lowered (if delta_t = 0.001 and R=0.02, then
        % 1.2566 is the max value)
        a1=0.0017/(models.agent_based_model.interaction_radius^2*pi);
        a2=0.3/(models.agent_based_model.interaction_radius^2*pi);
        a3=1.0/(models.agent_based_model.interaction_radius^2*pi);
        b1=9.0/(models.agent_based_model.interaction_radius^2*pi);
        g1=0.35/(models.agent_based_model.interaction_radius^2*pi);
        R_avg_per_um_sq = 60; % avg no. of molecules per um^2
        P_dnf_avg_per_um_sq = 80; % 60 would correspond to 10^5 molecules per 1650um^2 surface area (MCF7 cell)
        Rt;
        P_dnf_t;
        square_size = 3.5; % passed by abs
        time_steps = 0;  % passed by abs
    end
    
    properties (Constant)
        diffusion = 0.1; % um^2/s
        particle_radius = 0.01; % 10nm
        interaction_radius = 2*models.agent_based_model.particle_radius;
        delta_t = 0.1*models.agent_based_model.interaction_radius^2/(4*models.agent_based_model.diffusion); % s - maximal dt assuring sqrt(4*(2*D)*dt)=sqrt(msd)<=interaction_radius
    end
    
    methods
        function obj = agent_based_model(time_steps, square_size)
            if nargin>0
                obj.time_steps = time_steps;
            end
            if nargin>1
                obj.square_size = square_size;
            end
            
            obj.Rt = round(obj.square_size^2*obj.R_avg_per_um_sq);
            obj.P_dnf_t = round(obj.square_size^2*obj.P_dnf_avg_per_um_sq);
        end
        
        function reset_model(obj, par_name, par_value)
            pars = {'a1','a2','a3','b1','g1','k1','k2'};
            for i=1:length(pars)
                obj.(pars{i}) = 0.0;
            end
            obj.(par_name) = par_value/(models.agent_based_model.interaction_radius^2*pi);
            if strncmp(par_name,'k',1)
                obj.(par_name) = par_value;
            end
        end
    end
end