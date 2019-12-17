classdef multistep_experiment < experiments.experiment
    properties
        n_steps
        step_duration_min
        step_duration_frames
        step_amplitude
    end
        
    methods
        function obj = multistep_experiment(n_steps, model)
            obj@experiments.experiment();
            obj.n_steps = n_steps;
            
            obj.fr_per_sec = 1;
            obj.pre_stimulus_min = 5;
            obj.post_stimulus_min = 30;
            obj.step_duration_min = 20;
            
            obj.t_total_min = obj.pre_stimulus_min + (2*obj.n_steps-1)*obj.step_duration_min + obj.post_stimulus_min;
                        
            if nargin >=2; obj.set_up_model(model); end
            obj.set_up_time();
        end
                
        function step_duration_frames = get.step_duration_frames(obj)
            step_duration_frames = obj.step_duration_min*obj.minute_to_frames;
        end
        
        function input = set_up_input(obj)
            target_egf_egfr_p = 0.012;
            obj.step_amplitude = obj.model.tune_egf_free(target_egf_egfr_p); %pulse_amp = 1.5;
            input = zeros(size(obj.time));
            for i=1:(2*obj.n_steps-1)
                t_start = obj.pre_stimulus_frames + (i-1)*obj.step_duration_frames;
                t_end = t_start + obj.step_duration_frames;
                if i<=obj.n_steps
                    input(t_start:t_end) = obj.step_amplitude*i;
                else
                    input(t_start:t_end) = obj.step_amplitude*(obj.n_steps-mod(i,obj.n_steps));
                end
            end
            obj.input = input;
        end
    end
end