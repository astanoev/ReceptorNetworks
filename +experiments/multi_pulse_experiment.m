classdef multi_pulse_experiment < experiments.experiment
    properties
        n_pulses
        pulse_duration_min
        pulse_duration_frames
        pulse_interval_min
        pulse_amplitude
        between_pulses_min
        between_pulses_frames
    end
        
    methods
        function obj = multi_pulse_experiment(n_pulses, model)
            obj@experiments.experiment();
            obj.n_pulses = n_pulses;
            
            obj.fr_per_sec = 1;
            obj.pre_stimulus_min = 5;
            obj.post_stimulus_min = 15;
            obj.pulse_duration_min = 5;
            obj.pulse_interval_min = 20;
            
            obj.t_total_min = obj.pre_stimulus_min + obj.n_pulses*(obj.pulse_duration_min + obj.between_pulses_min) + obj.post_stimulus_min;
            
            if nargin >=2; obj.set_up_model(model); end
            obj.set_up_time();
        end
        
        function pulse_duration_frames = get.pulse_duration_frames (obj)
            pulse_duration_frames = obj.pulse_duration_min*obj.minute_to_frames;
        end
        
        function between_pulses_min = get.between_pulses_min(obj)
            between_pulses_min = obj.pulse_interval_min-obj.pulse_duration_min;
        end
        
        function between_pulses_frames = get.between_pulses_frames(obj)
            between_pulses_frames = obj.between_pulses_min*obj.minute_to_frames;
        end
        
        function set_pulse_interval_min(obj, pulse_interval_min)
            obj.pulse_interval_min = pulse_interval_min;
            obj.t_total_min = obj.pre_stimulus_min + obj.n_pulses*(obj.pulse_duration_min + obj.between_pulses_min) + obj.post_stimulus_min;
            obj.set_up_time();
        end
                
        function input = set_up_input(obj)
            if isempty(obj.pulse_amplitude)
                target_LRa = 0.15;
                obj.pulse_amplitude = obj.model.tune_Lt(target_LRa); %pulse_amp = 1.5;
            end
            input = zeros(size(obj.time));
            for i=1:obj.n_pulses
                t_start = obj.pre_stimulus_frames + (i-1)*(obj.pulse_duration_frames +obj.between_pulses_frames);
                t_end = obj.pre_stimulus_frames + (i-1)*(obj.pulse_duration_frames +obj.between_pulses_frames) + obj.pulse_duration_frames;
                input(t_start:t_end) = input(t_start:t_end) + obj.pulse_amplitude;
            end
            obj.input = input;
        end
    end
end