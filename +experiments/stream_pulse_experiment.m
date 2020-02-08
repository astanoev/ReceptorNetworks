classdef stream_pulse_experiment < experiments.experiment
    properties
        pulse_duration_min
        pulse_duration_frames
        pulse_amplitude
        interpulse_duration
        lambda % number of pulses per hour (60minute_to_frames)
        n_pulses
    end
        
    methods
        function obj = stream_pulse_experiment(n_pulses, lambda, t_total_min, model)
            obj@experiments.experiment();
            obj.n_pulses = n_pulses;
            obj.lambda = lambda;
            obj.t_total_min = t_total_min;
            
            obj.pre_stimulus_min = 2;
            obj.post_stimulus_min = 2.5;
            obj.pulse_duration_min = 5;
            
            if nargin >=4; obj.set_up_model(model); end
            obj.set_up_time();
        end
        
        function set_up_model(obj, model)
            obj.model = model;
        end
        
        function pulse_duration_frames = get.pulse_duration_frames(obj)
            pulse_duration_frames = obj.pulse_duration_min*obj.minute_to_frames;
        end
                
        function input = set_up_input(obj)
            if isempty(obj.pulse_amplitude)
                target_LRa = 0.15;
                obj.pulse_amplitude = obj.model.tune_Lt(target_LRa);
            end
            obj.set_up_time();
            obj.input = zeros(length(obj.time),1);
            times_pulses_frames = [0;randi(obj.t_total_frames - (obj.pre_stimulus_frames +obj.post_stimulus_frames) - (obj.n_pulses+1)*obj.pulse_duration_frames,obj.n_pulses-1,1)];
            times_pulses_frames = sort(times_pulses_frames) + (0:(obj.n_pulses-1))'*obj.pulse_duration_frames;
            for i=1:size(times_pulses_frames,1)
                t_start = obj.pre_stimulus_frames + times_pulses_frames(i);
                t_end = t_start + obj.pulse_duration_frames -1;
                obj.input(t_start:t_end,:) = repmat(obj.pulse_amplitude, obj.pulse_duration_frames, 1);
            end
            input = obj.input;
        end
    end
end