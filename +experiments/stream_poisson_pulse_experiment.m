classdef stream_poisson_pulse_experiment < experiments.experiment
    properties
        pulse_duration_min
        pulse_duration_frames
        interpulse_duration
        lambda % number of pulses per hour (60minute_to_frames)
        n_pulses
    end
        
    methods
        function obj = stream_poisson_pulse_experiment(n_pulses, lambda, t_total_min, model)
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
            target_egf_egfr_p = 0.15;
            pulse_amp = obj.model.tune_egf_free(target_egf_egfr_p);
            obj.set_up_time();
            obj.input = zeros(length(obj.time),1);
%             p = rand(obj.n_pulses,1);
%             r = -1/(obj.lambda/60)*log(1-p);
%             times = (cumsum(r)-r(1))*obj.minute_to_frames + [0:(obj.n_pulses-1)]'*obj.pulse_duration;
            %times_pulses_frames = [0;(obj.t_total_frames - (obj.pre_stimulus_frames +obj.post_stimulus_frames) - (obj.n_pulses+1)*obj.pulse_duration_frames)*rand(obj.n_pulses-1,1)];
            times_pulses_frames = [0;randi(obj.t_total_frames - (obj.pre_stimulus_frames +obj.post_stimulus_frames) - (obj.n_pulses+1)*obj.pulse_duration_frames,obj.n_pulses-1,1)];
            times_pulses_frames = sort(times_pulses_frames) + [0:(obj.n_pulses-1)]'*obj.pulse_duration_frames;
            for i=1:size(times_pulses_frames,1)
                t_start = obj.pre_stimulus_frames + times_pulses_frames(i);
                t_end = t_start + obj.pulse_duration_frames -1;
                obj.input(t_start:t_end,:) = repmat(pulse_amp, obj.pulse_duration_frames, 1);
            end
            input = obj.input;
        end
    end
end