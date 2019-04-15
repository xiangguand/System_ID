% This function is aimed to get amplitude and phase from the input signal.
% Make sure input signal can only solve one dimension.
function [amplitude, phase] = getSineAmplitudeAndPhase(signal)
    minimal_gap = 1e-9;
    max_s = max(signal);
    min_s = min(signal);
    data_size = size(signal);
    length = data_size(1);
    if length == 1
       length = data_size(2); 
    end
    amplitude = (abs(max_s) + abs(min_s)) / 2;
    
    p1 = 0;
    p2 = 0;
    temp_phase = 0;
    first_flag = true;
    for i=1:length
%         abs(abs(signal(i)) - max_s) 
        if abs(abs(signal(i)) - max_s) < minimal_gap
           if p1 == 0
                p1 = i; 
           elseif p1 ~= 0 && p2 == 0
                p2 = i;
           end
        end
        if (abs(abs(signal(i)) - max_s) < minimal_gap) && first_flag
            temp_phase = -i;
            first_flag = false;
        elseif (abs(abs(signal(i)) - min_s) < minimal_gap) && first_flag
            temp_phase = i;
            first_flag = false;
        end
    end
    % one sine period is p2 - p1
%     p1
%     p2
%     temp_phase
    full_period_dots = p2 - p1;
    phase =  temp_phase/full_period_dots*180;
    if phase < -90
        phase = phase + 180;
    elseif phase > 180
        phase = phase - 180;
    end
end

