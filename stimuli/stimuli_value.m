function val = stimuli_value(angle,stimulis_angle,type)
% Determine pixel value [0:1] for a stimulis pattern. 

% Type alterniatives are 'tophat' or 'cosine'.
% Angles in degrees. 

switch type
    case 'tophat'
        if angle <= stimulis_angle/4
            val = 1;
        else
            val = 0;
        end
        
    case 'cosine'
        if angle >= stimulis_angle/2
            val = 0;
        else
        val = cosd(angle*90/(stimulis_angle/2));
        end
        test = nan;
    case 'gaussian'
        disp('gaussian, not implemented. use cosine instead')
        val = nan;
    otherwise
        dip('Unsupported type of stimuli')
        val = nan;
end
end