function sample = highPassSample(sample,params)
% HIGHPASSSAMPLE  High pass filters a sample
%
%   SAMPLE = highPassSample(SAMPLE,PARAMS) filters out all frequencies
%   below a certain frequency, which is defined in the parameters as 
%   'highPassFq'.  A sampling rate must be defined in the
%   parameters as 'fs'.  

    if ~isfield(params, 'fs')
        error('highPassSample:Undefined','Undefined sampling frequency');
    end
    % 2nd order butterworth filter, supposedly...
    [bFilt,aFilt] = butter(2,params.highPassFq/(params.fs/2),'high');
    sample = filtfilt(bFilt,aFilt, sample);
end
