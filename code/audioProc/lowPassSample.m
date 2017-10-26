function sample = lowPassSample(sample,params)
    
    if params.fs/2 < params.lowPassFq
        warning('lowPassSample:nyquist', ...
            'Low pass frequency [%0.1f Hz] exceeds nyquist frequency', ...
            params.lowPassFq);
        return;
    end;
    [bFilt,aFilt] = butter(4,params.lowPassFq*2/params.fs,'low');
    sample = filtfilt(bFilt,aFilt, sample);
end
