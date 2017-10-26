function analyzeClip(clip, fs, varargin)
% function ANALYZECLIP
% simple, helpful script, just run whatever analysis we do on any sound vector with
% fs - can be used to run a much simpler version of plotForReview

params = processArgs(defaultParams, varargin{:});

params.fine.fs = fs;
spec = getMTSpectrumStats(clip, params.fine);

plotAllFigures(spec, [], params);
%title(sprintf('PSD from %s--%s with centroid freq.', ...
%    sampleToTimeString(sampleStart), ...
%    sampleToTimeString(sampleEnd)));
drawnow;

if ~params.quiet
    params.fs = fs;
    player = playSound(clip,params.fs,true);
end
end
