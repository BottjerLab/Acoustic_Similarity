function cl = plotDGFromRegion(songStruct, region, params, varargin)
% plots the derivative gram directly from a single region, if you never want to
% know the feature data
% 
% this function also uses what axes are preset, allowing gridded plotting
% of clips
% songStruct is an option if region has a reference to a disk-readable file

if nargin < 3
    params = defaultParams;
end
params = processArgs(params, varargin{:});
%{
if isempty(songStruct)
    if exist(region.file,'file') ~= 2 || ~strcmp(region.file(end-3:end), '.mat') 
        error('songStruct must be defined or valid region must be provided...'
    end
end
%}
[cl, fs] = getClipAndProcess(songStruct, region, params);

fP = params.fine; fP.fs = fs; fP.features = 'deriv';
exSpec = getMTSpectrumStats(cl, fP);       
plotDerivGram(exSpec, params);

% plot 10 ms scale bar
% starts 80% across the x axis and extends by the scaled length--JMA took
% this out
xx = xlim * [0.8 ; 0.2]; 
yy = ylim * [0.1 0.1; 0.9 0.9]; % 90% up the y axis
hold on; 
plot([xx xx-0.01],yy,'w-','LineWidth', 2);
hold off;

