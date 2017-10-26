
files = dir([matpath 'spike*']); files = {files.name};
shortNames = strrep(strrep(files, 'spikeRates-',''),'_times.mat','');

basalRate = []; spikeRate = []; isShell = []; pp = [];
for ii = 1:numel(shortNames)
    spikeData{ii} = load([matpath files{ii}]); 
    basalRate = [basalRate spikeData{ii}.basalRate];
    spikeRate = [spikeRate spikeData{ii}.spikeRate];
    isShell   = [isShell   spikeData{ii}.isShell  ];
    pp        = [pp        spikeData{ii}.rate_tp  ];
end
isShell = logical(isShell);

pSig = 0.005;
col = {'b','r'}; mSize = [20,5];
for ii = 0:1
    for jj = -1:2:1
        flags = isShell == ii & sign(pp - pSig) == jj;
        plot(basalRate(flags), spikeRate(flags), [col{ii+1} '.'],'MarkerSize',mSize((jj+3)/2));
        hold on;
    end
end
maxCoord = max(max(ylim), max(xlim));
plot([0 maxCoord], [0 maxCoord],'k-');
hold off;
xlabel('Basal rate (Hz)'); ylabel('Song spike rate (Hz)');