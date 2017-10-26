function clicks = findClicks(songStruct, events, params)

if nargin < 3
	params = defaultParams;
end

fs = 1/songStruct.interval;
   
for ii = 1:numel(events)
	thisEvent = events(ii);
	% high pass the clip to 7500Hz
	clipOnly = getClipAndProcess(songStruct, thisEvent, params);
	highCl = getClipAndProcess(songStruct, thisEvent, params, 'highPassFq', 9000);
    
    %playSound(clipOnly, fs, true);
    %beep;
    %playSound(highCl, fs, true);
	% in general, any fast clips are clicks
	specparams = params.(params.editSpecType);
    specparams.fs = fs;
	clipspec = getMTSpectrumStats(clipOnly, specparams);
	highspec = getMTSpectrumStats(highCl,   specparams);
	
	highPower = highspec.totalPower;
	highFrac = highPower ./ clipspec.totalPower;
	
	dHP = [0 diff(highPower)];
	dHF = [0 diff(highFrac )];
	
    figure(1)
    times = highspec.times;
    xx = [min(times) max(times)];
	subplot(511)
	semilogy(times, highPower,'r-'); xlim(xx);
	ylabel('power in high freq band'); 
	subplot(512)
	plot(times, highFrac, 'r-'); ylim([0 0.2]); xlim(xx);
	ylabel('fraction of total power in high freq band');
	subplot(514)
	plotSpectrogram(clipspec);
	subplot(515)
	plotSpectrogram(highspec);
%	plot(times, dHP); xlim(xx);
	ylabel('Diff power');
	subplot(513)
	plotWaveform(clipOnly, fs);% xlim(xx);
	ylabel('Diff fraction power');
    
    figure(2)
    params.fs = fs;
    plotAllFigures(clipspec, [], params, 'optGraphs', {'waveform' 'deriv' 'totalPower' 'rawAM'});
	%pause;
    % todo: filter by raw amplitude modulation, totalPower, and high power
    % fraction 
end
clicks = [];