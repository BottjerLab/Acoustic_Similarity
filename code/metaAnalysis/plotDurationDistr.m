function syllableLengths = plotDurationDistr(sylls)

syllableLengths = [sylls.stop]-[sylls.start];
% plotting
binWidth = 0.005;
syllBins = 0:binWidth:max(syllableLengths); % critical; describes syllable bins

binCenters = syllBins + binWidth/2;
syllPDF = hist(syllableLengths,syllBins) / numel(syllableLengths);
bar(binCenters, syllPDF, 1)