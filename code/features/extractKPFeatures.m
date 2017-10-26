function kpFeats = extractKPFeatures(spec, Nsamples, sampWidth)
	featNames = whichFeaturesKP(spec);
	specWidth = numel(spec.times);
	timePoints = floor(linspace(1,specWidth - sampWidth, Nsamples));
	for ii = 1:Nsamples
        idx = timePoints(ii);
		fV = cellfun(@(x) spec.(x)(:, idx:idx+sampWidth-1), featNames,...
'UniformOutput', false);
        fV = vertcat(fV{:});
		kpFeats(ii,:) = fV(:);
end
end