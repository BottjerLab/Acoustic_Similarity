function [candIdx, codebookCenters, callInds] = callRejector(songStruct, ROIs, possibleSyllables)

%possibleSyllables = candidateSyl;
% criteria: calls are unmodulated in power, short, stereotypical length (100 ms), and appear in sequences with breaks of 40 ms - 2s (so not in a motif)

% often have low frequency modulation, and low AM

%ROIs = bouts(11:5:31);
candDiffs = [];
candIdx = [];
durrDiff = [];
assignment = findParent(ROIs, possibleSyllables);
possibleSyllables(isnan(assignment)) = [];
assignment(isnan(assignment)) = [];

fprintf('%d syllables in %d ROIs...\n', numel(possibleSyllables), numel(ROIs));
% some programmable constants
Nsamples = 10;
sampLength = 3;

NKP = 75; % size of the vector from extractKPFeatures, TODO: replace with programmatic finding of size
KPfeats = zeros(Nsamples * numel(possibleSyllables), NKP);

for ii = 1:numel(ROIs)
	syllsInBout = possibleSyllables(assignment==ii);
	rIdxBout = find(assignment==ii);
    nSib = numel(syllsInBout);
    
    Nsamples = 10;
    sampWidth = 3;
	durations = zeros(1,nSib);
    
	for jj = 1:nSib
		[features(rIdxBout(jj)), paddedSpecs(jj)] = getFeatures(songStruct, syllsInBout(jj));
		durations(jj) = syllsInBout(jj).stop - syllsInBout(jj).start;
        % extract key features
        tryFeat = extractKPFeatures(paddedSpecs(jj), Nsamples, sampWidth);
        KPfeats((1:Nsamples) + Nsamples * (rIdxBout(jj) - 1),:) = tryFeat;
    end

    
	% find similar calls next to each other
    durDiffTol = 0.01; %s
    minExpectedPower = 5e-4; %raw power, normal 
	if any(abs(diff(durations)) < durDiffTol) 
		% compare the feature vectors
        simToNext = find(abs(diff(durations)) < durDiffTol);
		for jj = 1:numel(simToNext)           
			idx = simToNext(jj);		
            
            % syllables have to have a minimum power to qualify as same
            if features(rIdxBout(idx  )).totalPower_nanmean <= minExpectedPower || ...
               features(rIdxBout(idx+1)).totalPower_nanmean <= minExpectedPower,
                continue;
            end           
			% convert labeled features to vectors
            fV1 = cellfun(@(x) features(rIdxBout(idx)).(x), fieldnames(features))';
            fV2 = cellfun(@(x) features(rIdxBout(idx+1)).(x), fieldnames(features))';
            fprintf('Found potential match in syllable #%d...\n', rIdxBout(idx));
              
            calldiff = norm(fV1-fV2);
            candIdx = [candIdx rIdxBout(idx)];
            candDiffs = [candDiffs calldiff];
            durrDiff = [durrDiff durations(idx+1) - durations(idx)];
        end
    end
    clear features
end

minMatchScore = 1e3; % empirically found
isGoodCandMatch = (candDiffs < minMatchScore);
% now how to do some blind signal separation using these candidates as a
% seed

% construct clusters based on the time samples
Nclusters = 20;
[codeIdxs, codebookCenters] = kmeans(KPfeats, Nclusters, ...
    'emptyaction','drop', 'replicates', 30, ...
    'options',struct('Iterations', 150));

% now look at the histograms of the heuristically found calls
%callInds = codeIdxs((candIdx(ones(Nsamples,1),:))*10 + (1:10)' * ones(1,numel(candIdx)))
callInds = reshape(codeIdxs,10,numel(possibleSyllables));
%%note: scratches usually follow a exp decay like signal but this is hard
%%to find 

%input: ROIs/syllables
%goal: look for ROIs with repeated syllables of equal length

%each time a repeated 'syllable' appears, cast a vote for that syllable length as a call

%each call needs a specific fingerprint (more than just length)

%requires data structure:
%vote bins (syllable length), which vote 


%exclude ROIs that are exclusively calls

%trim ROIs around calls