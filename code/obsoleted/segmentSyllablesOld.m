function regions = segmentSyllables(spectrum, params, varargin)
% function regions = segmentSyllables(spectrum, params, varargin)

% similar to finding regions exceeding minimal power, but this function separates at the
% syllable level with lower thresholds, assuming a finer spectrogram,
% and accepting smaller syllable lengths
% assumes the spectrum is found with the .fine parameters and is 
% highpassed to remove low freq noise


%% process
if nargin < 2, params = defaultParams; end;
params = processArgs(params, varargin{:});

regions = initEvents;
%% smooth the power with simple binomial (approx gaussian) window
smoothedP = smoothSignal(spectrum.totalPower,17);


        %% use marking based on simple power threshold

    % TODO: explore derivative thresholding
    exceedsThresh = (smoothedP > params.syllablePowerThresh);
    crossesThresh = [diff(exceedsThresh) 0];
    risesThresh = find(crossesThresh == 1);
    fallsThresh = find(crossesThresh == -1);

    %% if regions are detected, check their boundaries
    if numel(risesThresh) + numel(fallsThresh) > 0
        % boundary conditions - will fix later
        if numel(fallsThresh) == numel(risesThresh) + 1
            risesThresh = [1 risesThresh];
        elseif numel(risesThresh) > numel(fallsThresh),
            fallsThresh = [fallsThresh numel(spectrum.times)];
        elseif numel(risesThresh) == numel(fallsThresh) && ...
                risesThresh(1) > fallsThresh(1)
            risesThresh = [1 risesThresh];
            fallsThresh = [fallsThresh numel(spectrum.times)];
        end
    end
    if params.verbose
        fprintf('%d syllables originally detected\n',numel(risesThresh));
    end

    %% extend all syllables to where the power is increasing
    dPower = diff(smoothedP) / (params.fine.windowSize - params.fine.nOverlap);
    for istat = 1:numel(risesThresh)
        newRise = find([ -1 dPower(1:risesThresh(istat))] < params.syllableMinRise, ...
            1, 'last');
        newFall = find([dPower(fallsThresh(istat):end) 1 ] > -params.syllableMinRise, ...
            1, 'first') + fallsThresh(istat) - 1;
        risesThresh(istat) = newRise;
        fallsThresh(istat) = newFall;
    end

    %% combine any segments that are close by
    doLink = (risesThresh(2:end) - fallsThresh(1:end-1) < ...
        params.syllableComboLength / (params.fine.windowSize - params.fine.nOverlap) );
    if any(doLink)
        fallsThresh([doLink false]) = [];
        risesThresh([false doLink]) = [];
    end
    if params.verbose
        fprintf('%d syllables linked into larger syllables\n',sum(doLink));
    end

    %% throw away any marks where window is too small
    minSamples = ceil(params.syllableMinLength / (params.fine.windowSize - params.fine.nOverlap));
    isTooShort = (fallsThresh - risesThresh < minSamples);
    risesThresh(isTooShort)=[];
    fallsThresh(isTooShort)=[];
    if params.verbose
        fprintf('%d syllables too short to include\n',sum(isTooShort));
    end
    %% convert to events
    regions = codeEvents([risesThresh' fallsThresh'], regions, 1, spectrum, params);
    if params.plot
        plotAllFigures(spectrum,regions,params);
    end
    end