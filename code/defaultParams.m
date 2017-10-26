function params = defaultParams
%DEFAULTPARAMS returns parameters for parsing functions
%
%    function params = defaultParams returns the parameter structure for
%    the parsing codebase.  Brief descriptions and places where the
%    parameters are used are within the .m file.

% %%%%%%%%%%%%%%%%%%%%% parameters to control overall program behavior 
% for many functions
params.quiet = true; % do we want to play the whole clip?
params.plot = false; % do we want to plot figures?
params.pause = false; % do we want to pause when we are plotting?
params.playsample = false; % do we want to play any selected song regions?
params.verbose = true; % do we want to print messages?
params.saveplot = false; % do we want to automatically save plots?

% %%%%%%%%%%%%%%%%%%%%% parameters to control plotting appearances
% used in plotAllFigures.m and plot directory
params.optGraphs = {'waveform','deriv','totalPower','wienerEntropy'};
params.showLabels = false; % do we want to show the region label?
% equivalent of contrastSlider in SAP
params.dgram.minContrast = 1e-8;
params.showDgramRegionStyle = 'none'; % {none|lines|fancy}
% sampling rate, passed for convenience
params.fs = NaN; % this should not need to be changed
% here the "songStructure" is the raw spike 6 output when we export 
% channels as .mat files.

% NB: we pass around a songStructure to most main functions which contains an
% interval field that tells us the sampling rate
% whenever we pass wave data without the songStructure, we need to define
% params.fs, either a priori or from the songStructure
% here we declare but do not fill the field

% %%%%%%%%%%%%%%%%%%%%% parameters for rough parsing
% (stepSpectrogram.m)
params.Nsplits = 100; % the number of parts to scan [#]
params.overlapSplit = 300; %ms, needed to make sure all the syllables are detected in their entirety

% initial event filtering - rough pass for findPossibleSounds.m
% these would be equivalent to sliding thresholds in SAP/evsonganaly

params.powerThresh = 1e-3;  % How loud does an initial sound need to be? (RMS power)
                             % a plausible range would be 10^-6 to 10^-3.
                             % NB: This should change often depending on
                             % the recording environment.
                             
params.riseThresh = 10^-1.6; % RMS power / s, if a sound increases in power 
                             % by this much, it is also counted.
                             % a plausible range would be 10^-3 to 10^-1.
params.minCenterFreq = 1000; % Hz, A sound needs to be at least this mean frequency to pass
params.minLength = 30;       % ms, A sound should be at least this long 
params.minRise = 1e-6;       % The threshold for setting boundaries of a sound - 
                             % at what minimum change in power / s do we put the boundary of a sound
                             % This should be small (1e-8 - 1e-5) to avoid
                             % cutting off sounds
                             
% %%%%%%%%%%%%%%%%%%%%% spectrogram parameters
% spectrogram parameters for a rough pass 
params.coarse.NfreqBands = 256; %#
params.coarse.windowSize = 30.0; %ms
params.coarse.nOverlap = 24.0; %ms
params.coarse.pitchValueIfNaN = 100; % Hz, used during getMTSpectrumStats
params.coarse.multitaper = false;
params.coarse.features = {'totalPower','wienerEntropy','deriv','centerFreq','fundamentalFreq','waveform'};

% spectrogram parameters for a rough pass 
params.rough.NfreqBands = 256; %#
params.rough.windowSize = 30.0; %ms
params.rough.nOverlap = 24.0; %ms
params.rough.pitchValueIfNaN = 100; % Hz, used during getMTSpectrumStats
params.rough.multitaper = true;
params.rough.features = {'totalPower','wienerEntropy','deriv','centerFreq','fundamentalFreq','waveform'};

% spectrogram parameters - intermediate pass (can't be too fine or we run out of memory)
params.inter.NfreqBands = 1024; %#
params.inter.freqBands = linspace(1,8192,params.inter.NfreqBands); % the frequencies to read
params.inter.windowSize = 9.27; %ms
params.inter.nOverlap = 9.27 - 1.36; %ms
params.inter.pitchValueIfNaN = 100; % Hz, used during getMTSpectrumStats
params.inter.multitaper = true;
params.inter.features = {'totalPower','wienerEntropy','deriv','centerFreq','fundamentalFreq','waveform'};

% spectrogram parameters for a noise filtering pass
params.noiseReduce.NfreqBands = 4096; %#
params.noiseReduce.windowSize = 12.0; %ms
params.noiseReduce.nOverlap=9.0; %ms
params.noiseReduce.pitchValueIfNaN = NaN; % Hz, used during getMTSpectrumStats
params.noiseReduce.multitaper = false;
params.noiseReduce.features = {'waveform'};

% spectrogram parameters - fine pass (can't be too fine or we run out of memory)
params.fine.NfreqBands = 2048; %#
params.fine.freqBands = linspace(1,8192,params.fine.NfreqBands); % the frequencies to read
params.fine.windowSize = 9.27; %ms
params.fine.nOverlap = 9.27 - 1.36; %ms
params.fine.pitchValueIfNaN = 100; % Hz, used during getMTSpectrumStats
params.fine.multitaper = true;
params.fine.features = {'totalPower','wienerEntropy','deriv','centerFreq','fundamentalFreq','waveform'};

% spectrogram parameters - fine pass (can't be too fine or we run out of memory)
params.ultrafine.NfreqBands = 2048; %#
params.ultrafine.freqBands = linspace(1,8192,params.fine.NfreqBands); % the frequencies to read
params.ultrafine.windowSize = (9.27)/1.5; %ms
params.ultrafine.nOverlap = (9.27 - 1.36)/1.5; %ms
params.ultrafine.pitchValueIfNaN = 100; % Hz, used during getMTSpectrumStats
params.ultrafine.multitaper = true;
params.ultrafine.features = {'totalPower','wienerEntropy','deriv','centerFreq','fundamentalFreq','waveform'};
     
params.best.NfreqBands = 2048 * 2; %#
params.best.freqBands = linspace(1,10240,params.fine.NfreqBands); % the frequencies to read
params.best.windowSize = (9.27)/1.5; %ms
params.best.nOverlap = (9.27 - 1.36)/1.5; %ms
params.best.pitchValueIfNaN = 100; % Hz, used during getMTSpectrumStats
params.best.multitaper = true;
params.best.features = {'totalPower','wienerEntropy','deriv','centerFreq','fundamentalFreq','waveform'};

% %%%%%%%%%%%%%%%%%%%%% syllable isolation - fine pass segmentation
% used in segmentSyllables.m
params.parseSpecType='inter';
params.syllable.smoothingWindow = 17; % # samples, width of smoothing window

% factors that change syllable inclusion
params.syllable.minRiseLength = 0;    % milliseconds;
params.syllable.minRiseAmp = 15;      % maximum has to be this much louder than minimum
params.syllable.minPower = 3.1e-6;    % power, the minimum acceptable volume for a syllable 
                                      % (will be frequently changed at runtime)

params.syllable.flatFactor = 2.5;     % the factor at which to call minima/maxima the same
                                      % larger flat factor means rejecting
                                      % more extrema
                                      % this shouldn't need to be changed unless syllables are being neglected
                                      % (then this value is too high)
                                      % or if syllables are cut in half (then this value is too low)
                                      % range should be 1-100

% factors that change the locations of the boundaries
params.syllable.borderRise = 1e-3;    % power / ms, this threshold set on the minimal change in amplitude across the syllable
                                      % this means that the slope of the power at the boundaries must be less
                                      % this can be changed if the boundaries are too tight/loose

params.syllable.comboLength = 15; % ms, minimum tolerated distance between syllables.
                                 % syllables closer than this length are
                                 % automatically linked
 
params.syllable.minLength = 20;  % ms, minimum allowed syllable length

% %%%%%%%%%%%%%%%%%%%%% preprocessing parameters
% playback/cutting parameters
% used by getClip, getClipAndProcess
params.preroll = 30; % ms
params.postroll = 30; % ms

params.highPassFq = 400; %signal filtering parameters for rough pass, Hz
params.lowPassFq = 12000; %signal filtering parameters for fine pass, Hz

params.doFilterNoise = true;
% params.noiseFilter = 1e-8 * ones(2049,1); % uniform, not helpful
params.noiseFilter = []; % no noise filter pre-specified

% %%%%% spectral parameters for applying noise gate
% used by analyzeNoise
params.nps.attack = 0.005; %s
params.nps.hold = 0.005; %s
params.nps.release = 0.005; %s
params.nps.reduction = -20.0; %dB
params.nps.hysteresis = 0; %dB, not applied yet
params.nps.freqSmooth = 120; %Hz, halfwidth
params.nps.sampleWindow = 9; %#, flexible parameter, the window over which to find the max
% generally smaller window means less smoothing, more noise rejection

% %%%%%%%%%%%%%%%%%%%%% syllable to cluster mapping
% used by clustering folder
% and because fewer clusters the better
params.nClusters = 20; % should not exceed 26 for now, just for alphanumeric representation
params.clusterMethod = 'ward'; % look at "help linkage" for the full list

% %%%%%%%%%%%%%%%%%%%%%
% used by arrangeBouts 
params.silentPeriod = 2.0; % seconds

% %%%%%%%%%%%%%%%%%%%%% parameters for autodetecting silence (i.e. noise)
%  used by autodetectNoise
params.noiseDynamicRange = 5e-4; % (RMS power)
params.noiseVariation = 2e-4; % (RMS power)
params.noiseLength = 2.0; % seconds

% %%%%%%%%%%%%%%%%%%%%% parameters for interactive labeling/adjustment, used by plotAndAdjust
params.editSpecType = 'inter'; % which set of spectrogram parameters to use
params.adjustLabels = false; % whether or not to allow label editing

% %%%%%%%%%%%%%%%%%%%%% parameters for distance metrics between dimensions
% log space for fundamentalFrequency
% params.featureCatalog = struct('name',{'fundamentalFreq','FM','AM','wienerEntropy','aperiodicity'},... % harmonicPitch? pitchgoodness?    
%     'doLog', {true, false, false, false, true},...
%     'median', {6.534, 44.10, 0.00, -1.810, 5.386},...
%     'MAD', {0.687, 22.60, 0.5, 0.937, 0.540});

% used for timeWarped and standardDistance
% linear space for fundamental frequency
params.featureCatalog = struct('name',{'fundamentalFreq','FM','AM','wienerEntropy','pitchGoodness'},... % harmonicPitch? pitchGoodness?    
    'doLog', {false, false, false, false, true},...
    'median', {688.2, 44.10, 0.00, -1.810, 5.386},...
    'MAD', {168.9, 22.60, 0.5, 0.937, 0.540},...
    'weights', {2 1 1 1 2});
params.warpingCost = 0;
params.maxWarpAllowed = 0;

% used for fixBoundaries/findSilence
params.boundaryAdvance = 0.05; % seconds
params.minimumSilentGap = 0.02; % second

% used for findSilenceBaseline
params.silenceThresh = 0.5; % between 0 and 1
% used for reportOnData
params.rejectSpikeFiles = true;
params.rejectNoNeuronSessions = false;

% used for getFeatures
% the features for summarization in the 
params.reduceFeatures = {'AM','FM','pitchGoodness','wienerEntropy','fundamentalFreq','times'}; 
params.removeFeatures = {'harmonicPitch'};

% used for mosaicDRspec
params.maxMosaicLength = 7.0; %seconds
params.mosaicRowLength = 1.0; %seconds

% used for makeDerivedMotifs
params.derivedMotifMaxGap = 0.3; % seconds

% used for drawClustersUpdated
params.manuallyUpdateClusters = false;
end

