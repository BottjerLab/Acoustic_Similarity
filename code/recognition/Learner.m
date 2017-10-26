classdef Learner < handle
    properties (Dependent)   
        threshold = 0;
        fastThreshold = 0;
    end
    properties (Dependent, SetAccess = private)
        typicalLength
        nExamples
    end
    properties (SetAccess = private, GetAccess = private) 
        itsSpectra = struct([])
        threshold_ = 0;
        fastThreshold_ = 0;
        localParams
    end
    methods
        % constructor
        function model = Learner(songStruct, regions, params, varargin)
            if nargin < 3, params = defaultParams; end
            model.localParams = processArgs(params, varargin{:});
            if nargin >= 2 && numel(regions)>0,
                model.addSpectra(songStruct, regions);
            end
            
            % for speed considerations
            model.localParams.editSpecType = 'rough';
            model.localParams.(model.localParams.editSpecType).fs = 1/songStruct.interval;
        end
        
        function ret = isInited(model) 
            ret = ~isempty(model.itsSpectra);
        end
        
        function adaptThreshold(model, estFrac)
         model.localParams = processArgs(defaultParams, 'Y231early');
            % can't bootstrap if we only have one example
            if model.nExamples == 1, return; end;
            
            if nargin < 2, estFrac = sqrt(model.nExamples); end;
            % sample with replacement a random fraction of the
            % intercategory similarity (can be greater than 1)
            % and take the 75th, e.g. percentile
            % value as the category boundary
            nBootstrap = ceil(estFrac * model.nExamples ^ 2);
            scores = zeros(1,nBootstrap);
            fastScores = zeros(1,nBootstrap);
            progressbar('Bootstrapping thresholds...');
            for iter = 1:nBootstrap
                % get two unequal numbers between 1 and nExamples;
                ii = randi(model.nExamples-1);
                jj = randi(model.nExamples); 
                if ii>=jj, ii=ii+1; end;
                scores(iter) = model.getDistance(model.itsSpectra(ii),model.itsSpectra(jj));
                fastScores(iter) = model.getFastDistance(model.itsSpectra(ii), model.itsSpectra(jj));
                progressbar(iter/nBootstrap);
            end
            % get the 90th percentile for inclusion (this is an arbitrary
            % cutoff)
            threshPer = 0.9;
            model.threshold_ = quantile(scores, threshPer);
            model.fastThreshold_ = quantile(fastScores, threshPer);
        end
        
        function narrowClassifier(model, narrowNumber)
            if nargin < 2, 
                narrowNumber = ceil(sqrt(model.nExamples));
            end
            distVector = zeros(1,model.nExamples*(model.nExamples-1)/2);
            ctr = 1;
            progressbar('Narrowing classifier...')
            for ii = 1:model.nExamples
                for jj = ii+1:model.nExamples
                    distVector(ctr) = model.getDistance(model.itsSpectra(ii), model.itsSpectra(jj));
                    ctr = ctr+1;
                    progressbar(((ii-1)*model.nExamples + jj)/(model.nExamples^2));
                end
            end
            
            linkTree = linkage(distVector, 'complete');
            clustIdx = cluster(linkTree,'maxclust', narrowNumber);
            [~,uniqueRep] = unique(clustIdx);
            model.itsSpectra = model.itsSpectra(uniqueRep);
        end            
        
        function dist = getDistance(model, spec1, spec2)
            dist = timeWarpedDistance(spec1,spec2, model.localParams, 'warpingCost', 0.05);
        end
        
        function dist = getFastDistance(model, spec1, spec2)
            revisedCatalog = model.localParams.featureCatalog;
            revisedCatalog = revisedCatalog(strcmp({revisedCatalog.name}, 'wienerEntropy'));
            dist = timeWarpedDistance(spec1, spec2, ...
                model.localParams, 'featureCatalog', revisedCatalog);
        end
        
        function [isMatch, averageScore, allScores] = score(model, songStruct, regions)
            nR = numel(regions);
            isMatch = false(nR,1);
            averageScore = zeros(nR,1);
            allScores = zeros(nR,model.nExamples);
            progressbar('Analyzing regions','Calculating distance');
            for ii = 1:nR
                thisClip = getClipAndProcess(songStruct, regions(ii), model.localParams);
                
                specParams = model.localParams.(model.localParams.editSpecType);
                specParams.fs = 1/songStruct.interval;
                
                spec = getMTSpectrumStats(thisClip,specParams);
                
                model.localParams.(model.localParams.editSpecType).fs = 1/songStruct.interval;
                
                for jj = 1:model.nExamples
                    allScores(ii,jj) = model.getDistance(model.itsSpectra(jj), spec);
                    progressbar([],jj/model.nExamples);
                end
                progressbar(ii/nR,[]);
                averageScore(ii) = mean(allScores(ii,:));
                isMatch(ii) = (averageScore(ii) < model.threshold);
            end
        end
        
        function [averageScore, allScores] = fastScore(model, songStruct, regions)
            nR = numel(regions);
            %isMatch = false(nR,1);
            averageScore = zeros(nR,1);
            allScores = zeros(nR,model.nExamples);
            progressbar('Analyzing regions','Calculating distance');
            for ii = 1:nR
                thisClip = getClipAndProcess(songStruct, regions(ii), model.localParams);
                
                specParams = model.localParams.(model.localParams.editSpecType);
                specParams.fs = 1/songStruct.interval;
                
                spec = getMTSpectrumStats(thisClip,specParams);
                
                model.localParams.(model.localParams.editSpecType).fs = 1/songStruct.interval;
                for jj = 1:model.nExamples
                    allScores(ii,jj) = model.getFastDistance(model.itsSpectra(jj), spec);
                    progressbar([],jj/model.nExamples);
                end
                progressbar(ii/nR,[]);
                averageScore(ii) = mean(allScores(ii,:));
                
                %isMatch(ii) = (averageScore(ii) < model.threshold);
            end
            progressbar(1,[]);
        end
        
        function [bestEvents, bestScores] = findBestRegions(model, songStruct, regions)
            % same as score, but read output from timeWarpedDistance
            % criteria is that this should already be a pretty good match,
            % not checked here
            nR = numel(regions);
            matchInterval = zeros(nR,2);
            bestScores = zeros(nR,1);
            progressbar('Finding best matches','Finding individual matches');
            for ii = 1:nR
                thisClip = getClipAndProcess(songStruct, regions(ii), model.localParams);
                
                specParams = model.localParams.(model.localParams.editSpecType);
                specParams.fs = 1/songStruct.interval;
                
                spec = getMTSpectrumStats(thisClip,specParams);
                
                model.localParams.(model.localParams.editSpecType).fs = 1/songStruct.interval;
                matchPoss = zeros(model.nExamples, 2);
                scores = zeros(model.nExamples,1);
                for jj = 1:model.nExamples
                    [scores(jj),~,matchPoss(jj,:)] = timeWarpedDistance(model.itsSpectra(jj), spec);
                    progressbar((ii-1 + jj/model.nExamples)/nR ,jj/model.nExamples);
                end
                matchInterval(ii,:) = mean(matchPoss,1);
                bestScores(ii) = mean(scores);
                progressbar(ii/nR,[]);
            end
            bestEvents = eventFromTimes(matchInterval(:,1), matchInterval(:,2), 1/songStruct.interval);
        end
        function addSpectra(model, songStruct, regions)
            firstClip = getClipAndProcess(songStruct, regions(1), model.localParams);
            
            fs = 1/songStruct.interval;
            specParams = model.localParams.(model.localParams.editSpecType);
            specParams.fs = fs; 
            
            progressbar('adding spectra')
            firstSpec = getMTSpectrumStats(firstClip, specParams);
            itsNewSpectra = initEvents(numel(regions), firstSpec);
            itsNewSpectra(1) = firstSpec;
            progressbar(1/numel(regions));
            for ii = 2:numel(regions)
                thisClip = getClipAndProcess(songStruct, regions(ii), model.localParams);
                itsNewSpectra(ii)= getMTSpectrumStats(thisClip, specParams);
                progressbar(ii/numel(regions));
            end
            model.itsSpectra = [model.itsSpectra itsNewSpectra];
            
        end
        
        function importParams(model, params)
            model.localParams = params;
        end
        %%%%%%%%%%%%%%%% accessors %%%%%%%%%%%%%%%%%%%%%%%%
        function set.threshold(model,t)
            if t < 0
               error('Learner:BadThreshold', ...
                   'Threshold must be a non-negative number');
            end
            if size(t) > 1
                error('Learner:BadThreshold', ...
                   'Threshold must be a scalar');
            end
            model.threshold_ = t;
        end
        
        function t = get.threshold(model)
            t = model.threshold_;
        end     
        
        function set.fastThreshold(model,t)
            if t < 0
               error('Learner:BadThreshold', ...
                   'Threshold must be a non-negative number');
            end
            if size(t) > 1
                error('Learner:BadThreshold', ...
                   'Threshold must be a scalar');
            end
            model.fastThreshold_ = t;
        end
        
        function t = get.fastThreshold(model)
            t = model.fastThreshold_;
        end 
        % length of model in seconds
        function n = get.typicalLength(model)
            n = zeros(model.nExamples,1);
            for ii = 1:model.nExamples
                n(ii) = model.itsSpectra(ii).times(end) - model.itsSpectra(ii).times(1);
            end
            n = mean(n);
        end
         
        function u = get.nExamples(model)
            u = numel(model.itsSpectra);
        end
    end
end
    