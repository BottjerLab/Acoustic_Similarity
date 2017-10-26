function [ageSylls, sessions, ageSpectra, ageFeatures] = loadAgeSylls(birdID, age)
    % get age and session ID for each syllable of a given bird/age
dataDir  = ['data' filesep birdID filesep];

load([dataDir 'allSpecs-' birdID '.mat']);
if strcmp(birdID,'Db113') %JMA added to fix mistake on ages from before clustering apparently
    for a = 1: length(DRsylls)
        b = DRsylls(a).age - 3;
        DRsylls(a).age = b;
    end
end
if strcmp(birdID,'Dg138') %JMA added to fix mistake on ages from before clustering apparently
    for a = 1: length(DRsylls)
        b = DRsylls(a).age - 4;
        DRsylls(a).age = b;
    end
end
if strcmp(birdID,'R247') %JMA added to fix mistake on ages from before clustering apparently
    for a = 1: length(DRsylls)
        b = DRsylls(a).age - 3;
        DRsylls(a).age = b;
    end
end
if strcmp(birdID,'R288') %JMA added to fix mistake on ages from before clustering apparently
    for a = 1: length(DRsylls)
        b = DRsylls(a).age - 2;
        DRsylls(a).age = b;
    end
end
uAges = unique([DRsylls.age]);

if ~any(age == uAges), error('Age %d not represented', age); end
seldAge = [DRsylls.age]==age;

ageSylls = DRsylls(seldAge);
% get the session ID for each syllable
if nargout >= 2
    sessions = cell(1,numel(ageSylls));
    for ii = 1:numel(ageSylls)
        [~,sessions{ii}] = fileparts(ageSylls(ii).file);
    end
end

if nargout >= 3
    ageSpectra = spectra(seldAge);
end

if nargout >= 4
    ageFeatures = featureTable(seldAge);
    %todo: censor features
end